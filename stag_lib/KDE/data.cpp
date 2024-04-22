/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/
#include <iostream>
#include <fstream>
#include <random>
#include <utility>
#include <deque>
#include <future>

#include "data.h"
#include "utility.h"
#include "Graph/random.h"
#include "KDE/kde.h"
#include "multithreading/ctpl_stl.h"

#define KDE_TREE_CUTOFF 1000

/*
 * Used to disable compiler warning for unused variable.
 */
template<class T> void ignore_warning(const T&){}

//------------------------------------------------------------------------------
// Implementation of the DataPoint class.
//------------------------------------------------------------------------------
stag::DataPoint::DataPoint(DenseMat& all_data, StagInt row_index) {
  dimension = all_data.cols();
  coordinates = all_data.row(row_index).data();
}

stag::DataPoint::DataPoint(std::vector<StagReal>& point_vector) {
  dimension = point_vector.size();
  coordinates = &point_vector[0];
}

//------------------------------------------------------------------------------
// Converting between data formats.
//------------------------------------------------------------------------------
std::vector<stag::DataPoint> stag::matrix_to_datapoints(DenseMat* data) {
  std::vector<stag::DataPoint> res;

  StagInt n = data->rows();
  StagInt d = data->cols();
  for (auto i = 0; i < n; i++) {
    res.emplace_back(d, data->row(i).data());
  }
  return res;
}

//------------------------------------------------------------------------------
// Loading and saving matrices to file.
//------------------------------------------------------------------------------
void stag::save_matrix(DenseMat& data, std::string& filename) {
  // Attempt to open the specified file
  std::ofstream os(filename);

  // If the file could not be opened, throw an exception
  if (!os.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Iterate through the entries in the matrix, and write the data file.
  for (auto i = 0; i < data.rows(); i++) {
    for (auto j = 0; j < data.cols(); j++) {
      os << data.coeffRef(i, j) << " ";
    }
    os << std::endl;
  }

  // Close the output file stream
  os.close();
}

/**
 * Parse a single content line of a 'csv'-style file containing numeric values.
 * This method assumes that the line is not a comment, and tries to parse either
 * by splitting on commas or whitespace.
 *
 * @param line the string to parse
 * @return a vector containing the values
 * @throw std::invalid_argument the line cannot be parsed
 */
std::vector<StagReal> parse_numeric_csv_content_line(std::string line) {
  // List the possible delimiters for the elements on the line
  std::vector<std::string> delimiters{",", " ", "\t"};

  // Split the line to extract the information
  StagUInt num_tokens_found = 0;
  std::vector<StagReal> tokens;

  // Try splitting by each delimiter in turn
  size_t split_pos = 0;
  std::string token;
  for (const std::string &delimiter: delimiters) {
    // If we have not found any delimiters of the previous types yet,
    // then try this one.
    if (num_tokens_found == 0) {
      while ((split_pos = line.find(delimiter)) != std::string::npos) {
        // Extract the portion of the line up to the delimiter
        token = line.substr(0, split_pos);
        line.erase(0, split_pos + delimiter.length());

        // If the token has length 0, then skip
        if (token.length() == 0) continue;

        // Parse the token as a real
        // throws an exception if the token cannot be parsed.
        tokens.push_back(std::stod(token));

        // Increase the counter of the number of tokens found
        num_tokens_found++;
      }
    }
  }

  // Extract the final token in the line
  try {
    // Try extracting another token from the rest of the line, but ignore any
    // errors - there might be none.
    StagReal last_token = std::stod(line);
    tokens.push_back(last_token);
    num_tokens_found++;
  } catch (std::exception &e) {
    // Ignore this
  }

  // Return the data
  return tokens;
}

DenseMat stag::load_matrix(std::string& filename) {
  // Attempt to open the provided file
  std::ifstream is(filename);

  // If the file could not be opened, throw an exception
  if (!is.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Initialise the data matrix to be 1 x 1.
  // Once we know the number of columns, we will update num_cols
  StagInt num_cols = -1;
  StagInt num_rows = 0;
  DenseMat data(1, 1);

  // Read in the file one line at a time
  std::string line;
  while (stag::safeGetline(is, line)) {
    if (line.length() > 0 && line[0] != '#' && line[0] != '/') {
      // Get all of the data from this line
      std::vector<StagReal> tokens = parse_numeric_csv_content_line(line);

      // If we don't yet know the number of columns in the data matrix, set it
      // now.
      if (num_cols < 0) {
        num_cols = (StagInt) tokens.size();
        data.conservativeResize(Eigen::NoChange_t::NoChange, num_cols);
      }

      // If the number of columns on this line does not match previous lines,
      // then throw an exception.
      if (num_cols != (StagInt) tokens.size()) {
        throw std::runtime_error("Rows of different size in matrix file.");
      }

      num_rows++;

      if (num_rows <= data.rows()) {
        for (auto i = 0; i < num_cols; i++) {
          data.coeffRef(num_rows-1, i) = tokens[i];
        }
      } else {
        // We need to resize the data matrix
        data.conservativeResize(ceil(1.1 * (StagReal) num_rows),
                                Eigen::NoChange_t::NoChange);

        for (auto i = 0; i < num_cols; i++) {
          data.coeffRef(num_rows-1, i) = tokens[i];
        }
      }
    }
  }

  // Finally, we resize the matrix to its final size.
  data.conservativeResize(num_rows, Eigen::NoChange_t::NoChange);

  return data;
}

//------------------------------------------------------------------------------
// Constructing approximate similarity graph.
//------------------------------------------------------------------------------
class KDETreeEntry {
public:
  KDETreeEntry(DenseMat* data, StagReal a, StagInt min_id, StagInt max_id)
    : sampling_dist(0.0, 1.0)
  {
     min_idx = min_id;
     max_idx = max_id;
    LOG_DEBUG("Initialising from " << min_idx << " to " << max_idx << std::endl);

     if (max_idx - min_idx >= KDE_TREE_CUTOFF) {
       below_cutoff = false;
       initialise_estimator(data, a);
     } else {
       below_cutoff = true;
       initialise_exact(data, a);
     }

     if (min_idx != max_idx) {
       StagInt midpoint = (min_idx + max_idx) / 2;
       assert(midpoint >= min_idx);
       assert(midpoint < max_idx);
       left_child = new KDETreeEntry(data, a, min_idx, midpoint);
       right_child = new KDETreeEntry(data, a, midpoint + 1, max_idx);
     }
  }

  void initialise_estimator(DenseMat* data, StagReal a) {
    StagInt n = data->rows();
    StagInt K1 = ceil(log((StagReal) n) / 0.25);
    StagReal K2_constant = 5 * log((StagReal) n);
    StagReal min_mu = 1.0 / (StagReal) n;
    this_estimator = stag::CKNSGaussianKDE(data, a, min_mu, K1, K2_constant,
                                           0, min_idx, max_idx + 1);
  }

  void initialise_exact(DenseMat* data, StagReal a) {
    exact_kde = stag::ExactGaussianKDE(data, a, min_idx, max_idx + 1);
  }

  StagReal estimate_weight(const stag::DataPoint& q, const StagInt q_id) {
    StagReal weight;
    cache_mutex.lock();
    if (!cached_weights.contains(q_id)) {
      cache_mutex.unlock();
      if (!below_cutoff) {
        weight = this_estimator.query(q);
      } else {
        weight = exact_kde.query(q);
      }

      cache_mutex.lock();
      cached_weights[q_id] = weight;
      cache_mutex.unlock();

    } else {
      weight = cached_weights[q_id];
      cache_mutex.unlock();
    }
    return weight;
  }

  std::vector<StagReal> estimate_weights(DenseMat* q) {
    if (!below_cutoff) {
      return this_estimator.query(q);
    } else {
      return exact_kde.query(q);
    }
  }

  StagInt sample_neighbor(const stag::DataPoint& q, const StagInt q_id) {
    if (min_idx == max_idx) {
      return min_idx;
    } else {
      StagReal left_est = left_child->estimate_weight(q, q_id);
      StagReal right_est = right_child->estimate_weight(q, q_id);
      StagReal my_est = left_est + right_est;

      if (sampling_dist(*stag::get_global_rng()) <= left_est / my_est) {
        return left_child->sample_neighbor(q, q_id);
      } else {
        return right_child->sample_neighbor(q, q_id);
      }
    }
  }

private:
  bool below_cutoff;
  stag::CKNSGaussianKDE this_estimator;
  stag::ExactGaussianKDE exact_kde;
  StagInt min_idx;
  StagInt max_idx;
  KDETreeEntry* left_child;
  KDETreeEntry* right_child;
  std::uniform_real_distribution<double> sampling_dist;
  std::unordered_map<StagInt, StagReal> cached_weights;
  std::mutex cache_mutex;
};

void sample_asg_edges(DenseMat* data,
                      KDETreeEntry& tree_root,
                      StagInt chunk_start,
                      StagInt chunk_end,
                      std::vector<EdgeTriplet>& edges,
                      StagInt edges_per_node) {
  for (auto i = chunk_start; i < chunk_end; i++) {
    stag::DataPoint q(*data, i);

    for (auto j = 0; j < edges_per_node; j++) {
      StagInt neighbor = tree_root.sample_neighbor(q, i);
      StagInt this_base_index = 2 * ((i * j) + j);
      edges.at(this_base_index) = EdgeTriplet(i, neighbor, 1);
      edges.at(this_base_index + 1) = EdgeTriplet(neighbor, i, 1);
    }
  }
}

stag::Graph stag::approximate_similarity_graph(DenseMat* data, StagReal a) {
  // Begin by creating the tree of kernel density estimators.
  // Creating each node is parallelized by the KDE code.
  KDETreeEntry tree_root(data, a, 0, data->rows() - 1);

  // Work out how many total edges we will sample
  auto edges_per_node = (StagInt) (3 * log((StagReal) data->rows()));
  StagInt num_edges = data->rows() * edges_per_node;

  // First, compute the degrees of all nodes
  std::vector<StagReal> degrees = tree_root.estimate_weights(data);

  StagInt num_threads = std::thread::hardware_concurrency();
  std::vector<EdgeTriplet> graph_edges(2 * num_edges);
  if (data->rows() <= num_threads * 2) {
    sample_asg_edges(data, tree_root, 0, data->rows(), graph_edges, edges_per_node);
  } else {
    // Start the thread pool
    ctpl::thread_pool pool((int) num_threads);
    std::vector<std::future<void>> futures;

    StagInt chunk_size = floor((StagReal) data->rows() / (StagReal) num_threads);

    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      StagInt this_chunk_start = chunk_id * chunk_size;
      StagInt this_chunk_end = this_chunk_start + chunk_size;
      if (chunk_id == num_threads - 1) this_chunk_end = data->rows();

      assert(this_chunk_start <= (StagInt) data->rows());
      assert(this_chunk_end <= (StagInt) data->rows());
      assert(this_chunk_end > this_chunk_start);

      futures.push_back(
        pool.push(
          [&, this_chunk_start, this_chunk_end, edges_per_node] (int id) {
            ignore_warning(id);
            sample_asg_edges(data, tree_root, this_chunk_start,
                             this_chunk_end, graph_edges, edges_per_node);
          }
        )
      );
    }

    // Join the futures
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      futures[chunk_id].get();
    }
  }

  // Return a graph
  SprsMat adj_mat(data->rows(), data->rows());
  adj_mat.setFromTriplets(graph_edges.begin(), graph_edges.end());
  return stag::Graph(adj_mat);
}
