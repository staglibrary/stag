/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/
#include <iostream>

#include "random.h"
#include "graph.h"

std::random_device dev_g;
std::mt19937_64 rng_g(dev_g());

std::mt19937_64* stag::get_global_rng() {
  return &rng_g;
}

std::mt19937_64 stag::create_rng() {
  std::random_device local_dev;
  std::mt19937_64 local_rng(local_dev());
  return local_rng;
}

/**
 * Get an upper estimate on the number of neighbours of each node in a
 * graph generated from the stochastic block model.
 *
 * This will be used to pre-allocate the memory of the adjacency matrix as we
 * construct the graph.
 *
 * @param cluster_sizes a vector of length \f$k\f$ with the number of vertices
 *                      in each cluster.
 * @param probabilities a \f$k \times k\f$ matrix with the inter-cluster
 *                      probabilities.
 * @return an Eigen vector with the neighbour estimates.
 */
Eigen::VectorXi estimate_sbm_neighbours(std::vector<StagInt>& cluster_sizes,
                                        DenseMat probabilities) {
  // Get the total number of vertices in the graph
  auto k = (StagInt) cluster_sizes.size();
  StagInt n = 0;
  for (StagInt s : cluster_sizes) n += s;

  // Create the vector which we'll return
  Eigen::VectorXi neighbours(n);

  // Compute the expected number of neighbours for a node in each cluster
  // This is basically the matrix product of probabilities with cluster_sizes
  Eigen::VectorXi cluster_neighbours(k);
  for (auto i = 0; i < k; i++) {
    StagInt this_cluster_neighbours = 0;
    for (auto j = 0; j < k; j++) {
      this_cluster_neighbours += cluster_sizes.at(j) * probabilities(i, j);
    }

    // Add a safety factor of 2
    cluster_neighbours.coeffRef(i) = 2 * this_cluster_neighbours;
  }

  // Now add the neighbour estimates for each actual node
  StagInt current_idx = 0;
  for (StagInt i = 0; i < k; i++) {
    for (StagInt j = 0; j < cluster_sizes.at(i); j++) {
      neighbours.coeffRef(current_idx) = cluster_neighbours(i);
      current_idx++;
    }
  }

  return neighbours;
}

/**
 * Sample edges between SBM clusters by directly iterating through each
 * edge and 'tossing a coin'. This technique should be used for 'large' values
 * of p.
 *
 * @param adj_mat the adjacency matrix to be updated - can pass null to not update anything
 * @param edgelist_os the edgelist file stream to be updated - can pass null to not update anything
 * @param cluster_idx the index of the 'source' cluster
 * @param other_cluster_idx the index of the 'target' cluster
 * @param this_cluster_vertices the number of vertices in the 'source' cluster
 * @param other_cluster_vertices the number of vertices in the 'target' cluster
 * @param this_cluster_start_idx the index of the first vertex in the 'source' cluster
 * @param other_cluster_start_idx the index of the first vertex in the 'target' cluster
 * @param p the probability of including each edge.
 */
void sample_edges_directly(SprsMat* adj_mat,
                           std::ostream* edgelist_os,
                           StagInt cluster_idx,
                           StagInt other_cluster_idx,
                           StagInt this_cluster_vertices,
                           StagInt other_cluster_vertices,
                           StagInt this_cluster_start_idx,
                           StagInt other_cluster_start_idx,
                           StagReal p) {
  // Validate the function inputs
  assert(0 <= p && p <= 1);

  // Prepare the random number generator
  std::bernoulli_distribution sampleDist(p);

  for (StagInt i = this_cluster_start_idx;
      i < this_cluster_start_idx + this_cluster_vertices; i++) {
    for (StagInt j = other_cluster_start_idx;
        j < other_cluster_start_idx + other_cluster_vertices; j++) {
      // If we are in the same cluster, then don't double sample
      if (cluster_idx == other_cluster_idx && j <= i) continue;

      // Toss a coin
      if (sampleDist(*stag::get_global_rng())) {
        if (adj_mat != nullptr) {
          adj_mat->insert(i, j) = 1;
          adj_mat->insert(j, i) = 1;
        }
        if (edgelist_os != nullptr) {
          *edgelist_os << i << " " << j << " " << 1 << std::endl;
        }
      }
    }
  }
}

/**
 * Sample edges between SBM clusters using the 'binomial trick'. This technique
 * should be used for 'small' values of p.
 *
 * @param adj_mat the adjacency matrix to be updated - can pass null to not update anything
 * @param edgelist_os the edgelist file stream to be updated - can pass null to not update anything
 * @param this_cluster_vertices the number of vertices in the 'source' cluster
 * @param other_cluster_vertices the number of vertices in the 'target' cluster
 * @param this_cluster_start_idx the index of the first vertex in the 'source' cluster
 * @param other_cluster_start_idx the index of the first vertex in the 'target' cluster
 * @param p
 */
void sample_edges_binomial(SprsMat* adj_mat,
                           std::ostream* edgelist_os,
                           StagInt this_cluster_vertices,
                           StagInt other_cluster_vertices,
                           StagInt this_cluster_start_idx,
                           StagInt other_cluster_start_idx,
                           StagReal p) {
  // Validate the function inputs
  assert(0 <= p && p <= 1);

  // Get the total number of possible edges and the expected number of edges
  StagInt max_edges = this_cluster_vertices * other_cluster_vertices;
  if (this_cluster_start_idx == other_cluster_start_idx) max_edges /= 2;
  StagReal expected_num_edges = p * ((StagReal) max_edges);

  // If the exppected number of edges is 0, don't do anything
  if (expected_num_edges < 1) {
    return;
  }

  // Prepare the random number generator. We will approximate the binomial
  // distribution with the normal distribution
  assert(sqrt(1 - p) * expected_num_edges > 0);
  std::normal_distribution<StagReal> numEdgesDist(expected_num_edges,
                                                  sqrt((1 - p) * expected_num_edges));
  std::uniform_int_distribution<StagInt> thisVertexDist(0, this_cluster_vertices - 1);
  std::uniform_int_distribution<StagInt> otherVertexDist(0, other_cluster_vertices - 1);

  // Decide how many edges to sample based on the 'binomial' distribution
  auto raw_sample = (StagInt) floor(numEdgesDist(*stag::get_global_rng()));
  StagInt numEdges = std::max((StagInt) 0, std::min(max_edges, raw_sample));

  // Sample the specific vertices
  StagInt randU = 0;
  StagInt randV = 0;
  for (StagInt i = 0; i < numEdges; i++) {
    // Choose two random vertices in the cluster
    randU = 0;
    randV = 0;
    while (randU == randV) {
      // Ignore the edge if u and v are identical
      randU = this_cluster_start_idx + thisVertexDist(*stag::get_global_rng());
      randV = other_cluster_start_idx + otherVertexDist(*stag::get_global_rng());
    }

    // Add this vertex to the adjacency matrix
    if (adj_mat != nullptr) {
      adj_mat->coeffRef(randU, randV)++;
      adj_mat->coeffRef(randV, randU)++;
    }
    if (edgelist_os != nullptr) {
      *edgelist_os << randU << " " << randV << " " << 1 << std::endl;
    }
  }
}

stag::Graph stag::sbm(StagInt n, StagInt k, StagReal p, StagReal q) {
  return stag::sbm(n, k, p, q, false);
}

stag::Graph stag::sbm(StagInt n, StagInt k, StagReal p, StagReal q, bool exact) {
  if (n < 1) throw std::invalid_argument("Number of vertices must be at least 1.");
  if (k < 1 || k > n/2) {
    throw std::invalid_argument("Number of clusters must be between 1 and n/2.");
  }
  if (p < 0 || p > 1) {
    throw std::invalid_argument("p must be between 0 and 1.");
  }
  if (q < 0 || q > 1) {
    throw std::invalid_argument("q must be between 0 and 1.");
  }

  // Create the cluster size vector and probabilities matrix
  std::vector<StagInt> cluster_sizes;
  DenseMat probabilities(k, k);
  for (auto i = 0; i < k; i++) {
    cluster_sizes.push_back(floor(((StagReal) n) / ((StagReal) k)));
    probabilities(i, i) = p;

    for (auto j = i + 1; j < k; j++) {
      probabilities(i, j) = q;
      probabilities(j, i) = q;
    }
  }

  return general_sbm(cluster_sizes, probabilities, exact);
}

void general_sbm_internal(SprsMat* adj_mat,
                          std::ostream* edgelist_os,
                          std::vector<StagInt>& cluster_sizes,
                          DenseMat& probabilities, bool exact) {
  // The number of clusters is the length of the cluster_sizes vector
  StagInt k = cluster_sizes.size();

  // Check that the input parameters make sense.
  for (StagInt size : cluster_sizes) {
    if (size < 1) throw std::invalid_argument("Number of vertices in each cluster must be at least 1.");
  }
  if (probabilities.rows() != k || probabilities.cols() != k) {
    throw std::invalid_argument("Probability matrix must be of size k * k.");
  }
  for (auto i = 0; i < k; i++) {
    for (auto j = 0; j < k; j++) {
      if (probabilities(i, j) < 0 || probabilities(i, j) > 1) {
        throw std::invalid_argument("All probabilities must be between 0 and 1.");
      }
    }
  }

  // Estimate the number of neighbours of each node and reserve
  // memory for the adjacency matrix.
  if (adj_mat != nullptr) {
    Eigen::VectorXi neighbour_estimates = estimate_sbm_neighbours(cluster_sizes,
                                                                  probabilities);
    adj_mat->reserve(neighbour_estimates);
  }

  // Iterate through the clusters
  StagInt this_cluster_start_idx = 0;
  for (StagInt cluster_idx = 0; cluster_idx < k; cluster_idx++) {
    // Get the number of vertices in the current cluster
    StagInt this_cluster_vertices = cluster_sizes.at(cluster_idx);

    StagInt other_cluster_start_idx = this_cluster_start_idx;
    for (StagInt other_cluster_idx = cluster_idx;
         other_cluster_idx < k; other_cluster_idx++){
      // Get the number of vertices in the other cluster
      StagInt other_cluster_vertices = cluster_sizes.at(other_cluster_idx);

      // Get the sampling probability between the two clusters
      StagReal prob = probabilities(cluster_idx, other_cluster_idx);

      // Sample the edges between this cluster and the other cluster
      if (this_cluster_vertices * other_cluster_vertices >= 10000 &&
            prob < 0.5 && !exact) {
        // For small probabilities, use the 'binomial trick' for sampling
        sample_edges_binomial(adj_mat,
                              edgelist_os,
                              this_cluster_vertices,
                              other_cluster_vertices,
                              this_cluster_start_idx,
                              other_cluster_start_idx,
                              prob);
      } else {
        // For large probabilities, we just iterate over every pair of vertices
        sample_edges_directly(adj_mat,
                              edgelist_os,
                              cluster_idx,
                              other_cluster_idx,
                              this_cluster_vertices,
                              other_cluster_vertices,
                              this_cluster_start_idx,
                              other_cluster_start_idx,
                              prob);
      }

      // Update the starting index for the other cluster
      other_cluster_start_idx += other_cluster_vertices;
    }

    // Update the starting index for this cluster
    this_cluster_start_idx += this_cluster_vertices;
  }

  // Finally, construct and return the graph.
  if (adj_mat != nullptr) adj_mat->makeCompressed();
}

stag::Graph stag::general_sbm(std::vector<StagInt> &cluster_sizes,
                              DenseMat &probabilities, bool exact) {
  StagInt n = 0;
  for (auto s : cluster_sizes) n += s;

  // Initialise a sparse adjacency matrix, and use a null pointer for the
  // edgelist file. This ensures that the general_sbm_internal method
  // writes the generated graph to the adjacency matrix and does not write
  // the graph to disk.
  SprsMat adj_mat(n, n);
  std::ostream* edgelist_os = nullptr;

  general_sbm_internal(&adj_mat, edgelist_os, cluster_sizes,
                       probabilities, exact);
  return stag::Graph(adj_mat);
}

stag::Graph stag::general_sbm(std::vector<StagInt>& cluster_sizes,
                              DenseMat& probabilities) {
  return stag::general_sbm(cluster_sizes, probabilities, false);
}

stag::Graph stag::erdos_renyi(StagInt n, StagReal p) {
  return stag::erdos_renyi(n, p, false);
}

stag::Graph stag::erdos_renyi(StagInt n, StagReal p, bool exact) {
  return stag::sbm(n, 1, p, 0, exact);
}

void stag::general_sbm_edgelist(std::string &filename,
                                std::vector<StagInt> &cluster_sizes,
                                DenseMat &probabilities,
                                bool exact) {
  SprsMat* adj_mat = nullptr;
  std::ofstream os(filename);

  // If the file could not be opened, throw an exception
  if (!os.is_open()) {
    throw std::runtime_error(std::strerror(errno));
  }

  // Write a header to the output stream giving the parameters of the
  // SBM model.
  StagInt k = cluster_sizes.size();
  StagInt n = 0;
  for (auto size: cluster_sizes) n += size;
  os << "# This graph was generated from a stochastic block model with the ";
  os << "following parameters." << std::endl;
  os << "#    n = " << n << std::endl;
  os << "#    k = " << k << std::endl;

  if (k <= 20) {
    os << "#    cluster sizes = ";
    for (StagInt size : cluster_sizes) os << size << " ";
    os << std::endl;
    os << "#    probability matrix = " << std::endl;
    for (auto i = 0; i < k; i++) {
      os << "#        ";
      for (auto j = 0; j < k; j++) {
        os << probabilities.coeffRef(i, j) << " ";
      }
      os << std::endl;
    }
  } else {
    os << "# (Probability matrix omitted as it is too large.)" << std::endl;
  }

  // Generate the graph.
  general_sbm_internal(adj_mat, &os, cluster_sizes, probabilities, exact);

  // Close the file output stream.
  os.close();
}

void stag::general_sbm_edgelist(std::string &filename,
                                std::vector<StagInt> &cluster_sizes,
                                DenseMat &probabilities) {
  stag::general_sbm_edgelist(filename, cluster_sizes, probabilities, false);
}

std::vector<StagInt> stag::sbm_gt_labels(StagInt n, StagInt k) {
  if (n < 1) throw std::invalid_argument("Number of vertices must be at least 1.");
  if (k < 1 || k > n/2) {
    throw std::invalid_argument("Number of clusters must be between 1 and n/2.");
  }

  // Create the cluster size vector
  std::vector<StagInt> cluster_sizes;
  for (auto i = 0; i < k; i++) {
    cluster_sizes.push_back(floor(((StagReal) n) / ((StagReal) k)));
  }

  return general_sbm_gt_labels(cluster_sizes);
}

std::vector<StagInt> stag::general_sbm_gt_labels(std::vector<StagInt>& cluster_sizes) {
  for (StagInt size : cluster_sizes) {
    if (size < 1) throw std::invalid_argument("Number of vertices in each cluster must be at least 1.");
  }

  std::vector<StagInt> labels;

  StagInt current_cluster = 0;
  for (auto this_size : cluster_sizes) {
    for (auto j = 0; j < this_size; j++) {
      labels.push_back(current_cluster);
    }
    current_cluster++;
  }

  return labels;
}