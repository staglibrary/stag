/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/
// Standard C++ libraries
#include <vector>
#include <set>
#include <deque>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <cmath>

// Additional libraries
#include <Eigen/Sparse>
#include "KMeansRex/KMeansRexCoreInterface.h"
#include "indicators/indicators.hpp"
#include "multithreading/ctpl_stl.h"

// STAG modules
#include "random.h"
#include "graph.h"
#include "cluster.h"
#include "utility.h"
#include "spectrum.h"
#include "data.h"
#include "kde.h"

#define ASG_TREE_CUTOFF 5000

/*
 * Used to disable compiler warning for unused variable.
 */
template<class T> void ignore_warning(const T&){}


std::vector<StagInt> stag::spectral_cluster(stag::Graph *graph, StagInt k) {
  // Check that the number of clusters is valid.
  if (k < 1 || k > graph->number_of_vertices() /2) {
    throw std::invalid_argument("Number of clusters must be between 1 and n/2.");
  }

  // Start by computing the 'first' k eigenvalues of the normalised graph
  // laplacian matrix.
  Eigen::MatrixXd eigvecs = stag::compute_eigenvectors(
      graph, stag::GraphMatrix::NormalisedLaplacian, k, stag::EigenSortRule::Smallest);

  // Run k-means clustering on the spectral embedding of the vertices
  Eigen::MatrixXd centres = Eigen::MatrixXd::Zero(k, k);
  Eigen::VectorXd clusters = Eigen::VectorXd::Zero(eigvecs.rows());
  char initialisation[9] = "plusplus";
  RunKMeans(eigvecs.data(),
            eigvecs.rows(),
            k,
            k,
            k * 100,
            42,
            initialisation,
            centres.data(),
            clusters.data()
            );

  // Move the eigen data into the vector to return
  return {clusters.data(), clusters.data() + clusters.rows()};
}

std::vector<StagInt> stag::cheeger_cut(stag::Graph* graph) {
  // First, compute the first 2 eigenvectors of the normalised graph Laplacian
  // matrix.
  stag::EigenSystem eigsys = stag::compute_eigensystem(
      graph, stag::GraphMatrix::NormalisedLaplacian, 2, stag::EigenSortRule::Smallest);

  // The vector to pass to the sweep set method is the second eigenvector,
  // with each entry normalised by the square root of the node degree.
  SprsMat vec;

  // Check the eigenvalues to ensure we are selecting the second eigenvector
  if (get<0>(eigsys).coeff(0) > get<0>(eigsys).coeff(1)) {
    vec = get<1>(eigsys).col(0).sparseView();
  } else {
    vec = get<1>(eigsys).col(1).sparseView();
  }

  // Normalise by the square root of vertex degrees.
  for (StagInt i = 0; i < graph->number_of_vertices(); i++) {
    vec.coeffRef(i, 0) = std::sqrt(1 / graph->degree(i)) * vec.coeff(i, 0);
  }

  // Perform the sweep
  std::vector<StagInt> clusters (graph->number_of_vertices(), 0);
  std::vector<StagInt> cut = stag::sweep_set_conductance(graph, vec);
  for (StagInt i : cut) clusters.at(i) = 1;
  return clusters;
}

std::vector<StagInt> stag::local_cluster(stag::LocalGraph *graph, StagInt seed_vertex, double target_volume) {
  if (target_volume <= 0) throw std::invalid_argument("Target volume must be positive.");

  // The 'locality' parameter should essentially be a constant if we expect the conductance
  // of the target cluster to be a constant. Set it to 0.01.
  // The error parameter should decrease as the inverse of the target volume.
  return stag::local_cluster_acl(graph,
                                 seed_vertex,
                                 0.01,
                                 1./ target_volume);
}

//------------------------------------------------------------------------------
// Implementation of ACL Local Clustering Algorithm based on PageRank vectors
//------------------------------------------------------------------------------
std::vector<StagInt> stag::local_cluster_acl(stag::LocalGraph* graph,
                                             StagInt seed_vertex,
                                             double locality) {
  return stag::local_cluster_acl(graph, seed_vertex, locality, 0.001);
}

std::vector<StagInt> stag::local_cluster_acl(stag::LocalGraph *graph,
                                             StagInt seed_vertex,
                                             double locality,
                                             double error) {
  // Check that the arguments are valid.
  if (!graph->vertex_exists(seed_vertex)) {
    throw std::invalid_argument("Seed vertex does not exist.");
  }
  if (locality < 0 || locality > 1) {
    throw std::invalid_argument("Locality parameter must be between 0 and 1.");
  }
  if (error <= 0) {
    throw std::invalid_argument("Error parameter must be greater than 0.");
  }

  // Compute the approximate pagerank vector
  SprsMat seedDist(seed_vertex + 1, 1);
  seedDist.coeffRef(seed_vertex, 0) = 1;
  std::tuple<SprsMat, SprsMat> apr = stag::approximate_pagerank(graph, seedDist, locality, error);
  SprsMat p = std::get<0>(apr);

  // Normalise the values in p by their degree. Make just one 'call' to the
  // LocalGraph object to get the degrees for all relevant vertices.
  std::vector<double> degrees = graph->degrees(stag::sprsMatInnerIndices(&p));
  StagInt degree_index = 0;
  for (SprsMat::InnerIterator it(p, 0); it; ++it) {
    it.valueRef() = it.value() / degrees.at(degree_index);
    degree_index++;
  }

  // Perform the sweep set operation on the approximate pagerank p to find the
  // sweep set with the optimal conductance.
  return stag::sweep_set_conductance(graph, p);
}

//------------------------------------------------------------------------------
// Pagerank and Approximate Pagerank Implementations
//------------------------------------------------------------------------------

/**
 * Perform the push operation on vertex u in the given graph.
 *
 * Updates the sparse vectors p and r, according to:
 *   p'(u) = p(u) + alpha * r(u)
 *   r'(v) = r(v) + (1 - alpha) * w(u, v) * r(u) / (2 * deg(u))
 *   r'(u) = (1 - alpha) * r(u) / 2
 * for all v which are neighbors of u.
 *
 * @param graph
 * @param p
 * @param r
 * @param alpha
 * @param u
 */
void push(stag::LocalGraph *graph, SprsMat *p, SprsMat *r, double alpha, StagInt u) {
  // First, make sure that the vector p is large enough.
  if (u >= p->rows()) {
    p->conservativeResize(u + 1, 1);
  }

  // Update p according to the push operation
  p->coeffRef(u, 0) = p->coeff(u, 0) + alpha * r->coeff(u, 0);

  // Iterate through the neighbors of u
  double deg = graph->degree(u);
  StagInt v;
  for (stag::edge e : graph->neighbors(u)) {
    v = e.v2;
    assert(v != u);

    // Make sure that r is large enough.
    if (v >= r->rows()) {
      r->conservativeResize(v + 1, 1);
    }

    // Perform the push operation on r for the vertex v
    r->coeffRef(v, 0) = r->coeff(v, 0) + (1 - alpha) * e.weight * r->coeff(u, 0) / (deg * 2);
  }

  // Update r(u) according to the push operation
  // This must happen after updating the neighbors of r
  r->coeffRef(u, 0) = (1 - alpha) * r->coeff(u, 0) / 2;
}

std::tuple<SprsMat, SprsMat> stag::approximate_pagerank(stag::LocalGraph *graph,
                                                        SprsMat &seed_vector,
                                                        double alpha,
                                                        double epsilon) {
  // Check that the arguments are valid
  if (seed_vector.cols() > 1) throw std::invalid_argument("Seed vector must be a column vector.");
  if (!graph->vertex_exists(seed_vector.rows() - 1)) {
    throw std::invalid_argument("Seed vector dimension must be less than the number of vertices in the graph");
  }
  if (alpha < 0 || alpha > 1) {
    throw std::invalid_argument("Alpha parameter must be between 0 and 1.");
  }
  if (epsilon <= 0) {
    throw std::invalid_argument("Epsilon parameter must be greater than 0.");
  }

  // Initialise p to be the all-zeros vector.
  SprsMat p(seed_vector.rows(), 1);

  // Initialise r to be equal to the seed vector
  SprsMat r(seed_vector);

  // Compress the seed vector
  seed_vector.makeCompressed();

  // We will maintain a queue of vertices satisfying
  //    r(u) >= epsilon * deg(u)
  // Along with the queue, maintain an unordered set which tracks the vertices
  // already in the queue.
  StagInt u;
  std::vector<double> degrees = graph->degrees(
      stag::sprsMatInnerIndices(&seed_vector));
  std::deque<StagInt> vertex_queue;
  std::unordered_set<StagInt> queue_members;
  StagInt degree_index = 0;
  for (SprsMat::InnerIterator it(seed_vector, 0); it; ++it) {
    u = it.row();
    if (r.coeff(u, 0) >= epsilon * degrees.at(degree_index)
          && degrees.at(degree_index) != 0) {
      vertex_queue.push_back(u);
      queue_members.insert(u);
    }
    degree_index++;
  }

  // While the queue is not empty, push an entry off the queue and apply the
  // push operation.
  while (!vertex_queue.empty()) {
    // Get the next vertex from the queue
    u = vertex_queue.front();
    vertex_queue.pop_front();
    queue_members.erase(u);

    // Perform the push operation on this vertex
    push(graph, &p, &r, alpha, u);

    // Check u to see if it should be added back to the queue
    // If so, then we add it back to the start of the queue in order to be as
    // cache-efficient as possible when accessing a graph from disk.
    if (r.coeff(u, 0) >= epsilon * graph->degree(u)) {
      vertex_queue.push_front(u);
      queue_members.insert(u);
    }

    // Check the neighbors of u to see if they should be added back to the queue
    // Skip any neighbors which are already in the queue.
    std::vector<StagInt> neighbors = graph->neighbors_unweighted(u);
    std::vector<double> neighbor_degrees = graph->degrees(neighbors);

    // The length of neighbors and neighbor_degrees should always be equal.
    // If they are not, there must be a bug in the implementation of graph->degrees.
    assert(neighbors.size() == neighbor_degrees.size());

    degree_index = 0;
    StagInt v;
    for (stag::edge e : graph->neighbors(u)) {
      v = e.v2;
      if (r.coeff(e.v2, 0) >= epsilon * neighbor_degrees.at(degree_index) &&
            !queue_members.contains(v)) {
        vertex_queue.push_back(v);
        queue_members.insert(v);
      }
      degree_index++;
    }
  }

  // Finally, return the approximate pagerank vector
  p.makeCompressed();
  r.makeCompressed();
  return {p, r};
}

//------------------------------------------------------------------------------
// Sweep set implementation
//------------------------------------------------------------------------------
std::vector<StagInt> sweep_set_conductance_inner(stag::LocalGraph* graph,
                                                 SprsMat& vec,
                                                 double total_volume) {
  // If the provided total volume of the graph is not 0, then we use the correct
  // min(vol(S), vol(V \ S)) on the denominator.

  // The given vector must be one dimensional
  assert(vec.cols() == 1);

  // Sort the indices according to the values in vec.
  std::vector<double> orig_values = stag::sprsMatValues(&vec);
  std::vector<StagInt> sorted_indices = stag::sprsMatInnerIndices(&vec);
  std::stable_sort(sorted_indices.begin(), sorted_indices.end(),
                   [&vec](StagInt i1, StagInt i2) {return vec.coeff(i1, 0) > vec.coeff(i2, 0);});

  // Get the degrees of every vertex in the vector
  std::vector<double> degrees = graph->degrees(sorted_indices);

  // We will iterate through the indices in order of their increasing value
  // and add them to the vertex_set
  std::unordered_set<StagInt> vertex_set;
  double cut_weight = 0;
  double set_volume = 0;
  double set_complement_volume = total_volume;
  double best_conductance = 2;
  StagInt best_idx = 0;
  StagInt current_idx = 0;
  for (StagInt v : sorted_indices) {
    // The index we store is 'one more' than the index we are looking at.
    // See the constructor in the return statement of this method.
    current_idx++;

    // Add the next vertex to the vertex set
    vertex_set.insert(v);

    // Update the vertex set volume
    set_volume += degrees.at(current_idx - 1);
    set_complement_volume = total_volume - set_volume;

    // Update the cut weight. We need to add the total degree of the node v,
    // and then remove any edges from v to the rest of the vertex set.
    cut_weight += degrees.at(current_idx - 1);
    for (stag::edge e : graph->neighbors(v)) {
      if (vertex_set.contains(e.v2)) cut_weight -= 2 * e.weight;
    }

    // Check whether the current conductance is the best we've seen
    double this_denominator;
    if (total_volume == 0) {
      this_denominator = set_volume;
    } else {
      this_denominator = std::min(set_volume, set_complement_volume);
    }
    if (this_denominator > 0 && cut_weight / this_denominator < best_conductance) {
      best_conductance = cut_weight / this_denominator ;
      best_idx = current_idx;
    }
  }

  // Return the best cut
  return {sorted_indices.begin(), sorted_indices.begin() + best_idx};
}

std::vector<StagInt> sweep_set_conductance_inner(stag::LocalGraph* graph,
                                                 SprsMat& vec) {
  return sweep_set_conductance_inner(graph, vec, 0);
}

std::vector<StagInt> stag::sweep_set_conductance(stag::Graph* graph,
                                                 SprsMat& vec) {
  return sweep_set_conductance_inner(graph, vec, graph->total_volume());
}

std::vector<StagInt> stag::sweep_set_conductance(stag::LocalGraph* graph,
                                                 SprsMat& vec) {
  return sweep_set_conductance_inner(graph, vec);
}

//------------------------------------------------------------------------------
// Clustering metrics
//------------------------------------------------------------------------------
StagInt nChoose2(StagInt n)
{
  if (n < 0) throw std::invalid_argument("n must be non-negative.");
  if (n < 2) return 0;

  return n * (n-1) / 2;
}

Eigen::MatrixXi contingency_table(std::vector<StagInt>& gt_labels,
                                  std::vector<StagInt>& labels) {
  auto n = (StagInt) gt_labels.size();
  if ((StagInt) labels.size() != n) {
    throw std::invalid_argument("Label vectors must be the same size.");
  }

  // Find the number of clusters
  StagInt k = 1;
  for (auto label : gt_labels) {
    if (label < 0) throw std::invalid_argument("Cluster labels must be non-negative.");
    if (label > k - 1) {
      k = label + 1;
    }
  }
  for (auto label: labels) {
    if (label < 0) throw std::invalid_argument("Cluster labels must be non-negative.");
    if (label > k - 1) {
      k = label + 1;
    }
  }

  // Construct the k by k contingency table
  Eigen::MatrixXi contingency(k, k);

  // Initialize everything to 0
  for (auto i = 0; i < k; i++) {
    for (auto j = 0; j < k; j++) {
      contingency(i, j) = 0;
    }
  }

  for (auto i = 0; i < n; i++) {
    contingency(gt_labels.at(i), labels.at(i))++;
  }

  return contingency;
}

std::unordered_map<StagInt, StagInt> compute_cluster_sizes(
    std::vector<StagInt>& labels){
  std::unordered_map<StagInt, StagInt> cluster_sizes;
  for (StagInt l : labels) {
    cluster_sizes[l]++;
  }
  return cluster_sizes;
}

double stag::adjusted_rand_index(std::vector<StagInt>& gt_labels,
                                 std::vector<StagInt>& labels) {
  auto n = (StagInt) gt_labels.size();

  // Start by constructing the k by k contingency table
  // and the sizes of every cluster
  std::unordered_map<StagInt, StagInt> gt_sizes = compute_cluster_sizes(gt_labels);
  std::unordered_map<StagInt, StagInt> label_sizes = compute_cluster_sizes(labels);
  Eigen::MatrixXi contingency = contingency_table(gt_labels, labels);
  StagInt k = contingency.rows();

  // Now compute three components of the ARI
  // See https://stats.stackexchange.com/questions/207366/calculating-the-adjusted-rand-index.
  StagInt c1 = 0;
  for (auto i = 0; i < k; i++) {
    for (auto j = 0; j < k; j++) {
      c1 += nChoose2(contingency(i, j));
    }
  }

  StagInt c2 = 0;
  StagInt c3 = 0;
  for (auto i = 0; i < k; i++) {
    c2 += nChoose2(gt_sizes[i]);
    c3 += nChoose2(label_sizes[i]);
  }

  StagInt nC2 = nChoose2(n);

  return (c1 - ((double) (c2 * c3) / nC2)) / (((double) (c2 + c3)/2) - ((double) (c2 * c3)/nC2));
}

double entropy(std::vector<StagInt>& labels){
  // The entropy of a clustering is defined to be
  //   H(U) = sum p_i log(p_i)
  // where p_i = |C_i| / N is the proportion of the items in the ith cluster.
  auto N = (double) labels.size();

  // We will build a map of cluster IDs to the number of elements in that
  // cluster
  std::unordered_map<StagInt, StagInt> cluster_sizes =
      compute_cluster_sizes(labels);

  // Compute the entropy
  double entropy = 0;
  for (auto& it : cluster_sizes) {
    if (it.second > 0) {
      entropy += ((double) it.second / N) * std::log2(N / (double) it.second);
    }
  }
  return entropy;
}

double stag::mutual_information(std::vector<StagInt> &gt_labels,
                                std::vector<StagInt> &labels) {
  // Start by constructing the k by k contingency table
  // and the sizes of every cluster
  std::unordered_map<StagInt, StagInt> gt_sizes = compute_cluster_sizes(gt_labels);
  std::unordered_map<StagInt, StagInt> label_sizes = compute_cluster_sizes(labels);
  Eigen::MatrixXi contingency = contingency_table(gt_labels, labels);
  StagInt k = contingency.rows();
  auto N = (double) labels.size();

  // Compute the MI
  double mi = 0;
  for (StagInt k1 = 0; k1 < k; k1++) {
    for (StagInt k2 = 0; k2 < k; k2++) {
      if (contingency(k1, k2) > 0) {
        mi += (contingency(k1, k2) / N) * std::log2((N * contingency(k1, k2)) /
                                                    ((double) gt_sizes[k1] * (double) label_sizes[k2]));
      }
    }
  }
  return mi;
}

double stag::normalised_mutual_information(std::vector<StagInt> &gt_labels,
                                           std::vector<StagInt> &labels) {
  double mi = mutual_information(gt_labels, labels);
  double gt_entropy = entropy(gt_labels);
  double label_entropy = entropy(labels);
  double mean_entropy = (gt_entropy + label_entropy) / 2;
  return mi / mean_entropy;
}

double stag::conductance(stag::LocalGraph* graph, std::vector<StagInt>& cluster) {
  // We need to efficiently check whether a given node is in the cluster.
  // For this purpose, convert the list of cluster nodes into a set.
  std::set<StagInt> cluster_set(cluster.begin(), cluster.end());

  // For each node in the cluster, we will query its neighbours, and update
  // the numerator and denominator of the conductance formula.
  double cut = 0;
  double volume = 0;
  for (auto v : cluster) {
    std::vector<stag::edge> neighbors = graph->neighbors(v);
    double this_deg = 0;
    for (stag::edge e : neighbors) {
      this_deg += e.weight;

      // Check whether the neighbor of v is in the set
      if (!cluster_set.count(e.v2)) {
        // v2 is not in the cluster - add this weight to the cut
        cut += e.weight;
      }
    }
    volume += this_deg;
  }

  // Check for a divide-by-zero: volume of a set with volume 0 is defined to be
  // 0.
  if (volume == 0) return 0;
  else return cut / volume;
}

std::vector<StagInt> stag::symmetric_difference(std::vector<StagInt> &S,
                                                std::vector<StagInt> &T) {
  // For the later steps, we assume that the vectors are sorted
  std::sort(S.begin(), S.end());
  std::sort(T.begin(), T.end());

  // Remove duplicates from the vectors
  auto last_idx = std::unique(S.begin(), S.end());
  S.erase(last_idx, S.end());
  last_idx = std::unique(T.begin(), T.end());
  T.erase(last_idx, T.end());

  // Compute and return the symmetric difference
  std::vector<StagInt> symmetric_difference;
  std::set_symmetric_difference(S.begin(), S.end(),
                                T.begin(), T.end(),
                                std::back_inserter(symmetric_difference));
  return symmetric_difference;
}

//------------------------------------------------------------------------------
// Find connected components
//------------------------------------------------------------------------------
std::vector<StagInt> stag::connected_component(stag::LocalGraph* g,
                                               StagInt v) {
  // We track the connected component as both a set and a vector for efficient
  // algorithm.
  std::unordered_set<StagInt> component_set = {v};
  std::vector<StagInt> component_vec = {v};
  std::vector<StagInt> frontier = {v};

  // Perform a breadth-first search from v, adding all discovered nodes to the
  // connected component to return.
  while (!frontier.empty()) {
    // Pop the next vertex to search off the back of the frontier.
    StagInt this_vertex = frontier.back();
    frontier.pop_back();

    // Iterate through the neighbours of this vertex
    for (auto n : g->neighbors_unweighted(this_vertex)) {
      // If we have not seen the neighbour before, add it to the connected
      // component, and the frontier of the search.
      if (component_set.find(n) == component_set.end()) {
        component_set.insert(n);
        component_vec.push_back(n);
        frontier.push_back(n);
      }
    }
  }

  // Return the component that we found
  return component_vec;
}

std::vector<std::vector<StagInt>> stag::connected_components(
    stag::Graph* g) {
  // The components to be returned
  std::vector<std::vector<StagInt>> components;

  // Keep track of the vertices which we've already found
  std::unordered_set<StagInt> found_vertices;

  for (StagInt v = 0; v < g->number_of_vertices(); v++) {
    // If this vertex has not been returned as part of one of the connected
    // components yet, then we find the connected component containing it.
    if (!found_vertices.contains(v)) {
      std::vector<StagInt> this_component = stag::connected_component(g, v);
      components.push_back(this_component);

      // Track that we've visited every node in the newly found component.
      for (auto idx : this_component) found_vertices.insert(idx);
    }
  }

  return components;
}

//------------------------------------------------------------------------------
// Constructing approximate similarity graph.
//------------------------------------------------------------------------------
class KDETreeEntry {
public:
  KDETreeEntry(DenseMat* data, StagReal a, StagInt min_id, StagInt max_id, StagInt dep, StagInt dep_cutoff)
      : sampling_dist(0.0, 1.0)
  {
    min_idx = min_id;
    max_idx = max_id;
    n_node = max_idx - min_idx;
    depth = dep;
    depth_cutoff = dep_cutoff;

    if (max_idx - min_idx >= ASG_TREE_CUTOFF) {
      below_cutoff = false;

      // We only initialise an estimator at this node if we are below the depth
      // depth cutoff. Above the depth cutoff and we just use the children to estimate.
      if (depth > depth_cutoff) initialise_estimator(data, a);
    } else {
      below_cutoff = true;
      initialise_exact(data, a);
    }

    if (!below_cutoff) {
      StagInt midpoint = (min_idx + max_idx) / 2;
      assert(midpoint >= min_idx);
      assert(midpoint < max_idx);
      left_child = new KDETreeEntry(data, a, min_idx, midpoint, depth + 1, depth_cutoff);
      right_child = new KDETreeEntry(data, a, midpoint + 1, max_idx, depth + 1, depth_cutoff);
    }
  }

  void initialise_estimator(DenseMat* data, StagReal a) {
    StagInt K1 = ceil(1 * log((StagReal) n_node));
    StagReal K2_constant = 0.1;
    StagReal min_mu = 1.0 / (StagReal) n_node;
    StagInt offset = 0;
    this_estimator = stag::CKNSGaussianKDE(data, a, min_mu, K1, K2_constant,
                                           offset, min_idx, max_idx + 1);
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
        if (depth > depth_cutoff) {
          // If we are below the cutoff depth, query the estimator.
          weight = (StagReal) n_node * this_estimator.query(q);
        } else {
          // If we are above the cutoff depth, get the estimate from the
          // children.
          weight = left_child->estimate_weight(q, q_id) + right_child->estimate_weight(q, q_id);
        }
      } else {
        weight = (StagReal) n_node * exact_kde.query(q);
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

  std::vector<StagReal> estimate_weights(DenseMat* q, indicators::ProgressBar* prog_bar,
                                         StagInt progress_start, StagInt progress_end) {
    std::vector<StagReal> weights;
    if (!below_cutoff) {
      if (depth > depth_cutoff) {
        weights = this_estimator.query(q);

        for (auto & weight : weights) {
          weight = (StagReal) n_node * weight;
        }
      } else {
        auto progress_middle = (StagInt) (progress_start + progress_end) / 2;
        std::vector<StagReal> left_weights = left_child->estimate_weights(q, prog_bar, progress_start, progress_middle);
        std::vector<StagReal> right_weights = right_child->estimate_weights(q, prog_bar, progress_middle, progress_end);
        assert(left_weights.size() == right_weights.size());
        for (StagInt i = 0; i < (StagInt) left_weights.size(); i++) {
          weights.push_back(left_weights.at(i) + right_weights.at(i));
        }
      }
    } else {
      weights = exact_kde.query(q);

      for (auto & weight : weights) {
        weight = (StagReal) n_node * weight;
      }
    }

    for (auto q_id = 0; q_id < (StagInt) weights.size(); q_id++) {
      cached_weights[q_id] = weights.at(q_id);
    }

    prog_bar->set_progress(progress_end);

    return weights;
  }

  std::vector<StagInt> sample_neighbors(const stag::DataPoint& q, const StagInt q_id, StagInt num_to_sample) {
    if (below_cutoff) {
      std::vector<StagReal> rs;
      for (auto i = 0; i < num_to_sample; i++) {
        rs.push_back(sampling_dist(*stag::get_global_rng()));
      }
      StagReal deg = estimate_weight(q, q_id); // should be cached
      return exact_kde.sample_neighbors(q, deg, rs);
    } else {
      StagReal left_est = left_child->estimate_weight(q, q_id);
      StagReal right_est = right_child->estimate_weight(q, q_id);
      StagReal my_est = left_est + right_est;

      // Get the number of left samples from a binomial distribution.
      std::binomial_distribution<> d(num_to_sample, left_est / my_est);
      StagInt num_left_samples = d(*stag::get_global_rng());

      assert(num_left_samples >= 0);
      std::vector<StagInt> left_samples;
      std::vector<StagInt> right_samples;
      if (num_left_samples > 0) left_samples = left_child->sample_neighbors(q, q_id, num_left_samples);
      if (num_left_samples < num_to_sample) right_samples = right_child->sample_neighbors(q, q_id, num_to_sample - num_left_samples);
      assert((StagInt) left_samples.size() + (StagInt) right_samples.size() == num_to_sample);

      std::vector<StagInt> samples;
      samples.reserve(num_to_sample);
      samples.insert(samples.end(), left_samples.begin(), left_samples.end());
      samples.insert(samples.end(), right_samples.begin(), right_samples.end());
      return samples;
    }
  }

private:
  bool below_cutoff;
  stag::CKNSGaussianKDE this_estimator;
  stag::ExactGaussianKDE exact_kde;
  StagInt min_idx;
  StagInt max_idx;
  StagInt n_node;
  KDETreeEntry* left_child;
  KDETreeEntry* right_child;
  std::uniform_real_distribution<double> sampling_dist;
  std::unordered_map<StagInt, StagReal> cached_weights;
  std::mutex cache_mutex;
  StagInt depth;
  StagInt depth_cutoff;
};

void sample_asg_edges(DenseMat* data,
                      KDETreeEntry& tree_root,
                      StagInt chunk_start,
                      StagInt chunk_end,
                      std::vector<EdgeTriplet>& edges,
                      StagInt edges_per_node,
                      std::vector<StagReal>& degrees,
                      indicators::ProgressBar* prog_bar) {
  auto nodes_per_tick = (StagInt) ceil((StagReal) data->rows() / 70);
  std::uniform_real_distribution<double> sampling_dist(0, 0.8);
  auto progress_offset = (StagInt) (sampling_dist(*stag::get_global_rng()) * nodes_per_tick);
  auto total_ticks = (StagInt) (chunk_end - chunk_start) / nodes_per_tick;
  StagInt ticks_made = 0;

  StagInt next_index = 2 * (chunk_start * edges_per_node);
  for (auto i = chunk_start; i < chunk_end; i++) {
    stag::DataPoint q(*data, i);

    std::vector<StagInt> node_samples = tree_root.sample_neighbors(q, i, edges_per_node);

    assert(i <= (StagInt) degrees.size());
    StagReal added_weight = degrees[i] / 5;

    for (auto j = 0; j < edges_per_node; j++) {
      assert(j < (StagInt) node_samples.size());
      StagInt neighbor = node_samples[j];
      assert(next_index < (StagInt) edges.size());
      edges[next_index] = EdgeTriplet(i, neighbor, added_weight);
      next_index++;
      edges[next_index] = EdgeTriplet(neighbor, i, added_weight);
      next_index++;
    }

    // Update the progress bar
    if ((i - chunk_start + progress_offset) % nodes_per_tick == 0) {
      if (ticks_made < total_ticks) {
        prog_bar->tick();
        ticks_made++;
      }
    }
  }
}

stag::Graph stag::approximate_similarity_graph(DenseMat* data, StagReal a) {
  // Work out how many total edges we will sample
  auto edges_per_node = (StagInt) (3 * log((StagReal) data->rows()));
  StagInt num_edges = data->rows() * edges_per_node;

  // Track the construction with a progress bar
  indicators::ProgressBar prog_bar{indicators::option::BarWidth{50},
                                   indicators::option::Start{"["},
                                   indicators::option::Fill{"■"},
                                   indicators::option::Lead{"■"},
                                   indicators::option::Remainder{"-"},
                                   indicators::option::End{"]"},
                                   indicators::option::ForegroundColor{indicators::Color::cyan},
                                   indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                                   indicators::option::PostfixText("Initialising KDE Estimators")};
  prog_bar.set_progress(1);

  // Begin by creating the tree of kernel density estimators.
  // Creating each node is parallelized by the KDE code.
  KDETreeEntry tree_root(data, a, 0, data->rows() - 1, 0, (StagInt) log((StagReal) edges_per_node));

  prog_bar.set_option(indicators::option::PostfixText("Estimating node degrees"));
  prog_bar.set_progress(10);

  // First, compute the degrees of all nodes
  std::vector<StagReal> degrees = tree_root.estimate_weights(data, &prog_bar, 10, 30);

  prog_bar.set_option(indicators::option::PostfixText("Sampling edges"));
  prog_bar.set_progress(30);

  StagInt num_threads = std::thread::hardware_concurrency();
  std::vector<EdgeTriplet> graph_edges(2 * num_edges);
  if (data->rows() <= num_threads * 2) {
    sample_asg_edges(data, tree_root, 0, data->rows(), graph_edges, edges_per_node, degrees, &prog_bar);
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
                                 this_chunk_end, graph_edges, edges_per_node, degrees, &prog_bar);
              }
          )
      );
    }

    // Join the futures
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      futures[chunk_id].get();
    }
  }

  prog_bar.set_option(indicators::option::PostfixText("Done!"));
  prog_bar.set_progress(100);

  // Return a graph
  SprsMat adj_mat(data->rows(), data->rows());
  adj_mat.setFromTriplets(graph_edges.begin(), graph_edges.end());
  return stag::Graph(adj_mat);
}

stag::Graph stag::similarity_graph(DenseMat* data, StagReal a) {
  StagInt n = data->rows();
  std::vector<EdgeTriplet> graph_edges;

  std::vector<stag::DataPoint> datapoints = stag::matrix_to_datapoints(data);

  StagInt num_threads = std::thread::hardware_concurrency();
  if (n <= num_threads * 2) {
    for (StagInt i = 0; i < n; i++) {
      for (StagInt j = i; j < n; j++) {
        StagReal weight = stag::gaussian_kernel(a, datapoints.at(i), datapoints.at(j));
        EdgeTriplet e1(i, j, weight);
        graph_edges.push_back(e1);
        if (i != j) {
          EdgeTriplet e2(j, i, weight);
          graph_edges.push_back(e2);
        }
      }
    }
  } else {
    // Start the thread pool
    ctpl::thread_pool pool((int) num_threads);
    std::vector<std::future<std::vector<EdgeTriplet>>> futures;

    StagInt chunk_size = floor((StagReal) n / (StagReal) num_threads);

    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      StagInt this_chunk_start = chunk_id * chunk_size;
      StagInt this_chunk_end = this_chunk_start + chunk_size;
      if (chunk_id == num_threads - 1) this_chunk_end = n;

      assert(this_chunk_start <= n);
      assert(this_chunk_end <= n);
      assert(this_chunk_end > this_chunk_start);

      futures.push_back(
          pool.push(
              [&, this_chunk_start, this_chunk_end] (int id) {
                ignore_warning(id);
                std::vector<EdgeTriplet> this_thread_edges;
                for (StagInt i = this_chunk_start; i < this_chunk_end; i++) {
                  for (StagInt j = i; j < n; j++) {
                    StagReal weight = stag::gaussian_kernel(a, datapoints.at(i), datapoints.at(j));
                    EdgeTriplet e1(i, j, weight);
                    this_thread_edges.push_back(e1);
                    if (i != j) {
                      EdgeTriplet e2(j, i, weight);
                      this_thread_edges.push_back(e2);
                    }
                  }
                }
                return this_thread_edges;
              }
          )
      );
    }

    // Join the futures
    for (auto chunk_id = 0; chunk_id < num_threads; chunk_id++) {
      std::vector<EdgeTriplet> this_thread_edges = futures[chunk_id].get();
      graph_edges.insert(graph_edges.end(), this_thread_edges.begin(), this_thread_edges.end());
    }

    assert((StagInt) graph_edges.size() == n * n);
  }

  // Return a graph
  SprsMat adj_mat(n, n);
  adj_mat.setFromTriplets(graph_edges.begin(), graph_edges.end());
  return stag::Graph(adj_mat);
}
