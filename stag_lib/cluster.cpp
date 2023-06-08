//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#include <vector>
#include <set>
#include <deque>
#include <unordered_set>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <Eigen/Sparse>
#include <graph.h>
#include <cluster.h>
#include <utility.h>
#include <spectrum.h>
#include <KMeansRex/KMeansRexCoreInterface.h>

std::vector<stag_int> stag::spectral_cluster(stag::Graph *graph, stag_int k) {
  // Check that the number of clusters is valid.
  if (k < 1 || k > graph->number_of_vertices() /2) {
    throw std::invalid_argument("Number of clusters must be between 1 and n/2.");
  }

  // Start by computing the 'first' k eigenvalues of the normalised graph
  // laplacian matrix.
  const SprsMat* lap = graph->normalised_laplacian();
  Eigen::MatrixXd eigvecs = stag::compute_eigenvectors(lap, k);

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

std::vector<stag_int> stag::cheeger_cut(stag::Graph* graph) {
  // First, compute the first 2 eigenvectors of the normalised graph Laplacian
  // matrix.
  const SprsMat* lap = graph->normalised_laplacian();
  stag::EigenSystem eigsys = stag::compute_eigensystem(lap, 2);

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
  for (stag_int i = 0; i < graph->number_of_vertices(); i++) {
    vec.coeffRef(i, 0) = std::sqrt(1 / graph->degree(i)) * vec.coeff(i, 0);
  }

  // Perform the sweep
  std::vector<stag_int> clusters (graph->number_of_vertices(), 0);
  std::vector<stag_int> cut = stag::sweep_set_conductance(graph, vec);
  for (stag_int i : cut) clusters.at(i) = 1;
  return clusters;
}

std::vector<stag_int> stag::local_cluster(stag::LocalGraph *graph, stag_int seed_vertex, double target_volume) {
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
std::vector<stag_int> stag::local_cluster_acl(stag::LocalGraph* graph,
                                              stag_int seed_vertex,
                                              double locality) {
  return stag::local_cluster_acl(graph, seed_vertex, locality, 0.001);
}

std::vector<stag_int> stag::local_cluster_acl(stag::LocalGraph *graph,
                                              stag_int seed_vertex,
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
  stag_int degree_index = 0;
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
void push(stag::LocalGraph *graph, SprsMat *p, SprsMat *r, double alpha, stag_int u) {
  // First, make sure that the vector p is large enough.
  if (u >= p->rows()) {
    p->conservativeResize(u + 1, 1);
  }

  // Update p according to the push operation
  p->coeffRef(u, 0) = p->coeff(u, 0) + alpha * r->coeff(u, 0);

  // Iterate through the neighbors of u
  double deg = graph->degree(u);
  stag_int v;
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
  stag_int u;
  std::vector<double> degrees = graph->degrees(
      stag::sprsMatInnerIndices(&seed_vector));
  std::deque<stag_int> vertex_queue;
  std::unordered_set<stag_int> queue_members;
  stag_int degree_index = 0;
  for (SprsMat::InnerIterator it(seed_vector, 0); it; ++it) {
    u = it.row();
    if (r.coeff(u, 0) >= epsilon * degrees.at(degree_index)) {
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
    std::vector<stag_int> neighbors = graph->neighbors_unweighted(u);
    std::vector<double> neighbor_degrees = graph->degrees(neighbors);

    // The length of neighbors and neighbor_degrees should always be equal.
    // If they are not, there must be a bug in the implementation of graph->degrees.
    assert(neighbors.size() == neighbor_degrees.size());

    degree_index = 0;
    stag_int v;
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
std::vector<stag_int> sweep_set_conductance_inner(stag::LocalGraph* graph,
                                                  SprsMat& vec,
                                                  double total_volume) {
  // If the provided total volume of the graph is not 0, then we use the correct
  // min(vol(S), vol(V \ S)) on the denominator.

  // The given vector must be one dimensional
  assert(vec.cols() == 1);

  // Sort the indices according to the values in vec.
  std::vector<double> orig_values = stag::sprsMatValues(&vec);
  std::vector<stag_int> sorted_indices = stag::sprsMatInnerIndices(&vec);
  std::stable_sort(sorted_indices.begin(), sorted_indices.end(),
                   [&vec](stag_int i1, stag_int i2) {return vec.coeff(i1, 0) > vec.coeff(i2, 0);});

  // Get the degrees of every vertex in the vector
  std::vector<double> degrees = graph->degrees(sorted_indices);

  // We will iterate through the indices in order of their increasing value
  // and add them to the vertex_set
  std::unordered_set<stag_int> vertex_set;
  double cut_weight = 0;
  double set_volume = 0;
  double set_complement_volume = total_volume;
  double best_conductance = 2;
  stag_int best_idx = 0;
  stag_int current_idx = 0;
  for (stag_int v : sorted_indices) {
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

std::vector<stag_int> sweep_set_conductance_inner(stag::LocalGraph* graph,
                                                  SprsMat& vec) {
  return sweep_set_conductance_inner(graph, vec, 0);
}

std::vector<stag_int> stag::sweep_set_conductance(stag::Graph* graph,
                                                  SprsMat& vec) {
  return sweep_set_conductance_inner(graph, vec, graph->total_volume());
}

std::vector<stag_int> stag::sweep_set_conductance(stag::LocalGraph* graph,
                                                  SprsMat& vec) {
  return sweep_set_conductance_inner(graph, vec);
}

//------------------------------------------------------------------------------
// Clustering metrics
//------------------------------------------------------------------------------
stag_int nChoose2(stag_int n)
{
  if (n < 0) throw std::invalid_argument("n must be non-negative.");
  if (n < 2) return 0;

  return n * (n-1) / 2;
}

Eigen::MatrixXi contingency_table(std::vector<stag_int>& gt_labels,
                                  std::vector<stag_int>& labels) {
  auto n = (stag_int) gt_labels.size();
  if ((stag_int) labels.size() != n) {
    throw std::invalid_argument("Label vectors must be the same size.");
  }

  // Find the number of clusters
  stag_int k = 1;
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

std::unordered_map<stag_int, stag_int> compute_cluster_sizes(
    std::vector<stag_int>& labels){
  std::unordered_map<stag_int, stag_int> cluster_sizes;
  for (stag_int l : labels) {
    cluster_sizes[l]++;
  }
  return cluster_sizes;
}

double stag::adjusted_rand_index(std::vector<stag_int>& gt_labels,
                                 std::vector<stag_int>& labels) {
  auto n = (stag_int) gt_labels.size();

  // Start by constructing the k by k contingency table
  // and the sizes of every cluster
  std::unordered_map<stag_int, stag_int> gt_sizes = compute_cluster_sizes(gt_labels);
  std::unordered_map<stag_int, stag_int> label_sizes = compute_cluster_sizes(labels);
  Eigen::MatrixXi contingency = contingency_table(gt_labels, labels);
  stag_int k = contingency.rows();

  // Now compute three components of the ARI
  // See https://stats.stackexchange.com/questions/207366/calculating-the-adjusted-rand-index.
  stag_int c1 = 0;
  for (auto i = 0; i < k; i++) {
    for (auto j = 0; j < k; j++) {
      c1 += nChoose2(contingency(i, j));
    }
  }

  stag_int c2 = 0;
  stag_int c3 = 0;
  for (auto i = 0; i < k; i++) {
    c2 += nChoose2(gt_sizes[i]);
    c3 += nChoose2(label_sizes[i]);
  }

  stag_int nC2 = nChoose2(n);

  return (c1 - ((double) (c2 * c3) / nC2)) / (((double) (c2 + c3)/2) - ((double) (c2 * c3)/nC2));
}

double entropy(std::vector<stag_int>& labels){
  // The entropy of a clustering is defined to be
  //   H(U) = sum p_i log(p_i)
  // where p_i = |C_i| / N is the proportion of the items in the ith cluster.
  auto N = (double) labels.size();

  // We will build a map of cluster IDs to the number of elements in that
  // cluster
  std::unordered_map<stag_int, stag_int> cluster_sizes =
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

double stag::mutual_information(std::vector<stag_int> &gt_labels,
                                std::vector<stag_int> &labels) {
  // Start by constructing the k by k contingency table
  // and the sizes of every cluster
  std::unordered_map<stag_int, stag_int> gt_sizes = compute_cluster_sizes(gt_labels);
  std::unordered_map<stag_int, stag_int> label_sizes = compute_cluster_sizes(labels);
  Eigen::MatrixXi contingency = contingency_table(gt_labels, labels);
  stag_int k = contingency.rows();
  auto N = (double) labels.size();

  // Compute the MI
  double mi = 0;
  for (stag_int k1 = 0; k1 < k; k1++) {
    for (stag_int k2 = 0; k2 < k; k2++) {
      if (contingency(k1, k2) > 0) {
        mi += (contingency(k1, k2) / N) * std::log2((N * contingency(k1, k2)) /
                                                    ((double) gt_sizes[k1] * (double) label_sizes[k2]));
      }
    }
  }
  return mi;
}

double stag::normalised_mutual_information(std::vector<stag_int> &gt_labels,
                                           std::vector<stag_int> &labels) {
  double mi = mutual_information(gt_labels, labels);
  double gt_entropy = entropy(gt_labels);
  double label_entropy = entropy(labels);
  double mean_entropy = (gt_entropy + label_entropy) / 2;
  return mi / mean_entropy;
}

double stag::conductance(stag::LocalGraph* graph, std::vector<stag_int>& cluster) {
  // We need to efficiently check whether a given node is in the cluster.
  // For this purpose, convert the list of cluster nodes into a set.
  std::set<stag_int> cluster_set(cluster.begin(), cluster.end());

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

std::vector<stag_int> stag::symmetric_difference(std::vector<stag_int> &S,
                                                 std::vector<stag_int> &T) {
  // For the later steps, we assume that the vectors are sorted
  std::sort(S.begin(), S.end());
  std::sort(T.begin(), T.end());

  // Remove duplicates from the vectors
  auto last_idx = std::unique(S.begin(), S.end());
  S.erase(last_idx, S.end());
  last_idx = std::unique(T.begin(), T.end());
  T.erase(last_idx, T.end());

  // Compute and return the symmetric difference
  std::vector<stag_int> symmetric_difference;
  std::set_symmetric_difference(S.begin(), S.end(),
                                T.begin(), T.end(),
                                std::back_inserter(symmetric_difference));
  return symmetric_difference;
}

//------------------------------------------------------------------------------
// Find connected components
//------------------------------------------------------------------------------
std::vector<stag_int> stag::connected_component(stag::LocalGraph* g,
                                                stag_int v) {
  // We track the connected component as both a set and a vector for efficient
  // algorithm.
  std::unordered_set<stag_int> component_set = {v};
  std::vector<stag_int> component_vec = {v};
  std::vector<stag_int> frontier = {v};

  // Perform a breadth-first search from v, adding all discovered nodes to the
  // connected component to return.
  while (!frontier.empty()) {
    // Pop the next vertex to search off the back of the frontier.
    stag_int this_vertex = frontier.back();
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

std::vector<std::vector<stag_int>> stag::connected_components(
    stag::Graph* g) {
  // The components to be returned
  std::vector<std::vector<stag_int>> components;

  // Keep track of the vertices which we've already found
  std::unordered_set<stag_int> found_vertices;

  for (stag_int v = 0; v < g->number_of_vertices(); v++) {
    // If this vertex has not been returned as part of one of the connected
    // components yet, then we find the connected component containing it.
    if (!found_vertices.contains(v)) {
      std::vector<stag_int> this_component = stag::connected_component(g, v);
      components.push_back(this_component);

      // Track that we've visited every node in the newly found component.
      for (auto idx : this_component) found_vertices.insert(idx);
    }
  }

  return components;
}
