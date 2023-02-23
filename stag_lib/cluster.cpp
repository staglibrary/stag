//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#include <vector>
#include <deque>
#include <unordered_set>
#include <stdexcept>
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
std::vector<stag_int> stag::sweep_set_conductance(stag::LocalGraph* graph,
                                                  SprsMat& vec) {
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

    // Update the cut weight. We need to add the total degree of the node v,
    // and then remove any edges from v to the rest of the vertex set.
    cut_weight += degrees.at(current_idx - 1);
    for (stag::edge e : graph->neighbors(v)) {
      if (vertex_set.contains(e.v2)) cut_weight -= 2 * e.weight;
    }

    // Check whether the current conductance is the best we've seen
    if (cut_weight / set_volume < best_conductance) {
      best_conductance = cut_weight / set_volume;
      best_idx = current_idx;
    }
  }

  // Return the best cut
  return {sorted_indices.begin(), sorted_indices.begin() + best_idx};
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

double stag::adjusted_rand_index(std::vector<stag_int>& gt_labels,
                                 std::vector<stag_int>& labels) {
  stag_int n = gt_labels.size();
  if (labels.size() != n) {
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

  // Start by constructing the k by k contingency table
  // and the sizes of every cluster
  Eigen::VectorXi gt_sizes(k);
  Eigen::VectorXi label_sizes(k);
  Eigen::MatrixXi contingency(k, k);

  // Initialize everything to 0
  for (auto i = 0; i < k; i++) {
    gt_sizes(i) = 0;
    label_sizes(i) = 0;
    for (auto j = 0; j < k; j++) {
      contingency(i, j) = 0;
    }
  }

  for (auto i = 0; i < n; i++) {
    contingency(gt_labels.at(i), labels.at(i))++;
    gt_sizes(gt_labels.at(i))++;
    label_sizes(labels.at(i))++;
  }

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
    c2 += nChoose2(gt_sizes(i));
    c3 += nChoose2(label_sizes(i));
  }

  stag_int nC2 = nChoose2(n);

  return (c1 - ((double) (c2 * c3) / nC2)) / (((double) (c2 + c3)/2) - ((double) (c2 * c3)/nC2));
}
