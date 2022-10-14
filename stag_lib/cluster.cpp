/**
* Implementations of the clustering algorithms defined in cluster.h.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
*/
#include <vector>
#include <queue>
#include <unordered_set>
#include <stdexcept>
#include <algorithm>
#include <Eigen/Sparse>
#include <graph.h>
#include <cluster.h>
#include <utility.h>

std::vector<stag_int> stag::local_cluster(stag::LocalGraph *graph, stag_int seed_vertex, stag_int target_volume) {
  // The 'locality' parameter should essentially be a constant if we expect the conductance
  // of the target cluster to be a constant. Set it to 0.01.
  // The error parameter should decrease as the inverse of the target volume.
  return stag::local_cluster_acl(graph,
                                 seed_vertex,
                                 0.01,
                                 1./ (double) target_volume);
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
  // Compute the approximate pagerank vector
  SprsMat seedDist(seed_vertex + 1, 1);
  seedDist.coeffRef(seed_vertex, 0) = 1;
  std::tuple<SprsMat, SprsMat> apr = stag::approximate_pagerank(graph, seedDist, locality, error);
  SprsMat p = std::get<0>(apr);

  // Normalise the values in p by their degree
  double deg;
  for (SprsMat::InnerIterator it(p, 0); it; ++it) {
    deg = graph->degree(it.row());
    it.valueRef() = it.value() / deg;
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
  // Check that the provided seed vector is a column vector
  if (seed_vector.cols() > 1) throw std::invalid_argument("Seed vector must be a column vector.");

  // Initialise p to be the all-zeros vector.
  SprsMat p(seed_vector.rows(), 1);

  // Initialise r to be equal to the seed vector
  SprsMat r(seed_vector);

  // We will maintain a queue of vertices satisfying
  //    r(u) >= epsilon * deg(u)
  // Along with the queue, maintain an unordered set which tracks the vertices
  // already in the queue.
  stag_int u;
  double deg;
  std::queue<stag_int> vertex_queue;
  std::unordered_set<stag_int> queue_members;
  for (SprsMat::InnerIterator it(seed_vector, 0); it; ++it) {
    u = it.row();
    deg = graph->degree(u);
    if (r.coeff(u, 0) >= epsilon * deg) {
      vertex_queue.push(u);
      queue_members.insert(u);
    }
  }

  // While the queue is not empty, push an entry off the queue and apply the
  // push operation.
  while (!vertex_queue.empty()) {
    // Get the next vertex from the queue
    u = vertex_queue.front();
    vertex_queue.pop();
    queue_members.erase(u);

    // Perform the push operation on this vertex
    push(graph, &p, &r, alpha, u);

    // Check u to see if it should be added back to the queue
    deg = graph->degree(u);
    if (r.coeff(u, 0) >= epsilon * deg) {
      vertex_queue.push(u);
      queue_members.insert(u);
    }

    // Check the neighbors of u to see if they should be added back to the queue
    // Skip any neighbors which are already in the queue.
    stag_int v;
    for (stag::edge e : graph->neighbors(u)) {
      v = e.v2;
      deg = graph->degree(v);
      if (r.coeff(e.v2, 0) >= epsilon * deg && !queue_members.contains(v)) {
        vertex_queue.push(v);
        queue_members.insert(v);
      }
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
    set_volume += graph->degree(v);

    // Update the cut weight. We need to add the total degree of the node v,
    // and then remove any edges from v to the rest of the vertex set.
    cut_weight += graph->degree(v);
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
