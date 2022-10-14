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
#include <iostream>
#include <Eigen/Sparse>
#include <graph.h>
#include <cluster.h>
#include <utility.h>

std::vector<stag_int> stag::local_cluster(stag::LocalGraph *graph, stag_int seed_vertex) {
  return stag::local_cluster_acl(graph, seed_vertex);
}

//------------------------------------------------------------------------------
// Implementation of ACL Local Clustering Algorithm based on PageRank vectors
//------------------------------------------------------------------------------
std::vector<stag_int> stag::local_cluster_acl(stag::LocalGraph *graph,
                                              stag_int seed_vertex) {
  // Compute the approximate pagerank vector
  SprsMat seedDist(seed_vertex + 1, 1);
  seedDist.coeffRef(seed_vertex, 0) = 1;
  std::tuple<SprsMat, SprsMat> apr = stag::approximate_pagerank(graph, seedDist, 0.1, 0.01);
  SprsMat p = std::get<0>(apr);
  return stag::sprsMatInnerIndices(&p);
}

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
    // The neighbor 'v' is the 'second' element in the returned edge
    v = e.u;
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
      v = e.u;
      deg = graph->degree(v);
      if (r.coeff(e.u, 0) >= epsilon * deg && !queue_members.contains(v)) {
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
