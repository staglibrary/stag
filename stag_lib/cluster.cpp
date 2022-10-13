/**
* Implementations of the clustering algorithms defined in cluster.h.
 *
 * Copyright 2022 Peter Macgregor
*/
#include <vector>
#include <Eigen/Sparse>
#include <graph.h>
#include <cluster.h>
#include <utility.h>

std::vector<int> stag::local_cluster(stag::LocalGraph *graph, int seed_vertex) {
  return stag::local_cluster_acl(graph, seed_vertex);
}

//------------------------------------------------------------------------------
// Implementation of ACL Local Clustering Algorithm based on PageRank vectors
//------------------------------------------------------------------------------
std::vector<int> stag::local_cluster_acl(stag::LocalGraph *graph,
                                        int seed_vertex) {
  // Compute the approximate pagerank vector
  SprsMat seedDist(seed_vertex + 1, 1);
  seedDist.coeffRef(seed_vertex, 0) = 1;
  SprsMat apr = stag::approximate_pagerank(graph, seedDist, 1, 1);
  apr.makeCompressed();
  return stag::sprsMatInnerIndices(&apr);
}

/**
 * Perform the push operation on vertex u in the given graph.
 *
 * Updates the sparse vectors p and r, according to:
 *   p'(u) = p(u) + alpha * r(u)
 *   r'(u) = (1 - alpha) * r(u) / 2
 *   r'(v) = r(v) + (1 - alpha) * w(u, v) * r(u) / (2 * deg(u))
 * for all v which are neighbors of u.
 *
 * @param graph
 * @param p
 * @param r
 * @param alpha
 * @param u
 */
void push(stag::LocalGraph *graph, SprsMat *p, SprsMat *r, double alpha, int u) {
  // First, make sure that the vector p is large enough.
  if (u >= p->rows()) {
    p->conservativeResize(u + 1, 1);
  }

  // Update p according to the push operation
  p->coeffRef(u, 0) = p->coeff(u, 0) + alpha * r->coeff(u, 0);

  // Update r according to the push operation
  r->coeffRef(u, 0) = (1 - alpha) * r->coeff(u, 0) / 2;

  // Iterate through the neighbors of u
  double deg = graph->degree(u);
  int v;
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
}

SprsMat stag::approximate_pagerank(stag::LocalGraph *graph,
                                   SprsMat &seed_vector,
                                   double alpha,
                                   double epsilon) {
  SprsMat p(seed_vector.rows(), 1);
  push(graph, &p, &seed_vector, alpha, 0);

  // Finally, return the approximate pagerank vector
  p.makeCompressed();
  return p;
}
