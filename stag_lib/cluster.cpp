/**
* Implementations of the clustering algorithms defined in cluster.h.
 *
 * Copyright 2022 Peter Macgregor
*/
#include <vector>
#include <graph.h>
#include <cluster.h>

std::vector<int> stag::local_cluster(stag::LocalGraph *graph, int seed_vertex) {
  return stag::local_cluster_acl(graph, seed_vertex);
}

//------------------------------------------------------------------------------
// Implementation of ACL Local Clustering Algorithm based on PageRank vectors
//------------------------------------------------------------------------------
std::vector<int> stag::local_cluster_acl(stag::LocalGraph *graph,
                                        int seed_vertex) {
  // Compute the approximate pagerank vector
  SprsMat apr = stag::approximate_pagerank(seed_vertex, 1, 1);

  return {graph->neighbors_unweighted(seed_vertex)[0]};
}

SprsMat stag::approximate_pagerank(stag::LocalGraph *graph,
                                   SprsMat &seed_vector,
                                   double alpha,
                                   double epsilon) {
  return seed_vector;
}
