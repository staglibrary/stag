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

std::vector<int> stag::local_cluster_acl(stag::LocalGraph *graph,
                                        int seed_vertex) {
  return {graph->neighbors_unweighted(seed_vertex)[0]};
}
