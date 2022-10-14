/**
 * Graph clustering algorithms based on spectral methods.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#ifndef STAG_TEST_CLUSTER_H
#define STAG_TEST_CLUSTER_H

#include <vector>

#include <graph.h>

namespace stag {
  /**
   * Default local clustering algorithm. Given a graph and starting vertex,
   * return a cluster which is close to the starting vertex.
   *
   * This method defaults to use the ACL local clustering algorithm.
   *
   * [ACL] Andersen, Reid, Fan Chung, and Kevin Lang.
   * "Local graph partitioning using pagerank vectors." 2006
   * 47th Annual IEEE Symposium on Foundations of Computer Science (FOCS'06). IEEE, 2006.
   *
   * @param graph - a graph object implementing the LocalGraph interface
   * @param seed_vertex - the starting vertex in the graph
   * @return a vector containing the indices of vectors considered to be in the
   *         same cluster as the seed_vertex.
   */
  std::vector<stag_int> local_cluster(stag::LocalGraph* graph, stag_int seed_vertex);

  /**
   * The ACL local clustering algorithm. Given a graph and starting vertex,
   * returns a cluster close to the starting vertex, constructed in a local way.
   *
   * [ACL] Andersen, Reid, Fan Chung, and Kevin Lang.
   * "Local graph partitioning using pagerank vectors." 2006
   * 47th Annual IEEE Symposium on Foundations of Computer Science (FOCS'06). IEEE, 2006.
   *
   * @param graph - a graph object implementing the LocalGraph interface
   * @param seed_vertex - the starting vertex in the graph
   * @return a vector containing the indices of vectors considered to be in the
   *         same cluster as the seed_vertex.
   */
  std::vector<stag_int> local_cluster_acl(stag::LocalGraph* graph, stag_int seed_vertex);

  /**
   * Compute the approximate pagerank vector as described in ACL:
   *
   * [ACL] Andersen, Reid, Fan Chung, and Kevin Lang.
   * "Local graph partitioning using pagerank vectors." 2006
   * 47th Annual IEEE Symposium on Foundations of Computer Science (FOCS'06). IEEE, 2006.
   *
   * The parameters s, alpha, and epsilon are used as described in the paper.
   *
   * Note that the dimension of the returned vectors may not match the true
   * number of vertices in the graph provided since the approximate
   * pagerank is computed locally.
   *
   * @param graph
   * @param seed_vector
   * @param alpha
   * @param epsilon
   * @return A tuple of sparse column vectors corresponding to
   *            p - the approximate pagerank vector
   *            r - the residual vector
   *         By the definition of approximate pagerank, it is the case that
   *            p + pr(r, alpha) = pr(s, alpha)
   *
   * @raises ArgumentError if the provided seed_vector is not a column vector.
   */
  std::tuple<SprsMat, SprsMat> approximate_pagerank(stag::LocalGraph* graph,
                                                    SprsMat &seed_vector,
                                                    double alpha,
                                                    double epsilon);
}

#endif //STAG_TEST_CLUSTER_H
