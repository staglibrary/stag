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
   * @param target_volume - the approximate volume of the cluster you would like to find
   * @return a vector containing the indices of vectors considered to be in the
   *         same cluster as the seed_vertex.
   */
  std::vector<stag_int> local_cluster(stag::LocalGraph* graph, stag_int seed_vertex, double target_volume);

  /**
   * The ACL local clustering algorithm. Given a graph and starting vertex,
   * returns a cluster close to the starting vertex, constructed in a local way.
   *
   * The locality parameter is passed as the alpha parameter in the personalised
   * pagerank calculation.
   *
   * [ACL] Andersen, Reid, Fan Chung, and Kevin Lang.
   * "Local graph partitioning using pagerank vectors." 2006
   * 47th Annual IEEE Symposium on Foundations of Computer Science (FOCS'06). IEEE, 2006.
   *
   * @param graph - a graph object implementing the LocalGraph interface
   * @param seed_vertex - the starting vertex in the graph
   * @param locality - a value in [0, 1] indicating how 'local' the cluster should
   *                   be. A value of '1' will return the return only the seed vertex
   *                   and a value of '0' will explore the whole graph.
   * @param error (optional) - the acceptable error in the calculation of the approximate
   *                           pagerank. Default 0.001.
   * @return a vector containing the indices of vectors considered to be in the
   *         same cluster as the seed_vertex.
   */
  std::vector<stag_int> local_cluster_acl(stag::LocalGraph* graph, stag_int seed_vertex, double locality);
  std::vector<stag_int> local_cluster_acl(stag::LocalGraph* graph, stag_int seed_vertex, double locality, double error);

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

  /**
   * Find the sweep set of the given vector with the minimum conductance.
   *
   * First, sort the vector, and then let
   *     S_i = {v_j : j <= i}
   * and return the set of original indices corresponding to
   *     argmin_i conducance(S_i)
   *
   * This method is expected to be run on vectors whose support is much less
   * than the total size of the graph. If the total volume of the support of vec
   * is larger than half of the volume of the total graph, then this method may
   * return unexpected results.
   *
   * Note that the caller is responsible for any required normalisation of the
   * input vector. In particular, this method does not normalise the vector by
   * the node degrees.
   *
   * @param graph
   * @param vec
   * @return a vector containing the indices of vec which give the minimum
   *         conductance in the given graph.
   */
  std::vector<stag_int> sweep_set_conductance(stag::LocalGraph* graph,
                                              SprsMat& vec);
}

#endif //STAG_TEST_CLUSTER_H
