//
// Graph clustering algorithms based on spectral methods.
//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file cluster.h
 * \brief Algorithms for finding clusters in graphs.
 *
 * The two key clustering methods provided by this module are stag::spectral_cluster and stag::local_cluster.
 */

#ifndef STAG_TEST_CLUSTER_H
#define STAG_TEST_CLUSTER_H

#include <vector>

#include <graph.h>

namespace stag {

  /**
   * Spectral clustering algorithm.
   *
   * This is a simple graph clustering method, which provides a clustering of the entire graph.
   * To use spectral clustering, simply pass a `stag::Graph` object
   * and the number of clusters you would like to find.
   *
   * \code{.cpp}
   *     #include <iostream>
   *     #include <stag/graph.h>
   *     #include <stag/cluster.h>
   *
   *     int main() {
   *       stag::Graph myGraph = stag::barbell_graph(10);
   *       std::vector<stag_int> clusters = stag::spectral_cluster(&myGraph, 2);
   *
   *       for (auto c : clusters) {
   *         std::cout << c << ", ";
   *       }
   *       std::cout << std::endl;
   *
   *       return 0;
   *     }
   * \endcode
   * 
   * The spectral clustering algorithm has the following steps.
   *   - Compute the \f$k\f$ smallest eigenvectors of the normalised Laplacian matrix.
   *   - Embed the vertices into \f$\mathbb{R}^k\f$ according to the eigenvectors.
   *   - Cluster the vertices into \f$k\f$ clusters using a \f$k\f$-means clustering algorithm.
   *
   * @param graph the graph object to be clustered
   * @param k the number of clusters to find. Should be less than \f$n/2\f$.
   * @return a vector giving the cluster membership for each vertex in the graph
   *
   * \par References
   * A. Ng, M. Jordan, Y. Weiss.
   * On spectral clustering: Analysis and an algorithm. NeurIPS'01
   */
  std::vector<stag_int> spectral_cluster(stag::Graph* graph, stag_int k);


  /**
   * Local clustering algorithm based on personalised Pagerank.
   *
   * Given a graph and starting vertex, return a cluster which is close to the
   * starting vertex.
   *
   * This method uses the ACL local clustering algorithm.
   *
   * @param graph a graph object implementing the LocalGraph interface
   * @param seed_vertex the starting vertex in the graph
   * @param target_volume the approximate volume of the cluster you would like to find
   * @return a vector containing the indices of vectors considered to be in the
   *         same cluster as the seed_vertex.
   *
   * \par References
   * R. Andersen, F. Chung, K. Lang.
   * Local graph partitioning using pagerank vectors. FOCS'06
   */
  std::vector<stag_int> local_cluster(stag::LocalGraph* graph, stag_int seed_vertex, double target_volume);

  /**
   * The ACL local clustering algorithm. Given a graph and starting vertex,
   * return a cluster close to the starting vertex, constructed in a local way.
   *
   * The locality parameter is passed as the alpha parameter in the personalised
   * Pagerank calculation.
   *
   * @param graph a graph object implementing the LocalGraph interface
   * @param seed_vertex the starting vertex in the graph
   * @param locality a value in \f$[0, 1]\f$ indicating how 'local' the cluster should
   *                 be. A value of \f$1\f$ will return only the seed vertex,
   *                 and a value of \f$0\f$ will explore the whole graph.
   * @param error (optional) - the acceptable error in the calculation of the approximate
   *                           pagerank. Default \f$0.001\f$.
   * @return a vector containing the indices of vectors considered to be in the
   *         same cluster as the seed_vertex.
   *
   * \par References
   * R. Andersen, F. Chung, K. Lang.
   * Local graph partitioning using pagerank vectors. FOCS'06
   */
  std::vector<stag_int> local_cluster_acl(stag::LocalGraph* graph, stag_int seed_vertex, double locality, double error);

  /**
   * \overload
   */
  std::vector<stag_int> local_cluster_acl(stag::LocalGraph* graph, stag_int seed_vertex, double locality);

  /**
   * Compute the approximate Pagerank vector.
   *
   * The parameters seed_vector, alpha, and epsilon are used as described in the ACL paper.
   *
   * Note that the dimension of the returned vectors may not match the correct
   * number of vertices in the graph provided since the approximate
   * Pagerank is computed locally.
   *
   * @param graph a stag::LocalGraph object
   * @param seed_vector the seed vector of the personalised Pagerank
   * @param alpha the locality parameter of the personalised Pagerank
   * @param epsilon the error parameter of the personalised Pagerank
   * @return A tuple of sparse column vectors corresponding to
   *          - p: the approximate Pagerank vector
   *          - r: the residual vector
   *
   *         By the definition of approximate Pagerank, it holds that
   *            p + ppr(r, alpha) = ppr(s, alpha).
   *
   * @throws std::invalid_argument if the provided seed_vector is not a column vector.
   *
   * \par References
   * R. Andersen, F. Chung, K. Lang.
   * Local graph partitioning using pagerank vectors. FOCS'06
   */
  std::tuple<SprsMat, SprsMat> approximate_pagerank(stag::LocalGraph* graph,
                                                    SprsMat &seed_vector,
                                                    double alpha,
                                                    double epsilon);

  /**
   * Find the sweep set of the given vector with the minimum conductance.
   *
   * First, sort the vector such that \f$v_1<= \ldots <= v_n\f$. Then let
   *
   * \f[
   *     S_i = \{v_j : j <= i\}
   * \f]
   *
   * and return the set of original indices corresponding to
   *
   * \f[
   *     \mathrm{argmin}_i \phi(S_i)
   * \f]
   *
   * where \f$\phi(S)\f$ is the conductance of \f$S\f$.
   *
   * This method is expected to be run on vectors whose support is much less
   * than the total size of the graph. If the total volume of the support of vec
   * is larger than half of the volume of an entire graph, then this method may
   * return unexpected results.
   *
   * Note that the caller is responsible for any required normalisation of the
   * input vector. In particular, this method does not normalise the vector by
   * the node degrees.
   *
   * @param graph a stag::LocalGraph object
   * @param vec the vector to sweep over
   * @return a vector containing the indices of vec which give the minimum
   *         conductance in the given graph
   */
  std::vector<stag_int> sweep_set_conductance(stag::LocalGraph* graph,
                                              SprsMat& vec);

  /**
   * Compute the Adjusted Rand Index between two label vectors.
   *
   * @param gt_labels the ground truth labels for the dataset
   * @param labels the candidate labels whose ARI should be calculated
   * @return the ARI between the two labels vectors
   *
   * \par References
   * W. M. Rand.
   * Objective criteria for the evaluation of clustering methods.
   * Journal of the American Statistical Association. 66 (336): 846–850. 1971.
   */
  double adjusted_rand_index(std::vector<stag_int>& gt_labels,
                             std::vector<stag_int>& labels);
}

#endif //STAG_TEST_CLUSTER_H
