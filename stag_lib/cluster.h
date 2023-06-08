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
   * Find the Cheeger cut in a graph.
   *
   * Let \f$G = (V, E)\f$ be a graph and \f$\mathcal{L}\f$ be its normalised Laplacian
   * matrix with eigenvalues \f$0 = \lambda_1 \leq \lambda_2 \leq \ldots \leq \lambda_n\f$.
   * Then, Cheeger's inequality states that
   *
   * \f[
   *   \frac{\lambda_2}{2} \leq \Phi_G \leq \sqrt{2 \lambda_2},
   * \f]
   *
   * where
   *
   * \f[
   *    \Phi_G = \min_{S \subset V} \phi(S)
   * \f]
   *
   * is the conductance of \f$G\f$. The proof of Cheeger's inequality is
   * constructive: by computing the eigenvector corresponding to \f$\lambda_2\f$,
   * and performing the sweep set operation, we are able to find a set \f$S\f$
   * with conductance close to the optimal. The partition returned by this
   * algorithm is called the 'Cheeger cut' of the graph.
   *
   * @param graph the graph object to be partitioned
   * @return A vector giving the cluster membership for each vertex in the graph.
   *         Each entry in the vector is either \f$0\f$ or \f$1\f$ to indicate
   *         which side of the cut the vertex belongs to.
   */
  std::vector<stag_int> cheeger_cut(stag::Graph* graph);

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
   * When the provided graph is a stag::LocalGraph, the volume of the support of
   * the provided vector should be less than half the total volume of the graph.
   * The method does not (and cannot) check this condition.
   *
   * When the provided graph is a stag::Graph, there is no restriction on the
   * volume of the support of the provided vector.
   *
   * Note that the caller is responsible for any required normalisation of the
   * input vector. In particular, this method does not normalise the vector by
   * the node degrees.
   *
   * @param graph a stag::LocalGraph or stag::Graph object
   * @param vec the vector to sweep over
   * @return a vector containing the indices of vec which give the minimum
   *         conductance in the given graph
   */
  std::vector<stag_int> sweep_set_conductance(stag::LocalGraph* graph,
                                              SprsMat& vec);

  /**
   * @overload
   */
  std::vector<stag_int> sweep_set_conductance(stag::Graph* graph,
                                              SprsMat& vec);

  /**
   * Return the vertex indices of every vertex in the same connected
   * component as the specified vertex.
   *
   * The running time of this method is proportional to the size of the returned
   * connected component.
   *
   * The returned vector is not sorted.
   *
   * @param g a stag::LocalGraph instance
   * @param v a vertex of the graph
   * @return a vector containing the vertex ids of every vertex in the
   *         connected component corresponding to v
   */
  std::vector<stag_int> connected_component(stag::LocalGraph* g, stag_int v);

  /**
   * Return a vector of the connected components in the specified graph.
   *
   * @param g a stag::Graph instance
   * @return a vector containing the connected components of the graph
   */
  std::vector<std::vector<stag_int>> connected_components(stag::Graph* g);

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

  /**
   * Compute the Mutual Information between two label vectors.
   *
   * @param gt_labels the ground truth labels for the dataset
   * @param labels the candidate labels whose MI should be calculated
   * @return the MI between the two labels vectors
   */
  double mutual_information(std::vector<stag_int>& gt_labels,
                            std::vector<stag_int>& labels);

  /**
   * Compute the Normalised Mutual Information between two label vectors.
   *
   * @param gt_labels the ground truth labels for the dataset
   * @param labels the candidate labels whose NMI should be calculated
   * @return the NMI between the two labels vectors
   *
   * \par References
   * Vinh, Epps, and Bailey, (2009).
   * Information theoretic measures for clusterings comparison.
   * 26th Annual International Conference on Machine Learning (ICML ‘09).
   */
  double normalised_mutual_information(std::vector<stag_int>& gt_labels,
                                       std::vector<stag_int>& labels);

  /**
   * Compute the conductance of the given cluster in a graph.
   *
   * Given a graph \f$G = (V, E)\f$, the conductance of \f$S \subseteq V\f$
   * is defined to be
   *
   * \f[
   *    \phi(S) = \frac{w(S, V \setminus S)}{\mathrm{vol}(S)},
   * \f]
   *
   * where \f$\mathrm{vol}(S) = \sum_{v \in S} \mathrm{deg}(v)\f$ is the volume
   * of \f$S\f$ and \f$w(S, V \setminus S)\f$ is the total weight of edges crossing
   * the cut between \f$S\f$ and \f$V \setminus S\f$.
   *
   * @param graph a stag::LocalGraph object representing \f$G\f$.
   * @param cluster a vector of node IDs in \f$S\f$.
   * @return the conductance \f$\phi_G(S)\f$.
   */
  double conductance(stag::LocalGraph* graph,
                     std::vector<stag_int>& cluster);

  /**
   * Compute the symmetric difference of two sets of integers.
   *
   * Given sets \f$S\f$ and \f$T\f$, the symmetric difference
   * \f$S \triangle T\f$ is defined to be
   *
   * \f[
   *    S \triangle T = \{S \setminus T\} \cup \{T \setminus S\}.
   * \f]
   *
   * Although \f$S\f$ and \f$T\f$ are provided as vectors, they are treated
   * as sets and any duplicates will be ignored.
   *
   * @param S a vector containing the first set of integers
   * @param T a vector containina the second set of integers
   * @return a vector containing the integers in the symmetric difference of S
   *         and T.
   */
  std::vector<stag_int> symmetric_difference(std::vector<stag_int>& S,
                                             std::vector<stag_int>& T);
}

#endif //STAG_TEST_CLUSTER_H
