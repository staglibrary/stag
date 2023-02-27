//
// Methods for generating graphs from random models.
//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file random.h
 * \brief Methods for generating graphs from random graph models
 */

#ifndef STAG_TEST_RANDOM_H
#define STAG_TEST_RANDOM_H

#include "graph.h"

namespace stag {
  /**
   * Generate a graph from the symmetric stochastic block model.
   *
   * Generates a graph with \f$n\f$ vertices, divided into \f$k\f$ evenly-sized
   * clusters.
   * For each pair of vertices \f$u\f$ and \f$v\f$, the probability of including
   * the edge \f$\{u, v\}\f$ in the graph is
   *  - \f$p\f$ if \f$u\f$ and \f$v\f$ are in the same cluster, and
   *  - \f$q\f$ otherwise.
   *
   * For large enough values of
   * \f$n\f$, this method samples from an approximate stochastic block model by
   * default which significantly speeds up the execution time. To sample
   * exactly from the stochastic block model, pass the optional 'exact'
   * parameter to the method.
   *
   * The approximate sampling method has running time \f$O(k^2 + \mathrm{nnz})\f$
   * where \f$\mathrm{nnz}\f$ is the number of non-zeros in the generated
   * graph's adjacency matrix,
   * and the exact
   * method has running time \f$O(n^2)\f$.
   *
   * @param n the number of vertices in the graph
   * @param k the number of clusters; vertices are split evenly between clusters
   * @param p the probability of including each edge inside a cluster
   * @param q the probability of including each edge between two clusters
   * @param exact (optional) whether to use the exact probability distribution. Default: false.
   * @return the randomly generated graph
   */
  Graph sbm(stag_int n, stag_int k, double p, double q, bool exact);

  /**
   * @overload
   */
  Graph sbm(stag_int n, stag_int k, double p, double q);

  /**
   * Generate a graph from the general stochastic block model.
   *
   * The `cluster_sizes` vector specifies the number of vertices in each
   * generated cluster.
   * Let \f$k\f$ be the length of the cluster_sizes vector.
   *
   * Then, `probabilities` should be a \f$k \times k\f$ matrix which specifies
   * the edge probability between every pair of vertices.
   * That is, for each pair of vertices \f$u\f$ and \f$v\f$, the probability of
   * including the edge \f$\{u, v\}\f$ in the graph is \f$P(u, v)\f$, where
   * \f$P\f$ is the `probabilities` matrix.
   *
   * The approximate sampling method has running time \f$O(k^2 + \mathrm{nnz})\f$
   * where \f$\mathrm{nnz}\f$ is the number of non-zeros in the generated
   * graph's adjacency matrix,
   * and the exact
   * method has running time \f$O(n^2)\f$.
   *
   * \par Example
   *
   * \code{cpp}
   * #include <stag/graph.h>
   * #include <stag/random.h>
   *
   * int main() {
   *   std::vector<stag_int> cluster_sizes = {100, 20, 10};
   *   DenseMat prob_mat {{0.4, 0.1, 0.1}, {0.1, 0.7, 0}, {0.1, 0, 1}};
   *   stag::Graph myGraph = stag::general_sbm(cluster_sizes, prob_mat);
   *   std::cout << *myGraph.adjacency() << std::endl;
   *   return 0;
   * }
   * \endcode
   *
   * @param cluster_sizes a vector of length \f$k\f$ with the number of vertices
   *                      in each cluster.
   * @param probabilities a \f$k \times k\f$ matrix with the inter-cluster
   *                      probabilities.
   * @param exact (optional) whether to use the exact probability distribution. Default: false.
   * @return the randomly generated graph
   */
  Graph general_sbm(std::vector<stag_int>& cluster_sizes,
                    DenseMat& probabilities, bool exact);

  /**
   * @overload
   */
  Graph general_sbm(std::vector<stag_int>& cluster_sizes,
                    DenseMat& probabilities);

  /**
   * Generate a graph from the Erdos-Renyi model.
   *
   * Generates a graph with \f$n\f$ vertices. For each pair of vertices \f$u\f$ and
   * \f$v\f$, the edge \f$\{u, v\}\f$ is included in the graph with probability \f$p\f$.
   *
   * For large values of n, this method will use an approximate version of the
   * random model with running time \f$O(\mathrm{nnz})\f$ where
   * \f$\mathrm{nnz}\f$ is the number of edges in the generated graph.
   *
   * If the 'exact' parameter is true, then the true Erdos-Renyi distribution
   * will be used, with running time \f$O(n^2)\f$.
   *
   * @param n the number of vertices in the graph
   * @param p the probability of including each edge
   * @param exact (optional) whether to sample from the exact model. Default: false.
   * @return the randomly generated graph
   */
  Graph erdos_renyi(stag_int n, double p, bool exact);

  /**
   * \overload
   */
  Graph erdos_renyi(stag_int n, double p);

  /**
   * Construct a vector with the ground truth labels for a graph drawn from the
   * symmetric stochastic block model.
   *
   * \par Example
   *
   * \code{cpp}
   * #include <stag/graph.h>
   * #include <stag/random.h>
   *
   * int main() {
   *   stag_int n = 6;
   *   stag_int k = 3;
   *   stag::Graph myGraph = stag::sbm(n, k, 0.8, 0.1);
   *
   *   std::vector<stag_int> gt_labels = stag::sbm_gt_labels(n, k);
   *
   *   // gt_labels is the vector {0, 0, 1, 1, 2, 2}.
   *
   *   return 0;
   * }
   * \endcode
   *
   *
   * @param n the number of vertices in the graph
   * @param k the number of clusters
   * @return a vector containing the ground truth labels for the vertices in the
   *         graph.
   */
  std::vector<stag_int> sbm_gt_labels(stag_int n, stag_int k);

  /**
   * Construct a vector with the ground truth labels for a graph drawn from the
   * general stochastic block model.
   *
   * \par Example
   *
   * \code{cpp}
   * #include <stag/graph.h>
   * #include <stag/random.h>
   *
   * int main() {
   *   std::vector<stag_int> cluster_sizes = {4, 2};
   *   DenseMat prob_mat {{0.4, 0.1}, {0.1, 0.7}};
   *   stag::Graph myGraph = stag::general_sbm(cluster_sizes, prob_mat);
   *
   *   std::vector<stag_int> gt_labels = stag::general_sbm_gt_labels(cluster_sizes);
   *
   *   // gt_labels is the vector {0, 0, 0, 0, 1, 1}.
   *
   *   return 0;
   * }
   * \endcode
   *
   * @param cluster_sizes a vector of length \f$k\f$ with the number of vertices
   *                      in each cluster.
   * @return a vector containing the ground truth labels for the vertices in the
   *         graph.
   */
  std::vector<stag_int> general_sbm_gt_labels(std::vector<stag_int>& cluster_sizes);
}

#endif //STAG_TEST_RANDOM_H
