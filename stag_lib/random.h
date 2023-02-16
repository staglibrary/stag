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
   * n, this method samples from an approximate stochastic block model by
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
}

#endif //STAG_TEST_RANDOM_H
