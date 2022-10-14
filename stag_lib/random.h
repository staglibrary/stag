/**
 * Methods for generating graphs from random models.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#ifndef STAG_TEST_RANDOM_H
#define STAG_TEST_RANDOM_H

#include "graph.h"

namespace stag {
  /**
   * Generate a graph from the symmetric stochastic block model.
   *
   * Every cluster has the same number of vertices. For large enough values of
   * n, this method samples from an approximate stochastic block model by
   * default which significantly speeds up the execution time. To sample
   * exactly from the stochastic block model, pass the optional 'exact'
   * parameter to the method.
   *
   * The approximate sampling method has running time O(k^2 + nnz) and the exact
   * method has running time O(n^2) where nnz is the number of non-zeros in the
   * generated graph's adjacency matrix.
   *
   * @param n - the number of vertices in the graph
   * @param k - the number of clusters; vertices are split evenly between clusters
   * @param p - the probability of including an edge inside a cluster
   * @param q - the probability of including an edge between two clusters
   * @param exact (optional) - whether to use the exact probability distribution
   *                           or an approximation.
   * @return the randomly generated graph
   */
  Graph sbm(stag_int n, stag_int k, double p, double q);
  Graph sbm(stag_int n, stag_int k, double p, double q, bool exact);

  /**
   * Generate a graph from the Erdos-Renyi model.
   *
   * For large values of n, this method will use an approximate version of the
   * random model with running time O(nnz) where nnz is the number of edges in
   * the sampled graph.
   *
   * If the 'exact' parameter is true, then the true Erdos-Renyi distribution
   * will be used, with running time O(n^2).
   *
   * @param n - the number of vertices in the graph
   * @param p - the probability of including each edge
   * @param exact (optional) - whether to sample from the exact model
   * @return the randomly generated graph
   */
  Graph erdos_renyi(stag_int n, double p);
  Graph erdos_renyi(stag_int n, double p, bool exact);
}

#endif //STAG_TEST_RANDOM_H
