/**
 * Methods for generating graphs from random models.
 *
 * Copyright 2022 Peter Macgregor
 */
#ifndef STAG_TEST_RANDOM_H
#define STAG_TEST_RANDOM_H

#include "graph.h"

namespace stag {
  /**
   * Generate a graph from the symmetric stochastic block model.
   *
   * Every cluster has the same number of vertices.
   *
   * @param n - the number of vertices in the graph
   * @param k - the number of clusters; vertices are split evenly between clusters
   * @param p - the probability of including an edge inside a cluster
   * @param q - the probability of including an edge between two clusters
   * @return the randomly generated graph
   */
  Graph sbm(int n, int k, double p, double q);

  /**
   * Generate a graph from the Erdos-Renyi model.
   *
   * @param n - the number of vertices in the graph
   * @param p - the probability of including each edge
   * @return the randomly generated graph
   */
  Graph erdos_renyi(int n, double p);
}

#endif //STAG_TEST_RANDOM_H
