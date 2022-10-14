/**
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <random>
#include <iostream>

#include "random.h"
#include "graph.h"

/**
 * Sample edges between SBM clusters by directly iterating through each
 * edge and 'tossing a coin'. This technique should be used for 'large' values
 * of p.
 *
 * @param cluster_idx
 * @param verticesInCluster
 * @param p
 * @return a list of the sampled edges
 */
std::vector<EdgeTriplet> sample_edges_directly(stag_int cluster_idx,
                                               stag_int other_cluster_idx,
                                               stag_int verticesInCluster,
                                               double p) {
  // Prepare the random number generator
  std::random_device dev;
  std::mt19937 prng(dev());
  std::bernoulli_distribution sampleDist(p);

  // Store the sampled edges
  std::vector<EdgeTriplet> sampledEdges;

  for (stag_int i = cluster_idx * verticesInCluster; i < (cluster_idx + 1) * verticesInCluster; i++) {
    for (stag_int j = other_cluster_idx * verticesInCluster; j < (other_cluster_idx + 1) * verticesInCluster; j++) {
      // If we are in the same cluster, then don't double sample
      if (cluster_idx == other_cluster_idx && j <= i) continue;

      // Toss a coin
      if (sampleDist(prng)) {
        sampledEdges.emplace_back(i, j, 1);
        sampledEdges.emplace_back(j, i, 1);
      }
    }
  }

  return sampledEdges;
}

/**
 * Sample edges between SBM clusters using the 'binomial trick'. This technique
 * should be used for 'small' values of p.
 *
 * @param cluster_idx
 * @param other_cluster_idx
 * @param verticesInCluster
 * @param p
 * @return a list of the sampled edges
 */
std::vector<EdgeTriplet> sample_edges_binomial(stag_int cluster_idx,
                                               stag_int other_cluster_idx,
                                               stag_int verticesInCluster,
                                               double p) {
  // Validate the function inputs
  assert(cluster_idx >= 0);
  assert(verticesInCluster >= 0);
  assert(0 <= p <= 1);

  // Prepare the random number generator. We will approximate the binomial
  // distribution with the normal distribution
  std::random_device dev;
  std::mt19937 prng(dev());
  std::normal_distribution<double> numEdgesDist(p * ((double) verticesInCluster) * ((double) verticesInCluster),
                                                (1 - p) * p * ((double) verticesInCluster) * ((double) verticesInCluster));
  std::uniform_int_distribution<stag_int> vertexDist(0, verticesInCluster - 1);

  // Decide how many edges to sample based on the 'binomial' distribution
  stag_int max_possible_edges = verticesInCluster * (verticesInCluster - 1);

  auto raw_sample = (stag_int) floor(numEdgesDist(prng));
  stag_int numEdges = std::max((stag_int) 0, std::min(max_possible_edges, raw_sample));

  // Store the sampled edges
  std::vector<EdgeTriplet> sampledEdges(numEdges);

  // Sample the specific vertices
  stag_int randU = 0;
  stag_int randV = 0;
  for (stag_int i = 0; i < numEdges; i++) {
    // Choose two random vertices in the cluster
    randU = 0;
    randV = 0;
    while (randU == randV) {
      // Ignore the edge if u and v are identical
      randU = verticesInCluster * cluster_idx + vertexDist(prng);
      randV = verticesInCluster * other_cluster_idx + vertexDist(prng);
    }

    // Add this vertex to the sampled edges
    sampledEdges.emplace_back(randU, randV, 1);
    sampledEdges.emplace_back(randV, randU, 1);
  }

  return sampledEdges;
}

stag::Graph stag::sbm(stag_int n, stag_int k, double p, double q) {
  return stag::sbm(n, k, p, q, false);
}

stag::Graph stag::sbm(stag_int n, stag_int k, double p, double q, bool exact) {
  // All arguments must be positive
  assert(n > 0);
  assert(k > 0);
  assert(p >= 0);
  assert(q >= 0);

  // Strictly speaking, we will generate a graph with k * floor(n/k) vertices
  stag_int verticesPerCluster = std::floor(n/k);
  stag_int totalVertices = k * verticesPerCluster;

  // We will build the adjacency matrix as we go along.
  std::vector<EdgeTriplet> allEdges;

  // Iterate through the clusters
  for (stag_int cluster_idx = 0; cluster_idx < k; cluster_idx++) {
    // First, sample the edges inside each cluster
    std::vector<EdgeTriplet> sampledEdges;
    if (verticesPerCluster >= 100 && p < 0.5 && !exact) {
      // For small values of p, use the 'binomial trick' for sampling
      sampledEdges = sample_edges_binomial(cluster_idx, cluster_idx, verticesPerCluster, p);
    } else {
      // For large values of p, we just iterate over every pair of vertices
      sampledEdges = sample_edges_directly(cluster_idx, cluster_idx, verticesPerCluster, p);
    }

    // Add the new edges to the full list
    allEdges.insert(allEdges.end(), sampledEdges.begin(), sampledEdges.end());

    // Now, sample edges to the other clusters
    for (stag_int other_cluster_idx = cluster_idx + 1; other_cluster_idx < k; other_cluster_idx++) {
      if (verticesPerCluster >= 100 && q < 0.5 && !exact) {
        // For small values of p, use the 'binomial trick' for sampling
        sampledEdges = sample_edges_binomial(cluster_idx, other_cluster_idx, verticesPerCluster, q);
      } else {
        // For large values of p, we just iterate over every pair of vertices
        sampledEdges = sample_edges_directly(cluster_idx, other_cluster_idx, verticesPerCluster, q);
      }

      // Add the new edges to the full list
      allEdges.insert(allEdges.end(), sampledEdges.begin(), sampledEdges.end());
    }
  }

  // Finally, construct the graph
  SprsMat adj_mat(totalVertices, totalVertices);
  adj_mat.setFromTriplets(allEdges.begin(), allEdges.end());
  return stag::Graph(adj_mat);
}

stag::Graph stag::erdos_renyi(stag_int n, double p) {
  return stag::erdos_renyi(n, p, false);
}

stag::Graph stag::erdos_renyi(stag_int n, double p, bool exact) {
  return stag::sbm(n, 1, p, 0, exact);
}
