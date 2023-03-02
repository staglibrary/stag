//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#include <random>
#include <iostream>

#include "random.h"
#include "graph.h"

/**
 * Get an upper estimate on the number of neighbours of each node in a
 * graph generated from the stochastic block model.
 *
 * This will be used to pre-allocate the memory of the adjacency matrix as we
 * construct the graph.
 *
 * @param cluster_sizes a vector of length \f$k\f$ with the number of vertices
 *                      in each cluster.
 * @param probabilities a \f$k \times k\f$ matrix with the inter-cluster
 *                      probabilities.
 * @return an Eigen vector with the neighbour estimates.
 */
Eigen::VectorXi estimate_sbm_neighbours(std::vector<stag_int>& cluster_sizes,
                                        DenseMat probabilities) {
  // Get the total number of vertices in the graph
  auto k = cluster_sizes.size();
  stag_int n = 0;
  for (stag_int s : cluster_sizes) n += s;

  // Create the vector which we'll return
  Eigen::VectorXi neighbours(n);

  // Compute the expected number of neighbours for a node in each cluster
  // This is basically the matrix product of probabilities with cluster_sizes
  Eigen::VectorXi cluster_neighbours(k);
  for (auto i = 0; i < k; i++) {
    stag_int this_cluster_neighbours = 0;
    for (auto j = 0; j < k; j++) {
      this_cluster_neighbours += cluster_sizes.at(j) * probabilities(i, j);
    }

    // Add a safety factor of 2
    cluster_neighbours.coeffRef(i) = 2 * this_cluster_neighbours;
  }

  // Now add the neighbour estimates for each actual node
  stag_int current_idx = 0;
  for (stag_int i = 0; i < k; i++) {
    for (stag_int j = 0; j < cluster_sizes.at(i); j++) {
      neighbours.coeffRef(current_idx) = cluster_neighbours(i);
      current_idx++;
    }
  }

  return neighbours;
}

/**
 * Sample edges between SBM clusters by directly iterating through each
 * edge and 'tossing a coin'. This technique should be used for 'large' values
 * of p.
 *
 * @param adj_mat the adjacency matrix to be updated
 * @param cluster_idx the index of the 'source' cluster
 * @param other_cluster_idx the index of the 'target' cluster
 * @param this_cluster_vertices the number of vertices in the 'source' cluster
 * @param other_cluster_vertices the number of vertices in the 'target' cluster
 * @param this_cluster_start_idx the index of the first vertex in the 'source' cluster
 * @param other_cluster_start_idx the index of the first vertex in the 'target' cluster
 * @param p the probability of including each edge.
 */
void sample_edges_directly(SprsMat& adj_mat,
                           stag_int cluster_idx,
                           stag_int other_cluster_idx,
                           stag_int this_cluster_vertices,
                           stag_int other_cluster_vertices,
                           stag_int this_cluster_start_idx,
                           stag_int other_cluster_start_idx,
                           double p) {
  // Prepare the random number generator
  std::random_device dev;
  std::mt19937 prng(dev());
  std::bernoulli_distribution sampleDist(p);

  for (stag_int i = this_cluster_start_idx;
      i < this_cluster_start_idx + this_cluster_vertices; i++) {
    for (stag_int j = other_cluster_start_idx;
        j < other_cluster_start_idx + other_cluster_vertices; j++) {
      // If we are in the same cluster, then don't double sample
      if (cluster_idx == other_cluster_idx && j <= i) continue;

      // Toss a coin
      if (sampleDist(prng)) {
        adj_mat.insert(i, j) = 1;
        adj_mat.insert(j, i) = 1;
      }
    }
  }
}

/**
 * Sample edges between SBM clusters using the 'binomial trick'. This technique
 * should be used for 'small' values of p.
 *
 * @param adj_mat the adjacency matrix to be updated
 * @param this_cluster_vertices the number of vertices in the 'source' cluster
 * @param other_cluster_vertices the number of vertices in the 'target' cluster
 * @param this_cluster_start_idx the index of the first vertex in the 'source' cluster
 * @param other_cluster_start_idx the index of the first vertex in the 'target' cluster
 * @param p
 */
void sample_edges_binomial(SprsMat& adj_mat,
                           stag_int this_cluster_vertices,
                           stag_int other_cluster_vertices,
                           stag_int this_cluster_start_idx,
                           stag_int other_cluster_start_idx,
                           double p) {
  // Validate the function inputs
  assert(0 <= p <= 1);

  // Get the total number of possible edges and the expected number of edges
  stag_int max_edges = this_cluster_vertices * other_cluster_vertices;
  if (this_cluster_start_idx == other_cluster_start_idx) max_edges /= 2;
  double expected_num_edges = p * ((double) max_edges);

  // Prepare the random number generator. We will approximate the binomial
  // distribution with the normal distribution
  std::random_device dev;
  std::mt19937 prng(dev());
  std::normal_distribution<double> numEdgesDist(expected_num_edges,
                                                sqrt((1 - p) * expected_num_edges));
  std::uniform_int_distribution<stag_int> thisVertexDist(0, this_cluster_vertices - 1);
  std::uniform_int_distribution<stag_int> otherVertexDist(0, other_cluster_vertices - 1);

  // Decide how many edges to sample based on the 'binomial' distribution
  auto raw_sample = (stag_int) floor(numEdgesDist(prng));
  stag_int numEdges = std::max((stag_int) 0, std::min(max_edges, raw_sample));

  // Sample the specific vertices
  stag_int randU = 0;
  stag_int randV = 0;
  for (stag_int i = 0; i < numEdges; i++) {
    // Choose two random vertices in the cluster
    randU = 0;
    randV = 0;
    while (randU == randV) {
      // Ignore the edge if u and v are identical
      randU = this_cluster_start_idx + thisVertexDist(prng);
      randV = other_cluster_start_idx + otherVertexDist(prng);
    }

    // Add this vertex to the adjacency matrix
    adj_mat.coeffRef(randU, randV)++;
    adj_mat.coeffRef(randV, randU)++;
  }
}

stag::Graph stag::sbm(stag_int n, stag_int k, double p, double q) {
  return stag::sbm(n, k, p, q, false);
}

stag::Graph stag::sbm(stag_int n, stag_int k, double p, double q, bool exact) {
  if (n < 1) throw std::invalid_argument("Number of vertices must be at least 1.");
  if (k < 1 || k > n/2) {
    throw std::invalid_argument("Number of clusters must be between 1 and n/2.");
  }
  if (p < 0 || p > 1) {
    throw std::invalid_argument("p must be between 0 and 1.");
  }
  if (q < 0 || q > 1) {
    throw std::invalid_argument("q must be between 0 and 1.");
  }

  // Create the cluster size vector and probabilities matrix
  std::vector<stag_int> cluster_sizes;
  DenseMat probabilities(k, k);
  for (auto i = 0; i < k; i++) {
    cluster_sizes.push_back(floor(((double) n) / ((double) k)));
    probabilities(i, i) = p;

    for (auto j = i + 1; j < k; j++) {
      probabilities(i, j) = q;
      probabilities(j, i) = q;
    }
  }

  return general_sbm(cluster_sizes, probabilities, exact);
}

stag::Graph stag::general_sbm(std::vector<stag_int>& cluster_sizes,
                              DenseMat& probabilities, bool exact) {
  // The number of clusters is the length of the cluster_sizes vector
  stag_int k = cluster_sizes.size();

  // Check that the input parameters make sense.
  for (stag_int size : cluster_sizes) {
    if (size < 1) throw std::invalid_argument("Number of vertices in each cluster must be at least 1.");
  }
  if (probabilities.rows() != k || probabilities.cols() != k) {
    throw std::invalid_argument("Probability matrix must be of size k * k.");
  }
  for (auto i = 0; i < k; i++) {
    for (auto j = 0; j < k; j++) {
      if (probabilities(i, j) < 0 || probabilities(i, j) > 1) {
        throw std::invalid_argument("All probabilities must be between 0 and 1.");
      }
    }
  }

  // Initialise the adjacency matrix
  stag_int n = 0;
  for (auto s : cluster_sizes) n += s;
  SprsMat adj_mat(n, n);

  // Estimate the number of neighbours of each node and reserve
  // memory for the adjacency matrix.
  Eigen::VectorXi neighbour_estimates = estimate_sbm_neighbours(cluster_sizes,
                                                                probabilities);
  adj_mat.reserve(neighbour_estimates);

  // Iterate through the clusters
  stag_int this_cluster_start_idx = 0;
  for (stag_int cluster_idx = 0; cluster_idx < k; cluster_idx++) {
    // Get the number of vertices in the current cluster
    stag_int this_cluster_vertices = cluster_sizes.at(cluster_idx);

    stag_int other_cluster_start_idx = this_cluster_start_idx;
    for (stag_int other_cluster_idx = cluster_idx;
         other_cluster_idx < k; other_cluster_idx++){
      // Get the number of vertices in the other cluster
      stag_int other_cluster_vertices = cluster_sizes.at(other_cluster_idx);

      // Get the sampling probability between the two clusters
      double prob = probabilities(cluster_idx, other_cluster_idx);

      // Sample the edges between this cluster and the other cluster
      if (this_cluster_vertices * other_cluster_vertices >= 10000 &&
            prob < 0.5 && !exact) {
        // For small probabilities, use the 'binomial trick' for sampling
        sample_edges_binomial(adj_mat,
                              this_cluster_vertices,
                              other_cluster_vertices,
                              this_cluster_start_idx,
                              other_cluster_start_idx,
                              prob);
      } else {
        // For large probabilities, we just iterate over every pair of vertices
        sample_edges_directly(adj_mat,
                              cluster_idx,
                              other_cluster_idx,
                              this_cluster_vertices,
                              other_cluster_vertices,
                              this_cluster_start_idx,
                              other_cluster_start_idx,
                              prob);
      }

      // Update the starting index for the other cluster
      other_cluster_start_idx += other_cluster_vertices;
    }

    // Update the starting index for this cluster
    this_cluster_start_idx += this_cluster_vertices;
  }

  // Finally, construct and return the graph.
  adj_mat.makeCompressed();
  return stag::Graph(adj_mat);
}


stag::Graph stag::general_sbm(std::vector<stag_int>& cluster_sizes,
                              DenseMat& probabilities) {
  return stag::general_sbm(cluster_sizes, probabilities, false);
}

stag::Graph stag::erdos_renyi(stag_int n, double p) {
  return stag::erdos_renyi(n, p, false);
}

stag::Graph stag::erdos_renyi(stag_int n, double p, bool exact) {
  return stag::sbm(n, 1, p, 0, exact);
}

std::vector<stag_int> stag::sbm_gt_labels(stag_int n, stag_int k) {
  if (n < 1) throw std::invalid_argument("Number of vertices must be at least 1.");
  if (k < 1 || k > n/2) {
    throw std::invalid_argument("Number of clusters must be between 1 and n/2.");
  }

  // Create the cluster size vector
  std::vector<stag_int> cluster_sizes;
  for (auto i = 0; i < k; i++) {
    cluster_sizes.push_back(floor(((double) n) / ((double) k)));
  }

  return general_sbm_gt_labels(cluster_sizes);
}

std::vector<stag_int> stag::general_sbm_gt_labels(std::vector<stag_int>& cluster_sizes) {
  for (stag_int size : cluster_sizes) {
    if (size < 1) throw std::invalid_argument("Number of vertices in each cluster must be at least 1.");
  }

  std::vector<stag_int> labels;

  stag_int current_cluster = 0;
  for (auto this_size : cluster_sizes) {
    for (auto j = 0; j < this_size; j++) {
      labels.push_back(current_cluster);
    }
    current_cluster++;
  }

  return labels;
}