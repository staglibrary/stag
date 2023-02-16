/**
 * Tests for the clustering methods in cluster.h.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <gtest/gtest.h>
#include <cluster.h>
#include <graph.h>
#include <random.h>
#include <utility.h>

// Define some helper test assertions.
#define EXPECT_FLOATS_NEARLY_EQ(expected, actual, thresh) \
        EXPECT_EQ(expected.size(), actual.size()) << "Array sizes differ.";\
        for (size_t idx = 0; idx < std::min(expected.size(), actual.size()); ++idx) \
        { \
            EXPECT_NEAR(expected[idx], actual[idx], thresh) << "at index: " << idx;\
        }

TEST(ClusterTest, SpectralCluster) {
  // Construct a test graph from the SBM
  stag_int n = 1000;
  stag_int k = 5;
  stag::Graph testGraph = stag::sbm(n, k, 0.6, 0.1);

  // Find the clusters
  auto clusters = stag::spectral_cluster(&testGraph, k);

  // There should be approximately the same number of each cluster
  stag_int c1 = 0;
  stag_int c2 = 0;

  for (auto c : clusters) {
    if (c == 1) c1++;
    if (c == 2) c2++;
  }
  EXPECT_NEAR(c1 / c2, 1, 0.01);
}

TEST(ClusterTest, SCArguments) {
  // Construct a test graph
  stag_int n = 100;
  stag::Graph testGraph = stag::erdos_renyi(n, 0.1);

  // Passing an invalid number of clusters to the spectral clustering
  // algorithm should return an error.
  stag_int k = -1;
  EXPECT_THROW(stag::spectral_cluster(&testGraph, k), std::invalid_argument);

  k = 0;
  EXPECT_THROW(stag::spectral_cluster(&testGraph, k), std::invalid_argument);

  k = (n / 2) + 1;
  EXPECT_THROW(stag::spectral_cluster(&testGraph, k), std::invalid_argument);

  k = n + 1;
  EXPECT_THROW(stag::spectral_cluster(&testGraph, k), std::invalid_argument);
}

TEST(ClusterTest, SpectralClusterSparse) {
  // Construct a very sparse test graph from the SBM
  stag_int n = 10000;
  stag_int k = 5;
  stag::Graph testGraph = stag::sbm(n, k, 0.06, 0.01);

  // Find the clusters
  auto clusters = stag::spectral_cluster(&testGraph, k);

  // There should be approximately the same number of each cluster
  stag_int c1 = 0;
  stag_int c2 = 0;

  for (auto c : clusters) {
    if (c == 1) c1++;
    if (c == 2) c2++;
  }
  EXPECT_NEAR(c1, c2, 0.8 * c1);
  EXPECT_NEAR(c1, n / k, 0.8 * c1);
  EXPECT_NEAR(c2, n / k, 0.8 * c2);
}

TEST(ClusterTest, SpectralClusterDisconnected) {
  // Construct a sparse disconnected graph from the SBM
  stag_int n = 1000;
  stag_int k = 2;
  stag::Graph testGraph = stag::sbm(n, k, 0.5, 0);

  // Find the clusters
  auto clusters = stag::spectral_cluster(&testGraph, k);

  // There should be exactly the same number of each cluster
  stag_int c0 = 0;
  stag_int c1 = 0;

  for (auto c : clusters) {
    if (c == 0) c0++;
    if (c == 1) c1++;
  }
  EXPECT_EQ(c0, n / k);
  EXPECT_EQ(c1, n / k);
}

TEST(ClusterTest, ApproxPageRank) {
  // Construct a test graph
  std::vector<stag_int> rowStarts = {0, 2, 4, 6, 8};
  std::vector<stag_int> colIndices = {1, 2, 0, 3, 0, 3, 1, 2};
  std::vector<double> values = {1./2, 1./2, 1./2, 1./2, 1./2, 1./2, 1./2, 1./2};
  stag::Graph testGraph(rowStarts, colIndices, values);

  // Run the approximate pagerank method
  SprsMat seed(1, 1);
  seed.coeffRef(0, 0) = 1;
  std::tuple<SprsMat, SprsMat> apr = stag::approximate_pagerank(&testGraph, seed, 1./3, 1./8);

  // Check that the returned vectors are the expected ones.
  std::vector<double> expected_p = {41./81, 2./27, 2./27};
  std::vector<double> expected_r = {5./81, 2./27 + 5./162, 2./27 + 5./162, 2./27};
  EXPECT_FLOATS_NEARLY_EQ(expected_p, stag::sprsMatValues(&std::get<0>(apr)), 0.00001);
  EXPECT_FLOATS_NEARLY_EQ(expected_r, stag::sprsMatValues(&std::get<1>(apr)), 0.00001);
}

TEST(ClusterTest, ACL) {
  // Construct a test barbell graph
  stag::Graph testGraph = stag::barbell_graph(10);

  // Run the acl clustering method
  std::vector<stag_int> cluster = stag::local_cluster_acl(&testGraph, 0, 0.8,
                                                          0.0001);
  std::stable_sort(cluster.begin(), cluster.end());

  // Check that we found one of the clusters.
  std::vector<stag_int> expected_cluster = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  EXPECT_EQ(cluster, expected_cluster);
}

TEST(ClusterTest, pagerankStart) {
  // Ensure that the approximate pagerank method works for a starting vertex
  // containing multiple starting vertices.
  stag::Graph testGraph = stag::star_graph(10);

  // Create the starting vector
  SprsMat seed_vector(10, 1);
  seed_vector.coeffRef(0, 0) = 0.5;
  seed_vector.coeffRef(4, 0) = 0.5;

  // Run the approximate pagerank method.
  // We set epsilon to be 0.1 so that the central vertex in the star graph
  // does not get added to the 'push queue' at the start of the APR algorithm.
  auto result = stag::approximate_pagerank(&testGraph,
                                           seed_vector,
                                           0.1,
                                           0.1);

  // Make sure that the starting vertex 4 has a non-zero
  // value on the resulting p vector.
  ASSERT_GT(std::get<0>(result).coeff(4, 0), 0);
}

TEST(ClusterTest, localSBM) {
  // Construct an SBM and check the clustering method
  // Note that there is some small probability that this will fail.
  stag::Graph testGraph = stag::sbm(2000, 4, 0.9, 0.01);

  // Find a cluster using the default local clustering algorithm
  std::vector<stag_int> cluster = stag::local_cluster(&testGraph, 0, testGraph.total_volume() / 4);
  std::stable_sort(cluster.begin(), cluster.end());

  // Check the number of returned vertices
  EXPECT_GE(cluster.size(), 100);
  EXPECT_LE(cluster.size(), 1000);

  // Let's say that 50% of the found cluster should lie inside the first cluster
  // in the SBM graph
  stag_int inside = 0;
  for (auto v : cluster) {
    if (v < 500) inside++;
  }
  EXPECT_GE(inside / cluster.size(), 0.5);
}

TEST(ClusterTest, LocalClusterArguments) {
  // Create a test graph
  stag_int n = 100;
  stag::Graph testGraph = stag::erdos_renyi(n, 0.1);

  // Check vertex is valid
  stag_int v = -1;
  double target_vol = 100;
  EXPECT_THROW(stag::local_cluster(&testGraph, v, target_vol), std::invalid_argument);

  v = n;
  EXPECT_THROW(stag::local_cluster(&testGraph, v, target_vol), std::invalid_argument);

  // Check volume is valid
  v = 0;
  target_vol = 0;
  EXPECT_THROW(stag::local_cluster(&testGraph, v, target_vol), std::invalid_argument);

  target_vol = -1;
  EXPECT_THROW(stag::local_cluster(&testGraph, v, target_vol), std::invalid_argument);
}

TEST(ClusterTest, sweepSet) {
  // Construct a test graph
  stag::Graph testGraph = stag::barbell_graph(4);

  // Construct the vector on which to sweep
  SprsMat vec(8, 1);
  vec.coeffRef(0, 0) = 0.1;
  vec.coeffRef(1, 0) = 0.25;
  vec.coeffRef(2, 0) = 0.2;
  vec.coeffRef(3, 0) = 0.15;
  vec.coeffRef(4, 0) = 0.05;
  vec.makeCompressed();

  // Compute the sweep set
  std::vector<stag_int> sweep_set = stag::sweep_set_conductance(&testGraph, vec);
  std::stable_sort(sweep_set.begin(), sweep_set.end());

  // Test the result
  std::vector<stag_int> expected_sweep = {0, 1, 2, 3};
  EXPECT_EQ(sweep_set, expected_sweep);
}
