/**
 * Tests for the methods in the random.h header file. Includes methods for
 * generating random graphs.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <iostream>
#include <gtest/gtest.h>
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

TEST(RandomTest, SBMComplete) {
  // Generate what should be a complete graph using the stochastic block model.
  stag::Graph testGraph = stag::sbm(10, 1, 1, 0);

  stag::Graph referenceGraph = stag::complete_graph(10);

  std::vector<stag_int> refStarts = stag::sprsMatOuterStarts(
      referenceGraph.laplacian());
  std::vector<stag_int> testStarts = stag::sprsMatOuterStarts(
      testGraph.laplacian());
  EXPECT_EQ(refStarts, testStarts);

  std::vector<stag_int> refIndices= stag::sprsMatInnerIndices(
      referenceGraph.laplacian());
  std::vector<stag_int> testIndices = stag::sprsMatInnerIndices(
      testGraph.laplacian());
  EXPECT_EQ(refIndices, testIndices);

  std::vector<double> refValues = stag::sprsMatValues(
      referenceGraph.laplacian());
  std::vector<double> testValues = stag::sprsMatValues(
      testGraph.laplacian());
  EXPECT_FLOATS_NEARLY_EQ(refValues, testValues, 0.000001);
}

TEST(RandomTest, GeneralSBM) {
  std::vector<stag_int> cluster_sizes = {1000, 100, 10};
  DenseMat prob_mat {{0.4, 0.1, 0.1}, {0.1, 0.7, 0}, {0.1, 0, 1}};
  stag::Graph testGraph = stag::general_sbm(cluster_sizes, prob_mat);
  EXPECT_EQ(testGraph.number_of_vertices(), 1110);
}

TEST(RandomTest, GeneralSBMArguments) {
  std::vector<stag_int> cluster_sizes = {1000, 100, 10};
  DenseMat prob_mat {{0.4, 0.1, 0.1}, {0.1, 1.7, 0}, {0.1, 0, 1}};
  EXPECT_THROW(stag::general_sbm(cluster_sizes, prob_mat), std::invalid_argument);

  prob_mat(1, 1) = 0.7;
  prob_mat(2, 0) = -0.1;
  EXPECT_THROW(stag::general_sbm(cluster_sizes, prob_mat), std::invalid_argument);

  prob_mat(2, 0) = 0.1;
  cluster_sizes.at(2) = 0;
  EXPECT_THROW(stag::general_sbm(cluster_sizes, prob_mat), std::invalid_argument);

  cluster_sizes.at(2) = -10;
  EXPECT_THROW(stag::general_sbm(cluster_sizes, prob_mat), std::invalid_argument);
}

TEST(RandomTest, SBMApprox) {
  stag::Graph testGraph = stag::sbm(1000, 2, 0.1, 0.01);
  EXPECT_EQ(testGraph.number_of_vertices(), 1000);
}

TEST(RandomTest, SBMExact) {
  stag::Graph testGraph = stag::sbm(1000, 2, 0.1, 0.01, true);
  EXPECT_EQ(testGraph.number_of_vertices(), 1000);
}

TEST(RandomTest, ErdosRenyi) {
  stag::Graph testGraph = stag::erdos_renyi(1000, 0.1);
  EXPECT_EQ(testGraph.number_of_vertices(), 1000);
}

TEST(RandomTest, SBMArguments) {
  stag_int n = -1;
  stag_int k = 2;
  double p = 0.5;
  double q = 0.5;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);

  n = 0;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);

  n = 100;
  k = -1;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);

  k = 0;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);

  k = 51;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);

  k = 2;
  p = -0.1;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);

  p = 1.1;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);

  p = 0.5;
  q = -0.1;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);

  q = 1.1;
  EXPECT_THROW(stag::sbm(n, k, p, q), std::invalid_argument);
}

TEST(RandomTest, ErdosRenyiArguments) {
  stag_int n = -1;
  double p = 0.5;
  EXPECT_THROW(stag::erdos_renyi(n, p), std::invalid_argument);

  n = 0;
  EXPECT_THROW(stag::erdos_renyi(n, p), std::invalid_argument);

  n = 100;
  p = -0.1;
  EXPECT_THROW(stag::erdos_renyi(n, p), std::invalid_argument);

  p = 1.1;
  EXPECT_THROW(stag::erdos_renyi(n, p), std::invalid_argument);
}
