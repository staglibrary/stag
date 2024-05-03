/**
 * Tests for the methods in the random.h header file. Includes methods for
 * generating random graphs.
 *
 * This file is provided as part of the STAG library and released under the GPL
 * license.
 */
#include <iostream>
#include <gtest/gtest.h>
#include "graph.h"
#include "random.h"
#include <utility.h>
#include "graphio.h"

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

  std::vector<StagInt> refStarts = stag::sprsMatOuterStarts(
      referenceGraph.laplacian());
  std::vector<StagInt> testStarts = stag::sprsMatOuterStarts(
      testGraph.laplacian());
  EXPECT_EQ(refStarts, testStarts);

  std::vector<StagInt> refIndices= stag::sprsMatInnerIndices(
      referenceGraph.laplacian());
  std::vector<StagInt> testIndices = stag::sprsMatInnerIndices(
      testGraph.laplacian());
  EXPECT_EQ(refIndices, testIndices);

  std::vector<double> refValues = stag::sprsMatValues(
      referenceGraph.laplacian());
  std::vector<double> testValues = stag::sprsMatValues(
      testGraph.laplacian());
  EXPECT_FLOATS_NEARLY_EQ(refValues, testValues, 0.000001);
}

TEST(RandomTest, GeneralSBM) {
  std::vector<StagInt> cluster_sizes = {1000, 100, 10};
  DenseMat prob_mat {{0.4, 0.1, 0.1}, {0.1, 0.7, 0}, {0.1, 0, 1}};
  stag::Graph testGraph = stag::general_sbm(cluster_sizes, prob_mat);
  EXPECT_EQ(testGraph.number_of_vertices(), 1110);
}

TEST(RandomTest, GeneralSBMEdgelist) {
  std::string filename = "output.edgelist";
  std::vector<StagInt> cluster_sizes = {100, 100, 10};
  DenseMat prob_mat {{0.4, 0.1, 0.1}, {0.1, 0.7, 0}, {0.1, 0, 1}};
  stag::general_sbm_edgelist(filename, cluster_sizes, prob_mat);
  stag::Graph testGraph = stag::load_edgelist(filename);
  EXPECT_EQ(testGraph.number_of_vertices(), 210);
}

TEST(RandomTest, GeneralSBMArguments) {
  std::vector<StagInt> cluster_sizes = {1000, 100, 10};
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
  StagInt n = -1;
  StagInt k = 2;
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
  StagInt n = -1;
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

TEST(RandomTest, SBMLabels) {
  StagInt n = 6;
  StagInt k = 3;
  stag::Graph myGraph = stag::sbm(n, k, 0.8, 0.1);

  std::vector<StagInt> gt_labels = stag::sbm_gt_labels(n, k);
  std::vector<StagInt> expected_labels {0, 0, 1, 1, 2, 2};

  EXPECT_EQ(gt_labels, expected_labels);
}

TEST(RandomTest, GeneralSBMLabels) {
  std::vector<StagInt> cluster_sizes = {4, 2};
  DenseMat prob_mat {{0.4, 0.1}, {0.1, 0.7}};
  stag::Graph myGraph = stag::general_sbm(cluster_sizes, prob_mat);

  std::vector<StagInt> gt_labels = stag::general_sbm_gt_labels(cluster_sizes);
  std::vector<StagInt> expected_labels {0, 0, 0, 0, 1, 1};

  EXPECT_EQ(gt_labels, expected_labels);
}

TEST(RandomTest, SBMLabelsArguments) {
  StagInt n = -1;
  StagInt k = 2;
  EXPECT_THROW(stag::sbm_gt_labels(n, k), std::invalid_argument);

  n = 0;
  EXPECT_THROW(stag::sbm_gt_labels(n, k), std::invalid_argument);

  n = 100;
  k = -1;
  EXPECT_THROW(stag::sbm_gt_labels(n, k), std::invalid_argument);

  k = 0;
  EXPECT_THROW(stag::sbm_gt_labels(n, k), std::invalid_argument);

  k = 51;
  EXPECT_THROW(stag::sbm_gt_labels(n, k), std::invalid_argument);
}

TEST(RandomTest, GeneralSBMLabelsArguments) {
  std::vector<StagInt> cluster_sizes = {1000, 100, 0};
  EXPECT_THROW(stag::general_sbm_gt_labels(cluster_sizes), std::invalid_argument);

  cluster_sizes.at(2) = -10;
  EXPECT_THROW(stag::general_sbm_gt_labels(cluster_sizes), std::invalid_argument);
}

