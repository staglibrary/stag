/**
 * Tests for the clustering methods in cluster.h.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <iostream>
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
