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
