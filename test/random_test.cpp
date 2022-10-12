/**
 * Tests for the methods in the random.h header file. Includes methods for
 * generating random graphs.
 *
 * Copyright 2022 Peter Macgregor
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
