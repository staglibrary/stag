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

TEST(RandomTest, SBM) {
  stag::Graph testGraph = stag::sbm(10000000, 10, 0.000001, 0.0000001);
  std::cerr << testGraph.number_of_edges() << std::endl;
}

TEST(RandomTest, ErdosRenyi) {
  stag::Graph testGraph = stag::erdos_renyi(20, 0.5);
  std::cerr << *testGraph.adjacency() << std::endl;
}
