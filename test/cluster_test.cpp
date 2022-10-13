/**
 * Tests for the clustering methods in cluster.h.
 *
 * Copyright 2022 Peter Macgregor
 */
#include <iostream>
#include <gtest/gtest.h>
#include <cluster.h>
#include <graph.h>
#include <random.h>

TEST(ClusterTest, ACL) {
  // Construct a test graph
  stag::Graph testGraph = stag::sbm(20, 2, 0.5, 0);

  // Find a cluster in the test graph
  std::vector<int> cluster = stag::local_cluster(&testGraph, 0);

  // Display the cluster
  for (auto i : cluster) {
    std::cerr << i << ", ";
  }
  std::cerr << std::endl;
}
