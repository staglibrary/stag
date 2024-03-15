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
#include <graphio.h>

// Define some helper test assertions.
#define EXPECT_FLOATS_NEARLY_EQ(expected, actual, thresh) \
        EXPECT_EQ(expected.size(), actual.size()) << "Array sizes differ.";\
        for (size_t idx = 0; idx < std::min(expected.size(), actual.size()); ++idx) \
        { \
            EXPECT_NEAR(expected[idx], actual[idx], thresh) << "at index: " << idx;\
        }

TEST(ClusterTest, SpectralCluster) {
  // Construct a test graph from the SBM
  StagInt n = 1000;
  StagInt k = 5;
  stag::Graph testGraph = stag::sbm(n, k, 0.6, 0.1);

  // Find the clusters
  auto clusters = stag::spectral_cluster(&testGraph, k);

  // There should be approximately the same number of each cluster
  StagInt c1 = 0;
  StagInt c2 = 0;

  for (auto c : clusters) {
    if (c == 1) c1++;
    if (c == 2) c2++;
  }
  EXPECT_NEAR(c1 / c2, 1, 0.01);
}

TEST(ClusterTest, SCArguments) {
  // Construct a test graph
  StagInt n = 100;
  stag::Graph testGraph = stag::erdos_renyi(n, 0.1);

  // Passing an invalid number of clusters to the spectral clustering
  // algorithm should return an error.
  StagInt k = -1;
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
  StagInt n = 10000;
  StagInt k = 5;
  stag::Graph testGraph = stag::sbm(n, k, 0.06, 0.01);

  // Find the clusters
  auto clusters = stag::spectral_cluster(&testGraph, k);

  // There should be approximately the same number of each cluster
  StagInt c1 = 0;
  StagInt c2 = 0;

  for (auto c : clusters) {
    if (c == 1) c1++;
    if (c == 2) c2++;
  }
  EXPECT_NEAR(c1, c2, 1.5 * c1);
  EXPECT_NEAR(c1, n / k, 1.5 * c1);
  EXPECT_NEAR(c2, n / k, 1.5 * c2);
}

TEST(ClusterTest, SpectralClusterDisconnected) {
  // Construct a sparse disconnected graph from the SBM
  StagInt n = 1000;
  StagInt k = 2;
  stag::Graph testGraph = stag::sbm(n, k, 0.5, 0);

  // Find the clusters
  auto clusters = stag::spectral_cluster(&testGraph, k);

  // There should be exactly the same number of each cluster
  StagInt c0 = 0;
  StagInt c1 = 0;

  for (auto c : clusters) {
    if (c == 0) c0++;
    if (c == 1) c1++;
  }
  EXPECT_EQ(c0, n / k);
  EXPECT_EQ(c1, n / k);
}

TEST(ClusterTest, ApproxPageRank) {
  // Construct a test graph
  std::vector<StagInt> rowStarts = {0, 2, 4, 6, 8};
  std::vector<StagInt> colIndices = {1, 2, 0, 3, 0, 3, 1, 2};
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

TEST(ClusterTest, ApproxPageRankNoPush) {
  // Test the approximate pagerank algorithm when there is no push operation.
  // This happens when the degree of the starting vertex is greater than 1/eps.

  // Construct a test graph
  std::vector<StagInt> rowStarts = {0, 3, 5, 7, 10};
  std::vector<StagInt> colIndices = {1, 2, 3, 0, 3, 0, 3, 0, 1, 2};
  std::vector<double> values = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  stag::Graph testGraph(rowStarts, colIndices, values);

  // Run the approximate pagerank method
  SprsMat seed(1, 1);
  seed.coeffRef(0, 0) = 1;
  std::tuple<SprsMat, SprsMat> apr = stag::approximate_pagerank(&testGraph, seed, 1./3, 1./2);

  // Define the expected pagerank vector - this is the empty sparse matrix.
  std::vector<StagInt> colStarts = {0, 0};
  std::vector<StagInt> rowIndices;
  std::vector<StagInt> valuesVec;

  // Check that the pagerank vector has the form that we expect
  std::vector<StagInt> newStarts = stag::sprsMatOuterStarts(&std::get<0>(apr));
  std::vector<StagInt> newIndices = stag::sprsMatInnerIndices(&std::get<0>(apr));
  std::vector<double> newValues = stag::sprsMatValues(&std::get<0>(apr));
  EXPECT_EQ(colStarts, newStarts);
  EXPECT_EQ(rowIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(valuesVec, newValues, 0.000001);
}

TEST(ClusterTest, ACL) {
  // Construct a test barbell graph
  stag::Graph testGraph = stag::barbell_graph(10);

  // Run the acl clustering method
  std::vector<StagInt> cluster = stag::local_cluster_acl(&testGraph, 0, 0.8,
                                                         0.0001);
  std::stable_sort(cluster.begin(), cluster.end());

  // Check that we found one of the clusters.
  std::vector<StagInt> expected_cluster = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
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
  std::vector<StagInt> cluster = stag::local_cluster(&testGraph, 0, testGraph.total_volume() / 4);
  std::stable_sort(cluster.begin(), cluster.end());

  // Check the number of returned vertices
  EXPECT_GE(cluster.size(), 100);
  EXPECT_LE(cluster.size(), 1000);

  // Let's say that 50% of the found cluster should lie inside the first cluster
  // in the SBM graph
  StagInt inside = 0;
  for (auto v : cluster) {
    if (v < 500) inside++;
  }
  EXPECT_GE(inside / cluster.size(), 0.5);
}

TEST(ClusterTest, LocalClusterArguments) {
  // Create a test graph
  StagInt n = 100;
  stag::Graph testGraph = stag::cycle_graph(n);

  // Check vertex is valid
  StagInt v = -1;
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

  // A very small (< 1) target_volume should work
  target_vol = 0.1;
  EXPECT_NO_THROW(stag::local_cluster(&testGraph, v, target_vol));
}

TEST(ClusterTest, ACLArguments) {
  // Create a test graph
  StagInt n = 100;
  stag::Graph testGraph = stag::cycle_graph(n);

  // Check vertex is valid
  StagInt v = -1;
  double alpha = 0.1;
  double eps = 0.1;
  EXPECT_THROW(stag::local_cluster_acl(&testGraph, v, alpha, eps), std::invalid_argument);

  v = n;
  EXPECT_THROW(stag::local_cluster_acl(&testGraph, v, alpha, eps), std::invalid_argument);

  // Check alpha is valid
  v = 0;
  alpha = -0.1;
  EXPECT_THROW(stag::local_cluster_acl(&testGraph, v, alpha, eps), std::invalid_argument);

  alpha = 2;
  EXPECT_THROW(stag::local_cluster_acl(&testGraph, v, alpha, eps), std::invalid_argument);

  // Check epsilon is valid
  alpha = 0.1;
  eps = -0.1;
  EXPECT_THROW(stag::local_cluster_acl(&testGraph, v, alpha, eps), std::invalid_argument);

  eps = 0;
  EXPECT_THROW(stag::local_cluster_acl(&testGraph, v, alpha, eps), std::invalid_argument);
}

TEST(ClusterTest, PageRankArguments) {
  // Create a test graph
  StagInt n = 100;
  stag::Graph testGraph = stag::cycle_graph(n);

  // Check seed vector is valid
  SprsMat v(101, 1);
  v.coeffRef(0, 0) = 0.5;
  v.coeffRef(100, 0) = 0.5;
  double alpha = 0.1;
  double eps = 0.1;
  EXPECT_THROW(stag::approximate_pagerank(&testGraph, v, alpha, eps), std::invalid_argument);

  // Check alpha is valid
  SprsMat v2(100, 1);
  v.coeffRef(0, 0) = 1;
  alpha = -0.1;
  EXPECT_THROW(stag::approximate_pagerank(&testGraph, v2, alpha, eps), std::invalid_argument);

  alpha = 2;
  EXPECT_THROW(stag::approximate_pagerank(&testGraph, v2, alpha, eps), std::invalid_argument);

  // Check epsilon is valid
  alpha = 0.1;
  eps = -0.1;
  EXPECT_THROW(stag::approximate_pagerank(&testGraph, v2, alpha, eps), std::invalid_argument);

  eps = 0;
  EXPECT_THROW(stag::approximate_pagerank(&testGraph, v2, alpha, eps), std::invalid_argument);
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
  std::vector<StagInt> sweep_set = stag::sweep_set_conductance(&testGraph, vec);
  std::stable_sort(sweep_set.begin(), sweep_set.end());

  // Test the result
  std::vector<StagInt> expected_sweep = {0, 1, 2, 3};
  EXPECT_EQ(sweep_set, expected_sweep);
}

TEST(ClusterTest, ARI) {
  std::vector<StagInt> gt_labels {0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
  std::vector<StagInt> labels    {0, 1, 0, 1, 1, 2, 2, 2, 2, 2};

  double expected_ari = 0.31257344;
  double actual_ari = stag::adjusted_rand_index(gt_labels, labels);
  EXPECT_NEAR(actual_ari, expected_ari, 0.0001);

  std::vector<StagInt> labels2   {1, 1, 2, 2, 2, 2, 0, 0, 0, 0};
  expected_ari = 1;
  actual_ari = stag::adjusted_rand_index(gt_labels, labels2);
  EXPECT_NEAR(actual_ari, expected_ari, 0.0001);
}

TEST(ClusterTest, ARIArguments) {
  std::vector<StagInt> gt_labels {-1, -1, 1, 1, 1, 1, 2, 2, 2, 2};
  std::vector<StagInt> labels    {0, 1, 0, 1, 1, 2, 2, 2, 2, 2};
  EXPECT_THROW(stag::adjusted_rand_index(gt_labels, labels), std::invalid_argument);

  gt_labels.at(0) = 0;
  gt_labels.at(1) = 0;
  labels.at(3) = -1;
  EXPECT_THROW(stag::adjusted_rand_index(gt_labels, labels), std::invalid_argument);
}

TEST(ClusterTest, NMI) {
  std::vector<StagInt> gt_labels {0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
  std::vector<StagInt> labels    {0, 1, 0, 1, 1, 2, 2, 2, 2, 2};

  double expected_nmi = 0.4558585;
  double actual_nmi = stag::normalised_mutual_information(gt_labels, labels);
  EXPECT_NEAR(actual_nmi, expected_nmi, 0.0001);

  std::vector<StagInt> labels2   {1, 1, 2, 2, 2, 2, 0, 0, 0, 0};
  expected_nmi = 1;
  actual_nmi = stag::adjusted_rand_index(gt_labels, labels2);
  EXPECT_NEAR(actual_nmi, expected_nmi, 0.0001);
}

TEST(ClusterTest, NMIArguments) {
  std::vector<StagInt> gt_labels {-1, -1, 1, 1, 1, 1, 2, 2, 2, 2};
  std::vector<StagInt> labels    {0, 1, 0, 1, 1, 2, 2, 2, 2, 2};
  EXPECT_THROW(stag::normalised_mutual_information(gt_labels, labels), std::invalid_argument);

  gt_labels.at(0) = 0;
  gt_labels.at(1) = 0;
  labels.at(3) = -1;
  EXPECT_THROW(stag::normalised_mutual_information(gt_labels, labels), std::invalid_argument);
}

TEST(ClusterTest, ALLGLocalClustering) {
  std::string filename = "test/data/hugegraph.adjacencylist";
  stag::AdjacencyListLocalGraph testGraph(filename);

  // Find a cluster using the default local clustering algorithm
  std::vector<StagInt> cluster = stag::local_cluster(&testGraph,
                                                     500000,
                                                     10000);

  // Check the number of returned vertices
  EXPECT_GE(cluster.size(), 100);
  EXPECT_LE(cluster.size(), 10000);
}

TEST(ClusterTest, ALLGLocalClusteringConsistent) {
  std::string adj_filename = "output.adjacencylist";
  std::string edge_filename = "output.edgelist";

  // Create a graph from the SBM
  StagInt n = 500;
  StagInt k = 5;
  double p = 0.1;
  double q = 0.01;

  // Create the general sbm parameters
  std::vector<StagInt> cluster_sizes;
  DenseMat probabilities(k, k);
  for (auto i = 0; i < k; i++) {
    cluster_sizes.push_back(floor(((double) n) / ((double) k)));
    probabilities(i, i) = p;

    for (auto j = i + 1; j < k; j++) {
      probabilities(i, j) = q;
      probabilities(j, i) = q;
    }
  }
  stag::general_sbm_edgelist(edge_filename, cluster_sizes, probabilities);

  // Convert it to an adjacencylist
  stag::edgelist_to_adjacencylist(edge_filename, adj_filename);

  // Get the clusters two ways: Graph and AdjacencyListLocalGraph
  stag::Graph fullGraph = stag::load_adjacencylist(adj_filename);
  stag::AdjacencyListLocalGraph adjGraph(adj_filename);

  // Clustering on the graph loaded into memory should give the same result.
  std::vector<StagInt> adj_cluster = stag::local_cluster(&adjGraph,
                                                         0,
                                                         1500);
  std::vector<StagInt> full_cluster = stag::local_cluster(&fullGraph,
                                                          0,
                                                          1500);
  std::stable_sort(adj_cluster.begin(), adj_cluster.end());
  std::stable_sort(full_cluster.begin(), full_cluster.end());

  EXPECT_EQ(adj_cluster, full_cluster);
}

TEST(ClusterTest, Conductance) {
  // Construct a test graph
  std::vector<StagInt> rowStarts = {0, 2, 4, 6, 8};
  std::vector<StagInt> colIndices = {1, 2, 0, 3, 0, 3, 1, 2};
  std::vector<double> values = {1./2, 1./2, 1./2, 1./2, 1./2, 1./2, 1./2, 1./2};
  stag::Graph testGraph(rowStarts, colIndices, values);

  // Check the conductance
  std::vector<StagInt> cluster = {0, 3};
  double cond = stag::conductance(&testGraph, cluster);
  EXPECT_NEAR(cond, 1, 0.00001);

  cluster = {1, 3};
  cond = stag::conductance(&testGraph, cluster);
  EXPECT_NEAR(cond, 1./2, 0.00001);

  // Check the conductance method for an AdjacencyListLocalGraph
  std::string filename = "test/data/test3.adjacencylist";
  stag::AdjacencyListLocalGraph adjGraph(filename);
  cluster = {0, 1};
  cond = stag::conductance(&adjGraph, cluster);
  EXPECT_NEAR(cond, 1.5/3.5, 0.00001);

  // Conductance of empty set is 0
  cluster = {};
  cond = stag::conductance(&adjGraph, cluster);
  EXPECT_EQ(cond, 0);
}

TEST(ClusterTest, ConductanceArguments) {
  // Construct a test graph
  std::vector<StagInt> rowStarts = {0, 2, 4, 6, 8};
  std::vector<StagInt> colIndices = {1, 2, 0, 3, 0, 3, 1, 2};
  std::vector<double> values = {1./2, 1./2, 1./2, 1./2, 1./2, 1./2, 1./2, 1./2};
  stag::Graph testGraph(rowStarts, colIndices, values);

  // Negative integers in the cluster should throw an argument exception
  std::vector<StagInt> cluster = {0, -1};
  EXPECT_THROW(stag::conductance(&testGraph, cluster), std::invalid_argument);
}

TEST(ClusterTest, SymmetricDifference) {
  // Construct some sets
  std::vector<StagInt> S = {0, 3, 2, 6};
  std::vector<StagInt> T = {3, 1, 5, 2};

  // Compute the symmetric difference
  std::vector<StagInt> calculated_difference = stag::symmetric_difference(S, T);
  std::vector<StagInt> expected_difference = {0, 1, 5, 6};
  EXPECT_EQ(calculated_difference, expected_difference);

  // If a value appears twice, the duplicates should be ignored
  S.push_back(3);
  calculated_difference = stag::symmetric_difference(S, T);
  EXPECT_EQ(calculated_difference, expected_difference);
  T.push_back(5);
  calculated_difference = stag::symmetric_difference(S, T);
  EXPECT_EQ(calculated_difference, expected_difference);

  // If one of the vectors is empty, then the symmetric difference is equal
  // to the non-empty vector.
  S.clear();
  calculated_difference = stag::symmetric_difference(S, T);
  EXPECT_EQ(calculated_difference, T);
}

TEST(ClusterTest, ConnComp) {
  // Create a disconnected graph
  StagInt n = 10;
  StagInt k = 2;
  stag::Graph testGraph = stag::sbm(n, k, 1, 0);

  // Check the connected components
  std::vector<StagInt> comp = stag::connected_component(&testGraph, 0);
  std::sort(comp.begin(), comp.end());
  std::vector<StagInt> expected_comp = {0, 1, 2, 3, 4};
  EXPECT_EQ(comp, expected_comp);

  comp = stag::connected_component(&testGraph, 7);
  std::sort(comp.begin(), comp.end());
  expected_comp = {5, 6, 7, 8, 9};
  EXPECT_EQ(comp, expected_comp);

  // Find the connected components by the full decomposition algorithm
  std::vector<std::vector<StagInt>> components = stag::connected_components(
      &testGraph);
  comp = components.at(0);
  std::sort(comp.begin(), comp.end());
  expected_comp = {0, 1, 2, 3, 4};
  EXPECT_EQ(comp, expected_comp);
  comp = components.at(1);
  std::sort(comp.begin(), comp.end());
  expected_comp = {5, 6, 7, 8, 9};
  EXPECT_EQ(comp, expected_comp);
}

TEST(ClusterTest, ConnCompALLG){
  std::string filename = "test/data/test8.adjacencylist";
  stag::AdjacencyListLocalGraph testGraph(filename);

  // Check the connected components
  std::vector<StagInt> comp = stag::connected_component(&testGraph, 2);
  std::sort(comp.begin(), comp.end());
  std::vector<StagInt> expected_comp = {1, 2, 3};
  EXPECT_EQ(comp, expected_comp);

  comp = stag::connected_component(&testGraph, 6);
  std::sort(comp.begin(), comp.end());
  expected_comp = {4, 5, 6};
  EXPECT_EQ(comp, expected_comp);
}

TEST(ClusterTest, CheegerCut){
  // Construct a test graph from the SBM
  StagInt n = 1000;
  StagInt k = 2;
  stag::Graph testGraph = stag::sbm(n, k, 0.6, 0.1);

  // Find the clusters
  auto clusters = stag::cheeger_cut(&testGraph);

  // There should be approximately the same number of each cluster
  StagInt c0 = 0;
  StagInt c1 = 0;

  for (auto c : clusters) {
    if (c == 0) c0++;
    if (c == 1) c1++;
  }
  EXPECT_NE(c0, 0);
  EXPECT_NE(c1, 0);
  EXPECT_NEAR(c0 / c1, 1, 0.01);
}

TEST(ClusterTest, CheegerCutComplete){
  // Construct a complete graph
  StagInt n = 200;
  stag::Graph testGraph = stag::complete_graph(n);

  // Find the cheeger cut of the graph
  std::vector<StagInt> clusters = stag::cheeger_cut(&testGraph);

  // There should be exactly n / 2 vertices in one cluster.
  StagInt c0 = 0;
  for (auto c : clusters) {
    if (c == 0) c0++;
  }
  EXPECT_EQ(c0, n / 2);
}
