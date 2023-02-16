/**
 * Tests for the methods in the spectrum.h header file. Includes methods for
 * computing eigenvalues and eigenvectors.
 */
#include <stdexcept>
#include <iostream>
#include <gtest/gtest.h>
#include <graph.h>
#include <utility.h>
#include <random.h>
#include <spectrum.h>
#include <graphio.h>

TEST(SpectrumTest, NormalisedLaplacianEigensystem) {
  // Create a small complete graph
  stag_int n = 10;
  stag::Graph testGraph = stag::complete_graph(n);

  // Extract the normalised laplacian matrix
  const SprsMat* lap = testGraph.normalised_laplacian();

  // Compute the first few eigenvalues and eigenvectors - there should be no exceptions
  // thrown.
  stag_int k = 4;
  stag::EigenSystem eigensystem = stag::compute_eigensystem(lap, k);

  // The eigenvalues should be equal to 0, and (n-1) copies of n/(n-1).
  Eigen::VectorXd eigenvalues = get<0>(eigensystem);
  std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

  EXPECT_NEAR(eigenvalues[0], 0, 0.000001);
  for (int i = 1; i < eigenvalues.size(); i++) {
      EXPECT_NEAR(eigenvalues[i], ((float) n) / (n-1), 0.000001);
  }
}

TEST(SpectrumTest, RandomGraphSpectrum) {
    // Create a graph from the SBM
    stag_int n = 100;
    stag::Graph testGraph = stag::sbm(n, 2, 0.5, 0.01);

    // Extract the normalised laplacian matrix
    const SprsMat* lap = testGraph.normalised_laplacian();

    // Compute the first few eigenvalues and eigenvectors - there should be no exceptions
    // thrown.
    stag_int k = 3;
    stag::EigenSystem eigensystem = stag::compute_eigensystem(lap, k);

    // The eigenvalues should be equal to 0, something 'small' and something 'large'.
    Eigen::VectorXd eigenvalues = get<0>(eigensystem);
    std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

    EXPECT_NEAR(eigenvalues[0], 0, 0.000001);
    EXPECT_LE(eigenvalues[1], 0.2);
    EXPECT_GE(eigenvalues[2], 0.5);
}

TEST(SpectrumTest, DisconnectedGraph) {
    // Create the data for the graph adjacency matrix of a disconnected graph
    //     0 2 0 0
    //     2 0 0 0
    //     0 0 0 1
    //     0 0 1 0
    std::vector<stag_int> rowStarts = {0, 1, 2, 3, 4};
    std::vector<stag_int> colIndices = {1, 0, 3, 2};
    std::vector<double> values = {2, 2, 1, 1};

    // Create the stag Graph object
    stag::Graph testGraph = stag::Graph(rowStarts, colIndices, values);

    // Extract the unnormalised laplacian matrix
    const SprsMat* lap = testGraph.laplacian();

    // Compute the first 3 eigenvalues and eigenvectors.
    stag_int k = 3;
    stag::EigenSystem eigensystem = stag::compute_eigensystem(lap, k);

    // The eigenvalues should be equal to 0, 0, and something else.
    Eigen::VectorXd eigenvalues = get<0>(eigensystem);
    std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

    EXPECT_NEAR(eigenvalues[0], 0, 0.000001);
    EXPECT_NEAR(eigenvalues[1], 0, 0.000001);
    EXPECT_GE(eigenvalues[2], 0.1);
}

TEST(SpectrumTest, HugeGraph) {
  // Load the graph from file
  std::string filename = "test/data/test6.edgelist";
  stag::Graph graph = stag::load_edgelist(filename);

  // Extract the unnormalised laplacian matrix
  const SprsMat* lap = graph.laplacian();

  // Compute the first 3 eigenvalues and eigenvectors.
  stag_int k = 3;
  stag::EigenSystem eigensystem = stag::compute_eigensystem(lap, k);
}

TEST(SpectrumTest, ArgumentChecking) {
  stag_int n = 10;
  stag::Graph testGraph = stag::complete_graph(n);
  const SprsMat* lap = testGraph.laplacian();

  stag_int k = -1;
  EXPECT_THROW(stag::compute_eigensystem(lap, k), std::invalid_argument);

  k = 0;
  EXPECT_THROW(stag::compute_eigensystem(lap, k), std::invalid_argument);

  k = n + 1;
  EXPECT_THROW(stag::compute_eigensystem(lap, k), std::invalid_argument);

  k = n;
  EXPECT_THROW(stag::compute_eigensystem(lap, k), std::invalid_argument);
}
