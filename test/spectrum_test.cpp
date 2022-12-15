/**
 * Tests for the methods in the spectrum.h header file. Includes methods for
 * computing eigenvalues and eigenvectors.
 */
#include <stdexcept>
#include <iostream>
#include <gtest/gtest.h>
#include <graph.h>
#include <utility.h>
#include <spectrum.h>

TEST(SpectrumTest, NormalisedLaplacianEigensystem) {
  // Create a small complete graph
  stag_int n = 10;
  stag::Graph testGraph = stag::complete_graph(n);

  // Extract the normalised laplacian matrix
  const SprsMat* lap = testGraph.normalised_laplacian();

  // Compute the first few eigenvalues and eigenvectors
  stag_int k = 4;
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd> eigensystem = stag::compute_eigensystem(
      lap, k);

  // Display the computed eigenvalues and eigenvectors
  std::cout << get<0>(eigensystem) << std::endl;
  std::cout << get<1>(eigensystem) << std::endl;
}
