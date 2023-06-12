/**
 * Tests for the methods in the solve.h header file. Includes methods for
 * solving Laplacian systems.
 */
#include <stdexcept>
#include <iostream>
#include <gtest/gtest.h>
#include <graph.h>
#include <solve.h>

TEST(SolveTest, JacobiIteration) {
  // Create a small cycle graph
  stag_int n = 3;
  stag::Graph testGraph = stag::cycle_graph(n);

  // Construct the target vector
  DenseVec b(n);
  b << 2, 4, 1; //, 6, 4, 5, 3, 2, 1, 3;

  // Solve a linear system
  DenseVec soln = stag::solve_laplacian(&testGraph, b, 0.1);
}

