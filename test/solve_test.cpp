/**
 * Tests for the methods in the solve.h header file. Includes methods for
 * solving Laplacian systems.
 */
#include <stdexcept>
#include <gtest/gtest.h>
#include <graph.h>
#include <random.h>
#include <solve.h>

TEST(SolveTest, JacobiIterationConvergenceFailure) {
  // Create a small cycle graph - Jacobi iteration cannot solve this system.
  stag_int n = 5;
  stag::Graph testGraph = stag::cycle_graph(n);

  // Construct the target vector
  DenseVec b(n);
  b << 2, 4, 1, 6, 4;

  // Solve a linear system
  EXPECT_THROW(stag::solve_laplacian_jacobi(&testGraph, b, 0.1),
               stag::ConvergenceError);
}

TEST(SolveTest, JacobiIteration) {
  // Create a small second difference graph
  stag_int n = 5;
  stag::Graph testGraph = stag::second_difference_graph(n);

  // Construct the target vector
  DenseVec b(n);
  b << 2, 4, 1, 6, 4;

  // Solve the linear system
  DenseVec x_hat = stag::solve_laplacian_jacobi(&testGraph, b, 0.1);

  // Compare with the expected answer
  DenseVec x(n);
  x << (15./2), 13, (29./2), 15, (19./2);
  EXPECT_TRUE(x_hat.isApprox(x, 0.01));
}

TEST(SolveTest, GaussSeidelConvergenceFailure) {
  // Create a small cycle graph - Gauss-Seidel iteration cannot solve this system.
  stag_int n = 5;
  stag::Graph testGraph = stag::cycle_graph(n);

  // Construct the target vector
  DenseVec b(n);
  b << 2, 4, 1, 6, 4;

  // Solve a linear system
  EXPECT_THROW(stag::solve_laplacian_gauss_seidel(&testGraph, b, 0.1),
               stag::ConvergenceError);
}

TEST(SolveTest, GaussSeidel) {
  // Create a small second difference graph
  stag_int n = 5;
  stag::Graph testGraph = stag::second_difference_graph(n);

  // Construct the target vector
  DenseVec b(n);
  b << 2, 4, 1, 6, 4;

  // Solve the linear system
  DenseVec x_hat = stag::solve_laplacian_gauss_seidel(&testGraph, b, 0.01);

  // Compare with the expected answer
  DenseVec x(n);
  x << (15./2), 13, (29./2), 15, (19./2);
  EXPECT_TRUE(x_hat.isApprox(x, 0.01));
}

TEST(SolveTest, ConjugateGradientExact) {
  // Create a small second difference graph
  stag_int n = 5;
  stag::Graph testGraph = stag::second_difference_graph(n);

  // Construct the target vector
  DenseVec b(n);
  b << 2, 4, 1, 6, 4;

  // Solve the linear system
  DenseVec x_hat = stag::solve_laplacian_exact_conjugate_gradient(&testGraph, b);

  // Check that we're close to the solution
  DenseVec b_hat = *testGraph.laplacian() * x_hat;

  EXPECT_TRUE(b_hat.isApprox(b, 0.0001));

  // Try the same with a random graph
  n = 10;
  testGraph = stag::erdos_renyi(n, 0.5);
  DenseVec b2(n);
  b2 << 4, 2, -4, -1, 3, -2, 3, -6, 1, 0; // Must be orthogonal to all-ones
  x_hat = stag::solve_laplacian_exact_conjugate_gradient(&testGraph, b2);
  b_hat = *testGraph.laplacian() * x_hat;
  EXPECT_TRUE(b_hat.isApprox(b2, 0.0001));
}
