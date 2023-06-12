//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#include "iostream"
#include "graph.h"
#include "utility.h"
#include "solve.h"

#include <utility>

double solution_error(const SprsMat* A, DenseVec& b, DenseVec& x) {
  DenseVec e_vec = *A * x - b;
  return e_vec.norm();
}

DenseVec solve_laplacian_jacobi(stag::Graph* g, DenseVec& b, double eps) {
  // Solve the system by Jacobi iteration. Since the graph Laplacian is not
  // strictly diagonally dominant, we add a small positive constant to the
  // diagonal.
  SprsMat I(g->number_of_vertices(), g->number_of_vertices());
  I.setIdentity();
  SprsMat A = *g->laplacian() + (eps) * I;
  DenseVec x = stag::jacobi_iteration(&A, b, eps);
  std::cout << "Final Error: " << solution_error(g->laplacian(), b, x) << std::endl;
  return x;
}

DenseVec stag::solve_laplacian(Graph* g, DenseVec& b, double eps) {
  return solve_laplacian_jacobi(g, b, eps);
}

/**
 * Test whether the provided sparse matrix is symmatric diagonally dominant.
 *
 * @param A the matrix to test
 * @return a boolean result
 */
bool symmetric_diagonally_dominant(const SprsMat* A) {
  if (!stag::isSymmetric(A)) {
    return false;
  }

  // Get the diagonal elements from A
  SprsMat D = ((DenseMat) A->diagonal().asDiagonal()).sparseView();

  // Construct a version of A without the diagonal and with absolute values
  SprsMat A_abs = *A - D;
  for (stag_int k = 0; k < A_abs.outerSize(); k++) {
    for (SprsMat::InnerIterator it(A_abs, k); it; ++it) {
      it.valueRef() = std::abs(it.value());
    }
  }

  // Construct a vector of the sum of absolute values of the off-diagonal entries
  Eigen::VectorXd degrees = A_abs * Eigen::VectorXd::Ones(A_abs.cols());

  // Check for dominance in each column
  for (auto i = 0; i < A_abs.cols(); i++) {
    if (D.coeff(i, i) <= degrees.coeff(i)) return false;
  }

  // If we get here, the matrix is diagonally dominant
  return true;
}

DenseVec stag::jacobi_iteration(const SprsMat* A, DenseVec& b, double eps) {
  // Check that A is strictly diagonally dominant
  assert(symmetric_diagonally_dominant(A));

  // The preconditioner is the diagonal matrix of A.
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> P = A->diagonal().asDiagonal();
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> P_inv = P.inverse();
  SprsMat I(A->cols(), A->cols());
  I.setIdentity();

  // We will iterate until the error is below the provided epsilon
  DenseVec x_t = b;
  double error = solution_error(A, b, x_t);
  while (error > eps) {
    x_t = x_t + P_inv * (b - *A * x_t);

    double new_error = solution_error(A, b, x_t);
    assert(new_error < error);
    error = new_error;

    std::cout << "Error: " << error << std::endl;
    std::cout << x_t << std::endl;
  }

  return x_t;
}
