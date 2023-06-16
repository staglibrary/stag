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

DenseVec stag::solve_laplacian(Graph* g, DenseVec& b, double eps) {
  return solve_laplacian_gauss_seidel(g, b, eps);
}

//------------------------------------------------------------------------------
// Jacobi iteration
//------------------------------------------------------------------------------

DenseVec stag::solve_laplacian_jacobi(stag::Graph* g, DenseVec& b, double eps,
                                      stag_int max_iterations) {
  DenseVec x = stag::jacobi_iteration(g->laplacian(), b, eps,
                                      max_iterations);
  return x;
}

DenseVec stag::solve_laplacian_jacobi(stag::Graph* g, DenseVec& b, double eps) {
  return stag::solve_laplacian_jacobi(g, b, eps, STAG_MAX_ITERATIONS);
}

/**
 * Test whether the provided sparse matrix is symmetric diagonally dominant.
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

DenseVec stag::jacobi_iteration(const SprsMat* A, DenseVec& b, double eps,
                                stag_int max_iterations) {
  // The preconditioner is the diagonal matrix of A.
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> P = A->diagonal().asDiagonal();
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> P_inv = P.inverse();
  SprsMat I(A->cols(), A->cols());
  I.setIdentity();

  // We will iterate until the error is below the provided epsilon
  DenseVec x_t = b;
  double error = solution_error(A, b, x_t);
  stag_int iteration = 0;
  while (error > eps) {
    x_t = x_t + P_inv * (b - *A * x_t);

    error = solution_error(A, b, x_t);

    // Check whether we have passed the maximum number of iterations.
    iteration++;
    if (iteration > max_iterations) throw stag::ConvergenceError();
  }

  return x_t;
}

//------------------------------------------------------------------------------
// Gauss-Seidel Method
//------------------------------------------------------------------------------

DenseVec stag::solve_laplacian_gauss_seidel(stag::Graph* g, DenseVec& b,
                                            double eps,
                                            stag_int max_iterations) {
  DenseVec x = stag::gauss_seidel_iteration(g->laplacian(), b, eps,
                                            max_iterations);
  return x;
}

DenseVec stag::solve_laplacian_gauss_seidel(stag::Graph* g, DenseVec& b,
                                            double eps) {
  return stag::solve_laplacian_gauss_seidel(g, b, eps, STAG_MAX_ITERATIONS);
}

DenseVec stag::gauss_seidel_iteration(const SprsMat *A,
                                      DenseVec &b,
                                      double eps,
                                      stag_int max_iterations) {
  // We will make use of the upper triangular part of A
  Eigen::TriangularView<const SprsMat, Eigen::StrictlyUpper> U =
      A->triangularView<Eigen::StrictlyUpper>();

  // In each iteration, we compute the solution to
  //    L x_{k+1} = b - U x_k,
  // where L is the lower triangular part of A.
  // We will iterate until the error is below the provided epsilon.
  DenseVec x_t = b;
  double error = solution_error(A, b, x_t);
  stag_int iteration = 0;
  while (error > eps) {
    // Solve the lower-triangular system
    //   L x = b - U x
    DenseVec bp = b - U * x_t;
    for (auto i = 0; i < A->rows(); i++) {
      double target = bp.coeff(i);
      for (auto j = 0; j < i; j++) {
        target -= A->coeff(i, j) * x_t.coeff(j);
      }
      x_t.coeffRef(i) = target / A->coeff(i,i);
    }

    error = solution_error(A, b, x_t);

    // Check whether we have passed the maximum number of iterations.
    iteration++;
    if (iteration > max_iterations) throw stag::ConvergenceError();
  }

  return x_t;
}
