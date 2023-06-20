//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#include "iostream"
#include <cmath>
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

/**
 * Compute the first \f$n\f$ Arnoldi vectors of \f$A\f$, beginning with \f$b\f$.
 *
 * The Arnoldi vectors are an orthogonalization of the vectors
 * \f[
 *  b, A b, A^2 b, \ldots, A^n b,
 * \f]
 * which are the basis vectors of the Krylov subspace.
 *
 * In conjugate Arnoldi, the Krylov basis vectors are orthogonalised with
 * respect to the \f$A\f$ inner product.
 *
 * @param A
 * @param b
 * @param dim
 * @return
 */
DenseMat conjugate_arnoldi(const SprsMat* A,
                           DenseVec &b,
                           stag_int n) {
  // Comments will give 'latex' type pseudocode for the algorithm as presented
  // in many sources on the Arnoldi method.
  DenseMat Q(b.rows(), n);

  // q_i = b
  Q.col(0) = b / std::sqrt(b.dot(*A * b));

  // FOR i := 2, ... n:
  for (stag_int i = 1; i < n; i++) {
    // x_i = A q_{i-1}
    Q.col(i) = *A * Q.col(i-1);

    // FOR j := 1, ..., i-1
    for (stag_int j = 0; j < i; j++) {
      // h_{i,j} = x_i^T A q_j
      double hij = Q.col(i).dot(*A * Q.col(j));

      // x_i = x_i - h_{i, j} q_j
      Q.col(i) = Q.col(i) - hij * Q.col(j);
      assert(std::abs(Q.col(i).dot(*A * Q.col(j))) < 0.000001);
    } // END FOR

    // h_{i, i} = \sqrt{x_i^T A x_i}
    double hii = std::sqrt(Q.col(i).dot(*A * Q.col(i)));

    // q_i = x_i / h_{i, i}
    Q.col(i) = (1 / hii) * Q.col(i);
    assert(std::abs(Q.col(i).dot(*A * Q.col(i)) - 1) < 0.000001);
  } // END FOR

  // RETURN q_1, ..., q_n
  return Q;
}

DenseMat conjugate_arnoldi(const SprsMat* A,
                           DenseVec &b) {
  return conjugate_arnoldi(A, b, A->cols());
}

DenseVec stag::solve_laplacian_exact_conjugate_gradient(Graph *g, DenseVec &b) {
  // Exactly solving a Laplacian system by the conjugate gradient method
  // involves the steps:
  //   1. Compute a set of n conjugate vectors p_1, ..., p_n
  //   2. Compute the coefficients a_k = (p_k^T b) / (p_k^T A p_k)
  //   3. Return x = sum_k a_k p_k
  DenseMat P = conjugate_arnoldi(g->laplacian(), b);

  DenseVec x(b.rows());
  x.setZero();

  for (stag_int k = 0; k < b.rows(); k++) {
    // We do not bother diving by p_k^T A p_k since this is guaranteed to be
    // 1 by the output of conjugate_arnoldi.
    double ak = P.col(k).dot(b);
    x += ak * P.col(k);

    // Check for convergence - if the Laplacian has an eigenvalue close to 0,
    // we find that sometimes the process converges one step early.
    if ((*g->laplacian() * x - b).norm() <= EPSILON) break;
  }

  return x;
}
