//
// Methods for solving Laplacian systems of equations.
//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file solve.h
 * \brief Methods for solving Laplacian systems of equations.
 */

#ifndef STAG_TEST_SOLVE_H
#define STAG_TEST_SOLVE_H

#include <stdexcept>
#include "graph.h"

/**
 * \cond do not document the STAG_MAX_ITERATIONS definition
 */
// Define a default maximum number of iterations for iterative methods.
#define STAG_MAX_ITERATIONS 1000

/**
 * \endcond
 */

namespace stag {

  /**
   * Solve the Laplacian system \f$L x = b\f$, where \f$L\f$ is the Laplacian
   * of the provided STAG graph object.
   *
   * The method for solving the system is chosen automatically by STAG.
   *
   * @param g
   * @param b
   * @param eps
   * @return
   */
  DenseVec solve_laplacian(Graph* g, DenseVec& b, double eps);

  /**
   * Solve the Laplacian system \f$L x = b\f$ by Jacobi iteration.
   *
   * At each iteration, we solve the system
   *
   * \f[
   *    P x_{k+1} = (P - L) x_k + b,
   * \f]
   *
   * where \f$P = \mathrm{diag}(L)\f$ is the diagonal of the Laplacian matrix
   * \f$L\f$.
   * The error at iteration \f$k\f$ is given by
   *
   * \f[
   *    e_k = \| L x_k - b \|_2,
   * \f]
   *
   * and we terminate the algorithm when \f$e_k\f$ is less than the provided
   * error parameter.
   *
   * \note
   * Jacobi iteration is guaranteed to converge if the Laplacian matrix is
   * Strictly Diagonally Dominant (SDD), and may also converge in other cases.
   *
   * @param g the graph representing the Laplacian matrix to be used
   * @param b the vector \f$b\f$
   * @param eps the error parameter \f$\epsilon\f$ controlling the permitted
   *            approximation error.
   * @param max_iterations (optional) the maximum number of iterations to perform.
   *                       If this parameter is omitted, STAG will automatically
   *                       set the maximum iterations.
   * @return the approximate solution \f$\hat{x}\f$ such that
   *         \f$\| L \hat{x} - b \|_2 \leq \epsilon\f$.
   * @throws stag::ConvergenceError if the algorithm does not converge
   */
  DenseVec solve_laplacian_jacobi(Graph* g, DenseVec& b, double eps,
                                  stag_int max_iterations);

  /**
   * @overload
   */
  DenseVec solve_laplacian_jacobi(Graph* g, DenseVec& b, double eps);

  /**
   * Solve a linear system \f$A x = b\f$ by Jacobi iteration.
   *
   * For more information about Jacobi iteration, see
   * stag::solve_laplacian_jacobi.
   *
   * @param A the matrix \f$A\f$
   * @param b the vector \f$b\f$
   * @param eps the error parameter \f$\epsilon\f$ controlling the permitted
   *            approximation error
   * @param max_iterations the maximum number of iterations to perform
   * @return the approximate solution \f$\hat{x}\f$ such that
   *         \f$\| A \hat{x} - b \|_2 \leq \epsilon\f$.
   * @throws stag::ConvergenceError if the algorithm does not converge
   */
  DenseVec jacobi_iteration(const SprsMat* A, DenseVec& b, double eps,
                            stag_int max_iterations);

  //----------------------------------------------------------------------------
  // Custom Exceptions for errors during the solve
  //----------------------------------------------------------------------------
  /**
   * @brief Exception thrown when an iterative solver fails to converge.
   */
  class ConvergenceError : public std::runtime_error {
  public:
    /**
     * Exception class constructor.
     */
    ConvergenceError()
      : std::runtime_error("Iterative solver failed to converge.") {}
  };
}

#endif //STAG_TEST_SOLVE_H
