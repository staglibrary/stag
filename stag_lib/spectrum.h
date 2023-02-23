//
// Methods relating to the spectrum of the graph matrices. Includes methods for
// computing the eigenvectors and eigenvalues of the graph matrices.
//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file spectrum.h
 * \brief Methods for computing eigenvalues and eigenvectors of sparse matrices.
 */

#ifndef STAG_TEST_SPECTRUM_H
#define STAG_TEST_SPECTRUM_H

#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <tuple>
#include "graph.h"


namespace stag {
  /**
   * \cond
   * The spectra operation for multiplying the SprsMat type
   *
   * Undocumented by doxygen.
   */
  typedef Spectra::SparseSymMatProd<double, Eigen::Upper, Eigen::ColMajor, stag_int> SprsMatOp;
  /**
   * \endcond
   */

  /**
   * A convenience type for returning eigenvalues and eigenvectors together.
   */
  typedef std::tuple<Eigen::VectorXd, Eigen::MatrixXd> EigenSystem;

  /**
   * Compute the eigenvalues and eigenvectors of a given matrix.
   *
   * By default, this will compute the eigenvalues of smallest magnitude.
   * This default can be overridden by the sort parameter which takes a
   * `Spectra::SortRule` object. This is likely to be one of the following.
   *   - `Spectra::SortRule::SmallestMagn` will return the eigenvalues with smallest magnitude
   *   - `Spectra::SortRule::LargestMagn` will return the eigenvalues with largest magnitude
   *
   * The following example demonstrates how to compute the 3 largest eigenvectors
   * and eigenvalues of a cycle graph.
   *
   * \code{.cpp}
   *     #include <Spectra/SymEigsSolver.h>
   *     #include <stag/graph.h>
   *     #include <stag/spectrum.h>
   *
   *     int main() {
   *       // Create a cycle graph
   *       stag::Graph myGraph = stag::cycle_graph(10);
   *
   *       // Extract the normalised laplacian matrix
   *       const SprsMat* lap = testGraph.normalised_laplacian();
   *
   *       // Compute a few eigenvalues and eigenvectors
   *       stag_int k = 3;
   *       stag::EigenSystem eigensystem = stag::compute_eigensystem(
   *           lap, k, Spectra::SortRule::LargestMagn);
   *
   *       Eigen::VectorXd eigenvalues = get<0>(eigensystem);
   *       Eigen::MatrixXd eigenvectors = get<1>(eigensystem);
   *
   *       return 0;
   *     }
   * \endcode
   *
   * @param mat the matrix on which to operate
   * @param num the number of eigenvalues and eigenvectors to compute
   * @param sort (optional) a Spectra::SortRule indicating which eigenvectors to calculate
   * @returns a stag::EigenSystem object containing the computed eigenvalues and eigenvectors
   * @throws std::runtime_error if the eigenvalue calculation does not converge
   */
  stag::EigenSystem compute_eigensystem(const SprsMat* mat, stag_int num,
                                        Spectra::SortRule sort);

  /**
   * \overload
   */
  stag::EigenSystem compute_eigensystem(const SprsMat* mat, stag_int num);

  /**
   * Compute the eigenvalues of a given matrix.
   *
   * By default, this will compute the eigenvalues of smallest magnitude.
   * This default can be overridden by the sort parameter which takes a `Spectra::SortRule`
   * object. This is likely to be one of the following.
   *   - `Spectra::SortRule::SmallestMagn` will return the eigenvalues with smallest magnitude
   *   - `Spectra::SortRule::LargestMagn` will return the eigenvalues with largest magnitude
   *
   * If you would like to calculate the eigenvectors and eigenvalues together, then
   * you should instead use stag::compute_eigensystem.
   *
   * @param mat the matrix on which to operate
   * @param num the number of eigenvalues to compute
   * @param sort (optional) a Spectra::SortRule indicating which eigenvalues to calculate
   * @returns an Eigen::VectorXd object containing the computed eigenvalues
   * @throws std::runtime_error if the eigenvalue calculation does not converge
   */
  Eigen::VectorXd compute_eigenvalues(const SprsMat* mat, stag_int num,
                                      Spectra::SortRule sort);

  /**
   * \overload
   */
  Eigen::VectorXd compute_eigenvalues(const SprsMat* mat, stag_int num);

  /**
   * Compute the eigenvectors of a given matrix.
   *
   * By default, this will compute the eigenvectors corresponding to the
   * eigenvalues of smallest magnitude.
   * This default can be overridden by the sort parameter which takes a `Spectra::SortRule`
   * object. This is likely to be one of the following.
   *   - `Spectra::SortRule::SmallestMagn` will return the eigenvalues with smallest magnitude
   *   - `Spectra::SortRule::LargestMagn` will return the eigenvalues with largest magnitude
   *
   * If you would like to calculate the eigenvectors and eigenvalues together, then
   * you should instead use stag::compute_eigensystem.
   *
   * @param mat the matrix on which to operate
   * @param num the number of eigenvectors to compute
   * @param sort (optional) a Spectra::SortRule indicating which eigenvectors to calculate
   * @returns an Eigen::MatrixXd object containing the computed eigenvectors
   * @throws std::runtime_error if the eigenvector calculation does not converge
   */
  Eigen::MatrixXd compute_eigenvectors(const SprsMat* mat, stag_int num,
                                      Spectra::SortRule sort);

  /**
   * \overload
   */
  Eigen::MatrixXd compute_eigenvectors(const SprsMat* mat, stag_int num);

  /**
   * Apply the power method to compute the dominant eigenvector of a matrix.
   *
   * Given a matrix \f$M\f$, an initial vector \f$v_0\f$, and a number of
   * iterations \f$t\f$, the power method calculates the vector
   *
   * \f[
   *    v_t = M^t v_0,
   * \f]
   *
   * which is close to the eigenvector of \f$M\f$ with largest eigenvalue.
   *
   * The running time of the power method is \f$O(t \cdot \mathrm{nnz}(M))\f$, where
   * \f$\mathrm{nnz}(M)\f$ is the number of non-zero elements in the matrix \f$M\f$.
   *
   * @param mat the matrix \f$M\f$ on which to operate.
   * @param num_iterations (optional) the number of iterations of the power
   *                       method to apply. It this argument is omitted,
   *                       \f$O(\log(n))\f$ iterations are used which results
   *                       in a vector whose Rayleigh quotient is a \f$(1 - \epsilon)\f$
   *                       approximation of the dominant eigenvalue.
   * @param initial_vector (optional) the initial vector to use for the power
   *                       iteration. If this argument is omitted, a random unit
   *                       vector will be used.
   * @return the vector \f$v_t\f$ computed by repeated multiplication with \f$M\f$.
   */
  Eigen::VectorXd power_method(const SprsMat* mat, stag_int num_iterations,
                               Eigen::VectorXd initial_vector);

  /**
   * \overload
   */
  Eigen::VectorXd power_method(const SprsMat* mat, stag_int num_iterations);

  /**
   * \overload
   */
  Eigen::VectorXd power_method(const SprsMat* mat,
                               Eigen::VectorXd initial_vector);

  /**
   * \overload
   */
  Eigen::VectorXd power_method(const SprsMat* mat);

  /**
   * Compute the Rayleigh quotient of the given vector and matrix.
   *
   * Given a matrix \f$M\f$, the Rayleigh quotient of vector \f$v\f$ is
   *
   * \f[
   *    R(M, v) = \frac{v^\top M v}{v^\top v}.
   * \f]
   *
   * @param mat a sparse matrix \f$M \in \mathbb{R}^{n \times n}\f$.
   * @param vec a vector \f$v \in \mathbb{R}^n\f$.
   * @return the Rayleigh quotient \f$R(M, v)\f$.
   */
  double rayleigh_quotient(const SprsMat* mat, Eigen::VectorXd& vec);
}


#endif //STAG_TEST_SPECTRUM_H
