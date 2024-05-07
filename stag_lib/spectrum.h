/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/

/**
 * @file spectrum.h
 * \brief Methods for computing eigenvalues and eigenvectors of graph matrices.
 *
 * For applications in spectral graph theory, we are interested in the spectrum
 * of the Graph adjacency and Laplacian matrix spectrums.
 *
 * Moreover, for each matrix, we are usually interested in the extreme
 * eigenvalues at one end of the spectrum or the other.
 *
 * As such, this module provides methods for computing the extreme eigenvalues
 * and eigenvectors of normalised or unnormalised adjacency or Laplacian matrix
 * of some graph.
 *
 * We limit ourselves to these cases as this is the most common application in
 * spectral Graph theory, and these matrices are 'well-behaved' for the solver
 * algorithms. If you'd like to compute the spectrum of a more general matrix,
 * you could consider using the Spectra C++ library directly.
 */

#ifndef STAG_TEST_SPECTRUM_H
#define STAG_TEST_SPECTRUM_H

// Standard C++ libraries
#include <tuple>

// Other libraries
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

// STAG modules
#include "graph.h"

namespace stag {
  /**
   * \cond
   * The spectra operation for multiplying the SprsMat type
   *
   * Undocumented by doxygen.
   */
  typedef Spectra::SparseSymMatProd<StagReal, Eigen::Upper, Eigen::ColMajor, StagInt> SprsMatProdOp;
  typedef Spectra::SparseSymShiftSolve<StagReal, Eigen::Upper, Eigen::ColMajor, StagInt> SprsMatShiftSolveOp;
  /**
   * \endcond
   */

  /**
   * A convenience type for returning eigenvalues and eigenvectors together.
   */
  typedef std::tuple<Eigen::VectorXd, Eigen::MatrixXd> EigenSystem;

  /**
   * When computing eigenvectors and eigenvalues, we always compute the values
   * at one end of the spectrum or the other.
   *
   * The EigenSortRule is used to specify whether to compute the largest or
   * smallest eigenvalues.
   */
  enum EigenSortRule {Largest, Smallest};

  /**
   * When computing eigenvectors and eigenvalues, these values are used to
   * specify which Graph matrix we are using to compute the spectrum.
   */
  enum GraphMatrix {Adjacency, Laplacian, NormalisedLaplacian};

  /**
   * Compute the eigenvalues and eigenvectors of a graph matrix.
   *
   * Computes a given number of eigenvectors at one end of the spectrum of a
   * graph matrix.
   *
   * The following example demonstrates how to compute the 3 largest eigenvectors
   * and eigenvalues of the normalised laplacian of a cycle graph.
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
   *       // Compute a few eigenvalues and eigenvectors
   *       StagInt k = 3;
   *       stag::EigenSystem eigensystem = stag::compute_eigensystem(
   *           myGraph, stag::GraphMatrix::NormalisedLaplacian,
   *           k, stag::EigenSortRule::Largest);
   *
   *       Eigen::VectorXd eigenvalues = get<0>(eigensystem);
   *       Eigen::MatrixXd eigenvectors = get<1>(eigensystem);
   *
   *       return 0;
   *     }
   * \endcode
   *
   * @param g the graph on which to operate
   * @param mat which graph matrix to use to compute the spectrum
   * @param num_eigs the number of eigenvalues and eigenvectors to compute
   * @param which a stag::EigenSortRule value indicating which eigenvectors to return
   * @returns a stag::EigenSystem object containing the computed eigenvalues and eigenvectors
   * @throws std::runtime_error if the eigenvalue calculation does not converge
   */
  stag::EigenSystem compute_eigensystem(stag::Graph* g,
                                        stag::GraphMatrix mat,
                                        StagInt num_eigs,
                                        stag::EigenSortRule which);

  /**
   * Compute the eigenvectors of a graph matrix.
   *
   * If you would like to calculate the eigenvectors and eigenvalues together, then
   * you should instead use stag::compute_eigensystem.
   *
   * @param g the graph on which to operate
   * @param mat which graph matrix to use to compute the spectrum
   * @param num_eigs the number of eigenvectors to compute
   * @param which a stag::EigenSortRule value indicating which eigenvectors to return
   * @returns an Eigen::MatrixXd object containing the computed eigenvectors
   * @throws std::runtime_error if the eigenvector calculation does not converge
   */
  Eigen::MatrixXd compute_eigenvectors(stag::Graph* g,
                                       stag::GraphMatrix mat,
                                       StagInt num_eigs,
                                       stag::EigenSortRule which);

  /**
   * Compute the eigenvalues of a graph matrix.
   *
   * Computes a given number of eigenvalues at one end of the spectrum of a
   * graph matrix.
   *
   * If you would like to calculate the eigenvectors and eigenvalues together, then
   * you should instead use stag::compute_eigensystem.
   *
   * @param g the graph on which to operate
   * @param mat which graph matrix to use to compute the spectrum
   * @param num_eigs the number of eigenvalues to compute
   * @param which a stag::EigenSortRule value indicating which eigenvalues to return
   * @returns an Eigen::VectorXd object containing the computed eigenvalues
   * @throws std::runtime_error if the eigenvalue calculation does not converge
   */
  Eigen::VectorXd compute_eigenvalues(stag::Graph* g,
                                      stag::GraphMatrix mat,
                                      StagInt num_eigs,
                                      stag::EigenSortRule which);

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
  Eigen::VectorXd power_method(const SprsMat* mat, StagInt num_iterations,
                               Eigen::VectorXd initial_vector);

  /**
   * \overload
   */
  Eigen::VectorXd power_method(const SprsMat* mat, StagInt num_iterations);

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
  StagReal rayleigh_quotient(const SprsMat* mat, Eigen::VectorXd& vec);
}


#endif //STAG_TEST_SPECTRUM_H
