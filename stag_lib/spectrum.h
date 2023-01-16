//
// Methods relating to the spectrum of the graph matrices. Includes methods for
// computing the eigenvectors and eigenvalues of the graph matrices.
//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#ifndef STAG_TEST_SPECTRUM_H
#define STAG_TEST_SPECTRUM_H

#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <tuple>
#include "graph.h"


namespace stag {
  /*
   * The spectra operation for multiplying the SprsMat type
   */
  typedef Spectra::SparseSymMatProd<double, Eigen::Upper, Eigen::ColMajor, stag_int> SprsMatOp;

  /*
   * An eigensystem is defined to be a tuple of eigenvalues and eigenvectors
   */
  typedef std::tuple<Eigen::VectorXd, Eigen::MatrixXd> EigenSystem;

  /*
   * Compute the eigenvalues and eigenvectors of a given matrix.
   * By default, this will compute the eigenvalues of smallest magnitude.
   * This default can be overriden by the sort parameter which takes a sortrule
   * from the Spectra library. Likely to be one of:
   *   - Spectra::SortRule::SmallestMagn (default),
   *   - Spectra::SortRule::LargestMagn.
   */
  stag::EigenSystem compute_eigensystem(const SprsMat* mat, stag_int num,
                                        Spectra::SortRule sort);
  stag::EigenSystem compute_eigensystem(const SprsMat* mat, stag_int num);

  /*
   * Compute the eigenvalues of the given matrix.
   * By default, this will compute the eigenvalues of smallest magnitude.
   * This default can be overriden by the sort parameter which takes a sortrule
   * from the Spectra library. Likely to be one of:
   *   - Spectra::SortRule::SmallestMagn (default),
   *   - Spectra::SortRule::LargestMagn.
   */
  Eigen::VectorXd compute_eigenvalues(const SprsMat* mat, stag_int num,
                                      Spectra::SortRule sort);
  Eigen::VectorXd compute_eigenvalues(const SprsMat* mat, stag_int num);

  /*
   * Compute the eigenvectors of the given matrix.
   * By default, this will compute the eigenvectors with eigenvalues
   * of smallest magnitude. This default can be overriden by the sort parameter
   * which takes a sortrule from the Spectra library. Likely to be one of:
   *   - Spectra::SortRule::SmallestMagn (default),
   *   - Spectra::SortRule::LargestMagn.
   */
  Eigen::MatrixXd compute_eigenvectors(const SprsMat* mat, stag_int num,
                                      Spectra::SortRule sort);
  Eigen::MatrixXd compute_eigenvectors(const SprsMat* mat, stag_int num);
}


#endif //STAG_TEST_SPECTRUM_H
