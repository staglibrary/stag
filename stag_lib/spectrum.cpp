//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#include "spectrum.h"

#include <Spectra/SymEigsSolver.h>
#include <stdexcept>
#include <algorithm>


stag::EigenSystem stag::compute_eigensystem(
    const SprsMat* mat, stag_int num, Spectra::SortRule sort) {
  if (num < 1 || num >= mat->rows()) {
    throw std::invalid_argument("Number of computed eigenvectors must be between 1 and n - 1.");
  }

  stag::SprsMatOp op(*mat);

  // Construct eigen solver object, requesting the smallest k eigenvalues
  long ncv = std::min<stag_int>(20 * num, mat->rows());
  Spectra::SymEigsSolver<SprsMatOp> eigs(op, num, ncv);

  // Initialize and compute
  eigs.init();
  eigs.compute(sort);

  // Ensure that the calculation has converged
  if (eigs.info() != Spectra::CompInfo::Successful) {
    throw std::runtime_error("Eigenvalue calculation failed to converge.");
  }

  // Retrieve results
  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd eigenvectors;
  eigenvalues = eigs.eigenvalues();
  eigenvectors = eigs.eigenvectors();

  return {eigenvalues, eigenvectors};
}

stag::EigenSystem stag::compute_eigensystem(const SprsMat *mat, stag_int num) {
  return stag::compute_eigensystem(mat, num, Spectra::SortRule::SmallestMagn);
}

Eigen::VectorXd stag::compute_eigenvalues(const SprsMat *mat, stag_int num,
                                          Spectra::SortRule sort) {
  return get<0>(stag::compute_eigensystem(mat, num, sort));
}

Eigen::VectorXd stag::compute_eigenvalues(const SprsMat *mat, stag_int num) {
  return stag::compute_eigenvalues(mat, num, Spectra::SortRule::SmallestMagn);
}

Eigen::MatrixXd stag::compute_eigenvectors(const SprsMat *mat, stag_int num,
                                           Spectra::SortRule sort) {
  return get<1>(stag::compute_eigensystem(mat, num, sort));
}

Eigen::MatrixXd stag::compute_eigenvectors(const SprsMat *mat, stag_int num) {
  return stag::compute_eigenvectors(mat, num, Spectra::SortRule::SmallestMagn);
}
