/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/
// Standard C++ libraries
#include <stdexcept>
#include <algorithm>
#include <random>
#include <utility>

// Other libraries
#include <Spectra/SymEigsSolver.h>

// STAG modules
#include "spectrum.h"


/**
 * Compute the eigensystem of a matrix, beginning with the largest eigenvalues.
 */
stag::EigenSystem compute_eigensystem_largestmag(const SprsMat* mat, StagInt num) {
  if (num < 1 || num >= mat->rows()) {
    throw std::invalid_argument("Number of computed eigenvectors must be between 1 and n - 1.");
  }

  stag::SprsMatProdOp op(*mat);

  // Construct eigen solver object, requesting the smallest k eigenvalues
  long ncv = std::min<StagInt>(10 * num, mat->rows());
  Spectra::SymEigsSolver<stag::SprsMatProdOp> eigs(op, num, ncv);

  // Initialize and compute
  eigs.init();
  eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10, Spectra::SortRule::LargestMagn);

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

/**
 * Compute the eigensystem of a matrix, beginning with the smallest eigenvalues.
 */
stag::EigenSystem compute_eigensystem_smallestmag(const SprsMat* mat, StagInt num) {
  if (num < 1 || num >= mat->rows()) {
    throw std::invalid_argument("Number of computed eigenvectors must be between 1 and n - 1.");
  }

  stag::SprsMatShiftSolveOp op(*mat);

  // Construct eigen solver object, requesting the smallest k eigenvalues
  long ncv = std::min<StagInt>(10 * num, mat->rows());
  Spectra::SymEigsShiftSolver<stag::SprsMatShiftSolveOp> eigs(op, num, ncv, -1e-6);

  // Initialize and compute
  // When computing eigenvalues with the smallest magnitude, we use the
  // shift-and-invert solver.
  eigs.init();
  eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10, Spectra::SortRule::SmallestMagn);

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

stag::EigenSystem stag::compute_eigensystem(
    const SprsMat* mat, StagInt num, Spectra::SortRule sort) {
  if (num < 1 || num >= mat->rows()) {
    throw std::invalid_argument("Number of computed eigenvectors must be between 1 and n - 1.");
  }

  if (sort == Spectra::SortRule::LargestMagn) {
    return compute_eigensystem_largestmag(mat, num);
  } else {
    return compute_eigensystem_smallestmag(mat, num);
  }
}

stag::EigenSystem stag::compute_eigensystem(const SprsMat *mat, StagInt num) {
  return stag::compute_eigensystem(mat, num, Spectra::SortRule::SmallestMagn);
}

Eigen::VectorXd stag::compute_eigenvalues(const SprsMat *mat, StagInt num,
                                          Spectra::SortRule sort) {
  return get<0>(stag::compute_eigensystem(mat, num, sort));
}

Eigen::VectorXd stag::compute_eigenvalues(const SprsMat *mat, StagInt num) {
  return stag::compute_eigenvalues(mat, num, Spectra::SortRule::SmallestMagn);
}

Eigen::MatrixXd stag::compute_eigenvectors(const SprsMat *mat, StagInt num,
                                           Spectra::SortRule sort) {
  return get<1>(stag::compute_eigensystem(mat, num, sort));
}

Eigen::MatrixXd stag::compute_eigenvectors(const SprsMat *mat, StagInt num) {
  return stag::compute_eigenvectors(mat, num, Spectra::SortRule::SmallestMagn);
}

/**
 * Generate a random unit vector with the given dimension.
 *
 * @param dimension the dimension of the random vector.
 * @return a random unit vector in \f$\mathbb{R}^d\f$.
 */
Eigen::VectorXd random_unit_vector(StagInt dimension) {
  // We can generate a random unit vector by generating a vector with
  // each entry drawn from the standard normal distribution and then
  // normalising the resulting vector.
  std::random_device dev;
  std::mt19937 prng(dev());
  std::normal_distribution<StagReal> gaussian_distribution(0, 1);

  Eigen::VectorXd random_vector(dimension);
  for (auto i = 0; i < dimension; i++) {
    random_vector(i) = gaussian_distribution(prng);
  }

  random_vector.normalize();
  return random_vector;
}

Eigen::VectorXd stag::power_method(const SprsMat *mat, StagInt num_iterations,
                                   Eigen::VectorXd initial_vector) {
  if (num_iterations < 0) throw std::invalid_argument("Number of iterations must be non-negative.");
  if (mat->rows() != mat->cols()) throw std::invalid_argument("Matrix must be square.");
  if (initial_vector.size() != mat->rows()) throw std::invalid_argument("Vector and matrix must have the same dimension");

  for (auto t = 0; t < num_iterations; t++) {
    initial_vector = *mat * initial_vector;
    initial_vector.normalize();
  }

  return initial_vector;
}

Eigen::VectorXd stag::power_method(const SprsMat *mat,
                                   Eigen::VectorXd initial_vector) {
  StagInt n = mat->rows();
  StagInt t = 10 * ((int) ceil(log((StagReal) n)));
  return stag::power_method(mat, t, std::move(initial_vector));
}

Eigen::VectorXd stag::power_method(const SprsMat *mat, StagInt num_iterations) {
  return stag::power_method(mat, num_iterations, random_unit_vector(mat->rows()));
}

Eigen::VectorXd stag::power_method(const SprsMat *mat) {
  return stag::power_method(mat, random_unit_vector(mat->rows()));
}

StagReal stag::rayleigh_quotient(const SprsMat *mat, Eigen::VectorXd& vec) {
  if (mat->rows() != mat->cols()) throw std::invalid_argument("Matrix must be square.");
  if (vec.size() != mat->rows()) throw std::invalid_argument("Vector and matrix must have the same dimension");
  if (vec.norm() == 0) throw std::invalid_argument("Vector with norm 0 had undefined Rayleigh quotient");

  StagReal numerator = vec.transpose() * *mat * vec;
  StagReal denominator = pow(vec.norm(), 2);
  return numerator / denominator;
}
