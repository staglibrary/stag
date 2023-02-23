//
// This file is provided as part of the STAG library and released under the MIT
// license.
//
#include "spectrum.h"

#include <Spectra/SymEigsSolver.h>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <utility>


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

/**
 * Generate a random unit vector with the given dimension.
 *
 * @param dimension the dimension of the random vector.
 * @return a random unit vector in \f$\mathbb{R}^d\f$.
 */
Eigen::VectorXd random_unit_vector(stag_int dimension) {
  // We can generate a random unit vector by generating a vector with
  // each entry drawn from the standard normal distribution and then
  // normalising the resulting vector.
  std::random_device dev;
  std::mt19937 prng(dev());
  std::normal_distribution<double> gaussian_distribution(0, 1);

  Eigen::VectorXd random_vector(dimension);
  for (auto i = 0; i < dimension; i++) {
    random_vector(i) = gaussian_distribution(prng);
  }

  random_vector.normalize();
  return random_vector;
}

Eigen::VectorXd stag::power_method(const SprsMat *mat, stag_int num_iterations,
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
  stag_int n = mat->rows();
  stag_int t = 10 * ((int) ceil(log((double) n)));
  return stag::power_method(mat, t, std::move(initial_vector));
}

Eigen::VectorXd stag::power_method(const SprsMat *mat, stag_int num_iterations) {
  return stag::power_method(mat, num_iterations, random_unit_vector(mat->rows()));
}

Eigen::VectorXd stag::power_method(const SprsMat *mat) {
  return stag::power_method(mat, random_unit_vector(mat->rows()));
}

double stag::rayleigh_quotient(const SprsMat *mat, Eigen::VectorXd& vec) {
  if (mat->rows() != mat->cols()) throw std::invalid_argument("Matrix must be square.");
  if (vec.size() != mat->rows()) throw std::invalid_argument("Vector and matrix must have the same dimension");
  if (vec.norm() == 0) throw std::invalid_argument("Vector with norm 0 had undefined Rayleigh quotient");

  double numerator = vec.transpose() * *mat * vec;
  double denominator = pow(vec.norm(), 2);
  return numerator / denominator;
}
