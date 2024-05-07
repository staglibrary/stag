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
 *
 * Add offset to the eigenvalues before returning them.
 */
stag::EigenSystem compute_eigensystem_largestmag(
    const SprsMat* mat, StagInt num, StagReal offset, bool invert) {
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

  // Add the offset to the eigenvalues
  if (offset != 0 || invert) {
    for (auto i = 0; i < eigenvalues.rows(); i++) {
      eigenvalues.coeffRef(i) += offset;

      if (invert) eigenvalues.coeffRef(i) *= -1;
    }
  }

  return {eigenvalues, eigenvectors};
}

stag::EigenSystem stag::compute_eigensystem(
    stag::Graph* g, stag::GraphMatrix mat, StagInt num, stag::EigenSortRule which) {
  if (num < 1 || num >= g->number_of_vertices()) {
    throw std::invalid_argument("Number of computed eigenvectors must be between 1 and n - 1.");
  }

  // Get the maximum degree of the graph
  StagReal max_degree = 0;
  for (auto i = 0; i < g->number_of_vertices(); i++) {
    max_degree = MAX(g->degree(i), max_degree);
  }

  switch (mat) {
    case stag::GraphMatrix::Adjacency:
      if (which == stag::EigenSortRule::Largest) {
        // We will find the maximum eigenvalues of A + d_max I.
        SprsMat adjusted_adjacency = *g->adjacency();
        for (auto i = 0; i < g->number_of_vertices(); i++) {
          adjusted_adjacency.coeffRef(i, i) += max_degree;
        }
        return compute_eigensystem_largestmag(
            &adjusted_adjacency, num, -max_degree, false);
      } else {
        // We will find the maximum eigenvalues of -A + d_max I.
        SprsMat adjusted_adjacency = - (*g->adjacency());
        for (auto i = 0; i < g->number_of_vertices(); i++) {
          adjusted_adjacency.coeffRef(i, i) += max_degree;
        }
        return compute_eigensystem_largestmag(
            &adjusted_adjacency, num, -max_degree, true);
      }
      break;
    case stag::GraphMatrix::Laplacian:
      if (which == stag::EigenSortRule::Largest) {
        // We can just compute the largest eigenvalues directly.
        return compute_eigensystem_largestmag(g->laplacian(), num, 0, false);
      } else {
        // We will find the maximum eigenvalues of (2 d_max I - L).
        SprsMat adjusted_laplacian = - (*g->laplacian());
        for (auto i = 0; i < g->number_of_vertices(); i++) {
          adjusted_laplacian.coeffRef(i, i) += 2 * max_degree;
        }
        return compute_eigensystem_largestmag(
            &adjusted_laplacian, num, -(2 * max_degree), true);
      }
      break;
    case stag::GraphMatrix::NormalisedLaplacian:
      if (which == stag::EigenSortRule::Largest) {
        // We can just compute the largest eigenvalues directly.
        return compute_eigensystem_largestmag(
            g->normalised_laplacian(), num, 0, false);
      } else {
        // We will find the maximum eigenvalues of (2 I - L).
        SprsMat adjusted_laplacian = - (*g->normalised_laplacian());
        for (auto i = 0; i < g->number_of_vertices(); i++) {
          adjusted_laplacian.coeffRef(i, i) += 2;
        }
        return compute_eigensystem_largestmag(
            &adjusted_laplacian, num, -2, true);
      }
      break;
    default:
      break;
  }

  // Should never get here
  assert(false);
  throw std::runtime_error("Failed to compute eigenvectors.");
}

Eigen::MatrixXd stag::compute_eigenvectors(
    stag::Graph* g, stag::GraphMatrix mat, StagInt num, stag::EigenSortRule which) {
  return get<1>(stag::compute_eigensystem(g, mat, num, which));
}

Eigen::VectorXd stag::compute_eigenvalues(
    stag::Graph* g, stag::GraphMatrix mat, StagInt num, stag::EigenSortRule which) {
  return get<0>(stag::compute_eigensystem(g, mat, num, which));
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
