/**
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <iterator>
#include "utility.h"

std::vector<stag_int> stag::sprsMatInnerIndices(const SprsMat *matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const stag_int *indexPtr = matrix->innerIndexPtr();
  stag_int nonZeros = matrix->nonZeros();
  return {indexPtr, indexPtr + nonZeros};
}

std::vector<stag_int> stag::sprsMatOuterStarts(const SprsMat* matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const stag_int *indexPtr = matrix->outerIndexPtr();
  stag_int outerSize = matrix->outerSize();
  return {indexPtr, indexPtr + outerSize + 1};
}

std::vector<double> stag::sprsMatValues(const SprsMat* matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const double *valuePtr = matrix->valuePtr();
  stag_int nonZeros = matrix->nonZeros();
  return {valuePtr, valuePtr + nonZeros};
}

std::vector<double> stag::sprsMatToVec(const SprsMat* matrix) {
  // If the number of dimensions is not given, use the dimension of the sparse
  // matrix.
  return stag::sprsMatToVec(matrix, matrix->rows());
}

std::vector<double> stag::sprsMatToVec(const SprsMat* matrix, stag_int n) {
  // Initialise the solution vector.
  std::vector<double> dense_vec;

  for (stag_int i = 0; i < n; i++) {
    if (i < matrix->rows()) {
      // Get the i-th entry of the sparse matrix
      dense_vec.push_back(matrix->coeff(i, 0));
    } else {
      // If the sparse matrix is not long enough, fill the vector with 0s.
      dense_vec.push_back(0);
    }
  }
  return dense_vec;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "ArgumentSelectionDefects"
bool stag::isSymmetric(const SprsMat *matrix) {
  // Iterate through the non-zero elements in the matrix
  for (int k = 0; k < matrix->outerSize(); ++k) {
    for (SprsMat::InnerIterator it(*matrix, k); it; ++it) {
      // If the value in the symmetrically opposite position is not the same,
      // then return false.
      if (it.value() != matrix->coeff(it.col(), it.row())) {
        return false;
      }
    }
  }

  // We didn't find any symmetrically opposite coefficients with different
  // values, and so this matrix is symmetric.
  return true;
}
#pragma clang diagnostic pop
