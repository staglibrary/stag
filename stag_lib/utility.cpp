/**
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <iterator>
#include "utility.h"

std::vector<int> stag::sprsMatInnerIndices(const SprsMat *matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const int *indexPtr = matrix->innerIndexPtr();
  long nonZeros = matrix->nonZeros();
  return {indexPtr, indexPtr + nonZeros};
}

std::vector<int> stag::sprsMatOuterStarts(const SprsMat *matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const int *indexPtr = matrix->outerIndexPtr();
  long outerSize = matrix->outerSize();
  return {indexPtr, indexPtr + outerSize + 1};
}

std::vector<double> stag::sprsMatValues(const SprsMat *matrix) {
  // Make sure that the given matrix is compressed
  assert(matrix->isCompressed());

  // Return the required indices vector
  const double *valuePtr = matrix->valuePtr();
  long nonZeros = matrix->nonZeros();
  return {valuePtr, valuePtr + nonZeros};
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
