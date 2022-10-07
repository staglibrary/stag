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
