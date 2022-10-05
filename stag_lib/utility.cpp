#include <iterator>
#include "utility.h"

std::vector<int> stag::sprsMatInnerIndices(SprsMat &matrix) {
  // Make sure that the given matrix is compressed
  matrix.makeCompressed();

  // Return the required indices vector
  int *indexPtr = matrix.innerIndexPtr();
  long nonZeros = matrix.nonZeros();
  return {indexPtr, indexPtr + nonZeros};
}

std::vector<int> stag::sprsMatOuterStarts(SprsMat &matrix) {
  // Make sure that the given matrix is compressed
  matrix.makeCompressed();

  // Return the required indices vector
  int *indexPtr = matrix.outerIndexPtr();
  long outerSize = matrix.outerSize();
  return {indexPtr, indexPtr + outerSize + 1};
}

std::vector<double> stag::sprsMatValues(SprsMat &matrix) {
  // Make sure that the given matrix is compressed
  matrix.makeCompressed();

  // Return the required indices vector
  double *valuePtr = matrix.valuePtr();
  long nonZeros = matrix.nonZeros();
  return {valuePtr, valuePtr + nonZeros};
}
