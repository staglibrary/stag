/**
 * Definitions of utility methods for dealing with sparse matrices.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#ifndef STAG_TEST_UTILITY_H
#define STAG_TEST_UTILITY_H

#include "graph.h"

namespace stag {

  /**
   * Given a sparse matrix, return the values vector, compatible with the CSR
   * format of other libraries.
   *
   * @param matrix
   * @return
   */
  std::vector<double> sprsMatValues(const SprsMat *matrix);

  /**
   * Given a sparse matrix, return the InnerIndices vector, compatible with the
   * CSR format of other libraries.
   *
   * @param matrix
   * @return
   */
  std::vector<stag_int> sprsMatInnerIndices(const SprsMat *matrix);

  /**
   * Given a sparse matrix, return the OuterStarts vector, compatible with the
   * CSR format of other libraries.
   *
   * @param matrix
   * @return
   */
  std::vector<stag_int> sprsMatOuterStarts(const SprsMat *matrix);

  /**
   * Given a sparse 'matrix' with only one column, convert it to a dense vector.
   *
   * @param matrix - the sparse vector to convert
   * @param n (optional) - the dimension of the dense vector to construct
   * @return a vector
   */
   std::vector<double> sprsMatToVec(const SprsMat *matrix, stag_int n);
   std::vector<double> sprsMatToVec(const SprsMat *matrix);

  /**
   * Check whether a sparse matrix is symmetric.
   */
  bool isSymmetric(const SprsMat *matrix);
}

#endif //STAG_TEST_UTILITY_H
