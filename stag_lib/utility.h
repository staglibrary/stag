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
  std::vector<int> sprsMatInnerIndices(const SprsMat *matrix);

  /**
   * Given a sparse matrix, return the OuterStarts vector, compatible with the
   * CSR format of other libraries.
   *
   * @param matrix
   * @return
   */
  std::vector<int> sprsMatOuterStarts(const SprsMat *matrix);
}

#endif //STAG_TEST_UTILITY_H
