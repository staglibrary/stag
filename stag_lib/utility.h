//
// Definitions of utility methods for dealing with sparse matrices.
//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file utility.h
 * \brief Helper methods for processing sparse matrices
 */

#ifndef STAG_TEST_UTILITY_H
#define STAG_TEST_UTILITY_H

#include <iostream>

#include "graph.h"

namespace stag {

  /**
   * Given a sparse matrix, return the values vector, compatible with the CSC
   * format of other libraries.
   */
  std::vector<double> sprsMatValues(const SprsMat *matrix);

  /**
   * Given a sparse matrix, return the InnerIndices vector, compatible with the
   * CSC format of other libraries.
   */
  std::vector<stag_int> sprsMatInnerIndices(const SprsMat *matrix);

  /**
   * Given a sparse matrix, return the OuterStarts vector, compatible with the
   * CSC format of other libraries.
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

   /**
    * \overload
    */
   std::vector<double> sprsMatToVec(const SprsMat *matrix);

   /**
    * Construct a sparse matrix from the CSC data vectors.
    *
    * For documentation on the format of the data vectors, please see the
    * documentation for the Eigen sparse matrix object.
    *
    * This method does not perform any error checking on the provided
    * vectors. The caller is responsible for ensuring that the provided data
    * vectors are well-formed.
    */
   SprsMat sprsMatFromVectors(std::vector<stag_int>& column_starts,
                              std::vector<stag_int>& row_indices,
                              std::vector<double>& values);

   /**
    * Add two vectors together element-wise.
    */
   template <typename T>
   std::vector<T> addVectors(std::vector<T>& v1, std::vector<T>& v2) {
     auto length = (stag_int) std::max(v1.size(), v2.size());
     std::vector<T> ans;
     T this_entry;

     for (stag_int i = 0; i < length; i++) {
       this_entry = 0;
       if (v1.size() > i) this_entry += v1.at(i);
       if (v2.size() > i) this_entry += v2.at(i);
       ans.push_back(this_entry);
     }

     return ans;
   }


  /**
   * Check whether a sparse matrix is symmetric.
   */
  bool isSymmetric(const SprsMat *matrix);

  /**
   * \cond
   * Do not document the stdErrVec or safeGetline methods
   */

  /**
   * Print a vector to stderr.
   */
  template <typename T>
  void stdErrVec(std::vector<T>& vec){
    for (auto i : vec) {
      std::cerr << i << ", ";
    }
    std::cerr << std::endl;
  }

  /**
   * Get the next line from an input stream, while safely handling all types of
   * line endings (CR, LF, CRLF).
   *
   * @param is the input stream to process
   * @param t the string variable in which to store the returned line
   */
  std::istream& safeGetline(std::istream& is, std::string& t);

  /**
   * \endcond
   */
}

#endif //STAG_TEST_UTILITY_H
