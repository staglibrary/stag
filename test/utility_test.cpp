/**
 * Tests for the utility.h header. Includes helper methods for dealing with
 * sparse matrices.
 *
 * This file is provided as part of the STAG library and released under the GPL
 * license.
 */
#include <gtest/gtest.h>
#include <utility.h>

TEST(UtilityTest, IsSymmetric) {
  // Construct a symmetric matrix.
  std::vector<StagInt> rowStarts = {0, 2, 4, 7, 8};
  std::vector<StagInt> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};
  SprsMat matrix = Eigen::Map<SprsMat>(4, 4, 8,
                                       rowStarts.data(),
                                       colIndices.data(),
                                       values.data());

  // Check that the isSymmetric method gives the correct answer.
  EXPECT_EQ(stag::isSymmetric(&matrix), true);

  // Construct a matrix which is not symmetric.
  rowStarts = {0, 2, 4, 7, 8};
  colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  values = {2, 3.3333, 2, 6, 3, 6, 1, 1};
  matrix = Eigen::Map<SprsMat>(4, 4, 8,
                                       rowStarts.data(),
                                       colIndices.data(),
                                       values.data());

  // Check that the isSymmetric method gives the correct answer.
  EXPECT_EQ(stag::isSymmetric(&matrix), false);
}

TEST(UtilityTest, SprsMatToVec) {
  // Construct a sparse matrix.
  std::vector<StagInt> colStarts = {0, 4};
  std::vector<StagInt> rowIndices = {0, 1, 2, 3};
  std::vector<double> values = {1, 2, 3, 4};
  SprsMat sparse_matrix = Eigen::Map<SprsMat>(4, 1, 4,
                                              colStarts.data(),
                                              rowIndices.data(),
                                              values.data());

  // Check that the dense vector gives the correct answer
  EXPECT_EQ(stag::sprsMatToVec(&sparse_matrix), values);

  // Construct a sparse vector with some zeros
  rowIndices = {0, 2, 4, 6};
  sparse_matrix = Eigen::Map<SprsMat>(7, 1, 4,
                                      colStarts.data(),
                                      rowIndices.data(),
                                      values.data());

  // The expected vector
  std::vector<double> expected_vec = {1, 0, 2, 0, 3, 0, 4, 0, 0};
  EXPECT_EQ(stag::sprsMatToVec(&sparse_matrix, 9), expected_vec);
}

TEST(UtilityTest, SprsMatToVecArguments) {
  // Construct a sparse matrix.
  std::vector<StagInt> colStarts = {0, 4};
  std::vector<StagInt> rowIndices = {0, 1, 2, 3};
  std::vector<double> values = {1, 2, 3, 4};
  SprsMat sparse_matrix = Eigen::Map<SprsMat>(4, 1, 4,
                                              colStarts.data(),
                                              rowIndices.data(),
                                              values.data());

  // Check the parameter checking of the sprsMatToVec method
  StagInt n = -1;
  EXPECT_THROW(stag::sprsMatToVec(&sparse_matrix, n), std::invalid_argument);

  n = 0;
  EXPECT_THROW(stag::sprsMatToVec(&sparse_matrix, n), std::invalid_argument);
}

TEST(UtilityTest, AddDoubleVectors) {
  std::vector<double> v1 = {1, 2, 3};
  std::vector<double> v2 = {1, 1};
  std::vector<double> expected_ans = {2, 3, 3};
  std::vector<double> ans = stag::addVectors(v1, v2);
  EXPECT_EQ(ans, expected_ans);
}

TEST(UtilityTest, AddIntVectors) {
  std::vector<StagInt> v1 = {1, 4};
  std::vector<StagInt> v2 = {9, 8, 9};
  std::vector<StagInt> expected_ans = {10, 12, 9};
  std::vector<StagInt> ans = stag::addVectors(v1, v2);
  EXPECT_EQ(ans, expected_ans);
}
