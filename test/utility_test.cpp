/**
 * Tests for the utility.h header. Includes helper methods for dealing with
 * sparse matrices.
 *
 * Copyright 2022 Peter Macgregor
 */
#include <gtest/gtest.h>
#include <utility.h>

TEST(UtilityTest, IsSymmetric) {
  // Construct a symmetric matrix.
  std::vector<int> rowStarts = {0, 2, 4, 7, 8};
  std::vector<int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
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
