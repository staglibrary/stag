/**
 * Tests for the methods in the lsh.h header file. Includes methods for
 * computing the Euclidean locality-sensitive hashing.
 *
 * This file is provided as part of the STAG library and released under the GPL
 * license.
 */
#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>
#include <lsh.h>

//------------------------------------------------------------------------------
// Tests for the DataPoint class.
//------------------------------------------------------------------------------
TEST(LSHTest, DataPointDirectInit) {
  // Create an underlying array
  StagUInt dim = 10;
  StagReal data[dim];

  // Initialise the data point
  stag::DataPoint data_point(dim, data);

  // Check that the coordinates are correct
  EXPECT_EQ(dim, data_point.dimension);
  for (StagUInt i = 0; i < dim; i++) {
    EXPECT_EQ(data[i], data_point.coordinates[i]);
  }
}

TEST(LSHTest, DataPointMatInit) {
  // Create a dense matrix
  DenseMat data_mat {{0, 1, 2}, {10, 11, 12}, {20, 21, 22}, {30, 31, 32}};

  // Initialise the data point
  StagInt row_index = 1;
  stag::DataPoint data_point(data_mat, row_index);

  // Check that the coordinates are correct
  EXPECT_EQ(3, data_point.dimension);
  for (StagUInt i = 0; i < 3; i++) {
    EXPECT_EQ(data_mat.coeff(row_index, i), data_point.coordinates[i]);
  }
}

TEST(LSHTest, DataPointVectorInit) {
  // Create a data vector
  StagUInt dim = 10;
  std::vector<StagReal> data;
  data.resize(dim);

  // Initialise the data point
  stag::DataPoint data_point(data);

  // Check that the coordinates are correct
  EXPECT_EQ(dim, data_point.dimension);
  for (StagUInt i = 0; i < dim; i++) {
    EXPECT_EQ(data[i], data_point.coordinates[i]);
  }
}

