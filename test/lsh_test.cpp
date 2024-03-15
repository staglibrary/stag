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

TEST(LSHTest, LSHFunction) {
  // Create two data vectors
  StagInt dim = 3;
  DenseMat data_mat {{0, 1, 1}, {1, 0, 0}};
  stag::DataPoint x1(data_mat, 0);
  stag::DataPoint x2(data_mat, 1);

  // The distance between these vectors is root 3.
  StagReal distance = sqrt(dim);

  // Get the collision probability under the LSH function.
  StagReal prob = stag::LSHFunction::collision_probability(distance);

  // Construct 1000 LSH functions
  StagInt num_functions = 1000;
  std::vector<stag::LSHFunction> functions;
  for (StagInt i = 0; i < num_functions; i++) {
    functions.emplace_back(dim);
  }

  // Compute the number for which the given vectors collide
  StagUInt num_collisions = 0;
  for (auto fun : functions) {
    if (fun.apply(x1) == fun.apply(x2)) {
      num_collisions++;
    }
  }

  // Check that the proportion of collisions is close to the expected number
  EXPECT_LE(num_collisions, 1.2 * prob * num_functions);
  EXPECT_GE(num_collisions, 0.8 * prob * num_functions);
}

