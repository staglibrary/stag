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
#include "lsh.h"

//------------------------------------------------------------------------------
// Tests for the DataPoint class.
//------------------------------------------------------------------------------
TEST(LSHTest, DataPointDirectInit) {
  // Create an underlying array
  StagUInt dim = 10;
  StagReal data[10];

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

TEST(LSHTest, CloseCollisionProb) {
  StagReal distance = 0.05;

  // Get the collision probability under the LSH function.
  StagReal prob = stag::LSHFunction::collision_probability(distance);

  // For this small value, the collision probability should be close to
  // 1 - 0.2 x.
  StagReal expected_prob = 1 - 0.2 * distance;
  EXPECT_NEAR(prob, expected_prob, 0.0001);
}

TEST(LSHTest, NegDistCollisionProb) {
  StagReal distance = 5.6;
  StagReal prob1 = stag::LSHFunction::collision_probability(distance);
  StagReal prob2 = stag::LSHFunction::collision_probability(-distance);
  EXPECT_EQ(prob1, prob2);

  distance = 0.2;
  prob1 = stag::LSHFunction::collision_probability(distance);
  prob2 = stag::LSHFunction::collision_probability(-distance);
  EXPECT_EQ(prob1, prob2);
}

TEST(LSHTest, E2LSH11) {
  // Check that an E2LSH table with K=1 and L=1 behaves like a single
  // LSHFunction.

  // Create two data vectors
  StagInt dim = 3;
  DenseMat data_mat {{0, 1, 1}, {1, 0, 0}};
  stag::DataPoint x1(data_mat, 0);
  stag::DataPoint x2(data_mat, 1);

  // The 'data' vector for the E2LSH function will be just one of the vectors.
  std::vector<stag::DataPoint> dataset = {x1};

  // Construct 1000 E2LSH tables
  StagInt num_tables = 1000;
  std::vector<stag::E2LSH> tables;
  for (StagInt i = 0; i < num_tables; i++) {
    tables.emplace_back(1, 1, dataset);
  }

  // Compute the number for which the given vectors collide
  StagUInt num_collisions = 0;
  for (auto tab : tables) {
    std::vector<stag::DataPoint> results = tab.get_near_neighbors(x2);
    if (!results.empty()) {
      num_collisions++;
    }
  }

  // Get the collision probability under the LSH function.
  StagReal distance = sqrt((StagReal) dim);
  StagReal prob = stag::LSHFunction::collision_probability(distance);

  EXPECT_EQ(prob, stag::E2LSH::collision_probability(1, 1, distance));

  // Check that the proportion of collisions is close to the expected number
  EXPECT_LE(num_collisions, 1.2 * prob * num_tables);
  EXPECT_GE(num_collisions, 0.8 * prob * num_tables);
}


TEST(LSHTest, E2LSH510) {
  // Check that an E2LSH table with K=5 and L=10 creates the correct number of
  // collisions.
  StagUInt K = 5;
  StagUInt L = 10;

  // Create two data vectors
  StagInt dim = 3;
  DenseMat data_mat {{0, 1, 1}, {1, 0, 0}};
  stag::DataPoint x1(data_mat, 0);
  stag::DataPoint x2(data_mat, 1);

  // The 'data' vector for the E2LSH function will be just one of the vectors.
  std::vector<stag::DataPoint> dataset = {x1};

  // Construct 1000 E2LSH tables
  StagInt num_tables = 1000;
  std::vector<stag::E2LSH> tables;
  for (StagInt i = 0; i < num_tables; i++) {
    tables.emplace_back(K, L, dataset);
  }

  // Compute the number for which the given vectors collide
  StagUInt num_collisions = 0;
  for (auto tab : tables) {
    std::vector<stag::DataPoint> results = tab.get_near_neighbors(x2);
    if (!results.empty()) {
      num_collisions++;
    }
  }

  // Get the collision probability under the LSH function.
  StagReal distance = sqrt((StagReal) dim);
  StagReal prob = tables[0].collision_probability(distance);

  // Check that the proportion of collisions is close to the expected number
  EXPECT_LE(num_collisions, 1.2 * prob * num_tables);
  EXPECT_GE(num_collisions, 0.8 * prob * num_tables);
}

TEST(LSHTest, E2LSHMoreData) {
  // Check that an E2LSH table with K=1 and L=1 creates the correct number of
  // collisions for different data points.
  StagUInt K = 1;
  StagUInt L = 1;

  // Create data vectors
  DenseMat data_mat {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}};
  stag::DataPoint x1(data_mat, 0);
  stag::DataPoint x2(data_mat, 1);
  stag::DataPoint x3(data_mat, 2);
  stag::DataPoint x4(data_mat, 3);

  // The 'data' vector for the E2LSH function.
  std::vector<stag::DataPoint> dataset = {x2, x3, x4};

  // Construct 1000 E2LSH tables
  StagInt num_tables = 1000;
  std::vector<stag::E2LSH> tables;
  for (StagInt i = 0; i < num_tables; i++) {
    tables.emplace_back(K, L, dataset);
  }

  // Compute the number for which the given vectors collide
  std::vector<StagUInt> num_collisions = {0, 0, 0};
  for (auto tab : tables) {
    std::vector<stag::DataPoint> results = tab.get_near_neighbors(x1);
    for (auto res : results) {
      if (res.coordinates[2] == 1) {
        num_collisions[2]++;
      } else if (res.coordinates[1] == 1) {
        num_collisions[1]++;
      } else if (res.coordinates[0] == 1) {
        num_collisions[0]++;
      } else {
        // Should never get here.
        assert(false);
      }
    }
  }

  for (auto i = 1; i <= 3; i++) {
    // Get the collision probability for each data vector
    StagReal distance = sqrt((StagReal) i);
    StagReal prob = stag::E2LSH::collision_probability(K, L, distance);

    // Check that the proportion of collisions is close to the expected number
    EXPECT_LE(num_collisions[i-1], 1.2 * prob * num_tables);
    EXPECT_GE(num_collisions[i-1], 0.8 * prob * num_tables);
  }
}


