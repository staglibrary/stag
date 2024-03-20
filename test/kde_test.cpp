/**
 * Tests for the methods in the kde.h header file. Includes methods for
 * computing kernel density estimates.
 *
 * This file is provided as part of the STAG library and released under the GPL
 * license.
 */
#include <stdexcept>
#include <gtest/gtest.h>
#include <kde.h>

//------------------------------------------------------------------------------
// Tests for the Gaussian kernel.
//------------------------------------------------------------------------------
TEST(KDETest, GaussianKernelDist) {
  // Create some test distances
  StagReal a = 1.2;
  std::vector<StagReal> distances = {0, 0.1, 0.5, 1, 1.5, 2, 10};

  // Define the expected Gaussian kernel values
  std::vector<StagReal> expected_kernel_values = {1,
                                                  0.8869204,
                                                  0.548812,
                                                  0.301194,
                                                  0.165299,
                                                  0.0907180,
                                                  0.00000614421};

  assert(distances.size() == expected_kernel_values.size());

  for (auto i = 0; i < distances.size(); i++) {
    StagReal value = stag::gaussian_kernel(a, distances[i]);
    EXPECT_LE(value, 1.01 * expected_kernel_values[i]);
    EXPECT_GE(value, 0.99 * expected_kernel_values[i]);
  }
}

TEST(KDETest, GaussianKernelPoint) {
  // Create some test points
  DenseMat data {{0, 0, 0}, {0, 1, 1}};
  stag::DataPoint x1(data, 0);
  stag::DataPoint x2(data, 1);

  // Define the expected kernel distance
  StagReal a = 1.5;
  StagReal expected_kernel_value = 0.0497871;

  // Check that the calculated kernel value is correct
  StagReal value = stag::gaussian_kernel(a, x1, x2);
  EXPECT_LE(value, 1.01 * expected_kernel_value);
  EXPECT_GE(value, 0.99 * expected_kernel_value);
}
