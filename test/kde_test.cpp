/**
 * Tests for the methods in the kde.h header file. Includes methods for
 * computing kernel density estimates.
 *
 * This file is provided as part of the STAG library and released under the GPL
 * license.
 */
#include <stdexcept>
#include <gtest/gtest.h>
#include "KDE/kde.h"

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

//------------------------------------------------------------------------------
// Tests for the CKNS Gaussian KDE.
//------------------------------------------------------------------------------

TEST(KDETest, CKNSTwoMoons) {
  // Load the two moons dataset
  std::string filename = "test/data/moons.txt";
  DenseMat data = stag::load_matrix(filename);

  // Create a CKNS KDE estimator
  StagReal a = 20;
  StagReal eps = 0.5;
  stag::CKNSGaussianKDE ckns_kde(&data, a, eps);

  // Create an exact Gaussian KDE
  stag::ExactGaussianKDE exact_kde(&data, a);

  // Check that the estimates are accurate
  std::vector<StagReal> kde_exact = exact_kde.query(&data);
  std::vector<StagReal> kde_estimates = ckns_kde.query(&data);

  StagReal total_error = 0;
  for (auto i = 0; i < kde_estimates.size(); i++) {
    total_error += abs(kde_estimates.at(i) - kde_exact.at(i)) / kde_exact.at(i);
  }
  StagReal avg_error = total_error / (StagReal) kde_estimates.size();
  EXPECT_LE(avg_error, 0.5 * eps);
}

TEST(KDETest, CKNSDefaultConstructor) {
  // Load the two moons dataset
  std::string filename = "test/data/moons.txt";
  DenseMat data = stag::load_matrix(filename);

  // Create a CKNS KDE estimator
  StagReal a = 20;
  stag::CKNSGaussianKDE ckns_kde(&data, a);

  // Create an exact Gaussian KDE
  stag::ExactGaussianKDE exact_kde(&data, a);

  // Check that the estimates are accurate
  std::vector<StagReal> kde_exact = exact_kde.query(&data);
  std::vector<StagReal> kde_estimates = ckns_kde.query(&data);

  StagReal total_error = 0;
  for (auto i = 0; i < kde_estimates.size(); i++) {
    total_error += abs(kde_estimates.at(i) - kde_exact.at(i)) / kde_exact.at(i);
  }
  StagReal avg_error = total_error / (StagReal) kde_estimates.size();
  EXPECT_LE(avg_error, 0.25);
}

TEST(KDETest, CKNSMinMu) {
  // Load the two moons dataset
  std::string filename = "test/data/moons.txt";
  DenseMat data = stag::load_matrix(filename);

  // Create a CKNS KDE estimator
  StagReal a = 10;
  StagReal eps = 0.5;
  StagReal min_mu = 0.005;
  stag::CKNSGaussianKDE ckns_kde(&data, a, eps, min_mu);

  // Create an exact Gaussian KDE
  stag::ExactGaussianKDE exact_kde(&data, a);

  // Check that the estimates are accurate
  std::vector<StagReal> kde_exact = exact_kde.query(&data);
  std::vector<StagReal> kde_estimates = ckns_kde.query(&data);

  StagReal total_error = 0;
  for (auto i = 0; i < kde_estimates.size(); i++) {
    total_error += abs(kde_estimates.at(i) - kde_exact.at(i)) / kde_exact.at(i);
  }
  StagReal avg_error = total_error / (StagReal) kde_estimates.size();
  EXPECT_LE(avg_error, 0.5 * eps);
}

TEST(KDETest,CKNSExplicitConstants){
  // Load the two moons dataset
  std::string filename = "test/data/moons.txt";
  DenseMat data = stag::load_matrix(filename);

  // Create a CKNS KDE estimator
  StagReal a = 10;
  StagInt K1 = 40;
  StagReal K2_constant = 50;
  StagReal min_mu = 0.005;
  StagInt offset = 1;
  stag::CKNSGaussianKDE ckns_kde(&data, a, min_mu, K1, K2_constant, offset);

  // Create an exact Gaussian KDE
  stag::ExactGaussianKDE exact_kde(&data, a);

  // Check that the estimates are accurate
  std::vector<StagReal> kde_exact = exact_kde.query(&data);
  std::vector<StagReal> kde_estimates = ckns_kde.query(&data);

  StagReal total_error = 0;
  for (auto i = 0; i < kde_estimates.size(); i++) {
    total_error += abs(kde_estimates.at(i) - kde_exact.at(i)) / kde_exact.at(i);
  }
  StagReal avg_error = total_error / (StagReal) kde_estimates.size();
  EXPECT_LE(avg_error, 0.25);
}

TEST(KDETest, CKNSMnist) {
  // Load the two MNIST dataset
  std::string filename = "test/data/mnist.txt";
  DenseMat data = stag::load_matrix(filename);
  StagReal a = 0.000001;

  // Create an exact Gaussian KDE
  stag::ExactGaussianKDE exact_kde(&data, a);
  std::vector<StagReal> kde_exact = exact_kde.query(&data);

  // Create a CKNS KDE estimator
  stag::CKNSGaussianKDE ckns_kde(&data, a);
  std::vector<StagReal> kde_estimates = ckns_kde.query(&data);

  // Check that the estimates are accurate
  StagReal total_error = 0;
  for (auto i = 0; i < kde_estimates.size(); i++) {
    total_error += abs(kde_estimates.at(i) - kde_exact.at(i)) / kde_exact.at(i);
  }
  StagReal avg_error = total_error / (StagReal) kde_estimates.size();
  EXPECT_LE(avg_error, 0.5);
}
