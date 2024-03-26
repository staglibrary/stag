
#include "definitions.h"
#include <string>
#include "KDE/data.h"
#include "KDE/kde.h"
#include <iostream>

void test_mnist() {
  // Load the mnist dataset
  std::string filename = "test/data/mnist.txt";
  DenseMat data = stag::load_matrix(filename);

  // Create a CKNS KDE estimator
  StagReal a = 0.000001;
  StagReal eps = 0.99;
  stag::CKNSGaussianKDE ckns_kde(&data, a, eps);
  std::vector<StagReal> kde_estimates = ckns_kde.query(&data);

  for (auto estimate : kde_estimates) {
    std::cout << estimate << std::endl;
  }
}

void test_moons() {
  // Load the two moons dataset
  std::string filename = "test/data/moons.txt";
  DenseMat data = stag::load_matrix(filename);

  // Create a CKNS KDE estimator
  StagReal a = 20;
  StagReal eps = 0.5;
  stag::CKNSGaussianKDE ckns_kde(&data, a, eps);
  std::vector<StagReal> kde_estimates = ckns_kde.query(&data);

  for (auto estimate : kde_estimates) {
    std::cout << estimate << std::endl;
  }
}


int main(int argc, char** args) {
  test_moons();
  return 0;
}
