
#include "definitions.h"
#include <string>
#include "KDE/data.h"
#include "KDE/kde.h"
#include <iostream>
#include "Graph/graphio.h"

void test_mnist() {
  // Load the mnist dataset
  std::string filename = "test/data/mnist.txt";
  DenseMat data = stag::load_matrix(filename);

  // Create tha approximate similarity graph from this matrix.
  StagReal a = 0.000001;
  stag::Graph asg = stag::approximate_similarity_graph(&data, a);
  std::cout << asg.number_of_vertices() << std::endl;
  std::cout << asg.number_of_edges() << std::endl;

  std::string fname = "mnist_similarity.el";
  stag::save_edgelist(asg, fname);
}

void test_moons() {
  // Load the two moons dataset
  std::string filename = "test/data/moons_short.txt";
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
  test_mnist();
  return 0;
}
