#include "stag.h"
#include <iostream>
#include <utility>

stag::Graph::Graph(Eigen::MatrixXd adjacency_matrix) {
  adjacency_matrix_ = std::move(adjacency_matrix);
}

Eigen::MatrixXd stag::Graph::laplacian() {
  // Construct the degree matrix
  Eigen::MatrixXd degree_matrix = adjacency_matrix_.rowwise().sum().asDiagonal();
  return degree_matrix - adjacency_matrix_;
};
