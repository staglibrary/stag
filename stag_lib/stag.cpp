#include "stag.h"

stag::Graph::Graph(const SprsMat& adjacency_matrix) {
  adjacency_matrix_ = adjacency_matrix;
}

stag::Graph::Graph(std::vector<int> &outerStarts, std::vector<int> &innerIndices,
                   std::vector<double> &values) {
  // Map the provided data vectors to the sparse matrix type.
  adjacency_matrix_ = Eigen::Map<SprsMat>(long(outerStarts.size()) - 1,
                                          long(outerStarts.size()) - 1,
                                          long(values.size()),
                                          outerStarts.data(),
                                          innerIndices.data(),
                                          values.data());
}

SprsMat stag::Graph::adjacency() {
  return adjacency_matrix_;
}

SprsMat stag::Graph::laplacian() {
  // Construct the vertex degrees.
  Eigen::VectorXd degrees = adjacency_matrix_ * Eigen::VectorXd::Ones(adjacency_matrix_.cols());
  SprsMat degree_matrix(adjacency_matrix_.cols(), adjacency_matrix_.cols());
  for (int i = 0; i < adjacency_matrix_.cols(); i++) {
    degree_matrix.insert(i, i) = degrees[i];
  }

  // Construct and return the laplacian matrix.
  return degree_matrix - adjacency_matrix_;
}

double stag::Graph::volume() {
  Eigen::VectorXd degrees = adjacency_matrix_ * Eigen::VectorXd::Ones(adjacency_matrix_.cols());
  return degrees.sum();
}
