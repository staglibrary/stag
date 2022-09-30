#include "stag.h"

stag::Graph::Graph(const SprsMat& adjacency_matrix) {
  adjacency_matrix_ = adjacency_matrix;
}

stag::Graph::Graph(int vertices, std::vector<int> rows, std::vector<int> cols,
                   std::vector<double> vals) {
  adjacency_matrix_ = SprsMat(vertices, vertices);

  // For now, do the lazy thing of iterating over the inputs and adding them to
  // the sparse matrix.
  for (int i = 0; i < rows.size(); i++) {
    adjacency_matrix_.insert(rows[i], cols[i]) = vals[i];
  }
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
