#include "graph.h"

stag::Graph::Graph(const SprsMat& adjacency_matrix) {
  adjacency_matrix_ = adjacency_matrix;
  adjacency_matrix_.makeCompressed();
  adj_init_ = true;
  lap_init_ = false;
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
  adjacency_matrix_.makeCompressed();
  adj_init_ = true;
  lap_init_ = false;
}

const SprsMat* stag::Graph::adjacency() {
  return &adjacency_matrix_;
}

void stag::Graph::initialise_laplacian_() {
  // If the laplacian matrix has already been initialised, then we do not
  // initialise it again.
  if (lap_init_) return;

  // Construct the vertex degrees.
  Eigen::VectorXd degrees = adjacency_matrix_ * Eigen::VectorXd::Ones(adjacency_matrix_.cols());
  SprsMat degree_matrix(adjacency_matrix_.cols(), adjacency_matrix_.cols());
  for (int i = 0; i < adjacency_matrix_.cols(); i++) {
    degree_matrix.insert(i, i) = degrees[i];
  }

  // Construct and return the laplacian matrix.
  laplacian_matrix_ = degree_matrix - adjacency_matrix_;
  laplacian_matrix_.makeCompressed();

  // We have now initialised the laplacian.
  lap_init_ = true;
}

const SprsMat* stag::Graph::laplacian() {
  initialise_laplacian_();
  return &laplacian_matrix_;
}

double stag::Graph::volume() {
  Eigen::VectorXd degrees = adjacency_matrix_ * Eigen::VectorXd::Ones(adjacency_matrix_.cols());
  return degrees.sum();
}

stag::Graph stag::cycle_graph(int n) {
  SprsMat adj_mat(n, n);
  std::vector<Eigen::Triplet<double>> non_zero_entries;
  for (int i = 0; i < n; i++) {
    non_zero_entries.emplace_back(i, (i + n + 1) % n, 1);
    non_zero_entries.emplace_back(i, (i + n - 1) % n, 1);
  }
  adj_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());
  return stag::Graph(adj_mat);
}

stag::Graph stag::complete_graph(int n) {
  SprsMat adj_mat(n, n);
  std::vector<Eigen::Triplet<double>> non_zero_entries;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) {
        non_zero_entries.emplace_back(i, j, 1);
      }
    }
  }
  adj_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());
  return stag::Graph(adj_mat);
}
