#include "graph.h"

//------------------------------------------------------------------------------
// Graph Object Constructors
//------------------------------------------------------------------------------

stag::Graph::Graph(const SprsMat& adjacency_matrix) {
  adjacency_matrix_ = adjacency_matrix;
  adjacency_matrix_.makeCompressed();
  adj_init_ = true;
  lap_init_ = false;
  deg_init_ = false;
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
  deg_init_ = false;
}

//------------------------------------------------------------------------------
// Graph Object Public Methods
//------------------------------------------------------------------------------

const SprsMat* stag::Graph::adjacency() {
  return &adjacency_matrix_;
}

const SprsMat* stag::Graph::laplacian() {
  initialise_laplacian_();
  return &laplacian_matrix_;
}

const SprsMat* stag::Graph::degree_matrix() {
  initialise_degree_matrix_();
  return &degree_matrix_;
}

double stag::Graph::total_volume() {
  Eigen::VectorXd degrees = adjacency_matrix_ * Eigen::VectorXd::Ones(adjacency_matrix_.cols());
  return degrees.sum();
}

//------------------------------------------------------------------------------
// Graph Object Private Methods
//------------------------------------------------------------------------------

void stag::Graph::initialise_laplacian_() {
  // If the laplacian matrix has already been initialised, then we do not
  // initialise it again.
  if (lap_init_) return;

  // Ensure that the degree matrix is initialised
  initialise_degree_matrix_();

  // Construct and return the laplacian matrix.
  laplacian_matrix_ = degree_matrix_ - adjacency_matrix_;
  laplacian_matrix_.makeCompressed();

  // We have now initialised the laplacian.
  lap_init_ = true;
}

void stag::Graph::initialise_degree_matrix_() {
  // If the degree matrix has already been initialised, then we do not
  // initialise it again.
  if (deg_init_) return;

  // Construct the vertex degrees.
  Eigen::VectorXd degrees = adjacency_matrix_ * Eigen::VectorXd::Ones(adjacency_matrix_.cols());
  degree_matrix_ = SprsMat(adjacency_matrix_.cols(), adjacency_matrix_.cols());
  for (int i = 0; i < adjacency_matrix_.cols(); i++) {
    degree_matrix_.insert(i, i) = degrees[i];
  }

  // Compress the degree matrix storage, and set the initialised flag
  degree_matrix_.makeCompressed();
  deg_init_ = true;
}

//------------------------------------------------------------------------------
// Standard Graph Constructors
//------------------------------------------------------------------------------

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
