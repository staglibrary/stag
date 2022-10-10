#include <stdexcept>
#include "graph.h"
#include "utility.h"

//------------------------------------------------------------------------------
// Graph Object Constructors
//------------------------------------------------------------------------------

stag::Graph::Graph(const SprsMat& adjacency_matrix) {
  // Load the adjacency matrix into this object.
  adjacency_matrix_ = adjacency_matrix;
  adjacency_matrix_.makeCompressed();

  // The number of vertices is the dimensions of the adjacency matrix
  number_of_vertices_ = adjacency_matrix_.outerSize();

  // Set the flags to indicate which matrices have been initialised.
  adj_init_ = true;
  lap_init_ = false;
  deg_init_ = false;
  norm_lap_init_ = false;

  // Check that the graph is configured correctly
  self_test_();
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

  // The number of vertices is the dimensions of the adjacency matrix
  number_of_vertices_ = adjacency_matrix_.outerSize();

  // Set the flags to indicate which matrices have been initialised.
  adj_init_ = true;
  lap_init_ = false;
  deg_init_ = false;
  norm_lap_init_ = false;

  // Check that the graph is configured correctly
  self_test_();
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

const SprsMat* stag::Graph::normalised_laplacian() {
  initialise_normalised_laplacian_();
  return &normalised_laplacian_matrix_;
}

const SprsMat* stag::Graph::degree_matrix() {
  initialise_degree_matrix_();
  return &degree_matrix_;
}

double stag::Graph::total_volume() {
  Eigen::VectorXd degrees = adjacency_matrix_ * Eigen::VectorXd::Ones(adjacency_matrix_.cols());
  return degrees.sum();
}

long stag::Graph::number_of_vertices() const {
  return number_of_vertices_;
}

long stag::Graph::number_of_edges() const {
  return adjacency_matrix_.nonZeros() / 2;
}

//------------------------------------------------------------------------------
// Graph Object Private Methods
//------------------------------------------------------------------------------

void stag::Graph::self_test_() {
  // Check that the adjacency matrix is symmetric.
  if (!stag::isSymmetric(&adjacency_matrix_)) {
    throw std::domain_error("Graph adjacency matrix must be symmetric.");
  }
}

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

void stag::Graph::initialise_normalised_laplacian_() {
  // If the normalised laplacian matrix has already been initialised, then we
  // do not initialise it again.
  if (norm_lap_init_) return;

  // Ensure that the degree matrix is initialised
  initialise_degree_matrix_();

  // Construct the inverse degree matrix
  SprsMat sqrt_inv_deg_mat(number_of_vertices_, number_of_vertices_);
  std::vector<Eigen::Triplet<double>> non_zero_entries;
  for (int i = 0; i < number_of_vertices_; i++) {
    non_zero_entries.emplace_back(i, i, 1 / sqrt(degree_matrix_.coeff(i, i)));
  }
  sqrt_inv_deg_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());

  // The normalised laplacian is defined by I - D^{-1/2} A D^{-1/2}
  SprsMat identity_matrix(number_of_vertices_, number_of_vertices_);
  identity_matrix.setIdentity();
  normalised_laplacian_matrix_ = identity_matrix - sqrt_inv_deg_mat * adjacency_matrix_ * sqrt_inv_deg_mat;
  normalised_laplacian_matrix_.makeCompressed();

  // We have now initialised the normalised laplacian matrix.
  norm_lap_init_ = true;
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
