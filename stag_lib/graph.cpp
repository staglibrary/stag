/**
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <stdexcept>
#include "graph.h"
#include "utility.h"


//------------------------------------------------------------------------------
// Local Graph Methods
//------------------------------------------------------------------------------

std::vector<double> stag::LocalGraph::degrees(std::vector<stag_int> vertices) {
  std::vector<double> degrees;

  for (stag_int v : vertices) {
    degrees.emplace_back(degree(v));
  }

  return degrees;
}

std::vector<stag_int> stag::LocalGraph::degrees_unweighted(
    std::vector<stag_int> vertices) {
  std::vector<stag_int> degrees;

  for (stag_int v : vertices) {
    degrees.emplace_back(degree_unweighted(v));
  }

  return degrees;
}

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
  lap_init_ = false;
  deg_init_ = false;
  inv_deg_init_ = false;
  norm_lap_init_ = false;
  lazy_rand_walk_init_ = false;

  // Check that the graph is configured correctly
  self_test_();
}

stag::Graph::Graph(std::vector<stag_int> &outerStarts, std::vector<stag_int> &innerIndices,
                   std::vector<double> &values) {
  // Map the provided data vectors to the sparse matrix type.
  adjacency_matrix_ = stag::sprsMatFromVectors(outerStarts, innerIndices, values);

  // The number of vertices is the dimensions of the adjacency matrix
  number_of_vertices_ = adjacency_matrix_.outerSize();

  // Set the flags to indicate which matrices have been initialised.
  lap_init_ = false;
  deg_init_ = false;
  inv_deg_init_ = false;
  norm_lap_init_ = false;
  lazy_rand_walk_init_ = false;

  // Check that the graph is configured correctly
  self_test_();
}

//------------------------------------------------------------------------------
// Graph Object Public Methods
//------------------------------------------------------------------------------

const SprsMat* stag::Graph::adjacency() const {
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

const SprsMat* stag::Graph::inverse_degree_matrix() {
  initialise_inverse_degree_matrix_();
  return &inverse_degree_matrix_;
}

const SprsMat* stag::Graph::lazy_random_walk_matrix() {
  initialise_lazy_random_walk_matrix_();
  return &lazy_random_walk_matrix_;
}

double stag::Graph::total_volume() {
  Eigen::VectorXd degrees = adjacency_matrix_ * Eigen::VectorXd::Ones(adjacency_matrix_.cols());
  return degrees.sum();
}

stag_int stag::Graph::number_of_vertices() const {
  return number_of_vertices_;
}

stag_int stag::Graph::number_of_edges() const {
  return adjacency_matrix_.nonZeros() / 2;
}

double stag::Graph::degree(stag_int v) {
  // For now, we can be a little lazy and use the degree matrix. Once this is
  // initialised, then checking degree is constant time.
  initialise_degree_matrix_();

  // If the vertex number is larger than the number of vertices in the graph,
  // then we say that the degree is 0.
  if (v < number_of_vertices_) {
    return degree_matrix_.coeff(v, v);
  } else {
    return 0;
  }
}

stag_int stag::Graph::degree_unweighted(stag_int v) {
  // If the requested vertex number is greater than the number of vertices, then
  // return a degree of 0.
  if (v >= number_of_vertices_) {
    return 0;
  }

  // The combinatorical degree of a vertex is equal to the number of non-zero
  // entries in its adjacency matrix row.
  const stag_int *indexPtr = adjacency_matrix_.outerIndexPtr();
  stag_int rowStart = *(indexPtr + v);
  stag_int nextRowStart = *(indexPtr + v + 1);
  return nextRowStart - rowStart;
}

std::vector<stag::edge> stag::Graph::neighbors(stag_int v) {
  // If the vertex v does not exist, return the empty vector.
  if (v >= number_of_vertices_) {
    return {};
  }

  std::vector<stag::edge> edges;

  // Iterate through the non-zero entries in the vth row of the adjacency matrix
  const double *weights = adjacency_matrix_.valuePtr();
  const stag_int *innerIndices = adjacency_matrix_.innerIndexPtr();
  const stag_int *rowStarts = adjacency_matrix_.outerIndexPtr();
  stag_int vRowStart = *(rowStarts + v);
  stag_int degree_unw = degree_unweighted(v);

  for (stag_int i = 0; i < degree_unw; i++) {
    if (*(weights + vRowStart + i) != 0) {
      edges.push_back({v, *(innerIndices + vRowStart + i), *(weights + vRowStart + i)});
    }
  }

  return edges;
}

std::vector<stag_int> stag::Graph::neighbors_unweighted(stag_int v) {
  // If the vertex v does not exist, return the empty vector.
  if (v >= number_of_vertices_) {
    return {};
  }

  // Return the non-zero indices in the vth row of the adjacency matrix
  const stag_int *innerIndices = adjacency_matrix_.innerIndexPtr();
  const stag_int *rowStarts = adjacency_matrix_.outerIndexPtr();
  stag_int vRowStart = *(rowStarts + v);
  stag_int degree = degree_unweighted(v);
  return {innerIndices + vRowStart, innerIndices + vRowStart + degree};
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
  std::vector<EdgeTriplet> non_zero_entries;
  for (stag_int i = 0; i < number_of_vertices_; i++) {
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
  for (stag_int i = 0; i < adjacency_matrix_.cols(); i++) {
    degree_matrix_.insert(i, i) = degrees[i];
  }

  // Compress the degree matrix storage, and set the initialised flag
  degree_matrix_.makeCompressed();
  deg_init_ = true;
}

void stag::Graph::initialise_inverse_degree_matrix_() {
  // If the inverse degree matrix has already been initialised, then we do not
  // initialise it again.
  if (inv_deg_init_) return;

  // We will construct the inverse degree matrix from the degree matrix itself
  initialise_degree_matrix_();
  inverse_degree_matrix_ = SprsMat(adjacency_matrix_.cols(), adjacency_matrix_.cols());
  for (stag_int i = 0; i < adjacency_matrix_.cols(); i++) {
    inverse_degree_matrix_.insert(i, i) = 1./degree_matrix_.coeff(i, i);
  }

  // Compress the degree matrix storage, and set the initialised flag
  inverse_degree_matrix_.makeCompressed();
  inv_deg_init_ = true;
}

void stag::Graph::initialise_lazy_random_walk_matrix_() {
  // If the lazy random walk matrix has already been initialised, then we do not
  // initialise it again.
  if (lazy_rand_walk_init_) return;

  // The lazy random walk matrix is defined to be
  //   (1/2) I + (1/2) A * D^{-1}
  initialise_inverse_degree_matrix_();
  SprsMat identityMatrix(number_of_vertices_, number_of_vertices_);
  identityMatrix.setIdentity();

  lazy_random_walk_matrix_ = SprsMat(number_of_vertices_, number_of_vertices_);
  lazy_random_walk_matrix_ = (1./2) * identityMatrix + (1./2) * adjacency_matrix_ * inverse_degree_matrix_;

  // Compress and set initialisation flag
  lazy_random_walk_matrix_.makeCompressed();
  lazy_rand_walk_init_ = true;
}

//------------------------------------------------------------------------------
// Equality Operators
//------------------------------------------------------------------------------
bool stag::operator==(const stag::Graph& lhs, const stag::Graph& rhs) {
  bool outerIndicesEqual = stag::sprsMatOuterStarts(lhs.adjacency()) == stag::sprsMatOuterStarts(rhs.adjacency());
  bool innerIndicesEqual = stag::sprsMatInnerIndices(lhs.adjacency()) == stag::sprsMatInnerIndices(rhs.adjacency());
  bool valuesEqual = stag::sprsMatValues(lhs.adjacency()) == stag::sprsMatValues(rhs.adjacency());
  return (outerIndicesEqual && innerIndicesEqual) && valuesEqual;
}

bool stag::operator!=(const stag::Graph &lhs, const stag::Graph &rhs) {
  return !(lhs == rhs);
}

bool stag::operator==(const stag::edge &lhs, const stag::edge &rhs) {
  return lhs.v1 == rhs.v1 && lhs.v2 == rhs.v2 && lhs.weight == rhs.weight;
}

bool stag::operator!=(const stag::edge &lhs, const stag::edge &rhs) {
  return !(lhs == rhs);
}

//------------------------------------------------------------------------------
// Standard Graph Constructors
//------------------------------------------------------------------------------
stag::Graph stag::cycle_graph(stag_int n) {
  SprsMat adj_mat(n, n);
  std::vector<EdgeTriplet> non_zero_entries;
  for (stag_int i = 0; i < n; i++) {
    non_zero_entries.emplace_back(i, (i + n + 1) % n, 1);
    non_zero_entries.emplace_back(i, (i + n - 1) % n, 1);
  }
  adj_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());
  return stag::Graph(adj_mat);
}

stag::Graph stag::complete_graph(stag_int n) {
  SprsMat adj_mat(n, n);
  std::vector<EdgeTriplet> non_zero_entries;
  for (stag_int i = 0; i < n; i++) {
    for (stag_int j = 0; j < n; j++) {
      if (i != j) {
        non_zero_entries.emplace_back(i, j, 1);
      }
    }
  }
  adj_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());
  return stag::Graph(adj_mat);
}

stag::Graph stag::barbell_graph(stag_int n) {
  // Construct the non-zero entries in the complete blocks of the adjacency
  // matrix
  std::vector<EdgeTriplet> non_zero_entries;
  for (stag_int i = 0; i < n; i++) {
    for (stag_int j = 0; j < n; j++) {
      if (i != j) {
        non_zero_entries.emplace_back(i, j, 1);
        non_zero_entries.emplace_back(n + i, n + j, 1);
      }
    }
  }

  // Add a single edge to connect the complete graphs
  non_zero_entries.emplace_back(n - 1, n, 1);
  non_zero_entries.emplace_back(n, n - 1, 1);

  // Construct the final adjacency matrix
  SprsMat adj_mat(2 * n, 2 * n);
  adj_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());
  return stag::Graph(adj_mat);
}

stag::Graph stag::star_graph(stag_int n) {
  // Construct the non-zero entries in the adjacency matrix
  std::vector<EdgeTriplet> non_zero_entries;
  for (stag_int i = 1; i < n; i++) {
    non_zero_entries.emplace_back(0, i, 1);
    non_zero_entries.emplace_back(i, 0, 1);
  }

  // Construct the final adjacency matrix
  SprsMat adj_mat(n, n);
  adj_mat.setFromTriplets(non_zero_entries.begin(), non_zero_entries.end());
  return stag::Graph(adj_mat);
}
