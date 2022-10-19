/**
 * Definitions related to the core Graph object used to represent graphs within
 * the library.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#ifndef STAG_LIBRARY_H
#define STAG_LIBRARY_H

// The fundamental datatype used in this library is the sparse matrix. For
// convenience, we define the sparse matrix type here.
#define stag_int long long
#define SprsMat Eigen::SparseMatrix<double, Eigen::ColMajor, stag_int>
#define EdgeTriplet Eigen::Triplet<double, stag_int>

// Redefine the eigen index type to be the same as stag_int
#undef EIGEN_DEFAULT_DENSE_INDEX_TYPE
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE stag_int

#include <Eigen/Sparse>
#include <vector>

namespace stag {
  /**
   * A structure representing a weighted edge in a graph.
   */
  struct edge {
    // The first vertex in the edge.
    stag_int v1;

    // The second vertex in the edge
    stag_int v2;

    // The weight of the edge.
    double weight;
  };

  /**
   * LocalGraph is an abstract class which defines methods for exploring the
   * local neighborhood of a graph.
   */
  class LocalGraph {
    public:
      /**
       * Given a vertex v, return its weighted degree.
       */
      virtual double degree(stag_int v) = 0;

      /**
       * Given a vertex v, return its unweighted degree. That is, the number
       * of neighbors of v, ignoring the edge weights.
       */
       virtual stag_int degree_unweighted(stag_int v) = 0;

      /**
       * Given a vertex v, return a vector of edges representing the
       * neighborhood of v.
       *
       * The returned edges will all have the ordering (v, x) such that
       * edge.v = v.
       *
       * @param v an int representing some vertex in the graph
       * @return an edge vector containing the neighborhood information
       */
      virtual std::vector<edge> neighbors(stag_int v) = 0;

      /**
       * Given a vertex v, return a vector containing the neighbors of v.
       *
       * The weights of edges to the neighbors are not returned by this method.
       *
       * @param v an int representing some vertex in the graph
       * @return an int vector giving the neighbors of v
       */
      virtual std::vector<stag_int> neighbors_unweighted(stag_int v) = 0;

      /**
       * Destructor for the LocalGraph object.
       */
      virtual ~LocalGraph() = default;
  };

  /**
   * The core object used to represent a graph for use with the library. Graphs
   * are always constructed from sparse matrices, and this is the internal
   * representation used as well.
   */
  class Graph : public LocalGraph {
    public:
      /**
       * Create a graph from an Eigen matrix.
       *
       * @param adjacency_matrix the sparse eigen matrix representing the adjacency matrix
       *               of the graph.
       */
      explicit Graph(const SprsMat& adjacency_matrix);

      /**
       * Create a graph from raw arrays describing a CSR sparse matrix.
       *
       * To use this constructor, you should understand the CSR sparse matrix
       * format. Note that this library uses the RowMajor format from the Eigen
       * library.
       *
       * @param outerStarts the indices of the start of each row in the CSR
       *                    matrix
       * @param innerIndices the column indices of each non-zero element in the
       *                     matrix
       * @param values the values of each non-zero element in the matrix
       */
      Graph(std::vector<stag_int> &outerStarts, std::vector<stag_int> &innerIndices,
            std::vector<double> &values);

      /**
       * Return the sparse adjacency matrix of the graph
       *
       * @return a sparse Eigen matrix representing the graph adjacency matrix.
       */
      const SprsMat* adjacency() const;

      /**
       * Construct the Laplacian matrix of the graph.
       *
       * The Laplacian matrix is defined by
       *   L = D - A
       * where D is the diagonal matrix of vertex degrees and A is the adjacency
       * matrix of the graph.
       *
       * @return a sparse Eigen matrix representing the graph Laplacian
       */
      const SprsMat* laplacian();

      /**
       * Construct the normalised Laplacian matrix of the graph.
       *
       * The normalised Laplacian matrix is defined by
       *   Ln = D^{-1/2} L D^{-1/2}
       * where D is the diagonal matrix of vertex degrees and L is the Laplacian
       * matrix of the graph.
       *
       * @return a sparse Eigen matrix representing the normalised Laplacian
       */
       const SprsMat* normalised_laplacian();

      /**
       * The degree matrix of the graph.
       *
       * The degree matrix is an n x n matrix such that each diagonal entry is
       * the degree of the corresponding node.
       *
       * @return a sparse Eigen matrix
       */
      const SprsMat* degree_matrix();

      /**
       * The inverse degree matrix of the graph.
       *
       * The inverse degree matrix is an n x n matrix such that each diagonal entry is
       * the inverse of the degree of the corresponding node, or 0 if the node
       * has degree 0.
       *
       * @return a sparse Eigen matrix
       */
      const SprsMat* inverse_degree_matrix();

      /**
       * The lazy random walk matrix of the graph.
       *
       * The lazy random walk matrix is defined to be
       *    1/2 I + 1/2 A D^{-1}
       * where I is the identity matrix, A is the graph adjacency matrix and
       * D is the degree matrix of the graph.
       *
       * @return a sparse Eigen matrix
       */
      const SprsMat* lazy_random_walk_matrix();

      /**
       * The total volume of the graph.
       *
       * The volume is defined as the sum of the node degrees.
       *
       * @return the graph's volume.
       */
      double total_volume();

      /**
       * The number of vertices in the graph.
       */
      stag_int number_of_vertices() const;

      /**
       * The number of edges in the graph.
       *
       * This is defined based on the number of non-zero elements in the
       * adjacency matrix, and ignores the weights of the edges.
       */
       stag_int number_of_edges() const;

       // Override the abstract methods in the LocalGraph base class.
       double degree(stag_int v) override;
       stag_int degree_unweighted(stag_int v) override;
       std::vector<edge> neighbors(stag_int v) override;
       std::vector<stag_int> neighbors_unweighted(stag_int v) override;
       ~Graph() override = default;

    private:
      /**
       * Initialise the laplacian matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_laplacian_();

      /**
       * Initialise the normalised Laplacian matrix of the graph is it has not
       * been initialised yet.
       */
      void initialise_normalised_laplacian_();

      /**
       * Initialise the degree matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_degree_matrix_();

      /**
       * Initialise the inverse degree matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_inverse_degree_matrix_();

      /**
       * Initialise the lazy random walk matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_lazy_random_walk_matrix_();

      /**
       * Check that the graph conforms to all assumptions that are currently
       * made within the library.
       *
       * @throws std::domain_error if the graph is not formatted correctly
       */
      void self_test_();

      // The number of vertices in the constructed graph.
      stag_int number_of_vertices_;

      // The ground truth definition of the graph object is the adjacency
      // matrix, stored in a sparse format. The adj_init_ variable is used to
      // indicate whether the matrix has been initialised yet.
      SprsMat adjacency_matrix_;

      // The laplacian matrix of the graph. The lap_init_ variable is used to
      // indicate whether the matrix has been initialised yet.
      bool lap_init_;
      SprsMat laplacian_matrix_;

      // The normalised Laplacian matrix of the graph. The norm_lap_init_
      // variable is used to indicate whether the matrix has been initialised
      // yet.
      bool norm_lap_init_;
      SprsMat normalised_laplacian_matrix_;

      // The degree matrix of the graph. The deg_init_ variable is used to
      // indicate whether the matrix has been initialised yet.
      bool deg_init_;
      SprsMat degree_matrix_;

      // The inverse degree matrix of the graph. The inv_deg_init_ variable is used to
      // indicate whether the matrix has been initialised yet.
      bool inv_deg_init_;
      SprsMat inverse_degree_matrix_;

      // The lazy random walk matrix of the graph. The lazy_rand_walk_init_ variable
      // is used to indicate whether the matrix has been initialised yet.
      bool lazy_rand_walk_init_;
      SprsMat lazy_random_walk_matrix_;
  };

  /**
   * Define equality for two graphs. Two graphs are equal iff their adjacency
   * matrices are equal
   */
  bool operator==(const Graph& lhs, const Graph& rhs);
  bool operator!=(const Graph& lhs, const Graph& rhs);

  /**
   * Define equality for edges.
   */
  bool operator==(const edge& lhs, const edge& rhs);
  bool operator!=(const edge& lhs, const edge& rhs);

  /**
   * Construct a cycle graph on n vertices.
   *
   * @param n
   * @return a graph object representing the n-cycle
   */
  Graph cycle_graph(stag_int n);

  /**
   * Construct a complete graph on n vertices.
   *
   * @param m
   * @return a graph object
   */
  Graph complete_graph(stag_int n);

  /**
   * Construct a barbell graph. The barbell graph consists of 2 cliques on n
   * vertices, connected by a single edge.
   *
   * @param n
   * @return
   */
  Graph barbell_graph(stag_int n);
}
#endif //STAG_LIBRARY_H
