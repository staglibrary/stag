/**
 * Definitions related to the core Graph object used to represent graphs within
 * the library.
 *
 * Copyright 2022 Peter Macgregor
 */
#ifndef STAG_LIBRARY_H
#define STAG_LIBRARY_H

#include <Eigen/Sparse>
#include <vector>

// The fundamental datatype used in this library is the sparse matrix. For
// convenience, we define the sparse matrix type here.
#define SprsMat Eigen::SparseMatrix<double>

namespace stag {

  /**
   * A structure representing a weighted edge in a graph.
   */
  struct edge {
    // The first vertex in the edge.
    int v;

    // The second vertex in the edge
    int u;

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
      virtual double degree(int v) = 0;

      /**
       * Given a vertex v, return its unweighted degree. That is, the number
       * of neighbors of v, ignoring the edge weights.
       */
       virtual int degree_unweighted(int v) = 0;

      /**
       * Given a vertex v, return a vector of edges representing the
       * neighborhood of v.
       *
       * @param v an int representing some vertex in the graph
       * @return an edge vector containing the neighborhood information
       */
      virtual std::vector<edge> neighbors(int v) = 0;

      /**
       * Given a vertex v, return a vector containing the neighbors of v.
       *
       * The weights of edges to the neighbors are not returned by this method.
       *
       * @param v an int representing some vertex in the graph
       * @return an int vector giving the neighbors of v
       */
      virtual std::vector<int> neighbors_unweighted(int v) = 0;
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
      Graph(std::vector<int> &outerStarts, std::vector<int> &innerIndices,
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
      long number_of_vertices() const;

      /**
       * The number of edges in the graph.
       *
       * This is defined based on the number of non-zero elements in the
       * adjacency matrix, and ignores the weights of the edges.
       */
       long number_of_edges() const;

       // Override the abstract methods in the LocalGraph base class.
       double degree(int v) override;
       int degree_unweighted(int v) override;
       std::vector<edge> neighbors(int v) override;
       std::vector<int> neighbors_unweighted(int v) override;

    private:
      /**
       * Initialise the laplacian matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_laplacian_();

      /**
       * Initialise the normalised Laplacian matrix of the grpah is it has not
       * been initialised yet.
       */
      void initialise_normalised_laplacian_();

      /**
       * Initialise the degree matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_degree_matrix_();

      /**
       * Check that the graph conforms to all assumptions that are currently
       * made within the library.
       *
       * @throws std::domain_error if the graph is not formatted correctly
       */
      void self_test_();

      // The number of vertices in the constructed graph.
      long number_of_vertices_;

      // The ground truth definition of the graph object is the adjacency
      // matrix, stored in a sparse format. The adj_init_ variable is used to
      // indicate whether the matrix has been initialised yet.
      bool adj_init_;
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
  Graph cycle_graph(int n);

  /**
   * Construct a complete graph on n vertices.
   *
   * @param m
   * @return a graph object
   */
  Graph complete_graph(int n);
}
#endif //STAG_LIBRARY_H
