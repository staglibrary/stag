#ifndef STAG_LIBRARY_H
#define STAG_LIBRARY_H

#include <Eigen/Sparse>
#include <vector>

// The fundamental datatype used in this library is the sparse matrix. For
// convenience, we define the sparse matrix type here.
#define SprsMat Eigen::SparseMatrix<double, Eigen::RowMajor>

namespace stag {

  /**
   * The core object used to represent a graph for use with the library. Graphs
   * are always constructed from sparse matrices, and this is the internal
   * representation used as well.
   */
  class Graph {
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
      const SprsMat* adjacency();

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
