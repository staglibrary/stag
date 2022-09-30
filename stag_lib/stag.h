#ifndef STAG_LIBRARY_H
#define STAG_LIBRARY_H

#include <Eigen/Sparse>
#include <vector>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SprsMat;

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
       * Create a graph from raw arrays describing a sparse matrix.
       *
       * @param vertices the number of vertices in the graph
       * @param rows the row indices of each non-zero element in the matrix
       * @param cols the column indices of each non-zero element in the matrix
       * @param vals the values of each non-zero element in the matrix
       */
       Graph(int vertices, std::vector<int> rows, std::vector<int> cols, std::vector<double> vals);

      /**
       * Construct the Laplacian matrix of the graph.
       *
       * The Laplacian matrix is defined by
       *   L = D - A
       * where D is the diagonal matrix of vertex degrees and A is the adjacency
       * matrix of the graph.
       *
       * @return a sparse Eigen matrix represengint the graph Laplacian
       */
      SprsMat laplacian();

      /*
       * The volume of the graph.
       *
       * The volume is defined as the sum of all of the node degrees.
       *
       * @return the graph's volume.
       */
      double volume();

    private:
      SprsMat adjacency_matrix_;
  };
}
#endif //STAG_LIBRARY_H
