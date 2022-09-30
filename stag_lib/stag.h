#ifndef STAG_LIBRARY_H
#define STAG_LIBRARY_H

#include <Eigen/Core>

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
      explicit Graph(Eigen::MatrixXd adjacency_matrix);

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
      Eigen::MatrixXd laplacian();

    private:
      Eigen::MatrixXd adjacency_matrix_;
  };
}
#endif //STAG_LIBRARY_H
