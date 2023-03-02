//
// Definitions related to the core Graph object used to represent graphs within
// the library.
//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file graph.h
 * \brief Core graph class definitions and constructors.
 */

#ifndef STAG_LIBRARY_H
#define STAG_LIBRARY_H

#include <Eigen/Sparse>
#include <vector>


/**
 * The integer type used throughout the library.
 */
typedef long long stag_int;

/**
 * The fundamental datatype used in this library is the sparse matrix.
 * We use the `Eigen::SparseMatrix` class in column-major format.
 */
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, stag_int> SprsMat;

/**
 *  Occasionally, it is more efficient to use a dense matrix, such as when
 *  the matrix is very small.
 *
 *  In this case, we use the `Eigen::MatrixXd` class.
 */
typedef Eigen::MatrixXd DenseMat;

/**
 * An Eigen::Triplet representing an edge in a graph. Stores the row, column, and value
 * of an entry in a graph adjacency matrix.
 */
typedef Eigen::Triplet<double, stag_int> EdgeTriplet;

/**
 * \cond
 */
// Redefine the eigen index type to be the same as stag_int
#undef EIGEN_DEFAULT_DENSE_INDEX_TYPE
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE stag_int
/**
 * \endcond
 */

namespace stag {
  /**
   * \brief A structure representing a weighted edge in a graph.
   */
  struct edge {
    /**
     * The first vertex in the edge.
     */
    stag_int v1;

    /**
     * The second vertex in the edge.
     */
    stag_int v2;

    /**
     * The weight of the edge.
     */
    double weight;
  };

  /**
   * \brief An abstract class which defines methods for exploring the
   * local neighborhood of vertices in a graph.
   *
   * To maximise the performance of the local algorithms using this class,
   * subclasses should cache the results of expensive queries. For example,
   * if querying the neighbors of a vertex requires accessing the disk, then
   * the result should be cached.
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
       * Given a list of vertices, return the degrees of each vertex in the
       * list.
       *
       * Providing a method for computing the degrees 'in bulk' increases the
       * efficiency of algorithms on graphs which are not stored directly in
       * memory.
       *
       * @param vertices a vector of ints representing the vertices to be
       *                 queried.
       * @return a vector of degrees
       */
      virtual std::vector<double> degrees(std::vector<stag_int> vertices) = 0;

      /**
       * Given a list of vertices, return their unweighted degrees.
       *
       * @param vertices a vector of ints representing the vertices to be
       *                 queried.
       * @return a vector of integer degrees
       */
      virtual std::vector<stag_int> degrees_unweighted(std::vector<stag_int> vertices) = 0;

      /**
       * Given a vertex ID, returns true or false to indicate whether the vertex exists
       * in the graph.
       *
       * @param v the vertex index to check
       * @return a boolean indicating whether there exists a vertex with the given index
       */
       virtual bool vertex_exists(stag_int v) = 0;

      /**
       * Destructor for the LocalGraph object.
       */
      virtual ~LocalGraph() = default;
  };

  /**
   * \brief The core object used to represent graphs for use with the library.
   *
   * Graphs are always constructed from sparse matrices, and this is the internal
   * representation used as well.
   * Vertices of the graph are always referred to by their unique integer index.
   * This index corresponds to the position of the vertex in the stored adjacency
   * matrix of the graph.
   */
  class Graph : public LocalGraph {
    public:
      /**
       * Create a graph from an Eigen matrix.
       *
       * \par Example
       *
       * \code{cpp}
       * #include <iostream>
       * #include <stag/graph.h>
       *
       * int main() {
       *   // Construct a sparse matrix representing the
       *   // triangle graph adjacency matrix.
       *   stag_int n = 3;
       *   SprsMat adj(n, n);
       *   adj.coeffRef(0, 1) = 1;
       *   adj.coeffRef(0, 2) = 1;
       *   adj.coeffRef(1, 0) = 1;
       *   adj.coeffRef(1, 2) = 1;
       *   adj.coeffRef(2, 0) = 1;
       *   adj.coeffRef(2, 1) = 1;
       *
       *   // Create a new STAG graph
       *   stag::Graph myGraph(adj);
       *
       *   // Display the adjacency matrix of the graph
       *   std::cout << *myGraph.adjacency() << std::endl;
       *
       *   return 0;
       * }
       * \endcode
       *
       * The provided adjacency matrix must be symmetric.
       *
       * @param adjacency_matrix the sparse eigen matrix representing the adjacency matrix
       *               of the graph.
       * @throws domain_error if the adjacency matrix is not symmetric
       */
      explicit Graph(const SprsMat& adjacency_matrix);

      /**
       * Create a graph from raw arrays describing a CSC sparse matrix.
       *
       * To use this constructor, you should understand the CSC sparse matrix
       * format. Note that this library uses the ColMajor format from the Eigen
       * library.
       * For more information, refer to the
       * [Eigen Documentation](https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html).
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
       * Return the sparse adjacency matrix of the graph.
       *
       * @return a sparse Eigen matrix representing the graph adjacency matrix.
       */
      const SprsMat* adjacency() const;

      /**
       * Return the Laplacian matrix of the graph.
       *
       * The Laplacian matrix is defined by
       *
       * \f[
       *   L = D - A
       * \f]
       *
       * where \f$D\f$ is the diagonal matrix of vertex degrees
       * (stag::Graph::degree_matrix)
       * and A is the adjacency matrix of the graph
       * (stag::Graph::adjacency).
       *
       * @return a sparse Eigen matrix representing the graph Laplacian
       */
      const SprsMat* laplacian();

      /**
       * Return the normalised Laplacian matrix of the graph.
       *
       * The normalised Laplacian matrix is defined by
       *
       * \f[
       *   \mathcal{L} = D^{-1/2} L D^{-1/2}
       * \f]
       *
       * where \f$D\f$ is the diagonal matrix of vertex degrees
       * (stag::Graph::degree_matrix)
       * and \f$L\f$ is the Laplacian matrix of the graph
       * (stag::Graph::laplacian).
       *
       * @return a sparse Eigen matrix representing the normalised Laplacian
       */
      const SprsMat* normalised_laplacian();

      /**
       * Return the signless Laplacian matrix of the graph.
       *
       * The signless Laplacian matrix is defined by
       *
       * \f[
       *   J = D + A
       * \f]
       *
       * where \f$D\f$ is the diagonal matrix of vertex degrees
       * (stag::Graph::degree_matrix)
       * and A is the adjacency matrix of the graph
       * (stag::Graph::adjacency).
       *
       * @return a sparse Eigen matrix representing the signless graph Laplacian
       */
      const SprsMat* signless_laplacian();

      /**
       * Return the normalised signless Laplacian matrix of the graph.
       *
       * The normalised signless Laplacian matrix is defined by
       *
       * \f[
       *   \mathcal{J} = D^{-1/2} J D^{-1/2}
       * \f]
       *
       * where \f$D\f$ is the diagonal matrix of vertex degrees
       * (stag::Graph::degree_matrix)
       * and \f$J\f$ is the signless Laplacian matrix of the graph
       * (stag::Graph::signless_laplacian).
       *
       * @return a sparse Eigen matrix representing the normalised Laplacian
       */
      const SprsMat* normalised_signless_laplacian();

      /**
       * The degree matrix of the graph.
       *
       * The degree matrix is an \f$n \times n\f$ matrix such that each diagonal entry is
       * the degree of the corresponding node.
       *
       * @return a sparse Eigen matrix
       */
      const SprsMat* degree_matrix();

      /**
       * The inverse degree matrix of the graph.
       *
       * The inverse degree matrix is an \f$n \times n\f$ matrix such that each diagonal entry is
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
       *
       * \f[
       *    \frac{1}{2} I + \frac{1}{2} A D^{-1}
       * \f]
       *
       * where \f$I\f$ is the identity matrix, \f$A\f$ is the graph adjacency matrix and
       * \f$D\f$ is the degree matrix of the graph.
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
       * The average degree of the graph.
       *
       * This is defined as the sum of the node degrees divided by the number of nodes.
       *
       * @return the graph's average degree.
       */
      double average_degree();

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
       std::vector<double> degrees(std::vector<stag_int> vertices) override;
       std::vector<stag_int> degrees_unweighted(std::vector<stag_int> vertices) override;
       bool vertex_exists(stag_int v) override;
       ~Graph() override = default;

    private:
      /**
       * Initialise the laplacian matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_laplacian_();

      /**
       * Initialise the signless Laplacian matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_signless_laplacian_();

    /**
       * Initialise the normalised Laplacian matrix of the graph is it has not
       * been initialised yet.
       */
      void initialise_normalised_laplacian_();

      /**
       * Initialise the signless Laplacian matrix of the graph if it has not been
       * initialised yet.
       */
      void initialise_normalised_signless_laplacian_();

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

      /**
       * \cond
       * Do not document the check vertex argument method.
       */

      /**
       * Check the validity of a method argument which is supposed to refer
       * to a vertex in the graph.
       *
       * @throws std::invalid_argument if the check does not pass
       */
       void check_vertex_argument(stag_int v);

       /**
        * \endcond
        */

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

      // The signless Laplacian matrix of the graph.
      bool signless_lap_init_;
      SprsMat signless_laplacian_matrix_;

      // The normalised Laplacian matrix of the graph. The norm_lap_init_
      // variable is used to indicate whether the matrix has been initialised
      // yet.
      bool norm_lap_init_;
      SprsMat normalised_laplacian_matrix_;

      // The normalised signless Laplacian matrix of the graph.
      bool signless_norm_lap_init_;
      SprsMat normalised_signless_laplacian_matrix_;

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
   * \cond
   * Do not generate documentation for operator definitions.
   */

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
   * \endcond
   */

  /**
   * Construct a cycle graph on n vertices.
   *
   * @param n the number of vertices in the constructed graph
   * @return a stag::Graph object representing a cycle graph
   */
  stag::Graph cycle_graph(stag_int n);

  /**
   * Construct a complete graph on n vertices.
   *
   * @param n the number of vertices in the constructed graph
   * @return a stag::Graph object representing a complete graph
   */
  stag::Graph complete_graph(stag_int n);

  /**
   * Construct a barbell graph. The barbell graph consists of 2 cliques on n
   * vertices, connected by a single edge.
   *
   * @param n the number of vertices in each of the two cliques.
   *          The returned graph will have \f$2n\f$ vertices.
   * @return a stag::Graph object representing the barbell graph
   */
  stag::Graph barbell_graph(stag_int n);

  /**
   * Construct a star graph. The star graph consists of one central vertex
   * connected by an edge to n-1 outer vertices.
   *
   * @param n the number of vertices in the constructed graph
   * @return a stag::Graph object representing the star graph
   */
   stag::Graph star_graph(stag_int n);

}
#endif //STAG_LIBRARY_H

