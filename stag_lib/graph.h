/*
   This file is provided as part of the STAG library and released under the GPL
   license.
*/

/**
 * @file graph.h
 * \brief Core graph class definitions and constructors.
 */

#ifndef STAG_LIBRARY_H
#define STAG_LIBRARY_H

#include <vector>
#include <fstream>
#include <unordered_map>

#include "definitions.h"


namespace stag {
  /**
   * \brief A structure representing a weighted edge in a graph.
   */
  struct edge {
    /**
     * The first vertex in the edge.
     */
    StagInt v1;

    /**
     * The second vertex in the edge.
     */
    StagInt v2;

    /**
     * The weight of the edge.
     */
    StagReal weight;
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
       *
       * A self-loop of weight \f$1\f$ contributes \f$2\f$ to the vertex's
       * degree.
       */
      virtual StagReal degree(StagInt v) = 0;

      /**
       * Given a vertex v, return its unweighted degree. That is, the number
       * of neighbors of v, ignoring the edge weights.
       */
      virtual StagInt degree_unweighted(StagInt v) = 0;

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
      virtual std::vector<edge> neighbors(StagInt v) = 0;

      /**
       * Given a vertex v, return a vector containing the neighbors of v.
       *
       * The weights of edges to the neighbors are not returned by this method.
       *
       * @param v an int representing some vertex in the graph
       * @return an int vector giving the neighbors of v
       */
      virtual std::vector<StagInt> neighbors_unweighted(StagInt v) = 0;

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
      virtual std::vector<StagReal> degrees(std::vector<StagInt> vertices) = 0;

      /**
       * Given a list of vertices, return their unweighted degrees.
       *
       * @param vertices a vector of ints representing the vertices to be
       *                 queried.
       * @return a vector of integer degrees
       */
      virtual std::vector<StagInt> degrees_unweighted(std::vector<StagInt> vertices) = 0;

      /**
       * Given a vertex ID, returns true or false to indicate whether the vertex exists
       * in the graph.
       *
       * @param v the vertex index to check
       * @return a boolean indicating whether there exists a vertex with the given index
       */
       virtual bool vertex_exists(StagInt v) = 0;

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
   *
   * This class supports graphs with positive edge weights. Self-loops are
   * permitted.
   */
  class Graph : public LocalGraph {
    public:
      /**
       * Create a graph from an Eigen matrix.
       *
       * The provided matrix should correspond either to the
       * adjacency matrix or Laplacian matrix of the graph. STAG will
       * automatically detect whether the provided matrix is an adjacency matrix
       * or a Laplacian matrix.
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
       *   StagInt n = 3;
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
       * The provided matrix must be symmetric, and may include
       * self-loops.
       *
       * @param matrix the sparse eigen matrix representing the adjacency matrix
       *               or Laplacian matrix of the graph.
       * @throws domain_error if the provided matrix is not symmetric
       */
      explicit Graph(const SprsMat& matrix);

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
      Graph(std::vector<StagInt> &outerStarts, std::vector<StagInt> &innerIndices,
            std::vector<StagReal> &values);

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
      StagReal total_volume();

      /**
       * The average degree of the graph.
       *
       * This is defined as the sum of the node degrees divided by the number of nodes.
       *
       * @return the graph's average degree.
       */
      StagReal average_degree();

      /**
       * The number of vertices in the graph.
       */
      StagInt number_of_vertices() const;

      /**
       * The number of edges in the graph.
       *
       * This is defined based on the number of non-zero elements in the
       * adjacency matrix, and ignores the weights of the edges.
       */
       StagInt number_of_edges() const;

       /**
        * Add an edge to the graph.
        *
        * The edge goes from node i to node j, and is added with weight w.
        * If there is already an edge from i to j, then w is added to
        * its weight.
        *
        * If either of \f$i\f$ or \f$j\f$ are larger than the number
        * of nodes in the graph, the graph is resized to have enough nodes.
        *
        * @param i
        * @param j
        * @param w
        */
       void add_edge(StagInt i, StagInt j, StagReal w);

       /**
        * Remove an edge from the graph.
        *
        * Remove any edge between nodes \f$i\f$ and \f$j\f$.
        *
        * @param i
        * @param j
        */
       void remove_edge(StagInt i, StagInt j);

       /**
        * Returns a boolean indicating whether this graph contains self loops.
        */
       bool has_self_loops() const;

       /**
        * Returns a boolean indicating whether the graph is connected.
        *
        * The running time of this method is \f$O(m)\f$, where \f$m\f$ is the
        * number of edges in the graph.
        */
       bool is_connected();

       /**
        * Construct and return a subgraph of this graph.
        *
        * Note that the vertex indices will be changed in the subgraph.
        *
        * @param vertices the vertices in the induced subgraph
        * @return a new stag::Graph object representing the subgraph induced by
        *         the given vertices
        */
       Graph subgraph(std::vector<StagInt>& vertices);

       /**
        * Construct and return the disjoint union of this graph and another.
        *
        * The disjoint union of two graphs \f$G\f$ and \f$H\f$ is a graph
        * containing \f$G\f$ and \f$H\f$ as disconnected subgraphs.
        *
        * @param other the other graph to be combined with this one
        * @return a new stag::Graph object representing the union of this graph
        *         with the other one
        */
       Graph disjoint_union(Graph& other);

       // Override the abstract methods in the LocalGraph base class.
       StagReal degree(StagInt v) override;
       StagInt degree_unweighted(StagInt v) override;
       std::vector<edge> neighbors(StagInt v) override;
       std::vector<StagInt> neighbors_unweighted(StagInt v) override;
       std::vector<StagReal> degrees(std::vector<StagInt> vertices) override;
       std::vector<StagInt> degrees_unweighted(std::vector<StagInt> vertices) override;
       bool vertex_exists(StagInt v) override;
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
       void check_vertex_argument(StagInt v) const;

       /**
        * \endcond
        */

      // The number of vertices in the constructed graph.
      StagInt number_of_vertices_;

      // The ground truth definition of the graph object is the adjacency
      // matrix, stored in a sparse format. The adj_init_ variable is used to
      // indicate whether the matrix has been initialised yet.
      SprsMat adjacency_matrix_;

      // Whether the graph has self loops
      bool has_self_loops_;

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
   * \brief A local graph backed by an adjacency list file on disk.
   *
   * The graph is loaded into memory in a local way only. That is, an adjacency
   * list data structure is constructed in memory as node neighbours are queried.
   * If a node is not found in the cached adjacency list, then the neighbours of
   * the node are queried from the adjacency list on disk.
   * This allows for local algorithms to be executed on very large graphs stored
   * on disk without loading the whole graph into memory.
   *
   * See [Graph File Formats](@ref file-formats) for more information
   * about the adjacency list file format.
   *
   * \note
   * It is important that the adjacency list on disk is stored with sorted
   * node indices. This allows us to query the neighbours of a given node in
   * \f$O(log(n))\f$ time using binary search.
   *
   */
  class AdjacencyListLocalGraph : public LocalGraph {
  public:
    /**
     * Construct a local graph backed by an adjacency list file.
     *
     * The adjacency list file must not be modified externally while it is in
     * use by this object.
     *
     * @param filename the name of the adjacencylist file which defines the graph
     */
    AdjacencyListLocalGraph(const std::string& filename);

    // Override the abstract methods in the LocalGraph base class.
    StagReal degree(StagInt v) override;
    StagInt degree_unweighted(StagInt v) override;
    std::vector<edge> neighbors(StagInt v) override;
    std::vector<StagInt> neighbors_unweighted(StagInt v) override;
    std::vector<StagReal> degrees(std::vector<StagInt> vertices) override;
    std::vector<StagInt> degrees_unweighted(std::vector<StagInt> vertices) override;
    bool vertex_exists(StagInt v) override;
    ~AdjacencyListLocalGraph() override;

  private:
    /**
     * Move the ifstream head to the start of the next content line, and return
     * the ID of the corresponding node.
     *
     * If there is no content line before the end of the file, return -1 and
     * set the ifstream head to the end of the file.
     *
     * @return the node ID of the next content line
     */
    StagInt goto_next_content_line();

    /**
     * Perform a binary search to find the given vertex in the adjacencylist
     * file.
     *
     * If the vertex v is defined, the input stream pointer will be pointing
     * at the start of the corresponding content line.
     *
     * If the vertex does not exist, will throw a runtime exception.
     *
     * @param v the vertex to search for
     */
    void find_vertex(StagInt v);

    // The input file stream corresponding to the adjacencylist file backing
    // this graph. The implementation makes random access to this file to
    // read the vertex adjacency information.
    std::ifstream is_;
    std::streampos end_of_file_;

    // In order to increase the efficiency of looking up neighbourhood
    // information in the graph, we cache the node ids corresponding to certain
    // locations in the file. The cached locations will correspond to the binary
    // search locations for the nodes we've queried.
    std::unordered_map<StagInt, StagInt> fileloc_to_node_id_;

    // We also store the full adjacency list of the graph queried so far.
    std::unordered_map<StagInt, std::vector<edge>> node_id_to_edgelist_;
  };

  /**
   * Construct a cycle graph on n vertices.
   *
   * @param n the number of vertices in the constructed graph
   * @return a stag::Graph object representing a cycle graph
   */
  stag::Graph cycle_graph(StagInt n);

  /**
   * Construct a complete graph on n vertices.
   *
   * @param n the number of vertices in the constructed graph
   * @return a stag::Graph object representing a complete graph
   */
  stag::Graph complete_graph(StagInt n);

  /**
   * Construct a barbell graph. The barbell graph consists of 2 cliques on n
   * vertices, connected by a single edge.
   *
   * @param n the number of vertices in each of the two cliques.
   *          The returned graph will have \f$2n\f$ vertices.
   * @return a stag::Graph object representing the barbell graph
   */
  stag::Graph barbell_graph(StagInt n);

  /**
   * Construct a star graph. The star graph consists of one central vertex
   * connected by an edge to n-1 outer vertices.
   *
   * @param n the number of vertices in the constructed graph
   * @return a stag::Graph object representing the star graph
   */
   stag::Graph star_graph(StagInt n);

  /**
   * Construct the 'identity graph'. The identity graph consists of \f$n\f$
   * vertices, each with a self-loop of weight 1.
   *
   * Both the adjacency matrix and Laplacian matrix of the identity graph are
   * equal to the identity matrix.
   *
   * @param n the number of vertices in the constructed graph
   * @return a stag::Graph object representing the identity graph
   */
  stag::Graph identity_graph(StagInt n);

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
   * Multiplying a graph by a scalar is equivalent to multiplying the weight
   * of each edge by the given value.
   *
   * For example, the following code creates a complete graph with edges of
   * weight 2.
   *
   * \code{.cpp}
   *     #include <stag/graph.h>
   *
   *     int main() {
   *       stag::Graph myGraph = 2 * stag::complete_graph(10);
   *
   *       return 0;
   *     }
   * \endcode
   *
   */
  template <typename Scalar>
  stag::Graph operator*(Scalar lhs, const stag::Graph& rhs) {
    const SprsMat new_adj = lhs * *rhs.adjacency();
    return stag::Graph(new_adj);
  }

  /**
   * \overload
   */
  template <typename Scalar>
  stag::Graph operator*(const stag::Graph& lhs, Scalar rhs) {
    const SprsMat new_adj = rhs * *lhs.adjacency();
    return stag::Graph(new_adj);
  }

  /**
   * Adding two graphs is equivalent to adding their adjacency matrices.
   *
   * The graphs must have the same number of vertices. For example, the following
   * code adds a complete graph and a cycle graph on \f$5\f$ vertices.
   *
   * \code{.cpp}
   *    #include <stag/graph.h>
   *
   *    int main() {
   *        stag::Graph myGraph = stag::complete_graph(5) + stag::cycle_graph(5);
   *
   *        return 0;
   *    }
   * \endcode
   *
   * @throws std::invalud_argument if the graphs have different sizes.
   */
  stag::Graph operator+(const stag::Graph& lhs, const stag::Graph& rhs);
}

#endif //STAG_LIBRARY_H
