/**
 * Tests for the methods in the graph.h header file. Includes the main Graph
 * object.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <stdexcept>
#include <gtest/gtest.h>
#include <graph.h>
#include <utility.h>

// Define some helper test assertions.
#define EXPECT_FLOATS_NEARLY_EQ(expected, actual, thresh) \
        EXPECT_EQ(expected.size(), actual.size()) << "Array sizes differ.";\
        for (size_t idx = 0; idx < std::min(expected.size(), actual.size()); ++idx) \
        { \
            EXPECT_NEAR(expected[idx], actual[idx], thresh) << "at index: " << idx;\
        }

/**
 * Create a useful test graph which we will use repeatedly.
 *
 * @return a stag::Graph object to be tested
 */
stag::Graph createTestGraph() {
  // Create the data for the graph adjacency matrix.
  //     0       2  3.3333 0
  //     2       0   6     0
  //     3.3333  6   0     1
  //     0       0   1     0
  std::vector<stag_int> rowStarts = {0, 2, 4, 7, 8};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};

  // Create the stag Graph object
  return {rowStarts, colIndices, values};
}

TEST(GraphTest, Volume) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.total_volume(), 24.6666);
}

TEST(GraphTest, NumberOfVertices) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.number_of_vertices(), 4);
}

TEST(GraphTest, NumberOfEdges) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.number_of_edges(), 4);
}

TEST(GraphTest, UnweightedDegree) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.degree_unweighted(0), 2);
  EXPECT_EQ(testGraph.degree_unweighted(1), 2);
  EXPECT_EQ(testGraph.degree_unweighted(2), 3);
  EXPECT_EQ(testGraph.degree_unweighted(3), 1);

  // The unweighted degree of a non-existent vertex is 0.
  EXPECT_EQ(testGraph.degree_unweighted(100), 0);
}

TEST(GraphTest, Degree) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_NEAR(testGraph.degree(0), 5.3333, 0.000001);
  EXPECT_NEAR(testGraph.degree(1), 8, 0.000001);
  EXPECT_NEAR(testGraph.degree(2), 10.3333, 0.000001);
  EXPECT_NEAR(testGraph.degree(3), 1, 0.000001);

  // The degree of a non-existent vertex is 0.
  EXPECT_EQ(testGraph.degree(100), 0);
}

TEST(GraphTest, UnweightedNeighbors) {
  stag::Graph testGraph = createTestGraph();

  std::vector<stag_int> expectedNeighbors = {1, 2};
  EXPECT_EQ(testGraph.neighbors_unweighted(0), expectedNeighbors);

  expectedNeighbors = {0, 2};
  EXPECT_EQ(testGraph.neighbors_unweighted(1), expectedNeighbors);

  expectedNeighbors = {0, 1, 3};
  EXPECT_EQ(testGraph.neighbors_unweighted(2), expectedNeighbors);

  expectedNeighbors = {2};
  EXPECT_EQ(testGraph.neighbors_unweighted(3), expectedNeighbors);
}

TEST(GraphTest, Neighbors) {
  stag::Graph testGraph = createTestGraph();

  std::vector<stag::edge> expectedNeighbors = {{0, 1, 2}, {0, 2, 3.3333}};
  EXPECT_EQ(testGraph.neighbors(0), expectedNeighbors);

  expectedNeighbors = {{1, 0, 2}, {1, 2, 6}};
  EXPECT_EQ(testGraph.neighbors(1), expectedNeighbors);

  expectedNeighbors = {{2, 0, 3.3333}, {2, 1, 6}, {2, 3, 1}};
  EXPECT_EQ(testGraph.neighbors(2), expectedNeighbors);

  expectedNeighbors = {{3, 2, 1}};
  EXPECT_EQ(testGraph.neighbors(3), expectedNeighbors);
}

TEST(GraphTest, AdjacencyMatrix) {
  // Create the test grpah object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 7, 8};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, AsymmetricAdjacency) {
  // Creating a graph with an asymetric adjacency matrix is not currently
  // supported and should throw an error.

  // Create the data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 7, 8};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3, 6, 1, 1};

  // Attempt to create the graph object. Should throw an exception.
  EXPECT_THROW({stag::Graph testGraph(rowStarts, colIndices, values);},
               std::domain_error);
}

TEST(GraphTest, LaplacianMatrix) {
  // Create the test graph object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the graph laplacian matrix.
  std::vector<stag_int> rowStarts = {0, 3, 6, 10, 12};
  std::vector<stag_int> colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> values = {5.3333, -2, -3.3333, -2, 8, -6, -3.3333, -6, 10.3333, -1, -1, 1};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, NormalisedLaplacianMatrix) {
  // Create the test graph object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the normalised graph laplacian matrix.
  std::vector<stag_int> rowStarts = {0, 3, 6, 10, 12};
  std::vector<stag_int> colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> values = {1, -2/(sqrt(5.3333) * sqrt(8)), -3.3333/(sqrt(5.3333) * sqrt(10.3333)),
                                -2/(sqrt(8) * sqrt(5.3333)), 1, -6/(sqrt(8) * sqrt(10.3333)),
                                -3.3333/(sqrt(10.3333) * sqrt(5.3333)), -6/(sqrt(10.3333) * sqrt(8)), 1, -1/sqrt(10.3333),
                                -1/sqrt(10.3333), 1};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.normalised_laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.normalised_laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.normalised_laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, DegreeMatrix) {
  // Create the test graph object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the graph degree matrix.
  std::vector<stag_int> rowStarts = {0, 1, 2, 3, 4};
  std::vector<stag_int> colIndices = {0, 1, 2, 3};
  std::vector<double> values = {5.3333, 8, 10.3333, 1};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.degree_matrix());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.degree_matrix());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.degree_matrix());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, InverseDegreeMatrix) {
  // Create the test graph object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the inverse degree matrix.
  std::vector<stag_int> rowStarts = {0, 1, 2, 3, 4};
  std::vector<stag_int> colIndices = {0, 1, 2, 3};
  std::vector<double> values = {1./5.3333, 1./8, 1./10.3333, 1./1};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.inverse_degree_matrix());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.inverse_degree_matrix());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.inverse_degree_matrix());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, LazyRandomWalkMatrix) {
  // Create the test graph object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the lazy random walk matrix.
  std::vector<stag_int> colStarts = {0, 3, 6, 10, 12};
  std::vector<stag_int> rowIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> values = {1./2, 1./5.3333, 3.3333/10.6666,
                                1./8, 1./2, 3./8,
                                3.3333/20.6666, 3./10.3333, 1./2, 1./20.6666,
                                1./2, 1./2};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.lazy_random_walk_matrix());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.lazy_random_walk_matrix());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.lazy_random_walk_matrix());

  EXPECT_EQ(colStarts, newStarts);
  EXPECT_EQ(rowIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, CycleGraphVolume) {
  // The volume of the cycle graph should be twice the number of vertices
  std::vector<stag_int> sizes = {3, 3, 5, 10, 20, 100};
  for (stag_int n: sizes) {
    stag::Graph testGraph = stag::cycle_graph(n);
    EXPECT_EQ(testGraph.total_volume(), 2 * n);
  }
}

TEST(GraphTest, CycleGraphLaplacian) {
  // Create a small cycle graph
  stag::Graph testGraph = stag::cycle_graph(4);

  // Define the expected laplacian matrix
  std::vector<stag_int> rowStarts = {0, 3, 6, 9, 12};
  std::vector<stag_int> colIndices = {0, 1, 3, 0, 1, 2, 1, 2, 3, 0, 2, 3};
  std::vector<double> values = {2, -1, -1, -1, 2, -1, -1, 2, -1, -1, -1, 2};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, CycleGraphDegrees) {
  // Create a small cycle graph
  stag::Graph testGraph = stag::cycle_graph(4);

  // Define the expected degree matrix
  std::vector<stag_int> rowStarts = {0, 1, 2, 3, 4};
  std::vector<stag_int> colIndices = {0, 1, 2, 3};
  std::vector<double> values = {2, 2, 2, 2};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.degree_matrix());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.degree_matrix());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.degree_matrix());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, CompleteGraphVolume) {
  // The volume of the complete graph should be n(n-1)
  std::vector<stag_int> sizes = {3, 3, 5, 10, 20, 100};
  for (stag_int n: sizes) {
    stag::Graph testGraph = stag::complete_graph(n);
    EXPECT_EQ(testGraph.total_volume(), n * (n - 1));
  }
}

TEST(GraphTest, CompleteGraphLaplacian) {
  // Create a small complete graph
  stag::Graph testGraph = stag::complete_graph(4);

  // Define the expected laplacian matrix
  std::vector<stag_int> rowStarts = {0, 4, 8, 12, 16};
  std::vector<stag_int> colIndices = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  std::vector<double> values = {3, -1, -1, -1, -1, 3, -1, -1, -1, -1, 3, -1, -1, -1, -1, 3};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, CompleteGraphNormalisedLaplacian) {
  // Create a small complete graph
  stag::Graph testGraph = stag::complete_graph(4);

  // Define the expected normalised laplacian matrix
  std::vector<stag_int> rowStarts = {0, 4, 8, 12, 16};
  std::vector<stag_int> colIndices = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  std::vector<double> values = {1, -1./3, -1./3, -1./3,
                                -1./3, 1., -1./3, -1./3,
                                -1./3, -1./3, 1, -1./3,
                                -1./3, -1./3, -1./3, 1};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.normalised_laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.normalised_laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.normalised_laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, BarbellGraph) {
  // Create a small complete graph
  stag::Graph testGraph = stag::barbell_graph(4);

  // Define the expected laplacian matrix
  std::vector<stag_int> colStarts = {0, 3, 6, 9, 13, 17, 20, 23, 26};
  std::vector<stag_int> rowIndices = {1, 2, 3,
                                      0, 2, 3,
                                      0, 1, 3,
                                      0, 1, 2, 4,
                                      3, 5, 6, 7,
                                      4, 6, 7,
                                      4, 5, 7,
                                      4, 5, 6};
  std::vector<double> values = {1, 1, 1,
                                1, 1, 1,
                                1, 1, 1,
                                1, 1, 1, 1,
                                1, 1, 1, 1,
                                1, 1, 1,
                                1, 1, 1,
                                1, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(colStarts, newStarts);
  EXPECT_EQ(rowIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, Equality) {
  // Create two identical graphs
  stag::Graph graph1 = stag::cycle_graph(10);
  stag::Graph graph2 = stag::cycle_graph(10);
  EXPECT_EQ(graph1, graph2);

  // Create another two identical graphs
  graph1 = createTestGraph();
  graph2 = createTestGraph();
  EXPECT_EQ(graph1, graph2);

  // Create two different graphs;
  graph1 = createTestGraph();
  graph2 = stag::complete_graph(4);
  EXPECT_NE(graph1, graph2);
}
