/**
 * Tests for the methods in the graph.h header file. Includes the main Graph
 * object.
 */
#include <stdexcept>
#include <gtest/gtest.h>
#include <graph.h>
#include <graphio.h>
#include <utility.h>
#include <random.h>
#include <cluster.h>
#include <stdio.h>

// Define some helper test assertions.
#define EXPECT_FLOATS_NEARLY_EQ(expected, actual, thresh) \
        EXPECT_EQ(expected.size(), actual.size()) << "Array sizes differ.";\
        for (size_t idx = 0; idx < std::min(expected.size(), actual.size()); ++idx) \
        { \
            EXPECT_NEAR(expected[idx], actual[idx], thresh) << "at index: " << idx;\
        }

#define EXPECT_GRAPHS_NEARLY_EQ(expected, actual) \
  EXPECT_FLOATS_NEARLY_EQ(stag::sprsMatInnerIndices(expected.adjacency()), \
                          stag::sprsMatInnerIndices(actual.adjacency()),   \
                          EPSILON)                                         \
  EXPECT_FLOATS_NEARLY_EQ(stag::sprsMatOuterStarts(expected.adjacency()), \
                          stag::sprsMatOuterStarts(actual.adjacency()),   \
                          EPSILON)                                         \
  EXPECT_FLOATS_NEARLY_EQ(stag::sprsMatValues(expected.adjacency()), \
                          stag::sprsMatValues(actual.adjacency()),   \
                          EPSILON)

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

/**
 * Create a simple test graph with self-loops.
 *
 * @return a stag::Graph object to be tested
 */
stag::Graph createTestGraphSelfLoops() {
  // Create the data for the graph adjacency matrix.
  //     1       2  3.3333 0
  //     2       0   6     0
  //     3.3333  6   2     1
  //     0       0   1     1
  std::vector<stag_int> rowStarts = {0, 3, 5, 9, 11};
  std::vector<stag_int> colIndices = {0, 1, 2, 0, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> values = {1, 2, 3.3333, 2, 6, 3.3333, 6, 2, 1, 1, 1};

  // Create the stag Graph object
  return {rowStarts, colIndices, values};
}

TEST(GraphTest, LaplacianInitialisation) {
  // Create a graph Laplacian matrix
  //     5.3333  -2  -3.3333   0
  //    -2        8  -6        0
  //    -3.3333  -6   10.3333 -1
  //     0        0  -1        1
  std::vector<stag_int> lapRowStarts = {0, 3, 6, 10, 12};
  std::vector<stag_int> lapColIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> lapValues = {5.3333, -2, -3.3333, -2, 8, -6, -3.3333, -6, 10.3333, -1, -1, 1};

  // Create the corresponding adjacency matrix
  //     0       2  3.3333 0
  //     2       0   6     0
  //     3.3333  6   0     1
  //     0       0   1     0
  std::vector<stag_int> adjRowStarts = {0, 2, 4, 7, 8};
  std::vector<stag_int> adjColIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> adjValues = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};

  // We will initialise the same graph in three ways:
  //    - SprsMat Laplacian matrix
  //    - Data vectors for sparse Laplacian matrix
  //    - Adjacency matrix
  // All three methods should result in the same graph.
  stag::Graph testGraph1 = stag::Graph(lapRowStarts, lapColIndices, lapValues);
  SprsMat lap = stag::sprsMatFromVectors(lapRowStarts, lapColIndices, lapValues);
  stag::Graph testGraph2 = stag::Graph(lap);
  stag::Graph testGraph3 = stag::Graph(adjRowStarts, adjColIndices, adjValues);

  // Check for graph equality
  EXPECT_GRAPHS_NEARLY_EQ(testGraph1, testGraph2);
  EXPECT_GRAPHS_NEARLY_EQ(testGraph2, testGraph3);
  EXPECT_GRAPHS_NEARLY_EQ(testGraph1, testGraph3);
}

TEST(GraphTest, LaplacianSelfLoopInitialisation) {
  // Create a graph Laplacian matrix, with self-loops
  //     6.3333  -2  -3.3333   0
  //    -2        8  -6        0
  //    -3.3333  -6   12.3333 -1
  //     0        0  -1        2
  std::vector<stag_int> lapRowStarts = {0, 3, 6, 10, 12};
  std::vector<stag_int> lapColIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> lapValues = {6.3333, -2, -3.3333, -2, 8, -6, -3.3333, -6, 12.3333, -1, -1, 2};

  // Create the corresponding adjacency matrix
  //     1       2   3.3333 0
  //     2       0   6      0
  //     3.3333  6   2      1
  //     0       0   1      1
  std::vector<stag_int> adjRowStarts = {0, 3, 5, 9, 11};
  std::vector<stag_int> adjColIndices = {0, 1, 2, 0, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> adjValues = {1, 2, 3.3333, 2, 6, 3.3333, 6, 2, 1, 1, 1};

  // We will initialise the same graph in three ways:
  //    - SprsMat Laplacian matrix
  //    - Data vectors for sparse Laplacian matrix
  //    - Adjacency matrix
  // All three methods should result in the same graph.
  stag::Graph testGraph1 = stag::Graph(lapRowStarts, lapColIndices, lapValues);
  SprsMat lap = stag::sprsMatFromVectors(lapRowStarts, lapColIndices, lapValues);
  stag::Graph testGraph2 = stag::Graph(lap);
  stag::Graph testGraph3 = stag::Graph(adjRowStarts, adjColIndices, adjValues);

  // Check for graph equality
  EXPECT_GRAPHS_NEARLY_EQ(testGraph1, testGraph2);
  EXPECT_GRAPHS_NEARLY_EQ(testGraph2, testGraph3);
  EXPECT_GRAPHS_NEARLY_EQ(testGraph1, testGraph3);
}

TEST(GraphTest, InvalidLaplacianInitialisation) {
  // Try to initialise a graph with a matrix which cannot be a Laplacian or
  // adjacency matrix.
  std::vector<stag_int> matRowStarts = {0, 2, 4, 6, 8};
  std::vector<stag_int> matColIndices = {1, 3, 0, 2, 1, 3, 0, 2};
  std::vector<double> matValues = {-1, 1, -1, 1, 1, 1, 1, 1};

  EXPECT_THROW(stag::Graph(matRowStarts, matColIndices, matValues),
               std::invalid_argument);
}

TEST(GraphTest, Volume) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.total_volume(), 24.6666);

  testGraph = createTestGraphSelfLoops();
  EXPECT_EQ(testGraph.total_volume(), 32.6666);
}

TEST(GraphTest, NumberOfVertices) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.number_of_vertices(), 4);

  testGraph = createTestGraphSelfLoops();
  EXPECT_EQ(testGraph.number_of_vertices(), 4);
}

TEST(GraphTest, NumberOfEdges) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.number_of_edges(), 4);

  testGraph = createTestGraphSelfLoops();
  EXPECT_EQ(testGraph.number_of_edges(), 7);
}

TEST(GraphTest, HasSelfLoops) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.has_self_loops(), false);

  testGraph = createTestGraphSelfLoops();
  EXPECT_EQ(testGraph.has_self_loops(), true);
}

TEST(GraphTest, UnweightedDegree) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_EQ(testGraph.degree_unweighted(0), 2);
  EXPECT_EQ(testGraph.degree_unweighted(1), 2);
  EXPECT_EQ(testGraph.degree_unweighted(2), 3);
  EXPECT_EQ(testGraph.degree_unweighted(3), 1);

  testGraph = createTestGraphSelfLoops();
  EXPECT_EQ(testGraph.degree_unweighted(0), 4);
  EXPECT_EQ(testGraph.degree_unweighted(1), 2);
  EXPECT_EQ(testGraph.degree_unweighted(2), 5);
  EXPECT_EQ(testGraph.degree_unweighted(3), 3);
}

TEST(GraphTest, Degree) {
  stag::Graph testGraph = createTestGraph();
  EXPECT_NEAR(testGraph.degree(0), 5.3333, 0.000001);
  EXPECT_NEAR(testGraph.degree(1), 8, 0.000001);
  EXPECT_NEAR(testGraph.degree(2), 10.3333, 0.000001);
  EXPECT_NEAR(testGraph.degree(3), 1, 0.000001);

  testGraph = createTestGraphSelfLoops();
  EXPECT_NEAR(testGraph.degree(0), 7.3333, 0.000001);
  EXPECT_NEAR(testGraph.degree(1), 8, 0.000001);
  EXPECT_NEAR(testGraph.degree(2), 14.3333, 0.000001);
  EXPECT_NEAR(testGraph.degree(3), 3, 0.000001);
}

TEST(GraphTest, OutOfRangeVertices) {
  // Check that we cannot request information for a vertex index that is out of range.
  stag::Graph testGraph = createTestGraph();

  EXPECT_THROW(testGraph.degree(4), std::invalid_argument);
  EXPECT_THROW(testGraph.degree(-1), std::invalid_argument);

  EXPECT_THROW(testGraph.degree_unweighted(4), std::invalid_argument);
  EXPECT_THROW(testGraph.degree_unweighted(-1), std::invalid_argument);

  EXPECT_THROW(testGraph.neighbors(4), std::invalid_argument);
  EXPECT_THROW(testGraph.neighbors(-1), std::invalid_argument);

  EXPECT_THROW(testGraph.neighbors_unweighted(4), std::invalid_argument);
  EXPECT_THROW(testGraph.neighbors_unweighted(-1), std::invalid_argument);

  // The degrees and degrees_unweighted methods take vectors of degrees
  std::vector<stag_int> vs1 = {0, 2, 4};
  EXPECT_THROW(testGraph.degrees(vs1), std::invalid_argument);

  std::vector<stag_int> vs2 = {0, -1, 3};
  EXPECT_THROW(testGraph.degrees_unweighted(vs2), std::invalid_argument);
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

  testGraph = createTestGraphSelfLoops();

  expectedNeighbors = {0, 1, 2};
  EXPECT_EQ(testGraph.neighbors_unweighted(0), expectedNeighbors);

  expectedNeighbors = {0, 2};
  EXPECT_EQ(testGraph.neighbors_unweighted(1), expectedNeighbors);

  expectedNeighbors = {0, 1, 2, 3};
  EXPECT_EQ(testGraph.neighbors_unweighted(2), expectedNeighbors);

  expectedNeighbors = {2, 3};
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

  testGraph = createTestGraphSelfLoops();

  expectedNeighbors = {{0, 0, 1}, {0, 1, 2}, {0, 2, 3.3333}};
  EXPECT_EQ(testGraph.neighbors(0), expectedNeighbors);

  expectedNeighbors = {{1, 0, 2}, {1, 2, 6}};
  EXPECT_EQ(testGraph.neighbors(1), expectedNeighbors);

  expectedNeighbors = {{2, 0, 3.3333}, {2, 1, 6}, {2, 2, 2}, {2, 3, 1}};
  EXPECT_EQ(testGraph.neighbors(2), expectedNeighbors);

  expectedNeighbors = {{3, 2, 1}, {3, 3, 1}};
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

  // Create the test grpah object
  testGraph = createTestGraphSelfLoops();

  // Create the expected data for the graph adjacency matrix.
  rowStarts = {0, 3, 5, 9, 11};
  colIndices = {0, 1, 2, 0, 2, 0, 1, 2, 3, 2, 3};
  values = {1, 2, 3.3333, 2, 6, 3.3333, 6, 2, 1, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  newValues = stag::sprsMatValues(testGraph.adjacency());

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

TEST(GraphTest, MalformedAdjacency1) {
  // Creating a graph with invalid sparse matrix data
  // should throw an error.

  // Create the incorrect data for the graph adjacency matrix.
  // Notice that there are not enough specified values for the length of the column
  // vector.
  std::vector<stag_int> rowStarts = {0, 2, 4, 7, 8};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6};

  // Attempt to create the graph object. Should throw an exception.
  EXPECT_THROW({stag::Graph testGraph(rowStarts, colIndices, values);},
               std::invalid_argument);
}

TEST(GraphTest, MalformedAdjacency2) {
  // Creating a graph with invalid sparse matrix data
  // should throw an error.

  // Create the incorrect data for the graph adjacency matrix.
  // Notice that the last value of the rowStarts vector is larger than the length
  // of the data vectors.
  std::vector<stag_int> rowStarts = {0, 2, 4, 7, 9};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};

  // Attempt to create the graph object. Should throw an exception.
  EXPECT_THROW({stag::Graph testGraph(rowStarts, colIndices, values);},
               std::invalid_argument);
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

  // Create the test graph object
  testGraph = createTestGraphSelfLoops();

  // Create the expected data for the graph laplacian matrix.
  rowStarts = {0, 3, 6, 10, 12};
  colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  values = {6.3333, -2, -3.3333, -2, 8, -6, -3.3333, -6, 12.3333, -1, -1, 2};

  // Check that the laplacian matrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(testGraph.laplacian());
  newIndices = stag::sprsMatInnerIndices(testGraph.laplacian());
  newValues = stag::sprsMatValues(testGraph.laplacian());

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

  // Create the test graph object
  testGraph = createTestGraphSelfLoops();

  // Create the expected data for the normalised graph laplacian matrix.
  rowStarts = {0, 3, 6, 10, 12};
  colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  values = {1 - 1/7.3333, -2/(sqrt(7.3333) * sqrt(8)), -3.3333/(sqrt(7.3333) * sqrt(14.3333)),
            -2/(sqrt(8) * sqrt(7.3333)), 1, -6/(sqrt(8) * sqrt(14.3333)),
            -3.3333/(sqrt(14.3333) * sqrt(7.3333)), -6/(sqrt(14.3333) * sqrt(8)), 1 - 2/14.3333, -1/(sqrt(14.3333) * sqrt(3)),
            -1/(sqrt(14.3333) * sqrt(3)), 1 - 1./3};

  // Check that the laplacian matrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(testGraph.normalised_laplacian());
  newIndices = stag::sprsMatInnerIndices(testGraph.normalised_laplacian());
  newValues = stag::sprsMatValues(testGraph.normalised_laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, SignlessLaplacianMatrix) {
  // Create the test graph object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the signless graph laplacian matrix.
  std::vector<stag_int> rowStarts = {0, 3, 6, 10, 12};
  std::vector<stag_int> colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> values = {5.3333, 2, 3.3333, 2, 8, 6, 3.3333, 6, 10.3333, 1, 1, 1};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.signless_laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.signless_laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.signless_laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);

  // Create the test graph object
  testGraph = createTestGraphSelfLoops();

  // Create the expected data for the signless graph Laplacian matrix.
  rowStarts = {0, 3, 6, 10, 12};
  colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  values = {8.3333, 2, 3.3333, 2, 8, 6, 3.3333, 6, 16.3333, 1, 1, 4};

  // Check that the laplacian matrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(testGraph.signless_laplacian());
  newIndices = stag::sprsMatInnerIndices(testGraph.signless_laplacian());
  newValues = stag::sprsMatValues(testGraph.signless_laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, NormalisedSignlessLaplacianMatrix) {
  // Create the test graph object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the normalised signless Laplacian matrix.
  std::vector<stag_int> rowStarts = {0, 3, 6, 10, 12};
  std::vector<stag_int> colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> values = {1, 2/(sqrt(5.3333) * sqrt(8)), 3.3333/(sqrt(5.3333) * sqrt(10.3333)),
                                2/(sqrt(8) * sqrt(5.3333)), 1, 6/(sqrt(8) * sqrt(10.3333)),
                                3.3333/(sqrt(10.3333) * sqrt(5.3333)), 6/(sqrt(10.3333) * sqrt(8)), 1, 1/sqrt(10.3333),
                                1/sqrt(10.3333), 1};

  // Check that the laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.normalised_signless_laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.normalised_signless_laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.normalised_signless_laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);

  // Create the test graph object
  testGraph = createTestGraphSelfLoops();

  // Create the expected data for the normalised signless Laplacian matrix.
  rowStarts = {0, 3, 6, 10, 12};
  colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  values = {1 + 1./7.3333, 2/(sqrt(7.3333) * sqrt(8)), 3.3333/(sqrt(7.3333) * sqrt(14.3333)),
            2/(sqrt(8) * sqrt(7.3333)), 1, 6/(sqrt(8) * sqrt(14.3333)),
            3.3333/(sqrt(14.3333) * sqrt(7.3333)), 6/(sqrt(14.3333) * sqrt(8)), 1 + 2./14.3333, 1/(sqrt(14.3333) * sqrt(3)),
            1/(sqrt(3) * sqrt(14.3333)), 1 + 1./3};

  // Check that the laplacian matrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(testGraph.normalised_signless_laplacian());
  newIndices = stag::sprsMatInnerIndices(testGraph.normalised_signless_laplacian());
  newValues = stag::sprsMatValues(testGraph.normalised_signless_laplacian());

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

  // Create the test graph object
  testGraph = createTestGraphSelfLoops();

  // Create the expected data for the graph degree matrix.
  rowStarts = {0, 1, 2, 3, 4};
  colIndices = {0, 1, 2, 3};
  values = {7.3333, 8, 14.3333, 3};

  // Check that the laplacian matrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(testGraph.degree_matrix());
  newIndices = stag::sprsMatInnerIndices(testGraph.degree_matrix());
  newValues = stag::sprsMatValues(testGraph.degree_matrix());

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

  // Create the test graph object
  testGraph = createTestGraphSelfLoops();

  // Create the expected data for the inverse degree matrix.
  rowStarts = {0, 1, 2, 3, 4};
  colIndices = {0, 1, 2, 3};
  values = {1./7.3333, 1./8, 1./14.3333, 1./3};

  // Check that the laplacian matrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(testGraph.inverse_degree_matrix());
  newIndices = stag::sprsMatInnerIndices(testGraph.inverse_degree_matrix());
  newValues = stag::sprsMatValues(testGraph.inverse_degree_matrix());

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

  // Create the test graph object
  testGraph = createTestGraphSelfLoops();

  // Create the expected data for the lazy random walk matrix.
  colStarts = {0, 3, 6, 10, 12};
  rowIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  values = {1./2 + 1./14.6666, 1./7.3333, 3.3333/14.6666,
            1./8, 1./2, 3./8,
            3.3333/28.6666, 3./14.3333, 1./2 + 1./14.3333, 1./28.6666,
            1./6, 1./2 + 1./6};

  // Check that the laplacian matrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(testGraph.lazy_random_walk_matrix());
  newIndices = stag::sprsMatInnerIndices(testGraph.lazy_random_walk_matrix());
  newValues = stag::sprsMatValues(testGraph.lazy_random_walk_matrix());

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

TEST(GraphTest, AverageDegree) {
  // The average degree of the cycle graph is always 2
  std::vector<stag_int> sizes = {3, 5, 10, 20, 100};
  for (stag_int n: sizes) {
    stag::Graph testGraph = stag::cycle_graph(n);
    EXPECT_EQ(testGraph.average_degree(), 2);
  }

  // The average degree of the star graph is 2 * (n-1) / n.
  for (stag_int n: sizes) {
    stag::Graph testGraph = stag::star_graph(n);
    EXPECT_EQ(testGraph.average_degree(), 2. * (n-1) / n);
  }

  stag::Graph testGraph = createTestGraphSelfLoops();
  EXPECT_EQ(testGraph.average_degree(), 32.6666/4);
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

TEST(GraphTest, StarGraph) {
  // Create a small star graph
  stag::Graph testGraph = stag::star_graph(5);

  // Define the expected adjacency matrix
  std::vector<stag_int> colStarts = {0, 4, 5, 6, 7, 8};
  std::vector<stag_int> rowIndices = {1, 2, 3, 4, 0, 0, 0, 0};
  std::vector<stag_int> values = {1, 1, 1, 1, 1, 1, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(colStarts, newStarts);
  EXPECT_EQ(rowIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, BarbellGraph) {
  // Create a small barbell graph
  stag::Graph testGraph = stag::barbell_graph(4);

  // Define the expected adjacency matrix
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

TEST(GraphTest, IdentityGraph) {
  // Create a small identity graph
  stag::Graph testGraph = stag::identity_graph(5);

  // Define the expected adjacency matrix
  std::vector<stag_int> colStarts = {0, 1, 2, 3, 4, 5};
  std::vector<stag_int> rowIndices = {0, 1, 2, 3, 4};
  std::vector<stag_int> values = {1, 1, 1, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(colStarts, newStarts);
  EXPECT_EQ(rowIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);

  // The Laplacian matrix should also be the identity
  newStarts = stag::sprsMatOuterStarts(testGraph.laplacian());
  newIndices = stag::sprsMatInnerIndices(testGraph.laplacian());
  newValues = stag::sprsMatValues(testGraph.laplacian());
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

  graph2 = createTestGraphSelfLoops();
  EXPECT_NE(graph1, graph2);
}

TEST(GraphTest, VertexExists) {
  stag_int n = 4;
  stag::Graph testGraph = stag::complete_graph(n);

  EXPECT_TRUE(testGraph.vertex_exists(0));
  EXPECT_TRUE(testGraph.vertex_exists(1));
  EXPECT_TRUE(testGraph.vertex_exists(2));
  EXPECT_TRUE(testGraph.vertex_exists(3));

  EXPECT_FALSE(testGraph.vertex_exists(4));
  EXPECT_FALSE(testGraph.vertex_exists(-1));
  EXPECT_FALSE(testGraph.vertex_exists(100));
}

TEST(GraphTest, ArgumentChecking) {
  stag_int n = -1;
  EXPECT_THROW(stag::cycle_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::complete_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::barbell_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::star_graph(n), std::invalid_argument);

  n = 0;
  EXPECT_THROW(stag::cycle_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::complete_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::barbell_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::star_graph(n), std::invalid_argument);

  n = 1;
  EXPECT_THROW(stag::cycle_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::complete_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::barbell_graph(n), std::invalid_argument);
  EXPECT_THROW(stag::star_graph(n), std::invalid_argument);
}

//------------------------------------------------------------------------------
// AdjacencyListLocalGraph tests
//------------------------------------------------------------------------------
TEST(GraphTest, ALLGVertexExists) {
  stag::AdjacencyListLocalGraph testGraph("test/data/test1.adjacencylist");
  EXPECT_TRUE(testGraph.vertex_exists(0));
  EXPECT_TRUE(testGraph.vertex_exists(1));
  EXPECT_TRUE(testGraph.vertex_exists(2));

  EXPECT_FALSE(testGraph.vertex_exists(3));
  EXPECT_FALSE(testGraph.vertex_exists(-1));
  EXPECT_FALSE(testGraph.vertex_exists(100));

  stag::AdjacencyListLocalGraph graph2("test/data/test7.adjacencylist");
  EXPECT_TRUE(graph2.vertex_exists(1));
  EXPECT_TRUE(graph2.vertex_exists(2));
  EXPECT_TRUE(graph2.vertex_exists(3));

  EXPECT_FALSE(graph2.vertex_exists(0));
  EXPECT_FALSE(graph2.vertex_exists(-1));
  EXPECT_FALSE(graph2.vertex_exists(100));
}

TEST(GraphTest, ALLGHugeGraph) {
  stag::AdjacencyListLocalGraph testGraph("test/data/test6.adjacencylist");
  EXPECT_NEAR(testGraph.degree(18950), 20.088781, 0.0001);
  EXPECT_EQ(testGraph.degree_unweighted(32107), 26);

  // Check the degree of the same node again. We should use the
  // cached adjacency matrix for this.
  EXPECT_NEAR(testGraph.degree(18950), 20.088781, 0.0001);
}

TEST(GraphTest, ALLGNeighboursOrder) {
  stag::AdjacencyListLocalGraph testGraph("test/data/test1.adjacencylist");
  stag_int node = 0;
  std::vector<stag_int> ns = testGraph.neighbors_unweighted(node);
  for (stag_int n : ns) {
    EXPECT_NE(n, node);
  }
}

TEST(GraphTest, Subgraph) {
  // Create a complete graph
  stag::Graph testGraph = stag::complete_graph(8);

  // Extract a subgraph (this will also be complete)
  std::vector<stag_int> vertices = {3, 5, 2, 1};
  stag::Graph subgraph = testGraph.subgraph(vertices);

  // Define the expected Laplacian matrix
  std::vector<stag_int> rowStarts = {0, 4, 8, 12, 16};
  std::vector<stag_int> colIndices = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  std::vector<double> values = {3, -1, -1, -1, -1, 3, -1, -1, -1, -1, 3, -1, -1, -1, -1, 3};

  // Check that the Laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(subgraph.laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(subgraph.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(subgraph.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);

  // Create a graph with self-loops
  testGraph = createTestGraphSelfLoops();

  // Extract a subgraph
  vertices = {1, 2, 3};
  subgraph = testGraph.subgraph(vertices);

  // Define the expected adjacency matrix
  rowStarts = {0, 1, 4, 6};
  colIndices = {1, 0, 1, 2, 1, 2};
  values = {6, 6, 2, 1, 1, 1};

  // Check that the adjacencymatrix has the form that we expect
  newStarts = stag::sprsMatOuterStarts(subgraph.adjacency());
  newIndices = stag::sprsMatInnerIndices(subgraph.adjacency());
  newValues = stag::sprsMatValues(subgraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, SubgraphCycle) {
  // Create a cycle graph
  stag::Graph testGraph = stag::cycle_graph(8);

  // Extract a subgraph - this will be a path on 3 vertices, and an isolated
  // vertex.
  std::vector<stag_int> vertices = {3, 4, 5, 7};
  stag::Graph subgraph = testGraph.subgraph(vertices);

  // Define the expected Laplacian matrix
  std::vector<stag_int> rowStarts = {0, 2, 5, 7, 8};
  std::vector<stag_int> colIndices = {0, 1, 0, 1, 2, 1, 2, 3};
  std::vector<double> values = {1, -1, -1, 2, -1, -1, 1, 0};

  // Check that the Laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(subgraph.laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(subgraph.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(subgraph.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, DisjointUnion) {
  // Construct two graphs
  stag::Graph g1 = stag::complete_graph(3);
  stag::Graph g2 = stag::cycle_graph(4);

  // Combine the two graphs into one
  stag::Graph g3 = g1.disjoint_union(g2);

  // Define the expected Laplacian matrix
  std::vector<stag_int> rowStarts = {0, 3, 6, 9, 12, 15, 18, 21};
  std::vector<stag_int> colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 4, 6, 3, 4, 5, 4, 5, 6, 3, 5, 6};
  std::vector<double> values = {2, -1, -1, -1, 2, -1, -1, -1, 2, 2, -1, -1, -1, 2, -1, -1, 2, -1, -1, -1, 2};

  // Check that the Laplacian matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(g3.laplacian());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(g3.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(g3.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);

  // Check that one of the graphs is not changed
  rowStarts = {0, 3, 6, 9};
  colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  values = {2, -1, -1, -1, 2, -1, -1, -1, 2};
  newStarts = stag::sprsMatOuterStarts(g1.laplacian());
  newIndices = stag::sprsMatInnerIndices(g1.laplacian());
  newValues = stag::sprsMatValues(g1.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, DisjoinUnionSelfLoops) {
  // Construct two graphs
  stag::Graph g1 = createTestGraph();
  stag::Graph g2 = createTestGraphSelfLoops();

  // One has self loops, the other doesn't
  EXPECT_EQ(g1.has_self_loops(), false);
  EXPECT_EQ(g2.has_self_loops(), true);

  // Combine the two graphs into one
  stag::Graph g3 = g1.disjoint_union(g2);

  // The new graph has self loops
  EXPECT_EQ(g3.has_self_loops(), true);

  // Define the expected adjacency matrix
  std::vector<stag_int> rowStarts = {0, 2, 4, 7, 8, 11, 13, 17, 19};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2, 4, 5, 6, 4, 6, 4, 5, 6, 7, 6, 7};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1, 1, 2, 3.3333, 2, 6, 3.3333, 6, 2, 1, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(g3.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(g3.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(g3.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, UnionComponents) {
  // Check that the connected components of a graph union are the original
  // graphs.

  // Create two random graphs.
  stag_int n = 100;
  stag_int k = 2;
  stag::Graph g1 = stag::sbm(n, k, 0.5, 0.5);
  stag::Graph g2 = stag::sbm(n, k, 0.5, 0.5);

  // Join the graphs
  stag::Graph union_graph = g1.disjoint_union(g2);

  // Find the first connected component of the union graph
  std::vector<stag_int> cc = stag::connected_component(&union_graph, 0);
  stag::Graph cc_graph = union_graph.subgraph(cc);

  // It's difficult to check equality of graphs, so let's just make sure that
  // the degree sequences are the same.
  EXPECT_EQ(cc_graph.number_of_vertices(), g1.number_of_vertices());
  std::vector<stag_int> cc_degrees;
  std::vector<stag_int> g1_degrees;
  for (auto v = 0; v < cc_graph.number_of_vertices(); v++) {
    cc_degrees.push_back(cc_graph.degree_unweighted(v));
    g1_degrees.push_back(g1.degree_unweighted(v));
  }
  std::sort(cc_degrees.begin(), cc_degrees.end());
  std::sort(g1_degrees.begin(), g1_degrees.end());
  EXPECT_TRUE(cc_degrees == g1_degrees);
}

TEST(GraphTest, GraphConnected) {
  // Create a random disconnected graph
  stag_int n = 100;
  stag::Graph g1 = stag::sbm(n, 2, 0.5, 0);
  EXPECT_FALSE(g1.is_connected());

  // The barbell graph is connected
  stag::Graph g2 = stag::barbell_graph(n);
  EXPECT_TRUE(g2.is_connected());

  // The identity graph is not connected
  stag::Graph g3 = stag::identity_graph(n);
  EXPECT_FALSE(g3.is_connected());
}

TEST(GraphTest, ALLGSelfLoopDegree) {
  stag::Graph testGraph = stag::identity_graph(5);
  std::string adjacencylist_filename = "output.adjacencylist";
  stag::save_adjacencylist(testGraph, adjacencylist_filename);
  stag::AdjacencyListLocalGraph allg{adjacencylist_filename};

  // The degree of a node in the identity graph is 2 - the self-loop counts
  // twice.
  for (stag_int i = 0; i < 5; i++) {
    EXPECT_EQ(allg.degree(i), 2);
    EXPECT_EQ(allg.degree_unweighted(i), 1);
  }
}

TEST(GraphTest, AddGraphs) {
  stag_int n = 4;
  stag::Graph testGraph = stag::complete_graph(n) + stag::cycle_graph(n);

  // Define the expected adjacency matrix
  std::vector<stag_int> rowStarts = {0, 3, 6, 9, 12};
  std::vector<stag_int> colIndices = {1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2};
  std::vector<double> values = {2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, AddGraphsBadSize) {
  stag::Graph g1 = stag::complete_graph(5);
  stag::Graph g2 = stag::cycle_graph(6);
  EXPECT_THROW(g1 + g2, std::invalid_argument);
}

TEST(GraphTest, ScalarMultiplication) {
  // Create a small star graph, multiplied by 3
  stag::Graph testGraph = 3 * stag::star_graph(5);

  // Define the expected adjacency matrix
  std::vector<stag_int> colStarts = {0, 4, 5, 6, 7, 8};
  std::vector<stag_int> rowIndices = {1, 2, 3, 4, 0, 0, 0, 0};
  std::vector<stag_int> values = {3, 3, 3, 3, 3, 3, 3, 3};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());
  EXPECT_EQ(colStarts, newStarts);
  EXPECT_EQ(rowIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);

  // Try post-multiplying by a double
  testGraph = stag::star_graph(5) * 3.0;
  newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(colStarts, newStarts);
  EXPECT_EQ(rowIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}
