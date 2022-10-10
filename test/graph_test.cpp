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
 * @return a stag::Graph object to be tested
 */
stag::Graph createTestGraph() {
  // Create the data for the graph adjacency matrix.
  std::vector<int> rowStarts = {0, 2, 4, 7, 8};
  std::vector<int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};

  // Create the stag Graph object
  return {rowStarts, colIndices, values};
}

TEST(GraphTest, Volume) {
  stag::Graph testGraph = createTestGraph();

  // The volume should be 24.6666
  EXPECT_EQ(testGraph.volume(), 24.6666);
}

TEST(GraphTest, AdjacencyMatrix) {
  // Create the data for the graph adjacency matrix.
  std::vector<int> rowStarts = {0, 2, 4, 7, 8};
  std::vector<int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};

  // Create the stag Graph object
  stag::Graph testGraph(rowStarts, colIndices, values);

  // Check that the adjacency matrix has the form that we expect
  std::vector<int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, LaplacianMatrix) {
  // Create the test graph object
  stag::Graph testGraph = createTestGraph();

  // Create the expected data for the graph laplacian matrix.
  std::vector<int> rowStarts = {0, 3, 6, 10, 12};
  std::vector<int> colIndices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
  std::vector<double> values = {5.3333, -2, -3.3333, -2, 8, -6, -3.3333, -6, 10.3333, -1, -1, 1};

  // Check that the laplacian matrix has the form that we expect
  std::vector<int> newStarts = stag::sprsMatOuterStarts(testGraph.laplacian());
  std::vector<int> newIndices = stag::sprsMatInnerIndices(testGraph.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, CycleGraphVolume) {
  // The volume of the cycle graph should be twice the number of vertices
  std::vector<int> sizes = {3, 3, 5, 10, 20, 100};
  for (int n: sizes) {
    stag::Graph testGraph = stag::cycle_graph(n);
    EXPECT_EQ(testGraph.volume(), 2 * n);
  }
}

TEST(GraphTest, CycleGraphLaplacian) {
  // Create a small cycle graph
  stag::Graph testGraph = stag::cycle_graph(4);

  // Define the expected laplacian matrix
  std::vector<int> rowStarts = {0, 3, 6, 9, 12};
  std::vector<int> colIndices = {0, 1, 3, 0, 1, 2, 1, 2, 3, 0, 2, 3};
  std::vector<double> values = {2, -1, -1, -1, 2, -1, -1, 2, -1, -1, -1, 2};

  // Check that the laplacian matrix has the form that we expect
  std::vector<int> newStarts = stag::sprsMatOuterStarts(testGraph.laplacian());
  std::vector<int> newIndices = stag::sprsMatInnerIndices(testGraph.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}

TEST(GraphTest, CompleteGraphVolume) {
  // The volume of the complete graph should be n(n-1)
  std::vector<int> sizes = {3, 3, 5, 10, 20, 100};
  for (int n: sizes) {
    stag::Graph testGraph = stag::complete_graph(n);
    EXPECT_EQ(testGraph.volume(), n * (n - 1));
  }
}

TEST(GraphTest, CompleteGraphLaplacian) {
  // Create a small complete graph
  stag::Graph testGraph = stag::complete_graph(4);

  // Define the expected laplacian matrix
  std::vector<int> rowStarts = {0, 4, 8, 12, 16};
  std::vector<int> colIndices = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  std::vector<double> values = {3, -1, -1, -1, -1, 3, -1, -1, -1, -1, 3, -1, -1, -1, -1, 3};

  // Check that the laplacian matrix has the form that we expect
  std::vector<int> newStarts = stag::sprsMatOuterStarts(testGraph.laplacian());
  std::vector<int> newIndices = stag::sprsMatInnerIndices(testGraph.laplacian());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.laplacian());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_FLOATS_NEARLY_EQ(values, newValues, 0.000001);
}
