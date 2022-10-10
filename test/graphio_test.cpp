/**
* Tests for the graphio.h header. Includes methods for reading and writing
 * graphs to disk.
 *
 * Copyright 2022 Peter Macgregor
*/
#include <gtest/gtest.h>
#include <stag.h>
#include "graphio.h"
#include "utility.h"

TEST(GraphioTest, FromEdgelistSimple) {
  std::string filename = "test/data/test1.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<int> rowStarts = {0, 2, 4, 6};
  std::vector<int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {1, 1, 1, 1, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, FromEdgelistWeights) {
  std::string filename = "test/data/test2.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<int> rowStarts = {0, 2, 4, 6};
  std::vector<int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {0.5, 0.5, 0.5, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, FromEdgelistExtraSpace) {
  std::string filename = "test/data/test3.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<int> rowStarts = {0, 2, 4, 6};
  std::vector<int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {1, 0.5, 1, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, FromEdgelistMissingWeights) {
  std::string filename = "test/data/test4.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<int> rowStarts = {0, 2, 4, 6};
  std::vector<int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {1, 0.5, 1, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, EdgelistBadFilename) {
  std::string badFilename = "thisfiledoesntexist.edgelist";
  EXPECT_THROW({stag::Graph testGraph = stag::load_edgelist(badFilename);},
               std::runtime_error);
}
