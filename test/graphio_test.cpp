/**
* Tests for the graphio.h header. Includes methods for reading and writing
 * graphs to disk.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
*/
#include <iostream>
#include <gtest/gtest.h>
#include "stag.h"
#include "graphio.h"
#include "utility.h"
#include "graph.h"

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

TEST(GraphioTest, SaveEdgelist) {
  // Save an edgelist file
  std::string filename = "output.edgelist";
  stag::Graph testGraph = stag::cycle_graph(10);
  stag::save_edgelist(testGraph, filename);

  // Reading the edgelist file back in should result in the same graph
  stag::Graph newGraph = stag::load_edgelist(filename);
  EXPECT_EQ(testGraph, newGraph);

  // Save another edgelist file
  testGraph = stag::complete_graph(10);
  stag::save_edgelist(testGraph, filename);

  // Reading the edgelist file back in should result in the same graph
  newGraph = stag::load_edgelist(filename);
  EXPECT_EQ(testGraph, newGraph);
}
