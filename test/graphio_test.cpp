/**
* Tests for the graphio.h header. Includes methods for reading and writing
 * graphs to disk.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
*/
#include <iostream>
#include <gtest/gtest.h>
#include "graphio.h"
#include "utility.h"
#include "graph.h"
#include "random.h"

TEST(GraphioTest, FromEdgelistSimple) {
  std::string filename = "test/data/test1.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 6};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {1, 1, 1, 1, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, FromEdgelistWeights) {
  std::string filename = "test/data/test2.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 6};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {0.5, 0.5, 0.5, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, FromEdgelistExtraSpace) {
  std::string filename = "test/data/test3.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 6};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {1, 0.5, 1, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, FromEdgelistMissingWeights) {
  std::string filename = "test/data/test4.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 6};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {1, 0.5, 1, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, FromEdgelistTabs) {
  std::string filename = "test/data/test5.edgelist";
  stag::Graph testGraph = stag::load_edgelist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 6};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {0.5, 0.5, 0.5, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, EdgelistBadFilename) {
  std::string badFilename = "thisfiledoesntexist.edgelist";
  EXPECT_THROW({stag::Graph testGraph = stag::load_edgelist(badFilename);},
               std::runtime_error);

  // Now try writing to a file which cannot be written to
  stag::Graph testGraph = stag::complete_graph(5);
  badFilename = "////cannotwritehere.edgelist";
  EXPECT_THROW(stag::save_edgelist(testGraph, badFilename),
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

TEST(GraphioTest, LoadAdjacencylistSimple) {
  std::string filename = "test/data/test1.adjacencylist";
  stag::Graph testGraph = stag::load_adjacencylist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 6};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {1, 1, 1, 1, 1, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, LoadAdjacencylistWeights) {
  std::string filename = "test/data/test2.adjacencylist";
  stag::Graph testGraph = stag::load_adjacencylist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 6};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {0.5, 0.5, 0.5, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, LoadAdjacencylistExtraSpace) {
  std::string filename = "test/data/test3.adjacencylist";
  stag::Graph testGraph = stag::load_adjacencylist(filename);

  // Create the expected data for the graph adjacency matrix.
  std::vector<stag_int> rowStarts = {0, 2, 4, 6};
  std::vector<stag_int> colIndices = {1, 2, 0, 2, 0, 1};
  std::vector<double> values = {1, 0.5, 1, 1, 0.5, 1};

  // Check that the adjacency matrix has the form that we expect
  std::vector<stag_int> newStarts = stag::sprsMatOuterStarts(testGraph.adjacency());
  std::vector<stag_int> newIndices = stag::sprsMatInnerIndices(testGraph.adjacency());
  std::vector<double> newValues = stag::sprsMatValues(testGraph.adjacency());

  EXPECT_EQ(rowStarts, newStarts);
  EXPECT_EQ(colIndices, newIndices);
  EXPECT_EQ(values, newValues);
}

TEST(GraphioTest, AdjacencyListErrors) {
  std::string badFilename = "thisfiledoesntexist.adjacencylist";
  EXPECT_THROW({stag::Graph testGraph = stag::load_adjacencylist(badFilename);},
               std::runtime_error);

  std::string filename = "test/data/badgraph.adjacencylist";
  EXPECT_THROW({stag::Graph testGraph = stag::load_adjacencylist(filename);},
               std::domain_error);
}

TEST(GraphioTest, SaveAdjacencylist) {
  // Save an edgelist file
  std::string filename = "output.adjacencylist";
  stag::Graph testGraph = stag::cycle_graph(10);
  stag::save_adjacencylist(testGraph, filename);

  // Reading the adjacencylist file back in should result in the same graph
  stag::Graph newGraph = stag::load_adjacencylist(filename);
  EXPECT_EQ(testGraph, newGraph);

  // Save another adjacencylist file
  testGraph = stag::complete_graph(10);
  stag::save_adjacencylist(testGraph, filename);

  // Reading the adjacencylist file back in should result in the same graph
  newGraph = stag::load_adjacencylist(filename);
  EXPECT_EQ(testGraph, newGraph);

  // Finally, try a random graph
  testGraph = stag::erdos_renyi(100, 0.05);
  stag::save_adjacencylist(testGraph, filename);
  newGraph = stag::load_adjacencylist(filename);
  EXPECT_EQ(testGraph, newGraph);
}

TEST(GraphioTest, Conversions) {
  std::string edgelist_filename = "output.edgelist";
  std::string adjacencylist_filename = "output.adjacencylist";
  stag::Graph testGraph = stag::erdos_renyi(100, 0.05);
  stag::save_edgelist(testGraph, edgelist_filename);
  stag::edgelist_to_adjacencylist(edgelist_filename, adjacencylist_filename);
  stag::Graph newGraph = stag::load_adjacencylist(adjacencylist_filename);
  EXPECT_EQ(testGraph, newGraph);

  // Now convert back
  testGraph = stag::erdos_renyi(100, 0.05);
  stag::save_adjacencylist(testGraph, adjacencylist_filename);
  stag::adjacencylist_to_edgelist(adjacencylist_filename, edgelist_filename);
  newGraph = stag::load_edgelist(edgelist_filename);
  EXPECT_EQ(testGraph, newGraph);

  edgelist_filename = "test/data/test5.edgelist";
  stag::stream_edgelist_to_adjacencylist(edgelist_filename, adjacencylist_filename);
  testGraph = stag::load_edgelist(edgelist_filename);
  newGraph = stag::load_adjacencylist(adjacencylist_filename);
  EXPECT_EQ(testGraph, newGraph);

  edgelist_filename = "output.edgelist";
  testGraph = stag::erdos_renyi(3000, 0.01);
  stag::save_edgelist(testGraph, edgelist_filename);
  stag::stream_edgelist_to_adjacencylist(edgelist_filename,
                                         adjacencylist_filename);
  newGraph = stag::load_adjacencylist(adjacencylist_filename);
  EXPECT_EQ(testGraph, newGraph);
}
