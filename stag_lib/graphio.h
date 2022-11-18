/**
 * Definitions related to reading and writing graphs from disk.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#ifndef STAG_TEST_GRAPHIO_H
#define STAG_TEST_GRAPHIO_H

#include <string>

#include "neo4j-client.h"

#include <graph.h>

namespace stag {
  /**
   * Load a graph from an edgelist file.
   *
   * We define an edgelist file in the following way.
   *   - Any lines beginning with '#' or '//' are ignored
   *   - Any blank lines are ignored
   *   - All other lines are of some format:
   *       <u>, <v>, <weight>
   *       <u>, <v>
   *       <u> <v> <weight>
   *       <u> <v>
   *     where <u> and <v> can be parsed as integers, and <weight> can be parsed
   *     as a double. If the weight is ommitted, it is taken to be 1. The
   *     data lines of the edgelist file should all have the same format.
   *
   * @return stag::Graph object
   * @throws std::runtime_error if the file cannot be parsed as an edgelist
   */
  stag::Graph load_edgelist(std::string &filename);

  /**
   * Save the given graph as an edgelist file.
   *
   * @param graph the graph object to be saved
   * @param filename the name of the file to save the edgelist data to
   */
  void save_edgelist(stag::Graph &graph, std::string &filename);

  /**
   * A local graph backed by a neo4j database.
   */
  class Neo4jGraph : public LocalGraph {
    public:
      Neo4jGraph(const std::string& uri,
                 const std::string& username,
                 const std::string& password);

      double degree(stag_int v) override;
      stag_int degree_unweighted(stag_int v) override;
      std::vector<edge> neighbors(stag_int v) override;
      std::vector<stag_int> neighbors_unweighted(stag_int v) override;
      ~Neo4jGraph() override;

    private:
      // The connection to the underlying neo4j database.
      neo4j_connection_t* connection;
  };
}

#endif //STAG_TEST_GRAPHIO_H
