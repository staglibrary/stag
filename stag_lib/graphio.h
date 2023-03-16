//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file graphio.h
 * \brief Methods for reading and writing graphs to disk.
 *
 * For an introduction to the file formats supported by STAG, see
 * [Graph File Formats](@ref file-formats).
 */

#ifndef STAG_TEST_GRAPHIO_H
#define STAG_TEST_GRAPHIO_H

#include <string>

#include <graph.h>

namespace stag {
  /**
   * Load a graph from an edgelist file.
   *
   * We define an edgelist file in the following way.
   *   - Any lines beginning with '#' or '//' are ignored
   *   - Any blank lines are ignored
   *   - All other lines are of one of the formats
   *       - `<u>, <v>, <weight>`
   *       - `<u>, <v>`
   *       - `<u> <v> <weight>`
   *       - `<u> <v>`
   *
   *     where `<u>` and `<v>` can be parsed as integers, and `<weight>` can be parsed
   *     as a double. If the weight is omitted, it is taken to be 1.
   *
   * Here is an example edgelist file.
   *
   *     # This line is ignored
   *     0, 1, 0.5
   *     1, 2, 1
   *     2, 0, 0.5
   *
   * @param filename the name of the edgelist file to be loaded
   * @return stag::Graph object
   * @throws std::runtime_error if the file doesn't exist or cannot be parsed as an edgelist
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
   * \cond
   * Parse a single content line of a STAG adjacency list file.
   */
  std::vector<stag::edge> parse_adjacencylist_content_line(std::string line);
  /**
   * \endcond
   */

  /**
   * \cond
   * Sort the edgelist file by the first vertex in each edge.
   *
   * All sorting happens on disk, keeping memory requirements low.
   */
  void sort_edgelist(std::string& filename);

  /**
   * Copy the edgelist file infile to outfile, while copying every edge to have
   * both directions.
   */
  void copy_edgelist_duplicate_edges(std::string& infile, std::string& outfile);
  /**
   * \endcond
   */

  /**
   * Load a graph from an adjacencylist file.
   *
   * The adjacency list file format is defined in the following way.
   *   - Any lines beginning with `#` or `//` are ignored
   *   - Any blank lines are ignored
   *   - All other lines have the format `<node_id>: <list of neighbours>`, where
   *     `<node_id>` is an integer and `<list of neighbours>` is either a
   *     space-separated list of integers or a space-separated list of
   *     `<id>:<weight>` where `<id>` gives the id of the neighbour and
   *     `<weight>` is the weight of the edge.
   *   - The `<node IDs>` of each line must be sorted in increasing order.
   *   - The graph should have no self-loops.
   *
   * Here is an example adjacencylist file.
   *
   *     # This line is ignored
   *     0: 1 2
   *     1: 0 2 3
   *     2: 0 1
   *     3: 1
   *
   * The following example includes weighted edges.
   *
   *     # This line is ignored
   *     0: 1:0.5 2:0.5
   *     1: 0:0.5 2:1
   *     2: 0:0.5 1:1
   *
   * Note that this method loads the entire graph into memory. For large graphs,
   * the stag::AdjacencyListLocalGraph object can be used to access the graph
   * in a 'local' way without reading the entire graph into memory.
   *
   * @param filename the name of the adjacency list file to be loaded
   * @return stag::Graph object
   * @throws std::runtime_error if the file doesn't exist or cannot be parsed as
   *         an adjacency list
   */
  stag::Graph load_adjacencylist(std::string& filename);

  /**
   * Save the given graph as an adjacencylist file.
   *
   * @param graph the graph object to be saved
   * @param filename the name of the file to save the adjacencylist data to
   */
  void save_adjacencylist(stag::Graph &graph, std::string& filename);

  /**
   * Convert an edgelist file to an adjacency list.
   *
   * @param edgelist_fname the name of the file containing the edgelist.
   * @param adjacencylist_fname the name of the file to write the adjacencylist.
   */
  void edgelist_to_adjacencylist(std::string& edgelist_fname,
                                 std::string& adjacencylist_fname);

  /**
   * Convert an adjacency list file to an edgelist.
   *
   * @param adjacencylist_fname the name of the file containing the adjacency list.
   * @param edgelist_fname the name of the file to write the edgelist.
   */
  void adjacencylist_to_edgelist(std::string& adjacencylist_fname,
                                 std::string& edgelist_fname);
}

#endif //STAG_TEST_GRAPHIO_H