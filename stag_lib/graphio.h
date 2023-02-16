//
// This file is provided as part of the STAG library and released under the MIT
// license.
//

/**
 * @file graphio.h
 * \brief Methods for reading and writing graphs to disk.
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
}

#endif //STAG_TEST_GRAPHIO_H