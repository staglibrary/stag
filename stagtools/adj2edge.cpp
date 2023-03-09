/**
 * This is a command-line tool for converting between edgelist and adjacency
 * lists.
 */
#include <iostream>
#include <cerrno>
#include "graphio.h"

int main(int argc, char** args) {
  // This program takes two arguments: the edgelist file and the adjacencylist
  // file to write to.
  if (argc != 3) {
    std::cout << "This program expects 2 command line arguments." << std::endl;
    return EINVAL;
  }

  // Extract the command lien arguments.
  std::string adj_fname(args[1]);
  std::string edge_fname(args[2]);

  stag::adjacencylist_to_edgelist(adj_fname, edge_fname);

  return 0;
}
