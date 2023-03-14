/**
 * This is a command-line tool for converting between edgelist and adjacency
 * lists.
 */
#include <iostream>
#include <cerrno>
#include "graphio.h"

void print_usage() {
  std::cout << "Usage: stag_adj2edge [adjacencylist] [edgelist]" << std::endl;
  std::cout << std::endl;
  std::cout << "Convert a STAG adjacency list file into an edgelist file." << std::endl;
  std::cout << std::endl;
  std::cout << "Arguments:" << std::endl;
  std::cout << "  [adjacencylist]   the name of the adjacencylist file to be converted" << std::endl;
  std::cout << "  [edgelist]        the name of the new edgelist file to be written" << std::endl;
}


int main(int argc, char** args) {
  // This program takes two arguments: the edgelist file and the adjacencylist
  // file to write to.
  if (argc != 3) {
    print_usage();
    return EINVAL;
  }

  // Extract the command lien arguments.
  std::string adj_fname;
  std::string edge_fname;
  try {
    adj_fname = std::string(args[1]);
    edge_fname = std::string(args[2]);
  } catch (...) {
    print_usage();
    return EINVAL;
  }

  stag::adjacencylist_to_edgelist(adj_fname, edge_fname);

  return 0;
}
