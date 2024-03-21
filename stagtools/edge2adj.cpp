/**
 * This is a command-line tool for converting between edgelist and adjacency
 * lists.
 */
#include <iostream>
#include <cerrno>
#include "Graph/graphio.h"


void print_usage() {
  std::cout << "Usage: stag_edge2adj [edgelist] [adjacencylist]" << std::endl;
  std::cout << std::endl;
  std::cout << "Convert an edgelist file to a STAG adjacency list file." << std::endl;
  std::cout << std::endl;
  std::cout << "Arguments:" << std::endl;
  std::cout << "  [edgelist]        the name of the edgelist file to be converted" << std::endl;
  std::cout << "  [adjacencylist]   the name of the new adjacencylist file to be written" << std::endl;
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
    edge_fname = std::string(args[1]);
    adj_fname = std::string(args[2]);
  } catch (...) {
    print_usage();
    return EINVAL;
  }

  stag::edgelist_to_adjacencylist(edge_fname, adj_fname);

  return 0;
}
