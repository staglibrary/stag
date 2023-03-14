/**
 * This is a command-line tool for generating graphs from the symmetric
 * stochastic block model.
 */
#include <iostream>
#include <cerrno>
#include "graphio.h"
#include "random.h"

void print_usage() {
  std::cout << "Usage: stag_sbm [filename] [n] [k] [p] [q]" << std::endl;
  std::cout << std::endl;
  std::cout << "This program generates a graph from the symmetric stochastic block model." << std::endl;
  std::cout << std::endl;
  std::cout << "Arguments:" << std::endl;
  std::cout << "  [filename]        the name of the edgelist file to be written" << std::endl;
  std::cout << "  [n]               the number of vertices in the graph" << std::endl;
  std::cout << "  [k]               the number of clusters in the graph" << std::endl;
  std::cout << "  [p]               each edge inside a cluster is included with probability p" << std::endl;
  std::cout << "  [q]               each edge between clusters is included with probability q" << std::endl;
}

int main(int argc, char** args) {
  // This program takes 5 arguments:
  //   - the name of the edgelist file to generate
  //   - the total number of vertices in the graph
  //   - the number of clusters in the graph
  //   - the probability of each internal edge
  //   - the probability of each external edge
  if (argc != 6) {
    print_usage();
    return EINVAL;
  }

  // Extract the command line arguments.
  std::string edge_fname;
  stag_int n, k;
  double p, q;
  try {
    edge_fname = std::string(args[1]);
    n = std::stoi(args[2]);
    k = std::stoi(args[3]);
    p = std::stod(args[4]);
    q = std::stod(args[5]);
  } catch (...) {
    print_usage();
    return EINVAL;
  }

  // Prepare the general SBM arguments
  stag_int cluster_size = floor(((double) n) / ((double) k));
  std::vector<stag_int> cluster_sizes;
  DenseMat probabilities(k, k);
  for (auto i = 0; i < k; i++) {
    cluster_sizes.push_back(cluster_size);
    probabilities(i, i) = p;

    for (auto j = i + 1; j < k; j++) {
      probabilities(i, j) = q;
      probabilities(j, i) = q;
    }
  }

  stag::general_sbm_edgelist(edge_fname, cluster_sizes, probabilities);

  return 0;
}
