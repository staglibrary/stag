/**
 * This is a command-line tool for generating graphs from the symmetric
 * stochastic block model.
 */
#include <iostream>
#include <cerrno>
#include <stag/graphio.h>
#include <stag/random.h>

int main(int argc, char** args) {
  // This program takes 5 arguments:
  //   - the name of the edgelist file to generate
  //   - the total number of vertices in the graph
  //   - the number of clusters in the graph
  //   - the probability of each internal edge
  //   - the probability of each external edge
  if (argc != 6) {
    std::cout << "This program expects 5 command line arguments:" << std::endl;
    std::cout << "  - the name of the edgelist file to write the graph to" << std::endl;
    std::cout << "  - the number of vertices in the graph" << std::endl;
    std::cout << "  - the number of clusters in the graph" << std::endl;
    std::cout << "  - the probability of each internal edge" << std::endl;
    std::cout << "  - the probability of each external edge" << std::endl;
    std::cout << "For example:" << std::endl;
    std::cout << "    stag_sbm mygraph.edgelist 1000 5 0.3 0.01" << std::endl;
    return EINVAL;
  }

  // Extract the command line arguments.
  std::string edge_fname(args[1]);
  stag_int n = std::stoi(args[2]);
  stag_int k = std::stoi(args[3]);
  double p = std::stod(args[4]);
  double q = std::stod(args[5]);

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
