#include <stag.h>
#include <Eigen/Core>

#include <iostream>

int main() {
  // Create the data for the graph adjacency matrix.
  int vertices = 4;
  std::vector<int> rows = {0, 0, 1, 1, 2, 2, 2, 3};
  std::vector<int> cols = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> vals = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};

  // Create the stag Graph object
  stag::Graph myGraph(vertices, rows, cols, vals);

  // Retrieve and print the Laplacian matrix
  SprsMat lap = myGraph.laplacian();
  std::cout << lap << std::endl;

  // Compute the graph's volume
  double vol = myGraph.volume();
  std::cout << "Volume: " << vol << std::endl;

  // Return
  return 0;
}
