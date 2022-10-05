#include <stag.h>
#include <utility.h>
#include <Eigen/Core>

#include <iostream>

int main() {
  // Create the data for the graph adjacency matrix.
  std::vector<int> rowIndices = {0, 2, 4, 7, 8};
  std::vector<int> colIndices = {1, 2, 0, 2, 0, 1, 3, 2};
  std::vector<double> values = {2, 3.3333, 2, 6, 3.3333, 6, 1, 1};

  // Create the stag Graph object
  stag::Graph myGraph(rowIndices, colIndices, values);

  // Print the adjacency matrix
  std::cout << "Adjacency Matrix" << std::endl;
  std::cout << myGraph.adjacency() << std::endl;

  // View the outer starts of the adjacency matrix
  std::cout << "Outer starts" << std::endl;
  SprsMat adj = myGraph.adjacency();
  std::vector<int> starts = stag::sprsMatOuterStarts(adj);
  for (int i: starts){
    std::cout << i << ' ';
  }
  std::cout << std::endl << std::endl;

  // Retrieve and print the Laplacian matrix
  SprsMat lap = myGraph.laplacian();
  std::cout << "Laplacian Matrix" << std::endl;
  std::cout << lap << std::endl;

  // Compute the graph's volume
  double vol = myGraph.volume();
  std::cout << "Volume: " << vol << std::endl;

  // Return
  return 0;
}
