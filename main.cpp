#include <stag.h>
#include <Eigen/Core>

#include <iostream>

int main()
{
  // Create a simple wee graph to experiment with
  Eigen::Matrix3d m1;
  m1 << 0, 2, 3.33333, 2, 0, 6, 3.33333, 6, 0;

  // Create the stag Graph object
  stag::Graph myGraph(m1);

  // Retrieve and print the Laplacian matrix
  Eigen::Matrix3d lap = myGraph.laplacian();
  std::cout << lap << std::endl;

  // Return
  return 0;
}