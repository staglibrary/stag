#include "stag.h"
#include <iostream>
#include <utility>

stag::Graph::Graph(Eigen::MatrixXd adjMat) {
   adjacencyMatrix = std::move(adjMat);
};

void stag::Graph::printAdjMat() {
  std::cout << adjacencyMatrix << std::endl;
}
