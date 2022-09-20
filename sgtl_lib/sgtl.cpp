#include "sgtl.h"
#include <iostream>
#include <utility>

sgtl::Graph::Graph(Eigen::MatrixXd adjMat) {
   adjacencyMatrix = std::move(adjMat);
};

void sgtl::Graph::printAdjMat() {
  std::cout << adjacencyMatrix << std::endl;
}
