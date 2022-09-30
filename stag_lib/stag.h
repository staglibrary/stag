#ifndef STAG_LIBRARY_H
#define STAG_LIBRARY_H

#include <Eigen/Core>

namespace stag {
  class Graph {
    public:
      Graph(Eigen::MatrixXd);
      void printAdjMat();

    private:
      Eigen::MatrixXd adjacencyMatrix;
  };
}
#endif //STAG_LIBRARY_H
