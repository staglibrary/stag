#ifndef SGTL_LIBRARY_H
#define SGTL_LIBRARY_H

#include <Eigen/Core>

namespace sgtl {
  class Graph {
    public:
      Graph(Eigen::MatrixXd);
      void printAdjMat();

    private:
      Eigen::MatrixXd adjacencyMatrix;
  };
}
#endif //SGTL_LIBRARY_H
