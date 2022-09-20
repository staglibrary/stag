#include <sgtl.h>
#include <Eigen/Core>

int main()
{
  Eigen::Matrix3d m1;
  m1 << 1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9;
  sgtl::Graph myGraph(m1);
  myGraph.printAdjMat();
  return 0;
}