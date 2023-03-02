//
// Created by macgr on 02/03/2023.
//

#include <chrono>
#include <thread>
#include <stdexcept>
#include <graph.h>
#include <utility.h>
#include <random.h>

int main() {
  stag_int n = 1000;
  stag_int k = 10;
  stag::Graph testGraph = stag::sbm(n * k, k, 0.3, 0.01);
  std::this_thread::sleep_for(std::chrono::milliseconds(5000));
  const SprsMat* lap_mat = testGraph.normalised_laplacian();
  std::this_thread::sleep_for(std::chrono::milliseconds(5000));
  return 0;
}
