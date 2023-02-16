# STAG: Spectral Toolkit of Algorithms for Graphs

[![Build](https://github.com/pmacg/stag/actions/workflows/github-actions-test.yml/badge.svg?branch=main)](https://github.com/pmacg/stag/actions/workflows/github-actions-test.yml)
[![Docs](https://github.com/staglibrary/stag/actions/workflows/github-actions-docs.yml/badge.svg)](https://staglibrary.io/docs/cpp)

STAG is a C++ library providing efficient spectral algorithms for the analysis of massive graphs.
There is also an accompanying [STAG Python library](https://github.com/staglibrary/stagpy).

The library provides a simple and convenient interface for developing and applying graph algorithms.
For example, the following will read a graph from file and apply the spectral clustering algorithm.

```c++
#include <iostream>
#include "graphio.h"
#include "graph.h"
#include "cluster.h"

int main() {
  // Read a graph from an edgelist file.
  std::string filename = "mygraph.edgelist";
  stag::Graph myGraph = stag::load_edgelist(filename);

  // Apply the spectral clustering algorithm to the graph
  int k = 4;
  auto labels = stag::spectral_cluster(&myGraph, k);

  // Display the cluster membership of each vertex in the graph
  for (auto c : labels) {
    std::cout << c << ", ";
  }
  std::cout << std::endl;

  return 0;
}
```

STAG provides several methods for constructing graphs, including generating random graphs
from the stochastic block model.

```c++
#include "random.h"
...
// Construct a graph with 5 clusters using the 
// stochastic block model.
stag_int n = 1000;
stag_int k = 5;
stag::Graph myGraph = stag::sbm(n, k, 0.6, 0.1);
```

STAG also provides an implementation of the popular ACL local clustering algorithm.

```c++
#include "cluster.h"
...
double target_volume = 1000;
int seed_vertex = 0;
auto cluster = stag::local_cluster(&myGraph, seed_vertex, target_volume);
```

The complete library documentation is available on the [STAG library website](https://staglibrary.io/docs/cpp).

## Quick Start

The STAG library has the following dependencies.

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) (version >= 3.1)
- [Spectra](https://spectralib.org/) (version >= 1.0.1)

In order to use STAG in your C++ project, clone the git repository and include the `stag_lib`
directory in your project. Adding the following to your `CMakeLists.txt` will compile STAG
and make it available to the code in your project.

```cmake
set(CMAKE_CXX_STANDARD 20) # STAG requires at least C++20.
 
include_directories(stag_lib)
add_subdirectory(stag_lib stag_lib)
 
target_link_libraries(
        YOUR_PROJECT
        stag
)
```

## Contact

If you have any questions about using the STAG library, please contact a member of the STAG
development team.

- He Sun ([h.sun@ed.ac.uk](mailto:h.sun@ed.ac.uk))
- Peter Macgregor ([peter.macgregor@ed.ac.uk](mailto:peter.macgregor@ed.ac.uk))
