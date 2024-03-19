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
StagInt n = 1000;
StagInt k = 5;
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

You should install these dependencies before installing STAG.

To install STAG, clone the git repository or download the source code for the release you
would like to use. Then you can install the library with the following commands from the
root directory of the repository.

```bash
mkdir build_dir
cd build_dir
cmake ..
make install
```

The final make command may require administrator privileges.

Once STAG is installed, you can use it in your C++ project by including the following in
your `CMakeLists.txt` file.

```cmake
set(CMAKE_CXX_STANDARD 20) # STAG requires at least C++20.
 
# Find and include the STAG library
find_package(stag REQUIRED)
message(STATUS "Found STAG!")
include_directories(${STAG_INCLUDE_DIRS})
 
target_link_libraries(
        YOUR_PROJECT
        stag
)
```

You may find it helpful to refer to the [example STAG project](https://github.com/staglibrary/example-stag-project).

## Citing STAG

If you find the STAG library useful in your research, please cite our [technical report](https://arxiv.org/abs/2304.03170).

```
@article{stag,
    author       = {Peter Macgregor and He Sun},
    title        = {Spectral Toolkit of Algorithms for Graphs: Technical Report {(1)}},
    journal      = {CoRR},
    volume       = {abs/2304.03170},
    year         = {2023},
    doi          = {10.48550/arXiv.2304.03170},
    eprinttype    = {arXiv},
    eprint       = {2304.03170},
}
```

## Licensing
The source code for the STAG library is available freely for use and modification.
The code is licenced under the GPL license, which requires that derivative works
are also released with the same license. 

## Acknowledgements

STAG includes open-source code developed by others. We are grateful to the following people for their contributions of
open-source code.

- Michael C. Hughes for the C++ implementation of the k-means algorithm available [here](https://github.com/michaelchughes/KMeansRex/).
- Alexandr Andoni and Piotr Indyk for their theoretical work and open source implementation of the Euclidean LSH algorithm, available [here](https://www.mit.edu/~andoni/LSH/).

## Contact

If you have any questions about using the STAG library, please contact a member of the STAG
development team.

- He Sun ([h.sun@ed.ac.uk](mailto:h.sun@ed.ac.uk))
- Peter Macgregor ([peter.macgregor@ed.ac.uk](mailto:peter.macgregor@ed.ac.uk))
