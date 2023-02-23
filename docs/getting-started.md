Getting Started {#getting-started}
===============

[TOC]

The STAG library is an easy-to-use C++ library providing several spectral algorithms for graphs.
STAG is built on the Eigen library for representing sparse matrices and so it may be useful to also refer to the [Eigen documentation](https://eigen.tuxfamily.org/index.php?title=Main_Page).

Installation
------------

### Dependencies
To include the STAG library in your C++ project, you should first install the following
dependencies.

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) (version >= 3.1)
- [Spectra](https://spectralib.org/) (version >= 1.0.1)


You should refer to their documentation for installation instruction,
although the following should work on a standard Linux system.

~~~~~~~~~~~~~~~~~~~~
# Create a directory to work in
mkdir libraries
cd libraries

# Install Eigen
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar xzvf eigen-3.4.0.tar.gz
cd eigen-3.4.0
mkdir build_dir
cd build_dir
cmake ..
sudo make install
cd ../..

# Install Spectra
wget https://github.com/yixuan/spectra/archive/v1.0.1.tar.gz
tar xzvf v1.0.1.tar.gz
cd spectra-1.0.1
mkdir build_dir
cd build_dir
cmake ..
sudo make install
cd ../..
~~~~~~~~~~~~~~~~~~~~

### Installing STAG

You should download the latest version of the STAG library. For example:

```bash
wget https://github.com/staglibrary/stag/archive/refs/tags/v1.1.0.tar.gz
tar xzvf v1.1.0.tar.gz
cd stag-1.1.0
```

Then, you can build and install STAG with cmake.

```bash
mkdir build_dir
cd build_dir
cmake ..
sudo make install
```

Adding the following to your ``CMakeLists.txt`` will make the STAG library available to your project.

~~~~~~~~~~~~~~~~~~~~{.cmake}
set(CMAKE_CXX_STANDARD 20) # STAG requires at least C++20.
 
# Find and include the STAG library
find_package(stag REQUIRED)
message(STATUS "Found STAG!")
include_directories(${STAG_INCLUDE_DIRS})
 
target_link_libraries(
        YOUR_PROJECT
        stag
)
~~~~~~~~~~~~~~~~~~~~

You may find it helpful to refer to the [example STAG project](https://github.com/staglibrary/example-stag-project).

Constructing Graphs
-------------------
There are four methods for creating a stag::Graph object:

1. Constructing a standard named graph.
2. Constructing a graph with a sparse adjacency matrix.
3. Loading a graph from an edgelist file.
4. Generating a random graph.

### Named Graphs

The STAG library can easily construct several standard graphs such as the
cycle graph, the complete graph, and the barbell graph.
See graph.h for the full list and their documentation.

You can construct a named graph with the following code.

~~~~~~~~~~~~~~~~~~~~{.cpp}
    #include <iostream>
    #include <stag/graph.h>

    int main() {
      stag::Graph myGraph = stag::cycle_graph(10);

      // Display the adjacency matrix of the graph
      std::cout << *myGraph.adjacency() << std::endl;

      return 0;
    }
~~~~~~~~~~~~~~~~~~~~

### Adjacency Matrices

You can construct a stag::Graph object with any sparse adjacency matrix using
the default constructor.

~~~~~~~~~~~~~~~~~~~~{.cpp}
    #include <iostream>
    #include <stag/graph.h>

    int main() {
        // Construct a sparse matrix representing the
        // triangle graph adjacency matrix.
        stag_int n = 3;
        SprsMat adj(n, n);
        adj.coeffRef(0, 1) = 1;
        adj.coeffRef(0, 2) = 1;
        adj.coeffRef(1, 0) = 1;
        adj.coeffRef(1, 2) = 1;
        adj.coeffRef(2, 0) = 1;
        adj.coeffRef(2, 1) = 1;

        // Create a new STAG graph
        stag::Graph myGraph(adj);

        // Display the adjacency matrix of the graph
        std::cout << *myGraph.adjacency() << std::endl;
  
        return 0;
    }
~~~~~~~~~~~~~~~~~~~~

There are many other ways to construct the sparse adjacency matrix.
Please refer to the [Eigen Documentation](https://eigen.tuxfamily.org/dox/group__TutorialSparse.html).

### Loading Graphs from Disk

The STAG library can read and write graphs to a simple edgelist file format.
In an edgelist file, each row corresponds to one edge in the graph and lists 
the indices of the edges endpoints, and optionally the weight of the edge.

For example, here is a simple edgelist file.

~~~~~~~~~~~~~
# This line is ignored
0, 1, 0.5
1, 2, 1
2, 0, 0.5
~~~~~~~~~~~~~

This represents a graph on three vertices with three edges.
For example, the first row corresponds to an edge between vertices 0 and 1 with weight
0.5.

To load a STAG graph from an edgelist file, you can use the following code.

~~~~~~~~~~~~~~~{.cpp}
    #include <iostream>
    #include <stag/graph.h>
    #include <stag/graphio.h>

    int main() {
        // Load the graph from the file
        std::string filename = "mygraph.edgelist";
        stag::Graph myGraph = stag::load_edgelist(filename);
        
        // Display the adjacency matrix of the graph
        std::cout << *myGraph.adjacency() << std::endl;
  
        return 0;
    }
~~~~~~~~~~~~~~~

The STAG library can also save a graph to an edgelist file.

~~~~~~~~~~~~~~~{.cpp}
    #include <iostream>
    #include <stag/graph.h>
    #include <stag/graphio.h>

    int main() {
        // Create a cycle graph
        stag::Graph myGraph = stag::cycle_graph(10);
        
        // Save the graph as an edgelist file
        std::string filename = "mygraph.edgelist";
        stag::save_edgelist(myGraph, filename);
  
        return 0;
    }
~~~~~~~~~~~~~~~

### Random Graphs

The random.h module provides several methods for generating
Erdos-Renyi graphs and random graphs drawn from the stochastic
block model.

~~~~~~~~~~~~~~~~{.cpp}
#include <iostream>
#include <stag/graph.h>
#include <stag/random.h>

int main() {
    // Create an Erdos-Renyi random graph with
    // 100 vertices and edge pobability 0.1.
    double p = 0.1;
    stag_int n = 100;
    stag::Graph myGraph = stag::erdos_renyi(n, p);
    
    // Display the adjacency matrix of the graph
    std::cout << *myGraph.adjacency() << std::endl;
    
    // Create a graph drawn from the stochastic block model
    // with two clusters of size 100, internal edge probability 0.1
    // and external edge probability 0.01
    n = 200;
    stag_int k = 2
    double q = 0.01;
    stag::Graph sbmGraph = stag::sbm(n, k, p, q);
    
    return 0;
}
~~~~~~~~~~~~~~~~

Clustering
----------

The cluster.h module provides two clustering methods which are applicable for a wide variety of applications.

### Spectral clustering

Spectral clustering is a popular and simple graph clustering algorithm.
STAG provides a simple spectral clustering interface which is demonstrated in the following example.

~~~~~~~~~~~~~~{.cpp}
#include <iostream>
#include <stag/graph.h>
#include <stag/random.h>
#include <stag/cluster.h>

int main() {
  // Construct a graph with 5 clusters using the 
  // stochastic block model.
  stag_int n = 1000;
  stag_int k = 5;
  stag::Graph testGraph = stag::sbm(n, k, 0.6, 0.1);

  // Find the clusters
  auto labels = stag::spectral_cluster(&testGraph, k);
  
  // Display the cluster membership of each vertex in the graph
  for (auto c : labels) {
    std::cout << c << ", ";
  }
  std::cout << std::endl;
  
  return 0;
}
~~~~~~~~~~~~~~

The stag::spectral_cluster method returns a vector giving the cluster label
of each vertex in the graph.

### Local clustering

Given a vertex of the underlying graph as input, a local clustering algorithm finds a
cluster around the input vertex. 
The algorithm runs in time proportional to the size of the returned 
cluster and independent of the size of the entire graph.

Local clustering is provided by the stag::local_cluster method which can be used as
in the following example.

~~~~~~~~~~~~~~{.cpp}
#include <iostream>
#include <stag/graph.h>
#include <stag/random.h>
#include <stag/cluster.h>

int main() {
  // Construct a graph with 4 clusters using the 
  // stochastic block model.
  stag::Graph myGraph = stag::sbm(2000, 4, 0.9, 0.01);

  // Find a cluster containing the vertex 0
  std::vector<stag_int> cluster = stag::local_cluster(
    &myGraph, 0, myGraph.total_volume() / 4);
 
  // Display the vertices in the found cluster 
  for (auto v : cluster) {
    std::cout << v << ", ";
  }
  std::cout << std::endl;
  
  return 0;
}
~~~~~~~~~~~~~~

Notice that the local clustering method requires an estimate of the volume of the target cluster.
This parameter does *not* impose a hard constraint on the algorithm and so an approximate volume is enough.

The stag::local_cluster method returns a vector containing the ID of each vertex in the local cluster.

