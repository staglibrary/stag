Getting Started
===============
The STAG library is an easy-to-use C++ library providing several spectral
algorithms for graphs.
STAG is built on the Eigen library for representing sparse matrices and
so it may be useful to also refer to the `Eigen documentation <https://eigen.tuxfamily.org/dox/>`_.

Installation
------------
To include the STAG library in your C++ project, you should first install the following
dependencies.

- Eigen (version >= 3.1)
- Spectra (version >= 1.0.1)

You should refer to their documentation for installation instructions,
although the following should work on a standard linux system.

.. code-block::

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

Then, download STAG from github and copy the ``stag_lib`` directory into you project.
Adding the following to your ``CMakeLists.txt`` will compile STAG and make it available
to the code in your project.

.. code-block::

   set(CMAKE_CXX_STANDARD 20) # STAG requires at least C++20.

   include_directories(stag_lib)
   add_subdirectory(stag_lib stag_lib)

   target_link_libraries(
           YOUR_PROJECT
           stag
   )

Then, you can construct a graph in your project with the following code.

.. code-block::

   #include <iostream>
   #include "graph.h"

   int main() {
     stag::Graph myGraph = stag::cycle_graph(10);

     // Display the adjacency matrix of the graph
     std::cout << *myGraph.adjacency() << std::endl;

     return 0;
   }


