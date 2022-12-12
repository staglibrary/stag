Getting Started
===============
The STAG library is an easy-to-use C++ library providing several spectral
algorithms for graphs.
STAG is built on the Eigen library for representing sparse matrices and
so it may be useful to also refer to the `Eigen documentation <https://eigen.tuxfamily.org/dox/>`_.

Installation
------------
To include the STAG library in your C++ project:

- download STAG from github
- copy the ``stag_lib`` and ``eigen-3.3.9`` directories into your project
- add the following to your ``CMakeLists.txt``.

.. code-block::

   set(CMAKE_CXX_STANDARD 20) # STAG requires at least C++20.

   include_directories(eigen-3.3.9) # Eigen - header-only
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


