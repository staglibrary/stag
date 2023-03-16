Graph File Formats {#file-formats}
=======================

[TOC]

The STAG library supports two simple file formats for storing graphs on disk:
EdgeList and AdjacencyList.

EdgeList File Format
--------------------
In an EdgeList file, each line in the file corresponds to one edge in the graph.
A line consists of two node IDs which should be integers and optionally an
edge weight which can be a real number.
The IDs and weight should be separated with spaces.

For example, here is a simple EdgeList file.

~~~~~~~~~~~~~
# This is a comment
0 1
1 2
2 0
~~~~~~~~~~~~~

This file defines a graph on three nodes (`0`, `1`, and `2`) with three edges:
`(0, 1)`, `(1, 2)`, and `(0, 2)`.
Here is another example, which specifies edge weights.

~~~~~~~~~~~~~
# The third number on each row is the weight of the edge.
0 1 0.5
1 2 1
2 0 3
~~~~~~~~~~~~~

AdjacencyList File Format
-------------------------
In an AdjacencyList file, each line in the file corresponds to one node in
the graph.
A line consists of the node ID, followed by a list of adjacent nodes.
The node IDs at the beginning of each line must be sorted in increasing order.
For example, here is a simple AdjacencyList file.

~~~~~~
# This is a comment
0: 1 2
1: 0 3 2
2: 0 1
3: 1
~~~~~~

In this example, the node with ID `1` has edges to nodes `0`, `2`, and `3`.
All edges are taken to have weight 1.
The following example shows how to specify the weight of each edge.

~~~~~~
# The number after each colon is the weight of the edge.
0: 1:0.5 2:0.5
1: 0:0.5 2:1
2: 0:0.5 1:1
~~~~~~

Here, the node with ID `1` has an edge of weight `0.5` to node `0` and an edge
of weight `1` to node `2`.


Working with Files
------------------

### Reading and Writing
The graphio.h module provides several methods for reading, writing and converting
between EdgeList and AdjacencyList files.
For example, the following code will read an AdjacencyList file and create a
stag::Graph object.

~~~~~~{.cpp}
    #include <iostream>
    #include <stag/graph.h>
    #include <stag/graphio.h>

    int main() {
        // Load the graph from the file
        std::string filename = "mygraph.adjacencylist";
        stag::Graph myGraph = stag::load_adjacencylist(filename);
        
        // Display the adjacency matrix of the graph
        std::cout << *myGraph.adjacency() << std::endl;
  
        return 0;
    }
~~~~~~

The following example will generate a random graph and save it as an edgelist.

~~~~~~{.cpp}
    #include <iostream>
    #include <stag/graph.h>
    #include <stag/graphio.h>
    #include <stag/random.h>

    int main() {
        // Create a random graph
        stag::Graph myGraph = stag::erdos_renyi(100, 0.1);
        
        // Save the graph
        std::string filename = "mygraph.edgelist";
        stag::save_edgelist(myGraph, filename);
  
        return 0;
    }
~~~~~~

The stag::edgelist_to_adjacencylist and stag::adjacencylist_to_edgelist methods
allow for conversion between the two file formats.

### Local Access
When working with huge graphs, the stag::AdjacencyListLocalGraph object
provides a way to query a graph locally without loading the entire graph into
memory.
For example, the following code will query the neighbours of the node with ID
`100`.

~~~~~~{.cpp}
    #include <iostream>
    #include <stag/graph.h>

    int main() {
        // Create the local graph object
        std::string filename = "mygraph.adjacencylist";
        std::AdjacencyListLocalGraph myGraph(filename);
        
        // Display the neighbours of node 100
        auto neighbors = myGraph.neighbors_unweighted(100);
        for (auto n : neighbors) {
            std::cout << n << ", ";
        }
        std::cout << std::endl;
  
        return 0;
    }
~~~~~~

One of the most useful applications of the stag::AdjacencyListLocalGraph object
is for local clustering.
The following example finds a cluster containing node `100`.

~~~~~~{.cpp}
    #include <iostream>
    #include <stag/graph.h>
    #include <stag/cluster.h>

    int main() {
        // Create the local graph object
        std::string filename = "mygraph.adjacencylist";
        std::AdjacencyListLocalGraph myGraph(filename);
        
        // Find a cluster
        auto node_id = 100;
        auto target_volume = 1000;
        auto cluster = stag::local_cluster(&myGraph, node_id, target_volume);
       
        // Display the cluster 
        for (auto v : cluster) {
            std::cout << v << ", ";
        }
        std::cout << std::endl;
  
        return 0;
    }
~~~~~~

### STAG Tools
Finally, the [STAG tools](@ref stag-tools) provide command line tools
for converting between graph file formats.
