Constructing Graphs
===================

There are two provided graph classes:

- :cpp:class:`stag::LocalGraph` is an abstract class providing methods for local
  access to the graph.
- :cpp:class:`stag::Graph` is a general-purpose graph object, defined by a sparse
  Eigen matrix. This class implements and extends the :cpp:class:`stag::LocalGraph` class.

.. doxygenclass :: stag::Graph
   :members:

.. doxygenclass :: stag::LocalGraph
   :members:

