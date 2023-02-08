Clustering Algorithms
=====================

The STAG library provides two graph clustering algorithms:
    - :cpp:func:`stag::spectral_cluster`
    - :cpp:func:`stag::local_cluster`

.. doxygenfunction:: stag::spectral_cluster

.. doxygenfunction:: stag::local_cluster

.. doxygenfunction:: local_cluster_acl(stag::LocalGraph *graph, stag_int seed_vertex, double locality, double error)

.. doxygenfunction:: local_cluster_acl(stag::LocalGraph *graph, stag_int seed_vertex, double locality)

.. doxygenfunction:: stag::approximate_pagerank

.. doxygenfunction:: stag::sweep_set_conductance


