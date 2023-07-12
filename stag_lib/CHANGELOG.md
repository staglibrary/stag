# Changelog
All notable changes to the library are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2023-07-12
### Added
- [Issue #92](https://github.com/staglibrary/stag/issues/92): Method for calculating the symmetric difference
- [Issue #94](https://github.com/staglibrary/stag/issues/94): Compute the connected components of a graph
- [Issue #94](https://github.com/staglibrary/stag/issues/94): Add a subgraph method to stag Graph class
- [Issue #93](https://github.com/staglibrary/stag/issues/93): Add methods to compute mutual information between clusters
- [Issue #98](https://github.com/staglibrary/stag/issues/98): Add method to construct the graph union
- [Issue #22](https://github.com/staglibrary/stag/issues/22): Add Cheeger cut method
- [Issue #107](https://github.com/staglibrary/stag/issues/107): Add support for self-loops in graphs
- [Issue #108](https://github.com/staglibrary/stag/issues/108): Allow graphs to be initialised with Laplacian matrix
- [Issue #111](https://github.com/staglibrary/stag/issues/111): Add `is_connected` method to Graph object
- [Issue #113](https://github.com/staglibrary/stag/issues/113): Construct the identity graph
- [Issue #114](https://github.com/staglibrary/stag/issues/114): Implement scalar multiplication of graphs
- [Issue #112](https://github.com/staglibrary/stag/issues/112): Add graph addition operator

### Changed
- [Issue #86](https://github.com/staglibrary/stag/issues/86): Remove test files from the release source archive

### Fixed
- [Issue #87](https://github.com/staglibrary/stag/issues/87): Occasional bug with sorting edgelist file

## [1.2.0] - 2023-03-16
### Added
- STAG Tools for creating graphs and converting between file formats
- [Issue #82](https://github.com/staglibrary/stag/issues/82): Methods for using AdjacencyList files
- The `AdjacencyListLocalGraph` class for using local access to a graph stored on disk
- [Issue #79](https://github.com/staglibrary/stag/issues/79): Conductance method

## [1.1.1] - 2023-03-02
### Fixed
- Reduced the memory usage of the SBM methods
- Correct the SBM approximate sampling distribution

### Added
- Added signless Laplacian methods to the graph object

## [1.1.0] - 2023-02-23
### Added
- [Issue #70](https://github.com/staglibrary/stag/issues/70): Add general stochastic block model method
- [Issue #69](https://github.com/staglibrary/stag/issues/69): Add methods to get the ground truth labels for SBM graphs
- [Issue #71](https://github.com/staglibrary/stag/issues/71): Add stag::adjusted_rand_index method
- [Issue #68](https://github.com/staglibrary/stag/issues/68): Add the power method function

## [1.0.0] - 2023-02-16
### Changed
- [Issue #28](https://github.com/staglibrary/stag/issues/28): Requesting the degree of a non-existant vertex now results
in an error.

### Fixed
- [Issue #56](https://github.com/staglibrary/stag/issues/56): Increase convergence speed of eigenvalue calculation
- [Issue #49](https://github.com/staglibrary/stag/issues/49): Add assertion to pagerank calculation that LocalGraph.degrees
returns a vector with the correct length.

### Added
- [Issue #28](https://github.com/staglibrary/stag/issues/28): Add argument error checking
- [Issue #50](https://github.com/staglibrary/stag/issues/50): Add stag::Graph::average_degree method on stag::Graph object.
- Add stag::LocalGraph::vertex_exists method to stag::LocalGraph object

## [0.4.0] - 2023-01-16
### Fixed
- [Issue #45](https://github.com/staglibrary/stag/issues/45): Make the stag::LocalGraph class pure virtual

### Added
- [Issue #51](https://github.com/staglibrary/stag/issues/51): Add a simple spectral clustering method

## [0.3.0] - 2022-11-18
### Fixed
- [Issue #37](https://github.com/staglibrary/stag/issues/37): Allow float type as target volume in local clustering method

### Added
- [Issue #43](https://github.com/staglibrary/stag/issues/43): Add batched queries for degrees in local graph class
- [Issue #36](https://github.com/staglibrary/stag/issues/36): Allow edgelist entries to be separated by tabs
- [Issue #42](https://github.com/staglibrary/stag/issues/42): Improve the cache-efficiency of the approximate pagerank algorithm
- Add constructor for star graph: `stag::star_graph(n)`

## [0.2.1] - 2022-10-19
### Changed
- [Issue #27](https://github.com/staglibrary/stag/issues/27): Use `long long` integer type throughout the library

### Added
- [Issue #29](https://github.com/staglibrary/stag/issues/29): Add `sprsMatFromVectors` helper method to construct sparse matrices

### Fixed
- [Issue #26](https://github.com/staglibrary/stag/issues/26): Add destructor to abstract `LocalGraph` class

## [0.2.0] - 2022-10-15
### Changed
- [Issue #5](https://github.com/staglibrary/stag/issues/5): Rename the `Graph::volume()` method to `Graph::total_volume()`
- [Issue #5](https://github.com/staglibrary/stag/issues/5): Attempting to construct a graph with an assymetric adjacency matrix
throws an exception
- Sparse matrices within the library are now stored in Column-Major format

### Added
- [Issue #5](https://github.com/staglibrary/stag/issues/5): `Graph::degree_matrix()` method
- [Issue #5](https://github.com/staglibrary/stag/issues/5): `Graph::normalised_laplacian()` method
- [Issue #5](https://github.com/staglibrary/stag/issues/5): `Graph::number_of_vertices()` method
- [Issue #5](https://github.com/staglibrary/stag/issues/5): `Graph::number_of_edges()` method
- [Issue #4](https://github.com/staglibrary/stag/issues/4): Add `load_edgelist()` method to read graphs from edgelist files
- [Issue #4](https://github.com/staglibrary/stag/issues/4): Add `save_edgelist()` method to save graphs to edgelist files
- Add graph equality operators `==` and `!=`
- [Issue #6](https://github.com/staglibrary/stag/issues/6): Add `LocalGraph` abstract class for providing local graph access
- [Issue #15](https://github.com/staglibrary/stag/issues/15): Add methods for generating random graphs from the stochastic
block model
- New `inverse_degree_matrix()` method on the `stag::Graph` object
- New `lazy_random_walk_matrix()` method on the `stag::Graph` object
- [Issue #10](https://github.com/staglibrary/stag/issues/10): Add approximate pagerank methods in `cluster.h`.
- New `addVectors(v1, v2)` utility method
- `barbell_graph(n)` graph constructor
- [Issue #10](https://github.com/staglibrary/stag/issues/10): Add ACL local clustering algorithm in `cluster.h`.

## [0.1.6] - 2022-10-10
### Added
- Graph class with a few basic methods
- `adjacency()` and `laplacian()` methods
- `cycle_graph(n)` and `complete_graph(n)` graph constructors