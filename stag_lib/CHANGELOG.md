# Changelog
All notable changes to the library are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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