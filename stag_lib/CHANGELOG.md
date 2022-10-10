# Changelog
All notable changes to the library are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- [Issue #5](https://github.com/pmacg/stag/issues/5): Rename the `Graph::volume()` method to `Graph::total_volume()`

### Added
- [Issue #5](https://github.com/pmacg/stag/issues/5): `Graph::degree_matrix()` method

## [0.1.6] - 2022-10-10
### Added
- Graph class with a few basic methods
- `adjacency()` and `laplacian()` methods
- `cycle_graph(n)` and `complete_graph(n)` graph constructors