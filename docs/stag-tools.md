STAG Tools {#stag-tools}
==========

The STAG library includes several command line tools for working with graph files.
These tools are automatically installed with the library.

Generating random graphs
------------------------
The `stag_sbm` command line tool can be used to generate random graphs from
the symmetric stochastic block model.

### Usage

```bash
stag_sbm [filename] [n] [k] [p] [q]
```

Creates a new EdgeList file named `[filename]` containing a graph drawn from the
symmetric stochastic block model with parameters `[n]`, `[k]`, `[p]`, and `[q]`.
Equivalent to calling the stag::sbm method and saving the resulting graph as
an EdgeList.

Converting between EdgeList and AdjacencyList
---------------------------------------------
STAG supports two file formats (see [Graph File Formats](@ref file-formats)),
and the `stag_edge2adj` and `stag_adj2edge` command line tools
can be used to convert between them.

### Usage

```bash
stag_edge2adj [edgelist] [adjacencylist]
```

Converts the EdgeList file to a new AdjacencyList file.
Equivalent to calling stag::edgelist_to_adjacencylist.

```bash
stag_adj2edge [adjacencylist] [edgelist]
```

Converts the AdjacencyList file to a new EdgeList file.
Equivalent to calling stag::adjacencylist_to_edgelist.
