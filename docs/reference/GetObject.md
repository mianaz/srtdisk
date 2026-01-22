# Figure out which objects to load from an h5Seurat file

Figure out which objects to load from an h5Seurat file

## Usage

``` r
GetAssays(assays, index)

GetCommands(index, assays = NULL)

GetDimReducs(reductions, index, assays = NULL)

GetGraphs(graphs, index, assays = NULL)

GetImages(images, index, assays = NULL)

GetNeighbors(neighbors, index)
```

## Arguments

- assays:

  One of:

  - A character vector with names of assays

  - A character vector with one or more of `counts`, `data`,
    `scale.data` describing which slots of **all assays** to load

  - A named list where each entry is either the name of an assay or a
    vector describing which slots (described above) to take from which
    assay

  - `NULL` for all assays

- index:

  An h5Seurat index
  ([`h5SI`](https://mianaz.github.io/srtdisk/reference/h5SI.md)) object

- reductions:

  One of:

  - A character vector with names of reductions

  - `NULL` for all reductions

  - `NA` for
    [global](https://satijalab.github.io/seurat-object/reference/IsGlobal.html)
    reductions

  - `FALSE` for no reductions

  **Note**: Only reductions associated with an assay loaded in `assays`
  or marked as
  [global](https://satijalab.github.io/seurat-object/reference/IsGlobal.html)
  will be loaded

- graphs:

  One of:

  - A character vector with names of graphs

  - `NULL` for all graphs

  - `FALSE` for no graphs

  **Note**: Only graphs associated with an assay loaded in `assays` will
  be loaded

- images:

  One of:

  - A character vector with names of images

  - `NULL` for all images

  - `NA` for
    [global](https://satijalab.github.io/seurat-object/reference/IsGlobal.html)
    images

  - `FALSE` for no images

- neighbors:

  One of:

  - A character vector with the names of neighbors

  - `NULL` for all neighbors

  - `FALSE` for no neighbors

## Value

`GetAssays`: A named list where each entry is a vector describing the
slots of an assay to load and the names are the assays to load

`GetCommands`: A vector of command log names that are derived from an
assay in `assay`

`GetDimReducs`: A vector of reduction names that are derived from an
assay in `assays` or global dimensional reductions

`GetGraphs`: A vector of graph names that are derived from an assay in
`assays`

`GetImages`: A vector of image names

`GetNeighbors`: A vector of neighbor names

## See also

[`LoadH5Seurat`](https://mianaz.github.io/srtdisk/reference/LoadH5Seurat.md)
