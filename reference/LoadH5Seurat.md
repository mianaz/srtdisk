# Load a saved `Seurat` object from an h5Seurat file

Load a saved `Seurat` object from an h5Seurat file

## Usage

``` r
LoadH5Seurat(file, ...)

# S3 method for class 'character'
LoadH5Seurat(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
)

# S3 method for class 'H5File'
LoadH5Seurat(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
)

# S3 method for class 'h5Seurat'
LoadH5Seurat(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
)

# S3 method for class 'h5Seurat'
as.Seurat(
  x,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = TRUE,
  tools = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- file, x:

  Name of h5Seurat or connected h5Seurat file to load

- ...:

  Arguments passed to other methods

- assays:

  One of:

  - A character vector with names of assays

  - A character vector with one or more of `counts`, `data`,
    `scale.data` describing which slots of **all assays** to load

  - A named list where each entry is either the name of an assay or a
    vector describing which slots (described above) to take from which
    assay

  - `NULL` for all assays

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

- neighbors:

  One of:

  - A character vector with the names of neighbors

  - `NULL` for all neighbors

  - `FALSE` for no neighbors

- images:

  One of:

  - A character vector with names of images

  - `NULL` for all images

  - `NA` for
    [global](https://satijalab.github.io/seurat-object/reference/IsGlobal.html)
    images

  - `FALSE` for no images

- meta.data:

  Load object metadata

- commands:

  Load command information  
  **Note**: only commands associated with an assay loaded in `assays`
  will be loaded

- misc:

  Load miscellaneous data

- tools:

  Load tool-specific information

- verbose:

  Show progress updates

## Value

A `Seurat` object with the data requested
