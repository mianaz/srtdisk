# Append data from an h5Seurat file to a preexisting [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html) object

Append data from an h5Seurat file to a preexisting
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
object

## Usage

``` r
AppendData(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = "commands",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'character'
AppendData(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = "commands",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'H5File'
AppendData(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = "commands",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'h5Seurat'
AppendData(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = "commands",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- file:

  Path to an h5Seurat file or an open h5Seurat object

- object:

  A
  [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
  object to append data to

- assays:

  One of:

  - A character vector with names of assays

  - A character vector with one or more of `counts`, `data`,
    `scale.data` describing which slots of **all assays** to load

  - A named list where each entry is either the name of an assay or a
    vector describing which slots (described above) to take from which
    assay

  - `NULL` for all assays

  - `FALSE` for no assays

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

- extras:

  Extra information to load; supports any combination of the following
  values:

  “commands”

  :   Load command logs. If `overwrite = TRUE`, replaces existing
      command logs

- overwrite:

  Overwrite existing data in `object` with data from `file`

- verbose:

  Show progress updates

- ...:

  Arguments passed to other methods

## Value

`object` with the extra data requested
