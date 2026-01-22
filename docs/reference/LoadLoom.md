# Loom-file Loading

Load data from a loom file into a
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
object

## Usage

``` r
LoadLoom(
  file,
  assay = NULL,
  cells = "CellID",
  features = "Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE,
  ...
)

# S3 method for class 'character'
LoadLoom(file, ...)

# S3 method for class 'H5File'
LoadLoom(file, ...)

# S3 method for class 'loom'
LoadLoom(file, ...)

# S3 method for class 'loom'
as.Seurat(
  x,
  assay = NULL,
  cells = "CellID",
  features = "Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE,
  ...
)
```

## Arguments

- file, x:

  Name of loom file or a `loom` object to load data from

- assay:

  Name of assay to store expression data as; if `NULL`, will search for
  an HDF5 attribute named `SEURAT_ASSAY` or an attribute dataset named
  `/attrs/SEURAT_ASSAY` for assay name. If not found, defaults to “RNA”

- cells:

  Name of dataset in `/col_attrs` with cell names

- features:

  Name of dataset in `/row_attrs` with feature names

- normalized:

  Name of matrix in `/layers` to store normalized data as; pass
  “/matrix” to store `/matrix` as normalized data instead of raw counts

- scaled:

  Name of dataset in `/layers` to store scaled data as

- filter:

  Keep only selected cells and/or features as specified by
  `/col_attrs/Valid` and `/row_attrs/Valid`, respectively

- verbose:

  Show progress updates

- ...:

  Arguments passed to other methods

## Value

A [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
object

## Details

`LoadLoom` will try to automatically fill slots of a `Seurat` object
based on data presence or absence in a given loom file. This method
varies by loom specification version. For version-specific details, see
sections below

## Loom 0.1 Loading

Loading data from loom files less than version 3.0.0 is not currently
supported

## Loom 3.0.0 Loading

blah

## See also

[Loom file
conventions](http://linnarssonlab.org/loompy/conventions/index.md)
