# Save a [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html) object to a loom file

Save a
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
object to a loom file

## Usage

``` r
SaveLoom(object, filename, overwrite = FALSE, verbose = TRUE, ...)

as.loom(x, ...)

# Default S3 method
SaveLoom(object, filename, overwrite = FALSE, verbose = TRUE, ...)

# S3 method for class 'Seurat'
SaveLoom(
  object,
  filename = paste0(Project(object = object), ".loom"),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# Default S3 method
as.loom(x, filename, overwrite = FALSE, verbose = TRUE, ...)

# S3 method for class 'H5File'
as.loom(x, ...)

# S3 method for class 'Seurat'
as.loom(
  x,
  filename = paste0(Project(object = x), ".loom"),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object, x:

  An object

- filename:

  Name of file to save the object to

- overwrite:

  Overwrite `filename` if present

- verbose:

  Show progress updates

- ...:

  Arguments passed to other methods

## Value

`SaveLoom`: Invisibly returns `filename`

`as.loom`: A
[`loom`](https://mianaz.github.io/srtdisk/reference/loom-class.md)
object
