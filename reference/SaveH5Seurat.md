# Save a `Seurat` object to an h5Seurat file

Save a `Seurat` object to an h5Seurat file

## Usage

``` r
SaveH5Seurat(object, filename, overwrite = FALSE, verbose = TRUE, ...)

as.h5Seurat(x, ...)

# Default S3 method
SaveH5Seurat(object, filename, overwrite = FALSE, verbose = TRUE, ...)

# S3 method for class 'Seurat'
SaveH5Seurat(
  object,
  filename = paste0(Project(object = object), ".h5Seurat"),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# Default S3 method
as.h5Seurat(x, filename, overwrite = FALSE, verbose = TRUE, ...)

# S3 method for class 'H5File'
as.h5Seurat(x, ...)

# S3 method for class 'Seurat'
as.h5Seurat(
  x,
  filename = paste0(Project(object = x), ".h5seurat"),
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

`SaveH5Seurat`: Invisbly returns `filename`

`as.h5Seurat`: An
[`h5Seurat`](https://mianaz.github.io/srtdisk/reference/h5Seurat-class.md)
object
