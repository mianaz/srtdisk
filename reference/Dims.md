# Get the dimensions of an HDF5 dataset or sparse matrix

Get the dimensions of an HDF5 dataset or sparse matrix

## Usage

``` r
Dims(x)
```

## Arguments

- x:

  An HDF5 dataset or sparse matrix

## Value

A vector with the dimensions of the dataset or sparse matrix. For sparse
matrices, if no dimensions are found in either the “dims” or “shape”
attribute, returns `c(NA_integer_, NA_integer_)`

## See also

[vignette("h5Seurat-spec")](https://mianaz.github.io/srtdisk/doc/h5Seurat-spec.md)
