# Check the datatype of an HDF5 dataset

Effectively, an implementation of
[`is`](https://rdrr.io/r/methods/is.html) for HDF5 datasets; useful to
ensure HDF5 validity for specific file structures

## Usage

``` r
IsDType(x, dtype)
```

## Arguments

- x:

  An HDF5 dataset (object of type
  [`H5D`](http://hhoeflin.github.io/hdf5r/reference/H5D-class.md))

- dtype:

  A character vector of HDF5 datatype names, must be present in
  [`h5types`](http://hhoeflin.github.io/hdf5r/reference/h5types.md)

## Value

A logical

## See also

[`h5types`](http://hhoeflin.github.io/hdf5r/reference/h5types.md)
