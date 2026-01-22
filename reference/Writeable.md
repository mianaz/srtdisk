# Is an HDF5 file or group writeable

Is an HDF5 file or group writeable

## Usage

``` r
Writeable(x, error = TRUE)
```

## Arguments

- x:

  An
  [`H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md)
  or
  [`H5Group`](http://hhoeflin.github.io/hdf5r/reference/H5Group-class.md)
  object

- error:

  Throw an error if `x` is not writeable

## Value

`TRUE` if `x` is writeable otherwise `FALSE`
