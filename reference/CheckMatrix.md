# Check that a dataset is a proper loom matrix

Check that a dataset is a proper loom matrix

## Usage

``` r
CheckMatrix(lfile, name, dims = NULL)
```

## Arguments

- lfile:

  A [`loom`](https://mianaz.github.io/srtdisk/reference/loom-class.md)
  or
  [`H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md)
  object

- name:

  Name of matrix to check

- dims:

  If provided, ensure `lfile[[name]]` has these dimensions; should be a
  two-dimensional numeric vector with ncells/ncol as the first value and
  nfeature/nrow as the second

## Value

If all checks pass successfully, invisibly returns `name`
