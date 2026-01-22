# Transpose a matrix

Transpose a matrix

## Usage

``` r
Transpose(x, ...)

# S3 method for class 'dgCMatrix'
Transpose(x, ...)

# S3 method for class 'H5D'
Transpose(
  x,
  dest = GetParent(x = x),
  dname = paste0("t_", basename(path = x$get_obj_name())),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'H5Group'
Transpose(
  x,
  dest = GetParent(x = x),
  dname = paste0("t_", basename(path = x$get_obj_name())),
  overwrite = FALSE,
  ...
)
```

## Arguments

- x:

  A matrix to transpose

- ...:

  Arguments passed to other methods

- dest:

  ...

- dname:

  ...

- overwrite:

  ...

- verbose:

  Show progress updates

## Value

`dgCMatrix` method: returns a `dgCMatrix` with the data of `x`
transposed

[`H5D`](http://hhoeflin.github.io/hdf5r/reference/H5D-class.md) and
[`H5Group`](http://hhoeflin.github.io/hdf5r/reference/H5Group-class.md)
methods: Invisibly returns `NULL`
