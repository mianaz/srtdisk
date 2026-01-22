# Load data from an HDF5 File

HDF5 allows storing data in an arbitrary fashion, which makes reading
data into memory a hassle. The methods here serve as convenience
functions for reading data stored in a certain format back into a
certain R object. For details regarding how data should be stored on
disk, please see the [h5Seurat file
specification](https://mianaz.github.io/srtdisk/doc/h5Seurat-spec.md).

## Usage

``` r
# S3 method for class 'H5D'
as.array(x, ...)

# S3 method for class 'H5D'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'H5Group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S4 method for class 'H5Group'
as.factor(x)

# S4 method for class 'H5Group'
as.list(x, which = NULL, ...)

# S3 method for class 'H5D'
as.matrix(x, transpose = FALSE, ...)

# S3 method for class 'H5Group'
as.matrix(x, ...)

# S3 method for class 'H5D'
as.sparse(x, ...)

# S3 method for class 'H5Group'
as.sparse(x, ...)
```

## Arguments

- x:

  An HDF5 dataset or group

- ...:

  Arguments passed to other methods

- which:

  Optional character vector specifying which elements to read from the
  H5Group

- row.names:

  `NULL` or a character vector giving the row names for the data frame.
  Missing values are not allowed.

- optional:

  logical. If `TRUE`, setting row names and converting column names (to
  syntactic names: see
  [`make.names`](https://rdrr.io/r/base/make.names.html)) is optional.
  Note that all of R's base package `as.data.frame()` methods use
  `optional` only for column names treatment, basically with the meaning
  of
  [`data.frame`](https://rdrr.io/r/base/data.frame.html)`(*, check.names = !optional)`.
  See also the `make.names` argument of the `matrix` method.

- transpose:

  Transpose the data upon reading it in, used when writing data in
  row-major order (eg. from C or Python)

## Value

`as.array`: returns an [`array`](https://rdrr.io/r/base/array.html) with
the data from the HDF5 dataset

`as.data.frame`: returns a
[`data.frame`](https://rdrr.io/r/base/data.frame.html) with the data
from the HDF5 dataset or group

`as.factor`: returns a [`factor`](https://rdrr.io/r/base/factor.html)
with the data from the HDF5 group

`as.list`: returns a [`list`](https://rdrr.io/r/base/list.html) with the
data from the HDF5 group

`as.logical`: returns a [`logical`](https://rdrr.io/r/base/logical.html)
with the data from the HDF5 dataset

`as.matrix`, `H5D` method: returns a
[`matrix`](https://rdrr.io/r/base/matrix.html) with the data from the
HDF5 dataset

`as.sparse`, `H5D` method: returns a sparse matrix with the data from
the HDF5 dataset

`as.sparse`, `as.matrix`, `H5Group` method: returns a
[`sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) with
the data from the HDF5 group

`dimnames`: returns a two-length list of character vectors for row and
column names. Row names should be in a column named `index`
