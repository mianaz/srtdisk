# Check to see if a matrix is empty

Determine if a matrix is empty or not. A matrix is considered empty if
it satisfies one of the following conditions:

- The dimensions of the matrix are 0-by-0 (`all(dim(x) == 0)`)

- The dimensions of the matrix are 1-by-1 (`all(dim(x) == 1)`) and the
  sole vlaue is `NA`

These two situations correspond to matrices generated with either
`new('matrix')` or [`matrix()`](https://rdrr.io/r/base/matrix.html)

## Usage

``` r
IsMatrixEmpty(x)
```

## Arguments

- x:

  A matrix

## Value

`TRUE` if the matrix is empty otherwise `FALSE`

## Examples

``` r
# \donttest{
srtdisk:::IsMatrixEmpty(new('matrix'))
#> [1] TRUE
srtdisk:::IsMatrixEmpty(matrix())
#> [1] TRUE
srtdisk:::IsMatrixEmpty(matrix(1:9, nrow = 3))
#> [1] FALSE
# }
```
