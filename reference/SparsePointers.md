# Convert sparse matrix pointers to indices and vice versa

Convert sparse matrix pointers to indices and vice versa

## Usage

``` r
IndexToPointer(j)

PointerToIndex(p)
```

## Source

`PointerToIndex` came from
[StackOverflow](https://stackoverflow.com/questions/20008200/r-constructing-sparse-matrix)

## Arguments

- j:

  A vector of sparse matrix colum indices

- p:

  A vector of sparse matrix pointers

## Value

`IndexToPointer`: A vector of index pointers (p)

`PointerToIndex`: A vector of column (j) indices

## Author

`PointerToIndex` was written by [Josh O'Brien on
StackOverflow](https://stackoverflow.com/users/980833/josh-obrien)

## Examples

``` r
# \donttest{
dat <- dat <- c(0, 0, 1, 4, 0, 2, 0, 9, 0)
smat <- Matrix::Matrix(data = dat, nrow = 3, sparse = TRUE)
j <- srtdisk:::PointerToIndex(p = smat@p)
Matrix::sparseMatrix(i = smat@i + 1, j = j, x = smat@x)
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>           
#> [1,] . 4 .
#> [2,] . . 9
#> [3,] 1 2 .
p <- srtdisk:::IndexToPointer(j = j)
Matrix::sparseMatrix(i = smat@i + 1, p = p, x= smat@x)
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>           
#> [1,] . 4 .
#> [2,] . . 9
#> [3,] 1 2 .
# }
```
