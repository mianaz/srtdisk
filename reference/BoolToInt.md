# Convert a logical to an integer

Unlike most programming languages, R has three possible
[logical](https://rdrr.io/r/base/logical.html) (boolean) values: `TRUE`,
`FALSE`, and [`NA`](https://rdrr.io/r/base/NA.html); moreover, the `NA`
value has representations in other data types, such as `NA_integer_`,
`NA_real_`, and `NA_character_`. Simply writing out the logical values
to an HDF5 file would cause issues when trying to read the data in to
another language, such as Python. To encode these three logical values
for other languages, we can encode the logicals as integers:

- `FALSE` becomes `0L`

- `TRUE` becomes `1L`

- `NA` becomes `2L`

This encoding scheme allows other languages to handle `NA`s in their own
manner while preserving all three logicals for R

## Usage

``` r
BoolToInt(x)
```

## Arguments

- x:

  A logical vector

## Value

An integer vector

## See also

[integer](https://rdrr.io/r/base/integer.html)
[logical](https://rdrr.io/r/base/logical.html)
[`NA`](https://rdrr.io/r/base/NA.html) `WriteH5Seurat`

## Examples

``` r
# \donttest{
srtdisk:::BoolToInt(x = c(TRUE, FALSE, NA))
#> [1] 1 0 2
# }
```
