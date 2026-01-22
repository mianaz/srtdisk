# Generate a random string of characters

Generate a random string of characters

## Usage

``` r
RandomName(length = 5L, ...)
```

## Arguments

- length:

  Length ([`nchar`](https://rdrr.io/r/base/nchar.html)) of string to
  generate

- ...:

  Extra parameters passed to
  [`sample`](https://rdrr.io/r/base/sample.html)

## Value

A random string of characters of length
([`nchar`](https://rdrr.io/r/base/nchar.html)) of `length`

## Examples

``` r
# \donttest{
srtdisk:::RandomName()
#> [1] "lbfdp"
# }
```
