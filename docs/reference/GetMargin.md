# Determine the margin to use for a dataset

Determine the margin to use for a dataset

## Usage

``` r
GetMargin(dims, MARGIN = getOption(x = "SeuratDisk.chunking.MARGIN"))
```

## Arguments

- dims:

  Dimensions of a dataset

- MARGIN:

  Either an integer value contained within `1:length(x = dims)` or one
  of the possible values of `SeuratDisk.chunking.MARGIN`

## Value

An integer value with the `MARGIN`

## See also

`SeuratDisk.chunking.MARGIN`

## Examples

``` r
# \donttest{
srtdisk:::GetMargin(c(4, 10))
#> [1] 2
# }
```
