# Does an R6 class inherit from scdisk

Does an R6 class inherit from scdisk

## Usage

``` r
IsSCDisk(r6class)
```

## Arguments

- r6class:

  An R6 class generator or a character name of an R6 class generator

## Value

If `r6class` inherits from scdisk, returns `TRUE`; otherwise, returns
`FALSE`

## Examples

``` r
# \donttest{
srtdisk:::IsSCDisk("H5File")
#> [1] FALSE
srtdisk:::IsSCDisk("scdisk")
#> [1] TRUE
srtdisk:::IsSCDisk("h5Seurat")
#> [1] TRUE
# }
```
