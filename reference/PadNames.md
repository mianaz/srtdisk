# Add names for unnamed or partially named objects

Add names for unnamed or partially named objects

## Usage

``` r
PadNames(x, prefix = "index")
```

## Arguments

- x:

  An object that can be named

- prefix:

  A prefix to be added to each name

## Value

`x` with unnamed values named

## Examples

``` r
# \donttest{
a <- list(1, b = 2, 3)
srtdisk:::PadNames(a)
#> $index1
#> [1] 1
#> 
#> $b
#> [1] 2
#> 
#> $index3
#> [1] 3
#> 
# }
```
