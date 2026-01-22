# Generate chunk points

Generate chunk points

## Usage

``` r
ChunkPoints(dsize, csize)
```

## Arguments

- dsize:

  Size of data being chunked

- csize:

  Size of chunk; if `NA`, assumes single chunk

## Value

A matrix where each row is a chunk, column 1 is start points, column 2
is end points

## Examples

``` r
# \donttest{
srtdisk:::ChunkPoints(100, 3)
#>       start end
#>  [1,]     1   3
#>  [2,]     4   6
#>  [3,]     7   9
#>  [4,]    10  12
#>  [5,]    13  15
#>  [6,]    16  18
#>  [7,]    19  21
#>  [8,]    22  24
#>  [9,]    25  27
#> [10,]    28  30
#> [11,]    31  33
#> [12,]    34  36
#> [13,]    37  39
#> [14,]    40  42
#> [15,]    43  45
#> [16,]    46  48
#> [17,]    49  51
#> [18,]    52  54
#> [19,]    55  57
#> [20,]    58  60
#> [21,]    61  63
#> [22,]    64  66
#> [23,]    67  69
#> [24,]    70  72
#> [25,]    73  75
#> [26,]    76  78
#> [27,]    79  81
#> [28,]    82  84
#> [29,]    85  87
#> [30,]    88  90
#> [31,]    91  93
#> [32,]    94  96
#> [33,]    97  99
#> [34,]   100 100
srtdisk:::ChunkPoints(100, NA)
#>      start end
#> [1,]     1 100
# }
```
