# Update a Seurat key

Attempts to validate a string to use as a Seurat key. Valid keys must
match the regular expression `^[[:alnum:]]+_$`; if `key` fails this
regular expression, an attempt to modify it to said key will be made by
removing all non-alphanumeric characters, collapsing the resulting
vector, and appending “\_”. If this stil fails, a random string of
lowercase characters will be generated, followed by “\_”, to be used as
the key

## Usage

``` r
UpdateKey(key)
```

## Arguments

- key:

  A key to validate and update

## Value

`key`, updated if invalid

## See also

[`Key`](https://satijalab.github.io/seurat-object/reference/Key.html)
[`RandomName`](https://mianaz.github.io/srtdisk/reference/RandomName.md)

## Examples

``` r
# \donttest{
srtdisk:::UpdateKey("RNA_")
#> [1] "RNA_"
srtdisk:::UpdateKey("potato")
#> [1] "potato_"
srtdisk:::UpdateKey("*@)")
#> [1] "vju_"
# }
```
