# Find the closest version

API changes happen at set versions, and knowing how a current running
version relates to versions introducing API changes is important.
`ClosestVersion` approximages both “rounding down” (eg. to determine
minimum version with new API addition) and “rounding up” (eg. to
determine maximum version before API deletion) for semantic versions.

## Usage

``` r
ClosestVersion(
  query,
  targets,
  direction = c("min", "max"),
  inclusive = direction == "min"
)
```

## Arguments

- query:

  A query version ([`character`](https://rdrr.io/r/base/character.html)
  or [`numeric_version`](https://rdrr.io/r/base/numeric_version.html))

- targets:

  A vector of target versions
  ([`character`](https://rdrr.io/r/base/character.html) or
  [`numeric_version`](https://rdrr.io/r/base/numeric_version.html))

- direction:

  Which way should we check for closest version? Choose from:

  min

  :   Closest version less than or equal to `query`

  max

  :   Closest version greater than or equal to `query`

- inclusive:

  Perform an inclusive comparison (eg. `>=` or `<=` versus to `>` or
  `<`) for “rounding”

## Value

The version from `targets` that is closest to `query` as a
[`character`](https://rdrr.io/r/base/character.html) vector

## See also

[`numeric_version`](https://rdrr.io/r/base/numeric_version.html)

## Examples

``` r
# \donttest{
srtdisk:::ClosestVersion('3.1.0', targets = c('3.0.0', '1.4.9', '4.3.2'))
#> [1] "3.0.0"
srtdisk:::ClosestVersion('3.1.0', targets = c('3.0.0', '1.4.9', '4.3.2'), direction = 'max')
#> [1] "4.3.2"
# }
```
