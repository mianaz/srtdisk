# Generate an HDF5 string dtype

Presets for encoding variations of
[`H5T_STRING`](http://hhoeflin.github.io/hdf5r/reference/H5T_STRING-class.md);
used to generate HDF5 datatype specifications with specific string
encodings

## Usage

``` r
StringType(stype = c("utf8", "ascii7"))
```

## Arguments

- stype:

  Type of string encoding to use, choose from:

  utf8

  :   Variable-width, UTF-8

  ascii7

  :   Fixed-width (7 bits), ASCII

## Value

An
[`H5T_STRING`](http://hhoeflin.github.io/hdf5r/reference/H5T_STRING-class.md)
object

## See also

[`H5T_STRING`](http://hhoeflin.github.io/hdf5r/reference/H5T_STRING-class.md)

## Examples

``` r
# \donttest{
srtdisk:::StringType()
srtdisk:::StringType('ascii7')
#> Class: H5T_STRING
#> Datatype: H5T_STRING {
#>       STRSIZE 7;
#>       STRPAD H5T_STR_NULLTERM;
#>       CSET H5T_CSET_ASCII;
#>       CTYPE H5T_C_S1;
#>    }
# }
```
