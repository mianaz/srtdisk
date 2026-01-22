# Guess an HDF5 Datatype

Wrapper around
[`hdf5r::guess_dtype`](http://hhoeflin.github.io/hdf5r/reference/guess_dtype.md),
allowing for the customization of string types rather than defaulting to
variable-length ASCII-encoded strings. Also encodes logicals as
[`H5T_INTEGER`](http://hhoeflin.github.io/hdf5r/reference/H5T_INTEGER-class.md)
instead of
[`H5T_LOGICAL`](http://hhoeflin.github.io/hdf5r/reference/H5T_LOGICAL-class.md)
to ensure cross-language compatibility (controlled via [package
options](https://mianaz.github.io/srtdisk/reference/SeuratDisk-package.md))

## Usage

``` r
GuessDType(x, stype = "utf8", ...)
```

## Arguments

- x:

  The object for which to guess the HDF5 datatype or the dimension or
  the number of elements

- stype:

  Type of string encoding to use, choose from:

  utf8

  :   Variable-width, UTF-8

  ascii7

  :   Fixed-width (7 bits), ASCII

- ...:

  Arguments passed on to
  [`hdf5r::guess_dtype`](http://hhoeflin.github.io/hdf5r/reference/guess_dtype.md)

  `ds_dim`

  :   Can explicitly set the dimension of the dataset object. For
      `scalar`, this is one. Otherwise, this can be used so that a
      multi-dimensional object can be represented so that some of its
      dimension are in the dataset, and some are inside an
      [`H5T_ARRAY`](http://hhoeflin.github.io/hdf5r/reference/H5T_ARRAY-class.md)

  `scalar`

  :   Should the datatype be created so that `x` can be represented as a
      scalar with that datatype? This is intended to know if a
      vector/array should be represented as an
      [`H5T_ARRAY`](http://hhoeflin.github.io/hdf5r/reference/H5T_ARRAY-class.md)
      or not.

  `string_len`

  :   If a string is in the R object, the length to which the
      corresponding HDF5 type should be set. If it is a positive
      integer, the string is of that length. If it is `Inf`, it is
      variable length. If it is set to `estimate`, it is set to the
      length of the longest string in the `x`.

## Value

An object of class
[`H5T`](http://hhoeflin.github.io/hdf5r/reference/H5T-class.md)

## See also

[`guess_dtype`](http://hhoeflin.github.io/hdf5r/reference/guess_dtype.md)
[`BoolToInt`](https://mianaz.github.io/srtdisk/reference/BoolToInt.md)
[`StringType`](https://mianaz.github.io/srtdisk/reference/StringType.md)

## Examples

``` r
# \donttest{
# Characters can either be variable-width UTF8-encoded or
# fixed-width ASCII-encoded
srtdisk:::GuessDType(x = 'hello')
#> Class: H5T_STRING
#> Datatype: H5T_STRING {
#>       STRSIZE H5T_VARIABLE;
#>       STRPAD H5T_STR_NULLTERM;
#>       CSET H5T_CSET_UTF8;
#>       CTYPE H5T_C_S1;
#>    }
srtdisk:::GuessDType(x = 'hello', stype = 'ascii7')
#> Class: H5T_STRING
#> Datatype: H5T_STRING {
#>       STRSIZE 7;
#>       STRPAD H5T_STR_NULLTERM;
#>       CSET H5T_CSET_ASCII;
#>       CTYPE H5T_C_S1;
#>    }

# Data frames are a compound type; character columns follow the same rules
# as character vectors
df <- data.frame(x = c('g1', 'g2', 'g3'), y = 1, 2, 3, stringsAsFactors = FALSE)
srtdisk:::GuessDType(x = df)
#> Class: H5T_COMPOUND
#> Datatype: H5T_COMPOUND {
#>       H5T_STRING {
#>          STRSIZE H5T_VARIABLE;
#>          STRPAD H5T_STR_NULLTERM;
#>          CSET H5T_CSET_UTF8;
#>          CTYPE H5T_C_S1;
#>       } "x" : 0;
#>       H5T_IEEE_F64LE "y" : 8;
#>       H5T_IEEE_F64LE "X2" : 16;
#>       H5T_IEEE_F64LE "X3" : 24;
#>    }
srtdisk:::GuessDType(x = df, stype = 'ascii7')
#> Class: H5T_COMPOUND
#> Datatype: H5T_COMPOUND {
#>       H5T_STRING {
#>          STRSIZE 7;
#>          STRPAD H5T_STR_NULLTERM;
#>          CSET H5T_CSET_ASCII;
#>          CTYPE H5T_C_S1;
#>       } "x" : 0;
#>       H5T_IEEE_F64LE "y" : 7;
#>       H5T_IEEE_F64LE "X2" : 15;
#>       H5T_IEEE_F64LE "X3" : 23;
#>    }

# Logicals are turned into integers to ensure compatibility with Python
# TRUE evaluates to 1, FALSE to 0, and NA to 2
srtdisk:::GuessDType(x = c(TRUE, FALSE, NA))
#> Class: H5T_INTEGER
#> Datatype: H5T_STD_I32LE
# }
```
