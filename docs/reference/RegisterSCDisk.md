# Get and Register [`scdisk`](https://mianaz.github.io/srtdisk/reference/scdisk-class.md) Subclasses

Mechanisms for registration of
[`scdisk`](https://mianaz.github.io/srtdisk/reference/scdisk-class.md)
subclass generators for use in functions that rely on the class
definition instead of an object.

## Usage

``` r
GetSCDisk(r6class = NULL)

RegisterSCDisk(r6class)
```

## Arguments

- r6class:

  An R6 class generator or a character name of an R6 class generator

## Value

`GetSCDisk`: if `r6class` is `NULL`, then a vector of all registered
`scdisk` subclasses; otherwise, a generator for the requested `scdisk`
subclass

`RegisterSCDisk`: adds `r6class` to the internal subclass registry and
invisibly returns `NULL`

## Details

While `scdisk`-subclassed objects (eg.
[`h5Seurat`](https://mianaz.github.io/srtdisk/reference/h5Seurat-class.md)
objects) follow traditional inheritance patterns (can be determined
through [`inherits`](https://rdrr.io/r/base/class.html)), the class
definitions and object generators do not. These functions provide a
simple mechanism for adding and getting the defintions of `scdisk`
subclasses for functions that utilize the object generators or other
aspects of the class definition (such as
[`Convert`](https://mianaz.github.io/srtdisk/reference/Convert.md))

To register a subclass of `scdisk`, simply add a call to
`RegisterSCDisk` in your load hook


    .onLoad <- function(libname, pkgname) {
      RegisterSCDisk(classgen)
      # Other code to be run on load
    }

## Examples

``` r
GetSCDisk()
#> [1] "h5Seurat" "loom"    
GetSCDisk("h5Seurat")
#> <h5Seurat> object generator
#>   Inherits from: <scdisk>
#>   Public:
#>     index: function () 
#>     set.version: function (version) 
#>     version: function () 
#>     detect.version: function () 
#>     is.v5: function () 
#>   Private:
#>     index.internal: list
#>     versions: 3.1.2 3.1.5.9900 5.0.0 5.2.1
#>     build.index: function (version) 
#>     create: function (version, verbose = TRUE) 
#>     validate: function (verbose = TRUE, ...) 
#>     v3.1.2: function () 
#>     v3.2.0: function () 
#>     v5.0.0: function () 
#>   Parent env: <environment: namespace:srtdisk>
#>   Locked objects: TRUE
#>   Locked class: TRUE
#>   Portable: TRUE

if (FALSE) { # \dontrun{
RegisterSCDisk(h5Seurat)
} # }
```
