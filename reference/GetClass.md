# Get a class string with package information

S4 classes are useful in the context of their defining package (benefits
of stricter typing). In order to ensure class information is properly
retained in HDF5 files, S4 class names are written as
“package:classname” with certain exceptions (eg. S4 classes defined by
[Seurat](https://satijalab.org/seurat/reference/Seurat-package.html))

## Usage

``` r
GetClass(class, packages = "Seurat")
```

## Arguments

- class:

  Class name

- packages:

  A vector of packages to exclude from resulting class information

## Value

A character vector with the class

## Examples

``` r
# \donttest{
srtdisk:::GetClass('Seurat')
#> [1] "SeuratObject:Seurat"
srtdisk:::GetClass('Matrix')
#> [1] "Matrix:Matrix"
# }
```
