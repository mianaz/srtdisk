# Write data to an HDF5 group

Writing data to HDF5 files can be done simply with usually sensible
defaults. However, when wanting any semblance of control over how an R
object is written out, the code constructs get complicated quickly.
`WriteH5Group` provides a wrapper with sensible defaults over some of
these complex code constructs to provide greater control over how data
are written to disk. These defaults were chosen to fit best with
[h5Seurat](https://mianaz.github.io/srtdisk/reference/h5Seurat-class.md)
files, see
[vignette("h5Seurat-spec")](https://mianaz.github.io/srtdisk/doc/h5Seurat-spec.md)
for more details

## Usage

``` r
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'ANY'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'array'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'Assay'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'Assay5'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'data.frame'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'dgCMatrix'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'DimReduc'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'factor'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'Graph'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'list'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'logical'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'Neighbor'
WriteH5Group(x, name, hgroup, verbose = TRUE)

# S4 method for class 'SeuratCommand'
WriteH5Group(x, name, hgroup, verbose = TRUE)
```

## Arguments

- x:

  An object

- name:

  Name to save data as

- hgroup:

  An HDF5 file or group (`H5File` or `H5Group` objects from hdf5r)

- verbose:

  Show progress updates

## Value

Invisibly returns `NULL`

## Examples

``` r
# \donttest{
# Setup an HDF5 file
hfile <- hdf5r::H5File$new(filename = tempfile(fileext = '.h5'), mode = 'a')
# }

# \donttest{
# Data frames are stored as either datasets or groups, depending on the
# presence of factor columns
df <- data.frame(
  x = c('g1', 'g1', 'g2', 'g1', 'g2'),
  y = 1:5,
  stringsAsFactors = FALSE
)

# When no factor columns are present, the data frame is written as a single
# HDF5 compound dataset
WriteH5Group(x = df, name = 'df', hgroup = hfile)
hfile[['df']]
#> Class: H5Group
#> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/Rtmpouj95k/file48373725e429.h5
#> Group: /df
#> Attributes: colnames, _index
#> Listing:
#>    name    obj_type dataset.dims dataset.type_class
#>  _index H5I_DATASET            5         H5T_STRING
#>       x H5I_DATASET            5         H5T_STRING
#>       y H5I_DATASET            5        H5T_INTEGER

# When factors are present, the data frame is written as a group
# This is because h5py does not implement HDF5 Enums, so factor level
# information would be lost
df$x <- factor(x = df$x)
WriteH5Group(x = df, name = 'df.factor', hgroup = hfile)
hfile[['df.factor']]
#> Class: H5Group
#> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/Rtmpouj95k/file48373725e429.h5
#> Group: /df.factor
#> Attributes: colnames, _index
#> Listing:
#>    name    obj_type dataset.dims dataset.type_class
#>  _index H5I_DATASET            5         H5T_STRING
#>       x   H5I_GROUP         <NA>               <NA>
#>       y H5I_DATASET            5        H5T_INTEGER
# }

# \donttest{
# Factors turn into a group with two components: values and levels
# This is to preserve level information for HDF5 APIs that don't implement
# the HDF5 Enum type (eg. h5py)
# values corresponds to the integer values of each member of a factor
# levels is a string dataset with one entry per level
fctr <- factor(x = c('g1', 'g1', 'g2', 'g1', 'g2'))
WriteH5Group(x = fctr, name = 'factor', hgroup = hfile)
hfile[['factor']]
#> Class: H5Group
#> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/Rtmpouj95k/file48373725e429.h5
#> Group: /factor
#> Listing:
#>    name    obj_type dataset.dims dataset.type_class
#>  levels H5I_DATASET            2         H5T_STRING
#>  values H5I_DATASET            5        H5T_INTEGER
# }

# \donttest{
# Logicals get encoded as integers with the following mapping
# FALSE becomes 0L
# TRUE becomes 1L
# NA becomes 2L
# These are stored as H5T_INTEGERS instead of H5T_LOGICALS
# Additionally, an attribute called "s3class" is written with the value of "logical"
WriteH5Group(c(TRUE, FALSE, NA), name = "logicals", hgroup = hfile)
hfile[["logicals"]]
#> Class: H5D
#> Dataset: /logicals
#> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/Rtmpouj95k/file48373725e429.h5
#> Access type: H5F_ACC_RDWR
#> Attributes: s3class
#> Datatype: H5T_STD_I32LE
#> Space: Type=Simple     Dims=3     Maxdims=Inf
#> Chunk: 2048
hfile[["logicals"]]$attr_open("s3class")$read()
#> [1] "logical"
# }

# \donttest{
# Close and remove the HDF5 file
hfile$close_all()
file.remove(hfile$filename)
#> [1] TRUE
# }
```
