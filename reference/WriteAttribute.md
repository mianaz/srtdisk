# Write an attribute to an HDF5 file, group, or dataset

Write an attribute to an HDF5 file, group, or dataset

## Usage

``` r
WriteAttribute(x, name, lfile, stype, transpose = TRUE, verbose = TRUE)
```

## Arguments

- x:

  An object to write out

- name:

  Name to store attribute as

- lfile:

  An HDF5
  [file](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md),
  [group](http://hhoeflin.github.io/hdf5r/reference/H5Group-class.md),
  or [dataset](http://hhoeflin.github.io/hdf5r/reference/H5D-class.md)

- stype:

  Data type of attribute

- transpose:

  Whether to transpose the data

- verbose:

  Whether to print messages

## Value

Invisibly returns `NULL`
