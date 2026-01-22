# Convert an HDF5 compound dataset to a group

Convert an HDF5 compound dataset to a group

## Usage

``` r
CompoundToGroup(src, dst, dname, order, index = NULL, name_map_fn = NULL)
```

## Arguments

- src:

  An HDF5 dataset
  ([`H5D`](http://hhoeflin.github.io/hdf5r/reference/H5D-class.md)) of
  type
  [`H5T_COMPOUND`](http://hhoeflin.github.io/hdf5r/reference/H5T_COMPOUND-class.md)

- dst:

  An HDF5 file
  ([`H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md))
  or group
  ([`H5Group`](http://hhoeflin.github.io/hdf5r/reference/H5Group-class.md))

- dname:

  Name of group in `dst`

- order:

  Column order specification

- index:

  Integer values of which values to pull; defaults to all values

- name_map_fn:

  Optional function to map compound field names

## Value

Invisibly returns `NULL`
