# Write a SpatialImage object to an HDF5 dataset

Write a SpatialImage object to an HDF5 dataset

## Usage

``` r
ImageWrite(x, name, hgroup, verbose = TRUE)
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
