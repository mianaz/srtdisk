# Connect to a single-cell HDF5 dataset

Connect to a single-cell HDF5 dataset

## Usage

``` r
Connect(filename, type = NULL, mode = c("r", "r+"), force = FALSE)
```

## Arguments

- filename:

  Name of on-disk file

- type:

  Type of single-cell dataset to connect as; choose from:

  - h5seurat

  Leave as `NULL` to guess type from file extension

- mode:

  Mode to connect to data as; choose from:

  r

  :   Open existing dataset in read-only mode

  r+

  :   Open existing dataset in read/write mode

- force:

  Force a connection if validation steps fail; returns a
  [`H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md)
  object

## Value

An object of class `type`, opened in mode `mode`
