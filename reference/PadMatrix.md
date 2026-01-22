# Pad a matrix

Pad a matrix

## Usage

``` r
PadMatrix(src, dest, dname, dims, index)

# S4 method for class 'H5D'
PadMatrix(src, dest, dname, dims, index)

# S4 method for class 'H5Group'
PadMatrix(src, dest, dname, dims, index)
```

## Arguments

- src:

  A source matrix

- dest:

  Destination HDF5 file or group for the padded matrix

- dname:

  Destination name for the padded matrix

- dims:

  A two-length integer vector with the number of rows and number of
  columns in the padded matrix

- index:

  A two-length list of integer vectors describing the rows and columns
  that `src` exists in

## Value

...
