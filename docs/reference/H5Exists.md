# Check to see if a dataset, group, or attribute exists in an HDF5 file, group, or dataset

Check to see if a dataset, group, or attribute exists in an HDF5 file,
group, or dataset

## Usage

``` r
AttrExists(x, name)

Exists(x, name)
```

## Arguments

- x:

  An HDF5
  [file](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md) or
  [group](http://hhoeflin.github.io/hdf5r/reference/H5Group-class.md);
  for `AttrExists`, may also be a
  [dataset](http://hhoeflin.github.io/hdf5r/reference/H5D-class.md)

- name:

  Name of dataset, group, or attribute to test for

## Value

`TRUE` if `name` exists in `x`, otherwise `FALSE`
