# Seurat bindings for h5Seurat files

Seurat bindings for h5Seurat files

## Usage

``` r
# S3 method for class 'h5Seurat'
Cells(x, ...)

# S3 method for class 'h5Seurat'
DefaultAssay(object, ...)

# S3 method for class 'h5Seurat'
DefaultAssay(object, ...) <- value

# S3 method for class 'h5Seurat'
Idents(object, ...)

# S3 method for class 'H5Group'
IsGlobal(object, ...)

# S3 method for class 'H5Group'
Key(object, ...)

# S3 method for class 'h5Seurat'
Project(object, ...)

# S3 method for class 'h5Seurat'
Project(object, ...) <- value

# S3 method for class 'h5Seurat'
Stdev(object, reduction = "pca", ...)
```

## Arguments

- x:

  An h5Seurat object

- object:

  An h5Seurat or H5Group object

- reduction:

  Name of the dimensional reduction (default: "pca")

- value:

  Value to set

- ...:

  Additional arguments passed to methods
