# Convert an on-disk single-cell dataset to another format

HDF5-based single-cell datasets can be converted from one format to
another using minimal memory. Details about conversion formats
implemented are provided below

## Usage

``` r
Convert(source, dest, assay, overwrite = FALSE, verbose = TRUE, standardize = FALSE, ...)

# S3 method for class 'character'
Convert(source, dest, assay, overwrite = FALSE, verbose = TRUE, standardize = FALSE, ...)

# S3 method for class 'H5File'
Convert(
  source,
  dest = "h5seurat",
  assay = "RNA",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'h5Seurat'
Convert(
  source,
  dest = "h5ad",
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
)
```

## Arguments

- source:

  Source dataset

- dest:

  Name of destination dataset

- assay:

  Converting from
  [`h5Seurat`](https://mianaz.github.io/srtdisk/reference/h5Seurat-class.md):
  name of assay to write out; converting to
  [`h5Seurat`](https://mianaz.github.io/srtdisk/reference/h5Seurat-class.md):
  name to store assay data as

- overwrite:

  Overwrite existing `dest`

- verbose:

  Show progress updates

- standardize:

  Logical; if `TRUE`, convert Seurat-style metadata column names to
  scanpy/AnnData conventions when converting to h5ad format. For
  example, `nCount_RNA` becomes `n_counts`, `nFeature_RNA` becomes
  `n_genes`. Only applicable for h5Seurat to h5ad conversions. Default
  is `FALSE`.

- ...:

  Arguments passed to other methods

## Value

If `source` is a `character`, invisibly returns `dest`; otherwise,
returns an
[`H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md),
or filetype-specific subclass of `H5File` (eg.
[`h5Seurat`](https://mianaz.github.io/srtdisk/reference/h5Seurat-class.md)),
connection to `dest`

## AnnData/H5AD to h5Seurat

The AnnData/H5AD to h5Seurat conversion will try to automatically fill
in datasets based on data presence. It works in the following manner:

### Expression data

The expression matrices `counts`, `data`, and `scale.data` are filled by
`/X` and `/raw/X` in the following manner:

- `counts` will be filled with `/raw/X` if present; otherwise, it will
  be filled with `/X`

- `data` will be filled with `/raw/X` if `/raw/X` is present and `/X` is
  dense; otherwise, it will be filled with `/X`

- `scale.data` will be filled with `/X` if it dense; otherwise, it will
  be empty

Feature names are taken from the feature-level metadata

### Feature-level metadata

Feature-level metadata is added to the `meta.features` datasets in each
assay. Feature names are taken from the dataset specified by the
“\_index” attribute, the “\_index” dataset, or the “index” dataset, in
that order. Metadata is populated with `/raw/var` if present, otherwise
with `/var`; if both `/raw/var` and `/var` are present, then
`meta.features` will be populated with `/raw/var` first, then `/var`
will be added to it. For columns present in both `/raw/var` and `/var`,
the values in `/var` will be used instead. **Note**: it is possible for
`/var` to have fewer features than `/raw/var`; if this is the case, then
only the features present in `/var` will be overwritten, with the
metadata for features *not* present in `/var` remaining as they were in
`/raw/var` or empty

### Cell-level metadata

Cell-level metadata is added to `meta.data`; the row names of the
metadata (as determined by the value of the “\_index” attribute, the
“\_index” dataset, or the “index” dataset, in that order) are added to
the “cell.names” dataset instead. If the “\_\_categories” dataset is
present, each dataset within “\_\_categories” will be stored as a factor
group. Cell-level metadata will be added as an HDF5 group unless factors
are **not** present and `srtdisk.dtype.dataframe_as_group` is `FALSE`

### Dimensional reduction information:

Cell embeddings are taken from `/obsm`; dimensional reductions are named
based on their names from `obsm` by removing the preceding “X\_”.For
example, if a dimensional reduction is named “X_pca” in `/obsm`, the
resulting dimensional reduction information will be named “pca”. The key
will be set to one of the following:

- “PC\_” if “pca” is present in the dimensional reduction name
  (`grepl("pca", reduction.name, ignore.case = TRUE)`)

- “tSNE\_” if “tsne” is present in the dimensional reduction name
  (`grepl("tsne", reduction.name, ignore.case = TRUE)`)

- `reduction.name_` for all other reductions

Remember that the preceding “X\_” will be removed from the reduction
name before converting to a key. Feature loadings are taken from `/varm`
and placed in the associated dimensional reduction. The dimensional
reduction is determine from the loadings name in `/varm`:

- “PCs” will be added to a dimensional reduction named “pca”

- All other loadings in `/varm` will be added to a dimensional reduction
  named `tolower(loading)` (eg. a loading named “ICA” will be added to a
  dimensional reduction named “ica”)

If a dimensional reduction cannot be found according to the rules above,
the loading will not be taken from the AnnData/H5AD file. Miscellaneous
information will be taken from `/uns/reduction` where `reduction` is the
name of the reduction in `/obsm` without the preceding “X\_”; if no
dimensional reduction information present, then miscellaneous
information will not be taken from the AnnData/H5AD file. Standard
deviations are taken from a dataset `/uns/reduction/variance`; the
variances will be converted to standard deviations and added to the
`stdev` dataset of a dimensional reduction

### Nearest-neighbor graph

If a nearest neighbor graph is present in `/uns/neighbors/distances`, it
will be added as a graph dataset in the h5Seurat file and associated
with `assay`; if a value is present in `/uns/neighbors/params/method`,
the name of the graph will be `assay_method`, otherwise, it will be
`assay_anndata`

### Layers

TODO: add this

### Miscellaneous information

All groups and datasets from `/uns` will be copied to `misc` in the
h5Seurat file except for the following:

- Any group or dataset named the same as a dimensional reduction (eg.
  `/uns/pca`)

- `/uns/neighbors`

## h5Seurat to AnnData/H5AD

The h5Seurat to AnnData/H5AD conversion will try to automatically fill
in datasets based on data presence. Data presense is determined by the
h5Seurat index (`source$index()`). It works in the following manner:

### Assay data

- `X` will be filled with `scale.data` if `scale.data` is present;
  otherwise, it will be filled with `data`

- `var` will be filled with `meta.features` **only** for the features
  present in `X`; for example, if `X` is filled with `scale.data`, then
  `var` will contain only features that have been scaled

- `raw.X` will be filled with `data` if `X` is filled with `scale.data`;
  otherwise, it will be filled with `counts`. If `counts` is not
  present, then `raw` will not be filled

- `raw.var` will be filled with `meta.features` with the features
  present in `raw.X`; if `raw.X` is not filled, then `raw.var` will not
  be filled

### Cell-level metadata

Cell-level metadata is added to `obs`

### Dimensional reduction information

Only dimensional reductions associated with `assay` or marked as
[global](https://satijalab.github.io/seurat-object/reference/IsGlobal.html)
will be transfered to the H5AD file. For every reduction `reduc`:

- cell embeddings are placed in `obsm` and renamed to `X_reduc`

- feature loadings, if present, are placed in `varm` and renamed to
  either “PCs” if `reduc` is “pca” otherwise `reduc` in all caps

For example, if `reduc` is “ica”, then cell embeddings will be “X_ica”
in `obsm` and feature loaodings, if present, will be “ICA” in `varm`

### Nearest-neighbor graphs

If a nearest-neighbor graph is associated with `assay`, it will be added
to `uns/neighbors/distances`; if more than one graph is present, then
**only** the last graph according to the index will be added.

### Layers

Data from other assays can be added to `layers` if they have the same
shape as `X` (same number of cells and features). To determine this, the
shape of each alternate assays's `scale.data` and `data` slots are
determined. If they are the same shape as `X`, then that slot
(`scale.data` is given priority over `data`) will be added as a layer
named the name of the assay (eg. “SCT”). In addition, the features names
will be added to `var` as `assay_features` (eg. “SCT_features”).
