# srtdisk: Interfaces for HDF5-Based Single Cell File Formats

The h5Seurat file format is specifically designed for the storage and
analysis of multi-modal single-cell and spatially-resolved expression
experiments, for example, from CITE-seq or 10X Visium technologies. It
holds all molecular information and associated metadata, including (for
example) nearest-neighbor graphs, dimensional reduction information,
spatial coordinates and image data, and cluster labels. We also support
rapid and on-disk conversion between h5Seurat and AnnData objects, with
the goal of enhancing interoperability between Seurat and Scanpy.

## Package options

srtdisk uses the following options to control behavior, users can
configure these with [`options`](https://rdrr.io/r/base/options.html):

- `srtdisk.dtypes.logical_to_int`:

  When writing [logical](https://rdrr.io/r/base/logical.html) vectors,
  coerce to integer types to ensure compatibility across languages (see
  [`BoolToInt`](https://mianaz.github.io/srtdisk/reference/BoolToInt.md)
  for more details)

- `srtdisk.dtypes.dataframe_as_group`:

  When writing [data.frame](https://rdrr.io/r/base/data.frame.html)s,
  always write out as a group regardless of factor presence

- `srtdisk.chunking.MARGIN`:

  Default direction for chunking datasets; choose from:

  largest

  :   Chunk along the largest dimension of a dataset

  smallest

  :   Chunk along the smallest dimension

  first

  :   Chunk along the first dimension

  last

  :   Chunk along the last dimension

- `srtdisk.dimreducs.allglobal`:

  Treat all DimReducs as global, regardless of actual global status

## See also

Useful links:

- <https://mojaveazure.github.io/seurat-disk/>

- <https://github.com/mojaveazure/seurat-disk>

- Report bugs at <https://github.com/mojaveazure/seurat-disk/issues>

## Author

**Maintainer**: Paul Hoffman <phoffman@nygenome.org>
([ORCID](https://orcid.org/0000-0002-7693-8957))

Other contributors:

- Rahul Satija <rsatija@nygenome.org>
  ([ORCID](https://orcid.org/0000-0001-9448-8833)) \[contributor\]
