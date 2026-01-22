# Changelog

## srtdisk 0.1.0

### New Features

#### Direct Seurat to H5AD Conversion

- Added `SeuratToH5AD()` function for direct conversion from Seurat
  objects to H5AD format
- Automatically handles intermediate h5Seurat file creation and cleanup
- Supports all
  [`Convert()`](https://mianaz.github.io/srtdisk/reference/Convert.md)
  parameters including `standardize` for scanpy-compatible naming

#### Improved Spatial Data Support

- Enhanced Visium spatial data conversion between h5Seurat and h5ad
  formats
- Better handling of multi-library spatial datasets
- Improved coordinate transformation and scale factor preservation

#### Seurat v5 Compatibility

- Full support for Seurat v5 Assay5 objects
- Fixed SeuratObject 5.0+ API compatibility (layer vs slot parameters)
- Proper handling of V5 layered data structure in conversions

### Bug Fixes

- Fixed spatial coordinate X/Y flip in h5Seurat to h5ad conversion
  (coordinates are now correctly stored as \[X, Y\] for scanpy/squidpy)
- Fixed “attribute already exists” error when converting h5ad graphs to
  h5seurat
- Fixed `GetAssayDataCompat()` and `SetAssayDataCompat()` to use `layer`
  parameter (SeuratObject 5.0+ requirement)
- Fixed vignette assay references (using correct RNA assay instead of
  non-existent SCT/Spatial)

### Code Improvements

- Removed unused utility functions (`WithAssayCompat`, `GetSlotMapping`,
  `GetSeuratSlotMapping`, `ValidateSlotMapping`)
- Removed deprecated `SafeExistsDeprecated` function
- Simplified verbose messaging in `AssembleObject.R`
- Extracted `CreateFakeCellNames` helper to reduce code duplication in
  `Convert.R`
- Removed dead commented-out code

### Documentation

- Updated vignettes with Python/reticulate integration for
  scanpy/squidpy examples
- Added conditional evaluation for Python chunks based on package
  availability
- Improved vignette examples with bundled test data
- Added cellxgene spatial dataset download examples

### Test Data

- Added `pbmc3k_final.h5ad` test file for round-trip conversion testing
