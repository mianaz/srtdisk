# srtdisk 0.2.1

> **Release Date:** 2026-02-05

## Highlights

### Correct Gene Set in h5Seurat to h5ad Conversion

Fixed a critical issue where `scale.data` (containing only ~2,000 variable features) was prioritized as `X` during h5Seurat to h5ad conversion, causing `var` to contain only variable genes instead of the full gene set. The conversion now correctly uses `data` (normalized, all genes) as `X` and `counts` as `raw/X`, matching standard AnnData conventions. `scale.data` is only used as a last-resort fallback.

| Seurat | h5ad | Notes |
|---|---|---|
| `data` (normalized) | `X` | All genes |
| `counts` (raw) | `raw/X` | All genes |
| `scale.data` | *(skipped)* | Can be recomputed with `sc.pp.scale()` |
| variable features | `var['highly_variable']` | Boolean column marking ~2,000 genes |

### Native AnnData Preprocessing

srtdisk now handles "messy" AnnData files natively without requiring Python preprocessing:

- **Column Name Sanitization**: Automatically sanitizes obs/var column names during h5ad to h5Seurat conversion
  - Replaces problematic characters (`/`, spaces, commas, semicolons, colons, backslashes) with underscores
  - Handles duplicate names by appending `__dup1`, `__dup2`, etc.
  - Ensures valid R names (prefixes with `X` if column starts with a number)
  - Example: `cell/type` -> `cell_type`, `sample name` -> `sample_name`

- **List/Dict Column Conversion**: Converts list and dict columns in compound datasets to R-compatible strings
  - Simple lists converted to semicolon-delimited strings (e.g., `["A", "B", "C"]` -> `"A;B;C"`)
  - Complex nested structures converted to JSON strings
  - Dict columns converted to JSON strings

- **Nullable Dtype Handling**: Properly handles pandas nullable extension dtypes stored as mask+values structures
  - Automatically flattens mask+values groups to simple vectors with NA values
  - Preserves missing value semantics during conversion

## Bug Fixes

- Fixed h5Seurat to h5ad conversion using `scale.data` (2,000 variable features) as `X` instead of `data` (all genes)
- Fixed handling of obs/var column names containing special characters
- Improved robustness of compound dataset reading with Python fallback

---

# srtdisk 0.2.0

> **Release Date:** 2026-01-29

## Highlights

### Full scRNA-seq Conversion Support with Seurat v5

srtdisk now provides complete bidirectional conversion between Seurat objects and h5ad format, with full compatibility for Seurat v5's new Assay5 architecture:

- **Seurat to h5ad**: Use `SeuratToH5AD()` for direct one-step conversion or the traditional two-step `SaveH5Seurat()` + `Convert()` workflow
- **h5ad to Seurat**: Use `Convert()` + `LoadH5Seurat()` to import scanpy/AnnData objects
- Properly handles V5 layered data structure (counts, data, scale.data layers)
- Supports multi-assay objects (e.g., CITE-seq with RNA + ADT)

### Metadata Preservation Fixes

- **Fixed categorical metadata loss**: Factor/categorical variables in `obs` are now correctly preserved during h5ad conversion instead of being dropped or converted to strings
- Improved handling of cell-level and feature-level metadata during round-trip conversions

### Spatial Data Support (Visium)

- **Seurat to h5ad**: Visium spatial data conversion is fully functional
  - Preserves spatial coordinates, scale factors, and tissue images
  - Compatible with scanpy/squidpy spatial analysis workflows
- **h5ad to Seurat**: Spatial conversion fully supported
- Fixed spatial coordinate X/Y orientation for scanpy compatibility

## New Features

- `SeuratToH5AD()` wrapper function for convenient direct conversion
- Enhanced Visium spatial data conversion with proper coordinate handling
- Support for SlideSeq and FOV-based spatial technologies (experimental)
- Improved scanpy-compatible naming conventions with `standardize = TRUE`

## Bug Fixes

- Fixed categorical metadata loss during h5ad conversion
- Fixed "attribute already exists" error when converting h5ad graphs to h5seurat
- Fixed `GetAssayDataCompat()` and `SetAssayDataCompat()` to use `layer` parameter (SeuratObject 5.0+ requirement)
- Fixed spatial coordinate X/Y flip in h5Seurat to h5ad conversion
- Improved UMAP display in h5ad to Seurat tutorial

## Documentation

- Comprehensive README with quick start conversion examples
- Comparison table: srtdisk vs SeuratDisk feature improvements
- Multi-assay (CITE-seq) conversion documentation
- Spatial data conversion examples (beta)

## Known Limitations

- Multi-assay conversion requires separate h5ad files per assay (AnnData limitation)
- Large datasets may require sufficient memory for in-memory conversion

# srtdisk 0.1.0

## New Features

### Direct Seurat to H5AD Conversion
- Added `SeuratToH5AD()` function for direct conversion from Seurat objects to H5AD format
- Automatically handles intermediate h5Seurat file creation and cleanup
- Supports all `Convert()` parameters including `standardize` for scanpy-compatible naming

### Improved Spatial Data Support
- Enhanced Visium spatial data conversion between h5Seurat and h5ad formats
- Better handling of multi-library spatial datasets
- Improved coordinate transformation and scale factor preservation

### Seurat v5 Compatibility
- Full support for Seurat v5 Assay5 objects
- Fixed SeuratObject 5.0+ API compatibility (layer vs slot parameters)
- Proper handling of V5 layered data structure in conversions

## Bug Fixes

- Fixed spatial coordinate X/Y flip in h5Seurat to h5ad conversion (coordinates are now correctly stored as [X, Y] for scanpy/squidpy)
- Fixed "attribute already exists" error when converting h5ad graphs to h5seurat
- Fixed `GetAssayDataCompat()` and `SetAssayDataCompat()` to use `layer` parameter (SeuratObject 5.0+ requirement)
- Fixed vignette assay references (using correct RNA assay instead of non-existent SCT/Spatial)

## Code Improvements

- Removed unused utility functions (`WithAssayCompat`, `GetSlotMapping`, `GetSeuratSlotMapping`, `ValidateSlotMapping`)
- Removed deprecated `SafeExistsDeprecated` function
- Simplified verbose messaging in `AssembleObject.R`
- Extracted `CreateFakeCellNames` helper to reduce code duplication in `Convert.R`
- Removed dead commented-out code

## Documentation

- Updated vignettes with Python/reticulate integration for scanpy/squidpy examples
- Added conditional evaluation for Python chunks based on package availability
- Improved vignette examples with bundled test data
- Added cellxgene spatial dataset download examples

## Test Infrastructure

- Simplified test suite: consolidated 4 test files into single focused `test-conversion.R`
- Tests now use real datasets:
  - CellxGene colorectal cancer sample (935 cells) for h5ad testing
  - `pbmc3k.final` from SeuratData for Seurat testing
  - `stxBrain` anterior1 from SeuratData for Visium spatial testing
- Removed synthetic test data in favor of real-world datasets
- Added `UpdateSeuratObject()` calls for compatibility with older SeuratData objects

## Test Data

- Replaced `pbmc_small.rds` with `crc_sample.h5ad` (CellxGene colorectal cancer, 935 cells x 25,344 genes)
- Vignettes now use SeuratData (`pbmc3k.final`, `stxBrain`) instead of bundled files
