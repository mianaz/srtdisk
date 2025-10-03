# SeuratDisk V5 Fork

## Version 0.1.0 (2024)

### Seurat V5 Compatibility
* Complete support for Seurat V5 Assay5 objects and layer-based data structures
* Native handling of V5 sparse matrix formats in HDF5
* Backward compatibility maintained with Seurat V3/V4 objects
* Support for mixed assay versions (V3/V4/V5) in same object

### AnnData Format Improvements
* Updated to anndata 0.8+ categorical encoding format (categories/codes)
* Proper encoding attributes for all datasets and groups
* Complete h5ad structure with required empty groups (obsm, obsp, varm, varp, layers, uns)
* Correct shape attributes for sparse matrices
* Support for counts-only objects (data layer now optional)

### H5Seurat to H5AD Conversion
* Improved metadata transfer with proper categorical encoding
* Fixed dimension handling for 1D and 2D feature datasets
* Enhanced layer mapping (data → X, counts → raw/X)
* Graph structures properly transferred to obsp
* All dimensional reductions preserved with encoding attributes

### HDF5 Method Patching
* Patched hdf5r H5D methods for improved multidimensional array access
* Automatic fallback to safe read methods
* Better handling of sparse matrix formats

### Bug Fixes
* Fixed "Number of arguments not equal to dimensions" error
* Resolved "Cannot find data slot in V5 h5Seurat file" error
* Fixed categorical variable encoding for anndata compatibility
* Corrected shape attribute missing error for sparse matrices
* Fixed deprecated GetAssayData slot/layer parameter warnings
* Improved metadata preservation during V5 conversions

### New Functions
* `LoadH5AD()`: Direct h5ad to Seurat loading (bypasses h5Seurat)
* `SafeH5DRead()`: Robust HDF5 dataset reading for any dimensionality
* `SafeGetAssayData()`: V3/V4/V5 compatible assay data access
* `SafeGetLayers()`: Helper for layer enumeration across Seurat versions

### Internal Improvements
* Added H5DAccess.R with safe dataset access utilities
* Enhanced V5Compatibility.R with comprehensive helper functions
* Improved error handling and validation
* Better support for environment-wrapped metadata (V5 pattern)

### Testing
* Added comprehensive test suite for save/load cycles
* Multi-assay functionality tests
* Data integrity validation

### Known Limitations
* Round-trip conversions (h5Seurat ↔ Seurat) may lose some V5 metadata
* Multi-modal assays with different feature dimensions have limited support
* Spatial transcriptomics support is basic
* Large datasets (>100K cells) require substantial memory