# CLAUDE.md - SeuratDisk Fork Analysis & Development Guide

## Project Overview

**Package**: srtdisk (fork of SeuratDisk)
**Author**: Miana Zeng
**Version**: 0.1.0
**Purpose**: Extended HDF5-based single cell file format support with Seurat V5 compatibility

## Current Architecture Analysis

### Core Components

#### 1. File Format Hierarchy
```
h5Seurat (central format)
    ├── Conversion Layer
    │   ├── h5ad (AnnData/Python)
    │   ├── Loom (limited support)
    │   └── h5mu (in development)
    ├── I/O Layer
    │   ├── SaveH5Seurat
    │   ├── LoadH5Seurat
    │   └── Convert
    └── Compatibility Layer
        ├── V3/V4 Seurat
        ├── V5 Seurat (Assay5)
        └── Spatial (Visium V1/V2)
```

#### 2. Key Function Categories

**Core Classes** (`R/h5Seurat.R`, `R/scdisk.R`)
- `h5Seurat`: R6 class for h5Seurat file connections
- `scdisk`: Base class for single-cell disk formats
- Version detection and management

**Conversion Functions** (`R/Convert.R`)
- `Convert()`: Main conversion dispatcher
- `H5ADToH5Seurat()`: AnnData → h5Seurat
- `H5SeuratToH5AD()`: h5Seurat → AnnData
- Memory-efficient on-disk operations

**V5 Compatibility** (`R/V5Compatibility.R`, `R/V5LayerSupport.R`)
- Layer-based data structure support
- Assay5 object handling
- BPCells matrix conversion (basic)

**Spatial Support** (`R/SpatialConversion.R`, `R/SpatialFOV.R`)
- Visium image reconstruction
- Spatial coordinate handling
- Scale factor management

**Utilities**
- Matrix operations (`R/UtilsMatrix.R`, `R/UtilsSparseMatrix.R`)
- HDF5 access (`R/UtilsH5Access.R`)
- Assay compatibility (`R/UtilsAssayCompat.R`)

### Dependencies Analysis

**Core Dependencies**:
- `hdf5r` (≥1.3.0): HDF5 file operations
- `Seurat` (≥3.2.0): Seurat object compatibility
- `SeuratObject` (≥4.0.0): Object structures
- `Matrix` (≥1.2.18): Sparse matrix support

**Utility Dependencies**:
- `cli`, `crayon`: Terminal output
- `R6`: Object-oriented programming
- `rlang`, `withr`: R programming utilities

**Not Using**:
- `anndataR`: Not integrated
- `MuDataSeurat`: Not integrated
- `BPCells`: Basic conversion only, no full integration

## Architecture Questions & Answers

### Q1: h5Seurat Intermediate Format Analysis

**Memory Advantages**: ✅ Yes
- On-disk operations without loading full datasets
- Partial data loading capability via HDF5
- Chunked reading/writing for large matrices

**Compression Support**: ✅ Yes
- HDF5 native compression (gzip levels 0-9)
- Chunk-based compression for sparse matrices
- Efficient storage of repetitive data

**Performance**: ⚠️ Good but improvable
- Fast for datasets <100k cells
- Some bottlenecks in very large datasets (>500k cells)
- Conversion typically takes 2-10 minutes for standard datasets

**Data Preservation**: ⚠️ Mostly complete
- ✅ Expression matrices (counts, data, scale.data)
- ✅ Metadata (cell and feature annotations)
- ✅ Reductions (PCA, UMAP, tSNE)
- ✅ Graphs (SNN, KNN)
- ✅ Spatial information (images, coordinates)
- ⚠️ Some custom slots may be lost
- ⚠️ Complex nested lists need special handling

### Q2: Bypassing h5Seurat Intermediate

**Recommendation**: Keep h5Seurat-centric approach

**Reasons**:
1. **Established architecture**: All conversions flow through h5Seurat
2. **Consistency**: Single source of truth for data structure
3. **Maintenance**: Easier to debug and maintain one central format
4. **Flexibility**: Can add new formats without rewriting everything

**Proposed Enhancement**:
```r
# Wrapper for direct conversion (internally uses h5Seurat)
SaveH5AD <- function(object, filename, cleanup = TRUE) {
  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  SaveH5Seurat(object, temp_h5seurat)
  Convert(temp_h5seurat, dest = filename)
  if (cleanup) unlink(temp_h5seurat)
  invisible(filename)
}
```

### Q3: External Package Recommendations

**Current Approach**: h5Seurat-centric ✅ Recommended to maintain

**anndataR**: ❌ Not recommended
- Would duplicate existing conversion logic
- Different approach (in-memory R AnnData objects)
- Additional dependency burden

**MuDataSeurat**: ⚠️ Consider selective integration
- Could help with h5mu format
- But may conflict with current architecture
- Suggest: Study their h5mu implementation, adapt key concepts

**Strategy**: Build native h5mu support using existing infrastructure

### Q4: BPCells Support Status

**Current State**: Basic conversion only
- `ConvertBPCellsMatrix()` in `R/zzz.R`
- Converts BPCells to dgCMatrix for compatibility
- No on-disk operations preserved

**Implementation Roadmap**:
```r
# Phase 1: Detect and preserve BPCells objects
IsBPCellsObject <- function(x) {
  inherits(x, c("IterableMatrix", "RenameDims"))
}

# Phase 2: Write BPCells-aware save function
SaveH5SeuratWithBPCells <- function(object, filename, preserve_on_disk = TRUE) {
  # Detect BPCells layers
  # Save references instead of materializing
  # Store BPCells metadata
}

# Phase 3: Load with BPCells reconstruction
LoadH5SeuratAsBPCells <- function(filename, layers = NULL) {
  # Reconstruct BPCells objects
  # Maintain on-disk properties
}
```

## Code Quality Analysis

### Function Complexity
Most functions are well-sized (<200 lines), but some need refactoring:
- `SaveH5Seurat.Seurat`: Could be modularized
- `Convert.h5Seurat`: Multiple responsibility areas
- `LoadH5Seurat`: Complex version handling

### Test Coverage
**Current**: Minimal (2 test files)
- `test-pbmc-small.R`: Basic roundtrip test
- `test-h5ad-conversion.R`: Conversion tests

**Needed**:
- Spatial data tests
- V5 compatibility tests
- Large dataset performance tests
- Edge case handling

### Documentation
- Good function documentation
- Vignettes present but need updates
- Missing: Architecture documentation, contribution guide

## Development Roadmap - Next Week Focus

### Priority 1: Complete h5mu/MuData Support (Days 1-2)
```r
# Target implementation
SaveH5MU <- function(object, filename, assays = NULL) {
  # Direct Seurat to h5mu
  # Handle multiple modalities
  # Preserve relationships
}

LoadH5MU <- function(filename, assays = NULL) {
  # Load as multimodal Seurat
  # Reconstruct assay relationships
}
```

### Priority 2: Enhance BPCells Integration (Days 3-4)
- Implement detection functions
- Add on-disk preservation options
- Create BPCells-aware save/load variants

### Priority 3: Performance Optimization (Day 5)
- Profile conversion bottlenecks
- Implement parallel processing where possible
- Optimize memory usage for large matrices

### Priority 4: Test Suite Expansion (Ongoing)
- Add multimodal tests
- Spatial data roundtrip tests
- Performance benchmarks
- CI/CD improvements

## File Organization Plan

```
srtdisk/
├── R/                      # Source code
│   ├── core/              # Core classes (future refactor)
│   ├── converters/        # Conversion functions (future refactor)
│   └── utils/            # Utility functions (future refactor)
├── tests/
│   ├── testthat/
│   └── benchmark/         # Performance tests (new)
├── vignettes/
├── devs/                  # Development notes (new, gitignored)
│   ├── performance/       # Profiling results
│   ├── debug/            # Debug outputs
│   └── notes/            # Development notes
└── inst/
    └── testdata/         # Test datasets
```

## Immediate Next Steps

1. **Today**: Create devs directory structure
2. **Today**: Document current function call hierarchy
3. **Tomorrow**: Begin h5mu implementation
4. **This Week**: Complete at least one major enhancement
5. **End of Week**: Update README with new features

## Performance Targets

- Conversion speed: <5 minutes for 100k cells
- Memory usage: <2x dataset size
- File size: Within 20% of source format
- Data integrity: 100% preservation of standard slots

## Known Issues to Address

1. Some V5 layer types not fully supported
2. Large spatial images may cause memory issues
3. Factor levels sometimes incorrectly indexed in h5ad
4. Missing progress bars for some long operations
5. Need better error messages for common failures

## Success Metrics

- [ ] All standard Seurat workflows supported
- [ ] Seamless Python interoperability
- [ ] Performance suitable for Atlas-scale data
- [ ] Comprehensive test coverage (>80%)
- [ ] Clear documentation and examples

## Notes for Miana

This fork has made excellent progress on V5 compatibility and spatial support. The architecture is sound but needs enhancement for modern workflows (BPCells, multimodal). The h5Seurat-centric approach is the right choice - it provides consistency and maintainability.

**Recommended focus**: Complete h5mu support first (highest user demand), then optimize for large-scale datasets with BPCells. Keep the codebase clean and avoid over-engineering - the current simplicity is a strength.

---
*Document created: 2025-10-31*
*Last updated: 2025-10-31*
*Analysis by: Claude*