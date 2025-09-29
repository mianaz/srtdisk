# SeuratDisk (Community V5 Fork)

<!-- badges: start -->
[![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)]()
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)]()
[![Seurat](https://img.shields.io/badge/Seurat-v5-green.svg)]()
<!-- badges: end -->

## Overview

SeuratDisk provides interfaces for HDF5-based single cell file formats. This is a community-maintained fork with Seurat V5 compatibility. The h5Seurat format enables efficient storage of large single-cell datasets and conversion between Seurat and AnnData/h5ad formats for interoperability with Scanpy.

## Features

- **Seurat V5 Compatible**: Full support for Seurat V5 objects and assay structures
- **Efficient Storage**: HDF5-based format for large single-cell datasets
- **Format Conversion**: Bidirectional conversion between Seurat, h5Seurat, and h5ad formats
- **Metadata Preservation**: Complete transfer of cell metadata, embeddings, and analysis results
- **Direct h5ad Loading**: Load h5ad files directly into Seurat (bypassing h5Seurat intermediate)

## Installation

```r
# Install from this repository
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# Replace 'username' with the actual GitHub username hosting this fork
remotes::install_github("username/seuratdisk-V5")

# Or install locally
R CMD INSTALL .
```

## Quick Start

### Save and Load Seurat Objects

```r
library(SeuratDisk)

# Save Seurat object to h5Seurat format
SaveH5Seurat(seurat_obj, filename = "data.h5Seurat")

# Load h5Seurat file back to Seurat
seurat_obj <- LoadH5Seurat("data.h5Seurat")
```

### Convert Between Formats

```r
# Convert Seurat to h5ad (for Scanpy/Python)
Convert("data.h5Seurat", dest = "data.h5ad")

# Convert h5ad to h5Seurat
Convert("data.h5ad", dest = "data.h5Seurat")

# Load h5ad directly into Seurat (recommended for V5)
seurat_obj <- LoadH5AD("data.h5ad")
```

### Working with AnnData Files

```r
# Direct h5ad to Seurat conversion (bypasses h5Seurat)
seurat_obj <- LoadH5AD("scanpy_output.h5ad")

# The LoadH5AD function preserves:
# - Expression matrices (X, layers)
# - Cell metadata (obs)
# - Gene metadata (var)
# - Dimensional reductions (obsm)
# - Graphs (obsp)
# - Unstructured data (uns)
```

## System Requirements

- **R**: >= 4.0.0
- **Seurat**: >= 5.0.0
- **SeuratObject**: >= 5.0.0
- **hdf5r**: >= 1.3.0

### Dependencies

Core dependencies are automatically installed:
- **Matrix**: Sparse matrix operations
- **hdf5r**: HDF5 file interface
- **R6**: Object system
- Additional utility packages (cli, crayon, rlang, stringi, withr)

## Supported Features

✅ **Working**:
- Seurat V5 → h5Seurat → h5ad conversion
- Direct h5ad → Seurat loading (LoadH5AD)
- All cell metadata preservation (V5 fix)
- Multiple assays (RNA, SCT)
- Dimensional reductions (PCA, UMAP, tSNE)
- Highly variable genes
- Graph structures

⚠️ **Limited Support**:
- Multi-modal assays with different dimensions (ADT, HTO)
- Spatial transcriptomics data
- Complex nested metadata structures

## Known Limitations

1. **Round-trip conversions**: h5Seurat → Seurat → h5Seurat may encounter issues with V5 objects
2. **Multi-modal data**: Assays with different feature dimensions have limited support
3. **Large files**: Datasets >100K cells may require increased memory
4. **Spatial data**: Basic structure preserved but not fully integrated

## Troubleshooting

### HDF5 Errors
If you encounter HDF5 errors, ensure you're using compatible library versions:
```r
# Check hdf5r version
packageVersion("hdf5r")  # Should be >= 1.3.0
```

### Memory Issues
For large datasets:
```r
# Increase memory limit (Windows)
memory.limit(size = 16000)

# Monitor memory usage
gc()
```

### V5 Compatibility
For best V5 compatibility, use `LoadH5AD()` for direct h5ad loading rather than the h5Seurat intermediate format.

## Contributing

This is a community-maintained fork. Contributions are welcome! Key areas for improvement:
- Enhanced multi-modal support
- Spatial data integration
- Performance optimization for large datasets

## Original Project

This fork builds upon the original SeuratDisk package. The original repository by mojaveazure may no longer be actively maintained.

## Citation

If you use SeuratDisk in your research, please cite the Seurat papers:

```
Hao and Hao et al. Dictionary learning for integrative, multimodal and scalable
single-cell analysis. Nature Biotechnology (2024) [Seurat V5]

Hao and Hao et al. Integrated analysis of multimodal single-cell data.
Cell (2021) [Seurat V4]
```

## License

GPL-3

## Support

For issues specific to this V5-compatible fork, please open an issue in this repository. For general Seurat questions, visit the [Seurat community](https://github.com/satijalab/seurat/discussions).