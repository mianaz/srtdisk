
<!-- README.md is generated from README.Rmd. Please edit that file -->

# srtdisk v0.1.0

<!-- badges: start -->

[![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://github.com/miana/seuratdisk-V5)
[![R](https://img.shields.io/badge/R-%3E%3D3.6-blue.svg)]()
[![Seurat](https://img.shields.io/badge/Seurat-v5-green.svg)]()
<!-- badges: end -->

## Overview

Extended implementation of SeuratDisk for Seurat v5 with enhanced
support for HDF5-based single cell file formats. The h5Seurat file
format is specifically designed for the storage and analysis of
multi-modal single-cell and spatially-resolved expression experiments,
for example, from CITE-seq or 10X Visium technologies. It holds all
molecular information and associated metadata, including (for example)
nearest-neighbor graphs, dimensional reduction information, spatial
coordinates and image data, and cluster labels.

We also support rapid and on-disk conversion between h5Seurat and
AnnData objects, with the goal of enhancing interoperability between
Seurat and Scanpy. This version includes **Seurat v5 Assay5
compatibility** and **improved spatial data handling**.

## Key Features

- **Seurat V5 Assay5 Support**: Full compatibility with Seurat v5
  layer-based data structure
- **h5Seurat ↔ h5ad Conversion**: Seamless conversion between Seurat and
  AnnData formats
- **Visium Spatial Data**: Complete support for 10X Visium V1/V2 spatial
  transcriptomics
- **On-Disk Operations**: Memory-efficient file format conversion
  without loading full datasets
- **Spatial Image Reconstruction**: Automatic reconstruction of Visium
  images from h5ad files

## Installation

srtdisk is not currently available on CRAN. You can install it from
GitHub with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("miana/seuratdisk-v5")
```

## Quick Start

### Seurat to AnnData (h5ad)

``` r
library(Seurat)
library(srtdisk) # changed package name

# Save Seurat object as h5Seurat
SaveH5Seurat(seurat_obj, filename = "data.h5Seurat")

# Convert to h5ad for Python/Scanpy
Convert("data.h5Seurat", dest = "h5ad")
```

### AnnData to Seurat

``` r
# Convert h5ad to h5Seurat
Convert("data.h5ad", dest = "h5seurat")

# Load as Seurat object
seurat_obj <- LoadH5Seurat("data.h5seurat")
```

### Spatial Visium Data

``` r
# Load Visium data from h5ad (created with Scanpy/Squidpy)
brain <- Convert("visium_brain.h5ad", dest = "h5seurat")
brain <- LoadH5Seurat("visium_brain.h5seurat")

# Spatial images are automatically reconstructed
Images(brain)
SpatialFeaturePlot(brain, features = "PTPRC")
```

## Supported Conversions

| From          | To            | Support        |
|---------------|---------------|----------------|
| Seurat Object | h5Seurat      | ✓ Full         |
| h5Seurat      | Seurat Object | ✓ Full         |
| h5Seurat      | h5ad          | ✓ Full         |
| h5ad          | h5Seurat      | ✓ Full         |
| h5ad          | Seurat Object | ✓ Via h5Seurat |
| Seurat Object | h5ad          | ✓ Via h5Seurat |

## Seurat V5 Compatibility

This package fully supports Seurat v5 Assay5 objects:

- **Layer-based data storage**: Automatic handling of
  counts/data/scale.data layers
- **Assay5 conversion**: Seamless conversion when saving to h5Seurat
  format
- **Backward compatibility**: Works with both legacy Assay and modern
  Assay5 objects
- **Spatial data**: Enhanced support for VisiumV1 and VisiumV2 formats

## Documentation

See the vignettes for detailed usage:

- **[convert-anndata.Rmd](vignettes/convert-anndata.Rmd)**: h5Seurat ↔
  h5ad conversion workflows
- **[h5Seurat-load.Rmd](vignettes/h5Seurat-load.Rmd)**: Saving and
  loading h5Seurat files
- **[h5Seurat-spec.Rmd](vignettes/h5Seurat-spec.Rmd)**: h5Seurat file
  format specification

## System Requirements

- R \>= 3.6.0
- Seurat \>= 3.2.0
- SeuratObject \>= 4.0.0
- hdf5r \>= 1.3.0
