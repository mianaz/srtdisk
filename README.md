
<!-- README.md is generated from README.Rmd. Please edit that file -->

# srtdisk v0.2.0

<!-- badges: start -->

[![CRAN/METACRAN](https://img.shields.io/cran/v/srtdisk)](https://cran.r-project.org/package=srtdisk)
[![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://github.com/mianaz/srtdisk)
<!-- badges: end -->

<!-- Interfaces for HDF5-based Single Cell File Formats -->

Extended implementation of SeuratDisk for Seurat v5 with enhanced
support for HDF5-based single cell file formats. The h5Seurat file
format is specifically designed for the storage and analysis of
multi-modal single-cell and spatially-resolved expression experiments,
for example, from CITE-seq or 10X Visium technologies. It holds all
molecular information and associated metadata, including (for example)
nearest-neighbor graphs, dimensional reduction information, spatial
coordinates and image data, and cluster labels. We also support rapid
and on-disk conversion between h5Seurat and AnnData objects, with the
goal of enhancing interoperability between Seurat and Scanpy. This
version includes Seurat v5 Assay5 compatibility and improved spatial
data handling.

## What’s New in srtdisk (vs SeuratDisk)

| Feature              | SeuratDisk             | srtdisk                     |
|----------------------|------------------------|-----------------------------|
| Seurat v5 Assay5     | Not supported          | Full support                |
| Spatial Technologies | Visium only            | Visium, SlideSeq, FOV       |
| Categorical Metadata | Lost during conversion | Preserved                   |
| Direct Conversion    | Two-step only          | `SeuratToH5AD()` wrapper    |
| h5ad Compatibility   | Basic                  | Improved scanpy conventions |

## Quick Start

### Seurat to h5ad (Two-step)

``` r
library(Seurat)
library(srtdisk)
library(SeuratData)

data("pbmc3k.final", package = "pbmc3k.SeuratData")
pbmc <- UpdateSeuratObject(pbmc3k.final)

# Step 1: Save as h5Seurat
SaveH5Seurat(pbmc, filename = "pbmc3k.h5Seurat", overwrite = TRUE)

# Step 2: Convert to h5ad
Convert("pbmc3k.h5Seurat", dest = "h5ad", overwrite = TRUE)
```

### Seurat to h5ad (One-step)

``` r
# Direct conversion using wrapper function
SeuratToH5AD(pbmc, filename = "pbmc3k_direct.h5ad", overwrite = TRUE)
```

### h5ad to Seurat

``` r
# Using bundled CRC sample
h5ad_path <- system.file("testdata", "crc_sample.h5ad", package = "srtdisk")
Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
crc <- LoadH5Seurat("crc_sample.h5seurat")
```

### Multi-assay (CITE-seq)

``` r
# Note: Converts one assay at a time
library(SeuratData)
InstallData("cbmc")
data("cbmc", package = "cbmc.SeuratData")
SeuratToH5AD(cbmc, filename = "cbmc_rna.h5ad", assay = "RNA", overwrite = TRUE)
SeuratToH5AD(cbmc, filename = "cbmc_adt.h5ad", assay = "ADT", overwrite = TRUE)
```

### Spatial Data Conversion (Beta)

``` r
library(SeuratData)
InstallData("stxBrain")
brain <- UpdateSeuratObject(LoadData("stxBrain", type = "anterior1"))
SeuratToH5AD(brain, filename = "brain_spatial.h5ad", overwrite = TRUE)
```

## Installation

srtdisk is not currently available on CRAN. You can install it from
[GitHub](https://github.com/mianaz/srtdisk) with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mianaz/srtdisk")
```

## Dependencies

srtdisk depends on the following non-standard packages:

| Package | CRAN Webpage | Source | Website |
|:--:|:--:|:--:|:--:|
| cli | [CRAN](https://cran.r-project.org/package=cli) | [GitHub](https://github.com/r-lib/cli) | [Website](https://cli.r-lib.org) |
| crayon | [CRAN](https://cran.r-project.org/package=crayon) | [GitHub](https://github.com/r-lib/crayon) | [Website](https://r-lib.github.io/crayon/) |
| hdf5r | [CRAN](https://cran.r-project.org/package=hdf5r) | [GitHub](https://github.com/hhoeflin/hdf5r/) | [Website](https://hhoeflin.github.io/hdf5r/) |
| Matrix | [CRAN](https://cran.r-project.org/package=Matrix) | – | [Website](https://Matrix.R-forge.R-project.org) |
| R6 | [CRAN](https://cran.r-project.org/package=R6) | [GitHub](https://github.com/r-lib/R6) | [Website](https://r6.r-lib.org) |
| rlang | [CRAN](https://cran.r-project.org/package=rlang) | [GitHub](https://github.com/r-lib/rlang) | [Website](https://rlang.r-lib.org) |
| Seurat | [CRAN](https://cran.r-project.org/package=Seurat) | [GitHub](https://github.com/satijalab/seurat) | [Website](https://satijalab.org/seurat) |
| SeuratObject | [CRAN](https://cran.r-project.org/package=SeuratObject) | [GitHub](https://github.com/satijalab/seurat-object) | [Website](https://satijalab.github.io/seurat-object/) |
| stringi | [CRAN](https://cran.r-project.org/package=stringi) | [GitHub](https://github.com/gagolews/stringi) | [Website](https://stringi.gagolewski.com/) |
