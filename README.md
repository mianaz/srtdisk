
<!-- README.md is generated from README.Rmd. Please edit that file -->

# srtdisk v0.1.0

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
| withr | [CRAN](https://cran.r-project.org/package=withr) | [GitHub](https://github.com/r-lib/withr#readme) | [Website](https://withr.r-lib.org) |
