#!/usr/bin/env Rscript

# Script to save SeuratData datasets as RDS files for testing without dependencies

library(Seurat)

# Save stxBrain from SeuratData
cat("=== Saving stxBrain dataset ===\n")
tryCatch({
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    cat("Installing SeuratData...\n")
    remotes::install_github("satijalab/seurat-data", upgrade = "never")
  }

  library(SeuratData)

  # Install and load stxBrain
  if (!"stxBrain" %in% InstalledData()$Dataset) {
    cat("Installing stxBrain dataset...\n")
    InstallData("stxBrain")
  }

  # Load using specific dataset name (posterior1 section)
  stxBrain <- LoadData("stxBrain", type = "posterior1")
  cat("Loaded stxBrain dataset (posterior1)\n")
  cat("  Cells:", ncol(stxBrain), "\n")
  cat("  Features:", nrow(stxBrain), "\n")

  # Save as RDS
  saveRDS(stxBrain, "inst/testdata/stxBrain.rds")
  cat("Saved to inst/testdata/stxBrain.rds\n")
  cat("  File size:", file.size("inst/testdata/stxBrain.rds") / 1024^2, "MB\n")

}, error = function(e) {
  cat("Error saving stxBrain:", conditionMessage(e), "\n")
})

# Save pbmc3k from SeuratData
cat("\n=== Saving pbmc3k dataset ===\n")
tryCatch({
  library(SeuratData)

  # Install and load pbmc3k
  if (!"pbmc3k" %in% InstalledData()$Dataset) {
    cat("Installing pbmc3k dataset...\n")
    InstallData("pbmc3k")
  }

  # Load using LoadData instead of data()
  pbmc3k <- LoadData("pbmc3k")
  cat("Loaded pbmc3k dataset\n")
  cat("  Cells:", ncol(pbmc3k), "\n")
  cat("  Features:", nrow(pbmc3k), "\n")

  # Save as RDS
  saveRDS(pbmc3k, "inst/testdata/pbmc3k.rds")
  cat("Saved to inst/testdata/pbmc3k.rds\n")
  cat("  File size:", file.size("inst/testdata/pbmc3k.rds") / 1024^2, "MB\n")

}, error = function(e) {
  cat("Error saving pbmc3k:", conditionMessage(e), "\n")
})

# Save pbmcMultiome from SeuratData (RNA + ATAC multimodal data)
cat("\n=== Saving pbmcMultiome dataset ===\n")
tryCatch({
  library(SeuratData)

  # Install and load pbmcMultiome
  if (!"pbmcMultiome" %in% InstalledData()$Dataset) {
    cat("Installing pbmcMultiome dataset...\n")
    InstallData("pbmcMultiome")
  }

  # Load the multiome data (has RNA and ATAC assays)
  pbmcMultiome <- LoadData("pbmcMultiome")
  cat("Loaded pbmcMultiome dataset\n")
  cat("  Cells:", ncol(pbmcMultiome), "\n")
  cat("  Assays:", paste(names(pbmcMultiome@assays), collapse = ", "), "\n")

  # Save as RDS
  saveRDS(pbmcMultiome, "inst/testdata/pbmcMultiome.rds")
  cat("Saved to inst/testdata/pbmcMultiome.rds\n")
  cat("  File size:", file.size("inst/testdata/pbmcMultiome.rds") / 1024^2, "MB\n")

}, error = function(e) {
  cat("Error saving pbmcMultiome:", conditionMessage(e), "\n")
})

cat("\nDataset saving complete.\n")
