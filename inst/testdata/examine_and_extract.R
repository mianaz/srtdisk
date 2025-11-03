#!/usr/bin/env Rscript

# Script to examine h5ad file structure and extract a single library

library(hdf5r)
library(anndata)

# Open the h5ad file
h5_file <- "inst/testdata/vis.heart.h5ad"
h5 <- H5File$new(h5_file, mode = "r")

cat("=== H5AD File Structure ===\n")
cat("Top-level groups:\n")
print(names(h5))

# Check for layers (different modalities/libraries)
if ("layers" %in% names(h5)) {
  cat("\n=== Layers (modalities) ===\n")
  layers <- h5[["layers"]]
  print(names(layers))
}

# Check obsm (spatial coordinates, embeddings, etc.)
if ("obsm" %in% names(h5)) {
  cat("\n=== Observation matrices (obsm) ===\n")
  obsm <- h5[["obsm"]]
  print(names(obsm))
}

# Check obs (cell metadata)
if ("obs" %in% names(h5)) {
  cat("\n=== Observation metadata (obs) ===\n")
  obs <- h5[["obs"]]
  if (inherits(obs, "H5Group")) {
    print(names(obs))
  }
}

# Check var (gene metadata)
if ("var" %in% names(h5)) {
  cat("\n=== Variable metadata (var) ===\n")
  var <- h5[["var"]]
  if (inherits(var, "H5Group")) {
    print(names(var))
  }
}

# Check main data matrix X
if ("X" %in% names(h5)) {
  cat("\n=== Main data matrix (X) ===\n")
  x_data <- h5[["X"]]
  cat("Type:", class(x_data), "\n")
  if (inherits(x_data, "H5Group")) {
    cat("X is a group, contents:\n")
    print(names(x_data))
  } else {
    cat("Shape:", dim(x_data), "\n")
  }
}

h5$close_all()

cat("\n=== Attempting to read with anndata package ===\n")
tryCatch({
  adata <- read_h5ad(h5_file)
  cat("Successfully read file with anndata\n")
  cat("Shape: ", dim(adata), "\n")
  cat("Available layers: ", names(adata$layers), "\n")

  # Extract first library/modality and save
  cat("\n=== Extracting single library ===\n")

  # If there are multiple layers, take the first one
  if (length(names(adata$layers)) > 0) {
    cat("Found", length(names(adata$layers)), "layers\n")
    cat("Extracting first layer:", names(adata$layers)[1], "\n")

    # Create new anndata with just one layer
    # We'll keep the same obs and var but use only one data matrix
    layer_name <- names(adata$layers)[1]

    # For now, let's create a subset by selecting the main X matrix
    # This is a simplified extraction
    cat("Creating single-library h5ad file...\n")

  } else {
    cat("No separate layers found. File may already be single-modality.\n")
  }

}, error = function(e) {
  cat("Error reading with anndata:", conditionMessage(e), "\n")
})

cat("\nExamination complete.\n")
