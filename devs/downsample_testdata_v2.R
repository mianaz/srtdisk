#!/usr/bin/env Rscript

# Downsample vis.heart.h5ad test data - simplified version
# Goal: Reduce spots to half and genes to 10k for faster vignette building

library(srtdisk)
library(Seurat)
library(SeuratObject)

# Path to original data
orig_path <- "inst/testdata/vis.heart.h5ad"

# Convert to h5seurat first
message("Converting original h5ad to h5seurat...")
Convert(orig_path, dest = "h5seurat", overwrite = TRUE)

# Load as Seurat object
message("Loading h5seurat file...")
heart_list <- LoadH5Seurat("inst/testdata/vis.heart.h5seurat")

if (!is.list(heart_list)) {
  stop("Expected a list of libraries, but got a single object")
}

message(sprintf("Found %d libraries", length(heart_list)))

# --- Process single library version ---
message("\n=== Creating single library version ===")
heart_single <- heart_list[[1]]
lib_name_single <- names(heart_list)[1]

message(sprintf("Using library: %s", lib_name_single))
message(sprintf("  Original: %d spots, %d genes", ncol(heart_single), nrow(heart_single)))

# Downsample spots (half)
n_spots <- ncol(heart_single)
n_spots_keep <- ceiling(n_spots / 2)
set.seed(42)
spots_keep <- sample(colnames(heart_single), n_spots_keep)

# Downsample genes (top 10k by variance)
n_genes <- nrow(heart_single)
if (n_genes > 10000) {
  message("  Selecting top 10k variable genes...")
  gene_vars <- apply(GetAssayData(heart_single, slot = "data"), 1, var)
  genes_keep <- names(sort(gene_vars, decreasing = TRUE)[1:10000])
} else {
  genes_keep <- rownames(heart_single)
}

# Subset
heart_single_sub <- subset(heart_single,
                           features = genes_keep,
                           cells = spots_keep)

message(sprintf("  Downsampled: %d spots, %d genes",
               ncol(heart_single_sub), nrow(heart_single_sub)))

# Save single library
message("  Saving single library dataset...")
SaveH5Seurat(heart_single_sub,
             filename = "inst/testdata/vis.heart.single.h5seurat",
             overwrite = TRUE)
Convert("inst/testdata/vis.heart.single.h5seurat",
        dest = "h5ad",
        overwrite = TRUE)

message("  Single library saved: vis.heart.single.h5ad")

# --- Process multi-library version ---
message("\n=== Creating multi-library version ===")

# Keep only first 2 libraries to reduce size
n_libs_to_keep <- min(2, length(heart_list))
heart_list_subset <- heart_list[1:n_libs_to_keep]

message(sprintf("Keeping %d libraries", n_libs_to_keep))

# Downsample each library
for (i in seq_along(heart_list_subset)) {
  lib_name <- names(heart_list_subset)[i]
  heart <- heart_list_subset[[i]]

  message(sprintf("\nLibrary %d: %s", i, lib_name))
  message(sprintf("  Original: %d spots, %d genes", ncol(heart), nrow(heart)))

  # Downsample spots (half)
  n_spots <- ncol(heart)
  n_spots_keep <- ceiling(n_spots / 2)
  set.seed(42 + i)  # Different seed for each library
  spots_keep <- sample(colnames(heart), n_spots_keep)

  # Downsample genes (top 10k by variance)
  n_genes <- nrow(heart)
  if (n_genes > 10000) {
    gene_vars <- apply(GetAssayData(heart, slot = "data"), 1, var)
    genes_keep <- names(sort(gene_vars, decreasing = TRUE)[1:10000])
  } else {
    genes_keep <- rownames(heart)
  }

  # Subset
  heart_sub <- subset(heart,
                     features = genes_keep,
                     cells = spots_keep)

  message(sprintf("  Downsampled: %d spots, %d genes",
                 ncol(heart_sub), nrow(heart_sub)))

  # Update in list
  heart_list_subset[[i]] <- heart_sub
}

# Now we need to save this as a multi-library h5ad
# The approach: add library_id metadata and save separately, then manually combine

message("\n=== Saving multi-library version ===")

# For multi-library h5ad, we need to ensure proper structure
# Let's use a different approach: save them individually with library_id,
# then use Python to combine them properly

# Actually, let's just join them in a simple way
# Add library_id to metadata first
for (i in seq_along(heart_list_subset)) {
  lib_name <- names(heart_list_subset)[i]
  heart_list_subset[[i]]$library_id <- lib_name

  # Make cell names unique across libraries
  new_cell_names <- paste0(lib_name, "_", colnames(heart_list_subset[[i]]))
  colnames(heart_list_subset[[i]]) <- new_cell_names
}

# Now join the layers properly using JoinLayers
message("Merging libraries using JoinLayers...")

# First merge the Seurat objects
heart_multi <- heart_list_subset[[1]]
if (length(heart_list_subset) > 1) {
  for (i in 2:length(heart_list_subset)) {
    heart_multi <- merge(heart_multi, heart_list_subset[[i]])
  }
}

message(sprintf("Merged: %d spots, %d genes", ncol(heart_multi), nrow(heart_multi)))

# Now join the layers to create a single counts/data matrix
message("Joining layers...")
heart_multi <- JoinLayers(heart_multi)

message("Layers joined successfully")
message(sprintf("Final: %d spots, %d genes", ncol(heart_multi), nrow(heart_multi)))

# Save multi-library version
message("Saving multi-library dataset...")
SaveH5Seurat(heart_multi,
             filename = "inst/testdata/vis.heart.h5seurat",
             overwrite = TRUE)
Convert("inst/testdata/vis.heart.h5seurat",
        dest = "h5ad",
        overwrite = TRUE)

message("  Multi-library saved: vis.heart.h5ad")

# Clean up temporary h5seurat files
message("\n=== Cleaning up ===")
temp_files <- c(
  "inst/testdata/vis.heart.h5seurat",
  "inst/testdata/vis.heart.single.h5seurat"
)

for (f in temp_files) {
  if (file.exists(f)) {
    file.remove(f)
    message(sprintf("  Removed: %s", f))
  }
}

# Check final file sizes
message("\n=== Final file sizes ===")
files_to_check <- c(
  "inst/testdata/vis.heart.h5ad",
  "inst/testdata/vis.heart.single.h5ad"
)

for (f in files_to_check) {
  if (file.exists(f)) {
    size_mb <- file.size(f) / 1024^2
    message(sprintf("  %s: %.1f MB", basename(f), size_mb))
  }
}

message("\nAll done!")
