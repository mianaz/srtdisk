#!/usr/bin/env Rscript

# Downsample vis.heart.h5ad test data
# Goal: Reduce spots to half and genes to 10k for faster vignette building

library(srtdisk)
library(Seurat)

# Path to original data
orig_path <- "inst/testdata/vis.heart.h5ad"

# Convert to h5seurat first
message("Converting original h5ad to h5seurat...")
Convert(orig_path, dest = "h5seurat", overwrite = TRUE)

# Load as Seurat object
message("Loading h5seurat file...")
heart_list <- LoadH5Seurat("inst/testdata/vis.heart.h5seurat")

# Check if it's a list (multiple libraries) or single object
if (is.list(heart_list)) {
  message(sprintf("Found %d libraries", length(heart_list)))

  # Process each library
  for (i in seq_along(heart_list)) {
    lib_name <- names(heart_list)[i]
    heart <- heart_list[[i]]

    message(sprintf("\nProcessing library %d: %s", i, lib_name))
    message(sprintf("  Original: %d spots, %d genes", ncol(heart), nrow(heart)))

    # Downsample spots (half)
    n_spots <- ncol(heart)
    n_spots_keep <- ceiling(n_spots / 2)
    set.seed(42)
    spots_keep <- sample(colnames(heart), n_spots_keep)

    # Downsample genes (top 10k by variance)
    n_genes <- nrow(heart)
    if (n_genes > 10000) {
      gene_vars <- apply(GetAssayData(heart, slot = "data"), 1, var)
      genes_keep <- names(sort(gene_vars, decreasing = TRUE)[1:10000])
    } else {
      genes_keep <- rownames(heart)
    }

    # Subset the object
    heart_sub <- subset(heart,
                       features = genes_keep,
                       cells = spots_keep)

    message(sprintf("  Downsampled: %d spots, %d genes",
                   ncol(heart_sub), nrow(heart_sub)))

    # Update the list
    heart_list[[i]] <- heart_sub
  }

  # Save downsampled multi-library h5ad
  message("\nSaving downsampled multi-library dataset...")

  # For multi-library, we need to merge them back or save separately
  # Let's save the first library as single, and create a merged version

  # Save first library as single
  message("Saving single library (first library)...")
  SaveH5Seurat(heart_list[[1]],
               filename = "inst/testdata/vis.heart.single.h5seurat",
               overwrite = TRUE)
  Convert("inst/testdata/vis.heart.single.h5seurat",
          dest = "h5ad",
          overwrite = TRUE)

  # For multi-library, we'll merge them
  message("Merging libraries for multi-library dataset...")
  # Add library_id to each object's metadata
  for (i in seq_along(heart_list)) {
    lib_name <- names(heart_list)[i]
    heart_list[[i]]$library_id <- lib_name
  }

  # Merge all libraries
  heart_merged <- heart_list[[1]]
  if (length(heart_list) > 1) {
    for (i in 2:length(heart_list)) {
      heart_merged <- merge(heart_merged, heart_list[[i]])
    }
  }

  message(sprintf("Merged object: %d spots, %d genes",
                 ncol(heart_merged), nrow(heart_merged)))

  # Save merged multi-library dataset
  message("Saving downsampled multi-library dataset...")
  SaveH5Seurat(heart_merged,
               filename = "inst/testdata/vis.heart.h5seurat",
               overwrite = TRUE)
  Convert("inst/testdata/vis.heart.h5seurat",
          dest = "h5ad",
          overwrite = TRUE)

} else {
  # Single library case
  heart <- heart_list
  message(sprintf("Single library: %d spots, %d genes", ncol(heart), nrow(heart)))

  # Downsample spots (half)
  n_spots <- ncol(heart)
  n_spots_keep <- ceiling(n_spots / 2)
  set.seed(42)
  spots_keep <- sample(colnames(heart), n_spots_keep)

  # Downsample genes (top 10k by variance)
  n_genes <- nrow(heart)
  if (n_genes > 10000) {
    gene_vars <- apply(GetAssayData(heart, slot = "data"), 1, var)
    genes_keep <- names(sort(gene_vars, decreasing = TRUE)[1:10000])
  } else {
    genes_keep <- rownames(heart)
  }

  # Subset the object
  heart_sub <- subset(heart,
                     features = genes_keep,
                     cells = spots_keep)

  message(sprintf("Downsampled: %d spots, %d genes",
                 ncol(heart_sub), nrow(heart_sub)))

  # Save both versions
  message("Saving downsampled datasets...")

  # Save as both single and multi (same in this case)
  SaveH5Seurat(heart_sub,
               filename = "inst/testdata/vis.heart.single.h5seurat",
               overwrite = TRUE)
  Convert("inst/testdata/vis.heart.single.h5seurat",
          dest = "h5ad",
          overwrite = TRUE)

  SaveH5Seurat(heart_sub,
               filename = "inst/testdata/vis.heart.h5seurat",
               overwrite = TRUE)
  Convert("inst/testdata/vis.heart.h5seurat",
          dest = "h5ad",
          overwrite = TRUE)
}

message("\nDone! Test data downsampled successfully.")
message("Files created:")
message("  - inst/testdata/vis.heart.h5ad (downsampled multi-library)")
message("  - inst/testdata/vis.heart.single.h5ad (single library)")

# Clean up temporary files
message("\nCleaning up temporary files...")
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

message("\nAll done!")
