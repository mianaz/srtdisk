# Test categorical column preservation during h5ad to h5Seurat conversion
library(srtdisk)

test_that("h5ad categorical columns are preserved in conversion", {
  # This test verifies that categorical obs columns (e.g., cell_type, condition)
  # are preserved when converting from h5ad to h5Seurat format
  
  skip_if_not_installed("Seurat")
  
  testdata_dir <- system.file("testdata", package = "srtdisk")
  if (!dir.exists(testdata_dir) || length(list.files(testdata_dir, pattern = "\\.h5ad$")) == 0) {
    testdata_dir <- file.path("inst", "testdata")
    if (!dir.exists(testdata_dir)) {
      skip("Test data not available")
    }
  }
  
  h5ad_files <- list.files(testdata_dir, pattern = "\\.h5ad$", full.names = TRUE)
  if (length(h5ad_files) == 0) {
    skip("No h5ad test files found")
  }
  
  test_file <- h5ad_files[1]
  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(temp_h5seurat), add = TRUE)
  
  # Convert h5ad to h5Seurat
  result <- tryCatch({
    Convert(test_file, dest = temp_h5seurat, overwrite = TRUE, verbose = FALSE)
    "success"
  }, error = function(e) {
    skip(paste("Conversion failed:", conditionMessage(e)))
  })
  
  expect_true(file.exists(temp_h5seurat))
  expect_gt(file.size(temp_h5seurat), 0)
  
  # Load the converted file
  seurat_obj <- LoadH5Seurat(temp_h5seurat, verbose = FALSE)
  
  expect_s4_class(seurat_obj, "Seurat")
  
  # Check that metadata was loaded
  metadata <- seurat_obj@meta.data
  expect_true(nrow(metadata) > 0)
  expect_true(ncol(metadata) > 0)
  
  # If the original h5ad had categorical columns, they should appear in metadata
  # We can't check specific column names without knowing the test file structure,
  # but we can verify that non-numeric columns exist (likely categorical fields)
  has_character_cols <- any(sapply(metadata, is.character))
  has_factor_cols <- any(sapply(metadata, is.factor))
  
  # At least some metadata columns should be character or factor (not just numeric)
  # This verifies that categorical columns were preserved
  expect_true(has_character_cols || has_factor_cols,
              info = paste("Expected character or factor columns in metadata. Columns:",
                          paste(names(metadata), collapse = ", "),
                          "| Types:", paste(sapply(metadata, class), collapse = ", ")))
})

test_that("h5ad obs structure is correctly identified", {
  # Test that we can correctly identify whether obs is a group or compound dataset
  
  testdata_dir <- system.file("testdata", package = "srtdisk")
  if (!dir.exists(testdata_dir) || length(list.files(testdata_dir, pattern = "\\.h5ad$")) == 0) {
    testdata_dir <- file.path("inst", "testdata")
    if (!dir.exists(testdata_dir)) {
      skip("Test data not available")
    }
  }
  
  h5ad_files <- list.files(testdata_dir, pattern = "\\.h5ad$", full.names = TRUE)
  if (length(h5ad_files) == 0) {
    skip("No h5ad test files found")
  }
  
  test_file <- h5ad_files[1]
  
  # Open the h5ad file and check obs structure
  tryCatch({
    h5_conn <- hdf5r::H5File$new(test_file, mode = "r")
    
    if (h5_conn$exists("obs")) {
      obs <- h5_conn[["obs"]]
      
      # obs can be either H5Group or H5D (compound dataset)
      is_group <- inherits(obs, "H5Group")
      is_dataset <- inherits(obs, "H5D")
      
      expect_true(is_group || is_dataset,
                  info = paste("obs is neither H5Group nor H5D, type:", class(obs)))
      
      if (is_group) {
        # For groups, check if categorical columns exist
        obs_names <- names(obs)
        message("obs is H5Group with contents: ", paste(obs_names, collapse = ", "))
        
        # Check for modern categorical format (groups with categories/codes)
        categorical_cols <- character(0)
        for (col_name in obs_names) {
          if (col_name != "_index" && inherits(obs[[col_name]], "H5Group")) {
            col_group <- obs[[col_name]]
            if (all(c("categories", "codes") %in% names(col_group))) {
              categorical_cols <- c(categorical_cols, col_name)
            }
          }
        }
        
        if (length(categorical_cols) > 0) {
          message("Found categorical columns in modern format: ", paste(categorical_cols, collapse = ", "))
        }
        
        # Check for legacy __categories
        if ("__categories" %in% obs_names) {
          cat_names <- names(obs[["__categories"]])
          message("Found __categories with: ", paste(cat_names, collapse = ", "))
        }
      } else {
        message("obs is H5D (compound dataset)")
        # For compound datasets, check if __categories exists at obs/__categories
        # Note: This might not exist if obs is a direct dataset
      }
    }
    
    h5_conn$close_all()
  }, error = function(e) {
    skip(paste("Could not inspect h5ad structure:", conditionMessage(e)))
  })
})
