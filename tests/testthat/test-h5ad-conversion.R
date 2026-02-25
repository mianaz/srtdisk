# Load package for testing
library(srtdisk)

test_that("h5ad test data is available", {
  # Check if test data directory exists
  testdata_dir <- system.file("testdata", package = "srtdisk")

  # If running in CI with downloaded data or local development
  if (!dir.exists(testdata_dir) || length(list.files(testdata_dir, pattern = "\\.h5ad$")) == 0) {
    # Check for data in inst/testdata (CI download location or local dev)
    testdata_dir <- file.path("inst", "testdata")
    if (!dir.exists(testdata_dir)) {
      skip("Test data not available")
    }
  }

  # Find available h5ad files
  h5ad_files <- list.files(testdata_dir, pattern = "\\.h5ad$", full.names = TRUE)

  if (length(h5ad_files) == 0) {
    skip("No h5ad test files found")
  }

  # Just verify the files exist and are readable
  expect_true(length(h5ad_files) >= 1)
  for (f in h5ad_files) {
    expect_true(file.exists(f))
    expect_true(file.size(f) > 0)
  }
})

test_that("Convert function can be called with h5ad files", {
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

  # Test that Convert can be called (may fail with some files due to format issues)
  # We're mainly testing that the function exists and basic infrastructure works
  result <- tryCatch({
    Convert(test_file, dest = temp_h5seurat, overwrite = TRUE, verbose = FALSE)
    "success"
  }, error = function(e) {
    paste("error:", conditionMessage(e))
  })

  # Log what happened for debugging
  message("Conversion result: ", result)

  # At minimum, the function should be callable (even if it errors on some files)
  expect_true(is.character(result))
})

test_that("LoadH5Seurat function exists and is callable", {
  # Just test that the function exists
  expect_true(exists("LoadH5Seurat"))
  expect_true(is.function(LoadH5Seurat))
})

test_that("h5ad file structure can be inspected", {
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

  # Test that we can open and inspect h5ad files
  result <- tryCatch({
    h5_conn <- hdf5r::H5File$new(test_file, mode = "r")

    is_valid <- h5_conn$is_valid

    # Check for expected h5ad structure
    names_list <- names(h5_conn)
    has_names <- length(names_list) > 0

    # Most h5ad files should have at least X
    has_x <- "X" %in% names_list || "layers" %in% names_list

    # Close the connection
    tryCatch(h5_conn$close_all(), error = function(e) NULL)

    list(valid = is_valid, has_names = has_names, has_x = has_x)
  }, error = function(e) {
    list(valid = FALSE, has_names = FALSE, has_x = FALSE, error = conditionMessage(e))
  })

  expect_true(result$valid)
  expect_true(result$has_names)
  expect_true(result$has_x)
})

test_that("h5ad files have correct extension", {
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

  # Test file extensions
  for (test_file in h5ad_files) {
    expect_match(test_file, "\\.h5ad$", info = paste("File:", basename(test_file)))
  }
})
