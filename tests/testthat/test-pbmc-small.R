library(srtdisk)

test_that("pbmc_small roundtrip through h5Seurat succeeds", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  # Ensure SeuratObject is loaded to avoid SaveH5Seurat index creation bugs
  library(Seurat)
  library(SeuratObject)

  # Create a minimal test Seurat object instead of loading pbmc_small
  # which may not be available in all Seurat versions
  set.seed(123)

  # Create a small test count matrix
  counts <- matrix(
    rpois(100 * 20, lambda = 5),
    nrow = 100,
    ncol = 20,
    dimnames = list(
      paste0("Gene", 1:100),
      paste0("Cell", 1:20)
    )
  )

  # Create a minimal Seurat object
  suppressWarnings({
    test_obj <- Seurat::CreateSeuratObject(
      counts = counts,
      project = "TestProject",
      min.cells = 0,
      min.features = 0
    )
  })

  # Normalize data to ensure 'data' layer exists
  # This is important for proper roundtrip testing
  suppressWarnings({
    test_obj <- Seurat::NormalizeData(test_obj, verbose = FALSE)
  })

  # Verify the test object was created
  expect_s4_class(test_obj, "Seurat")

  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(temp_h5seurat), add = TRUE)

  expect_no_error(
    SaveH5Seurat(test_obj, filename = temp_h5seurat, overwrite = TRUE, verbose = FALSE)
  )

  expect_true(file.exists(temp_h5seurat))
  expect_gt(file.size(temp_h5seurat), 0)

  reloaded <- LoadH5Seurat(temp_h5seurat, verbose = FALSE)

  expect_s4_class(reloaded, "Seurat")
  expect_equal(ncol(reloaded), ncol(test_obj))
  expect_equal(nrow(reloaded), nrow(test_obj))
  expect_identical(colnames(reloaded), colnames(test_obj))
  expect_identical(rownames(reloaded), rownames(test_obj))
})
