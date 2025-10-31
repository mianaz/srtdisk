library(srtdisk)

test_that("pbmc_small roundtrip through h5Seurat succeeds", {
  skip_if_not_installed("Seurat")

  data("pbmc_small", package = "Seurat")
  expect_true(exists("pbmc_small"), info = "pbmc_small dataset failed to load")

  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(temp_h5seurat), add = TRUE)

  expect_no_error(
    SaveH5Seurat(pbmc_small, filename = temp_h5seurat, overwrite = TRUE, verbose = FALSE)
  )

  expect_true(file.exists(temp_h5seurat))
  expect_gt(file.size(temp_h5seurat), 0)

  reloaded <- LoadH5Seurat(temp_h5seurat, verbose = FALSE)

  expect_s4_class(reloaded, "Seurat")
  expect_equal(ncol(reloaded), ncol(pbmc_small))
  expect_equal(nrow(reloaded), nrow(pbmc_small))
  expect_identical(colnames(reloaded), colnames(pbmc_small))
  expect_identical(rownames(reloaded), rownames(pbmc_small))
})
