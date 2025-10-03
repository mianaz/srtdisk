# Test save/load functionality

test_that("Basic save and load works", {
  skip_if_not_installed("Seurat")

  # Create simple test object
  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- CreateSeuratObject(counts = counts)

  # Test save and load
  tmp <- tempfile(fileext = ".h5seurat")
  SaveH5Seurat(obj, tmp, overwrite = TRUE, verbose = FALSE)
  loaded <- LoadH5Seurat(tmp, verbose = FALSE)

  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), 4)
  expect_equal(nrow(loaded), 3)

  unlink(tmp)
})

test_that("Multi-assay objects work correctly", {
  skip_if_not_installed("Seurat")

  # Create multi-assay object
  rna <- matrix(rpois(100, 1), 10, 10)
  adt <- matrix(rpois(50, 5), 5, 10)

  obj <- CreateSeuratObject(counts = rna)
  obj[["ADT"]] <- CreateAssayObject(counts = adt)

  # Test save and load
  tmp <- tempfile(fileext = ".h5seurat")
  SaveH5Seurat(obj, tmp, overwrite = TRUE, verbose = FALSE)
  loaded <- LoadH5Seurat(tmp, verbose = FALSE)

  expect_equal(Assays(loaded), c("RNA", "ADT"))
  expect_equal(nrow(loaded[["RNA"]]), 10)
  expect_equal(nrow(loaded[["ADT"]]), 5)

  unlink(tmp)
})

test_that("V5 objects are handled correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not(packageVersion("Seurat") >= "5.0.0")

  # Create V5 object
  counts <- matrix(rpois(100, 1), 10, 10)
  obj <- CreateSeuratObject(counts = counts)

  # Check if it's V5
  is_v5 <- inherits(obj[["RNA"]], "Assay5")

  if (is_v5) {
    # Test V5 specific features
    tmp <- tempfile(fileext = ".h5seurat")
    SaveH5Seurat(obj, tmp, overwrite = TRUE, verbose = FALSE)
    loaded <- LoadH5Seurat(tmp, verbose = FALSE)

    expect_s4_class(loaded[["RNA"]], "Assay5")
    unlink(tmp)
  }
})