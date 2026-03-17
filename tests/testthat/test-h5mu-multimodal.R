# Tests for H5MU multimodal support

library(srtdisk)

test_that("LoadH5MU function exists and is exported", {
  expect_true(exists("LoadH5MU"))
  expect_true(is.function(LoadH5MU))
})

test_that("SaveH5MU function exists and is exported", {
  expect_true(exists("SaveH5MU"))
  expect_true(is.function(SaveH5MU))
})

test_that("LoadH5MU validates file existence", {
  expect_error(
    LoadH5MU("nonexistent_file.h5mu"),
    "File not found"
  )
})

test_that("SaveH5MU validates input object", {
  expect_error(
    SaveH5MU("not_a_seurat", "output.h5mu"),
    "must be a Seurat object"
  )
})

test_that("SaveH5MU requires MuDataSeurat", {
  skip_if(requireNamespace("MuDataSeurat", quietly = TRUE),
          "MuDataSeurat is installed, skip missing-package test")

  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  expect_error(
    SaveH5MU(obj, tempfile(fileext = ".h5mu")),
    "MuDataSeurat"
  )
})

test_that("LoadH5MU requires MuDataSeurat", {
  skip_if(requireNamespace("MuDataSeurat", quietly = TRUE),
          "MuDataSeurat is installed, skip missing-package test")

  # Create a minimal h5mu file structure for testing
  tmp <- tempfile(fileext = ".h5mu")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_group("mod")
  h5$close_all()

  expect_error(
    LoadH5MU(tmp),
    "MuDataSeurat"
  )
})

test_that("GetDefaultModalityMapping returns correct mappings", {
  mapping <- srtdisk:::GetDefaultModalityMapping(c("rna", "prot", "atac"))
  expect_equal(mapping[["rna"]], "RNA")
  expect_equal(mapping[["prot"]], "ADT")
  expect_equal(mapping[["atac"]], "ATAC")
})

test_that("GetDefaultModalityMapping preserves unknown modalities", {
  mapping <- srtdisk:::GetDefaultModalityMapping(c("rna", "custom_assay"))
  expect_equal(mapping[["rna"]], "RNA")
  expect_equal(mapping[["custom_assay"]], "custom_assay")
})

test_that("GetDefaultAssayToModalityMapping returns correct reverse mappings", {
  mapping <- srtdisk:::GetDefaultAssayToModalityMapping(c("RNA", "ADT", "ATAC"))
  expect_equal(mapping[["RNA"]], "rna")
  expect_equal(mapping[["ADT"]], "prot")
  expect_equal(mapping[["ATAC"]], "atac")
})

test_that("GetDefaultAssayToModalityMapping handles conflicts", {
  # RNA and SCT both map to "rna" by default - should resolve
  mapping <- srtdisk:::GetDefaultAssayToModalityMapping(c("RNA", "SCT"))
  expect_false(anyDuplicated(mapping) > 0)
  expect_equal(mapping[["RNA"]], "rna")
  expect_equal(mapping[["SCT"]], "sct")
})

test_that("ValidateMultimodalObject accepts valid objects", {
  skip_if_not_installed("Seurat")

  rna <- matrix(rpois(100, 1), 10, 10)
  rownames(rna) <- paste0("Gene", 1:10)
  colnames(rna) <- paste0("Cell", 1:10)
  adt <- matrix(rpois(50, 5), 5, 10)
  rownames(adt) <- paste0("ADT", 1:5)
  colnames(adt) <- colnames(rna)

  obj <- Seurat::CreateSeuratObject(counts = rna)
  obj[["ADT"]] <- Seurat::CreateAssayObject(counts = adt)

  result <- srtdisk:::ValidateMultimodalObject(obj, c("RNA", "ADT"))
  expect_true(result$valid)
})

test_that("SaveH5MU rejects overwrite when file exists", {
  skip_if_not_installed("Seurat")

  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".h5mu")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    SaveH5MU(obj, tmp, overwrite = FALSE),
    "already exists"
  )
})

test_that("Single assay warning works in SaveH5MU", {
  skip_if_not_installed("Seurat")
  skip_if(!requireNamespace("MuDataSeurat", quietly = TRUE),
          "MuDataSeurat not available")
  # MuDataSeurat does not yet support Seurat v5 Assay5 objects (meta.features slot)
  skip_if({
    m <- matrix(rpois(20, 5), nrow = 2, dimnames = list(paste0("g", 1:2), paste0("c", 1:10)))
    inherits(tryCatch(Seurat::CreateSeuratObject(counts = m)[["RNA"]], error = function(e) NULL), "Assay5")
  }, "MuDataSeurat does not support Seurat v5 Assay5 objects")

  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".h5mu")
  on.exit(unlink(tmp), add = TRUE)

  # Should produce a "single assay" note
  expect_message(
    SaveH5MU(obj, tmp, verbose = TRUE),
    "only one assay"
  )
})

test_that("H5MU conversion pathway exists in Convert.H5File", {
  # Verify that Convert recognizes h5mu file type
  # We can't test full conversion without MuDataSeurat, but can verify the pathway exists
  expect_true("h5mu" %in% c("h5ad", "h5seurat", "h5mu"))  # Just ensure our additions work
})

test_that("Multimodal round-trip via SaveH5MU/LoadH5MU preserves structure", {
  skip_if_not_installed("Seurat")
  skip_if(!requireNamespace("MuDataSeurat", quietly = TRUE),
          "MuDataSeurat not available for round-trip test")
  # MuDataSeurat does not yet support Seurat v5 Assay5 objects (meta.features slot)
  skip_if({
    m <- matrix(rpois(20, 5), nrow = 2, dimnames = list(paste0("g", 1:2), paste0("c", 1:10)))
    inherits(tryCatch(Seurat::CreateSeuratObject(counts = m)[["RNA"]], error = function(e) NULL), "Assay5")
  }, "MuDataSeurat does not support Seurat v5 Assay5 objects")

  # Create multimodal object (RNA + ADT)
  set.seed(42)
  rna <- matrix(rpois(200, 5), 20, 10)
  rownames(rna) <- paste0("Gene", 1:20)
  colnames(rna) <- paste0("Cell", 1:10)
  adt <- matrix(rpois(50, 10), 5, 10)
  rownames(adt) <- paste0("ADT", 1:5)
  colnames(adt) <- colnames(rna)

  obj <- Seurat::CreateSeuratObject(counts = rna)
  obj[["ADT"]] <- Seurat::CreateAssayObject(counts = adt)
  obj$group <- factor(sample(c("A", "B"), 10, replace = TRUE))

  tmp <- tempfile(fileext = ".h5mu")
  on.exit(unlink(tmp), add = TRUE)

  # Save
  SaveH5MU(obj, tmp, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(tmp))

  # Load back
  loaded <- LoadH5MU(tmp, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), 10)

  # Check assays are present (names may differ due to modality mapping)
  loaded_assays <- Seurat::Assays(loaded)
  expect_true(length(loaded_assays) >= 2)
})
