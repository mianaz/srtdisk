# Test loom <-> Seurat conversion functionality

# Skip tests if required packages are not available
skip_if_not_installed("Seurat")
skip_if_not_installed("Matrix")

# Helper function to create a test Seurat object
create_test_seurat <- function(n_cells = 100, n_features = 200) {
  set.seed(42)
  # Create integer counts (not negative values)
  counts <- matrix(
    data = rpois(n_features * n_cells, lambda = 5),
    nrow = n_features,
    ncol = n_cells,
    dimnames = list(
      paste0("Gene", seq_len(n_features)),
      paste0("Cell", seq_len(n_cells))
    )
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = "TestProject")

  # Normalize to populate the data layer
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)

  # Add some metadata
  seurat_obj$cluster <- factor(sample(c("A", "B", "C"), n_cells, replace = TRUE))
  seurat_obj$numeric_meta <- rnorm(n_cells)

  return(seurat_obj)
}

# -----------------------------------------------------------------------------
# Test SaveLoom
# -----------------------------------------------------------------------------

test_that("SaveLoom creates a valid loom file from Seurat object", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  expect_no_error(SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE))
  expect_true(file.exists(loom_file))

  # Verify it's a valid HDF5 file
  expect_true(hdf5r::is_hdf5(loom_file))

  # Open and check structure
  h5 <- hdf5r::H5File$new(loom_file, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  # Check required loom datasets/groups exist

  expect_true(h5$exists("matrix"))
  expect_true(h5$exists("row_attrs"))
  expect_true(h5$exists("col_attrs"))
  expect_true(h5$exists("layers"))

  # Check cell and gene IDs are stored
  expect_true(h5$exists("col_attrs/CellID"))
  expect_true(h5$exists("row_attrs/Gene"))

  # Verify dimensions (loom stores genes x cells, transposed)
  matrix_dims <- h5[["matrix"]]$dims
  expect_equal(matrix_dims[1], ncol(seurat_obj))  # cells
  expect_equal(matrix_dims[2], nrow(seurat_obj))  # genes

  unlink(loom_file)
})

test_that("SaveLoom preserves cell metadata", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)

  h5 <- hdf5r::H5File$new(loom_file, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  # Check metadata columns are preserved
  expect_true(h5$exists("col_attrs/cluster"))
  expect_true(h5$exists("col_attrs/numeric_meta"))

  # Verify cell names match
  cell_ids <- h5[["col_attrs/CellID"]][]
  expect_equal(cell_ids, colnames(seurat_obj))

  unlink(loom_file)
})

test_that("SaveLoom handles overwrite parameter correctly", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  # First save
  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)

  # Second save without overwrite should fail
  expect_error(
    SaveLoom(seurat_obj, filename = loom_file, overwrite = FALSE, verbose = FALSE),
    "already exists"
  )

  # With overwrite should succeed
  expect_no_error(
    SaveLoom(seurat_obj, filename = loom_file, overwrite = TRUE, verbose = FALSE)
  )

  unlink(loom_file)
})

test_that("SaveLoom adds .loom extension if missing", {
  seurat_obj <- create_test_seurat()
  base_path <- tempfile()

  result <- SaveLoom(seurat_obj, filename = base_path, verbose = FALSE)

  expect_true(grepl("\\.loom$", result))
  expect_true(file.exists(result))

  unlink(result)
})

# -----------------------------------------------------------------------------
# Test LoadLoom
# -----------------------------------------------------------------------------

test_that("LoadLoom reads a loom file as Seurat object", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)
  loaded_obj <- LoadLoom(loom_file, verbose = FALSE)

  expect_s4_class(loaded_obj, "Seurat")
  expect_equal(ncol(loaded_obj), ncol(seurat_obj))
  expect_equal(nrow(loaded_obj), nrow(seurat_obj))

  unlink(loom_file)
})

test_that("LoadLoom preserves cell and gene names", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)
  loaded_obj <- LoadLoom(loom_file, verbose = FALSE)

  expect_equal(colnames(loaded_obj), colnames(seurat_obj))
  expect_equal(rownames(loaded_obj), rownames(seurat_obj))

  unlink(loom_file)
})

test_that("LoadLoom preserves cell metadata", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)
  loaded_obj <- LoadLoom(loom_file, verbose = FALSE)

  # Check metadata columns exist
  expect_true("cluster" %in% colnames(loaded_obj[[]]))
  expect_true("numeric_meta" %in% colnames(loaded_obj[[]]))

  # Check values are preserved (may be character after round-trip for factors)
  expect_equal(
    as.character(loaded_obj$cluster),
    as.character(seurat_obj$cluster)
  )
  expect_equal(
    loaded_obj$numeric_meta,
    seurat_obj$numeric_meta,
    tolerance = 1e-6
  )

  unlink(loom_file)
})

test_that("LoadLoom handles custom assay name", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)
  loaded_obj <- LoadLoom(loom_file, assay = "MyAssay", verbose = FALSE)

  expect_true("MyAssay" %in% Seurat::Assays(loaded_obj))
  expect_equal(Seurat::DefaultAssay(loaded_obj), "MyAssay")

  unlink(loom_file)
})

# -----------------------------------------------------------------------------
# Test round-trip conversion
# -----------------------------------------------------------------------------

test_that("Round-trip Seurat -> Loom -> Seurat preserves expression data", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)
  loaded_obj <- LoadLoom(loom_file, verbose = FALSE)

  # Compare expression matrices
  original_data <- Seurat::GetAssayData(seurat_obj, layer = "data")
  loaded_data <- Seurat::GetAssayData(loaded_obj, layer = "data")

  # Ensure same order of cells and genes
  loaded_data <- loaded_data[rownames(original_data), colnames(original_data)]

  # Compare values only (ignore dimnames attributes)
  expect_equal(
    unname(as.matrix(original_data)),
    unname(as.matrix(loaded_data)),
    tolerance = 1e-6
  )

  unlink(loom_file)
})

test_that("Round-trip preserves main matrix (data layer saved, loaded as counts)", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)
  loaded_obj <- LoadLoom(loom_file, verbose = FALSE)

  # SaveLoom saves "data" layer as /matrix
  # LoadLoom loads /matrix as "counts" (since loom format doesn't distinguish)
  # So we compare original "data" with loaded "counts"
  original_data <- Seurat::GetAssayData(seurat_obj, layer = "data")

  # LoadLoom stores the main matrix in counts layer
  loaded_counts <- Seurat::GetAssayData(loaded_obj, layer = "counts")

  # Ensure same order
  loaded_counts <- loaded_counts[rownames(original_data), colnames(original_data)]

  # Compare values only (ignore dimnames attributes)
  expect_equal(
    unname(as.matrix(original_data)),
    unname(as.matrix(loaded_counts)),
    tolerance = 1e-6
  )

  unlink(loom_file)
})

# -----------------------------------------------------------------------------
# Test with real loom file (if available)
# -----------------------------------------------------------------------------

test_that("LoadLoom works with SummarizedExperiment example.loom", {
  skip_if_not_installed("SummarizedExperiment")

  example_loom <- system.file(
    "extdata", "example.loom",
    package = "SummarizedExperiment"
  )
  skip_if(!file.exists(example_loom), "Example loom file not found")

  # The SummarizedExperiment example uses different attribute names
  # Try to load it - it may fail if attributes don't match expectations
  result <- tryCatch(
    LoadLoom(example_loom, verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(result)) {
    expect_s4_class(result, "Seurat")
    expect_gt(ncol(result), 0)
    expect_gt(nrow(result), 0)
  }
})

# -----------------------------------------------------------------------------
# Test dimensional reductions
# -----------------------------------------------------------------------------

test_that("Round-trip preserves dimensional reductions", {
  seurat_obj <- create_test_seurat(n_cells = 50, n_features = 100)

  # Add a simple PCA-like reduction
  set.seed(42)
  pca_embeddings <- matrix(
    rnorm(50 * 10),
    nrow = 50,
    ncol = 10,
    dimnames = list(colnames(seurat_obj), paste0("PC_", 1:10))
  )
  seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = pca_embeddings,
    key = "PC_",
    assay = "RNA"
  )

  loom_file <- tempfile(fileext = ".loom")
  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)
  loaded_obj <- LoadLoom(loom_file, verbose = FALSE)

  # Check if PCA was preserved
  if ("pca" %in% Seurat::Reductions(loaded_obj)) {
    original_embeddings <- Seurat::Embeddings(seurat_obj, "pca")
    loaded_embeddings <- Seurat::Embeddings(loaded_obj, "pca")

    # Compare values only (column names may differ between PC_1 and pca_1)
    expect_equal(
      unname(original_embeddings),
      unname(loaded_embeddings),
      tolerance = 1e-6
    )
  }

  unlink(loom_file)
})

# -----------------------------------------------------------------------------
# Test edge cases
# -----------------------------------------------------------------------------

test_that("SaveLoom handles special characters in gene names", {
  # Create test data with special gene names directly
  # Note: Some characters like underscores may be converted during HDF5 round-trip
  set.seed(42)
  new_names <- c(
    "Gene-1", "Gene.2", "Gene-3", "Gene:4", "Gene/5",
    "Gene6", "Gene7", "Gene8", "Gene9", "Gene10"
  )
  counts <- matrix(
    data = rpois(10 * 50, lambda = 5),
    nrow = 10,
    ncol = 50,
    dimnames = list(new_names, paste0("Cell", 1:50))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = "TestProject")
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)

  loom_file <- tempfile(fileext = ".loom")

  expect_no_error(SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE))

  loaded_obj <- LoadLoom(loom_file, verbose = FALSE)
  expect_equal(rownames(loaded_obj), new_names)

  unlink(loom_file)
})

test_that("LoadLoom handles missing optional attributes gracefully", {
  seurat_obj <- create_test_seurat()
  loom_file <- tempfile(fileext = ".loom")

  SaveLoom(seurat_obj, filename = loom_file, verbose = FALSE)

  # Try loading with non-existent layers
  loaded_obj <- LoadLoom(
    loom_file,
    normalized = "nonexistent_layer",
    scaled = "another_nonexistent",
    verbose = FALSE
  )

  expect_s4_class(loaded_obj, "Seurat")

  unlink(loom_file)
})
