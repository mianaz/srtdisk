# Test native AnnData preprocessing features
# These tests verify column name sanitization, list/dict conversion, and nullable dtype handling

# Skip tests if required packages are not available
skip_if_not_installed("Seurat")
skip_if_not_installed("Matrix")
skip_if_not_installed("hdf5r")

# -----------------------------------------------------------------------------
# Test SanitizeColumnName helper function
# -----------------------------------------------------------------------------

test_that("SanitizeColumnName replaces problematic characters", {
  # Access the internal function from the package namespace
  SanitizeColumnName <- srtdisk:::SanitizeColumnName

  # Test forward slash replacement

  expect_equal(SanitizeColumnName("cell/type"), "cell_type")

  # Test space replacement
  expect_equal(SanitizeColumnName("cell type"), "cell_type")

  # Test comma replacement
  expect_equal(SanitizeColumnName("cell,type"), "cell_type")

  # Test semicolon replacement
  expect_equal(SanitizeColumnName("cell;type"), "cell_type")

  # Test colon replacement
  expect_equal(SanitizeColumnName("cell:type"), "cell_type")

  # Test backslash replacement
  expect_equal(SanitizeColumnName("cell\\type"), "cell_type")

  # Test multiple consecutive problematic characters become single underscore
  expect_equal(SanitizeColumnName("cell//type"), "cell_type")
  expect_equal(SanitizeColumnName("cell  type"), "cell_type")

  # Test leading/trailing underscores are removed
  expect_equal(SanitizeColumnName("/celltype"), "celltype")
  expect_equal(SanitizeColumnName("celltype/"), "celltype")
  expect_equal(SanitizeColumnName("/celltype/"), "celltype")

  # Test names starting with numbers get X prefix
  expect_equal(SanitizeColumnName("1_celltype"), "X1_celltype")
  expect_equal(SanitizeColumnName("123abc"), "X123abc")

  # Test normal names are unchanged
  expect_equal(SanitizeColumnName("cell_type"), "cell_type")
  expect_equal(SanitizeColumnName("CellType"), "CellType")
})

test_that("SanitizeColumnNames handles duplicate names", {
  SanitizeColumnNames <- srtdisk:::SanitizeColumnNames

  # Test that duplicates after sanitization get __dupN suffix
  result <- SanitizeColumnNames(c("cell/type", "cell type", "cell:type"))
  expect_equal(unname(result), c("cell_type", "cell_type__dup1", "cell_type__dup2"))

  # Test that original names are preserved in the names attribute
  expect_equal(names(result), c("cell/type", "cell type", "cell:type"))

  # Test that non-duplicates remain unchanged
  result2 <- SanitizeColumnNames(c("col_a", "col_b", "col_c"))
  expect_equal(unname(result2), c("col_a", "col_b", "col_c"))
})

# -----------------------------------------------------------------------------
# Test FlattenNullable helper function
# -----------------------------------------------------------------------------

test_that("FlattenNullable handles mask+values structures", {
  FlattenNullable <- srtdisk:::FlattenNullable

  # Create a temporary h5 file with mask+values structure
  temp_h5 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5), add = TRUE)

  h5 <- hdf5r::H5File$new(temp_h5, mode = "w")
  on.exit(h5$close_all(), add = TRUE, after = FALSE)

  # Create a nullable column group
  h5$create_group("test_nullable")
  values <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  mask <- c(FALSE, TRUE, FALSE, FALSE, TRUE)  # TRUE = missing in h5ad
  h5[["test_nullable"]]$create_dataset("values", robj = values)
  h5[["test_nullable"]]$create_dataset("mask", robj = mask)

  # Test flattening
  result <- FlattenNullable(h5[["test_nullable"]])

  expect_false(is.null(result))
  expect_equal(length(result), 5)
  expect_equal(result[1], 1.0)
  expect_true(is.na(result[2]))  # mask = TRUE -> NA
  expect_equal(result[3], 3.0)
  expect_equal(result[4], 4.0)
  expect_true(is.na(result[5]))  # mask = TRUE -> NA
})

test_that("FlattenNullable returns NULL for non-mask structures", {
  FlattenNullable <- srtdisk:::FlattenNullable

  temp_h5 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5), add = TRUE)

  h5 <- hdf5r::H5File$new(temp_h5, mode = "w")
  on.exit(h5$close_all(), add = TRUE, after = FALSE)

  # Create a categorical group (should not be flattened)
  h5$create_group("test_categorical")
  h5[["test_categorical"]]$create_dataset("categories", robj = c("A", "B", "C"))
  h5[["test_categorical"]]$create_dataset("codes", robj = c(0L, 1L, 2L))

  result <- FlattenNullable(h5[["test_categorical"]])
  expect_null(result)

  # Create a regular group without mask/values
  h5$create_group("test_regular")
  h5[["test_regular"]]$create_dataset("data", robj = c(1, 2, 3))

  result2 <- FlattenNullable(h5[["test_regular"]])
  expect_null(result2)
})

# -----------------------------------------------------------------------------
# Test integrated preprocessing in h5ad conversion
# -----------------------------------------------------------------------------

test_that("h5ad conversion handles special column names", {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("anndata"))
  skip_if_not(reticulate::py_module_available("numpy"))

  # Create a test h5ad file with problematic column names
  temp_h5ad <- tempfile(fileext = ".h5ad")
  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit({
    unlink(temp_h5ad)
    unlink(temp_h5seurat)
  }, add = TRUE)

  # Create h5ad with problematic obs column names using Python
  reticulate::py_run_string(sprintf('
import anndata as ad
import numpy as np
import pandas as pd

# Create simple expression matrix
X = np.random.poisson(5, (50, 100)).astype(np.float32)

# Create obs with problematic column names
obs = pd.DataFrame({
    "cell/type": ["A"] * 25 + ["B"] * 25,
    "sample name": ["S1"] * 50,
    "gene:count": np.random.randn(50),
    "1_numeric_start": np.arange(50)
}, index=[f"Cell{i}" for i in range(50)])

# Create var
var = pd.DataFrame({
    "gene/symbol": [f"Gene{i}" for i in range(100)]
}, index=[f"Gene{i}" for i in range(100)])

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.write("%s")
', temp_h5ad))

  # Convert h5ad to h5seurat
  expect_no_error(
    Convert(temp_h5ad, dest = temp_h5seurat, verbose = FALSE, overwrite = TRUE)
  )

  # Load and check the converted file
  h5 <- hdf5r::H5File$new(temp_h5seurat, mode = "r")
  on.exit(h5$close_all(), add = TRUE, after = FALSE)

  # Check that problematic column names were sanitized
  meta_cols <- names(h5[["meta.data"]])
  meta_cols <- meta_cols[!meta_cols %in% c("_index", "__categories")]

  # Expect sanitized names
  expect_true("cell_type" %in% meta_cols)
  expect_true("sample_name" %in% meta_cols)
  expect_true("gene_count" %in% meta_cols)
  expect_true("X1_numeric_start" %in% meta_cols)

  # Original problematic names should NOT be present
  expect_false("cell/type" %in% meta_cols)
  expect_false("sample name" %in% meta_cols)
})

test_that("h5ad conversion handles list columns", {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("anndata"))
  skip_if_not(reticulate::py_module_available("numpy"))
  skip_if_not(reticulate::py_module_available("h5py"))
  skip_if_not(reticulate::py_module_available("pandas"))

  # Create a test h5ad with list-type columns
  temp_h5ad <- tempfile(fileext = ".h5ad")
  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit({
    unlink(temp_h5ad)
    unlink(temp_h5seurat)
  }, add = TRUE)

  # Create h5ad with list columns using Python
  # Note: This creates a compound dataset obs with list-valued columns
  reticulate::py_run_string(sprintf('
import h5py
import numpy as np

# Create a minimal h5ad file structure
with h5py.File("%s", "w") as f:
    # Create X matrix (simple dense for testing)
    X = np.random.poisson(5, (10, 20)).astype(np.float32)
    f.create_dataset("X", data=X)

    # Create compound obs with a list-like column
    # Store as compound dataset which triggers Python fallback
    obs_dtype = np.dtype([
        ("index", "S20"),
        ("cluster", "i4"),
        ("tags", "S50")  # List will be serialized as string
    ])
    obs_data = np.array([
        (f"Cell{i}".encode(), i %% 3, "A;B;C".encode())
        for i in range(10)
    ], dtype=obs_dtype)
    f.create_dataset("obs", data=obs_data)

    # Create compound var
    var_dtype = np.dtype([
        ("index", "S20"),
        ("gene_name", "S20")
    ])
    var_data = np.array([
        (f"Gene{i}".encode(), f"Gene{i}".encode())
        for i in range(20)
    ], dtype=var_dtype)
    f.create_dataset("var", data=var_data)

    # Add required attributes
    f.attrs["encoding-type"] = "anndata"
    f.attrs["encoding-version"] = "0.1.0"
', temp_h5ad))

  # Convert should succeed
  expect_no_error(
    Convert(temp_h5ad, dest = temp_h5seurat, verbose = FALSE, overwrite = TRUE)
  )

  expect_true(file.exists(temp_h5seurat))
})

# -----------------------------------------------------------------------------
# Test with real datasets (if available)
# -----------------------------------------------------------------------------

test_that("Conversion preserves metadata from real h5ad files", {
  # Skip if test data not available
  test_h5ad <- system.file("testdata", "test_sample.h5ad", package = "srtdisk")
  skip_if(!file.exists(test_h5ad), "Test h5ad file not available")

  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(temp_h5seurat), add = TRUE)

  expect_no_error(
    Convert(test_h5ad, dest = temp_h5seurat, verbose = FALSE, overwrite = TRUE)
  )

  # Load converted file
  seurat_obj <- LoadH5Seurat(temp_h5seurat, verbose = FALSE)

  # Basic checks
  expect_s4_class(seurat_obj, "Seurat")
  expect_gt(ncol(seurat_obj), 0)
  expect_gt(nrow(seurat_obj), 0)

  # Check that metadata columns exist and are accessible
  meta_cols <- colnames(seurat_obj@meta.data)
  expect_true(length(meta_cols) > 0)
})

# -----------------------------------------------------------------------------
# Test round-trip preservation
# -----------------------------------------------------------------------------

test_that("Round-trip conversion preserves sanitized column names", {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("anndata"))

  # Create Seurat object with normal column names
  set.seed(42)
  counts <- matrix(
    data = rpois(100 * 50, lambda = 5),
    nrow = 100,
    ncol = 50,
    dimnames = list(
      paste0("Gene", seq_len(100)),
      paste0("Cell", seq_len(50))
    )
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = "TestProject")

  # Add metadata with valid names
  seurat_obj$cell_type <- factor(sample(c("A", "B", "C"), 50, replace = TRUE))
  seurat_obj$sample_id <- factor(rep(c("S1", "S2"), 25))
  seurat_obj$score <- rnorm(50)

  temp_h5seurat1 <- tempfile(fileext = ".h5seurat")
  temp_h5ad <- tempfile(fileext = ".h5ad")
  temp_h5seurat2 <- tempfile(fileext = ".h5seurat")
  on.exit({
    unlink(temp_h5seurat1)
    unlink(temp_h5ad)
    unlink(temp_h5seurat2)
  }, add = TRUE)

  # Save -> Convert to h5ad -> Convert back
  SaveH5Seurat(seurat_obj, filename = temp_h5seurat1, verbose = FALSE)
  Convert(temp_h5seurat1, dest = temp_h5ad, verbose = FALSE, overwrite = TRUE)
  Convert(temp_h5ad, dest = temp_h5seurat2, verbose = FALSE, overwrite = TRUE)

  # Load and compare
  loaded_obj <- LoadH5Seurat(temp_h5seurat2, verbose = FALSE)

  # Check metadata columns are preserved
  expect_true("cell_type" %in% colnames(loaded_obj@meta.data))
  expect_true("sample_id" %in% colnames(loaded_obj@meta.data))
  expect_true("score" %in% colnames(loaded_obj@meta.data))

  # Check values are preserved
  expect_equal(
    as.character(loaded_obj$cell_type),
    as.character(seurat_obj$cell_type)
  )
  expect_equal(
    loaded_obj$score,
    seurat_obj$score,
    tolerance = 1e-6
  )
})
