# Comprehensive tests for LoadH5AD and Seurat <-> h5ad roundtrip
#
# Tests cover: CSR/CSC/dense matrices, categorical metadata (modern + legacy),
# numeric/string metadata, feature metadata, deduplication, dimensional
# reductions, neighbor graphs, raw/X semantics, layers, _index attribute
# convention, dimension mismatch warnings, roundtrip conversion, real data,
# empty sparse matrices, multi-assay, and spatial coordinates.

library(srtdisk)

skip_if_not_installed("hdf5r")
skip_if_not_installed("Seurat")
skip_if_not_installed("Matrix")

# ---------------------------------------------------------------------------
# Shared helper: write a minimal h5ad frame (obs + var groups with _index)
# ---------------------------------------------------------------------------

write_h5ad_frame <- function(h5, group_name, index_values, columns = list()) {
  grp <- h5$create_group(group_name)
  grp$create_dataset("_index", robj = index_values)
  grp$create_attr(
    "encoding-type", robj = "dataframe",
    dtype = hdf5r::H5T_STRING$new(size = Inf),
    space = hdf5r::H5S$new(type = "scalar")
  )
  grp$create_attr(
    "encoding-version", robj = "0.2.0",
    dtype = hdf5r::H5T_STRING$new(size = Inf),
    space = hdf5r::H5S$new(type = "scalar")
  )
  if (length(columns) > 0) {
    col_order <- names(columns)
    grp$create_attr(
      "column-order", robj = col_order,
      dtype = hdf5r::H5T_STRING$new(size = Inf)
    )
    for (nm in col_order) {
      val <- columns[[nm]]
      if (is.numeric(val) || is.logical(val) || is.character(val)) {
        dtype <- if (is.double(val)) {
          hdf5r::h5types$H5T_NATIVE_DOUBLE
        } else if (is.integer(val)) {
          hdf5r::h5types$H5T_NATIVE_INT
        } else if (is.logical(val)) {
          hdf5r::h5types$H5T_NATIVE_INT
        } else {
          hdf5r::H5T_STRING$new(size = Inf)
        }
        grp$create_dataset(nm, robj = val, dtype = dtype)
      }
    }
  }
  invisible(grp)
}

# Helper: write a CSR sparse matrix group to h5ad
write_csr_group <- function(parent, name, mat_csr, n_rows, n_cols) {
  grp <- parent$create_group(name)
  grp$create_dataset("data", robj = mat_csr@x, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  grp$create_dataset("indices", robj = mat_csr@j, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  grp$create_dataset("indptr", robj = mat_csr@p, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  grp$create_attr("shape", robj = as.integer(c(n_rows, n_cols)),
                   dtype = hdf5r::h5types$H5T_NATIVE_INT)
  grp$create_attr("encoding-type", robj = "csr_matrix",
                   dtype = hdf5r::H5T_STRING$new(size = Inf),
                   space = hdf5r::H5S$new(type = "scalar"))
  grp$create_attr("encoding-version", robj = "0.1.0",
                   dtype = hdf5r::H5T_STRING$new(size = Inf),
                   space = hdf5r::H5S$new(type = "scalar"))
  invisible(grp)
}

# ---------------------------------------------------------------------------
# 1. CSR sparse matrix
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads CSR sparse matrix correctly", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(1)
  n_cells <- 20
  n_genes <- 15
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))

  # Build dense, convert to dgCMatrix (CSC), then to RsparseMatrix (CSR)
  dense <- matrix(rpois(n_cells * n_genes, lambda = 2), nrow = n_cells, ncol = n_genes)
  csc <- Matrix::Matrix(dense, sparse = TRUE)  # dgCMatrix
  csr <- as(csc, "RsparseMatrix")

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  write_csr_group(h5, "X", csr, n_cells, n_genes)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes)
  expect_identical(colnames(obj), cell_names)

  # Verify non-zero values survive the round-trip
  loaded_mat <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
  # h5ad stores cells-x-genes, LoadH5AD transposes to genes-x-cells
  expect_equal(sum(loaded_mat), sum(dense))
})

# ---------------------------------------------------------------------------
# 2. CSC sparse matrix
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads CSC sparse matrix with correct transposition", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(2)
  n_cells <- 12
  n_genes <- 8
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- paste0("g", seq_len(n_genes))

  dense <- matrix(rpois(n_cells * n_genes, lambda = 3), nrow = n_cells, ncol = n_genes)
  csc <- as(Matrix::Matrix(dense, sparse = TRUE), "dgCMatrix")  # cells x genes CSC

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  x_grp <- h5$create_group("X")
  x_grp$create_dataset("data", robj = csc@x, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  x_grp$create_dataset("indices", robj = csc@i, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_dataset("indptr", robj = csc@p, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_attr("shape", robj = as.integer(c(n_cells, n_genes)),
                     dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_attr("encoding-type", robj = "csc_matrix",
                     dtype = hdf5r::H5T_STRING$new(size = Inf),
                     space = hdf5r::H5S$new(type = "scalar"))
  x_grp$create_attr("encoding-version", robj = "0.1.0",
                     dtype = hdf5r::H5T_STRING$new(size = Inf),
                     space = hdf5r::H5S$new(type = "scalar"))

  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes)

  loaded <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
  # genes-x-cells in Seurat => compare to transpose of dense (cells-x-genes)
  expect_equal(loaded[, 1], dense[1, ], ignore_attr = TRUE)
  expect_equal(sum(loaded), sum(dense))
})

# ---------------------------------------------------------------------------
# 3. Dense matrix
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads dense matrix with auto-transpose detection", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(3)
  n_cells <- 10
  n_genes <- 6
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, lambda = 4), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  # h5ad stores cells-x-genes in row-major; hdf5r writes column-major,
  # so we store the matrix as-is and LoadH5AD should handle the transpose
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes)
  expect_equal(sum(as.matrix(Seurat::GetAssayData(obj, layer = "counts"))), sum(dense))
})

# ---------------------------------------------------------------------------
# 4. Modern categorical metadata (codes/categories groups in obs)
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads modern categorical metadata (codes/categories groups)", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(4)
  n_cells <- 15
  n_genes <- 5
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 3), nrow = n_cells, ncol = n_genes)

  categories <- c("T-cell", "B-cell", "Monocyte")
  codes <- as.integer(c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2))

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "var", gene_names)

  # obs with modern categorical

  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  obs$create_attr("encoding-version", robj = "0.2.0",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # Create categorical group for cell_type
  ct <- obs$create_group("cell_type")
  ct$create_dataset("codes", robj = codes, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  ct$create_dataset("categories", robj = categories)
  ct$create_attr("encoding-type", robj = "categorical",
                 dtype = hdf5r::H5T_STRING$new(size = Inf),
                 space = hdf5r::H5S$new(type = "scalar"))
  ct$create_attr("encoding-version", robj = "0.2.0",
                 dtype = hdf5r::H5T_STRING$new(size = Inf),
                 space = hdf5r::H5S$new(type = "scalar"))

  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_true("cell_type" %in% colnames(obj@meta.data))
  expect_s3_class(obj$cell_type, "factor")
  expect_equal(levels(obj$cell_type), categories)
  expect_equal(as.character(obj$cell_type), categories[codes + 1L])
})

# ---------------------------------------------------------------------------
# 5. Modern categorical with NA values (codes = -1)
# ---------------------------------------------------------------------------

test_that("LoadH5AD handles categorical columns with NA codes (-1)", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(50)
  n_cells <- 8
  n_genes <- 4
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- paste0("g", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  categories <- c("A", "B")
  codes <- as.integer(c(0, 1, -1, 0, -1, 1, 0, 1))  # -1 means NA

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "var", gene_names)

  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  grp <- obs$create_group("status")
  grp$create_dataset("codes", robj = codes, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  grp$create_dataset("categories", robj = categories)
  grp$create_attr("encoding-type", robj = "categorical",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  grp$create_attr("encoding-version", robj = "0.2.0",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_true("status" %in% colnames(obj@meta.data))
  expect_s3_class(obj$status, "factor")
  expect_true(is.na(obj$status[3]))
  expect_true(is.na(obj$status[5]))
  expect_equal(as.character(obj$status[1]), "A")
  expect_equal(as.character(obj$status[2]), "B")
})

# ---------------------------------------------------------------------------
# 5b. Legacy categorical metadata (__categories group)
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads legacy categorical metadata (__categories)", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(51)
  n_cells <- 10
  n_genes <- 5
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 3), nrow = n_cells, ncol = n_genes)

  categories <- c("TypeA", "TypeB", "TypeC")
  codes <- as.integer(c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0))

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "var", gene_names)

  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # Legacy: codes stored as dataset, categories in __categories group
  obs$create_dataset("cluster", robj = codes, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  cats <- obs$create_group("__categories")
  cats$create_dataset("cluster", robj = categories)

  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_true("cluster" %in% colnames(obj@meta.data))
  expect_s3_class(obj$cluster, "factor")
  expect_equal(levels(obj$cluster), categories)
  expect_equal(as.character(obj$cluster), categories[codes + 1L])
})

# ---------------------------------------------------------------------------
# 6. Numeric and string metadata columns
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads numeric and string metadata columns", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(6)
  n_cells <- 10
  n_genes <- 5
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  n_counts <- as.double(rowSums(dense))
  batch_labels <- rep(c("batch1", "batch2"), each = 5)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "var", gene_names)
  write_h5ad_frame(h5, "obs", cell_names, columns = list(
    n_counts = n_counts,
    batch = batch_labels
  ))
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_true("n_counts" %in% colnames(obj@meta.data))
  expect_true("batch" %in% colnames(obj@meta.data))
  expect_equal(as.numeric(obj$n_counts), n_counts, ignore_attr = TRUE)
  expect_equal(as.character(obj$batch), batch_labels, ignore_attr = TRUE)
  expect_true(is.numeric(obj$n_counts))
  expect_true(is.character(obj$batch))
})

# ---------------------------------------------------------------------------
# 7. Feature metadata (var) including highly_variable boolean
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads feature metadata and sets VariableFeatures from highly_variable", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(7)
  n_cells <- 10
  n_genes <- 20
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 3), nrow = n_cells, ncol = n_genes)

  # 5 highly variable genes
  hv <- c(rep(1L, 5), rep(0L, 15))
  total_counts <- as.double(colSums(dense))

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)

  var_grp <- h5$create_group("var")
  var_grp$create_dataset("_index", robj = gene_names)
  var_grp$create_attr("encoding-type", robj = "dataframe",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))
  var_grp$create_attr("encoding-version", robj = "0.2.0",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))
  var_grp$create_dataset("highly_variable", robj = hv,
                         dtype = hdf5r::h5types$H5T_NATIVE_INT)
  var_grp$create_dataset("total_counts", robj = total_counts,
                         dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  # Check variable features were set
  hvf <- Seurat::VariableFeatures(obj)
  expect_equal(length(hvf), 5)
  expect_true(all(hvf %in% gene_names[1:5]))
})

# ---------------------------------------------------------------------------
# 8. Feature name deduplication
# ---------------------------------------------------------------------------

test_that("LoadH5AD deduplicates feature names", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(8)
  n_cells <- 6
  n_genes <- 4
  cell_names <- paste0("c", seq_len(n_cells))
  # Duplicate gene names
  gene_names <- c("TP53", "BRCA1", "TP53", "EGFR")
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  # make.unique should produce TP53, BRCA1, TP53.1, EGFR
  expect_equal(nrow(obj), n_genes)
  rn <- rownames(obj)
  expect_true(all(c("TP53", "BRCA1", "EGFR") %in% rn))
  expect_equal(length(unique(rn)), n_genes)
  # The second TP53 should be TP53.1
  expect_true("TP53.1" %in% rn)
})

# ---------------------------------------------------------------------------
# 9. Dimensional reductions (PCA 5D + UMAP 2D) in obsm
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads PCA and UMAP from obsm with transpose detection", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(9)
  n_cells <- 15
  n_genes <- 8
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 3), nrow = n_cells, ncol = n_genes)

  pca_embed <- matrix(rnorm(n_cells * 5), nrow = n_cells, ncol = 5)
  umap_embed <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)

  obsm <- h5$create_group("obsm")
  obsm$create_dataset("X_pca", robj = pca_embed, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  obsm$create_dataset("X_umap", robj = umap_embed, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_true("pca" %in% names(obj@reductions))
  expect_true("umap" %in% names(obj@reductions))

  pca_loaded <- Seurat::Embeddings(obj, "pca")
  expect_equal(nrow(pca_loaded), n_cells)
  expect_equal(ncol(pca_loaded), 5)

  umap_loaded <- Seurat::Embeddings(obj, "umap")
  expect_equal(nrow(umap_loaded), n_cells)
  expect_equal(ncol(umap_loaded), 2)

  # Key naming
  expect_equal(Seurat::Key(obj[["pca"]]), "PC_")
  expect_equal(Seurat::Key(obj[["umap"]]), "UMAP_")
})

# ---------------------------------------------------------------------------
# 10. Neighbor graphs (connectivities in obsp)
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads neighbor graphs from obsp", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(10)
  n_cells <- 12
  n_genes <- 5
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 3), nrow = n_cells, ncol = n_genes)

  # Create a sparse connectivity matrix (cells x cells)
  conn_dense <- matrix(0, nrow = n_cells, ncol = n_cells)
  for (i in seq_len(n_cells)) {
    # Each cell connected to 3 nearest neighbors
    neighbors <- sample(setdiff(seq_len(n_cells), i), min(3, n_cells - 1))
    conn_dense[i, neighbors] <- runif(length(neighbors))
  }
  conn_csr <- as(Matrix::Matrix(conn_dense, sparse = TRUE), "RsparseMatrix")

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)

  obsp <- h5$create_group("obsp")
  write_csr_group(obsp, "connectivities", conn_csr, n_cells, n_cells)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_true("RNA_snn" %in% names(obj@graphs))
  graph <- obj@graphs[["RNA_snn"]]
  expect_equal(nrow(graph), n_cells)
  expect_equal(ncol(graph), n_cells)
})

# ---------------------------------------------------------------------------
# 11. raw/X layer semantics
# ---------------------------------------------------------------------------

test_that("LoadH5AD maps X to data and raw/X to counts when both exist", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(11)
  n_cells <- 10
  n_genes <- 6
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))

  raw_counts <- matrix(rpois(n_cells * n_genes, lambda = 5), nrow = n_cells, ncol = n_genes)
  # "Normalized" data: log1p of raw counts
  normalized <- log1p(raw_counts)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")

  # X = normalized data (cells x genes, dense)
  h5$create_dataset("X", robj = normalized, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  # raw/X = raw counts (CSR)
  raw_grp <- h5$create_group("raw")
  raw_csr <- as(Matrix::Matrix(raw_counts, sparse = TRUE), "RsparseMatrix")
  write_csr_group(raw_grp, "X", raw_csr, n_cells, n_genes)

  # raw/var with _index (feature names)
  raw_var <- raw_grp$create_group("var")
  raw_var$create_dataset("_index", robj = gene_names)
  raw_var$create_attr("encoding-type", robj = "dataframe",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))
  raw_var$create_attr("encoding-version", robj = "0.2.0",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))

  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  # Verify counts layer has raw counts
  counts_mat <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
  expect_equal(sum(counts_mat), sum(raw_counts))

  # Verify data layer has normalized values
  data_mat <- as.matrix(Seurat::GetAssayData(obj, layer = "data"))
  expect_equal(sum(data_mat), sum(normalized), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 12. Multiple layers (layers/counts, layers/log_normalized)
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads multiple named layers", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(12)
  n_cells <- 10
  n_genes <- 6
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  raw_counts <- matrix(rpois(n_cells * n_genes, lambda = 3), nrow = n_cells, ncol = n_genes)
  log_norm <- log1p(raw_counts / rowSums(raw_counts) * 1e4)

  # X = raw counts
  raw_csr <- as(Matrix::Matrix(raw_counts, sparse = TRUE), "RsparseMatrix")

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  write_csr_group(h5, "X", raw_csr, n_cells, n_genes)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)

  # layers
  layers <- h5$create_group("layers")
  log_norm_csr <- as(Matrix::Matrix(log_norm, sparse = TRUE), "RsparseMatrix")
  write_csr_group(layers, "log_normalized", log_norm_csr, n_cells, n_genes)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  # The main X goes to counts; log_normalized maps to data layer
  expect_true("data" %in% SeuratObject::Layers(obj[["RNA"]]))
})

# ---------------------------------------------------------------------------
# 13. _index attribute convention (attribute pointing to a column name)
# ---------------------------------------------------------------------------

test_that("LoadH5AD resolves _index attribute naming the index column", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(13)
  n_cells <- 8
  n_genes <- 4
  cell_names <- paste0("barcode_", seq_len(n_cells))
  gene_names <- paste0("symbol_", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  # obs: _index attribute points to "barcodes" column
  obs <- h5$create_group("obs")
  obs$create_attr("_index", robj = "barcodes",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  obs$create_dataset("barcodes", robj = cell_names)

  # var: _index attribute points to "gene_symbols" column
  var_grp <- h5$create_group("var")
  var_grp$create_attr("_index", robj = "gene_symbols",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))
  var_grp$create_attr("encoding-type", robj = "dataframe",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))
  var_grp$create_dataset("gene_symbols", robj = gene_names)

  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes)
  expect_identical(colnames(obj), cell_names)
})

# ---------------------------------------------------------------------------
# 14. Legacy compound dataset format (obs as H5D compound)
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads legacy compound dataset obs", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(14)
  n_cells <- 8
  n_genes <- 4
  cell_names <- paste0("cell", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  # Build a compound data.frame for obs
  obs_df <- data.frame(
    index = cell_names,
    n_counts = as.double(rowSums(dense)),
    batch = rep(c("A", "B"), each = 4),
    stringsAsFactors = FALSE
  )

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  h5$create_dataset("obs", robj = obs_df)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_equal(ncol(obj), n_cells)
  expect_true("n_counts" %in% colnames(obj@meta.data))
  expect_true("batch" %in% colnames(obj@meta.data))
  expect_equal(as.numeric(obj$n_counts), as.double(rowSums(dense)), ignore_attr = TRUE)
})

# ---------------------------------------------------------------------------
# 15. Dimension mismatch warning
# ---------------------------------------------------------------------------

test_that("LoadH5AD warns on feature count mismatch and truncates names", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(15)
  n_cells <- 8
  n_genes_matrix <- 4  # actual matrix columns
  n_genes_names <- 6   # more names than matrix columns => triggers warning
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- paste0("g", seq_len(n_genes_names))
  dense <- matrix(rpois(n_cells * n_genes_matrix, 2), nrow = n_cells, ncol = n_genes_matrix)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  expect_warning(
    obj <- LoadH5AD(tmp, verbose = FALSE),
    "mismatch"
  )
  # Feature names should be truncated to match matrix dimensions
  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes_matrix)
})

# ---------------------------------------------------------------------------
# 16. Seurat -> h5ad -> Seurat roundtrip
# ---------------------------------------------------------------------------

test_that("Seurat -> h5ad -> Seurat roundtrip preserves data", {
  set.seed(16)
  n_cells <- 50
  n_genes <- 30

  counts <- Matrix::Matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells,
    sparse = TRUE,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Cell", seq_len(n_cells))
    )
  )
  srt <- Seurat::CreateSeuratObject(counts = counts, project = "RoundTrip",
                                     min.cells = 0, min.features = 0)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$cluster <- factor(sample(c("A", "B", "C"), n_cells, replace = TRUE))
  srt$score <- rnorm(n_cells)

  dest <- tempfile(fileext = ".h5ad")
  on.exit(unlink(dest), add = TRUE)

  Convert(srt, dest = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(dest))

  # Read back
  loaded <- LoadH5AD(dest, verbose = FALSE)

  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), n_cells)
  expect_equal(nrow(loaded), n_genes)

  # Check expression sum is preserved (within tolerance for float conversion)
  orig_sum <- sum(Seurat::GetAssayData(srt, layer = "counts"))
  loaded_sum <- sum(Seurat::GetAssayData(loaded, layer = "counts"))
  expect_equal(loaded_sum, orig_sum, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 17. Real data: crc_sample.h5ad
# ---------------------------------------------------------------------------

test_that("LoadH5AD loads crc_sample.h5ad with correct structure", {
  crc_path <- system.file("testdata", "crc_sample.h5ad", package = "srtdisk")
  skip_if(crc_path == "" || !file.exists(crc_path), "crc_sample.h5ad not available")

  obj <- LoadH5AD(crc_path, verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  # 935 cells, 25344 genes (from shape attribute)
  expect_equal(ncol(obj), 935)
  expect_equal(nrow(obj), 25344)

  # Check some known obs columns exist
  md_cols <- colnames(obj@meta.data)
  expected_cols <- c("cell_type", "disease", "tissue", "sex", "donor_id")
  for (col in expected_cols) {
    expect_true(col %in% md_cols, info = paste("Missing metadata column:", col))
  }

  # cell_type should be a factor (categorical in h5ad)
  expect_s3_class(obj$cell_type, "factor")

  # UMAP should be present
  expect_true("umap" %in% names(obj@reductions))
  umap_embed <- Seurat::Embeddings(obj, "umap")
  expect_equal(nrow(umap_embed), 935)
  expect_equal(ncol(umap_embed), 2)
})

# ---------------------------------------------------------------------------
# 18. Empty sparse matrix edge case
# ---------------------------------------------------------------------------

test_that("LoadH5AD handles empty sparse matrix (0 non-zero elements)", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  n_cells <- 5
  n_genes <- 3
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- paste0("g", seq_len(n_genes))

  h5 <- hdf5r::H5File$new(tmp, mode = "w")

  # Empty CSR: no data, all-zero indptr
  x_grp <- h5$create_group("X")
  x_grp$create_dataset("data", robj = numeric(0), dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  x_grp$create_dataset("indices", robj = integer(0), dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_dataset("indptr", robj = rep(0L, n_cells + 1L), dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_attr("shape", robj = as.integer(c(n_cells, n_genes)),
                     dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_attr("encoding-type", robj = "csr_matrix",
                     dtype = hdf5r::H5T_STRING$new(size = Inf),
                     space = hdf5r::H5S$new(type = "scalar"))
  x_grp$create_attr("encoding-version", robj = "0.1.0",
                     dtype = hdf5r::H5T_STRING$new(size = Inf),
                     space = hdf5r::H5S$new(type = "scalar"))

  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes)
  expect_equal(sum(Seurat::GetAssayData(obj, layer = "counts")), 0)
})

# ---------------------------------------------------------------------------
# 19. CITE-seq-like multi-assay roundtrip (RNA + ADT -> h5ad of RNA assay)
# ---------------------------------------------------------------------------

test_that("Multi-assay Seurat converts RNA assay to h5ad and loads back", {
  set.seed(19)
  n_cells <- 20
  n_rna <- 50
  n_adt <- 10

  rna_counts <- Matrix::Matrix(
    rpois(n_rna * n_cells, lambda = 3),
    nrow = n_rna, ncol = n_cells, sparse = TRUE,
    dimnames = list(paste0("Gene", seq_len(n_rna)), paste0("Cell", seq_len(n_cells)))
  )
  adt_counts <- Matrix::Matrix(
    rpois(n_adt * n_cells, lambda = 100),
    nrow = n_adt, ncol = n_cells, sparse = TRUE,
    dimnames = list(paste0("ADT", seq_len(n_adt)), paste0("Cell", seq_len(n_cells)))
  )

  srt <- Seurat::CreateSeuratObject(counts = rna_counts, assay = "RNA",
                                     min.cells = 0, min.features = 0)
  srt[["ADT"]] <- Seurat::CreateAssayObject(counts = adt_counts)

  dest <- tempfile(fileext = ".h5ad")
  on.exit(unlink(dest), add = TRUE)

  # Convert RNA assay
  SeuratToH5AD(srt, filename = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(dest))

  loaded <- LoadH5AD(dest, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), n_cells)
  # Should have RNA features
  expect_equal(nrow(loaded), n_rna)

  # Verify RNA expression values
  orig_sum <- sum(Seurat::GetAssayData(srt, assay = "RNA", layer = "counts"))
  loaded_sum <- sum(Seurat::GetAssayData(loaded, layer = "counts"))
  expect_equal(loaded_sum, orig_sum, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 20. Spatial coordinates in obsm
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads spatial coordinates from obsm/spatial", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(20)
  n_cells <- 12
  n_genes <- 5
  cell_names <- paste0("spot", seq_len(n_cells))
  gene_names <- paste0("gene", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 3), nrow = n_cells, ncol = n_genes)

  spatial_coords <- matrix(runif(n_cells * 2, 0, 1000), nrow = n_cells, ncol = 2)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)

  obsm <- h5$create_group("obsm")
  obsm$create_dataset("spatial", robj = spatial_coords,
                       dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  # spatial should be available either as a reduction or via the spatial handler
  # LoadH5AD stores obsm entries as reductions (stripping X_ prefix if any)
  # "spatial" doesn't have X_ prefix, so it's stored as-is
  expect_true("spatial" %in% names(obj@reductions))
  spatial_embed <- Seurat::Embeddings(obj, "spatial")
  expect_equal(nrow(spatial_embed), n_cells)
  expect_equal(ncol(spatial_embed), 2)
})

# ---------------------------------------------------------------------------
# Additional tests: assay name parameter
# ---------------------------------------------------------------------------

test_that("LoadH5AD respects custom assay.name parameter", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(99)
  n_cells <- 8
  n_genes <- 4
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- paste0("g", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, assay.name = "Spatial", verbose = FALSE)

  expect_true("Spatial" %in% names(obj@assays))
  expect_equal(Seurat::DefaultAssay(obj), "Spatial")
})

# ---------------------------------------------------------------------------
# Additional: multiple categoricals in obs
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads multiple categorical columns correctly", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(100)
  n_cells <- 12
  n_genes <- 4
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- paste0("g", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "var", gene_names)

  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # First categorical: cell_type
  ct_cats <- c("T-cell", "B-cell", "NK")
  ct_codes <- as.integer(rep(0:2, 4))
  ct <- obs$create_group("cell_type")
  ct$create_dataset("codes", robj = ct_codes, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  ct$create_dataset("categories", robj = ct_cats)
  ct$create_attr("encoding-type", robj = "categorical",
                 dtype = hdf5r::H5T_STRING$new(size = Inf),
                 space = hdf5r::H5S$new(type = "scalar"))
  ct$create_attr("encoding-version", robj = "0.2.0",
                 dtype = hdf5r::H5T_STRING$new(size = Inf),
                 space = hdf5r::H5S$new(type = "scalar"))

  # Second categorical: tissue
  tis_cats <- c("blood", "spleen")
  tis_codes <- as.integer(rep(0:1, 6))
  tis <- obs$create_group("tissue")
  tis$create_dataset("codes", robj = tis_codes, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  tis$create_dataset("categories", robj = tis_cats)
  tis$create_attr("encoding-type", robj = "categorical",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  tis$create_attr("encoding-version", robj = "0.2.0",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # Also add a numeric column alongside categoricals
  obs$create_dataset("score", robj = rnorm(n_cells),
                     dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_s3_class(obj$cell_type, "factor")
  expect_equal(levels(obj$cell_type), ct_cats)

  expect_s3_class(obj$tissue, "factor")
  expect_equal(levels(obj$tissue), tis_cats)

  expect_true(is.numeric(obj$score))
})

# ---------------------------------------------------------------------------
# Additional: SeuratToH5AD convenience wrapper roundtrip
# ---------------------------------------------------------------------------

test_that("SeuratToH5AD convenience wrapper produces loadable h5ad", {
  set.seed(200)
  n_cells <- 25
  n_genes <- 15

  counts <- Matrix::Matrix(
    rpois(n_genes * n_cells, lambda = 4),
    nrow = n_genes, ncol = n_cells, sparse = TRUE,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Cell", seq_len(n_cells))
    )
  )
  srt <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)

  dest <- tempfile(fileext = ".h5ad")
  on.exit(unlink(dest), add = TRUE)

  result <- SeuratToH5AD(srt, filename = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(dest))
  expect_equal(result, dest)

  loaded <- LoadH5AD(dest, verbose = FALSE)
  expect_equal(ncol(loaded), n_cells)
  expect_equal(nrow(loaded), n_genes)
})

# ---------------------------------------------------------------------------
# Additional: var with _index attribute (not _index dataset)
# ---------------------------------------------------------------------------

test_that("LoadH5AD resolves var _index attribute correctly", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(300)
  n_cells <- 6
  n_genes <- 3
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- c("TP53", "EGFR", "KRAS")
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)

  # var with _index attribute pointing to "gene_symbols"
  var_grp <- h5$create_group("var")
  var_grp$create_attr("_index", robj = "gene_symbols",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))
  var_grp$create_attr("encoding-type", robj = "dataframe",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))
  var_grp$create_dataset("gene_symbols", robj = gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_equal(nrow(obj), n_genes)
  expect_true(all(c("TP53", "EGFR", "KRAS") %in% rownames(obj)))
})

# ---------------------------------------------------------------------------
# Additional: file not found error
# ---------------------------------------------------------------------------

test_that("LoadH5AD errors on non-existent file", {
  expect_error(LoadH5AD("does_not_exist.h5ad"), "File not found")
})

# ---------------------------------------------------------------------------
# Additional: no X matrix error
# ---------------------------------------------------------------------------

test_that("LoadH5AD errors when no X matrix is present", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  write_h5ad_frame(h5, "obs", c("c1", "c2"))
  h5$close_all()

  expect_error(LoadH5AD(tmp, verbose = FALSE))
})

# ---------------------------------------------------------------------------
# Additional: generated cell/feature names when obs/var are missing
# ---------------------------------------------------------------------------

test_that("LoadH5AD generates cell and feature names when obs/var are absent", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(400)
  n_cells <- 5
  n_genes <- 3
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes)
  # Auto-generated names
  expect_true(all(grepl("^Cell", colnames(obj))))
  expect_true(all(grepl("^Gene", rownames(obj))))
})

# ---------------------------------------------------------------------------
# Additional: roundtrip with metadata and reductions
# ---------------------------------------------------------------------------

test_that("Seurat -> h5ad -> Seurat roundtrip preserves metadata and structure", {
  set.seed(500)
  n_cells <- 40
  n_genes <- 20

  counts <- Matrix::Matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells, sparse = TRUE,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Cell", seq_len(n_cells))
    )
  )
  srt <- Seurat::CreateSeuratObject(counts = counts, project = "TestRT",
                                     min.cells = 0, min.features = 0)
  srt$condition <- factor(sample(c("Ctrl", "Treat"), n_cells, replace = TRUE))
  srt$nUMI <- as.double(Matrix::colSums(counts))

  dest <- tempfile(fileext = ".h5ad")
  on.exit(unlink(dest), add = TRUE)

  Convert(srt, dest = dest, overwrite = TRUE, verbose = FALSE)
  loaded <- LoadH5AD(dest, verbose = FALSE)

  # Cell count preserved
  expect_equal(ncol(loaded), n_cells)
  # Feature count preserved
  expect_equal(nrow(loaded), n_genes)
  # Metadata columns should include condition
  expect_true("condition" %in% colnames(loaded@meta.data))
})

# ---------------------------------------------------------------------------
# Additional: CSR with explicit encoding-version attribute
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads CSR with all standard attributes", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(600)
  n_cells <- 8
  n_genes <- 5
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- paste0("g", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)
  csr <- as(Matrix::Matrix(dense, sparse = TRUE), "RsparseMatrix")

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  write_csr_group(h5, "X", csr, n_cells, n_genes)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)
  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes)
  expect_equal(sum(as.matrix(Seurat::GetAssayData(obj, layer = "counts"))), sum(dense))
})

# ---------------------------------------------------------------------------
# Additional: uns data goes to misc
# ---------------------------------------------------------------------------

test_that("LoadH5AD stores uns scalar datasets in misc slot", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(700)
  n_cells <- 6
  n_genes <- 3
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- paste0("g", seq_len(n_genes))
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  write_h5ad_frame(h5, "obs", cell_names)
  write_h5ad_frame(h5, "var", gene_names)

  uns <- h5$create_group("uns")
  uns$create_dataset("schema_version", robj = "3.0.0")
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_true("schema_version" %in% names(obj@misc))
  expect_equal(as.character(obj@misc$schema_version), "3.0.0")
})

# ---------------------------------------------------------------------------
# Additional: var with "index" dataset (no underscore prefix)
# ---------------------------------------------------------------------------

test_that("LoadH5AD reads var with 'index' dataset (no underscore prefix)", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(800)
  n_cells <- 6
  n_genes <- 3
  cell_names <- paste0("c", seq_len(n_cells))
  gene_names <- c("GeneA", "GeneB", "GeneC")
  dense <- matrix(rpois(n_cells * n_genes, 2), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("X", robj = dense, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  # obs with "index" (no underscore)
  obs <- h5$create_group("obs")
  obs$create_dataset("index", robj = cell_names)
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # var with "index" (no underscore)
  var_grp <- h5$create_group("var")
  var_grp$create_dataset("index", robj = gene_names)
  var_grp$create_attr("encoding-type", robj = "dataframe",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))
  h5$close_all()

  obj <- LoadH5AD(tmp, verbose = FALSE)

  expect_equal(ncol(obj), n_cells)
  expect_equal(nrow(obj), n_genes)
  expect_identical(colnames(obj), cell_names)
  expect_true(all(c("GeneA", "GeneB", "GeneC") %in% rownames(obj)))
})
