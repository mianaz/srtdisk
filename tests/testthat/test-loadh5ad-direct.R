# Tests for direct LoadH5AD functionality

library(srtdisk)

test_that("LoadH5AD function exists and is exported", {
  expect_true(exists("LoadH5AD"))
  expect_true(is.function(LoadH5AD))
})

test_that("LoadH5AD validates file existence", {
  expect_error(
    LoadH5AD("nonexistent_file.h5ad"),
    "File not found"
  )
})

test_that("LoadH5AD creates Seurat object from synthetic dense h5ad", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  # Create a synthetic h5ad file with dense matrix
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(42)
  n_cells <- 20
  n_genes <- 10
  counts <- matrix(rpois(n_cells * n_genes, lambda = 5), nrow = n_cells, ncol = n_genes)
  cell_names <- paste0("Cell_", seq_len(n_cells))
  gene_names <- paste0("Gene_", seq_len(n_genes))

  h5 <- hdf5r::H5File$new(tmp, mode = "w")

  # Write X matrix (cells x genes in h5ad)
  h5$create_dataset("X", robj = counts, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  # Write obs (cell metadata)
  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  obs$create_dataset("n_counts", robj = rowSums(counts),
                     dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  # Add encoding attributes for obs
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  obs$create_attr("encoding-version", robj = "0.2.0",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # Write var (feature metadata)
  var <- h5$create_group("var")
  var$create_dataset("_index", robj = gene_names)

  # Add encoding attributes for var
  var$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  var$create_attr("encoding-version", robj = "0.2.0",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  h5$close_all()

  # Load with LoadH5AD
  result <- LoadH5AD(tmp, verbose = FALSE)

  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), n_cells)
  expect_equal(nrow(result), n_genes)
  expect_identical(colnames(result), cell_names)
  # Note: Seurat v5 may replace underscores with dashes in feature names
  expect_equal(length(rownames(result)), length(gene_names))
})

test_that("LoadH5AD reads sparse h5ad correctly", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(123)
  n_cells <- 15
  n_genes <- 8

  # Create sparse data in CSR format (cells x genes)
  cell_names <- paste0("Cell_", seq_len(n_cells))
  gene_names <- paste0("Gene_", seq_len(n_genes))

  # Build a simple sparse matrix
  mat <- Matrix::rsparsematrix(nrow = n_cells, ncol = n_genes, density = 0.3)
  mat@x <- abs(mat@x) * 10  # Make values positive integers-ish

  # Convert to CSR for h5ad storage
  csr <- as(mat, "RsparseMatrix")

  h5 <- hdf5r::H5File$new(tmp, mode = "w")

  # Write sparse X in CSR format
  x_group <- h5$create_group("X")
  x_group$create_dataset("data", robj = csr@x, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  x_group$create_dataset("indices", robj = csr@j, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_group$create_dataset("indptr", robj = csr@p, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_group$create_attr("shape", robj = as.integer(c(n_cells, n_genes)),
                       dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_group$create_attr("encoding-type", robj = "csr_matrix",
                       dtype = hdf5r::H5T_STRING$new(size = Inf),
                       space = hdf5r::H5S$new(type = "scalar"))

  # Write obs
  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # Write var
  var <- h5$create_group("var")
  var$create_dataset("_index", robj = gene_names)
  var$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  h5$close_all()

  result <- LoadH5AD(tmp, verbose = FALSE)

  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), n_cells)
  expect_equal(nrow(result), n_genes)
  expect_identical(colnames(result), cell_names)
  # Note: Seurat v5 may replace underscores with dashes in feature names
  expect_equal(length(rownames(result)), length(gene_names))
})

test_that("LoadH5AD preserves categorical metadata", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  n_cells <- 10
  n_genes <- 5
  cell_names <- paste0("Cell_", seq_len(n_cells))
  gene_names <- paste0("Gene_", seq_len(n_genes))
  counts <- matrix(rpois(n_cells * n_genes, 3), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")

  # Write X
  h5$create_dataset("X", robj = counts, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  # Write obs with categorical column
  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)

  # Add categorical column (cluster)
  categories <- c("TypeA", "TypeB", "TypeC")
  codes <- as.integer(c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
  obs$create_dataset("cluster", robj = codes, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  cats_group <- obs$create_group("__categories")
  cats_group$create_dataset("cluster", robj = categories)

  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # Write var
  var <- h5$create_group("var")
  var$create_dataset("_index", robj = gene_names)
  var$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  h5$close_all()

  result <- LoadH5AD(tmp, verbose = FALSE)

  # Check categorical was preserved as factor
  expect_true("cluster" %in% colnames(result@meta.data))
  expect_s3_class(result$cluster, "factor")
  expect_equal(levels(result$cluster), categories)
  expect_equal(as.character(result$cluster), categories[codes + 1])
})

test_that("LoadH5AD preserves dimensional reductions", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  n_cells <- 10
  n_genes <- 5
  cell_names <- paste0("Cell_", seq_len(n_cells))
  gene_names <- paste0("Gene_", seq_len(n_genes))
  counts <- matrix(rpois(n_cells * n_genes, 3), nrow = n_cells, ncol = n_genes)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")

  # Write X
  h5$create_dataset("X", robj = counts, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  # Write obs
  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # Write var
  var <- h5$create_group("var")
  var$create_dataset("_index", robj = gene_names)
  var$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  # Write obsm with PCA and UMAP
  obsm <- h5$create_group("obsm")
  pca_embed <- matrix(rnorm(n_cells * 5), nrow = n_cells, ncol = 5)
  obsm$create_dataset("X_pca", robj = pca_embed, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  umap_embed <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)
  obsm$create_dataset("X_umap", robj = umap_embed, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  h5$close_all()

  result <- LoadH5AD(tmp, verbose = FALSE)

  # Check reductions
  expect_true("pca" %in% names(result@reductions))
  expect_true("umap" %in% names(result@reductions))

  # Check PCA dimensions
  pca <- Seurat::Embeddings(result, "pca")
  expect_equal(nrow(pca), n_cells)
  expect_equal(ncol(pca), 5)

  # Check UMAP dimensions
  umap <- Seurat::Embeddings(result, "umap")
  expect_equal(nrow(umap), n_cells)
  expect_equal(ncol(umap), 2)
})

test_that("LoadH5AD handles missing X matrix gracefully", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = c("Cell1", "Cell2"))
  h5$close_all()

  # Should error because there's no X matrix (may fail at different points)
  expect_error(LoadH5AD(tmp, verbose = FALSE))
})
