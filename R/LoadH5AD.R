#' Load an AnnData H5AD file as a Seurat object
#'
#' Direct conversion from H5AD format to Seurat object without intermediate h5Seurat.
#' Supports optional BPCells on-disk matrix loading for large datasets that exceed
#' available memory.
#'
#' @param file Path to H5AD file
#' @param assay.name Name for the primary assay (default: "RNA")
#' @param use.bpcells If not NULL, a directory path where BPCells will store the
#'   expression matrix on disk. Requires the BPCells package. The resulting Seurat
#'   object will reference the on-disk matrix instead of loading it into memory,
#'   enabling analysis of datasets larger than available RAM.
#' @param verbose Show progress messages
#'
#' @return A \code{Seurat} object. If \code{use.bpcells} is set, the count matrix
#'   is stored on disk in BPCells format and the object uses minimal memory.
#'
#' @importFrom hdf5r H5File h5attr h5attr_names
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat CreateSeuratObject SetAssayData CreateDimReducObject
#' @importFrom SeuratObject Cells
#'
#' @export
#'
LoadH5AD <- function(file, assay.name = "RNA", use.bpcells = NULL, verbose = TRUE) {
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  h5ad <- H5File$new(file, mode = "r")
  on.exit(h5ad$close_all())

  if (verbose) {
    message("Loading H5AD file: ", file)
  }

  # Helper function to read H5AD sparse or dense matrix
  ReadH5ADMatrix <- function(h5_obj, transpose = TRUE) {
    if (inherits(h5_obj, "H5Group")) {
      # Sparse matrix (CSR or CSC format in h5ad)
      if (h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
        data_vals <- h5_obj[["data"]][]
        indices <- h5_obj[["indices"]][]   # 0-based
        indptr <- h5_obj[["indptr"]][]

        # Detect encoding type: CSR vs CSC
        encoding <- tryCatch(h5attr(h5_obj, "encoding-type"), error = function(e) "csr_matrix")
        is_csc <- identical(encoding, "csc_matrix")

        # Get dimensions from shape attribute
        if (h5_obj$attr_exists("shape")) {
          shape <- h5attr(h5_obj, "shape")
          n_rows <- shape[1]
          n_cols <- shape[2]
        } else if (is_csc) {
          n_cols <- length(indptr) - 1L
          n_rows <- if (length(indices) > 0) max(indices) + 1L else 0L
        } else {
          n_rows <- length(indptr) - 1L
          n_cols <- if (length(indices) > 0) max(indices) + 1L else 0L
        }

        if (is_csc) {
          # CSC format: indptr = column pointers, indices = row indices (0-based)
          # sparseMatrix with i + p form constructs dgCMatrix directly
          mat <- sparseMatrix(
            i = as.integer(indices) + 1L,
            p = as.integer(indptr),
            x = as.numeric(data_vals),
            dims = c(n_rows, n_cols)
          )
          if (transpose) mat <- t(mat)
        } else {
          # CSR format: indptr = row pointers, indices = column indices (0-based)
          indices_1based <- indices + 1L
          row_indices <- rep(seq_len(n_rows), diff(indptr))
          mat <- sparseMatrix(
            i = row_indices,
            j = indices_1based,
            x = data_vals,
            dims = c(n_rows, n_cols),
            index1 = TRUE
          )
          if (transpose) mat <- t(mat)
        }
        return(mat)
      }
    } else if (inherits(h5_obj, "H5D")) {
      # Dense matrix
      mat <- h5_obj[,]
      if (transpose) mat <- t(mat)
      return(mat)
    }
    stop("Unknown matrix format", call. = FALSE)
  }

  # 1. Read cell names
  if (verbose) message("Reading cell names...")
  cell.names <- NULL
  obs_compound_df <- NULL
  if (h5ad$exists("obs")) {
    obs_obj <- h5ad[["obs"]]
    if (inherits(obs_obj, "H5Group")) {
      # Modern h5ad format: obs is a group
      if (obs_obj$exists("_index")) {
        cell.names <- as.character(obs_obj[["_index"]][])
      } else if (obs_obj$exists("index")) {
        cell.names <- as.character(obs_obj[["index"]][])
      } else if (obs_obj$attr_exists("_index")) {
        # AnnData convention: _index attribute names the index column
        idx_col <- h5attr(obs_obj, "_index")
        if (obs_obj$exists(idx_col)) {
          cell.names <- as.character(obs_obj[[idx_col]][])
        }
      }
    } else if (inherits(obs_obj, "H5D")) {
      # Legacy h5ad format: obs is a compound HDF5 dataset
      obs_compound_df <- obs_obj$read()
      if ("index" %in% names(obs_compound_df)) {
        cell.names <- as.character(obs_compound_df$index)
      } else if ("_index" %in% names(obs_compound_df)) {
        cell.names <- as.character(obs_compound_df[["_index"]])
      } else if (!is.null(rownames(obs_compound_df))) {
        cell.names <- rownames(obs_compound_df)
      }
    }
  }

  # Generate cell names if not found
  if (is.null(cell.names)) {
    n_cells <- if (h5ad$exists("X")) {
      if (inherits(h5ad[["X"]], "H5Group") && h5ad[["X"]]$attr_exists("shape")) {
        h5attr(h5ad[["X"]], "shape")[1]
      } else if (inherits(h5ad[["X"]], "H5D")) {
        h5ad[["X"]]$dims[1]
      }
    } else {
      stop("Cannot determine number of cells", call. = FALSE)
    }
    cell.names <- paste0("Cell", seq_len(n_cells))
  }

  # 2. Read feature names
  if (verbose) message("Reading feature names...")
  feature.names <- NULL
  var_compound_df <- NULL
  if (h5ad$exists("var")) {
    var_obj <- h5ad[["var"]]
    if (inherits(var_obj, "H5Group")) {
      if (var_obj$exists("_index")) {
        feature.names <- as.character(var_obj[["_index"]][])
      } else if (var_obj$exists("index")) {
        feature.names <- as.character(var_obj[["index"]][])
      } else if (var_obj$attr_exists("_index")) {
        idx_col <- h5attr(var_obj, "_index")
        if (var_obj$exists(idx_col)) {
          feature.names <- as.character(var_obj[[idx_col]][])
        }
      }
    } else if (inherits(var_obj, "H5D")) {
      # Legacy compound dataset
      var_compound_df <- var_obj$read()
      if ("index" %in% names(var_compound_df)) {
        feature.names <- as.character(var_compound_df$index)
      } else if ("_index" %in% names(var_compound_df)) {
        feature.names <- as.character(var_compound_df[["_index"]])
      }
    }
  }

  # Generate feature names if not found
  if (is.null(feature.names)) {
    n_features <- if (h5ad$exists("X")) {
      if (inherits(h5ad[["X"]], "H5Group") && h5ad[["X"]]$attr_exists("shape")) {
        h5attr(h5ad[["X"]], "shape")[2]
      } else if (inherits(h5ad[["X"]], "H5D")) {
        h5ad[["X"]]$dims[2]
      }
    } else {
      stop("Cannot determine number of features", call. = FALSE)
    }
    feature.names <- paste0("Gene", seq_len(n_features))
  }

  # Deduplicate feature names (some datasets have duplicates)
  if (anyDuplicated(feature.names)) {
    feature.names <- make.unique(feature.names)
  }

  # 3. Read main expression matrix
  if (verbose) message("Reading expression matrix...")
  if (!h5ad$exists("X")) {
    stop("No expression matrix (X) found in h5ad file", call. = FALSE)
  }

  use_bpcells <- !is.null(use.bpcells)
  if (use_bpcells) {
    if (!requireNamespace("BPCells", quietly = TRUE)) {
      stop("BPCells package is required for on-disk loading. ",
           "Install with: remotes::install_github('bnprks/BPCells')", call. = FALSE)
    }

    # Prefer raw/X (raw counts) over X (often normalized) — matches non-BPCells path
    has_raw <- h5ad$exists("raw") && h5ad[["raw"]]$exists("X")
    bp_group <- if (has_raw) "raw/X" else "X"
    if (verbose) {
      if (has_raw) {
        message("Loading raw counts (raw/X) via BPCells (on-disk)...")
      } else {
        message("Loading expression matrix (X) via BPCells (on-disk)...")
      }
    }

    # Close hdf5r handle before BPCells opens the file (avoid lock conflicts)
    h5ad$close_all()

    if (isTRUE(use.bpcells)) {
      # use.bpcells = TRUE: read directly from h5ad (backed by HDF5)
      expr_matrix <- BPCells::open_matrix_anndata_hdf5(path = file, group = bp_group)
    } else {
      # use.bpcells = directory path: write to BPCells dir for faster repeated access
      bpcells_dir <- file.path(use.bpcells, "counts")
      dir.create(bpcells_dir, recursive = TRUE, showWarnings = FALSE)
      bpcells_mat <- BPCells::open_matrix_anndata_hdf5(path = file, group = bp_group)
      BPCells::write_matrix_dir(mat = bpcells_mat, dir = bpcells_dir, overwrite = TRUE)
      expr_matrix <- BPCells::open_matrix_dir(dir = bpcells_dir)
    }
    # Reopen h5ad for metadata reading
    h5ad <- H5File$new(file, mode = "r")
  } else {
    expr_matrix <- ReadH5ADMatrix(h5ad[["X"]], transpose = TRUE)
  }

  # Handle dimension mismatches: check if matrix needs transposing
  # This can happen with dense matrices where R/Python storage conventions differ
  n_features <- length(feature.names)
  n_cells <- length(cell.names)
  if (nrow(expr_matrix) == n_cells && ncol(expr_matrix) == n_features &&
      nrow(expr_matrix) != n_features) {
    expr_matrix <- t(expr_matrix)
  }

  if (nrow(expr_matrix) != n_features) {
    warning(sprintf(
      "Feature count mismatch: matrix has %d rows but %d feature names. Truncating.",
      nrow(expr_matrix), n_features
    ), immediate. = TRUE)
    feature.names <- feature.names[seq_len(nrow(expr_matrix))]
  }
  if (ncol(expr_matrix) != n_cells) {
    warning(sprintf(
      "Cell count mismatch: matrix has %d columns but %d cell names. Truncating.",
      ncol(expr_matrix), n_cells
    ), immediate. = TRUE)
    cell.names <- cell.names[seq_len(ncol(expr_matrix))]
  }

  rownames(expr_matrix) <- feature.names
  colnames(expr_matrix) <- cell.names

  # 4. Create Seurat object
  if (verbose) message("Creating Seurat object...")
  seurat_obj <- CreateSeuratObject(
    counts = expr_matrix,
    project = "H5AD",
    assay = assay.name,
    min.cells = 0,
    min.features = 0
  )

  # 5. Add raw counts if present
  if (h5ad$exists("raw") && h5ad[["raw"]]$exists("X")) {
    if (use_bpcells) {
      # In BPCells mode, skip raw/X to preserve on-disk matrix.
      # Loading raw/X in-memory would defeat the purpose of on-disk mode.
      if (verbose) message("Skipping raw counts (BPCells on-disk mode preserves X as counts)")
    } else {
      if (verbose) message("Adding raw counts...")

      raw_features <- NULL
      if (h5ad[["raw"]]$exists("var")) {
        raw_var <- h5ad[["raw/var"]]
        if (raw_var$exists("_index")) {
          raw_features <- as.character(raw_var[["_index"]][])
        } else if (raw_var$exists("index")) {
          raw_features <- as.character(raw_var[["index"]][])
        }
      }

      if (!is.null(raw_features)) {
        raw_matrix <- ReadH5ADMatrix(h5ad[["raw/X"]], transpose = TRUE)

        # Handle dense matrix transpose
        n_raw_features <- length(raw_features)
        if (nrow(raw_matrix) == n_cells && ncol(raw_matrix) == n_raw_features &&
            nrow(raw_matrix) != n_raw_features) {
          raw_matrix <- t(raw_matrix)
        }

        # Match dimensions
        raw_features <- raw_features[seq_len(min(length(raw_features), nrow(raw_matrix)))]
        rownames(raw_matrix) <- raw_features
        colnames(raw_matrix) <- cell.names

        # Find common features
        common_features <- intersect(feature.names, raw_features)
        if (length(common_features) > 0) {
          # raw/X -> counts layer (actual raw counts)
          raw_subset <- raw_matrix[common_features, , drop = FALSE]
          seurat_obj[[assay.name]] <- SetAssayData(
            object = seurat_obj[[assay.name]],
            layer = "counts",
            new.data = raw_subset
          )
          # X -> data layer (normalized) when raw/X exists
          x_subset <- expr_matrix[common_features, , drop = FALSE]
          seurat_obj[[assay.name]] <- SetAssayData(
            object = seurat_obj[[assay.name]],
            layer = "data",
            new.data = x_subset
          )
          if (verbose) message("  Set raw/X as counts, X as data (normalized)")
        }
      }
    }
  }

  # 6. Add layers if present
  if (h5ad$exists("layers")) {
    if (use_bpcells) {
      if (verbose) message("Skipping additional layers (BPCells on-disk mode)")
    } else {
      if (verbose) message("Adding layers...")
      layer_names <- names(h5ad[["layers"]])

      for (layer_name in layer_names) {
        if (verbose) message("  Adding layer: ", layer_name)

        layer_matrix <- ReadH5ADMatrix(h5ad[["layers"]][[layer_name]], transpose = TRUE)

        # Ensure dimensions match
        if (nrow(layer_matrix) == nrow(expr_matrix) && ncol(layer_matrix) == ncol(expr_matrix)) {
          rownames(layer_matrix) <- feature.names
          colnames(layer_matrix) <- cell.names

          # Map layer names to Seurat slots
          seurat_slot <- switch(layer_name,
            "counts" = "counts",
            "data" = "data",
            "log_normalized" = "data",
            "scale.data" = "scale.data",
            "scaled" = "scale.data",
            layer_name
          )

          tryCatch({
            seurat_obj[[assay.name]] <- SetAssayData(
              object = seurat_obj[[assay.name]],
              layer = seurat_slot,
              new.data = layer_matrix
            )
          }, error = function(e) {
            if (verbose) warning("Could not add layer ", layer_name, ": ", e$message, immediate. = TRUE)
          })
        }
      }
    }
  }

  # 7. Add cell metadata
  if (h5ad$exists("obs")) {
    if (verbose) message("Adding cell metadata...")

    if (!is.null(obs_compound_df)) {
      # Legacy compound dataset: obs was already read as a data.frame
      meta_cols <- setdiff(names(obs_compound_df), c("_index", "index"))
      for (col in meta_cols) {
        tryCatch({
          seurat_obj[[col]] <- obs_compound_df[[col]]
        }, error = function(e) {
          if (verbose) message("Could not add metadata column '", col, "': ", e$message)
        })
      }
    } else {
      obs_group <- h5ad[["obs"]]
      # Exclude index columns
      obs_exclude <- c("_index", "index", "__categories")
      if (obs_group$attr_exists("_index")) {
        obs_exclude <- c(obs_exclude, h5attr(obs_group, "_index"))
      }
      obs_cols <- setdiff(names(obs_group), obs_exclude)

      # Cache __categories group reference if it exists (legacy format)
      has_legacy_cats <- obs_group$exists("__categories")
      legacy_cats <- if (has_legacy_cats) obs_group[["__categories"]] else NULL
      legacy_cat_names <- if (has_legacy_cats) names(legacy_cats) else character(0)

      for (col in obs_cols) {
        if (verbose) message("  Adding metadata: ", col)

        tryCatch({
          meta_values <- NULL
          col_obj <- obs_group[[col]]

          if (inherits(col_obj, "H5Group")) {
            # Modern h5ad categorical format: group with categories/codes sub-datasets
            encoding_type <- tryCatch(h5attr(col_obj, "encoding-type"), error = function(e) "")
            if (encoding_type == "categorical" && col_obj$exists("categories") && col_obj$exists("codes")) {
              codes <- col_obj[["codes"]]$read()
              categories <- as.character(col_obj[["categories"]]$read())
              codes[codes == -1L] <- NA_integer_
              meta_values <- factor(categories[codes + 1L], levels = categories)
            }
          } else if (inherits(col_obj, "H5D")) {
            # Legacy categorical format: __categories group with per-column category lists
            if (has_legacy_cats && col %in% legacy_cat_names) {
              codes <- col_obj$read()
              categories <- as.character(legacy_cats[[col]]$read())
              codes[codes == -1L] <- NA_integer_
              meta_values <- factor(categories[codes + 1L], levels = categories)
            } else {
              # Numeric or string dataset
              meta_values <- col_obj$read()
              if (is.character(meta_values)) {
                meta_values <- as.character(meta_values)
              }
            }
          }

          # Add to Seurat object if we got values
          if (!is.null(meta_values)) {
            seurat_obj[[col]] <- meta_values
          }
        }, error = function(e) {
          if (verbose) warning("Could not add metadata column '", col, "': ", e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 8. Add dimensional reductions
  if (h5ad$exists("obsm")) {
    if (verbose) message("Adding dimensional reductions...")

    for (reduc_name in names(h5ad[["obsm"]])) {
      clean_name <- gsub("^X_", "", reduc_name)
      if (verbose) message("  Adding reduction: ", clean_name)

      embeddings <- h5ad[["obsm"]][[reduc_name]][,]

      # HDF5 stores in row-major (C) order, R reads in column-major (Fortran) order
      # h5ad obsm shape is (n_obs, n_dims) but hdf5r may return (n_dims, n_obs)
      if (nrow(embeddings) != length(cell.names) && ncol(embeddings) == length(cell.names)) {
        embeddings <- t(embeddings)
      }

      # Check dimensions
      if (nrow(embeddings) != length(cell.names)) {
        warning("Skipping reduction ", clean_name, " - dimension mismatch", immediate. = TRUE)
        next
      }

      rownames(embeddings) <- cell.names
      colnames(embeddings) <- paste0(toupper(clean_name), "_", seq_len(ncol(embeddings)))

      # Determine key
      key <- switch(clean_name,
        "pca" = "PC_",
        "tsne" = "tSNE_",
        "umap" = "UMAP_",
        paste0(toupper(clean_name), "_")
      )

      reduc_obj <- CreateDimReducObject(
        embeddings = embeddings,
        key = key,
        assay = assay.name
      )

      seurat_obj[[clean_name]] <- reduc_obj
    }
  }

  # 9. Add feature metadata
  if (h5ad$exists("var")) {
    if (verbose) message("Adding feature metadata...")

    if (!is.null(var_compound_df)) {
      # Legacy compound dataset: var was already read as a data.frame
      meta_cols <- setdiff(names(var_compound_df), c("_index", "index"))
      for (col in meta_cols) {
        tryCatch({
          meta_values <- var_compound_df[[col]]
          if (col == "highly_variable") {
            if (is.character(meta_values)) meta_values <- meta_values == "True"
            else if (is.numeric(meta_values)) meta_values <- as.logical(meta_values)
          }
          if (length(meta_values) == nrow(seurat_obj)) {
            names(meta_values) <- feature.names
            seurat_obj[[assay.name]][[col]] <- meta_values
            if (col == "highly_variable" && is.logical(meta_values)) {
              VariableFeatures(seurat_obj) <- feature.names[meta_values]
            }
          }
        }, error = function(e) {
          if (verbose) warning("Could not add feature metadata ", col, ": ", e$message, immediate. = TRUE)
        })
      }
    } else {
      var_group <- h5ad[["var"]]
      var_exclude <- c("_index", "index", "__categories")
      if (var_group$attr_exists("_index")) {
        var_exclude <- c(var_exclude, h5attr(var_group, "_index"))
      }
      var_cols <- setdiff(names(var_group), var_exclude)

      # Cache __categories if present (legacy format)
      has_var_cats <- var_group$exists("__categories")
      var_cats <- if (has_var_cats) var_group[["__categories"]] else NULL
      var_cat_names <- if (has_var_cats) names(var_cats) else character(0)

      for (col in var_cols) {
        if (verbose) message("  Adding feature metadata: ", col)

        tryCatch({
          meta_values <- NULL
          col_obj <- var_group[[col]]

          if (inherits(col_obj, "H5Group")) {
            # Modern h5ad categorical format
            encoding_type <- tryCatch(h5attr(col_obj, "encoding-type"), error = function(e) "")
            if (encoding_type == "categorical" && col_obj$exists("categories") && col_obj$exists("codes")) {
              codes <- col_obj[["codes"]]$read()
              categories <- as.character(col_obj[["categories"]]$read())
              codes[codes == -1L] <- NA_integer_
              meta_values <- factor(categories[codes + 1L], levels = categories)
            }
          } else if (inherits(col_obj, "H5D")) {
            # Check legacy categorical format
            if (has_var_cats && col %in% var_cat_names) {
              codes <- col_obj$read()
              categories <- as.character(var_cats[[col]]$read())
              codes[codes == -1L] <- NA_integer_
              meta_values <- factor(categories[codes + 1L], levels = categories)
            } else {
              meta_values <- col_obj$read()
            }
          }

          if (is.null(meta_values)) next

          # Special handling for highly_variable: convert to logical
          if (col == "highly_variable") {
            if (is.factor(meta_values)) {
              meta_values <- as.character(meta_values) == "True"
            } else if (is.numeric(meta_values)) {
              meta_values <- as.logical(meta_values)
            }
          }

          # Ensure length matches and name with feature names for Seurat v5 compatibility
          if (length(meta_values) == nrow(seurat_obj)) {
            names(meta_values) <- feature.names
            seurat_obj[[assay.name]][[col]] <- meta_values

            # Set variable features if highly_variable column exists
            if (col == "highly_variable" && is.logical(meta_values)) {
              VariableFeatures(seurat_obj) <- feature.names[meta_values]
            }
          }
        }, error = function(e) {
          if (verbose) warning("Could not add feature metadata ", col, ": ", e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 10. Add neighbor graphs from obsp
  if (h5ad$exists("obsp")) {
    if (verbose) message("Adding neighbor graphs...")

    for (graph_name in names(h5ad[["obsp"]])) {
      if (verbose) message("  Adding graph: ", graph_name)

      graph_matrix <- ReadH5ADMatrix(h5ad[["obsp"]][[graph_name]], transpose = FALSE)

      # Ensure square matrix with correct dimensions
      if (nrow(graph_matrix) == ncol(graph_matrix) && nrow(graph_matrix) == ncol(seurat_obj)) {
        rownames(graph_matrix) <- cell.names
        colnames(graph_matrix) <- cell.names

        # Map to Seurat graph names
        seurat_graph_name <- switch(graph_name,
          "connectivities" = paste0(assay.name, "_snn"),
          "distances" = paste0(assay.name, "_nn"),
          graph_name
        )

        seurat_obj@graphs[[seurat_graph_name]] <- as.Graph(graph_matrix)
      }
    }
  }

  # 11. Add uns (unstructured) data to misc
  if (h5ad$exists("uns")) {
    if (verbose) message("Adding unstructured data...")
    uns_group <- h5ad[["uns"]]

    for (item in names(uns_group)) {
      tryCatch({
        if (inherits(uns_group[[item]], "H5D")) {
          # Simple dataset
          seurat_obj@misc[[item]] <- uns_group[[item]][]
        } else if (inherits(uns_group[[item]], "H5Group")) {
          # Complex group - store as list
          if (verbose) message("  Storing complex uns item: ", item)
          # For now, just note it exists
          seurat_obj@misc[[paste0(item, "_present")]] <- TRUE
        }
      }, error = function(e) {
        if (verbose) warning("Could not add uns item ", item, ": ", e$message, immediate. = TRUE)
      })
    }
  }

  # 12. Add spatial data support
  if (h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
    seurat_obj <- ConvertH5ADSpatialToSeurat(
      h5ad_file = h5ad,
      seurat_obj = seurat_obj,
      assay_name = assay.name,
      verbose = verbose
    )
  }

  if (verbose) {
    message("\nSuccessfully loaded H5AD file")
    message("  Cells: ", ncol(seurat_obj))
    message("  Features: ", nrow(seurat_obj))
    message("  Assays: ", paste(names(seurat_obj@assays), collapse = ", "))
    if (length(seurat_obj@reductions) > 0) {
      message("  Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))
    }
    if (length(seurat_obj@graphs) > 0) {
      message("  Graphs: ", paste(names(seurat_obj@graphs), collapse = ", "))
    }
    message("  Metadata columns: ", ncol(seurat_obj@meta.data))
  }

  return(seurat_obj)
}
