#' @include zzz.R
#' @importFrom hdf5r H5Group H5File h5attr
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#'
NULL

#' V5 Layer Support Functions
#'
#' Internal functions for reading and writing V5 layer structures in h5Seurat files
#'
#' @name V5LayerSupport
#' @rdname V5LayerSupport
#'
#' @keywords internal
#'
NULL

#' Write V5 Layer Data
#'
#' Write matrix data in V5 layer format (sparse or dense)
#'
#' @param h5_group The HDF5 group to write to
#' @param layer_name Name of the layer (e.g., "counts", "data", "scale.data")
#' @param matrix_data The matrix data to write
#' @param features Character vector of feature names
#' @param verbose Show progress messages
#'
#' @return Invisible NULL
#'
#' @keywords internal
#'
WriteV5Layer <- function(h5_group, layer_name, matrix_data, features = NULL, verbose = FALSE) {
  if (verbose) message(paste0("Writing V5 layer '", layer_name, "'"))

  if (!h5_group$exists("layers")) {
    h5_group$create_group("layers")
  }
  layers_group <- h5_group[["layers"]]

  if (layers_group$exists(layer_name)) {
    if (verbose) message(paste0("Removing existing layer '", layer_name, "'"))
    layers_group$link_delete(layer_name)
  }

  is_sparse <- inherits(matrix_data, "sparseMatrix") || inherits(matrix_data, "dgCMatrix")

  if (is_sparse) {
    if (!inherits(matrix_data, "dgCMatrix")) {
      matrix_data <- as(matrix_data, "dgCMatrix")
    }

    layer_group <- layers_group$create_group(layer_name)

    # Write sparse matrix in CSC format (V5 standard)
    layer_group[["data"]] <- matrix_data@x
    layer_group[["indices"]] <- as.integer(matrix_data@i)
    layer_group[["indptr"]] <- as.integer(matrix_data@p)

    layer_group$create_attr("shape", c(nrow(matrix_data), ncol(matrix_data)))
    layer_group$create_attr("encoding-type", "csc_matrix")
    layer_group$create_attr("encoding-version", "0.1.0")

    if (verbose) message(paste0("  Wrote sparse matrix: ", length(matrix_data@x), " non-zero values, ",
            "shape ", nrow(matrix_data), " x ", ncol(matrix_data)))
  } else {
    layers_group[[layer_name]] <- as.matrix(matrix_data)

    if (verbose) message(paste0("  Wrote dense matrix: shape ", nrow(matrix_data), " x ", ncol(matrix_data)))
  }

  if (!is.null(features) && !h5_group$exists("meta.data")) {
    h5_group$create_group("meta.data")
  }
  if (!is.null(features) && !h5_group[["meta.data"]]$exists("_index")) {
    h5_group[["meta.data"]][["_index"]] <- features
    if (verbose) message(paste0("  Stored ", length(features), " feature names"))
  }

  return(invisible(NULL))
}

#' Read V5 Layer Data
#'
#' Read matrix data from V5 layer format
#'
#' @param h5_group The HDF5 assay group
#' @param layer_name Name of the layer to read
#' @param features Expected feature names (for dimension validation)
#' @param cells Expected cell names (for dimension validation)
#' @param verbose Show progress messages
#'
#' @return A matrix object, or NULL if layer doesn't exist
#'
#' @keywords internal
#'
ReadV5Layer <- function(h5_group, layer_name, features = NULL, cells = NULL, verbose = FALSE) {
  if (!h5_group$exists("layers")) {
    if (verbose) message("No layers group found")
    return(NULL)
  }

  layers_group <- h5_group[["layers"]]

  if (!layers_group$exists(layer_name)) {
    if (verbose) message(paste0("Layer '", layer_name, "' not found"))
    return(NULL)
  }

  layer_obj <- layers_group[[layer_name]]

  if (inherits(layer_obj, "H5Group")) {
    if (verbose) message(paste0("Reading sparse V5 layer '", layer_name, "'"))

    if (exists("ReadSparseMatrix", envir = asNamespace("SeuratDisk"))) {
      mat <- ReadSparseMatrix(layer_obj, verbose = verbose)
    } else {
      if (layer_obj$exists("data") && layer_obj$exists("indices") && layer_obj$exists("indptr")) {
        data_vals <- layer_obj[["data"]][]
        indices <- layer_obj[["indices"]][] + 1L
        indptr <- layer_obj[["indptr"]][]

        if (layer_obj$attr_exists("shape")) {
          shape <- h5attr(x = layer_obj, which = "shape")
          nrows <- shape[1]
          ncols <- shape[2]
        } else {
          ncols <- length(indptr) - 1L
          nrows <- max(indices)
        }

        sparse_mat <- Matrix::sparseMatrix(
          i = indices,
          p = indptr,
          x = data_vals,
          dims = c(nrows, ncols),
          index1 = TRUE
        )

        mat <- as.matrix(sparse_mat)
      } else {
        warning("Unknown sparse matrix format in layer '", layer_name, "'")
        return(NULL)
      }
    }
  } else {
    if (verbose) message(paste0("Reading dense V5 layer '", layer_name, "'"))
    mat <- as.matrix(layer_obj[,])
  }

  if (!is.null(features) && !is.null(cells)) {
    expected_nrows <- length(features)
    expected_ncols <- length(cells)

    if (nrow(mat) != expected_nrows || ncol(mat) != expected_ncols) {
      if (nrow(mat) == expected_ncols && ncol(mat) == expected_nrows) {
        if (verbose) message("Transposing matrix to match expected dimensions")
        mat <- t(mat)
      } else if (nrow(mat) < expected_nrows && ncol(mat) == expected_ncols) {
        if (verbose) message("Expanding matrix to full feature space")
        full_mat <- matrix(0, nrow = expected_nrows, ncol = expected_ncols)
        full_mat[1:nrow(mat), ] <- mat
        mat <- full_mat
      }
    }

    rownames(mat) <- features
    colnames(mat) <- cells
  }

  return(mat)
}

#' Get Layer Path
#'
#' Get the correct path for a layer/slot in V4 vs V5 format
#'
#' @param h5_group The HDF5 assay group
#' @param slot_name Name of the slot/layer
#' @param verbose Show progress messages
#'
#' @return The path to the data, or NULL if not found
#'
#' @keywords internal
#'
GetLayerPath <- function(h5_group, slot_name, verbose = FALSE) {
  # First check V5 structure
  v5_path <- paste0("layers/", slot_name)
  if (h5_group$exists(v5_path)) {
    if (verbose) {
      message("Found V5 path: ", v5_path)
    }
    return(v5_path)
  }

  # Check V4 structure (direct slot)
  if (h5_group$exists(slot_name)) {
    if (verbose) {
      message("Found V4 path: ", slot_name)
    }
    return(slot_name)
  }

  if (verbose) {
    message("No path found for slot '", slot_name, "'")
  }
  return(NULL)
}

#' Migrate V4 to V5 Structure
#'
#' Convert V4 h5Seurat structure to V5 format
#'
#' @param h5_file Path to the h5Seurat file or h5Seurat object
#' @param output_file Path for the converted file (optional, overwrites if NULL)
#' @param verbose Show progress messages
#'
#' @return Path to the converted file
#'
#' @keywords internal
#'
MigrateV4ToV5 <- function(h5_file, output_file = NULL, verbose = FALSE) {
  # Open the file
  if (is.character(h5_file)) {
    h5f <- h5Seurat$new(h5_file, mode = "r+")
    on.exit(h5f$close_all())
  } else {
    h5f <- h5_file
  }

  # Check version
  current_version <- h5f$version()
  if (verbose) {
    message("Current version: ", current_version)
  }

  # If output file specified, copy first
  if (!is.null(output_file)) {
    file.copy(h5_file, output_file, overwrite = TRUE)
    h5f$close_all()
    h5f <- h5Seurat$new(output_file, mode = "r+")
  }

  # Process each assay
  assays <- names(h5f[["assays"]])
  for (assay in assays) {
    if (verbose) {
      message("Migrating assay: ", assay)
    }

    assay_group <- h5f[["assays"]][[assay]]

    # Check if already has layers structure
    if (assay_group$exists("layers")) {
      if (verbose) {
        message("  Already has layers structure, skipping")
      }
      next
    }

    # Create layers group
    layers_group <- assay_group$create_group("layers")

    # Move direct slots to layers
    slots_to_move <- c("counts", "data", "scale.data")
    for (slot in slots_to_move) {
      if (assay_group$exists(slot)) {
        if (verbose) {
          message("  Moving ", slot, " to layers/", slot)
        }

        # Read the data
        slot_obj <- assay_group[[slot]]

        # Check if it's sparse
        if (inherits(slot_obj, "H5Group")) {
          # Already a sparse matrix group, just move it
          # Note: HDF5 doesn't support direct move, so we need to copy
          # This is a simplified version - in production, would need proper copying
          warning("Cannot directly move sparse matrix groups - manual intervention needed")
        } else {
          # Read and rewrite in V5 format
          mat <- as.matrix(slot_obj[,])

          # Convert to sparse if beneficial
          if (sum(mat == 0) / length(mat) > .SPARSITY_THRESHOLD) {
            mat <- as(mat, "dgCMatrix")
          }

          # Write in V5 format using WriteV5Layer
          WriteV5Layer(
            h5_group = assay_group,
            layer_name = slot,
            matrix_data = mat,
            verbose = verbose
          )

          # Delete old slot
          assay_group$link_delete(slot)
        }
      }
    }
  }

  # Update version
  h5f$set.version("5.2.1")

  if (verbose) {
    message("Migration complete")
  }

  return(ifelse(is.null(output_file), h5_file, output_file))
}

#' Check if h5Seurat file uses V5 structure
#'
#' @param h5_file Path to h5Seurat file or h5Seurat object
#'
#' @return Logical indicating if file uses V5 structure
#'
#' @keywords internal
#'
IsV5Structure <- function(h5_file) {
  if (is.character(h5_file)) {
    h5f <- h5Seurat$new(h5_file, mode = "r")
    on.exit(h5f$close_all())
  } else {
    h5f <- h5_file
  }

  # Check for V5 indicators
  if (h5f$exists("assays")) {
    assays <- names(h5f[["assays"]])
    if (length(assays) > 0) {
      first_assay <- h5f[["assays"]][[assays[1]]]
      return(first_assay$exists("layers"))
    }
  }

  return(FALSE)
}