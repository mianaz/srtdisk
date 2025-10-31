#' @keywords internal
ReadSparseMatrixDims <- function(h5_group, indices = NULL, indptr = NULL,
                                  attr_names = c("shape", "dims", "Dim"), verbose = FALSE) {
  for (attr_name in attr_names) {
    if (h5_group$attr_exists(attr_name)) {
      dims <- as.integer(h5attr(h5_group, attr_name))
      if (length(dims) >= 2) {
        if (verbose) message("Read dimensions from '", attr_name, "': ", dims[1], " x ", dims[2])
        return(list(nrows = dims[1], ncols = dims[2]))
      }
    }
  }

  if (!is.null(indices) && !is.null(indptr) && length(indices) > 0) {
    nrows <- max(indices) + 1L
    ncols <- length(indptr) - 1L
    if (verbose) message("Inferred dimensions: ", nrows, " x ", ncols)
    return(list(nrows = nrows, ncols = ncols))
  }

  if (verbose) message("Cannot determine sparse matrix dimensions")
  list(nrows = NA, ncols = NA)
}

#' @keywords internal
ValidateSparseMatrixDims <- function(nrows, ncols, verbose = FALSE) {
  if (is.na(nrows) || is.na(ncols) || nrows <= 0 || ncols <= 0) {
    if (verbose) message("Invalid dimensions: ", nrows, " x ", ncols)
    return(FALSE)
  }
  TRUE
}

#' @keywords internal
CreateSparseMatrixSafe <- function(indices, indptr, data, nrows, ncols, verbose = FALSE) {
  if (!ValidateSparseMatrixDims(nrows, ncols, verbose)) return(NULL)

  tryCatch(
    Matrix::sparseMatrix(i = indices + 1L, p = indptr, x = data,
                        dims = c(nrows, ncols), index1 = TRUE),
    error = function(e) {
      if (verbose) {
        message("Error creating sparse matrix: ", conditionMessage(e))
        message("  dims: ", nrows, "x", ncols, ", indptr: ", length(indptr),
                ", indices: ", length(indices), ", data: ", length(data))
      }
      NULL
    }
  )
}

