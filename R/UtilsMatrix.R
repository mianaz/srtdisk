#' Matrix Utility Functions
#'
#' This file contains utility functions for matrix operations including
#' setting names, validation, and transformation.
#'
#' @name UtilsMatrix
#' @keywords internal
NULL

# Removed: SetMatrixNames() - unused wrapper function
# Use direct rownames()/colnames() assignment or dimnames<- instead

#' Validate matrix dimensions match expected values
#'
#' Checks if a matrix has expected dimensions and optionally attempts to
#' transpose or expand if dimensions don't match.
#'
#' @param mat Matrix to validate
#' @param expected_nrows Expected number of rows
#' @param expected_ncols Expected number of columns
#' @param allow_transpose If TRUE, transpose matrix if dimensions are swapped
#' @param verbose Print diagnostic messages
#'
#' @return List with elements:
#'   \item{valid}{Logical, TRUE if dimensions match (after any corrections)}
#'   \item{matrix}{Corrected matrix (may be transposed)}
#'   \item{transposed}{Logical, TRUE if matrix was transposed}
#'
#' @keywords internal
#'
ValidateMatrixDimensions <- function(
  mat,
  expected_nrows,
  expected_ncols,
  allow_transpose = TRUE,
  verbose = FALSE
) {
  actual_nrows <- nrow(mat)
  actual_ncols <- ncol(mat)

  if (actual_nrows == expected_nrows && actual_ncols == expected_ncols) {
    if (verbose) {
      message("Matrix dimensions match: ", actual_nrows, " x ", actual_ncols)
    }
    return(list(valid = TRUE, matrix = mat, transposed = FALSE))
  }

  if (allow_transpose &&
      actual_nrows == expected_ncols &&
      actual_ncols == expected_nrows) {
    if (verbose) {
      message(
        "Matrix appears transposed (",
        actual_nrows, " x ", actual_ncols,
        "), transposing to match expected (",
        expected_nrows, " x ", expected_ncols, ")"
      )
    }
    return(list(valid = TRUE, matrix = t(mat), transposed = TRUE))
  }

  if (verbose) {
    message(
      "Matrix dimensions mismatch: got ",
      actual_nrows, " x ", actual_ncols,
      ", expected ",
      expected_nrows, " x ", expected_ncols
    )
  }
  return(list(valid = FALSE, matrix = mat, transposed = FALSE))
}

#' Expand matrix to full feature space
#'
#' Expands a matrix that contains a subset of features to the full feature
#' space by adding zero rows for missing features.
#'
#' @param mat Matrix to expand (features x cells)
#' @param current_features Character vector of features in mat
#' @param full_features Character vector of all features in full space
#' @param verbose Print diagnostic messages
#'
#' @return Expanded matrix with all features
#' @keywords internal
#'
ExpandMatrixToFullFeatures <- function(
  mat,
  current_features,
  full_features,
  verbose = FALSE
) {
  if (length(current_features) == length(full_features) &&
      all(current_features == full_features)) {
    if (verbose) {
      message("Matrix already contains all features")
    }
    return(mat)
  }

  if (verbose) {
    message(
      "Expanding matrix from ",
      length(current_features),
      " to ",
      length(full_features),
      " features"
    )
  }

  if (inherits(mat, "dgCMatrix") || inherits(mat, "Matrix")) {
    full_mat <- Matrix::Matrix(
      0,
      nrow = length(full_features),
      ncol = ncol(mat),
      sparse = TRUE
    )
  } else {
    full_mat <- matrix(
      0,
      nrow = length(full_features),
      ncol = ncol(mat)
    )
  }

  rownames(full_mat) <- full_features
  colnames(full_mat) <- colnames(mat)

  feature_idx <- match(current_features, full_features)
  if (any(is.na(feature_idx))) {
    warning(
      "Some features in matrix not found in full feature space: ",
      paste(current_features[is.na(feature_idx)], collapse = ", ")
    )
    feature_idx <- feature_idx[!is.na(feature_idx)]
    current_features <- current_features[!is.na(match(current_features, full_features))]
  }

  full_mat[feature_idx, ] <- mat[current_features, , drop = FALSE]

  return(full_mat)
}

#' Check if matrix is empty or all zeros
#'
#' @param mat Matrix to check
#' @return Logical indicating if matrix is empty
#' @keywords internal
#'
IsMatrixEmptyUtil <- function(mat) {
  if (is.null(mat)) {
    return(TRUE)
  }

  if (length(mat) == 0) {
    return(TRUE)
  }

  if (inherits(mat, "dgCMatrix") || inherits(mat, "Matrix")) {
    return(Matrix::nnzero(mat) == 0)
  }

  return(all(mat == 0))
}

#' Ensure matrix has unique row and column names
#'
#' Checks for duplicate row/column names and makes them unique if necessary.
#'
#' @param mat Matrix to check
#' @param fix_rownames If TRUE, make duplicate rownames unique
#' @param fix_colnames If TRUE, make duplicate colnames unique
#' @param verbose Print warnings for duplicates
#'
#' @return Matrix with unique names
#' @keywords internal
#'
EnsureUniqueMatrixNames <- function(
  mat,
  fix_rownames = TRUE,
  fix_colnames = TRUE,
  verbose = FALSE
) {
  if (fix_rownames && !is.null(rownames(mat))) {
    if (anyDuplicated(rownames(mat))) {
      if (verbose) {
        n_dup <- sum(duplicated(rownames(mat)))
        message(
          "Found ", n_dup,
          " duplicate rownames, making unique"
        )
      }
      rownames(mat) <- make.unique(rownames(mat))
    }
  }

  if (fix_colnames && !is.null(colnames(mat))) {
    if (anyDuplicated(colnames(mat))) {
      if (verbose) {
        n_dup <- sum(duplicated(colnames(mat)))
        message(
          "Found ", n_dup,
          " duplicate colnames, making unique"
        )
      }
      colnames(mat) <- make.unique(colnames(mat))
    }
  }

  return(mat)
}

#' Convert dense matrix to sparse if beneficial
#'
#' Converts a dense matrix to sparse format if the sparsity level exceeds
#' a threshold.
#'
#' @param mat Dense matrix
#' @param sparsity_threshold Proportion of zeros required to convert (default 0.5)
#' @param verbose Print diagnostic messages
#'
#' @return Sparse matrix if conversion beneficial, otherwise original matrix
#' @keywords internal
#'
ConvertToSparseIfBeneficial <- function(
  mat,
  sparsity_threshold = 0.5,
  verbose = FALSE
) {
  if (inherits(mat, "dgCMatrix") || inherits(mat, "Matrix")) {
    return(mat)
  }

  n_zeros <- sum(mat == 0)
  sparsity <- n_zeros / length(mat)

  if (sparsity >= sparsity_threshold) {
    if (verbose) {
      message(
        "Converting to sparse matrix (sparsity: ",
        round(sparsity * 100, 1), "%)"
      )
    }
    return(Matrix::Matrix(mat, sparse = TRUE))
  }

  if (verbose) {
    message(
      "Keeping as dense matrix (sparsity: ",
      round(sparsity * 100, 1),
      "% < threshold ",
      sparsity_threshold * 100, "%)"
    )
  }
  return(mat)
}
