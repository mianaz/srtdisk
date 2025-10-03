#' @title Safe H5D Dataset Access
#' @description Functions for safely accessing HDF5 datasets regardless of dimensionality
#' @keywords internal
#'
NULL

#' Safely read an H5D dataset
#'
#' @param dataset An H5D dataset object
#' @param subset Optional subset indices
#' @return The data from the dataset
#' @keywords internal
#'
SafeH5DRead <- function(dataset, subset = NULL) {
  if (!inherits(x = dataset, what = 'H5D')) {
    stop("Input must be an H5D dataset", call. = FALSE)
  }

  ndims <- length(dataset$dims)

  # If subset is provided, use it
  if (!is.null(subset)) {
    if (ndims == 1) {
      return(dataset[subset])
    } else if (ndims == 2) {
      if (is.list(subset)) {
        return(dataset[subset[[1]], subset[[2]]])
      } else {
        # Assume subset is for first dimension, take all of second
        return(dataset[subset, ])
      }
    } else {
      stop("Cannot subset datasets with >2 dimensions", call. = FALSE)
    }
  }

  # No subset - read all data
  if (ndims == 0) {
    # Scalar
    return(dataset$read())
  } else if (ndims == 1) {
    # 1D array
    return(dataset[])
  } else if (ndims == 2) {
    # 2D array - use appropriate access method
    # Try different methods to find one that works
    result <- tryCatch(
      expr = dataset[,],
      error = function(e) {
        tryCatch(
          expr = dataset$read(),
          error = function(e2) {
            # Last resort: explicit indexing
            dataset[1:dataset$dims[1], 1:dataset$dims[2]]
          }
        )
      }
    )

    # Handle special case where 2nd dimension is 1 (should be 1D vector)
    if (!is.null(dim(result)) && ncol(result) == 1) {
      result <- drop(result)
    }

    return(result)
  } else {
    # Higher dimensional - use read()
    return(dataset$read())
  }
}

#' Check if H5D dataset is effectively 1D
#'
#' @param dataset An H5D dataset object
#' @return Logical indicating if dataset is effectively 1D
#' @keywords internal
#'
IsEffectively1D <- function(dataset) {
  if (!inherits(x = dataset, what = 'H5D')) {
    return(FALSE)
  }

  ndims <- length(dataset$dims)

  if (ndims == 1) {
    return(TRUE)
  } else if (ndims == 2) {
    # Check if one dimension is 1
    return(dataset$dims[1] == 1 || dataset$dims[2] == 1)
  } else {
    return(FALSE)
  }
}

