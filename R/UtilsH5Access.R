#' HDF5 Access Utility Functions
#'
#' This file contains utility functions for safely accessing HDF5 files,
#' including existence checking with caching, dataset reading with V5
#' compatibility, and feature retrieval.
#'
#' @name UtilsH5Access
#' @keywords internal
NULL

#' Create a cached HDF5 existence checker
#'
#' Creates a closure that checks for HDF5 path existence with result caching
#' to improve performance for repeated checks.
#'
#' @return A function with signature \code{function(group, path)} that checks
#'   if \code{path} exists in \code{group}, using a cache to avoid repeated
#'   HDF5 calls.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' safe_exists <- CreateCachedExistsChecker()
#' if (safe_exists(h5_group, "data/matrix")) {
#'   # path exists
#' }
#' }
#'
CreateCachedExistsChecker <- function() {
  cache <- new.env(hash = TRUE, parent = emptyenv())

  function(group, path) {
    cache_key <- paste0(group$get_obj_name(), "::", path)

    if (exists(cache_key, envir = cache, inherits = FALSE)) {
      return(get(cache_key, envir = cache))
    }

    result <- tryCatch(
      expr = {
        group$exists(name = path)
      },
      error = function(e) {
        path %in% names(group)
      }
    )

    assign(cache_key, result, envir = cache)
    return(result)
  }
}

#' Safely read HDF5 dataset with V5 compatibility
#'
#' Reads character data from an HDF5 dataset, handling V5 2D dummy datasets
#' and providing fallback to special index locations.
#'
#' @param dataset H5D dataset object to read
#' @param h5_group Optional H5Group for V5 feature index lookup
#' @param dataset_name Optional dataset name to identify special cases
#'   (e.g., "features" for V5 feature lookup)
#' @param verbose Print diagnostic messages
#'
#' @return Character vector of dataset values
#' @keywords internal
#'
SafeReadDatasetV5 <- function(
  dataset,
  h5_group = NULL,
  dataset_name = "",
  verbose = FALSE
) {
  if (dataset_name == "features" && !is.null(h5_group)) {
    if (h5_group$exists("meta.data/_index")) {
      if (verbose) {
        message("Using V5 feature index from meta.data/_index")
      }
      return(as.character(h5_group[["meta.data/_index"]][]))
    }
  }

  if (length(dataset$dims) > 1) {
    if (verbose) {
      message("Reading 2D dataset, extracting first column")
    }
    result <- tryCatch(
      expr = {
        if (inherits(dataset, "H5D")) {
          dataset[, 1]
        } else {
          dataset[, 1]
        }
      },
      error = function(e) {
        SafeH5DRead(dataset)
      }
    )
    return(as.character(result))
  }

  return(as.character(dataset[]))
}

#' Get feature names with V5 fallback
#'
#' Retrieves feature names from an HDF5 group, trying V5 location first
#' (meta.data/_index), then falling back to direct features dataset.
#'
#' @param h5_group HDF5 group containing feature information
#' @param verbose Print diagnostic messages
#'
#' @return Character vector of feature names, or NULL if not found
#' @keywords internal
#'
GetFeaturesV5Safe <- function(h5_group, verbose = FALSE) {
  v5_index_exists <- tryCatch(
    h5_group$exists("meta.data/_index"),
    error = function(e) FALSE
  )

  if (v5_index_exists) {
    if (verbose) {
      message("Using V5 feature index from meta.data/_index")
    }
    return(as.character(h5_group[["meta.data/_index"]][]))
  }

  features_exists <- tryCatch(
    h5_group$exists("features"),
    error = function(e) FALSE
  )

  if (features_exists) {
    if (verbose) {
      message("Using direct features dataset")
    }
    result <- tryCatch(
      expr = as.character(h5_group[["features"]][]),
      error = function(e) {
        if (length(h5_group[["features"]]$dims) > 1) {
          return(as.character(h5_group[["features"]][, 1]))
        }
        return(NULL)
      }
    )
    return(result)
  }

  if (verbose) {
    message("No features found in h5_group")
  }
  return(NULL)
}

#' Set V5 feature index in HDF5 group
#'
#' Stores feature names in the V5 location (meta.data/_index) if not already
#' present.
#'
#' @param h5_group HDF5 group to store features in
#' @param features Character vector of feature names
#' @param verbose Print diagnostic messages
#'
#' @return NULL (invisible)
#' @keywords internal
#'
SetFeaturesV5 <- function(h5_group, features, verbose = FALSE) {
  if (is.null(features) || length(features) == 0) {
    if (verbose) {
      message("No features provided, skipping V5 index creation")
    }
    return(invisible(NULL))
  }

  if (!h5_group$exists("meta.data")) {
    h5_group$create_group("meta.data")
  }

  if (!h5_group[["meta.data"]]$exists("_index")) {
    h5_group[["meta.data"]][["_index"]] <- features
    if (verbose) {
      message("Stored ", length(features), " feature names in V5 index")
    }
  }

  return(invisible(NULL))
}

#' Check if HDF5 path exists with caching (deprecated - use CreateCachedExistsChecker)
#'
#' This function is deprecated. Use \code{CreateCachedExistsChecker()} to create
#' a cached existence checker instead.
#'
#' @param group H5Group object
#' @param path Path to check
#' @param cache_env Environment to use for caching
#'
#' @return Logical indicating if path exists
#' @keywords internal
#'
#' @note Deprecated: Use CreateCachedExistsChecker instead
#'
SafeExistsDeprecated <- function(group, path, cache_env = NULL) {
  .Deprecated("CreateCachedExistsChecker")

  if (is.null(cache_env)) {
    return(tryCatch(
      expr = group$exists(name = path),
      error = function(e) path %in% names(group)
    ))
  }

  cache_key <- paste0(group$get_obj_name(), "::", path)

  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    return(get(cache_key, envir = cache_env))
  }

  result <- tryCatch(
    expr = group$exists(name = path),
    error = function(e) path %in% names(group)
  )

  assign(cache_key, result, envir = cache_env)
  return(result)
}
