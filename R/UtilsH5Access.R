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

#' Resolve a nested HDF5 path to its object
#'
#' @param base_group The starting H5Group
#' @param path A path string, possibly containing "/" for nested paths
#' @return The H5 object at the resolved path, or NULL if not found
#' @keywords internal
ResolveNestedH5Path <- function(base_group, path) {
  if (!inherits(x = base_group, what = c('H5Group', 'H5File'))) {
    stop("base_group must be an H5Group or H5File", call. = FALSE)
  }
  if (!grepl("/", path)) {
    if (base_group$exists(name = path)) {
      return(base_group[[path]])
    }
    return(NULL)
  }
  parts <- strsplit(path, "/")[[1]]
  current_obj <- base_group
  for (part in parts) {
    if (!inherits(x = current_obj, what = c('H5Group', 'H5File'))) return(NULL)
    if (!current_obj$exists(name = part)) return(NULL)
    current_obj <- current_obj[[part]]
  }
  return(current_obj)
}

#' Copy HDF5 matrix data handling nested paths
#'
#' @param src_group The source H5Group containing the data
#' @param src_path Path to the source data (can be nested)
#' @param dst_loc Destination H5File or H5Group
#' @param dst_name Name for the copied data at destination
#' @param verbose Show progress messages
#' @return Invisible NULL
#' @keywords internal
CopyH5MatrixData <- function(src_group, src_path, dst_loc, dst_name, verbose = FALSE) {
  src_obj <- ResolveNestedH5Path(base_group = src_group, path = src_path)
  if (is.null(src_obj)) {
    stop("Source path '", src_path, "' not found", call. = FALSE)
  }
  if (dst_loc$exists(name = dst_name)) {
    dst_loc$link_delete(name = dst_name)
  }

  if (inherits(x = src_obj, what = 'H5Group')) {
    # Sparse matrix: copy group structure
    dst_group <- dst_loc$create_group(name = dst_name)
    for (comp in c('data', 'indices', 'indptr')) {
      if (src_obj$exists(name = comp)) {
        src_obj$obj_copy_to(dst_loc = dst_group, dst_name = comp, src_name = comp)
      }
    }
    for (attr_name in c('dims', 'shape', 'encoding-type', 'encoding-version',
                        'h5sparse_format', 'h5sparse_shape')) {
      if (src_obj$attr_exists(attr_name = attr_name)) {
        attr_val <- h5attr(x = src_obj, which = attr_name)
        dst_group$create_attr(attr_name = attr_name, robj = attr_val,
                              dtype = GuessDType(x = attr_val))
      }
    }
  } else if (inherits(x = src_obj, what = 'H5D')) {
    # Dense matrix: copy from parent
    parent_obj <- GetParent(x = src_obj)
    dataset_name <- basename(path = src_obj$get_obj_name())
    parent_obj$obj_copy_to(dst_loc = dst_loc, dst_name = dst_name, src_name = dataset_name)
  }
  invisible(NULL)
}

