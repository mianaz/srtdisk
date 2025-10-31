#' Assay Compatibility Utility Functions
#'
#' This file contains utility functions for handling compatibility between
#' Seurat V4 Assay objects and V5 Assay5 objects. These functions provide
#' unified interfaces for operations that differ between versions.
#'
#' @name UtilsAssayCompat
#' @keywords internal
NULL

#' Check if object is Assay5 (V5 format)
#'
#' @param object Object to check
#' @return Logical indicating if object inherits from Assay5
#' @keywords internal
#'
IsAssay5 <- function(object) {
  inherits(object, "Assay5")
}

#' Get assay data with V4/V5 compatibility
#'
#' Retrieves data from an Assay or Assay5 object, automatically using the
#' appropriate parameter name (layer for V5, slot for V4).
#'
#' @param object Assay or Assay5 object
#' @param layer_or_slot Name of layer (V5) or slot (V4) to retrieve
#' @param verbose Print diagnostic messages
#'
#' @return Matrix or sparse matrix of assay data
#' @keywords internal
#'
#' @seealso \code{\link{SafeGetAssayData}} in V5Compatibility.R for the
#'   original implementation
#'
GetAssayDataCompat <- function(object, layer_or_slot, verbose = FALSE) {
  if (IsAssay5(object)) {
    if (verbose) {
      message("Using layer parameter for Assay5: ", layer_or_slot)
    }
    return(tryCatch(
      expr = GetAssayData(object = object, layer = layer_or_slot),
      error = function(e) {
        if (verbose) {
          message("Error getting layer '", layer_or_slot, "': ", conditionMessage(e))
        }
        return(NULL)
      }
    ))
  } else {
    if (verbose) {
      message("Using slot parameter for Assay: ", layer_or_slot)
    }
    return(tryCatch(
      expr = GetAssayData(object = object, slot = layer_or_slot),
      error = function(e) {
        if (verbose) {
          message("Error getting slot '", layer_or_slot, "': ", conditionMessage(e))
        }
        return(NULL)
      }
    ))
  }
}

#' Set assay data with V4/V5 compatibility
#'
#' Sets data in an Assay or Assay5 object, automatically using the appropriate
#' parameter name (layer for V5, slot for V4).
#'
#' @param object Assay or Assay5 object
#' @param layer_or_slot Name of layer (V5) or slot (V4) to set
#' @param new.data New data to set
#' @param verbose Print diagnostic messages
#'
#' @return Modified Assay or Assay5 object
#' @keywords internal
#'
SetAssayDataCompat <- function(object, layer_or_slot, new.data, verbose = FALSE) {
  if (IsAssay5(object)) {
    if (verbose) {
      message("Setting layer for Assay5: ", layer_or_slot)
    }
    return(SetAssayData(
      object = object,
      layer = layer_or_slot,
      new.data = new.data
    ))
  } else {
    if (verbose) {
      message("Setting slot for Assay: ", layer_or_slot)
    }
    return(SetAssayData(
      object = object,
      slot = layer_or_slot,
      new.data = new.data
    ))
  }
}

#' Get available layers or slots from assay
#'
#' Returns list of available data layers (V5) or slots (V4) from an Assay
#' or Assay5 object.
#'
#' @param object Assay or Assay5 object
#' @param verbose Print diagnostic messages
#'
#' @return Character vector of layer/slot names
#' @keywords internal
#'
GetAssayLayersCompat <- function(object, verbose = FALSE) {
  if (IsAssay5(object)) {
    if (verbose) {
      message("Getting layers from Assay5")
    }
    layers <- tryCatch(
      expr = {
        if (exists("Layers", where = "package:SeuratObject", mode = "function")) {
          SeuratObject::Layers(object)
        } else {
          slotNames(object)
        }
      },
      error = function(e) {
        if (verbose) {
          message("Error getting layers: ", conditionMessage(e))
        }
        return(character(0))
      }
    )
    return(layers)
  } else {
    if (verbose) {
      message("Getting slots from Assay (V4)")
    }
    return(c("counts", "data", "scale.data"))
  }
}

#' Execute function with V4/V5 compatibility wrapper
#'
#' Executes a function with appropriate parameter names based on whether the
#' object is an Assay (V4) or Assay5 object.
#'
#' @param object Assay or Assay5 object
#' @param func Function to execute
#' @param layer_or_slot Layer/slot name to pass to function
#' @param use_layer_param Logical, if TRUE uses 'layer' parameter name,
#'   if FALSE uses 'slot' parameter name. If NULL (default), automatically
#'   determined from object type.
#' @param ... Additional arguments to pass to func
#'
#' @return Result from function call
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Automatically detects V4 vs V5 and uses appropriate parameter
#' result <- WithAssayCompat(
#'   object = assay_obj,
#'   func = GetAssayData,
#'   layer_or_slot = "counts"
#' )
#' }
#'
WithAssayCompat <- function(
  object,
  func,
  layer_or_slot,
  use_layer_param = NULL,
  ...
) {
  if (is.null(use_layer_param)) {
    use_layer_param <- IsAssay5(object)
  }

  params <- list(object = object, ...)

  if (use_layer_param) {
    params$layer <- layer_or_slot
  } else {
    params$slot <- layer_or_slot
  }

  do.call(func, params)
}

#' Check if assay requires S4 slot reconstruction
#'
#' Determines if an assay object requires S4 slot reconstruction (V4 and below)
#' or should be left intact (V5).
#'
#' @param object Assay or Assay5 object
#' @param h5_group Optional HDF5 group to check for s4class attribute
#' @param verbose Print diagnostic messages
#'
#' @return Logical indicating if S4 reconstruction is needed
#' @keywords internal
#'
RequiresS4Reconstruction <- function(object, h5_group = NULL, verbose = FALSE) {
  if (IsAssay5(object)) {
    if (verbose) {
      message("Assay5 detected - skipping S4 reconstruction")
    }
    return(FALSE)
  }

  if (!is.null(h5_group)) {
    if (h5_group$attr_exists("s4class")) {
      if (verbose) {
        message("Legacy Assay with s4class attribute - needs reconstruction")
      }
      return(TRUE)
    }
  }

  if (verbose) {
    message("No S4 reconstruction required")
  }
  return(FALSE)
}

#' Get slot mapping for Seurat version
#'
#' Returns a mapping of slot names appropriate for different Seurat versions
#' and formats (h5Seurat, h5ad, etc.)
#'
#' @param format Character, format type ("h5seurat", "h5ad")
#' @param seurat_version Character, Seurat version (e.g., "4.0", "5.0")
#'
#' @return Named list mapping generic slot names to format-specific names
#' @keywords internal
#'
GetSlotMapping <- function(format = "h5seurat", seurat_version = "5.0") {
  format <- tolower(format)

  if (format == "h5seurat") {
    if (seurat_version >= "5.0") {
      return(list(
        counts = "layers/counts",
        data = "layers/data",
        scale.data = "layers/scale.data",
        features = "meta.data/_index",
        cells = "obs_names"
      ))
    } else {
      return(list(
        counts = "counts",
        data = "data",
        scale.data = "scale.data",
        features = "features",
        cells = "cell.names"
      ))
    }
  } else if (format == "h5ad") {
    return(list(
      counts = "layers/counts",
      data = "X",
      features = "var/_index",
      cells = "obs/_index"
    ))
  }

  return(list())
}
