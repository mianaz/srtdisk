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
  # In SeuratObject 5.0+, the 'slot' argument is defunct - always use 'layer'
  # This works for both Assay and Assay5 objects in modern SeuratObject
  if (verbose) {
    message("Getting layer/slot '", layer_or_slot, "' from ", class(object)[1])
  }
  return(tryCatch(
    expr = GetAssayData(object = object, layer = layer_or_slot),
    error = function(e) {
      if (verbose) {
        message("Error getting '", layer_or_slot, "': ", conditionMessage(e))
      }
      return(NULL)
    }
  ))
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
  # In SeuratObject 5.0+, the 'slot' argument is defunct - always use 'layer'
  # This works for both Assay and Assay5 objects in modern SeuratObject
  if (verbose) {
    message("Setting layer/slot '", layer_or_slot, "' in ", class(object)[1])
  }
  return(SetAssayData(
    object = object,
    layer = layer_or_slot,
    new.data = new.data
  ))
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

