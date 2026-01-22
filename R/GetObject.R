#' @include zzz.R
#' @importFrom stats setNames
#'
NULL

#' Figure out which objects to load from an h5Seurat file
#'
#' @inheritParams LoadH5Seurat
#' @param index An h5Seurat index (\code{\link{h5SI}}) object
#'
#' @seealso \code{\link{LoadH5Seurat}}
#'
#' @rdname GetObject
#' @name GetObject
#'
#' @keywords internal
#'
NULL

#' @return \code{GetAssays}: A named list where each entry is a vector
#' describing the slots of an assay to load and the names are the assays to load
#'
#' @rdname GetObject
#'
# Robust GetAssays implementation
GetAssays <- function(assays = NULL, index) {
  # index is expected to contain element `assays` which is a named list of assay metadata
  all_assays <- setdiff(names(index), c('global', 'no.assay'))
  if (is.null(all_assays) || length(all_assays) == 0) {
    stop("No assays found in the h5Seurat index.")
  }

  # helper: return named list of length n with value val
  fill_named <- function(names_vec, val = NULL) {
    setNames(rep(list(val), length(names_vec)), names_vec)
  }

  # NULL -> load all assays with default slots (NULL)
  if (is.null(assays)) {
    return(fill_named(all_assays, NULL))
  }

  # If assays is already a named list, validate names and normalize values
  if (is.list(assays) && !is.null(names(assays)) && any(names(assays) != "")) {
    requested <- intersect(names(assays), all_assays)
    if (length(requested) == 0) {
      stop(
        "No assays found matching the requested assay names: ",
        paste(names(assays), collapse = ", "),
        ". Available assays: ", paste(all_assays, collapse = ", ")
      )
    }
    out <- assays[requested]
    out <- lapply(out, function(x) {
      if (is.logical(x) && isTRUE(x)) return(NULL)
      if (identical(x, "")) return(NULL)
      x
    })
    return(out)
  }

  # If character vector:
  if (is.character(assays)) {
    # Special sentinel: "layers-only" => select assays that have any layers
    if (length(assays) == 1 && identical(assays, "layers-only")) {
      has_layers <- vapply(all_assays, function(a) {
        meta <- index[[a]]
        layers_meta <- meta[['layers']]
        !is.null(layers_meta) && length(layers_meta) > 0
      }, logical(1))
      selected <- all_assays[has_layers]
      if (length(selected) == 0) {
        stop(
          "Requested 'layers-only' but no assays in the file contain layers. Available assays: ",
          paste(all_assays, collapse = ", ")
        )
      }
      # Return named list where each assay requests its available layers (caller may expect "layers" or explicit names)
      # Prefer to return the literal layer names if available
      return(setNames(lapply(selected, function(a) {
        layers_meta <- index[[a]][['layers']]
        if (!is.null(layers_meta) && length(layers_meta) > 0) {
          return(as.character(layers_meta))
        }
        # fallback to sentinel if meta doesn't list them explicitly
        "layers"
      }), selected))
    }

    # If any of the provided strings match assay names -> treat as assay names to load (NULL slots)
    if (any(assays %in% all_assays)) {
      selected <- intersect(as.character(assays), all_assays)
      if (length(selected) == 0) {
        stop("No assay names found matching: ", paste(as.character(assays), collapse = ", "))
      }
      return(fill_named(selected, NULL))
    }

    # Otherwise treat the character vector as slot names to load for every assay
    requested_slots <- as.character(assays)
    return(setNames(rep(list(requested_slots), length(all_assays)), all_assays))
  }

  stop(
    "Unsupported 'assays' argument type. Expect NULL, character vector, or named list. ",
    "Available assays: ", paste(all_assays, collapse = ", ")
  )
}

#' @return \code{GetCommands}: A vector of command log names that are derived
#' from an assay in \code{assay}
#'
#' @rdname GetObject
#'
GetCommands <- function(index, assays = NULL) {
  assays <- GetAssays(assays = assays, index = index)
  unique(unlist(lapply(names(assays), function(x) index[[x]]$commands)))
}

#' @return \code{GetDimReducs}: A vector of reduction names that are derived
#' from an assay in \code{assays} or global dimensional reductions
#'
#' @rdname GetObject
#'
GetDimReducs <- function(reductions, index, assays = NULL) {
  if (isFALSE(reductions)) {
    return(NULL)
  }
  if (!is.null(reductions) && all(is.na(reductions))) {
    return(index$global$reductions)
  }
  if (is.null(reductions)) {
    reductions <- unique(unlist(lapply(
      setdiff(names(index), 'no.assay'),
      function(x) {
        if (x == 'global') index[[x]]$reductions else names(index[[x]]$reductions)
      }
    )))
  }
  if (isTRUE(getOption('SeuratDisk.dimreducs.allglobal', default = FALSE))) {
    return(reductions)
  }
  assays <- GetAssays(assays = assays, index = index)
  assays.reducs <- unique(c(
    unlist(lapply(names(assays), function(x) names(index[[x]]$reductions))),
    index$global$reductions
  ))
  intersect(reductions, assays.reducs)
}

#' @return \code{GetGraphs}: A vector of graph names that are derived from an
#' assay in \code{assays}
#'
#' @rdname GetObject
#'
GetGraphs <- function(graphs, index, assays = NULL) {
  if (isFALSE(graphs)) {
    return(NULL)
  }
  if (is.null(graphs)) {
    graphs <- unique(unlist(lapply(
      setdiff(names(index), c('global', 'no.assay')),
      function(x) index[[x]]$graphs
    )))
  }
  assays <- GetAssays(assays = assays, index = index)
  assays.graphs <- unique(unlist(lapply(names(assays), function(x) index[[x]]$graphs)))
  intersect(graphs, assays.graphs)
}

#' @return \code{GetImages}: A vector of image names
#'
#' @rdname GetObject
#'
GetImages <- function(images, index, assays = NULL) {
  if (isFALSE(images)) {
    return(NULL)
  }
  if (!is.null(images) && all(is.na(images))) {
    return(index$global$images)
  }
  if (is.null(images)) {
    images <- unique(unlist(lapply(names(index), function(x) index[[x]]$images)))
  }
  assays <- GetAssays(assays = assays, index = index)
  assays.images <- unique(c(
    unlist(lapply(names(assays), function(x) index[[x]]$images)),
    index$global$images
  ))
  intersect(images, assays.images)
}

#' @return \code{GetNeighbors}: A vector of neighbor names
#'
#' @rdname GetObject
#'
GetNeighbors <- function(neighbors, index) {
  if (isFALSE(neighbors)) {
    return(NULL)
  }
  if (is.null(neighbors)) {
    return(index$global$neighbors)
  }
  intersect(neighbors, index$global$neighbors)
}
