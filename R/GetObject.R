#' @include zzz.R
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
GetAssays <- function(assays, index) {
  index.assays <- setdiff(x = names(x = index), y = c('global', 'no.assay'))
  assay.slots <- c('counts', 'data', 'scale.data')
  assay.msg <- paste(
    'Assay specification must include either the name of an assay',
    'or one or more assay slots'
  )
  assays <- assays %||% index.assays
  if (!is.null(x = names(x = assays))) {
    assays <- as.list(x = assays)
  }
  if (!is.list(x = assays)) {
    if (any(assays %in% index.assays) && any(assays %in% assay.slots)) {
      stop("Ambiguous assays", call. = FALSE)
    } else if (any(assays %in% index.assays)) {
      assays <- assays[assays %in% index.assays]
      assays <- sapply(
        X = assays,
        FUN = function(...) {
          return(assay.slots)
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
    } else {
      assays <- assays[assays %in% assay.slots]
      if (length(x = assays) < 1) {
        stop(assay.msg, call. = FALSE)
      }
      assays <- list(assays)
      assays <- rep_len(x = assays, length.out = length(x = index.assays))
      names(x = assays) <- index.assays
    }
  } else {
    for (i in 1:length(x = assays)) {
      assay.name <- names(x = assays)[i] %||% index.assays[i] %||% ''
      if (!assay.name %in% index.assays) {
        if (assays[[i]][1] %in% index.assays) {
          assay.name <- assays[[i]][1]
        } else if (any(assays[[i]] %in% assay.slots)) {
          assay.name <- hdf5r::h5attr(x = file, which = 'active.assay')
        }
      }
      if (nchar(x = assay.name) < 0 || !assay.name %in% index.assays) {
        stop(assay.msg, call. = FALSE)
      }
      assay.content <- assays[[i]]
      if (assay.content[1] %in% index.assays) {
        assay.content <- assay.slots
      } else {
        assay.content <- assay.content[assay.content %in% assay.slots]
        if (length(x = assay.content) < 1) {
          stop(assay.msg, call. = FALSE)
        }
      }
      assays[i] <- list(assay.content)
      names(x = assays)[i] <- assay.name
    }
  }
  assays.checked <- assays
  unique.assays <- unique(x = names(x = assays.checked))
  assays <- vector(mode = 'list', length = length(x = unique.assays))
  names(x = assays) <- unique.assays
  for (i in unique.assays) {
    assays.use <- which(x = names(x = assays.checked) == i)
    slots.use <- unique(x = unlist(
      x = assays.checked[assays.use],
      use.names = FALSE
    ))
    slots.use <- slots.use[match(x = names(x = index[[i]]$slots), table = slots.use)]

    # Filter to only slots that are TRUE in the index
    # Use which() to get indices where index[[i]]$slots is TRUE
    slots.use <- as.character(x = na.omit(object = slots.use[which(index[[i]]$slots)]))

    # Fallback: if index incorrectly marks all slots as FALSE, use all requested slots
    # This handles cases where SaveH5Seurat creates an invalid index
    if (length(slots.use) == 0 && any(names(index[[i]]$slots) %in% assay.slots)) {
      # Use the originally requested slots that exist in the index
      slots.use <- names(index[[i]]$slots)[names(index[[i]]$slots) %in% assay.slots]
    }

    assays[[i]] <- slots.use
  }
  # Remove assays that have no matching slots (instead of throwing error)
  # This allows loading specific slots from only the assays that have them
  assays_before_filter <- assays
  assays <- Filter(function(x) length(x) > 0, assays)

  # Only error if NO assays remain
  if (length(assays) == 0) {
    # Add detailed error message for debugging
    msg <- "No assays found with the requested slots/layers\n"
    msg <- paste0(msg, "  Available assays in index: ", paste(index.assays, collapse = ", "), "\n")
    msg <- paste0(msg, "  Assays before filtering: ", paste(names(assays_before_filter), collapse = ", "), "\n")
    for (a in names(assays_before_filter)) {
      msg <- paste0(msg, "    ", a, ": ", paste(assays_before_filter[[a]], collapse = ", "), "\n")
    }
    stop(msg, call. = FALSE)
  }

  return(assays)
}

#' @return \code{GetCommands}: A vector of command log names that are derived
#' from an assay in \code{assay}
#'
#' @rdname GetObject
#'
GetCommands <- function(index, assays = NULL) {
  assays <- GetAssays(assays = assays, index = index)
  return(unique(x = unlist(x = lapply(
    X = names(x = assays),
    FUN = function(x) {
      return(index[[x]]$commands)
    }
  ))))
}

#' @return \code{GetDimReducs}: A vector of reduction names that are derived
#' from an assay in \code{assays} or global dimensional reductions
#'
#' @rdname GetObject
#'
GetDimReducs <- function(reductions, index, assays = NULL) {
  if (isFALSE(x = reductions)) {
    return(NULL)
  } else if (!is.null(x = reductions) && all(is.na(x = reductions))) {
    return(index$global$reductions)
  }
  if (is.null(x = reductions)) {
    reductions <- unique(x = unlist(x = lapply(
      X = setdiff(x = names(x = index), y = 'no.assay'),
      FUN = function(x) {
        if (x == 'global') {
          return(index[[x]]$reductions)
        }
        return(names(x = index[[x]]$reductions))
      }
    )))
  }
  if (isTRUE(x = getOption(x = 'SeuratDisk.dimreducs.allglobal', default = FALSE))) {
    return(reductions)
  }
  assays <- GetAssays(assays = assays, index = index)
  assays.reducs <- lapply(
    X = names(x = assays),
    FUN = function(x) {
      return(names(x = index[[x]]$reductions))
    }
  )
  assays.reducs <- unique(x = c(unlist(x = assays.reducs), index$global$reductions))
  return(intersect(x = reductions, y = assays.reducs))
}

#' @return \code{GetGraphs}: A vector of graph names that are derived from an
#' assay in \code{assays}
#'
#' @rdname GetObject
#'
GetGraphs <- function(graphs, index, assays = NULL) {
  if (isFALSE(x = graphs)) {
    return(NULL)
  } else if (is.null(x = graphs)) {
    graphs <- unique(x = unlist(x = lapply(
      X = setdiff(x = names(x = index), y = c('global', 'no.assay')),
      FUN = function(x) {
        return(index[[x]]$graphs)
      }
    )))
  }
  assays <- GetAssays(assays = assays, index = index)
  assays.graphs <- unique(x = unlist(x = lapply(
    X = names(x = assays),
    FUN = function(x) {
      return(index[[x]]$graphs)
    }
  )))
  return(intersect(x = graphs, y = assays.graphs))
}

#' @return \code{GetImages}: A vector of image names
#'
#' @rdname GetObject
#'
GetImages <- function(images, index, assays = NULL) {
  if (isFALSE(x = images)) {
    return(NULL)
  } else if (!is.null(x = images) && all(is.na(x = images))) {
    return(index$global$images)
  } else if (is.null(x = images)) {
    images <- unique(x = unlist(x = lapply(
      X = names(x = index),
      FUN = function(x) {
        return(index[[x]]$images)
      }
    )))
  }
  assays <- GetAssays(assays = assays, index = index)
  assays.images <- lapply(
    X = names(x = assays),
    FUN = function(x) {
      return(index[[x]]$images)
    }
  )
  assays.images <- unique(x = c(unlist(x = assays.images), index$global$images))
  return(intersect(x = images, y = assays.images))
}

#' @return \code{GetNeighbors}: A vector of neighbor names
#'
#' @rdname GetObject
#'
GetNeighbors <- function(neighbors, index) {
  if (isFALSE(x = neighbors)) {
    return(NULL)
  } else if (is.null(x = neighbors)) {
    return(index[["global"]]$neighbors)
  } else {
    return(intersect(x = neighbors, y = index[["global"]]$neighbors))
  }
}
