#' @include zzz.R
#' @include h5Seurat.R
#'
NULL

#' Connect to a single-cell HDF5 dataset
#'
#' @param filename Name of on-disk file
#' @param type Type of single-cell dataset to connect as; choose from:
#' \itemize{
#'  \item h5seurat
#' }
#' Leave as \code{NULL} to guess type from file extension
#' @param mode Mode to connect to data as; choose from:
#' \describe{
#'  \item{r}{Open existing dataset in read-only mode}
#'  \item{r+}{Open existing dataset in read/write mode}
#' }
#' @param force Force a connection if validation steps fail; returns a
#' \code{\link[hdf5r]{H5File}} object
#'
#' @return An object of class \code{type}, opened in mode \code{mode}
#'
#' @importFrom hdf5r H5File
#'
#' @export
#'
Connect <- function(
  filename,
  type = NULL,
  mode = c('r', 'r+'),
  force = FALSE
) {
  type <- type %||% FileType(file = filename)
  mode <- match.arg(arg = mode)
  if (!file.exists(filename)) {
    stop("Cannot find ", type, " file ", filename, call. = FALSE)
  }
  return(tryCatch(
    expr = {
      cls <- GetSCDisk(r6class = type)
      cls$new(filename = filename, mode = mode)
    },
    error = function(err) {
      if (!isTRUE(x = force)) {
        stop(err$message, call. = FALSE)
      }
      # Silently fall back for HDF5-based formats without R6 classes (h5ad, h5mu)
      # These are opened as raw H5File objects and handled by format-specific code
      return(H5File$new(filename = filename, mode = mode))
    }
  ))
}
