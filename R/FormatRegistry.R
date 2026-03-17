#' @include zzz.R
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Format Registry
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Internal storage for format loaders, savers, and direct conversion paths.
# Populated via RegisterFormat() and RegisterDirectPath() calls in .onLoad().
.format_registry <- new.env(parent = emptyenv())
.format_registry$loaders <- list()
.format_registry$savers <- list()
.format_registry$direct_paths <- list()

#' Register a file format with its load and save functions
#'
#' @param ext Lowercase file extension (e.g., "h5seurat", "h5ad", "loom", "rds")
#' @param loader Function(file, assay, verbose, ...) -> Seurat object
#' @param saver Function(object, dest, overwrite, verbose, ...) -> invisible(dest)
#'
#' @keywords internal
#'
RegisterFormat <- function(ext, loader = NULL, saver = NULL) {
  if (!is.null(loader)) {
    .format_registry$loaders[[ext]] <- loader
  }
  if (!is.null(saver)) {
    .format_registry$savers[[ext]] <- saver
  }
}

#' Register a direct conversion path between two file formats
#'
#' Direct paths bypass the Seurat hub for efficiency (e.g., h5ad <-> h5seurat
#' can be converted via direct HDF5 operations without loading into memory).
#'
#' @param stype Source format extension string
#' @param dtype Destination format extension string
#' @param fn Conversion function(source, dest, assay, overwrite, verbose, ...)
#'
#' @keywords internal
#'
RegisterDirectPath <- function(stype, dtype, fn) {
  key <- paste0(stype, "->", dtype)
  .format_registry$direct_paths[[key]] <- fn
}

#' @keywords internal
GetLoader <- function(ext) {
  .format_registry$loaders[[ext]]
}

#' @keywords internal
GetSaver <- function(ext) {
  .format_registry$savers[[ext]]
}

#' @keywords internal
GetDirectPath <- function(stype, dtype) {
  .format_registry$direct_paths[[paste0(stype, "->", dtype)]]
}

#' List all registered format conversions
#'
#' @return A list with loaders, savers, and direct_paths names
#'
#' @keywords internal
#'
ListFormats <- function() {
  list(
    loaders = names(.format_registry$loaders),
    savers = names(.format_registry$savers),
    direct_paths = names(.format_registry$direct_paths)
  )
}
