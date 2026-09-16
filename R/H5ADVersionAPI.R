#' @include H5ADVersion.R
NULL

#' Inspect, upgrade or downgrade the on-disk layout of an h5ad file
#'
#' anndata changed its h5ad layout in 0.8.0 (every element became tagged with
#' \code{encoding-type}/\code{encoding-version}, categoricals became groups,
#' pandas nullable dtypes gained a values + mask layout), extended it in 0.11
#' (\code{nullable-string-array}) and made those extensions the default in
#' 0.13 when used with pandas 3 (string columns and the obs/var index are
#' written as nullable string arrays, \code{None} as \code{null}). Older
#' anndata releases, scanpy pinned to them, and every converter derived from
#' SeuratDisk cannot open files written in the newer layouts, while modern
#' anndata still reads the legacy layout.
#'
#' These functions rewrite an h5ad file between layouts \emph{losslessly}:
#' matrices are copied at the HDF5 level, dataframes and \code{uns} are
#' re-encoded, and anything the target layout cannot represent (a nullable
#' integer with missing values, a \code{null} entry, a nullable string index)
#' is written in the closest representable form and recorded in
#' \code{uns/__h5ad_compat_manifest__}. Running the opposite conversion on the
#' result restores the original encodings exactly.
#'
#' @param source Path to the input h5ad file
#' @param dest Path to the output h5ad file
#' @param overwrite Overwrite \code{dest} if it exists
#' @param strings How to write nullable string arrays when upgrading:
#'   \code{"keep"} preserves the \code{nullable-string-array} encoding
#'   (anndata >= 0.11); \code{"categorical"} rewrites them as categoricals
#'   (when they hold missing values) or plain string arrays, which every
#'   anndata >= 0.8 can read.
#' @param verbose Print progress
#'
#' @return \code{UpgradeH5AD} and \code{DowngradeH5AD} return \code{dest}
#'   invisibly. \code{H5ADLayout} returns an object of class
#'   \code{h5ad_layout} describing the layout (\code{"encoded"},
#'   \code{"legacy"} or \code{"compound"}), the oldest anndata release able to
#'   read the file, and which encodings it contains.
#'
#' @section Layouts:
#' \describe{
#'   \item{\code{encoded}}{anndata >= 0.8. Written by \code{UpgradeH5AD}.}
#'   \item{\code{legacy}}{anndata 0.7 (\code{dataframe} 0.1.0 with a
#'     \code{__categories} group referenced by HDF5 object references, no
#'     encoding attributes on arrays, dicts and scalars). Written by
#'     \code{DowngradeH5AD}; readable by anndata 0.7 through 0.13 and by
#'     SeuratDisk.}
#'   \item{\code{compound}}{anndata < 0.7 (obs/var stored as compound
#'     datasets). Read-only; \code{UpgradeH5AD} converts it.}
#' }
#'
#' @examples
#' \dontrun{
#' H5ADLayout("modern.h5ad")
#' DowngradeH5AD("modern.h5ad", "legacy.h5ad")   # for anndata 0.7 / SeuratDisk
#' UpgradeH5AD("legacy.h5ad", "restored.h5ad")   # exact original encodings back
#' }
#'
#' @name H5ADVersion
#' @rdname H5ADVersion
#' @export
UpgradeH5AD <- function(source, dest, overwrite = FALSE,
                        strings = c("keep", "categorical"), verbose = TRUE) {
  strings <- match.arg(strings)
  .h5ad_rewrite(source, dest, target = "encoded", overwrite = overwrite,
                strings = strings, manifest = TRUE, verbose = verbose)
}

#' @rdname H5ADVersion
#' @export
DowngradeH5AD <- function(source, dest, overwrite = FALSE, verbose = TRUE) {
  .h5ad_rewrite(source, dest, target = "legacy", overwrite = overwrite,
                manifest = TRUE, verbose = verbose)
}

#' @param file Path to an h5ad file
#' @rdname H5ADVersion
#' @export
H5ADLayout <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  .h5ad_layout_info(file)
}
