#' @include SeuratVersion.R
NULL

#' Upgrade or downgrade a Seurat object between the v3/v4 and v5 generations
#'
#' Seurat v5 (SeuratObject >= 5.0) stores expression data in \code{Assay5}
#' objects (named layers, per-layer cell and feature sets, on-disk matrices)
#' and, since Seurat 5.1, Visium images in the FOV-based \code{VisiumV2}
#' class. Objects built that way cannot be loaded by Seurat v3/v4, and the
#' coercion SeuratObject ships (\code{as(x, "Assay")}) is lossy: split
#' layers are joined, layers other than \code{counts}/\code{data}/
#' \code{scale.data} are dropped, the default layer is forgotten and on-disk
#' matrices are pulled into memory.
#'
#' \code{DowngradeSeurat} converts every \code{Assay5} to a v3-style
#' \code{Assay} and every \code{VisiumV2} image to \code{VisiumV1}, and stores
#' everything the older classes cannot hold in a sidecar list under
#' \code{misc$.seurat_version_sidecar} (of the assay, or of the object for
#' images): split-layer cell membership, extra layers, the default layer, the
#' original version stamp, image boundaries. \code{UpgradeSeurat} converts back
#' and consumes the sidecar, so that
#' \code{UpgradeSeurat(DowngradeSeurat(x))} reproduces \code{x}. Assay
#' subclasses (\code{SCTAssay}, \code{ChromatinAssay}) are upgraded to
#' \code{Assay5} with their extra slots kept in the sidecar and restored on
#' downgrade.
#'
#' Both functions also accept a file path (\code{.rds} or \code{.h5seurat}):
#' the object is read, converted and written to \code{dest}.
#'
#' @param object A \code{Seurat} object, or the path to an \code{.rds} /
#'   \code{.h5seurat} file holding one
#' @param dest When \code{object} is a path: where to write the converted
#'   object (\code{.rds} or \code{.h5seurat})
#' @param to Target generation: \code{"v5"} (\code{Assay5}), \code{"v4"} or
#'   \code{"v3"} (\code{Assay}; the two differ only in the version stamp)
#' @param overwrite Overwrite \code{dest} if it exists
#' @param verbose Print what was converted and what went into the sidecar
#'
#' @return The converted \code{Seurat} object, or \code{dest} (invisibly)
#'   when converting a file. \code{SeuratGeneration} returns an object of
#'   class \code{seurat_generation} describing the assay and image classes
#'   in use.
#'
#' @section What is lossless:
#' Layers (including split layers such as \code{counts.sample1} and
#' non-standard names), the default layer, per-layer cell/feature sets,
#' feature-level metadata, variable features, keys, misc, assay subclasses,
#' Visium image boundaries and the version stamp all round-trip exactly.
#' On-disk (BPCells) layers are loaded into memory on downgrade and stay in
#' memory after the upgrade; the sidecar records which layers were on disk.
#'
#' @examples
#' \dontrun{
#' v4 <- DowngradeSeurat(pbmc_v5)          # loadable by Seurat v4
#' saveRDS(v4, "pbmc_v4.rds")
#' v5 <- UpgradeSeurat(v4)                 # identical layers to pbmc_v5
#' UpgradeSeurat("old_v3.rds", dest = "new_v5.h5seurat")
#' SeuratGeneration(v5)
#' }
#'
#' @name SeuratVersion
#' @rdname SeuratVersion
#' @export
UpgradeSeurat <- function(object, dest = NULL, to = "v5", overwrite = FALSE, verbose = TRUE) {
  to <- match.arg(to, c("v5"))
  if (is.character(object)) {
    if (is.null(dest)) stop("`dest` is required when `object` is a file path", call. = FALSE)
    return(.srt_convert_file(object, dest, to = to, overwrite = overwrite, verbose = verbose,
                             reader = function(f) LoadH5Seurat(f, verbose = verbose),
                             writer = function(x, f) SaveH5Seurat(x, filename = f, overwrite = overwrite, verbose = verbose)))
  }
  .srt_convert_object(object, to = to, verbose = verbose)
}

#' @rdname SeuratVersion
#' @export
DowngradeSeurat <- function(object, dest = NULL, to = c("v4", "v3"), overwrite = FALSE, verbose = TRUE) {
  to <- match.arg(to)
  if (is.character(object)) {
    if (is.null(dest)) stop("`dest` is required when `object` is a file path", call. = FALSE)
    return(.srt_convert_file(object, dest, to = to, overwrite = overwrite, verbose = verbose,
                             reader = function(f) LoadH5Seurat(f, verbose = verbose),
                             writer = function(x, f) SaveH5Seurat(x, filename = f, overwrite = overwrite, verbose = verbose)))
  }
  .srt_convert_object(object, to = to, verbose = verbose)
}

#' @rdname SeuratVersion
#' @export
SeuratGeneration <- function(object) {
  .srt_object_version(object)
}
