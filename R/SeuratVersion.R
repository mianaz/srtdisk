#' @include zzz.R
#' @importFrom methods as is new slot slot<- slotNames validObject
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lossless Seurat object upgrade / downgrade
#
# Seurat v5 (SeuratObject >= 5.0.0) replaced the v3/v4 `Assay` class with
# `Assay5` (arbitrary named layers, per-layer cell/feature sets, on-disk
# matrices) and Seurat 5.1 replaced `VisiumV1` with the FOV-based `VisiumV2`
# for 10x Visium images. A v5 object cannot be loaded by Seurat v4, and
# SeuratObject's own `as(x, "Assay")` coercion is lossy: it joins split
# layers, drops every layer that is not counts/data/scale.data, forgets the
# default layer and converts on-disk matrices to memory.
#
# The functions below convert objects between the two generations and keep
# everything the target generation cannot hold in a *sidecar* stored in the
# `misc` slot of the assay (or of the object, for images). Converting back
# consumes the sidecar and rebuilds the original structure exactly. The
# sidecar is a plain list, so it survives RDS / h5Seurat round trips through
# either package.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.SRT_SIDECAR <- ".seurat_version_sidecar"
.SRT_STANDARD_LAYERS <- c("counts", "data", "scale.data")

# ---- helpers ----------------------------------------------------------------

.srt_layer_base <- function(layer) {
  # "counts.sample1" -> "counts"; "scale.data" -> "scale.data"; "foo" -> NA
  for (std in .SRT_STANDARD_LAYERS) {
    if (identical(layer, std) || startsWith(layer, paste0(std, "."))) return(std)
  }
  NA_character_
}

.srt_is_ondisk <- function(x) {
  inherits(x, c("IterableMatrix", "DelayedMatrix", "DelayedArray", "HDF5Matrix"))
}

.srt_get_sidecar <- function(x) {
  m <- tryCatch(slot(x, "misc"), error = function(e) NULL)
  if (is.list(m)) m[[.SRT_SIDECAR]] else NULL
}

.srt_set_sidecar <- function(x, sidecar) {
  m <- tryCatch(slot(x, "misc"), error = function(e) list())
  if (!is.list(m)) m <- list()
  if (is.null(sidecar)) {
    m[[.SRT_SIDECAR]] <- NULL
  } else {
    m[[.SRT_SIDECAR]] <- sidecar
  }
  slot(x, "misc") <- m
  x
}

# The default layer SeuratObject's own Assay -> Assay5 coercion would pick for
# a given layer set; only a different default needs recording in the sidecar
.srt_expected_default <- function(layers) {
  if ("counts" %in% layers && !"scale.data" %in% layers) "counts" else "data"
}

.srt_matrix_identical <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  if (!identical(dim(a), dim(b))) return(FALSE)
  if (inherits(a, "dgCMatrix") && inherits(b, "dgCMatrix")) {
    return(identical(a@p, b@p) && identical(a@i, b@i) && isTRUE(all.equal(a@x, b@x)))
  }
  isTRUE(all.equal(as.matrix(a), as.matrix(b), check.attributes = FALSE))
}

# ---- assay: v3/v4 Assay -> v5 Assay5 --------------------------------------------

#' Convert a v3/v4 \code{Assay} to a v5 \code{Assay5}, restoring a sidecar
#'
#' @keywords internal
#' @noRd
.srt_assay_to_v5 <- function(assay, verbose = TRUE) {
  if (inherits(assay, "Assay5")) return(assay)
  if (!inherits(assay, "Assay")) {
    stop("Cannot upgrade an assay of class ", class(assay)[1], call. = FALSE)
  }
  sidecar <- .srt_get_sidecar(assay)
  cls <- class(assay)[1]
  extra_slots <- setdiff(slotNames(assay), slotNames("Assay"))

  # Start from SeuratObject's own coercion (counts/data/scale.data,
  # feature-level meta data, variable features, key, misc)
  base <- assay
  if (cls != "Assay") {
    # Subclasses (SCTAssay, ChromatinAssay, ...) coerce through plain Assay;
    # their extra slots are kept in the sidecar so a downgrade restores them.
    base <- methods::as(assay, "Assay", strict = TRUE)
  }
  to <- suppressWarnings(methods::as(base, "Assay5"))

  if (!is.null(sidecar) && identical(sidecar$kind, "assay5")) {
    to <- .srt_restore_assay5_from_sidecar(to, assay, sidecar, verbose = verbose)
  } else {
    if (cls != "Assay" || length(extra_slots)) {
      sc <- list(kind = "assay", class = cls, package = attr(class(assay), "package"),
                 extra_slots = lapply(setNames(extra_slots, extra_slots),
                                      function(nm) slot(assay, nm)))
      to <- .srt_set_sidecar(to, sc)
      if (verbose) {
        message("Assay class ", cls, " has slots Assay5 cannot hold (",
                paste(extra_slots, collapse = ", "), "); kept in the sidecar")
      }
    }
  }
  methods::validObject(to)
  to
}

.srt_restore_assay5_from_sidecar <- function(to, from_v3, sidecar, verbose = TRUE) {
  # Original layer set, in order, with per-layer cells and features
  layers <- sidecar$layers
  joined <- list(
    counts = tryCatch(SeuratObject::GetAssayData(from_v3, layer = "counts"), error = function(e) NULL),
    data = tryCatch(SeuratObject::GetAssayData(from_v3, layer = "data"), error = function(e) NULL),
    scale.data = tryCatch(SeuratObject::GetAssayData(from_v3, layer = "scale.data"), error = function(e) NULL)
  )
  # If `data` was synthesised at downgrade time and has not been changed
  # since, it must not survive the upgrade
  if (isTRUE(sidecar$synthesized_data) &&
      .srt_matrix_identical(joined$data, joined$counts)) {
    joined$data <- NULL
  }
  # Rebuild every original layer from the joined matrices or the sidecar
  rebuilt <- list()
  for (nm in names(layers)) {
    rec <- layers[[nm]]
    mat <- if (identical(rec$source, "extra")) {
      sidecar$extra_layers[[nm]]
    } else {
      src <- joined[[rec$source]]
      if (is.null(src)) next
      feats <- intersect(rec$features, rownames(src))
      cells <- intersect(rec$cells, colnames(src))
      src[feats, cells, drop = FALSE]
    }
    if (is.null(mat)) next
    rebuilt[[nm]] <- mat
  }
  # Replace the coerced layers by the original ones (SeuratObject warns
  # about the default layer moving while layers are swapped; the default is
  # set explicitly below)
  suppressWarnings({
    for (nm in names(rebuilt)) {
      SeuratObject::LayerData(to, layer = nm) <- rebuilt[[nm]]
    }
    for (nm in setdiff(SeuratObject::Layers(to), names(rebuilt))) {
      SeuratObject::LayerData(to, layer = nm) <- NULL
    }
  })
  # Restore the original layer order (the layers slot and the cell/feature
  # LogMaps are all keyed by layer name)
  ord <- intersect(names(layers), SeuratObject::Layers(to))
  if (length(ord) == length(SeuratObject::Layers(to)) && !identical(ord, SeuratObject::Layers(to))) {
    slot(to, "layers") <- slot(to, "layers")[ord]
    slot(to, "cells") <- slot(to, "cells")[, ord]
    slot(to, "features") <- slot(to, "features")[, ord]
  }
  dl <- intersect(sidecar$default_layer, SeuratObject::Layers(to))
  if (length(dl)) {
    suppressWarnings(SeuratObject::DefaultLayer(to) <- dl)
  }
  if (!is.null(sidecar$assay_orig)) slot(to, "assay.orig") <- sidecar$assay_orig
  if (!is.null(sidecar$misc)) slot(to, "misc") <- sidecar$misc
  to <- .srt_set_sidecar(to, NULL)
  if (isTRUE(verbose) && length(sidecar$on_disk)) {
    message("Layers ", paste(sidecar$on_disk, collapse = ", "),
            " were on-disk (", paste(unique(sidecar$on_disk_class), collapse = ", "),
            ") before the downgrade and are now in memory")
  }
  to
}

# ---- assay: v5 Assay5 -> v3/v4 Assay --------------------------------------------

#' Convert a v5 \code{Assay5} to a v3/v4 \code{Assay}, recording a sidecar
#'
#' @keywords internal
#' @noRd
.srt_assay_to_v3 <- function(assay, verbose = TRUE) {
  if (inherits(assay, "Assay")) {
    # Already v3-style; restore a subclass recorded by a previous upgrade
    return(assay)
  }
  if (!inherits(assay, "Assay5")) {
    stop("Cannot downgrade an assay of class ", class(assay)[1], call. = FALSE)
  }
  prior <- .srt_get_sidecar(assay)
  layers <- SeuratObject::Layers(assay)
  sidecar <- list(kind = "assay5", layers = list(), extra_layers = list(),
                  default_layer = tryCatch(SeuratObject::DefaultLayer(assay), error = function(e) NULL),
                  assay_orig = slot(assay, "assay.orig"),
                  on_disk = character(0), on_disk_class = character(0),
                  synthesized_data = FALSE)
  misc <- slot(assay, "misc")
  if (is.list(misc)) misc[[.SRT_SIDECAR]] <- NULL
  sidecar$misc <- misc

  # Classify every layer
  for (nm in layers) {
    base <- .srt_layer_base(nm)
    mat <- SeuratObject::LayerData(assay, layer = nm)
    if (.srt_is_ondisk(mat)) {
      sidecar$on_disk <- c(sidecar$on_disk, nm)
      sidecar$on_disk_class <- c(sidecar$on_disk_class, class(mat)[1])
    }
    rec <- list(
      cells = SeuratObject::Cells(assay, layer = nm),
      features = SeuratObject::Features(assay, layer = nm),
      source = if (is.na(base)) "extra" else base,
      on_disk = .srt_is_ondisk(mat)
    )
    sidecar$layers[[nm]] <- rec
    if (is.na(base)) {
      if (.srt_is_ondisk(mat)) mat <- methods::as(mat, "dgCMatrix")
      sidecar$extra_layers[[nm]] <- mat
    }
  }
  # SeuratObject's coercion joins split layers and drops the rest
  to <- suppressWarnings(methods::as(assay, "Assay"))
  # `data` gets a copy of counts when the v5 assay had no data layer
  has_data <- any(vapply(names(sidecar$layers), function(nm) {
    identical(sidecar$layers[[nm]]$source, "data")
  }, logical(1)))
  if (!has_data) sidecar$synthesized_data <- TRUE

  # Restore a subclass recorded by a previous upgrade (SCTAssay, ...)
  if (!is.null(prior) && identical(prior$kind, "assay") && !identical(prior$class, "Assay")) {
    restored <- tryCatch({
      cls <- prior$class
      if (!is.null(prior$package) && !methods::isClass(cls)) {
        requireNamespace(prior$package, quietly = TRUE)
      }
      obj <- methods::as(to, cls)
      for (nm in names(prior$extra_slots)) {
        slot(obj, nm) <- prior$extra_slots[[nm]]
      }
      obj
    }, error = function(e) {
      warning("Could not restore assay class ", prior$class, ": ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
    if (!is.null(restored)) to <- restored
  }
  # misc: the assay's own misc plus the sidecar
  slot(to, "misc") <- misc
  dropped <- names(sidecar$extra_layers)
  split_layers <- setdiff(names(sidecar$layers)[!names(sidecar$layers) %in% .SRT_STANDARD_LAYERS], dropped)
  default_differs <- !identical(sidecar$default_layer,
                                .srt_expected_default(names(sidecar$layers)))
  if (length(dropped) || length(split_layers) || sidecar$synthesized_data ||
      length(sidecar$on_disk) || default_differs) {
    to <- .srt_set_sidecar(to, sidecar)
  }
  if (isTRUE(verbose)) {
    if (length(split_layers)) {
      message("Joined split layers ", paste(split_layers, collapse = ", "),
              " (cell membership kept in the sidecar)")
    }
    if (length(dropped)) {
      message("Layers ", paste(dropped, collapse = ", "),
              " have no v3 slot; kept in the sidecar")
    }
    if (length(sidecar$on_disk)) {
      message("On-disk layers ", paste(sidecar$on_disk, collapse = ", "), " loaded into memory")
    }
  }
  methods::validObject(to)
  to
}

# ---- images: VisiumV2 <-> VisiumV1 ----------------------------------------------

# spot.radius as Read10X_Image(image.type = "VisiumV1") computes it:
# spot diameter (full-res pixels) * lowres scale factor / image size
.srt_visium_spot_radius <- function(img) {
  sf <- slot(img, "scale.factors")
  dims <- dim(slot(img, "image"))
  r <- tryCatch({
    if (!is.null(sf[["spot"]]) && !is.null(sf[["lowres"]]) && length(dims) >= 2L && max(dims[1:2]) > 0) {
      (sf[["spot"]] * sf[["lowres"]]) / max(dims[1:2])
    } else NULL
  }, error = function(e) NULL)
  if (is.null(r) || !is.finite(r)) {
    r <- tryCatch(as.numeric(Seurat::Radius(img)), error = function(e) NULL)
  }
  if (is.null(r) || !length(r) || !is.finite(r)) r <- 0
  as.numeric(r)
}

.srt_visium_v2_to_v1 <- function(img, verbose = TRUE) {
  if (!inherits(img, "VisiumV2")) return(list(image = img, sidecar = NULL))
  coords <- Seurat::GetTissueCoordinates(img)
  cells <- as.character(coords$cell %||% rownames(coords))
  coordinates <- data.frame(
    tissue = 1L,
    row = NA_integer_,
    col = NA_integer_,
    imagerow = as.numeric(coords$y),
    imagecol = as.numeric(coords$x),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  # Space Ranger tissue_positions columns, when a previous conversion kept them
  prior <- if (is.list(slot(img, "misc"))) slot(img, "misc")[["tissue_positions"]] else NULL
  if (is.data.frame(prior) && all(c("row", "col") %in% colnames(prior))) {
    common <- intersect(rownames(prior), cells)
    coordinates[common, "row"] <- prior[common, "row"]
    coordinates[common, "col"] <- prior[common, "col"]
    if ("tissue" %in% colnames(prior)) coordinates[common, "tissue"] <- prior[common, "tissue"]
  }
  v1 <- methods::new(
    Class = "VisiumV1",
    image = slot(img, "image"),
    scale.factors = slot(img, "scale.factors"),
    coordinates = coordinates,
    spot.radius = .srt_visium_spot_radius(img),
    assay = SeuratObject::DefaultAssay(img),
    key = SeuratObject::Key(img)
  )
  if ("misc" %in% slotNames(v1)) {
    m <- slot(img, "misc"); if (!is.list(m)) m <- list()
    slot(v1, "misc") <- m
  }
  # Everything VisiumV1 cannot hold: boundaries (with radius / segmentations),
  # molecules, orientation. The image array itself stays in the V1 object.
  stub <- img
  slot(stub, "image") <- array(numeric(0), dim = c(0L, 0L, 0L))
  sidecar <- list(kind = "VisiumV2", object = stub)
  if (verbose) message("Converted VisiumV2 image to VisiumV1 (boundaries kept in the sidecar)")
  list(image = v1, sidecar = sidecar)
}

.srt_visium_v1_to_v2 <- function(img, sidecar = NULL, verbose = TRUE) {
  if (!inherits(img, "VisiumV1")) return(img)
  if (!is.null(sidecar) && identical(sidecar$kind, "VisiumV2")) {
    v2 <- sidecar$object
    slot(v2, "image") <- slot(img, "image")
    slot(v2, "scale.factors") <- slot(img, "scale.factors")
    if (verbose) message("Restored VisiumV2 image from the sidecar")
    return(v2)
  }
  # No sidecar: build a VisiumV2 the way Read10X_Image() does
  coordinates <- slot(img, "coordinates")
  if (!all(c("imagecol", "imagerow") %in% colnames(coordinates))) return(img)
  sf <- slot(img, "scale.factors")
  fov <- SeuratObject::CreateFOV(
    coords = coordinates[, c("imagecol", "imagerow"), drop = FALSE],
    type = "centroids",
    radius = sf[["spot"]] %||% slot(img, "spot.radius"),
    assay = SeuratObject::DefaultAssay(img),
    key = SeuratObject::Key(img)
  )
  v2 <- methods::new(
    Class = "VisiumV2",
    boundaries = slot(fov, "boundaries"),
    molecules = slot(fov, "molecules"),
    assay = SeuratObject::DefaultAssay(img),
    key = SeuratObject::Key(img),
    image = slot(img, "image"),
    scale.factors = sf,
    coords_x_orientation = "horizontal"
  )
  if ("misc" %in% slotNames(v2)) {
    m <- if ("misc" %in% slotNames(img)) slot(img, "misc") else list()
    if (!is.list(m)) m <- list()
    # keep the Space Ranger row/col/tissue columns so a later downgrade can
    # restore them
    m[["tissue_positions"]] <- coordinates
    slot(v2, "misc") <- m
  }
  if (verbose) message("Converted VisiumV1 image to VisiumV2")
  v2
}

# ---- object level --------------------------------------------------------------

#' Describe which Seurat generation an object belongs to
#'
#' @keywords internal
#' @noRd
.srt_object_version <- function(object) {
  stopifnot(inherits(object, "Seurat"))
  assays <- SeuratObject::Assays(object)
  classes <- vapply(assays, function(a) class(object[[a]])[1], character(1))
  images <- tryCatch(SeuratObject::Images(object), error = function(e) character(0))
  img_classes <- vapply(images, function(i) class(object[[i]])[1], character(1))
  has_v5 <- any(vapply(assays, function(a) inherits(object[[a]], "Assay5"), logical(1)))
  on_disk <- character(0)
  for (a in assays[classes == "Assay5"]) {
    for (l in SeuratObject::Layers(object[[a]])) {
      if (.srt_is_ondisk(SeuratObject::LayerData(object[[a]], layer = l))) on_disk <- c(on_disk, paste0(a, "/", l))
    }
  }
  generation <- if (has_v5 || "VisiumV2" %in% img_classes) "v5" else "v4"
  structure(list(
    version = as.character(slot(object, "version")),
    generation = generation,
    assays = classes,
    images = img_classes,
    has_assay5 = has_v5,
    has_visium_v2 = "VisiumV2" %in% img_classes,
    on_disk_layers = on_disk,
    has_sidecar = !is.null(.srt_get_sidecar(object)) ||
      any(vapply(assays, function(a) !is.null(.srt_get_sidecar(object[[a]])), logical(1)))
  ), class = "seurat_generation")
}

#' @export
#' @noRd
print.seurat_generation <- function(x, ...) {
  cat("Seurat object generation:", x$generation, "(version slot", x$version, ")\n")
  for (a in names(x$assays)) cat(sprintf("  assay %-14s %s\n", a, x$assays[[a]]))
  for (i in names(x$images)) cat(sprintf("  image %-14s %s\n", i, x$images[[i]]))
  if (length(x$on_disk_layers)) cat("  on-disk layers:", paste(x$on_disk_layers, collapse = ", "), "\n")
  if (x$has_sidecar) cat("  carries a version-conversion sidecar\n")
  invisible(x)
}

#' Convert a Seurat object between the v3/v4 and v5 generations
#'
#' @keywords internal
#' @noRd
.srt_convert_object <- function(object, to = c("v5", "v4", "v3"), verbose = TRUE) {
  to <- match.arg(to)
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object", call. = FALSE)
  target_v5 <- identical(to, "v5")
  obj_sidecar <- .srt_get_sidecar(object) %||% list(kind = "seurat", images = list())
  if (is.null(obj_sidecar$images)) obj_sidecar$images <- list()

  # Assays
  for (a in SeuratObject::Assays(object)) {
    assay <- object[[a]]
    converted <- if (target_v5) {
      .srt_assay_to_v5(assay, verbose = verbose)
    } else {
      .srt_assay_to_v3(assay, verbose = verbose)
    }
    if (!identical(class(converted), class(assay)) || !identical(converted, assay)) {
      # Replace the assay in place without triggering cell-name validation
      # on partial assays
      slot(object, "assays")[[a]] <- converted
    }
  }

  # Images
  images <- tryCatch(SeuratObject::Images(object), error = function(e) character(0))
  for (i in images) {
    img <- object[[i]]
    if (target_v5) {
      restored <- .srt_visium_v1_to_v2(img, sidecar = obj_sidecar$images[[i]], verbose = verbose)
      if (!identical(restored, img)) {
        slot(object, "images")[[i]] <- restored
        obj_sidecar$images[[i]] <- NULL
      }
    } else {
      res <- .srt_visium_v2_to_v1(img, verbose = verbose)
      if (!is.null(res$sidecar)) {
        slot(object, "images")[[i]] <- res$image
        obj_sidecar$images[[i]] <- res$sidecar
      }
    }
  }

  # Version stamp
  if (target_v5) {
    if (!is.null(obj_sidecar$original_version)) {
      slot(object, "version") <- package_version(obj_sidecar$original_version)
      obj_sidecar$original_version <- NULL
    } else {
      slot(object, "version") <- utils::packageVersion("SeuratObject")
    }
  } else {
    obj_sidecar$original_version <- as.character(slot(object, "version"))
    slot(object, "version") <- package_version(if (identical(to, "v3")) "3.2.3" else "4.4.0")
  }

  keep_sidecar <- length(obj_sidecar$images) > 0L || !is.null(obj_sidecar$original_version)
  object <- .srt_set_sidecar(object, if (keep_sidecar) obj_sidecar else NULL)
  methods::validObject(object)
  object
}

#' Convert a Seurat object held in a file
#'
#' @keywords internal
#' @noRd
.srt_convert_file <- function(source, dest, to, overwrite = FALSE, verbose = TRUE,
                              reader, writer) {
  if (!file.exists(source)) stop("File not found: ", source, call. = FALSE)
  if (file.exists(dest) && !overwrite) stop("Destination exists: ", dest, call. = FALSE)
  ext <- tolower(tools::file_ext(source))
  object <- if (ext == "rds") readRDS(source) else reader(source)
  if (!inherits(object, "Seurat")) {
    stop("File does not contain a Seurat object: ", source, call. = FALSE)
  }
  object <- .srt_convert_object(object, to = to, verbose = verbose)
  dext <- tolower(tools::file_ext(dest))
  if (dext == "rds") {
    saveRDS(object, dest)
  } else {
    if (file.exists(dest)) file.remove(dest)
    writer(object, dest)
  }
  invisible(dest)
}
