#' @include zzz.R
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AnnData on-disk encoding compatibility helpers (hdf5r backend)
#
# One place that knows how to read every element encoding that anndata has
# written to h5ad files from 0.6 through 0.13, so the readers in this package
# do not each carry their own partial copy of the rules:
#
#   anndata < 0.8  ("legacy")   dataframe 0.1.0: obs/var groups with a
#                               `_index` attribute, categoricals stored as
#                               integer-code datasets whose categories live in
#                               a sibling `__categories/<col>` dataset (the
#                               `ordered` flag sits on that dataset). Even
#                               older files store obs/var as a single compound
#                               dataset. No `encoding-type` attributes on
#                               columns.
#   anndata >= 0.8 ("encoded")  every element carries `encoding-type` /
#                               `encoding-version`; categoricals are groups
#                               (codes + categories, `ordered` attribute);
#                               pandas nullable dtypes are groups holding
#                               `values` + `mask` (nullable-integer 0.1.0,
#                               nullable-boolean 0.1.0). anndata 0.11 added
#                               nullable-string-array 0.1.0 (same layout, with
#                               an `na-value` attribute), which anndata 0.13 +
#                               pandas 3 write for *every* string column and
#                               index by default, and the `null` encoding
#                               (empty dataspace) for `None` values in uns.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read an HDF5 attribute, returning a default when it is absent
#'
#' @param obj hdf5r object (H5File, H5Group or H5D)
#' @param name Attribute name
#' @param default Value returned when the attribute does not exist or cannot
#'   be read
#'
#' @return The attribute value, or \code{default}
#'
#' @keywords internal
#' @noRd
.h5ad_attr <- function(obj, name, default = NULL) {
  if (is.null(obj)) return(default)
  ok <- tryCatch(isTRUE(obj$attr_exists(attr_name = name)), error = function(e) FALSE)
  if (!ok) return(default)
  tryCatch(hdf5r::h5attr(x = obj, which = name), error = function(e) default)
}

#' AnnData encoding type of an element (\code{""} when absent)
#'
#' @keywords internal
#' @noRd
.h5ad_encoding <- function(obj) {
  enc <- .h5ad_attr(obj, "encoding-type", default = "")
  if (length(enc) != 1L || is.na(enc)) return("")
  as.character(enc)
}

#' Coerce what hdf5r hands back for an HDF5 boolean to an R logical
#'
#' anndata writes booleans as an HDF5 enum (\code{FALSE}/\code{TRUE}); hdf5r
#' usually decodes that to a logical, but older builds return a factor, and
#' legacy files store booleans as uint8.
#'
#' @keywords internal
#' @noRd
#' Close an hdf5r object handle without failing
#'
#' Handles are closed as soon as they are no longer needed: hdf5r's
#' \code{close_all()} and its GC finalizers misbehave on HDF5 1.12 (the
#' library bundled by the CRAN Windows binary) when a file is closed while
#' many child handles are still open.
#'
#' @keywords internal
#' @noRd
.h5_close_quietly <- function(obj) {
  if (!is.null(obj) && inherits(obj, "H5RefClass") && !inherits(obj, "H5File")) {
    tryCatch(obj$close(), error = function(e) NULL)
  }
  invisible(NULL)
}

#' Open a child of an HDF5 group, apply a function and close the handle
#'
#' @keywords internal
#' @noRd
.h5_with_child <- function(parent, name, fn) {
  child <- parent[[name]]
  on.exit(.h5_close_quietly(child), add = TRUE)
  fn(child)
}

.h5ad_as_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.factor(x)) return(as.logical(toupper(as.character(x))))
  if (is.character(x)) return(as.logical(toupper(x)))
  as.logical(x)
}

#' Decode a nullable (\code{values} + \code{mask}) group
#'
#' Covers \code{nullable-integer}, \code{nullable-boolean} and
#' \code{nullable-string-array}. In every case \code{mask == TRUE} marks a
#' missing value. Groups without encoding attributes are decoded the same way
#' as long as they have both datasets.
#'
#' @param grp H5Group
#'
#' @return integer, logical or character vector with \code{NA} where masked,
#'   or \code{NULL} if the group is not a nullable layout
#'
#' @keywords internal
#' @noRd
.h5ad_read_nullable <- function(grp) {
  if (!inherits(grp, "H5Group")) return(NULL)
  if (!(grp$exists("values") && grp$exists("mask"))) return(NULL)
  values <- .h5_with_child(grp, "values", function(d) d$read())
  mask <- .h5ad_as_logical(.h5_with_child(grp, "mask", function(d) d$read()))
  enc <- .h5ad_encoding(grp)
  if (identical(enc, "nullable-boolean")) {
    values <- .h5ad_as_logical(values)
  } else if (identical(enc, "nullable-string-array")) {
    values <- as.character(values)
  } else if (identical(enc, "nullable-integer")) {
    # hdf5r returns int64 as integer64/double; keep R integer where it fits
    if (!is.integer(values)) {
      dv <- suppressWarnings(as.double(values))
      values <- if (all(is.na(dv) | abs(dv) < .Machine$integer.max)) as.integer(dv) else dv
    }
  } else if (is.factor(values) && all(levels(values) %in% c("FALSE", "TRUE"))) {
    values <- .h5ad_as_logical(values)
  }
  if (length(mask) == length(values)) {
    values[mask] <- NA
  }
  values
}

#' Decode a categorical group (codes + categories + ordered)
#'
#' @param grp H5Group with \code{codes} and \code{categories}
#'
#' @return factor (ordered when the AnnData \code{ordered} flag is set), or
#'   \code{NULL} if the group is not a categorical layout
#'
#' @keywords internal
#' @noRd
.h5ad_read_categorical_group <- function(grp) {
  if (!inherits(grp, "H5Group")) return(NULL)
  if (!(grp$exists("codes") && grp$exists("categories"))) return(NULL)
  codes <- .h5_with_child(grp, "codes", function(d) d$read())
  categories <- .h5ad_read_string_like(grp, "categories")
  if (is.null(categories)) {
    categories <- as.character(.h5_with_child(grp, "categories", function(d) d$read()))
  }
  is_ordered <- isTRUE(.h5ad_as_logical(.h5ad_attr(grp, "ordered", FALSE))[1])
  .h5ad_decode_codes(as.integer(codes), as.character(categories), ordered = is_ordered)
}

#' Turn 0-based category codes into an R factor
#'
#' Level order follows the stored category order verbatim; \code{-1} codes
#' (and any out-of-range code) become \code{NA}.
#'
#' @keywords internal
#' @noRd
.h5ad_decode_codes <- function(codes, categories, ordered = FALSE) {
  codes[is.na(codes) | codes < 0L | codes >= length(categories)] <- NA_integer_
  factor(categories[codes + 1L], levels = categories, ordered = isTRUE(ordered))
}

#' Read a string-like child: a plain dataset or a nullable-string-array group
#'
#' @param parent H5Group
#' @param name Child name
#'
#' @return character vector (\code{NA} for masked entries), or \code{NULL}
#'   when the child does not exist or is not string-like
#'
#' @keywords internal
#' @noRd
.h5ad_read_string_like <- function(parent, name) {
  if (is.null(parent) || !parent$exists(name)) return(NULL)
  child <- parent[[name]]
  on.exit(.h5_close_quietly(child), add = TRUE)
  if (inherits(child, "H5D")) {
    return(as.character(child$read()))
  }
  if (inherits(child, "H5Group")) {
    vals <- .h5ad_read_nullable(child)
    if (!is.null(vals)) return(as.character(vals))
  }
  NULL
}

#' Name of the index column of an AnnData dataframe group
#'
#' The \code{_index} attribute names it (anndata >= 0.7); older files use a
#' literal \code{_index} or \code{index} child.
#'
#' @keywords internal
#' @noRd
.h5ad_index_name <- function(grp) {
  idx <- .h5ad_attr(grp, "_index", default = NULL)
  if (!is.null(idx) && length(idx) == 1L && nzchar(idx) && grp$exists(idx)) {
    return(as.character(idx))
  }
  for (cand in c("_index", "index")) {
    if (grp$exists(cand)) return(cand)
  }
  NULL
}

#' Read the index (row names) of an AnnData dataframe group
#'
#' @param grp H5Group (obs, var, raw/var, or a dataframe in uns)
#'
#' @return character vector, or \code{NULL} when no index can be found
#'
#' @keywords internal
#' @noRd
.h5ad_read_index <- function(grp) {
  if (is.null(grp) || !inherits(grp, "H5Group")) return(NULL)
  idx <- .h5ad_index_name(grp)
  if (is.null(idx)) return(NULL)
  .h5ad_read_string_like(grp, idx)
}

#' Column names of an AnnData dataframe group, in stored order
#'
#' Uses the \code{column-order} attribute when present (that is the pandas
#' column order), otherwise the HDF5 link order; the index column and the
#' legacy \code{__categories} group are never columns.
#'
#' @keywords internal
#' @noRd
.h5ad_dataframe_columns <- function(grp) {
  if (is.null(grp) || !inherits(grp, "H5Group")) return(character(0))
  present <- names(grp)
  drop <- c("__categories", "_index", "index", .h5ad_index_name(grp))
  cols <- .h5ad_attr(grp, "column-order", default = NULL)
  if (!is.null(cols) && length(cols) > 0L) {
    cols <- as.character(cols)
    cols <- c(cols[cols %in% present], setdiff(present, cols))
  } else {
    cols <- present
  }
  setdiff(cols, drop)
}

#' Read one column of an AnnData dataframe group in any encoding
#'
#' Understands, in this order: nullable groups (\code{values}/\code{mask}),
#' categorical groups (codes/categories/ordered), legacy code datasets whose
#' categories live in \code{__categories/<col>} (with that dataset's
#' \code{ordered} attribute), HDF5 boolean enums, and plain datasets.
#'
#' @param grp H5Group (obs / var)
#' @param col Column name
#'
#' @return An R vector, or \code{NULL} when the column cannot be decoded
#'
#' @keywords internal
#' @noRd
.h5ad_read_column <- function(grp, col) {
  if (!grp$exists(col)) return(NULL)
  obj <- grp[[col]]
  on.exit(.h5_close_quietly(obj), add = TRUE)
  if (inherits(obj, "H5Group")) {
    vals <- .h5ad_read_nullable(obj)
    if (!is.null(vals)) return(vals)
    vals <- .h5ad_read_categorical_group(obj)
    return(vals)
  }
  if (!inherits(obj, "H5D")) return(NULL)
  # Legacy (< 0.8) categorical: integer codes + __categories/<col>
  if (grp$exists("__categories")) {
    cats_grp <- grp[["__categories"]]
    on.exit(.h5_close_quietly(cats_grp), add = TRUE)
    if (cats_grp$exists(col)) {
      codes <- obj$read()
      cats_dset <- cats_grp[[col]]
      on.exit(.h5_close_quietly(cats_dset), add = TRUE)
      categories <- as.character(cats_dset$read())
      # anndata 0.7 stored booleans as categoricals with "False"/"True" categories
      if (length(categories) <= 2L && all(categories %in% c("False", "True"))) {
        codes <- as.integer(codes)
        codes[codes < 0L] <- NA_integer_
        return(categories[codes + 1L] == "True")
      }
      is_ordered <- isTRUE(.h5ad_as_logical(.h5ad_attr(cats_dset, "ordered", FALSE))[1])
      return(.h5ad_decode_codes(as.integer(codes), categories, ordered = is_ordered))
    }
  }
  vals <- obj$read()
  if (is.factor(vals) && all(levels(vals) %in% c("FALSE", "TRUE"))) {
    vals <- .h5ad_as_logical(vals)
  }
  vals
}

#' Read a whole AnnData dataframe group into a data.frame
#'
#' @param grp H5Group with the \code{dataframe} layout (any version)
#' @param row.names Optional row names to use instead of the stored index
#'
#' @return data.frame
#'
#' @keywords internal
#' @noRd
.h5ad_read_dataframe <- function(grp, row.names = NULL) {
  idx <- row.names %||% .h5ad_read_index(grp)
  cols <- .h5ad_dataframe_columns(grp)
  out <- list()
  for (col in cols) {
    v <- tryCatch(.h5ad_read_column(grp, col), error = function(e) NULL)
    if (is.null(v)) next
    if (!is.null(idx) && length(v) != length(idx)) next
    out[[col]] <- v
  }
  n <- if (!is.null(idx)) length(idx) else if (length(out)) length(out[[1]]) else 0L
  df <- if (length(out)) as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE,
                                       check.names = FALSE)
        else data.frame(matrix(nrow = n, ncol = 0))
  if (!is.null(idx) && length(idx) == nrow(df)) {
    rownames(df) <- make.unique(as.character(idx))
  }
  df
}

#' Read a scalar / array element from uns, tolerating every encoding
#'
#' \code{null} (empty dataspace) becomes \code{NULL}; nullable groups,
#' categoricals and dataframes decode through the helpers above; other
#' groups become named lists.
#'
#' @keywords internal
#' @noRd
.h5ad_read_element <- function(obj) {
  if (inherits(obj, "H5D")) {
    enc <- .h5ad_encoding(obj)
    if (identical(enc, "null")) return(NULL)
    simple <- tryCatch({
      sp <- obj$get_space()
      on.exit(.h5_close_quietly(sp), add = TRUE)
      sp$is_simple()
    }, error = function(e) TRUE)
    if (!isTRUE(simple)) return(NULL)
    vals <- obj$read()
    if (is.factor(vals) && all(levels(vals) %in% c("FALSE", "TRUE"))) {
      vals <- .h5ad_as_logical(vals)
    }
    return(vals)
  }
  if (!inherits(obj, "H5Group")) return(NULL)
  enc <- .h5ad_encoding(obj)
  if (identical(enc, "dataframe") || !is.null(.h5ad_index_name(obj))) {
    return(.h5ad_read_dataframe(obj))
  }
  nn <- .h5ad_read_nullable(obj)
  if (!is.null(nn)) return(nn)
  cg <- .h5ad_read_categorical_group(obj)
  if (!is.null(cg)) return(cg)
  out <- list()
  for (nm in names(obj)) {
    out[nm] <- list(tryCatch(.h5_with_child(obj, nm, .h5ad_read_element), error = function(e) NULL))
  }
  out
}

#' Seurat's feature-name normalisation, applied to a character vector
#'
#' \code{CreateSeuratObject()} replaces underscores with dashes in feature
#' names. Anything that later has to match those names (variable features,
#' scaled features, raw/X rows, layer rows) must go through the same
#' transformation, otherwise joins silently fail.
#'
#' @keywords internal
#' @noRd
.seurat_feature_names <- function(x) {
  x <- as.character(x)
  gsub("_", "-", x, fixed = TRUE)
}

#' Fill in columns the compiled reader does not decode
#'
#' The C reader (\code{C_read_h5ad}) handles plain datasets and categorical
#' groups. Nullable groups (\code{nullable-integer}, \code{nullable-boolean},
#' \code{nullable-string-array}) come back as \code{NULL} entries or are
#' missing; decode those with hdf5r so the result matches the R reader.
#'
#' @param cols Named list returned by the C reader (index already removed)
#' @param grp H5Group for obs / var
#' @param n Expected column length
#'
#' @keywords internal
#' @noRd
.h5ad_patch_c_columns <- function(cols, grp, n) {
  if (is.null(grp) || !inherits(grp, "H5Group")) return(cols)
  cols <- cols[!vapply(cols, is.null, logical(1))]
  wanted <- .h5ad_dataframe_columns(grp)
  # Legacy (< 0.8) categoricals are integer-code datasets; the C reader hands
  # them back as bare integers, so re-decode every column that has an entry
  # under __categories (this also restores the `ordered` flag).
  redo <- character(0)
  if (grp$exists("__categories")) {
    redo <- intersect(.h5_with_child(grp, "__categories", names), wanted)
    cols[redo] <- NULL
  }
  for (col in setdiff(wanted, names(cols))) {
    v <- tryCatch(.h5ad_read_column(grp, col), error = function(e) NULL)
    if (!is.null(v) && length(v) == n) cols[[col]] <- v
  }
  # keep pandas column order
  cols[c(intersect(wanted, names(cols)), setdiff(names(cols), wanted))]
}
