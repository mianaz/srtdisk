#' @include AnnDataCompat.R
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# On-disk h5ad layout upgrade / downgrade
#
# anndata changed its h5ad layout in 0.8.0 (March 2022) and again, more
# subtly, in 0.11 (nullable strings) and 0.13 (pandas 3 string inference
# makes nullable-string-array the default for every string column and the
# obs/var index, and `None` is written with the `null` encoding). Files
# written by a newer anndata cannot be opened by an older one, and the
# original SeuratDisk converters only understand the pre-0.8 layout.
#
# The rewriter below walks an h5ad file element by element and re-encodes it
# for a target layout:
#
#   "encoded"  the current anndata layout (>= 0.8): every element carries
#              encoding-type / encoding-version; categoricals are groups;
#              nullable pandas dtypes are values + mask groups.
#   "legacy"   the anndata 0.7 layout: dataframe 0.1.0 with `__categories`
#              sibling datasets referenced by an HDF5 object reference, no
#              encoding attributes on plain arrays, dict groups and scalars.
#
# Anything the target layout cannot represent (a nullable integer with
# missing values, a `None` in uns, ...) is written in the closest
# representable form *and recorded* in a manifest inside `uns`, so that
# converting back restores the original encodings exactly. Matrices are
# copied at the HDF5 level and never decoded.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.H5AD_MANIFEST <- "__h5ad_compat_manifest__"
.H5AD_MATRIX_ENCODINGS <- c("csr_matrix", "csc_matrix")
.H5AD_OPAQUE_ENCODINGS <- c("awkward-array", "rec-array")

# ---- low level hdf5r writing helpers ----------------------------------------

.h5w_string_type <- function() {
  hdf5r::H5T_STRING$new(size = Inf)$set_cset("UTF-8")
}

.h5w_bool_type <- function() {
  # h5py reads a FALSE/TRUE enum over int8 back as numpy bool
  hdf5r::H5T_ENUM$new(labels = c("FALSE", "TRUE"), values = c(0L, 1L))
}

.h5w_scalar_space <- function() hdf5r::H5S$new(type = "scalar")

.h5w_attr <- function(obj, name, value, force_array = FALSE) {
  if (isTRUE(obj$attr_exists(attr_name = name))) obj$attr_delete(attr_name = name)
  scalar <- length(value) == 1L && !force_array
  if (force_array && is.character(value)) {
    # anndata iterates `column-order`; a length-1 (or empty) value must still
    # be a 1-D array attribute, not a scalar string
    obj$create_attr(attr_name = name, robj = value, dtype = .h5w_string_type(),
                    space = hdf5r::H5S$new(type = "simple", dims = length(value),
                                           maxdims = length(value)))
    return(invisible(NULL))
  }
  if (is.character(value)) {
    if (scalar) {
      obj$create_attr(attr_name = name, robj = value, dtype = .h5w_string_type(),
                      space = .h5w_scalar_space())
    } else {
      obj$create_attr(attr_name = name, robj = value, dtype = .h5w_string_type())
    }
  } else if (is.logical(value)) {
    if (scalar) {
      obj$create_attr(attr_name = name, robj = value, dtype = .h5w_bool_type(),
                      space = .h5w_scalar_space())
    } else {
      obj$create_attr(attr_name = name, robj = value, dtype = .h5w_bool_type())
    }
  } else if (scalar) {
    obj$create_attr(attr_name = name, robj = value, space = .h5w_scalar_space())
  } else {
    obj$create_attr(attr_name = name, robj = value)
  }
  invisible(NULL)
}

.h5w_drop_attr <- function(obj, name) {
  if (isTRUE(tryCatch(obj$attr_exists(attr_name = name), error = function(e) FALSE))) {
    obj$attr_delete(attr_name = name)
  }
  invisible(NULL)
}

.h5w_encoding <- function(obj, type, version) {
  .h5w_attr(obj, "encoding-type", type)
  .h5w_attr(obj, "encoding-version", version)
}

.h5w_strip_encoding <- function(obj) {
  .h5w_drop_attr(obj, "encoding-type")
  .h5w_drop_attr(obj, "encoding-version")
}

#' Write an R vector as a 1-D (or scalar) HDF5 dataset
#'
#' Strings become variable-length UTF-8, logicals an h5py-compatible bool
#' enum (no NA allowed), integers/doubles their native types unless
#' \code{dtype} is given.
#'
#' @keywords internal
#' @noRd
.h5w_dataset <- function(grp, name, x, dtype = NULL, scalar = FALSE) {
  if (grp$exists(name)) grp$link_delete(name)
  if (is.null(dtype)) {
    dtype <- if (is.character(x)) {
      .h5w_string_type()
    } else if (is.logical(x)) {
      .h5w_bool_type()
    } else if (is.integer(x)) {
      hdf5r::h5types$H5T_STD_I64LE
    } else {
      NULL
    }
  }
  if (is.factor(x)) x <- as.character(x)
  if (scalar) {
    return(invisible(grp$create_dataset(
      name = name, robj = x, dtype = dtype,
      space = .h5w_scalar_space(), chunk_dims = NULL
    )))
  }
  n <- length(x)
  if (n == 0L) {
    return(invisible(grp$create_dataset(name = name, robj = x, dtype = dtype,
                                        chunk_dims = NULL)))
  }
  gz <- tryCatch(GetCompressionLevel(), error = function(e) 4L)
  if (n < 1024L || is.character(x)) gz <- 0L
  invisible(grp$create_dataset(
    name = name, robj = x, dtype = dtype,
    chunk_dims = if (gz > 0L) min(n, 65536L) else NULL,
    gzip_level = gz
  ))
}

.h5w_codes_dtype <- function(n_levels) {
  if (n_levels <= 127L) hdf5r::h5types$H5T_STD_I8LE
  else if (n_levels <= 32767L) hdf5r::h5types$H5T_STD_I16LE
  else hdf5r::h5types$H5T_STD_I32LE
}

.h5w_copy <- function(src_parent, name, dst_parent, dst_name = name) {
  if (dst_parent$exists(dst_name)) dst_parent$link_delete(dst_name)
  dst_parent$obj_copy_from(src_loc = src_parent, src_name = name, dst_name = dst_name)
  dst_parent[[dst_name]]
}

# Object reference to `obj` (an H5D / H5Group of `file`) for a
# H5T_STD_REF_OBJ attribute. hdf5r's own `create_reference()` opens a
# temporary file handle whose garbage-collection finalizer corrupts hdf5r's
# reference bookkeeping ("r_count can never be more than 1 larger than
# h5_count" after enough references, at GC-dependent moments). A classic
# object reference is the object's header address, so it is assembled here
# from `obj_info()` and attached to the file handle we already hold; the
# caller closes the returned object right after writing the attribute.
.h5w_object_reference <- function(file, obj) {
  addr <- as.numeric(obj$obj_info()$addr)
  bytes <- integer(8)
  for (i in seq_len(8)) {
    bytes[i] <- addr %% 256
    addr <- addr %/% 256
  }
  ref <- hdf5r::H5R_OBJECT$new(1, file)
  ref$ref <- as.raw(bytes)
  ref
}

.h5_is_scalar <- function(dset) {
  isTRUE(tryCatch(dset$get_space()$get_simple_extent_type() == "H5S_SCALAR",
                  error = function(e) FALSE)) ||
    isTRUE(tryCatch(length(dset$dims) == 0L, error = function(e) FALSE))
}

.h5_is_null_space <- function(dset) {
  isTRUE(tryCatch(dset$get_space()$get_simple_extent_type() == "H5S_NULL",
                  error = function(e) FALSE))
}

.h5_is_compound <- function(dset) {
  isTRUE(tryCatch(dset$get_type()$get_class() == "H5T_COMPOUND",
                  error = function(e) FALSE))
}

.h5_is_string <- function(dset) {
  isTRUE(tryCatch(dset$get_type()$get_class() == "H5T_STRING",
                  error = function(e) FALSE))
}

# ---- layout inspection ------------------------------------------------------

#' Inspect the on-disk layout of an h5ad file
#'
#' Reports which anndata layout a file uses and which encodings it contains,
#' which is what decides whether an older or newer anndata (or a converter
#' written against one of them) can read it.
#'
#' @param file Path to an h5ad file
#'
#' @return A list of class \code{h5ad_layout} with elements \code{layout}
#'   (\code{"encoded"}, \code{"legacy"} or \code{"compound"}),
#'   \code{min_anndata} (the oldest anndata release able to read the file),
#'   \code{encodings} (a table of the encoding types found), and logical
#'   flags \code{has_nullable}, \code{has_nullable_strings},
#'   \code{has_null}, \code{has_manifest}.
#'
#' @keywords internal
#' @noRd
.h5ad_layout_info <- function(file) {
  h5 <- hdf5r::H5File$new(file, mode = "r")
  on.exit(tryCatch(h5$close_all(), error = function(e) NULL), add = TRUE)
  encs <- character(0)
  legacy_cats <- FALSE
  compound <- FALSE
  walk <- function(grp, depth) {
    for (nm in names(grp)) {
      obj <- tryCatch(grp[[nm]], error = function(e) NULL)
      if (is.null(obj)) next
      enc <- .h5ad_encoding(obj)
      if (nzchar(enc)) encs <<- c(encs, enc)
      if (identical(nm, "__categories")) legacy_cats <<- TRUE
      if (inherits(obj, "H5D") && .h5_is_compound(obj) && depth == 0L &&
          nm %in% c("obs", "var")) compound <<- TRUE
      if (inherits(obj, "H5Group") && depth < 3L && !identical(nm, "__categories") &&
          !enc %in% .H5AD_MATRIX_ENCODINGS) {
        walk(obj, depth + 1L)
      }
      .h5_close_quietly(obj)
    }
  }
  walk(h5, 0L)
  root_enc <- .h5ad_encoding(h5)
  obs_is_group <- h5$exists("obs") && .h5_with_child(h5, "obs", function(o) inherits(o, "H5Group"))
  obs_df_version <- if (obs_is_group) {
    as.character(.h5_with_child(h5, "obs", function(o) .h5ad_attr(o, "encoding-version", "")))
  } else ""
  layout <- if (compound) {
    "compound"
  } else if (legacy_cats || identical(obs_df_version, "0.1.0") ||
             (obs_is_group && !nzchar(root_enc) && !"categorical" %in% encs)) {
    "legacy"
  } else {
    "encoded"
  }
  has_manifest <- h5$exists("uns") && .h5_with_child(h5, "uns", function(u) u$exists(.H5AD_MANIFEST))
  has_ns <- "nullable-string-array" %in% encs
  has_nb <- any(c("nullable-integer", "nullable-boolean") %in% encs)
  has_null <- "null" %in% encs
  min_anndata <- if (layout == "compound") {
    "0.6"
  } else if (layout == "legacy") {
    "0.7"
  } else if (has_ns) {
    "0.11"
  } else if (has_null) {
    "0.9"
  } else {
    "0.8"
  }
  structure(list(
    file = file,
    layout = layout,
    min_anndata = min_anndata,
    root_encoding = root_enc,
    encodings = table(encs),
    has_nullable = has_nb,
    has_nullable_strings = has_ns,
    has_null = has_null,
    has_manifest = has_manifest,
    needs_string_materialization = has_ns || has_null
  ), class = "h5ad_layout")
}

#' @export
#' @noRd
print.h5ad_layout <- function(x, ...) {
  cat("h5ad layout:", x$layout, "\n")
  cat("  file:                ", x$file, "\n")
  cat("  readable by anndata >=", x$min_anndata, "\n")
  cat("  nullable columns:    ", x$has_nullable, "\n")
  cat("  nullable strings:    ", x$has_nullable_strings, "\n")
  cat("  null (None) entries: ", x$has_null, "\n")
  cat("  compat manifest:     ", x$has_manifest, "\n")
  if (length(x$encodings)) {
    cat("  encodings:\n")
    for (nm in names(x$encodings)) cat(sprintf("    %-22s %d\n", nm, x$encodings[[nm]]))
  }
  invisible(x)
}

# ---- manifest ----------------------------------------------------------------

.h5ad_manifest_key <- function(path) gsub("/", "|", path, fixed = TRUE)
.h5ad_manifest_path <- function(key) gsub("|", "/", key, fixed = TRUE)

.h5ad_manifest_read <- function(h5) {
  out <- list(elements = list(), source_layout = NULL)
  if (!(h5$exists("uns") && h5[["uns"]]$exists(.H5AD_MANIFEST))) return(out)
  mg <- h5[["uns"]][[.H5AD_MANIFEST]]
  if (mg$exists("source_layout")) {
    out$source_layout <- as.character(mg[["source_layout"]]$read())
  }
  if (mg$exists("elements")) {
    eg <- mg[["elements"]]
    for (key in names(eg)) {
      rec <- strsplit(as.character(eg[[key]]$read()), "\t", fixed = TRUE)[[1]]
      out$elements[[.h5ad_manifest_path(key)]] <- list(
        type = rec[1], version = if (length(rec) > 1L) rec[2] else "",
        note = if (length(rec) > 2L) rec[3] else ""
      )
    }
  }
  out
}

.h5ad_manifest_write <- function(h5, ctx) {
  if (!length(ctx$manifest)) return(invisible(NULL))
  if (!h5$exists("uns")) {
    ug <- h5$create_group("uns")
    if (ctx$target == "encoded") .h5w_encoding(ug, "dict", "0.1.0")
  }
  ug <- h5[["uns"]]
  if (ug$exists(.H5AD_MANIFEST)) ug$link_delete(.H5AD_MANIFEST)
  mg <- ug$create_group(.H5AD_MANIFEST)
  eg <- mg$create_group("elements")
  if (ctx$target == "encoded") {
    .h5w_encoding(mg, "dict", "0.1.0")
    .h5w_encoding(eg, "dict", "0.1.0")
  }
  wr <- function(name, value) {
    d <- .h5w_dataset(mg, name, value, scalar = TRUE)
    if (ctx$target == "encoded") .h5w_encoding(d, "string", "0.2.0")
  }
  wr("writer", ctx$writer)
  wr("writer_version", ctx$writer_version)
  wr("source_layout", ctx$source_layout)
  wr("target_layout", ctx$target)
  for (path in names(ctx$manifest)) {
    rec <- ctx$manifest[[path]]
    d <- .h5w_dataset(eg, .h5ad_manifest_key(path),
                      paste(rec$type, rec$version, rec$note %||% "", sep = "\t"),
                      scalar = TRUE)
    if (ctx$target == "encoded") .h5w_encoding(d, "string", "0.2.0")
  }
  invisible(NULL)
}

.ctx_record <- function(ctx, path, type, version, note = "") {
  ctx$manifest[[path]] <- list(type = type, version = version, note = note)
  invisible(NULL)
}

.ctx_restore <- function(ctx, path) {
  ctx$restore[[path]]
}

# ---- column classification ---------------------------------------------------

# What a decoded R column would naturally be written as in the encoded layout
.h5ad_natural_encoding <- function(x) {
  if (is.factor(x)) return("categorical")
  if (is.logical(x)) return(if (anyNA(x)) "nullable-boolean" else "array")
  if (is.character(x)) return(if (anyNA(x)) "nullable-string-array" else "string-array")
  if (is.integer(x) || inherits(x, "integer64")) return(if (anyNA(x)) "nullable-integer" else "array")
  "array"
}

# Encoding actually present in the source for a column (before decoding)
.h5ad_source_column_encoding <- function(grp, col) {
  if (!grp$exists(col)) return("")
  obj <- grp[[col]]
  enc <- .h5ad_encoding(obj)
  if (nzchar(enc)) return(enc)
  if (inherits(obj, "H5Group")) {
    if (obj$exists("codes") && obj$exists("categories")) return("categorical")
    if (obj$exists("values") && obj$exists("mask")) return("nullable")
    return("")
  }
  if (grp$exists("__categories") && grp[["__categories"]]$exists(col)) {
    return("legacy-categorical")
  }
  if (.h5_is_string(obj)) return("string-array")
  "array"
}

# ---- dataframe writers -------------------------------------------------------

.h5ad_write_categorical_encoded <- function(parent, name, x) {
  g <- parent$create_group(name)
  codes <- as.integer(x) - 1L
  codes[is.na(codes)] <- -1L
  lv <- levels(x)
  d <- .h5w_dataset(g, "codes", codes, dtype = .h5w_codes_dtype(length(lv)))
  .h5w_encoding(d, "array", "0.2.0")
  d <- .h5w_dataset(g, "categories", lv)
  .h5w_encoding(d, "string-array", "0.2.0")
  .h5w_attr(g, "ordered", is.ordered(x))
  .h5w_encoding(g, "categorical", "0.2.0")
  invisible(g)
}

.h5ad_write_nullable_encoded <- function(parent, name, x, type, na_value = "NaN") {
  g <- parent$create_group(name)
  mask <- is.na(x)
  values <- x
  if (type == "nullable-boolean") {
    values[mask] <- FALSE
    d <- .h5w_dataset(g, "values", as.logical(values))
  } else if (type == "nullable-string-array") {
    values[mask] <- ""
    d <- .h5w_dataset(g, "values", as.character(values))
    .h5w_encoding(d, "string-array", "0.2.0")
    # anndata >= 0.11 stores the pandas missing-value flavour: "NA" (pd.NA,
    # dtype "string") or "NaN" (np.nan, the pandas 3 "str" dtype)
    .h5w_attr(g, "na-value", if (identical(na_value, "NA")) "NA" else "NaN")
  } else {
    values[mask] <- 0L
    d <- .h5w_dataset(g, "values", as.integer(values), dtype = hdf5r::h5types$H5T_STD_I64LE)
  }
  if (type != "nullable-string-array") .h5w_encoding(d, "array", "0.2.0")
  d <- .h5w_dataset(g, "mask", as.logical(mask))
  .h5w_encoding(d, "array", "0.2.0")
  .h5w_encoding(g, type, "0.1.0")
  invisible(g)
}

.h5ad_write_array_encoded <- function(parent, name, x, force_string = FALSE) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x) && anyNA(x)) x[is.na(x)] <- FALSE
  if (is.character(x) && anyNA(x)) x[is.na(x)] <- ""
  dtype <- if (is.integer(x)) hdf5r::h5types$H5T_STD_I64LE else NULL
  d <- .h5w_dataset(parent, name, x, dtype = dtype)
  .h5w_encoding(d, if (is.character(x) || force_string) "string-array" else "array", "0.2.0")
  invisible(d)
}

# Write one decoded column into an encoded (>= 0.8) dataframe group.
# `want` is the encoding to produce (from the manifest when upgrading,
# otherwise the natural one); `strings` controls nullable strings.
.h5ad_write_column_encoded <- function(g, col, x, want, ctx, na_value = "NaN") {
  if (want == "categorical" || (is.factor(x) && want != "string-array" && want != "nullable-string-array")) {
    if (!is.factor(x)) x <- factor(x)
    return(.h5ad_write_categorical_encoded(g, col, x))
  }
  if (want %in% c("nullable-integer", "nullable-boolean")) {
    return(.h5ad_write_nullable_encoded(g, col, x, want))
  }
  if (want == "nullable-string-array") {
    if (identical(ctx$strings, "categorical")) {
      if (anyNA(x)) {
        return(.h5ad_write_categorical_encoded(g, col, factor(x)))
      }
      return(.h5ad_write_array_encoded(g, col, as.character(x), force_string = TRUE))
    }
    return(.h5ad_write_nullable_encoded(g, col, as.character(x), want,
                                        na_value = na_value))
  }
  if (is.factor(x)) x <- as.character(x)
  if ((is.logical(x) || is.integer(x) || is.character(x)) && anyNA(x)) {
    # target says plain array but values have NA: fall back to nullable
    return(.h5ad_write_column_encoded(g, col, x, .h5ad_natural_encoding(x), ctx))
  }
  .h5ad_write_array_encoded(g, col, x)
}

# Write one decoded column into a legacy (0.7) dataframe group. Anything the
# layout cannot hold is written as its closest representable form and
# recorded in the manifest so an upgrade can restore it.
.h5ad_write_column_legacy <- function(g, col, x, src_enc, path, ctx, na_note = "") {
  write_cat <- function(f, note = "") {
    if (!g$exists("__categories")) .h5_close_quietly(g$create_group("__categories"))
    cg <- g[["__categories"]]
    on.exit(.h5_close_quietly(cg), add = TRUE)
    lv <- levels(f)
    cd <- .h5w_dataset(cg, col, lv)
    on.exit(.h5_close_quietly(cd), add = TRUE)
    .h5w_attr(cd, "ordered", is.ordered(f))
    codes <- as.integer(f) - 1L
    codes[is.na(codes)] <- -1L
    d <- .h5w_dataset(g, col, codes, dtype = .h5w_codes_dtype(length(lv)))
    on.exit(.h5_close_quietly(d), add = TRUE)
    ref <- .h5w_object_reference(ctx$dst_file, cd)
    on.exit(.h5_close_quietly(ref), add = TRUE)
    d$create_attr(attr_name = "categories", robj = ref,
                  dtype = hdf5r::h5types$H5T_STD_REF_OBJ, space = .h5w_scalar_space())
    invisible(NULL)
  }
  if (is.factor(x)) {
    if (src_enc %in% c("nullable-string-array", "string-array")) {
      .ctx_record(ctx, path, src_enc, if (src_enc == "string-array") "0.2.0" else "0.1.0",
                  note = na_note)
    }
    return(write_cat(x))
  }
  if (is.logical(x)) {
    if (anyNA(x)) {
      .ctx_record(ctx, path, "nullable-boolean", "0.1.0")
      return(write_cat(factor(ifelse(x, "True", "False"), levels = c("False", "True"))))
    }
    return(.h5w_dataset(g, col, x))
  }
  if (is.character(x)) {
    if (anyNA(x)) {
      .ctx_record(ctx, path, if (nzchar(src_enc)) src_enc else "nullable-string-array", "0.1.0",
                  note = na_note)
      return(write_cat(factor(x)))
    }
    if (src_enc == "nullable-string-array") .ctx_record(ctx, path, src_enc, "0.1.0", note = na_note)
    return(.h5w_dataset(g, col, x))
  }
  if (is.integer(x) || inherits(x, "integer64")) {
    if (anyNA(x)) {
      .ctx_record(ctx, path, "nullable-integer", "0.1.0")
      return(.h5w_dataset(g, col, as.double(x)))
    }
    if (src_enc == "nullable-integer") .ctx_record(ctx, path, src_enc, "0.1.0")
    return(.h5w_dataset(g, col, as.integer(x), dtype = hdf5r::h5types$H5T_STD_I64LE))
  }
  .h5w_dataset(g, col, as.double(x))
}

#' Rewrite an AnnData dataframe (obs, var, raw/var, or a dataframe in uns)
#'
#' @keywords internal
#' @noRd
.h5ad_write_dataframe <- function(src, name, dst_parent, path, ctx) {
  obj <- src[[name]]
  on.exit(.h5_close_quietly(obj), add = TRUE)
  child_attr <- function(cn, attr, default) {
    .h5_with_child(obj, cn, function(o) as.character(.h5ad_attr(o, attr, default)))
  }
  # ---- decode
  if (inherits(obj, "H5D")) {
    # anndata < 0.7: single compound dataset. Categories (if any) live in
    # uns/<col>_categories.
    df <- obj$read()
    idx_name <- if ("index" %in% names(df)) "index" else if ("_index" %in% names(df)) "_index" else NULL
    index <- if (!is.null(idx_name)) as.character(df[[idx_name]]) else rownames(df)
    cols <- setdiff(names(df), idx_name)
    columns <- lapply(cols, function(cn) {
      v <- df[[cn]]
      if (is.numeric(v) && ctx$src_uns_exists(paste0(cn, "_categories"))) {
        cats <- as.character(ctx$src_uns_read(paste0(cn, "_categories")))
        v <- .h5ad_decode_codes(as.integer(v), cats)
      }
      v
    })
    names(columns) <- cols
    src_encs <- setNames(rep("", length(cols)), cols)
    idx_enc <- "string-array"
    index_name <- "_index"
  } else {
    index_name <- .h5ad_index_name(obj) %||% "_index"
    index <- .h5ad_read_index(obj)
    idx_enc <- if (obj$exists(index_name)) .h5ad_source_column_encoding(obj, index_name) else "string-array"
    cols <- .h5ad_dataframe_columns(obj)
    columns <- list()
    src_encs <- character(0)
    for (cn in cols) {
      v <- tryCatch(.h5ad_read_column(obj, cn), error = function(e) NULL)
      if (is.null(v)) {
        # undecodable column (e.g. awkward): copy verbatim later
        columns[cn] <- list(NULL)
      } else {
        columns[[cn]] <- v
      }
      src_encs[cn] <- .h5ad_source_column_encoding(obj, cn)
    }
  }
  n <- length(index)
  if (is.null(index)) {
    n <- if (length(columns)) length(columns[[1]]) else 0L
    index <- as.character(seq_len(n) - 1L)
  }
  # ---- write
  if (dst_parent$exists(name)) dst_parent$link_delete(name)
  g <- dst_parent$create_group(name)
  keep <- character(0)
  for (cn in names(columns)) {
    v <- columns[[cn]]
    cpath <- paste0(path, "/", cn)
    if (is.null(v)) {
      if (inherits(obj, "H5Group") && obj$exists(cn)) {
        .h5w_copy(obj, cn, g)
        keep <- c(keep, cn)
      }
      next
    }
    if (length(v) != n) next
    keep <- c(keep, cn)
    if (ctx$target == "encoded") {
      rec <- .ctx_restore(ctx, cpath)
      want <- rec$type %||% NULL
      if (is.null(want)) {
        want <- switch(src_encs[[cn]],
          "nullable-integer" = , "nullable-boolean" = , "nullable-string-array" = src_encs[[cn]],
          "categorical" = , "legacy-categorical" = "categorical",
          .h5ad_natural_encoding(v))
      }
      na_value <- rec$note %||% ""
      if (!nzchar(na_value) && inherits(obj, "H5Group") && obj$exists(cn)) {
        na_value <- child_attr(cn, "na-value", "NaN")
      }
      .h5ad_write_column_encoded(g, cn, v, want, ctx, na_value = na_value)
    } else {
      na_note <- if (inherits(obj, "H5Group") && obj$exists(cn)) {
        child_attr(cn, "na-value", "")
      } else ""
      .h5ad_write_column_legacy(g, cn, v, src_encs[[cn]], cpath, ctx, na_note = na_note)
    }
  }
  # index
  idx <- as.character(index)
  ipath <- paste0(path, "/", index_name)
  if (ctx$target == "encoded") {
    irec <- .ctx_restore(ctx, ipath)
    want <- irec$type %||% idx_enc
    if (want == "nullable-string-array" && !identical(ctx$strings, "categorical")) {
      na_value <- irec$note %||% ""
      if (!nzchar(na_value) && inherits(obj, "H5Group") && obj$exists(index_name)) {
        na_value <- child_attr(index_name, "na-value", "NaN")
      }
      .h5ad_write_nullable_encoded(g, index_name, idx, "nullable-string-array",
                                   na_value = na_value)
    } else {
      if (anyNA(idx)) idx[is.na(idx)] <- ""
      .h5ad_write_array_encoded(g, index_name, idx, force_string = TRUE)
    }
    .h5w_encoding(g, "dataframe", "0.2.0")
  } else {
    if (idx_enc == "nullable-string-array") {
      .ctx_record(ctx, ipath, idx_enc, "0.1.0",
                  note = child_attr(index_name, "na-value", ""))
    }
    if (anyNA(idx)) idx[is.na(idx)] <- ""
    .h5w_dataset(g, index_name, idx)
    .h5w_encoding(g, "dataframe", "0.1.0")
  }
  .h5w_attr(g, "_index", index_name)
  .h5w_attr(g, "column-order", keep, force_array = TRUE)
  .h5_close_quietly(g)
  invisible(NULL)
}

# ---- generic element rewriter -----------------------------------------------

.h5ad_rewrite_dataset <- function(src, name, dst, path, ctx) {
  obj <- src[[name]]
  on.exit(.h5_close_quietly(obj), add = TRUE)
  enc <- .h5ad_encoding(obj)
  if (identical(enc, "null") || .h5_is_null_space(obj)) {
    # anndata < 0.9 has no `null` encoding: drop it for the legacy layout and
    # for the compatibility ("categorical" strings) flavour of the encoded one
    if (ctx$target == "legacy" || isTRUE(ctx$drop_null)) {
      .ctx_record(ctx, path, "null", "0.1.0")
      return(invisible(NULL))
    }
    d <- .h5w_copy(src, name, dst)
    .h5w_encoding(d, "null", "0.1.0")
    .h5_close_quietly(d)
    return(invisible(NULL))
  }
  d <- .h5w_copy(src, name, dst)
  on.exit(.h5_close_quietly(d), add = TRUE)
  if (ctx$target == "encoded") {
    if (!nzchar(enc) || enc == "array" || enc == "string-array" ||
        enc == "numeric-scalar" || enc == "string" || enc == "bytes") {
      if (.h5_is_scalar(d)) {
        .h5w_encoding(d, if (.h5_is_string(d)) "string" else "numeric-scalar", "0.2.0")
      } else if (.h5_is_compound(d)) {
        .h5w_encoding(d, "rec-array", "0.2.0")
      } else {
        .h5w_encoding(d, if (.h5_is_string(d)) "string-array" else "array", "0.2.0")
      }
    }
  } else {
    .h5w_strip_encoding(d)
  }
  invisible(NULL)
}

.h5ad_rewrite_group <- function(src, name, dst, path, ctx, role = "generic") {
  obj <- src[[name]]
  on.exit(.h5_close_quietly(obj), add = TRUE)
  enc <- .h5ad_encoding(obj)
  if (identical(name, .H5AD_MANIFEST) && role == "uns") {
    return(invisible(NULL))  # never carry a stale manifest forward
  }
  # Sparse matrices and opaque encodings: copy at the HDF5 level
  if (enc %in% .H5AD_MATRIX_ENCODINGS ||
      (obj$exists("data") && obj$exists("indices") && obj$exists("indptr"))) {
    g <- .h5w_copy(src, name, dst)
    on.exit(.h5_close_quietly(g), add = TRUE)
    if (!nzchar(enc)) enc <- "csr_matrix"
    .h5w_encoding(g, enc, "0.1.0")
    if (!isTRUE(g$attr_exists("shape")) && ctx$target == "encoded") {
      # shape is required by anndata >= 0.8; derive it from indptr/indices
      indptr <- .h5_with_child(obj, "indptr", function(d) d$read())
      indices <- .h5_with_child(obj, "indices", function(d) d$read())
      major <- length(indptr) - 1L
      minor <- if (length(indices)) max(indices) + 1L else 0L
      .h5w_attr(g, "shape", if (enc == "csc_matrix") c(minor, major) else c(major, minor))
    }
    return(invisible(NULL))
  }
  if (enc %in% .H5AD_OPAQUE_ENCODINGS) {
    .h5_close_quietly(.h5w_copy(src, name, dst))
    return(invisible(NULL))
  }
  # Dataframes
  is_df <- identical(enc, "dataframe") || role %in% c("obs", "var") ||
    !is.null(.h5ad_index_name(obj)) && (obj$exists("__categories") || nzchar(as.character(.h5ad_attr(obj, "_index", ""))))
  if (is_df) {
    return(invisible(.h5ad_write_dataframe(src, name, dst, path, ctx)))
  }
  # Standalone categorical / nullable groups (uns, obsm)
  if (identical(enc, "categorical") || (obj$exists("codes") && obj$exists("categories"))) {
    v <- .h5ad_read_categorical_group(obj)
    if (ctx$target == "encoded") {
      return(invisible(.h5ad_write_categorical_encoded(dst, name, v)))
    }
    .ctx_record(ctx, path, "categorical", "0.2.0", note = if (is.ordered(v)) "ordered" else "")
    return(invisible(.h5w_dataset(dst, name, as.character(v))))
  }
  if (obj$exists("values") && obj$exists("mask")) {
    v <- .h5ad_read_nullable(obj)
    type <- if (nzchar(enc)) enc else .h5ad_natural_encoding(v)
    na_value <- .ctx_restore(ctx, path)$note %||% as.character(.h5ad_attr(obj, "na-value", "NaN"))
    if (ctx$target == "encoded") {
      return(invisible(.h5ad_write_column_encoded(dst, name, v, type, ctx, na_value = na_value)))
    }
    .ctx_record(ctx, path, type, "0.1.0", note = na_value)
    if (is.character(v) || is.logical(v)) {
      f <- if (is.logical(v)) factor(ifelse(v, "True", "False"), levels = c("False", "True")) else factor(v)
      codes <- as.integer(f) - 1L; codes[is.na(codes)] <- -1L
      g <- dst$create_group(name)
      .h5w_dataset(g, "codes", codes, dtype = .h5w_codes_dtype(nlevels(f)))
      .h5w_dataset(g, "categories", levels(f))
      return(invisible(g))
    }
    return(invisible(.h5w_dataset(dst, name, as.double(v))))
  }
  # Plain (dict) group: recurse
  g <- dst$create_group(name)
  on.exit(.h5_close_quietly(g), add = TRUE)
  child_role <- if (role == "raw") "raw-child" else if (role == "uns") "uns" else "generic"
  for (child in names(obj)) {
    cpath <- paste0(path, "/", child)
    child_is_dataset <- .h5_with_child(obj, child, function(o) inherits(o, "H5D"))
    r <- child_role
    if (role == "raw" && child == "var") r <- "var"
    if (child_is_dataset) {
      .h5ad_rewrite_dataset(obj, child, g, cpath, ctx)
    } else {
      .h5ad_rewrite_group(obj, child, g, cpath, ctx, role = r)
    }
  }
  if (ctx$target == "encoded") {
    .h5w_encoding(g, if (role == "raw") "raw" else "dict", "0.1.0")
  } else {
    .h5w_strip_encoding(g)
  }
  invisible(NULL)
}

#' Rewrite an h5ad file into another on-disk layout
#'
#' @param source Path to the input h5ad file
#' @param dest Path to the output h5ad file
#' @param target \code{"encoded"} (anndata >= 0.8 layout) or
#'   \code{"legacy"} (anndata 0.7 layout)
#' @param overwrite Overwrite \code{dest} if it exists
#' @param strings For \code{target = "encoded"}: \code{"keep"} writes
#'   nullable-string-array elements as such; \code{"categorical"} rewrites
#'   them as categoricals (with missing values) or plain string arrays
#'   (without), which every anndata >= 0.8 and every SeuratDisk-derived
#'   converter can read.
#' @param manifest Record every lossy re-encoding in
#'   \code{uns/__h5ad_compat_manifest__} so the reverse conversion can restore
#'   the original encodings
#' @param verbose Print progress
#'
#' @return \code{dest}, invisibly
#'
#' @keywords internal
#' @noRd
.h5ad_rewrite <- function(source, dest, target = c("encoded", "legacy"),
                          overwrite = FALSE, strings = c("keep", "categorical"),
                          manifest = TRUE, verbose = TRUE) {
  target <- match.arg(target)
  strings <- match.arg(strings)
  if (!file.exists(source)) stop("File not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite) stop("Destination exists: ", dest, call. = FALSE)
    file.remove(dest)
  }
  info <- .h5ad_layout_info(source)
  src <- hdf5r::H5File$new(source, mode = "r")
  on.exit(tryCatch(src$close_all(), error = function(e) NULL), add = TRUE)
  dst <- hdf5r::H5File$new(dest, mode = "w")
  on.exit(tryCatch(dst$close_all(), error = function(e) NULL), add = TRUE)

  prior <- .h5ad_manifest_read(src)
  pkg <- utils::packageName() %||% "scConvert"
  ctx <- new.env(parent = emptyenv())
  ctx$target <- target
  ctx$strings <- strings
  ctx$drop_null <- identical(strings, "categorical")
  ctx$manifest <- list()
  ctx$restore <- if (manifest) prior$elements else list()
  ctx$source_layout <- info$layout
  ctx$writer <- pkg
  ctx$writer_version <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "")
  ctx$dst_file <- dst
  ctx$src_uns_exists <- function(key) src$exists("uns") && src[["uns"]]$exists(key)
  ctx$src_uns_read <- function(key) src[["uns"]][[key]]$read()

  if (verbose) {
    message("Rewriting ", basename(source), " (", info$layout, " layout) as ",
            target, " layout -> ", basename(dest))
  }
  roles <- c(obs = "obs", var = "var", raw = "raw", uns = "uns")
  for (nm in names(src)) {
    role <- roles[nm] %||% "generic"
    if (is.na(role)) role <- "generic"
    if (verbose) message("  ", nm)
    if (.h5_with_child(src, nm, function(o) inherits(o, "H5D"))) {
      if (nm %in% c("obs", "var")) {
        .h5ad_write_dataframe(src, nm, dst, nm, ctx)
      } else {
        .h5ad_rewrite_dataset(src, nm, dst, nm, ctx)
      }
    } else {
      .h5ad_rewrite_group(src, nm, dst, nm, ctx, role = role)
    }
  }
  # Restore elements the previous downgrade had to drop (null values)
  if (target == "encoded" && length(ctx$restore)) {
    for (path in names(ctx$restore)) {
      rec <- ctx$restore[[path]]
      if (identical(rec$type, "null")) {
        parts <- strsplit(path, "/", fixed = TRUE)[[1]]
        parent <- dst
        ok <- TRUE
        for (p in parts[-length(parts)]) {
          if (!parent$exists(p)) { ok <- FALSE; break }
          parent <- parent[[p]]
        }
        if (ok && !parent$exists(parts[length(parts)])) {
          d <- parent$create_dataset(parts[length(parts)], dtype = hdf5r::h5types$H5T_IEEE_F32LE,
                                     space = hdf5r::H5S$new(type = "null"), chunk_dims = NULL)
          .h5w_encoding(d, "null", "0.1.0")
        }
      }
    }
  }
  if (target == "encoded") {
    .h5w_encoding(dst, "anndata", "0.1.0")
  } else {
    .h5w_strip_encoding(dst)
  }
  if (manifest) .h5ad_manifest_write(dst, ctx)
  if (verbose) {
    n_rec <- length(ctx$manifest)
    message("Done. ", if (n_rec) paste0(n_rec, " element(s) re-encoded and recorded in uns/",
                                        .H5AD_MANIFEST) else "No lossy re-encodings were needed.")
  }
  invisible(dest)
}

#' Materialise an h5ad file so that SeuratDisk-derived converters can read it
#'
#' Returns \code{file} unchanged when nothing in it needs rewriting;
#' otherwise writes a temporary copy with nullable strings turned into
#' categoricals / string arrays and \code{null} entries dropped, and returns
#' that path. Callers should delete the temporary file when done.
#'
#' @keywords internal
#' @noRd
.h5ad_normalize_for_convert <- function(file, verbose = FALSE) {
  info <- tryCatch(.h5ad_layout_info(file), error = function(e) NULL)
  if (is.null(info) || !isTRUE(info$needs_string_materialization)) {
    return(list(path = file, temporary = FALSE))
  }
  tmp <- tempfile(fileext = ".h5ad")
  if (verbose) {
    message("Input uses nullable-string / null encodings (anndata >= 0.11); ",
            "materialising a converter-compatible copy")
  }
  .h5ad_rewrite(file, tmp, target = "encoded", strings = "categorical",
                manifest = FALSE, verbose = FALSE)
  list(path = tmp, temporary = TRUE)
}
