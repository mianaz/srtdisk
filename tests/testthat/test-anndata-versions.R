# Cross-version AnnData compatibility: fixtures written by anndata 0.7.8
# (legacy layout: __categories, no encoding attributes) and by anndata
# 0.13.3 with pandas 3 (nullable-string-array index and columns,
# nullable-integer / nullable-boolean, `null` in uns), plus the on-disk
# upgrade / downgrade round trip between the two layouts.
#
# The fixtures were generated with scripts in the scConvert development
# repository (dev/anndata-fixtures); srtdisk ships the same files; each holds a 20 x 8 AnnData with
# sparse X, a counts layer and a dense layer, categorical (ordered and
# unordered), boolean, integer and float obs columns, obsm/varm/obsp/varp and
# a nested uns (including a DataFrame).

skip_if_not_installed("hdf5r")
skip_if_not_installed("Seurat")
suppressPackageStartupMessages({library(Matrix); library(Seurat)})

legacy_file <- system.file("testdata", "anndata", "anndata_0.7.8.h5ad", package = "srtdisk")
modern_file <- system.file("testdata", "anndata", "anndata_0.13.3_pandas3.h5ad", package = "srtdisk")
zarr_dir <- system.file("testdata", "anndata", "anndata_0.13.3_pandas3.zarr", package = "srtdisk")

expect_fixture_object <- function(obj, modern = TRUE) {
  md <- obj@meta.data
  expect_equal(dim(obj), c(8L, 20L))
  expect_true(is.factor(md$cluster))
  expect_equal(levels(md$cluster), c("b", "a", "c"))   # stored order, not alphabetical
  expect_true(is.ordered(md$grade))
  expect_equal(levels(md$grade), c("low", "mid", "high"))
  expect_true(is.logical(md$is_doublet))
  expect_true(is.numeric(md$n_counts))
  if (modern) {
    expect_true("batch_id" %in% colnames(md))
    expect_equal(sum(is.na(md$batch_id)), 1L)
    expect_true(is.logical(md$qc_pass))
    expect_equal(sum(is.na(md$qc_pass)), 1L)
    expect_true("note" %in% colnames(md))
    expect_equal(sum(is.na(md$note)), 1L)
  }
  expect_true(all(c("pca", "umap") %in% names(obj@reductions)))
  expect_true(all(c("RNA_snn", "RNA_nn") %in% names(obj@graphs)))
  # underscores in the file's gene names become dashes
  expect_true("gene-3-with-underscore" %in% rownames(obj))
}

test_that("LoadH5AD reads the anndata 0.7 legacy layout", {
  skip_if(!file.exists(legacy_file))
  obj <- suppressWarnings(LoadH5AD(legacy_file, verbose = FALSE))
  expect_fixture_object(obj, modern = FALSE)
  expect_true(all(c("counts", "data") %in% Layers(obj[["RNA"]])))
})

test_that("LoadH5AD reads anndata 0.13 + pandas 3 encodings", {
  skip_if(!file.exists(modern_file))
  obj_r <- suppressWarnings(LoadH5AD(modern_file, verbose = FALSE))
  expect_fixture_object(obj_r)
  expect_null(obj_r@misc$none_value)
  expect_true(is.data.frame(obj_r@misc$df))
})

test_that("h5ad -> h5Seurat -> Seurat handles nullable strings and underscore genes", {
  skip_if(!file.exists(modern_file))
  tmp <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(tmp), add = TRUE)
  Convert(modern_file, dest = tmp, verbose = FALSE)
  obj <- suppressWarnings(LoadH5Seurat(tmp, verbose = FALSE))
  md <- obj@meta.data
  expect_equal(dim(obj), c(8L, 20L))
  expect_true(is.ordered(md$grade))
  expect_true("batch_id" %in% colnames(md))
  expect_equal(sum(is.na(md$batch_id)), 1L)
  expect_true("qc_pass" %in% colnames(md))
  expect_true("note" %in% colnames(md))
  expect_equal(sum(is.na(md$note)), 1L)
  expect_true("gene-3-with-underscore" %in% rownames(obj))
})

test_that("H5ADLayout reports layout and minimum anndata version", {
  skip_if(!file.exists(legacy_file) || !file.exists(modern_file))
  old <- H5ADLayout(legacy_file)
  expect_s3_class(old, "h5ad_layout")
  expect_equal(old$layout, "legacy")
  expect_equal(old$min_anndata, "0.7")
  new <- H5ADLayout(modern_file)
  expect_equal(new$layout, "encoded")
  expect_equal(new$min_anndata, "0.11")
  expect_true(new$has_nullable_strings)
  expect_true(new$has_null)
  expect_output(print(new), "nullable strings")
})

test_that("DowngradeH5AD / UpgradeH5AD round-trip the modern layout losslessly", {
  skip_if(!file.exists(modern_file))
  legacy <- tempfile(fileext = ".h5ad")
  restored <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(legacy, restored)), add = TRUE)
  DowngradeH5AD(modern_file, legacy, verbose = FALSE)
  info <- H5ADLayout(legacy)
  expect_equal(info$layout, "legacy")
  expect_true(info$has_manifest)
  expect_false(info$has_nullable_strings)
  expect_false(info$has_null)
  # the legacy file is readable and complete
  obj <- suppressWarnings(LoadH5AD(legacy, verbose = FALSE))
  expect_equal(dim(obj), c(8L, 20L))
  # (checked on the slot: Seurat's `$` accessor goes through FetchData, which
  # drops the ordered class)
  expect_true(is.ordered(obj@meta.data$grade))
  expect_true("qc_pass" %in% colnames(obj@meta.data))
  # and the upgrade restores every original encoding
  UpgradeH5AD(legacy, restored, verbose = FALSE)
  back <- H5ADLayout(restored)
  expect_equal(back$layout, "encoded")
  expect_false(back$has_manifest)
  expect_equal(as.list(back$encodings), as.list(H5ADLayout(modern_file)$encodings))
  h5a <- hdf5r::H5File$new(modern_file, "r"); h5b <- hdf5r::H5File$new(restored, "r")
  expect_equal(hdf5r::h5attr(h5b[["obs"]][["_index"]], "encoding-type"), "nullable-string-array")
  expect_equal(hdf5r::h5attr(h5b[["obs"]][["_index"]], "na-value"),
               hdf5r::h5attr(h5a[["obs"]][["_index"]], "na-value"))
  expect_equal(hdf5r::h5attr(h5b[["obs"]][["batch_id"]], "encoding-type"), "nullable-integer")
  expect_equal(h5b[["obs"]][["batch_id"]][["mask"]]$read(), h5a[["obs"]][["batch_id"]][["mask"]]$read())
  expect_true(h5b[["uns"]]$exists("none_value"))
  expect_equal(h5b[["X"]][["data"]]$read(), h5a[["X"]][["data"]]$read())
  h5a$close_all(); h5b$close_all()
})

test_that("UpgradeH5AD converts the legacy layout to the encoded layout", {
  skip_if(!file.exists(legacy_file))
  out <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out), add = TRUE)
  UpgradeH5AD(legacy_file, out, verbose = FALSE)
  info <- H5ADLayout(out)
  expect_equal(info$layout, "encoded")
  expect_equal(info$min_anndata, "0.8")
  h5 <- hdf5r::H5File$new(out, "r")
  expect_equal(hdf5r::h5attr(h5, "encoding-type"), "anndata")
  expect_equal(hdf5r::h5attr(h5[["obs"]], "encoding-version"), "0.2.0")
  expect_equal(hdf5r::h5attr(h5[["obs"]][["grade"]], "encoding-type"), "categorical")
  expect_true(isTRUE(as.logical(hdf5r::h5attr(h5[["obs"]][["grade"]], "ordered"))))
  expect_false(h5[["obs"]]$exists("__categories"))
  h5$close_all()
  obj <- suppressWarnings(LoadH5AD(out, verbose = FALSE))
  expect_fixture_object(obj, modern = FALSE)
})
