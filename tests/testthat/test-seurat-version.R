# Lossless upgrade / downgrade between Seurat v3/v4 (Assay) and v5 (Assay5)

skip_if_not_installed("Seurat")
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

make_v5_object <- function() {
  set.seed(1)
  counts <- as(matrix(rpois(2000, 2), 50, 40,
                      dimnames = list(paste0("g", 1:50), paste0("c", 1:40))), "dgCMatrix")
  withr::with_options(list(Seurat.object.assay.version = "v5"), {
    obj <- CreateSeuratObject(counts)
    obj$batch <- rep(c("a", "b"), each = 20)
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 20, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    LayerData(obj[["RNA"]], layer = "ambient") <- counts * 0.5
    Misc(obj[["RNA"]], "note") <- "hello"
    obj
  })
}

layers_identical <- function(a, b) {
  identical(Layers(a), Layers(b)) && all(vapply(Layers(a), function(l) {
    x <- LayerData(a, layer = l); y <- LayerData(b, layer = l)
    identical(dim(x), dim(y)) && identical(dimnames(x), dimnames(y)) &&
      isTRUE(all.equal(as.matrix(x), as.matrix(y)))
  }, logical(1)))
}

test_that("v3 Assay upgrades to Assay5 and back without a sidecar", {
  data("pbmc_small", package = "SeuratObject")
  v3 <- pbmc_small
  expect_equal(SeuratGeneration(v3)$generation, "v4")
  v5 <- UpgradeSeurat(v3, verbose = FALSE)
  expect_s4_class(v5[["RNA"]], "Assay5")
  expect_setequal(Layers(v5[["RNA"]]), c("counts", "data", "scale.data"))
  expect_identical(VariableFeatures(v5), VariableFeatures(v3))
  expect_equal(SeuratGeneration(v5)$generation, "v5")
  back <- DowngradeSeurat(v5, verbose = FALSE)
  expect_identical(class(back[["RNA"]])[1], "Assay")
  expect_equal(as.matrix(GetAssayData(back, layer = "counts")), as.matrix(GetAssayData(v3, layer = "counts")))
  expect_equal(GetAssayData(back, layer = "scale.data"), GetAssayData(v3, layer = "scale.data"))
  expect_null(back[["RNA"]]@misc[[".seurat_version_sidecar"]])
})

test_that("Assay5 with split and extra layers round-trips exactly", {
  skip_if_not_installed("withr")
  orig <- make_v5_object()
  expect_identical(Layers(orig[["RNA"]]), c("counts.a", "counts.b", "data.a", "data.b", "scale.data", "ambient"))
  down <- DowngradeSeurat(orig, verbose = FALSE)
  expect_identical(class(down[["RNA"]])[1], "Assay")
  expect_setequal(Layers(down[["RNA"]]), c("counts", "data", "scale.data"))
  expect_equal(as.character(down@version), "4.4.0")
  sidecar <- down[["RNA"]]@misc[[".seurat_version_sidecar"]]
  expect_false(is.null(sidecar))
  expect_named(sidecar$extra_layers, "ambient")
  # joined counts carry every cell
  counts <- GetAssayData(down, layer = "counts")
  expect_equal(sort(colnames(counts)), sort(colnames(orig)))
  up <- UpgradeSeurat(down, verbose = FALSE)
  expect_true(layers_identical(orig[["RNA"]], up[["RNA"]]))
  expect_identical(DefaultLayer(up[["RNA"]]), DefaultLayer(orig[["RNA"]]))
  expect_identical(VariableFeatures(up), VariableFeatures(orig))
  expect_identical(Misc(up[["RNA"]], "note"), "hello")
  expect_null(up[["RNA"]]@misc[[".seurat_version_sidecar"]])
  expect_equal(as.character(up@version), as.character(orig@version))
})

test_that("SCTAssay survives an upgrade / downgrade cycle", {
  set.seed(2)
  counts <- as(matrix(rpois(4000, 3), 80, 50,
                      dimnames = list(paste0("g", 1:80), paste0("c", 1:50))), "dgCMatrix")
  sct <- suppressWarnings(SCTransform(CreateSeuratObject(counts), verbose = FALSE))
  expect_s4_class(sct[["SCT"]], "SCTAssay")
  up <- UpgradeSeurat(sct, verbose = FALSE)
  expect_s4_class(up[["SCT"]], "Assay5")
  expect_false(is.null(up[["SCT"]]@misc[[".seurat_version_sidecar"]]))
  down <- DowngradeSeurat(up, verbose = FALSE)
  expect_s4_class(down[["SCT"]], "SCTAssay")
  expect_identical(names(slot(down[["SCT"]], "SCTModel.list")), names(slot(sct[["SCT"]], "SCTModel.list")))
  expect_equal(slot(down[["SCT"]], "SCTModel.list")[[1]]@feature.attributes,
               slot(sct[["SCT"]], "SCTModel.list")[[1]]@feature.attributes)
})

test_that("file-based conversion works for RDS and h5Seurat", {
  skip_if_not_installed("withr")
  skip_if_not_installed("hdf5r")
  orig <- make_v5_object()
  rds <- tempfile(fileext = ".rds"); saveRDS(orig, rds)
  v4 <- tempfile(fileext = ".rds")
  v5 <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(c(rds, v4, v5)), add = TRUE)
  DowngradeSeurat(rds, dest = v4, verbose = FALSE)
  expect_identical(class(readRDS(v4)[["RNA"]])[1], "Assay")
  UpgradeSeurat(v4, dest = v5, verbose = FALSE)
  back <- suppressWarnings(LoadH5Seurat(v5, verbose = FALSE))
  expect_s4_class(back[["RNA"]], "Assay5")
  expect_true(all(Layers(orig[["RNA"]]) %in% Layers(back[["RNA"]])))
})

test_that("VisiumV2 images downgrade to VisiumV1 and are restored", {
  demo <- system.file("extdata", "spatial_demo.rds", package = "srtdisk")
  skip_if(!nzchar(demo))
  sp <- readRDS(demo)
  img <- Images(sp)[1]
  skip_if(!inherits(sp[[img]], "VisiumV2"))
  sp4 <- DowngradeSeurat(sp, verbose = FALSE)
  expect_s4_class(sp4[[img]], "VisiumV1")
  expect_gt(sp4[[img]]@spot.radius, 0)
  c0 <- GetTissueCoordinates(sp[[img]])
  c1 <- sp4[[img]]@coordinates
  expect_equal(unname(as.matrix(c1[rownames(c0), c("imagecol", "imagerow")])),
               unname(as.matrix(c0[, c("x", "y")])))
  sp5 <- UpgradeSeurat(sp4, verbose = FALSE)
  expect_s4_class(sp5[[img]], "VisiumV2")
  expect_equal(GetTissueCoordinates(sp5[[img]]), c0)
})
