library(srtdisk)

test_that("assays argument layers-only expands to all assays", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  # Create a minimal test Seurat object
  set.seed(123)

  # Create a small test count matrix
  counts <- matrix(
    rpois(100 * 20, lambda = 5),
    nrow = 100,
    ncol = 20,
    dimnames = list(
      paste0("Gene", 1:100),
      paste0("Cell", 1:20)
    )
  )

  # Create a minimal Seurat object
  suppressWarnings({
    test_obj <- Seurat::CreateSeuratObject(
      counts = counts,
      project = "TestProject",
      min.cells = 0,
      min.features = 0
    )
  })

  # Normalize data to ensure data slot exists
  suppressWarnings({
    test_obj <- Seurat::NormalizeData(test_obj, verbose = FALSE)
  })

  # Verify the test object was created
  expect_s4_class(test_obj, "Seurat")

  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(temp_h5seurat), add = TRUE)

  # Save the object
  expect_no_error(
    SaveH5Seurat(test_obj, filename = temp_h5seurat, overwrite = TRUE, verbose = FALSE)
  )

  expect_true(file.exists(temp_h5seurat))
  expect_gt(file.size(temp_h5seurat), 0)

  # Test loading with assays = c("data") - this should expand to all assays
  obj <- LoadH5Seurat(temp_h5seurat, assays = c("data"), verbose = FALSE)

  # obj may be a Seurat object or a list, handle both
  if (is.list(obj) && !inherits(obj, "Seurat")) {
    obj <- obj[[1]]
  }

  expect_s4_class(obj, "Seurat")
  expect_true(length(names(obj@assays)) >= 1)
  
  # Check that each assay has a data layer/slot
  for (a in names(obj@assays)) {
    assay_obj <- obj@assays[[a]]
    # Check if 'data' slot exists
    has_data_slot <- "data" %in% slotNames(assay_obj)
    
    # For V5 assays, check if data layer exists in assay storage
    has_data_layer <- FALSE
    if ("layers" %in% slotNames(assay_obj)) {
      has_data_layer <- "data" %in% names(assay_obj@layers)
    } else if (methods::hasArg("assay.storage") && "assay.storage" %in% slotNames(assay_obj)) {
      has_data_layer <- "data" %in% names(assay_obj@assay.storage)
    }
    
    expect_true(has_data_slot || has_data_layer, 
                info = paste("missing data layer in assay", a))
  }
})

test_that("assays argument with counts and data layers works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  # Create a minimal test Seurat object
  set.seed(456)

  counts <- matrix(
    rpois(50 * 10, lambda = 5),
    nrow = 50,
    ncol = 10,
    dimnames = list(
      paste0("Gene", 1:50),
      paste0("Cell", 1:10)
    )
  )

  suppressWarnings({
    test_obj <- Seurat::CreateSeuratObject(
      counts = counts,
      project = "TestProject2",
      min.cells = 0,
      min.features = 0
    )
  })

  # Normalize data
  suppressWarnings({
    test_obj <- Seurat::NormalizeData(test_obj, verbose = FALSE)
  })

  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(temp_h5seurat), add = TRUE)

  SaveH5Seurat(test_obj, filename = temp_h5seurat, overwrite = TRUE, verbose = FALSE)

  # Test loading with multiple layers
  obj <- LoadH5Seurat(temp_h5seurat, assays = c("counts", "data"), verbose = FALSE)

  if (is.list(obj) && !inherits(obj, "Seurat")) {
    obj <- obj[[1]]
  }

  expect_s4_class(obj, "Seurat")
  expect_true(length(names(obj@assays)) >= 1)
})
