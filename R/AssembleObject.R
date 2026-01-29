#' @include zzz.R
#' @include UtilsH5Access.R
#' @importFrom hdf5r h5attr
#' @importFrom methods slot<- new
#' @importFrom Seurat Cells Key<- Key Cells scalefactors
#' @importFrom SeuratObject CreateFOV CreateCentroids
#'
NULL

#' Read sparse matrix from HDF5
#'
#' Internal function to read sparse matrix data from h5Seurat files,
#' handling both V4 and V5 formats
#'
#' @param h5_group The HDF5 group containing sparse matrix data
#' @param verbose Show progress updates
#'
#' @return A matrix object, or NULL if reading fails
#'
#' @keywords internal
#'
ReadSparseMatrix <- function(h5_group, verbose = FALSE) {
  safe_exists <- CreateCachedExistsChecker()

  if (safe_exists(h5_group, "data") && safe_exists(h5_group, "indices") && safe_exists(h5_group, "indptr")) {
    if (verbose) {
      message("Reading sparse matrix in CSR/CSC format")
    }

    data_vals <- h5_group[["data"]][]
    indices <- h5_group[["indices"]][]
    indptr <- h5_group[["indptr"]][]

    nrows <- NA
    ncols <- NA

    if (h5_group$attr_exists("shape")) {
      shape <- h5attr(x = h5_group, which = "shape")
      nrows <- as.integer(shape[1])
      ncols <- as.integer(shape[2])
      if (verbose) {
        message("Read dimensions from 'shape' attribute: ", nrows, " x ", ncols)
      }
    } else if (h5_group$attr_exists("dims")) {
      dims <- h5attr(x = h5_group, which = "dims")
      nrows <- as.integer(dims[1])
      ncols <- as.integer(dims[2])
      if (verbose) {
        message("Read dimensions from 'dims' attribute: ", nrows, " x ", ncols)
      }
    } else {
      if (length(indices) == 0) {
        if (verbose) {
          message("Warning: Empty indices array, cannot infer dimensions")
        }
        return(NULL)
      }
      ncols <- length(indptr) - 1L
      nrows <- max(indices) + 1L
      if (verbose) {
        message("Warning: No shape/dims attribute found, inferring dimensions: ", nrows, " x ", ncols)
      }
    }

    if (is.na(nrows) || is.na(ncols) || nrows <= 0 || ncols <= 0) {
      if (verbose) {
        message("Error: Invalid dimensions (", nrows, " x ", ncols, "), cannot create sparse matrix")
      }
      return(NULL)
    }

    tryCatch({
      if (length(indptr) - 1 == ncols) {
        if (verbose) {
          message("Detected CSC format")
        }
        sparse_mat <- Matrix::sparseMatrix(
          i = indices + 1L,
          p = indptr,
          x = data_vals,
          dims = c(nrows, ncols),
          index1 = TRUE
        )
      } else if (length(indptr) - 1 == nrows) {
        if (verbose) {
          message("Detected CSR format, will transpose")
        }
        sparse_mat <- Matrix::sparseMatrix(
          j = indices + 1L,
          p = indptr,
          x = data_vals,
          dims = c(ncols, nrows),
          index1 = TRUE
        )
        sparse_mat <- Matrix::t(sparse_mat)
      } else {
        if (verbose) {
          message("Error: Cannot determine sparse matrix format. indptr length: ",
                  length(indptr) - 1, ", nrows: ", nrows, ", ncols: ", ncols)
        }
        return(NULL)
      }

      return(sparse_mat)
    }, error = function(e) {
      if (verbose) {
        message("Error creating sparse matrix: ", conditionMessage(e))
        message("Dimensions: ", nrows, " x ", ncols)
        message("indptr length: ", length(indptr))
        message("indices length: ", length(indices))
        message("data length: ", length(data_vals))
      }
      return(NULL)
    })
  }

  if (safe_exists(h5_group, "i") && safe_exists(h5_group, "p") && safe_exists(h5_group, "x")) {
    if (verbose) {
      message("Reading sparse matrix in dgCMatrix format")
    }

    i_vals <- h5_group[["i"]][]
    p_vals <- h5_group[["p"]][]
    x_vals <- h5_group[["x"]][]

    if (h5_group$attr_exists("Dim")) {
      dims <- as.integer(h5attr(x = h5_group, which = "Dim"))
      if (verbose) {
        message("Read dimensions from 'Dim' attribute: ", dims[1], " x ", dims[2])
      }
    } else {
      if (length(i_vals) == 0) {
        if (verbose) {
          message("Warning: Empty i_vals array, cannot infer dimensions")
        }
        return(NULL)
      }
      dims <- c(max(i_vals) + 1L, length(p_vals) - 1L)
      if (verbose) {
        message("Inferred dimensions: ", dims[1], " x ", dims[2])
      }
    }

    if (length(dims) < 2 || any(is.na(dims)) || any(dims <= 0)) {
      if (verbose) {
        message("Error: Invalid dimensions (", paste(dims, collapse = " x "), "), cannot create sparse matrix")
      }
      return(NULL)
    }

    tryCatch({
      sparse_mat <- Matrix::sparseMatrix(
        i = i_vals + 1L,
        p = p_vals,
        x = x_vals,
        dims = dims,
        index1 = TRUE
      )

      return(sparse_mat)
    }, error = function(e) {
      if (verbose) {
        message("Error creating dgCMatrix sparse matrix: ", conditionMessage(e))
        message("Dimensions: ", paste(dims, collapse = " x "))
        message("i_vals length: ", length(i_vals))
        message("p_vals length: ", length(p_vals))
        message("x_vals length: ", length(x_vals))
      }
      return(NULL)
    })
  }

  # If it's not a sparse matrix group, try reading as dense matrix
  if (verbose) {
    message("No sparse matrix components found, attempting dense matrix read")
  }

  tryCatch({
    # Read the matrix data
    mat_data <- h5_group[,]

    # Check if we should convert dense to sparse
    if (is.matrix(mat_data)) {
      zero_proportion <- sum(mat_data == 0) / length(mat_data)
      if (zero_proportion > 0.5) {
        if (verbose) {
          message("Converting dense to sparse (", round(zero_proportion * 100, 1), "% zeros)")
        }
        return(as(mat_data, "dgCMatrix"))
      }
    }

    return(mat_data)
  }, error = function(e) {
    if (verbose) {
      message("Failed to read matrix: ", conditionMessage(e))
    }
    return(NULL)
  })
}

#' Assemble an object from an h5Seurat file
#'
#' @param assay,reduction,graph,image,neighbor,cmd Name of assay, reduction,
#' graph, image, neighbor, or command to load
#' @param file A connected h5Seurat file to pull the data from
#' @param verbose Show progress updates
#'
#' @name AssembleObject
#' @rdname AssembleObject
#'
#' @keywords internal
#'
NULL

#' @param slots Optional vector of assay slots to load, defaults to all slots
#' present in assay
#'
#' @importFrom Seurat CreateAssayObject GetAssayData SetAssayData AddMetaData
#' VariableFeatures<-
#'
#' @return \code{AssembleAssay}: An \code{Assay} object
#'
#' @rdname AssembleObject
#'
#' @aliases AssembleAssay
#'
AssembleAssay <- function(assay, file, slots = NULL, verbose = TRUE) {
  index <- file$index()
  if (!assay %in% names(x = index)) {
    stop("Cannot find assay ", assay, " in this h5Seurat file", call. = FALSE)
  }
  slots.assay <- names(x = Filter(f = isTRUE, x = index[[assay]]$slots))
  
  # Handle NULL slots (means load all available slots)
  if (is.null(slots)) {
    slots <- slots.assay
  }
  
  # Check that slots is a non-empty character vector and contains valid values before calling match.arg
  if (is.null(slots) || length(slots) == 0 || !any(slots %in% slots.assay)) {
    stop("'slots' must be a non-empty character vector. Available slots for assay '", 
         assay, "': ", paste(slots.assay, collapse = ", "), call. = FALSE)
  }
  
  slots <- match.arg(arg = slots, choices = slots.assay, several.ok = TRUE)
  if (!any(c('counts', 'data') %in% slots)) {
    stop("At least one of 'counts' or 'data' must be loaded", call. = FALSE)
  }
  assay.group <- file[['assays']][[assay]]

  safe_exists <- CreateCachedExistsChecker()

  safe_read_dataset <- function(dataset) {
    as.character(SafeH5DRead(dataset))
  }

  features <- GetFeaturesV5Safe(h5_group = assay.group, verbose = verbose)
  if (is.null(features)) {
    features <- safe_read_dataset(assay.group[['features']])
  }
  features <- FixFeatures(features = features)
  if (verbose) {
    message("Loaded ", length(features), " features for assay")
  }

  read_matrix_data <- function(slot_name) {
    if (exists("ReadV5Layer", envir = asNamespace("srtdisk"))) {
      mat <- ReadV5Layer(
        h5_group = assay.group,
        layer_name = slot_name,
        features = features,
        cells = Cells(x = file),
        verbose = verbose
      )
      if (!is.null(mat)) {
        return(mat)
      }
    }

    # Fallback to direct reading if ReadV5Layer not available or returns NULL
    # First try direct matrix (V4 style)
    if (safe_exists(assay.group, slot_name)) {
      if (verbose) {
        message("Reading direct matrix for slot '", slot_name, "'")
      }
      # Use ReadSparseMatrix which handles both sparse and dense formats
      mat <- ReadSparseMatrix(h5_group = assay.group[[slot_name]], verbose = verbose)
      if (!is.null(mat)) {
        # Ensure proper dimensions for V5 compatibility
        if (ncol(mat) != length(Cells(x = file)) || nrow(mat) != length(features)) {
          # Check if transposed
          if (nrow(mat) == length(Cells(x = file)) && ncol(mat) == length(features)) {
            if (verbose) {
              message("Transposing matrix to match expected dimensions")
            }
            mat <- t(mat)
          }
        }
        return(mat)
      } else {
        # ReadSparseMatrix returned NULL, try V5 sparse layers instead
        if (verbose) {
          message("Direct matrix read returned NULL, trying V5 layers")
        }
      }
    }

    # Try V5 sparse layers structure
    layers_path <- paste0("layers/", slot_name)
    if (safe_exists(assay.group, layers_path)) {
      if (verbose) {
        message("Reading V5 sparse layers for slot '", slot_name, "'")
      }
      layer_group <- assay.group[[layers_path]]

      # Use ReadSparseMatrix for V5 layers
      mat <- ReadSparseMatrix(h5_group = layer_group, verbose = verbose)

      # V5 sparse matrices might need special dimension handling
      if (!is.null(mat) && length(features) > 0) {
        expected_nrows <- length(features)
        expected_ncols <- length(Cells(x = file))

        # Check dimensions and adjust if needed
        if (nrow(mat) != expected_nrows || ncol(mat) != expected_ncols) {
          # If transposed
          if (nrow(mat) == expected_ncols && ncol(mat) == expected_nrows) {
            if (verbose) {
              message("Transposing V5 matrix to match expected dimensions")
            }
            mat <- t(mat)
          } else if (nrow(mat) < expected_nrows) {
            # V5 might have filtered features, need to expand
            if (verbose) {
              message("Expanding sparse matrix to full feature space: ", nrow(mat), " -> ", expected_nrows)
            }
            # Create full matrix with zeros for missing features
            full_mat <- matrix(0, nrow = expected_nrows, ncol = ncol(mat))
            # Copy existing data (assuming first features match)
            full_mat[1:nrow(mat), ] <- mat
            mat <- full_mat
          }
        }
      }

      return(mat)
    }

    # If neither direct nor sparse layers exist, return NULL instead of error
    # This allows more graceful fallback behavior
    if (verbose) {
      message("Warning: Cannot find matrix data for slot '", slot_name, "' in assay '", assay, "'")
    }
    return(NULL)
  }
  # Load the appropriate slots - prioritize counts for initialization
  # We need to load counts first if available, then add other layers
  # Initialize assay with counts or data
  init_with_counts <- 'counts' %in% slots
  init_slot_name <- if (init_with_counts) 'counts' else 'data'

  if (verbose) {
    message("Initializing ", assay, " with ", init_slot_name)
  }

  init_data <- read_matrix_data(init_slot_name)
  if (is.null(init_data)) {
    stop("Failed to read ", init_slot_name, " matrix for assay '", assay, "'", call. = FALSE)
  }

  rownames(x = init_data) <- features
  colnames(x = init_data) <- Cells(x = file)

  # Create V5-compatible assay
  # CreateSeuratObject always puts data in 'counts' layer initially
  temp_seurat <- CreateSeuratObject(counts = init_data, min.cells = -1, min.features = -1)
  obj <- temp_seurat[['RNA']]

  # If we initialized with 'data' (not counts), we need to move the data to the correct layer
  # The data is currently in 'counts' layer but should be in 'data' layer
  if (!init_with_counts) {
    # Move data from counts to data layer, clear counts
    obj <- SetAssayDataCompat(object = obj, layer_or_slot = "data", new.data = init_data)
    # Remove the incorrectly placed counts layer by setting empty matrix
    empty_counts <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(nrow(init_data), ncol(init_data)),
      dimnames = list(features, Cells(x = file))
    )
    obj <- SetAssayDataCompat(object = obj, layer_or_slot = "counts", new.data = empty_counts)
  }

  tryCatch(
    expr = Key(object = obj) <- Key(object = assay.group),
    error = function(e) NULL
  )
  # Add remaining slots/layers (V5 compatibility)
  # Determine which slot was used for initialization
  init_slot <- if ('counts' %in% slots) 'counts' else 'data'

  for (slot in slots) {
    # Skip the slot used for initialization
    if (slot == init_slot) {
      next
    }

    slot_empty <- tryCatch({
      GetAssayDataCompat(object = obj, layer_or_slot = slot) |>
        IsMatrixEmpty()
    }, error = function(e) {
      TRUE
    })
    
    if (slot_empty) {
      if (verbose) {
        message("Adding ", slot, " for ", assay)
      }
      tryCatch({
        dat <- read_matrix_data(slot)

        if (is.null(dat)) {
          if (verbose) {
            message("Skipping slot '", slot, "' - not found in file")
          }
          next
        }

        colnames(x = dat) <- Cells(x = file)
        rownames(x = dat) <- if (slot == 'scale.data') {
          FixFeatures(features = safe_read_dataset(assay.group[['scaled.features']]))
        } else {
          features
        }

        obj <- SetAssayDataCompat(object = obj, layer_or_slot = slot, new.data = dat)

      }, error = function(e) {
        if (verbose) {
          message("Skipping slot '", slot, "' - failed to set: ", conditionMessage(e))
        }
      })
    } else {
      if (verbose) {
        message("Slot '", slot, "' already populated, skipping")
      }
    }
  }
  # Add meta features
  if (safe_exists(assay.group, 'meta.features')) {
    if (verbose) {
      message("Adding feature-level metadata for ", assay)
    }
    meta.data <- as.data.frame(
      x = assay.group[['meta.features']],
      row.names = features
    )
    if (ncol(x = meta.data)) {
      obj <- AddMetaData(
        object = obj,
        metadata = meta.data
      )
    }
  }
  # Add variable feature information
  if (safe_exists(assay.group, 'variable.features')) {
    if (verbose) {
      message("Adding variable feature information for ", assay)
    }
    VariableFeatures(object = obj) <- safe_read_dataset(assay.group[['variable.features']])
  }
  # Add miscellaneous information
  if (safe_exists(assay.group, 'misc')) {
    if (verbose) {
      message("Adding miscellaneous information for ", assay)
    }
    tryCatch({
      if ('misc' %in% slotNames(obj)) {
        slot(object = obj, name = 'misc') <- SafeH5GroupToList(h5obj = assay.group[['misc']], recursive = TRUE)
      } else {
        if (verbose) {
          message("Misc slot not available in this assay type (Assay5), skipping")
        }
      }
    }, error = function(e) {
      if (verbose) {
        message("Failed to set misc information: ", conditionMessage(e))
      }
    })
  }
  # Handle S4 class reconstruction - skip for Assay5 objects
  if (!inherits(obj, "Assay5") && assay.group$attr_exists(attr_name = 's4class')) {
    tryCatch({
      classdef <- unlist(x = strsplit(
        x = h5attr(x = assay.group, which = 's4class'),
        split = ':'
      ))
      pkg <- classdef[1]
      cls <- classdef[2]

      formal <- methods::getClassDef(Class = cls, package = pkg, inherits = FALSE)
      missing <- setdiff(
        x = slotNames(x = formal),
        y = slotNames(x = methods::getClass(Class = 'Assay'))
      )
      missing <- intersect(x = missing, y = names(x = assay.group))
      missing <- sapply(
        X = missing,
        FUN = function(x) as.list(x = assay.group[[x]], recursive = TRUE),
        simplify = FALSE
      )
      obj <- c(SeuratObject::S4ToList(object = obj), missing)
      attr(x = obj, which = 'classDef') <- paste(classdef, collapse = ':')
      obj <- SeuratObject::ListToS4(x = obj)
    }, error = function(e) {
      if (verbose) {
        message("S4 class reconstruction failed: ", conditionMessage(e))
      }
    })
  }
  return(obj)
}

#' @importClassesFrom Seurat JackStrawData
#' @importFrom Seurat CreateDimReducObject Stdev IsGlobal JS<-
#'
#' @rdname AssembleObject
#'
AssembleDimReduc <- function(reduction, file, verbose = TRUE) {
  safe_read_dataset <- function(dataset) {
    as.character(SafeH5DRead(dataset))
  }

  index <- file$index()
  index.check <- vapply(
    X = setdiff(x = names(x = index), y = c('global', 'no.assay')),
    FUN = function(x) {
      return(reduction %in% names(x = index[[x]]$reductions))
    },
    FUN.VALUE = logical(length = 1L)
  )
  if (!any(index.check)) {
    stop(
      "Cannot find reduction ",
      reduction,
      " in this h5Seurat file",
      call. = FALSE
    )
  } else if (sum(index.check) > 1) {
    stop("Multiple reductions named ", reduction, call. = FALSE)
  }
  assay <- names(x = which(x = index.check))
  reduc.group <- file[['reductions']][[reduction]]
  key <- Key(object = reduc.group)
  # Pull cell embeddings
  if (index[[assay]]$reductions[[reduction]][['cell.embeddings']]) {
    if (verbose) {
      message("Adding cell embeddings for ", reduction)
    }
    embeddings <- as.matrix(x = reduc.group[['cell.embeddings']])
    rownames(x = embeddings) <- Cells(x = file)
    colnames(x = embeddings) <- paste0(key, 1:ncol(x = embeddings))
  } else {
    if (verbose) {
      warning(
        "No cell embeddings for ",
        reduction,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    embeddings <- new(Class = 'matrix')
  }
  # Pull feature loadings
  if (index[[assay]]$reductions[[reduction]][['feature.loadings']]) {
    if (verbose) {
      message("Adding feature loadings for ", reduction)
    }
    loadings <- as.matrix(x = reduc.group[['feature.loadings']])
    rownames(x = loadings) <- safe_read_dataset(reduc.group[['features']])
    colnames(x = loadings) <- paste0(key, 1:ncol(x = loadings))
  } else {
    loadings <- new(Class = 'matrix')
  }
  # Pull projected loadings
  if (index[[assay]]$reductions[[reduction]][['feature.loadings.projected']]) {
    if (verbose) {
      message("Adding projected loadings for ", reduction)
    }
    projected <- as.matrix(x = reduc.group[['feature.loadings.projected']])
    rownames(x = projected) <- safe_read_dataset(reduc.group[['projected.features']])
    colnames(x = projected) <- paste0(key, 1:ncol(x = projected))
  } else {
    projected <- new(Class = 'matrix')
  }
  # Build the object
  obj <- CreateDimReducObject(
    embeddings = embeddings,
    loadings = loadings,
    projected = projected,
    assay = assay,
    stdev = Stdev(object = file, reduction = reduction),
    key = key,
    global = IsGlobal(object = reduc.group)
  )
  # Add misc
  if (reduc.group$exists(name = 'misc')) {
    if (verbose) {
      message("Adding miscellaneous information for ", reduction)
    }
    # Use SafeH5GroupToList to handle 3D+ arrays (e.g., from Squidpy UMAP)
    slot(object = obj, name = 'misc') <- SafeH5GroupToList(h5obj = reduc.group[['misc']], recursive = TRUE)
  }
  # Add jackstraw
  if (index[[assay]]$reductions[[reduction]][['jackstraw']]) {
    if (verbose) {
      message("Loading JackStraw data for ", reduction)
    }
    js <- new(Class = 'JackStrawData')
    for (slot in names(x = reduc.group[['jackstraw']])) {
      JS(object = js, slot = slot) <- as.matrix(x = reduc.group[['jackstraw']][[slot]])
    }
    JS(object = obj) <- js
  }
  return(obj)
}

#' @importFrom Seurat as.sparse as.Graph
#'
#' @rdname AssembleObject
#'
AssembleGraph <- function(graph, file, verbose = TRUE) {
  index <- file$index()
  obj <- as.sparse(x = file[['graphs']][[graph]])
  rownames(x = obj) <- colnames(x = obj) <- Cells(x = file)
  obj <- as.Graph(x = obj)
  if (file[['graphs']][[graph]]$attr_exists(attr_name = 'assay.used')) {
    assay <- h5attr(x = file[['graphs']][[graph]], which = 'assay.used')
    if (graph %in% index[[assay]]$graphs) {
      DefaultAssay(object = obj) <- assay
    }
  }
  return(obj)
}

#' @rdname AssembleObject
#'
AssembleImage <- function(image, file, verbose = TRUE) {
  index <- file$index()
  img_group <- file[['images']][[image]]

  # Get assay information
  assay <- if (img_group$attr_exists(attr_name = 'assay')) {
    h5attr(x = img_group, which = 'assay')
  } else {
    'Spatial'  # Default assay name for spatial data
  }

  # Get s4class to determine object type
  s4class <- if (img_group$attr_exists(attr_name = 's4class')) {
    h5attr(x = img_group, which = 's4class')
  } else {
    NULL
  }

  # Construct spatial image object (VisiumV1/VisiumV2/SliceImage)
  if (!is.null(s4class) && s4class %in% c('VisiumV1', 'VisiumV2', 'SliceImage')) {
    tryCatch({
      # Get image data
      image_data <- NULL
      if (img_group$exists(name = 'image')) {
        image_data <- img_group[['image']]$read()
        # Image data should be in Seurat format: height x width x channels
        if (length(dim(image_data)) == 3L) {
          # Normalize to 0-1 range if needed
          if (max(image_data) > 1) {
            image_data <- image_data / 255
          }
          storage.mode(image_data) <- 'double'
        }
      }

      # Get scale factors
      scale_factors <- NULL
      if (img_group$exists(name = 'scale.factors')) {
        sf_group <- img_group[['scale.factors']]

        # Read scale factor values
        spot <- if (sf_group$exists('spot')) sf_group[['spot']][] else NA_real_
        fiducial <- if (sf_group$exists('fiducial')) sf_group[['fiducial']][] else NA_real_
        hires <- if (sf_group$exists('hires')) sf_group[['hires']][] else NA_real_
        lowres <- if (sf_group$exists('lowres')) sf_group[['lowres']][] else NA_real_

        # Create scalefactors object
        scale_factors <- tryCatch({
          scalefactors(
            spot = as.numeric(spot),
            fiducial = as.numeric(fiducial),
            hires = as.numeric(hires),
            lowres = as.numeric(lowres)
          )
        }, error = function(e) {
          if (verbose) {
            message("Could not create scalefactors object: ", conditionMessage(e))
          }
          NULL
        })
      }

      # Get spatial coordinates and boundaries - different for VisiumV1 vs VisiumV2
      coordinates <- NULL
      boundaries_list <- list()

      if (verbose) {
        message("Assembling spatial image: ", image)
        message("  s4class: ", s4class)
        message("  has image_data: ", !is.null(image_data))
        message("  has scale_factors: ", !is.null(scale_factors))
      }

      # Try boundaries/centroids (VisiumV2 style)
      if (img_group$exists(name = 'boundaries')) {
        boundaries_group <- img_group[['boundaries']]
        if (boundaries_group$exists(name = 'centroids')) {
          centroids_group <- boundaries_group[['centroids']]

          # Read all centroids data
          if (centroids_group$exists(name = 'coords')) {
            if (verbose) {
              message("Reading spatial coordinates from boundaries/centroids")
            }

            coords_mat <- centroids_group[['coords']][,]
            cell_names <- if (centroids_group$exists('cells')) {
              as.character(centroids_group[['cells']][])
            } else {
              Cells(x = file)
            }

            # Note: Centroids are already filtered per library during conversion
            # No need to filter again here

            # Read centroid parameters
            radius_val <- if (centroids_group$exists('radius')) {
              centroids_group[['radius']][]
            } else {
              as.numeric(scale_factors[['spot']])
            }

            theta_val <- if (centroids_group$exists('theta')) {
              centroids_group[['theta']][]
            } else {
              0
            }

            nsides_val <- if (centroids_group$exists('nsides')) {
              centroids_group[['nsides']][]
            } else {
              0L  # 0 = infinite sides (circle), returns centroids not polygon vertices
            }

            # Create Centroids object using SeuratObject function
            tryCatch({
              # Ensure coords has proper column names and rownames
              if (is.null(colnames(coords_mat))) {
                colnames(coords_mat) <- c('x', 'y')
              }

              # Set rownames to cell names - CreateCentroids uses rownames for cells
              if (length(cell_names) == nrow(coords_mat)) {
                rownames(coords_mat) <- cell_names
              }

              centroids_obj <- SeuratObject::CreateCentroids(
                coords = coords_mat,
                nsides = as.integer(nsides_val),
                radius = as.numeric(radius_val),
                theta = as.numeric(theta_val)
              )

              boundaries_list[['centroids']] <- centroids_obj

              if (verbose) {
                message("Created Centroids object with ", nrow(coords_mat), " cells")
              }
            }, error = function(e) {
              if (verbose) {
                message("Could not create Centroids object: ", conditionMessage(e))
              }
            })
          }
        }
      }

      if (verbose) {
        message("  boundaries_list length after boundaries check: ", length(boundaries_list))
      }

      # Fallback: If no boundaries/centroids found, try to create from spatial reduction
      # This handles native h5ad files from scanpy/squidpy that store coords in obsm/spatial
      if (length(boundaries_list) == 0 && file$exists(name = 'reductions/spatial')) {
        if (verbose) {
          message("  No boundaries found, attempting to create centroids from spatial reduction")
        }

        tryCatch({
          spatial_reduc <- file[['reductions/spatial']]
          if (spatial_reduc$exists(name = 'cell.embeddings')) {
            coords_mat <- as.matrix(spatial_reduc[['cell.embeddings']])
            all_cells <- Cells(x = file)

            # Filter coordinates to cells from this library if applicable
            cells_to_keep <- all_cells
            if (file$exists(name = 'meta.data')) {
              meta_group <- file[['meta.data']]
              lib_col_names <- c('sangerID', 'library_id', 'sample', 'batch')

              for (col_name in lib_col_names) {
                if (meta_group$exists(col_name)) {
                  col_obj <- meta_group[[col_name]]

                  if (inherits(col_obj, 'H5Group')) {
                    if (col_obj$exists('values') && col_obj$exists('levels')) {
                      values_int <- col_obj[['values']]$read()
                      levels_str <- col_obj[['levels']]$read()
                      lib_ids <- levels_str[values_int + 1]
                    } else {
                      lib_ids <- NULL
                    }
                  } else {
                    lib_ids <- as.character(col_obj$read())
                  }

                  if (!is.null(lib_ids)) {
                    cells_to_keep <- all_cells[lib_ids == image]
                    if (verbose) {
                      message("  Filtering to ", length(cells_to_keep), " cells from library ", image)
                    }
                    break
                  }
                }
              }
            }

            # Filter and prepare coordinates
            cell_indices <- which(all_cells %in% cells_to_keep)
            coords_mat_filtered <- coords_mat[cell_indices, , drop = FALSE]

            # Ensure proper column names
            if (is.null(colnames(coords_mat_filtered))) {
              colnames(coords_mat_filtered) <- c('x', 'y')
            }
            rownames(coords_mat_filtered) <- cells_to_keep

            # Create Centroids object for VisiumV2 compatibility
            radius_val <- if (!is.null(scale_factors)) {
              as.numeric(scale_factors[['spot']])
            } else {
              1  # Default radius
            }

            centroids_obj <- SeuratObject::CreateCentroids(
              coords = coords_mat_filtered,
              nsides = 0L,  # Circle
              radius = radius_val,
              theta = 0
            )

            boundaries_list[['centroids']] <- centroids_obj

            if (verbose) {
              message("  Created Centroids from spatial reduction with ", nrow(coords_mat_filtered), " cells")
            }
          }
        }, error = function(e) {
          if (verbose) {
            message("  Could not create centroids from spatial reduction: ", conditionMessage(e))
          }
        })
      }

      if (verbose) {
        message("  Final boundaries_list length: ", length(boundaries_list))
      }

      # For VisiumV1 and SliceImage, read coordinates from reductions or image group
      if (s4class %in% c('VisiumV1', 'SliceImage')) {
        if (file$exists(name = 'reductions/spatial')) {
          if (verbose) {
            message("Reading spatial coordinates from reduction for VisiumV1/SliceImage")
          }
          spatial_reduc <- file[['reductions/spatial']]
          if (spatial_reduc$exists(name = 'cell.embeddings')) {
            coords_mat <- as.matrix(spatial_reduc[['cell.embeddings']])
            all_cells <- Cells(x = file)

            # Filter coordinates to only cells from this library
            # Check for library ID in meta.data (could be sangerID, library_id, etc.)
            cells_to_keep <- all_cells
            if (file$exists(name = 'meta.data')) {
              meta_group <- file[['meta.data']]
              # Try common library ID column names
              lib_col_names <- c('sangerID', 'library_id', 'sample', 'batch')
              lib_col_found <- FALSE

              for (col_name in lib_col_names) {
                if (meta_group$exists(col_name)) {
                  col_obj <- meta_group[[col_name]]

                  # Handle both simple datasets and factor/categorical structures
                  if (inherits(col_obj, 'H5Group')) {
                    # It's a factor with values and levels
                    if (col_obj$exists('values') && col_obj$exists('levels')) {
                      values_int <- col_obj[['values']]$read()  # 0-indexed integers
                      levels_str <- col_obj[['levels']]$read()   # level strings
                      # Convert to 1-indexed for R and map to levels
                      lib_ids <- levels_str[values_int + 1]
                    } else {
                      lib_ids <- NULL
                    }
                  } else {
                    # Simple dataset
                    lib_ids <- as.character(col_obj$read())
                  }

                  if (!is.null(lib_ids)) {
                    # Filter to cells matching this image name
                    cells_to_keep <- all_cells[lib_ids == image]
                    lib_col_found <- TRUE
                    if (verbose) {
                      message("  Filtering to ", length(cells_to_keep), " cells from library ", image)
                    }
                    break
                  }
                }
              }

              if (!lib_col_found && verbose) {
                message("  Warning: No library ID column found, using all cells")
              }
            }

            # Filter coords_mat to only include matching cells
            cell_indices <- which(all_cells %in% cells_to_keep)
            coords_mat_filtered <- coords_mat[cell_indices, , drop = FALSE]

            # Create coordinate dataframe for VisiumV1
            # reductions/spatial/cell.embeddings is transposed from h5ad obsm/X_spatial
            # h5ad stores [X, Y] = [imagecol, imagerow] format (verified with scanpy)
            # After transpose: column 1 = imagecol (X), column 2 = imagerow (Y)
            coordinates <- data.frame(
              imagerow = coords_mat_filtered[, 2],  # Y from column 2
              imagecol = coords_mat_filtered[, 1],  # X from column 1
              row.names = cells_to_keep,
              stringsAsFactors = FALSE
            )
          }
        } else if (img_group$exists(name = 'coordinates')) {
          if (verbose) {
            message("Reading spatial coordinates from image group (SliceImage format)")
          }
          tryCatch({
            coords_group <- img_group[['coordinates']]
            coordinates <- as.data.frame(x = coords_group, row.names = Cells(x = file))
          }, error = function(e) {
            if (verbose) {
              message("Could not read coordinates from image group: ", conditionMessage(e))
            }
          })
        }
      }

      # Get key
      image_key <- if (img_group$exists(name = 'key')) {
        img_group[['key']][]
      } else {
        paste0(image, '_')
      }

      # Only attempt to create a full Visium object if we have the required components
      can_create_v2 <- !is.null(image_data) && !is.null(scale_factors) && length(boundaries_list) > 0
      can_create_v1 <- !is.null(image_data) && !is.null(scale_factors) && !is.null(coordinates)

      if (verbose) {
        message("  Can create VisiumV2: ", can_create_v2,
                " (image: ", !is.null(image_data),
                ", scale_factors: ", !is.null(scale_factors),
                ", boundaries: ", length(boundaries_list), ")")
        message("  Can create VisiumV1: ", can_create_v1,
                " (image: ", !is.null(image_data),
                ", scale_factors: ", !is.null(scale_factors),
                ", coordinates: ", !is.null(coordinates), ")")
      }

      if (s4class == 'VisiumV2' && can_create_v2) {
        # VisiumV2 object
        if (verbose) {
          n_cells <- if (!is.null(boundaries_list$centroids)) {
            length(slot(boundaries_list$centroids, 'cells'))
          } else {
            0
          }
          message("Creating VisiumV2 object with ", n_cells, " spots")
        }

        tryCatch({
          obj <- new(
            Class = 'VisiumV2',
            image = image_data,
            scale.factors = scale_factors,
            molecules = list(),
            boundaries = boundaries_list,
            coords_x_orientation = "horizontal",
            assay = assay,
            key = image_key
          )
          return(obj)
        }, error = function(e) {
          if (verbose) {
            message("Error creating VisiumV2 object: ", conditionMessage(e))
          }
        })
      } else if (s4class %in% c('VisiumV1', 'SliceImage') && can_create_v1) {
        if (verbose) {
          message("Creating ", s4class, " object with ", nrow(coordinates), " spots")
        }

        tryCatch({
          obj <- new(
            Class = s4class,
            image = image_data,
            scale.factors = scale_factors,
            coordinates = coordinates,
            spot.radius = as.numeric(scale_factors[['spot']]) / 2,
            assay = assay,
            key = image_key
          )
          return(obj)
        }, error = function(e) {
          if (verbose) {
            message("Error creating ", s4class, " object: ", conditionMessage(e))
          }
        })
      }

      # If we get here, we couldn't create a full object
      if (verbose) {
        message("Could not create full ", s4class, " object, falling back to list")
      }
    }, error = function(e) {
      if (verbose) {
        message("Error creating Visium object: ", conditionMessage(e))
        message("Falling back to simple list structure")
      }
    })
  }

  # Fallback: Use simple list structure (original behavior)
  obj <- tryCatch({
    as.list(x = img_group, recursive = TRUE, row.names = Cells(x = file))
  }, error = function(e) {
    if (verbose) {
      message("Standard as.list failed for image, using safe conversion: ", conditionMessage(e))
    }
    SafeH5GroupToList(h5obj = img_group, recursive = TRUE)
  })

  # Try to set assay attribute
  tryCatch({
    # Check if image is in the index - wrap in tryCatch to handle structured objects
    if (!is.null(index[[assay]]$images) && image %in% index[[assay]]$images) {
      tryCatch({
        DefaultAssay(object = obj) <- assay
      }, error = function(e) {
        if (verbose) {
          message("Note: Setting assay as attribute instead of DefaultAssay")
        }
        attr(obj, 'assay') <- assay
      })
    } else {
      # If not in index or index check failed, just set assay attribute
      if (verbose) {
        message("Setting assay as attribute for image")
      }
      attr(obj, 'assay') <- assay
    }
  }, error = function(e) {
    # If the index check itself fails, just set the assay attribute
    if (verbose) {
      message("Could not check image index, setting assay as attribute: ", conditionMessage(e))
    }
    attr(obj, 'assay') <- assay
  })

  return(obj)
}

#' @importClassesFrom Seurat Neighbor
#'
#' @rdname AssembleObject
#'
AssembleNeighbor <- function(neighbor, file, verbose = TRUE) {
  neighbor.group <- file[['neighbors']][[neighbor]]
  obj <- new(
    Class = 'Neighbor',
    nn.idx =  as.matrix(x = neighbor.group[["nn.idx"]]),
    nn.dist = as.matrix(x = neighbor.group[["nn.dist"]]),
    cell.names =  as.matrix(x = neighbor.group[["cell.names"]])[,1]
  )
  return(obj)
}

#' @importClassesFrom Seurat SeuratCommand
#'
#' @rdname AssembleObject
#'
AssembleSeuratCommand <- function(cmd, file, verbose = TRUE) {
  index <- file$index()
  index.check <- vapply(
    X = setdiff(x = names(x = index), y = 'global'),
    FUN = function(x) {
      return(cmd %in% index[[x]]$commands)
    },
    FUN.VALUE = logical(length = 1L)
  )
  if (!any(index.check)) {
    stop("Cannot find command ", cmd, " in this h5Seurat file", call. = FALSE)
  } else if (sum(index.check) > 1) {
    stop("Multiple commands named ", cmd, call. = FALSE)
  }
  cmd.group <- file[['commands']][[cmd]]
  cmdlog <- new(
    Class = 'SeuratCommand',
    name = h5attr(x = cmd.group, which = 'name'),
    time.stamp = as.POSIXct(x = h5attr(x = cmd.group, which = 'time.stamp')),
    call.string = h5attr(x = cmd.group, 'call.string'),
    params = as.list(x = cmd.group, recursive = TRUE)
  )
  assay <- names(x = which(x = index.check))
  if (assay != 'no.assay') {
    slot(object = cmdlog, name = 'assay.used') <- assay
  }
  return(cmdlog)
}
