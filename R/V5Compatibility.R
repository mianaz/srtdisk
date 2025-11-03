#' @include UtilsAssayCompat.R
#' @title Seurat V5 Compatibility Functions
#' @description Functions for handling Seurat V5 object structures
#' @name V5Compatibility
#' @keywords internal
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# V5-specific helper functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Removed: SafeGetAssayData() - duplicate of GetAssayDataCompat() in UtilsAssayCompat.R
# Use GetAssayDataCompat(object, layer_or_slot) instead

# Removed: SafeGetLayers() - duplicate of GetAssayLayersCompat() in UtilsAssayCompat.R
# Use GetAssayLayersCompat(object) instead

#' Transfer metadata from V5 h5Seurat to h5ad
#'
#' @param source Source H5 group (potentially wrapped in environment)
#' @param dfile Destination H5 file
#' @param dname Destination group name
#' @param index Cell index for subsetting
#' @param verbose Logical, print messages
#' @return NULL invisibly
#' @keywords internal
#'
TransferMetadataV5 <- function(source, dfile, dname = "obs", index = NULL, verbose = TRUE) {
  # Handle environment wrapper (common in V5)
  if (is.environment(x = source)) {
    if (exists("hgroup", envir = source, inherits = FALSE)) {
      source <- get("hgroup", envir = source, inherits = FALSE)
    } else {
      # Direct metadata transfer from parent file
      parent_file <- dfile$get_file_id()
      if (parent_file$exists(name = 'meta.data')) {
        source <- parent_file[['meta.data']]
      } else {
        warning("Could not extract metadata from environment wrapper", immediate. = TRUE)
        return(invisible(x = NULL))
      }
    }
  }

  if (!inherits(x = source, what = c('H5D', 'H5Group'))) {
    warning("Invalid source type for metadata transfer", immediate. = TRUE)
    return(invisible(x = NULL))
  }

  if (verbose) {
    message("Transferring metadata to ", dname)
  }

  # Create obs group
  if (!dfile$exists(name = dname)) {
    dfile$create_group(name = dname)
  }

  if (inherits(x = source, what = 'H5Group')) {
    # Transfer each metadata column
    for (col in names(x = source)) {
      if (col == '__categories' || col == '_index') next

      if (IsFactor(x = source[[col]])) {
        # Handle categorical data
        if (!dfile[[dname]]$exists(name = '__categories')) {
          dfile[[dname]]$create_group(name = '__categories')
        }

        # Convert to 0-based indexing for anndata
        values <- if (is.null(index)) {
          source[[col]][['values']][] - 1L
        } else {
          source[[col]][['values']][index] - 1L
        }

        dfile[[dname]]$create_dataset(
          name = col,
          robj = values,
          dtype = source[[col]][['values']]$get_type()
        )

        dfile[[dname]][['__categories']]$create_dataset(
          name = col,
          robj = source[[col]][['levels']][],
          dtype = source[[col]][['levels']]$get_type()
        )
      } else {
        # Regular column
        data <- if (is.null(index)) {
          source[[col]][]
        } else {
          source[[col]][index]
        }

        dfile[[dname]]$create_dataset(
          name = col,
          robj = data,
          dtype = source[[col]]$get_type()
        )
      }
    }

    # Add column order
    if (source$attr_exists(attr_name = 'colnames')) {
      dfile[[dname]]$create_attr(
        attr_name = 'column-order',
        robj = h5attr(x = source, which = 'colnames'),
        dtype = GuessDType(x = h5attr(x = source, which = 'colnames'))
      )
    }
  }

  # Add encoding attributes for anndata
  encoding.info <- c('type' = 'dataframe', 'version' = '0.1.0')
  names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))

  for (i in seq_along(along.with = encoding.info)) {
    attr.name <- names(x = encoding.info)[i]
    if (!dfile[[dname]]$attr_exists(attr_name = attr.name)) {
      dfile[[dname]]$create_attr(
        attr_name = attr.name,
        robj = encoding.info[i],
        dtype = GuessDType(x = encoding.info[i]),
        space = H5S$new(type = "scalar")
      )
    }
  }

  return(invisible(x = NULL))
}

# Note: WriteH5GroupAssay5 is now defined in WriteH5Group.R to avoid duplication

#' Handle spatial data for V5 objects
#'
#' @param object A Seurat object
#' @param hgroup An H5Group to write to
#' @param verbose Logical, print messages
#' @return NULL invisibly
#' @keywords internal
#'
WriteH5SpatialV5 <- function(object, hgroup, verbose = TRUE) {
  # Check for spatial data
  if (length(x = Images(object = object)) > 0) {
    if (verbose) {
      message("Adding spatial information")
    }

    spatial.group <- hgroup$create_group(name = 'images')

    for (img in Images(object = object)) {
      if (verbose) {
        message("Adding image: ", img)
      }

      img.data <- object[[img]]
      img.group <- spatial.group$create_group(name = img)

      # Write coordinates if available
      coords <- GetTissueCoordinates(object = object[[img]])
      if (!is.null(coords)) {
        WriteH5Group(
          x = coords,
          name = 'coordinates',
          hgroup = img.group,
          verbose = verbose
        )
      }

      # Write scale factors
      scale.factors <- ScaleFactors(object = object[[img]])
      if (!is.null(scale.factors)) {
        for (sf.name in names(scale.factors)) {
          img.group$create_attr(
            attr_name = sf.name,
            robj = scale.factors[[sf.name]],
            dtype = GuessDType(x = scale.factors[[sf.name]])
          )
        }
      }

      # Write image data if available and not too large
      # Note: Large images should be handled with care
      if (inherits(x = img.data, what = "VisiumV1")) {
        # Handle Visium-specific data
        img.group$create_attr(
          attr_name = 'assay',
          robj = DefaultAssay(object = img.data),
          dtype = GuessDType(x = DefaultAssay(object = img.data))
        )
      }
    }
  }

  return(invisible(x = NULL))
}

#' Map all Seurat object slots to h5Seurat/h5ad
#'
#' @param format Target format ("h5seurat" or "h5ad")
#' @return A list describing the slot mapping
#' @export
#'
GetSeuratSlotMapping <- function(format = c("h5seurat", "h5ad")) {
  format <- match.arg(format)

  # Comprehensive slot mapping for Seurat objects
  mapping <- list(
    # Core data slots
    assays = list(
      h5seurat = "/assays",
      h5ad = list(
        X = "/X",  # Default assay data
        layers = "/layers",  # Additional layers
        var = "/var",  # Feature metadata
        varm = "/varm"  # Feature embeddings
      )
    ),

    # Cell metadata
    meta.data = list(
      h5seurat = "/meta.data",
      h5ad = "/obs"
    ),

    # Dimensionality reductions
    reductions = list(
      h5seurat = "/reductions",
      h5ad = "/obsm"  # Cell embeddings
    ),

    # Graphs (nearest neighbor, SNN)
    graphs = list(
      h5seurat = "/graphs",
      h5ad = "/obsp"  # Pairwise cell annotations
    ),

    # Spatial data
    images = list(
      h5seurat = "/images",
      h5ad = "/uns/spatial"  # Unstructured spatial data
    ),

    # Cell identities/clusters
    active.ident = list(
      h5seurat = "/active.ident",
      h5ad = "/obs/ident"  # Stored as metadata column
    ),

    # Project information
    project.name = list(
      h5seurat = "/project",
      h5ad = "/uns/project"
    ),

    # Version information
    version = list(
      h5seurat = "/version",
      h5ad = "/uns/version"
    ),

    # Commands/history
    commands = list(
      h5seurat = "/commands",
      h5ad = "/uns/commands"
    ),

    # Miscellaneous data
    misc = list(
      h5seurat = "/misc",
      h5ad = "/uns/misc"
    ),

    # Tools results (e.g., from specialized analyses)
    tools = list(
      h5seurat = "/tools",
      h5ad = "/uns/tools"
    )
  )

  if (format == "h5seurat") {
    return(lapply(mapping, function(x) x$h5seurat))
  } else {
    return(lapply(mapping, function(x) x$h5ad))
  }
}

#' Validate slot mapping completeness
#'
#' @param object A Seurat object
#' @param verbose Logical, print validation results
#' @return Logical indicating if all slots are mappable
#' @keywords internal
#'
ValidateSlotMapping <- function(object, verbose = TRUE) {
  if (!inherits(x = object, what = "Seurat")) {
    stop("Object must be a Seurat object", call. = FALSE)
  }

  # Get all slot names from the object
  all.slots <- slotNames(x = object)

  # Get our mapping
  h5seurat.mapping <- GetSeuratSlotMapping(format = "h5seurat")
  h5ad.mapping <- GetSeuratSlotMapping(format = "h5ad")

  # Check coverage
  mapped.slots <- names(h5seurat.mapping)
  unmapped <- setdiff(all.slots, mapped.slots)

  if (length(unmapped) > 0 && verbose) {
    warning("Unmapped slots: ", paste(unmapped, collapse = ", "))
  }

  # Check for data in each slot
  slot.status <- list()
  for (slot in all.slots) {
    has.data <- !is.null(slot(object = object, name = slot))
    if (slot == "assays") {
      has.data <- length(Assays(object = object)) > 0
    } else if (slot == "reductions") {
      has.data <- length(Reductions(object = object)) > 0
    } else if (slot == "graphs") {
      has.data <- length(Graphs(object = object)) > 0
    } else if (slot == "images") {
      has.data <- length(Images(object = object)) > 0
    }

    slot.status[[slot]] <- list(
      has.data = has.data,
      mapped = slot %in% mapped.slots
    )
  }

  if (verbose) {
    message("Slot mapping validation:")
    for (slot in names(slot.status)) {
      status <- slot.status[[slot]]
      if (status$has.data) {
        symbol <- if (status$mapped) "[+]" else "[-]"
        message("  ", symbol, " ", slot,
                if (!status$mapped) " (unmapped)" else "")
      }
    }
  }

  return(all(sapply(slot.status, function(x) !x$has.data || x$mapped)))
}