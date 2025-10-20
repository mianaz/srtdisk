#' @include zzz.R
#' @include h5Seurat.R
#' @include GetObject.R
#' @include AssembleObject.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load a Seurat object from an h5Seurat file
#'
#' Load a previously saved Seurat object from an h5Seurat file. This function supports
#' flexible loading options, allowing you to load only the components you need (e.g.,
#' specific assays, reductions) to minimize memory usage on large datasets.
#'
#' @param file,x Name of an h5Seurat file path (character) or connected h5Seurat file to load
#' @param assays One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all assays
#'  \item A character vector with names of assays to load (e.g., \code{c("RNA", "ADT")})
#'  \item A character vector specifying which data layers to load for all assays:
#'    \code{c("counts", "data")} loads only the counts and data layers, skipping scale.data
#'  \item A named list for fine-grained control, e.g., \code{list(RNA = "data", ADT = c("data", "scale.data"))}
#' }
#' @param reductions One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all reductions (PCA, UMAP, etc.)
#'  \item A character vector with names of specific reductions (e.g., \code{c("pca", "umap")})
#'  \item \code{NA}: Load only global (assay-independent) reductions
#'  \item \code{FALSE}: Skip loading all reductions
#' }
#' \strong{Note}: Only reductions associated with a loaded assay or marked as global will be loaded.
#' @param graphs One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all graphs
#'  \item A character vector with specific graph names (e.g., \code{c("RNA_snn", "ADT_snn")})
#'  \item \code{FALSE}: Skip loading graphs
#' }
#' \strong{Note}: Only graphs associated with loaded assays will be available.
#' @param neighbors One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all neighbor information
#'  \item A character vector with neighbor names
#'  \item \code{FALSE}: Skip neighbors
#' }
#' @param images One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all images (for spatial experiments)
#'  \item A character vector with image names
#'  \item \code{NA}: Load only global images
#'  \item \code{FALSE}: Skip images
#' }
#' @param meta.data Logical; if \code{TRUE} (default), load cell metadata
#' @param commands Logical; if \code{TRUE} (default), load command history.
#'   Commands are only loaded if their associated assays are loaded.
#' @param misc Logical; if \code{TRUE} (default when all assays loaded), load miscellaneous data
#' @param tools Logical; if \code{TRUE} (default when all assays loaded), load tool-specific information
#' @param verbose Logical; if \code{TRUE} (default), show progress messages
#' @param ... Arguments passed to other methods
#'
#' @return A \code{Seurat} object containing the requested components
#'
#' @details
#' The h5Seurat format is highly flexible for selective loading. This is particularly useful when:
#' \itemize{
#'   \item Working with very large datasets where loading everything would exceed memory
#'   \item You only need specific assays or reductions for downstream analysis
#'   \item You want to quickly inspect object structure without full data loading
#' }
#'
#' @section Seurat V5 Layer Support:
#' For Seurat V5 objects with multiple layers, you can selectively load layers per assay.
#' For example, use \code{assays = list(RNA = "data")} to load only the normalized expression layer,
#' skipping raw counts and scaled data.
#'
#' @seealso
#' \code{\link{SaveH5Seurat}} to save a Seurat object to h5Seurat format
#' \code{\link{Convert}} to convert to other formats
#'
#' @examples
#' \dontrun{
#' library(SeuratDisk)
#'
#' # Load entire h5Seurat file
#' seurat_obj <- LoadH5Seurat("data.h5seurat")
#'
#' # Load only specific assays
#' seurat_obj <- LoadH5Seurat("data.h5seurat", assays = c("RNA", "ADT"))
#'
#' # Load only specific data layers (memory-efficient for large files)
#' seurat_obj <- LoadH5Seurat("data.h5seurat", assays = c("data"))  # Only normalized expression
#'
#' # Load specific assays with different layers
#' seurat_obj <- LoadH5Seurat(
#'   "data.h5seurat",
#'   assays = list(RNA = c("data", "scale.data"), ADT = "data")
#' )
#'
#' # Load without reductions (faster)
#' seurat_obj <- LoadH5Seurat("data.h5seurat", reductions = FALSE)
#'
#' # Load UMAP and PCA reductions only
#' seurat_obj <- LoadH5Seurat("data.h5seurat", reductions = c("umap", "pca"))
#'
#' # Load spatial data without graphs (for Visium experiments)
#' seurat_obj <- LoadH5Seurat("visium.h5seurat", images = TRUE, graphs = FALSE)
#' }
#'
#' @export
#'
LoadH5Seurat <- function(file, ...) {
  UseMethod(generic = 'LoadH5Seurat', object = file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat character
#' @export
#'
LoadH5Seurat.character <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  hfile <- h5Seurat$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(LoadH5Seurat(
    file = hfile,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    neighbors = neighbors,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat H5File
#' @export
#'
LoadH5Seurat.H5File <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  return(LoadH5Seurat(
    file = as.h5Seurat(x = file),
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    neighbors = neighbors,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @importFrom Seurat as.Seurat
#'
#' @rdname LoadH5Seurat
#' @method LoadH5Seurat h5Seurat
#' @export
#'
LoadH5Seurat.h5Seurat <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  return(as.Seurat(
    x = file,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    neighbors = neighbors,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @importClassesFrom Seurat Seurat
#' @importFrom methods slot<-
#' @importFrom Seurat as.Seurat DefaultAssay Cells
#' Idents<- Idents Project<- Project
#' AddMetaData
#'
#' @aliases as.Seurat
#'
#' @rdname LoadH5Seurat
#' @method as.Seurat h5Seurat
#' @export
#'
as.Seurat.h5Seurat <- function(
  x,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = TRUE,
  tools = TRUE,
  verbose = TRUE,
  ...
) {
  index <- x$index()
  obj.all <- all(vapply(
    X = c(assays, reductions, graphs),
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  ))
  # Load Assays
  assays <- GetAssays(assays = assays, index = index)
  if (!DefaultAssay(object = index) %in% names(x = assays)) {
    active.assay <- names(x = assays)[1]
    warning(
      "Default assay not requested, using ",
      active.assay,
      " instead",
      call. = FALSE,
      immediate. = TRUE
    )
  } else {
    active.assay <- DefaultAssay(object = index)
  }
  assay.objects <- vector(mode = 'list', length = length(x = assays))
  names(x = assay.objects) <- names(x = assays)
  for (assay in names(x = assays)) {
    assay.objects[[assay]] <- AssembleAssay(
      assay = assay,
      file = x,
      slots = assays[[assay]],
      verbose = verbose
    )
  }
  default.assay <- list(assay.objects[[active.assay]])
  names(x = default.assay) <- active.assay
  object <- new(
    Class = 'Seurat',
    assays = default.assay,
    active.assay = active.assay,
    meta.data = data.frame(row.names = Cells(x = x)),
    version = package_version(x = x$version())
  )
  for (assay in names(x = assay.objects)) {
    if (assay != active.assay) {
      object[[assay]] <- assay.objects[[assay]]
    }
  }
  # Load DimReducs
  reductions <- GetDimReducs(
    reductions = reductions,
    index = index,
    assays = assays
  )
  for (reduc in reductions) {
    if (verbose) {
      message("Adding reduction ", reduc)
    }
    reduction <- AssembleDimReduc(
      reduction = reduc,
      file = x,
      verbose = verbose
    )
    if (isTRUE(x = getOption(x = 'SeuratDisk.dimreducs.allglobal', default = FALSE))) {
      slot(object = reduction, name = 'global') <- TRUE
    }
    object[[reduc]] <- reduction
  }
  # Load Graphs
  graphs <- GetGraphs(graphs = graphs, index = index, assays = assays)
  for (graph in graphs) {
    if (verbose) {
      message("Adding graph ", graph)
    }
    object[[graph]] <- AssembleGraph(graph = graph, file = x, verbose = verbose)
  }
  # Load Neighbors
  neighbors <- GetNeighbors(neighbors = neighbors, index = index)
  for (neighbor in neighbors) {
    if (verbose) {
      message("Adding neighbors ", neighbor)
    }
    object[[neighbor]] <- AssembleNeighbor(
      neighbor = neighbor,
      file = x,
      verbose = verbose
      )
  }
  # Load SpatialImages
  if (packageVersion(pkg = 'Seurat') >= numeric_version(x = spatial.version)) {
    images <- GetImages(images = images, index = index, assays = assays)
    for (image in images) {
      if (verbose) {
        message("Adding image ", image)
      }
      tryCatch({
        img_obj <- AssembleImage(
          image = image,
          file = x,
          verbose = verbose
        )
        object[[image]] <- img_obj
      }, error = function(e) {
        if (verbose) {
          message("Could not add image ", image, " to Seurat object: ", conditionMessage(e))
          message("Storing image in misc slot instead")
        }
        # Store in misc slot as fallback
        tryCatch({
          img_obj <- AssembleImage(
            image = image,
            file = x,
            verbose = verbose
          )
          slot(object = object, name = 'misc')[[paste0('image_', image)]] <<- img_obj
        }, error = function(e2) {
          if (verbose) {
            message("Failed to store image in misc: ", conditionMessage(e2))
          }
        })
      })
    }
  }
  # Load SeuratCommands
  if (commands) {
    if (verbose) {
      message("Adding command information")
    }
    cmds <- GetCommands(index = index, assays = assays)
    cmdlogs <- vector(mode = 'list', length = length(x = cmds))
    names(x = cmdlogs) <- cmds
    for (cmd in cmds) {
      cmdlogs[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = x,
        verbose = verbose
      )
    }
    slot(object = object, name = 'commands') <- cmdlogs
  }
  # Load meta.data
  if (meta.data && x$exists('meta.data')) {
    if (verbose) {
      message("Adding cell-level metadata")
    }

    # Safely read metadata with error handling
    tryCatch({
      md <- as.data.frame(x = x[['meta.data']], row.names = Cells(x = x))

      # Validate columns before subsetting
      if (ncol(x = md) > 0) {
        # Remove any columns that might cause issues
        problematic_cols <- sapply(md, function(col) {
          inherits(col, "environment") || is.null(col)
        })
        if (any(problematic_cols)) {
          md <- md[, !problematic_cols, drop = FALSE]
        }

        if (ncol(md) > 0) {
          object <- AddMetaData(object = object, metadata = md)
        }
      }
    }, error = function(e) {
      if (verbose) {
        warning("Could not load metadata: ", conditionMessage(e),
                call. = FALSE, immediate. = TRUE)
      }
    })
  }
  # Set cell identities and object project
  Idents(object = object) <- Idents(object = x)
  Project(object = object) <- Project(object = x)
  # Load misc
  if (misc) {
    if (verbose) {
      message("Adding miscellaneous information")
    }
    # Use SafeH5GroupToList to handle 3D+ arrays
    slot(object = object, name = 'misc') <- SafeH5GroupToList(h5obj = x[['misc']], recursive = TRUE)
  }
  # Load tools
  if (tools) {
    if (verbose) {
      message("Adding tool-specific results")
    }
    # Use SafeH5GroupToList to handle 3D+ arrays
    slot(object = object, name = 'tools') <- SafeH5GroupToList(h5obj = x[['tools']], recursive = TRUE)
  }
  # Load no.assay information
  if (obj.all && !is.null(x = index$no.assay)) {
    if (verbose) {
      message("Adding data that was not associated with an assay")
    }
    for (graph in index$no.assay$graphs) {
      object[[graph]] <- AssembleGraph(
        graph = graph,
        file = x,
        verbose = verbose
      )
    }
    for (cmd in index$no.assay$commands) {
      object[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = x,
        verbose = verbose
      )
    }
  }
  return(object)
}
