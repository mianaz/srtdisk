#' @include scdisk.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Version constants should now be defined in zzz.R for global availability

#' A class for connections to h5Seurat files
#'
#' @docType class
#' @name h5Seurat-class
#' @rdname h5Seurat-class
#' @aliases h5Seurat
#' @format An \code{\link[R6]{R6Class}} object
#' @seealso \code{\link[hdf5r]{H5File}}
#'
#' @importFrom R6 R6Class
#' @importFrom hdf5r H5File h5attr
#' @importFrom utils packageVersion
#'
#' @export
#'
h5Seurat <- R6Class(
  classname = 'h5Seurat',
  inherit = scdisk,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
    #' @description Get the index for this h5Seurat file
    index = function() {
      if (!length(x = private$index.internal)) {
        private$build.index(
          version = ClosestVersion(
            query = self$version(),
            targets = private$versions
          )
        )
      }
      return(private$index.internal)
    },
    #' @description Set the version attribute
    #' @param version A version number matching the regex
    #' \code{^\\d+(\\.\\d+){2}(\\.9\\d{3})?$}
    set.version = function(version) {
      version <- as.character(x = version)
      if (!grepl(pattern = version.regex, x = version)) {
        stop("Invalid version specification: ", version, call. = FALSE)
      }
      ClosestVersion(query = version, targets = private$versions)
      if (self$attr_exists(attr_name = 'version')) {
        self$attr_delete(attr_name = 'version')
      }
      self$create_attr(
        attr_name = 'version',
        robj = version,
        dtype = GuessDType(x = version)
      )
      if (!self$mode %in% modes$new) {
        private$validate()
      }
      return(invisible(x = self))
    },
    #' @description Get the version attribute
    version = function() {
      return(h5attr(x = self, which = 'version'))
    },
    #' @description Detect the Seurat version of the stored object
    #' @return A string indicating the detected version ("V3", "V4", or "V5")
    detect.version = function() {
      # Check if we have V5-specific structure
      if (self$exists("assays")) {
        # Get first assay for inspection
        assays <- names(self[["assays"]])
        if (length(assays) > 0) {
          first_assay <- self[["assays"]][[assays[1]]]

          # V5 indicator: layers subdirectory with sparse matrix components
          if (first_assay$exists("layers")) {
            layers_group <- first_assay[["layers"]]
            if (inherits(layers_group, "H5Group")) {
              # Check if layers contain sparse matrix structures
              layer_names <- layers_group$names
              if (length(layer_names) > 0) {
                first_layer <- layers_group[[layer_names[1]]]
                if (inherits(first_layer, "H5Group") &&
                    first_layer$exists("data") &&
                    first_layer$exists("indices") &&
                    first_layer$exists("indptr")) {
                  return("V5")
                }
              }
            }
          }

          # V4 indicator: direct counts/data without layers
          if (first_assay$exists("counts") || first_assay$exists("data")) {
            # Check if these are direct datasets (not groups)
            if (first_assay$exists("counts")) {
              counts_obj <- first_assay[["counts"]]
              if (!inherits(counts_obj, "H5Group")) {
                return("V4")
              }
            }
          }
        }
      }

      # Fall back to version attribute if structure is ambiguous
      version_attr <- self$version()
      if (!is.null(version_attr)) {
        version_num <- numeric_version(version_attr)
        if (version_num >= numeric_version("5.0.0")) {
          return("V5")
        } else if (version_num >= numeric_version("4.0.0")) {
          return("V4")
        } else {
          return("V3")
        }
      }

      # Default to V4 if we can't determine
      return("V4")
    },
    #' @description Check if this is a V5 format file
    #' @return Logical indicating if file uses V5 format
    is.v5 = function() {
      return(self$detect.version() == "V5")
    }
  ),
  private = list(
    # Fields
    index.internal = list(),
    versions = c('3.1.2', spatial.version, v5.version, '5.2.1'),
    # Methods
    build.index = function(version) {
      version <- match.arg(arg = version, choices = private$versions)
      version <- numeric_version(x = version)
      # Get Assay information
      index <- sapply(
        X = names(x = self[['assays']]),
        FUN = function(x) {
          # Handle V5 structure if applicable
          if (version >= numeric_version(x = v5.version) && 
              self[['assays']][[x]]$attr_exists(attr_name = 's4class')) {
            # Check for V5 Assay structure
            assay_class <- self[['assays']][[x]]$attr_open(attr_name = 's4class')$read()
            
            # Check if we have a V5-style "layers" subdirectory
            potential_layers <- c()
            if (self[['assays']][[x]]$exists("layers") && 
                inherits(self[['assays']][[x]][['layers']], "H5Group")) {
              # V5 structure with separate layers group
              potential_layers <- self[['assays']][[x]][['layers']]$names
            } else {
              # Legacy structure, look for direct children
              all_datasets <- self[['assays']][[x]]$names
              non_layer_names <- c("features", "scaled.features", "variable.features", "meta.features", "misc", "layers")
              potential_layers <- setdiff(all_datasets, non_layer_names)
            }
            
            # Create a check vector for all potential layers
            check <- rep(TRUE, length(potential_layers))
            names(check) <- potential_layers
            
            # Ensure standard slots are included for backward compatibility
            std_slots <- c('counts', 'data', 'scale.data')
            for (slot in std_slots) {
              if (!(slot %in% names(check))) {
                check[slot] <- FALSE
              }
            }
            
            # Special handling for scale.data
            if ('scale.data' %in% names(check)) {
              check[['scale.data']] <- check[['scale.data']] && 
                (self[['assays']][[x]]$exists(name = 'scaled.features') || 
                 (self[['assays']][[x]]$exists(name = 'scale.data') && 
                  self[['assays']][[x]][['scale.data']]$attr_exists(attr_name = 'scaled.features')))
            }
            
            # Add V5-specific information
            assay_result <- list(
              slots = check,
              v5class = assay_class,
              layers = potential_layers  # Include all potential layers
            )
            return(assay_result)
          } else {
            # Handle pre-V5 structure
            slots <- c('counts', 'data', 'scale.data')
            check <- slots %in% names(x = self[['assays']][[x]])
            names(x = check) <- slots
            check[['scale.data']] <- check[['scale.data']] && self[['assays']][[x]]$exists(name = 'scaled.features')
            check <- list(check)
            names(x = check) <- 'slots'
            return(check)
          }
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
      # Get DimReduc information
      reduc.slots <- c(
        'cell.embeddings',
        'feature.loadings',
        'feature.loadings.projected',
        'jackstraw'
      )
      for (reduc in names(x = self[['reductions']])) {
        reduc.assay <- h5attr(
          x = self[['reductions']][[reduc]],
          which = 'active.assay'
        )
        if (!reduc.assay %in% names(x = index)) {
          warning(
            "Cannot find assay ",
            reduc.assay,
            " in the H5Seurat file",
            call. = FALSE,
            immediate. = TRUE
          )
          next
        }
        check <- reduc.slots %in% names(x = self[['reductions']][[reduc]])
        names(x = check) <- reduc.slots
        if (check[['feature.loadings']]) {
          check[['feature.loadings']] <- self[['reductions']][[reduc]]$exists(name = 'features')
        }
        if (check[['feature.loadings.projected']]) {
          check[['feature.loadings.projected']] <- self[['reductions']][[reduc]]$exists(name = 'projected.features')
        }
        index[[reduc.assay]][['reductions']][[reduc]] <- check
        if (IsGlobal(object = self[['reductions']][[reduc]])) {
          index$global$reductions <- c(index$global$reductions, reduc)
        }
      }
      # Get graph information
      for (graph in names(x = self[['graphs']])) {
        if (self[['graphs']][[graph]]$attr_exists(attr_name = 'assay.used')) {
          graph.assay <- h5attr(x = self[['graphs']][[graph]], which = 'assay.used')
          if (graph.assay %in% names(x = index)) {
            index[[graph.assay]]$graphs <- c(index[[graph.assay]]$graphs, graph)
          } else {
            warning(
              "Cannot find assay ",
              graph.assay,
              " in the h5Seurat file",
              call. = FALSE,
              immediate. = TRUE
            )
          }
        } else {
          index$no.assay$graphs <- c(index$no.assay$graphs, graph)
        }
      }
      # Get images
      if (version >= numeric_version(x = spatial.version)) {
        for (image in names(x = self[['images']])) {
          img.group <- self[['images']][[image]]
          if (!img.group$attr_exists(attr_name = 'assay') || !img.group$attr_exists(attr_name = 's4class')) {
            next
          }
          img.assay <- h5attr(x = img.group, which = 'assay')
          if (!img.assay %in% names(x = index)) {
            warning(
              "Cannot find assay ",
              img.assay,
              " in the H5Seurat file",
              call. = FALSE,
              immediate. = TRUE
            )
            index$no.assay$images <- c(index$no.assay$images, image)
          } else {
            index[[img.assay]]$images <- c(index[[img.assay]]$images, image)
          }
          if (IsGlobal(object = img.group)) {
            index$global$images <- c(index$global$images, image)
          }
        }
      }
      # Get commands
      for (cmd in names(x = self[['commands']])) {
        assay <- ifelse(
          test = self[['commands']][[cmd]]$attr_exists(attr_name = 'assay.used'),
          yes = h5attr(x = self[['commands']][[cmd]], which = 'assay.used'),
          no = NA_character_
        )
        if (assay %in% setdiff(x = names(x = index), y = c('global', 'no.assay'))) {
          index[[assay]]$commands <- c(index[[assay]]$commands, cmd)
        } else if (!is.na(x = assay)) {
          warning(
            "Cannot find assay",
            assay,
            " in the h5Seurat file",
            call. = FALSE,
            immediate. = TRUE
          )
        } else {
          index$no.assay$commands <- c(index$no.assay$commands, cmd)
        }
      }
      # Get neighbors information
      if (Exists(x = self, name = 'neighbors')) {
        index$global$neighbors <- names(x = self[['neighbors']])
      }
      # TODO: Get metadata
      # TODO: Get miscellaneous data
      # TODO: Get tool-specific results
      # Finalize the index
      private$index.internal <- structure(
        .Data = index,
        class = c('h5SI', 'list'),
        active.assay = DefaultAssay(object = self)
      )
      return(invisible(x = NULL))
    },
    create = function(version, verbose = TRUE) {
      if (self$mode == 'r') {
        stop(private$errors(type = 'mode'), call. = FALSE)
      }
      version <- ClosestVersion(query = version, targets = private$versions)
      if (verbose) {
        message("Creating h5Seurat file for version ", version)
      }
      self$set.version(version = version)
      if (numeric_version(x = version) >= numeric_version(x = '3.1.2')) {
        for (group in c('assays', 'commands', 'neighbors', 'graphs', 'misc', 'reductions', 'tools')) {
          if (!private$is.data(name = group, type = 'H5Group')) {
            self$create_group(name = group)
          }
        }
        attrs <- c(
          'active.assay' = '',
          'project' = 'SeuratDiskProject'
        )
        for (i in seq_along(along.with = attrs)) {
          if (!self$attr_exists(attr_name = names(x = attrs)[i])) {
            self$create_attr(
              attr_name = names(x = attrs)[i],
              robj = attrs[i],
              dtype = GuessDType(x = attrs[i])
            )
          }
        }
      }
      if (numeric_version(x = version) >= numeric_version(x = spatial.version)) {
        self$create_group(name = 'images')
      }
      if (numeric_version(x = version) >= numeric_version(x = v5.version)) {
        # V5 specific initialization
        # Currently, we're keeping the same basic structure but can add
        # V5-specific requirements here as needed
        if (verbose) {
          message("Creating h5Seurat file for Seurat V5")
        }
      }
      return(invisible(x = self))
    },
    validate = function(verbose = TRUE, ...) {
      if (self$mode %in% modes$new) {
        private$create(
          version = packageVersion(pkg = 'Seurat'),
          verbose = verbose
        )
        return(invisible(x = NULL))
      }
      if (verbose) {
        message("Validating h5Seurat file")
      }
      if (!self$attr_exists(attr_name = 'version')) {
        stop(
          "Invalid h5Seurat file: cannot find attribute 'version'",
          call. = FALSE
        )
      }
      version <- h5attr(x = self, which = 'version')
      version <- ClosestVersion(query = version, targets = private$versions)
      version <- numeric_version(x = version)
      if (version >= numeric_version(x = '3.1.2')) {
        private$v3.1.2()
      }
      if (version >= numeric_version(x = '3.1.3.9900')) {
        private$v3.2.0()
      }
      if (version >= numeric_version(x = v5.version)) {
        private$v5.0.0()
      }
      private$build.index(version = as.character(x = version))
      return(invisible(x = NULL))
    },
    v3.1.2 = function() {
      # TODO: Check top-level attributes
      attrs <- c('project', 'active.assay', 'version')
      for (attr in attrs) {
        if (!self$attr_exists(attr_name = attr)) {
          stop("Missing attribute ", attr, call. = FALSE)
        }
      }
      # TODO: Check cell.names and meta.data
      if (!private$is.data(name = 'cell.names')) {
        stop("Cannot find dataset with cell names", call. = FALSE)
      }
      ncells <- self[['cell.names']]$dims
      if (length(x = ncells) != 1) {
        stop("Cell names must be one-dimensional", call. = FALSE)
      }
      if (private$is.data(name = 'meta.data')) {
        if (length(x = self[['meta.data']]$dims) != 1) {
          stop("Cell-level metadata must be one-dimensional")
        } else if (self[['meta.data']]$dims != ncells) {
          stop(
            "Cell number mismatch between cell names and cell-level metadata",
            call. = FALSE
          )
        }
        if (!inherits(x = self[['meta.data']]$get_type(), what = 'H5T_COMPOUND')) {
          stop("Cell-level metadata must be a data frame", call. = FALSE)
        }
      } else if (private$is.data(name = 'meta.data', type = 'H5Group')) {
        # warning("Validation for group meta data not yet implemented", call. = FALSE, immediate. = TRUE)
        ''
      } else {
        stop("Cannot find cell-level metadata")
      }
      # TODO: Check Assays
      if (!private$is.data(name = 'assays', type = 'H5Group')) {
        stop("Cannot find assay expression data", call. = FALSE)
      }
      if (!DefaultAssay(object = self) %in% names(x = self[['assays']])) {
        stop("Default assay not present", call. = FALSE)
      }
      for (assay in names(x = self[['assays']])) {
        if (!private$is.data(name = file.path('assays', assay), type = 'H5Group')) {
          stop(
            "Assay representations must be HDF5 groups, offending entry: ",
            assay,
            call. = FALSE
          )
        }
      }
      # TODO: Check DimReducs
      # TODO: Check Graphs
      # TODO: Check SeuratCommands
      # TODO: Check miscellaneous data
      # TODO: Check tool-specific results
      return(invisible(x = NULL))
    },
    v3.2.0 = function() {
      return(invisible(x = NULL))
      .NotYetImplemented()
    },
    v5.0.0 = function() {
      # Validate V5-specific structures
      # For now, we'll keep basic validation and build upon it as we implement more V5 features
      
      # Check for required groups for V5
      req.groups <- c('assays', 'reductions', 'graphs', 'misc', 'tools')
      for (group in req.groups) {
        if (!private$is.data(name = group, type = 'H5Group')) {
          stop("Missing required group ", group, " for Seurat V5", call. = FALSE)
        }
      }
      
      # Additional V5-specific validation will be added as we implement more features
      
      return(invisible(x = NULL))
    }
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tools for handling h5Seurat indexes (h5SI objects)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Tools for handling h5Seurat indexes
#'
#' @param x,object An h5Seurat index (\code{h5SI})
#'
#' @name h5SI
#' @rdname h5SI
#'
#' @seealso \code{\link[Seurat]{DefaultAssay}} \code{\link[base]{print}}
#'
#' @keywords internal
#'
NULL

#' @importFrom Seurat DefaultAssay
#'
#' @inheritParams Seurat::DefaultAssay
#'
#' @return \code{DefaultAssay}: the assay set as the default assay
#'
#' @rdname h5SI
#' @method DefaultAssay h5SI
#' @export
#'
DefaultAssay.h5SI <- function(object) {
  return(attr(x = object, which = 'active.assay'))
}

#' @return \code{print}: Invisibly returns \code{x}
#'
#' @importFrom cli symbol
#' @importFrom tools toTitleCase
#' @importFrom crayon red green yellow col_align
#'
#' @rdname h5SI
#' @method print h5SI
#' @export
#'
print.h5SI <- function(x, ...) {
  # Some constants
  catn <- function(...) {
    cat(..., '\n', sep = '')
  }
  reduc.header <- c(
    'Embeddings' = 'cell.embeddings',
    'Loadings' = 'feature.loadings',
    'Projected' = 'feature.loadings.projected',
    'JackStraw' = 'jackstraw'
  )
  symbols <- c(red(symbol$cross), green(symbol$tick))
  # Get the assays
  assays <- setdiff(x = names(x = x), y = c('global', 'no.assay'))
  assays <- assays[order(assays == DefaultAssay(object = x), decreasing = TRUE)]
  for (assay in assays) {
    header <- paste("Data for assay", assay)
    if (assay == DefaultAssay(object = x)) {
      header <- paste0(header, yellow(symbol$star), ' (default assay)')
    }
    catn(header)
    # Show slot information
    catn(col_align(
      text = c('counts', 'data', 'scale.data'),
      width = nchar(x = 'scale.data') + 1,
      align = 'center'
    ))
    catn(col_align(
      text = symbols[x[[assay]]$slots + 1],
      width = nchar(x = 'scale.data') + 1,
      align = 'center'
    ))
    # Show dimensional reduction information
    if (!is.null(x = x[[assay]]$reductions)) {
      catn("Dimensional reductions:")
      reductions <- names(x = x[[assay]]$reductions)
      reductions <- paste0(' ', reductions, ': ')
      reductions <- col_align(text = reductions, width = max(nchar(x = reductions)))
      catn(
        MakeSpace(n = max(nchar(x = reductions))),
        col_align(
          text = names(x = reduc.header),
          width = max(nchar(x = names(x = reduc.header))) + 1,
          align = 'center'
        )
      )
      for (i in seq_along(along.with = reductions)) {
        reduc <- names(x = x[[assay]]$reductions)[i]
        catn(
          reductions[i],
          col_align(
            text = symbols[x[[assay]]$reductions[[reduc]] + 1],
            width = max(nchar(x = names(x = reduc.header))) + 1,
            align = 'center'
          )
        )
      }
    }
    # Show graph information
    if (!is.null(x = x[[assay]]$graphs)) {
      catn("Graphs:")
      catn(paste0(' ', symbol$line, ' ', x[[assay]]$graphs, collapse = '\n'))
    }
    # Show image information
    if (!is.null(x = x[[assay]]$images)) {
      catn("Images:")
      catn(paste0(' ', symbol$line, ' ', x[[assay]]$images, collapse = '\n'))
    }
    # TODO: Show command information
    # TODO: Show globals
    # if (!is.null(x = x$global)) {
    #   catn("Globally available information:")
    # }
    # Show no assay
  }
  return(invisible(x = x))
}
