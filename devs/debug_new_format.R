library(srtdisk)
library(hdf5r)

# Create test h5ad with NEW categorical format (categories/codes)
temp_h5ad <- tempfile(fileext = ".h5ad")
temp_h5seurat <- tempfile(fileext = ".h5seurat")

h5ad <- H5File$new(temp_h5ad, mode = "w")

n_cells <- 100
n_genes <- 50
h5ad[["X"]] <- matrix(rpois(n_cells * n_genes, 5), nrow = n_cells, ncol = n_genes)

obs_group <- h5ad$create_group("obs")
obs_group[["_index"]] <- paste0("Cell_", seq_len(n_cells))

# Add categorical using NEW format (categories/codes)
cell_type_group <- obs_group$create_group("cell_type")
cell_types <- rep(c("T cell", "B cell", "Monocyte"), length.out = n_cells)
cell_type_codes <- as.integer(factor(cell_types)) - 1L

cell_type_group[["codes"]] <- cell_type_codes
cell_type_group[["categories"]] <- unique(cell_types)

cat("Created h5ad with NEW categorical format\n")
cat("obs/cell_type/codes:", paste(head(cell_type_codes, 10), collapse = ", "), "\n")
cat("obs/cell_type/categories:", paste(unique(cell_types), collapse = ", "), "\n")

obs_group[["n_genes"]] <- rpois(n_cells, 1000)

var_group <- h5ad$create_group("var")
var_group[["_index"]] <- paste0("Gene_", seq_len(n_genes))

uns_group <- h5ad$create_group("uns")

h5ad$close_all()

cat("\n=== Converting h5ad to h5seurat ===\n")
Convert(temp_h5ad, temp_h5seurat, overwrite = TRUE, verbose = FALSE)

cat("\n=== Checking h5seurat structure ===\n")
h5s <- H5File$new(temp_h5seurat, mode = "r")
cat("meta.data columns: ", paste(names(h5s[["meta.data"]]), collapse = ", "), "\n")

if (h5s[["meta.data"]]$exists("cell_type")) {
  cat("\n✓ cell_type exists in h5seurat!\n")
  ct <- h5s[["meta.data"]][["cell_type"]]
  cat("Class:", class(ct), "\n")
  if (inherits(ct, "H5Group")) {
    cat("Group contents:", paste(names(ct), collapse = ", "), "\n")
    if (ct$exists("levels")) {
      cat("Levels:", paste(ct[["levels"]][], collapse = ", "), "\n")
    }
    if (ct$exists("values")) {
      cat("Values (first 10):", paste(head(ct[["values"]][], 10), collapse = ", "), "\n")
    }
  }
} else {
  cat("\n❌ cell_type does NOT exist in h5seurat!\n")
  cat("This means NormalizeH5ADCategorical failed to process it.\n")
}

h5s$close_all()

cat("\n=== Loading into Seurat ===\n")
seurat_obj <- LoadH5Seurat(temp_h5seurat, verbose = FALSE)
cat("Seurat meta.data columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
cat("Has cell_type?", "cell_type" %in% colnames(seurat_obj@meta.data), "\n")

if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
  cat("✓ SUCCESS! cell_type values (first 10):", paste(head(seurat_obj@meta.data$cell_type, 10), collapse = ", "), "\n")
} else {
  cat("❌ FAILURE: cell_type was lost during loading\n")
}

unlink(c(temp_h5ad, temp_h5seurat))
