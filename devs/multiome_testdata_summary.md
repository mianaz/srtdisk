# Multiome Test Data Setup Summary

**Created**: 2025-11-12
**Purpose**: Demonstrate multimodal (RNA + ATAC) h5mu conversion

---

## What Was Created

### 1. Combined Multiome Dataset

**File**: `inst/testdata/pbmc_multiome_1k.rds`

- **Size**: 52 MB (manageable for testing)
- **Cells**: 1,000 (subsampled from 11,909)
- **Modalities**:
  - RNA: 36,601 genes
  - ATAC: 108,377 peaks
- **Metadata**: 6 columns including cell type annotations

**Source**:
- Combined from `pbmc_multiome_rna.rds` + `pbmc_multiome_atac.rds`
- Subset to 1,000 common cells (seed = 42 for reproducibility)
- Includes metadata from both modalities

### 2. Preparation Script

**File**: `inst/scripts/prepare_multiome_testdata.R`

**Features**:
- Loads RNA and ATAC data separately
- Finds common cells
- Subsamples to 1,000 cells
- Combines into single Seurat object
- Merges metadata from both modalities
- Adds ATAC reductions with prefix

**Usage**:
```r
source("inst/scripts/prepare_multiome_testdata.R")
```

### 3. Test Suite

**File**: `tests/testthat/test-multiome-h5mu.R`

**Tests**:
1. ✅ Save multiome to h5mu format
2. ✅ Verify h5mu file structure
3. ✅ Load h5mu back to Seurat (roundtrip)
4. ✅ Custom modality names
5. ✅ Selective assay saving (e.g., RNA only)

### 4. Example Script

**File**: `examples/multiome_conversion_example.R`

**Demonstrates**:
- Loading multiome data
- Saving to h5mu with default settings
- Custom modality naming
- Selective modality saving
- Reloading into Seurat
- Python interoperability workflow
- h5mu file structure inspection

### 5. Documentation

**Updated**:
- `inst/testdata/README.md`: Added multiome dataset documentation

---

## Quick Start

### Create the Multiome Test Data

```r
# Run the preparation script
source("inst/scripts/prepare_multiome_testdata.R")
```

**Output**:
```
=== Multiome object summary ===
Cells: 1000
Assays: RNA, ATAC
  RNA : 36601 features
  ATAC : 108377 features
Metadata columns: 6
  orig.ident, nCount_RNA, nFeature_RNA, seurat_annotations, nCount_ATAC, nFeature_ATAC

✓ Multiome test data created successfully!
```

### Use the Data

```r
library(srtdisk)

# Load multiome object
multiome_obj <- readRDS("inst/testdata/pbmc_multiome_1k.rds")

# Save to h5mu format
SaveH5MU(multiome_obj, "test_multiome.h5mu")

# Load back
multiome_reloaded <- LoadH5MU("test_multiome.h5mu")
```

### Run Tests

```r
# Run multiome tests
testthat::test_file("tests/testthat/test-multiome-h5mu.R")

# Or run all tests
devtools::test()
```

### Run Example

```r
source("examples/multiome_conversion_example.R")
```

---

## File Sizes

| File | Size | Description |
|------|------|-------------|
| `pbmc_multiome_1k.rds` | 52 MB | Combined RNA+ATAC (1k cells) |
| `pbmc_multiome_rna.rds` | 105 MB | RNA only (11.9k cells) |
| `pbmc_multiome_atac.rds` | 481 MB | ATAC only (11.9k cells) |
| `pbmc3k_final.h5ad` | 39 MB | Single-cell RNA (2.7k cells) |

---

## h5mu Format Overview

**h5mu** is the standard format for multimodal single-cell data:

### Structure

```
.h5mu file
├── mod/                    # Modalities
│   ├── rna/               # RNA modality
│   │   ├── X              # Expression matrix
│   │   ├── obs/           # Cell metadata
│   │   ├── var/           # Gene metadata
│   │   ├── obsm/          # Reductions (PCA, UMAP)
│   │   └── obsp/          # Graphs
│   └── atac/              # ATAC modality
│       ├── X              # Peak matrix
│       ├── obs/           # Cell metadata
│       ├── var/           # Peak metadata
│       └── ...
├── obs/                   # Global cell metadata
└── var/                   # (optional) Global features
```

### Python Compatibility

h5mu files can be used in Python with `muon`:

```python
import muon as mu

# Load multiome data
mdata = mu.read_h5mu('test_multiome.h5mu')

# Access modalities
rna = mdata.mod['rna']
atac = mdata.mod['atac']

# Process in Python...
# Then save and reload in R
```

---

## Testing Checklist

- [x] ✅ Multiome data created (1000 cells, RNA + ATAC)
- [x] ✅ SaveH5MU() function works
- [x] ✅ LoadH5MU() function works
- [x] ✅ Roundtrip conversion preserves data
- [x] ✅ Custom modality naming works
- [x] ✅ Selective assay saving works
- [x] ✅ h5mu file structure is valid
- [x] ✅ Documentation updated
- [x] ✅ Example script created
- [x] ✅ Test suite created

---

## Key Functions

### SaveH5MU()

Save Seurat multimodal object to h5mu format:

```r
SaveH5MU(
  object,                  # Seurat object with multiple assays
  filename,                # Output .h5mu file
  assays = NULL,          # Which assays to save (NULL = all)
  modality.names = NULL,  # Custom naming (e.g., c(RNA = "rna"))
  overwrite = FALSE,
  verbose = TRUE
)
```

### LoadH5MU()

Load h5mu file into Seurat:

```r
LoadH5MU(
  filename,                # Input .h5mu file
  modalities = NULL,      # Which modalities to load (NULL = all)
  assay.names = NULL,     # Custom Seurat assay names
  verbose = TRUE
)
```

---

## Use Cases

### 1. CITE-seq (RNA + Protein)

```r
# RNA + ADT assays
SaveH5MU(
  citeseq_obj,
  "citeseq.h5mu",
  modality.names = c(RNA = "rna", ADT = "prot")
)
```

### 2. Multiome (RNA + ATAC)

```r
# RNA + ATAC assays
SaveH5MU(
  multiome_obj,
  "multiome.h5mu",
  modality.names = c(RNA = "rna", ATAC = "atac")
)
```

### 3. Spatial + RNA

```r
# Spatial coordinates in metadata
SaveH5MU(
  spatial_obj,
  "spatial.h5mu"
)
```

### 4. Triple-omics

```r
# RNA + ATAC + Protein
SaveH5MU(
  triple_obj,
  "triple.h5mu",
  modality.names = c(RNA = "rna", ATAC = "atac", ADT = "prot")
)
```

---

## Next Steps

1. **Validate with real Python workflow**:
   - Save h5mu in R
   - Process in Python with muon
   - Reload in R

2. **Add more test cases**:
   - CITE-seq data
   - Spatial + multimodal
   - Large-scale multiome (>10k cells)

3. **Performance optimization**:
   - Chunked writing for large matrices
   - Compression options

4. **Documentation**:
   - Add vignette for multimodal workflows
   - Python interop examples

---

## Summary

✅ **Complete multiome test infrastructure** is now in place:
- Real multimodal data (RNA + ATAC)
- Preparation scripts
- Comprehensive tests
- Usage examples
- Full documentation

The test data is ready for validating h5mu conversion functionality and demonstrating multimodal workflows! 🎉
