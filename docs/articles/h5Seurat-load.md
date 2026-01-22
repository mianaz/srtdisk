# Saving and Loading Data from an h5Seurat File

The h5Seurat file format, based on
[HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format), is
specifically designed for the storage and analysis of multi-modal
single-cell and spatially-resolved expression experiments, for example,
from CITE-seq or 10X Visium technologies. It holds all molecular
information and associated metadata, including (for example)
nearest-neighbor graphs, dimensional reduction information, spatial
coordinates and image data, and cluster labels.

This vignette serves as a guide to saving and loading `Seurat` objects
to h5Seurat files. **This implementation supports Seurat v5 with full
Assay5 compatibility**, automatically handling both legacy Assay and
modern Assay5 objects. For more details about h5Seurat files, please see
the h5Seurat-spec vignette.

## Saving a dataset

Saving a `Seurat` object to an h5Seurat file is a fairly painless
process. All assays, dimensional reductions, spatial images, and
nearest-neighbor graphs are automatically saved as well as extra
metadata such as miscellaneous data, command logs, or cell identity
classes from a `Seurat` object. To save a `Seurat` object, we need the
[Seurat](https://satijalab.org/seurat) and srtdisk R packages.

``` r

library(Seurat)
library(srtdisk)
```

For this vignette, we’ll use the `pbmc3k.final` dataset from
SeuratData - a fully processed PBMC dataset with clustering, dimensional
reductions, and nearest-neighbor graphs.

``` r

library(SeuratData)

# Install pbmc3k if not already installed
if (!"pbmc3k.final" %in% rownames(InstalledData())) {
  InstallData("pbmc3k")
}

# Load the processed pbmc3k dataset
data("pbmc3k.final", package = "pbmc3k.SeuratData")
pbmc <- UpdateSeuratObject(pbmc3k.final)
pbmc
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  3 layers present: counts, data, scale.data
#>  2 dimensional reductions calculated: pca, umap
```

As seen, we have a dataset with multiple components to it. Saving an
object is as simple as calling
[`SaveH5Seurat()`](https://mianaz.github.io/srtdisk/reference/SaveH5Seurat.md);
minimally, this function takes a `Seurat` object and nothing else.
Optional arguments are present for specifying a filename and whether or
not you want to overwrite a preexisting file.

``` r

SaveH5Seurat(pbmc, filename = "pbmc3k.h5Seurat", overwrite = TRUE)
```

This process is quick and results in a compact on-disk file.

``` r

size <- file.size("pbmc3k.h5Seurat")
print(structure(size, class = "object_size"), units = "Mb")
#> 52.1 Mb
```

h5Seurat files are HDF5-based and can be read in other languages, such
as Python, providing better interoperability than Rds files.

## Connecting to and querying h5Seurat files

Unlike most data formats, HDF5 files can be connected to and explored
without loading the data into memory. To facilitate this, we’ve built an
`h5Seurat` object to serve as an interface to h5Seurat files in R.
`h5Seurat` objects are built off of the `H5File` object from the hdf5r
package.

`h5Seurat` objects and R6 classes

One thing to note, `h5Seurat` and `H5File` objects are R6 objects.
Unlike most R objects (called S3 and S4), and more like objects in
Python, R6 objects are *encapsulated* objects; this means that methods
are attached directly to the object instead of to a generic function.  
In Seurat, most functions take an object as input and return an object
as output. These functions actually run differently depending on the
class of the object passed to them. For example, `RunUMAP` has 4
different modes of operation, depending on the type of object that’s
passed to it. One can see which objects trigger different routines by
using the
[`methods`](https://stat.ethz.ch/R-manual/R-patched/library/utils/html/methods.html)
function. Functions that change behavior are known as “generics” and the
exact implementations are known as “methods”; these methods are
associated with the generic instead of with the object itself.  
R6 objects, however, have their methods attached directly to the object.
Calling an R6 method is done similarly to data access in R’s S3 and S4
object system: using the `$` operator. For example, creating a new
Seurat object is done with `CreateSeuratObject` or
`new(Class = "Seurat")` (for advanced users), while initializing a new
`h5Seurat` object is done with `h5Seurat$new()`  
For more details about R6 objects, please see the [R6 website and
documentation](https://r6.r-lib.org/)

Connecting to an h5Seurat file is as simple as instantiating an
`h5Seurat` object.

``` r

hfile <- Connect("pbmc3k.h5Seurat")
hfile
#> Class: h5Seurat
#> Filename: /Users/miana/Documents/GitHub/srtdisk/vignettes/pbmc3k.h5Seurat
#> Access type: H5F_ACC_RDONLY
#> Attributes: version, project, active.assay
#> Listing:
#>          name    obj_type dataset.dims dataset.type_class
#>  active.ident   H5I_GROUP         <NA>               <NA>
#>        assays   H5I_GROUP         <NA>               <NA>
#>    cell.names H5I_DATASET         2638         H5T_STRING
#>      commands   H5I_GROUP         <NA>               <NA>
#>        graphs   H5I_GROUP         <NA>               <NA>
#>        images   H5I_GROUP         <NA>               <NA>
#>     meta.data   H5I_GROUP         <NA>               <NA>
#>          misc   H5I_GROUP         <NA>               <NA>
#>     neighbors   H5I_GROUP         <NA>               <NA>
#>    reductions   H5I_GROUP         <NA>               <NA>
#> < Printed 10, out of 11>
```

As seen, the h5Seurat file is structured similarly to a `Seurat` object,
with different HDF5 groups sharing the names of slots in a `Seurat`
object. However, it’s difficult to glean what data is present in this
dataset similar to calling a `Seurat` object in the R console. To get
around this, we’ve created an `index` method for `h5Seurat` objects;
this method creates a summary of the data stored within the h5Seurat
object. As `Seurat` objects are organized around the assay data, this
h5Seurat index showcases the data grouped by assay.

``` r

hfile$index()
#> Data for assay RNA★ (default assay)
#>    counts      data    scale.data
#>      ✔          ✔          ✔     
#> Dimensional reductions:
#>         Embeddings  Loadings  Projected  JackStraw 
#>  pca:       ✔          ✔          ✖          ✔     
#>  umap:      ✔          ✖          ✖          ✖     
#> Graphs:
#>  ─ RNA_nn
#>  ─ RNA_snn
```

First we get a breakdown of what slots are filled within each assay,
followed by a table of dimensional reduction information. This table
shows which bits of information (eg. cell embeddings, feature loadings,
JackStraw data) are present. these tables, we get a list of
nearest-neighbor graphs and spatial image data. This way, we can see
what data gets loaded on a per-assay basis as is required by Seurat.

To explore an h5Seurat file deeper, we can use the double bracket `[[`
operator to explore various aspects of the dataset. The double bracket
`[[` operator takes a UNIX-style path comprised of dataset names.

``` r

hfile[["assays"]]
#> Class: H5Group
#> Filename: /Users/miana/Documents/GitHub/srtdisk/vignettes/pbmc3k.h5Seurat
#> Group: /assays
#> Listing:
#>  name  obj_type dataset.dims dataset.type_class
#>   RNA H5I_GROUP         <NA>               <NA>
hfile[["assays/RNA"]]
#> Class: H5Group
#> Filename: /Users/miana/Documents/GitHub/srtdisk/vignettes/pbmc3k.h5Seurat
#> Group: /assays/RNA
#> Attributes: key
#> Listing:
#>               name    obj_type dataset.dims dataset.type_class
#>             counts   H5I_GROUP         <NA>               <NA>
#>               data   H5I_GROUP         <NA>               <NA>
#>           features H5I_DATASET        13714         H5T_STRING
#>      meta.features   H5I_GROUP         <NA>               <NA>
#>         scale.data H5I_DATASET 13714 x 2638          H5T_FLOAT
#>    scaled.features H5I_DATASET        13714         H5T_STRING
#>  variable.features H5I_DATASET         2000         H5T_STRING
hfile[["reductions"]]
#> Class: H5Group
#> Filename: /Users/miana/Documents/GitHub/srtdisk/vignettes/pbmc3k.h5Seurat
#> Group: /reductions
#> Listing:
#>  name  obj_type dataset.dims dataset.type_class
#>   pca H5I_GROUP         <NA>               <NA>
#>  umap H5I_GROUP         <NA>               <NA>
hfile[["reductions/umap"]]
#> Class: H5Group
#> Filename: /Users/miana/Documents/GitHub/srtdisk/vignettes/pbmc3k.h5Seurat
#> Group: /reductions/umap
#> Attributes: active.assay, key, global
#> Listing:
#>             name    obj_type dataset.dims dataset.type_class
#>  cell.embeddings H5I_DATASET     2638 x 2          H5T_FLOAT
#>             misc   H5I_GROUP         <NA>               <NA>
```

When finished exploring an h5Seurat file, remember to close the
connection. Because we’re working with file on disk directly rather than
loading it into memory, we need to close it to prevent file corruption.
You can also open the file in read-only mode (`mode = "r"`) to help
alleviate file corruption, though it’s still a good habit to close the
h5Seurat file when done working with it.

``` r

hfile$close_all()
```

## Loading datasets

Reading data from an h5Seurat file is as simple as calling
[`LoadH5Seurat()`](https://mianaz.github.io/srtdisk/reference/LoadH5Seurat.md);
by default, it loads the entire object into memory.

``` r

pbmc2 <- LoadH5Seurat("pbmc3k.h5Seurat")
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  1 layer present: counts
#>  2 dimensional reductions calculated: pca, umap
```

However, there are situations in which loading an entire `Seurat` object
is not desirable. As such, we can leverage the HDF5 format and load only
parts of a dataset at a time. `LoadH5Seurat` makes use of assay
association to limit the data loaded. In `Seurat` objects, all
dimensional reduction information, nearest-neighbor graphs, and spatial
image data have an assay they “belong” to (see
[`DefaultAssay()`](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html)
in the Seurat package for more details). If only certain assays are
requested, then only the object associated with those assays are loaded.

There are four main parameters for controlling data loading. The first
is the `assays` parameter; this parameter controls which assays are
loaded and which slots of each assay are loaded. The simplest level of
control is specifying the assays to load. For our dataset, passing an
assay name will load the entire assay object for the assay specified.

``` r

pbmc2 <- LoadH5Seurat("pbmc3k.h5Seurat", assays = "RNA")
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  1 layer present: counts
#>  2 dimensional reductions calculated: pca, umap
```

We can also choose the slots to load; the slots available are “counts”
for the raw expression data, “data” for the normalized expression data,
or “scale.data” for the scaled expression data. Specifying slots instead
of assays will load the desired slots from all assays that have the
requested slots. When specifying slots, one of either “counts” or “data”
**must** be specified as the `Seurat` object uses these slots to control
dataset dimensionality information.

``` r

pbmc2 <- LoadH5Seurat("pbmc3k.h5Seurat", assays = "data")
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  1 layer present: counts
#>  2 dimensional reductions calculated: pca, umap
```

For more fine-tuned control, the `assays` parameter can also take a
named list or vector, where the names are the names of the assays to
load and the values are the slots to load.

``` r

pbmc2 <- LoadH5Seurat("pbmc3k.h5Seurat", assays = list("RNA" = c("data", "scale.data")))
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  2 layers present: counts, scale.data
#>  2 dimensional reductions calculated: pca, umap
```

Finally, passing `NULL` to `assays` (the default behavior) loads all
assays and all slots.

The second of the main parameters is the `reductions` parameter; this
parameter controls which dimensional reductions are loaded. As
dimensional reductions are tied to assays, the data request needs to be
either associated to a loaded assay or marked as global to be loaded
(see details below).

``` r

pbmc2 <- LoadH5Seurat("pbmc3k.h5Seurat", assays = "RNA", reductions = "pca")
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  1 layer present: counts
#>  1 dimensional reduction calculated: pca
```

There are three special values the `reductions` parameter can take:
`NULL` for all dimensional reductions that can be loaded (the default
behavior), `NA` for *global* dimensional reductions only, or `FALSE` for
no dimensional reduction information.

The `graphs` parameter is the third main parameter; this parameter
controls which nearest-neighbor graphs to load. Just like dimensional
reduction information, nearest-neighbor graphs are tied to assays, and
thus are only loaded when their associated assay is loaded as well.
There are two special values the `graphs` parameter can take: `NULL` for
all graphs that can be loaded (the default behavior) or `FALSE` for no
nearest-neighbor graphs.

The final main parameter is the `images` parameter; this parameter
controls which spatial image data is loaded. All spatial image data are
marked global by default, so they are loaded whether or not their
associated assays are loaded as well. The `images` parameter has three
special values: `NULL` for all spatial image data (the default), `NA`
for *global* spatial image data (typically the same as `NULL`), or
`FALSE` for no spatial image data.

With these four parameters, there is a lot of customization for loading
`Seurat` objects from h5Seurat files. For example, the following will
load the “data” and “scale.data” slots from the “RNA” assay, *global*
dimensional reductions, none of the nearest neighbor graphs, and no
spatial images.

``` r

pbmc2 <- LoadH5Seurat("pbmc3k.h5Seurat", assays = list("RNA" = c("data", "scale.data")), reductions = NA, graphs = FALSE, images = FALSE)
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  2 layers present: counts, scale.data
#>  1 dimensional reduction calculated: umap
```

In addition, there are four secondary parameters to `LoadH5Seurat`:
`meta.data`, `commands`, `misc`, and `tools`; these all take simple
`TRUE`/`FALSE` values to control the loading of cell-level metadata,
command logs, miscellaneous information, or tool-specific results,
respectively.

Global objects

The concept of global objects in Seurat is designed as an extension to
the assay-centric nature of `Seurat` objects. In Seurat, each assay is
considered to be one experiment or measurement of data for a common
group of cells. These assays are then used to generate working
summaries, such as reduced dimension space or nearest-neighbor graphs.
Generally, if an assay is removed from an object, the working summaries
are of little use so they get removed as well.  
However, there are some instances in which these working summaries are
useful outside the context of their assay. For example, some reduced
representations of the data such as tSNE or UMAP are useful for
visualization regardless of the the assay. As such, certain reduced
representations and all spatial image data are marked as *global*,
allowing them to persist as useful visualization contexts without their
associated assay and large expression matrices being present.  
For more details about global objects, please see the documentation for
[globality in Seurat](https://rdrr.io/cran/Seurat/man/IsGlobal.html)

Partial loading of datasets is an excellent way to limit memory usage
and prevent the loading of massive datasets into memory. However, there
can be instances in which a partial dataset was loaded, but then needs
to be expanded with additional data from the h5Seurat file. Instead of
redoing the partial load, we can make use of
[`AppendData()`](https://mianaz.github.io/srtdisk/reference/AppendData.md)
to add additional objects from an h5Seurat file to an already-loaded
`Seurat` object. To show how this works, we’ll start off by loading just
the “data” slot from the “RNA” assay, but not load any dimensional
reduction information or nearest-neighbor graphs.

``` r

pbmc2 <- LoadH5Seurat("pbmc3k.h5Seurat", assays = c("RNA" = "data"), reductions = FALSE, graphs = FALSE, images = FALSE)
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  1 layer present: counts
```

`AppendData` takes the h5Seurat file, the `Seurat` object generated from
`LoadH5Seurat` and uses the four main paramters from `LoadH5Seurat`
(`assays`, `reductions`, `graphs`, and `images`). These parameters are
used in the same way as `LoadH5Seurat` with one exception: `assays` can
now take `FALSE` as a value. By passing `FALSE`, we prevent other assay
information from being loaded; this is useful if we only want to add
other bits of data to our `Seurat` object. For example, we can choose to
add only *global* dimensional reductions to the already existing
`Seurat` object.

``` r

pbmc2 <- AppendData("pbmc3k.h5Seurat", pbmc2, assays = FALSE, reductions = NA, graphs = FALSE, images = FALSE)
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  1 layer present: counts
#>  1 dimensional reduction calculated: umap
```

The only limits to the number of times `AppendData` can be run is when
the `h5Seurat` file has run out of data not present in the `Seurat`
object. Otherwise, it can be run multiple times, adding new bits of data
to our `Seurat` object. Here, we fill out the rest of the “RNA” assay,
but load no other information

``` r

pbmc2 <- AppendData("pbmc3k.h5Seurat", pbmc2, assays = "RNA", reductions = FALSE, graphs = FALSE, images = FALSE)
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  1 layer present: counts
#>  1 dimensional reduction calculated: umap
```

If we want to perform a “full append” (loading all bits of data of a
`Seurat` object from an h5Seurat file), we can set the four parameters
to `NULL`, which happens to be the default values for these parmaters.
This loads the rest of the `Seurat` object from the h5Seurat file into
memory.

``` r

pbmc2 <- AppendData("pbmc3k.h5Seurat", pbmc2)
pbmc2
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  1 layer present: counts
#>  2 dimensional reductions calculated: umap, pca
```

## Working with Spatial Data

For spatial transcriptomics data, we can use the stxBrain dataset from
SeuratData:

``` r

library(SeuratData)

# Install stxBrain if not already installed
if (!"stxBrain" %in% rownames(InstalledData())) {
  InstallData("stxBrain")
}

# Load anterior1 section
brain <- UpdateSeuratObject(LoadData("stxBrain", type = "anterior1"))
brain
#> An object of class Seurat 
#> 31053 features across 2696 samples within 1 assay 
#> Active assay: Spatial (31053 features, 0 variable features)
#>  1 layer present: counts
#>  1 spatial field of view present: anterior1
```

Save and reload to demonstrate h5Seurat with spatial data:

``` r

SaveH5Seurat(brain, filename = "stxBrain.h5Seurat", overwrite = TRUE)

# Reload
brain2 <- LoadH5Seurat("stxBrain.h5Seurat")

# Verify spatial images are preserved
cat("Original images:", Images(brain), "\n")
#> Original images: anterior1
cat("Reloaded images:", Images(brain2), "\n")
#> Reloaded images: anterior1
```

## Session Info

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.2
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/Indiana/Indianapolis
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] stxKidney.SeuratData_0.1.0    stxBrain.SeuratData_0.1.2    
#>  [3] ssHippo.SeuratData_3.1.4      pbmcref.SeuratData_1.0.0     
#>  [5] pbmcMultiome.SeuratData_0.1.4 pbmc3k.SeuratData_3.1.4      
#>  [7] panc8.SeuratData_3.0.2        cbmc.SeuratData_3.1.4        
#>  [9] SeuratData_0.2.2.9002         srtdisk_0.1.0                
#> [11] Seurat_5.4.0                  SeuratObject_5.3.0           
#> [13] sp_2.2-0                     
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.4        
#>   [4] spatstat.utils_3.2-1   farver_2.1.2           rmarkdown_2.30        
#>   [7] fs_1.6.6               ragg_1.5.0             vctrs_0.7.0           
#>  [10] ROCR_1.0-11            spatstat.explore_3.6-0 htmltools_0.5.9       
#>  [13] sass_0.4.10            sctransform_0.4.3      parallelly_1.46.1     
#>  [16] KernSmooth_2.23-26     bslib_0.9.0            htmlwidgets_1.6.4     
#>  [19] desc_1.4.3             ica_1.0-3              plyr_1.8.9            
#>  [22] plotly_4.11.0          zoo_1.8-15             cachem_1.1.0          
#>  [25] igraph_2.2.1           mime_0.13              lifecycle_1.0.5       
#>  [28] pkgconfig_2.0.3        Matrix_1.7-4           R6_2.6.1              
#>  [31] fastmap_1.2.0          fitdistrplus_1.2-4     future_1.69.0         
#>  [34] shiny_1.12.1           digest_0.6.39          patchwork_1.3.2       
#>  [37] tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.5.1         
#>  [40] textshaping_1.0.4      progressr_0.18.0       spatstat.sparse_3.1-0 
#>  [43] httr_1.4.7             polyclip_1.10-7        abind_1.4-8           
#>  [46] compiler_4.5.2         bit64_4.6.0-1          S7_0.2.1              
#>  [49] fastDummies_1.7.5      MASS_7.3-65            rappdirs_0.3.4        
#>  [52] tools_4.5.2            lmtest_0.9-40          otel_0.2.0            
#>  [55] httpuv_1.6.16          future.apply_1.20.1    goftest_1.2-3         
#>  [58] glue_1.8.0             nlme_3.1-168           promises_1.5.0        
#>  [61] grid_4.5.2             Rtsne_0.17             cluster_2.1.8.1       
#>  [64] reshape2_1.4.5         generics_0.1.4         hdf5r_1.3.12          
#>  [67] gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2           
#>  [70] data.table_1.18.0      spatstat.geom_3.7-0    RcppAnnoy_0.0.23      
#>  [73] ggrepel_0.9.6          RANN_2.6.2             pillar_1.11.1         
#>  [76] stringr_1.6.0          spam_2.11-3            RcppHNSW_0.6.0        
#>  [79] later_1.4.5            splines_4.5.2          dplyr_1.1.4           
#>  [82] lattice_0.22-7         survival_3.8-6         bit_4.6.0             
#>  [85] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2          
#>  [88] pbapply_1.7-4          knitr_1.51             gridExtra_2.3         
#>  [91] scattermore_1.2        xfun_0.56              matrixStats_1.5.0     
#>  [94] stringi_1.8.7          lazyeval_0.2.2         yaml_2.3.12           
#>  [97] evaluate_1.0.5         codetools_0.2-20       tibble_3.3.1          
#> [100] cli_3.6.5              uwot_0.2.4             xtable_1.8-4          
#> [103] reticulate_1.44.1      systemfonts_1.3.1      jquerylib_0.1.4       
#> [106] dichromat_2.0-0.1      Rcpp_1.1.1             globals_0.18.0        
#> [109] spatstat.random_3.4-4  png_0.1-8              spatstat.univar_3.1-6 
#> [112] parallel_4.5.2         pkgdown_2.2.0          ggplot2_4.0.1         
#> [115] dotCall64_1.2          listenv_0.10.0         viridisLite_0.4.2     
#> [118] scales_1.4.0           ggridges_0.5.7         purrr_1.2.1           
#> [121] crayon_1.5.3           rlang_1.1.7            cowplot_1.2.0
```
