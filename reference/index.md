# Package index

## Package Information

Overview and utilities

- [`srtdisk-package`](https://mianaz.github.io/srtdisk/reference/SeuratDisk-package.md)
  [`_PACKAGE`](https://mianaz.github.io/srtdisk/reference/SeuratDisk-package.md)
  [`srtdisk`](https://mianaz.github.io/srtdisk/reference/SeuratDisk-package.md)
  : srtdisk: Interfaces for HDF5-Based Single Cell File Formats

- [`scdisk-class`](https://mianaz.github.io/srtdisk/reference/scdisk-class.md)
  [`scdisk`](https://mianaz.github.io/srtdisk/reference/scdisk-class.md)
  : A disk-based object for single-cell analysis

- [`IsSCDisk()`](https://mianaz.github.io/srtdisk/reference/IsSCDisk.md)
  : Does an R6 class inherit from scdisk

- [`GetSCDisk()`](https://mianaz.github.io/srtdisk/reference/RegisterSCDisk.md)
  [`RegisterSCDisk()`](https://mianaz.github.io/srtdisk/reference/RegisterSCDisk.md)
  :

  Get and Register `scdisk` Subclasses

## h5Seurat Format

Working with h5Seurat files - the native HDF5-based R/Seurat format

- [`SaveH5Seurat()`](https://mianaz.github.io/srtdisk/reference/SaveH5Seurat.md)
  [`as.h5Seurat()`](https://mianaz.github.io/srtdisk/reference/SaveH5Seurat.md)
  :

  Save a `Seurat` object to an h5Seurat file

- [`LoadH5Seurat()`](https://mianaz.github.io/srtdisk/reference/LoadH5Seurat.md)
  [`as.Seurat(`*`<h5Seurat>`*`)`](https://mianaz.github.io/srtdisk/reference/LoadH5Seurat.md)
  :

  Load a saved `Seurat` object from an h5Seurat file

- [`h5Seurat-class`](https://mianaz.github.io/srtdisk/reference/h5Seurat-class.md)
  [`h5Seurat`](https://mianaz.github.io/srtdisk/reference/h5Seurat-class.md)
  : A class for connections to h5Seurat files

- [`Cells(`*`<h5Seurat>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  [`DefaultAssay(`*`<h5Seurat>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  [`` `DefaultAssay<-`( ``*`<h5Seurat>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  [`Idents(`*`<h5Seurat>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  [`IsGlobal(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  [`Key(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  [`Project(`*`<h5Seurat>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  [`` `Project<-`( ``*`<h5Seurat>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  [`Stdev(`*`<h5Seurat>`*`)`](https://mianaz.github.io/srtdisk/reference/h5Seurat-bindings.md)
  : Seurat bindings for h5Seurat files

- [`Connect()`](https://mianaz.github.io/srtdisk/reference/Connect.md) :
  Connect to a single-cell HDF5 dataset

- [`DefaultAssay(`*`<h5SI>`*`)`](https://mianaz.github.io/srtdisk/reference/h5SI.md)
  [`print(`*`<h5SI>`*`)`](https://mianaz.github.io/srtdisk/reference/h5SI.md)
  : Tools for handling h5Seurat indexes

## AnnData (h5ad) Format

Converting between Seurat and Python AnnData format

- [`H5ADToH5Seurat()`](https://mianaz.github.io/srtdisk/reference/H5ADToH5Seurat.md)
  : Convert AnnData/H5AD files to h5Seurat files
- [`H5SeuratToH5AD()`](https://mianaz.github.io/srtdisk/reference/H5SeuratToH5AD.md)
  : Convert h5Seurat files to H5AD files

## Loom Format

Loading and saving Loom files for interoperability with loompy

- [`LoadLoom()`](https://mianaz.github.io/srtdisk/reference/LoadLoom.md)
  [`as.Seurat(`*`<loom>`*`)`](https://mianaz.github.io/srtdisk/reference/LoadLoom.md)
  : Loom-file Loading

- [`SaveLoom()`](https://mianaz.github.io/srtdisk/reference/SaveLoom.md)
  [`as.loom()`](https://mianaz.github.io/srtdisk/reference/SaveLoom.md)
  :

  Save a
  [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
  object to a loom file

- [`loom-class`](https://mianaz.github.io/srtdisk/reference/loom-class.md)
  [`loom`](https://mianaz.github.io/srtdisk/reference/loom-class.md) : A
  class for connections to loom files

- [`DefaultAssay(`*`<loom>`*`)`](https://mianaz.github.io/srtdisk/reference/loom-bindings.md)
  [`dim(`*`<loom>`*`)`](https://mianaz.github.io/srtdisk/reference/loom-bindings.md)
  : Seurat binding for loom files

- [`LoadLoom0.1()`](https://mianaz.github.io/srtdisk/reference/LoomLoading.md)
  [`LoadLoom3.0()`](https://mianaz.github.io/srtdisk/reference/LoomLoading.md)
  : Loom-file Loading

- [`LoomValidate0.1()`](https://mianaz.github.io/srtdisk/reference/ValidateLoom.md)
  [`LoomValidate3.0.0()`](https://mianaz.github.io/srtdisk/reference/ValidateLoom.md)
  : Validate Loom Files

## Format Conversion

General conversion utilities between all supported formats

- [`Convert()`](https://mianaz.github.io/srtdisk/reference/Convert.md) :
  Convert an on-disk single-cell dataset to another format
- [`FileType()`](https://mianaz.github.io/srtdisk/reference/FileType.md)
  : Determine a filetype based on its extension

## Reading & Writing

Low-level HDF5 reading and writing functions

- [`as.array(`*`<H5D>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  [`as.data.frame(`*`<H5D>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  [`as.data.frame(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  [`as.factor(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  [`as.list(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  [`as.matrix(`*`<H5D>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  [`as.matrix(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  [`as.sparse(`*`<H5D>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  [`as.sparse(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/ReadH5.md)
  : Load data from an HDF5 File
- [`WriteH5Group()`](https://mianaz.github.io/srtdisk/reference/WriteH5Group.md)
  : Write data to an HDF5 group
- [`BasicWrite()`](https://mianaz.github.io/srtdisk/reference/BasicWrite.md)
  : Write lists and other data to an HDF5 dataset
- [`ImageWrite()`](https://mianaz.github.io/srtdisk/reference/ImageWrite.md)
  : Write a SpatialImage object to an HDF5 dataset
- [`SparseWrite()`](https://mianaz.github.io/srtdisk/reference/SparseWrite.md)
  : Write a sparse matrix to an HDF5 dataset
- [`WriteAttribute()`](https://mianaz.github.io/srtdisk/reference/WriteAttribute.md)
  : Write an attribute to an HDF5 file, group, or dataset
- [`AttrExists()`](https://mianaz.github.io/srtdisk/reference/H5Exists.md)
  [`Exists()`](https://mianaz.github.io/srtdisk/reference/H5Exists.md) :
  Check to see if a dataset, group, or attribute exists in an HDF5 file,
  group, or dataset
- [`H5Path()`](https://mianaz.github.io/srtdisk/reference/H5Path.md) :
  Create an HDF5 object path

## Object Operations

Utilities for working with HDF5 objects and attributes

- [`GetAssays()`](https://mianaz.github.io/srtdisk/reference/GetObject.md)
  [`GetCommands()`](https://mianaz.github.io/srtdisk/reference/GetObject.md)
  [`GetDimReducs()`](https://mianaz.github.io/srtdisk/reference/GetObject.md)
  [`GetGraphs()`](https://mianaz.github.io/srtdisk/reference/GetObject.md)
  [`GetImages()`](https://mianaz.github.io/srtdisk/reference/GetObject.md)
  [`GetNeighbors()`](https://mianaz.github.io/srtdisk/reference/GetObject.md)
  : Figure out which objects to load from an h5Seurat file

- [`AssembleAssay()`](https://mianaz.github.io/srtdisk/reference/AssembleObject.md)
  [`AssembleDimReduc()`](https://mianaz.github.io/srtdisk/reference/AssembleObject.md)
  [`AssembleGraph()`](https://mianaz.github.io/srtdisk/reference/AssembleObject.md)
  [`AssembleImage()`](https://mianaz.github.io/srtdisk/reference/AssembleObject.md)
  [`AssembleNeighbor()`](https://mianaz.github.io/srtdisk/reference/AssembleObject.md)
  [`AssembleSeuratCommand()`](https://mianaz.github.io/srtdisk/reference/AssembleObject.md)
  : Assemble an object from an h5Seurat file

- [`AppendData()`](https://mianaz.github.io/srtdisk/reference/AppendData.md)
  :

  Append data from an h5Seurat file to a preexisting
  [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
  object

- [`Transpose()`](https://mianaz.github.io/srtdisk/reference/Transpose.md)
  : Transpose a matrix

- [`PadMatrix()`](https://mianaz.github.io/srtdisk/reference/PadMatrix.md)
  : Pad a matrix

## Testing & Validation

Functions for testing and validating HDF5 files and Seurat objects

- [`IsDataFrame(`*`<H5D>`*`)`](https://mianaz.github.io/srtdisk/reference/TestH5.md)
  [`IsDataFrame(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/TestH5.md)
  [`IsFactor(`*`<H5D>`*`)`](https://mianaz.github.io/srtdisk/reference/TestH5.md)
  [`IsFactor(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/TestH5.md)
  [`IsList(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/TestH5.md)
  [`IsLogical(`*`<H5D>`*`)`](https://mianaz.github.io/srtdisk/reference/TestH5.md)
  [`IsMatrix(`*`<H5D>`*`)`](https://mianaz.github.io/srtdisk/reference/TestH5.md)
  [`IsMatrix(`*`<H5Group>`*`)`](https://mianaz.github.io/srtdisk/reference/TestH5.md)
  : Test HDF5 datasets and groups to see what kind of data they are
- [`IsDataFrame()`](https://mianaz.github.io/srtdisk/reference/TestObject.md)
  [`IsFactor()`](https://mianaz.github.io/srtdisk/reference/TestObject.md)
  [`IsList()`](https://mianaz.github.io/srtdisk/reference/TestObject.md)
  [`IsLogical()`](https://mianaz.github.io/srtdisk/reference/TestObject.md)
  [`IsMatrix()`](https://mianaz.github.io/srtdisk/reference/TestObject.md)
  : Test an object's class
- [`CheckMatrix()`](https://mianaz.github.io/srtdisk/reference/CheckMatrix.md)
  : Check that a dataset is a proper loom matrix
- [`IsMatrixEmpty()`](https://mianaz.github.io/srtdisk/reference/IsMatrixEmpty.md)
  : Check to see if a matrix is empty

## Data Type Utilities

Utilities for handling data types and format conversions

- [`GuessDType()`](https://mianaz.github.io/srtdisk/reference/GuessDType.md)
  : Guess an HDF5 Datatype
- [`IsDType()`](https://mianaz.github.io/srtdisk/reference/IsDType.md) :
  Check the datatype of an HDF5 dataset
- [`StringType()`](https://mianaz.github.io/srtdisk/reference/StringType.md)
  : Generate an HDF5 string dtype
- [`BoolToInt()`](https://mianaz.github.io/srtdisk/reference/BoolToInt.md)
  : Convert a logical to an integer
- [`Dims()`](https://mianaz.github.io/srtdisk/reference/Dims.md) : Get
  the dimensions of an HDF5 dataset or sparse matrix
- [`GetClass()`](https://mianaz.github.io/srtdisk/reference/GetClass.md)
  : Get a class string with package information
- [`GetMargin()`](https://mianaz.github.io/srtdisk/reference/GetMargin.md)
  : Determine the margin to use for a dataset
- [`GetParent()`](https://mianaz.github.io/srtdisk/reference/GetParent.md)
  : Get the parent of an HDF5 dataset or group
- [`IndexToPointer()`](https://mianaz.github.io/srtdisk/reference/SparsePointers.md)
  [`PointerToIndex()`](https://mianaz.github.io/srtdisk/reference/SparsePointers.md)
  : Convert sparse matrix pointers to indices and vice versa
- [`ClosestVersion()`](https://mianaz.github.io/srtdisk/reference/ClosestVersion.md)
  : Find the closest version
- [`ChunkPoints()`](https://mianaz.github.io/srtdisk/reference/ChunkPoints.md)
  : Generate chunk points
- [`MakeSpace()`](https://mianaz.github.io/srtdisk/reference/MakeSpace.md)
  : Make a space

## Metadata & Attributes

Working with HDF5 attributes and metadata

- [`PadNames()`](https://mianaz.github.io/srtdisk/reference/PadNames.md)
  : Add names for unnamed or partially named objects
- [`RandomName()`](https://mianaz.github.io/srtdisk/reference/RandomName.md)
  : Generate a random string of characters
- [`FormatTime()`](https://mianaz.github.io/srtdisk/reference/Timestamp.md)
  [`Timestamp()`](https://mianaz.github.io/srtdisk/reference/Timestamp.md)
  [`TSFormats()`](https://mianaz.github.io/srtdisk/reference/Timestamp.md)
  : Create and work with timestamps
- [`UpdateKey()`](https://mianaz.github.io/srtdisk/reference/UpdateKey.md)
  : Update a Seurat key
- [`UpdateSlots()`](https://mianaz.github.io/srtdisk/reference/UpdateSlots.md)
  : Update slots in an object
- [`FixFeatures()`](https://mianaz.github.io/srtdisk/reference/FixFeatures.md)
  : Fix Feature Names
- [`CompoundToGroup()`](https://mianaz.github.io/srtdisk/reference/CompoundToGroup.md)
  : Convert an HDF5 compound dataset to a group
- [`Writeable()`](https://mianaz.github.io/srtdisk/reference/Writeable.md)
  : Is an HDF5 file or group writeable
- [`WriteMode()`](https://mianaz.github.io/srtdisk/reference/WriteMode.md)
  : Get the proper HDF5 connection mode for writing depending on
  overwrite status
