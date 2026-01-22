# h5Seurat File Format Specification

## Overall File Structure

### Required Attributes

There are three required top-level HDF5 attributes: “project”,
“active.assay”, and “version”. Each of these must be a single
[character](#character-representation). The “project” attribute
corresponds to the project value of a `Seurat` object; the
“active.assay” attribute is the name of the default assay and must be
present in the [“assays” group](#assay-expression-data). The “version”
corresponds to the version of Seurat that the h5Seurat file is based on.

### Top-Level Datasets and Groups

There are two required top-level HDF5 datasets: “cell.names” and
“meta.data”. The “cell.names” dataset should be a one-dimensional
[character](#character-representation) dataset, with a length equal to
the number of cells present in the data. Cell names are not stored
anywhere else in the h5Seurat file.

The “meta.data” dataset contains cell-level metadata. It should be
stored as either an [HDF5 dataset or group](#data-frame-representation),
depending on the contents of the meta data. See the [data frame
representation](#data-frame-representation) for more details.

## Assay Expression Data

[`Assay`](https://rdrr.io/cran/SeuratObject/man/Assay-class.html)
objects are stored in the top-level group “assays”; each assay is stored
as its own group within the “assays” group. Within each assay group,
there must be a dataset named “features” and either a dataset or group
named “data”; the “features” dataset must be a one-dimensional
[character](#character-representation) dataset with a length equal to
the number of total features within the assay. The “data” entry is a
matrix, with dimensions of $`m_{features} x n_{cells}`$; this entry may
be either a dataset, if “data” is a [dense
matrix](#dense-matrix-representation), or a group, if “data” is a
[sparse matrix](#sparse-matrix-representation). Assay groups must also
have an attribute named “key”; this is a single
[character](#character-representation) value.

Assay groups may also have the following optional groups and datasets:

- “counts”: either a [dense](#dense-matrix-representation) or
  [sparse](#sparse-matrix-representation) matrix; must have the same
  dimesions as “data”
- “scale.data”: a [dense](#dense-matrix-representation) matrix; if
  “scale.data” is present, a one-dimensional
  [character](#character-representation) dataset must also be present.
  The “scale.data” matrix must be of dimensions
  $`m_{scaled features} x n_{cells}`$
- “meta.features”: a [data frame](#data-frame-representation) with the
  same number of rows as values present in “features”
- “variable.features”: a one-dimensional
  [character](#character-representation) dataset
- “misc”: a [list](#list-and-custom-class-representation)

Subclasses of
[`Assay`](https://rdrr.io/cran/SeuratObject/man/Assay-class.html)
objects must also follow the same rules as [custom S4
classes](#list-and-custom-class-representation).

**Seurat V5 Assay5 Support**: This implementation fully supports Seurat
v5 `Assay5` objects. When saving to h5Seurat, `Assay5` objects are
automatically converted to the storage format described above, with
layers (counts, data, scale.data) preserved as separate matrices. When
loading, the format is compatible with both legacy `Assay` and modern
`Assay5` objects.

## Dimensional Reductions

[Dimensional reduction
information](https://rdrr.io/cran/SeuratObject/man/DimReduc-class.html)
is stored in the top-level group “reductions”; each dimensional
reduction is stored as its own group within the “reductions” group.
Within each dimensional reduction group, there are three required
attributes: “active.assay”, “key”, and “global”; “active.assay” must be
one or more [character](#character-representation) values where each
value is a name of an [assay](#assay-expression-data), “key” must be a
single a [character](#character-representation) value, and “global” must
be a single [logical](#logical-representation) value. In addition, there
must also be a dataset named “cell.embeddings” representing a [dense
matrix](#dense-matrix-representation). This matrix must have the same
number of rows as cells present in the h5Seurat file.

Dimensional reduction groups may also have the following optional groups
and datasets:

- “feature.loadings”: …
- “feature.loadings.projected”: …
- “misc”: a [list](#list-and-custom-class-representation)
- “jackstraw”: a [custom S4
  group](#list-and-custom-class-representation)

## Nearest-Neighbor Graphs

[Nearest-neighbor
graphs](https://rdrr.io/cran/SeuratObject/man/Graph-class.html) are
stored in the top-level group “graphs”; each graph is stored as its own
group within the “graphs” group. Graph names become graph group names.
Graphs are stored as [sparse matrices](#sparse-matrix-representation)
with an additional HDF5 attribute: “assay.used”. This HDF5 attribute
should be a single [character](#character-representation) value.

## Spatial Image Data

Spatial image data is stored in the top-level group “images”; each image
is stored as its own group within the “images” group. Actual structure
of the image group is dependent on the structure of the spatial image
data. However, it follows the same rules as [custom S4
classes](#list-and-custom-class-representation).

**Note**: spatial images are only supported in objects that were
generated by a version of Seurat that has spatial support. Currently,
this is restricted to version 3.1.5.9900 or higher.

## Command Logs

## Miscellaneous Information and Tool-Specific Results

Miscellaneous information is stored in the top-level group “misc”; this
group follows the same runs as
[lists](#list-and-custom-class-representation). The “misc” group is
required to be present, but not required to be filled.

Tool-specific results are stored in the top-level group “tools”; this
group follows the same runs as
[lists](#list-and-custom-class-representation). The “tools” group is
required to be present, but not required to be filled.

## Common Data Structures

Some data types are found commonly throughout `Seurat` objects

### Character Representation

All
[character](https://stat.ethz.ch/R-manual/R-devel/library/base/html/character.html)
values (strings in other languages) should be encoding as
variable-length UTF-8 strings; this applies to HDF5 datasets (including
standalone
[string](https://hhoeflin.github.io/hdf5r/reference/H5T_STRING-class.html)
datasets as well as parts of
[compound](https://hhoeflin.github.io/hdf5r/reference/H5T_COMPOUND-class.html)
datasets) and HDF5 attributes.

### Dense Matrix Representation

Dense matrices should be stored as a two-dimensional dataset of any
type. Datasets should be written in a
[column-major](https://en.wikipedia.org/wiki/Row-_and_column-major_order)
order. For column-major implementations (eg. R, Fortran), dataset
dimensions on-disk should be the same as dimensions in-memory (eg.
$`m_{diskrow} x n_{diskcol} \sim m_{memrow} x n_{memcol}`$). For
row-major implmentations (eg. C/C++, Python), dataset dimensions on-disk
should appear *transposed* to dimensions in-memory (eg.
$`m_{diskrow} x n_{diskcol} \sim n_{memrow} x m_{memcol}`$); row-major
implmemetnations transpose datasets prior to reading and writing data.

    #>      [,1] [,2] [,3] [,4]
    #> [1,]    0    0    0    1
    #> [2,]    0    1    1    1
    #> [3,]    1    1    0    0
    #> Class: H5D
    #> Dataset: /densemat
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Access type: H5F_ACC_RDWR
    #> Datatype: H5T_IEEE_F64LE
    #> Space: Type=Simple     Dims=3 x 4     Maxdims=Inf x Inf
    #> Chunk: 32 x 32
    #>      [,1] [,2] [,3] [,4]
    #> [1,]    0    0    0    1
    #> [2,]    0    1    1    1
    #> [3,]    1    1    0    0

### Sparse Matrix Representation

Sparse matrices are stored as an HDF5 group with three datasets:
“indices”, “indptr”, and “data”; the “indices” and “data” datasets must
be the same length. “data” represents each non-zero element of the
matrix. “indices” represents the $`0`$-based row numbers for each value
in “data”

    #> 3 x 4 sparse Matrix of class "dgCMatrix"
    #>             
    #> [1,] . . . 1
    #> [2,] . 1 1 1
    #> [3,] 1 1 . .
    #> Class: H5Group
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Group: /sparsemat
    #> Attributes: dims
    #> Listing:
    #>     name    obj_type dataset.dims dataset.type_class
    #>     data H5I_DATASET            6          H5T_FLOAT
    #>  indices H5I_DATASET            6        H5T_INTEGER
    #>   indptr H5I_DATASET            5        H5T_INTEGER
    #> Class: H5D
    #> Dataset: /sparsemat/indices
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Access type: H5F_ACC_RDWR
    #> Datatype: H5T_STD_I32LE
    #> Space: Type=Simple     Dims=6     Maxdims=Inf
    #> Chunk: 2048
    #> [1] 2 1 2 1 0 1
    #> Class: H5D
    #> Dataset: /sparsemat/data
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Access type: H5F_ACC_RDWR
    #> Datatype: H5T_IEEE_F64LE
    #> Space: Type=Simple     Dims=6     Maxdims=Inf
    #> Chunk: 1024
    #> [1] 1 1 1 1 1 1

“indptr” represents the points in “data” at which a new column is
started. This dataset is $`0`$-based and should be $`n_{columns} + 1`$
in length.

    #> Class: H5D
    #> Dataset: /sparsemat/indptr
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Access type: H5F_ACC_RDWR
    #> Datatype: H5T_STD_I32LE
    #> Space: Type=Simple     Dims=5     Maxdims=Inf
    #> Chunk: 2048
    #> [1] 0 1 3 4 6

The “indices”, “indptr”, and “data” datasets correspond to the “i”, “p”,
and “x” slots in a
[`dgCMatrix`](https://rdrr.io/cran/Matrix/man/dgCMatrix-class.html),
respectively.

There may optionally be an HDF5 attribute called “dims”; this attribute
should be a two
[integer](https://hhoeflin.github.io/hdf5r/reference/H5T_INTEGER-class.html)
values corresponding to the number of rows and number of columns, in
that order, in the sparse matrix.

    #> [1] 3 4

### Factor Representation

[Factors](https://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html)
should be stored as an HDF5 group with two datasets: “levels” and
“values”

    #> [1] g1 g2 g1 g1 g2
    #> Levels: g1 g2
    #> Class: H5Group
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Group: /factor
    #> Listing:
    #>    name    obj_type dataset.dims dataset.type_class
    #>  levels H5I_DATASET            2         H5T_STRING
    #>  values H5I_DATASET            5        H5T_INTEGER

The “levels” dataset should be a [character](#character-representation)
dataset with one entry per
[level](https://stat.ethz.ch/R-manual/R-devel/library/base/html/levels.html)

    #> Class: H5D
    #> Dataset: /factor/levels
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Access type: H5F_ACC_RDWR
    #> Datatype: H5T_STRING {
    #>       STRSIZE H5T_VARIABLE;
    #>       STRPAD H5T_STR_NULLTERM;
    #>       CSET H5T_CSET_UTF8;
    #>       CTYPE H5T_C_S1;
    #>    }
    #> Space: Type=Simple     Dims=2     Maxdims=Inf
    #> Chunk: 1024
    #> [1] "g1" "g2"

The “values” dataset should be an integer dataset with one entry per
value in the original factor. These integers should correspond to the
factor level they had in R

    #> Class: H5D
    #> Dataset: /factor/values
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Access type: H5F_ACC_RDWR
    #> Datatype: H5T_STD_I32LE
    #> Space: Type=Simple     Dims=5     Maxdims=Inf
    #> Chunk: 2048
    #> [1] 1 2 1 1 2

The number of unique entries in “values” should not exceed the number of
unique entries in “levels”

Rationale

Storing factors in this manner seems excessive from an R perspective.
HDF5 has the concept of [enumerated datasets
(enums)](https://hhoeflin.github.io/hdf5r/reference/H5T_ENUM-class.html)
which are an efficient way to store R factors on-disk. However, [some
implementations of HDF5 do not support
enums](http://docs.h5py.org/en/stable/special.html#enumerated-types) and
thus lose factor level information. In order to make h5Seurat as
cross-language as possible, we’ve opted to store factors as HDF5 groups
instead of HDF5 enums.

### Logical Representation

[Logical](https://stat.ethz.ch/R-manual/R-devel/library/base/html/logical.html)
values (booleans in other languages) are encoded as integers in the
following manner: `FALSE` is encoded as `0`, `TRUE` is encoded as `1`,
and `NA` is encoded as `2`

    #> [1]  TRUE FALSE    NA
    #> Class: H5D
    #> Dataset: /logicals
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Access type: H5F_ACC_RDWR
    #> Attributes: s3class
    #> Datatype: H5T_STD_I32LE
    #> Space: Type=Simple     Dims=3     Maxdims=Inf
    #> Chunk: 2048
    #> [1] 1 0 2

An optional HDF5 attribute named “s3class” with the value “logical” is
allowed to enforce reading in the dataset as logical values. This HDF5
attribute is a single [character](#character-representation) value.

Rationale

Unlike most languages,
[logicals](https://stat.ethz.ch/R-manual/R-devel/library/base/html/logical.html)
(or booleans) can take one of *three* values: `TRUE`, `FALSE`, or
[`NA`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/NA.html);
as such, an extra integer value is needed to handle the additional
logical value. Typically, these values are stored as
[enums](https://hhoeflin.github.io/hdf5r/reference/H5T_ENUM-class.html)
with mappings between the logical representation and integer value.
However, [some implementations of HDF5 do not support
enums](http://docs.h5py.org/en/stable/special.html#enumerated-types) and
thus lose the logical representation. Since the mappings are lost, all
logicals are stored as integers instead.

### Data Frame Representation

There are two ways of storing [data
frames](https://stat.ethz.ch/R-manual/R-devel/library/base/html/data.frame.html)
in h5Seurat files: as [datasets](#data-frame-datasets) or as
[groups](#data-frame-groups). Data frame groups are required when data
frames contain
[factors](https://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html);
if no factors are present, data frames can be stored in either type.

Rationale

[Storing factors](#factor-representation) on-disk in an HDF5 file poses
a unique set of challenges. Namely, [enumerated datasets
(enums)](https://hhoeflin.github.io/hdf5r/reference/H5T_ENUM-class.html),
which are an ideal method of storing mapping values, [are not supported
in some implementations of
HDF5](http://docs.h5py.org/en/stable/special.html#enumerated-types). As
factor level information is lost under implementations that do not
support HDF5 enums, we need a method of storing factors in a
cross-language manner. Two options were presented: groups or [compound
datasets](https://hhoeflin.github.io/hdf5r/reference/H5T_COMPOUND-class.html).
While the former seems excessive, the latter presents issues with
unequal dataset length. Therefore, to accomodate factor level
information in data frames, we utilize HDF5 groups to store data frames
when one or more columns are
[factors](https://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html).

#### Data Frame Datasets

Data frames stored as datasets should be a one-dimensional [compound
dataset](https://hhoeflin.github.io/hdf5r/reference/H5T_COMPOUND-class.html).
The single dimension should be equal to the number of observations
(number of rows) in the data frame. Each data type in the compound
dataset must adhere to the same requirements as standard datasets (eg.
[character encodings](#character-representation), [logical
mapping](#logical-representation), etc). The compound labels should
correspond to the data frame column names.

    #>    x y
    #> 1 g1 1
    #> 2 g1 2
    #> 3 g2 3
    #> 4 g1 4
    #> 5 g2 5
    #> Class: H5D
    #> Dataset: /dfdataset
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Access type: H5F_ACC_RDWR
    #> Datatype: H5T_COMPOUND {
    #>       H5T_STRING {
    #>          STRSIZE H5T_VARIABLE;
    #>          STRPAD H5T_STR_NULLTERM;
    #>          CSET H5T_CSET_UTF8;
    #>          CTYPE H5T_C_S1;
    #>       } "x" : 0;
    #>       H5T_STD_I32LE "y" : 8;
    #>    }
    #> Space: Type=Simple     Dims=5     Maxdims=Inf
    #> Chunk: 682
    #>    x y
    #> 1 g1 1
    #> 2 g1 2
    #> 3 g2 3
    #> 4 g1 4
    #> 5 g2 5

Row names are not stored with the dataset itself, but may be stored
elsewhere in the h5Seurat file, typically named
`dataset_name.row.names`; an optional HDF5 attribute called “logicals”
containing the names of
[logical](https://stat.ethz.ch/R-manual/R-devel/library/base/html/logical.html)
columns is allowed. This attribute consists of [character
values](#character-representation).

#### Data Frame Groups

Data frames stored as groups are used when
[factors](https://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html)
are present in the data frame. Within the data frame group, there should
be one dataset or group per column. Columns that are factors are [stored
as groups](#factor-representation) while all other columns are stored as
one-dimensional datasets. Each dataset within the group must adhere to
the same requirements as standard datasets (eg. [character
encodings](#character-representation), [logical
mapping](#logical-representation), etc). The names of the datasets
within the group correspond to the data frame column names

    #> Class: H5Group
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Group: /dfgroup
    #> Attributes: colnames, _index
    #> Listing:
    #>    name    obj_type dataset.dims dataset.type_class
    #>  _index H5I_DATASET            5         H5T_STRING
    #>       x   H5I_GROUP         <NA>               <NA>
    #>       y H5I_DATASET            5        H5T_INTEGER

Data frame row names may be stored in a dataset called “row.names”
within the group; this dataset should be a one-dimensional
[character](#character-representation) dataset. There are two optional
attributes allowed: “colnames” and “logicals”; the “colnames” attribute
contains the column names in the same order as was present in the
in-memory [data
frame](https://stat.ethz.ch/R-manual/R-devel/library/base/html/data.frame.html)
as [character values](#character-representation). This is used to
control column order when reading the data frame back into memory. Note,
the “colnames” attribute does not need to contain the name of every
dataset.

The “logicals” attribute contains the names of
[logical](https://stat.ethz.ch/R-manual/R-devel/library/base/html/logical.html)
columns; this attribute should consist of [character
values](#character-representation).

### List and Custom Class Representation

[Lists](https://stat.ethz.ch/R-manual/R-devel/library/base/html/list.html)
are stored as HDF5 groups. Each entry in a list must be named; the names
serve as the names of datasets and groups within the list group. List
values are stored as HDF5 datasets or groups, depending on their R
object type. For example, a list within a list would be stored as a
group within the first group.

    #> $a
    #> [1] 1 2 3
    #> 
    #> $b
    #> $b$b1
    #> [1] "hello"
    #> 
    #> $b$b2
    #> [1] "tomato" "potato"
    #> Class: H5Group
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Group: /list
    #> Attributes: names
    #> Listing:
    #>  name    obj_type dataset.dims dataset.type_class
    #>     a H5I_DATASET            3        H5T_INTEGER
    #>     b   H5I_GROUP         <NA>               <NA>
    #> Class: H5Group
    #> Filename: /private/var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T/RtmpwezU36/file55f71184be5.h5
    #> Group: /list/b
    #> Attributes: names
    #> Listing:
    #>  name    obj_type dataset.dims dataset.type_class
    #>    b1 H5I_DATASET            1         H5T_STRING
    #>    b2 H5I_DATASET            2         H5T_STRING

Custom classes are stored as
[lists](#list-and-custom-class-representation) with a special HDF5
attribute. S3 classes are stored with the attribute “s3class”, which has
a value equal to the class of the object. This attribute is a
variable-length [character](#character-representation) value.

Custom S4 classes are also stored as
[lists](#list-and-custom-class-representation) where each entry is a
slot in the S4 class. S4 class groups have an attribute named “s4class”;
this attribute should be a single-length
[character](#character-representation) storing the name of the class and
the package that defines the class in the form of `package:class` (eg.
`Signac:ChromatinAssay`); custom S4 classes defined in the
[Seurat](https://satijalab.org/seurat) package can be named with just
the class of the object.
