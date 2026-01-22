# A class for connections to h5Seurat files

A class for connections to h5Seurat files

A class for connections to h5Seurat files

## Format

An [`R6Class`](https://r6.r-lib.org/reference/R6Class.html) object

## See also

[`H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md)

## Super classes

[`hdf5r::H5RefClass`](http://hhoeflin.github.io/hdf5r/reference/H5RefClass-class.md)
-\>
[`hdf5r::H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md)
-\>
[`srtdisk::scdisk`](https://mianaz.github.io/srtdisk/reference/scdisk-class.md)
-\> `h5Seurat`

## Methods

### Public methods

- [`h5Seurat$index()`](#method-h5Seurat-index)

- [`h5Seurat$set.version()`](#method-h5Seurat-set.version)

- [`h5Seurat$version()`](#method-h5Seurat-version)

Inherited methods

- [`hdf5r::H5RefClass$close()`](http://hhoeflin.github.io/hdf5r/reference/H5RefClass-class.html#method-close)
- [`hdf5r::H5RefClass$dec_ref()`](http://hhoeflin.github.io/hdf5r/reference/H5RefClass-class.html#method-dec_ref)
- [`hdf5r::H5RefClass$get_file_id()`](http://hhoeflin.github.io/hdf5r/reference/H5RefClass-class.html#method-get_file_id)
- [`hdf5r::H5RefClass$get_obj_type()`](http://hhoeflin.github.io/hdf5r/reference/H5RefClass-class.html#method-get_obj_type)
- [`hdf5r::H5RefClass$get_ref()`](http://hhoeflin.github.io/hdf5r/reference/H5RefClass-class.html#method-get_ref)
- [`hdf5r::H5RefClass$inc_ref()`](http://hhoeflin.github.io/hdf5r/reference/H5RefClass-class.html#method-inc_ref)
- [`hdf5r::H5RefClass$methods()`](http://hhoeflin.github.io/hdf5r/reference/H5RefClass-class.html#method-methods)
- [`hdf5r::H5File$attr_delete()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_delete)
- [`hdf5r::H5File$attr_delete_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_delete_by_idx)
- [`hdf5r::H5File$attr_delete_by_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_delete_by_name)
- [`hdf5r::H5File$attr_exists()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_exists)
- [`hdf5r::H5File$attr_exists_by_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_exists_by_name)
- [`hdf5r::H5File$attr_get_number()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_get_number)
- [`hdf5r::H5File$attr_info_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_info_by_idx)
- [`hdf5r::H5File$attr_info_by_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_info_by_name)
- [`hdf5r::H5File$attr_name_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_name_by_idx)
- [`hdf5r::H5File$attr_open()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_open)
- [`hdf5r::H5File$attr_open_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_open_by_idx)
- [`hdf5r::H5File$attr_open_by_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_open_by_name)
- [`hdf5r::H5File$attr_rename()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_rename)
- [`hdf5r::H5File$attr_rename_by_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-attr_rename_by_name)
- [`hdf5r::H5File$close_all()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-close_all)
- [`hdf5r::H5File$commit()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-commit)
- [`hdf5r::H5File$create_attr()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-create_attr)
- [`hdf5r::H5File$create_attr_by_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-create_attr_by_name)
- [`hdf5r::H5File$create_dataset()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-create_dataset)
- [`hdf5r::H5File$create_group()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-create_group)
- [`hdf5r::H5File$create_reference()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-create_reference)
- [`hdf5r::H5File$exists()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-exists)
- [`hdf5r::H5File$file_info()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-file_info)
- [`hdf5r::H5File$flush()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-flush)
- [`hdf5r::H5File$get_filename()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-get_filename)
- [`hdf5r::H5File$get_filesize()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-get_filesize)
- [`hdf5r::H5File$get_intent()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-get_intent)
- [`hdf5r::H5File$get_obj_count()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-get_obj_count)
- [`hdf5r::H5File$get_obj_ids()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-get_obj_ids)
- [`hdf5r::H5File$get_obj_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-get_obj_name)
- [`hdf5r::H5File$group_info()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-group_info)
- [`hdf5r::H5File$group_info_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-group_info_by_idx)
- [`hdf5r::H5File$group_info_by_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-group_info_by_name)
- [`hdf5r::H5File$link()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link)
- [`hdf5r::H5File$link_copy_from()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_copy_from)
- [`hdf5r::H5File$link_copy_to()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_copy_to)
- [`hdf5r::H5File$link_create_external()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_create_external)
- [`hdf5r::H5File$link_create_hard()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_create_hard)
- [`hdf5r::H5File$link_create_soft()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_create_soft)
- [`hdf5r::H5File$link_delete()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_delete)
- [`hdf5r::H5File$link_delete_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_delete_by_idx)
- [`hdf5r::H5File$link_exists()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_exists)
- [`hdf5r::H5File$link_info()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_info)
- [`hdf5r::H5File$link_info_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_info_by_idx)
- [`hdf5r::H5File$link_move_from()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_move_from)
- [`hdf5r::H5File$link_move_to()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_move_to)
- [`hdf5r::H5File$link_name_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_name_by_idx)
- [`hdf5r::H5File$link_value()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_value)
- [`hdf5r::H5File$link_value_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-link_value_by_idx)
- [`hdf5r::H5File$ls()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-ls)
- [`hdf5r::H5File$mount()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-mount)
- [`hdf5r::H5File$obj_copy_from()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-obj_copy_from)
- [`hdf5r::H5File$obj_copy_to()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-obj_copy_to)
- [`hdf5r::H5File$obj_info()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-obj_info)
- [`hdf5r::H5File$obj_info_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-obj_info_by_idx)
- [`hdf5r::H5File$obj_info_by_name()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-obj_info_by_name)
- [`hdf5r::H5File$open()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-open)
- [`hdf5r::H5File$open_by_idx()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-open_by_idx)
- [`hdf5r::H5File$path_valid()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-path_valid)
- [`hdf5r::H5File$print()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-print)
- [`hdf5r::H5File$unmount()`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.html#method-unmount)
- [`srtdisk::scdisk$chunk.points()`](https://mianaz.github.io/srtdisk/reference/scdisk.html#method-chunk.points)
- [`srtdisk::scdisk$finalizer()`](https://mianaz.github.io/srtdisk/reference/scdisk.html#method-finalizer)
- [`srtdisk::scdisk$initialize()`](https://mianaz.github.io/srtdisk/reference/scdisk.html#method-initialize)
- [`srtdisk::scdisk$last.modified()`](https://mianaz.github.io/srtdisk/reference/scdisk.html#method-last.modified)
- [`srtdisk::scdisk$timestamp()`](https://mianaz.github.io/srtdisk/reference/scdisk.html#method-timestamp)

------------------------------------------------------------------------

### Method `index()`

Get the index for this h5Seurat file

#### Usage

    h5Seurat$index()

------------------------------------------------------------------------

### Method `set.version()`

Set the version attribute

#### Usage

    h5Seurat$set.version(version)

#### Arguments

- `version`:

  A version number matching the regex `^\d+(\.\d+){2}(\.9\d{3})?$`

------------------------------------------------------------------------

### Method [`version()`](https://rdrr.io/r/base/Version.html)

Get the version attribute

#### Usage

    h5Seurat$version()
