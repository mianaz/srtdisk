# Bug诊断：AnnData obs分类变量在转换时丢失

## 问题描述

用户报告：从 h5ad (AnnData) 转换到 h5seurat，再加载到 Seurat 对象时，obs 中的分类变量（如细胞类型）丢失了。只有数值列能可靠地出现在 Seurat 的 meta.data 中。

## 根本原因分析

### 数据流转路径

```
h5ad file
  ↓ (H5ADToH5Seurat in Convert.R)
h5seurat file (meta.data group)
  ↓ (LoadH5Seurat in LoadH5Seurat.R)
  ↓ (as.data.frame.H5Group in ReadH5.R)
Seurat object (meta.data slot)
```

### 关键代码位置

#### 1. Convert.R:586-620 - H5AD到h5Seurat的转换

```r
# Add cell-level metadata
if (source$exists(name = 'obs') && inherits(x = source[['obs']], what = 'H5Group')) {
  dfile$obj_copy_from(
    src_loc = source,
    src_name = 'obs',
    dst_name = 'meta.data'
  )
  # Normalize h5ad categorical format (categories/codes) to SeuratDisk format (levels/values)
  NormalizeH5ADCategorical <- function(dfgroup) {
    for (col_name in names(dfgroup)) {
      col_obj <- dfgroup[[col_name]]
      # Check if this is an h5ad categorical group (has 'categories' and 'codes')
      if (inherits(col_obj, 'H5Group') &&
          all(c('categories', 'codes') %in% names(col_obj))) {
        # Rename 'categories' to 'levels'
        # Rename 'codes' to 'values'
        # ...
      }
    }
  }
  NormalizeH5ADCategorical(dfgroup = dfile[['meta.data']])
  ColToFactor(dfgroup = dfile[['meta.data']])
}
```

**问题1**：`NormalizeH5ADCategorical` 只处理新格式（categories/codes），而 `ColToFactor` 只处理旧格式（__categories）。这两个函数可能存在兼容性问题。

**问题2**：`ColToFactor` (Convert.R:337-393) 期望的是 `__categories` 组加上各列的编码，但实际h5ad文件可能使用不同的格式。

#### 2. ReadH5.R:104-212 - 读取H5Group为data.frame

**核心bug在这里！**

```r
as.data.frame.H5Group <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- NULL
  idx <- NULL
  colnames <- NULL

  for (i in names(x = x)) {
    if (i == 'index') {
      next
    }

    # Handle H5Groups (which might be factors)
    if (inherits(x[[i]], "H5Group")) {
      if (IsFactor(x = x[[i]])) {
        tryCatch({
          values <- as.integer(x = x[[i]][['values']][])
          levels <- x[[i]][['levels']][]
          values <- values + 1L  # Convert to 1-based
          # ...
          df[[i]] <- factor(x = levels[values], levels = levels)
        }, error = function(e) {
          # Skip if can't read as factor
          NULL  # ⚠️ 问题在这里！静默跳过！
        })
      }
      next
    }

    # Safe reading with error handling
    tryCatch({
      dset <- as.vector(x = x[[i]][])
      # ...
      df[[i]] <- dset
    }, error = function(e) {
      NULL  # ⚠️ 问题在这里！静默跳过！
    })
  }

  return(df)
}
```

**关键问题**：
1. **Line 150-153**: 如果读取factor时出错，静默跳过（不报错，不记录）
2. **Line 192-195**: 如果读取任何列时出错，静默跳过
3. **没有日志或警告**，用户完全不知道发生了什么

### 可能的失败场景

#### 场景1：values索引越界
```r
values <- as.integer(x = x[[i]][['values']][])  # 可能是 [0, 1, 2, 3]
values <- values + 1L  # 变成 [1, 2, 3, 4]
levels <- x[[i]][['levels']][]  # 可能只有 ["Type1", "Type2", "Type3"] (3个)
# 如果 values 中有4，levels[4] 会越界！
```

#### 场景2：levels或values不存在
- h5ad可能使用不同的命名：`categories` vs `levels`
- 转换过程中可能没有正确重命名

#### 场景3：数据类型不匹配
- values可能不是整数类型
- levels可能不是字符串类型

## 验证步骤

### 1. 检查h5ad文件中obs的实际结构

```r
library(hdf5r)
h5ad <- H5File$new("your_file.h5ad", mode = "r")

# 查看obs结构
print(names(h5ad[["obs"]]))

# 检查某个分类列
if (h5ad[["obs"]]$exists("cell_type")) {
  cat_col <- h5ad[["obs"]][["cell_type"]]
  print(class(cat_col))
  print(names(cat_col))  # 应该看到是什么结构
}

# 检查__categories
if (h5ad[["obs"]]$exists("__categories")) {
  print(names(h5ad[["obs"]][["__categories"]]))
}

h5ad$close_all()
```

### 2. 检查h5seurat文件中meta.data的结构

```r
library(hdf5r)
h5s <- H5File$new("converted.h5seurat", mode = "r")

# 查看meta.data结构
print(names(h5s[["meta.data"]]))

# 检查某个列的结构
if (h5s[["meta.data"]]$exists("cell_type")) {
  col <- h5s[["meta.data"]][["cell_type"]]
  print(class(col))
  print(names(col))  # 应该有 levels 和 values

  if (inherits(col, "H5Group")) {
    if (col$exists("levels")) {
      print("Levels:")
      print(col[["levels"]][])
    }
    if (col$exists("values")) {
      print("Values (first 10):")
      print(head(col[["values"]][], 10))
    }
  }
}

h5s$close_all()
```

## 修复方案

### 方案1：改进error handling（推荐）

在 `ReadH5.R` 的 `as.data.frame.H5Group` 中添加详细的error logging：

```r
tryCatch({
  values <- as.integer(x = x[[i]][['values']][])
  levels <- x[[i]][['levels']][]
  values <- values + 1L

  # 验证索引范围
  if (any(values > length(levels) | values < 1, na.rm = TRUE)) {
    warning(sprintf(
      "Column '%s': values out of range [%d-%d], levels length: %d",
      i, min(values, na.rm = TRUE), max(values, na.rm = TRUE), length(levels)
    ), call. = FALSE, immediate. = TRUE)
    # 截断超出范围的值为NA
    values[values > length(levels) | values < 1] <- NA
  }

  df[[i]] <- factor(x = levels[values], levels = levels)
}, error = function(e) {
  warning(sprintf(
    "Failed to read factor column '%s': %s",
    i, conditionMessage(e)
  ), call. = FALSE, immediate. = TRUE)
})
```

### 方案2：修复ColToFactor和NormalizeH5ADCategorical的协调

确保 `NormalizeH5ADCategorical` 正确处理所有h5ad分类格式。

### 方案3：添加验证函数

在转换完成后验证数据完整性：

```r
ValidateMetadataConversion <- function(source, dest, verbose = TRUE) {
  # 检查源文件obs列数
  source_cols <- setdiff(names(source[["obs"]]), c("_index", "index", "__categories"))

  # 检查目标文件meta.data列数
  dest_cols <- setdiff(names(dest[["meta.data"]]), c("_index", "index"))

  # 比较
  missing_cols <- setdiff(source_cols, dest_cols)

  if (length(missing_cols) > 0 && verbose) {
    warning(sprintf(
      "The following columns were not converted: %s",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE, immediate. = TRUE)
  }

  return(list(
    source_cols = source_cols,
    dest_cols = dest_cols,
    missing_cols = missing_cols
  ))
}
```

## 测试计划

### 1. 创建最小可复现示例

创建一个简单的h5ad文件，包含：
- 数值列
- 字符串列
- 分类列（使用__categories格式）
- 分类列（使用categories/codes格式）

### 2. 单元测试

```r
test_that("Categorical obs columns are preserved in h5ad->h5seurat conversion", {
  # 创建测试h5ad文件
  temp_h5ad <- create_test_h5ad_with_categorical()
  temp_h5seurat <- tempfile(fileext = ".h5seurat")

  # 转换
  Convert(temp_h5ad, temp_h5seurat)

  # 加载
  seurat_obj <- LoadH5Seurat(temp_h5seurat)

  # 验证
  expect_true("cell_type" %in% colnames(seurat_obj@meta.data))
  expect_true(is.factor(seurat_obj@meta.data$cell_type))
  expect_gt(length(levels(seurat_obj@meta.data$cell_type)), 0)
})
```

## 下一步行动

1. ✅ 诊断完成
2. ⏳ 创建测试用例重现问题
3. ⏳ 实现修复（改进error handling）
4. ⏳ 验证修复
5. ⏳ 更新文档

---
*诊断日期: 2025-11-12*
*相关issue: GitHub issue报告的分类变量丢失问题*
