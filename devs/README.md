# Development Tools

This directory contains development utilities and scripts for srtdisk package maintenance.

## roxygenize_workaround.R

**Purpose:** Workaround for roxygen2 bug with S4 generics (parser_setGeneric error)

**Issue:** Roxygen2 versions 7.0.0 through 7.3.3 have a bug where `methods::getGeneric()` returns NULL during documentation generation, causing the error:
```
Error in value@generic :
  no applicable method for `@` applied to an object of class "NULL"
```

**Solution:** This script temporarily patches roxygen2's `parser_setGeneric` function to work around the bug.

### Usage

From the package root directory:

```r
# From R console
source("devs/roxygenize_workaround.R")

# Or from command line
Rscript devs/roxygenize_workaround.R
```

### What it does

1. Patches roxygen2's `parser_setGeneric` function in memory
2. Creates placeholder objects for S4 generics instead of calling `methods::getGeneric()`
3. Runs `roxygen2::roxygenize()` normally
4. Preserves all @export directives and documentation

### Notes

- **Safe:** Only modifies roxygen2 in memory during execution
- **No permanent changes** to R installation
- **R CMD build/check/test work normally** - they don't use roxygen2
- **Remove this workaround** once roxygen2 bug is fixed upstream

### Bug Report

This bug should be reported to: https://github.com/r-lib/roxygen2/issues

The bug was introduced in roxygen2 version 7.0.0.
