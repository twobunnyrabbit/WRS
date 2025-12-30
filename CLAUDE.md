# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

WRS is an R package containing robust statistics functions developed by R.R. Wilcox. The package is **no longer maintained** (v0.45 is the final version as of 2025-10-29).

## Package Structure

This is a standard R package with an unconventional source organization:

- `pkg/` - The actual R package directory
  - `pkg/R/Rallfun-v45.R` - **Single monolithic file containing ~1,971 functions (~97K lines)**
  - `pkg/DESCRIPTION` - Package metadata
  - `pkg/NAMESPACE` - Exports all functions with `exportPattern("^[^\\.]")`
  - `pkg/man/` - Limited documentation (only a few man files)
  - `pkg/inst/` - Contains CITATION file
- `old R files/` - Historical versions (v20-v44)
- `old Rd files/` - Archived documentation
- `check.R` - Package check and test script

All package functionality is in the single `Rallfun-v45.R` file. There is no modular organization.

## Development Commands

### Install Dependencies
```r
# Required dependencies
install.packages(c("MASS", "akima", "robustbase"))

# Suggested packages (used by various functions)
install.packages(c("akima", "cobs", "robust", "mgcv", "scatterplot3d",
                   "quantreg", "rrcov", "lars", "pwr", "trimcluster",
                   "mc2d", "psych", "Rfit", "DepthProc", "class", "fda",
                   "rankFD", "parallel", "e1071"))
```

### Build and Check Package
```r
library(tools)
library(devtools)

# Check for non-ASCII characters
for (nm in list.files("pkg/R", full.names=TRUE)) {
  cat(nm, ":")
  showNonASCIIfile(nm)
  cat("\n")
}

# Run package check (without rebuilding documentation)
devtools::check("pkg", document=FALSE)

# Install the package locally
devtools::install("pkg")
```

### Install from GitHub
```r
# Install from the GitHub repository
devtools::install_github("nicebread/WRS", subdir="pkg")
```

### Testing Functions
After installation, test individual functions:
```r
library(WRS)

# Example tests from check.R
timevec <- c(77, 87, 88, 114, 151, 210, 219, 246, 253, 262, 296, 299,
             306, 376, 428, 515, 666, 1310, 2611)
trimse(timevec, 0.1)
winmean(timevec, 0.1)
msmedse(timevec)

x <- rnorm(100)
y <- .3*x + rnorm(100)
m <- matrix(cbind(x, y), ncol=2)
outproMC(m)
wmw.RZR(x, y)
runbin.CI(x, y)
```

## Code Architecture

### Monolithic Design
The entire package is a single R source file with no internal module structure. Functions are defined sequentially without explicit organization into sections or categories.

### Function Categories
The package contains functions for:

1. **Robust Location/Dispersion Measures**: trimmed means (`trimse`, `winmean`), Harrell-Davis estimator (`hd`), M-estimators (`mest`), winsorized variance (`winvar`)

2. **Group Comparisons**:
   - Independent groups: `yuen` (trimmed means), `pb2gen`, `qcomhdMC` (quantile comparisons)
   - Dependent groups: `yuend`, `DqdifMC`
   - Bootstrap methods: most functions ending in `b` or `MC`

3. **Regression Methods**:
   - Robust regression: Theil-Sen (`tsreg`), modifications for tied values (`tshdreg`)
   - Quantile regression smoothers: `qhdsm`, `qhdsm2g`
   - Regression plotting: `rplot`, `rplotCI`, `rplotpbCI`, `reg2plot`

4. **ANCOVA Methods**:
   - Independent groups: `ancova`, `ancJN`, `anctspb`, `ancGLOB`, `ancGpar`
   - Dependent groups: `Dancova`, `Dancovamp`, `Danctspb`
   - Visualization: `ancdifplot`

5. **ANOVA Methods**: `t1way`, `MEDanova`, functions with naming pattern `*anova*`

6. **Bootstrap Utilities**: Functions with `MC` suffix use parallel processing via `mclapply`

7. **Outlier Detection**: `outpro`, `outproMC`, `outfun` parameter in many functions

### Common Patterns

- **Trimming parameter `tr`**: Default is 0.2 (20% trimming on each tail)
- **Quantile parameter `q`**: Used in quantile-based methods
- **Bootstrap resampling `nboot`**: Number of bootstrap samples (varies by function)
- **Random seed `SEED`**: Most bootstrap functions accept `SEED=TRUE/FALSE` parameter
- **Multicore `MC` suffix**: Functions ending in `MC` use parallel processing
- **Outlier handling**: Many functions have `xout`/`eout` parameters and `outfun` argument

### Dependencies on External Packages
Functions conditionally load packages as needed (not all are hard dependencies):
- `parallel` - for MC functions using `mclapply`
- Various suggested packages loaded within specific functions via `library()`

## Important Notes

1. **No documentation**: The package has minimal Rd files. Function behavior must be understood from comments in the source code and the referenced textbook.

2. **Academic license**: The package uses USC-RL v1.0 license (academic/non-commercial use only). See header of `Rallfun-v45.R`.

3. **Version history**: Old versions are preserved in `old R files/` directory. All development happens in the single main file.

4. **No active maintenance**: As of v0.45 (2025-10-29), this package is no longer maintained. No future updates are planned.

5. **WRScpp**: A C++ companion package was mentioned for performance but is no longer functional/maintained.

6. **Related textbook**: Functions implement methods from R.R. Wilcox's robust statistics textbook (5th edition). Function comments sometimes reference specific sections or papers.

## Git Workflow

- Main branch: `master`
- The repository contains the full package plus historical files
- Recent commits show maintenance updates and PR merges (last significant update: 2025-10-29)

## Active Refactoring Project

**Status**: Phase 1 COMPLETED as of 2025-12-30 ✅

This package is undergoing a major refactoring to modernize its structure while maintaining 100% backward compatibility:

- **Goal**: Transform from single 97K-line file to 20 modular files with full documentation
- **Target Version**: v0.46
- **Progress Tracking**: See `REFACTORING-PROGRESS.md` in the root directory
- **Implementation Plan**: See `.claude/plans/curious-questing-clock.md`

**Phase 1 Complete - Module Extraction** ✅:
1. ✅ Split `pkg/R/Rallfun-v45.R` → 20 focused module files in `pkg/R-new/`
2. ✅ All 1,828 unique functions successfully extracted
3. ✅ All modules source without errors
4. ✅ 100% backward compatibility maintained

**Next Steps - Phase 2 (Optimization)**:
1. Remove ~600 redundant library() calls (increased due to duplication)
2. Resolve 1,126 duplicate function definitions across modules
3. Extract common patterns to utility functions
4. Optimize and clean up code

**Remaining Phases**:
- Phase 3-4: Add full roxygen2 documentation (~1,500 user-facing functions)
- Phase 5: Final testing and validation

**Critical Files**:
- `Rallfun-v45.R.ORIGINAL` - Backup of original file (DO NOT MODIFY)
- `pkg/R-new/` - 20 new modular files (all functional) ✅
- `all-functions.txt` - Inventory of all functions
- `reference-outputs.rds` - Baseline test data for validation
- `pkg/tests/test-backward-compat.R` - Backward compatibility test suite

**Module Files** (in `pkg/R-new/`):
- Core: `00-utils-core.R`, `location.R`, `outliers.R`, `bootstrap.R`
- Analysis: `two-sample.R`, `anova.R`, `ancova.R`, `correlation.R`, `mcp.R`
- Regression: `regression.R`, `regression-advanced.R`, `covariance.R`
- Other: `medians.R`, `plotting.R`, `effect-size.R`, `power.R`, `winsorize.R`
- Specialized: `classification.R`, `special.R`, `zzz-internal.R`

**If contributing**: Always run `source("pkg/tests/test-backward-compat.R"); test_backward_compatibility()` after making changes to ensure no breaking changes.
