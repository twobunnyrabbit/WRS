# WRS Package Refactoring Progress

**Project**: Transform WRS from monolithic 97K-line file to modular, documented package
**Version**: v0.45 → v0.46
**Timeline**: 12 weeks
**Started**: 2025-12-30
**Last Updated**: 2025-12-30

---

## Quick Reference

- **Full Plan**: See `.claude/plans/curious-questing-clock.md` (approved plan with all details)
- **Original Source**: `pkg/R/Rallfun-v45.R` (97,199 lines, 1,971 functions)
- **Backup**: `Rallfun-v45.R.ORIGINAL` (2.6 MB - DO NOT MODIFY)
- **Function Inventory**: `all-functions.txt` (sorted list of all functions)
- **Reference Tests**: `reference-outputs.rds` (baseline outputs from v0.45)

---

## Current Status

**Phase 1, Week 3 IN PROGRESS - 2025-12-30**

### Progress Summary
- ✅ **Week 1 COMPLETED**: Foundation modules (4 files, 219 functions)
  - 00-utils-core.R, location.R, outliers.R, bootstrap.R
- ✅ **Week 2 COMPLETED**: Main analysis modules (3 files, 256 functions)
  - two-sample.R, anova.R, correlation.R
- 🔄 **Week 3 IN PROGRESS**: Specialized modules
  - ✅ ancova.R extracted (149 functions, 11,014 lines, 286.1 KB)
  - ✅ regression.R extracted (98 functions, 5,051 lines, 139.8 KB)
  - ✅ regression-advanced.R extracted (75 functions, 3,433 lines, 95.2 KB)
  - ✅ covariance.R extracted (43 functions, 1,034 lines, 27 KB)
  - ✅ mcp.R extracted (106 functions, 9,046 lines, 255 KB)
  - ✅ medians.R extracted (42 functions, 2,474 lines, 69 KB)
  - 🔄 Next: plotting.R and remaining modules

### Overall Metrics
- **Modules completed**: 13 of 20 (65%)
- **Functions extracted**: 1,003 of 1,971 (50.9%)
- **Lines extracted**: 55,601 of 97,199 (57.2%)
- **Total size**: 1.6 MB of ~2.6 MB (61.5%)
- **Status**: All 13 modules source successfully ✅

### Recently Completed (2025-12-30)
**medians.R module** - Median-based Methods:
- 42 functions, 2,474 lines, 69 KB
- Marginal medians: msmed, msmedse, msmedci, msmedsub
- Median SE/CI: bpmed, bpmedse, medcipb
- Two-group comparisons: med2g, medhd2g, medpb, medpb2, dmedpb
- ANOVA: med1way, med2way, MEDanova, med2mcp
- Median effects: med.effect, MED.ES, medpb.es
- Specialized: medind, dlinmed, wwmed, wwwmed, medcurve, dmedian
- Ophthalmology median functions: oph.*.commedian, oph.*.comMedAE
- Utilities: WMW2med, runstest.med, exmed
- Successfully sourced and validated with all dependencies
- **Note**: Added missing pb2gen and pb2genMC to two-sample.R (dependency fix)

### Next Steps
1. Continue with remaining Week 3 modules (plotting.R, effect-size.R, etc.)
2. Extract remaining modules (power.R, winsorize.R, parallel.R, classification.R, special.R)
3. End-of-phase validation when all 20 modules extracted

### ✅ Phase 0: Preparation (COMPLETED - 2025-12-30)

All preparation tasks completed successfully:

1. **Safety Backup Created**
   - Original file backed up to: `/home/mando/coding/R-Projects/WRS/Rallfun-v45.R.ORIGINAL`
   - Size: 2.6 MB
   - Contains all 1,971 functions

2. **Working Directories Created**
   - `pkg/R-new/` - Staging area for new modular files
   - `pkg/tests/` - Test scripts for validation

3. **Function Inventory Generated**
   - File: `/home/mando/coding/R-Projects/WRS/all-functions.txt`
   - Total functions: **1,971**
   - Format: Sorted alphabetically, one function per line
   - Note: Duplicates visible (e.g., `adpchk` appears multiple times)

4. **Test Infrastructure Created**
   - `pkg/tests/create-reference-outputs.R` - Generates v0.45 baseline
   - `pkg/tests/test-backward-compat.R` - Validates compatibility
   - Reference data will be saved to: `reference-outputs.rds`

### 🔄 Current Phase: Phase 1, Week 3

**In Progress**: Specialized module extraction
**Next**: Extract ancova.R, regression.R, and remaining specialized modules

---

## Project Goals

### Critical Constraints
- ✅ **ZERO breaking changes** - All 1,971 functions remain available
- ✅ **Maintain signatures** - No changes to function parameters or defaults
- ✅ **Backward compatible** - All outputs must match v0.45 exactly
- ✅ **Keep exports** - Maintain `exportPattern("^[^\\.]")` in NAMESPACE

### Target State (v0.46)
- [ ] 20 focused module files (from 1 monolithic file)
- [ ] Full roxygen2 documentation (all ~1,500 user-facing functions)
- [ ] 0 redundant library() calls (currently 331)
- [ ] 0 duplicate functions (currently 62 duplicates)
- [ ] Clean namespace (471 internal helpers marked)

---

## Key Statistics

### Current Codebase Analysis

| Metric | Count | Notes |
|--------|-------|-------|
| Total Functions | 1,971 | Confirmed from inventory |
| User-facing Functions | ~1,500 | Estimated, need documentation |
| Internal Helpers | ~471 | Functions with `_sub`, `*MC`, etc. |
| Code Lines | 97,199 | Single file |
| File Size | 2.6 MB | Rallfun-v45.R |
| Redundant library() Calls | 331 | 84× parallel, 75× MASS, 21× akima |
| Duplicate Functions | 62 | Need resolution |
| Current Documentation | 1 | Only WRS-package.Rd |
| Old .Rd Files | 15 | In `old Rd files/` directory |

### Extracted Modules (12 of 20 Complete)

| Module | Status | Functions | Lines | Size | Key Functions |
|--------|--------|-----------|-------|------|---------------|
| 00-utils-core.R | ✅ | 54 | 2,836 | 72.7 KB | elimna, listm, matl, winvar, pbos, pbvar |
| location.R | ✅ | 74 | 3,242 | 86.5 KB | mest, mom, hd, tmean, onestep |
| outliers.R | ✅ | 64 | 2,713 | 67.4 KB | outpro, out, outbox, depth |
| bootstrap.R | ✅ | 30 | 915 | 21.9 KB | bootdpci, onesampb, trimcibt |
| two-sample.R | ✅ | 103 | 5,251 | 142 KB | yuen, yuend, wmw, pb2gen, cid |
| anova.R | ✅ | 57 | 3,547 | 99.8 KB | t1way, t2way, t3way, bwtrim |
| correlation.R | ✅ | 108 | 5,045 | 162.0 KB | pbcor, scor, tau, wincor, mscor |
| ancova.R | ✅ | 149 | 11,014 | 286.1 KB | ancova, Dancova, ancES, ancGLOB |
| regression.R | ✅ | 98 | 5,051 | 139.8 KB | tsreg, opreg, ltsreg, regci, reg1way |
| regression-advanced.R | ✅ | 75 | 3,433 | 95.2 KB | qhdsm, smean, logreg, mlrreg, KNNreg |
| covariance.R | ✅ | 43 | 1,034 | 27 KB | covogk, wincov, skipcov, dcov, covmba |
| mcp.R | ✅ | 106 | 9,046 | 255 KB | con1way, linconb, pairdepb, rmmcp, mcppb |
| medians.R | ✅ | 42 | 2,474 | 69 KB | msmed, med2g, medpb, MEDanova, med.effect |
| **TOTAL** | **13/20** | **1,003** | **55,601** | **1.6 MB** | **50.9% of functions** |

### Target Module Structure (20 Files)

| Module | Functions | Purpose |
|--------|-----------|---------|
| 00-utils-core.R | ~50 | Foundation utilities ✅ |
| location.R | ~80 | Robust location estimators ✅ |
| outliers.R | ~70 | Outlier detection ✅ |
| bootstrap.R | ~60 | Bootstrap infrastructure ✅ |
| two-sample.R | ~60 | Two-group comparisons ✅ |
| anova.R | ~90 | ANOVA methods ✅ |
| correlation.R | ~45 | Correlation methods ✅ |
| ancova.R | ~98 | ANCOVA methods (largest) ✅ |
| regression.R | ~85 | Regression methods ✅ |
| regression-advanced.R | ~60 | Quantile regression, etc. ✅ |
| covariance.R | ~50 | Covariance estimation ✅ |
| mcp.R | ~55 | Multiple comparisons ✅ |
| medians.R | ~40 | Median-based methods ✅ |
| plotting.R | ~50 | Visualization |
| effect-size.R | ~35 | Effect sizes |
| power.R | ~25 | Power analysis |
| winsorize.R | ~30 | Winsorization methods |
| parallel.R | ~80 | Multicore functions (*MC) |
| classification.R | ~40 | Classification/discriminant |
| special.R | ~80 | Specialized methods |
| zzz-internal.R | ~471 | Internal helpers |

---

## Phase Checklist

### ✅ Phase 0: Preparation (Week 0) - COMPLETED 2025-12-30
- [x] Create safety backup
- [x] Create working directories (R-new, tests)
- [x] Extract function inventory
- [x] Create reference test suite script
- [x] Create backward compatibility test script
- [x] **COMPLETED**: Run reference test suite to establish baseline
  - Created `reference-outputs.rds` (2.4KB, 23 test functions)
  - Baseline p-values and outputs captured successfully
- [x] **COMPLETED**: Analyze function dependencies
  - Created `dependency-analysis.rds` with full dependency graph
  - **Key finding**: Top 50 core utilities identified
  - **Critical**: `elimna` called by 928 functions (47%!)
  - **Must extract first**: elimna, listm, matl, near, hd, winvar
  - **Zero isolated functions**: All 1,971 functions are interconnected

### 🔄 Phase 1: Module Extraction (Weeks 1-3)
**Status**: In Progress - Starting Week 1
**Goal**: Split monolithic file into 20 modules WITHOUT changing code

#### ✅ Week 1: Foundation Modules (COMPLETED 2025-12-30)
- [x] Extract 00-utils-core.R (~50 functions)
  - **COMPLETED**: 50 functions extracted (2,765 lines, 72KB)
  - Top utilities: elimna (928 calls), listm (327 calls), matl (219 calls)
  - Successfully sourced without errors
- [x] Extract location.R (~80 functions)
  - **COMPLETED**: 75 functions extracted (3,230 lines, 87KB)
  - Includes: mest, mom, onestep, hd variants, trimmed means
  - Added dependencies: hpsi, lloc, tmean
- [x] Extract outliers.R (~70 functions)
  - **COMPLETED**: 64 functions extracted (2,713 lines, 68KB)
  - Includes: outpro variants, depth methods, bagging functions
- [x] Extract bootstrap.R (~60 functions)
  - **COMPLETED**: 30 functions extracted (915 lines, 22KB)
  - Includes: BCA methods, permutation tests, bootstrap infrastructure
- [x] Validation: All modules sourced successfully, core functions tested

#### ✅ Week 2: Main Analysis Modules (COMPLETED 2025-12-30)
- [x] Extract two-sample.R (~60 functions)
  - **COMPLETED**: 103 functions extracted (5,251 lines, 142 KB)
  - Includes: yuen, yuend, wmw, pb2gen, pb2genMC, cid, qcomhd, trimpb2
  - Bootstrap methods, quantile comparisons, effect sizes
  - Successfully sourced and tested
  - **NOTE**: Excluded mulwmwv2 (brace mismatch in original source)
  - **UPDATE 2025-12-30**: Added pb2gen and pb2genMC (were initially missed)
- [x] Extract anova.R (~90 functions)
  - **COMPLETED**: 57 functions extracted (3,547 lines, 99.8 KB)
  - Includes: t1way, t2way, t3way, bwtrim, rmanova, pbanova
  - One-way, two-way, three-way ANOVA methods
  - Bootstrap and robust variants
- [x] Extract correlation.R (~45 functions)
  - **COMPLETED**: 96 functions extracted (5,045 lines, 162.0 KB)
  - Includes: pbcor, pcor, scor, wincor, tau, mscor
  - Pearson, Spearman, Kendall correlations
  - Percentage bend, skipped correlations
  - **Added**: pbos function to 00-utils-core.R (dependency)
- [x] Validation: All modules sourced successfully, all tests passed

#### Week 3: Specialized Modules (IN PROGRESS)
- [x] Extract ancova.R (~98 functions)
  - **COMPLETED**: 149 functions extracted (11,014 lines, 287 KB)
  - Includes: ancova, Dancova, ancES, ancGLOB, ancJN, anctspb
  - Core ANCOVA, dependent ANCOVA, bootstrap methods, effect sizes
  - Successfully sourced with dependencies
  - **Note**: Some advanced functions have dependencies that will be resolved as more modules are extracted
- [x] Extract regression.R (~85 functions)
  - **COMPLETED**: 98 functions extracted (5,051 lines, 139.8 KB)
  - Includes: tsreg, tshdreg, opreg, ltsreg, qreg, regci, regtest, reg1way
  - Theil-Sen, LTS, M-regression, outlier-pruned, inference methods
  - Two-group comparisons and one-way regression ANOVA
  - Correlation-based regression: scorreg, correg, taureg
  - Successfully sourced and validated
- [x] Extract regression-advanced.R (~60 functions)
  - **COMPLETED**: 75 functions extracted (3,433 lines, 95.2 KB)
  - Includes: qhdsm, qhdsm2g, smean, smeancr, logreg, mlrreg, mulgreg
  - Quantile regression smoothers, smoothing methods, logistic regression
  - Multivariate/multilevel regression, KNN, random forest
  - GAM-related methods, regression inference (regYci, regYband)
  - Mediation, PCA, instrumental variables regression
  - Successfully sourced and validated
- [x] Extract covariance.R (~50 functions)
  - **COMPLETED**: 43 functions extracted (972 lines, 26.5 KB)
  - Includes: covogk, wincov, skipcov, dcov, covmba, covmtrim
  - Robust covariance: OGK, MVE, MCD, MBA, S-estimators
  - Winsorized/trimmed covariance, skipped covariance
  - Median-based, distance, and ROC-based covariance
  - Mixed design covariance (bwwcovm, bbwcovm)
  - Successfully sourced and validated
- [x] Extract mcp.R (~55 functions)
  - **COMPLETED**: 106 functions extracted (9,046 lines, 255 KB)
  - Includes: con1way, con2way, con3way (contrast generators)
  - Linear contrasts: linconb, linconpb, linconbt, linconEP, linconES, linconQS
  - Dependent contrasts: lindep, lindepbt, pairdepb
  - Factorial MCP: mcp2a, mcp2atm, mcp3atm, mcp3med, rm3mcp
  - Repeated measures: rmmcp, wmcp, rmmcppb, rmmcpES, rmmcpQS, rmmismcp
  - Between-within: bwmcp, bwwmcp, bbwmcp, bwrmcp, bwimcp, bwbmcp
  - Bootstrap MCP: mcppb, tmcppb, bmcppb, pbmcp, bbmcppb, wwmcppb
  - Split-plot: spmcpa, spmcpi, spmcpb, sppba, sppbb, sppbi
  - Quantile-based: qdmcp, qdmcpdif, mwwmcp, twwmcp, tkmcp
  - Specialized: stepmcp, signmcp, discmcp, sintmcp, anctsmcp, skmcp
  - P-value adjustment: mcpPV, mcpKadjp
  - Successfully sourced and validated
- [x] Extract medians.R
  - **COMPLETED**: 42 functions extracted (2,474 lines, 69 KB)
  - Includes: msmed, msmedse, msmedci, msmedsub
  - Two-group comparisons: med2g, medhd2g, medpb, medpb2, dmedpb
  - ANOVA: med1way, med2way, MEDanova, med2mcp
  - Effect sizes: med.effect, MED.ES, medpb.es
  - Specialized: medind, dlinmed, wwmed, medcurve, dmedian
  - Ophthalmology: oph.*.commedian, oph.*.comMedAE
  - Successfully sourced and validated
  - **Note**: Fixed missing pb2gen/pb2genMC in two-sample.R
- [ ] Extract plotting.R, effect-size.R, power.R, winsorize.R, parallel.R
- [ ] Extract classification.R, special.R
- [ ] Extract zzz-internal.R (last)
- [ ] **End-of-phase validation**: Replace R/ with R-new/, build & check package

### 📋 Phase 2: Optimization (Weeks 4-5)
**Status**: Not started
**Goal**: Remove redundant code, fix duplicates, optimize patterns

#### Week 4: Eliminate library() Calls
- [ ] Update NAMESPACE (add imports)
- [ ] Update DESCRIPTION (move to Imports)
- [ ] Remove library(parallel) calls (84 instances)
- [ ] Remove library(MASS) calls (75 instances)
- [ ] Remove library(akima) calls (21 instances)
- [ ] Handle optional packages (conditional loading)
- [ ] Validate after each module

#### Week 5: Fix Duplicates & Extract Patterns
- [ ] Resolve 62 duplicate functions
  - Document in `pkg/duplicate-functions-resolution.md`
  - Example: lintestMC at lines 426 & 33604
- [ ] Extract .bootstrap_setup() pattern (456 occurrences)
- [ ] Extract .bootstrap_ci() pattern (90+ occurrences)
- [ ] Extract .bootstrap_pvalue() pattern (70+ occurrences)
- [ ] Extract .prepare_regression_data() pattern
- [ ] Extract .apply_outlier_removal() pattern
- [ ] Refactor top 20 functions to use utilities
- [ ] Validate outputs match exactly

### 📋 Phase 3: Documentation - Core (Weeks 6-8)
**Status**: Not started
**Goal**: Document all user-facing functions in core modules

- [ ] Setup: Add Roxygen to DESCRIPTION
- [ ] Create common-params.R (shared param docs)
- [ ] Week 6: Document utils-core, location, outliers (~180 functions)
- [ ] Week 7: Document two-sample, anova, regression, correlation (~280 functions)
- [ ] Week 8: Document ancova, bootstrap, mcp, covariance (~280 functions)

### 📋 Phase 4: Documentation - Advanced (Weeks 9-11)
**Status**: Not started

- [ ] Week 9: Document advanced modules (~370 functions)
- [ ] Week 10: Document internal helpers (~471 functions)
- [ ] Week 11: Review, polish, create vignettes

### 📋 Phase 5: Final Testing (Week 12)
**Status**: Not started

- [ ] Comprehensive testing
- [ ] Performance validation
- [ ] Export validation
- [ ] Update package metadata (v0.46)
- [ ] Final checklist (see below)

---

## Important Files & Locations

### Source Files
```
/home/mando/coding/R-Projects/WRS/
├── Rallfun-v45.R.ORIGINAL          # Backup - DO NOT MODIFY
├── REFACTORING-PROGRESS.md         # This file
├── all-functions.txt               # Function inventory (1,971 functions)
├── reference-outputs.rds           # Baseline test data (to be created)
├── .claude/plans/
│   └── curious-questing-clock.md   # Approved implementation plan
├── pkg/
│   ├── R/
│   │   └── Rallfun-v45.R          # Current monolithic file
│   ├── R-new/                      # Staging area for new modules
│   ├── tests/
│   │   ├── create-reference-outputs.R
│   │   └── test-backward-compat.R
│   ├── NAMESPACE                   # To be updated
│   └── DESCRIPTION                 # To be updated
└── old Rd files/                   # 15 reference .Rd files
```

### Files to Be Created (Phase 1)
```
pkg/R-new/
├── 00-utils-core.R
├── location.R
├── outliers.R
├── bootstrap.R
├── two-sample.R
├── anova.R
├── ancova.R
├── regression.R
├── regression-advanced.R
├── correlation.R
├── covariance.R
├── mcp.R
├── medians.R
├── plotting.R
├── effect-size.R
├── power.R
├── winsorize.R
├── parallel.R
├── classification.R
├── special.R
├── common-params.R                 # Phase 3
└── zzz-internal.R
```

---

## Test Strategy

### Before Making Changes
1. **Run reference test suite**:
   ```r
   source("pkg/tests/create-reference-outputs.R")
   # Creates: reference-outputs.rds
   ```

### After Each Change
2. **Run compatibility tests**:
   ```r
   source("pkg/tests/test-backward-compat.R")
   test_backward_compatibility()
   # Must pass with 0 failures
   ```

### After Each Module
3. **Build validation**:
   ```bash
   R CMD build pkg/
   R CMD check --as-cran WRS_*.tar.gz
   # Target: 0 errors, 0 warnings
   ```

---

## Key Technical Changes

### NAMESPACE Updates (Phase 2)
**Current**:
```r
exportPattern("^[^\\.]")
importFrom("grDevices", "chull")
importFrom("graphics", "abline", ...)
importFrom("stats", "TukeyHSD", ...)
```

**Target** (maintain exportPattern, add imports):
```r
exportPattern("^[^\\.]")  # KEEP - backward compatibility

# Add package-level imports
import(parallel)
importFrom(MASS, ginv, cov.rob, cov.mcd, cov.mve, lqs)
importFrom(akima, interp, aspline)
importFrom(quantreg, rq, rq.fit)
# ... more imports
```

### DESCRIPTION Updates (Phase 2)
**Move to Imports**:
- parallel (currently Suggested)
- MASS (currently Suggested)
- akima (currently Suggested)

**Add**:
```
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.2.3
```

---

## Known Issues & Findings

### Duplicate Functions (62 total)
- `lintestMC`: Lines 426 and 33604
- `adpchk`: Multiple definitions
- Need systematic comparison and resolution
- Document in `pkg/duplicate-functions-resolution.md`

### Redundant Library Calls (331 total)
- `library(parallel)`: 84 occurrences
- `library(MASS)`: 75 occurrences
- `library(akima)`: 21 occurrences
- 40+ other packages

### Common Patterns to Extract
1. Bootstrap setup: ~456 occurrences
2. CI computation: ~90 occurrences
3. P-value computation: ~70 occurrences
4. Data preparation: ~100 occurrences
5. Outlier removal: ~80 occurrences

---

## Success Criteria

### Must Achieve (Final Validation)
- [ ] 0 errors, 0 warnings from `R CMD check --as-cran`
- [ ] All 1,971 functions present and functional
- [ ] 100% reference test pass rate
- [ ] Same namespace exports as v0.45
- [ ] All function signatures unchanged
- [ ] Performance within 10% of v0.45

### Quality Targets
- [ ] ~1,500 functions with full roxygen2 docs
- [ ] ~471 internal helpers marked `@keywords internal`
- [ ] 0 library() calls in function bodies
- [ ] 0 duplicate functions
- [ ] 20 focused modules with clear responsibilities

---

## Daily Workflow Reminders

### During Work
1. Work on one module/section at a time
2. After significant changes, run compatibility tests
3. Commit frequently with clear messages

### End of Day
1. Run full package check:
   ```bash
   R CMD build pkg/
   R CMD check --as-cran WRS_*.tar.gz
   ```
2. Update this progress file
3. Commit and document any issues

### Weekly
1. Review all changes from the week
2. Performance benchmark on key functions
3. Update documentation if patterns emerged
4. Plan next week's targets

---

## Commands Reference

### Function Inventory
```bash
# List all functions
grep "^[a-zA-Z][a-zA-Z0-9_\.]*[[:space:]]*<-[[:space:]]*function" \
  pkg/R/Rallfun-v45.R | awk -F'<-' '{print $1}' | sort

# Count functions
grep -c "^[a-zA-Z][a-zA-Z0-9_\.]*[[:space:]]*<-[[:space:]]*function" \
  pkg/R/Rallfun-v45.R
```

### Find Duplicates
```bash
# Find duplicate function names
grep "^[a-zA-Z][a-zA-Z0-9_\.]*[[:space:]]*<-[[:space:]]*function" \
  pkg/R/Rallfun-v45.R | awk -F'<-' '{print $1}' | sort | uniq -d
```

### Find Library Calls
```bash
# Count library() calls
grep -c "library(" pkg/R/Rallfun-v45.R

# Find specific library calls
grep -n "library(parallel)" pkg/R/Rallfun-v45.R | wc -l
```

### Compare Function Definitions
```bash
# Find line numbers for a function
grep -n "^lintestMC<-function" pkg/R/Rallfun-v45.R

# Compare two definitions
diff <(sed -n '426,500p' pkg/R/Rallfun-v45.R) \
     <(sed -n '33604,33678p' pkg/R/Rallfun-v45.R)
```

---

## Notes & Observations

### Textbook Reference
- User has access to Wilcox's "Introduction to Robust Estimation and Hypothesis Testing" (5th Ed.)
- Can extract examples for documentation
- Map functions to chapters for accurate examples

### Old Documentation
- 15 .Rd files exist in `old Rd files/` directory
- Shows documentation structure to emulate
- Examples: yuen.Rd, pbcor.Rd, tsreg.Rd, etc.

### Package Status
- Currently at v0.45 (2025-10-29)
- Marked as "no longer maintained"
- Our refactoring creates maintained v0.46 fork

---

## Questions & Decisions Log

### 2025-12-30: Initial Planning
**Q**: Scope and approach?
**A**: Improve in-place, split into 15-20 modules, document all ~1,500 functions, has textbook

**Q**: Before starting implementation?
**A**: Create progress tracking document (this file)

---

## Next Session Checklist

When resuming work:
1. ✅ Read this file (REFACTORING-PROGRESS.md)
2. ✅ Check Current Status section above
3. ✅ Review most recent Phase checklist
4. ✅ Check Questions & Decisions Log
5. ✅ Review any uncommitted changes: `git status`
6. ✅ Continue from last incomplete task

---

## Contact & Resources

- **Full Plan**: `.claude/plans/curious-questing-clock.md`
- **Original README**: `README.md`
- **Package Info**: `pkg/DESCRIPTION`
- **Textbook**: Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing, 5th Ed.
- **GitHub**: https://github.com/nicebread/WRS/

---

*Last updated: 2025-12-30 - Phase 1, Week 3 In Progress (1,003 functions extracted across 13 modules, medians.R completed)*
