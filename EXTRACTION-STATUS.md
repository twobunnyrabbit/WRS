# WRS Package Refactoring - Extraction Status Report

**Generated:** 2025-12-30
**Phase:** Phase 1, Week 3 - Module Extraction
**Status:** IN PROGRESS

---

## Executive Summary

- **12 of 20 modules** completed (60%)
- **961 of 1,971 functions** extracted (48.8%)
- **53,062 of 97,199 lines** extracted (54.6%)
- **1.5 MB of 2.6 MB** extracted (57.7%)
- **All 12 modules source successfully** ✅

---

## Extracted Modules

| # | Module | Functions | Lines | Size | Status |
|---|--------|-----------|-------|------|--------|
| 1 | 00-utils-core.R | 54 | 2,836 | 72.7 KB | ✅ |
| 2 | location.R | 74 | 3,242 | 86.5 KB | ✅ |
| 3 | outliers.R | 64 | 2,713 | 67.4 KB | ✅ |
| 4 | bootstrap.R | 30 | 915 | 21.9 KB | ✅ |
| 5 | two-sample.R | 103 | 5,186 | 138.5 KB | ✅ |
| 6 | anova.R | 57 | 3,547 | 99.8 KB | ✅ |
| 7 | correlation.R | 108 | 5,045 | 162.0 KB | ✅ |
| 8 | ancova.R | 149 | 11,014 | 286.1 KB | ✅ |
| 9 | regression.R | 98 | 5,051 | 139.8 KB | ✅ |
| 10 | regression-advanced.R | 75 | 3,433 | 95.2 KB | ✅ |
| 11 | covariance.R | 43 | 1,034 | 27 KB | ✅ |
| 12 | mcp.R | 106 | 9,046 | 255 KB | ✅ |
| **TOTAL** | **12 modules** | **961** | **53,062** | **1.5 MB** | **60%** |

---

## Key Functions by Module

### 00-utils-core.R (54 functions)
Core utilities used throughout the package:
- `elimna`, `listm`, `matl` (called by 928, 327, 219 functions respectively)
- `near`, `winvar`, `pbvar`, `pbos`
- `winval`, `rdepth`, `fdepth`

### location.R (74 functions)
Robust location estimators:
- `mest`, `mom`, `onestep`, `tmean`
- `hd` (Harrell-Davis), `tauloc`, `trimpb`
- Variants: `hdpb`, `hdci`, `momci`, `mestci`

### outliers.R (64 functions)
Outlier detection methods:
- `outpro`, `outproMC`, `out`, `outbox`
- `depth`, `fdepth`, `rdepth`
- Projection-based methods, bagging

### bootstrap.R (30 functions)
Bootstrap infrastructure:
- `bootdpci`, `onesampb`, `trimcibt`
- `percentile_boot`, `Qhd`, `sintv2`

### two-sample.R (103 functions)
Two-group comparisons:
- `yuen`, `yuend`, `wmw`, `pb2gen`
- `cid`, `qcomhd`, `trimpb2`
- Bootstrap variants, effect sizes

### anova.R (57 functions)
ANOVA methods:
- `t1way`, `t2way`, `t3way`, `bwtrim`
- `rmanova`, `pbanova`, `sppba`, `sppbi`
- Robust one-way, two-way, three-way ANOVA

### correlation.R (108 functions)
Correlation methods:
- `pbcor`, `pcor`, `scor`, `wincor`, `tau`
- `mscor`, `indt`, `twocor`
- Percentage bend, skipped, Winsorized correlations

### ancova.R (149 functions)
ANCOVA methods:
- `ancova`, `Dancova`, `ancES`, `ancGLOB`
- `ancJN`, `anctspb`, `ancdet`, `ancpb`
- Independent and dependent group ANCOVA variants

### regression.R (98 functions)
Core regression methods:
- `tsreg`, `tshdreg`, `opreg`, `ltsreg`, `qreg`
- `regci`, `regtest`, `lintest`, `lintestMC`
- `difreg`, `reg2ci`, `DregG`, `reg1way`
- Theil-Sen, LTS, M-regression, outlier-pruned methods

### regression-advanced.R (75 functions)
Advanced regression methods:
- `qhdsm`, `qhdsm2g`, `smean`, `smeancr` (quantile smoothers, smoothing)
- `logreg`, `logreg.P.ci`, `logreg.pred` (logistic regression)
- `mlrreg`, `mulgreg`, `MULMreg` (multivariate/multilevel)
- `KNNreg`, `quantregForest`, `regR.Forest` (machine learning)
- `gamindt`, `gamplot` (GAM-related methods)
- `regYci`, `regYband`, `regmediate`, `regpca` (specialized)

### covariance.R (43 functions)
Covariance and scatter matrix estimation:
- `covogk`, `outogk`, `skipogk`, `cov.ogk` (OGK robust estimators)
- `covmve`, `mvecov`, `covmcd`, `mcdcov`, `DETMCD` (MVE/MCD estimators)
- `covmba`, `cov.mba`, `covmba2`, `rmba` (MBA estimators)
- `wincov`, `wincovN`, `covmtrim` (winsorized/trimmed)
- `skipcov`, `skip.cov` (skipped covariance)
- `dcov`, `Scov`, `tbscov` (distance, S-estimator, tau-scale)
- `cov2med`, `covmmed`, `covroc`, `cov.roc` (median-based, ROC)
- `bwwcovm`, `bbwcovm` (mixed designs)
- `cov2cor`, `longcov2mat`, `covloc` (utilities)

### mcp.R (106 functions)
Multiple Comparisons Procedures:
- `con1way`, `con2way`, `con3way` (contrast matrix generators)
- `linconb`, `linconpb`, `linconbt`, `linconEP`, `linconES`, `linconQS` (linear contrasts)
- `lindep`, `lindepbt`, `pairdepb` (dependent contrasts)
- `mcppb`, `tmcppb`, `bmcppb`, `pbmcp` (one-way MCP)
- `mcp2a`, `mcp2atm`, `mcp3atm`, `mcp3med` (factorial MCP)
- `rmmcp`, `wmcp`, `rmmcppb`, `rmmcpES`, `rmmcpQS` (repeated measures)
- `bwmcp`, `bwwmcp`, `bbwmcp`, `bwrmcp`, `bwimcp` (between-within designs)
- `spmcpa`, `spmcpi`, `spmcpb`, `sppba`, `sppbb`, `sppbi` (split-plot)
- `qdmcp`, `mwwmcp`, `twwmcp`, `tkmcp` (quantile-based)
- `mcpPV`, `mcpKadjp` (p-value adjustment)

---

## Remaining Work

### Modules to Extract (8 remaining)

| # | Module | Est. Functions | Purpose |
|---|--------|----------------|---------|
| 13 | medians.R | ~40 | Median-based methods |
| 14 | plotting.R | ~50 | Visualization functions |
| 15 | effect-size.R | ~35 | Effect size calculations |
| 16 | power.R | ~25 | Power analysis |
| 17 | winsorize.R | ~30 | Winsorization methods |
| 18 | parallel.R | ~80 | Multicore functions (*MC) |
| 19 | classification.R | ~40 | Classification/discriminant |
| 20 | special.R | ~80 | Specialized methods |
| 21 | zzz-internal.R | ~471 | Internal helper functions |

**Total estimated:** ~1,010 functions (note: actual may vary)

### Remaining Statistics
- **Functions:** 1,010 remaining (51.2%)
- **Lines:** 44,137 remaining (45.4%)
- **Size:** 1.1 MB remaining (42.3%)

---

## Validation Status

### ✅ Completed
- All 12 modules source without errors
- No syntax errors detected
- Key functions tested and available
- Dependencies properly resolved

### ⚠️ Pending
- Full backward compatibility testing
- Package build (`R CMD build`)
- Package check (`R CMD check --as-cran`)
- Performance benchmarking

---

## Known Issues

### Duplicate Functions Found
1. **lintestMC** - Lines 426 and 33604 (extracted first definition)
2. **Qreghat** - Lines 86686 and 93871 (extracted first definition)
3. **scorreg.sub** - Lines 65557 and 67269 (extracted first definition)

**Resolution:** Using first definition in each case. Need to investigate and resolve during Phase 2 (Optimization).

### Missing Functions
- **mulwmwv2** - Excluded from two-sample.R due to brace mismatch in original source
  - Location: Line ~unknown in original file
  - Action: Will investigate and fix separately

---

## Files Generated

### Core Refactoring Files
- `pkg/R-new/*.R` - 10 extracted module files
- `Rallfun-v45.R.ORIGINAL` - Backup of original file (2.6 MB)
- `all-functions.txt` - Inventory of all 1,971 functions
- `reference-outputs.rds` - Baseline test data (2.4 KB)
- `dependency-analysis.rds` - Dependency graph

### Progress Tracking
- `REFACTORING-PROGRESS.md` - Main progress document
- `EXTRACTION-STATUS.md` - This file
- `CLAUDE.md` - Project instructions for Claude

### Test Infrastructure
- `pkg/tests/create-reference-outputs.R` - Reference test generator
- `pkg/tests/test-backward-compat.R` - Compatibility validator

### Temporary/Analysis Files
- `/tmp/regression-functions.txt` - Function list for extraction
- `/tmp/regression-functions-lines.csv` - Line numbers for extraction
- `/tmp/regression-extraction-stats.csv` - Detailed extraction statistics
- `/tmp/current-module-stats.csv` - Current module statistics

---

## Next Steps

### Immediate (Week 3 continuation)
1. ✅ regression.R extracted - COMPLETED
2. ✅ regression-advanced.R extracted - COMPLETED (75 functions)
3. ✅ covariance.R extracted - COMPLETED (43 functions)
4. ✅ mcp.R extracted - COMPLETED (106 functions)
5. ⬜ Extract medians.R (~40 functions)

### Week 3 Completion
- Extract remaining specialized modules (plotting, effect-size, power, etc.)
- End-of-week validation: source all modules together

### Phase 1 Completion (End of Week 3)
- Extract all 20 modules
- Replace `pkg/R/` with `pkg/R-new/` contents
- Run full package build and check
- Validate 100% backward compatibility

### Phase 2 (Weeks 4-5)
- Eliminate redundant library() calls (331 instances)
- Resolve duplicate functions (62 duplicates)
- Extract common patterns (bootstrap, CI, p-value)

### Phase 3-4 (Weeks 6-11)
- Add roxygen2 documentation (~1,500 functions)
- Create vignettes and examples

### Phase 5 (Week 12)
- Final testing and validation
- Release v0.46

---

## Success Criteria

### Must Achieve
- ✅ 0 errors when sourcing all modules
- ⬜ 0 errors from `R CMD check --as-cran`
- ⬜ All 1,971 functions present and functional
- ⬜ 100% backward compatibility
- ⬜ Same namespace exports as v0.45
- ⬜ Performance within 10% of v0.45

### On Track
- ✅ 60% of modules extracted (OVER HALFWAY!)
- ✅ 48.8% of functions extracted
- ✅ 54.6% of lines extracted
- ✅ All modules source successfully
- ✅ Week 1 and Week 2 completed on schedule
- ✅ Week 3 in progress, on track

---

*Last updated: 2025-12-30*
*See REFACTORING-PROGRESS.md for detailed implementation notes*
