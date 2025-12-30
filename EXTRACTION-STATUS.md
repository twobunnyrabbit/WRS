# WRS Package Refactoring - Extraction Status Report

**Generated:** 2025-12-30
**Phase:** Phase 1, Week 3 - Module Extraction
**Status:** IN PROGRESS

---

## Executive Summary

- **9 of 20 modules** completed (45%)
- **737 of 1,971 functions** extracted (37.4%)
- **39,549 of 97,199 lines** extracted (40.7%)
- **1,074.7 KB of 2,663 KB** extracted (40.4%)
- **All 9 modules source successfully** ✅

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
| **TOTAL** | **9 modules** | **737** | **39,549** | **1,074.7 KB** | **45%** |

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

---

## Remaining Work

### Modules to Extract (11 remaining)

| # | Module | Est. Functions | Purpose |
|---|--------|----------------|---------|
| 10 | regression-advanced.R | ~60 | Advanced regression (mediation, logistic, Y-value) |
| 11 | covariance.R | ~50 | Covariance estimation methods |
| 12 | mcp.R | ~55 | Multiple comparison procedures |
| 13 | medians.R | ~40 | Median-based methods |
| 14 | plotting.R | ~50 | Visualization functions |
| 15 | effect-size.R | ~35 | Effect size calculations |
| 16 | power.R | ~25 | Power analysis |
| 17 | winsorize.R | ~30 | Winsorization methods |
| 18 | parallel.R | ~80 | Multicore functions (*MC) |
| 19 | classification.R | ~40 | Classification/discriminant |
| 20 | special.R | ~80 | Specialized methods |
| 21 | zzz-internal.R | ~471 | Internal helper functions |

**Total estimated:** ~1,016 functions (note: actual may vary)

### Remaining Statistics
- **Functions:** 1,234 remaining (62.6%)
- **Lines:** 57,650 remaining (59.3%)
- **Size:** 1,588.3 KB remaining (59.6%)

---

## Validation Status

### ✅ Completed
- All 9 modules source without errors
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
- `pkg/R-new/*.R` - 9 extracted module files
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
2. ⬜ Extract regression-advanced.R (~60 functions)
3. ⬜ Extract covariance.R (~50 functions)
4. ⬜ Extract mcp.R (~55 functions)
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
- ✅ 45% of modules extracted
- ✅ 37.4% of functions extracted
- ✅ 40.7% of lines extracted
- ✅ All modules source successfully
- ✅ Week 1 and Week 2 completed on schedule
- ✅ Week 3 in progress, on track

---

*Last updated: 2025-12-30*
*See REFACTORING-PROGRESS.md for detailed implementation notes*
