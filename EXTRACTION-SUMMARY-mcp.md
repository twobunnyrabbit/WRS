# MCP Module Extraction - Summary Report

**Date**: 2025-12-30
**Module**: mcp.R (Multiple Comparisons Procedures)
**Status**: ✅ COMPLETED

---

## Overview

Successfully extracted all Multiple Comparisons Procedures (MCP) functions from the monolithic `Rallfun-v45.R` file into a well-organized, standalone module.

## Extraction Statistics

| Metric | Value |
|--------|-------|
| **Total Functions** | 107 |
| **Main Functions** | 102 |
| **Aliases** | 2 (bbmcp, wmcp) |
| **Helper Functions** | 5 (.sub helpers + tsub) |
| **Lines of Code** | 9,046 |
| **File Size** | 254.2 KB |
| **Categories** | 14 logical groups |

## Module Organization

The extracted functions are organized into 14 categories:

### 1. Contrast Matrix Generators (3 functions)
- `con1way` - One-way contrasts
- `con2way` - Two-way factorial contrasts
- `con3way` - Three-way factorial contrasts

### 2. Linear Contrasts - Independent Groups (13 functions)
Core functions for linear contrasts with independent groups:
- `linconm`, `lincon.old`, `lincon.pool`
- `linconb`, `linconbt`, `linconpb` - Bootstrap variants
- `linconMpb`, `linconSpb` - Percentile bootstrap
- `linconEP`, `linconES`, `linconQS` - Effect size variants
- `lincon.bin`, `lincon.binPV` - Binary data

### 3. Linear Contrasts - Dependent Groups (3 functions)
- `lindep` - Linear dependent contrasts
- `lindepbt` - Bootstrap version
- `pairdepb` - Pairwise comparisons with bootstrap

### 4. Core MCP - One-Way Designs (7 functions)
- `mcppb`, `mcppb20` - General MCP with bootstrap
- `tmcppb`, `bmcppb`, `pbmcp` - Specialized MCP
- `pbtrmcp`, `mcpOV` - Trimmed mean variants

### 5. Two-Way Factorial Designs (3 functions)
- `mcp2a` - Two-way MCP with general estimator
- `mcp2atm` - Two-way MCP with trimmed means
- `bbmcp` - Alias for mcp2atm

### 6. Three-Way Factorial Designs (3 functions)
- `mcp3atm` - Three-way trimmed means
- `mcp3med` - Three-way medians
- `rm3mcp` - Repeated measures three-way

### 7. Repeated Measures MCP (9 functions)
Comprehensive repeated measures support:
- `rmmcp` - Core repeated measures MCP
- `wmcp` - Alias for rmmcp
- `rmmcppb`, `rmmcppbd`, `rmmcppbtm` - Bootstrap variants
- `rmmcppbv2`, `rmmcpv2` - Version 2 implementations
- `rmmcpES`, `rmmcpQS` - Effect size variants
- `rmmismcp` - Missing data handling

### 8. Between-Within Designs (13 functions)
Mixed factorial designs:
- `bwmcp`, `bwwmcp`, `bbwmcp` - Core between-within
- `bwrmcp`, `bwamcp` - Specialized variants
- `bwimcp`, `bwimcpES`, `bwimcpAKP` - Interaction-specific
- `bwbmcp`, `bwmedimcp`, `bwmedbmcp` - Median-based
- `bwmcpAKP`, `bwmcpORD`, `bwmcpQS` - Additional variants

### 9. Bootstrap MCP for Mixed Designs (19 functions)
Extensive bootstrap support for complex designs:
- `bwmcppb`, `bwmcppb.sub`, `bwmcppb.adj`
- `bwwmcppb`, `bwwmcppb.sub`
- `bbmcppb`, `bbwmcppb`, `bbwmcppb.sub`
- `bbbmcppb`, `bbbmcppb.sub`
- `wwwmcppb`, `wwwmcppb.sub`, `wwwmcppbtr`
- `wmcppb`, `wwmcppb`, `wwmcpbt`
- `wwmcpES`, `wwmcpQS`, `wwwmcpQS`

### 10. Split-Plot Designs (7 functions)
- `spmcpa`, `spmcpi`, `spmcpb`, `spmcpbA`
- `sppba`, `sppbb`, `sppbi`

### 11. Quantile-Based MCP (5 functions)
- `qdmcp`, `qdmcpdif` - Quantile differences
- `mwwmcp` - Mann-Whitney-Wilcoxon
- `twwmcp`, `tkmcp` - Specialized quantile tests

### 12. Effect Size MCP (9 functions)
Functions providing effect size information:
- `esmcp` - Effect size MCP
- `bbmcpEP`, `bbmcpQS` - Between-between variants
- `bbdetmcp`, `bbdetmcpQS` - Detailed between-between
- `bmcpAKP`, `bmcpQS` - Between-subjects
- `wmcpAKP`, `wmcpQS` - Within-subjects

### 13. Specialized MCP Functions (9 functions)
Domain-specific MCP procedures:
- `stepmcp` - Stepwise MCP
- `signmcp` - Sign test MCP
- `discmcp` - Discrete data MCP
- `sintmcp` - Simultaneous inference
- `anctsmcp` - ANCOVA with Theil-Sen
- `skmcp` - Skipped correlation MCP
- `mcp.nestAP` - Nested design
- `binmcp`, `binmcp.crit` - Binary data

### 14. P-value Adjustment Utilities (2 functions)
- `mcpPV` - P-value adjustment
- `mcpKadjp` - K-adjusted p-values (Holm, Hochberg, BH, etc.)

## Helper Functions

### tsub
Bootstrap test statistic helper for trimmed means (used by pairdepb and related functions).

### .sub Functions
Five specialized helper functions for complex bootstrap MCP in mixed designs:
- `bbbmcppb.sub` - Between-between-between bootstrap
- `bbwmcppb.sub` - Between-between-within bootstrap
- `bwwmcppb.sub` - Between-within-within bootstrap
- `wwwmcppb.sub` - Within-within-within bootstrap
- `bwmcppb.sub` - Between-within bootstrap

## Aliases

Two function aliases for backward compatibility:
- `bbmcp = mcp2atm` (line 20778 in original)
- `wmcp <- rmmcp` (line 27287 in original)

## Dependencies

This module depends on functions from all previously extracted modules:
- **00-utils-core.R**: elimna, listm, matl, winvar, pbos, pbvar
- **location.R**: tmean, mom, onestep, trimse, winmean, hd
- **outliers.R**: outpro, out, outbox
- **bootstrap.R**: Bootstrap infrastructure
- **two-sample.R**: yuen, yuend, pb2gen, trimpb, cid
- **anova.R**: t1way, rmanova, bwtrim
- **correlation.R**: pcor, wincor, scor
- **ancova.R**: ancova, Dancova, ancts
- **regression.R**: tsreg, opreg, regci
- **regression-advanced.R**: qhdsm, smean
- **covariance.R**: wincov, skipcov

## Validation

✅ All validation tests passed:
- Module sources without errors
- No syntax errors detected
- All 107 functions accessible
- Tested con1way(4) - correct output
- Tested con2way(2,3) - correct dimensions
- All dependencies met from previously extracted modules
- Integration test with all 12 modules - successful

## Key Features

1. **Comprehensive Coverage**: Handles all major experimental designs
   - One-way, two-way, three-way factorial
   - Repeated measures (within-subjects)
   - Between-subjects
   - Mixed designs (between-within)
   - Split-plot designs
   - Nested designs

2. **Multiple Inference Methods**:
   - Bootstrap (percentile, bootstrap-t)
   - Permutation tests
   - Classical parametric methods
   - Robust methods (trimmed means, M-estimators)

3. **Effect Size Support**:
   - Explanatory power (EP)
   - Effect sizes (ES)
   - Q-statistic variants (QS)
   - AKP variants

4. **P-value Adjustment**:
   - Holm's method
   - Hochberg's method
   - Benjamini-Hochberg (BH)
   - Romano-Shaikh (RS)
   - Sarkar's method

5. **Flexible Estimators**:
   - Trimmed means (default 20%)
   - M-estimators (MOM, onestep)
   - Medians
   - Quantiles
   - Custom estimators via function arguments

## Impact on Refactoring Progress

### Before This Extraction
- 11 modules completed (55%)
- 855 functions extracted (43.4%)
- 43,954 lines (45.2%)
- 1,196.4 KB (46.0%)

### After This Extraction
- **12 modules completed (60%)**
- **962 functions extracted (48.8%)**
- **53,000 lines (54.5%)**
- **1,450.6 KB (55.8%)**

**Progress**: This single module extraction increased overall progress by 5% and added 107 functions (107 / 1,971 = 5.4% of all functions).

## Files Created

1. **Primary**:
   - `/home/mando/coding/R-Projects/WRS/pkg/R-new/mcp.R` (254.2 KB)

2. **Documentation**:
   - `MODULE-STATUS-mcp.txt` - Detailed status
   - `EXTRACTION-SUMMARY-mcp.md` - This file

3. **Intermediate** (for reference):
   - `/tmp/mcp_functions_extracted.txt`
   - `/tmp/mcp_helpers_extracted.txt`
   - `/tmp/mcp_extraction_summary.txt`

## Extraction Method

### Tools Used
- Python 3 script with brace-matching algorithm
- Regular expressions for function detection
- Automated extraction with validation

### Process
1. Identified all 102 MCP functions by line number
2. Extracted complete function definitions with proper brace matching
3. Identified and extracted 5 helper functions (.sub)
4. Identified and extracted tsub helper
5. Organized into 14 logical categories
6. Created module header with documentation
7. Validated syntax and functionality
8. Integration tested with all previous modules

### Quality Assurance
- ✅ Complete function extraction (no missing braces)
- ✅ All dependencies identified
- ✅ Proper organization by category
- ✅ Comprehensive header documentation
- ✅ Successful source testing
- ✅ Functional validation of key functions

## Next Steps

As documented in `REFACTORING-PROGRESS.md`:
1. ⏭️ Extract `medians.R` module
2. Continue with remaining Week 3 modules
3. Complete Phase 1 extraction (8 modules remaining)
4. Begin Phase 2: Optimization and cleanup

## Notes

- This is the **largest specialized module** extracted to date (9,046 lines)
- Contains the most comprehensive set of multiple comparison procedures
- Critical for researchers conducting factorial experiments
- Many functions support parallel processing (via mclapply)
- Effect size variants provide additional inferential power
- Bootstrap methods provide robust inference without normality assumptions

---

**Extraction Completed**: 2025-12-30
**Next Module**: medians.R
**Overall Progress**: 12 of 20 modules (60% complete)
