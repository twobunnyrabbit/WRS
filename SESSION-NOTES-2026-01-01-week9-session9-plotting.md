# Session Notes: 2026-01-01 - Week 9, Session 9 - plotting.R Documentation Started

**Date**: 2026-01-01
**Phase**: Phase 3 - Documentation
**Week**: Week 9 - Advanced Modules Documentation
**Module**: plotting.R (Session 1 of ~8 estimated)
**Status**: IN PROGRESS (11/80 functions, 13.75%)

---

## Session Summary

Started documentation of plotting.R, the visualization and plotting module containing 80 functions for robust statistical plots. Documented 11 high-priority regression plotting functions focusing on the rplot family and variants.

---

## Functions Documented (11/80)

### Regression Plotting Functions (7)

1. **reg2plot** - Plot two regression lines for group comparison
   - Creates scatter plot with fitted regression lines for two independent groups
   - Default: Theil-Sen regression, group 1 solid line, group 2 dashed

2. **reg2g.p2plot** - 3D plot of two regression surfaces for group comparison
   - Creates 3D scatter plot with regression planes for two groups (2 predictors)
   - Requires scatterplot3d package
   - Color-coded: blue (group 1), red (group 2)

3. **regp2plot** - 3D plot of regression surface
   - Creates 3D scatter plot with fitted regression plane (2 predictors)
   - Single group version of reg2g.p2plot

4. **regplot** - General regression plot (1D or 2D)
   - Adapts automatically to number of predictors
   - 1 predictor: 2D scatter plot with regression line
   - 2 predictors: 3D regression surface

5. **riplot** - Regression interaction plot
   - Investigates regression interaction via additive model
   - Plots residuals to assess if additive model is adequate
   - Requires exactly 2 predictors

6. **rplot.res** - Regression residual plot with running interval smoother
   - Fits smooth using all predictors except one, plots residuals vs excluded
   - Useful for examining individual predictor contributions

### Running Interval Smoother Confidence Band Functions (4)

7. **rplotCI** - Running interval smoother with confidence band (trimmed mean)
   - Approximate simultaneous confidence band (1-alpha coverage)
   - Pre-computed critical values for npts=10 or 25 with alpha=0.05
   - Based on trimmed means

8. **rplotCIS** - Simple running interval smoother confidence band
   - Simpler alternative to rplotCI
   - Pointwise confidence intervals (no FWE adjustment)
   - Uses Student's t-distribution

9. **rplotpbCI** - Bootstrap running interval smoother confidence band
   - Uses percentile bootstrap for confidence intervals
   - Works with any location estimator (not just trimmed means)
   - More flexible but computationally intensive

10. **rplotCIM** - Running interval smoother with confidence band (median-based)
    - Uses Harrell-Davis median estimator
    - Robust alternative for heavy-tailed distributions
    - Bonferroni adjustment for simultaneous coverage

11. **rplotCIsmm** - Running interval smoother with flexible confidence band
    - Uses Studentized maximum modulus distribution
    - Works with any alpha level and any number of points
    - More flexible than rplotCI, slightly wider intervals

---

## Key Patterns Identified

### Common Parameters
- **fr**: Span parameter controlling smoothing (0.2 to 1+)
- **tr**: Trimming proportion (default: 0.2)
- **xout**: Remove outliers from predictors
- **outfun**: Outlier detection function (default: outpro)
- **plotit**: Create plot (default: TRUE)
- **pyhat**: Return fitted values (default: FALSE)
- **LP/LPCI**: Apply lowess smoothing for smoother appearance

### Function Families

1. **Basic Regression Plots**: reg2plot, regplot, regp2plot, reg2g.p2plot
2. **Diagnostic Plots**: riplot, rplot.res
3. **CI Plots (parametric)**: rplotCI, rplotCIS, rplotCIsmm
4. **CI Plots (bootstrap)**: rplotpbCI
5. **CI Plots (robust)**: rplotCIM

---

## Remaining Work

### Functions Remaining: 69/80

**Next Priorities**:

1. **Helper/Internal Functions** (~8 functions)
   - rplotCITAP.pv, rplotCITAP.sub (internal helpers for rplotCI)
   - rplotCIv2.pv, rplotCIv2.sub (internal helpers)
   - rplotCV, rplotsm, rplotN (specialized variants)

2. **Binary Regression & Comparison Plots** (~6 functions)
   - rplot.bin, rplot.binCI (binary response)
   - reg.vs.rplot, reg.vs.lplot (model comparison)
   - reg2difplot (two-group difference plot)
   - qplotreg, qregplots (quantile regression plots)

3. **Lowess/Loess Plotting Family** (~15 functions)
   - lplotv2, lplot2g, lplotCI, lplotse, lplotPV
   - lplotcom2, lplotcom2v2, lplotcomBCI, lplotcomBCIv2, etc.

4. **Group Comparison Plots** (~10 functions)
   - g2plot, g2plotdifxy, g5plot, g5.cen.plot, gplot
   - l2plot, linplot, lin2plot, sumplot2g, difQplot

5. **Box/Error Bar Plots** (~6 functions)
   - bplot, ebarplot, ebarplot.med, box_plot1, STRIPchart

6. **Functional Data Plots** (~8 functions)
   - fbplot, Flplot, FQplot, Flplot2g, func.plot, spag.plot

7. **Interaction & Special Plots** (~10 functions)
   - interplot, Qinterplot, plot.inter, reg.plot.inter, ols.plot.inter
   - logreg.plot, longreg.plot

8. **Misc/Utility Plots** (~6 functions)
   - Bagplot, rd2plot, splot, splotg5, kdplot, piplot
   - plot3D, bwiJ2plot, dlinplot, plot_robpca, testplot, plotDAP
   - prplot, MLRplot

---

## Documentation Strategy

### Approach
- **User-facing functions**: Full roxygen2 documentation with details, examples
- **Internal helpers**: Mark with @keywords internal, brief documentation
- **Similar functions**: Group documentation to maintain consistency

### Estimated Sessions
- **Session 1** (DONE): 11 high-priority regression plotting functions
- **Session 2-3**: Helper functions + binary/comparison plots (~14 functions)
- **Session 4-5**: Lowess/loess family (~15 functions)
- **Session 6**: Group comparison plots (~10 functions)
- **Session 7**: Box/error bar + functional data plots (~14 functions)
- **Session 8**: Interaction + misc/utility plots (~16 functions)

**Total**: ~8 sessions to complete all 80 functions

---

## Technical Notes

### Dependencies
- **scatterplot3d**: Required for 3D plotting functions
- **quantreg**: Required for quantile regression plots (qregplots)
- **akima**: Required for surface interpolation in some loess plots

### Parameter Inheritance
- Using `@inheritParams common-params` for standard parameters
- Using `@inheritParams rplotCI` for rplot family shared parameters

### Examples
- Most examples wrapped in `\dontrun{}` due to computational time
- Focus on demonstrating parameter usage and interpretation

---

## Statistics

### Overall Progress
- **Total functions in WRS**: ~1,500 user-facing
- **Total documented so far**: 903/~1,500 (60.2%)
- **Week 9 progress**: 112/181 target (61.9%)

### Module Progress
- **plotting.R**: 11/80 (13.75%)
- **Estimated completion**: 7 more sessions (~69 functions remaining)

---

## Next Session Plan

**Priority**: Continue with helper functions and binary/comparison plots

**Target functions** (12-15 functions):
1. rplotCITAP.pv (mark @keywords internal)
2. rplotCITAP.sub (mark @keywords internal)
3. rplotCIv2.pv (mark @keywords internal)
4. rplotCIv2.sub (mark @keywords internal)
5. rplotCV (cross-validation for running smoother)
6. rplotsm (smoothed running interval plot)
7. rplotN (large sample running interval plot)
8. rplot.bin (binary response running smoother)
9. rplot.binCI (binary response with confidence intervals)
10. reg.vs.rplot (regression vs running smoother comparison)
11. reg.vs.lplot (regression vs lowess comparison)
12. reg2difplot (two-group regression difference plot)
13. qplotreg (quantile regression lines plot)
14. qregplots (multiple quantile regression lines)
15. prplot (partial residual plot)

**Estimated time**: 1-2 hours for this batch

---

## Files Modified

1. **pkg/R-new/plotting.R**
   - Added roxygen2 documentation for 11 functions
   - Lines updated: ~31-1000 (added ~400 lines of documentation)

2. **REFACTORING-PROGRESS.md**
   - Updated Week 9 progress: 112/181 (61.9%)
   - Updated plotting.R: 11/80 (13.75%)
   - Updated overall documentation: 903/~1,500 (60.2%)
   - Updated status line

---

## Issues & Decisions

### None encountered

All functions documented successfully with consistent formatting and appropriate detail level.

---

*Session completed: 2026-01-01*
*Next session: Continue plotting.R documentation (Session 2)*
