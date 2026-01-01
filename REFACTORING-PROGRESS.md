# WRS Package Refactoring Progress

**Project**: Transform WRS from monolithic 97K-line file to modular, documented package
**Version**: v0.45 → v0.46
**Timeline**: 12 weeks
**Started**: 2025-12-30
**Last Updated**: 2026-01-01
**Phase 1 Status**: ✅ COMPLETED (All 20 modules extracted)
**Phase 2 Status**: ✅ COMPLETED (Weeks 4-5 completed, Week 6 skipped)
**Phase 3 Status**: 🔄 IN PROGRESS (Weeks 6-8 complete, Week 9 in progress: 122/181 functions, 67.4%)

---

## Quick Reference

- **Full Plan**: See `.claude/plans/curious-questing-clock.md` (approved plan with all details)
- **Original Source**: `pkg/R/Rallfun-v45.R` (97,199 lines, 1,869 function definitions, 1,828 unique)
- **Backup**: `Rallfun-v45.R.ORIGINAL` (2.6 MB - DO NOT MODIFY)
- **Refactored Modules**: `pkg/R-new/` (21 files, 98,599 lines, 1,909 functions, 0 duplicates) ✅
- **Common Parameters**: `pkg/R-new/common-params.R` (shared roxygen2 parameter docs)
- **Function Inventory**: `all-functions.txt` (sorted list of all functions)
- **Reference Tests**: `reference-outputs.rds` (baseline outputs from v0.45)
- **Deduplication Backup**: `pkg/R-new.BEFORE-DEDUP` (backup before duplicate removal)

---

## Current Status

**✅ PHASE 1 COMPLETED - 2025-12-30**
**✅ PHASE 2 COMPLETED - 2025-12-30**
**🔄 PHASE 3 IN PROGRESS - 2026-01-02** (Weeks 6-10: 1,058/~1,500 functions, 70.5%)

### Phase 1: Module Extraction (COMPLETED)
- ✅ **Week 1 COMPLETED**: Foundation modules (4 files, 219 functions)
  - 00-utils-core.R, location.R, outliers.R, bootstrap.R
- ✅ **Week 2 COMPLETED**: Main analysis modules (3 files, 256 functions)
  - two-sample.R, anova.R, correlation.R
- ✅ **Week 3 COMPLETED**: Specialized modules (13 files, 1,359 functions)
  - ✅ ancova.R extracted (149 functions, 11,014 lines, 287 KB)
  - ✅ regression.R extracted (98 functions, 5,051 lines, 140 KB)
  - ✅ regression-advanced.R extracted (74 functions, 3,433 lines, 96 KB)
  - ✅ covariance.R extracted (41 functions, 1,034 lines, 28 KB)
  - ✅ mcp.R extracted (104 functions, 9,046 lines, 256 KB)
  - ✅ medians.R extracted (42 functions, 2,474 lines, 72 KB)
  - ✅ plotting.R extracted (97 functions, 4,519 lines, 120 KB)
  - ✅ effect-size.R extracted (41 functions, 1,516 lines, 44 KB)
  - ✅ power.R extracted (10 functions, 292 lines, 8.0 KB)
  - ✅ winsorize.R extracted (9 functions, 138 lines, 4.0 KB)
  - ✅ classification.R extracted (27 functions, 1,668 lines, 52 KB)
  - ✅ special.R extracted (1,797 functions, 100,014 lines, 2.6 MB)
  - ✅ zzz-internal.R extracted (3 functions, 215 lines, 8.0 KB)

### Phase 2: Optimization (COMPLETED)
- ✅ **Week 4 COMPLETED** (2025-12-30): Library call elimination
  - Removed 325 of 558 library() calls (58% reduction)
  - Updated NAMESPACE and DESCRIPTION
  - All modules source successfully
  - All backward compatibility tests passed
- ✅ **Week 5 COMPLETED** (2025-12-30): Duplicate function resolution
  - Removed 1,171 duplicate function definitions
  - Reduced from 3,079 to 1,908 total definitions
  - special.R reduced from 1,886 to 859 functions (58% reduction, 57,881 lines removed)
  - All 20 modules source successfully
  - All 23 backward compatibility tests passed
- ⏭️ **Week 6 SKIPPED**: Pattern extraction (optional, deferred)
  - Decided to proceed directly to Phase 3 documentation (higher priority)

### Phase 3: Documentation - Core Modules (COMPLETED)
- ✅ **Setup COMPLETED** (2025-12-30): Roxygen2 infrastructure
  - Updated DESCRIPTION: Version 0.46, added Roxygen configuration
  - Created common-params.R with shared parameter documentation
- ✅ **00-utils-core.R COMPLETED** (2025-12-30): 53/53 functions documented
  - All core utility functions fully documented with roxygen2
  - File sources successfully
- ✅ **location.R COMPLETED** (2025-12-31): 71/71 functions documented
  - All robust location estimators fully documented
  - M-estimators, Harrell-Davis, group comparisons, ANOVA, multivariate location
- ✅ **outliers.R COMPLETED** (2025-12-31): 64/64 functions documented
  - All outlier detection and data depth methods fully documented
  - Projection-based, classical, robust methods, depth functions, classification/bagging

### Overall Metrics
- **Modules completed**: 20 of 20 (100%) ✅
- **Unique functions**: 1,828 of 1,828 (100%) ✅
- **Total function definitions**: 1,908 (was 3,079, removed 1,171 duplicates) ✅
- **Duplicate functions**: 0 (was 1,171, all resolved) ✅
- **Total size**: ~2.4 MB across 20 modular files (was 4.4 MB before deduplication)
- **Library calls optimized**: 325 removed, 233 remain (58% reduction) ✅
- **Roxygen2 documentation**: 1,058 of ~1,500 functions (70.5%) 🔄
  - common-params.R created ✅
  - 00-utils-core.R: 53/53 (100%) ✅
  - location.R: 71/71 (100%) ✅
  - outliers.R: 64/64 (100%) ✅
  - bootstrap.R: 27/27 (100%) ✅
  - two-sample.R: 88/88 (100%) ✅
  - anova.R: 52/52 (100%) ✅
  - correlation.R: 82/82 (100%) ✅
  - ancova.R: 125/125 (100%) ✅
  - regression.R: 84/84 (100%) ✅
  - mcp.R: 98/98 (100%) ✅
  - covariance.R: 37/37 (100%) ✅
  - regression-advanced.R: 69/69 (100%) ✅
  - medians.R: 32/32 (100%) ✅
  - plotting.R: 80/80 (100%) ✅
  - effect-size.R: 17/39 (43.6%) 🔄
- **Status**: All 20 modules source successfully ✅
- **Backward compatibility**: 100% maintained (23/23 tests pass) ✅

### Recently Completed (2026-01-01)

**Phase 3 Week 9 - Advanced Modules Documentation - COMPLETE** ✅ (181/181, 100%):
- ✅ **Completed regression-advanced.R documentation (69/69 functions, 100%)** - **Module Complete!**
  - **Session 1 COMPLETE**: Added 8 high-priority user-facing functions
    - ✅ khomreg: Test for homoscedasticity in linear regression
    - ✅ qhdsm: Quantile regression smoother using Harrell-Davis estimator
    - ✅ smean: Multivariate skipped measure of location
    - ✅ logreg: Logistic regression with outlier detection
    - ✅ gamplot: GAM regression surface plotting
    - ✅ mlrreg: Robust multivariate linear regression (Rousseeuw et al.)
    - ✅ KNNreg: K-nearest neighbors regression with robust estimation
    - ✅ regYci: Confidence intervals for predicted Y values
  - **Session 2 COMPLETE**: Added 2 additional functions
    - ✅ regYband: Confidence bands for regression predictions
    - ✅ regmediate: Mediation analysis via bootstrap
  - **Session 3 COMPLETE**: Added 8 advanced regression functions
    - ✅ mulgreg: Multivariate regression via robust covariance (Rousseeuw method)
    - ✅ gamindt: Test for association using GAM with bootstrap permutation
    - ✅ gamplotv2: Enhanced GAM plot with strength of association measures
    - ✅ logreg.P.ci: Confidence intervals for logistic regression probabilities
    - ✅ logreg.pred: Predict probabilities from logistic regression
    - ✅ mlrregCI: Bootstrap confidence intervals for multivariate regression
    - ✅ mlrregWtest: Test all slopes equal zero in multivariate regression
    - ✅ regmed2: Test mediation pathways (predictor→mediator, mediator→outcome)
  - **Session 4 COMPLETE**: Added 10 advanced regression and utility functions
    - ✅ regpca: Principal component analysis for regression
    - ✅ qhdsm2g: Two-group quantile regression smoother
    - ✅ qhdsm.pred: Predictions from quantile regression smoother
    - ✅ regse: Bootstrap standard errors for regression parameters
    - ✅ regbtci: Bootstrap-t confidence intervals for regression
    - ✅ MULMreg: Multivariate multiple regression
    - ✅ regIVcom: Compare strength of association for predictor subsets
    - ✅ regIVstr: Estimate strength of each predictor
    - ✅ regR.Forest: Random forest regression with robust location estimation
    - ✅ regGmcp: Global multiple comparisons for regression parameters
  - **Session 5 COMPLETE**: Added 10 robust multivariate and specialized regression functions
    - ✅ COVreg: Regression estimation via robust covariance matrix
    - ✅ Mreglde: Multivariate regression via least distance estimator
    - ✅ mlrreg.Stest: Test all slopes=0 in multivariate regression
    - ✅ mgvreg: Regression with MGV outlier detection
    - ✅ qhdplotsm: Plot smoothed quantile regression lines
    - ✅ longreg: Regression analysis for longitudinal data
    - ✅ regIVcommcp: Multiple comparisons of predictor strength for all pairs
    - ✅ gamplotINT: GAM regression surface with interaction term
    - ✅ regstr: Compute explanatory strength of association
    - ✅ multireg.prob: Multinomial logistic regression probabilities
  - **Session 6 COMPLETE**: Added 10 multivariate and utility regression functions
    - ✅ mgvfreg: Regression with MGV outlier detection (fast inward method)
    - ✅ smeancr: Test hypothesis about multivariate skipped mean
    - ✅ regvarp: Measure importance of predictors via generalized variance
    - ✅ smean2v2: Two-group comparison of multivariate skipped means
    - ✅ regpord: Compare strength of association for standardized predictors
    - ✅ smeancrv2: Test hypothesis about skipped mean (improved version)
    - ✅ regYci.sum: Summarize confidence intervals for regression predictions
    - ✅ regYci2Gv2: Two-group comparison of regression predictions (ANCOVA-style)
    - ✅ quantregForest: Quantile regression forests for robust prediction
    - ✅ regcon.out: Detect outliers in regression allowing heteroscedasticity
  - **Session 7 COMPLETE**: Added 10 specialized regression and diagnostic functions
    - ✅ regYciCV: Simulation-based critical values for regYci
    - ✅ regYciCV2G: Critical values for two-group regression comparison
    - ✅ reg.hyp.split: Hypothesis testing via recursive hyperplane splitting
    - ✅ regbin.hyp.split: Binary response hyperplane split testing
    - ✅ reg.resid: Compute regression residuals with any estimator
    - ✅ regIQR: Conditional interquartile range via quantile regression
    - ✅ qinvreg: Inverse quantile regression (find q for target prediction)
    - ✅ reg.con.dist: Estimate conditional distribution via quantile regression
    - ✅ reghet.blp: Regression removing bad leverage points (heteroscedastic)
    - ✅ reghet.blp.ci: Bootstrap CIs after removing bad leverage points
  - **Session 8 COMPLETE**: Added 11 remaining functions (9 internal helpers + 2 user-facing) - **Module Complete!**
    - ✅ regpord.sub: Bootstrap helper for regpord (internal)
    - ✅ Mreglde.sub: Optimization objective function for Mreglde (internal)
    - ✅ mlrreg.est: Bootstrap helper for mlrregCI (internal)
    - ✅ mlrreg.subest: Bootstrap helper for mlrregWtest (internal)
    - ✅ regunstack: Unstack data by groups for regression analysis (utility)
    - ✅ regY2G.sub: Bootstrap helper for two-group comparisons (internal)
    - ✅ regIVcom_sub: Bootstrap helper for regIVcom (internal)
    - ✅ regIVbinv2_sub: Multinomial logistic regression helper (internal)
    - ✅ qinvreg.sub: Objective function for qinvreg optimization (internal)
    - ✅ regHH: Regression after removing bad leverage points (HH method)
    - ✅ reg.break: Estimate regression breakpoint (piecewise linear regression)
- ✅ **Completed medians.R documentation (32/32 functions, 100%)** - **Module Complete!**
  - **Session 1 COMPLETE**: Added 10 high-priority median comparison functions
    - ✅ msmed: Test linear contrasts using medians with McKean-Schrader SE
    - ✅ med2g: Two-group median comparison using bootstrap
    - ✅ medpb: Multiple comparisons for medians with bootstrap (FWER control)
    - ✅ msmedci: Confidence interval for median using McKean-Schrader method
    - ✅ medcipb: Bootstrap confidence interval for median
    - ✅ med1way: Heteroscedastic one-way ANOVA for medians
    - ✅ med2way: Two-way ANOVA for medians (independent groups)
    - ✅ med.effect: Robust effect size based on median and robust variance
    - ✅ MED.ES: One-sample effect size based on median and MAD/Winsorized SD
    - ✅ dmedpb: Bootstrap comparisons for dependent groups medians
  - **Session 2 COMPLETE**: Added 2 additional user-facing functions
    - ✅ med2mcp: Multiple comparisons for two-way ANOVA using medians
    - ✅ medpb.es: Multiple comparisons for medians with shift-type effect sizes (Q)
  - **Session 3 COMPLETE**: Added 20 remaining functions (12 user-facing + 8 helpers)
    - ✅ User-facing: bpmed, bpmedse, exmed, medr, medind, medcurve, dmedian, dlinmed, wwmed, wwwmed, runstest.med, oph.astig.medianconvexpoly
    - ✅ Helpers: med2way.sub, med2way.crit, med1way.crit, msmedse, msmedsub, medpb.old, med.effect.sub, medindsub

**Phase 3 Week 8 - Advanced Analysis Module Documentation - COMPLETE** ✅ (344/344, 100%):
- ✅ **Completed ancova.R documentation (125/125 functions, 100%)** - **Module Complete!**
  - **Session 1**: Added 9 high-priority user-facing functions
    - Core ANCOVA variants: ancovaG, ancovam, ancovamp, ancovampG, ancovaV2
    - KMS effect size methods: ancova.KMS, ancova.KMSci, ancova.KMS.plot
    - Effect size function: ancova.ES
  - **Session 2**: Added 13 additional high-priority functions
    - ✅ Two-covariate ANCOVA (3): ancovap2.KMS, ancovap2.KMSci, ancovap2.KMS.plot
    - ✅ Theil-Sen methods (4): ancts, anctsmcp, anctsmp, anctspb
    - ✅ Dependent effect sizes (2): ancovad.ES, ancovad.ESci
    - ✅ Other priority functions (4): ancpar, anc.plot.es, ancsm, anctgen
  - **Session 3**: Added 13 dependent measures ANCOVA functions
    - ✅ Effect size summaries (1): Dancova.ES.sum
    - ✅ Multivariate dependent ANCOVA (2): Dancovamp, Dancovapb
    - ✅ Specified points & improved methods (3): Dancovapts, DancovaV2, DancGLOBv2
    - ✅ Theil-Sen dependent ANCOVA (3): Dancts, Danctspb, DanctspbMC
    - ✅ OLS & specialized methods (4): Dancols, DancCR, Dancdet, Danc.grid
  - **Session 4**: Added 26 helper/internal functions
    - ✅ Standard error bootstrap helpers (3): ancova.ES.SEpb, ancovad.ES.SEpb, ancovap2.KMS.SEpb
    - ✅ Internal p-value helpers (2): ancovaV2.pv, ancovaV2pv.sub
    - ✅ DancGLOB helpers (2): DancGLOB_pv, DancGLOB_sub (2 duplicate definitions)
    - ✅ Dancols helpers (3): Dancols_sub, Dancols_sub1, Dancols_sub2
    - ✅ DEP helpers (2): DEPanc, DEPancpb
    - ✅ Quantile ANCOVA helpers (4): Qancsm, QSanc, QS.ancbse, QS.ancbse.sub
    - ✅ Repeated measures helpers (7): rmanc.best, rmanc.best.crit, rmanc.best.crit.det, rmanc.best.DO, rmanc.best.ex, rmanc.bestPB, rmanc.best.PV
    - ✅ Other helpers (2): oancpb (deprecated), CLASSanc (non-robust)
- ✅ **Completed regression.R documentation (84/84 functions, 100%)** - **Module Complete!**
  - **Session 1**: Added 5 high-priority Theil-Sen and robust regression functions
    - ✅ Main Theil-Sen: tsreg, tshdreg, tsp1reg
    - ✅ LTS regression: ltsreg
    - ✅ Outlier-pruned: opreg
  - **Session 2**: Added 6 LTS variant and M-regression functions
    - ✅ M-regression methods (4): bireg, chreg, bmreg, winreg
    - ✅ LTS variants (2): MMreg, LMSreg
  - **Session 3**: Added 5 regression inference and comparison functions
    - ✅ Confidence intervals: regci, rregci, lsfitci
    - ✅ Hypothesis tests: regtest, lintest
    - ✅ Group comparisons: reg1way, reg2ci
  - **Session 4**: Added 4 high-priority functions
    - ✅ Bootstrap methods: lsfitci (LS-specific), rregci (robust), opregpb (outlier-pruned)
    - ✅ Quantile regression: Qreg
  - ✅ **Session 5**: Added 8 comparison and MC variant functions
    - ✅ Dependent groups: difreg, DregG, difregMC
    - ✅ MC variants: lintestMC, reg1wayMC, reg2ciMC
    - ✅ Quantile/other: qreg, chregF
  - ✅ **Session 6**: Added 9 diagnostic and specialized regression functions
    - ✅ Diagnostics: reglev, hratio
    - ✅ Robust methods: mbmreg, depreg, tsgreg, wreg
    - ✅ Advanced: snmreg, mdepreg, mopreg
  - ✅ **Session 7**: Added 10 user-facing regression functions
    - ✅ Global/dependent tests: DregGMC, difregOLS
    - ✅ Outlier-pruned MC: opregMC
    - ✅ Kernel regression: bkreg, kerreg
    - ✅ Correlation-based: taureg
    - ✅ Robust variants: ltsgreg, gyreg, tsregNW
    - ✅ Heteroscedastic inference: lsfitNci
  - ✅ **Session 8**: Added 10 remaining functions (82 user-facing + 2 internal helpers)
    - ✅ Conditional regression: regi (split by third variable)
    - ✅ Linearity checks: linchk (split-point test), lintests1 (internal helper)
    - ✅ Outlier detection: regout (boxplot rule on residuals)
    - ✅ S-type Theil-Sen: stsregp1 (single predictor), stsreg (general)
    - ✅ Modified Theil-Sen: tstsreg (with outlier removal), tssnmreg (TS + skipped)
    - ✅ Isotonic regression: reg1wayISOMC (one-way ANOVA, parallel)
    - ✅ Internal helpers: snmreg.sub (skipped regression optimization)
- ⏳ **mcp.R documentation IN PROGRESS: 65/98 functions (66.3%)**
  - **Session 1 COMPLETE** (2026-01-01): Documented 8 high-priority MCP functions
    - ✅ Contrast generators (1): con1way (all pairwise comparisons matrix)
    - ✅ Independent group contrasts (4): linconm (M-estimators), linconb (bootstrap-t), linconbt (alternative bootstrap-t), linconpb (percentile bootstrap with Rom's method)
    - ✅ Dependent group contrasts (2): lindep (with covariance matrix), pairdepb (pairwise bootstrap-t)
    - ✅ General MCP (1): mcppb (percentile bootstrap with Winsorization option)
  - **Session 2 COMPLETE** (2026-01-01): Documented 10 high-priority MCP functions
    - ✅ Effect size variants (3): linconEP (explanatory power), linconES (quantile shift), linconQS (quantile shift alternative)
    - ✅ Dependent contrasts bootstrap (1): lindepbt (bootstrap-t with Rom's method)
    - ✅ Factorial MCP (3): mcp2a (two-way with M-estimators), mcp2atm (two-way trimmed means), mcp3atm (three-way trimmed means)
    - ✅ Bootstrap MCP variants (3): tmcppb (trimmed means percentile bootstrap), bmcppb (alternative name), pbmcp (M-estimator percentile bootstrap)
  - **Session 3 COMPLETE** (2026-01-01): Documented 9 repeated measures MCP functions
    - ✅ Core repeated measures (2): rmmcp (trimmed means with Hochberg/Rom FWE control), wmcp (alias for rmmcp)
    - ✅ Bootstrap variants (3): rmmcppbd (percentile bootstrap with M-estimators), rmmcppbtm (percentile bootstrap for trimmed means), rmmcppbv2 (improved missing value handling)
    - ✅ Parametric v2 (1): rmmcpv2 (pairwise deletion for missing values)
    - ✅ Effect size variants (2): rmmcpES (with AKP robust Cohen's d), rmmcpQS (with quantile shift effect sizes)
    - ✅ Improved missing value support (1): rmmismcp (marginal comparisons, bootstrap)
  - **Session 4 COMPLETE** (2026-01-01): Documented 8 between-within design MCP functions
    - ✅ Two-way between-within bootstrap-t (3): bwmcp (main effects & interactions), bwwmcp (three-way B-W-W), bbwmcp (three-way B-B-W)
    - ✅ Rank-based between-within (1): bwrmcp (distribution-free method)
    - ✅ Interaction/Factor B tests (2): bwimcp (interaction contrasts only), bwbmcp (Factor B comparisons with Rom's method)
    - ✅ Percentile bootstrap variants (2): bwmcppb (two-way, any estimator), bwwmcppb (three-way B-W-W, any estimator)
  - **Session 5 COMPLETE** (2026-01-01): Documented 10 split-plot and quantile-based MCP functions
    - ✅ Split-plot designs (5): spmcpa (Factor A main effects), spmcpi (interactions), spmcpb (within-subjects pairwise), spmcpbA (by Factor A level), bwamcp (Factor A with pooling)
    - ✅ Quantile-based MCP (5): qdmcp (dependent groups), qdmcpdif (difference scores), mwwmcp (within-within quantile), twwmcp (within-within trimmed means), tkmcp (Tukey-Kramer)
  - **Session 6 COMPLETE** (2026-01-01): Documented 10 specialized & factorial bootstrap MCP functions
    - ✅ Multivariate contrasts (2): linconMpb (marginal estimators with Mahalanobis/projection distances), linconSpb (multivariate skipped estimators accounting for data structure)
    - ✅ Specialized MCP (2): mcpOV (O-type multivariate location estimator), mcp3med (three-way median-based ANOVA comparisons)
    - ✅ Between-subjects factorial (2): bbmcppb (2-way between-between), bbbmcppb (3-way between-between-between)
    - ✅ Within-subjects factorial (3): wmcppb (1-way repeated measures), wwmcppb (2-way within-within), wwwmcppb (3-way within-within-within)
    - ✅ Deprecated (1): pbtrmcp (deprecated function, redirects to bmcppb)
  - **Session 7 COMPLETE** (2026-01-01): Documented 10 between-within and within-within MCP functions
    - ✅ Three-way repeated measures (1): rm3mcp (3-way within-within-within using rmmcp)
    - ✅ Between-within interaction effect sizes (3): bwimcpES (general effect sizes), bwimcpAKP (AKP robust Cohen's d analog), bwmedimcp (median-based interactions)
    - ✅ Between-within Factor B comparisons (1): bwmedbmcp (median-based Factor B, pooled or separate)
    - ✅ Between-within comprehensive effect sizes (3): bwmcpAKP (AKP for all contrasts), bwmcpORD (ordinal/distribution-free), bwmcpQS (quantile shift)
    - ✅ Within-within two-way (2): wwmcpbt (bootstrap-t method), wwmcpQS (quantile shift effect sizes)
  - **Session 8 COMPLETE** (2026-01-01): Documented 33 remaining MCP functions - **Module Complete!**
    - ✅ Three-way within designs (2): wwwmcpQS (3-way within quantile shift), wwmcpES (deprecated, use ww.es)
    - ✅ Effect size MCP (7): esmcp (general effect sizes), bbmcpEP (between-between explanatory power), bbmcpQS (between-between quantile shift), bbdetmcp (deterministic between-between), bbdetmcpQS (deterministic QS), bmcpAKP (AKP for independent), bmcpQS (QS for independent), wmcpAKP (AKP for dependent), wmcpQS (QS for dependent)
    - ✅ Specialized MCP (5): stepmcp (step-down method), signmcp (sign test), discmcp (discrete distributions), sintmcp (median-based simultaneous), skmcp (Storer-Kim for binary)
    - ✅ Bootstrap helpers (7): bwmcppb.sub, bwmcppb.adj, bwwmcppb.sub, bbwmcppb (3-way B-B-W bootstrap), bbwmcppb.sub, wwwmcppb.sub, wwwmcppbtr (3-way within trimmed)
    - ✅ Nested and binomial (3): mcp.nestAP (nested designs), binmcp (binomial proportions), binmcp.crit (critical values)
    - ✅ Old/deprecated (4): lincon.old, lincon.pool, lincon.bin, lincon.binPV
    - ✅ P-value adjustment utilities (3): mcpPV (combined p-values), mcpKadjp (k-FWER adjustment), D1 (Romano-Shaikh helper)
    - ✅ Internal helper (1): tsub (bootstrap test statistic helper)
- ✅ **Completed covariance.R documentation (37/37 functions, 100%)** - **Module Complete!**
  - **Session 8 COMPLETE** (2026-01-01): Documented all 37 covariance estimation functions
    - ✅ OGK methods (4): gkcov, covogk, cov.ogk, skipogk (Orthogonalized Gnanadesikan-Kettenring)
    - ✅ MCD methods (3): covmcd, mcdcov, DETMCD (Minimum Covariance Determinant)
    - ✅ MVE methods (2): covmve, mvecov (Minimum Volume Ellipsoid)
    - ✅ MBA methods (5): covmba, rmba, cov.mba, covmba2, mgvcov (Median Ball Algorithm)
    - ✅ Skipped covariance (5): mscov, skipcov, covout, skip.cov (with multiple outlier detection options)
    - ✅ Biweight methods (3): bicov, bicovm, bicovM (Biweight midcovariance)
    - ✅ Winsorized methods (3): wincov, wincovN, wmean.cov
    - ✅ S-estimators (3): tbscov (Rocke's TBS), Scov (Davies), covroc
    - ✅ Mixed design covariance (2): bwwcovm (B-W-W), bbwcovm (B-B-W)
    - ✅ Median-based (2): cov2med, covmmed
    - ✅ Utilities (7): dcov (Donoho-Gasko depth), cov.roc, longcov2mat, cov.funl, covl, cov2cor
- 📊 Total Week 8 Progress: 344/344 functions documented (100%) ✅

**Phase 3 Week 7 - Main Analysis Module Documentation - COMPLETE** ✅:
- ✅ Completed bootstrap.R documentation (27/27 functions, 100%)
  - Bootstrap infrastructure, BCA methods, permutation tests
  - Bootstrap helpers for regression, ANCOVA, and correlations
- ✅ Completed two-sample.R documentation (88/88 functions, 100%)
  - Yuen's test, WMW tests, Cliff's analog methods
  - Quantile comparisons, bootstrap methods, effect sizes
  - Linear contrasts, interaction tests, ANCOVA variants
- ✅ Completed anova.R documentation (52/52 functions, 100%)
  - One-way, two-way, three-way ANOVA methods
  - Bootstrap and robust variants, trimmed mean ANOVA
  - Rank-based, median-based, and quantile-based methods
- ✅ Completed correlation.R documentation (82/82 functions, 100%)
  - Pearson, Spearman, Kendall correlations
  - Percentage bend, skipped, winsorized correlations
  - Robust correlation methods, multiple comparisons
- ✅ Total: 249 functions documented in Week 7 (plus 188 from Week 6 = 437 total, 29% of ~1,500)

**Phase 3 Week 6 - Core Module Documentation - COMPLETE** ✅:
- ✅ Completed location.R documentation (71/71 functions, 100%)
- ✅ Completed outliers.R documentation (64/64 functions, 100%)
- ✅ Completed 00-utils-core.R documentation (53/53 functions, 100%)
- ✅ Total: 188 functions documented across 3 foundation modules

**Phase 2, Week 5 - Duplicate Resolution - COMPLETE** ✅:
- ✅ Analyzed all 20 modules and identified 1,171 duplicate function definitions
- ✅ Created automated resolution strategy using R parser
- ✅ Removed 1,171 duplicates (1,027 from special.R, 144 from other modules)
- ✅ Reduced total functions from 3,079 to 1,908 (38% reduction)
- ✅ Reduced special.R from 1,886 to 859 functions (58% reduction, 57,881 lines removed)
- ✅ All 20 modules source successfully after changes
- ✅ All 23 backward compatibility tests PASSED
- ✅ 100% backward compatibility maintained
- ✅ Created backup at `pkg/R-new.BEFORE-DEDUP` for safety
- ✅ Documented resolution strategy in `pkg/duplicate-functions-resolution.md`

**Phase 2, Week 4 - Library Call Elimination - COMPLETE** ✅:
- ✅ Updated NAMESPACE with imports for parallel, MASS, akima
- ✅ Updated DESCRIPTION (moved 3 packages from Suggests to Imports)
- ✅ Removed 325 library() calls (58% reduction: 558 → 233)
  - library(parallel): 124 → 1 comment
  - library(MASS): 114 → 0
  - library(akima): 39 → 0
  - library(stats): 12 → 3 comments
- ✅ All 20 modules source successfully after changes
- ✅ All 23 backward compatibility tests PASSED
- ✅ 100% backward compatibility maintained

### Next Steps (For Next Session)

**Phase 3, Week 8: COMPLETED** ✅

**Current Status**: 903/~1,500 functions documented (60.2%)

**Completed in Week 8**:
1. ✅ **ancova.R** (125/125 functions, 100%) - 4 sessions
2. ✅ **regression.R** (84/84 functions, 100%) - 8 sessions
3. ✅ **mcp.R** (98/98 functions, 100%) - 8 sessions
4. ✅ **covariance.R** (37/37 functions, 100%) - 1 session

**Phase 3, Week 9: Advanced Modules Documentation** ✅ COMPLETE

**Week 9 Progress** (181/181 target, 100%):
- ✅ **regression-advanced.R** COMPLETE: 69/69 functions (100%)
- ✅ **medians.R** COMPLETE: 32/32 functions (100%)
- ✅ **plotting.R** COMPLETE: 80/80 functions (100%)
  - **Session 1 COMPLETE** (2026-01-01): Documented 11 high-priority regression plotting functions
    - ✅ reg2plot (two-group regression lines)
    - ✅ reg2g.p2plot (3D two-group regression surfaces)
    - ✅ regp2plot (3D regression surface)
    - ✅ regplot (general 1D/2D regression plot)
    - ✅ riplot (regression interaction plot)
    - ✅ rplot.res (regression residual plot)
    - ✅ rplotCI (running smoother with confidence band, trimmed mean)
    - ✅ rplotCIS (simple running smoother confidence band)
    - ✅ rplotpbCI (bootstrap running smoother confidence band)
    - ✅ rplotCIM (median-based running smoother confidence band)
    - ✅ rplotCIsmm (flexible running smoother confidence band)
  - **Session 2 COMPLETE** (2026-01-01): Documented 10 additional plotting functions
    - ✅ rplotCV (cross-validation prediction error for running smoother)
    - ✅ rplotsm (running smoother with strength of association)
    - ✅ rplotN (running smoother for large datasets)
    - ✅ rplot.bin (running smoother for binary outcomes)
    - ✅ rplot.binCI (confidence intervals for binary outcome probabilities)
    - ✅ reg.vs.rplot (compare regression vs running smoother estimates)
    - ✅ reg.vs.lplot (compare regression vs lowess estimates)
    - ✅ reg2difplot (3D plot of regression surface differences)
    - ✅ qplotreg (plot multiple quantile regression lines)
    - ✅ qregplots (plot quantile regression lines using quantreg package)
  - **Session 3 COMPLETE** (2026-01-01): Documented 10 lowess and group comparison plotting functions
    - ✅ lplotv2 (LOESS plot version 2 with strength of association)
    - ✅ lplot2g (LOESS plot for two groups)
    - ✅ lplotCI (LOESS confidence band with simultaneous coverage)
    - ✅ g2plot (two-group density plot)
    - ✅ g2plotdifxy (plot distribution of pairwise differences)
    - ✅ g5plot (five-group density plot)
    - ✅ g5.cen.plot (five-group centered density plot)
    - ✅ gplot (group scatterplot)
    - ✅ sumplot2g (summary plot panel for two-group comparison)
    - ✅ difQplot (quantile sums to assess symmetry about zero)
  - **Session 4 COMPLETE** (2026-01-01): Documented 10 diagnostic and LOESS comparison functions
    - ✅ prplot (partial residual plot for checking curvature)
    - ✅ MLRplot (forward response and residual plots for OLS)
    - ✅ lplotse (bootstrap SE for LOESS predictions - internal)
    - ✅ lplotPV (p-value for LOESS strength of association)
    - ✅ lplotcom2 (compare predictor importance with LOESS - Method 1)
    - ✅ lplotcom2v2 (compare predictor importance with LOESS - Method 2)
    - ✅ lplotcomBCI (bootstrap CI for comparing predictor importance)
    - ✅ lplotCIMC (bootstrap helper for lplotcom2 - internal)
    - ✅ lplotCIMCv2 (bootstrap helper for lplotcom2v2 - internal)
    - ✅ lplotbsepvv3 (critical p-value determination - internal)
  - **Session 5 COMPLETE** (2026-01-01): Documented 10 functional data and error bar plotting functions
    - ✅ lplotN (LOESS smoother for large datasets)
    - ✅ linplot (plot distribution of linear contrast)
    - ✅ lin2plot (plot two distributions based on positive/negative coefficients)
    - ✅ l2plot (LOESS smoothers for two groups)
    - ✅ bplot (split data by quantiles and compare groups)
    - ✅ ebarplot (error bar plot for group comparisons)
    - ✅ ebarplot.med (error bar plot for medians with distribution-free CIs)
    - ✅ box_plot1 (enhanced boxplot with ggplot2)
    - ✅ STRIPchart (strip chart for matrix data)
    - ✅ fbplot (functional boxplot for functional data)
  - **Session 6 COMPLETE** (2026-01-01): Documented 10 interaction and spaghetti plotting functions
    - ✅ Flplot (plot average of multiple curves)
    - ✅ FQplot (plot median and quartiles of multiple curves)
    - ✅ Flplot2g (plot average curves for two groups)
    - ✅ func.plot (functional boxplot for matrix input)
    - ✅ spag.plot (spaghetti plot for longitudinal data)
    - ✅ interplot (interaction plot for two-way design)
    - ✅ Qinterplot (quantile-based interaction plot for 2x2 design)
    - ✅ plot.inter (plot distribution of interaction contrasts)
    - ✅ reg.plot.inter (robust regression surface with interaction)
    - ✅ ols.plot.inter (OLS regression surface with interaction)
  - **Session 7 COMPLETE** (2026-01-01): Documented 10 logistic/specialized plotting functions and internal helpers
    - ✅ rplotCITAP.pv (internal: critical p-value simulation for rplotCI TAP method)
    - ✅ rplotCITAP.sub (internal: bootstrap helper for rplotCITAP.pv)
    - ✅ rplotCIv2.pv (internal: critical p-value simulation version 2)
    - ✅ rplotCIv2.sub (internal: bootstrap helper for rplotCIv2.pv)
    - ✅ lplotcomBCI9 (compare predictor importance at quartiles)
    - ✅ logreg.plot (logistic regression with MLE and robust estimates)
    - ✅ longreg.plot (individual regression lines for longitudinal data)
    - ✅ Bagplot (bivariate bagplot using data depth)
    - ✅ rd2plot (expected frequency curves for two groups)
    - ✅ splot (relative frequency plot for discrete data)
  - **Session 8 COMPLETE** (2026-01-01): Documented 9 distribution and utility plotting functions - **Module Complete!**
    - ✅ splotg5 (frequency plots for up to 5 groups)
    - ✅ kdplot (kernel density estimate plot)
    - ✅ piplot (FY plot with prediction intervals)
    - ✅ plot3D (3D scatter/surface plot)
    - ✅ bwiJ2plot (difference score distributions for J-by-2 mixed design)
    - ✅ dlinplot (linear contrast distribution for dependent groups)
    - ✅ plot_robpca (robust PCA diagnostic plot)
    - ✅ testplot (simple test utility plot)
    - ✅ plotDAP (double angle plot for ophthalmology data)
- ✅ **medians.R**: 32/32 functions (100%) - **Module Complete!**
  - **Session 1 COMPLETE** (2026-01-01): Documented 10 high-priority median comparison functions
  - **Session 2 COMPLETE** (2026-01-01): Documented 2 additional user-facing functions
  - **Session 3 COMPLETE** (2026-01-01): Documented 20 remaining functions (12 user-facing + 8 helpers)

**Phase 3 Week 10 - Utility Modules Documentation - IN PROGRESS** 🔄 (17/84, 20.2%):
- ⏳ **effect-size.R documentation IN PROGRESS: 17/39 functions (43.6%)**
  - **Session 1 COMPLETE** (2026-01-02): Documented 10 high-priority effect size functions
    - ✅ Main effect size summaries (4): ES.summary, ES.summary.CI, dep.ES.summary, dep.ES.summary.CI
    - ✅ Q-statistics (3): qhat (independent groups), qhatd (dependent groups), qhatds1 (helper)
    - ✅ AKP robust effect sizes (2): akp.effect, akp.effect.ci
    - ✅ Shift function (1): shiftdhd (decile differences for dependent groups)
  - **Session 2 COMPLETE** (2026-01-02): Documented 7 additional effect size functions
    - ✅ Pairwise comparisons (2): IND.PAIR.ES (independent), DEP.PAIR.ES (dependent)
    - ✅ Depth-based Q-statistics (2): qhatDEP, qhatdepPB
    - ✅ Multivariate/factorial (3): MUL.ES.sum, RCES, inter.ES
  - **Remaining**: 22 functions (interES.2by2, lin.ES, LCES, factorial design effect sizes, helpers)
- ⏳ **power.R**: 0/8 functions (0%)
- ⏳ **winsorize.R**: 0/10 functions (0%)
- ⏳ **classification.R**: 0/27 functions (0%)

**Phase 3, Week 9-10: Continue Documentation** 🔄

**Remaining modules to document** (~709 functions remaining to reach ~1,500 target):

1. **regression-advanced.R** (69 functions estimated)
   - Quantile regression smoothers, GAM methods, multivariate regression
   - Logistic regression, KNN, random forest methods
   - Mediation, PCA, instrumental variables

2. **medians.R** (32 functions estimated)
   - Median-based group comparisons, ANOVA, effect sizes
   - Specialized median methods

3. **plotting.R** (80 functions estimated)
   - Regression plots, ANCOVA plots, functional data plots
   - GAM plots, interaction plots, error bar plots

4. **effect-size.R** (39 functions estimated)
   - Q statistics, AKP effect sizes, ES summaries
   - Factorial/ANOVA effect sizes

5. **power.R** (8 functions estimated)
   - Power analysis for various tests

6. **winsorize.R** (10 functions estimated)
   - Winsorization utilities

7. **classification.R** (27 functions estimated)
   - KNN, K-means, ridge, LASSO, classifiers

8. **special.R** (859 functions - large module)
   - Ophthalmology methods (oph.*)
   - Binomial methods (bin.*)
   - Run tests, sign tests
   - Selection methods (selby*)
   - Many specialized domain methods

9. **zzz-internal.R** (4 functions)
   - Internal utility functions

**Recommended Next Steps**:
- **Week 9**: Document regression-advanced.R, medians.R, plotting.R (181 functions)
- **Week 10**: Document effect-size.R, power.R, winsorize.R, classification.R (84 functions)
- **Week 11-12**: Document special.R (859 functions - will require multiple sessions)
- **Week 12**: Document zzz-internal.R, review and polish all documentation

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

### 🔄 Current Phase: Phase 3 - Documentation (Week 8)

**Completed**: Phase 2 fully complete - Library optimization & duplicate resolution
**Completed**: Phase 3 Week 6 - Foundation modules documentation ✅
**Completed**: Phase 3 Week 7 - Main analysis modules documentation ✅
**Current**: Phase 3 Week 8 - Advanced analysis modules documentation (76.7% complete)
**Status**: 701/~1,500 functions documented (46.7%)
  - ✅ common-params.R created
  - ✅ 00-utils-core.R: 53/53 functions (100%)
  - ✅ location.R: 71/71 functions (100%)
  - ✅ outliers.R: 64/64 functions (100%)
  - ✅ bootstrap.R: 27/27 functions (100%)
  - ✅ two-sample.R: 88/88 functions (100%)
  - ✅ anova.R: 52/52 functions (100%)
  - ✅ correlation.R: 82/82 functions (100%)
  - ✅ ancova.R: 125/125 functions (100%) - **COMPLETE!**
  - ✅ regression.R: 84/84 functions (100%) - **COMPLETE!**
  - ⏳ mcp.R: 55/98 functions (56.1%) 🔄
  - ⏳ covariance.R: 0/37 functions (0%)
**Next**: Document mcp.R (43 remaining functions) and covariance.R (37 functions) to complete Week 8

---

## Project Goals

### Critical Constraints
- ✅ **ZERO breaking changes** - All 1,971 functions remain available
- ✅ **Maintain signatures** - No changes to function parameters or defaults
- ✅ **Backward compatible** - All outputs must match v0.45 exactly
- ✅ **Keep exports** - Maintain `exportPattern("^[^\\.]")` in NAMESPACE

### Target State (v0.46)
- [x] 20 focused module files (from 1 monolithic file) ✅
- [ ] Full roxygen2 documentation (all ~1,500 user-facing functions)
- [x] Minimal redundant library() calls (558 → 233, 58% reduction) ✅
- [x] 0 duplicate functions (1,171 duplicates resolved) ✅
- [ ] Clean namespace (471 internal helpers to mark)
- [x] Proper package imports (parallel, MASS, akima in Imports) ✅

---

## Key Statistics

### Original Codebase (Before Refactoring)

| Metric | Count | Notes |
|--------|-------|-------|
| Function Definitions | 1,869 | In Rallfun-v45.R |
| Unique Functions | 1,828 | After removing duplicates |
| Duplicate Definitions | 41 | In original file |
| Code Lines | 97,199 | Single monolithic file |
| File Size | 2.6 MB | Rallfun-v45.R |
| Redundant library() Calls | 331 | 84× parallel, 75× MASS, 21× akima |
| Documentation | 1 | Only WRS-package.Rd |
| Old .Rd Files | 15 | In `old Rd files/` directory |

### Refactored Codebase (After Phase 2, Week 5)

| Metric | Count | Notes |
|--------|-------|-------|
| Modules Created | 20 | All functional ✅ |
| Function Definitions | 1,909 | Duplicates removed ✅ |
| Unique Functions | 1,828 | 100% match with original ✅ |
| Duplicate Definitions | 0 | All 1,171 duplicates resolved ✅ |
| Code Lines | 98,599 | Across 20 modular files |
| Total Size | 2.6 MB | Reduced from 4.4 MB (41% reduction) ✅ |
| library() Calls | 233 | Reduced from 558 (58% reduction) ✅ |
| Package Imports | 3 | parallel, MASS, akima ✅ |
| Documentation | 1 | Phase 3 will add ~1,500 roxygen2 docs |
| Backward Compatibility | 100% | All 23 tests pass ✅ |

### Extracted Modules (20 of 20 Complete, Deduplicated) ✅

| Module | Status | Functions | Lines | Size | Key Functions |
|--------|--------|-----------|-------|------|---------------|
| 00-utils-core.R | ✅ | 53 | 2,809 | 72 KB | elimna, listm, matl, winvar, pbos, pbvar, yuen, yuend |
| location.R | ✅ | 71 | 3,182 | 85 KB | mest, mom, hd, tmean, onestep |
| outliers.R | ✅ | 64 | 2,707 | 67 KB | outpro, out, outbox, depth |
| bootstrap.R | ✅ | 26 | 835 | 20 KB | bootdpci, onesampb, trimcibt |
| two-sample.R | ✅ | 88 | 4,334 | 117 KB | wmw, pb2gen, cid, trimpb |
| anova.R | ✅ | 52 | 3,080 | 86 KB | t1way, t2way, t3way |
| correlation.R | ✅ | 83 | 4,669 | 152 KB | pbcor, scor, tau, wincor, mscor |
| ancova.R | ✅ | 125 | 9,333 | 246 KB | ancova, Dancova, ancES, ancGLOB |
| regression.R | ✅ | 84 | 4,376 | 119 KB | tsreg, opreg, ltsreg, regci, reg1way |
| regression-advanced.R | ✅ | 69 | 3,215 | 89 KB | qhdsm, smean, logreg, mlrreg, KNNreg |
| covariance.R | ✅ | 37 | 833 | 21 KB | covogk, wincov, skipcov, dcov, covmba |
| mcp.R | ✅ | 98 | 8,168 | 232 KB | con1way, linconb, pairdepb, rmmcp, mcppb |
| medians.R | ✅ | 32 | 2,076 | 58 KB | msmed, med2g, medpb, MEDanova, med.effect |
| plotting.R | ✅ | 80 | 3,313 | 90 KB | g2plot, gamplot, Bagplot, fbplot |
| effect-size.R | ✅ | 39 | 1,459 | 41 KB | qhat, ES.summary, akp.effect, dep.ES.summary |
| power.R | ✅ | 8 | 239 | 6 KB | pow1, powt1est, pow2an |
| winsorize.R | ✅ | 10 | 136 | 3 KB | win, winmean, winse, winci, winsd, winsdN |
| classification.R | ✅ | 27 | 1,667 | 48 KB | KNN, Kmeans, ridge.test, lasso.est, class.* |
| special.R | ✅ | 859 | 41,954 | 1.1 MB | oph.*, bin.*, run.*, selby*, specialized methods |
| zzz-internal.R | ✅ | 4 | 214 | 6 KB | wlogregv2, best.cell.crit, bestPB.DO |
| **TOTAL** | **20/20** | **1,909** | **98,599** | **2.6 MB** | **100% complete, 0 duplicates** |

Note: All 1,171 duplicates removed in Phase 2, Week 5. Unique functions: 1,828 (100% of original)

### Final Module Structure (20 Files, Deduplicated) ✅

| Module | Functions (After Dedup) | Original Estimate | Purpose |
|--------|------------------------|-------------------|---------|
| 00-utils-core.R | 53 | ~50 | Foundation utilities ✅ |
| location.R | 71 | ~80 | Robust location estimators ✅ |
| outliers.R | 64 | ~70 | Outlier detection ✅ |
| bootstrap.R | 26 | ~60 | Bootstrap infrastructure ✅ |
| two-sample.R | 88 | ~60 | Two-group comparisons ✅ |
| anova.R | 52 | ~90 | ANOVA methods ✅ |
| correlation.R | 83 | ~45 | Correlation methods ✅ |
| ancova.R | 125 | ~98 | ANCOVA methods ✅ |
| regression.R | 84 | ~85 | Regression methods ✅ |
| regression-advanced.R | 69 | ~60 | Quantile regression, etc. ✅ |
| covariance.R | 37 | ~50 | Covariance estimation ✅ |
| mcp.R | 98 | ~55 | Multiple comparisons ✅ |
| medians.R | 32 | ~40 | Median-based methods ✅ |
| plotting.R | 80 | ~50 | Visualization ✅ |
| effect-size.R | 39 | ~35 | Effect sizes ✅ |
| power.R | 8 | ~25 | Power analysis ✅ |
| winsorize.R | 10 | ~30 | Winsorization methods ✅ |
| classification.R | 27 | ~40 | Classification/ML methods ✅ |
| special.R | 859 | ~80 | Specialized & domain methods ✅ |
| zzz-internal.R | 4 | ~471 | Internal utilities ✅ |
| **TOTAL** | **1,909** | **~1,500** | **All functions extracted & deduplicated** |

**Notes**:
- All 1,171 duplicate function definitions resolved in Phase 2, Week 5 ✅
- Unique functions: 1,828 (100% of original) ✅
- special.R reduced from 1,886 to 859 functions (58% reduction)
- Most MC functions integrated into their domain modules rather than separate parallel.R

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

### ✅ Phase 1: Module Extraction (Weeks 1-3) - COMPLETED 2025-12-30
**Status**: COMPLETED ✅
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

#### ✅ Week 3: Specialized Modules (COMPLETED 2025-12-30)
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
- [x] Extract plotting.R
  - **COMPLETED**: 97 functions extracted (4,519 lines, 119 KB)
  - Includes: rplot, lplot, g2plot, gamplot, Bagplot, fbplot, ebarplot
  - Regression plots, lowess/loess, ANCOVA plots, group comparisons
  - GAM plots, functional data plots, interaction plots
  - Successfully sourced and validated
  - **Note**: Duplicate linplot at line 68059 excluded (kept line 28226 version)
- [x] Extract effect-size.R
  - **COMPLETED**: 41 functions extracted (1,516 lines, 43 KB)
  - Includes: qhat, qhatd, ES.summary, ES.summary.CI, akp.effect, dep.ES.summary
  - Q statistics, AKP robust effect sizes, general ES summaries
  - Factorial/ANOVA ES, interaction ES, linear combination ES
  - Successfully sourced and validated
- [x] Extract power.R
  - **COMPLETED**: 10 functions extracted (292 lines, 7.3 KB)
  - Includes: pow1, powt1est, powt1an, pow2an, powest, anova_power, ancmg1.power
  - One-sample power, two-sample power, ANOVA power, regression power
  - Successfully sourced and validated
- [x] Extract winsorize.R
  - **COMPLETED**: 10 functions extracted (138 lines, 3.2 KB)
  - Includes: win, winmean, winse, winci, winsd, winsd05, winsdN, winvarN, winsorized, WINCOR
  - Core winsorization, standard errors/CIs, standard deviations, normalized variance
  - Successfully sourced and validated
  - **Note**: Core utilities (winvar, winval, winall) remain in 00-utils-core.R
- [x] Extract classification.R
  - **COMPLETED**: 27 functions extracted (1,668 lines, 52 KB)
  - Includes: KNN, Kmeans, ridge regression, LASSO, various classifiers
  - Machine learning and classification methods
  - Successfully sourced and validated
- [x] Extract special.R
  - **COMPLETED**: 1,797 functions extracted (100,014 lines, 2.6 MB)
  - Includes: ophthalmology (oph.*), binomial (bin.*), run tests, sign tests
  - Selection methods (selby*), smoothing, ANOVA extensions, specialized methods
  - Successfully sourced and validated
- [x] Extract zzz-internal.R
  - **COMPLETED**: 3 functions extracted (215 lines, 8.0 KB)
  - Includes: wlogregv2, best.cell.crit, bestPB.DO
  - Internal utility functions
  - Successfully sourced and validated
- [x] **End-of-phase validation**: All 20 modules source successfully ✅
  - All modules tested and pass without errors
  - 1,828 unique functions extracted (100% of original)
  - Ready for Phase 2

### ✅ Phase 2: Optimization (Weeks 4-5)
**Status**: Weeks 4-5 COMPLETED - 2025-12-30 ✅
**Goal**: Remove redundant code, fix duplicates, optimize patterns

#### ✅ Week 4: Eliminate library() Calls (COMPLETED - 2025-12-30)
- [x] Update NAMESPACE (add imports)
  - Added importFrom() for parallel (mclapply, detectCores)
  - Added importFrom() for MASS (ginv, cov.rob, cov.mcd, cov.mve, lqs, rlm)
  - Added importFrom() for akima (interp, aspline)
- [x] Update DESCRIPTION (move to Imports)
  - Moved parallel, MASS, akima from Suggests to Imports
  - Added all discovered optional packages to Suggests
- [x] Remove library(parallel) calls (124 instances → 1 comment only)
- [x] Remove library(MASS) calls (114 instances → 0 remaining)
- [x] Remove library(akima) calls (39 instances → 0 remaining)
- [x] Remove library(stats) calls (12 instances → 3 comments only)
- [x] Handle optional packages (233 calls remain for packages in Suggests)
- [x] Validate after each module (all 20 modules source successfully)
- [x] **End-of-week validation**: All 23 backward compatibility tests PASSED ✅

**Week 4 Results**:
- **Library calls removed**: 325 of 558 (58% reduction)
- **Remaining calls**: 233 (all for optional packages in Suggests)
- **All modules**: Source successfully without errors
- **Backward compatibility**: 100% maintained (23/23 tests passed)

#### ✅ Week 5: Fix Duplicates (COMPLETED - 2025-12-30)
- [x] Analyze all modules for duplicates (found 1,171 total duplicates)
- [x] Create automated resolution strategy using R parser
- [x] Document resolution strategy in `pkg/duplicate-functions-resolution.md`
- [x] Remove duplicates from special.R (1,027 removed, 1,886 → 859 functions)
- [x] Remove duplicates from other 15 modules (144 removed)
- [x] Create backup at `pkg/R-new.BEFORE-DEDUP` before changes
- [x] Validate all 20 modules source successfully
- [x] Run backward compatibility tests (23/23 passed)
- [x] **End-of-week validation**: All tests PASSED, 100% backward compatibility ✅

**Week 5 Results**:
- **Duplicates removed**: 1,171 (1,027 from special.R, 144 from other modules)
- **Total functions**: Reduced from 3,079 to 1,908 (38% reduction)
- **special.R size**: Reduced by 58% (57,881 lines removed)
- **All modules**: Source successfully without errors
- **Backward compatibility**: 100% maintained (23/23 tests passed)

#### Week 6: Extract Common Patterns (SKIPPED)
- [x] **DECISION**: Skip optional pattern extraction, proceed to Phase 3 documentation
  - Pattern extraction can be done later as optimization
  - Documentation is higher priority for package usability

### 📋 Phase 3: Documentation - Core (Weeks 6-8)
**Status**: 🔄 IN PROGRESS (Started 2025-12-30)
**Goal**: Document all user-facing functions in core modules

#### Setup (COMPLETED - 2025-12-30)
- [x] Add Roxygen to DESCRIPTION (Version updated to 0.46, RoxygenNote: 7.3.2)
- [x] Create common-params.R (comprehensive shared parameter documentation)

#### Week 6: Document Foundation Modules (COMPLETED - 2025-12-31) ✅
- [x] **00-utils-core.R COMPLETED**: 53/53 functions documented (100%)
  - All core utilities: elimna, listm, matl, near, hd, winvar, etc.
  - File sources successfully
  - Complete roxygen2 documentation with examples
- [x] **location.R COMPLETED**: 71/71 functions documented (100%)
  - M-estimators, Harrell-Davis estimators, group comparisons
  - ANOVA functions, location comparisons, multivariate location
  - Mean estimators and utility functions
- [x] **outliers.R COMPLETED**: 64/64 functions documented (100%)
  - Projection-based, classical, and robust outlier detection
  - Data depth methods and classification/bagging functions
  - All functions fully documented with roxygen2

**Week 6 Results**: 188/188 foundation functions documented (100%) ✅
- [x] Roxygen2 syntax validated - generates valid .Rd files ✅
- [x] All 23 backward compatibility tests PASSED ✅
- [x] All 3 foundation modules source successfully ✅

#### Week 7: Document Main Analysis Modules (COMPLETED) ✅
- [x] **bootstrap.R COMPLETED**: 27/27 functions documented (100%)
  - Bootstrap infrastructure, BCA methods, permutation tests
  - Bootstrap helpers for regression, ANCOVA, and correlations
- [x] **two-sample.R COMPLETED**: 88/88 functions documented (100%)
  - Yuen's test, WMW tests, Cliff's analog methods
  - Quantile comparisons, bootstrap methods, effect sizes
  - Linear contrasts, interaction tests, ANCOVA variants
- [x] **anova.R COMPLETED**: 52/52 functions documented (100%)
  - One-way, two-way, three-way ANOVA methods
  - Bootstrap and robust variants, trimmed mean ANOVA
  - Rank-based, median-based, and quantile-based methods
- [x] **correlation.R COMPLETED**: 82/82 functions documented (100%)
  - Pearson, Spearman, Kendall correlations
  - Percentage bend, skipped, winsorized correlations
  - Robust correlation methods, multiple comparisons

**Week 7 Results**: 249/249 main analysis functions documented (100%) ✅
- [x] All 4 main analysis modules fully documented ✅
- [x] All modules source successfully ✅
- [x] Professional-grade CRAN-level documentation ✅

#### Week 8: Document Advanced Analysis Modules (PENDING)
- [ ] ancova.R: 125 functions
- [ ] regression.R: 84 functions
- [ ] mcp.R: 98 functions
- [ ] covariance.R: 37 functions

**Target**: ~344 additional functions documented

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

### 2025-12-30: Phase 2, Week 4 - Library Call Optimization
**Q**: Which packages should go to Imports vs Suggests?
**A**: Moved parallel, MASS, akima to Imports (heavily used: 124, 114, 39 calls respectively). Kept remaining packages in Suggests for optional functionality.

**Q**: Should we remove all library() calls?
**A**: Removed calls for imported packages (325 removed). Kept 233 calls for optional packages in Suggests - this is acceptable R practice.

**Results**: 58% reduction in library() calls, all tests pass, zero breaking changes.

### 2025-12-30: Phase 2, Week 5 - Duplicate Function Resolution
**Q**: How many duplicates and where?
**A**: Found 1,171 duplicate function definitions across 20 modules. special.R had 1,027 duplicates (58% of its functions were duplicates from other modules).

**Q**: What resolution strategy?
**A**: Automated resolution using R parser:
- Domain-specific modules take priority over special.R
- Core utilities stay in 00-utils-core.R
- Specialized functions (oph.*, bin.*, run.*, selby*) remain in special.R
- Used parser-based extraction to avoid cutting functions mid-definition

**Q**: How to ensure no breaking changes?
**A**: Created backup (`pkg/R-new.BEFORE-DEDUP`), validated each file sources successfully, ran full backward compatibility test suite (23/23 tests pass).

**Results**: 1,171 duplicates removed (38% reduction in total functions), special.R reduced by 58%, all tests pass, zero breaking changes.

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

*Last updated: 2026-01-02 - Phase 3 Week 10 IN PROGRESS 🔄 (Phase 2 complete: 325 library() calls removed, 1,171 duplicates resolved. Phase 3: 1,058/~1,500 functions documented [70.5%]. **Weeks 6-9 complete**: 14 modules fully documented. **Week 10 IN PROGRESS**: Utility modules (17/84 functions, 20.2%). **effect-size.R IN PROGRESS** [17/39, 43.6%]: Sessions 1-2 complete - documented main effect size summaries (ES.summary, ES.summary.CI, dep.ES.summary, dep.ES.summary.CI), Q-statistics (qhat, qhatd, qhatDEP, qhatdepPB), AKP robust effect sizes (akp.effect, akp.effect.ci), pairwise comparisons (IND.PAIR.ES, DEP.PAIR.ES), multivariate/factorial effect sizes (MUL.ES.sum, RCES, inter.ES), and shift function (shiftdhd). 22 functions remaining. All 23 backward compatibility tests passing, 100% compatible.)*
