# WRS Package Function Reference

A comprehensive guide to the robust statistics functions in the WRS package, organized by statistical application. This reference is designed to accompany Rand Wilcox's robust statistics textbook.

**Package Version:** 0.45 (final version)
**Total Functions:** ~1,755 user-facing functions
**Last Updated:** 2026-01-06

---

## Quick Navigation

1. [Location Estimation & Hypothesis Tests](#1-location-estimation--hypothesis-tests)
2. [Two-Sample Comparisons](#2-two-sample-comparisons)
3. [ANOVA & Group Comparisons](#3-anova--group-comparisons)
4. [ANCOVA (Analysis of Covariance)](#4-ancova-analysis-of-covariance)
5. [Multiple Comparisons & Contrasts](#5-multiple-comparisons--contrasts)
6. [Median-Based Methods](#6-median-based-methods)
7. [Correlation & Association](#7-correlation--association)
8. [Regression Methods](#8-regression-methods)
9. [Advanced Regression](#9-advanced-regression)
10. [Covariance & Scatter Estimation](#10-covariance--scatter-estimation)
11. [Outlier Detection & Data Depth](#11-outlier-detection--data-depth)
12. [Winsorization](#12-winsorization)
13. [Effect Size Estimation](#13-effect-size-estimation)
14. [Classification & Machine Learning](#14-classification--machine-learning)
15. [Plotting & Visualization](#15-plotting--visualization)
16. [Power Analysis](#16-power-analysis)
17. [Bootstrap & Resampling](#17-bootstrap--resampling)
18. [Core Utilities](#18-core-utilities)

---

## 1. Location Estimation & Hypothesis Tests
**71 functions** | Module: `location.R`

Robust measures of central tendency and location-based hypothesis tests.

### Core Location Estimators

**Trimmed Mean Methods**
- `tmean` - Calculate trimmed mean (removes extreme values)
- `trimse` - Standard error of trimmed mean
- `trimci` - Confidence interval for trimmed mean
- `trimvar` - Variance of trimmed mean

**M-Estimators**
- `mest` - Huber's M-estimator of location (iterative robust estimator)
- `mestci` - Confidence interval for M-estimator
- `mestse` - Standard error of M-estimator
- `mestseb` - Bootstrap standard error for M-estimator
- `onestep` - One-step M-estimator (single iteration)

**Modified One-Step (MOM)**
- `mom` - Modified one-step M-estimator (removes values beyond k*MAD)
- `momci` - Confidence interval for MOM
- `momse` - Standard error of MOM

**Harrell-Davis Quantile Estimator**
- `hd` - Harrell-Davis estimator for any quantile (default: median)
- `hdci` - Confidence interval using Harrell-Davis
- `hdpb` - Percentile bootstrap CI for Harrell-Davis
- `hdseb` - Bootstrap standard error for Harrell-Davis
- `hdmq` - Multiple quantiles (deciles) using Harrell-Davis

**Other Estimators**
- `idealf` - Ideal fourths (robust quartile estimation)
- `tauloc` - Tau estimate of location
- `tbs` - TBS (tau-based) location estimator

### Two-Sample Location Tests

**Independent Groups**
- `yuen` - Two-sample test for trimmed means (Yuen's method)
- `yuenbt` - Bootstrap version of Yuen's test
- `yuenv2` - Alternative version of Yuen's test
- `pb2gen` - Two-group percentile bootstrap for general estimators

**Dependent (Paired) Groups**
- `yuend` - Paired samples test for trimmed means
- `ydbt` - Bootstrap version for dependent groups
- `dif` - Simple paired difference test
- `difQpci` - Paired quantile differences with CI

**Multivariate Location**
- `mul.loc2g` - Multivariate two-group location test
- `Dmul.loc2g` - Dependent groups multivariate location test
- `L1median` - L1 (spatial) median estimator

### Location Comparisons & Effect Sizes

- `loc.dif.summary` - Summary of location differences across quantiles
- `dep.loc.summary` - Dependent groups location summary
- `loc2dif` - Shift function showing location changes across distribution
- `loc2difpb` - Bootstrap version of shift function
- `loc2gmulpb` - Ratio of location estimates between groups
- `loc2plot` - Graphical display of location differences

### One-Sample Tests

- `trimcibt` - Bootstrap CI for one-sample trimmed mean
- `sintv2` - One-sample CI using various estimators
- `momset` - One-sample test using MOM estimator
- `hdpbci` - One-sample percentile bootstrap CI using HD

---

## 2. Two-Sample Comparisons
**88 functions** | Module: `two-sample.R`

Methods for comparing two independent or dependent groups.

### Trimmed Mean Tests

**Independent Groups**
- `yuen` - Yuen's test for trimmed means (handles unequal variances)
- `yuenv2` - Alternative implementation
- `trimpb` - Percentile bootstrap for trimmed means
- `trimpb2` - Two-sample bootstrap with various estimators

**Dependent Groups**
- `yuend` - Dependent groups trimmed mean test
- `ydbt` - Bootstrap version for paired data
- `yuendbt` - Alternative bootstrap implementation

### Wilcoxon-Mann-Whitney Variants

- `wmw` - Robust Wilcoxon-Mann-Whitney test
- `wmwloc` - WMW with location parameter estimate
- `wmwloc2` - Extended WMW location test
- `wmwpb` - Percentile bootstrap version of WMW
- `qwmwhd` - Quantile-based WMW using Harrell-Davis
- `wmwaov` - WMW for multiple groups (ANOVA-style)
- `mwmw` - Multivariate WMW
- `mulwmw` - Multiple variable WMW comparison
- `linWMW` - Linear contrasts using WMW

### Quantile-Based Comparisons

- `qcomhd` - Compare specific quantiles using Harrell-Davis
- `qcomhdMC` - Parallel version for multiple quantiles
- `Dqcomhd` - Dependent groups quantile comparison
- `DqdifMC` - Dependent groups quantile differences (parallel)
- `qhat` - Q-hat measure (classification accuracy)
- `qhatd` - Q-hat for dependent groups

### Binomial & Discrete Data

- `twobinom` - Two-sample binomial test
- `twobici` - Binomial confidence interval for two proportions
- `disc2com` - Compare two discrete distributions
- `disc2comSK` - Skipped version for discrete data

### Effect Size Measures

- `cid` - Cliff's delta (robust effect size)
- `cidM` - Multiple Cliff's delta comparisons
- `cidv2` - Alternative Cliff's delta calculation
- `cidmul` - Multivariate Cliff's delta
- `yuen.effect` - Effect size for Yuen's test
- `yuen.effect.ci` - CI for Yuen effect size
- `akp.effect` - AKP robust effect size

### Correlation & Association Tests

- `tworhobt` - Compare two correlations (dependent groups)
- `twopcor` - Two percentage bend correlations
- `twohc4cor` - HC4 heteroscedasticity-consistent correlation comparison

### Other Specialized Tests

- `bmp` - Brunner-Munzel test (heteroscedastic ranks)
- `lplotg2` - Compare regression lines for two groups
- `ancovaWMW` - ANCOVA using WMW approach
- `RunS2` - Runs test for two samples

---

## 3. ANOVA & Group Comparisons
**52 functions** | Module: `anova.R`

Analysis of variance methods for comparing multiple groups.

### One-Way ANOVA

**Trimmed Means**
- `t1way` - One-way ANOVA for trimmed means
- `t1wayv2` - Alternative implementation with improved performance
- `t1wayF` - F-statistic version
- `t1waybt` - Bootstrap version
- `t1waybtv2` - Improved bootstrap implementation

**Alternative Methods**
- `rfanova` - Repeated measures one-way ANOVA
- `apanova` - Approximate ANOVA
- `pbanova` - Percentile bootstrap ANOVA
- `bdanova1` - Bounded differences ANOVA
- `med1way` - Median-based one-way ANOVA

**Effect Size**
- `t1way.effect` - Effect size for one-way trimmed means
- `t1way.EXES.ci` - Explanatory effect size with CI
- `KS.ANOVA.ES` - Kolmogorov-Smirnov ANOVA effect size

### Two-Way ANOVA

- `t2way` - Two-way ANOVA for trimmed means
- `t2wayv2` - Alternative implementation
- `t2waybt` - Bootstrap two-way ANOVA
- `t2way.no.p` - Two-way without p-values (diagnostics only)
- `med2way` - Median-based two-way ANOVA
- `Q2anova` - Projection-based depth two-way ANOVA

### Three-Way ANOVA

- `t3way` - Three-way ANOVA for trimmed means
- `t3wayv2` - Alternative implementation
- `Q3anova` - Projection-based three-way ANOVA

### Nested Designs

- `anova.nestA` - Nested ANOVA (one factor nested within another)
- `anova.nestAP` - Nested with random effects

### Multivariate ANOVA

- `YYmanova` - Multivariate ANOVA
- `MULtr.anova` - Multiple response trimmed ANOVA
- `mulrank` - Multivariate rank-based ANOVA

### Repeated Measures ANOVA

- `rmanova` - Repeated measures ANOVA
- `rmanovab` - Bootstrap repeated measures
- `rmanovab1` - One-way repeated measures bootstrap
- `sppba` - Split-plot design ANOVA (Type A)
- `sppbi` - Split-plot design ANOVA (Type I)
- `spp2fact` - Two-factor split-plot design

### Special Methods

- `discANOVA` - Discriminant-based ANOVA
- `BFANOVA` - Brown-Forsythe robust ANOVA
- `BFBANOVA` - Brown-Forsythe bootstrap ANOVA
- `MEDanova` - Comprehensive median-based ANOVA
- `Qanova` - Quantile-based ANOVA

---

## 4. ANCOVA (Analysis of Covariance)
**125 functions** | Module: `ancova.R`

Analysis of covariance methods for adjusting group comparisons by covariates.

### Basic Independent Groups ANCOVA

- `ancova` - Two-group ANCOVA with single covariate
- `anc2g` - Alternative two-group ANCOVA
- `ancdet` - Deterministic ANCOVA
- `ancdet2C` - Two-covariate deterministic ANCOVA

### Johnson-Neyman Technique

- `ancJN` - Johnson-Neyman regions of significance
- `ancJNtest` - Hypothesis test version

### Global Tests

- `ancGLOB` - Global test comparing all regression parameters
- `ancGLOB_pv` - P-value version of global test
- `ancGpar` - Test for parallel regression lines
- `ancGparMC` - Parallel version using multicore

### Dependent Groups ANCOVA

- `DEPanc` - Basic dependent groups ANCOVA
- `DEPancpb` - Percentile bootstrap version
- `Dancova` - Dependent ANCOVA alternative
- `Dancovamp` - Multiple predictors
- `Danctspb` - Theil-Sen based dependent ANCOVA
- `ancovad.sub` - Subsampling approach

### Advanced Methods

**Best Predictor Selection**
- `anc.best` - Find "best" predictor combination
- `anc.bestpb` - Bootstrap best predictor

**Multiple Covariates**
- `anc.grid` - Grid-based testing with multiple covariates
- `anc.grid.bin` - Binary outcome version
- `ancmg` - Multiple groups ANCOVA
- `ancmpb` - Multiple groups percentile bootstrap

**Effect Size & Visualization**
- `ancES` - Effect size for ANCOVA
- `ancES.sum` - Summary of ANCOVA effect sizes
- `ancdifplot` - Plot differences adjusted for covariate
- `reg2plot` - Compare two regression lines

### Quantile-Based ANCOVA

- `QSanc` - Quantile shift ANCOVA
- `QS.ancbse` - Bootstrap standard error for quantile shift
- `anclin.QS` - Linear ANCOVA with quantile shifts
- `anclin.QSM` - Multicore version

### Interaction Testing

- `t2way.KMS.inter` - Interaction effect size (Kraemer-Kupfer)
- `t2way.KMS.inter.ES` - Extended interaction effect size
- `anc.inter` - ANCOVA interaction test

### Regression Parameter Comparisons

- `reg2ci` - CI for difference in regression slopes
- `reg2ciMC` - Multicore version
- `ancpb` - Bootstrap ANCOVA
- `ancts` - Theil-Sen ANCOVA
- `anctspb` - Bootstrap Theil-Sen ANCOVA

---

## 5. Multiple Comparisons & Contrasts
**98 functions** | Module: `mcp.R`

Post-hoc tests and multiple comparison procedures.

### Linear Contrasts

**General Linear Contrasts**
- `lincon` - General linear contrasts among means
- `linconb` - Bootstrap version
- `linconbt` - Bootstrap-t method
- `linconpb` - Percentile bootstrap
- `linconM` - Multiple contrast testing
- `linconES` - Effect size for contrasts
- `linconQS` - Quantile shift contrasts

**Dependent Groups**
- `dlintrim` - Dependent groups linear contrasts (trimmed means)
- `dlintrimMC` - Multicore version

### Pairwise Comparisons

- `con1way` - All pairwise comparisons (one-way)
- `pairdepb` - Pairwise dependent groups bootstrap
- `Dpairdifpb` - Pairwise differences for dependent data
- `allpdiff` - All pairwise differences

### Post-Hoc Tests

**General MCP**
- `mcppb` - Multiple comparisons percentile bootstrap
- `tmcppb` - Trimmed means MCP bootstrap
- `bmcppb` - Between-groups MCP bootstrap
- `pbmcp` - Percentage bend MCP

**WMW-Based**
- `wmcpAKP` - Wilcoxon multiple comparisons (AKP method)
- `wmcpQS` - WMW quantile shift MCP
- `wmwmcp` - General WMW MCP

### Repeated Measures MCP

- `rmmcp` - Repeated measures multiple comparisons
- `rmmcppb` - Bootstrap version
- `rmmcppbd` - Pairwise differences
- `rmmcpES` - Effect size version
- `rmmcpQS` - Quantile shift version

### Complex Designs

**Between-Within Designs**
- `bwmcp` - Between-within MCP
- `bwmcppb` - Bootstrap version
- `bwmcppbES` - Effect size version

**Three-Factor Designs**
- `bbwmcp` - Between-between-within MCP
- `bwrmcp` - Between-within-repeated MCP
- `bwimcp` - Between-within interaction MCP

### Nested Designs

- `mcp.nestAP` - Nested design contrasts
- `lincon.nest` - Linear contrasts for nested designs

### Special Purpose

- `binmcp` - Binomial data multiple comparisons
- `discmcp` - Discriminant-based MCP
- `stepmcp` - Step-down multiple comparison procedure
- `Dqdif` - Quantile differences (dependent groups)

---

## 6. Median-Based Methods
**32 functions** | Module: `medians.R`

Statistical inference based on the median.

### Sample Median Estimation

- `msmed` - Sample median
- `msmedse` - Standard error of median
- `msmedci` - Confidence interval for median
- `sintv2` - Confidence interval using various estimators (includes median)

### Two-Sample Median Tests

- `med2g` - Two-group median comparison
- `medpb` - Percentile bootstrap median test
- `bpmed` - Bootstrap-t median test
- `medhd2g` - Median comparison using Harrell-Davis

### Multiple Group Median Tests

- `med1way` - One-way median ANOVA
- `med2way` - Two-way median ANOVA
- `med2mcp` - Multiple comparisons for medians
- `medanova` - General median-based ANOVA

### Effect Size

- `med.effect` - Median-based effect size
- `MED.ES` - Median explanatory effect size
- `medpb.es` - Effect size with bootstrap

### Specialized Median Methods

- `exmed` - Simulation-based median contrasts
- `medr` - Regression-based median tests
- `medind` - Median independence test
- `qsum` - Sum of quantiles (including median)

### Functional & Multivariate

- `medcurve` - Median curve for functional data
- `dmedian` - Depth-based multivariate median (deepest point)
- `dlinmed` - Linear combinations of dependent medians
- `L1median` - L1 (spatial) median

---

## 7. Correlation & Association
**83 functions** | Module: `correlation.R`

Robust correlation and association measures.

### Pearson-Type Robust Correlations

**Percentage Bend Correlation**
- `pbcor` - Percentage bend correlation (down-weights outliers)
- `pbcor.mc` - Multicore version for large datasets
- `pcorhc4` - HC4 heteroscedasticity-consistent CI
- `pcorall` - All pairwise percentage bend correlations
- `pbcorEP` - EP (explanatory power) version

**Skipped Correlation**
- `pcor` - Skipped correlation (removes outliers via projection)
- `scorall` - All pairwise skipped correlations
- `pcorb` - Bootstrap version

### Winsorized Correlation

- `wincor` - Winsorized correlation
- `wincor.pval` - P-value for winsorized correlation
- `winall` - All pairwise winsorized correlations

### Spearman Rank Correlation

- `spear` - Spearman's rho with robust inference
- `scor` - Skipped Spearman (removes outliers)
- `scorall` - All pairwise skipped Spearman
- `scorci` - CI for Spearman correlation

### Kendall's Tau

- `tau` - Kendall's tau with CI
- `tauall` - All pairwise Kendall's tau
- `tauci` - CI for tau
- `taureg` - Regression based on Kendall's tau

### Multiple Correlation

- `mscor` - Multiple correlation coefficient
- `mscorci` - CI for multiple correlation
- `pmcor` - Percentage bend multiple correlation
- `pmcorci` - CI for percentage bend multiple correlation

### Partial Correlation

- `pcor.test` - Partial correlation test
- `pbal` - Partial percentage bend correlation
- `indt` - Independence test via correlation

### Comparing Correlations

- `twopcor` - Compare two percentage bend correlations
- `twocor` - Compare two correlations (general)
- `tworhobt` - Bootstrap test for two correlations

### Ball Correlation (Distance-Based)

- `bcor` - Ball correlation (tests any type of dependence)
- `bcorci` - CI for ball correlation

---

## 8. Regression Methods
**84 functions** | Module: `regression.R`

Robust regression techniques resistant to outliers and violations of assumptions.

### Nonparametric Regression

**Theil-Sen Estimator**
- `tsreg` - Theil-Sen regression (median of all pairwise slopes)
- `tshdreg` - Theil-Sen using Harrell-Davis quantile
- `tstsreg` - Theil-Sen for time series
- `tssnmreg` - Theil-Sen smoothed non-monotonic regression

**Kernel Regression**
- `kerreg` - Kernel regression smoother
- `kercon` - Kernel regression with contrasts

### M-Regression

**Percentage Bend Regression**
- `opreg` - OP (Olive-Hawkins) percentage bend regression
- `opregpb` - Bootstrap version
- `mopreg` - Modified OP regression

**Winsorized Regression**
- `wreg` - Winsorized regression
- `wregCI` - CI for winsorized regression parameters

### Quantile Regression

- `qreg` - Quantile regression (fit regression at any quantile)
- `qregci` - CI for quantile regression parameters
- `rqfit` - Wrapper for quantreg package

### Robust Least Squares Variants

- `ltsreg` - Least trimmed squares regression
- `lmsreg` - Least median of squares regression
- `bireg` - Biweight regression
- `chreg` - Choi's regression method

### Bootstrap & CI Methods

- `regci` - CI for regression parameters (various methods)
- `reg2ci` - CI for difference in slopes (two groups)
- `reg2ciMC` - Multicore version
- `rregci` - Robust regression CI

### Linearity Testing

- `lintest` - Test for linearity
- `linchk` - Check linearity assumption
- `linmod.check` - Linear model diagnostics

### Regression Comparison

- `difreg` - Compare two regression lines
- `difregMC` - Multicore version
- `reg2plot` - Plot two regression lines
- `regtest` - Omnibus test for regression

### Dependent Groups Regression

- `DregG` - Compare dependent regression lines
- `DregGMC` - Multicore version
- `mdepreg` - Multiple dependent regressions

### Special Regression Methods

- `gyreg` - Gy regression
- `snmreg` - Smooth non-monotonic regression
- `rlplot` - Robust LOESS plot
- `rplot` - Running interval smoother plot

---

## 9. Advanced Regression
**69 functions** | Module: `regression-advanced.R`

Specialized and advanced regression techniques.

### Quantile Regression Smoothers

- `qhdsm` - Quantile smoother using Harrell-Davis
- `qhdsm2g` - Two-group quantile smoother comparison
- `qsmoothed` - General quantile smoothing
- `qsm.plot` - Plot quantile regression smoother

### Polynomial & Spline Regression

- `polyreg` - Polynomial regression
- `polyreg.plot` - Plot polynomial fits
- `adreg` - Additive regression
- `adrun` - Running additive regression

### Instrumental Variables

- `regIV` - Instrumental variable regression
- `regIVcom` - Compare IV regressions
- `tsregIV` - Theil-Sen IV regression

### Specialized Estimation

**Multiple Linear Robust**
- `mlrreg` - Multiple linear robust regression
- `mlrregES` - Effect size version

**Poisson-Inspired**
- `poireg` - Poisson-inspired robust regression
- `poisreg` - Poisson regression wrapper

**Order-Restricted**
- `regpord` - Parameter order-restricted regression
- `regpordMC` - Multicore version

### Regression Depth Methods

- `regdepth` - Regression depth
- `rdplot` - Regression depth plot

### Smoothing Techniques

- `rungen` - Running general smoother
- `runmean` - Running mean smoother
- `runhat` - Running hat matrix smoother
- `lplot` - LOESS plot
- `lplot.pred` - LOESS with prediction bands

---

## 10. Covariance & Scatter Estimation
**37 functions** | Module: `covariance.R`

Robust covariance matrix and scatter estimation.

### Covariance Matrix Estimation

**Trimmed & Winsorized**
- `covmtrim` - Trimmed covariance matrix
- `winmean` - Winsorized mean vector
- `covmve` - Minimum volume ellipsoid covariance

**Skipped Covariance**
- `skip.cov` - Covariance with outliers removed
- `skip.cor` - Correlation with outliers removed

**MVE/MCD Methods**
- `cov.mve` - Minimum volume ellipsoid
- `cov.mcd` - Minimum covariance determinant
- `covmcd` - MCD wrapper

**OGK Method**
- `cov.ogk` - Orthogonalized Gnanadesikan-Kettenring
- `covogk` - OGK wrapper

**Tau-Based**
- `cov.tbs` - Tau-based covariance estimator

### Correlation Matrices

- `corAKP` - AKP robust correlation matrix
- `corNPest` - Nonparametric correlation estimate
- `correg` - Regression-based correlation matrix
- `scorreg` - Skipped regression correlation

### Testing Covariance Structures

- `cov.test` - Test for covariance matrix equality
- `cov.mah` - Mahalanobis-based covariance test
- `homoscedasticity.test` - Test for equal covariances

---

## 11. Outlier Detection & Data Depth
**64 functions** | Module: `outliers.R`

Methods for identifying outliers and measuring data depth.

### Projection-Based Outlier Detection

**Stahel-Donoho Method**
- `outpro` - Projection-based outlier detection (Stahel-Donoho)
- `outproMC` - Multicore parallel version
- `outpro.depth` - Depth-based variant
- `outproad` - Adjusted distance variant

**Skipped Estimator**
- `out` - Basic outlier detection using MAD
- `outskip` - Skipped estimator outlier removal

### Multivariate Outlier Detection

**Distance-Based**
- `outmah` - Mahalanobis distance outliers
- `outmve` - MVE-based outliers
- `outmcd` - MCD-based outliers
- `outms` - Minimum spanning tree method

**Generalized Variance**
- `outmgv` - Generalized variance method
- `outmgvf` - Fast generalized variance
- `outmgvad` - Adjusted generalized variance

### Specialized Methods

- `outbag` - Bagplot (bivariate boxplot) outliers
- `out3d` - 3D outlier detection
- `outbox` - Univariate boxplot rule
- `outcov` - Covariance-based outliers
- `outICS` - ICS (invariant coordinate selection) outliers

### Data Depth Functions

**Projection Depth**
- `fdepth` - Projection-based depth (Stahel-Donoho)
- `fdepthv2` - Alternative implementation
- `fdepth.plot` - Plot depth contours

**Half-Space Depth**
- `depth` - Tukey half-space depth
- `depth2` - Bivariate depth
- `depthg2` - Grouped bivariate depth

**Specialized Depth**
- `indepth` - Induced depth
- `rdepth` - Regression depth
- `pdepth` - Projection depth
- `Depth.reg` - Regression depth wrapper

### Depth-Based Methods

- `dmedian` - Depth-based median (deepest point)
- `dmean` - Depth-based mean
- `dcenter` - Depth-based center
- `qchi.pval` - Depth-based chi-square test

---

## 12. Winsorization
**10 functions** | Module: `winsorize.R`

Winsorization replaces extreme values with less extreme ones.

### Core Winsorization

- `win` - Winsorize data to specified percentiles
- `winval` - Get winsorizing values (cutoff points)

### Winsorized Statistics

- `winmean` - Winsorized mean
- `winsd` - Winsorized standard deviation
- `winsd05` - Winsorized SD with different default
- `winvar` - Winsorized variance
- `winci` - CI for winsorized mean
- `winse` - Standard error of winsorized mean

### Testing with Winsorized Statistics

- `winmcp` - Multiple comparisons with winsorized means
- `win2.test` - Two-sample test using winsorized statistics

---

## 13. Effect Size Estimation
**41 functions** | Module: `effect-size.R`

Measures of effect size for various analyses.

### General Effect Size

- `ES.summary` - Summary of effect sizes across conditions
- `ES.summary.CI` - Effect size summary with CIs
- `akp.effect` - AKP robust effect size (explanatory measure)
- `akp.effectMC` - Multicore version

### Classification-Based Effect Size

- `qhat` - Q-hat measure (classification accuracy as effect size)
- `qhatd` - Q-hat for dependent groups
- `qhatsub` - Q-hat with subsampling

### Quantile-Based Effect Sizes

**Shift Functions**
- `shiftdhd` - Shift function using deciles (Harrell-Davis)
- `shifthd` - Shift function for two groups
- `shifthdMC` - Multicore version

**Quantile Comparisons**
- `qdif` - Quantile differences
- `qhatQS` - Q-hat using quantile shifts

### Test-Specific Effect Sizes

- `yuen.effect` - Effect size for Yuen's test (trimmed means)
- `yuen.effect.ci` - CI for Yuen effect size
- `t1way.effect` - Effect size for one-way ANOVA
- `med.effect` - Median-based effect size
- `reg.effect` - Regression effect size

### Comparative Effect Sizes

- `depeffect` - Dependent groups effect size
- `cidM` - Multiple Cliff's delta
- `MED.ES` - Median explanatory effect size

---

## 14. Classification & Machine Learning
**27 functions** | Module: `classification.R`

Robust classification and machine learning methods.

### Clustering

- `Kmeans` - Robust K-means clustering (uses medoid-type centers)
- `Kmeans.grp` - K-means with group membership output
- `kmeans.plot` - Plot K-means results

### K-Nearest Neighbors

- `KNN` - K-nearest neighbors classification
- `KNNv2` - Alternative implementation
- `KNNdist` - KNN using custom distance metric
- `NN.class` - Nearest neighbor classification

### Regression-Based Classification

**Ridge Regression**
- `ridge.test` - Ridge regression classification
- `ridge.est.k` - Ridge estimator with k selection
- `ridge.Liu` - Liu-type ridge estimator

**Lasso**
- `lasso.est` - Lasso regression classification
- `lasso.rep` - Lasso with replication

**Logistic Regression**
- `class.logR` - Logistic regression classification
- `log.class` - Alternative logistic classification

### Tree-Based Methods

- `class.forest` - Random forest wrapper for robust data
- `class.gbm` - Gradient boosting machine wrapper
- `class.ada` - Adaboost wrapper

### Support Vector Machines

- `class.svm` - SVM classification wrapper

### Error Assessment

- `class.error` - Classification error rate
- `class.error.cv` - Cross-validated error rate
- `class.error.boot` - Bootstrap error estimate

### Depth-Based Classification

- `Depth.class` - Classification using data depth
- `DDclass` - DD-plot based classification

---

## 15. Plotting & Visualization
**81 functions** | Module: `plotting.R`

Visualization functions for robust statistical analyses.

### Basic Regression Plots

- `rplot` - Regression plot with running interval smoother
- `lplot` - LOESS regression plot
- `lplot.pred` - LOESS with prediction bands
- `rdplot` - Robust density plot
- `akerd` - Adaptive kernel density estimate

### Regression with Confidence Intervals

- `rplotCI` - Regression plot with CI
- `rplotpbCI` - Regression plot with percentile bootstrap CI
- `regplotCI` - General regression plot with CI
- `reg2plot` - Compare two regression lines graphically

### 3D Visualization

- `ols` - 3D OLS regression surface
- `out3d` - 3D outlier detection plot
- `fdepth.plot` - Depth contour plot
- `rdepth.plot` - Robust density 3D plot

### Comparison Plots

- `loc2plot` - Compare locations between two groups
- `loc2dif` - Shift function plot (quantile differences)
- `g2plot` - General two-group comparison plot
- `difQpci` - Paired quantile difference plot

### Interaction & Effect Plots

- `interplot` - Interaction plot for two factors
- `Qinterplot` - Quantile-based interaction plot
- `plot.inter` - Interaction visualization
- `reg.plot.inter` - Regression interaction plot

### ANCOVA Plots

- `ancdifplot` - ANCOVA difference plot
- `reg2plot` - Compare regression lines (ANCOVA context)

### Repeated Measures Plots

- `spag.plot` - Spaghetti plot for repeated measures
- `bwplot` - Between-within design plot
- `bwbplot` - Bootstrap between-within plot

### Density & Distribution Plots

- `akerd` - Adaptive kernel density
- `akerd.bw` - Bandwidth selector for adaptive kernel density
- `rdplot` - Robust density estimation plot
- `kdplot` - Kernel density plot

### Boxplot Variants

- `bplot` - Bivariate boxplot (bagplot)
- `boxplot.grp` - Grouped boxplots
- `out3d` - 3D boxplot-style outlier detection

### Specialized Plots

- `runmean2g` - Running mean comparison for two groups
- `ancJNplot` - Johnson-Neyman plot (ANCOVA)
- `qsm.plot` - Quantile smoother plot
- `shiftdhd` - Shift function with deciles

---

## 16. Power Analysis
**7 functions** | Module: `power.R`

Statistical power and sample size calculations.

### One-Sample Power

- `pow1` - One-sample t-test power (Kraemer-Paik approximation)
- `powt1est` - Power estimation for one-sample tests

### Two-Sample Power

- `pow2an` - Two-sample ANOVA power
- `powt2est` - Power estimation for two-sample tests
- `powd` - Power for detecting differences

### Regression Power

- `epow` - Explanatory power (analog of R-squared)
- `power.reg` - Power analysis for regression

---

## 17. Bootstrap & Resampling
**27 functions** | Module: `bootstrap.R`

General bootstrap and resampling utilities.

### Generic Bootstrap

- `bootstrap` - General bootstrap function
- `bootse` - Bootstrap standard error for any estimator
- `bootpval` - Bootstrap p-value

### Specialized Bootstrap

- `bootdep` - Bootstrap for dependent data (time series)
- `bootcov` - Bootstrap covariance matrix
- `bootTT` - Two-sample bootstrap test
- `bootdpci` - Bootstrap CI for differences in parameters

### Parallel Bootstrap

Functions ending in `MC` use parallel processing:
- `bootMC` - Parallel bootstrap
- `outproMC` - Parallel outlier detection bootstrap

### Percentile Bootstrap

- `trimpb` - Trimmed mean percentile bootstrap
- `pb2gen` - General two-group percentile bootstrap
- `medpb` - Median percentile bootstrap

---

## 18. Core Utilities
**55 functions** | Module: `00-utils-core.R`

Data manipulation, contrasts, and helper functions.

### Data Handling

**Missing Data**
- `elimna` - Remove rows with NA values
- `NAmat` - Identify NA positions in matrix

**Data Structures**
- `listm` - Convert list to matrix
- `matl` - Convert matrix to list
- `fac2list` - Split data by factor levels
- `pool.a.list` - Concatenate list elements

### Contrast Matrices

**Factorial Contrasts**
- `con2way` - Contrasts for two-way design
- `con3way` - Contrasts for three-way design
- `con.all.pairs` - All pairwise comparison contrasts
- `conCOV` - Contrasts with covariate

**Custom Contrasts**
- `con1way` - One-way design contrasts
- `Jcon` - Johnson-Neyman contrasts

### Quantile & Location Utilities

- `qest` - Quantile estimation
- `idealf` - Ideal fourths (quartiles)
- `near` - Find nearest values
- `near3d` - Find nearest in 3D space

### Statistical Utilities

**Confidence Intervals**
- `binom.conf` - Binomial confidence intervals
- `sint` - Confidence interval for various estimators
- `trimci` - CI for trimmed mean

**Matrix Operations**
- `binmat` - Subset matrix by binary condition
- `standm` - Standardize matrix columns
- `covmtrim` - Trimmed covariance matrix

### Visualization Support

- `runmean2g` - Two-group running mean
- `plot.pval` - P-value plotting
- `ci.plot` - Confidence interval plots

---

## Appendix: Function Naming Conventions

### Common Suffixes

- `*MC` - Parallel/multicore version (uses `parallel::mclapply`)
- `*pb` - Percentile bootstrap method
- `*bt` - Bootstrap-t method
- `*b` - General bootstrap version
- `*ci` - Confidence interval function
- `*se` - Standard error function
- `*ES` - Effect size version
- `*QS` - Quantile shift version
- `*v2` - Alternative/improved implementation

### Common Prefixes

- `out*` - Outlier detection functions
- `cov*` - Covariance-related functions
- `win*` - Winsorization functions
- `lin*` - Linear contrast/model functions
- `qcom*` - Quantile comparison functions
- `med*` - Median-based methods
- `reg*` - Regression functions

### Parameter Conventions

- `tr` - Trim proportion (default 0.2 = 20% trimming)
- `q` - Quantile (default 0.5 = median)
- `nboot` - Number of bootstrap samples
- `SEED` - Random seed control (TRUE/FALSE)
- `alpha` - Significance level (default 0.05)
- `xout` - Remove outliers before analysis (TRUE/FALSE)
- `outfun` - Outlier detection function to use

---

## References

This function reference is designed to accompany:

**Wilcox, R.R.** (2017). *Introduction to Robust Estimation and Hypothesis Testing* (4th ed.). Academic Press.

**Package Citation:**
```r
citation("WRS")
```

For detailed mathematical descriptions and applications, refer to the textbook chapters corresponding to each function category.

---

**Document Version:** 1.0
**Generated:** 2026-01-06
**Package Version:** WRS 0.45 (final)
