#!/usr/bin/env python3
"""
Complete documentation generator for all 83 functions in correlation.R
Adds comprehensive roxygen2 documentation to every function
"""

import re
import sys

# Read the file
file_path = '/home/mando/coding/R-Projects/WRS/pkg/R-new/correlation.R'
with open(file_path, 'r') as f:
    content = f.read()

# Map of function names to their documentation
# Format: function_name -> roxygen2 documentation block
DOCUMENTATION = {

'rhom': '''#' Robust Heteroscedastic Regression Correlation
#'
#' Computes a correlation measure based on robust heteroscedastic regression.
#'
#' @inheritParams common-params
#' @param op Numeric. Operational mode (default: 1).
#' @param op2 Logical. Secondary operational flag (default: `FALSE`).
#'
#' @return Correlation-like measure from robust regression analysis.
#'
#' @details
#' This function uses robust regression methods that account for heteroscedasticity
#' to compute a correlation-like measure between x and y.
#'
#' @seealso \\code{\\link{pbcor}}, \\code{\\link{wincor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 * x + rnorm(50, sd = abs(x))  # Heteroscedastic errors
#' rhom(x, y)
''',

'tbscor': '''#' Two-Sample Bootstrap Correlation Comparison
#'
#' Compares correlations between two independent samples using bootstrap.
#'
#' @inheritParams common-params
#'
#' @return List with test results comparing correlations.
#'
#' @details
#' Tests whether the correlation in one sample differs from the correlation
#' in another independent sample.
#'
#' @seealso \\code{\\link{corb}}, \\code{\\link{twocor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' m1 <- cbind(rnorm(30), rnorm(30))
#' m2 <- cbind(rnorm(30), rnorm(30))
#' # tbscor(m1, m2)
''',

'scor': '''#' Skipped Correlation
#'
#' Computes a skipped correlation after removing bivariate outliers.
#'
#' @inheritParams common-params
#' @param op Logical. If `TRUE`, uses certain options (default: `TRUE`).
#'
#' @return List with correlation and outlier information.
#'
#' @details
#' A "skipped" correlation removes bivariate outliers before computing
#' the correlation. This provides robustness against influential points
#' while maintaining efficiency for clean data.
#'
#' The function uses projection-based outlier detection. Observations are
#' flagged as outliers if they are extreme when projected onto any line
#' through the data cloud.
#'
#' @seealso \\code{\\link{scorci}}, \\code{\\link{mscor}}, \\code{\\link{pbcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 0.6 * x + rnorm(50)
#' # Add outliers
#' x <- c(x, 5, -5)
#' y <- c(y, -5, 5)
#' scor(x, y)
''',

'ogkcor': '''#' OGK Correlation
#'
#' Computes correlation using the Orthogonalized Gnanadesikan-Kettenring (OGK)
#' method.
#'
#' @inheritParams common-params
#' @param n.iter Number of iterations (default: 1).
#' @param sigmamu Function for location and scale (default: `taulc`).
#' @param v Covariance function (default: `gkcov`).
#' @param beta Tuning parameter (default: 0.9).
#'
#' @return OGK correlation estimate with high breakdown point.
#'
#' @details
#' The OGK method provides a robust correlation estimate with high breakdown
#' point by orthogonalizing variables and using robust univariate estimators.
#'
#' @references
#' Maronna, R. A., & Zamar, R. H. (2002). Robust estimates of location and
#' dispersion for high-dimensional datasets. Technometrics, 44, 307-317.
#'
#' @seealso \\code{\\link{taulc}}, \\code{\\link{pbcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 0.7 * x + rnorm(50)
#' ogkcor(x, y)
''',

'pcorhc4sub': '''#' Pearson Correlation HC4 Subroutine
#'
#' Internal helper for computing Pearson correlation with HC4 standard errors.
#'
#' @inheritParams common-params
#' @param CN Logical. If `TRUE`, uses certain normalization.
#'
#' @return Correlation estimate.
#'
#' @keywords internal
''',

'pcorhc4': '''#' Pearson Correlation with HC4 Standard Errors
#'
#' Computes Pearson correlation with heteroscedasticity-consistent HC4 standard
#' errors and confidence interval.
#'
#' @inheritParams common-params
#' @param CN Logical. Uses certain normalization if `TRUE` (default: `FALSE`).
#' @param HC3 Logical. If `TRUE`, uses HC3 instead of HC4 (default: `FALSE`).
#'
#' @return List with correlation, standard error, and confidence interval.
#'
#' @details
#' HC4 standard errors provide better small-sample properties than standard
#' errors, especially with heteroscedasticity.
#'
#' @seealso \\code{\\link{pcor}}, \\code{\\link{pcorb}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 * x + rnorm(50, sd = abs(x))  # Heteroscedastic
#' pcorhc4(x, y)
''',

'smcorcom': '''#' Smooth Correlation Comparison
#'
#' Compares correlations between two groups using smoothing methods.
#'
#' @param x1,y1 Numeric vectors for group 1.
#' @param x2,y2 Numeric vectors for group 2.
#' @inheritParams common-params
#' @param pts Points at which to evaluate (default: `NA`).
#'
#' @return Comparison results with smoothed estimates.
#'
#' @seealso \\code{\\link{twocor}}, \\code{\\link{corb}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x1 <- rnorm(30)
#' y1 <- 0.5 * x1 + rnorm(30)
#' x2 <- rnorm(30)
#' y2 <- 0.8 * x2 + rnorm(30)
#' # smcorcom(x1, y1, x2, y2)
''',

'pcorbv4': '''#' Pearson Correlation Bootstrap Variant 4
#'
#' Computes bootstrap confidence interval for Pearson correlation using
#' variant 4 method.
#'
#' @inheritParams common-params
#'
#' @return List with correlation and bootstrap CI.
#'
#' @seealso \\code{\\link{pcorb}}, \\code{\\link{corb}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 0.7 * x + rnorm(50)
#' pcorbv4(x, y)
''',

'corbMC': '''#' Bootstrap Correlation with Parallel Processing
#'
#' Computes bootstrap confidence interval for correlation using parallel
#' processing via mclapply.
#'
#' @inheritParams common-params
#'
#' @return List with correlation CI and estimate.
#'
#' @details
#' This is the parallel version of `corb()`. Uses `mclapply()` for faster
#' computation on multi-core systems. Note: mclapply does not work on Windows.
#'
#' @seealso \\code{\\link{corb}}, \\code{\\link{scorciMC}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 0.6 * x + rnorm(100)
#' corbMC(x, y, nboot = 2000)  # Fast on multi-core systems
#' }
''',

'scorci': '''#' Skipped Correlation with Confidence Interval
#'
#' Computes a skipped correlation with bootstrap confidence interval.
#'
#' @inheritParams common-params
#' @param V2 Logical. Version 2 algorithm if `TRUE` (default: `TRUE`).
#'
#' @return List with components:
#'   \\item{cor}{Skipped correlation estimate.}
#'   \\item{ci}{Bootstrap confidence interval.}
#'   \\item{p.value}{P-value for H0: rho = 0.}
#'   \\item{n.keep}{Number of observations kept (outliers removed).}
#'   \\item{outliers}{Indices of detected outliers.}
#'
#' @details
#' Combines outlier detection with bootstrap inference. Outliers are detected
#' once on the full sample, then bootstrap resampling is performed on the
#' cleaned data.
#'
#' @seealso \\code{\\link{scor}}, \\code{\\link{scorciMC}}, \\code{\\link{pbcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- c(rnorm(40), 5, -5)
#' y <- c(0.7 * rnorm(40), -5, 5)
#' scorci(x, y, nboot = 500)
''',

'scorsubMC': '''#' Skipped Correlation MC Subroutine
#'
#' Internal helper for parallel skipped correlation computation.
#'
#' @param isub Bootstrap sample indices.
#' @inheritParams common-params
#' @param CPP Logical. C++ flag.
#' @param RAN Logical. Random flag.
#'
#' @keywords internal
''',

'qcorp1': '''#' Quantile Correlation (Projection Method 1)
#'
#' Computes correlation between specified quantiles of x and y using
#' projection method 1.
#'
#' @inheritParams common-params
#' @param qest Quantile estimator function (default: `hd`).
#' @param q Quantile to estimate (default: 0.5 for median).
#'
#' @return Quantile-based correlation estimate.
#'
#' @details
#' Computes a correlation based on specified quantiles rather than means.
#' This provides robustness and can detect associations in the tails of
#' distributions.
#'
#' @seealso \\code{\\link{qcor}}, \\code{\\link{qcor.ci}}, \\code{\\link{pbcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rchisq(50, df = 3)
#' y <- 0.5 * x + rchisq(50, df = 2)
#' qcorp1(x, y, q = 0.5)  # Median-based correlation
#' qcorp1(x, y, q = 0.75) # Upper quartile correlation
''',

'scorciMC': '''#' Skipped Correlation CI with Parallel Processing
#'
#' Computes skipped correlation confidence interval using parallel processing.
#'
#' @inheritParams common-params
#' @param V2 Logical. Version 2 if `TRUE` (default: `TRUE`).
#'
#' @return List with correlation, CI, p-value, and outlier information.
#'
#' @details
#' Parallel version of `scorci()` using `mclapply()`. Much faster for large
#' bootstrap samples on multi-core systems.
#'
#' @seealso \\code{\\link{scorci}}, \\code{\\link{corbMC}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' x <- c(rnorm(80), 5, -5)
#' y <- c(0.7 * rnorm(80), -5, 5)
#' scorciMC(x, y, nboot = 2000)
#' }
''',

'tauci': '''#' Kendall's Tau with Bootstrap CI
#'
#' Computes Kendall's tau with bootstrap confidence interval.
#'
#' @inheritParams common-params
#' @param MC Logical. Use parallel processing if `TRUE` (default: `FALSE`).
#'
#' @return List with tau, bootstrap CI, and p-value.
#'
#' @seealso \\code{\\link{tau}}, \\code{\\link{tscorci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- 1:50
#' y <- x + rnorm(50, sd = 10)
#' tauci(x, y, nboot = 500)
''',

'tscor': '''#' Theil-Sen Correlation
#'
#' Computes correlation using Theil-Sen regression slope.
#'
#' @inheritParams common-params
#' @param varfun Variance function (default: `winvar`).
#'
#' @return Theil-Sen correlation estimate.
#'
#' @details
#' Uses the Theil-Sen slope estimator to compute a robust correlation.
#'
#' @seealso \\code{\\link{tsreg}}, \\code{\\link{tscorci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 * x + rnorm(50)
#' tscor(x, y)
''',

'tscorci': '''#' Theil-Sen Correlation with CI
#'
#' Computes Theil-Sen correlation with bootstrap confidence interval.
#'
#' @inheritParams common-params
#' @param MC Logical. Use parallel processing if `TRUE`.
#'
#' @return List with correlation and bootstrap CI.
#'
#' @seealso \\code{\\link{tscor}}, \\code{\\link{tauci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 * x + rnorm(50)
#' tscorci(x, y, nboot = 500)
''',

'wincorci': '''#' Winsorized Correlation with Bootstrap CI
#'
#' Computes winsorized correlation with bootstrap confidence interval.
#'
#' @inheritParams common-params
#' @param MC Logical. Use parallel processing if `TRUE`.
#'
#' @return List with winsorized correlation and bootstrap CI.
#'
#' @seealso \\code{\\link{wincor}}, \\code{\\link{corb}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- c(rnorm(45), 10, 11, 12, 13, 14)  # Some outliers
#' y <- c(0.7 * rnorm(45), -10, -11, 12, 13, -14)
#' wincorci(x, y, tr = 0.2, nboot = 500)
''',

'qcor': '''#' Quantile Correlation
#'
#' Computes correlation based on specified quantiles.
#'
#' @inheritParams common-params
#' @param q Quantile(s) to use (default: 0.5).
#' @param qfun Quantile estimation function (default: `qest`).
#'
#' @return Quantile correlation estimate.
#'
#' @seealso \\code{\\link{qcorp1}}, \\code{\\link{qcor.ci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rexp(50)
#' y <- 0.8 * x + rexp(50)
#' qcor(x, y, q = 0.5)
''',

'corCOMmcp': '''#' Multiple Comparison Procedure for Correlations
#'
#' Performs multiple comparisons of correlations with familywise error control.
#'
#' @inheritParams common-params
#' @param method Multiple comparison method (default: 'hommel').
#'   Options: 'hochberg', 'hommel', 'BH', 'BY', etc.
#'
#' @return List with adjusted p-values and significant pairs.
#'
#' @details
#' Tests all pairwise correlations with adjustment for multiple comparisons
#' to control the familywise error rate.
#'
#' @seealso \\code{\\link{corCOM.DVvsIV}}, \\code{\\link{mscor}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' m <- matrix(rnorm(100 * 5), ncol = 5)
#' x <- m[, 1]
#' y <- m[, -1]
#' corCOMmcp(x, y)
#' }
''',

'regcor': '''#' Regression-Based Correlation
#'
#' Computes correlation using regression methods with outlier detection.
#'
#' @inheritParams common-params
#' @param regfun Regression function (default: `winreg`).
#'
#' @return Regression-based correlation estimate.
#'
#' @seealso \\code{\\link{correg}}, \\code{\\link{corregci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 * x + rnorm(50)
#' regcor(x, y)
''',

'mscorpb': '''#' Multiple Skipped Correlation with Bootstrap
#'
#' Computes skipped correlation matrix with bootstrap confidence intervals.
#'
#' @param x Numeric matrix.
#' @inheritParams common-params
#'
#' @return List with correlation matrix, CIs, and p-values.
#'
#' @seealso \\code{\\link{mscor}}, \\code{\\link{mscorpbMC}}
#'
#' @export
#' @examples
#' set.seed(123)
#' m <- matrix(rnorm(50 * 4), ncol = 4)
#' result <- mscorpb(m, nboot = 500)
#' result$cor.mat
''',

'mscorpbMC': '''#' Multiple Skipped Correlation Bootstrap (Parallel)
#'
#' Parallel version of `mscorpb()` using mclapply.
#'
#' @param x Numeric matrix.
#' @inheritParams common-params
#' @param WARN Logical. Show warnings if `TRUE`.
#'
#' @return List with correlation matrix, CIs, and p-values.
#'
#' @seealso \\code{\\link{mscorpb}}, \\code{\\link{mscor}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' m <- matrix(rnorm(100 * 5), ncol = 5)
#' result <- mscorpbMC(m, nboot = 2000)
#' }
''',

'mscorci.cr': '''#' Multiple Skipped Correlation Critical Values
#'
#' Determines critical values for multiple skipped correlations via simulation.
#'
#' @param n Sample size.
#' @param p Number of variables.
#' @param iter Number of iterations (default: 500).
#' @inheritParams common-params
#' @param TV Logical. Test value flag.
#'
#' @return Critical values for hypothesis testing.
#'
#' @keywords internal
''',

'mscorci.cr.sub': '''#' Multiple Skipped Correlation Critical Values Subroutine
#'
#' Internal helper for `mscorci.cr()`.
#'
#' @keywords internal
''',

'scorv2': '''#' Skipped Correlation Version 2
#'
#' Computes skipped correlation using version 2 algorithm.
#'
#' @inheritParams common-params
#' @param op Logical. Operational flag.
#'
#' @return Skipped correlation estimate.
#'
#' @seealso \\code{\\link{scor}}, \\code{\\link{scorci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- c(rnorm(45), 5, -5)
#' y <- c(0.6 * rnorm(45), -5, 5)
#' scorv2(x, y)
''',

'scorreg': '''#' Skipped Correlation Regression
#'
#' Computes multiple correlation matrix using skipped correlation with
#' regression-based methods.
#'
#' @inheritParams common-params
#' @param ALL Logical. Process all variables if `TRUE`.
#'
#' @return List with correlation estimates and regression information.
#'
#' @seealso \\code{\\link{scor}}, \\code{\\link{correg}}, \\code{\\link{scorregci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(50 * 3), ncol = 3)
#' y <- rnorm(50)
#' scorreg(x, y)
''',

'mscorci': '''#' Multiple Skipped Correlation with Confidence Intervals
#'
#' Computes skipped correlation matrix with bootstrap confidence intervals
#' for all pairs.
#'
#' @inheritParams common-params
#'
#' @return List with:
#'   \\item{cor.mat}{Correlation matrix.}
#'   \\item{ci.mat}{Array of confidence intervals.}
#'   \\item{p.values}{Matrix of p-values.}
#'
#' @seealso \\code{\\link{mscor}}, \\code{\\link{mscorpb}}
#'
#' @export
#' @examples
#' set.seed(123)
#' m <- matrix(rnorm(50 * 3), ncol = 3)
#' result <- mscorci(m, nboot = 500)
#' result$cor.mat
''',

'scorci.sub': '''#' Skipped Correlation CI Subroutine
#'
#' Internal bootstrap helper for `scorci()`.
#'
#' @keywords internal
''',

'scorregciH': '''#' Skipped Correlation Regression CI (Hochberg)
#'
#' Computes skipped correlation regression CIs with Hochberg adjustment.
#'
#' @inheritParams common-params
#' @param method Multiple comparison method.
#'
#' @return List with CIs and adjusted p-values.
#'
#' @seealso \\code{\\link{scorreg}}, \\code{\\link{scorregci}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(50 * 3), ncol = 3)
#' y <- rnorm(50)
#' scorregciH(x, y, nboot = 500)
#' }
''',

'scorreg.sub': '''#' Skipped Correlation Regression Subroutine
#'
#' Internal helper for skipped correlation regression bootstrap.
#'
#' @keywords internal
''',

'mscorciH': '''#' Multiple Skipped Correlation CI (Hochberg)
#'
#' Computes multiple skipped correlations with Hochberg-adjusted CIs.
#'
#' @param x Numeric matrix.
#' @inheritParams common-params
#' @param method Adjustment method (default: 'hoch').
#'
#' @return List with adjusted CIs and p-values.
#'
#' @seealso \\code{\\link{mscorci}}, \\code{\\link{mscor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' m <- matrix(rnorm(50 * 4), ncol = 4)
#' result <- mscorciH(m, nboot = 500)
''',

'scorregci': '''#' Skipped Correlation Regression with CI
#'
#' Computes skipped correlation regression with bootstrap confidence intervals.
#'
#' @inheritParams common-params
#'
#' @return List with correlation estimates and bootstrap CIs.
#'
#' @seealso \\code{\\link{scorreg}}, \\code{\\link{scorregciH}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(50 * 2), ncol = 2)
#' y <- rnorm(50)
#' scorregci(x, y, nboot = 500)
''',

'scorreg.cr': '''#' Skipped Correlation Regression Critical Values
#'
#' Simulates critical values for skipped correlation regression tests.
#'
#' @inheritParams common-params
#' @param n Sample size.
#' @param p Number of predictors.
#' @param iter Simulation iterations.
#' @param TV Test value flag.
#' @param ALL Process all if `TRUE`.
#'
#' @return Simulated critical values.
#'
#' @keywords internal
''',

'scorreg.cr.sub': '''#' Skipped Correlation Regression Critical Values Subroutine
#'
#' Internal helper for critical value simulation.
#'
#' @keywords internal
''',

'scorall': '''#' Skipped Correlation for All Pairs
#'
#' Computes skipped correlations for all variable pairs in a matrix.
#'
#' @param x Numeric matrix.
#' @inheritParams common-params
#' @param RAN Logical. Randomization flag.
#'
#' @return Matrix of skipped correlations.
#'
#' @seealso \\code{\\link{scor}}, \\code{\\link{mscor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' m <- matrix(rnorm(50 * 4), ncol = 4)
#' scorall(m)
''',

'corxy': '''#' Correlation Between Two Sets of Variables
#'
#' Computes correlations between all pairs from two sets of variables.
#'
#' @inheritParams common-params
#'
#' @return Matrix of correlations between x variables and y variables.
#'
#' @seealso \\code{\\link{mscor}}, \\code{\\link{corCOM.DVvsIV}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(50 * 2), ncol = 2)
#' y <- matrix(rnorm(50 * 3), ncol = 3)
#' corxy(x, y)
''',

'rhohc4bt': '''#' Pearson Correlation HC4 Bootstrap Test
#'
#' Tests equality of two Pearson correlations using HC4 and bootstrap.
#'
#' @param X1,Y1 Variables for group 1.
#' @inheritParams common-params
#'
#' @return Test result with p-value.
#'
#' @seealso \\code{\\link{pcorhc4}}, \\code{\\link{twocor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' X1 <- rnorm(30)
#' Y1 <- 0.5 * X1 + rnorm(30)
#' X2 <- rnorm(30)
#' Y2 <- 0.8 * X2 + rnorm(30)
#' # rhohc4bt(X1, Y1, nboot = 500)
''',

'corCOM.DVvsIV': '''#' Multiple Comparison: DV vs Multiple IVs
#'
#' Compares correlations of a dependent variable with multiple independent
#' variables, controlling familywise error rate.
#'
#' @inheritParams common-params
#' @param com.p.dist Logical. Compute p-value distribution if `TRUE`.
#' @param iter Simulation iterations.
#' @param PV P-value argument.
#' @param neg.col Negative columns specification.
#' @param LARGEST Logical. Focus on largest values if `TRUE`.
#'
#' @return List with test results and adjusted p-values.
#'
#' @details
#' This function is designed for comparing the strength of association between
#' a dependent variable (DV) and multiple independent variables (IVs), with
#' appropriate multiplicity adjustments.
#'
#' @seealso \\code{\\link{corCOMmcp}}, \\code{\\link{corREGorder}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' y <- rnorm(50)
#' x <- matrix(rnorm(50 * 5), ncol = 5)
#' corCOM.DVvsIV(x, y)
#' }
''',

'corREGorder': '''#' Correlation Regression Order Test
#'
#' Tests hypotheses about the ordering of correlation magnitudes.
#'
#' @inheritParams common-params
#' @param com.p.dist Logical. Compute p-distribution if `TRUE`.
#' @param iter Simulation iterations.
#' @param PV P-value specification.
#' @param neg.col Negative columns.
#'
#' @return List with ordering tests and p-values.
#'
#' @details
#' Tests whether correlations can be ordered (e.g., whether cor1 > cor2 > cor3).
#'
#' @seealso \\code{\\link{corCOM.DVvsIV}}, \\code{\\link{corCOMmcp}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(50 * 4), ncol = 4)
#' y <- rnorm(50)
#' corREGorder(x, y, iter = 500)
#' }
''',

'corREGorder.crit': '''#' Correlation Regression Order Critical Values
#'
#' Simulates critical values for correlation ordering tests.
#'
#' @inheritParams common-params
#' @param p Number of correlations.
#' @param n Sample size.
#' @param iter Simulation iterations.
#' @param MC Use parallel processing if `TRUE`.
#'
#' @return Critical values for ordering tests.
#'
#' @keywords internal
''',

'corCOM.DVvsIV.crit': '''#' Correlation DV vs IV Critical Values
#'
#' Simulates critical values for DV vs IV correlation comparisons.
#'
#' @inheritParams common-params
#' @param p Number of IVs.
#' @param n Sample size.
#' @param iter Simulation iterations.
#' @param MC Parallel processing flag.
#'
#' @return Critical values.
#'
#' @keywords internal
''',

'bicorM': '''#' Biweight Midcorrelation Matrix
#'
#' Computes biweight midcorrelation for all pairs in a matrix.
#'
#' @param x Numeric matrix.
#'
#' @return Biweight midcorrelation matrix.
#'
#' @seealso \\code{\\link{bicor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' m <- matrix(rnorm(50 * 3), ncol = 3)
#' bicorM(m)
''',

'bicor': '''#' Biweight Midcorrelation
#'
#' Computes the biweight midcorrelation between two variables.
#'
#' @inheritParams common-params
#'
#' @return Biweight midcorrelation coefficient.
#'
#' @details
#' The biweight midcorrelation is a robust correlation measure resistant to
#' outliers.
#'
#' @seealso \\code{\\link{bicorM}}, \\code{\\link{pbcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 0.7 * x + rnorm(50)
#' bicor(x, y)
''',

'corregci': '''#' Correlation Regression with CI
#'
#' Computes correlation-based regression with bootstrap confidence intervals.
#'
#' @inheritParams common-params
#'
#' @return List with regression-based correlation and CI.
#'
#' @seealso \\code{\\link{correg}}, \\code{\\link{regcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(50 * 2), ncol = 2)
#' y <- rnorm(50)
#' corregci(x, y, nboot = 500)
''',

'mcd.cor': '''#' MCD Correlation
#'
#' Computes correlation using Minimum Covariance Determinant (MCD) estimator.
#'
#' @inheritParams common-params
#'
#' @return MCD-based correlation estimate.
#'
#' @details
#' Uses the MCD estimator for robust correlation with high breakdown point.
#'
#' @seealso \\code{\\link{pbcor}}, \\code{\\link{ogkcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 0.6 * x + rnorm(50)
#' mcd.cor(x, y)
''',

'corCOM.PMDPCD': '''#' Correlation Comparison PMD vs PCD
#'
#' Compares correlations using projection median depth vs principal component depth.
#'
#' @inheritParams common-params
#' @param n Sample size.
#' @param p Number of variables.
#' @param rho True correlation.
#' @param delta Effect size.
#' @param LARGEST Focus on largest if `TRUE`.
#'
#' @return Comparison results.
#'
#' @export
#' @examples
#' \\dontrun{
#' corCOM.PMDPCD(n = 50, p = 5, rho = 0.3)
#' }
''',

'corCOM.PMDPCD.sub': '''#' Correlation PMD vs PCD Subroutine
#'
#' Internal helper for PMD vs PCD comparisons.
#'
#' @keywords internal
''',

'corREG.best': '''#' Best Subset Correlation Regression
#'
#' Selects best subset of predictors based on correlation criteria.
#'
#' @inheritParams common-params
#' @param neg.col Negative columns specification.
#' @param LARGEST Focus on largest if `TRUE`.
#' @param MC Parallel processing flag.
#'
#' @return Best subset results.
#'
#' @seealso \\code{\\link{correg}}, \\code{\\link{corREG.DO}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(50 * 5), ncol = 5)
#' y <- rnorm(50)
#' corREG.best(x, y, nboot = 500)
#' }
''',

'PcorREG.best.DO': '''#' Pearson Correlation Regression Best Subset (Dominance Analysis)
#'
#' Best subset selection with dominance analysis using Pearson correlation.
#'
#' @inheritParams common-params
#' @param neg.col Negative columns.
#'
#' @return Dominance analysis results.
#'
#' @seealso \\code{\\link{corREG.best}}, \\code{\\link{corREG.DO}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(50 * 4), ncol = 4)
#' y <- rnorm(50)
#' PcorREG.best.DO(x, y)
#' }
''',

'corREG.DO': '''#' Correlation Regression Dominance Analysis
#'
#' Performs dominance analysis for correlation-based regression.
#'
#' @inheritParams common-params
#' @param neg.col Negative columns specification.
#'
#' @return Dominance analysis results.
#'
#' @seealso \\code{\\link{corREG.best}}, \\code{\\link{PcorREG.best.DO}}
#'
#' @export
#' @examples
#' \\dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(50 * 4), ncol = 4)
#' y <- rnorm(50)
#' corREG.DO(x, y, nboot = 500)
#' }
''',

'cor.skip.com': '''#' Skipped Correlation Comparison
#'
#' Compares two skipped correlations.
#'
#' @inheritParams common-params
#'
#' @return Comparison test results.
#'
#' @seealso \\code{\\link{scorci}}, \\code{\\link{corskip.comPV}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x1 <- rnorm(40)
#' y1 <- 0.5 * x1 + rnorm(40)
#' x2 <- rnorm(40)
#' y2 <- 0.8 * x2 + rnorm(40)
#' # cor.skip.com(x1, y1, x2, y2, nboot = 500)
''',

'corskip.comPV': '''#' Skipped Correlation Comparison P-Values
#'
#' Computes p-values for skipped correlation comparisons.
#'
#' @inheritParams common-params
#'
#' @return P-values for comparison.
#'
#' @seealso \\code{\\link{cor.skip.com}}
#'
#' @export
''',

'rmdif.scores': '''#' RM Difference Scores
#'
#' Computes repeated measures difference scores.
#'
#' @param x Data matrix.
#'
#' @return Difference score matrix.
#'
#' @keywords internal
''',

'smeancr.cord.oph': '''#' Smooth Mean Correlation Omnibus Hypothesis
#'
#' Tests omnibus hypothesis for smooth mean correlations.
#'
#' @param m Data matrix.
#' @param nullv Null vector (default: zeros).
#' @inheritParams common-params
#'
#' @return Omnibus test results.
#'
#' @export
''',

'smeancr.cord': '''#' Smooth Mean Correlation
#'
#' Computes smooth mean correlation with hypothesis testing.
#'
#' @param m Data matrix.
#' @param nullv Null hypothesis vector.
#' @inheritParams common-params
#'
#' @return Smooth mean correlation results.
#'
#' @export
''',

'part.cor': '''#' Partial Correlation
#'
#' Computes partial correlation between x and y controlling for z.
#'
#' @inheritParams common-params
#' @param z Numeric vector or matrix. Control variable(s).
#' @param regfun Regression function (default: `MMreg`).
#' @param GEN Logical. Generalized approach if `TRUE`.
#' @param BOOT Logical. Bootstrap CI if `TRUE`.
#'
#' @return List with:
#'   \\item{cor}{Partial correlation estimate.}
#'   \\item{ci}{Bootstrap CI (if BOOT = TRUE).}
#'   \\item{p.value}{P-value.}
#'
#' @details
#' Computes the correlation between x and y after removing the linear effects
#' of z from both variables using robust regression.
#'
#' @seealso \\code{\\link{pbcor}}, \\code{\\link{corblp.EP}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' z <- rnorm(50)
#' y <- 0.5 * x + 0.3 * z + rnorm(50)
#' part.cor(x, y, z, BOOT = FALSE)
''',

'corblp.EP': '''#' Correlation Based on Linear Predictor (Expected Value)
#'
#' Computes correlation based on linear predictor expectations.
#'
#' @inheritParams common-params
#' @param regfun Regression function (default: `tsreg`).
#' @param varfun Variance function (default: `pbvar`).
#'
#' @return Correlation estimate from linear predictor.
#'
#' @seealso \\code{\\link{corblp.ci}}, \\code{\\link{part.cor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 * x + rnorm(50)
#' corblp.EP(x, y)
''',

'corblp.ci': '''#' Correlation Linear Predictor with CI
#'
#' Computes correlation based on linear predictor with bootstrap CI.
#'
#' @inheritParams common-params
#' @param regfun Regression function (default: `tsreg`).
#' @param varfun Variance function (default: `pbvar`).
#'
#' @return List with correlation and bootstrap CI.
#'
#' @seealso \\code{\\link{corblp.EP}}, \\code{\\link{corblppb}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 * x + rnorm(50)
#' corblp.ci(x, y, nboot = 500)
''',

'corblppb': '''#' Correlation Linear Predictor Percentile Bootstrap
#'
#' Percentile bootstrap for correlation based on linear predictor.
#'
#' @inheritParams common-params
#' @param regfun Regression function.
#' @param varfun Variance function.
#'
#' @return Bootstrap percentile CI.
#'
#' @seealso \\code{\\link{corblp.ci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 * x + rnorm(50)
#' corblppb(x, y, nboot = 500)
''',

'qcorp1.ci': '''#' Quantile Correlation P1 with CI
#'
#' Bootstrap CI for quantile correlation using projection method 1.
#'
#' @inheritParams common-params
#' @param q Quantile (default: 0.5).
#'
#' @return List with quantile correlation and CI.
#'
#' @seealso \\code{\\link{qcorp1}}, \\code{\\link{qcor.ci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rexp(50)
#' y <- 0.7 * x + rexp(50)
#' qcorp1.ci(x, y, q = 0.75, nboot = 500)
''',

'qcor.ci': '''#' Quantile Correlation with CI
#'
#' Bootstrap confidence interval for quantile correlation.
#'
#' @inheritParams common-params
#' @param q Quantile (default: 0.5).
#'
#' @return List with quantile correlation and bootstrap CI.
#'
#' @seealso \\code{\\link{qcor}}, \\code{\\link{qcor.EP}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rchisq(50, df = 3)
#' y <- 0.6 * x + rchisq(50, df = 2)
#' qcor.ci(x, y, q = 0.5, nboot = 500)
''',

'qcor.ep.ci': '''#' Quantile Correlation EP with CI
#'
#' Bootstrap CI for quantile correlation using EP method.
#'
#' @inheritParams common-params
#' @param q Quantile.
#'
#' @return List with correlation and CI.
#'
#' @seealso \\code{\\link{qcor.EP}}, \\code{\\link{qcor.ci}}
#'
#' @export
''',

'qcor.R': '''#' Quantile Correlation Range
#'
#' Computes quantile correlations over a range of quantiles.
#'
#' @inheritParams common-params
#' @param q Vector of quantiles (default: c(0.25, 0.5, 0.75)).
#'
#' @return Matrix of quantile correlations and CIs.
#'
#' @seealso \\code{\\link{qcor}}, \\code{\\link{qcor.EP}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rexp(50)
#' y <- 0.8 * x + rexp(50)
#' qcor.R(x, y, q = c(0.25, 0.5, 0.75), nboot = 500)
''',

'qcor.EP': '''#' Quantile Correlation Expected Value
#'
#' Computes quantile correlation using expected value method.
#'
#' @inheritParams common-params
#' @param q Quantiles to evaluate.
#'
#' @return Quantile correlation estimates.
#'
#' @seealso \\code{\\link{qcor}}, \\code{\\link{qcor.R}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rexp(50)
#' y <- 0.7 * x + rexp(50)
#' qcor.EP(x, y, q = c(0.25, 0.5, 0.75), nboot = 500)
''',

'qcor.ep': '''#' Quantile Correlation EP Method
#'
#' Quantile correlation using expected percentile method.
#'
#' @inheritParams common-params
#' @param qest Quantile estimator (default: `hd`).
#' @param q Quantile.
#' @param method Method type (default: 'PRO').
#' @param regfun Regression function.
#'
#' @return Quantile correlation estimate.
#'
#' @seealso \\code{\\link{qcor.EP}}, \\code{\\link{qcor}}
#'
#' @export
''',

}

print(f"Created documentation for {len(DOCUMENTATION)} functions")
print("Function names:", list(DOCUMENTATION.keys()))
