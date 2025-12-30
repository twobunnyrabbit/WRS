#' Common parameter documentation for WRS package functions
#'
#' This file contains documentation for parameters that are used across
#' multiple functions in the WRS package. Other functions can inherit
#' these parameter descriptions using @inheritParams.
#'
#' @name common-params
#' @keywords internal
NULL

#' @param x A numeric vector, matrix, or data frame containing the data.
#' @param y A numeric vector, matrix, or data frame containing the data for
#'   the second group (for two-sample tests) or the second variable (for
#'   correlation/regression).
#' @param tr Proportion of observations to trim from each tail of the
#'   distribution (default: 0.2 for 20% trimming). Must be between 0 and 0.5.
#' @param alpha Significance level for confidence intervals and hypothesis tests
#'   (default: 0.05 for 95% confidence).
#' @param nboot Number of bootstrap samples to generate (default varies by
#'   function, typically 599 or 2000).
#' @param SEED Logical value or integer. If `TRUE`, uses a fixed seed for
#'   reproducibility. If `FALSE`, uses random seed. If an integer, uses that
#'   value as the seed.
#' @param q Quantile to estimate. For example, `q = 0.5` estimates the median,
#'   `q = 0.25` estimates the first quartile (default: 0.5).
#' @param plotit Logical value. If `TRUE`, creates a plot of the results
#'   (default: `TRUE` for most plotting functions, `FALSE` for computational
#'   functions).
#' @param MC Logical value. If `TRUE`, uses parallel processing via `mclapply`
#'   to speed up bootstrap computations. Note: parallel processing may not work
#'   on all systems (default: `FALSE`).
#' @param grp A vector indicating group membership. Can be numeric, factor, or
#'   character. Used for grouping observations in multi-group comparisons.
#' @param est Estimator function to use for computing location or other
#'   statistics. Common choices include `mean`, `median`, `tmean` (trimmed mean),
#'   `mom` (modified one-step M-estimator), `mest` (M-estimator).
#' @param null.value Null hypothesis value for testing (default: 0 for most
#'   tests, indicating no difference or no correlation).
#' @param side Character string specifying the alternative hypothesis. Options:
#'   \itemize{
#'     \item `"both"` or `"two.sided"`: two-sided test (default)
#'     \item `"less"`: one-sided test for negative differences
#'     \item `"greater"`: one-sided test for positive differences
#'   }
#' @param op Number of outliers to remove from each tail when using outlier
#'   detection methods. Default varies by function.
#' @param outfun Outlier detection function to use. Common choices:
#'   \itemize{
#'     \item `outpro`: Projection-based outlier detection
#'     \item `out`: MAD-based outlier detection
#'     \item `outbox`: Boxplot-based outlier detection
#'     \item `outmgv`: MGV outlier detection
#'   }
#'   Default is typically `outpro`.
#' @param xout Logical value. If `TRUE`, outliers are removed from `x` before
#'   analysis (default: `FALSE`).
#' @param yout Logical value. If `TRUE`, outliers are removed from `y` before
#'   analysis (default: `FALSE`).
#' @param eout Logical value. If `TRUE`, outliers are removed based on
#'   eliminating points with large residuals or leverage (default: `FALSE`).
#' @param cop Type of correlation or covariance estimator to use:
#'   \itemize{
#'     \item `1`: Pearson correlation
#'     \item `2`: Winsorized correlation
#'     \item `3`: Percentage bend correlation
#'     \item `4`: Spearman correlation
#'   }
#' @param pr Logical value. If `TRUE`, prints progress messages during
#'   computation (default: `TRUE` for computationally intensive functions).
#' @param MM Logical value. If `TRUE`, uses the MM-estimator for robust
#'   regression (default: `FALSE`).
#' @param bend Bending constant for percentage bend correlation (default: 0.2).
#' @param WIN Logical value. If `TRUE`, uses winsorized values instead of
#'   trimmed values (default: `FALSE`).
#' @param JR Logical value. If `TRUE`, uses the Johansen-Rousseeuw estimator
#'   (default: `TRUE` for some functions).
#' @param na.rm Logical value. If `TRUE`, removes `NA` values before computation
#'   (default: `TRUE` for most functions).
#' @param ... Additional arguments passed to other functions.
#'
#' @keywords internal
#' @name common-params
NULL
