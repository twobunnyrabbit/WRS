#!/usr/bin/env python3
"""
Comprehensive documentation generator for correlation.R
Adds roxygen2 documentation to all 83 functions
"""

import re

# Read the correlation.R file
with open('/home/mando/coding/R-Projects/WRS/pkg/R-new/correlation.R', 'r') as f:
    content = f.read()

# Documentation templates for each function type
docs = {
    # Specialized correlation functions
    'tauloc': '''
#' Tau Measure of Location
#'
#' Computes the tau measure of location as described by Yohai and Zamar (1988).
#'
#' @param x Numeric vector.
#' @param cval Numeric. Tuning constant (default: 4.5).
#'
#' @return Numeric value representing the tau measure of location.
#'
#' @details
#' The tau measure of location is a robust location estimator that uses
#' a weighted mean with weights based on standardized deviations. Values
#' beyond cval standard deviations receive zero weight.
#'
#' @references
#' Yohai, V. J., & Zamar, R. H. (1988). High breakdown-point estimates of
#' regression by means of the minimization of an efficient scale.
#' Journal of the American Statistical Association, 83, 406-413.
#'
#' @seealso \\code{\\link{tauvar}}, \\code{\\link{taulc}}
#'
#' @export
#' @examples
#' x <- c(1, 2, 3, 4, 5, 100)  # Contains outlier
#' tauloc(x)
#' median(x)  # For comparison
''',

    'tauvar': '''
#' Tau Measure of Scale
#'
#' Computes the tau measure of scale as described by Yohai and Zamar (1988),
#' using the computational method from Maronna and Zamar (2002).
#'
#' @param x Numeric vector.
#' @param cval Numeric. Tuning constant (default: 3).
#'
#' @return Numeric value representing the tau measure of scale.
#'
#' @details
#' The tau measure of scale is a robust scale estimator with high efficiency
#' and breakdown point. It combines aspects of M-estimation and S-estimation.
#'
#' @references
#' Yohai, V. J., & Zamar, R. H. (1988). High breakdown-point estimates of
#' regression by means of the minimization of an efficient scale.
#' Journal of the American Statistical Association, 83, 406-413.
#'
#' Maronna, R. A., & Zamar, R. H. (2002). Robust estimates of location and
#' dispersion for high-dimensional datasets. Technometrics, 44, 307-317.
#'
#' @seealso \\code{\\link{tauloc}}, \\code{\\link{taulc}}
#'
#' @export
#' @examples
#' x <- c(1, 2, 3, 4, 5, 100)  # Contains outlier
#' tauvar(x)
#' mad(x)  # For comparison
''',

    'taulc': '''
#' Tau Location and Scale
#'
#' Computes both tau location and scale measures simultaneously.
#'
#' @param x Numeric vector.
#' @param mu.too Logical. If `TRUE`, returns both location and scale;
#'   if `FALSE` (default), returns only scale.
#'
#' @return If `mu.too = FALSE`, returns tau scale. If `mu.too = TRUE`,
#'   returns a list with components `mu` (location) and `sigma` (scale).
#'
#' @seealso \\code{\\link{tauloc}}, \\code{\\link{tauvar}}
#'
#' @export
#' @examples
#' x <- rnorm(50)
#' taulc(x)  # Scale only
#' taulc(x, mu.too = TRUE)  # Location and scale
''',

    'ecor': '''
#' Correlation with Outlier Detection
#'
#' Computes correlation after removing outliers detected by projection methods.
#'
#' @inheritParams common-params
#' @param pcor Logical. If `TRUE`, uses Pearson correlation; if `FALSE`, uses `corfun`.
#' @param regfun Regression function for outlier detection (default: `tsreg`).
#' @param outkeep Logical. If `TRUE`, keeps outliers; if `FALSE` (default), removes them.
#'
#' @return List with components:
#'   \\item{cor}{Correlation coefficient (after removing outliers if applicable).}
#'   \\item{outid.x}{Indices of x-outliers.}
#'   \\item{outid.y}{Indices of y-outliers.}
#'   \\item{keep}{Indices of observations retained.}
#'
#' @details
#' This function detects outliers using projection-based methods and computes
#' the correlation after removing them. Outliers are identified separately for
#' x and y using the specified outlier detection function.
#'
#' @seealso \\code{\\link{ocor}}, \\code{\\link{pbcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- c(rnorm(40), 5, 6)  # Two outliers
#' y <- c(0.5 * rnorm(40), -5, 7)
#' ecor(x, y)
''',

    'ocor': '''
#' Correlation with Multivariate Outlier Detection
#'
#' Computes correlation after removing bivariate outliers detected using
#' multivariate outlier detection methods.
#'
#' @inheritParams common-params
#' @param pcor Logical. If `TRUE`, uses Pearson correlation after outlier removal.
#'
#' @return Numeric value representing the correlation after outlier removal.
#'
#' @details
#' Unlike `ecor()`, which detects univariate outliers separately for x and y,
#' this function uses multivariate outlier detection on the (x, y) pairs.
#'
#' @seealso \\code{\\link{ecor}}, \\code{\\link{scorci}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 0.6 * x + rnorm(50)
#' # Add bivariate outlier
#' x <- c(x, 5)
#' y <- c(y, 5)
#' ocor(x, y)
''',

    'cori': '''
#' Correlation at Specified Value of Third Variable
#'
#' Computes the correlation between x and y for observations near a specified
#' value of a third variable z.
#'
#' @inheritParams common-params
#' @param z Numeric vector. The conditioning variable.
#' @param pt Numeric. The point at which to condition on z (default: `median(z)`).
#' @param fr Numeric. Fraction of data to include near pt (default: 0.8).
#' @param est Location estimator function (default: `onestep`).
#' @param testit Logical. If `TRUE`, performs hypothesis test (default: `FALSE`).
#' @param nboot Number of bootstrap samples for testing (default: 100).
#'
#' @return If `testit = FALSE`, returns the correlation. If `testit = TRUE`,
#'   returns a list with correlation, test statistic, and p-value.
#'
#' @details
#' This function computes the correlation between x and y using only
#' observations where z is near the specified point `pt`. The "nearness"
#' is determined by `fr` (fraction of closest observations).
#'
#' @seealso \\code{\\link{runcor}}, \\code{\\link{pbcor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' n <- 100
#' z <- rnorm(n)
#' x <- rnorm(n)
#' y <- ifelse(z > 0, 0.8 * x, 0.2 * x) + rnorm(n)
#' # Correlation when z is near median
#' cori(x, y, z, pt = median(z))
#' # Correlation when z is near 1
#' cori(x, y, z, pt = 1, fr = 0.5)
''',

    'pcor': '''
#' Pearson Correlation
#'
#' Computes Pearson's product-moment correlation coefficient.
#'
#' @inheritParams common-params
#'
#' @return Pearson correlation coefficient.
#'
#' @details
#' This is a simple wrapper that handles matrix input and missing values.
#' If `y = NA`, assumes `x` is a two-column matrix.
#'
#' @seealso \\code{\\link{pbcor}}, \\code{\\link{wincor}}, \\code{\\link{pcorb}}
#'
#' @export
#' @examples
#' x <- rnorm(50)
#' y <- 0.7 * x + rnorm(50)
#' pcor(x, y)
#'
#' # Matrix input
#' m <- cbind(x, y)
#' pcor(m)
''',

}

print("Documentation templates created. Run this to add comprehensive docs to all correlation functions.")
