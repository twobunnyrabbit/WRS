# =============================================================================
# WRS: Winsorization Methods
# =============================================================================
#
# Core winsorization functions and utilities
#
# Functions in this module:
#   - win: Core winsorization function
#   - winmean: Winsorized mean
#   - winsd, winsd05, winsdN: Winsorized standard deviations
#   - winse, winci: Standard error and confidence interval for winsorized mean
#   - winvarN: Normalized winsorized variance
#   - winsorized: General winsorization with threshold
#   - WINCOR: Winsorized correlation matrix convenience function
#
# Note: Core utilities like winvar, winval, winall are in 00-utils-core.R
#       Winsorized correlation (wincor) is in correlation.R
#       Winsorized covariance (wincov) is in covariance.R
#       Winsorized regression (winreg) is in regression.R
#
# =============================================================================

#' Winsorized Standard Deviation
#'
#' Computes the standard deviation based on Winsorized values. The data are
#' Winsorized by replacing the lowest and highest values with specified quantiles,
#' then computing the standard deviation.
#'
#' @inheritParams common-params
#'
#' @return A numeric value: the Winsorized standard deviation.
#'
#' @details
#' This function computes `sqrt(winvar(x, tr))` where `winvar` is the Winsorized
#' variance. The Winsorized variance is more robust to outliers than the classical
#' variance because extreme values are replaced with less extreme values before
#' computing the variance.
#'
#' The amount of Winsorization is controlled by `tr`. For example, with `tr = 0.2`,
#' the lowest 20% of values are replaced with the value at the 20th percentile,
#' and the highest 20% are replaced with the value at the 80th percentile.
#'
#' @seealso \code{\link{winsdN}} for normalized version, \code{\link{winvar}} for
#' Winsorized variance, \code{\link{winmean}} for Winsorized mean
#'
#' @examples
#' # Sample data with outliers
#' x <- c(1, 2, 3, 4, 5, 100)
#'
#' # Standard deviation (affected by outlier)
#' sd(x)
#'
#' # Winsorized standard deviation (robust to outlier)
#' winsd(x, tr = 0.2)
#'
#' # With 10% Winsorization
#' winsd(x, tr = 0.1)
#'
#' @export
winsd<-function(x,tr=.2,na.rm=FALSE){
val=sqrt(winvar(x,tr=tr,na.rm=na.rm))
val
}

#' Winsorized Standard Deviation (Alternative Name)
#'
#' Computes the standard deviation based on Winsorized values. This function is
#' identical to \code{\link{winsd}} and is provided for backward compatibility.
#'
#' @inheritParams common-params
#'
#' @return A numeric value: the Winsorized standard deviation.
#'
#' @details
#' This function is identical to \code{\link{winsd}}. It computes the square root
#' of the Winsorized variance. Use \code{\link{winsd}} for new code.
#'
#' @seealso \code{\link{winsd}} for the preferred function name
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5, 100)
#' winsd05(x, tr = 0.2)
#'
#' @export
winsd05<-function(x,tr=.2,na.rm=FALSE){
val=sqrt(winvar(x,tr=tr,na.rm=na.rm))
val
}

#' Gamma Winsorized Mean
#'
#' Computes the gamma Winsorized mean, a robust measure of central tendency that
#' replaces extreme values with less extreme values before computing the mean.
#'
#' @inheritParams common-params
#'
#' @return A numeric value: the gamma Winsorized mean.
#'
#' @details
#' The gamma Winsorized mean works as follows:
#' \enumerate{
#'   \item Sort the data in ascending order
#'   \item Identify the lower and upper cutoff points based on the trim proportion `tr`
#'   \item Replace all values below the lower cutoff with the lower cutoff value
#'   \item Replace all values above the upper cutoff with the upper cutoff value
#'   \item Compute the mean of the Winsorized data
#' }
#'
#' For example, with `tr = 0.2` and n = 10 observations, the lowest value is
#' replaced with the 2nd smallest value, and the highest value is replaced with
#' the 9th value (2nd largest). The mean is then computed on these Winsorized values.
#'
#' This is a simple implementation. For more sophisticated Winsorized mean estimation
#' with missing value handling, use \code{\link{winmean}}.
#'
#' @note This function does not handle missing values. Use \code{\link{winmean}}
#' for automatic NA removal.
#'
#' @seealso \code{\link{winmean}} for Winsorized mean with NA handling,
#' \code{\link{winval}} for Winsorized values without computing mean
#'
#' @examples
#' # Sample data with outliers
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)
#'
#' # Regular mean (affected by outlier)
#' mean(x)
#'
#' # Gamma Winsorized mean (robust to outlier)
#' win(x, tr = 0.2)
#'
#' # Compare with different trim levels
#' win(x, tr = 0.1)
#' win(x, tr = 0.3)
#'
#' @export
win<-function(x,tr=.2){
#
#  Compute the gamma Winsorized mean for the data in the vector x.
#
#  tr is the amount of Winsorization
#
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
y<-ifelse(y<=xbot,xbot,y)
y<-ifelse(y>=xtop,xtop,y)
win<-mean(y)
win
}

#' Winsorized Mean
#'
#' Computes the Winsorized mean, a robust measure of central tendency that replaces
#' extreme values with less extreme values before computing the mean.
#'
#' @inheritParams common-params
#'
#' @return A numeric value: the Winsorized mean.
#'
#' @details
#' The Winsorized mean is computed by:
#' \enumerate{
#'   \item Optionally removing NA values if `na.rm = TRUE`
#'   \item Winsorizing the data using \code{\link{winval}}
#'   \item Computing the mean of the Winsorized values
#' }
#'
#' Winsorization replaces the most extreme values with less extreme values. With
#' `tr = 0.2`, the lowest 20% of values are replaced with the value at the 20th
#' percentile, and the highest 20% are replaced with the value at the 80th percentile.
#'
#' The Winsorized mean provides a balance between the arithmetic mean (sensitive to
#' outliers) and the trimmed mean (which discards extreme values). All observations
#' contribute to the Winsorized mean, but extreme values have their influence limited.
#'
#' @seealso \code{\link{win}} for gamma Winsorized mean, \code{\link{tmean}} for
#' trimmed mean, \code{\link{winse}} for standard error of Winsorized mean,
#' \code{\link{winci}} for confidence interval
#'
#' @examples
#' # Sample data with outliers
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)
#'
#' # Regular mean (affected by outlier)
#' mean(x)
#'
#' # Winsorized mean (robust to outlier)
#' winmean(x, tr = 0.2)
#'
#' # Compare with trimmed mean (which discards extremes)
#' mean(x, trim = 0.2)
#'
#' # Data with missing values
#' y <- c(1, 2, 3, NA, 5, 100)
#' winmean(y, tr = 0.2, na.rm = TRUE)
#'
#' @export
winmean<-function(x,tr=.2,na.rm=TRUE){
if(na.rm)x=elimna(x)
winmean<-mean(winval(x,tr))
winmean
}

#' Normalized Winsorized Variance
#'
#' Computes the Winsorized variance rescaled so that it equals 1 for data from a
#' standard normal distribution. This provides a variance estimator that is more
#' comparable to the classical variance under normality.
#'
#' @inheritParams common-params
#'
#' @return A numeric value: the normalized Winsorized variance.
#'
#' @details
#' The normalized Winsorized variance is computed by dividing the Winsorized variance
#' by a correction term that depends on the trim proportion `tr`. The correction term
#' is chosen so that for data from a standard normal distribution (mean 0, variance 1),
#' the normalized Winsorized variance equals 1 on average.
#'
#' For common trim proportions, pre-computed correction factors are used:
#' \itemize{
#'   \item `tr = 0`: correction = 1 (no Winsorization)
#'   \item `tr = 0.1`: correction = 0.6786546
#'   \item `tr = 0.2`: correction = 0.4120867
#' }
#'
#' For other trim proportions, the correction factor is computed numerically as:
#' \deqn{c = \int_{\Phi^{-1}(tr)}^{\Phi^{-1}(1-tr)} \phi(x)^2 dx + 2[\Phi^{-1}(tr)]^2 \cdot tr}
#' where \eqn{\phi} is the standard normal density and \eqn{\Phi} is the standard
#' normal CDF.
#'
#' @note Missing values are automatically removed before computation.
#'
#' @seealso \code{\link{winvar}} for Winsorized variance, \code{\link{winsdN}} for
#' normalized Winsorized standard deviation
#'
#' @examples
#' # Standard normal data
#' set.seed(123)
#' x <- rnorm(100)
#'
#' # Regular variance (close to 1)
#' var(x)
#'
#' # Winsorized variance (less than 1 due to Winsorization)
#' winvar(x, tr = 0.2)
#'
#' # Normalized Winsorized variance (rescaled to be close to 1)
#' winvarN(x, tr = 0.2)
#'
#' # With 10% Winsorization
#' winvarN(x, tr = 0.1)
#'
#' @export
winvarN<-function(x,tr=.2){
#
# rescale the Winsorized variance so that it equals one for the standard
# normal distribution
#
x=elimna(x)
cterm=NULL
if(tr==0)cterm=1
if(tr==0.1)cterm=0.6786546
if(tr==0.2)cterm=0.4120867
if(is.null(cterm))cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
bot=winvar(x,tr=tr)/cterm
bot
}

#' Standard Error of Winsorized Mean
#'
#' Estimates the standard error of the Winsorized mean using the Winsorized variance.
#'
#' @inheritParams common-params
#'
#' @return A numeric value: the estimated standard error of the Winsorized mean.
#'
#' @details
#' The standard error is computed using the formula:
#' \deqn{SE = \frac{(n-1) \sqrt{s_w^2}}{(h-1) \sqrt{n}}}
#' where:
#' \itemize{
#'   \item \eqn{n} is the sample size
#'   \item \eqn{s_w^2} is the Winsorized variance
#'   \item \eqn{h = n - 2 \lfloor tr \cdot n \rfloor} is the effective sample size
#'     after accounting for Winsorization
#' }
#'
#' This standard error can be used for constructing confidence intervals and
#' hypothesis tests for the Winsorized mean. See \code{\link{winci}} for computing
#' confidence intervals using this standard error.
#'
#' @note Missing values are automatically removed before computation.
#'
#' @seealso \code{\link{winmean}} for Winsorized mean, \code{\link{winci}} for
#' confidence interval, \code{\link{winvar}} for Winsorized variance
#'
#' @examples
#' # Sample data
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)
#'
#' # Winsorized mean and its standard error
#' winmean(x, tr = 0.2)
#' winse(x, tr = 0.2)
#'
#' # Compare with standard error of the mean
#' sd(x) / sqrt(length(x))
#'
#' # Construct approximate 95% CI manually
#' m <- winmean(x, tr = 0.2)
#' se <- winse(x, tr = 0.2)
#' c(m - 2*se, m + 2*se)
#'
#' @export
winse<-function(x,tr=.2){
#
# Estimate the standard error of the Winsorized mean
#
x=elimna(x)
n=length(x)
h=n-2*floor(tr*n)
top=(n-1)*sqrt(winvar(x,tr=tr))
bot=(h-1)*sqrt(n)
se=top/bot
se
}

#' Confidence Interval for Winsorized Mean
#'
#' Computes a confidence interval for the population Winsorized mean and performs
#' a one-sample hypothesis test.
#'
#' @inheritParams common-params
#'
#' @return A list with components:
#' \item{ci}{A vector of length 2 containing the lower and upper confidence limits}
#' \item{test.stat}{The t-test statistic for testing H0: Winsorized mean = null.value}
#' \item{p.value}{The two-sided p-value for the hypothesis test}
#'
#' @details
#' The confidence interval is computed using a Student's t distribution with
#' \eqn{df = n - 2 \lfloor tr \cdot n \rfloor - 1} degrees of freedom, where
#' \eqn{n} is the sample size.
#'
#' The confidence interval is:
#' \deqn{[\bar{X}_w - t_{1-\alpha/2, df} \cdot SE, \bar{X}_w + t_{1-\alpha/2, df} \cdot SE]}
#' where \eqn{\bar{X}_w} is the Winsorized mean and SE is the standard error from
#' \code{\link{winse}}.
#'
#' The test statistic for H0: \eqn{\mu_w = \mu_0} is:
#' \deqn{t = \frac{\bar{X}_w - \mu_0}{SE}}
#'
#' @note
#' By default, this function prints a reminder that the p-value tests against the
#' null value specified in `null.value` (default 0). Set `pr = FALSE` to suppress
#' this message.
#'
#' @seealso \code{\link{winmean}} for Winsorized mean, \code{\link{winse}} for
#' standard error, \code{\link{trimci}} for confidence interval for trimmed mean
#'
#' @examples
#' # Sample data
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)
#'
#' # 95% confidence interval for Winsorized mean
#' winci(x, tr = 0.2, pr = FALSE)
#'
#' # Test if Winsorized mean equals 5
#' winci(x, tr = 0.2, null.value = 5, pr = FALSE)
#'
#' # 90% confidence interval
#' winci(x, tr = 0.2, alpha = 0.10, pr = FALSE)
#'
#' # With different amount of Winsorization
#' winci(x, tr = 0.1, pr = FALSE)
#'
#' @export
winci<-function(x,tr=.2,alpha=.05,null.value=0,pr=TRUE){
#
#  Compute a 1-alpha confidence interval for the Winsorized mean
#
#  The default amount of  Winsorizing is tr=.2
#
if(pr){
print("The p-value returned by the this function is based on the")
print("null value specified by the argument null.value, which defaults to 0")
}
x<-elimna(x)
se<-winse(x,tr=tr)
df<-length(x)-2*floor(tr*length(x))-1
trimci<-winmean(x,tr)-qt(1-alpha/2,df)*se
trimci[2]<-winmean(x,tr)+qt(1-alpha/2,df)*se
test<-(winmean(x,tr)-null.value)/se
sig<-2*(1-pt(abs(test),df))
list(ci=trimci,test.stat=test,p.value=sig)
}

#' Winsorize Data with Fixed Threshold
#'
#' Winsorizes data by replacing values outside a fixed threshold with the threshold
#' values. Unlike quantile-based Winsorization, this uses absolute thresholds.
#'
#' @param x A numeric vector to be Winsorized.
#' @param a Threshold multiplier (default: 1.5). Values outside \eqn{\pm a \cdot \sigma}
#'   are replaced with \eqn{\pm a \cdot \sigma}.
#' @param sigma Scale parameter for the threshold (default: 1). Typically the standard
#'   deviation or a robust scale estimator like MAD.
#'
#' @return A numeric vector of the same length as `x` with extreme values Winsorized.
#'
#' @details
#' This function implements threshold-based Winsorization:
#' \itemize{
#'   \item Values greater than \eqn{a \cdot \sigma} are replaced with \eqn{a \cdot \sigma}
#'   \item Values less than \eqn{-a \cdot \sigma} are replaced with \eqn{-a \cdot \sigma}
#'   \item Other values remain unchanged
#' }
#'
#' This differs from quantile-based Winsorization (as in \code{\link{win}} and
#' \code{\link{winmean}}) which replaces a fixed proportion of extreme values.
#' Threshold-based Winsorization is useful when you have a natural scale for what
#' constitutes an "extreme" value.
#'
#' Common choices for `sigma`:
#' \itemize{
#'   \item Standard deviation: `sd(x)`
#'   \item MAD (median absolute deviation): `mad(x)`
#'   \item A domain-specific scale (e.g., measurement precision)
#' }
#'
#' @note This function does not remove missing values. Ensure `x` contains no NAs
#' before calling this function.
#'
#' @seealso \code{\link{win}} for quantile-based Winsorization, \code{\link{winval}}
#' for Winsorized values
#'
#' @examples
#' # Sample data
#' x <- c(-10, -2, -1, 0, 1, 2, 3, 15)
#'
#' # Winsorize values outside ±1.5
#' winsorized(x, a = 1.5, sigma = 1)
#'
#' # Using standard deviation as scale
#' winsorized(x, a = 2, sigma = sd(x))
#'
#' # Using MAD (median absolute deviation) as robust scale
#' winsorized(x, a = 3, sigma = mad(x))
#'
#' @export
winsorized<- function(x,a=1.5,sigma=1) {
s<-sigma
newx<-x
indp<-x>(a*s)
newx[indp]<-(a*s)
indn<- x<(a*-s)
newx[indn]<- (-a*s)
newx
}

#' Winsorized Correlation Matrix
#'
#' Convenience function to compute only the Winsorized correlation matrix from
#' multivariate data.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or data frame. Each column represents a variable.
#'
#' @return A correlation matrix computed from Winsorized values.
#'
#' @details
#' This is a convenience wrapper that calls \code{\link{winall}} and extracts only
#' the correlation matrix component. The function \code{\link{winall}} computes
#' both the Winsorized covariance and correlation matrices; this function returns
#' only the correlation matrix.
#'
#' The Winsorized correlation is more robust to outliers than the Pearson correlation
#' because extreme values are replaced with less extreme values before computing
#' correlations.
#'
#' For more detailed information and additional output (including covariance matrix),
#' use \code{\link{winall}} directly.
#'
#' @seealso \code{\link{winall}} for Winsorized covariance and correlation,
#' \code{\link{wincor}} for Winsorized correlation between two vectors,
#' \code{\link{pbcor}} for percentage bend correlation
#'
#' @examples
#' # Multivariate data
#' set.seed(123)
#' X <- matrix(rnorm(100), ncol = 5)
#'
#' # Pearson correlation matrix
#' cor(X)
#'
#' # Winsorized correlation matrix (robust to outliers)
#' WINCOR(X, tr = 0.2)
#'
#' # With data frame
#' df <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
#' WINCOR(df, tr = 0.1)
#'
#' @export
WINCOR<-function(x,tr=.2){
#
# For convenience, compute Winsorized correlation matrix only.
#
a=winall(x,tr=tr)$cor
a
}

#' Normalized Winsorized Standard Deviation
#'
#' Computes the Winsorized standard deviation rescaled so that it estimates the
#' population standard deviation under normality. This provides a standard deviation
#' estimator that is more comparable to the classical standard deviation.
#'
#' @inheritParams common-params
#'
#' @return A numeric value: the normalized Winsorized standard deviation.
#'
#' @details
#' The normalized Winsorized standard deviation is computed by dividing the
#' Winsorized standard deviation by a correction factor. The correction factor
#' is the square root of the correction term used in \code{\link{winvarN}}.
#'
#' For data from a standard normal distribution (mean 0, variance 1), the normalized
#' Winsorized standard deviation equals 1 on average, making it comparable to the
#' classical standard deviation under normality while being more robust to outliers.
#'
#' The correction factor for `tr = 0` (no Winsorization) is 1. For other values of
#' `tr`, the correction factor is computed as:
#' \deqn{c = \sqrt{\int_{\Phi^{-1}(tr)}^{\Phi^{-1}(1-tr)} \phi(x)^2 dx + 2[\Phi^{-1}(tr)]^2 \cdot tr}}
#' where \eqn{\phi} is the standard normal density and \eqn{\Phi} is the standard
#' normal CDF.
#'
#' @note Missing values are automatically removed before computation.
#'
#' @seealso \code{\link{winsd}} for Winsorized standard deviation without normalization,
#' \code{\link{winvarN}} for normalized Winsorized variance
#'
#' @examples
#' # Standard normal data
#' set.seed(123)
#' x <- rnorm(100)
#'
#' # Regular standard deviation (close to 1)
#' sd(x)
#'
#' # Winsorized standard deviation (less than 1 due to Winsorization)
#' winsd(x, tr = 0.2)
#'
#' # Normalized Winsorized standard deviation (rescaled to be close to 1)
#' winsdN(x, tr = 0.2)
#'
#' # With contaminated data (10% outliers)
#' y <- c(rnorm(90), rnorm(10, mean = 5, sd = 3))
#' sd(y)        # Affected by outliers
#' winsdN(y, tr = 0.2)  # Robust to outliers
#'
#' @export
winsdN<-function(x,tr=.2){
#
# Rescale a Winsorized standard deviation so that it estimates
# the population standard deviation under normality.
#
x=elimna(x)
e=winsd(x,tr=tr)
if(tr==0)cterm=1
else
cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
cterm=sqrt(cterm)
e=e/cterm
e
}
