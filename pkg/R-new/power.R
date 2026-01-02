# power.R - Power Analysis Functions
# Part of WRS package refactoring
#
# Functions for statistical power analysis and sample size estimation.
# These functions estimate statistical power for various robust methods
# including trimmed means, ANOVA, ANCOVA, and chi-square tests.
#
# Dependencies: trimse, winvar, yuen, elimna, tsreg, pbcor, outmgvf, pbvar, ancmg1
# External packages: pwr (for power.chisq.test)

# ============================================================================
# Power Analysis Functions (10 functions)
# ============================================================================

#' Estimate Explanatory Power
#'
#' Estimates the explanatory power of predictor(s) x for predicting y, measuring
#' how much variability in y is explained by the regression on x.
#'
#' @param x Predictor variable(s). Can be a vector or matrix.
#' @param y Response variable (numeric vector).
#' @param pcor Logical. If \code{TRUE}, uses Pearson correlation squared. If \code{FALSE},
#'   uses robust correlation squared (default: \code{FALSE}).
#' @param regfun Regression function to use (default: \code{\link{tsreg}} for Theil-Sen).
#' @param corfun Correlation function to use when \code{pcor=FALSE} (default: \code{\link{pbcor}}
#'   for percentage bend correlation).
#' @param outkeep Logical. If \code{FALSE}, removes outliers before computing explanatory
#'   power (default: \code{FALSE}).
#' @param outfun Function for detecting outliers (default: \code{\link{outmgvf}}).
#' @param varfun Function for computing variance (default: \code{\link{pbvar}} for
#'   percentage bend variance).
#' @param op Logical. If \code{TRUE}, uses correlation-based measure. If \code{FALSE},
#'   uses variance ratio (default: \code{TRUE}).
#'
#' @return Estimated explanatory power (scalar between 0 and 1).
#'
#' @details
#' This function estimates explanatory power using two possible approaches:
#' \itemize{
#'   \item \strong{op = TRUE}: Computes correlation between fitted values and observed y,
#'     then squares it (R-squared analog)
#'   \item \strong{op = FALSE}: Computes variance of fitted values divided by variance of y
#' }
#'
#' The procedure:
#' \enumerate{
#'   \item Fits regression model using \code{regfun}
#'   \item Optionally removes outliers using \code{outfun}
#'   \item Computes explanatory power using fitted values
#' }
#'
#' Missing values are automatically removed.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{tsreg}}, \code{\link{pbcor}}, \code{\link{pbvar}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 0.5 * x + rnorm(50)
#' epow(x, y)
#'
#' # With multiple predictors
#' x2 <- cbind(x, rnorm(50))
#' epow(x2, y)
epow<-function(x,y,pcor=FALSE,regfun=tsreg,corfun=pbcor,outkeep=FALSE,outfun=outmgvf,varfun=pbvar,op=TRUE){
#
# Estimate the explanatory power between x and y
#
xx<-elimna(cbind(x,y))
pval<-1
if(is.matrix(x))pval<-ncol(x)
pp<-pval+1
x<-xx[,1:pval]
y<-xx[,pp]
x<-as.matrix(x)
flag<-rep(TRUE,nrow(x))
temp<-regfun(x,y)
ip<-ncol(x)+1
yhat<-y-temp$res
if(!outkeep){
temp<-outfun(cbind(x,y),plotit=FALSE)$out.id
flag[temp]<-FALSE
}
epow1<-varfun(yhat[flag])/varfun(y[flag])
if(pcor)epow2<-cor(yhat[flag],y[flag])^2
if(!pcor)epow2<-corfun(yhat[flag],y[flag])$cor^2
if(op)est<-epow2
if(!op)est<-epow1
est
}

#' Power for One-Sample Student's T-Test (One-Sided)
#'
#' Computes the statistical power of a one-sided, one-sample Student's t-test
#' using the Kraemer-Paik approximation for the noncentral t-distribution.
#'
#' @param n Sample size (positive integer).
#' @param Del Standardized effect size: \eqn{(\mu_0 - \mu_1) / \sigma}, where
#'   \eqn{\mu_0} is the hypothesized value and \eqn{\mu_1} is the true mean.
#' @param alpha Significance level / Type I error probability (must be between 0 and 1).
#'
#' @return Estimated statistical power (scalar between 0 and 1).
#'
#' @details
#' This function uses the Kraemer-Paik (1979) approximation to compute power for
#' a one-sided, one-sample t-test. The approximation provides accurate power estimates
#' without requiring direct computation with the noncentral t-distribution.
#'
#' The effect size \code{Del} is computed as the difference between the null hypothesis
#' mean and true mean, divided by the population standard deviation.
#'
#' @references
#' Kraemer, H. C., & Paik, M. (1979). A central t approximation to the noncentral
#' t distribution. Technometrics, 21(3), 357-360.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{powt1est}}, \code{\link{pow2an}}
#'
#' @export
#' @examples
#' # Power for n=30, effect size = 0.5, alpha = 0.05
#' pow1(n = 30, Del = 0.5, alpha = 0.05)
#'
#' # Compare different sample sizes
#' sapply(c(20, 30, 50, 100), function(n) pow1(n, Del = 0.5, alpha = 0.05))
pow1<-function(n,Del,alpha){
#
#  Determine power of Student's T in the
#  one-sided, one-sample case where
#
#  n=sample size
#  Del=(mu0-mu1)/sigma
#  alpha=Type I error probability
#  mu0 is hypothesized value
#  mu1 is some non-null value for the mean.
#
Del<-abs(Del)
if(alpha<=0 || alpha>=1)stop("alpha must be between 0 and 1")
K11<-1-alpha
K5<-sqrt(n)*Del
#  Next, use the Kraemer-Paik (1979, Technometrics, 21, 357-360)
#  approximation of the noncentral T.
K6<-n-1
K14<-qt(K11,K6)
K7<-K14*sqrt(1+K5*K5/K6)
K8<-K5*sqrt(1+K14*K14/K6)
K9<-K7-K8
pow1<-1-pt(K9,K6)
pow1
}

#' Power Analysis for Two-Sample Trimmed Means
#'
#' Performs power analysis for comparing 20\% trimmed means of two independent groups
#' using the percentile bootstrap method. Plots power as a function of effect size (delta).
#'
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param ci Logical. If \code{TRUE}, computes bootstrap confidence band for the power
#'   curve (default: \code{FALSE}).
#' @param plotit Logical. If \code{TRUE}, creates a plot of power vs. delta (default: \code{TRUE}).
#' @param nboot Number of bootstrap samples for confidence band computation (default: 800).
#'
#' @return A list with components:
#'   \item{delta}{Vector of effect sizes (delta values) evaluated.}
#'   \item{power}{Estimated power at each delta value.}
#'   \item{lowp}{Lower 95\% confidence bound for power at each delta (only if \code{ci=TRUE}).}
#'
#' @details
#' This function estimates statistical power for Yuen's test (\code{\link{yuen}}) comparing
#' two independent groups based on 20\% trimmed means.
#'
#' The procedure:
#' \enumerate{
#'   \item Computes standard error from the data using \code{\link{yuen}}
#'   \item Evaluates power at 15 delta values ranging from 0 to 3.5 times the SE
#'   \item Optionally bootstraps to obtain confidence bounds for the power estimates
#' }
#'
#' The plot shows power on the y-axis and effect size (delta) on the x-axis. If \code{ci=TRUE},
#' a dashed line shows the lower 95\% confidence bound.
#'
#' Missing values are removed before analysis.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{yuen}}, \code{\link{powest}}, \code{\link{powt1an}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 0, sd = 1)
#' y <- rnorm(30, mean = 0.3, sd = 1)
#'
#' # Power curve without confidence band
#' result <- pow2an(x, y, ci = FALSE)
#'
#' # With bootstrap confidence band (slower)
#' \dontrun{
#' result_ci <- pow2an(x, y, ci = TRUE, nboot = 500)
#' }
pow2an<-function(x,y,ci=FALSE,plotit=TRUE,nboot=800){
#
# Do a power analysis when comparing the 20% trimmed
# means of two independent groups with the percentile
# bootstrap method.
#
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
lp<-NA
se<-yuen(x,y)$se
gval<-NA
dv<-seq(0,3.5*se,length=15)
for(i in 1:length(dv)){
gval[i]<-powest(x,y,dv[i],se)
}
if(!ci){
if(plotit){
plot(dv,gval,type="n",xlab="delta",ylab="power")
lines(dv,gval)
}}
if(ci){
print("Taking bootstrap samples. Please wait.")
datax <- matrix(sample(x, size = length(x) * nboot, replace = TRUE),
                nrow = nboot)
datay <- matrix(sample(y, size = length(y) * nboot, replace = TRUE),
                nrow = nboot)
pboot<-matrix(NA,ncol=15,nrow=nboot)
for(i in 1:nboot){
se<-yuen(datax[i,],datay[i,])$se
for(j in 1:length(dv)){
pboot[i,j]<-powest(x,y,dv[j],se)
}}
ll<-floor(.05*nboot+.5)
for(i in 1:15){
temp<-sort(pboot[,i])
lp[i]<-temp[ll]
}
plot(c(dv,dv),c(gval,lp),type="n",xlab="delta",ylab="power")
lines(dv,gval)
lines(dv,lp,lty=2)
}
list(delta=dv,power=gval,lowp=lp)
}

#' Power Analysis for Chi-Square Test
#'
#' Computes power, sample size, effect size, or significance level for a chi-square
#' test of association. This is a wrapper for \code{pwr.chisq.test} from the pwr package.
#'
#' @param w Effect size (Cohen's w). \code{NULL} if solving for effect size.
#' @param N Total sample size. \code{NULL} if solving for sample size.
#' @param df Degrees of freedom: \code{(r - 1) * (c - 1)} for an r × c contingency table.
#'   \code{NULL} if solving for df.
#' @param sig.level Significance level (Type I error probability) (default: 0.05).
#'   \code{NULL} if solving for significance level.
#' @param power Statistical power (1 - Type II error probability). \code{NULL} if
#'   solving for power.
#'
#' @return An object of class \code{power.htest} containing the computed parameters.
#'
#' @details
#' This function requires the \code{pwr} package. Exactly one of \code{w}, \code{N},
#' \code{df}, \code{sig.level}, and \code{power} must be \code{NULL}; that parameter
#' will be computed from the others.
#'
#' Cohen's w effect size conventions:
#' \itemize{
#'   \item Small: w = 0.1
#'   \item Medium: w = 0.3
#'   \item Large: w = 0.5
#' }
#'
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.).
#' Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link[pwr]{pwr.chisq.test}}
#'
#' @export
#' @examples
#' # Power for a 2x3 table with w=0.3, N=100, alpha=0.05
#' power.chisq.test(w = 0.3, N = 100, df = (2-1)*(3-1))
#'
#' # Sample size needed for power=0.80
#' power.chisq.test(w = 0.3, power = 0.80, df = 2, sig.level = 0.05)
power.chisq.test<-function(w = NULL, N = NULL, df = NULL,
 sig.level = 0.05, power = NULL){
library(pwr)
res=pwr.chisq.test(w=w,N=N,df=df,sig.level=sig.level, power = power)
res
}

#' @keywords internal
powest<-function(x=NA,y=NA,delta=0,se=NA,wv1=NA,wv2=NA,n1=NA,n2=NA){
#
# wv1 = Winsorized variance for group 1
# wv2 = Winsorized variance for group 2
#
# Only 20% trimming is allowed.
#
tr<-.2
if(is.na(se)){
if(is.na(wv1)){
h1 <- length(x) - 2 * floor(tr * length(x))
h2 <- length(y) - 2 * floor(tr * length(y))
q1 <- ((length(x) - 1) * winvar(x, tr))/(h1 * (h1 - 1))
q2 <- ((length(y) - 1) * winvar(y, tr))/(h2 * (h2 - 1))
}
if(!is.na(wv1)){
if(is.na(n1))stop("Need to specify sample size for group 1")
if(is.na(n2))stop("Need to specify sample size for group 2")
h1<-n1-2*floor(tr*n1)
h2<-n2-2*floor(tr*n2)
q1<-(n1-1)*wv1/(h1*(h1-1))
q2<-(n2-1)*wv2/(h2*(h2-1))
}
se<-sqrt(q1+q2)
}
ygam<-sqrt(2*.01155)*c(0:35)/8
pow<-c(500.0,540.0,607.0, 706.0, 804.0,981.0,1176.0,1402.0,1681.0, 2008.0,
   2353.0, 2769.0, 3191.0, 3646.0, 4124.0, 4617.0, 5101.0, 5630.0,
   6117.0, 6602.0, 7058.0, 7459.0, 7812.0, 8150.0, 8479.0, 8743.0,
   8984.0, 9168.0, 9332.0, 9490.0, 9607.0, 9700.0, 9782.0, 9839.0,
   9868.0)/10000
flag<-(delta==0 & se==0)
if(flag)powest<-.05
else{
chk<-floor(8*delta/se)+1
chk1<-chk+1
gval<-delta/se
d1<-(gval-(chk-1)/8)*8
if(chk > length(pow))powest<-1
if(chk == length(pow))pow[chk1]<-1
if(chk <= length(pow))
powest<-pow[chk]+d1*(pow[chk1]-pow[chk])
}
powest
}

#' Power Analysis for One-Sample Trimmed Mean
#'
#' Performs power analysis for testing a one-sample 20\% trimmed mean using the
#' percentile bootstrap method. Plots power as a function of effect size (delta).
#'
#' @param x Numeric vector containing the sample data.
#' @param ci Logical. If \code{TRUE}, computes bootstrap confidence band for the power
#'   curve (default: \code{FALSE}).
#' @param plotit Logical. If \code{TRUE}, creates a plot of power vs. delta (default: \code{TRUE}).
#' @param nboot Number of bootstrap samples for confidence band computation (default: 800).
#'
#' @return A list with components:
#'   \item{delta}{Vector of effect sizes (delta values) evaluated.}
#'   \item{power}{Estimated power at each delta value.}
#'   \item{lowp}{Lower 95\% confidence bound for power at each delta (only if \code{ci=TRUE}).}
#'
#' @details
#' This function estimates statistical power for a one-sample test of the 20\% trimmed mean
#' using the percentile bootstrap method.
#'
#' The procedure:
#' \enumerate{
#'   \item Computes trimmed standard error using \code{\link{trimse}}
#'   \item Evaluates power at 15 delta values ranging from 0 to 3.5 times the SE
#'   \item Optionally bootstraps to obtain confidence bounds for the power estimates
#' }
#'
#' The plot shows power on the y-axis and effect size (delta) on the x-axis. If \code{ci=TRUE},
#' a dashed line shows the lower 95\% confidence bound.
#'
#' Missing values are removed before analysis.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{trimse}}, \code{\link{powt1est}}, \code{\link{pow2an}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 0.5, sd = 1)
#'
#' # Power curve without confidence band
#' result <- powt1an(x, ci = FALSE)
#'
#' # With bootstrap confidence band (slower)
#' \dontrun{
#' result_ci <- powt1an(x, ci = TRUE, nboot = 500)
#' }
powt1an<-function(x,ci=FALSE,plotit=TRUE,nboot=800){
#
# Do a power analysis for the one-sample case with 20% trimmed
# mean and when the percentile bootstrap is to be used to test
# hypoltheses.
#
x<-x[!is.na(x)]
lp<-NA
se<-trimse(x)
gval<-NA
dv<-seq(0,3.5*se,length=15)
for(i in 1:length(dv)){
gval[i]<-powest(x,rep(0,5),dv[i],se)
}
if(!ci){
if(plotit){
plot(dv,gval,type="n",xlab="delta",ylab="power")
lines(dv,gval)
}}
if(ci){
set.seed(2)
print("Taking bootstrap samples. Please wait.")
datax <- matrix(sample(x, size = length(x) * nboot, replace = TRUE),
                nrow = nboot)
pboot<-matrix(NA,nrow=nboot,ncol=length(dv))
for(i in 1:nboot){
se<-trimse(datax[i,])
for(j in 1:length(dv)){
pboot[i,j]<-powest(x,rep(0,5),dv[j],se)
}}
ll<-floor(.05*nboot+.5)
for(i in 1:15){
temp<-sort(pboot[,i])
lp[i]<-temp[ll]
}
plot(c(dv,dv),c(gval,lp),type="n",xlab="delta",ylab="power")
lines(dv,gval)
lines(dv,lp,lty=2)
}
list(delta=dv,power=gval,lowp=lp)
}

#' Estimate Power for One-Sample Trimmed Mean at Specific Delta
#'
#' Estimates statistical power for testing a one-sample 20\% trimmed mean at a
#' specified effect size (delta), with optional bootstrap confidence interval.
#'
#' @param x Numeric vector containing the sample data.
#' @param delta Effect size (difference from null hypothesis in original units) (default: 0).
#' @param ci Logical. If \code{TRUE}, computes a bootstrap confidence interval for the
#'   power estimate (default: \code{FALSE}).
#' @param nboot Number of bootstrap samples for confidence interval (default: 800).
#'
#' @return A list with components:
#'   \item{est.power}{Estimated power at the specified delta.}
#'   \item{ci}{Lower 95\% confidence bound for power (only if \code{ci=TRUE}).}
#'
#' @details
#' This function estimates power at a single specified effect size for a one-sample
#' test of the 20\% trimmed mean using the percentile bootstrap method. Only 20\%
#' trimming is allowed (hardcoded).
#'
#' The procedure:
#' \enumerate{
#'   \item Computes trimmed standard error using \code{\link{trimse}}
#'   \item Estimates power using \code{\link{powest}}
#'   \item If \code{ci=TRUE}, bootstraps to obtain lower 95\% confidence bound
#' }
#'
#' Unlike \code{\link{powt1an}} which evaluates power across multiple delta values,
#' this function focuses on a single delta and optionally provides a confidence interval.
#'
#' Missing values are removed before analysis.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{powt1an}}, \code{\link{powest}}, \code{\link{trimse}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 0.5, sd = 1)
#'
#' # Power at delta = 0.5
#' powt1est(x, delta = 0.5)
#'
#' # With bootstrap confidence interval (slower)
#' \dontrun{
#' powt1est(x, delta = 0.5, ci = TRUE, nboot = 500)
#' }
powt1est<-function(x,delta=0,ci=FALSE,nboot=800){
#
# Estimate power for a given value of delta
#
# Only 20% trimming is allowed.
#
temp1<-powest(x,rep(0,5),delta,se=trimse(x))
if(ci){
set.seed(2)
pboot<-NA
datay<-rep(0,5)
print("Taking bootstrap samples. Please wait.")
datax <- matrix(sample(x, size = length(x) * nboot, replace = TRUE
                        ), nrow = nboot)
for(i in 1:nboot) {
se <- trimse(datax[i,  ])
pboot[i] <- powest(x, rep(0,5), delta, se)
}
temp <- sort(pboot)
}
ll<-floor(0.05 * nboot + 0.5)
list(est.power=temp1,ci=temp[ll])
}

#' Power Analysis Using Normal Distribution Approximation
#'
#' Computes statistical power for a two-sided test using a normal distribution
#' approximation, given sample size, effect size, and variance.
#'
#' @param n Sample size (positive integer).
#' @param alpha Significance level (Type I error probability) (default: 0.05).
#' @param del Effect size: difference between null and alternative hypothesis
#'   values (in original units). \code{NULL} if not specified.
#' @param var Population variance. \code{NULL} if not specified.
#'
#' @return A list with component:
#'   \item{power}{Estimated statistical power (scalar between 0 and 1).}
#'
#' @details
#' This function uses a normal distribution approximation to compute power for a
#' two-sided hypothesis test. It assumes the test statistic follows a normal
#' distribution under both null and alternative hypotheses.
#'
#' The power is computed as:
#' \deqn{P(Z < -z_{\alpha/2} - \sqrt{n} \cdot \delta / \sigma) + P(Z > z_{\alpha/2} - \sqrt{n} \cdot \delta / \sigma)}
#'
#' where \eqn{z_{\alpha/2}} is the critical value, \eqn{\delta} is the effect size,
#' and \eqn{\sigma} is the population standard deviation.
#'
#' This is a general-purpose power function suitable when normality assumptions
#' hold. For robust methods with trimmed means, see \code{\link{powt1an}} or
#' \code{\link{pow2an}}.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{pow1}}, \code{\link{pow2an}}, \code{\link{powt1an}}
#'
#' @export
#' @examples
#' # Power for n=50, effect size=0.5, variance=1, alpha=0.05
#' z.power(n = 50, del = 0.5, var = 1, alpha = 0.05)
#'
#' # Compare different sample sizes
#' sapply(c(30, 50, 100, 200), function(n) z.power(n, del = 0.5, var = 1)$power)
z.power<-function(n,alpha=.05,del=NULL,var=NULL){
 q=qnorm(1-alpha/2)
 sig=sqrt(var)
 p1=pnorm(0-q-(sqrt(n)*del)/sig)
 p2=1-pnorm(q-(sqrt(n)*del)/sig)
 p=p1+p2
 list(power=p)
 }
