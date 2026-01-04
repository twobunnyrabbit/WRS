# WRS Package - Regression Methods
# Extracted from Rallfun-v45.R
#
# This module contains robust regression methods including:
# - Theil-Sen regression (tsreg, tshdreg, etc.)
# - Least trimmed squares (ltsreg, LMSreg, MMreg)
# - M-regression (chreg, bmreg, bireg, winreg)
# - Outlier-pruned regression (opreg, mopreg)
# - Skipped regression (snmreg)
# - Depth-based regression (depreg, mdepreg)
# - Quantile regression basics (qreg, Qreg)
# - Regression inference (regci, regtest, lintest)
# - Two-group comparisons (difreg, reg2ci)
# - One-way regression ANOVA (reg1way)
#
# Total functions: 98
# Extraction date: 2025-12-30


# ============================================================================
# lintestMC
# ============================================================================

#' Test Linearity of Regression Surface (Parallel Processing)
#'
#' Tests the hypothesis that the regression surface is a plane using parallel
#' processing via \code{mclapply}. This is the multicore version of
#' \code{\link{lintest}}.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param regfun Robust regression function to use (default: \code{tsreg}).
#'   If using \code{tshdreg}, be sure to include \code{RES=TRUE}.
#' @inheritParams common-params
#'
#' @return A list with components:
#' \describe{
#'   \item{dstat}{The Kolmogorov-Smirnov type test statistic}
#'   \item{wstat}{The Cramér-von Mises type test statistic}
#'   \item{p.value.d}{Bootstrap p-value for the KS-type statistic}
#'   \item{p.value.w}{Bootstrap p-value for the CvM-type statistic}
#' }
#'
#' @details
#' This function implements the linearity test of Stute, Thies, and Zhu (1998)
#' using parallel processing for the bootstrap. The test assesses whether the
#' regression function is a plane (linear in all predictors).
#'
#' The method constructs an empirical process based on residuals and compares
#' it to its expected value under linearity. Two test statistics are computed:
#' \itemize{
#'   \item \code{dstat}: Maximum absolute deviation (Kolmogorov-Smirnov type)
#'   \item \code{wstat}: Mean squared deviation (Cramér-von Mises type)
#' }
#'
#' Bootstrap is used to approximate the null distribution. The function uses
#' \code{mclapply} for parallel processing, making it faster than
#' \code{\link{lintest}} for large datasets or many bootstrap samples.
#'
#' **Important**: When using \code{tshdreg} as the regression function, you
#' must include \code{RES=TRUE} to ensure residuals are available.
#'
#' @references
#' Stute, W., Thies, S., & Zhu, L.-X. (1998). Model checks for regression: An
#' innovation process approach. \emph{Annals of Statistics}, 26, 1916-1934.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{lintest}} for non-parallel version,
#'   \code{\link{regtest}} for testing specific parameters
#'
#' @examples
#' \dontrun{
#' # Linear relationship (should not reject)
#' set.seed(123)
#' x <- matrix(rnorm(200), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(100)
#' lintestMC(x, y, nboot=500)
#'
#' # Nonlinear relationship (should reject)
#' y2 <- x[,1]^2 + rnorm(100)
#' lintestMC(x, y2, nboot=500)
#'
#' # Using tshdreg
#' lintestMC(x, y, regfun=tshdreg, RES=TRUE, nboot=500)
#' }
#'
#' @export
lintestMC<-function(x,y,regfun=tsreg,nboot=500,alpha=.05,xout=FALSE,outfun=out,...){
#
# Test the hypothesis that the regression surface is a plane.
# Stute et al. (1998, JASA, 93, 141-149).
#
set.seed(2)
if(identical(regfun,tshdreg))print('When using tshdreg, be sure to include RES=TRUE')
#if(identical(regfun,Qreg))print('When using Qreg, be sure to include res.vals=TRUE')
x<-as.matrix(x)
d<-ncol(x)
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
x<-as.matrix(x)
y<-temp[,d+1]
if(xout){
flag<-outfun(x)$keep
x<-x[flag,]
x<-as.matrix(x)
y<-y[flag]
}
mflag<-matrix(NA,nrow=length(y),ncol=length(y))
for (j in 1:length(y)){
for (k in 1:length(y)){
mflag[j,k]<-(sum(x[j,]<=x[k,])==ncol(x))
}
}
reg<-regfun(x,y,...)
yhat<-y-reg$residuals
print("Taking bootstrap samples, please wait.")
data<-matrix(runif(length(y)*nboot),nrow=nboot)
data<-sqrt(12)*(data-.5) # standardize the random numbers.
data=listm(t(data))
rvalb<-mclapply(data,lintests1,yhat,reg$residuals,mflag,x,regfun,mc.preschedule=TRUE,...)
# An n x nboot matrix of R values
rvalb=matl(rvalb)
rvalb<-rvalb/sqrt(length(y))
dstatb<-apply(abs(rvalb),2,max)
wstatb<-apply(rvalb^2,2,mean)
# compute test statistic
v<-c(rep(1,length(y)))
rval<-lintests1(v,yhat,reg$residuals,mflag,x,regfun,...)
rval<-rval/sqrt(length(y))
dstat<-max(abs(rval))
wstat<-mean(rval^2)
ib<-round(nboot*(1-alpha))
p.value.d<-1-sum(dstat>=dstatb)/nboot
p.value.w<-1-sum(wstat>=wstatb)/nboot
list(dstat=dstat,wstat=wstat,p.value.d=p.value.d,p.value.w=p.value.w)
}

# ============================================================================
# bireg
# ============================================================================

#' Biweight Midregression
#'
#' Computes a robust regression using biweight M-estimation (midregression). This
#' method combines biweight covariances of predictors with the dependent variable
#' to produce robust slope estimates, then iteratively refines them.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param iter Maximum number of iterations for convergence (default: 20).
#' @param bend Bending constant for biweight estimator (default: 1.28).
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{residuals}{Residuals from the final fit.}
#' }
#'
#' @details
#' The biweight midregression is computed by:
#' 1. Computing biweight M-estimates of location for each predictor
#' 2. Computing biweight covariances between predictors and outcome
#' 3. Solving for initial slopes
#' 4. Iteratively refining estimates until convergence (change < 0.0001)
#'
#' If convergence fails within \code{iter} iterations, a warning message is displayed.
#'
#' @family robust regression functions
#' @export
#' @examples
#' \dontrun{
#' # Simple biweight regression
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2*x + rnorm(50, sd=0.5)
#' bireg(x, y)
#'
#' # Multiple regression with outliers
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50, sd=0.3)
#' y[1:5] <- y[1:5] + 5  # add outliers
#' bireg(x, y)
#' }
bireg<-function(x,y,iter=20,bend=1.28){
#
# Compute a biweight midregression equation
# The predictors are assumed to be stored in the n by p matrix x.
#
x<-as.matrix(x)
ma<-matrix(0,ncol(x),1)
m<-matrix(0,ncol(x),ncol(x))
mvals<-apply(x,2,mest,bend)
for (i in 1:ncol(x)){
ma[i,1]<-bicov(x[,i],y)
for (j in 1:ncol(x))m[i,j]<-bicov(x[,i],x[,j])
}
slope<-solve(m,ma)
b0<-mest(y,bend)-sum(slope%*%mvals)
for(it in 1:iter){
res<-y-x%*%slope-b0
for (i in 1:ncol(x))ma[i,1]<-bicov(x[,i],res)
slopeadd<-solve(m,ma)
b0add<-mest(res,bend)-sum(slopeadd%*%mvals)
if(max(abs(slopeadd),abs(b0add)) <.0001)break
slope<-slope+slopeadd
b0<-b0+b0add
}
if(max(abs(slopeadd),abs(b0add)) >=.0001)
paste("failed to converge in",iter,"iterations")
list(coef=c(b0,slope),residuals=res)
}

# ============================================================================
# chreg
# ============================================================================

#' Coakley-Hettmansperger Robust Regression
#'
#' Computes Coakley and Hettmansperger's (1993) robust regression estimators
#' using minimum volume ellipsoid (MVE) for detecting leverage points and
#' least trimmed squares (LTS) for initial estimates, with iterative refinement
#' using Huber's Psi function.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param bend Bending constant for Huber's Psi function (default: 1.345).
#' @param SEED Logical; if TRUE (default), sets random seed for reproducibility.
#' @inheritParams opreg
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{residuals}{Residuals from the final fit.}
#' }
#'
#' @details
#' The Coakley-Hettmansperger method provides robust regression estimates by:
#' 1. Using MVE (minimum volume ellipsoid) to compute Mallows weights for leverage points
#' 2. Computing initial estimates via least trimmed squares (LTS)
#' 3. Computing standardized residuals with robust scale estimate
#' 4. Applying Huber's Psi function with the specified bending constant
#' 5. Updating coefficient estimates using weighted least squares
#'
#' **Important**: When using \code{chreg} with bootstrap methods, use \code{chregF}
#' instead to avoid the warning message that appears by default.
#'
#' The method downweights both leverage points (via Mallows weights) and residual
#' outliers (via Huber's Psi), providing excellent resistance to various types of
#' contamination.
#'
#' @references
#' Coakley, C.W. and Hettmansperger, T.P. (1993). A bounded influence, high breakdown,
#' efficient regression estimator. Journal of the American Statistical Association,
#' 88, 872-880.
#'
#' @family robust regression functions
#' @seealso \code{\link{chregF}} for bootstrap-compatible version
#' @export
#' @examples
#' \dontrun{
#' # Simple robust regression
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2*x + rnorm(50, sd=0.5)
#' chreg(x, y)
#'
#' # Multiple regression with leverage points
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50, sd=0.3)
#' x[1:3,] <- x[1:3,] + 5  # add leverage points
#' chreg(x, y, xout=TRUE)  # remove detected outliers
#' }
chreg<-function(x,y,bend=1.345,SEED=TRUE,xout=FALSE,outfun=outpro,pr=TRUE,...){
#
# Compute Coakley Hettmansperger robust regression estimators
# JASA, 1993, 88, 872-880
#
# x is a n by p matrix containing the predictor values.
#
# No missing values are allowed
#
#  Comments in this function follow the notation used
#  by Coakley and Hettmansperger
#
# with old version of R, need library(lqs) when using ltsreg
# as the initial estimate.
#
if(pr)print('If using chreg with a bootstrap method, use chregF instead')
if(SEED)set.seed(12) # Set seed so that results are always duplicated.
x<-as.matrix(x)
p<-ncol(x)
m<-elimna(cbind(x,y))
x<-m[,1:p]
p1<-p+1
y<-m[,p1]
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
x<-as.matrix(x)
cutoff<-bend
mve<-vector("list")
if(ncol(x)==1){
mve$center<-median(x)
mve$cov<-mad(x)^2
}
if(ncol(x)>=2)mve<-cov.mve(x)  # compute minimum volume ellipsoid measures of
                 # location and scale and store in mve.
reg0<-ltsReg(x,y) # compute initial regression est using least trimmed
                 # squares.
# Next, compute the rob-md2(i) values and store in rob
rob<-1  # Initialize vector rob
mx<-mve$center
rob<-mahalanobis(x,mx,mve$cov)
k21<-qchisq(.95,p)
c62<-k21/rob
vecone<-c(rep(1,length(y))) # Initialize vector vecone to 1
c30<-pmin(vecone,c62)  # mallows weights put in c30
k81<-median(abs(reg0$residuals)) # median of absolute residuals
k72<-1.4826*(1+(5/(length(y)-p-1)))*k81 # lms scale
c60<-reg0$residuals/(k72*c30) # standardized residuals
#  compute psi and store in c27
cvec<-c(rep(cutoff,length(y))) # Initialize vector cvec to cutoff
c27<-pmin(cvec,c60)
c27<-pmax(-1*cutoff,c27)  #c27 contains psi values
#
# compute B matrix and put in c66.
#  Also, transform B so that i th diag elem = 0 if c27[i] is
#  between -cutoff and cutoff, 1 otherwise.
#
c66<-ifelse(abs(c27)<=bend,1,0) # Have derivative of psi in c66
m1<-cbind(1,x)  # X matrix with col of 1's added
m2<-t(m1)   #X transpose
m5<-diag(c30) # matrix W, diagonal contains weights
m4<-diag(c66) # B matrix
m6<-m4%*%m1   # BX
m7<-m2%*%m6   # X'BX (nD=X'BX)
m8<-solve(m7)  #m8 = (X'-B-X)inverse
m9<-m8%*%m2 #m9=X prime-B-X inverse X'
m9<-m9%*%m5 # m9=X prime-B-X inverse X'W
m10<-m9%*%c27
c20<-m10*k72
c21<-reg0$coef+c20 #update initial estimate of parameters.
res<-y-m1%*%c21
list(coef=t(c21),residuals=res)
}

# ============================================================================
# bmreg
# ============================================================================

#' Bounded M-Regression with Schweppe Weights
#'
#' Computes a bounded M-regression estimator using Huber's Psi function with
#' Schweppe weights. The method iteratively reweights observations based on
#' standardized residuals and leverage, providing robustness to both outliers
#' and leverage points.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param iter Maximum number of iterations for convergence (default: 20).
#' @param bend Bending constant for Huber's Psi function. Default is
#'   \code{2*sqrt((ncol(x)+1)/nrow(x))}, an adaptive choice based on problem dimensions.
#' @inheritParams opreg
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{residuals}{Residuals from the final fit.}
#'   \item{w}{Final weights for each observation.}
#' }
#'
#' @details
#' The bounded M-regression with Schweppe weights:
#' 1. Starts with least squares initial estimates
#' 2. Computes Schweppe weights based on leverage (\code{sqrt(1-h)})
#' 3. Computes robust scale using median of middle residuals
#' 4. Applies Huber's Psi function to weighted standardized residuals
#' 5. Updates coefficients via weighted least squares
#' 6. Iterates until convergence (coefficient change < 0.0001)
#'
#' The Schweppe weighting scheme (\code{nu = sqrt(1-h)}) downweights high-leverage
#' points, while the Huber Psi function provides bounded influence for residual outliers.
#' Together, these provide protection against both types of contamination.
#'
#' The default bending constant is adaptive and increases with sample size, providing
#' more robust estimates for larger datasets.
#'
#' @family robust regression functions
#' @export
#' @examples
#' \dontrun{
#' # Simple bounded M-regression
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2*x + rnorm(50, sd=0.5)
#' bmreg(x, y)
#'
#' # Multiple regression with custom bending constant
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50, sd=0.3)
#' y[1:5] <- y[1:5] + 10  # add outliers
#' fit <- bmreg(x, y, bend=1.5)
#' print(fit$coef)
#' }
bmreg<-function(x,y,iter=20,bend=2*sqrt((ncol(x)+1)/nrow(x)),xout=FALSE,outfun=outpro,...){
# compute a bounded M regression using Huber Psi and Schweppe weights.
# The predictors are assumed to be stored in the n by p matrix x.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x<-as.matrix(x)
init<-lsfit(x,y)
resid<-init$residuals
x1<-cbind(1,x)
nu<-sqrt(1-hat(x1))
low<-ncol(x)+1
for(it in 1:iter){
ev<-sort(abs(resid))
scale<-median(ev[c(low:length(y))])/qnorm(.75)
rov<-(resid/scale)/nu
psi<-ifelse(abs(rov)<=bend,rov,bend*sign(rov))  # Huber Psi
wt<-nu*psi/(resid/scale)
new<-lsfit(x,y,wt)
if(max(abs(new$coef-init$coef))<.0001)break
init$coef<-new$coef
resid<-new$residuals
}
resid<-y-x1%*%new$coef
if(max(abs(new$coef-init$coef))>=.0001)
paste("failed to converge in",iter,"steps")
list(coef=new$coef,residuals=resid,w=wt)
}

# ============================================================================
# reglev
# ============================================================================

#' Detect Regression Leverage Points
#'
#' Identifies good and bad leverage points in regression using the Rousseeuw
#' and van Zomeren (1990) method. This diagnostic tool helps detect influential
#' observations and regression outliers.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param plotit Logical; if \code{TRUE} (default), creates a diagnostic plot
#'   of robust distances vs. standardized residuals.
#' @param SEED Logical; if \code{TRUE} (default), sets random seed for
#'   reproducibility.
#' @param DIS Logical; if \code{TRUE}, returns robust Mahalanobis distances
#'   in output (default: \code{FALSE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{levpoints}{Row numbers of leverage points (outliers in predictor space)}
#'   \item{regout}{Row numbers of regression outliers (large standardized residuals)}
#'   \item{bad.lev.points}{Row numbers of bad leverage points (both leverage and
#'     regression outliers - high influence points)}
#'   \item{keep}{Row numbers of observations to keep (excluding bad leverage points)}
#'   \item{dis}{Robust Mahalanobis distances (if \code{DIS=TRUE}, otherwise \code{NULL})}
#'   \item{stanres}{Standardized residuals from LMS regression}
#'   \item{crit}{Critical value for leverage (97.5th percentile of chi-square)}
#' }
#'
#' @details
#' The function uses least median of squares (LMS) regression to identify outliers:
#'
#' **Leverage Points**: Observations with robust Mahalanobis distance exceeding
#' \eqn{\sqrt{\chi^2_{0.975, p}}}, where p is the number of predictors. These are
#' outliers in the predictor space.
#'
#' **Regression Outliers**: Observations with standardized residuals exceeding 2.5
#' in absolute value. These don't fit the regression model well.
#'
#' **Bad Leverage Points**: Observations that are both leverage points AND
#' regression outliers. These are the most problematic as they have high influence
#' on the regression fit and don't follow the overall pattern.
#'
#' The standardized residuals use a robust scale estimator based on the median
#' absolute residual with a finite-sample correction factor.
#'
#' If a plot is requested, it shows:
#' - X-axis: Robust Mahalanobis distance
#' - Y-axis: Standardized residuals
#' - Vertical line at the leverage cutoff
#' - Horizontal lines at ±2.5 for regression outlier cutoffs
#'
#' @references
#' Rousseeuw, P.J. & van Zomeren, B.C. (1990). Unmasking multivariate outliers
#' and leverage points. \emph{Journal of the American Statistical Association},
#' 85, 633-639.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @family robust regression diagnostics
#' @seealso \code{\link{lmsreg}} for LMS regression, \code{\link{reg.reglev}}
#'   for interactive version
#' @export
#' @examples
#' \dontrun{
#' # Simple example
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50, sd=0.5)
#'
#' # Add a bad leverage point
#' x[1,] <- c(5, 5)  # leverage point
#' y[1] <- 20        # also regression outlier
#'
#' # Detect leverage points
#' result <- reglev(x, y)
#' result$bad.lev.points  # Should identify row 1
#'
#' # Without plot
#' reglev(x, y, plotit=FALSE)
#' }
reglev<-function(x,y,plotit=TRUE,SEED=TRUE,DIS=FALSE){
#
#  Search for good and bad leverage points using the
#  Rousseuw and van Zomeren method.
#
#  x is an n by p matrix
#
#  The function returns the number of the rows in x that are identified
#  as outliers. (The row numbers are stored in outliers.)
#  It also returns the distance of the points identified as outliers
#  in the variable dis.
#
xy=elimna(cbind(x,y))
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
plotit<-as.logical(plotit)
if(SEED)set.seed(12)
x<-as.matrix(x)
res<-lmsreg(x,y)$resid
sighat<-sqrt(median(res^2))
sighat<-1.4826*(1+(5/(length(y)-ncol(x)-1)))*sighat
stanres<-res/sighat
if(ncol(x)>=2)mve<-cov.mve(x)
if(ncol(x)==1){
mve<-vector("list")
mve$center<-median(x)
mve$cov<-mad(x)^2
}
dis<-mahalanobis(x,mve$center,mve$cov)
dis<-sqrt(dis)
crit<-sqrt(qchisq(.975,ncol(x)))
chk<-ifelse(dis>crit,1,0)
vec<-c(1:nrow(x))
id<-vec[chk==1]
chkreg<-ifelse(abs(stanres)>2.5,1,0)
idreg<-vec[chkreg==1]
if(plotit){
plot(dis,stanres,xlab="Robust distances",ylab="standardized residuals")
abline(-2.5,0)
abline(2.5,0)
abline(v=crit)
}
all=c(id,idreg)
ID=duplicated(all)
blp=all[ID]
vec=c(1:length(y))
nkeep=vec
if(length(blp)>0)nkeep=vec[-blp]
if(!DIS)dis=NULL
list(levpoints=id,regout=idreg,bad.lev.points=blp,keep=nkeep,dis=dis,stanres=stanres,crit=crit)
}

# ============================================================================
# winreg
# ============================================================================

#' Winsorized Regression
#'
#' Computes a robust regression estimator based on Winsorized covariances.
#' The method uses Winsorized means and Winsorized covariances to compute
#' regression coefficients, then iteratively refines them.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param iter Maximum number of iterations for convergence (default: 20).
#' @inheritParams opreg
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{resid}{Residuals from the final fit.}
#' }
#'
#' @details
#' Winsorized regression provides robustness by:
#' 1. Computing Winsorized means for each predictor
#' 2. Computing Winsorized covariances between predictors and outcome
#' 3. Solving for initial slopes using these robust covariances
#' 4. Computing intercept from Winsorized mean of outcome
#' 5. Iteratively refining estimates by Winsorizing residuals
#' 6. Converges when coefficient changes are less than 0.0001
#'
#' Winsorization replaces extreme values (beyond the \code{tr} quantile)
#' with less extreme values, providing bounded influence while retaining
#' more information than trimming would.
#'
#' The method is particularly effective when you want robustness without
#' completely discarding outlying observations.
#'
#' @family robust regression functions
#' @seealso \code{\link{win}}, \code{\link{wincor}}
#' @export
#' @examples
#' \dontrun{
#' # Simple Winsorized regression
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2*x + rnorm(50, sd=0.5)
#' winreg(x, y)
#'
#' # Multiple regression with 10% Winsorization
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50, sd=0.3)
#' y[1:5] <- y[1:5] + 8  # add outliers
#' winreg(x, y, tr=0.1)
#' }
winreg<-function(x,y,iter=20,tr=.2,xout=FALSE,outfun=outpro,...){
#
# Compute a Winsorized regression estimator
# The predictors are assumed to be stored in the n by p matrix x.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x=as.matrix(x)
ma<-matrix(0,ncol(x),1)
m<-matrix(0,ncol(x),ncol(x))
mvals<-apply(x,2,win,tr)
for (i in 1:ncol(x)){
ma[i,1]<-wincor(x[,i],y,tr=tr)$cov
for (j in 1:ncol(x))m[i,j]<-wincor(x[,i],x[,j],tr=tr)$cov
}
slope<-solve(m,ma)
b0<-win(y,tr)-sum(slope%*%mvals)
for(it in 1:iter){
res<-y-x%*%slope-b0
for (i in 1:ncol(x))ma[i,1]<-wincor(x[,i],res,tr=tr)$cov
slopeadd<-solve(m,ma)
b0add<-win(res,tr)-sum(slopeadd%*%mvals)
if(max(abs(slopeadd),abs(b0add)) <.0001)break
slope<-slope+slopeadd
b0<-b0+b0add
}
if(max(abs(slopeadd),abs(b0add)) >=.0001)
paste("failed to converge in",iter,"iterations")
list(coef=c(b0,slope),resid=res)
}

# ============================================================================
# regpres1
# ============================================================================
regpres1<-function(isub,x,y,regfun,mval){
#
#  Perform regression using x[isub] to predict y[isub]
#  isub is a vector of length n,
#  a bootstrap sample from the sequence of integers
#  1, 2, 3, ..., n
#
#  This function is used by other functions when computing
#  bootstrap estimates.
#
#  regfun is some regression method already stored in R
#  It is assumed that regfun$coef contains the  intercept and slope
#  estimates produced by regfun.  The regression methods written for
#  this  book, plus regression functions in R, have this property.
#
#  x is assumed to be a matrix containing values of the predictors.
#
xmat<-matrix(x[isub,],mval,ncol(x))
regboot<-regfun(xmat,y[isub])
regboot<-regboot$coef
regboot
}

# ============================================================================
# hratio
# ============================================================================

#' Compute Half-Slope Ratios
#'
#' Computes a matrix of half-slope ratios for regression diagnostics. This
#' method divides the data based on each predictor and compares regression
#' slopes in the two halves, helping detect interactions and non-constant
#' regression relationships.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param regfun Regression function to use (default: \code{bmreg}). Must return
#'   coefficients in \code{regfun$coef} with intercept as first element.
#'
#' @return A p by p matrix of half-slope ratios where:
#' - Row i: half-slope ratios when data is split based on the i-th predictor
#' - Column j: ratio of slopes for the j-th predictor
#' - Element [i,j]: ratio of slope for predictor j in upper half vs. lower half
#'   when split by predictor i
#'
#' @details
#' For each predictor variable:
#' 1. Order observations by that predictor
#' 2. Split data into two halves (lower and upper)
#' 3. Fit regression model to each half
#' 4. Compute ratio of slopes: (upper half slope) / (lower half slope)
#'
#' **Interpretation**:
#' - Ratio near 1: slope is similar in both halves (good)
#' - Ratio far from 1: slope differs between halves, suggesting:
#'   - Non-constant relationship
#'   - Interaction effects
#'   - Possible need for transformation
#'
#' The diagonal elements [i,i] show whether predictor i's relationship with y
#' changes across the range of predictor i. Off-diagonal elements [i,j] show
#' whether predictor j's relationship with y changes across the range of
#' predictor i.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @family robust regression diagnostics
#' @seealso \code{\link{bmreg}}, \code{\link{reglev}}
#' @export
#' @examples
#' \dontrun{
#' # Example with constant relationship
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50, sd=0.5)
#' hratio(x, y)  # Should have ratios near 1
#'
#' # Example with interaction
#' y2 <- x[,1] + 2*x[,2] + 3*x[,1]*x[,2] + rnorm(50, sd=0.5)
#' hratio(x, y2)  # Ratios will deviate from 1
#' }
hratio<-function(x,y,regfun=bmreg){
#
#   Compute a p by p matrix of half-slope ratios
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
#  OUTPUT:
#The first row reports the half-slope
#ratios when the data are divided into two groups using the first predictor.
#The first column is the half-slope ratio for the first predictor, the
#second column is the half-slope ratio for the second predictor, and so forth.
#The second row contains the half-slope ratios when the data are divided
#into two groups using the second predictor, and so on.
#
x<-as.matrix(x)
xmat<-matrix(0,nrow(x),ncol(x))
mval<-floor(length(y)/2)
mr<-length(y)-mval
xmatl<-matrix(0,mval,ncol(x))
xmatr<-matrix(0,mr,ncol(x))
hmat<-matrix(NA,ncol(x),ncol(x))
isub<-c(1:length(y))
ksub<-c(1:ncol(x))+1
for (k in 1:ncol(x)){
xord<-order(x[,k])
yord<-y[xord]
yl<-yord[isub<=mval]
yr<-yord[isub>mval]
for (j in 1:ncol(x)){
xmat[,j]<-x[xord,j]
xmatl[,j]<-xmat[isub<=mval,j]
xmatr[,j]<-xmat[isub>mval,j]
}
coefl<-regfun(xmatl,yl)$coef
coefr<-regfun(xmatr,yr)$coef
hmat[k,]<-coefr[ksub[ksub>=2]]/coefl[ksub[ksub>=2]]
}
hmat
}

# ============================================================================
# mbmreg
# ============================================================================

#' Modified Bounded M-Regression
#'
#' Computes a modified bounded M-regression estimator using Huber's psi function
#' with Schweppe weights. Regression outliers (identified via LMS) are given
#' zero weight, providing strong protection against outliers.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param iter Maximum number of iterations for convergence (default: 20).
#' @param bend Bending constant for Huber's psi (default: \code{2*sqrt((p+1)/n)}).
#' @inheritParams opreg
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{residuals}{Residuals from the final fit.}
#'   \item{w}{Final weights assigned to each observation.}
#' }
#'
#' @details
#' This is the modified M-regression estimator described in Wilcox (2022), Chapter 8.
#' The procedure:
#'
#' 1. **Initial outlier detection**: Uses LMS regression to identify regression
#'    outliers (standardized residuals > 2.5). These get weight 0 throughout.
#'
#' 2. **Iterative reweighting**:
#'    - Computes Schweppe weights: \eqn{w_i = \nu_i \psi(r_i / \sigma \nu_i)}
#'    - Where \eqn{\nu_i = \sqrt{1 - h_{ii}}} (leverages from hat matrix)
#'    - \eqn{\psi} is Huber's psi function with tuning constant \code{bend}
#'    - Scale \eqn{\sigma} estimated robustly from residuals
#'
#' 3. **Convergence**: Iterates until coefficient changes < 0.0001
#'
#' The Schweppe weighting down-weights high leverage points more than low
#' leverage points, providing better finite-sample performance than standard
#' M-estimation. Regression outliers identified initially maintain zero weight,
#' protecting against severe outliers.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press. Chapter 8.
#'
#' @family robust regression functions
#' @seealso \code{\link{bmreg}} for standard bounded M-regression,
#'   \code{\link{lmsreg}} for least median of squares
#' @export
#' @examples
#' \dontrun{
#' # Simple example
#' set.seed(123)
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2*x + rnorm(50, sd=0.5)
#' y[1:3] <- y[1:3] + 10  # add outliers
#' mbmreg(x, y)
#'
#' # Multiple regression
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50)
#' mbmreg(x, y)
#' }
mbmreg<-function(x,y,iter=20,bend=2*sqrt(ncol(x)+1)/nrow(x),xout=FALSE,outfun=outpro,...){
#
# Compute a bounded M regression estimator using
# Huber Psi and Schweppe weights with
# regression outliers getting a weight of zero.
#
# This is the modified M-regression estimator in Chapter 8
#
# The predictors are assumed to be stored in the n by p matrix x.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x<-as.matrix(x)
if(is.matrix(y)){
if(ncol(y)==1)y=as.vector(y)
}
x1<-cbind(1,x)
reslms<-lmsreg(x,y)$resid
sighat<-sqrt(median(reslms^2))
sighat<-1.4826*(1+(5/(length(y)-ncol(x)-1)))*sighat
if(sighat==0)warning("The estimated measure of scale, based on the residuals using lms regression, is zero")
temp<-ifelse(sighat*reslms>0,abs(reslms)/sighat,0*reslms)
wt<-ifelse(temp<=2.5,1,0)
init<-lsfit(x,y,wt)
resid<-init$residuals
nu<-sqrt(1-hat(x1))
low<-ncol(x)+1
for(it in 1:iter){
ev<-sort(abs(resid))
scale<-median(ev[c(low:length(y))])/qnorm(.75)
rov<-(resid/scale)/nu
psi<-ifelse(abs(rov)<=bend,rov,bend*sign(rov))  # Huber Psi
wt<-nu*psi/(resid/scale)
wt<-ifelse(temp<=2.5,wt,0)
new<-lsfit(x,y,wt)
if(abs(max(new$coef-init$coef)<.0001))break
init$coef<-new$coef
resid<-new$residuals
}
resid<-y-x1%*%new$coef
if(abs(max(new$coef-init$coef)>=.0001))
paste("failed to converge in",iter,"steps")
list(coef=new$coef,residuals=resid,w=wt)
}

# ============================================================================
# regts1
# ============================================================================
regts1<-function(vstar,yhat,res,mflag,x,tr){
ystar<-yhat+res*vstar
bres<-ystar-mean(ystar,tr)
rval<-0
for (i in 1:nrow(x)){
rval[i]<-sum(bres[mflag[,i]])
}
rval
}

# ============================================================================
# depreg
# ============================================================================

#' Depth-Based Regression (Univariate)
#'
#' Computes the depth regression estimator for simple linear regression (one
#' predictor). The method finds the regression line with maximum regression depth,
#' providing a robust alternative to least squares.
#'
#' @param x A numeric vector of the predictor variable (or single-column matrix).
#' @param y A numeric vector of the dependent variable (length n).
#' @inheritParams opreg
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients: \code{c(intercept, slope)}.}
#'   \item{residuals}{Residuals from the fit.}
#' }
#'
#' @details
#' The depth regression estimator finds the line with maximum regression depth.
#' For each pair of points (i,j):
#' 1. Compute the line passing through both points
#' 2. Calculate the regression depth of this line
#' 3. Select the line(s) with maximum depth
#' 4. Average coefficients if multiple lines achieve maximum depth
#'
#' **Regression depth** of a line measures how "central" it is relative to the
#' data cloud. A line with high depth is surrounded by many data points,
#' making it resistant to outliers.
#'
#' This method provides:
#' - High breakdown point (up to 33%)
#' - Robustness to outliers in both x and y
#' - Equivariance under affine transformations
#'
#' **Note**: This function only supports simple linear regression (one predictor).
#' For multiple predictors, use \code{\link{mdepreg}}.
#'
#' @references
#' Rousseeuw, P.J. & Hubert, M. (1999). Regression depth. \emph{Journal of the
#' American Statistical Association}, 94, 388-433.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @family robust regression functions
#' @seealso \code{\link{mdepreg}} for multiple regression, \code{\link{tsreg}}
#'   for Theil-Sen regression
#' @export
#' @examples
#' \dontrun{
#' # Simple linear regression
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + 1 + rnorm(50, sd=0.5)
#' depreg(x, y)
#'
#' # With outliers
#' y[1:3] <- y[1:3] + 10
#' depreg(x, y)  # Robust to outliers
#'
#' # Compare with least squares
#' lm(y ~ x)$coef  # Affected by outliers
#' }
depreg<-function(x,y,xout=FALSE,outfun=out,...){
#
# Compute the depth regression estimator.
# Only a single predictor is allowed in this version
#  Perhaps use instead
#
if(is.matrix(x)){
if(ncol(x)>=2)stop("Only a single predicor is allowed")
x<-as.vector(x)
}
xy=cbind(x,y)
xy=elimna(xy)
if(xout){
flag<-outfun(xy[,1],plotit=FALSE,...)$keep
xy<-xy[flag,]
}
x=xy[,1]
y=xy[,2]
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
slope<-v1/v2
vec3<-outer(ys,ys,"+")
vec4<-outer(xs,xs,"+")
v3<-vec3[vec2>0]
v4<-vec4[vec2>0]
deep<-NA
inter<-v3/2-slope*v4/2
temp<-matrix(c(inter,slope),ncol=2)
deep<-apply(temp,1,rdepth.orig,x,y)
best<-max(deep)
coef<-NA
coef[2]<-mean(slope[deep==best])
coef[1]<-mean(inter[deep==best])
res<-y-coef[2]*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# tsgreg
# ============================================================================

#' Theil-Sen Regression for Groups
#'
#' Computes the Theil-Sen regression estimator using all pairwise slopes,
#' which is particularly useful when data are divided into groups or when
#' you want the most robust version of Theil-Sen estimation.
#'
#' @param x A numeric vector of the predictor variable.
#' @param y A numeric vector of the dependent variable (same length as x).
#' @param tries Number of slope pairs to compute (default: all possible pairs
#'   = \code{n*(n-1)/2}).
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients: \code{c(intercept, slope)}.}
#'   \item{residuals}{Residuals from the fit.}
#' }
#'
#' @details
#' This function computes the Theil-Sen estimator by:
#' 1. Computing slopes for all (or \code{tries}) pairs of distinct points
#' 2. Taking the median of these slopes as the robust slope estimate
#' 3. Computing intercept as median of \eqn{y_i - \hat{\beta} x_i}
#'
#' The Theil-Sen estimator:
#' - Has a breakdown point of approximately 29%
#' - Is highly resistant to outliers
#' - Provides good efficiency for normal data
#' - Is distribution-free
#'
#' **Note**: This version uses simple median of all pairwise slopes. For
#' handling tied x-values or when speed is important, consider
#' \code{\link{tsreg}} or \code{\link{tshdreg}}.
#'
#' @references
#' Sen, P.K. (1968). Estimates of the regression coefficient based on Kendall's tau.
#' \emph{Journal of the American Statistical Association}, 63, 1379-1389.
#'
#' Theil, H. (1950). A rank-invariant method of linear and polynomial regression
#' analysis. \emph{Indagationes Mathematicae}, 12, 85-91.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @family robust regression functions
#' @seealso \code{\link{tsreg}} for faster Theil-Sen with iteration,
#'   \code{\link{tshdreg}} for Harrell-Davis version
#' @export
#' @examples
#' \dontrun{
#' # Simple example
#' set.seed(123)
#' x <- 1:50
#' y <- 2*x + rnorm(50, sd=2)
#' tsgreg(x, y)
#'
#' # With outliers
#' y[1:3] <- y[1:3] + 20
#' tsgreg(x, y)  # Robust result
#' }
tsgreg<-function(x,y,tries=(length(y)^2-length(y))/2){
#
#
x<-as.matrix(x)
if(nrow(x)!=length(y))stop("Length of y must match the number of rows of x")
# eliminate any rows with missing values.
m1<-cbind(x,y)
m1<-elimna(m1)
x<-m1[,1:ncol(x)]
y<-m1[,ncol(x)+1]
set.seed(2)
data<-matrix(NA,ncol=ncol(x)+1,nrow=tries)
for(i in 1:tries){
data[i,]<-sample(length(y),size=ncol(x)+1,replace=FALSE)
}
bvec <- apply(data, 1,tsgregs1,x,y)
coef<-0
numzero<-0
loc<-0
for (i in 1:ncol(x)){
ip<-i+1
temp<-bvec[ip,]
loc[i]<-median(x[,i])
coef[i+1]<-median(temp[temp!=0])
numzero[i]<-length(temp[temp==0])
}
ip<-ncol(x)+1
coef[1]<-median(y)-sum(coef[2:ip]*loc)
res<-y-x %*% coef[2:ip] - coef[1]
list(coef=coef,residuals=res,numzero=numzero)
}

# ============================================================================
# tsgregs1
# ============================================================================
tsgregs1<-function(isub,x,y){
#
#  This function is used by tsgreg
#
#  Perform regression using x[isub,] to predict y[isub]
#  isub is a vector of length nsub, determined by tsgreg
#
tsgregs1<-lsfit(x[isub,],y[isub])$coef
}

# ============================================================================
# lts1reg
# ============================================================================
lts1reg<-function(x,y,tr=.2,h=NA){
#
# Compute the least trimmed squares regression estimator.
# Only a single predictor is allowed in this version
#
if(is.na(h))h<-length(x)-floor(tr * length(x))
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
slope<-v1/v2
vec3<-outer(ys,ys,"+")
vec4<-outer(xs,xs,"+")
v3<-vec3[vec2>0]
v4<-vec4[vec2>0]
val<-NA
inter<-v3/2-slope*v4/2
for(i in 1:length(slope)){
#risk<-(y[vec2>0]-slope[i]*x[vec2>0]-inter[i])^2
risk<-(y-slope[i]*x-inter[i])^2
risk<-sort(risk)
val[i]<-sum(risk[1:h])
}
best<-min(val)
coef<-NA
coef[2]<-mean(slope[val==best])
coef[1]<-mean(inter[val==best])
res<-y-coef[2]*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# twolsreg
# ============================================================================
twolsregsub<-function(isub, x, y)
{
        #
        #  Compute least squares estimate of the
        #  slope using x[isub] and y[isub]
        #  isub is a vector of length n,
        #  a bootstrap sample from the sequence of integers
        #  1, 2, 3, ..., n
        #
        twolsregsub<-lsfit(x[isub],y[isub])$coef[2]
        twolsregsub
}

# ============================================================================
# regi
# ============================================================================

#' Conditional Regression Comparison Based on Third Variable
#'
#' Splits data into two groups based on whether a third variable \code{z} is
#' less than or greater than a specified point, then plots smoothed regression
#' lines for each group. Optionally tests whether the regression coefficients
#' differ between groups.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values (same length as \code{x}).
#' @param z Numeric vector of conditioning variable (same length as \code{x}).
#' @param pt Split point for \code{z} (default: \code{median(z)}).
#' @param fr Span parameter for running mean smoother (default: 0.8).
#' @param est Location estimator for smoothing (default: \code{onestep}).
#' @param regfun Regression function for testing (default: \code{tsreg}).
#'   Only used if \code{testit=TRUE}.
#' @param testit Logical; if TRUE, performs formal test of regression
#'   equality via \code{\link{reg2ci}} and plots fitted lines (default: FALSE).
#' @param ... Additional arguments passed to \code{runmean2g} or \code{regfun}.
#'
#' @return If \code{testit=FALSE}, returns "Done" after plotting. If
#'   \code{testit=TRUE}, returns the output from \code{\link{reg2ci}}, which
#'   includes confidence intervals for the difference in regression
#'   coefficients.
#'
#' @details
#' This function is useful for exploring whether the relationship between
#' \code{x} and \code{y} changes depending on the value of a third variable
#' \code{z}. The data are split at \code{pt} (default: median of \code{z}),
#' and smoothed regression curves are plotted for each group.
#'
#' When \code{testit=TRUE}, the function additionally:
#' \itemize{
#'   \item Fits regression lines using \code{regfun} to each group
#'   \item Overlays these fitted lines on the plot (solid for z < pt, dashed for z >= pt)
#'   \item Computes bootstrap confidence intervals for coefficient differences via \code{\link{reg2ci}}
#' }
#'
#' The smoothing uses \code{\link{runmean2g}}, which applies a running mean
#' smoother separately to each group and plots both on the same axes.
#'
#' @seealso \code{\link{reg2ci}} for comparing regression coefficients,
#'   \code{\link{runmean2g}} for two-group running means
#'
#' @examples
#' \dontrun{
#' # Example: Does the x-y relationship differ for low vs high z?
#' set.seed(123)
#' x <- rnorm(100)
#' z <- rnorm(100)
#' # Create different slopes for low vs high z
#' y <- ifelse(z < median(z), 2*x + rnorm(100), 0.5*x + rnorm(100))
#'
#' # Visual comparison with smooths
#' regi(x, y, z)
#'
#' # Formal test with fitted lines
#' result <- regi(x, y, z, testit=TRUE)
#' print(result)
#' }
#'
#' @export
regi<-function(x,y,z,pt=median(z),fr=.8,est=onestep,regfun=tsreg,testit=FALSE,...){
#
# split the data according to whether z is < or > pt, then
# use runmean2g to plot a smooth of the regression
# lines corresponding to these two groups.
#
m<-cbind(x,y,z)
m<-elimna(m)
x<-m[,1]
y<-m[,2]
z<-m[,3]
flag<-(z<pt)
runmean2g(x[flag],y[flag],x[!flag],y[!flag],fr=fr,est=est,...)
output<-"Done"
if(testit){
abline(regfun(x[flag],y[flag])$coef)
abline(regfun(x[!flag],y[!flag])$coef,lty=2)
output<-reg2ci(x[flag],y[flag],x[!flag],y[!flag],regfun=regfun,plotit=FALSE)
}
output
}

# ============================================================================
# linchk
# ============================================================================

#' Test Linearity by Splitting Data on Predictor Value
#'
#' Splits the data into two groups based on whether a specified predictor
#' variable is less than or equal to a split point, then tests whether the
#' regression coefficients (particularly slopes) are equal across the two groups.
#'
#' @param x Numeric matrix of predictor variables (n by p).
#' @param y Numeric vector of response values (length n).
#' @param sp Split point value for the predictor specified by \code{pv}.
#' @param pv Column number of the predictor to use for splitting (default: 1).
#' @param regfun Robust regression function to use (default: \code{tsreg}).
#' @inheritParams common-params
#' @param pr Logical; if TRUE, prints which predictor is being used for
#'   splitting (default: TRUE).
#'
#' @return A list (output from \code{\link{reg2ci}}) containing:
#' \describe{
#'   \item{ci.dif}{Confidence intervals for differences in coefficients}
#'   \item{p.value}{P-value for test of equal coefficients}
#'   \item{est.1}{Coefficient estimates for group where predictor <= sp}
#'   \item{est.2}{Coefficient estimates for group where predictor > sp}
#' }
#'
#' @details
#' This function tests for nonlinearity or interaction effects by checking
#' whether the regression relationship changes when a predictor crosses a
#' threshold value. The data are split into two groups:
#' \itemize{
#'   \item Group 1: Observations where \code{x[, pv] <= sp}
#'   \item Group 2: Observations where \code{x[, pv] > sp}
#' }
#'
#' The function then calls \code{\link{reg2ci}} to:
#' 1. Estimate regression coefficients separately for each group using \code{regfun}
#' 2. Compute bootstrap confidence intervals for the difference in coefficients
#' 3. Test the null hypothesis that coefficients are equal across groups
#'
#' If the confidence intervals for slope differences exclude zero, this suggests
#' that the relationship between predictors and response changes across the
#' split point, indicating nonlinearity or an interaction with the splitting
#' variable.
#'
#' @seealso \code{\link{reg2ci}} for two-group regression comparison,
#'   \code{\link{lintest}} for general linearity testing,
#'   \code{\link{regi}} for conditional regression with smoothing
#'
#' @examples
#' \dontrun{
#' # Generate data with different slopes above/below median
#' set.seed(123)
#' x <- matrix(rnorm(200), ncol=2)
#' # Slope changes when x[,1] > 0
#' y <- ifelse(x[,1] <= 0, 2*x[,1] + x[,2], 0.5*x[,1] + x[,2]) + rnorm(100)
#'
#' # Test if slopes differ when splitting x[,1] at 0
#' result <- linchk(x, y, sp=0, pv=1)
#' print(result$p.value)
#'
#' # Try different split point
#' linchk(x, y, sp=median(x[,1]), pv=1)
#' }
#'
#' @export
linchk<-function(x,y,sp,pv=1,regfun=tsreg,plotit=TRUE,nboot=599,alpha=.05,pr=TRUE,xout=FALSE){
#
# Split the data into two groups according to whether
# predictor variable pv has a value less than sp.
# Then test the hypothesis that slope coefficients,
# based on the regression method regfun, are equal.
#
x<-as.matrix(x)
if(pr)print(paste("Splitting data using predictor", pv))
xx<-x[,pv]
flag<-(xx<=sp)
temp<-reg2ci(x[flag,],y[flag],x[!flag,],y[!flag],regfun=regfun,plotit=plotit,nboot=nboot,alpha=alpha,xout=xout)
temp
}

# ============================================================================
# lintests1
# ============================================================================

#' Internal Helper for Linearity Test Bootstrap
#'
#' Computes R-values for a single bootstrap sample in the linearity test.
#' This function is called by \code{\link{lintestMC}} and \code{\link{lintest}}.
#'
#' @param vstar Vector of standardized random values for bootstrap.
#' @param yhat Fitted values from the regression.
#' @param res Residuals from the regression.
#' @param mflag Matrix indicating ordering relationships between observations.
#' @param x Matrix of predictor variables.
#' @param regfun Regression function to use.
#' @param ... Additional arguments passed to \code{regfun}.
#'
#' @return Vector of R-values for the bootstrap sample.
#'
#' @details
#' This is an internal function used in the bootstrap implementation of the
#' Stute et al. (1998) linearity test. It should not be called directly by users.
#'
#' @keywords internal
#' @noRd
lintests1<-function(vstar,yhat,res,mflag,x,regfun,...){
ystar<-yhat+res*vstar
bres<-regfun(x,ystar,...)$residuals
rval<-0
for (i in 1:nrow(x)){
rval[i]<-sum(bres[mflag[,i]])
}
rval
}

# ============================================================================
# lsfitci
# ============================================================================

#' Bootstrap Confidence Intervals for Least Squares Regression
#'
#' Computes bootstrap confidence intervals for intercept and slopes in ordinary
#' least squares regression. For a single predictor, uses an adjusted percentile
#' method that handles heteroscedasticity well.
#'
#' @inheritParams common_params
#' @param nboot Number of bootstrap samples (default: 599).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param ... Additional arguments passed to \code{outfun}.
#'
#' @return A list with components:
#' \item{intercept.ci}{Vector of length 2 giving the confidence interval for the intercept.}
#' \item{slope.ci}{Matrix with one row per predictor. Columns are lower and upper CI bounds.}
#' \item{crit.level}{Bonferroni-adjusted significance level when p > 1 (NULL when p=1).}
#' \item{p.values}{Matrix of p-values for each slope when p > 1 (NA when p=1).}
#'
#' @details
#' This function is designed specifically for ordinary least squares regression
#' and provides better performance than \code{regci} with \code{regfun=lsfit} when n < 250.
#'
#' For a single predictor (\code{p=1}), the function uses an adjusted percentile
#' bootstrap method that provides accurate coverage even under heteroscedasticity.
#' The adjustment varies with sample size to maintain nominal coverage.
#'
#' For multiple predictors (\code{p>1}), uses standard percentile bootstrap with
#' Bonferroni correction to control family-wise error rate.
#'
#' \strong{Important}: The adjusted method for p=1 has been calibrated only for
#' \code{alpha=0.05}. If a different alpha is specified with p=1, it will be
#' reset to 0.05 with a warning.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' @seealso \code{\link{regci}}, \code{\link{lsfitNci}}, \code{\link{rregci}}
#'
#' @examples
#' # Single predictor with heteroscedastic errors
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2 + 3*x + rnorm(100, sd=abs(x))  # Heteroscedastic errors
#' lsfitci(x, y, nboot=500)
#'
#' # Multiple predictors
#' x <- matrix(rnorm(200), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(100)
#' lsfitci(x, y, nboot=500)
#'
#' @export
lsfitci<-function(x,y,nboot=599,alpha=.05,SEED=TRUE,xout=FALSE,outfun=out){
#
#   Compute a confidence interval for the slope parameters of
#   a linear regression equation when using the least squares estimator.
#
#   For p=1 predictor,
#   this function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
#   For p>1, a standard percentile bootstrap method is used
#   with FWE (the probability of at least one type I error)
#   controlled via the Bonferroni inequality.
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=599
#
#   SEED=T causes the seed of the random number generator to be set to 2,
#   otherwise the seed is not set.
#
#   Warning: probability coverage has been studied only when alpha=.05
#
x<-as.matrix(x)
p<-ncol(x)
pp<-p+1
temp<-elimna(cbind(x,y)) # Remove any missing values.
x<-temp[,1:p]
y<-temp[,p+1]
if(xout){
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,pp]
}
x<-as.matrix(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples; please wait")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,lsfit) # A p+1 by n matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
if(p==1){
if(alpha != .05){print("Resetting alpha to .05")
print("With p=1, unknown how to adjust confidence interval")
print("when alpha is not equal to .05.")
}
ilow<-15
ihi<-584
if(length(y) < 250){
ilow<-13
ihi<-586
}
if(length(y) < 180){
ilow<-10
ihi<-589
}
if(length(y) < 80){
ilow<-7
ihi<-592
}
if(length(y) < 40){
ilow<-6
ihi<-593
}
ilow<-round((ilow/599)*nboot)
ihi<-round((ihi/599)*nboot)
}
if(p>1){
ilow<-round(alpha*nboot/2)+1
ihi<-nboot-ilow
}
lsfitci<-matrix(0,ncol(x),2)
for(i in 1:ncol(x)){
ip<-i+1
bsort<-sort(bvec[ip,])
lsfitci[i,1]<-bsort[ilow+1]
lsfitci[i,2]<-bsort[ihi]
}
bsort<-sort(bvec[1,])
interceptci<-c(bsort[15],bsort[584])
crit.level<-NA
pmat<-NA
if(p>1){
crit.level<-alpha/p
pmat<-matrix(NA,nrow=p,ncol=2)
dimnames(pmat) <- list(NULL, c("Slope","p-value"))
for(pv in 1:p){
pmat[pv,1]<-pv
pp<-pv+1
pmat[pv,2]<-(sum(bvec[pp,]<0)+.5*sum(bvec[pp,]==0))/nboot
temp3<-1-pmat[pv,2]
pmat[pv,2]<-2*min(pmat[pv,2],temp3)
}}
list(intercept.ci=interceptci,slope.ci=lsfitci,crit.level=crit.level,
p.values=pmat)
}

# ============================================================================
# lsfitNci
# ============================================================================

#' Heteroscedasticity-Consistent Confidence Intervals for OLS Regression
#'
#' Computes confidence intervals for ordinary least squares (OLS) regression
#' coefficients using the HC3 heteroscedasticity-consistent standard errors
#' recommended by Long and Ervin (2000).
#'
#' @param x Numeric vector or matrix of predictor variable(s).
#' @param y Numeric vector of the response variable.
#' @param alpha Significance level for confidence intervals (default: 0.05 for 95% CIs).
#'
#' @return A list with components:
#' \describe{
#'   \item{ci}{Matrix with columns: parameter index, lower CI limit, upper CI limit}
#'   \item{stand.errors}{Vector of HC3 heteroscedasticity-consistent standard errors}
#' }
#'
#' @details
#' This function performs ordinary least squares regression but uses
#' heteroscedasticity-consistent (robust) standard errors for inference,
#' specifically the HC3 estimator of Long and Ervin (2000).
#'
#' **HC3 Standard Errors**:
#' The HC3 estimator adjusts for heteroscedasticity (non-constant error variance)
#' and leverage points. It modifies residuals as:
#' \deqn{r_i^* = r_i / (1 - h_i)^2}
#' where \eqn{r_i} are OLS residuals and \eqn{h_i} are leverage values (diagonal
#' elements of the hat matrix).
#'
#' **When to use**:
#' - Ordinary least squares regression (not robust regression)
#' - Suspect heteroscedasticity (non-constant variance)
#' - Want valid inference without assuming homoscedasticity
#' - Alternative to weighted least squares when weights unknown
#'
#' **Advantages over classical OLS inference**:
#' - Valid even when errors have non-constant variance
#' - Accounts for leverage points
#' - No need to specify variance structure
#' - Better small-sample properties than HC0, HC1, HC2
#'
#' The function prints confidence intervals and returns both intervals and
#' standard errors.
#'
#' @references
#' Long, J.S., & Ervin, L.H. (2000). Using heteroscedasticity consistent standard
#' errors in the linear regression model. \emph{The American Statistician}, 54, 217-224.
#'
#' MacKinnon, J.G., & White, H. (1985). Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics}, 29, 305-325.
#'
#' @seealso \code{\link{lsfitci}} for bootstrap-based OLS inference,
#'   \code{\link{regci}} for robust regression with bootstrap inference
#'
#' @examples
#' # OLS with heteroscedastic errors
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2 + 3*x + rnorm(100, sd=abs(x))  # Variance increases with |x|
#'
#' # HC3 confidence intervals
#' result <- lsfitNci(x, y)
#' print(result$ci)
#'
#' # Multiple predictors
#' x <- matrix(rnorm(200), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(100, sd=abs(x[,1]))
#' lsfitNci(x, y, alpha=0.01)  # 99% confidence intervals
#'
#' @export
lsfitNci<-function(x,y,alpha=.05){
#
# Compute confidence interval for least squares
# regression using heteroscedastic method
# recommended by Long and Ervin (2000).
#
x<-as.matrix(x)
if(nrow(x) != length(y))stop("Length of y does not match number of x values")
m<-cbind(x,y)
m<-elimna(m)
y<-m[,ncol(x)+1]
temp<-lsfit(x,y)
x<-cbind(rep(1,nrow(x)),m[,1:ncol(x)])
xtx<-solve(t(x)%*%x)
h<-diag(x%*%xtx%*%t(x))
hc3<-xtx%*%t(x)%*%diag(temp$res^2/(1-h)^2)%*%x%*%xtx
df<-nrow(x)-ncol(x)
crit<-qt(1-alpha/2,df)
al<-ncol(x)
ci<-matrix(NA,nrow=al,ncol=3)
for(j in 1:al){
ci[j,1]<-j
ci[j,2]<-temp$coef[j]-crit*sqrt(hc3[j,j])
ci[j,3]<-temp$coef[j]+crit*sqrt(hc3[j,j])
}
print("Confidence intervals for intercept followed by slopes:")
list(ci=ci,stand.errors=sqrt(diag(hc3)))
}

# ============================================================================
# taureg
# ============================================================================

#' Correlation-Based Regression using Kendall's Tau
#'
#' Computes correlations between a response variable and each predictor using
#' Kendall's tau (or another specified correlation function). Useful for screening
#' predictors and identifying potentially important variables for regression.
#'
#' @param m Numeric matrix of predictor variables (n by p).
#' @param y Numeric vector of the response variable (length n).
#' @param corfun Correlation function to use (default: \code{tau} for Kendall's tau).
#'   Must return a list with components \code{$cor} and \code{$p.value}.
#'   Alternatives: \code{pbcor}, \code{wincor}, \code{scor}.
#' @param ... Additional arguments passed to \code{corfun}.
#'
#' @return A list with components:
#' \describe{
#'   \item{cor}{Vector of correlations between y and each column of m}
#'   \item{p.value}{Vector of two-sided p-values for testing each correlation = 0}
#' }
#'
#' @details
#' This function is useful for exploring relationships between a response
#' variable and multiple predictors using robust correlation measures. It
#' computes the correlation between y and each column of m separately.
#'
#' **Kendall's Tau** (default):
#' - Robust to outliers and heavy-tailed distributions
#' - Based on concordance/discordance of pairs
#' - Range: [-1, 1]
#' - Interpretation similar to Pearson correlation but more robust
#'
#' **Alternative correlation functions**:
#' - \code{pbcor}: Percentage bend correlation (robust alternative to Pearson)
#' - \code{wincor}: Winsorized correlation
#' - \code{scor}: Spearman's rho (rank-based)
#' - Any function returning \code{$cor} and \code{$p.value}
#'
#' **Use cases**:
#' - Variable screening in high-dimensional regression
#' - Identifying important predictors before model building
#' - Exploratory data analysis
#' - When relationship may be monotonic but not linear
#'
#' **Note**: This is not a regression fit but rather a correlation analysis.
#' For actual regression using correlation-based methods, see \code{\link{correg}}
#' or \code{\link{scorreg}}.
#'
#' @references
#' Kendall, M.G. (1938). A new measure of rank correlation. \emph{Biometrika},
#' 30, 81-93.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{tau}} for Kendall's tau,
#'   \code{\link{pbcor}} for percentage bend correlation,
#'   \code{\link{correg}} for correlation-based regression fitting
#'
#' @examples
#' # Correlation analysis with multiple predictors
#' set.seed(123)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' x3 <- rnorm(100)
#' y <- 2*x1 - x2 + rnorm(100)
#' m <- cbind(x1, x2, x3)
#'
#' # Kendall's tau for each predictor
#' result <- taureg(m, y)
#' print(result$cor)
#' print(result$p.value)
#'
#' # Using percentage bend correlation
#' taureg(m, y, corfun=pbcor)
#'
#' @export
taureg<-function(m,y,corfun=tau,...){
#
#    Compute Kendall's tau between y and each of the
#    p variables stored  in the n by p matrix m.
#
#    Alternative measures of correlation can be used via the
#    argument corfun. The only requirement is that the function
#    corfun returns the correlation in corfun$cor and the p-value
#    in corfun$p.value.
#
#    This function also returns the two-sided significance level
#    for all pairs of variables, plus a test of zero correlations
#    among all pairs. (See chapter 9 of Wilcox, 2005, for details.)
#
m<-as.matrix(m)
tauvec<-NA
siglevel<-NA
for (i in 1:ncol(m)){
pbc<-corfun(m[,i],y,...)
tauvec[i]<-pbc$cor
siglevel[i]<-pbc$p.value
}
list(cor=tauvec,p.value=siglevel)
}

# ============================================================================
# lintest
# ============================================================================

#' Test Linearity of Regression Surface
#'
#' Tests the hypothesis that the regression surface is a plane (i.e., that the
#' true regression function is linear) using the method of Stute et al. (1998).
#'
#' @inheritParams common_params
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#'   Must return residuals in \code{$residuals}.
#' @param nboot Number of bootstrap samples (default: 500).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param ... Additional arguments passed to \code{regfun} and \code{outfun}.
#'
#' @return A list with components:
#' \item{dstat}{Kolmogorov-Smirnov-type test statistic (maximum absolute deviation).}
#' \item{wstat}{Cramér-von Mises-type test statistic (mean squared deviation).}
#' \item{p.value.d}{P-value for the KS-type statistic.}
#' \item{p.value.w}{P-value for the CvM-type statistic.}
#'
#' @details
#' This function implements the linearity test of Stute, Thies, and Zhu (1998),
#' which tests whether the conditional expectation E(Y|X) is a linear function of X.
#' Under the null hypothesis of linearity, the regression function is a plane.
#'
#' The test constructs two statistics:
#' \itemize{
#'   \item{\code{dstat}: }{Kolmogorov-Smirnov-type statistic based on maximum absolute deviation}
#'   \item{\code{wstat}: }{Cramér-von Mises-type statistic based on mean squared deviation}
#' }
#'
#' Both statistics compare the empirical process based on residuals to its expected
#' value under linearity. Bootstrap is used to determine critical values.
#'
#' @references
#' Stute, W., Thies, S., & Zhu, L.-X. (1998). Model checks for regression: An innovation
#' process approach. \emph{Journal of the American Statistical Association}, 93, 141-149.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' @seealso \code{\link{regtest}}, \code{\link{tsreg}}, \code{\link{linchk}}
#'
#' @examples
#' # Test linearity with Theil-Sen regression
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2*x[,1] - x[,2] + rnorm(50)  # Linear relationship
#' lintest(x, y, nboot=500)
#'
#' # Test with nonlinear relationship
#' y2 <- x[,1]^2 + rnorm(50)  # Nonlinear
#' lintest(x, y2, nboot=500)
#'
#' @export
lintest<-function(x,y,regfun=tsreg,nboot=500,alpha=.05,xout=FALSE,SEED=TRUE,
outfun=out,...){
#
# Test the hypothesis that the regression surface is a plane.
# Stute et al. (1998, JASA, 93, 141-149).
#
if(SEED)set.seed(2)
#if(identical(regfun,Qreg))print('When using Qreg, be sure to include res.vals=TRUE')
#if(identical(regfun,tshdreg))print('When using tshdreg, be sure to include RES=TRUE')
#if(identical(regfun,MMreg))print('When using MMreg, be sure to include RES=TRUE') # no longer necessary
x<-as.matrix(x)
d<-ncol(x)
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
x<-as.matrix(x)
y<-temp[,d+1]
if(xout){
flag<-outfun(x,...)$keep
x<-x[flag,]
x<-as.matrix(x)
y<-y[flag]
}
mflag<-matrix(NA,nrow=length(y),ncol=length(y))
for (j in 1:length(y)){
for (k in 1:length(y)){
mflag[j,k]<-(sum(x[j,]<=x[k,])==ncol(x))
}
}
reg<-regfun(x,y,...)
yhat<-y-reg$residuals
#print("Taking bootstrap samples, please wait.")
data<-matrix(runif(length(y)*nboot),nrow=nboot)
data<-sqrt(12)*(data-.5) # standardize the random numbers.
rvalb<-apply(data,1,lintests1,yhat,reg$residuals,mflag,x,regfun,...)
# An n x nboot matrix of R values
rvalb<-rvalb/sqrt(length(y))
dstatb<-apply(abs(rvalb),2,max)
wstatb<-apply(rvalb^2,2,mean)
# compute test statistic
v<-c(rep(1,length(y)))
rval<-lintests1(v,yhat,reg$residuals,mflag,x,regfun,...)
rval<-rval/sqrt(length(y))
dstat<-max(abs(rval))
wstat<-mean(rval^2)
ib<-round(nboot*(1-alpha))
p.value.d<-1-sum(dstat>=dstatb)/nboot
p.value.w<-1-sum(wstat>=wstatb)/nboot
list(dstat=dstat,wstat=wstat,p.value.d=p.value.d,p.value.w=p.value.w)
}

# ============================================================================
# regtest
# ============================================================================

#' Test Hypotheses About Regression Parameters
#'
#' Tests the hypothesis that a subset of regression parameters (slopes) equals
#' specified constants using a bootstrap-based confidence ellipsoid approach with
#' Mahalanobis distance.
#'
#' @inheritParams common_params
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#' @param nboot Number of bootstrap samples (default: 600).
#' @param plotit Logical; if TRUE and testing 2 parameters, plots the confidence
#'   ellipsoid with null hypothesis point (default: TRUE).
#' @param grp Vector of predictor indices to test (default: all predictors).
#'   For example, \code{grp=c(1,3)} tests the 1st and 3rd slopes.
#' @param nullvec Vector of null hypothesis values corresponding to parameters
#'   in \code{grp} (default: all zeros).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param pr Logical; if TRUE, prints progress messages (default: TRUE).
#' @param ... Additional arguments passed to \code{regfun} and \code{outfun}.
#'
#' @return A list with components:
#' \item{test}{Test statistic (Mahalanobis distance).}
#' \item{crit}{Bootstrap critical value at level \code{alpha}.}
#' \item{p.value}{Bootstrap p-value.}
#' \item{nullvec}{Vector of null hypothesis values tested.}
#' \item{est}{Estimated values of the tested parameters.}
#' \item{n}{Sample size.}
#'
#' @details
#' By default, tests whether all slope parameters equal zero (overall regression test).
#' Use \code{grp} to test specific subsets of parameters. The method constructs a
#' bootstrap confidence ellipsoid using Mahalanobis distance and compares the null
#' hypothesis point to this region.
#'
#' When \code{plotit=TRUE} and exactly 2 parameters are tested, produces a scatterplot
#' of bootstrap slope estimates with the confidence region boundary and null hypothesis
#' point marked with a square.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' @seealso \code{\link{regci}}, \code{\link{reg1way}}, \code{\link{tsreg}}
#'
#' @examples
#' # Test overall regression (all slopes = 0)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2*x[,1] - x[,2] + rnorm(50)
#' regtest(x, y, nboot=500)
#'
#' # Test first slope only
#' regtest(x, y, grp=1, nullvec=0, nboot=500)
#'
#' @export
regtest<-function(x,y,regfun=tsreg,nboot=600,alpha=.05,plotit=TRUE,
grp=c(1:ncol(x)),nullvec=c(rep(0,length(grp))),xout=FALSE,outfun=outpro,SEED=TRUE,pr=TRUE,...){
#
#  Test the hypothesis that q of the p predictors are equal to
#  some specified constants. By default, the hypothesis is that all
#  p predictors have a coefficient equal to zero.
#  The method is based on a confidence ellipsoid.
#  The critical value is determined with the percentile bootstrap method
#  in conjunction with Mahalanobis distance.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
if(pr)print("Default for outfun is now outpro")
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,regfun=regfun,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x<-as.matrix(x)
if(length(grp)!=length(nullvec))stop("The arguments grp and nullvec must have the same length.")
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun) # A p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
grp<-grp+1  #Ignore the intercept.
est<-regfun(x,y)$coef
estsub<-est[grp]
bsub<-t(bvec[grp,])
if(length(grp)==1){
m1<-sum((bvec[grp,]-est)^2)/(length(y)-1)
dis<-(bsub-estsub)^2/m1
}
if(length(grp)>1){
mvec<-apply(bsub,2,FUN=mean)
m1<-var(t(t(bsub)-mvec+estsub))
dis<-mahalanobis(bsub,estsub,m1)
}
dis2<-order(dis)
dis<-sort(dis)
critn<-floor((1-alpha)*nboot)
crit<-dis[critn]
test<-mahalanobis(t(estsub),nullvec,m1)
sig.level<-1-sum(test>dis)/nboot
if(length(grp)==2 && plotit){
plot(bsub,xlab="Parameter 1",ylab="Parameter 2")
points(nullvec[1],nullvec[2],pch=0)
xx<-bsub[dis2[1:critn],]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(test=test,crit=crit,p.value=sig.level,nullvec=nullvec,est=estsub,n=length(y))
}

# ============================================================================
# reg2ci
# ============================================================================

#' Compare Two Independent Regression Lines
#'
#' Computes bootstrap confidence intervals for the differences in regression parameters
#' (intercepts and slopes) between two independent groups.
#'
#' @param x Predictor matrix or vector for group 1 (n1 x p).
#' @param y Response vector for group 1.
#' @param x1 Predictor matrix or vector for group 2 (n2 x p).
#' @param y1 Response vector for group 2.
#' @inheritParams common_params
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#'   Any function that returns coefficients in \code{$coef}.
#' @param nboot Number of bootstrap samples (default: 599).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param plotit Logical; if TRUE and p=1, plots both regression lines (default: TRUE).
#' @param xlab X-axis label for plot (default: "X").
#' @param ylab Y-axis label for plot (default: "Y").
#' @param pr Logical; if TRUE, prints progress messages (default: FALSE).
#' @param ... Additional arguments passed to \code{regfun} and \code{outfun}.
#'
#' @return A list with components:
#' \item{n}{Vector of sample sizes for the two groups.}
#' \item{output}{Matrix with rows for each parameter difference. Columns are:
#'   \code{Parameter} (parameter index, 0=intercept), \code{ci.lower} (lower CI),
#'   \code{ci.upper} (upper CI), \code{p.value} (two-sided test),
#'   \code{Group 1} (group 1 estimate), \code{Group 2} (group 2 estimate).}
#'
#' @details
#' This function compares regression lines from two independent samples by constructing
#' bootstrap confidence intervals for the difference in each parameter. The percentile
#' bootstrap method is used.
#'
#' P-values test whether each parameter difference equals zero (i.e., whether the
#' corresponding intercepts or slopes are equal across groups).
#'
#' When \code{plotit=TRUE} and there is a single predictor, the function plots both
#' regression lines: group 1 as a solid line with points, group 2 as a dashed line
#' with "+" symbols.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' @seealso \code{\link{reg1way}}, \code{\link{regci}}, \code{\link{tsreg}}
#'
#' @examples
#' # Compare regression lines for two groups
#' set.seed(123)
#' x1 <- rnorm(50)
#' y1 <- 2 + 3*x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- 1 + 2*x2 + rnorm(50)
#' reg2ci(x1, y1, x2, y2, nboot=500)
#'
#' @export
reg2ci<-function(x,y,x1,y1,regfun=tsreg,nboot=599,alpha=.05,plotit=TRUE,SEED=TRUE,
xout=FALSE,outfun=outpro,xlab="X",ylab="Y",pr=FALSE,...){
#
#   Compute a .95 confidence interval for the difference between the
#   the intercepts and slopes corresponding to two independent groups.
#   The default regression method is Theil-Sen.
#
#   The predictor values for the first group are
#   assumed to be in the n by p matrix x.
#   The predictors for the second group are in x1
#
#   The default number of bootstrap samples is nboot=599
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
x1<-as.matrix(x1)
xx1<-cbind(x1,y1)
xx1<-elimna(xx1)
x1<-xx1[,1:ncol(x1)]
x1<-as.matrix(x1)
y1<-xx1[,ncol(x1)+1]
x=as.matrix(x)
x1=as.matrix(x1)
if(xout){
if(pr)print("outfun now defaults to outpro rather than out")
if(identical(outfun,outblp)){
flag1=outblp(x,y,plotit=FALSE)$keep
flag2=outblp(x1,y2,plotit=FALSE)$keep
}
if(!identical(outfun,outblp)){
flag1=outfun(x,plotit=FALSE)$keep
flag2=outfun(x1,plotit=FALSE)$keep
}
x=x[flag1,]
y=y[flag1]
x1=x1[flag2,]
y1=y1[flag2]
}
n=length(y)
n[2]=length(y1)
x<-as.matrix(x)
x1<-as.matrix(x1)
est1=regfun(x,y,...)$coef
est2=regfun(x1,y1,...)$coef
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun,xout=FALSE,...) # A p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data,1,regboot,x1,y1,regfun,xout=FALSE,...)
bvec<-bvec-bvec1
p1<-ncol(x)+1
regci<-matrix(0,p1,6)
dimnames(regci)<-list(NULL,
c("Parameter","ci.lower","ci.upper","p.value","Group 1","Group 2"))
ilow<-round((alpha/2)*nboot)+1
ihi<-nboot-(ilow-1)
for(i in 1:p1){
temp<-sum(bvec[i,]<0)/nboot+sum(bvec[i,]==0)/(2*nboot)
regci[i,4]<-2*min(temp,1-temp)
bsort<-sort(bvec[i,])
regci[i,2]<-bsort[ilow]
regci[i,3]<-bsort[ihi]
regci[,1]<-c(0:ncol(x))
}
regci[,5]=est1
regci[,6]=est2
if(ncol(x)==1 && plotit){
plot(c(x,x1),c(y,y1),type="n",xlab=xlab,ylab=ylab)
points(x,y)
points(x1,y1,pch="+")
abline(regfun(x,y,...)$coef)
abline(regfun(x1,y1,...)$coef,lty=2)
}
list(n=n,output=regci)
}

# ============================================================================
# regpre
# ============================================================================
regpre<-function(x,y,regfun=lsfit,error=absfun,nboot=100,adz=TRUE,
mval=round(5*log(length(y))),model=NULL,locfun=mean,pr=FALSE,
xout=FALSE,outfun=out,STAND=TRUE,
plotit=TRUE,xlab="Model Number",ylab="Prediction Error",SEED=TRUE,...){
#
#   Estimate prediction error using the regression method
#   regfun. The .632 method is used.
#   (See Efron and Tibshirani, 1993, pp. 252--254)
#
#   The predictor values are assumed to be in the n-by-p matrix x.
#   The default number of bootstrap samples is nboot=100
#
#   Prediction error is the expected value of the function error.
#   The argument error defaults to squared error.
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
#   The default value for mval, the number of observations to resample
#   for each of the B bootstrap samples is based on results by
#   Shao (JASA, 1996, 655-665). (Resampling n vectors of observations
#   model selection may not lead to the correct model as n->infinity.
#
#   The argument model should have list mode, model[[1]] indicates
#   which predictors are used in the first model. For example, storing
#   1,4 in model[[1]] means predictors 1 and 4 are being considered.
#   If model is not specified, and number of predictors is at most 5,
#   then all models are considered.
#
#   If adz=T, added to the models to be considered is where
#   all regression slopes are zero. That is, use measure of location only
#   corresponding to
#   locfun.
#
if(pr){
print("By default, least squares regression is used, ")
print("But from Wilcox, R. R. 2008, Journal of Applied Statistics, 35, 1-8")
print("Setting regfun=tsreg appears to be a better choice for general use.")
print("That is, replace least squares with the Theil-Sen estimator")
print("Note: Default for the argument error is now absfun")
print(" meaning absolute error is used")
print("To use squared error, set error=sqfun")
}
x<-as.matrix(x)
d<-ncol(x)
p1<-d+1
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
y<-temp[,d+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(!STAND)flag<-outfun(x,plotit=FALSE,...)$keep
if(STAND)flag<-outpro(x,STAND=TRUE,plotit=FALSE)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(is.null(model)){
if(d<=5)model<-modgen(d,adz=adz)
if(d>5)model[[1]]<-c(1:ncol(x))
}
mout<-matrix(NA,length(model),5,dimnames=list(NULL,c("apparent.error",
"boot.est","err.632","var.used","rank")))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=mval*nboot,replace=TRUE),nrow=nboot)
bid<-apply(data,1,idb,length(y))
#  bid is an n by nboot matrix. If the jth bootstrap sample from
#  1, ..., mval contains the value i, bid[i,j]=0; otherwise bid[i,j]=1
for (imod in 1:length(model)){
nmod=length(model[[imod]])-1
temp=c(nmod:0)
mout[imod,4]=sum(model[[imod]]*10^temp)
if(sum(model[[imod]]==0)!=1){
xx<-x[,model[[imod]]]
xx<-as.matrix(xx)
if(sum(model[[imod]]==0)!=1)bvec<-apply(data,1,regpres1,xx,y,regfun,mval,...)
# bvec is a p+1 by nboot matrix. The first row
# contains the bootstrap intercepts, the second row
# contains the bootstrap values for first predictor, etc.
if(sum(model[[imod]]==0)!=1)yhat<-cbind(1,xx)%*%bvec
if(sum(model[[imod]]==0)==1){
bvec0<-matrix(0,nrow=p1,ncol=nboot)
for(it in 1:nboot){
bvec0[1,it]<-locfun(y[data[it,]])
}
yhat<-cbind(1,x)%*%bvec0
}
# yhat is n by nboot matrix of predicted values based on
                           # bootstrap regressions.
bi<-apply(bid,1,sum) # B sub i in notation of Efron and Tibshirani, p. 253
temp<-(bid*(yhat-y))
diff<-apply(temp,1,error)
ep0<-sum(diff/bi)/length(y)
aperror<-error(regfun(xx,y,...)$resid)/length(y) # apparent error
regpre<-.368*aperror+.632*ep0
mout[imod,1]<-aperror
mout[imod,3]<-regpre
temp<-yhat-y
diff<-apply(temp,1,error)
mout[imod,2]<-sum(diff)/(nboot*length(y))
}
if(sum(model[[imod]]==0)==1){
mout[imod,3]<-locpre(y,error=error,est=locfun,SEED=SEED,mval=mval)
}}
mout[,5]=rank(mout[,3])
if(plotit)plot(c(1:nrow(mout)),mout[,3],xlab=xlab,ylab=ylab)
list(estimates=mout)
}

# ============================================================================
# regout
# ============================================================================

#' Detect Regression Outliers Using Boxplot Rule on Residuals
#'
#' Identifies regression outliers by fitting a robust regression line and
#' applying a boxplot rule to the residuals. Points with outlying residuals
#' are flagged and optionally plotted.
#'
#' @param x Numeric vector or matrix of predictor values.
#' @param y Numeric vector of response values.
#' @param regest Robust regression function to use for fitting (default: \code{stsreg}).
#' @param plotit Logical; if TRUE, creates a scatterplot with the regression
#'   line and marks outliers with circles (default: TRUE).
#' @param mbox Logical; if TRUE, uses Carling's modification of the boxplot
#'   rule; if FALSE, uses ideal fourths with conventional boxplot rule
#'   (default: TRUE).
#'
#' @return A list with component:
#' \describe{
#'   \item{out.id}{Vector of indices for outlying observations.}
#' }
#'
#' @details
#' This function:
#' 1. Fits a robust regression line using \code{regest}
#' 2. Computes residuals from the fit
#' 3. Applies \code{\link{outbox}} to identify outlying residuals
#' 4. Returns indices of outliers and optionally plots the data
#'
#' The boxplot rule identifies outliers as observations with residuals beyond
#' the "fences" (typically 1.5 × IQR from the quartiles). When \code{mbox=TRUE},
#' Carling's modification adjusts the rule to better control the false positive
#' rate under normality.
#'
#' Outliers are marked with "o" symbols in the plot, making it easy to visually
#' identify observations that deviate from the fitted robust regression line.
#'
#' @seealso \code{\link{outbox}} for boxplot-based outlier detection,
#'   \code{\link{stsreg}} for S-type Theil-Sen regression,
#'   \code{\link{reglev}} for leverage point detection
#'
#' @examples
#' \dontrun{
#' # Generate data with outliers
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' y[1:3] <- y[1:3] + 5  # Add outliers
#'
#' # Detect and plot outliers
#' result <- regout(x, y)
#' print(result$out.id)  # Indices of outliers
#'
#' # Without plotting
#' regout(x, y, plotit=FALSE)
#'
#' # Using different regression estimator
#' regout(x, y, regest=tsreg)
#' }
#'
#' @export
regout<-function(x,y,regest=stsreg,plotit=TRUE,mbox=TRUE){
#
# Check for regression outliers by fitting a
# a line to data using regest and then applying
# a boxplot rule to the residuals.
# mbox=T uses Carling's method
# mbox=F uses ideal fourths with conventional boxplot rules.
#
chk<-regest(x,y)
flag<-outbox(chk$residuals,mbox=mbox)$out.id
if(plotit){
plot(x,y)
points(x[flag],y[flag],pch="o")
abline(chk$coef)
}
list(out.id=flag)
}

# ============================================================================
# stsregp1
# ============================================================================

#' S-Type Theil-Sen Regression (Single Predictor)
#'
#' Computes an S-type modification of the Theil-Sen regression estimator for
#' a single predictor variable. The modification chooses the slope that minimizes
#' a robust scale measure of the residuals.
#'
#' @param x Numeric vector of predictor values (single predictor only).
#' @param y Numeric vector of response values.
#' @param sc Robust scale estimator to minimize (default: \code{pbvar}).
#' @inheritParams common-params
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of coefficients (intercept, slope).}
#'   \item{residuals}{Residuals from the fit.}
#' }
#'
#' @details
#' This function is a specialized version of S-type regression for the single
#' predictor case. It:
#' 1. Computes all pairwise slopes (as in standard Theil-Sen)
#' 2. For each candidate slope, calculates the robust scale (\code{sc}) of residuals
#' 3. Selects the slope that minimizes this scale
#' 4. Determines the intercept using medians
#'
#' The S-type modification provides better efficiency than median-of-slopes
#' while maintaining high breakdown point. The scale function \code{sc} is
#' typically a robust variance estimator like percentage bend variance
#' (\code{\link{pbvar}}).
#'
#' **Note**: This version is limited to a single predictor. For multiple predictors,
#' use \code{\link{stsreg}} which employs a Gauss-Seidel algorithm.
#'
#' @seealso \code{\link{stsreg}} for general S-type Theil-Sen regression,
#'   \code{\link{tsreg}} for standard Theil-Sen,
#'   \code{\link{pbvar}} for percentage bend variance
#'
#' @examples
#' \dontrun{
#' # Simple linear regression
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50, sd=0.5)
#' stsregp1(x, y)
#'
#' # With outliers - S-type is more efficient than median-of-slopes
#' y[1:5] <- y[1:5] + 5
#' fit <- stsregp1(x, y)
#' print(fit$coef)
#'
#' # Using different scale estimator
#' stsregp1(x, y, sc=winvar)
#' }
#'
#' @export
stsregp1<-function(x,y,sc=pbvar,xout=FALSE,outfun=out,...){
#
# Compute the S-type modification of
# the Theil-Sen regression estimator.
# Only a single predictor is allowed in this version
#
xy=elimna(cbind(x,y))
p=ncol(as.matrix(x))
if(p!=1)stop("Current version is limited to one predictor")
p1=p+1
x=xy[,1:p]
y=xy[,p1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
slope<-v1/v2
allvar<-NA
for(i in 1:length(slope))allvar[i]<-sc(y-slope[i]*x,...)
temp<-order(allvar)
coef<-0
coef[2]<-slope[temp[1]]
coef[1]<-median(y)-coef[2]*median(x)
res<-y-coef[2]*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# stsreg
# ============================================================================

#' S-Type Theil-Sen Regression (General)
#'
#' Computes S-type regression coefficients using a modification of the
#' Theil-Sen estimator. For multiple predictors, uses a Gauss-Seidel
#' algorithm to iteratively refine slope estimates.
#'
#' @param x Numeric vector or matrix of predictor values (n by p).
#' @param y Numeric vector of response values (length n).
#' @inheritParams common-params
#' @param iter Maximum number of Gauss-Seidel iterations for multiple
#'   predictors (default: 10).
#' @param sc Robust scale function for selecting slopes (default: \code{pbvar}).
#' @param varfun Robust variance function for final residuals (default: \code{pbvar}).
#' @param corfun Robust correlation function used in iterations (default: \code{pbcor}).
#' @param plotit Logical; currently unused (default: FALSE).
#' @param ... Additional arguments passed to \code{outfun}.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{residuals}{Residuals from the final fit.}
#'   \item{Strength.Assoc}{Measures of association strength (correlation-like).}
#' }
#'
#' @details
#' The S-type Theil-Sen estimator improves upon the standard Theil-Sen
#' (median-of-pairwise-slopes) by selecting slopes that minimize a robust
#' scale measure of residuals. This provides:
#' \itemize{
#'   \item Higher efficiency than standard Theil-Sen
#'   \item Maintained high breakdown point (up to 50%)
#'   \item Better performance with skewed error distributions
#' }
#'
#' **Algorithm**:
#' - **Single predictor**: Uses \code{\link{stsregp1}} to compute all pairwise
#'   slopes and select the one minimizing \code{sc(residuals)}.
#' - **Multiple predictors**: Employs Gauss-Seidel iterations:
#'   1. Initialize slopes using \code{\link{tsp1reg}} for each predictor
#'   2. Iteratively update each slope holding others fixed
#'   3. Continue until convergence or \code{iter} iterations reached
#'
#' The correlation function \code{corfun} is used to compute association
#' strength measures returned in \code{Strength.Assoc}.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{stsregp1}} for single predictor version,
#'   \code{\link{tsreg}} for standard Theil-Sen,
#'   \code{\link{tstsreg}} for S-type with outlier removal
#'
#' @examples
#' \dontrun{
#' # Single predictor
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' stsreg(x, y)
#'
#' # Multiple predictors
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50)
#' fit <- stsreg(x, y)
#' print(fit$coef)
#'
#' # With outlier removal
#' stsreg(x, y, xout=TRUE)
#' }
#'
#' @export
stsreg<-function(x,y,xout=FALSE,outfun=outpro,iter=10,sc=pbvar,varfun=pbvar,
corfun=pbcor,plotit=FALSE,...){
#
#  Compute Theil-Sen regression estimator
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(ncol(x)==1){
temp1<-stsregp1(x,y,sc=sc)
coef<-temp1$coef
res<-temp1$res
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-tsp1reg(x[,p],y)$coef[2]
}
res<-y-x%*%temp
alpha<-median(res)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-stsregp1(x[,p],r[,p],sc=sc)$coef[2]
}
alpha<-median(y-x%*%temp)
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
yhat<-y-res
stre=NULL
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}
list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

# ============================================================================
# rregci
# ============================================================================

#' Bootstrap Confidence Intervals for Robust Regression
#'
#' Computes bootstrap confidence intervals for regression parameters using
#' robust M-regression methods. Default is the Coakley-Hettmansperger estimator
#' (bounded influence M-regression with Schweppe weights).
#'
#' @inheritParams common_params
#' @param regfun Robust regression function to use (default: \code{chreg} for
#'   Coakley-Hettmansperger). Any function that returns coefficients in \code{$coef}.
#' @param nboot Number of bootstrap samples (default: 599).
#' @param ... Additional arguments passed to \code{regfun}.
#'
#' @return A matrix with one row per parameter (intercept, then slopes) and columns:
#' \item{ci.lower}{Lower confidence interval bound.}
#' \item{ci.upper}{Upper confidence interval bound.}
#' \item{Estimate}{Point estimate of the parameter.}
#' \item{S.E.}{Bootstrap standard error.}
#' \item{p-value}{Two-sided p-value testing if parameter equals zero.}
#'
#' @details
#' This function provides bootstrap confidence intervals for robust regression
#' estimators. The default Coakley-Hettmansperger estimator (\code{chreg}) uses
#' bounded influence M-regression with Schweppe weights, providing good robustness
#' to outliers and leverage points.
#'
#' Other robust regression functions can be specified via \code{regfun}, such as
#' \code{bireg} (bisquare M-regression) or \code{bmreg} (modified M-regression).
#'
#' For least squares regression with n < 250, use \code{lsfitci} instead for
#' better performance.
#'
#' @references
#' Coakley, C.W., & Hettmansperger, T.P. (1993). A bounded influence, high breakdown,
#' efficient regression estimator. \emph{Journal of the American Statistical Association},
#' 88, 872-880.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' @seealso \code{\link{chreg}}, \code{\link{bireg}}, \code{\link{bmreg}},
#' \code{\link{regci}}, \code{\link{lsfitci}}
#'
#' @examples
#' # Robust regression CIs with outliers
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50)
#' y[1:3] <- y[1:3] + 10  # Add outliers
#' rregci(x, y, nboot=500)
#'
#' # Multiple predictors with bisquare regression
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(50)
#' rregci(x, y, regfun=bireg, nboot=500)
#'
#' @export
rregci<-function(x,y,regfun=chreg,nboot=599,alpha=.05, ...){
#
#   Compute a .95 confidence interval for each of the parameters of
#   a linear regression equation. The default regression method is
#   a bounded influence M-regression with Schweppe weights
#   (the Coakley-Hettmansperger estimator).
#
#   When using the least squares estimator, and when n<250, use
#   lsfitci instead.
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=599
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimated of
#   the first predictor, etc.
#
x<-as.matrix(x)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun,...) # A p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
p1<-ncol(x)+1
regci<-matrix(0,p1,2)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
se<-NA
pvec<-NA
for(i in 1:p1){
bsort<-sort(bvec[i,])
pvec[i]<-sum(bvec[i,]<0)/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
regci[i,1]<-bsort[ilow]
regci[i,2]<-bsort[ihi]
se[i]<-sqrt(var(bvec[i,]))
}
pvec<-2*pvec
list(regci=regci,p.value=pvec,se=se)
}

# ============================================================================
# wsp1reg
# ============================================================================
wsp1reg<-function(x,y,plotit=FALSE){
#
# Compute the Wilcoxon R estimate of the slope
# Only a single predictor is allowed in this version
#
temp<-matrix(c(x,y),ncol=2)
temp<-elimna(temp)     # Remove any pairs with missing values
x<-temp[,1]
y<-temp[,2]
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
slope<-v1/v2
tmin<-wrregfun(slope[1],x,y)
ikeep<-1
for(i in 2:length(slope)){
tryit<-wrregfun(slope[i],x,y)
if(tryit<tmin){
tmin<-tryit
ikeep<-i
}}
coef<-NA
coef[2]<-slope[ikeep]
coef[1]<-median(y-coef[2]*x)
if(plotit){
plot(x,y,xlab="X",ylab="Y")
abline(coef)
}
res<-y-coef[2]*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# wreg
# ============================================================================

#' Wilcoxon Rank-Based Regression
#'
#' Computes the Wilcoxon R-estimator for regression, a rank-based robust
#' regression method. For multiple predictors, uses the Gauss-Seidel algorithm
#' to iterate to a solution.
#'
#' @param x A numeric matrix of predictor variables (n by p), or a vector for
#'   simple linear regression.
#' @param y A numeric vector of the dependent variable (length n).
#' @param iter Maximum number of iterations for the Gauss-Seidel algorithm when
#'   p > 1 (default: 10).
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients: \code{c(intercept, slopes)}.}
#'   \item{residuals}{Residuals from the fit.}
#' }
#'
#' @details
#' The Wilcoxon R-estimator minimizes the dispersion of residuals based on
#' Wilcoxon rank scores, providing a robust alternative to least squares.
#'
#' **For simple regression (p=1)**:
#' - Uses \code{\link{wsp1reg}} to find the slope that minimizes
#'   \eqn{\sum_{i<j} |r_i - r_j|} where \eqn{r_i} are residuals
#' - Intercept is the median of residuals
#'
#' **For multiple regression (p>1)**:
#' - Uses Gauss-Seidel iterative algorithm
#' - In each iteration, updates one slope at a time holding others fixed
#' - Continues until convergence (change < 0.0001) or \code{iter} reached
#' - Initial values from univariate Wilcoxon regressions
#'
#' Properties:
#' - Distribution-free (no normality assumption)
#' - Resistant to outliers in y
#' - Asymptotically as efficient as least squares for normal data
#' - Related to Kendall's tau
#'
#' @references
#' Hettmansperger, T.P. & McKean, J.W. (2011). \emph{Robust Nonparametric
#' Statistical Methods} (2nd ed.). CRC Press.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @family robust regression functions
#' @seealso \code{\link{wsp1reg}} for simple Wilcoxon regression,
#'   \code{\link{tsreg}} for Theil-Sen regression
#' @export
#' @examples
#' \dontrun{
#' # Simple Wilcoxon regression
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' wreg(x, y)
#'
#' # Multiple regression
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50)
#' wreg(x, y)
#'
#' # With outliers
#' y[1:3] <- y[1:3] + 10
#' wreg(x, y)  # Robust result
#' }
wreg<-function(x,y,iter=10){
#
#  Compute Wilcoxon R estimate
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#
#  Argument iter is used when there is more than one
#  predictor and indicates maximum number
#  of iterations used by Gauss-Seidel algoritm.
#
temp<-NA
x<-as.matrix(x)
if(ncol(x)==1){
temp1<-wsp1reg(x,y)
coef<-temp1$coef
res<-temp1$res
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-wsp1reg(x[,p],y)$coef[2]
}
res<-y-x%*%temp
alpha<-median(res)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-wsp1reg(x[,p],r[,p],plotit=FALSE)$coef[2]
}
alpha<-median(y-x%*%temp)
if(max(abs(tempold-temp))<.0001)break
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
list(coef=coef,residuals=res)
}

# ============================================================================
# opreg
# ============================================================================
#' Outlier-Pruned Regression
#'
#' Performs regression after removing outliers detected using projection-based
#' outlier detection. This two-stage approach combines robust outlier detection
#' with any chosen regression method for improved performance.
#'
#' @param x A numeric vector or matrix containing the predictor variable(s). For
#'   multiple predictors, rows represent observations and columns represent variables.
#' @inheritParams common-params
#' @param regfun Regression function to use after outlier removal (default: `tsreg`).
#'   Can be any function that returns coefficients in `$coef` and residuals in
#'   `$residuals`. Common choices: `tsreg`, `ltsreg`, `chreg`, `lsfit`.
#' @param cop Type of correlation/covariance measure used in projection-based
#'   outlier detection (default: 3 for Winsorized correlation). See \code{\link{outpro}}
#'   for details. Options:
#'   \itemize{
#'     \item 1: Pearson correlation
#'     \item 2: Median-based correlation
#'     \item 3: Winsorized correlation
#'     \item 4: MVE (Minimum Volume Ellipsoid)
#'   }
#' @param MC Logical. If `TRUE`, uses parallel processing for outlier detection
#'   via \code{\link{outproMC}} (default: `FALSE`).
#' @param varfun Function used to compute variance for strength of association
#'   (default: `pbvar`).
#' @param corfun Correlation function used for computing explanatory power when
#'   the variance ratio exceeds 1 (default: `pbcor`).
#' @param STAND Logical. If `TRUE`, standardizes the data before applying outlier
#'   detection (default: `TRUE`). Recommended to keep as `TRUE`.
#' @param xout Logical. Not used in this function but included to avoid conflicts
#'   when using with other functions like `regci` (default: `FALSE`).
#'
#' @return A list with components:
#'   \item{coef}{Vector of regression coefficients from `regfun`}
#'   \item{residuals}{Vector of residuals}
#'   \item{Strength.Assoc}{Strength of association (square root of explanatory power)}
#'   \item{Explanatory.Power}{Proportion of variance explained}
#'
#' @details
#' This function implements a two-stage robust regression approach:
#' \enumerate{
#'   \item Identify and remove outliers using projection-based outlier detection
#'         (\code{\link{outpro}} or \code{\link{outproMC}})
#'   \item Apply the specified regression method to the cleaned data
#' }
#'
#' The projection-based outlier detection examines many random projections of
#' the multivariate data and flags points that appear as outliers in multiple
#' projections. This is particularly effective for detecting outliers in
#' high-dimensional data.
#'
#' Strength of association is computed as the square root of the ratio of
#' variance of fitted values to variance of y. If this ratio exceeds 1, the
#' squared correlation is used instead.
#'
#' @references
#' Wilcox, R. R. (2022). \emph{Introduction to Robust Estimation and Hypothesis
#' Testing} (5th ed.). Academic Press. Chapter 10.
#'
#' Stahel, W. A. (1981). Robuste Schätzungen: infinitesimale Optimalität und
#' Schätzungen von Kovarianzmatrizen. PhD thesis, ETH Zürich.
#'
#' @seealso \code{\link{outpro}} for outlier detection details,
#'   \code{\link{opregMC}} for parallel version with different output,
#'   \code{\link{mopreg}} for multiple outcome regression,
#'   \code{\link{tsreg}} for Theil-Sen regression
#'
#' @examples
#' # Simple regression with outlier removal
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50)
#' y[1:3] <- y[1:3] + 10  # Add some outliers
#' result <- opreg(x, y)
#' result$coef
#'
#' # Multiple regression with different regression method
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - 3*x[,2] + rnorm(50)
#' opreg(x, y, regfun=ltsreg)
#'
#' # Using parallel processing for large datasets
#' opreg(x, y, MC=TRUE)
#'
#' @export
opreg<-function(x,y,regfun=tsreg,cop=3,MC=FALSE,varfun=pbvar,corfun=pbcor,STAND=TRUE,xout=FALSE){
#
# Do regression on points not labled outliers
# using projection-type outlier detection method
#
#  Note: argument xout is not relevant here, but is included to avoid conflicts when using regci.
#
x<-as.matrix(x)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
if(!MC)ivec<-outpro(m,plotit=FALSE,cop=cop,STAND=STAND)$keep
if(MC)ivec<-outproMC(m,plotit=FALSE,cop=cop,STAND=STAND)$keep
np1<-ncol(x)+1
coef<-regfun(m[ivec,1:ncol(x)],m[ivec,np1])$coef
vec<-rep(1,length(y))
residuals<-y-cbind(vec,x)%*%coef
stre=NULL
yhat<-y-residuals
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
stre=sqrt(e.pow)
}
list(coef=coef,residuals=residuals,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

# ============================================================================
# wrregfun
# ============================================================================
wrregfun<-function(slope,x=x,y=y){
x<-as.matrix(x)
res<-y-x%*%slope
v1<-rank(res)
v2<-sqrt(12)*(v1/(length(y)+1)-.5)
wrregfun<-sum(v2*res)
wrregfun
}

# ============================================================================
# opregpb
# ============================================================================

#' Bootstrap Inference for Outlier-Pruned Regression
#'
#' Performs bootstrap hypothesis testing and confidence intervals for regression
#' parameters after removing outliers using projection-type outlier detection.
#' Combines Theil-Sen regression with robust outlier removal.
#'
#' @inheritParams common_params
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param om Logical; if TRUE and p > 1, performs omnibus test of all parameters
#'   (default: TRUE). If FALSE, only individual parameter tests.
#' @param ADJ Logical; if TRUE, adjusts p-values as described in Wilcox (2022)
#'   Section 11.1.5 (default: TRUE).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param nullvec Vector of null hypothesis values for each parameter, including
#'   intercept (default: all zeros).
#' @param plotit Logical; if TRUE, creates diagnostic plots (default: TRUE).
#' @param opdis Outlier distance type: 1 for projection distance, 2 for Mahalanobis
#'   distance (default: 2).
#' @param gval Critical value for projection-type outlier detection
#'   (default: sqrt(qchisq(0.95, p+1))).
#'
#' @return A list with components:
#' \item{n}{Sample size.}
#' \item{n.keep}{Number of observations retained after outlier removal.}
#' \item{regci}{Matrix of CIs and p-values for each parameter.}
#' \item{test}{Omnibus test statistic (if \code{om=TRUE} and p > 1).}
#' \item{crit}{Critical value for omnibus test.}
#' \item{p.value}{P-value for omnibus test.}
#' \item{adjusted.p.value}{Adjusted omnibus p-value (if \code{ADJ=TRUE}).}
#'
#' @details
#' This function implements a two-stage procedure:
#' \enumerate{
#'   \item{Detect and remove outliers using projection-based methods (default: Mahalanobis distance)}
#'   \item{Apply Theil-Sen regression to the remaining data}
#'   \item{Use bootstrap to construct CIs and test hypotheses}
#' }
#'
#' When \code{om=TRUE} and there are multiple predictors, an omnibus test examines
#' whether all parameters simultaneously equal their null values. Individual parameter
#' tests are always provided.
#'
#' The \code{ADJ=TRUE} option applies corrections described in Wilcox (2022) to improve
#' Type I error control when outliers are present.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press. (See Section 11.1.5)
#'
#' @seealso \code{\link{opreg}}, \code{\link{opregpbMC}} (parallel version),
#' \code{\link{regci}}, \code{\link{outpro}}
#'
#' @examples
#' # Bootstrap inference with outlier removal
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2 + 3*x + rnorm(100)
#' y[1:5] <- y[1:5] + 20  # Add outliers
#' opregpb(x, y, nboot=500)
#'
#' @export
opregpb<-function(x,y,nboot=1000,alpha=.05,om=TRUE,ADJ=TRUE,SEED=TRUE,
nullvec=rep(0,ncol(x)+1),plotit=TRUE,opdis=2,gval=sqrt(qchisq(.95,ncol(x)+1))){
#
# generate bootstrap estimates
# use projection-type outlier detection method followed by
# TS regression.
#
# om=T and ncol(x)>1, means an omnibus test is performed,
# otherwise only individual tests of parameters are performed.
#
# opdis=2, means that Mahalanobis distance is used
# opdis=1, means projection-type distance is used
#
# gval is critical value for projection-type outlier detection
# method
#
# ADJ=T, Adjust p-values as described in Section 11.1.5 of the text.
#
if(SEED)set.seed(2)
x<-as.matrix(x)
m<-cbind(x,y)
p1<-ncol(x)+1
m<-elimna(m) # eliminate any rows with missing data
x<-m[,1:ncol(x)]
x<-as.matrix(x)
y<-m[,p1]
if(nrow(x)!=length(y))stop("Sample size of x differs from sample size of y")
if(!is.matrix(x))stop("Data should be stored in a matrix")
print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun=opreg)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
# using Hochberg method
bvec<-t(bvec)
dvec<-alpha/(c(1:ncol(x)))
test<-NA
icl0<-round(alpha*nboot/2)
icl<-round(alpha*nboot/(2*ncol(x)))
icu0<-nboot-icl0
icu<-nboot-icl
output<-matrix(0,p1,6)
vlabs="Intercept"
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(output)<-list(vlabs,c("Param.","p.value","p.crit",
"ci.lower","ci.upper","s.e."))
pval<-NA
for(i in 1:p1){
output[i,1]<-i-1
se.val<-var(bvec[,i])
temp<-sort(bvec[,i])
output[i,6]<-sqrt(se.val)
if(i==1){
output[i,4]<-temp[icl0+1]
output[i,5]<-temp[icu0]
}
if(i>1){
output[i,4]<-temp[icl+1]
output[i,5]<-temp[icu]
}
pval[i]<-sum((temp>nullvec[i]))/length(temp)
if(pval[i]>.5)pval[i]<-1-pval[i]
}
fac<-2
if(ADJ){
# Adjust p-value if n<60
nval<-length(y)
if(nval<20)nval<-20
if(nval>60)nval<-60
fac<-2-(60-nval)/40
}
pval[1]<-2*pval[1]
pval[2:p1]<-fac*pval[2:p1]
output[,2]<-pval
temp2<-order(0-pval[2:p1])
zvec<-dvec[1:ncol(x)]
sigvec<-(test[temp2]>=zvec)
output[temp2+1,3]<-zvec
output[1,3]<-NA
output[,2]<-pval
om.pval<-NA
temp<-opreg(x,y)$coef
if(om && ncol(x)>1){
temp2<-rbind(bvec[,2:p1],nullvec[2:p1])
if(opdis==1)dis<-pdis(temp2,center=temp[2:p1])
if(opdis==2){
cmat<-var(bvec[,2:p1]-apply(bvec[,2:p1],2,mean)+temp[2:p1])
dis<-mahalanobis(temp2,temp[2:p1],cmat)
}
om.pval<-sum((dis[nboot+1]<=dis[1:nboot]))/nboot
}
# do adjusted p-value
nval<-length(y)
if(nval<20)nval<-20
if(nval>60)nval<-60
adj.pval<-om.pval/2+(om.pval-om.pval/2)*(nval-20)/40
if(ncol(x)==2 && plotit){
plot(bvec[,2],bvec[,3],xlab="Slope 1",ylab="Slope 2")
temp.dis<-order(dis[1:nboot])
ic<-round((1-alpha)*nboot)
xx<-bvec[temp.dis[1:ic],2:3]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(output=output,om.pval=om.pval,adj.om.pval=adj.pval)
}

# ============================================================================
# ltsgreg
# ============================================================================

#' Least Trimmed Absolute Value Regression
#'
#' Computes a robust regression estimator that minimizes the sum of the smallest
#' absolute residuals. This is a variant of least trimmed squares (LTS) regression
#' using absolute values instead of squared residuals.
#'
#' @param x Numeric vector or matrix of predictor variable(s).
#' @param y Numeric vector of the response variable.
#' @param tr Trimming proportion (default: 0.2). The fraction of largest absolute
#'   residuals to trim. Must be between 0 and 0.5.
#' @param h Number of observations to use in the fit (default: \code{NA}, computed
#'   as \code{n - floor(tr*n)}). If specified, overrides \code{tr}.
#' @inheritParams common-params
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes)}
#'   \item{residuals}{Vector of residuals from the fitted model}
#' }
#'
#' @details
#' This function minimizes the sum of the \code{h} smallest absolute residuals,
#' where \code{h = n - floor(tr*n)}. This provides robustness to outliers while
#' using absolute deviations instead of squared deviations.
#'
#' **Method**:
#' 1. Uses \code{ltsReg} from the \code{robustbase} package for initial estimates
#' 2. Refines estimates using Nelder-Mead optimization
#' 3. Minimizes: sum of h smallest |residual| values
#'
#' **Properties**:
#' - **Breakdown point**: Approximately (1-h/n) ≈ \code{tr}
#' - **Robustness**: High resistance to outliers
#' - **L1 criterion**: Uses absolute values (more robust than L2/squared)
#'
#' **When to use**:
#' - Heavy-tailed error distributions
#' - Presence of outliers in Y direction
#' - Want breakdown protection with L1-type estimator
#' - Alternative to \code{\link{ltsreg}} (which uses L2)
#'
#' **Trimming proportion**:
#' - Default (0.2): 20% trimming, good balance of robustness and efficiency
#' - Larger \code{tr}: more robustness, lower efficiency
#' - Smaller \code{tr}: higher efficiency, less robustness
#' - Maximum \code{tr}: 0.5 (50% breakdown point)
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{ltsreg}} for least trimmed squares (L2 version),
#'   \code{\link{tsreg}} for Theil-Sen regression,
#'   \code{\link{opreg}} for outlier-pruned regression
#'
#' @examples
#' \dontrun{
#' # Simple robust regression
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50)
#'
#' # Add outliers
#' y[1:5] <- y[1:5] + 10
#'
#' # Fit with 20% trimming
#' fit <- ltsgreg(x, y, tr=0.2)
#' print(fit$coef)
#'
#' # More aggressive trimming (30%)
#' fit2 <- ltsgreg(x, y, tr=0.3)
#' print(fit2$coef)
#'
#' # Multiple predictors
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rt(50, df=3)
#' ltsgreg(x, y)
#' }
#'
#' @export
ltsgreg<-function(x, y, tr = 0.2, h = NA,xout=FALSE,outfun=outpro,...)
{
        #
        # Compute the least trimmed absolute value regression estimator.
        # The default amount of trimming is .2
        #
        #  Can simply use ltsreg with tr=amount of trimming.
        #
x<-as.matrix(x)
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
p<-np-1
x<-X[,1:p]
x<-as.matrix(x)
y<-X[,np]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
X<-cbind(x,y)
}
if(is.na(h)) h <- length(y) - floor(tr * length(y))
temp<-ltsReg(y~x)$coef
START<-temp
coef<-nelderv2(X,np,FN=lts.sub,START=START,h=h,p=p)
        res <- y - x%*%coef[2:np] - coef[1]
        list(coef = coef, residuals = res)
}

# ============================================================================
# locreg
# ============================================================================
qreg.sub<-function(X,theta,qval=.5){
np<-ncol(X)
p<-np-1
x<-X[,1:p]
y<-X[,np]
temp<-t(t(x)*theta[2:np])
yhat<-apply(temp,1,sum)+theta[1]
res<-y-yhat
flag<-(res<=0)
rval<-(qval-flag)*res
val<-sum(rval)
val
}

# ============================================================================
# qreg
# ============================================================================

#' Quantile Regression using Koenker-Bassett Method
#'
#' Computes quantile regression using the Koenker-Bassett approach. This function
#' estimates the \code{qval}-th conditional quantile of Y given X. Default uses
#' the \code{rq} function from the \code{quantreg} package.
#'
#' @param x Numeric vector or matrix of predictor variable(s). For multiple
#'   predictors, rows represent observations and columns represent variables.
#' @param y Numeric vector of the dependent variable.
#' @param qval Quantile to estimate (default: 0.5 for median regression).
#'   Must be between 0 and 1.
#' @param q Alternative parameter name for \code{qval} (default: \code{NULL}).
#'   If specified, overrides \code{qval}.
#' @param pr Logical. If \code{TRUE}, prints progress messages (default: \code{FALSE}).
#' @inheritParams common-params
#' @param plotit Logical. If \code{TRUE} and there is a single predictor, creates
#'   a scatterplot with the quantile regression line (default: \code{FALSE}).
#' @param xlab Label for x-axis when \code{plotit=TRUE} (default: "X").
#' @param ylab Label for y-axis when \code{plotit=TRUE} (default: "Y").
#' @param op Operation code for the old version (default: 1). Only used when \code{v2=FALSE}.
#' @param v2 Logical. If \code{TRUE}, uses the \code{rq} function from
#'   \code{quantreg} package (default: \code{TRUE}). If \code{FALSE}, uses
#'   an older, slower implementation.
#' @param method Fitting method passed to \code{rq} (default: 'br' for Barrodale-Roberts).
#'   Can also use 'scad' for variable selection (Wu & Liu, 2009).
#' @param WARN Logical. If \code{FALSE}, suppresses warnings from \code{rq}
#'   (default: \code{FALSE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes)}
#'   \item{residuals}{Vector of residuals from the quantile regression fit}
#' }
#'
#' @details
#' Quantile regression estimates conditional quantiles rather than conditional
#' means. For \code{qval=0.5}, this gives median regression, which is robust to
#' outliers in Y. Other quantiles can be used to model different parts of the
#' conditional distribution.
#'
#' The Koenker-Bassett approach minimizes a weighted sum of absolute deviations:
#' \deqn{\sum_{i} \rho_{\tau}(y_i - x_i'\beta)}
#' where \eqn{\rho_{\tau}} is the check function that weights positive and
#' negative residuals differently based on the quantile \eqn{\tau}.
#'
#' When \code{v2=TRUE} (default), the function uses the highly optimized \code{rq}
#' function from the \code{quantreg} package. The \code{method} parameter controls
#' the algorithm used. The 'scad' method implements variable selection via SCAD
#' penalty (Wu & Liu, 2009).
#'
#' @references
#' Koenker, R. & Bassett, G. (1978). Regression quantiles. \emph{Econometrica}, 46, 33-50.
#'
#' Wu, Y. & Liu, Y. (2009). Variable selection in quantile regression.
#' \emph{Statistica Sinica}, 19, 801-817.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{Qreg}} for an alternative quantile regression implementation,
#'   \code{\link[quantreg]{rq}} for the underlying quantreg function
#'
#' @examples
#' # Median regression (50th percentile)
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2 + 3*x + rt(100, df=3)  # Heavy-tailed errors
#' result <- qreg(x, y, qval=0.5)
#' result$coef
#'
#' # Compare different quantiles
#' q10 <- qreg(x, y, qval=0.1)
#' q50 <- qreg(x, y, qval=0.5)
#' q90 <- qreg(x, y, qval=0.9)
#'
#' # Multiple predictors with 75th percentile
#' x <- matrix(rnorm(200), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(100)
#' qreg(x, y, qval=0.75)
#'
#' @export
qreg<-function(x, y,qval=.5, q=NULL,pr=FALSE,xout=FALSE, outfun=outpro,plotit=FALSE,xlab="X",ylab="Y",op=1,v2=TRUE,method='br',WARN=FALSE,...)
{
#
# Compute the quantile regression line. That is, the goal is to
# determine the qth (qval) quantile of Y given X using the
#  the Koenker-Bassett approach.
#
#  v2=T, uses the function rq in the R library quantreg
#  v2=F, uses an older and slower version
#  op=1 has to do with the old version.
#
#  method=scad, see Wu and  Liu (2009). VARIABLE SELECTION IN QUANTILE REGRESSION, Statistica Sinica 19, 801-817.
#
if(!is.null(q))qval=q
x<-as.matrix(x)
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
p<-np-1
x<-X[,1:p]
x<-as.matrix(x)
y<-X[,np]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(!v2){
temp<-ltareg(x,y,0,op=op)
if(qval==.5){
coef<-temp$coef
res<-temp$res
}
if(qval!=.5){
START<-temp$coef
coef<-nelderv2(X,np,FN=qreg.sub,START=START,qval=qval)
}}
if(v2){
library(quantreg)
x<-as.matrix(x)
if(!WARN)options(warn=-1)
temp<-rq(y~x,tau=qval,method=method)
coef<-temp[1]$coefficients
if(!WARN)options(warn=0)
}
if(ncol(x)==1){
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab)
abline(coef)
}}
res <- y - x%*%coef[2:np] - coef[1]
list(coef = coef, residuals = res)
}

# ============================================================================
# kerreg
# ============================================================================

#' Kernel Regression with Epanechnikov Kernel
#'
#' Computes local weighted regression using the Epanechnikov kernel. For one
#' predictor, the function calls \code{\link{locreg}}; for multiple predictors,
#' it implements kernel regression with automatic bandwidth selection.
#'
#' @param x A numeric vector or matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param pyhat Logical; if \code{TRUE}, returns the predicted values (default: \code{FALSE}).
#' @param pts For univariate case: points at which to compute fitted values. If
#'   \code{NA} (default), uses a grid of \code{np} equally spaced points.
#' @param plotit Logical; if \code{TRUE} (default), creates a plot of the fitted
#'   regression surface. For bivariate predictors, creates a 3D perspective plot;
#'   for univariate predictors, creates a scatterplot with fitted line.
#' @param theta Rotation angle for 3D perspective plot (default: 50). Only used
#'   when \code{plotit=TRUE} and there are 2 predictors.
#' @param phi Colatitude angle for 3D perspective plot (default: 25). Only used
#'   when \code{plotit=TRUE} and there are 2 predictors.
#' @param expand Expansion factor for z-axis in perspective plot (default: 0.5).
#' @param scale Logical; if \code{TRUE}, scales the 3D plot to fit in the viewing box
#'   (default: \code{FALSE}).
#' @param zscale Logical; if \code{TRUE}, standardizes all variables (x and y) using
#'   median and MAD before fitting (default: \code{FALSE}).
#' @param eout Logical; if \code{TRUE}, removes outliers from the combined (x,y)
#'   data before fitting (default: \code{FALSE}).
#' @param xout Logical; if \code{TRUE}, removes outliers from predictor space
#'   before fitting (default: \code{FALSE}). Cannot be used with \code{eout=TRUE}.
#' @param outfun Function for detecting outliers (default: \code{out}). Only used
#'   if \code{xout=TRUE} or \code{eout=TRUE}.
#' @param np For univariate case: number of points for plotting (default: 100).
#' @param xlab Label for x-axis (default: "X").
#' @param ylab Label for y-axis (default: "Y").
#' @param zlab Label for z-axis in 3D plot (default: "Z").
#' @param varfun Function for computing variance (default: \code{pbvar}). Currently
#'   not used in returned output.
#' @param e.pow Logical; if \code{TRUE}, computes explanatory power (default: \code{TRUE}).
#'   Currently not used in returned output.
#' @param pr Logical; if \code{TRUE} (default), prints informational messages.
#' @param ticktype Type of tick marks for perspective plot (default: "simple").
#' @param pch Plotting character for scatterplot points (default: '.').
#' @param ... Additional arguments passed to the outlier detection function.
#'
#' @return
#' If \code{pyhat=TRUE}, returns a vector of fitted values. If \code{pyhat=FALSE}
#' (default), returns \code{NULL} but produces a plot if \code{plotit=TRUE}.
#'
#' @details
#' The \code{kerreg} function implements local weighted regression using the
#' Epanechnikov kernel, which is optimal in a mean squared error sense among
#' kernels of bounded support.
#'
#' **For univariate predictors (p=1):**
#' The function calls \code{\link{locreg}} internally. See \code{\link{locreg}}
#' for details on bandwidth selection and estimation.
#'
#' **For multivariate predictors (p>1):**
#' The method performs the following steps:
#' \enumerate{
#'   \item Standardizes each predictor by dividing by \code{min(SD, IQR/1.34)}
#'   \item Selects bandwidth using the formula: \code{h = A * n^(-1/(p+4))} where
#'         A is a dimension-dependent constant
#'   \item For each observation, fits a local linear regression using weighted
#'         least squares with Epanechnikov kernel weights
#'   \item The Epanechnikov kernel is: \code{K(u) = 0.5*(p+2)*(1-u^2)/c_p} for
#'         \code{u < 1}, and 0 otherwise, where \code{c_p} is a normalizing constant
#' }
#'
#' The bandwidth is chosen to balance bias and variance according to theory
#' developed by Fan (1993). The method requires at least \code{p+1} observations
#' within the kernel support to compute each local fit.
#'
#' **Outlier Handling:**
#' \itemize{
#'   \item \code{xout=TRUE}: Removes outliers in predictor space before fitting
#'   \item \code{eout=TRUE}: Removes outliers in the joint (x,y) space
#'   \item \code{zscale=TRUE}: Robust standardization before outlier detection
#' }
#'
#' **Note**: Cannot specify both \code{xout=TRUE} and \code{eout=TRUE} simultaneously.
#'
#' @references
#' Fan, J. (1993). Local linear regression smoothers and their minimax efficiencies.
#' \emph{Annals of Statistics}, 21, 196-217.
#'
#' Bjerve, S. and Doksum, K.A. (1993). Correlation curves: Measures of association
#' as functions of covariate values. \emph{Annals of Statistics}, 21, 890-902.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{locreg}} for univariate local regression,
#'   \code{\link{rplot}} for robust regression plotting,
#'   \code{\link{lplot}} for LOWESS smoothing
#'
#' @examples
#' \dontrun{
#' # Univariate kernel regression
#' set.seed(123)
#' x <- seq(0, 2*pi, length.out=100)
#' y <- sin(x) + rnorm(100, sd=0.2)
#' kerreg(x, y, plotit=TRUE)
#'
#' # Get fitted values
#' yhat <- kerreg(x, y, pyhat=TRUE, plotit=FALSE)
#'
#' # Bivariate predictors with 3D visualization
#' x <- matrix(runif(200, -2, 2), ncol=2)
#' y <- x[,1]^2 + x[,2]^2 + rnorm(100, sd=0.3)
#' kerreg(x, y, theta=45, phi=20)
#'
#' # With outlier removal
#' x <- rnorm(100)
#' y <- 2 + 3*x + rnorm(100, sd=0.5)
#' y[1:5] <- y[1:5] + 5  # Add outliers
#' kerreg(x, y, eout=TRUE)
#'
#' # Robust standardization before fitting
#' kerreg(x, y, zscale=TRUE, plotit=TRUE)
#' }
#'
#' @export
kerreg<-function(x,y,pyhat=FALSE,pts=NA,plotit=TRUE,theta=50,phi=25,expand=.5,
scale=FALSE,zscale=FALSE,eout=FALSE,xout=FALSE,outfun=out,np=100,xlab="X",ylab="Y",zlab="Z",
varfun=pbvar,e.pow=TRUE,pr=TRUE,ticktype="simple",pch='.',...){
#
# Compute local weighted regression with Epanechnikov kernel
#
# See Fan, Annals of Statistics, 1993, 21, 196-217.
# cf. Bjerve and Doksum, Annals of Statistics, 1993, 21, 890-902
#
# With a single predictor, this function calls locreg
# See locreg for information about np and plotting
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
d<-ncol(x)
np1<-d+1
m<-elimna(cbind(x,y))
if(xout && eout)stop("Can't have eout=xout=T")
if(eout){
flag<-outfun(m,plotit=FALSE,...)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
}
if(zscale){
for(j in 1:np1){
m[,j]<-(m[,j]-median(m[,j]))/mad(m[,j])
}}
x<-m[,1:d]
x<-as.matrix(x)
y<-m[,np1]
n<-nrow(x)
if(d>1){
xrem<-x
pi<-gamma(.5)^2
cd<-c(2,pi)
if(d==2)A<-1.77
if(d==3)A<-2.78
if(d>2){
for(j in 3:d)cd[j]<-2*pi*cd[j-2]/j  # p. 76
}
if(d>3)A<-(8*d*(d+2)*(d+4)*(2*sqrt(pi))^d)/((2*d+1)*cd[d])  # p. 87
hval<-A*(1/n)^(1/(d+4))  # p. 86
for(j in 1:d){
sig<-sqrt(var(x[,j]))
temp<-idealf(x[,j])
iqr<-(temp$qu-temp$ql)/1.34
A<-min(c(sig,iqr))
x[,j]<-x[,j]/A
}
xx<-cbind(rep(1,nrow(x)),x)
yhat<-NA
for(j in 1:n){
yhat[j]<-NA
temp1<-t(t(x)-x[j,])/(hval)
temp1<-temp1^2
temp1<-apply(temp1,1,FUN="sum")
temp<-.5*(d+2)*(1-temp1)/cd[d]
epan<-ifelse(temp1<1,temp,0) # Epanechnikov kernel, p. 76
chkit<-sum(epan!=0)
if(chkit >= np1){
vals<-lsfit(x,y,wt=epan)$coef
yhat[j]<-xx[j,]%*%vals
}}
if(plotit  && d==2){
if(pr){
if(!scale){
print("scale=F is specified")
print("If there is dependence, might use scale=T")
}}
m<-elimna(cbind(xrem,yhat))
xrem<-m[,1:d]
yhat<-m[,np1]
fitr<-yhat
iout<-c(1:length(fitr))
nm1<-length(fitr)-1
for(i in 1:nm1){
ip1<-i+1
for(k in ip1:length(fitr))if(sum(xrem[i,]==xrem[k,])==2)iout[k]<-0
}
fitr<-fitr[iout>=1]
mkeep<-xrem[iout>=1,]
fit<-interp(mkeep[,1],mkeep[,2],fitr)
persp(fit,theta=theta,phi=phi,expand=expand,xlab=xlab,ylab=ylab,zlab=zlab,
scale=scale,ticktype=ticktype)
}}
if(d==1){
yhat<-locreg(x[,1],y,pyhat=TRUE,np=np,plotit=plotit,pts=pts,
xlab=xlab,ylab=ylab,pch=pch)
yhat2<-locreg(x[,1],y,pyhat=TRUE,np=0,plotit=FALSE)
}
if(d>1)yhat2<-yhat
m<-NULL
#E.pow<-varfun(yhat2[!is.na(yhat2)])/varfun(y)
# Estimate of explanatory power performs poorly.
if(pyhat)m<-yhat
#list(Strength.Assoc=sqrt(E.pow),Explanatory.Power=E.pow,yhat=m)
m
}

# ============================================================================
# snmreg.sub
# ============================================================================

#' Internal Helper for Skipped Regression Optimization
#'
#' Computes the variance of residuals for a given set of regression slopes
#' in the skipped regression algorithm. Used internally by \code{\link{snmreg}}.
#'
#' @param X Matrix with predictors in first p columns and response in last column.
#' @param theta Vector of slope parameters to evaluate.
#'
#' @return Robust variance of residuals (via \code{\link{pbvar}}) for the
#'   given slopes.
#'
#' @details
#' This function is called repeatedly by the optimization routine in
#' \code{\link{snmreg}} to find slopes that minimize a robust residual variance.
#' It should not be called directly by users.
#'
#' @keywords internal
#' @noRd
snmreg.sub<-function(X,theta){
np<-ncol(X)
p<-np-1
x<-X[,1:p]
y<-X[,np]
temp<-t(t(x)*theta)
yhat<-apply(temp,1,sum)
yhat<-yhat
res<-y-yhat
val<-pbvar(res)
val
}

# ============================================================================
# tstsreg
# ============================================================================

#' Modified Theil-Sen Regression with Outlier Removal
#'
#' Computes a modified Theil-Sen regression estimator by first using an S-type
#' initial estimate to identify and remove outlying residuals, then applying
#' standard Theil-Sen regression to the remaining data.
#'
#' @param x Numeric vector or matrix of predictor values (n by p).
#' @param y Numeric vector of response values (length n).
#' @param sc Robust scale function for S-type estimation (default: \code{pbvar}).
#' @inheritParams common-params
#' @param iter Number of iterations for S-type regression (default: 5).
#' @param plotit Logical; currently unused (default: FALSE).
#' @param ... Additional arguments passed to \code{outfun}.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{residuals}{Residuals from the final Theil-Sen fit.}
#' }
#'
#' @details
#' This two-stage robust regression procedure:
#' 1. **Stage 1 (Outlier identification)**: Fits an S-type Theil-Sen regression
#'    using \code{\link{stsreg}} to get initial residuals.
#' 2. **Stage 2 (Outlier removal)**: Identifies observations with standardized
#'    residuals exceeding 2 (in MAD units: |residual - median| / MAD > 2).
#' 3. **Stage 3 (Final fit)**: Applies standard Theil-Sen regression
#'    (\code{\link{tsreg}}) to the remaining observations.
#'
#' This approach combines:
#' - The efficiency of S-type estimation for initial outlier detection
#' - The simplicity and robustness of standard Theil-Sen for the final fit
#' - Improved performance when a small fraction of outliers is present
#'
#' The method is particularly useful when you expect some contamination but
#' want to ensure the final estimates are computed on relatively clean data.
#'
#' @seealso \code{\link{stsreg}} for S-type Theil-Sen,
#'   \code{\link{tsreg}} for standard Theil-Sen,
#'   \code{\link{tssnmreg}} for combining Theil-Sen with skipped regression
#'
#' @examples
#' \dontrun{
#' # Data with a few outliers
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' y[1:3] <- y[1:3] + 10  # Add outliers
#'
#' # Standard Theil-Sen (may be affected by outliers)
#' tsreg(x, y)
#'
#' # Modified version with outlier removal
#' tstsreg(x, y)
#'
#' # Multiple predictors
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50)
#' y[1:5] <- y[1:5] + 8
#' tstsreg(x, y)
#' }
#'
#' @export
tstsreg<-function(x,y,sc=pbvar,xout=FALSE,outfun=outpro,iter=5,plotit=FALSE,...){
#
# Compute a modified Theil-Sen regression estimator.
# Use s-type initial estimate, eliminate points with
# outlying residuals, then do regular Theil-Sen
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
res=stsreg(x,y,sc=sc,iter=iter)$res
chk<-abs(res-median(res))/mad(res)
xx<-x[chk<=2,]
yy<-y[chk<=2]
temp<-tsreg(xx,yy)
list(coef=temp$coef,residuals=temp$res)
}

# ============================================================================
# tssnmreg
# ============================================================================

#' Hybrid Theil-Sen and Skipped Regression with Outlier Removal
#'
#' Combines skipped regression with Theil-Sen estimation by first using
#' skipped regression to identify outliers, then applying standard Theil-Sen
#' to the remaining clean data.
#'
#' @param x Numeric vector or matrix of predictor values (n by p).
#' @param y Numeric vector of response values (length n).
#' @param sc Robust scale function for skipped regression (default: \code{pbvar}).
#' @inheritParams common-params
#' @param plotit Logical; currently unused (default: FALSE).
#' @param ... Additional arguments passed to \code{outfun}.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{residuals}{Residuals from the final Theil-Sen fit.}
#' }
#'
#' @details
#' This hybrid approach proceeds in three stages:
#' 1. **Stage 1**: Fits skipped regression (\code{\link{snmreg}}) which uses
#'    Nelder-Mead optimization with a "skipped" correlation approach.
#' 2. **Stage 2**: Identifies outliers as observations with standardized
#'    residuals exceeding 2 (in MAD units).
#' 3. **Stage 3**: Applies standard Theil-Sen regression (\code{\link{tsreg}})
#'    to the remaining observations after outlier removal.
#'
#' **Why this hybrid?**
#' - Skipped regression (\code{snmreg}) is effective at handling outliers
#'   in both x and y directions
#' - Theil-Sen provides simple, distribution-free inference once outliers
#'   are removed
#' - The combination leverages strengths of both methods
#'
#' The method is particularly useful when:
#' - You suspect leverage points (outliers in predictor space)
#' - You want the simplicity of Theil-Sen for final inference
#' - Initial outlier identification needs to handle multivariate outliers
#'
#' @seealso \code{\link{snmreg}} for skipped regression,
#'   \code{\link{tsreg}} for standard Theil-Sen,
#'   \code{\link{tstsreg}} for S-type Theil-Sen with outlier removal
#'
#' @examples
#' \dontrun{
#' # Data with leverage points
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' # Add leverage points
#' x[1:3] <- x[1:3] + 5
#' y[1:3] <- y[1:3] - 5
#'
#' # Standard Theil-Sen (affected by leverage)
#' tsreg(x, y)
#'
#' # Hybrid approach (more robust)
#' tssnmreg(x, y)
#'
#' # Multiple predictors
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50)
#' x[1:5,] <- x[1:5,] + 4  # Leverage points
#' tssnmreg(x, y)
#' }
#'
#' @export
tssnmreg<-function(x,y,sc=pbvar,xout=FALSE,outfun=out,plotit=FALSE,...){
#
# Compute a modified Theil--Sen regression estimator.
# Use s-type initial estimate, eliminate points with
# outlying residuals, then do regular Theil--Sen
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
res=snmreg(x,y)$res
chk<-abs(res-median(res))/mad(res)
xx<-x[chk<=2,]
yy<-y[chk<=2]
temp<-tsreg(xx,yy)
list(coef=temp$coef,residuals=temp$res)
}

# ============================================================================
# gyreg
# ============================================================================

#' Robust Regression with Outlier Detection Based on Residual Scale
#'
#' Computes robust regression estimates by identifying and removing observations
#' with large scaled residuals using a data-adaptive threshold based on the
#' comparison of empirical and theoretical quantile functions.
#'
#' @param x A numeric vector or matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param rinit Initial robust regression function to use (default: \code{lmsreg}).
#'   Must return residuals in \code{rinit$res}.
#' @param K Threshold multiplier for flagging outliers based on scaled residuals
#'   (default: 2.5). Observations with \code{|residual|/MAD > K} are candidates
#'   for removal.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients from the refit (intercept, slopes).}
#'   \item{res}{Residuals from the final least squares fit after outlier removal.}
#' }
#'
#' @details
#' The \code{gyreg} function implements a two-stage robust regression procedure:
#'
#' **Stage 1 - Initial Fit and Outlier Detection:**
#' \enumerate{
#'   \item Fits an initial robust regression using \code{rinit} (default: least median of squares)
#'   \item Computes scaled residuals: \code{res.scale = |residual| / MAD(residuals)}
#'   \item Flags observations where \code{res.scale >= K}
#' }
#'
#' **Stage 2 - Adaptive Outlier Removal:**
#' \enumerate{
#'   \item For flagged observations, compares empirical quantiles to theoretical
#'         normal quantiles
#'   \item Computes the maximum deviation: \code{dval = max(pnorm(sorted.residuals[i]) - i/n)}
#'   \item Removes the top \code{floor(n * dval)} most extreme observations
#'   \item Refits using ordinary least squares on the retained observations
#' }
#'
#' This adaptive approach determines the number of outliers to remove based on
#' the data rather than using a fixed threshold. The method is particularly
#' effective when outliers cause substantial deviations from normality in the
#' residual distribution.
#'
#' The final fit uses ordinary least squares rather than a robust method,
#' which is appropriate after outlier removal but may be sensitive to any
#' remaining leverage points.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{lmsreg}} for least median of squares regression,
#'   \code{\link{opreg}} for outlier-pruned regression with different detection method,
#'   \code{\link{chreg}} for Coakley-Hettmansperger robust regression
#'
#' @examples
#' \dontrun{
#' # Simple regression with outliers
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50, sd=0.5)
#' y[1:3] <- y[1:3] + 5  # Add outliers
#'
#' # Standard robust regression
#' result <- gyreg(x, y)
#' print(result$coef)
#'
#' # Multiple predictors
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(50, sd=0.3)
#' y[1:5] <- y[1:5] + 4  # Add outliers
#' gyreg(x, y)
#'
#' # Using different initial estimator and threshold
#' gyreg(x, y, rinit=ltsreg, K=3.0)
#' }
#'
#' @export
gyreg<-function(x,y,rinit=lmsreg,K=2.5){
xy=elimna(cbind(x,y))
p=ncol(as.matrix(x))
p1=p+1
x=xy[,1:p]
y=xy[,p1]
res<-rinit(y~x)$res
res.scale<-abs(res)/mad(res)
flag<-(res.scale >=K)
i0<-sum(flag)
il<-length(y)-i0+1
res.sort<-sort(res.scale)
if(i0>0){
dval<-pnorm(res.sort[il:length(y)])-c(il:length(y))/length(y)
}
if(i0<=0)dval<-0
dval<-max(dval)
ndval<-floor(length(y)*dval)
if(ndval<0)ndval<-0
iup<-length(y)-ndval
rord<-order(res.scale)
flag<-rord[1:iup]
x=as.matrix(x)
temp<-lsfit(x[flag,],y[flag])
list(coef=temp$coef,res=temp$residual)
}

# ============================================================================
# regpreS
# ============================================================================
regpreS<-function(x,y,regfun=lsfit,error=absfun,nboot=100,
mval=round(5*log(length(y))),locfun=mean,pr=TRUE,
xout=FALSE,outfun=out,
plotit=TRUE,xlab="Model Number",ylab="Prediction Error",SEED=TRUE,...){
#
#   Stepwise selection of predictors based on
#   estimates of  prediction error using the regression method
#   regfun,
#   which defaults to least squares.  Prediction error
#   is estimated with .632 method.
#   (See Efron and Tibshirani, 1993, pp. 252--254)
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=100
#
#   Prediction error is the expected value of the function error.
#   The argument error defaults to absolute  error. To use
#   squared error, set error=sqfun.
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
#   The default value for mval, the number of observations to resample
#   for each of the B bootstrap samples is based on results by
#   Shao (JASA, 1996, 655-665). (Resampling n vectors of observations,
#   model selection may not lead to the correct model as n->infinity.
#
if(SEED)set.seed(2)
q=ncol(x)
qm1=q-1
x<-as.matrix(x)
d<-ncol(x)
p1<-d+1
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
y<-temp[,d+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,SEED=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
adit=NULL
pval=c(1:ncol(x))
#pval=c(1:q)
allp=pval
for(ip in 1:qm1){
model=list()
for(j  in 1:length(pval))model[[j]]=c(adit,pval[j])
temp=regpre(x,y,model=model,pr=FALSE,plotit=FALSE,adz=FALSE,regfun=regfun,
SEED=SEED)$estimates
pbest=order(temp[,5])
adit=model[[pbest[1]]]
pval=allp[-adit]
}
output=model[[pbest[1]]]
output=c(output,allp[-output])
output
}

# ============================================================================
# tsp1reg
# ============================================================================
#' Single Predictor Theil-Sen Regression
#'
#' Computes the Theil-Sen robust regression estimator for a single predictor.
#' This function is used internally by \code{\link{tsreg}} and \code{\link{tshdreg}}
#' but can also be called directly for simple linear regression.
#'
#' @param x A numeric vector containing the predictor variable.
#' @param y A numeric vector containing the response variable (must have same
#'   length as `x`).
#' @param plotit Logical. If `TRUE`, creates a scatterplot with the fitted
#'   regression line (default: `FALSE`).
#' @param HD Logical. If `TRUE`, uses Harrell-Davis estimator for computing
#'   the slope and intercept; if `FALSE`, uses sample median (default: `FALSE`).
#' @param OPT Logical. If `TRUE`, computes intercept as `median(y) - slope*median(X)`;
#'   if `FALSE`, computes intercept as `median(y - slope*X)` (default: `TRUE`).
#' @param tr Logical or numeric. Trimming parameter for Harrell-Davis estimator
#'   when `HD=TRUE` (default: `FALSE` for no trimming).
#'
#' @return A list with components:
#'   \item{coef}{Numeric vector of length 2 containing intercept and slope}
#'   \item{residuals}{Vector of residuals (observed - fitted values)}
#'
#' @details
#' For a single predictor, the Theil-Sen slope is computed as the median of all
#' pairwise slopes (y[j] - y[i]) / (x[j] - x[i]) for all pairs where x[j] > x[i].
#' When `HD=TRUE`, the Harrell-Davis estimator (a weighted average of order statistics)
#' is used instead of the sample median, which can provide improved efficiency.
#'
#' The intercept can be computed in two ways controlled by `OPT`:
#' \itemize{
#'   \item `OPT=TRUE`: intercept = median(y) - slope * median(x)
#'   \item `OPT=FALSE`: intercept = median(y - slope * x)
#' }
#'
#' Missing values are automatically removed before computation.
#'
#' @references
#' Theil, H. (1950). A rank-invariant method of linear and polynomial regression
#' analysis. \emph{Indagationes Mathematicae}, 12, 85-91.
#'
#' Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau.
#' \emph{Journal of the American Statistical Association}, 63, 1379-1389.
#'
#' @seealso \code{\link{tsreg}} for multiple predictor version,
#'   \code{\link{tshdreg}} for Harrell-Davis variant, \code{\link{hd}}
#'   for Harrell-Davis estimator
#'
#' @examples
#' # Simple linear regression
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50)
#' result <- tsp1reg(x, y)
#' result$coef  # Intercept and slope
#'
#' # With plot
#' tsp1reg(x, y, plotit=TRUE)
#'
#' # Using Harrell-Davis estimator
#' tsp1reg(x, y, HD=TRUE)
#'
#' # Different intercept computation
#' tsp1reg(x, y, OPT=FALSE)
#'
#' @export
tsp1reg<-function(x,y,plotit=FALSE,HD=FALSE,OPT=TRUE,tr=FALSE){
#
# Compute the Theil-Sen regression estimator.
# Only a single predictor is allowed in this version
#
#  OPT=TRUE, compute the intercept using median(y)-beta_1median(X)
#  OPT=FALSE compute the intercept using median of y-beta_1X
#
temp<-matrix(c(x,y),ncol=2)
temp<-elimna(temp)     # Remove any pairs with missing values
x<-temp[,1]
y<-temp[,2]
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
if(!HD)slope<-median(v1/v2,na.rm=TRUE)
if(HD)slope<-hd(v1/v2,na.rm=TRUE,tr=tr)
if(OPT){
if(!HD)coef<-median(y,na.rm=TRUE)-slope*median(x,na.rm=TRUE)
if(HD)coef<-hd(y,na.rm=TRUE)-slope*hd(x,na.rm=TRUE,tr=tr)
}
if(!OPT){
if(!HD)coef<-median(y-slope*x,na.rm=TRUE)
if(HD)coef<-hd(y-slope*x,na.rm=TRUE,tr=tr)
}
names(coef)<-"Intercept"
coef<-c(coef,slope)
if(plotit){
plot(x,y,xlab="X",ylab="Y")
abline(coef)
}
res<-y-slope*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# mdepreg.sub
# ============================================================================
mdepreg.sub<-function(X,theta){
np<-ncol(X)
p<-np-1
x<-X[,1:p]
y<-X[,np]
temp<-t(t(x)*theta[2:np])
yhat<-apply(temp,1,sum)+theta[1]
res<-y-yhat
val<-0-mregdepth(x,res)
val
}

# ============================================================================
# poireg
# ============================================================================
poireg<-function(x,y,xout=FALSE,outfun=outpro,plotit=FALSE,xlab="X",ylab="Y",
varfun=var,YHAT=FALSE,STAND=TRUE,...){
#
# Perform Poisson regression.
# The predictors are assumed to be stored in the n by p matrix x.
# The y values are typically count data (integers).
#
# xout=T will remove outliers from among the x values and then fit
# the regression line.
#  Default:
# One predictor, a mad-median rule is used.
# With more than one, projection method is used.
#
# outfun=out will use MVE method
#
xy=elimna(cbind(x,y))
x<-as.matrix(x)
x=xy[,1:ncol(x)]
y=xy[,ncol(xy)]
x<-as.matrix(x)
if(xout){
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
temp=glm(formula=y~x,family=poisson)
init=summary(temp)
yhat=temp$coef[1]
for(j in 1:ncol(x)){
j1=j+1
yhat=yhat+temp$coef[j1]*x[,j]
}
yhat=exp(yhat)
if(plotit){
x=as.matrix(x)
if(ncol(x)>1)stop("Cannot plot with more than one predictor")
plot(x,y,xlab=xlab,ylab=ylab)
#points(x,yhat,pch=".")
xord=order(x)
lines(x[xord],yhat[xord])
init$coef
}
ex=varfun(yhat)/varfun(y)
str=sqrt(ex)
hatv=NULL
if(YHAT)hatv=yhat
list(results=init,Explanatory.Power=ex,Strength.Assoc=str,yhat=hatv)
}

# ============================================================================
# tsreg
# ============================================================================
#' Theil-Sen Regression Estimator
#'
#' Computes the Theil-Sen robust regression estimator, a non-parametric alternative
#' to least squares regression that is resistant to outliers. For multiple predictors,
#' uses a Gauss-Seidel algorithm for estimation.
#'
#' @param x A numeric vector or matrix containing the predictor variable(s). For
#'   multiple predictors, rows represent observations and columns represent variables.
#' @inheritParams common-params
#' @param iter Number of iterations for the Gauss-Seidel algorithm when there are
#'   multiple predictors (default: 5). Values of 1 can be unsatisfactory; values
#'   of 10 may result in excessive execution time without improving Type I error rates.
#' @param varfun Function used to compute variance for strength of association
#'   (default: `pbvar`, percentage bend midvariance). Can be any function that
#'   computes a measure of variation.
#' @param tr Logical or numeric. If `TRUE` or numeric, uses Harrell-Davis estimator
#'   for the intercept (default: `FALSE` for median).
#' @param do.stre Logical. If `TRUE`, computes strength of association measure
#'   (default: `TRUE`). Set to `FALSE` when `varfun` cannot be computed.
#' @param corfun Correlation function used for computing explanatory power when
#'   the variance ratio exceeds 1 (default: `pbcor`).
#' @param WARN Logical. If `TRUE`, prints warning when measure of variation equals
#'   zero (default: `TRUE`).
#' @param HD Logical. If `TRUE`, uses Harrell-Davis estimator for the intercept;
#'   if `FALSE`, uses median (default: `FALSE`).
#' @param OPT Logical. If `TRUE`, computes intercept as `median(y) - b1*median(X)`;
#'   if `FALSE`, computes intercept as `median(y - b1*X)` (default: `FALSE`).
#'   The default is consistent with other R implementations of Theil-Sen.
#' @param xlab Label for x-axis when `plotit=TRUE` (default: "X").
#' @param ylab Label for y-axis when `plotit=TRUE` (default: "Y").
#'
#' @return A list with components:
#'   \item{n}{Original sample size before outlier removal}
#'   \item{n.keep}{Sample size after outlier removal (if `xout=TRUE`)}
#'   \item{coef}{Vector of regression coefficients (intercept followed by slopes)}
#'   \item{residuals}{Vector of residuals}
#'   \item{Strength.Assoc}{Strength of association (square root of explanatory power)}
#'   \item{Explanatory.Power}{Proportion of variance explained}
#'
#' @details
#' The Theil-Sen estimator computes the slope as the median of all pairwise slopes
#' between points. For a single predictor, this is computed directly. For multiple
#' predictors, a Gauss-Seidel algorithm is used for iterative refinement.
#'
#' The Theil-Sen estimator has a breakdown point of approximately 29%, meaning
#' it remains robust even when up to 29% of the data are outliers. This makes it
#' substantially more robust than ordinary least squares regression.
#'
#' Strength of association is computed as the square root of the ratio of the
#' variance of fitted values to the variance of y (using `varfun`). If this
#' ratio exceeds 1, the squared correlation between fitted and observed values
#' is used instead.
#'
#' @references
#' Theil, H. (1950). A rank-invariant method of linear and polynomial regression
#' analysis. \emph{Indagationes Mathematicae}, 12, 85-91.
#'
#' Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau.
#' \emph{Journal of the American Statistical Association}, 63, 1379-1389.
#'
#' Wilcox, R. R. (2022). \emph{Introduction to Robust Estimation and Hypothesis
#' Testing} (5th ed.). Academic Press.
#'
#' @seealso \code{\link{tshdreg}} for Theil-Sen with Harrell-Davis estimator,
#'   \code{\link{tsp1reg}} for single predictor version, \code{\link{tsregNW}}
#'   for neighborhood-weighted version, \code{\link{ltsreg}} for least trimmed
#'   squares, \code{\link{opreg}} for outlier-pruned regression
#'
#' @examples
#' # Simple regression with one predictor
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50)
#' result <- tsreg(x, y)
#' result$coef  # Intercept and slope
#'
#' # Multiple regression
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - 3*x[,2] + rnorm(50)
#' tsreg(x, y)
#'
#' # With outlier removal
#' x <- c(rnorm(45), rep(10, 5))  # 5 outliers
#' y <- 2 + 3*x + rnorm(50)
#' tsreg(x, y, xout=TRUE)  # Remove outliers before fitting
#'
#' @export
tsreg<-function(x,y,xout=FALSE,outfun=outpro,iter=5,varfun=pbvar,tr=FALSE,do.stre=TRUE,
corfun=pbcor,plotit=FALSE,WARN=TRUE,HD=FALSE,OPT=FALSE,xlab='X',ylab='Y',...){
#
#  Compute Theil-Sen regression estimator
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#  By default, the number of iterations is
#  iter=5
#  starting with Rallfun-v35. It was 10 but this can result in high execution time and
#  does not appear to be necessary in terms of Type I errors.
# iter=1 can be unsatisfactory
#
#  OPT=TRUE, compute the intercept using median(y)-b_1median(X)
#  OPT=FALSE compute the intercept using median of y-b_1X
#
#  Starting with version Rallfun-v29, OPT=F is the default, which is consistent with
#  other R functions that have been supplied for computing the Theil--Sen estimator.
#
#do.stre=FALSE, strength not computed, can be useful when varfun can't be computed
#
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
n=nrow(x)
n.keep=n
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else flag<-outfun(x,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
n.keep=nrow(x)
}
if(ncol(x)==1){
temp1<-tsp1reg(x,y,HD=HD,OPT=OPT,tr=tr)
coef<-temp1$coef
res<-temp1$res
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-tsp1reg(x[,p],y,OPT=OPT,tr=tr)$coef[2]
}
res<-y-x%*%temp
if(!HD)alpha<-median(res)
if(HD)alpha<-hd(res,tr=tr)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-tsp1reg(x[,p],r[,p],plotit=FALSE,OPT=OPT,tr=tr)$coef[2]
}
if(!HD)alpha<-median(y-x%*%temp)
if(HD)alpha<-hd(y-x%*%temp,tr=tr)
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
yhat<-y-res
stre=e.pow=NULL
if(do.stre){
temp=varfun(y)
if(temp==0){
if(WARN)print("Warning: When computing strength of association, measure of variation=0")
}
e.pow=NULL
if(temp>0){
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}}}
if(plotit){
if(ncol(x)==1){
plot(x,y,xlab=xlab,ylab=ylab)
abline(coef)
}}
list(n=n,n.keep=n.keep,
coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

# ============================================================================
# snmreg
# ============================================================================

#' Skipped Nelder-Mead Regression S-Estimator
#'
#' Computes a regression S-estimator using the Nelder-Mead optimization method.
#' The S-estimator minimizes a robust scale measure (percentage bend midvariance)
#' of residuals, providing high breakdown point protection.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param SEED Logical; if \code{TRUE} (default), sets random seed for reproducibility
#'   (currently not used, retained for compatibility).
#' @inheritParams opreg
#' @param initreg Initial regression function for starting values (default: \code{MMreg}).
#'   Must return coefficients in \code{$coef}.
#' @param ... Additional arguments passed to \code{outfun} if \code{xout=TRUE}.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients: \code{c(intercept, slopes)}.}
#'   \item{residuals}{Residuals from the fit.}
#' }
#'
#' @details
#' S-estimators minimize a robust scale of residuals:
#' \deqn{\min_{\beta} s(r_1, \ldots, r_n)}
#' where \eqn{r_i = y_i - x_i'\beta} are residuals and \eqn{s} is a robust
#' scale estimator (percentage bend midvariance).
#'
#' This implementation:
#' 1. Optionally removes outliers using projection-based detection
#' 2. Gets initial values from \code{initreg} (default: MM-regression)
#' 3. Uses Nelder-Mead simplex method to minimize scale of residuals
#' 4. Computes intercept as median of residuals
#'
#' **Properties**:
#' - High breakdown point (up to 50%)
#' - Robust to outliers in both x and y
#' - Affine equivariant
#' - Computationally intensive due to optimization
#'
#' The "skipped" aspect refers to optional outlier removal before fitting.
#'
#' @references
#' Rousseeuw, P.J. & Yohai, V.J. (1984). Robust regression by means of
#' S-estimators. \emph{Robust and Nonlinear Time Series Analysis}, 256-272.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @family robust regression functions
#' @seealso \code{\link{MMreg}} for MM-regression, \code{\link{ltsreg}} for LTS,
#'   \code{\link{outpro}} for outlier detection
#' @export
#' @examples
#' \dontrun{
#' # Simple S-regression
#' set.seed(123)
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2*x + rnorm(50, sd=0.5)
#' snmreg(x, y)
#'
#' # With outliers
#' y[1:5] <- y[1:5] + 10
#' snmreg(x, y)
#'
#' # Multiple regression with outlier removal
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50)
#' snmreg(x, y, xout=TRUE)
#' }
snmreg<-function(x,y,SEED=TRUE,xout=FALSE,outfun=outpro,initreg=MMreg,...){
#
# Compute regression S-estimator via Nelder-Mead method
# The measure of scale is taken to be the percentage bend midvariance
#
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
x <- as.matrix(x)
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
npm1=np-1
x=X[,1:npm1]
x=as.matrix(x)
y=X[,np]
N<-np-1
#temp<-initreg(x,y,SEED=SEED)$coef
temp<-initreg(x,y)$coef
START<-temp[2:np]
temp<-nelder(X,N,FN=snmreg.sub,START=START)
alpha <- median(y - x %*% temp)
coef <- c(alpha,temp)
res <- y - x %*% temp - alpha
list(coef = coef, residuals = res)
}

# ============================================================================
# mopreg
# ============================================================================

#' Multiple Outcome Outlier-Pruned Regression
#'
#' Performs regression for multiple outcome variables, optionally removing
#' multivariate outliers using projection-based outlier detection before fitting.
#' This extends outlier-pruned regression to handle multiple dependent variables
#' simultaneously.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric matrix of outcome variables (n by q), where q is the number
#'   of outcomes.
#' @param regfun Regression function to apply to each outcome (default: \code{tsreg}).
#'   Must return coefficients in \code{$coef}.
#' @param cop Type of correlation/covariance for outlier detection (default: 3).
#'   See \code{\link{outpro}} for options.
#' @param KEEP Logical; if \code{TRUE} (default), keeps all observations (no outlier
#'   removal). If \code{FALSE}, removes outliers before regression.
#' @param MC Logical; if \code{TRUE}, uses parallel processing via \code{\link{outproMC}}
#'   for outlier detection (default: \code{FALSE}).
#' @param STAND Logical; if \code{TRUE} (default), standardizes data before outlier
#'   detection.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Matrix of regression coefficients (p+1 by q). Each column contains
#'     coefficients for one outcome: \code{c(intercept, slopes)}.}
#'   \item{residuals}{Matrix of residuals (n by q). Each column contains residuals
#'     for one outcome.}
#' }
#'
#' @details
#' This function extends outlier-pruned regression to multiple outcomes:
#'
#' 1. **If \code{KEEP=FALSE}**: Detects multivariate outliers in the combined
#'    (X, Y) space using projection-based methods. Outliers are removed before
#'    any regression.
#'
#' 2. **For each outcome**: Applies \code{regfun} to predict that outcome from X,
#'    using only the retained observations.
#'
#' The outlier detection considers all variables (predictors and all outcomes)
#' jointly, which can be more powerful than detecting outliers separately for
#' each outcome.
#'
#' **Use cases**:
#' - Multivariate regression with multiple correlated outcomes
#' - Protecting all analyses from the same set of outlying cases
#' - Ensuring consistency across multiple outcome models
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @family robust regression functions
#' @seealso \code{\link{opreg}} for single outcome version,
#'   \code{\link{outpro}} for outlier detection details,
#'   \code{\link{tsreg}} for Theil-Sen regression
#' @export
#' @examples
#' \dontrun{
#' # Multiple outcomes
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y1 <- x[,1] + 2*x[,2] + rnorm(50)
#' y2 <- 2*x[,1] - x[,2] + rnorm(50)
#' y <- cbind(y1, y2)
#'
#' # Without outlier removal
#' result1 <- mopreg(x, y, KEEP=TRUE)
#'
#' # With outlier removal
#' y[1:3,] <- y[1:3,] + 10  # Add outliers
#' result2 <- mopreg(x, y, KEEP=FALSE)
#' }
mopreg<-function(x,y,regfun=tsreg,cop=3,KEEP=TRUE,MC=FALSE,STAND=TRUE){
#
# Do multiple (outcomes) regression on points not labled outliers
# using projection-type outlier detection method
# Arg=regfun determines regression method;
# by default, Theil-Sen is used.
#
#  KEEP=F, outliers will be eliminated
#  KEEP=T, outliers are not eliminated
# cop: see function outpro
x<-as.matrix(x)
y<-as.matrix(y)
px<-ncol(x)
py<-ncol(y)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
if(KEEP)ivec<-c(1:nrow(x))
if(!KEEP){
if(!MC)ivec<-outpro(m,plotit=FALSE,cop=cop,STAND=STAND)$keep
if(MC)ivec<-outproMC(m,plotit=FALSE,cop=cop,STAND=STAND)$keep
}
np1<-ncol(x)+1
vec<-rep(1,nrow(m))
pxpy<-px+py
coef<-matrix(ncol=py,nrow=np1)
res<-matrix(ncol=py,nrow=nrow(m))
for(i in 1:py){
pv<-px+i
coef[,i]<-regfun(m[ivec,1:ncol(x)],m[ivec,pv])$coef
vec<-as.matrix(vec)
res[,i]<-m[,pv]-cbind(vec,m[,1:ncol(x)])%*%coef[,i]
}
list(coef=coef,residuals=res)
}

# ============================================================================
# mdepreg.coef
# ============================================================================
mdepreg.coef<-function(x,y,xout=FALSE,outfun=out,...){
#
# multiple depth regression
#
X<-cbind(x,y)
X<-elimna(X)
p1=ncol(X)
p=p1-1
if(xout){
flag<-outfun(X[,1:p],plotit=FALSE,...)$keep
X<-X[flag,]
}
library(mrfDepth)
a=rdepthmedian(X)$deep
list(coef=a)
}

# ============================================================================
# mdepreg
# ============================================================================

#' Multivariate Depth-Based Regression
#'
#' Computes regression coefficients using multivariate depth methods via the
#' \code{mrfDepth} package. This extends depth-based regression to handle
#' multiple predictors by finding the deepest regression hyperplane.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @inheritParams opreg
#' @param RES Logical; if \code{TRUE}, computes and returns residuals
#'   (default: \code{FALSE}).
#' @param ... Additional arguments passed to \code{outfun} if \code{xout=TRUE}.
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Original sample size before any outlier removal.}
#'   \item{n.keep}{Number of observations retained after outlier removal
#'     (same as \code{n} if \code{xout=FALSE}).}
#'   \item{coef}{Vector of regression coefficients: \code{c(intercept, slopes)}.}
#'   \item{residuals}{Residuals (if \code{RES=TRUE}, otherwise \code{NA}).}
#' }
#'
#' @details
#' This function uses regression depth via the \code{mrfDepth} package to find
#' the "deepest" regression hyperplane in the (X, y) space.
#'
#' The method:
#' 1. Combines predictors and outcome: \eqn{(X_1, \ldots, X_p, y)}
#' 2. Computes the regression depth median using \code{rdepthmedian()}
#' 3. Returns the deepest hyperplane coefficients
#'
#' **Regression depth** measures how "central" a hyperplane is relative to the
#' data cloud. The deepest hyperplane is most representative of the overall
#' data pattern and is highly resistant to outliers.
#'
#' **Properties**:
#' - High breakdown point (up to 33%)
#' - Affine equivariant
#' - Robust to outliers in both predictors and outcome
#' - No distributional assumptions
#'
#' **Note**: Requires the \code{mrfDepth} package. For simple regression
#' (one predictor), \code{\link{depreg}} may be faster.
#'
#' @references
#' Rousseeuw, P.J. & Hubert, M. (1999). Regression depth. \emph{Journal of the
#' American Statistical Association}, 94, 388-433.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @family robust regression functions
#' @seealso \code{\link{depreg}} for simple regression, \code{\link{tsreg}}
#'   for Theil-Sen regression
#' @export
#' @examples
#' \dontrun{
#' # Requires mrfDepth package
#' if (require(mrfDepth)) {
#'   # Multiple regression
#'   set.seed(123)
#'   x <- matrix(rnorm(100), ncol=2)
#'   y <- x[,1] + 2*x[,2] + rnorm(50)
#'   mdepreg(x, y, RES=TRUE)
#'
#'   # With outliers
#'   y[1:5] <- y[1:5] + 15
#'   mdepreg(x, y, RES=TRUE)
#' }
#' }
mdepreg<-function(x,y,xout=FALSE,outfun=out,RES=FALSE,...){
#
# multiple depth regression
#
X<-cbind(x,y)
X<-elimna(X)
n=nrow(X)
p1=ncol(X)
p=p1-1
if(xout){
flag<-outfun(X[,1:p],plotit=FALSE,...)$keep
X<-X[flag,]
}
n.keep=nrow(X)
library(mrfDepth)
a=rdepthmedian(X)$deepest
res=NA
if(RES)res=X[,p1]-X[,1:p]%*%a[2:p1]-a[1]
list(n=n,n.keep=n.keep,coef=a,residuals=res)
}

# ============================================================================
# mdepreg.orig
# ============================================================================
mdepreg.orig<-function(x,y,xout=FALSE,outfun=outpro){
#
# multiple depth regression
#
X<-cbind(x,y)
X<-elimna(X)
np=n.keep=ncol(X)
p=np-1
if(xout){
id=outfun(X[,1:p],plotit=FALSE)$keep
X=X[id,]
n.keep=nrow(X)
}
if(np==2){
temp=depreg(X[,1],X[,2])
coef=temp$coef
res=temp$residuals
}
if(np>2){
N<-np-1
x=X[,1:N]
y=X[,np]
START<-tsreg(x,y)$coef
coef<-nelderv2(X,np,FN=mdepreg.sub,START=START)
x <- as.matrix(x)
res <- y - x %*% coef[2:np] - coef[1]
}
list(n=n,n.keep=n.keep,coef = coef, residuals = res)
}

# ============================================================================
# MMreg
# ============================================================================

#' MM-Regression Estimator
#'
#' Computes Yohai's (1987) MM-regression estimator, which combines high breakdown
#' point with high efficiency. This is a wrapper around \code{robustbase::lmrob}
#' with optional outlier removal and association strength assessment.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param RES Logical; if TRUE (default), returns residuals.
#' @param STAND Logical; if TRUE (default), standardizes predictors for outlier detection.
#' @param varfun Function for computing variance (default: \code{pbvar}), used for
#'   strength of association.
#' @param corfun Function for computing correlation (default: \code{pbcor}), used as
#'   fallback for strength of association.
#' @param WARN Logical; if FALSE (default), suppresses warnings from \code{lmrob}.
#' @inheritParams opreg
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#'   \item{residuals}{Residuals from the fit (if \code{RES=TRUE}).}
#'   \item{Strength.Assoc}{Measure of association strength (explanatory power).}
#' }
#'
#' @details
#' MM-regression achieves:
#' - **High breakdown point** (50%) from the initial S-estimate
#' - **High efficiency** (95% under normality) from the final M-estimate
#' - **Equivariance** under affine transformations
#'
#' The method works in three stages:
#' 1. Compute a high-breakdown S-estimate of regression
#' 2. Compute M-estimate of scale using S-estimate residuals
#' 3. Compute final M-estimate of regression using this scale
#'
#' **Association Strength**: The function computes an explanatory power measure
#' as the ratio of variance of fitted values to variance of observed values,
#' using robust variance estimators. If this exceeds 1 (possible with robust methods),
#' it falls back to the squared robust correlation.
#'
#' **Requires**: The \code{robustbase} package must be installed.
#'
#' @references
#' Yohai, V.J. (1987). High breakdown-point and high efficiency robust estimates
#' for regression. The Annals of Statistics, 15, 642-656.
#'
#' @family robust regression functions
#' @seealso \code{\link{ltsreg}}, \code{\link{LMSreg}}
#' @export
#' @examples
#' \dontrun{
#' # Simple MM-regression
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2*x + rnorm(50, sd=0.5)
#' MMreg(x, y)
#'
#' # Multiple regression with outlier removal
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50, sd=0.3)
#' y[1:5] <- y[1:5] + 10  # add outliers
#' fit <- MMreg(x, y, xout=TRUE)
#' print(fit$Strength.Assoc)
#' }
MMreg<-function(x,y,RES=TRUE,xout=FALSE,outfun=outpro,STAND=TRUE,varfun=pbvar,corfun=pbcor,WARN=FALSE,...){
#
#  Compute MM regression estimate derived by Yohai (1987)
#  simply by calling the R function lmrob
#  This function will remove leverage points when
#  xout=T
#  using the outlier detection method indicated by
#  outfun, which defaults to the projection method.
#
library('robustbase')
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else flag<-outpro(x,STAND=STAND,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(!WARN)options(warn=-1)
temp=lmrob(y~x)
if(!WARN)options(warn=0)
coef=temp$coefficients
p1=ncol(x)+1
res<-y-x%*%coef[2:p1]-coef[1]
yhat<-y-res
stre=NULL
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}
if(!RES)res=NULL
list(coef=coef,residuals=res,Strength.Assoc=stre)
}

# ============================================================================
# opregpbMC
# ============================================================================
opregpbMC<-function(x,y,nboot=1000,alpha=.05,om=TRUE,ADJ=TRUE,cop=3,SEED=TRUE,
nullvec=rep(0,ncol(x)+1),plotit=TRUE,opdis=2,gval=sqrt(qchisq(.95,ncol(x)+1))){
#
#  Same as opregpb, only this function takes advantage of a multi-core
#  processor assuming one is availabe and that the R package
#  multicore has been installed.
#
# generate bootstrap estimates
# use projection-type outlier detection method followed by
# TS regression.
#
# om=T and ncol(x)>1, means an omnibus test is performed,
# otherwise only individual tests of parameters are performed.
#
# opdis=2, means that Mahalanobis distance is used
# opdis=1, means projection-type distance is used
#
# gval is critical value for projection-type outlier detection
# method
#
# ADJ=T, Adjust p-values as described in Section 11.1.5 of the text.
#
if(SEED)set.seed(2)
x<-as.matrix(x)
m<-cbind(x,y)
p1<-ncol(x)+1
m<-elimna(m) # eliminate any rows with missing data
x<-m[,1:ncol(x)]
x<-as.matrix(x)
y<-m[,p1]
if(nrow(x)!=length(y))stop("Sample size of x differs from sample size of y")
if(!is.matrix(x))stop("Data should be stored in a matrix")
print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun=opregMC)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
# using Hochberg method
bvec<-t(bvec)
dvec<-alpha/(c(1:ncol(x)))
test<-NA
icl0<-round(alpha*nboot/2)
icl<-round(alpha*nboot/(2*ncol(x)))
icu0<-nboot-icl0
icu<-nboot-icl
output<-matrix(0,p1,6)
dimnames(output)<-list(NULL,c("Param.","p.value","crit.p.value",
"ci.lower","ci.upper","s.e."))
pval<-NA
for(i in 1:p1){
output[i,1]<-i-1
se.val<-var(bvec[,i])
temp<-sort(bvec[,i])
output[i,6]<-sqrt(se.val)
if(i==1){
output[i,4]<-temp[icl0+1]
output[i,5]<-temp[icu0]
}
if(i>1){
output[i,4]<-temp[icl+1]
output[i,5]<-temp[icu]
}
pval[i]<-sum((temp>nullvec[i]))/length(temp)
if(pval[i]>.5)pval[i]<-1-pval[i]
}
fac<-2
if(ADJ){
# Adjust p-value if n<60
nval<-length(y)
if(nval<20)nval<-20
if(nval>60)nval<-60
fac<-2-(60-nval)/40
}
pval[1]<-2*pval[1]
pval[2:p1]<-fac*pval[2:p1]
output[,2]<-pval
temp2<-order(0-pval[2:p1])
zvec<-dvec[1:ncol(x)]
sigvec<-(test[temp2]>=zvec)
output[temp2+1,3]<-zvec
output[1,3]<-NA
output[,2]<-pval
om.pval<-NA
temp<-opregMC(x,y)$coef
if(om && ncol(x)>1){
temp2<-rbind(bvec[,2:p1],nullvec[2:p1])
if(opdis==1)dis<-pdisMC(temp2,center=temp[2:p1])
if(opdis==2){
cmat<-var(bvec[,2:p1]-apply(bvec[,2:p1],2,mean)+temp[2:p1])
dis<-mahalanobis(temp2,temp[2:p1],cmat)
}
om.pval<-sum((dis[nboot+1]<=dis[1:nboot]))/nboot
}
# do adjusted p-value
nval<-length(y)
if(nval<20)nval<-20
if(nval>60)nval<-60
adj.pval<-om.pval/2+(om.pval-om.pval/2)*(nval-20)/40
if(ncol(x)==2 && plotit){
plot(bvec[,2],bvec[,3],xlab="Slope 1",ylab="Slope 2")
temp.dis<-order(dis[1:nboot])
ic<-round((1-alpha)*nboot)
xx<-bvec[temp.dis[1:ic],2:3]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(output=output,om.pval=om.pval,adj.om.pval=adj.pval)
}

# ============================================================================
# opregMC
# ============================================================================

#' Outlier-Pruned Regression (Parallel Processing)
#'
#' Performs regression after removing outliers detected by the projection-based
#' outlier detection method \code{\link{outproMC}}. Uses parallel processing for
#' outlier detection.
#'
#' @param x Numeric vector or matrix of predictor variable(s).
#' @param y Numeric vector of the dependent variable.
#' @param regfun Regression function to use (default: \code{tsreg}).
#'   Can be any function that returns coefficients in \code{$coef}.
#' @param cop Critical outlier detection parameter passed to \code{outproMC}.
#'   Larger values make detection less sensitive (default: 3).
#' @param fast Logical. Currently not used (default: \code{FALSE}).
#' @param pr Logical. If \code{TRUE} and \code{fast=TRUE}, prints coefficients
#'   (default: \code{TRUE}).
#' @param prres Logical. If \code{TRUE} and \code{fast=TRUE}, prints residuals
#'   (default: \code{FALSE}).
#' @param STAND Logical. If \code{TRUE}, standardizes data before outlier
#'   detection (default: \code{TRUE}).
#' @param xout Logical. Included for compatibility but not used; outlier
#'   detection is always performed (default: \code{FALSE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept, slopes)}
#'   \item{residuals}{Vector of residuals for all observations (including outliers)}
#' }
#'
#' @details
#' This function first uses \code{\link{outproMC}} (projection-based outlier
#' detection with multicore processing) to identify outliers in the joint space
#' of predictors and response. Then it fits the regression using only the
#' non-outlier points.
#'
#' The residuals are computed for all observations, including those flagged as
#' outliers, so you can see how well the robust fit predicts the outlying cases.
#'
#' **Outlier Detection**: Uses projection pursuit methods to detect outliers
#' in high-dimensional space. The \code{cop} parameter controls sensitivity:
#' - Larger \code{cop}: fewer points flagged as outliers (more conservative)
#' - Smaller \code{cop}: more points flagged as outliers (more aggressive)
#' - Default (3): reasonable balance for most applications
#'
#' **When to use**:
#' - High-dimensional regression with potential outliers
#' - When you want automatic outlier removal
#' - As a more robust alternative to \code{\link{opreg}}
#'
#' This is the parallel processing version, faster than \code{\link{opreg}}
#' for large datasets.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{opreg}} for non-parallel version,
#'   \code{\link{outproMC}} for the outlier detection method,
#'   \code{\link{opregpb}} for bootstrap inference
#'
#' @examples
#' \dontrun{
#' # Regression with outliers
#' set.seed(123)
#' x <- matrix(rnorm(200), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(100)
#'
#' # Add some outliers
#' y[1:5] <- y[1:5] + 10
#' x[1:5,] <- x[1:5,] + 3
#'
#' # Fit with automatic outlier removal
#' fit <- opregMC(x, y)
#' print(fit$coef)
#'
#' # Compare with regular tsreg (influenced by outliers)
#' fit.regular <- tsreg(x, y)
#' print(fit.regular$coef)
#' }
#'
#' @export
opregMC<-function(x,y,regfun=tsreg,cop=3,fast=FALSE,pr=TRUE,prres=FALSE,STAND=TRUE,xout=FALSE){
#
# Do regression on points not labled outliers
# using projection-type outlier detection method
#
#  Note: argument xout is not relevant here, but is included to avoid conflicts when using regci.
#
x<-as.matrix(x)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
ivec<-outproMC(m,plotit=FALSE,cop=cop,fast=FALSE,pr=FALSE,STAND=STAND)$keep
np1<-ncol(x)+1
coef<-regfun(m[ivec,1:ncol(x)],m[ivec,np1])$coef
vec<-rep(1,length(y))
residuals<-y-cbind(vec,x)%*%coef
if(fast && pr){
print("Intercept, followed by slopes:")
print(coef)
if(prres){
print("Residuals:")
print(residuals)
}}
list(coef=coef,residuals=residuals)
}

# ============================================================================
# bkreg
# ============================================================================

#' Binary Kernel Regression
#'
#' Computes a kernel-based estimator for binary (0/1) response regression.
#' Estimates P(Y=1|X=x) nonparametrically using kernel density estimation.
#'
#' @param x Numeric vector or matrix of predictor variable(s).
#' @param y Numeric vector of binary responses (should contain only 0 and 1).
#' @param kerfun Kernel density estimation function (default: \code{akerd}).
#'   Must support \code{pyhat=TRUE} and \code{pts} arguments.
#' @param pyhat Logical. If \code{TRUE}, returns estimated probabilities;
#'   if \code{FALSE}, returns "Done" message (default: \code{FALSE}).
#' @inheritParams common-params
#' @param plotit Logical. If \code{TRUE}, creates a plot of the fitted
#'   probabilities (default: \code{TRUE}).
#' @param xlab Label for x-axis in plots (default: "X").
#' @param ylab Label for y-axis in plots (default: "Y").
#' @param zlab Label for z-axis in 3D plots (default: "Z").
#' @param pr Logical. If \code{TRUE}, prints warnings about scaling for
#'   2-predictor case (default: \code{TRUE}).
#' @param theta Viewing angle for 3D plot (default: 50).
#' @param phi Colatitude angle for 3D plot (default: 25).
#' @param duplicate How to handle duplicate points in \code{interp}
#'   (default: "error"). See \code{\link[akima]{interp}}.
#' @param expand Expansion factor for 3D plot (default: 0.5).
#' @param scale Logical. If \code{TRUE}, scales z-axis in 3D plot (default: \code{FALSE}).
#' @param ticktype Type of axis ticks for 3D plot (default: "simple").
#'
#' @return If \code{pyhat=TRUE}, returns a numeric vector of estimated probabilities
#'   P(Y=1|X=x) for each observation. If \code{pyhat=FALSE}, returns "Done".
#'
#' @details
#' This function implements a nonparametric estimator for binary regression
#' based on the method of Signorini and Jones (2004). The estimator uses
#' kernel density estimation to estimate the conditional probability:
#'
#' \deqn{P(Y=1|X=x) = \frac{n_1 f(x|Y=1)}{n_1 f(x|Y=1) + n_0 f(x|Y=0)}}
#'
#' where \eqn{f(x|Y=1)} and \eqn{f(x|Y=0)} are kernel density estimates
#' for the predictor distributions in each response group, and \eqn{n_1},
#' \eqn{n_0} are the sample sizes in each group.
#'
#' **Visualization**:
#' - For 1 predictor: creates a 2D plot showing data points and fitted curve
#' - For 2 predictors: creates a 3D surface plot (requires \code{akima} package)
#'
#' **When to use**:
#' - Binary response variable (logistic-type problems)
#' - Want nonparametric alternative to logistic regression
#' - Relationship between predictors and log-odds may be nonlinear
#' - Small to moderate sample sizes (kernel methods can be unstable with sparse data)
#'
#' @references
#' Signorini, D.F., & Jones, M.C. (2004). Kernel estimators for univariate
#' binary regression. \emph{Journal of the American Statistical Association},
#' 99, 119-126.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{akerd}} for the default kernel density estimator,
#'   \code{\link{kerreg}} for continuous response kernel regression
#'
#' @examples
#' \dontrun{
#' # Simulated binary response data
#' set.seed(123)
#' x <- rnorm(100)
#' p <- plogis(1 + 2*x)  # True probabilities
#' y <- rbinom(100, 1, p)  # Binary outcomes
#'
#' # Fit and visualize
#' probs <- bkreg(x, y, pyhat=TRUE, plotit=TRUE)
#'
#' # Two predictors with 3D plot
#' x <- matrix(rnorm(200), ncol=2)
#' p <- plogis(1 + x[,1] - 0.5*x[,2])
#' y <- rbinom(100, 1, p)
#' bkreg(x, y, pyhat=FALSE, plotit=TRUE, scale=TRUE)
#' }
#'
#' @export
bkreg<-function(x,y,kerfun=akerd,pyhat=FALSE,plotit=TRUE,xlab="X",ylab="Y",
zlab="Z",xout=FALSE,outfun=outpro,pr=TRUE,theta=50,phi=25,duplicate="error",
expand=.5,scale=FALSE,ticktype="simple",...){
#
# Kernel estimator for binary regression.
# (See Signorini and Jones, JASA, 2004, 119-)
#
x=as.matrix(x)
p=ncol(x)
p1=p+1
xx<-elimna(cbind(x,y))
x<-xx[,1:p]
y<-xx[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x=as.matrix(x)
flag<-(y==1)
mv=sum(flag)
nv=sum(!flag)
phat<-NA
fhat<-kerfun(x[flag,],pyhat=TRUE,plotit=FALSE,pts=x)
ghat<-kerfun(x[!flag,],pyhat=TRUE,plotit=FALSE,pts=x)
phat<-mv*fhat/(mv*fhat+nv*ghat)
if(p==1){
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab)
flag2<-order(x)
#lines(x[flag2],phat[flag2])
lines(x[flag2],phat)
}}
if(p==2){
if(plotit){
if(pr){
if(!scale)print("With dependence, suggest using scale=T")
}
fitr<-phat
iout<-c(1:length(fitr))
nm1<-length(fitr)-1
for(i in 1:nm1){
ip1<-i+1
for(k in ip1:length(fitr))if(sum(x[i,]==x[k,])==2)iout[k]<-0
}
fitr<-fitr[iout>=1] # Eliminate duplicate points in the x-y plane
#                 This is necessary when doing three dimensional plots
#                 with the R function interp
mkeep<-x[iout>=1,]
fit<-interp(mkeep[,1],mkeep[,2],fitr,duplicate=duplicate)
persp(fit,theta=theta,phi=phi,expand=expand,
scale=scale,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype)
}}
if(!pyhat)phat<-"Done"
phat
}

# ============================================================================
# regpreCV
# ============================================================================
regpreCV<-function(x,y,regfun=tsreg,varfun=pbvar,adz=TRUE,model=NULL,locfun=mean,
xout=FALSE,outfun=out,
plotit=TRUE,xlab="Model Number",ylab="Prediction Error",...){
#
# Estimate the prediction error using the regression method
#   regfun in conjunction with leave-one-out cross-validation
#
#   The argument model should have list mode, model[[1]] indicates
#   which predictors are used in the first model. For example, storing
#   1,4 in model[[1]] means predictors 1 and 4 are being considered.
#   If model is not specified, and number of predictors is at most 5,
#   then all models are considered.
#
#   If adz=T, added to the models to be considered is where
#   all regression slopes are zero. That is, use measure of location only
#   corresponding to
#   locfun.
#
x<-as.matrix(x)
d<-ncol(x)
p1<-d+1
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
y<-temp[,d+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(is.null(model)){
if(d<=5)model<-modgen(d,adz=adz)
if(d>5)model[[1]]<-c(1:ncol(x))
}
mout<-matrix(NA,length(model),3,dimnames=list(NULL,c("est.error",
"var.used","rank")))
for (imod in 1:length(model)){
nmod=length(model[[imod]])-1
temp=c(nmod:0)
mout[imod,2]=sum(model[[imod]]*10^temp)
#
if(sum(model[[imod]]==0)!=1){
xx<-x[,model[[imod]]]
xx<-as.matrix(xx)
mout[imod,1]<-regpecv(xx,y,regfun=regfun,varfun=varfun,...)
}
#
if(sum(model[[imod]]==0)==1){
mout[imod,1]<-locCV(y,varfun=varfun,locfun=locfun)
}}
mout[,3]=rank(mout[,1])
if(plotit)plot(c(1:nrow(mout)),mout[,1],xlab=xlab,ylab=ylab)
mout
}

# ============================================================================
# regci
# ============================================================================

#' Bootstrap Confidence Intervals for Regression Parameters
#'
#' Computes bootstrap confidence intervals for each parameter (intercept and slopes)
#' in a linear regression equation using the percentile bootstrap method.
#'
#' @inheritParams common_params
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#'   Any function that returns coefficients in \code{$coef} (intercept first, then slopes).
#' @param nboot Number of bootstrap samples (default: 599).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param pr Logical; if TRUE, prints progress messages (default: TRUE).
#' @param null.val Vector of null hypothesis values for each parameter (default: all zeros).
#'   Length should equal number of parameters (p+1).
#' @param method P-value adjustment method for slopes (default: 'hoch' for Hochberg).
#'   Passed to \code{p.adjust()}.
#' @param plotit Logical; if TRUE and there are exactly 2 predictors, plots the
#'   bootstrap confidence region (default: FALSE).
#' @param xlab X-axis label for plot (default: "Predictor 1").
#' @param ylab Y-axis label for plot (default: "Predictor 2").
#' @param WARNS Logical; if FALSE, suppresses warnings during bootstrap (default: FALSE).
#' @param LABELS Logical; if TRUE, uses variable names from x as labels (default: FALSE).
#' @param ... Additional arguments passed to \code{regfun} and \code{outfun}.
#'
#' @return A list with components:
#' \item{regci}{Matrix with rows for each parameter and columns:
#'   \code{ci.low} (lower CI bound), \code{ci.up} (upper CI bound),
#'   \code{Estimate} (point estimate), \code{S.E.} (bootstrap standard error),
#'   \code{p-value} (two-sided p-value), \code{p.adj} (adjusted p-values for slopes).}
#' \item{n}{Sample size after removing missing values.}
#' \item{n.keep}{Sample size after removing leverage points (if \code{xout=TRUE}).}
#'
#' @details
#' The function uses percentile bootstrap to construct confidence intervals for
#' regression parameters. P-values test whether each parameter differs from the
#' corresponding value in \code{null.val} (default 0). For slopes, adjusted p-values
#' control family-wise error rate using the specified method.
#'
#' When using least squares regression with n < 250, consider using \code{lsfitci}
#' instead for better performance.
#'
#' If duplicate y-values are detected with Theil-Sen regression, the function suggests
#' using \code{tshdreg} which may have better power.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' @seealso \code{\link{tsreg}}, \code{\link{lsfitci}}, \code{\link{regtest}},
#' \code{\link{reg2ci}}
#'
#' @examples
#' # Theil-Sen regression with bootstrap CIs
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2 + 3*x + rnorm(50)
#' regci(x, y)
#'
#' # Multiple predictors with leverage point removal
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(50)
#' regci(x, y, xout=TRUE, nboot=500)
#'
#' @export
regci<-function(x,y,regfun=tsreg,nboot=599,alpha=.05,SEED=TRUE,pr=TRUE,null.val=NULL,method='hoch',
xout=FALSE,outfun=outpro,plotit=FALSE,xlab="Predictor 1",ylab="Predictor 2",WARNS=FALSE,LABELS=FALSE,...){
#
#   Compute a .95 confidence interval for each of the parameters of
#   a linear regression equation. The default regression method is
#   the Theil-Sen estimator.
#
#   When using the least squares estimator, and when n<250, use
#   lsfitci instead.
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=599
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimated of
#   the first predictor, etc.
#
#   plotit=TRUE: If there are two predictors, plot 1-alpha confidence region based
#  on the bootstrap samples.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
nrem=length(y)
if(xout){
if(pr)print("Default for outfun is now outpro, not out")
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
estit=regfun(x,y,...)$coef
if(is.null(null.val))null.val=rep(0,p1)
flagF=FALSE
flagF=identical(regfun,tsreg)
if(flagF){if(pr){
if(sum(duplicated(y)>0))print("Duplicate values detected; tshdreg might have more power than tsreg")
}}
nv=length(y)
x<-as.matrix(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
if(!WARNS)options(warn=-1)
bvec<-apply(data,1,regboot,x,y,regfun,xout=FALSE,...)
options(warn=0)
#Leverage points already removed.
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
regci<-matrix(0,p1,6)
vlabs="Intercept"
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
if(LABELS)vlabs[2:p1]=labels(x)[[2]]
dimnames(regci)<-list(vlabs,c("ci.low","ci.up","Estimate","S.E.","p-value",'p.adj'))
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
se<-NA
pvec<-NA
for(i in 1:p1){
bsort<-sort(bvec[i,])
#pvec[i]<-(sum(bvec[i,]<0)+.5*sum(bvec[i,]==0))/nboot
pvec[i]<-(sum(bvec[i,]<null.val[i])+.5*sum(bvec[i,]==null.val[i]))/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
regci[i,1]<-bsort[ilow]
regci[i,2]<-bsort[ihi]
se[i]<-sqrt(var(bvec[i,]))
}
if(p1==3){
if(plotit){
plot(bvec[2,],bvec[3,],xlab=xlab,ylab=ylab)
}}
regci[,3]=estit
pvec<-2*pvec
regci[,4]=se
regci[,5]=regci[,6]=pvec
regci[2:p1,6]=p.adjust(pvec[2:p1],method=method)
list(regci=regci,n=nrem,n.keep=nv)
}

# ============================================================================
# regpecv
# ============================================================================
regpecv<-function(x,y,regfun=tsreg,varfun=pbvar,...){
#
# Estimate prediction error via leave-one-out cross-validation
#
# regfun defaults to Theil-Sen estimator
# function returns measure of prediction error: robust measure of variation
# applied to the n differences y_i-y_{-i}, i=1,...,n
# where y_{-1} is estimate of y when ith vector of observations is omitted.
#
xy=elimna(cbind(x,y))
x=as.matrix(x)
px=ncol(x)
px1=px+1
n=nrow(xy)
vals=NA
for(i in 1:n){
est=regfun(xy[-i,1:px],xy[-i,px1])$coef
vals[i]=xy[i,px1]-(est[1]+sum(est[2:px1]*xy[i,1:px]))
}
pe=varfun(vals)
pe
}

# ============================================================================
# reg1way
# ============================================================================

#' One-Way ANOVA for Independent Regression Lines
#'
#' Tests the hypothesis that regression parameters (intercepts and slopes) are equal
#' across two or more independent groups using a bootstrap-based Johansen-type MANOVA
#' approach.
#'
#' @param x List of predictor matrices, one per group. \code{x[[j]]} contains the
#'   n_j x p predictor matrix for group j.
#' @param y List of response vectors, one per group. \code{y[[j]]} contains the
#'   response values for group j.
#' @inheritParams common_params
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#' @param nboot Number of bootstrap samples per group (default: 100).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param STAND Logical; if TRUE, standardizes predictors before outlier detection (default: TRUE).
#' @param AD Logical; if TRUE, computes adjusted critical value and p-value using
#'   Johansen correction (default: FALSE).
#' @param pr Logical; if TRUE, prints progress messages and warnings (default: TRUE).
#' @param ... Additional arguments passed to \code{regfun} and \code{outfun}.
#'
#' @return A list with components:
#' \item{n}{Vector of sample sizes per group after removing missing values.}
#' \item{n.keep}{Vector of sample sizes per group after removing leverage points.}
#' \item{test.stat}{Johansen-type test statistic (follows chi-squared distribution).}
#' \item{crit.value}{Critical value at level \code{alpha}.}
#' \item{adjusted.crit}{Johansen-adjusted critical value (if \code{!xout} or \code{AD=TRUE}).}
#' \item{p.value}{Unadjusted p-value.}
#' \item{adjusted.p.value}{Johansen-adjusted p-value (if \code{!xout} or \code{AD=TRUE}).}
#' \item{est}{Data frame of estimated regression parameters for each group.}
#'
#' @details
#' This function performs a one-way ANOVA comparing regression lines across J independent
#' groups. It tests whether all intercepts and slopes are equal across groups. The test
#' uses bootstrap estimates of standard errors combined with a Johansen MANOVA-type
#' test statistic.
#'
#' If \code{xout=TRUE}, leverage points are removed using \code{outfun} (default: MVE method).
#' The function automatically suggests using \code{xout=TRUE} if it is not specified.
#'
#' The Johansen adjustment (controlled by \code{AD} parameter or automatically applied when
#' \code{xout=FALSE}) provides improved Type I error control when samples are small.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' Johansen, S. (1980). The Welch-James approximation to the distribution of the
#' residual sum of squares in a weighted linear regression. \emph{Biometrika}, 67, 85-92.
#'
#' @seealso \code{\link{reg1wayMC}} (parallel version), \code{\link{regci}},
#' \code{\link{reg2ci}}, \code{\link{tsreg}}
#'
#' @examples
#' # Compare regression lines for two groups
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol=1)
#' y1 <- 2 + 3*x1 + rnorm(50)
#' x2 <- matrix(rnorm(50), ncol=1)
#' y2 <- 1 + 2*x2 + rnorm(50)
#' reg1way(list(x1, x2), list(y1, y2), nboot=100)
#'
#' @export
reg1way<-function(x,y,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,outfun=outpro,STAND=TRUE,AD=FALSE,alpha=.05,pr=TRUE,...){
#
#  Test hypothesis that for two or more independent groups, all regression parameters
#  (the intercepts and slopes) are equal
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen MANOVA type test statistic.
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun, which defaults to the MVE method.
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#
if(pr){
if(!xout)print("Might want to consider xout=T to  remove leverage points")
}
if(SEED)set.seed(2)
if(!is.list(x))stop("Argument x should have list mode")
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop("Something is wrong. Number of covariates differs among the groups being compared")
nv=NULL
p=ncol(x[[1]])
p1=p+1
for(j in 1:J){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1:p]
y[[j]]=xy[,p1]
x[[j]]=as.matrix(x[[j]])
nv=c(nv,nrow(x[[j]]))
}
nv.keep=nv
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
}}
x=lapply(x,as.matrix)
K=p1
est=matrix(NA,nrow=J,ncol=p1)
grpnum=NULL
for(j in 1:J)grpnum[j]=paste("Group",j)
vlabs="Intercept"
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)=list(grpnum,vlabs)
ecov=list()
ecovinv=list()
W=rep(0,p1)
gmean=rep(0,p1)
for(j in 1:J){
est[j,]=regfun(x[[j]],y[[j]],xout=FALSE,...)$coef
nv.keep[j]=nrow(x[[j]])
vals=matrix(NA,nrow=nboot,ncol=p1)
data<-matrix(sample(length(y[[j]]),size=length(y[[j]])*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
bvec<-lapply(data,regbootMC,x[[j]],y[[j]],regfun,...)
# bvec is a p+1 by nboot matrix.
vals=t(matl(bvec))
ecov[[j]]=var(vals)
ecovinv[[j]]=solve(ecov[[j]])  #W_j
gmean=gmean+ecovinv[[j]]%*%est[j,]
W=W+ecovinv[[j]]
}
estall=solve(W)%*%gmean
F=0
for(k in 1:K){
for(m in 1:K){
for(j in 1:J){
F=F+ecovinv[[j]][k,m]*(est[j,k]-estall[k])*(est[j,m]-estall[m])
}}}
pvalad=NULL
# if xout=F or AD=T, compute corrected critical value, stemming from Johansen
df=K*(J-1)
if(!xout || AD){
iden=diag(p1)
Aw=0
for(j in 1:J){
temp=iden-solve(W)%*%ecovinv[[j]]
tempsq=temp%*%temp
Aw=Aw+(sum(diag(tempsq))+(sum(diag(temp)))^2)/(nv[j]-1)
}
Aw=Aw/2
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
crit=qchisq(alval[i],df)
critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
if(F<critad)break
}
pvalad=1-irem/1000
}
#
pval=1-pchisq(F,df)
crit=qchisq(1-alpha,df)
critad=NULL
if(!xout || AD)critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
est=data.frame(est)
hlabs=NULL
vlabs="Intercept"
for(j in 1:J)hlabs[j]=paste("Group",j)
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)<-list(hlabs,vlabs)
list(n=nv,n.keep=nv.keep,test.stat=F,crit.value=crit,adjusted.crit=critad,p.value=pval,adjusted.p.value=pvalad,est=est)
}

# ============================================================================
# reg1wayMC
# ============================================================================
#' Compare Regression Lines Across Multiple Independent Groups (Parallel Processing)
#'
#' Tests the hypothesis that two or more independent groups have identical
#' regression parameters using bootstrap estimation and a Johansen MANOVA-type
#' test statistic. This is the multicore version of \code{\link{reg1way}}.
#'
#' @param x A list of length J, where each element is a numeric matrix of predictors
#'   for one group. All groups must have the same number of predictors.
#' @param y A list of length J, where each element is a numeric vector of the
#'   dependent variable for one group.
#' @param regfun Robust regression function to use (default: \code{tsreg}).
#'   Common alternatives: \code{tshdreg}, \code{chreg}, \code{ltsreg}, \code{opreg}.
#' @param STAND Logical. If \code{TRUE}, standardizes predictors when detecting
#'   outliers (default: \code{TRUE}).
#' @param AD Logical. If \code{TRUE}, applies an adjusted critical value to control
#'   Type I error when leverage points are present (default: \code{FALSE}).
#'   Ignored when \code{xout=TRUE}.
#' @param pr Logical. If \code{TRUE}, prints warnings and suggestions (default: \code{TRUE}).
#' @inheritParams common-params
#'
#' @return A list with components:
#' \describe{
#'   \item{test}{The test statistic value}
#'   \item{p.value}{Bootstrap p-value for the global test}
#'   \item{n.all}{Original sample sizes before any outlier removal}
#'   \item{n.keep}{Sample sizes after outlier removal (if \code{xout=TRUE})}
#'   \item{param.est}{Matrix of regression coefficient estimates (rows = groups, columns = parameters)}
#'   \item{Var.Explained}{Matrix of explanatory power for each group}
#'   \item{crit}{Critical value (unadjusted or adjusted if \code{AD=TRUE})}
#' }
#'
#' @details
#' This function extends robust regression comparison to multiple independent groups.
#' It uses bootstrap to estimate the standard errors of regression coefficients and
#' constructs a Johansen-type MANOVA test statistic.
#'
#' The procedure:
#' 1. Fits robust regression for each group
#' 2. Uses bootstrap to estimate standard errors
#' 3. Constructs a test statistic based on the weighted sum of squared differences
#' 4. Compares to bootstrap null distribution
#'
#' When \code{AD=TRUE} and \code{xout=FALSE}, an adjusted critical value is used
#' to protect against Type I error inflation due to leverage points. The adjustment
#' is based on simulation results.
#'
#' The function uses \code{mclapply} for parallel processing, making it faster than
#' \code{\link{reg1way}} for large datasets or many bootstrap samples.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{reg1way}} for non-parallel version,
#'   \code{\link{reg2ci}} for two-group comparisons,
#'   \code{\link{reg1wayISOMC}} for isotonic version
#'
#' @examples
#' \dontrun{
#' # Compare regression lines across three groups
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol=1)
#' y1 <- 1 + 2*x1 + rnorm(50)
#'
#' x2 <- matrix(rnorm(50), ncol=1)
#' y2 <- 1.5 + 2*x2 + rnorm(50)  # Same slope, different intercept
#'
#' x3 <- matrix(rnorm(50), ncol=1)
#' y3 <- 1 + 3*x3 + rnorm(50)  # Different slope
#'
#' result <- reg1wayMC(list(x1, x2, x3), list(y1, y2, y3), nboot=100)
#' print(result$p.value)
#' }
#'
#' @export
reg1wayMC<-function(x,y,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,outfun=outpro,
STAND=TRUE,alpha=.05,pr=TRUE,AD=FALSE,...){
#
#  Test hypothesis that for two or more independent groups, all regression parameters are equal
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen MANOVA type test statistic
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#
if(pr){
if(!xout)print("Might want to consider xout=T to  remove leverage points")
}
if(SEED)set.seed(2)
if(!is.list(x))stop("Argument x should have list mode")
if(pr){
if(xout)print("xout=T, so an adjusted critical is not computed and apparently not needed")
}
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop("Something is wrong. Number of covariates differs among the groups being compared")
nv=NULL
nv.keep=NULL
nv.all=NULL
p=ncol(x[[1]])
p1=p+1
for(j in 1:J){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1:p]
y[[j]]=xy[,p1]
x[[j]]=as.matrix(x[[j]])
nv.all[j]=c(nv,nrow(x[[j]]))
}
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
}}
x=lapply(x,as.matrix)
p1=ncol(x[[1]])+1
K=p1
est=matrix(NA,nrow=J,ncol=p1)
hlabs=NULL
vlabs="Intercept"
for(j in 1:J)hlabs[j]=paste("Group",j)
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)<-list(hlabs,vlabs)
nv=NA
ecov=list()
ecovinv=list()
W=rep(0,p1)
gmean=rep(0,p1)
for(j in 1:J){
est[j,]=regfun(x[[j]],y[[j]])$coef
nv.keep[j]=nrow(x[[j]])
nv[j]=nv.keep[j]
vals=matrix(NA,nrow=nboot,ncol=p1)
data<-matrix(sample(nv[j],size=nv[j]*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
bvec<-mclapply(data,regbootMC,x[[j]],y[[j]],regfun,...)
vals=t(matl(bvec))
ecov[[j]]=var(vals)
ecovinv[[j]]=solve(ecov[[j]])  #W_j
gmean=gmean+ecovinv[[j]]%*%est[j,]
W=W+ecovinv[[j]]
}
estall=solve(W)%*%gmean
F=0
for(k in 1:K){
for(m in 1:K){
for(j in 1:J){
F=F+ecovinv[[j]][k,m]*(est[j,k]-estall[k])*(est[j,m]-estall[m])
}}}
df=K*(J-1)
pvalad=NULL
# if xout=F, compute corrected critical value, stemming from Johansen
df=K*(J-1)
if(!xout || AD){
iden=diag(p1)
Aw=0
for(j in 1:J){
temp=iden-solve(W)%*%ecovinv[[j]]
tempsq=temp%*%temp
Aw=Aw+(sum(diag(tempsq))+(sum(diag(temp)))^2)/(nv[j]-1)
}
Aw=Aw/2
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
crit=qchisq(alval[i],df)
critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
if(F<critad)break
}
pavida=1-irem/1000
}
#
pval=1-pchisq(F,df)
crit=qchisq(1-alpha,df)
critad=NULL
if(AD)critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
est=data.frame(est)
list(n=nv.all,n.keep=nv.keep,test.stat=F,crit.value=crit,adjusted.crit=critad,p.value=pval,adjusted.p.value=pvalad,est=est)
}

# ============================================================================
# reg2ciMC
# ============================================================================
#' Compare Two Independent Regression Lines (Parallel Processing)
#'
#' Computes bootstrap confidence intervals for the difference in regression
#' parameters between two independent groups using parallel processing.
#' This is the multicore version of \code{\link{reg2ci}}.
#'
#' @param x Numeric vector or matrix of predictor variable(s) for the first group.
#' @param y Numeric vector of the dependent variable for the first group.
#' @param x1 Numeric vector or matrix of predictor variable(s) for the second group.
#' @param y1 Numeric vector of the dependent variable for the second group.
#' @param regfun Robust regression function to use (default: \code{tsreg}).
#'   Must return coefficients in \code{regfun$coef} with intercept as first element.
#' @param plotit Logical. If \code{TRUE}, creates a scatterplot with both
#'   regression lines (default: \code{TRUE}).
#' @param pr Logical. If \code{TRUE}, prints progress messages (default: \code{FALSE}).
#' @param xlab Label for x-axis when \code{plotit=TRUE} (default: "X").
#' @param ylab Label for y-axis when \code{plotit=TRUE} (default: "Y").
#' @inheritParams common-params
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Sample sizes for both groups}
#'   \item{ci}{Matrix with columns: Param, ci.low, ci.hi for each regression parameter}
#'   \item{p.value}{P-values for testing difference = 0 for each parameter}
#'   \item{Est.1}{Regression coefficient estimates for group 1}
#'   \item{Est.2}{Regression coefficient estimates for group 2}
#'   \item{Est.dif}{Difference in estimates (group 1 - group 2)}
#' }
#'
#' @details
#' This function compares regression lines between two independent groups using
#' a percentile bootstrap to construct confidence intervals for the differences
#' in parameters. The function uses \code{mclapply} for parallel processing.
#'
#' The bootstrap procedure:
#' 1. Resamples each group independently with replacement
#' 2. Refits the regression models
#' 3. Computes the difference in parameters
#' 4. Constructs percentile confidence intervals from the bootstrap distribution
#'
#' P-values are computed from the proportion of bootstrap samples where the
#' difference has opposite sign from zero.
#'
#' When \code{plotit=TRUE}, creates a scatterplot with both regression lines
#' overlaid. Points from different groups can be distinguished visually.
#'
#' The function uses \code{mclapply} for parallel processing, making it faster than
#' \code{\link{reg2ci}} for large datasets or many bootstrap samples.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{reg2ci}} for non-parallel version,
#'   \code{\link{difreg}} for dependent groups,
#'   \code{\link{reg1wayMC}} for comparing more than two groups
#'
#' @examples
#' \dontrun{
#' # Compare two independent groups
#' set.seed(123)
#' x1 <- rnorm(50)
#' y1 <- 1 + 2*x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- 2 + 3*x2 + rnorm(50)  # Different intercept and slope
#'
#' result <- reg2ciMC(x1, y1, x2, y2, nboot=500)
#' print(result$ci)
#'
#' # Multiple predictors
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- 1 + 2*x1[,1] - x1[,2] + rnorm(50)
#' x2 <- matrix(rnorm(100), ncol=2)
#' y2 <- 1.5 + 2.5*x2[,1] - 1.5*x2[,2] + rnorm(50)
#' reg2ciMC(x1, y1, x2, y2, nboot=500)
#' }
#'
#' @export
reg2ciMC<-function(x,y,x1,y1,regfun=tsreg,nboot=599,alpha=.05,plotit=TRUE,SEED=TRUE,
xout=FALSE,outfun=outpro,pr=FALSE,xlab='X',ylab='Y',...){
#
#   Compute a .95 confidence interval for the difference between the
#   the intercepts and slopes corresponding to two independent groups.
#   The default regression method is Theil-Sen.
#
#   Same as reg2ci, only takes advantage of a multi-core processor
#
#   The predictor values for the first group are
#   assumed to be in the n by p matrix x.
#   The predictors for the second group are in x1
#
#   The default number of bootstrap samples is nboot=599
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
x1<-as.matrix(x1)
xx1<-cbind(x1,y1)
xx1<-elimna(xx1)
x1<-xx1[,1:ncol(x1)]
x1<-as.matrix(x1)
y1<-xx1[,ncol(x1)+1]
x=as.matrix(x)
x1=as.matrix(x1)
if(xout){
if(identical(outfun,outblp)){
flag1=outblp(x,y,plotit=FALSE)$keep
flag2=outblp(x1,y2,plotit=FALSE)$keep
}
if(!identical(outfun,outblp)){
flag1=outfun(x,plotit=FALSE)$keep
flag2=outfun(x1,plotit=FALSE)$keep
}
x=x[flag1,]
y=y[flag1]
x1=x1[flag2,]
y1=y1[flag2]
}
n=length(y)
n[2]=length(y1)
x<-as.matrix(x)
x1<-as.matrix(x1)
est1=regfun(x,y)$coef
est2=regfun(x1,y1)$coef
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec<-mclapply(data,regbootMC,x,y,regfun,mc.preschedule=TRUE,xout=FALSE,...)
bvec=matl(bvec)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec1<-mclapply(data,regbootMC,x1,y1,regfun,mc.preschedule=TRUE,xout=FALSE,...)
bvec1=matl(bvec1)
bvec<-bvec-bvec1
p1<-ncol(x)+1
regci<-matrix(0,p1,6)
dimnames(regci)<-list(NULL,
c("Parameter","ci.lower","ci.upper","p.value","Group 1","Group 2"))
ilow<-round((alpha/2)*nboot)+1
ihi<-nboot-(ilow-1)
for(i in 1:p1){
temp<-sum(bvec[i,]<0)/nboot+sum(bvec[i,]==0)/(2*nboot)
regci[i,4]<-2*min(temp,1-temp)
bsort<-sort(bvec[i,])
regci[i,2]<-bsort[ilow]
regci[i,3]<-bsort[ihi]
regci[,1]<-c(0:ncol(x))
}
regci[,5]=est1
regci[,6]=est2
if(ncol(x)==1 && plotit){
plot(c(x,x1),c(y,y1),type="n",xlab=xlab,ylab=ylab)
points(x,y)
points(x1,y1,pch="+")
abline(regfun(x,y,...)$coef)
abline(regfun(x1,y1,...)$coef,lty=2)
}
list(n=n,output=regci)
}

# ============================================================================
# reg1wayISO
# ============================================================================
reg1wayISO<-function(x,y,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,outfun=outpro,STAND=TRUE,alpha=.05,pr=TRUE,...){
#
#  Test hypothesis that for two or more independent groups, all slope parameters
#  are equal.
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen MANOVA type test statistic.
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun
#
if(SEED)set.seed(2)
if(pr){
if(!xout)print("Might want to consider xout=T to  remove leverage points")
}
if(!is.list(x))stop("Argument x should have list mode")
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop("Something is wrong. Number of covariates differs among the groups being compared")
nv=NULL
p=ncol(x[[1]])
p1=p+1
for(j in 1:J){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1:p]
y[[j]]=xy[,p1]
x[[j]]=as.matrix(x[[j]])
nv=c(nv,nrow(x[[j]]))
}
nv.keep=nv
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
}}
x=lapply(x,as.matrix)
K=p1
est=matrix(NA,nrow=J,ncol=p1)
nv.keep=NULL
ecov=list()
ecovinv=list()
W=rep(0,p1)
gmean=rep(0,p)
for(j in 1:J){
est[j,]=regfun(x[[j]],y[[j]],xout=FALSE,...)$coef
nv.keep[j]=nrow(x[[j]])
vals=matrix(NA,nrow=nboot,ncol=p1)
data<-matrix(sample(length(y[[j]]),size=length(y[[j]])*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
bvec<-lapply(data,regbootMC,x[[j]],y[[j]],regfun,...)
# bvec is a p+1 by nboot matrix.
vals=t(matl(bvec))
ecov[[j]]=var(vals)
ecovinv[[j]]=solve(ecov[[j]])  #W_j
gmean=gmean+ecovinv[[j]][2:K,2:K]%*%est[j,2:K]
W=W+ecovinv[[j]]
}
estall=solve(W[2:K,2:K])%*%gmean
estall=c(0,estall)
F=0
for(k in 2:K){
for(m in 2:K){
for(j in 1:J){
F=F+ecovinv[[j]][k,m]*(est[j,k]-estall[k])*(est[j,m]-estall[m])
}}}
df=p*(J-1)
pvalad=NULL
AD=FALSE # Seems adjusted critical is not needed
if(AD){
iden=diag(p1)
Aw=0
for(j in 1:J){
temp=iden-solve(W)%*%ecovinv[[j]]
tempsq=temp%*%temp
Aw=Aw+(sum(diag(tempsq))+(sum(diag(temp)))^2)/(nv[j]-1)
}
Aw=Aw/2
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
crit=qchisq(alval[i],df)
critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
if(F<critad)break
}
pvalad=1-irem/1000
}
#
pval=1-pchisq(F,df)
crit=qchisq(1-alpha,df)
critad=NULL
if(AD)critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
est=data.frame(est)
hlabs=NULL
vlabs="Intercept"
for(j in 1:J)hlabs[j]=paste("Group",j)
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)<-list(hlabs,vlabs)
list(n=nv,n.keep=nv.keep,test.stat=F,crit.value=crit,p.value=pval,est=est)
}

# ============================================================================
# reg1wayISOMC
# ============================================================================

#' One-Way Regression ANOVA with Isotonic Bootstrap (Parallel Processing)
#'
#' Tests the hypothesis that regression parameters are equal across J independent
#' groups using bootstrap standard errors and a Johansen-type MANOVA statistic.
#' This is the parallel processing (multicore) version using \code{mclapply}.
#'
#' @param x List of length J, where \code{x[[j]]} is an n_j by p matrix of
#'   predictors for group j.
#' @param y List of length J, where \code{y[[j]]} is a vector of length n_j
#'   containing the response values for group j.
#' @param regfun Regression function to use (default: \code{tsreg}).
#' @inheritParams common-params
#' @param STAND Logical; if TRUE, standardizes the test statistic (default: TRUE).
#' @param pr Logical; if TRUE, prints a reminder about using \code{xout=TRUE}
#'   (default: TRUE).
#' @param ... Additional arguments passed to \code{regfun}.
#'
#' @return A list with components:
#' \describe{
#'   \item{test}{The test statistic value.}
#'   \item{crit.val}{Critical value for the test at level \code{alpha}.}
#'   \item{p.value}{Bootstrap p-value for the test.}
#'   \item{n}{Vector of sample sizes for each group (before outlier removal).}
#'   \item{n.keep}{Vector of sample sizes for each group (after outlier removal).}
#' }
#'
#' @details
#' This function extends one-way regression ANOVA to multiple groups using
#' a robust approach:
#'
#' **Null Hypothesis**: All J groups have identical regression coefficients
#' (both intercepts and slopes).
#'
#' **Algorithm**:
#' 1. For each group j = 1, ..., J:
#'    - Optionally remove leverage points using \code{outfun} if \code{xout=TRUE}
#'    - Estimate regression coefficients using \code{regfun}
#'    - Bootstrap to estimate standard errors (using parallel processing)
#' 2. Compute Johansen-type test statistic comparing coefficient estimates
#' 3. Determine significance via bootstrap critical value
#'
#' The isotonic bootstrap is used to ensure better finite-sample performance.
#' Parallel processing via \code{mclapply} makes this feasible even with
#' large \code{nboot}.
#'
#' **Important**: This is a global test. If rejected, use post-hoc methods
#' like \code{\link{reg2ciMC}} for pairwise comparisons.
#'
#' @seealso \code{\link{reg1wayISO}} for non-parallel version,
#'   \code{\link{reg1way}} for standard one-way regression test,
#'   \code{\link{reg2ciMC}} for pairwise comparisons
#'
#' @examples
#' \dontrun{
#' # Three groups with different regression relationships
#' set.seed(123)
#' x <- list()
#' y <- list()
#' for(j in 1:3) {
#'   x[[j]] <- matrix(rnorm(100), ncol=2)
#'   # Different slopes for each group
#'   y[[j]] <- j * x[[j]][,1] + 2*x[[j]][,2] + rnorm(50)
#' }
#'
#' # Test for equal regression coefficients
#' result <- reg1wayISOMC(x, y, nboot=500)
#' print(result$p.value)
#'
#' # If significant, do pairwise comparisons
#' if(result$p.value < 0.05) {
#'   reg2ciMC(x[[1]], y[[1]], x[[2]], y[[2]])
#' }
#' }
#'
#' @export
reg1wayISOMC<-function(x,y,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,outfun=outpro,
STAND=TRUE,alpha=.05,pr=TRUE,...){
#
#  Test hypothesis that for two or more independent groups, all regression parameters are equal
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen MANOVA type test statistic
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#
if(SEED)set.seed(2)
if(pr){
if(!xout)print("Might want to consider xout=T to  remove leverage points")
}
if(!is.list(x))stop("Argument x should have list mode")
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop("Something is wrong. Number of covariates differs among the groups being compared")
nv=NULL
nv.keep=NULL
nv.all=NULL
p=ncol(x[[1]])
p1=p+1
for(j in 1:J){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1:p]
y[[j]]=xy[,p1]
x[[j]]=as.matrix(x[[j]])
nv.all[j]=c(nv,nrow(x[[j]]))
}
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
}}
x=lapply(x,as.matrix)
p1=ncol(x[[1]])+1
K=p1
est=matrix(NA,nrow=J,ncol=p1)
hlabs=NULL
vlabs="Intercept"
for(j in 1:J)hlabs[j]=paste("Group",j)
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)<-list(hlabs,vlabs)
nv=NA
ecov=list()
ecovinv=list()
W=rep(0,p1)
gmean=rep(0,p)
for(j in 1:J){
est[j,]=regfun(x[[j]],y[[j]])$coef
nv.keep[j]=nrow(x[[j]])
nv[j]=nv.keep[j]
vals=matrix(NA,nrow=nboot,ncol=p1)
data<-matrix(sample(nv[j],size=nv[j]*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
bvec<-mclapply(data,regbootMC,x[[j]],y[[j]],regfun,...)
vals=t(matl(bvec))
ecov[[j]]=var(vals)
ecovinv[[j]]=solve(ecov[[j]])  #W_j
gmean=gmean+ecovinv[[j]][2:K,2:K]%*%est[j,2:K]
W=W+ecovinv[[j]]
}
estall=solve(W[2:K,2:K])%*%gmean
estall=c(0,estall)
F=0
for(k in 2:K){
for(m in 2:K){
for(j in 1:J){
F=F+ecovinv[[j]][k,m]*(est[j,k]-estall[k])*(est[j,m]-estall[m])
}}}
df=p*(J-1)
#
df=p*(J-1)
#
pval=1-pchisq(F,df)
crit=qchisq(1-alpha,df)
critad=NULL
#if(!xout || AD)critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
est=data.frame(est)
list(n=nv.all,n.keep=nv.keep,test.stat=F,crit.value=crit,p.value=pval,est=est)
}

# ============================================================================
# tsregNW
# ============================================================================

#' Theil-Sen Regression Using Gauss-Seidel Algorithm
#'
#' Computes the Theil-Sen regression estimator using the Gauss-Seidel algorithm
#' for multiple predictors. This iterative method provides robust regression
#' estimates with high breakdown point and efficiency.
#'
#' @param x A numeric vector or matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param xout Logical; if \code{TRUE}, removes outliers from predictor space
#'   before fitting (default: \code{FALSE}).
#' @param outfun Function for detecting outliers in predictor space (default: \code{out}).
#'   Only used if \code{xout=TRUE}.
#' @param iter Maximum number of iterations for Gauss-Seidel algorithm (default: 10).
#' @param varfun Function for computing variance in explanatory power calculation
#'   (default: \code{pbvar}).
#' @param corfun Function for computing correlation in explanatory power calculation
#'   (default: \code{pbcor}).
#' @param plotit Logical; if \code{TRUE}, creates plots during outlier detection
#'   (default: \code{FALSE}). Only used if \code{xout=TRUE}.
#' @param tol Convergence tolerance for the Gauss-Seidel algorithm (default: 0.0001).
#'   Iteration stops when the maximum absolute change in coefficients is less than \code{tol}.
#' @param ... Additional arguments passed to the outlier detection function.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients: \code{c(intercept, slope_1, ..., slope_p)}.}
#'   \item{residuals}{Always \code{NULL} in current implementation.}
#'   \item{Strength.Assoc}{Strength of association, computed as \code{sqrt(Explanatory.Power)}.}
#'   \item{Explanatory.Power}{Measure of explanatory power, computed as the ratio
#'     \code{varfun(fitted.values) / varfun(y)}, or the squared correlation if
#'     this ratio exceeds 1. Returns \code{NULL} if \code{varfun(y) = 0}.}
#' }
#'
#' @details
#' The Theil-Sen estimator is a robust regression method based on the median of
#' slopes between all pairs of points. For multiple predictors, this function
#' uses the Gauss-Seidel algorithm to iteratively refine the estimates.
#'
#' **Algorithm:**
#' \enumerate{
#'   \item **Initialization**: For each predictor \code{j}, compute the univariate
#'         Theil-Sen slope between \code{x[,j]} and \code{y}
#'   \item **Gauss-Seidel Iteration**: For each predictor \code{j} in turn:
#'         \itemize{
#'           \item Compute partial residuals: \code{r[,j] = y - X %*% beta - alpha + beta[j] * x[,j]}
#'           \item Update \code{beta[j]} using univariate Theil-Sen on \code{(x[,j], r[,j])}
#'         }
#'   \item Update the intercept: \code{alpha = median(y - X %*% beta)}
#'   \item Repeat step 2-3 until convergence or \code{iter} iterations reached
#' }
#'
#' **Convergence**: The algorithm stops when the maximum absolute change in any
#' coefficient is less than \code{tol}, or after \code{iter} iterations.
#'
#' **For univariate regression (p=1):**
#' The function calls \code{\link{tsp1reg}} directly, which computes the exact
#' Theil-Sen estimator as the median of all pairwise slopes.
#'
#' **Explanatory Power:**
#' The function computes a robust measure of fit using the specified variance
#' function (default: percentage bend variance). If the ratio of fitted value
#' variance to response variance exceeds 1 (which can occur with robust variance
#' estimators), it falls back to the squared correlation.
#'
#' **Advantages:**
#' - High breakdown point (up to 29.3% for large samples)
#' - Robust to outliers in both x and y
#' - Asymptotically normal with good efficiency
#' - Resistant to leverage points when combined with \code{xout=TRUE}
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' Theil, H. (1950). A rank-invariant method of linear and polynomial regression
#' analysis, I, II, III. \emph{Proceedings of the Koninklijke Nederlandse Akademie
#' van Wetenschappen A}, 53, 386-392, 521-525, 1397-1412.
#'
#' Sen, P.K. (1968). Estimates of the regression coefficient based on Kendall's tau.
#' \emph{Journal of the American Statistical Association}, 63, 1379-1389.
#'
#' @seealso \code{\link{tsreg}} for standard Theil-Sen implementation,
#'   \code{\link{tsp1reg}} for univariate Theil-Sen,
#'   \code{\link{tshdreg}} for Theil-Sen with Harrell-Davis quantile estimator,
#'   \code{\link{regci}} for confidence intervals
#'
#' @examples
#' \dontrun{
#' # Simple linear regression
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50, sd=0.5)
#' result <- tsregNW(x, y)
#' print(result$coef)
#' print(result$Explanatory.Power)
#'
#' # Multiple regression
#' x <- matrix(rnorm(150), ncol=3)
#' y <- 1 + 2*x[,1] - x[,2] + 0.5*x[,3] + rnorm(50, sd=0.3)
#' tsregNW(x, y)
#'
#' # With outliers - use outlier removal
#' y[1:5] <- y[1:5] + 5
#' tsregNW(x, y, xout=TRUE)
#'
#' # Adjust convergence criteria
#' tsregNW(x, y, iter=20, tol=0.00001)
#'
#' # Using different robust variance estimator
#' tsregNW(x, y, varfun=winvar)
#' }
#'
#' @export
tsregNW<-function(x,y,xout=FALSE,outfun=out,iter=10,varfun=pbvar,
corfun=pbcor,plotit=FALSE,tol=.0001,...){
#
#  Compute Theil-Sen regression estimator
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(ncol(x)==1){
temp1<-tsp1reg(x,y)
coef<-temp1$coef
res<-temp1$res
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-tsp1reg(x[,p],y)$coef[2]
}
res<-y-x%*%temp
alpha<-median(res)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-tsp1reg(x[,p],r[,p],plotit=FALSE)$coef[2]
}
if(max(abs(temp-tempold))<tol)break
alpha<-median(y-x%*%temp)
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
yhat<-y-res
stre=NULL
temp=varfun(y)
#if(temp==0)print('Warning: When computing strength of association, measure of variation=0')
e.pow=NULL
if(temp>0){
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}}
res=NULL
list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

# ============================================================================
# reg2cimcp
# ============================================================================
reg2cimcp<-function(x,y,regfun=tsreg,nboot=599,alpha=0.05,
SEED=TRUE,xout=FALSE,outfun=out,...){
#
# Like reg2ci only x1 etc have list mode containing
# data for J>1 groups. For all pairs of groups are compared via a
# call to reg2ci.
#
#  x list mode contain a matrix of predictors.
#  x[[1]] contains predictors for first group
#  y[[1]] dependent variable for first group.
#
#
if(!is.list(x))stop('x and y should have list mode')
J=length(x) # number of groups
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
res=reg2ci(x[[j]],y[[j]],x[[k]],y[[k]],regfun=regfun,nboot=nboot,alpha=alpha,
plotit=FALSE,xout=xout,outfun=outfun,WARN=FALSE,...)
print(paste('Group', j,'Group', k))
print(res)
}}}
}

# ============================================================================
# reg1mcp
# ============================================================================
reg1mcp<-function(x,y,regfun=tsreg,SEED=TRUE,nboot=100,xout=FALSE,outfun=outpro,STAND=TRUE,alpha=.05,
pr=TRUE,MC=FALSE,...){
#
#  Perform all pairwise comparisons of intercepts among J independent groups
#  Do the same of the first slope, followed by the 2nd slope, etc.
#
#  Control FWE via Hochberg's methods for each set of
#  (J^2-J)/2 parameters. That is, control FWE for the intercepts
#  Do the same for the first slope, etc.
#
#  #  x and y are assumed to have list mode having
#  length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun,
#   which defaults to the projection method.
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#   output contains all pairwise comparisons
#   For each parameter, FWE is controlled using Hochberg's method
#   So by default, for the intercepts,
#   all pairwise comparisons are performed with FWE=.05
#   For the first slope, all pairwise comparisons are performed
#   with FWE=.05, etc.
#
if(SEED)set.seed(2)
if(!is.list(x))stop('Argument x should have list mode')
if(!is.list(y))stop('Argument y should have list mode')
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop('Something is wrong. Number of covariates differs among the groups being compared')
nv=NULL
p=ncol(x[[1]])
p1=p+1
for(j in 1:J){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1:p]
y[[j]]=xy[,p1]
x[[j]]=as.matrix(x[[j]])
nv=c(nv,nrow(x[[j]]))
}
nv.keep=nv
critrad=NULL
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
nv.keep[j]=length(y[[j]])
}}
tot=(J^2-J)/2
dvec<-alpha/c(1:tot)
outl=list()
nr=tot*p1
outp=matrix(NA,ncol=7,nrow=nr)
x=lapply(x,as.matrix)
rlab=rep('Intercept',tot)
xx=list()
yy=list()
iall=0
ivp=c(1,tot)-tot
for(ip in 1:p){
#iv=ip-1
rlab=c(rlab,rep(paste('slope',ip),tot))
}
i=0
sk=1+tot*p
st=seq(1,sk,tot)
st=st-1
for(j in 1:J){
for(k in 1:J){
if(j < k){
i=i+1
st=st+1
xx[[1]]=x[[j]][,1:p]
xx[[2]]=x[[k]][,1:p]
yy[[1]]=y[[j]]
yy[[2]]=y[[k]]
if(!MC)temp=reg2ci(xx[[1]],yy[[1]],xx[[2]],yy[[2]],regfun=regfun)$output
if(MC)temp=reg2ci(xx[[1]],yy[[1]],xx[[2]],yy[[2]],regfun=regfun)$output
iall=iall+1
outp[iall,1]=j
outp[iall,2]=k
outp[st,3]=temp[,4]
outp[st,5]=temp[,2]
outp[st,6]=temp[,3]
}}}
for(i in 1:p1){
ivp=ivp+tot
temp2<-order(0-outp[ivp[1]:ivp[2],3])
icc=c(ivp[1]:ivp[2])
icc[temp2]=dvec
outp[ivp[1]:ivp[2],4]=icc
}
flag=(outp[,3]<=outp[,4])
outp[,7]=rep(0,nr)
outp[flag,7]=1
v=outp[1:tot,1]
vall=rep(v,p1)
outp[,1]=vall
v=outp[1:tot,2]
vall=rep(v,p1)
outp[,2]=vall
#outp[,7]=p.adjust(outp[,3],method=method)
dimnames(outp)=list(rlab,c('Group','Group','p.value','p.crit','ci.low','ci.hi','Sig'))
list(n=nv,n.keep=nv.keep,output=outp)
}

# ============================================================================
# chregF
# ============================================================================

#' Fast Coakley-Hettmansperger Robust Regression
#'
#' Computes Coakley and Hettmansperger's robust regression estimators with
#' a streamlined implementation. This is a faster variant of \code{\link{chreg}}
#' that performs a single iteration rather than iterating to convergence.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @param bend Bending constant for Huber's Psi function (default: 1.345).
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility
#'   (default: \code{FALSE}). The seed is set to 12 when \code{TRUE}.
#' @inheritParams common-params
#' @param outfun Outlier detection function (default: \code{out}). See
#'   \code{\link{chreg}} for other options.
#'
#' @return A list with components:
#' \describe{
#'   \item{coef}{Vector of regression coefficients (intercept followed by slopes)}
#'   \item{residuals}{Vector of residuals from the fit}
#' }
#'
#' @details
#' This function implements the Coakley-Hettmansperger (1993) robust regression
#' method with a single-iteration approximation. The method provides protection
#' against both leverage points and outliers in Y.
#'
#' The algorithm:
#' 1. Computes an initial LTS (Least Trimmed Squares) estimate
#' 2. Computes Mallows weights based on Mahalanobis distance of predictors
#' 3. Computes Huber Psi weights based on standardized residuals
#' 4. Updates the coefficient estimates using weighted least squares
#'
#' Unlike \code{\link{chreg}} which iterates to convergence, \code{chregF}
#' performs only a single update step, making it faster but potentially less
#' accurate. Use \code{chregF} when speed is important and you're willing to
#' accept a less refined estimate.
#'
#' The Mallows weighting downweights high-leverage points, while the Huber
#' Psi function provides robustness to outliers in the dependent variable.
#'
#' @references
#' Coakley, C.W. and Hettmansperger, T.P. (1993). A bounded influence, high
#' breakdown, efficient regression estimator. \emph{Journal of the American
#' Statistical Association}, 88, 872-880.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{chreg}} for the fully iterated version,
#'   \code{\link{bmreg}} for bounded M-regression,
#'   \code{\link{ltsreg}} for LTS regression
#'
#' @examples
#' # Simple regression - fast approximation
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2 + 3*x + rnorm(100)
#' y[1:5] <- y[1:5] + 10  # Add outliers
#' result <- chregF(x, y, SEED=TRUE)
#' result$coef
#'
#' # Multiple regression with leverage points
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - x[,2] + rnorm(50)
#' x[1:3,] <- x[1:3,] + 5  # Add leverage points
#' chregF(x, y, xout=TRUE)
#'
#' @export
chregF<-function(x,y,bend=1.345,SEED=FALSE,xout=FALSE,outfun=out,...){
#
# Compute Coakley Hettmansperger robust regression estimators
# JASA, 1993, 88, 872-880
#
# x is a n by p matrix containing the predictor values.
#
# No missing values are allowed
#
#  Comments in this function follow the notation used
#  by Coakley and Hettmansperger
#
# with old version of R, need library(lqs) when using ltsreg
# as the initial estimate.
#
if(SEED)set.seed(12) # Set seed so that results are always duplicated.
x<-as.matrix(x)
p<-ncol(x)
m<-elimna(cbind(x,y))
x<-m[,1:p]
p1<-p+1
y<-m[,p1]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
x<-as.matrix(x)
cutoff<-bend
mve<-vector('list')
if(ncol(x)==1){
mve$center<-median(x)
mve$cov<-mad(x)^2
}
if(ncol(x)>=2)mve<-cov.mve(x)  # compute minimum volume ellipsoid measures of
                 # location and scale and store in mve.
reg0<-ltsReg(x,y) # compute initial regression est using least trimmed
                 # squares.
# Next, compute the rob-md2(i) values and store in rob
rob<-1  # Initialize vector rob
mx<-mve$center
rob<-mahalanobis(x,mx,mve$cov)
k21<-qchisq(.95,p)
c62<-k21/rob
vecone<-c(rep(1,length(y))) # Initialize vector vecone to 1
c30<-pmin(vecone,c62)  # mallows weights put in c30
k81<-median(abs(reg0$residuals)) # median of absolute residuals
k72<-1.4826*(1+(5/(length(y)-p-1)))*k81 # lms scale
c60<-reg0$residuals/(k72*c30) # standardized residuals
#  compute psi and store in c27
cvec<-c(rep(cutoff,length(y))) # Initialize vector cvec to cutoff
c27<-pmin(cvec,c60)
c27<-pmax(-1*cutoff,c27)  #c27 contains psi values
#
# compute B matrix and put in c66.
#  Also, transform B so that i th diag elem = 0 if c27[i] is
#  between -cutoff and cutoff, 1 otherwise.
#
c66<-ifelse(abs(c27)<=bend,1,0) # Have derivative of psi in c66
m1<-cbind(1,x)  # X matrix with col of 1's added
m2<-t(m1)   #X transpose
m5<-diag(c30) # matrix W, diagonal contains weights
m4<-diag(c66) # B matrix
m6<-m4%*%m1   # BX
m7<-m2%*%m6   # X'BX (nD=X'BX)
m8<-solve(m7)  #m8 = (X'-B-X)inverse
m9<-m8%*%m2 #m9=X prime-B-X inverse X'
m9<-m9%*%m5 # m9=X prime-B-X inverse X'W
m10<-m9%*%c27
c20<-m10*k72
c21<-reg0$coef+c20 #update initial estimate of parameters.
res<-y-m1%*%c21
list(coef=t(c21),residuals=res)
}

# ============================================================================
# DregGOLS
# ============================================================================
DregGOLS<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,SEED=TRUE,nboot=200,
STAND=TRUE,...){
#
#  Global test that two dependent (time 1 and time 2)
#  OLS regression lines are identical
#
#  Use a variation of Hotelling's test coupled with a bootstrap
#  estimate of the relevant covariance matrix associated with the differences
#  in the estimates of the parameters.
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
n=length(y1)
if(xout){
opf=identical(outfun,outpro)
if(!opf){
flag1=outfun(x1)$out.id
flag2=outfun(x2)$out.id
}
if(opf){
flag1=outpro(x1,STAND=STAND)$out.id
flag2=outfun(x2,STAND=STAND)$out.id
}
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
nk=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data,1,regboot,x1,y1,regfun=lsfit,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-apply(data,1,regboot,x2,y2,regfun=lsfit,...)
dif=t(bvec1-bvec2)
S=cov(dif)
est1=lsfit(x1,y1)$coef
est2=lsfit(x2,y2)$coef
est=est1-est2
k <- (nk-p1)/((nk - 1)*p1)
        stat <- k * crossprod(est, solve(S, est))[1, ]
        pvalue <- 1 - pf(stat, p1, nk - p1)
list(test.statistic = stat, degrees_of_freedom = c(p1, nk - p1), p.value =
pvalue,est.1=est1,est.2=est2,estimate.dif = est)
}

# ============================================================================
# difregOLS
# ============================================================================

#' Compare OLS Regression Parameters for Dependent Groups
#'
#' Computes bootstrap confidence intervals and tests for comparing ordinary
#' least squares (OLS) regression parameters between two dependent (paired) groups.
#' This is the OLS-specific version of \code{\link{difreg}}.
#'
#' @param x1 Numeric vector or matrix of predictor variable(s) for the first group.
#' @param y1 Numeric vector of the dependent variable for the first group.
#' @param x2 Numeric vector or matrix of predictor variable(s) for the second group.
#' @param y2 Numeric vector of the dependent variable for the second group.
#' @param regfun Regression function (default: \code{lsfit} for OLS).
#'   Should not be changed unless using a specific OLS variant.
#' @inheritParams common-params
#' @param plotit Logical. If \code{TRUE}, creates a scatterplot with both
#'   regression lines (default: \code{FALSE}).
#' @param xlab Label for x-axis when \code{plotit=TRUE} (default: "X").
#' @param ylab Label for y-axis when \code{plotit=TRUE} (default: "Y").
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Original sample size before any outlier removal}
#'   \item{n.keep}{Sample size after outlier removal (if \code{xout=TRUE})}
#'   \item{est.dif}{Difference in estimates (group 1 - group 2)}
#'   \item{est.1}{OLS coefficient estimates for group 1}
#'   \item{est.2}{OLS coefficient estimates for group 2}
#'   \item{test.stat}{t-statistics for each parameter difference}
#'   \item{standard.error}{Bootstrap standard errors}
#'   \item{p.values}{Two-sided p-values for each parameter difference}
#'   \item{conf.intervals}{Matrix with columns: Param, ci.low, ci.hi}
#' }
#'
#' @details
#' This function compares OLS regression parameters between two dependent (paired)
#' groups using bootstrap methods. It's specifically designed for OLS regression
#' with paired/repeated measures data.
#'
#' **Method**:
#' 1. Bootstrap resamples preserve pairing by resampling entire rows
#' 2. Fits OLS to each group in each bootstrap sample
#' 3. Computes differences in parameters
#' 4. Uses bootstrap SE and t-distribution for inference
#'
#' **Inference**: Unlike \code{\link{difreg}} which uses percentile bootstrap,
#' this function uses bootstrap standard errors with t-distribution critical values.
#' This provides both confidence intervals and p-values.
#'
#' **When to use**:
#' - OLS regression (not robust regression)
#' - Paired/repeated measures design (pre-post, matched pairs, etc.)
#' - Want parametric-style inference with bootstrap SEs
#' - Comparing specific parameters (intercepts, slopes)
#'
#' **Comparison with alternatives**:
#' - \code{\link{difreg}}: For robust regression (tsreg, chreg, etc.)
#' - \code{\link{DregGOLS}}: For global test of all parameters
#' - \code{\link{reg2ci}}: For independent groups
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press. Chapter 10.
#'
#' @seealso \code{\link{difreg}} for robust regression version,
#'   \code{\link{DregGOLS}} for global test,
#'   \code{\link{reg2ci}} for independent groups
#'
#' @examples
#' \dontrun{
#' # Pre-post intervention data
#' set.seed(123)
#' n <- 50
#' x1 <- rnorm(n)
#' y1 <- 1 + 2*x1 + rnorm(n)
#' x2 <- x1 + rnorm(n, sd=0.3)
#' y2 <- 1.5 + 2.5*x2 + rnorm(n)
#'
#' # Compare OLS regression parameters
#' result <- difregOLS(x1, y1, x2, y2, nboot=500)
#' print(result$p.values)
#' print(result$conf.intervals)
#'
#' # Multiple predictors
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- 1 + 2*x1[,1] - x1[,2] + rnorm(50)
#' x2 <- x1 + matrix(rnorm(100, sd=0.2), ncol=2)
#' y2 <- 1.5 + 2.2*x2[,1] - 1.1*x2[,2] + rnorm(50)
#' difregOLS(x1, y1, x2, y2, nboot=500, plotit=FALSE)
#' }
#'
#' @export
difregOLS<-function(x1,y1,x2,y2,regfun=lsfit,xout=FALSE,outfun=outpro,nboot=200,
alpha=.05,SEED=TRUE,plotit=FALSE,xlab='X',ylab='Y',...){
#
# OLS regression data from two different times i.e., two dependent groups
#
#  compute confidence interval for the difference between intercepts
#  and the slopes
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
n=length(y1)
if(xout){
flag1=outfun(x1)$out.id
flag2=outfun(x2)$out.id
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
nk=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data,1,regboot,x1,y1,regfun,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-apply(data,1,regboot,x2,y2,regfun,...)
dif=t(bvec1)-t(bvec2)
est1=lsfit(x1,y1)$coef
est2=lsfit(x2,y2)$coef
estdif=est1-est2
se=apply(dif,2,sd)
pvec=NA
test=NA
test=estdif/se
df=nk-1
pvec=2*(1-pt(abs(test),df))
if(plotit){
reg2plot(x1,y1,x2,y2,xlab=xlab,ylab=ylab)
}
lvec='Intercept'
ci=matrix(NA,nrow=p1,ncol=3)
ci[,1]=c(0:p)
ci[,2]=estdif+qt(alpha/2,df)*se
ci[,3]=estdif-qt(alpha/2,df)*se
dimnames(ci)=list(NULL,c('Param','ci.low','ci.hi'))
for(j in 2:p1)lvec=c(lvec,paste('slope',j-1))
pvec=array(pvec,dimnames=lvec)
list(n=n,n.keep=nk,est.dif=estdif,est.1=est1,est.2=est2,
test.stat=test,standard.error=se,p.values=pvec,conf.intervals=ci)
}

# ============================================================================
# difregYvar
# ============================================================================
difregYvar<-function(x1,y1,x2,y2,regfun=tsreg,pts=NULL,
nboot=100,xout=FALSE,outfun=out,SEED=TRUE,...){
#
#  Estimate standard error of difference between the predicted value of Y
#  corresponding to two dependent groups using regression estimator indicated by
#  the argument
#  regfun
#  corresponding to the points in
#  pts
#  regfun defaults to tsreg, the Theil--Sen estimator
#  pts default is to use all unique points among x1 and x2
#
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
x1<-as.matrix(x1)
x2=as.matrix(x2)
if(is.null(pts)){
pts=rbind(x1,x2)
pts=unique(pts)
}
pts=as.matrix(pts)
nvpts=nrow(pts)
bvec1=matrix(NA,nrow=nboot,ncol=nvpts)
bvec2=matrix(NA,nrow=nboot,ncol=nvpts)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
for(ib in 1:nboot){
bvec1[ib,]=regYsub(x1[data[ib,],],y1[data[ib,]],pts,p1=p1,regfun=regfun,...)
bvec2[ib,]=regYsub(x2[data[ib,],],y2[data[ib,]],pts,p1=p1,regfun=regfun,...)
}
bvec=bvec1-bvec2
sqsd=apply(bvec,2,var)
sqsd
}

# ============================================================================
# difreg
# ============================================================================

#' Compare Two Dependent Regression Lines
#'
#' Computes bootstrap confidence intervals for the difference in regression
#' parameters between two dependent groups (e.g., repeated measures at two time
#' points on the same subjects).
#'
#' @param x1 Numeric vector or matrix of predictor variable(s) for the first group.
#' @param y1 Numeric vector of the dependent variable for the first group.
#' @param x2 Numeric vector or matrix of predictor variable(s) for the second group.
#' @param y2 Numeric vector of the dependent variable for the second group.
#' @inheritParams common-params
#' @param regfun Robust regression function to use (default: \code{tsreg}).
#'   Common alternatives: \code{tshdreg}, \code{chreg}, \code{ltsreg}, \code{opreg}.
#' @param plotit Logical. If \code{TRUE}, creates a scatterplot with both
#'   regression lines (default: \code{FALSE}).
#' @param xlab Label for x-axis when \code{plotit=TRUE} (default: "X").
#' @param ylab Label for y-axis when \code{plotit=TRUE} (default: "Y").
#' @param pr Logical. If \code{TRUE}, prints warnings about duplicate Y values
#'   when using \code{tsreg} (default: \code{TRUE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Original sample size before any outlier removal}
#'   \item{n.keep}{Sample size after outlier removal (if \code{xout=TRUE})}
#'   \item{param}{Character vector labeling parameters (Intercept, slope1, ...)}
#'   \item{p.values}{P-values for testing difference = 0 for each parameter}
#'   \item{est.grp1}{Regression coefficient estimates for group 1}
#'   \item{est.grp2}{Regression coefficient estimates for group 2}
#'   \item{conf.intervals}{Matrix with columns: Param, ci.low, ci.hi}
#' }
#'
#' @details
#' This function analyzes regression data from dependent (paired) groups, such as
#' measurements taken at two different time points on the same subjects. The
#' function uses a percentile bootstrap to construct confidence intervals for
#' the difference in regression parameters.
#'
#' The bootstrap procedure resamples entire rows (preserving the dependency
#' between observations) and refits the regression models. P-values are computed
#' using the proportion of bootstrap samples where the difference has the opposite
#' sign from zero.
#'
#' When \code{xout=TRUE}, outliers detected by \code{outfun} in either group
#' are removed from both groups (to preserve pairing).
#'
#' If duplicate Y values are detected and \code{regfun=tsreg}, a warning suggests
#' using \code{tshdreg} which may have better power with tied values.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{DregG}} for an alternative approach,
#'   \code{\link{difregMC}} for parallel processing version,
#'   \code{\link{reg2ci}} for independent groups
#'
#' @examples
#' # Simulated pre-post intervention data
#' set.seed(123)
#' n <- 50
#' x1 <- rnorm(n)  # Pre-intervention predictor
#' y1 <- 1 + 2*x1 + rnorm(n)  # Pre-intervention outcome
#' x2 <- x1 + rnorm(n, sd=0.3)  # Post-intervention predictor
#' y2 <- 1.5 + 2.5*x2 + rnorm(n)  # Post-intervention (changed slope/intercept)
#'
#' # Test for difference in regression parameters
#' result <- difreg(x1, y1, x2, y2, nboot=500)
#' print(result$conf.intervals)
#'
#' # With plot
#' difreg(x1, y1, x2, y2, plotit=TRUE, nboot=500)
#'
#' @export
difreg<-function(x1,y1,x2,y2,regfun=tsreg,xout=FALSE,outfun=outpro,nboot=599,
alpha=.05,SEED=TRUE,plotit=FALSE,xlab='X',ylab='Y',pr=TRUE,...){
#
# regression data from two different times i.e., two dependent groups
#
#  compute confidence interval for the difference in the slopes
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
n=length(y1)
if(xout){
flag1=outfun(x1,...)$out.id
flag2=outfun(x2,...)$out.id
flag=unique(c(flag1,flag2))
X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
flagF=identical(regfun,tsreg)
if(flagF){
if(pr){
if(sum(duplicated(y1)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
pr=FALSE
}
if(pr){
if(sum(duplicated(y2)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
}}
nk=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec1<-lapply(data,regboot,x1,y1,regfun,xout=FALSE,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-lapply(data,regboot,x2,y2,regfun,xout=FALSE,...)
bvec1=matl(bvec1)
bvec2=matl(bvec2)
dif=t(bvec1)-t(bvec2)
dif.sort=apply(dif,2,sort)
pvec=NA
for(i in 1:p1){
pvec[i]<-(sum(dif[,i]<0)+.5*sum(dif[,i]==0))/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
}
pvec<-2*pvec
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=matrix(NA,nrow=p1,ncol=3)
ci[,1]=c(0:p)
for(i in 1:p1){
ci[i,2]=dif.sort[ilow,i]
ci[i,3]=dif.sort[ihi,i]
}
dimnames(ci)=list(NULL,c('Param','ci.low','ci.hi'))
if(plotit){
reg2plot(x1,y1,x2,y2,xlab=xlab,ylab=ylab,regfun=regfun,...)
}
lvec='Intercept'
for(j in 2:p1)lvec=c(lvec,paste('slope',j-1))
#pvec=array(pvec,dimnames=lvec)
est1=regfun(x1,y1,xout=FALSE,...)$coef
est2=regfun(x2,y2,xout=FALSE,...)$coef
list(n=n,n.keep=nk,param=lvec,p.values=pvec,est.grp1=est1,est.grp2=est2,conf.intervals=ci)
}

# ============================================================================
# tshdreg
# ============================================================================
#' Theil-Sen Regression with Harrell-Davis Estimator
#'
#' Computes the Theil-Sen robust regression estimator using the Harrell-Davis
#' estimator for the intercept. This variant can provide improved efficiency
#' compared to the standard median-based approach while maintaining robustness.
#'
#' @param x A numeric vector or matrix containing the predictor variable(s). For
#'   multiple predictors, rows represent observations and columns represent variables.
#' @inheritParams common-params
#' @param HD Logical. If `TRUE`, uses Harrell-Davis estimator for the intercept;
#'   if `FALSE`, uses median (default: `TRUE`).
#' @param outfun Outlier detection function (default: `out`). See \code{\link{tsreg}}
#'   for other options.
#' @param iter Number of iterations for the back-fitting algorithm when there are
#'   multiple predictors (default: 5).
#' @param varfun Function used to compute variance for strength of association
#'   (default: `pbvar`).
#' @param tr Logical or numeric. Trimming parameter for Harrell-Davis estimator
#'   (default: `FALSE` for no trimming).
#' @param do.stre Logical. If `TRUE`, computes strength of association measure
#'   (default: `TRUE`).
#' @param corfun Correlation function used for computing explanatory power when
#'   the variance ratio exceeds 1 (default: `pbcor`).
#' @param tol Convergence tolerance for the back-fitting algorithm (default: 0.0001).
#'   Iteration stops when the maximum absolute change in coefficients is less than `tol`.
#' @param RES Logical. If `TRUE`, returns residuals in the output (default: `TRUE`).
#'   Set to `FALSE` to save memory when residuals are not needed.
#' @param OPT Logical. If `TRUE`, computes intercept as `HD(y) - b1*HD(X)`;
#'   if `FALSE`, computes intercept as `HD(y - b1*X)` (default: `FALSE`).
#' @param xlab Label for x-axis when `plotit=TRUE` (default: "X").
#' @param ylab Label for y-axis when `plotit=TRUE` (default: "Y").
#'
#' @return A list with components:
#'   \item{coef}{Vector of regression coefficients (intercept followed by slopes)}
#'   \item{residuals}{Vector of residuals (NULL if `RES=FALSE`)}
#'   \item{Strength.Assoc}{Strength of association (square root of explanatory power)}
#'   \item{Explanatory.Power}{Proportion of variance explained}
#'
#' @details
#' This function implements the Theil-Sen estimator with the Harrell-Davis estimator
#' used for computing the intercept and for the back-fitting algorithm (when multiple
#' predictors are present). The Harrell-Davis estimator is a weighted average of
#' order statistics that can be more efficient than the sample median while
#' maintaining robustness.
#'
#' For multiple predictors, a back-fitting algorithm is used, which iteratively
#' updates each coefficient while holding the others fixed. The algorithm continues
#' until convergence (based on `tol`) or until `iter` iterations are completed.
#'
#' When using this function with bootstrap methods, the residuals must be available.
#' Ensure `RES=TRUE` when using `tshdreg` as the `regfun` argument in bootstrap
#' functions.
#'
#' @references
#' Harrell, F. E. & Davis, C. E. (1982). A new distribution-free quantile estimator.
#' \emph{Biometrika}, 69, 635-640.
#'
#' Wilcox, R. R. (2022). \emph{Introduction to Robust Estimation and Hypothesis
#' Testing} (5th ed.). Academic Press.
#'
#' @seealso \code{\link{tsreg}} for standard Theil-Sen regression,
#'   \code{\link{tsp1reg}} for single predictor version, \code{\link{hd}}
#'   for Harrell-Davis estimator
#'
#' @examples
#' # Simple regression with Harrell-Davis estimator
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50)
#' result <- tshdreg(x, y)
#' result$coef
#'
#' # Multiple regression with outlier removal
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - 3*x[,2] + rnorm(50)
#' tshdreg(x, y, xout=TRUE)
#'
#' # Ensure residuals are available for bootstrap use
#' result <- tshdreg(x, y, RES=TRUE)
#'
#' @export
tshdreg<-function(x,y,HD=TRUE,xout=FALSE,outfun=out,iter=5,varfun=pbvar,tr=FALSE,do.stre=TRUE,
corfun=pbcor,plotit=FALSE,tol=.0001,RES=TRUE,OPT=FALSE,xlab='X',ylab='Y',...){
#
#  Compute Theil-Sen regression estimator
#
#  Use back-fitting
#  when there is more than one predictor
#  and estimate intercept using Harrel-Davis estimator
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(ncol(x)==1){
temp1<-tshd(x,y,HD=HD,plotit=plotit,xlab=xlab,ylab=ylab,OPT=OPT,tr=tr)
coef<-temp1$coef
res<-y-coef[2]*x-coef[1]
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-tshd(x[,p],y)$coef[2]
}
res<-y-x%*%temp
alpha<-hd(res)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-tshd(x[,p],r[,p],plotit=FALSE,tr=tr)$coef[2]
}
if(max(abs(temp-tempold))<tol)break
alpha<-hd(y-x%*%temp,tr=tr)
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
yhat<-y-res
stre=NULL
e.pow=NULL
if(do.stre){
temp=varfun(y)
if(temp==0)print('Warning: When computing strength of association, measure of variation=0')
if(temp>0){
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}}}
if(!RES)res=NULL
list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow,residuals=res)
}

# ============================================================================
# ltsreg
# ============================================================================
#' Least Trimmed Squares Regression
#'
#' Computes the Least Trimmed Squares (LTS) regression estimator, a highly
#' robust regression method that minimizes the sum of the h smallest squared
#' residuals. Uses the \code{ltsReg} function from the \pkg{robustbase} package.
#'
#' @param x A numeric vector or matrix containing the predictor variable(s). For
#'   multiple predictors, rows represent observations and columns represent variables.
#' @inheritParams common-params
#' @param tr Proportion of observations to trim (default: 0.5). The function
#'   uses the (1-tr) proportion with smallest squared residuals. Higher values
#'   provide more robustness but lower efficiency.
#' @param outfun Outlier detection function for removing outliers from predictors
#'   when `xout=TRUE` (default: `outpro`).
#' @param STAND Logical. Passed to outlier detection function (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{coef}{Vector of regression coefficients (intercept followed by slopes)}
#'   \item{residuals}{Vector of raw residuals}
#'
#' @details
#' The Least Trimmed Squares estimator minimizes the sum of the h smallest
#' squared residuals, where h = n*(1-tr). This makes LTS highly robust to
#' outliers with a breakdown point of approximately 50% when tr=0.5.
#'
#' The default trimming of tr=0.5 provides maximum robustness, using only
#' the best-fitting 50% of the data. This is much more robust than ordinary
#' least squares (breakdown point of 0%) or M-estimators (breakdown point
#' typically around 10-20%).
#'
#' This function requires the \pkg{robustbase} package to be installed.
#'
#' @references
#' Rousseeuw, P. J. (1984). Least median of squares regression.
#' \emph{Journal of the American Statistical Association}, 79, 871-880.
#'
#' Rousseeuw, P. J. & Van Driessen, K. (2006). Computing LTS regression for
#' large data sets. \emph{Data Mining and Knowledge Discovery}, 12, 29-45.
#'
#' Wilcox, R. R. (2022). \emph{Introduction to Robust Estimation and Hypothesis
#' Testing} (5th ed.). Academic Press.
#'
#' @seealso \code{\link{tsreg}} for Theil-Sen regression, \code{\link{LMSreg}}
#'   for Least Median of Squares, \code{\link{MMreg}} for MM-estimator,
#'   \code{\link[robustbase]{ltsReg}} in package \pkg{robustbase}
#'
#' @examples
#' # Simple regression
#' x <- rnorm(50)
#' y <- 2 + 3*x + rnorm(50)
#' result <- ltsreg(x, y)
#' result$coef
#'
#' # With outliers - LTS is highly robust
#' y[1:5] <- y[1:5] + 10  # Add outliers
#' ltsreg(x, y)
#'
#' # Multiple regression with less trimming
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 1 + 2*x[,1] - 3*x[,2] + rnorm(50)
#' ltsreg(x, y, tr=0.25)  # Use 75% of data
#'
#' @export
ltsreg<-function(x,y,tr=.5,xout=FALSE,outfun=outpro,STAND=TRUE,...){
#
# Leasts trimmed squares regression via the function ltsReg in the
# R package robustbase
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
library(robustbase)
temp=ltsReg(y~x,alpha=1-tr)
#coef=ltsReg(y~x)[8]$coefficients
coef=temp[8]$coefficients
res=temp[7]$raw.resid
list(coef=coef,residuals=res)
}

# ============================================================================
# ltsreg.2
# ============================================================================
ltsreg.2<-function(x,y,tr=.2,xout=FALSE,outfun=outpro,STAND=TRUE,...){
#
# Leasts trimmed squares regression via the function ltsReg in the
# R package robustbase
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
library(robustbase)
temp=ltsReg(y~x,alpha=1-tr)
coef=temp[8]$coefficients
res=temp[7]$raw.resid
list(coef=coef,residuals=res)
}

# ============================================================================
# DregG
# ============================================================================

#' Global Test for Comparing Two Dependent Regression Lines
#'
#' Tests the global null hypothesis that two dependent groups (e.g., repeated
#' measures) have identical regression parameters using a robust analog of
#' Hotelling's T-squared test with bootstrap covariance estimation.
#'
#' @param x1 Numeric vector or matrix of predictor variable(s) for the first group.
#' @param y1 Numeric vector of the dependent variable for the first group.
#' @param x2 Numeric vector or matrix of predictor variable(s) for the second group.
#' @param y2 Numeric vector of the dependent variable for the second group.
#' @param nullv Numeric vector specifying the null hypothesis value for the
#'   difference in parameters (default: \code{NULL}, which tests all differences = 0).
#' @param regfun Robust regression function to use (default: \code{tshdreg}).
#'   Common alternatives: \code{tsreg}, \code{chreg}, \code{ltsreg}, \code{opreg}.
#' @inheritParams common-params
#' @param plotit Logical. Included for compatibility with outlier detection but
#'   not used for plotting (default: \code{FALSE}).
#' @param pr Logical. If \code{TRUE}, prints warnings about duplicate Y values
#'   when using \code{tsreg} (default: \code{TRUE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{p.value}{Bootstrap p-value for the global test}
#'   \item{est.1}{Regression coefficient estimates for group 1}
#'   \item{est.2}{Regression coefficient estimates for group 2}
#'   \item{estimate.dif}{Difference in estimates (group 1 - group 2)}
#' }
#'
#' @details
#' This function implements a global test for comparing all regression parameters
#' simultaneously between two dependent groups. The test statistic is based on
#' a projection distance (similar to Mahalanobis distance) computed from the
#' differences in parameter estimates.
#'
#' The procedure:
#' 1. Estimates regression parameters for both groups
#' 2. Computes the difference in parameter vectors
#' 3. Uses bootstrap to estimate the covariance of the differences
#' 4. Computes a projection distance from the differences to the null value
#' 5. Compares observed distance to bootstrap distribution
#'
#' The bootstrap resamples rows (preserving pairing between measurements) and
#' refits both regressions. This provides a robust estimate of the sampling
#' distribution under dependence.
#'
#' When \code{xout=TRUE}, outliers detected in either group are removed from
#' both groups to preserve pairing.
#'
#' For ordinary least squares regression, use \code{\link{DregGOLS}} instead,
#' which provides a specialized implementation.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{difreg}} for parameter-by-parameter comparisons,
#'   \code{\link{DregGMC}} for parallel processing version,
#'   \code{\link{DregGOLS}} for OLS-specific version
#'
#' @examples
#' # Simulated pre-post intervention data
#' set.seed(123)
#' n <- 50
#' x1 <- rnorm(n)
#' y1 <- 1 + 2*x1 + rnorm(n)
#' x2 <- x1 + rnorm(n, sd=0.3)
#' y2 <- 1.5 + 2.5*x2 + rnorm(n)  # Changed parameters
#'
#' # Global test
#' result <- DregG(x1, y1, x2, y2, nboot=500)
#' print(result$p.value)
#'
#' # Multiple predictors
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- 1 + 2*x1[,1] - x1[,2] + rnorm(50)
#' x2 <- x1 + matrix(rnorm(100, sd=0.2), ncol=2)
#' y2 <- 1.5 + 2.2*x2[,1] - 1.1*x2[,2] + rnorm(50)
#' DregG(x1, y1, x2, y2, nboot=500)
#'
#' @export
DregG<-function(x1,y1,x2,y2,nullv=NULL,regfun=tshdreg,nboot=500,xout=FALSE,outfun=outpro,
SEED=TRUE,plotit=FALSE,pr=TRUE,...){
#
#  Global test that two dependent groups have identical
#  regression parameters.
#
#  Use a variation of Hotelling's test coupled with a bootstrap
#  estimate of the relevant covariance matrix associated with the differences
#  in the estimates of the parameters.#  For OLS, use DregGOLS
#
#  (plotit=F is used so that in simulations, if xout=T, the seed is not
#  set everytime outpro is called.)
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
if(is.null(nullv))nullv=rep(0,p1)
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
n=length(y1)
if(xout){
flag1=outfun(x1,...)$out.id
flag2=outfun(x2,...)$out.id
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
x1=as.matrix(x1)
x2=as.matrix(x2)
flagF=FALSE
flagF1=identical(regfun,tsreg)
flagF1[2]=identical(regfun,tshdreg)
#flagF1[3]=identical(regfun,tshdreg_C) obsolete,now it causes an error
if(sum(flagF1)>0)flagF=TRUE
if(!flagF){if(pr){
if(sum(duplicated(y1)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
pr=FALSE
}
if(pr){
if(sum(duplicated(y2)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
}}
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data,1,regboot,x1,y1,regfun=regfun,xout=FALSE,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-apply(data,1,regboot,x2,y2,regfun=regfun,xout=FALSE,...)
dif=t(bvec1-bvec2)
temp<-pdis(rbind(dif,nullv))
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
est1=regfun(x1,y1)$coef
est2=regfun(x2,y2)$coef
est=est1-est2
list(p.value=sig.level,est.1=est1,est.2=est2,estimate.dif = est)
}

# ============================================================================
# DregGMC
# ============================================================================

#' Global Test for Dependent Regression Lines (Parallel Processing)
#'
#' Tests the hypothesis that two dependent groups have identical regression
#' parameters using a variation of Hotelling's test with bootstrap covariance
#' estimation and parallel processing.
#'
#' @param x1 Numeric vector or matrix of predictor variable(s) for the first group.
#' @param y1 Numeric vector of the dependent variable for the first group.
#' @param x2 Numeric vector or matrix of predictor variable(s) for the second group.
#' @param y2 Numeric vector of the dependent variable for the second group.
#' @param nullv Numeric vector specifying the null hypothesis values for the
#'   difference in parameters (default: \code{NULL}, uses zeros for all parameters).
#' @param regfun Robust regression function to use (default: \code{tsreg}).
#'   Common alternatives: \code{tshdreg}, \code{chreg}, \code{ltsreg}, \code{opreg}.
#' @inheritParams common-params
#' @param plotit Logical; currently not used but included for compatibility (default: \code{FALSE}).
#' @param pr Logical. If \code{TRUE}, prints warnings about duplicate Y values
#'   when using \code{tsreg} (default: \code{TRUE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{p.value}{P-value from the global test}
#'   \item{est.1}{Regression coefficient estimates for group 1}
#'   \item{est.2}{Regression coefficient estimates for group 2}
#'   \item{estimate.dif}{Difference in estimates (group 1 - group 2)}
#' }
#'
#' @details
#' This function performs a global test that all regression parameters
#' (intercept and slopes) are equal between two dependent (paired) groups.
#' It is the multicore version of \code{\link{DregG}}.
#'
#' The method uses a multivariate extension of the percentile bootstrap approach
#' with a variation of Hotelling's T-squared statistic. The test statistic is
#' based on the Mahalanobis distance of the observed difference from the null
#' hypothesis difference, using a bootstrap estimate of the covariance matrix.
#'
#' **Procedure**:
#' 1. Bootstrap resampling preserves the pairing by resampling entire rows
#' 2. Computes regression estimates for both groups in each bootstrap sample
#' 3. Estimates the covariance matrix of the parameter differences
#' 4. Computes a projection distance-based test statistic
#' 5. P-value from comparing observed statistic to bootstrap distribution
#'
#' **When to use**: This is a global test that considers all parameters
#' jointly. If the global test is significant, use \code{\link{difregMC}}
#' to determine which specific parameters differ.
#'
#' The function uses \code{mclapply} for parallel processing, making it
#' faster than \code{\link{DregG}} for large datasets or many bootstrap samples.
#'
#' **Note**: For OLS regression with dependent groups, use \code{\link{DregGOLS}}.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press. Chapter 10.
#'
#' @seealso \code{\link{DregG}} for non-parallel version,
#'   \code{\link{difregMC}} for parameter-specific tests,
#'   \code{\link{DregGOLS}} for OLS version
#'
#' @examples
#' \dontrun{
#' # Simulated pre-post intervention data
#' set.seed(123)
#' n <- 50
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- 1 + 2*x1[,1] - x1[,2] + rnorm(n)
#' x2 <- x1 + matrix(rnorm(100, sd=0.3), ncol=2)
#' y2 <- 1.2 + 2.3*x2[,1] - 0.9*x2[,2] + rnorm(n)
#'
#' # Global test
#' result <- DregGMC(x1, y1, x2, y2, nboot=500)
#' print(result$p.value)
#' print(result$estimate.dif)
#'
#' # If significant, follow up with parameter-specific tests
#' if(result$p.value < 0.05) {
#'   difregMC(x1, y1, x2, y2, nboot=500)
#' }
#' }
#'
#' @export
DregGMC<-function(x1,y1,x2,y2,nullv=NULL,regfun=tsreg,nboot=500,xout=FALSE,outfun=outpro,
SEED=TRUE,plotit=FALSE,pr=TRUE,...){
#
#  Global test that two dependent groups have identical
#  regression parameters.
#
#  Use a variation of Hotelling's test coupled with a bootstrap
#  estimate of the relevant covariance matrix associated with the differences
#  in the estimates of the parameters.#  For OLS, use DregGOLS
#
#  (plotit=F is used so that in simulations, if xout=T, the seed is not
#  set everytime outpro is called.)
#
flag=FALSE
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
if(is.null(nullv))nullv=rep(0,p1)
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
n=length(y1)
if(xout){
flag1=outfun(x1,...)$out.id
flag2=outfun(x2,...)$out.id
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
flagF=FALSE
flagF1=identical(regfun,tsreg)
flagF1[2]=identical(regfun,tshdreg)
#flagF1[3]=identical(regfun,tshdreg_C)
if(sum(flagF1)>0)flagF=TRUE
if(!flagF){
if(pr){
if(sum(duplicated(y1)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
pr=FALSE
}
if(pr){
if(sum(duplicated(y2)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
}}
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec1=mclapply(data,regbootMC,x1,y1,regfun,xout=FALSE,...)
bvec1=matl(bvec1)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2=mclapply(data,regbootMC,x2,y2,regfun,xout=FALSE,...)
bvec2=matl(bvec2)
dif=t(bvec1-bvec2)
temp<-pdisMC(rbind(dif,nullv))
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
est1=regfun(x1,y1)$coef
est2=regfun(x2,y2)$coef
est=est1-est2
list(p.value=sig.level,est.1=est1,est.2=est2,estimate.dif = est)
}

# ============================================================================
# difregMC
# ============================================================================

#' Compare Two Dependent Regression Lines (Parallel Processing)
#'
#' Computes bootstrap confidence intervals for the difference in regression
#' parameters between two dependent groups using parallel processing.
#' This is the multicore version of \code{\link{difreg}}.
#'
#' @param x1 Numeric vector or matrix of predictor variable(s) for the first group.
#' @param y1 Numeric vector of the dependent variable for the first group.
#' @param x2 Numeric vector or matrix of predictor variable(s) for the second group.
#' @param y2 Numeric vector of the dependent variable for the second group.
#' @param regfun Robust regression function to use (default: \code{tsreg}).
#'   Common alternatives: \code{tshdreg}, \code{chreg}, \code{ltsreg}, \code{opreg}.
#' @inheritParams common-params
#' @param plotit Logical. If \code{TRUE}, creates a scatterplot with both
#'   regression lines (default: \code{FALSE}).
#' @param xlab Label for x-axis when \code{plotit=TRUE} (default: "X").
#' @param ylab Label for y-axis when \code{plotit=TRUE} (default: "Y").
#' @param pr Logical. If \code{TRUE}, prints warnings about duplicate Y values
#'   when using \code{tsreg} (default: \code{TRUE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Original sample size before any outlier removal}
#'   \item{n.keep}{Sample size after outlier removal (if \code{xout=TRUE})}
#'   \item{param}{Character vector labeling parameters (Intercept, slope1, ...)}
#'   \item{p.values}{P-values for testing difference = 0 for each parameter}
#'   \item{est.grp1}{Regression coefficient estimates for group 1}
#'   \item{est.grp2}{Regression coefficient estimates for group 2}
#'   \item{conf.intervals}{Matrix with columns: Param, ci.low, ci.hi}
#' }
#'
#' @details
#' This function compares regression parameters between two dependent (paired) groups
#' using a percentile bootstrap with parallel processing via \code{mclapply}.
#'
#' The bootstrap procedure resamples entire rows (preserving the dependency between
#' observations) and refits both regression models. This provides valid inference
#' for paired/repeated measures designs where the same subjects are measured at
#' two time points or under two conditions.
#'
#' P-values are computed from the proportion of bootstrap samples where the
#' difference has opposite sign from zero.
#'
#' When \code{xout=TRUE}, outliers detected by \code{outfun} in either group
#' are removed from both groups to preserve pairing.
#'
#' The function uses \code{mclapply} for parallel processing, making it faster than
#' \code{\link{difreg}} for large datasets or many bootstrap samples.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{difreg}} for non-parallel version,
#'   \code{\link{DregG}} for global test,
#'   \code{\link{reg2ciMC}} for independent groups
#'
#' @examples
#' \dontrun{
#' # Simulated pre-post intervention data
#' set.seed(123)
#' n <- 50
#' x1 <- rnorm(n)
#' y1 <- 1 + 2*x1 + rnorm(n)
#' x2 <- x1 + rnorm(n, sd=0.3)
#' y2 <- 1.5 + 2.5*x2 + rnorm(n)
#'
#' # Compare with parallel processing
#' result <- difregMC(x1, y1, x2, y2, nboot=500)
#' print(result$conf.intervals)
#'
#' # Multiple predictors
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- 1 + 2*x1[,1] - x1[,2] + rnorm(50)
#' x2 <- x1 + matrix(rnorm(100, sd=0.2), ncol=2)
#' y2 <- 1.5 + 2.2*x2[,1] - 1.1*x2[,2] + rnorm(50)
#' difregMC(x1, y1, x2, y2, nboot=500)
#' }
#'
#' @export
difregMC<-function(x1,y1,x2,y2,regfun=tsreg,xout=FALSE,outfun=outpro,nboot=599,
alpha=.05,SEED=TRUE,plotit=FALSE,xlab='X',ylab='Y',pr=TRUE,...){
#
# regression data from two different times i.e., two dependent groups
#
#  compute confidence interval for the difference in the slopes and intercepts
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
n=length(y1)
if(xout){
flag1=outfun(x1,...)$out.id
flag2=outfun(x2,...)$out.id
flag=unique(c(flag1,flag2))
X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
flagF=identical(regfun,tsreg)
if(flagF){
if(pr){
if(sum(duplicated(y1)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
pr=FALSE
}
if(pr){
if(sum(duplicated(y2)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
}}
nk=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec1<-mclapply(data,regboot,x1,y1,regfun,mc.preschedule=TRUE,xout=FALSE,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-mclapply(data,regboot,x2,y2,regfun,mc.preschedule=TRUE,xout=FALSE,...)
bvec1=matl(bvec1)
bvec2=matl(bvec2)
dif=t(bvec1)-t(bvec2)
dif.sort=apply(dif,2,sort)
pvec=NA
for(i in 1:p1){
pvec[i]<-(sum(dif[,i]<0)+.5*sum(dif[,i]==0))/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
}
pvec<-2*pvec
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=matrix(NA,nrow=p1,ncol=3)
ci[,1]=c(0:p)
for(i in 1:p1){
ci[i,2]=dif.sort[ilow,i]
ci[i,3]=dif.sort[ihi,i]
}
dimnames(ci)=list(NULL,c('Param','ci.low','ci.hi'))
if(plotit){
reg2plot(x1,y1,x2,y2,xlab=xlab,ylab=ylab,regfun=regfun,...)
}
lvec='Intercept'
for(j in 2:p1)lvec=c(lvec,paste('slope',j-1))
#pvec=array(pvec,dimnames=lvec)
est1=regfun(x1,y1,xout=FALSE,...)$coef
est2=regfun(x2,y2,xout=FALSE,...)$coef
list(n=n,n.keep=nk,param=lvec,p.values=pvec,est.grp1=est1,est.grp2=est2,conf.intervals=ci)
}

# ============================================================================
# Qreg
# ============================================================================

#' Quantile Regression
#'
#' Performs quantile regression using numerical optimization. Handles tied values
#' in the dependent variable better than \code{qreg}. Estimates the qth conditional
#' quantile of Y given X.
#'
#' @inheritParams common_params
#' @param q Quantile to estimate (default: 0.5 for median regression).
#'   Must be between 0 and 1.
#' @param res.vals Logical; if TRUE, returns residuals (default: TRUE).
#' @param plotit Logical; if TRUE and p=1, plots the data with the fitted
#'   quantile regression line (default: FALSE).
#' @param xlab X-axis label for plot (default: 'X').
#' @param ylab Y-axis label for plot (default: 'Y').
#' @param pch Plotting character (default: '*').
#' @param ... Additional arguments passed to \code{outfun}.
#'
#' @return A list with components:
#' \item{coef}{Vector of estimated coefficients (intercept, then slopes).}
#' \item{residuals}{Vector of residuals (if \code{res.vals=TRUE}), otherwise NULL.}
#'
#' @details
#' This function estimates the qth conditional quantile of the response variable
#' given the predictors using numerical optimization (BFGS method). It is more
#' robust to tied values in the dependent variable than \code{qreg}, which can
#' encounter computational issues with ties.
#'
#' When \code{q=0.5}, performs median regression (least absolute deviations).
#' Other quantiles allow estimation of the conditional distribution of Y|X.
#'
#' The function uses ordinary least squares estimates as starting values for the
#' optimization routine.
#'
#' @references
#' Koenker, R., & Bassett, G. (1978). Regression quantiles. \emph{Econometrica},
#' 46, 33-50.
#'
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' @seealso \code{\link{qreg}}, \code{\link{Qreghat}}, \code{\link{tsreg}}
#'
#' @examples
#' # Median regression (50th percentile)
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2 + 3*x + rnorm(100)
#' Qreg(x, y, q=0.5)
#'
#' # 90th percentile regression
#' Qreg(x, y, q=0.9)
#'
#' # With outlier removal
#' y[1:5] <- y[1:5] + 20
#' Qreg(x, y, q=0.5, xout=TRUE)
#'
#' @export
Qreg<-function(x,y,q=.5,xout=FALSE,outfun=outpro,res.vals=TRUE,plotit=FALSE,xlab='X',ylab='Y',pch='*',...){
#
# Quantile regression. Like the function qreg, but avoids computational
# problems that can arise when there are tied values among the dependent
# variable
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
xx=as.matrix(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
init=ols(x,y)$coef
v=optim(init,qfun,x=x,y=y,q=q,method='BFGS')$par
p1=ncol(x)+1
res=NULL
if(res.vals)res<-y-x%*%v[2:p1]-v[1]
if(ncol(x)==1){
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab,pch=pch)
abline(v)
}}
list(coef=v,residuals=res)
}

# ============================================================================
# regcits
# ============================================================================
regcits<-function(x,y,regfun=tshdreg,nboot=599,alpha=.05,SEED=TRUE,pr=TRUE,
xout=FALSE,outfun=outpro,plotit=FALSE,xlab='Predictor 1',ylab='Predictor 2',
MC=TRUE,...){
if(MC)v=regciMC(x,y,regfun=regfun,nboot=nboot,alpha=alpha,SEED=SEED,pr=pr,xout=xout,
outfun=outfun,plotit=plotit,xlab=xlab,ylab=ylab,...)
if(!MC)v=regci(x,y,regfun=regfun,nboot=nboot,alpha=alpha,SEED=SEED,pr=pr,xout=xout,
outfun=outfun,plotit=plotit,xlab=xlab,ylab=ylab,...)
v
}

# ============================================================================
# scorreg
# ============================================================================
corregci.sub<-function(isub,x,y,corfun){
p=ncol(x)
xmat<-matrix(x[isub,],nrow(x),ncol(x))
e=NA
for(j in 1:p)e[j]=corfun(xmat[,j],y[isub])$cor
e
}

# ============================================================================
# LMSreg
# ============================================================================

#' Least Median of Squares Regression
#'
#' Computes the least median of squares (LMS) regression estimator, which has
#' the highest possible breakdown point (50%) but lower efficiency than LTS.
#' The LMS estimator minimizes the median of squared residuals.
#'
#' @param x A numeric matrix of predictor variables (n by p).
#' @param y A numeric vector of the dependent variable (length n).
#' @inheritParams opreg
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Original sample size.}
#'   \item{n.keep}{Sample size after outlier removal (if \code{xout=TRUE}).}
#'   \item{coef}{Vector of regression coefficients (intercept, slopes).}
#' }
#'
#' @details
#' **Least Median of Squares (LMS)**:
#'
#' The LMS estimator finds regression coefficients that minimize:
#' \deqn{median(r_i^2)}
#' where \eqn{r_i} are the residuals.
#'
#' **Properties**:
#' - **Breakdown point**: 50% (highest possible)
#' - **Efficiency**: Lower than LTS or MM-regression
#' - **Resistance**: Excellent resistance to outliers and leverage points
#' - **Computation**: Uses \code{MASS::lmsreg} internally
#'
#' **When to use LMS**:
#' - Maximum robustness is paramount
#' - Heavy contamination suspected (up to 50%)
#' - Initial screening of data
#' - As a starting point for other robust methods
#'
#' **Comparison with LTS**: LTS (least trimmed squares) typically has better
#' efficiency while maintaining high breakdown. LMS is more extreme in its
#' resistance to outliers.
#'
#' @references
#' Rousseeuw, P.J. (1984). Least median of squares regression. Journal of the
#' American Statistical Association, 79, 871-880.
#'
#' @family robust regression functions
#' @seealso \code{\link{ltsreg}}, \code{\link{MMreg}}
#' @export
#' @examples
#' \dontrun{
#' # Simple LMS regression
#' x <- matrix(rnorm(50), ncol=1)
#' y <- 2*x + rnorm(50, sd=0.5)
#' LMSreg(x, y)
#'
#' # Heavily contaminated data
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50, sd=0.3)
#' y[1:20] <- y[1:20] + 10  # 20% contamination
#' fit <- LMSreg(x, y)
#' print(fit$coef)
#' }
LMSreg<-function(x,y,xout=FALSE,outfun=outpro,...){
#
#  Least median of squares
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
n=nrow(x)
n.keep=n
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
n.keep=nrow(x)
}
a=lmsreg(x,y)
list(n=n,n.keep=n.keep,coef=a[3]$coefficients)
}

# ============================================================================
# Qreghat
# ============================================================================
Qreghat<-function(x,y,xr=x,q=.5,xout=FALSE,outfun=outpro,plotit.pts=FALSE){
#
#
#
xy=elimna(cbind(x,y))
xr=as.matrix(xr)
x=as.matrix(x)
p=ncol(x)
p1=p+1
vals=NA
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
xr=as.matrix(xr)
est=Qreg(x,y,q=q)$coef
if(ncol(xr)!=p)xr=t(xr)   # for a single point, need to transpose.
yhat=est[1]+xr%*%est[2:p1]
if(plotit.pts)points(xr,yhat)
yhat
}

# ============================================================================
# reg.reglev
# ============================================================================
reg.reglev<-function(x,y,plotit=TRUE,xlab='X',ylab='Y',GEN=TRUE,regfun=tsreg,outfun=outpro,pr=TRUE,...){

#
#  Remove any bad leverage points detected by
#  the fit using the estimator indicated by regun
#
# GEN=TRUE: use a generalization of the Rousseeuw van Zomeren method
# GEN=FALSE: usw the Rousseeuw van Zomeren method. Unknown when if ever this older approach
#     offers an advantage.
#
xy=elimna(cbind(x,y))
n=nrow(xy)
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
x<-as.matrix(x)
keep=c(1:n)
if(!GEN)a=reglev(x,y,plotit=FALSE,SEED=FALSE)$bad.lev.points
else a=reglev.gen(x,y,plotit=FALSE,regfun=regfun,outfun=outfun)$bad.lev
if(length(a)>0)keep=keep[-a]
nk=length(y[keep])
e=regfun(x[keep,],y[keep],...)
list(n=n,n.keep=nk,coef=e$coef)
}

