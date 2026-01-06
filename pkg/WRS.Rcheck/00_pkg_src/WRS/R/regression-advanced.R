# WRS Package - Advanced Regression Methods
# Extracted from Rallfun-v45.R
#
# This module contains advanced regression methods including:
#  - Quantile regression smoothers (qhdsm, qhdsm2g, etc.)
#  - Smoothing methods (smean, smeancr, etc.)
#  - Logistic regression (logreg, logreg.P.ci, etc.)
#  - Multivariate/multilevel regression (mlrreg, mulgreg, MULMreg)
#  - K-nearest neighbors regression (KNNreg)
#  - GAM-related methods (gamindt, gamplot, etc.)
#  - Regression inference methods (regYci, regYband, etc.)
#  - Mediation, PCA, random forest regression
#  - Instrumental variables regression (regIV*)
#  - Specialized regression utilities
#
# Total functions: 75
# Extraction date: 2025-12-30


# ============================================================================
# khomreg
# ============================================================================

#' Test for Homoscedasticity in Linear Regression
#'
#' Tests the hypothesis that the error term in a linear regression model is
#' homoscedastic using a modification of the Cook-Weisberg statistic derived
#' by Koenker. This is a robust test for constant variance of residuals.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values. Can have
#'   multiple columns for multiple predictors.
#' @param y A numeric vector containing the response variable.
#'
#' @return A list with components:
#'   \item{test}{The test statistic (chi-square distributed)}
#'   \item{p.value}{The p-value for the test of homoscedasticity}
#'
#' @details
#' This function tests H0: homoscedastic errors using Koenker's modification
#' of the Cook-Weisberg statistic. The test statistic follows a chi-square
#' distribution with 1 degree of freedom under the null hypothesis.
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' testing using the function specified by `outfun`.
#'
#' @references
#' Lyon, J.D., & Tsai, C.L. (1996). A comparison of tests for heteroscedasticity.
#' The Statistician, 45, 337-349.
#'
#' @seealso \code{\link{tsreg}}, \code{\link{ltsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test for homoscedasticity
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 2)
#' y <- 2 + 3*x[,1] - x[,2] + rnorm(50)
#' khomreg(x, y)
#' }
khomreg<-function(x,y,xout=FALSE,outfun=out,...){
#
# Test hypothesis that error term in a linear regression model
# is homoscedastic using modification of Cook-Weisberg
# statistic derived by Koenker;
# See Lyon and Tsai, 1996, Statistician, 45, 337-349
#
x<-as.matrix(x)
if(xout){
flag<-outfun(x,...)$keep
x<-as.matrix(x)
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
pv<-ncol(x)
pv1<-pv+1
m<-cbind(x,y)
m<-elimna(m)
x<-m[,1:pv]
x<-as.matrix(x)
y<-m[,pv1]
dvec<-NA
mat<-matrix(nrow=nrow(x),ncol=pv)
temp<-lsfit(x,y)
sigest<-mean(temp$res^2)
dvec<-y-temp$res
dbar<-dvec-mean(dvec)
uval<-temp$res^2
uval<-as.matrix(uval)
test<-t(uval)%*%dbar%*%solve(t(dbar)%*%dbar)%*%t(dbar)%*%uval
psihat<-mean((temp$res^2-sigest)^2)
test<-test/psihat
p.value<-1-pchisq(test,1)
list(test=test,p.value=p.value)
}


# ============================================================================
# mgvfreg
# ============================================================================

#' Regression with MGV Outlier Detection
#'
#' Performs regression after removing outliers detected using the faster
#' inward Minimum Generalized Variance (MGV) method. Fits regression only
#' on points not flagged as outliers.
#'
#' @inheritParams common-params
#' @param x A numeric vector or matrix containing predictor values.
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use on non-outlier points
#'   (default: `tsreg` for Theil-Sen regression). Can be any regression
#'   function that returns `coef` and `residuals` components.
#' @param outfun Outlier detection function to apply to the MGV method
#'   (default: `outbox`). Used within `outmgvf()`.
#' @param plotit Logical. If `TRUE` and data is bivariate, plots the data
#'   and overlays the regression line (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{coef}{Regression coefficients from the regression on non-outlier points.}
#'   \item{residuals}{Residuals from the regression fit.}
#'
#' @details
#' This function combines outlier detection with robust regression:
#' \enumerate{
#'   \item Identifies outliers using the faster inward MGV method (`outmgvf`)
#'   \item Removes flagged outliers from the dataset
#'   \item Fits the specified regression function to the remaining points
#' }
#'
#' The MGV (Minimum Generalized Variance) method provides a computationally
#' efficient approach to multivariate outlier detection, making it suitable
#' for larger datasets or when speed is important.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note Requires the `SEED` variable to be defined in the calling environment
#'   for reproducible outlier detection.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{mgvreg}}, \code{\link{outmgvf}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Regression with MGV outlier detection
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' # Add outliers
#' x[1:3] <- c(5, 5, -5)
#' y[1:3] <- c(-5, 5, 5)
#' result <- mgvfreg(x, y)
#' result$coef
#' }
mgvfreg<-function(x,y,regfun=tsreg,outfun=outbox,plotit=TRUE){
#
# Do regression on points not labled outliers
# Use the faster inward mgv method
#
m<-cbind(x,y)
m<-elimna(m) # elminate any rows with missing data
flag<-outmgvf(m,outfun=outfun,plotit,SEED=SEED)$out.id
ivec<-rep(TRUE,nrow(m))
ivec[flag]<-FALSE
x<-as.matrix(x)
temp<-regfun(x[ivec,],y[ivec])
coef<-temp$coef
if(plotit && ncol(m)==2)abline(coef)
residuals<-temp$residuals
list(coef=coef,residuals=residuals)
}


# ============================================================================
# smeancr
# ============================================================================

#' Test Hypothesis About Multivariate Skipped Mean
#'
#' Tests the hypothesis that the multivariate skipped mean equals a specified
#' null value using bootstrap methods. The skipped mean is computed after
#' removing outliers via projection-based detection.
#'
#' @param m An n-by-p numeric matrix containing the multivariate data.
#' @param nullv Numeric vector of length p specifying the null hypothesis values
#'   (default: vector of zeros).
#' @param cop Integer (1-4) specifying the center for projection-based outlier detection:
#'   \itemize{
#'     \item 1 = Donoho-Gasko median
#'     \item 2 = MCD (Minimum Covariance Determinant)
#'     \item 3 = Marginal medians (default)
#'     \item 4 = MVE (Minimum Volume Ellipsoid)
#'   }
#' @inheritParams common-params
#' @param FAST Logical. If `TRUE`, uses faster depth-based outlier detection;
#'   if `FALSE`, uses standard projection method (default: `FALSE`).
#' @param nboot Number of bootstrap samples (default: 500).
#' @param plotit Logical. If `TRUE` and p=2, creates a plot showing bootstrap
#'   distribution and confidence region (default: `TRUE`).
#' @param xlab,ylab Axis labels for bivariate plot.
#' @param STAND Logical. If `TRUE`, standardizes data before projection
#'   (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{p.value}{Bootstrap p-value for testing H0: skipped mean = nullv.}
#'   \item{crit.level}{Critical level used (adjusted based on sample size).}
#'
#' @details
#' This function tests whether the multivariate skipped mean equals a specified
#' null value using a bootstrap approach:
#' \enumerate{
#'   \item For each bootstrap sample, compute the skipped mean (mean after
#'     outlier removal via projection method)
#'   \item Compute projection distance of each bootstrap estimate from the null value
#'   \item p-value is the proportion of bootstrap distances less than the distance
#'     of the null value
#' }
#'
#' For bivariate data (p=2) with `plotit=TRUE`, creates a scatterplot showing:
#' \itemize{
#'   \item Bootstrap estimates of the skipped mean
#'   \item The estimated skipped mean (marked with +)
#'   \item Confidence region (convex hull of central bootstrap estimates)
#' }
#'
#' The critical level is automatically adjusted based on sample size to maintain
#' approximately 0.05 significance level.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{smean}}, \code{\link{smeancrv2}}, \code{\link{smean2v2}},
#'   \code{\link{outpro}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test if bivariate skipped mean equals (0,0)
#' set.seed(123)
#' m <- matrix(rnorm(100), ncol=2)
#' m <- rbind(m, matrix(c(5,5,-5,-5), ncol=2))  # Add outliers
#' smeancr(m)
#'
#' # Test different null value
#' smeancr(m, nullv=c(0.5, 0.5))
#' }
smeancr<-function(m,nullv=rep(0,ncol(m)),cop=3,MM=FALSE,SEED=TRUE,FAST=FALSE,
nboot=500,plotit=TRUE,xlab="VAR 1",ylab="VAR 2",STAND=TRUE){
#
# m is an n by p matrix
#
# Test hypothesis that multivariate skipped estimators
# are all equal to the null value, which defaults to zero.
# The level of the test is .05.
#
# Eliminate outliers using a projection method
# That is, determine center of data using:
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
# cop=4 MVE
#
# For each point
# consider the line between it and the center
# project all points onto this line, and
# check for outliers using
#
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# Eliminate any outliers and compute means
#  using remaining data.
#
if(SEED)set.seed(2)
#if(!is.na(SEED))set.seed(SEED)
m<-elimna(m)
n<-nrow(m)
crit.level<-.05
if(n<=120)crit.level<-.045
if(n<=80)crit.level<-.04
if(n<=60)crit.level<-.035
if(n<=40)crit.level<-.03
if(n<=30)crit.level<-.025
if(n<=20)crit.level<-.02
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
val<-matrix(NA,ncol=ncol(m),nrow=nboot)
for(j in 1: nboot){
mm<-m[data[j,],]
if(FAST)temp<-outpro.depth(mm,plotit=FALSE,SEED=FALSE)$keep
if(!FAST)temp<-outpro(mm,plotit=FALSE,cop=cop,STAND=STAND)$keep
val[j,]<-apply(mm[temp,],2,mean)
}
temp<-pdis(rbind(val,nullv))
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
if(ncol(m)==2 && plotit){
plot(val[,1],val[,2],xlab=xlab,ylab=ylab)
temp3<-smean(m,cop=cop,STAND=STAND)
points(temp3[1],temp3[2],pch="+")
ic<-round((1-crit.level)*nboot)
temp<-pdis(val)
temp.dis<-order(temp)
xx<-val[temp.dis[1:ic],]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(p.value=sig.level,crit.level=crit.level)
}


# ============================================================================
# gamplot
# ============================================================================

#' Plot Regression Surface Using Generalized Additive Model
#'
#' Fits and plots a regression surface using a generalized additive model (GAM)
#' with optional spline smoothing. Supports up to 4 predictors with automatic
#' visualization for 1D and 2D cases.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values (maximum 4 columns).
#' @param y A numeric vector containing the response variable.
#' @param sop Logical. If `TRUE`, uses spline smoothing (s()) in the GAM;
#'   if `FALSE`, uses linear terms (default: `TRUE`).
#' @param pyhat Logical. If `TRUE`, returns predicted values; if `FALSE`,
#'   returns "Done" message (default: `FALSE`).
#' @param theta,phi Viewing angles for 3D perspective plot (defaults: 50, 25).
#' @param expand Expansion factor for 3D plot (default: 0.5).
#' @param scale Logical. If `TRUE`, scales the 3D plot (default: `TRUE`).
#' @param ticktype Type of tick marks for 3D plot (default: "simple").
#' @param xlab,ylab,zlab Axis labels for plots.
#'
#' @return If `pyhat = TRUE`, returns fitted values from the GAM. If
#'   `pyhat = FALSE`, returns "Done". Creates plots as a side effect.
#'
#' @details
#' This function fits a generalized additive model using the `mgcv` package:
#' - With `sop = TRUE`: Uses smoothing splines s(x) for each predictor
#' - With `sop = FALSE`: Uses linear terms for each predictor
#'
#' The function automatically handles visualization:
#' - 1 predictor: Line plot of fitted values vs. x
#' - 2 predictors: 3D perspective plot of fitted surface
#' - 3-4 predictors: No automatic plot, but returns fitted values
#'
#' Outlier removal options:
#' - `xout = TRUE`: Remove outliers from predictor space
#' - `eout = TRUE`: Remove outliers from combined (x,y) space
#' - Cannot have both `xout` and `eout` set to `TRUE`
#'
#' @note Requires the `mgcv` package for GAM fitting and `akima` package for
#'   3D surface interpolation.
#'
#' @references
#' Wood, S.N. (2017). Generalized Additive Models: An Introduction with R (2nd ed.).
#' Chapman and Hall/CRC.
#'
#' @seealso \code{\link{gamplotv2}}, \code{\link{gamplotINT}}, \code{\link{gamindt}},
#'   \code{\link[mgcv]{gam}}
#'
#' @export
#' @examples
#' \dontrun{
#' # GAM with single predictor
#' set.seed(123)
#' x <- seq(-3, 3, length=100)
#' y <- sin(x) + rnorm(100, sd=0.2)
#' gamplot(x, y)
#'
#' # GAM with two predictors
#' x <- matrix(rnorm(200), ncol=2)
#' y <- x[,1]^2 + x[,2]^2 + rnorm(100)
#' gamplot(x, y)
#' }
gamplot<-function(x,y,sop=TRUE,pyhat=FALSE,eout=FALSE,xout=FALSE,outfun=out,plotit=TRUE,
xlab="X",ylab="",zlab="",theta=50,phi=25,expand=.5,scale=TRUE,ticktype="simple"){
#
# Plot regression surface using generalized additive model
#
# sop=F, use usual linear model y~x1+x2...
# sop=T, use splines
#
library(mgcv)
x<-as.matrix(x)
np<-ncol(x)
np1<-np+1
if(ncol(x)>4)stop("x should have at most four columns of data")
m<-elimna(cbind(x,y))
x<-m[,1:np]
x<-as.matrix(x)
y<-m[,np1]
if(xout && eout)stop("Can't have xout=eout=T")
if(eout){
flag<-outfun(m)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
}
x<-m[,1:np]
x<-as.matrix(x)
y<-m[,np1]
if(!sop){
if(ncol(x)==1)fitr<-fitted(gam(y~x[,1]))
if(ncol(x)==2)fitr<-fitted(gam(y~x[,1]+x[,2]))
if(ncol(x)==3)fitr<-fitted(gam(y~x[,1]+x[,2]+x[,3]))
if(ncol(x)==4)fitr<-fitted(gam(y~x[,1]+x[,2]+x[,3]+x[,4]))
}
if(sop){
if(ncol(x)==1)fitr<-fitted(gam(y~s(x[,1])))
if(ncol(x)==2)fitr<-fitted(gam(y~s(x[,1])+s(x[,2])))
if(ncol(x)==3)fitr<-fitted(gam(y~s(x[,1])+s(x[,2])+s(x[,3])))
if(ncol(x)==4)fitr<-fitted(gam(y~s(x[,1])+s(x[,2])+s(x[,3])+s(x[,4])))
}
last<-fitr
if(plotit){
if(ncol(x)==1){
plot(x,fitr,xlab=xlab,ylab=ylab)
}
if(ncol(x)==2){
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
fitr<-interp(mkeep[,1],mkeep[,2],fitr)
persp(fitr,theta=theta,phi=phi,expand=expand,xlab=xlab,ylab=ylab,zlab=zlab,
scale=scale,ticktype=ticktype)
}
}
if(!pyhat)last <- "Done"
last
}


# ============================================================================
# mulgreg
# ============================================================================

#' Multivariate Regression via Robust Covariance
#'
#' Performs multivariate linear regression using robust covariance-based
#' estimation following Rousseeuw et al. (2004). Estimates regression
#' coefficients from robust location and scatter estimates.
#'
#' @param x A numeric matrix containing predictor values (independent variables).
#' @param y A numeric matrix containing response values (dependent variables).
#'   Must be multivariate (matrix, not vector).
#' @param cov.fun Covariance function for robust estimation (default: `rmba`
#'   for Median Ball Algorithm). Must return components `$center` and `$cov`.
#'
#' @return A list with components:
#'   \item{coef}{Matrix of regression coefficients. First row contains intercepts
#'     (one per response variable), subsequent rows contain slopes.}
#'   \item{residuals}{Matrix of residuals (observed - predicted) with same
#'     dimensions as y.}
#'
#' @details
#' This function implements a simpler variant of the Rousseeuw et al. (2004)
#' multivariate regression method. Unlike `mlrreg()`, this function:
#' \itemize{
#'   \item Computes a single robust covariance matrix for (X, Y)
#'   \item Derives regression coefficients directly from the covariance structure
#'   \item Does not perform iterative outlier removal
#' }
#'
#' The regression coefficients are estimated as:
#' \deqn{B = \Sigma_{XX}^{-1} \Sigma_{XY}}
#' where \eqn{\Sigma_{XX}} and \eqn{\Sigma_{XY}} are robust estimates of the
#' covariance matrices.
#'
#' The choice of `cov.fun` determines the robustness properties:
#' - `rmba`: Median Ball Algorithm (default, good breakdown point)
#' - `covmcd`: Minimum Covariance Determinant (high breakdown)
#' - `covogk`: Orthogonalized Gnanadesikan-Kettenring
#' - `skipcov`: Skipped covariance (projection-based outlier removal)
#'
#' @references
#' Rousseeuw, P.J., Van Aelst, S., Van Driessen, K., & Agulló, J. (2004).
#' Robust multivariate regression. Technometrics, 46, 293-305.
#'
#' @seealso \code{\link{mlrreg}}, \code{\link{MULMreg}}, \code{\link{COVreg}},
#'   \code{\link{Mreglde}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Multivariate regression with 2 predictors and 2 responses
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- matrix(NA, nrow=50, ncol=2)
#' y[,1] <- 1 + 2*x[,1] - x[,2] + rnorm(50)
#' y[,2] <- -1 + x[,1] + 3*x[,2] + rnorm(50)
#' result <- mulgreg(x, y)
#' result$coef  # First row = intercepts, remaining rows = slopes
#' }
mulgreg<-function(x,y,cov.fun=rmba){
#
# Do Multivariate regression in Rousseeuw, Van Aelst, Van Driessen Agullo
# (2004) Technometrics, 46, 293-305
#
# (y can be multivariate)
#
if(!is.matrix(y))stop("y is not a matrix")
X<-cbind(x,y)
X<-elimna(X)
qy<-ncol(y)
qx<-ncol(x)
qxp1<-qx+1
tqyqx<-qy+qx
y<-X[,qxp1:tqyqx]
# compute initial estimate of slopes and intercept:
locscat<-cov.fun(X)
sig<-locscat$cov
mu<-locscat$center
sigxx<-sig[1:qx,1:qx]
sigxy<-sig[1:qx,qxp1:tqyqx]
sigyy<-sig[qxp1:tqyqx,qxp1:tqyqx]
Bhat<-solve(sigxx)%*%sigxy
sige<-sigyy-t(Bhat)%*%sigxx%*%Bhat
sige.inv<-solve(sige)
Ahat<-t(mu[qxp1:tqyqx]-t(Bhat)%*%mu[1:qx])
resL<-matrix(nrow=nrow(X),ncol=qy)
for(i in 1:nrow(X))resL[i,]<-y[i,]-t(Bhat)%*%X[i,1:qx]
for(j in 1:qy)resL[,j]<-resL[,j]-Ahat[j]
list(coef=rbind(Ahat,Bhat),residuals=resL)
}


# ============================================================================
# regpord.sub
# ============================================================================

#' Bootstrap Helper for regpord
#'
#' Internal helper function used by \code{\link{regpord}} to compute generalized
#' variance for bootstrap samples when assessing predictor importance.
#'
#' @param isub Bootstrap sample indices.
#' @param x Matrix of predictors.
#' @param y Response vector.
#' @param cov.fun Covariance function for generalized variance.
#'
#' @return Numeric vector of generalized variance values from \code{regvarp()}.
#'
#' @keywords internal
#' @export
regpord.sub<-function(isub,x,y,cov.fun){
xmat<-matrix(x[isub,],nrow(x),ncol(x))
vals<-regvarp(xmat,y[isub],cov.fun=cov.fun)
vals
}


# ============================================================================
# regvarp
# ============================================================================

#' Measure Importance of Predictors in Regression
#'
#' Assesses the importance of individual predictors or predictor combinations
#' in a regression problem by computing generalized variance for each predictor
#' or p-predictor subset.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values (must have multiple columns).
#' @param y A numeric vector containing the response variable.
#' @param p Integer specifying the number of predictors to consider simultaneously
#'   (default: 1 for individual predictors). If p > 1, evaluates all p-way
#'   predictor combinations.
#' @param locfun Location function for standardizing predictors (default: `lloc`).
#' @param scat Scatter function for standardizing predictors (default: `var`).
#' @param est Estimator function for standardizing predictors (default: `mean`).
#' @param cov.fun Covariance function to use for generalized variance computation
#'   (default: `cov.mba` for Median Ball Algorithm covariance).
#'
#' @return A numeric vector containing generalized variance values:
#'   \itemize{
#'     \item If `p=1`: One value per predictor
#'     \item If `p>1`: One value per p-way predictor combination
#'   }
#'   Higher values indicate stronger association with the response.
#'
#' @details
#' This function evaluates predictor importance by:
#' \enumerate{
#'   \item Standardizing all predictors using specified location, scatter, and
#'     estimator functions
#'   \item For each predictor or predictor combination, computing the generalized
#'     variance (determinant of robust covariance matrix) of the response with
#'     the predictor(s)
#'   \item Returning these variance measures as importance metrics
#' }
#'
#' The generalized variance provides a multivariate measure of association that
#' accounts for the covariance structure. Lower generalized variance indicates
#' stronger association (less unexplained variability).
#'
#' Missing values are automatically removed before analysis.
#'
#' @note Used internally by `regpord()` and `regpord.sub()` for comparing
#'   predictor strengths.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regpord}}, \code{\link{regIVstr}}, \code{\link{gvarg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Measure importance of individual predictors
#' set.seed(123)
#' x <- matrix(rnorm(100*3), ncol=3)
#' y <- 2*x[,1] + 0.5*x[,2] + rnorm(100)
#' regvarp(x, y, p=1)  # Evaluate each predictor
#'
#' # Evaluate all 2-way combinations
#' regvarp(x, y, p=2)
#' }
regvarp<-function(x,y,p=1,locfun=lloc,scat=var,est=mean,cov.fun=cov.mba){
#
# Measure the importance of each of p variables in a regression
# problem, p>1
#
xy=cbind(x,y)
xy<-elimna(xy)
m<-ncol(x)
x=xy[,1:m]
n<-nrow(x)
m1=m+1
y=xy[,m1]
x=standm(x,locfun=locfun,est=est,scat=scat)
vals=NA
if(p==1)for(j in 1:m){
vals[j]=gvarg(cbind(y,x[,j]),cov.fun)
}
if(p>1){
temp=modgen(m)
ic=0
for(j in 1:length(temp)){
if(length(temp[[j]])==p){
ic=ic+1
vals[ic]=gvarg(cbind(y,x[,temp[[j]]]),cov.fun)
z=cbind(y,x[,temp[[j]]])
}}}
vals
}


# ============================================================================
# smean2v2
# ============================================================================

#' Two-Group Comparison of Multivariate Skipped Means
#'
#' Tests the hypothesis that two independent groups have equal multivariate
#' skipped means using bootstrap methods. The skipped mean is computed after
#' removing outliers via projection-based detection.
#'
#' @param m1 An n1-by-p numeric matrix containing data for group 1.
#' @param m2 An n2-by-p numeric matrix containing data for group 2 (must have
#'   same number of columns as m1).
#' @param nullv Numeric vector of length p specifying the null hypothesis
#'   difference (default: vector of zeros, testing equality).
#' @param cop Integer (1-4) specifying the center for projection-based outlier detection:
#'   \itemize{
#'     \item 1 = Donoho-Gasko median
#'     \item 2 = MCD (Minimum Covariance Determinant)
#'     \item 3 = Marginal medians (default)
#'     \item 4 = MVE (Minimum Volume Ellipsoid)
#'   }
#' @inheritParams common-params
#' @param nboot Number of bootstrap samples (default: 500).
#' @param plotit Logical. If `TRUE` and p=2, creates a plot showing bootstrap
#'   distribution of the difference (default: `TRUE`).
#' @param STAND Logical. If `TRUE`, standardizes data before projection
#'   (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{Est.1}{Skipped mean for group 1.}
#'   \item{Est.2}{Skipped mean for group 2.}
#'   \item{dif}{Difference in skipped means (group 1 - group 2).}
#'   \item{p.value}{Bootstrap p-value for testing H0: difference = nullv.}
#'   \item{crit.level}{Critical level used (adjusted based on sample size).}
#'
#' @details
#' This function compares two independent groups using skipped means:
#' \enumerate{
#'   \item Compute skipped means for each group (mean after outlier removal)
#'   \item For each bootstrap sample, resample from pooled data and compute
#'     the difference in skipped means
#'   \item Compute projection distance of each bootstrap difference from null value
#'   \item p-value is based on comparing observed difference to bootstrap distribution
#' }
#'
#' For bivariate data (p=2) with `plotit=TRUE`, creates a scatterplot showing:
#' \itemize{
#'   \item Bootstrap estimates of the difference
#'   \item Confidence region (convex hull of central bootstrap estimates)
#' }
#'
#' The critical level is automatically adjusted based on sample size to maintain
#' approximately 0.05 significance level.
#'
#' @note Groups must have the same number of variables (columns).
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{smean}}, \code{\link{smeancr}}, \code{\link{smeancrv2}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two groups on bivariate data
#' set.seed(123)
#' m1 <- matrix(rnorm(50*2, mean=0), ncol=2)
#' m2 <- matrix(rnorm(50*2, mean=0.5), ncol=2)
#' # Add outliers
#' m1 <- rbind(m1, matrix(c(5,5), ncol=2))
#' m2 <- rbind(m2, matrix(c(-5,-5), ncol=2))
#' smean2v2(m1, m2)
#' }
smean2v2<-function(m1,m2,nullv=rep(0,ncol(m1)),cop=3,MM=FALSE,SEED=TRUE,
nboot=500,plotit=TRUE,MC=FALSE,STAND=TRUE){
#
# m is an n by p matrix
#
# For two independent groups,
# test hypothesis that multivariate skipped estimators
# are all equal.
#
# The level of the test is .05.
#
# Skipped estimator is used, i.e.,
# eliminate outliers using a projection method.
# That is, determine center of data using:
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
# cop=4 MVE
#
# For each point
# consider the line between it and the center,
# project all points onto this line, and
# check for outliers using
#
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection is an outlier
# using a modification of the usual boxplot rule.
#
# Eliminate any outliers and compute means
#  using remaining data.
#
if(ncol(m1) != ncol(m2)){
stop("Number of variables in group 1 does not equal the number in group 2.")
}
if(SEED)set.seed(2)
m1<-elimna(m1)
m2<-elimna(m2)
n1<-nrow(m1)
n2<-nrow(m2)
n<-min(c(n1,n2))
crit.level<-.05
if(n<=120)crit.level<-.045
if(n<=80)crit.level<-.04
if(n<=60)crit.level<-.035
if(n<=40)crit.level<-.03
if(n<=30)crit.level<-.025
if(n<=20)crit.level<-.02
val<-matrix(NA,ncol=ncol(m1),nrow=nboot)
est1=smean(m1)
est2=smean(m2)
#est=smean(m1)-smean(m2)
est=est1-est2
for(j in 1: nboot){
data1<-sample(n1,size=n1,replace=TRUE)
data2<-sample(n2,size=n2,replace=TRUE)
mm1<-m1[data1,]
temp<-outpro(mm1,plotit=FALSE,cop=cop,STAND=STAND)$keep
v1<-apply(mm1[temp,],2,mean)
mm2<-m2[data2,]
temp<-outpro(mm2,plotit=FALSE,cop=cop,STAND=STAND)$keep
v2<-apply(mm2[temp,],2,mean)
val[j,]<-v1-v2
}
if(!MC)temp<-pdis(rbind(val,nullv))
if(MC)temp<-pdisMC(rbind(val,nullv))
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
if(ncol(m1)==2 && plotit){
plot(val[,1],val[,2],xlab="VAR 1",ylab="VAR 2")
if(!MC)temp3<-smean(m1,cop=cop)-smean(m2,cop=cop)
if(MC)temp3<-smeanMC(m1,cop=cop)-smeanMC(m2,cop=cop)
points(temp3[1],temp3[2],pch="+")
ic<-round((1-crit.level)*nboot)
if(!MC)temp<-pdis(val)
if(MC)temp<-pdisMC(val)
temp.dis<-order(temp)
xx<-val[temp.dis[1:ic],]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(Est.1=est1,Est.2=est2,dif=est,p.value=sig.level,crit.level=crit.level)
}


# ============================================================================
# regpord
# ============================================================================

#' Compare Strength of Association for Standardized Predictors
#'
#' Compares the strength of association between the response and each of
#' multiple predictors using robust covariance-based generalized variance,
#' with predictors standardized to enable fair comparison.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values (must have 2 or more columns).
#' @param y A numeric vector containing the response variable.
#' @param nboot Number of bootstrap samples for standard error estimation (default: 100).
#' @param cov.fun Covariance function to use for generalized variance computation
#'   (default: `cov.mba` for Median Ball Algorithm).
#' @param pr Logical. If `TRUE`, prints progress messages (default: `TRUE`).
#' @param plotit Logical. If `TRUE` and there are exactly 2 predictors, creates
#'   a scatterplot distinguishing the two predictors (default: `TRUE`).
#' @param xlab,ylab Axis labels for the plot.
#' @param est Estimator for standardizing predictors (default: `mean`).
#' @param scat Scatter measure for standardizing predictors (default: `var`).
#'
#' @return A list with components:
#'   \item{crit.value}{Critical value for significance testing (approximately 2.06 - 5.596/sqrt(n)).}
#'   \item{est}{Vector of estimated generalized variance values for each predictor
#'     (lower values indicate stronger association).}
#'   \item{results}{Data frame with columns:
#'     \itemize{
#'       \item Pred.: Index of first predictor in comparison
#'       \item Pred: Index of second predictor in comparison
#'       \item test.stat: Test statistic for comparing the two predictors
#'       \item Decision: "reject" if predictors differ significantly, "fail to reject" otherwise
#'     }}
#'
#' @details
#' This function evaluates and compares predictor importance:
#' \enumerate{
#'   \item Standardizes all predictors to have comparable scales
#'   \item For each predictor, computes the generalized variance (via robust
#'     covariance) of the response with that predictor
#'   \item Uses bootstrap to estimate standard errors for pairwise comparisons
#'   \item Tests whether predictors differ significantly in their association strength
#' }
#'
#' The generalized variance provides a robust measure of association strength.
#' Lower values indicate stronger association (less residual variability).
#'
#' For visualization when p=2, creates a scatterplot with different symbols
#' for each predictor plotted against Y (both on standardized scale).
#'
#' @note Outlier removal with `xout=TRUE` removes cases from predictor space
#'   before analysis using the function specified by `outfun`.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regvarp}}, \code{\link{regIVstr}}, \code{\link{regIVcom}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare strength of association for 3 predictors
#' set.seed(123)
#' x <- matrix(rnorm(100*3), ncol=3)
#' y <- 2*x[,1] + 0.5*x[,2] + rnorm(100)  # x[,1] is strongest
#' regpord(x, y)
#' }
regpord<-function(x,y,nboot=100,alpha=.05,SEED=TRUE,xout=FALSE,cov.fun=cov.mba,pr=TRUE,
plotit=TRUE,xlab="Standardized Predictors",ylab="Y",est=mean,scat=var,...){
#
# Compare strength of association of two predictors via
# some robust covariance matrix, with predictors standardized.
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
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
n=nrow(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regpord.sub,x,y,cov.fun)
ptot=(p^2-p)/2
# bvec is a p by nboot matrix.
est=regvarp(x,y,est=est,scat=scat)
regci<-matrix(0,ptot,4)
dimnames(regci)<-list(NULL,c("Pred.","Pred","test.stat","Decision"))
ic=0
crit05=2.06-5.596/sqrt(n)
if(pr){
print("est is the estimated generalized variance")
}
if(p==2){
if(plotit){
z=standm(x,locfun=lloc,est=mean,scat=var)
z1=cbind(z[,1],y)
z2=cbind(z[,2],y)
plot(rbind(z1,z2),type="n",xlab=xlab,ylab=ylab)
points(z1,pch="*")
points(z2,pch="+")
}}
for(j in 1:p){
for(k in 1:p){
if(j<k){
sqse<-mean((bvec[j,]-est[j]-bvec[k,]+est[k])^2)*nboot/(nboot-1)
test=(est[j]-est[k])/sqrt(sqse)
ic=ic+1
regci[ic,1]<-j
regci[ic,2]<-k
regci[ic,3]<-test
regci[ic,4]<-0
if(abs(test)>=crit05)regci[ic,4]<-1
}}}
regci=data.frame(regci)
flag=(regci[,4]==0)
regci[flag,4]="fail to reject"
regci[!flag,4]="reject"
list(crit.value=crit05,est=est,results=regci)
}


# ============================================================================
# smean
# ============================================================================

#' Multivariate Skipped Measure of Location
#'
#' Computes a multivariate skipped mean - the mean of the data after removing
#' outliers detected using projection-based or other multivariate outlier
#' detection methods.
#'
#' @param m An n-by-p numeric matrix containing the multivariate data.
#' @param cop Integer (1-6) specifying the center to use for outlier detection:
#'   \itemize{
#'     \item 1 = Donoho-Gasko median
#'     \item 2 = MCD (Minimum Covariance Determinant)
#'     \item 3 = Marginal medians (default)
#'     \item 4 = MVE (Minimum Volume Ellipsoid)
#'     \item 5 = TBS (Rocke's S-estimator)
#'     \item 6 = MBA (Median Ball Algorithm)
#'   }
#' @inheritParams common-params
#' @param op Integer (1-3) specifying outlier detection method:
#'   \itemize{
#'     \item 1 = Projection method (default)
#'     \item 2 = MGV (Minimum Generalized Variance) method
#'     \item 3 = Use method specified by `outfun`
#'   }
#' @param outfun Outlier detection function to use when `op = 3`
#'   (default: `outogk`).
#' @param cov.fun Covariance function for MGV method when `op = 2`
#'   (default: `rmba`).
#' @param STAND Logical. If `TRUE`, standardizes data before projection
#'   (default: `TRUE`).
#'
#' @return A numeric vector of length p containing the multivariate skipped
#'   mean (mean of the data after outliers are removed).
#'
#' @details
#' The skipped mean provides a robust multivariate location estimator by:
#' \enumerate{
#'   \item Detecting outliers using one of several multivariate methods
#'   \item Computing the usual mean of the remaining (non-outlier) points
#' }
#'
#' For projection-based detection (op=1), each point is projected onto the line
#' connecting it to the center. A point is flagged as an outlier if it is an
#' outlier in any projection using a boxplot rule (MM=FALSE) or MAD-based rule
#' (MM=TRUE).
#'
#' If n < 14, the function automatically switches to MGV method (op=2).
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{smeancr}}, \code{\link{smean2v2}}, \code{\link{outpro}},
#'   \code{\link{outmgv}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compute skipped mean for bivariate data
#' set.seed(123)
#' m <- matrix(rnorm(100), ncol = 2)
#' # Add some outliers
#' m <- rbind(m, matrix(c(5, 5, -5, -5), ncol = 2, byrow = TRUE))
#' smean(m)
#' }
smean<-function(m,cop=3,MM=FALSE,op=1,outfun=outogk,cov.fun=rmba,MC=FALSE,STAND=TRUE,...){
#
# m is an n by p matrix
#
# Compute a multivariate skipped measure of location
#
# op=1:
# Eliminate outliers using a projection method
# If in addition, MC=T, a multicore processor is used
# assuming your computer has multiple cores and the package
# multicore has been installed.
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)
#
# For each point
# consider the line between it and the center,
# project all points onto this line, and
# check for outliers using
#
# MM=FALSE, a boxplot rule.
# MM=TRUE, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# op=2 use mgv (function outmgv) method to eliminate outliers
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# op=3 use outlier method indicated by outfun
#
# Eliminate any outliers and compute means
#  using remaining data.
#
m<-elimna(m)
m=as.matrix(m)
if(nrow(m)<14)op=2
if(op==1){
if(!MC)temp<-outpro(m,plotit=FALSE,cop=cop,MM=MM,STAND=STAND)$keep
if(MC)temp<-outproMC(m,plotit=FALSE,cop=cop,MM=MM,STAND=STAND)$keep
}
if(op==2)temp<-outmgv(m,plotit=FALSE,cov.fun=cov.fun)$keep
if(op==3)temp<-outfun(m,plotit=FALSE,...)$keep
val<-apply(m[temp,],2,mean)
val
}


# ============================================================================
# smeancrv2
# ============================================================================

#' Test Hypothesis About Multivariate Skipped Mean (Version 2)
#'
#' Tests the hypothesis that the multivariate skipped mean equals a specified
#' null value using bootstrap methods. This is an improved version of `smeancr()`
#' with enhanced distance computation and optional multicore support.
#'
#' @param m An n-by-p numeric matrix containing the multivariate data.
#' @param nullv Numeric vector of length p specifying the null hypothesis values
#'   (default: vector of zeros).
#' @param cop Integer (1-4) specifying the center for projection-based outlier detection:
#'   \itemize{
#'     \item 1 = Donoho-Gasko median
#'     \item 2 = MCD (Minimum Covariance Determinant)
#'     \item 3 = Marginal medians (default)
#'     \item 4 = MVE (Minimum Volume Ellipsoid)
#'   }
#' @inheritParams common-params
#' @param nboot Number of bootstrap samples (default: 500).
#' @param plotit Logical. If `TRUE` and p=2, creates a plot showing bootstrap
#'   distribution and confidence region (default: `TRUE`).
#' @param xlab,ylab Axis labels for bivariate plot.
#' @param STAND Logical. If `TRUE`, standardizes data before projection
#'   (default: `TRUE`).
#'
#' @return A list with component:
#'   \item{p.value}{Bootstrap p-value for testing H0: skipped mean = nullv.}
#'
#' @details
#' This function is an improved version of `smeancr()` that:
#' \itemize{
#'   \item Uses `smean()` directly to compute skipped means (more efficient)
#'   \item Computes projection distances centered at the estimated skipped mean
#'     rather than at the origin
#'   \item Supports multicore processing via `MC` parameter
#' }
#'
#' The test procedure:
#' \enumerate{
#'   \item Compute the skipped mean estimate for the full sample
#'   \item For each bootstrap sample, compute the skipped mean
#'   \item Calculate projection distances of bootstrap estimates from the null value,
#'     centered at the sample estimate
#'   \item p-value is the proportion of bootstrap distances less than the distance
#'     of the null value from the estimate
#' }
#'
#' For bivariate data (p=2) with `plotit=TRUE`, creates a scatterplot showing:
#' \itemize{
#'   \item Bootstrap estimates of the skipped mean
#'   \item The estimated skipped mean (marked with +)
#'   \item Confidence region (convex hull of central bootstrap estimates)
#' }
#'
#' @note This version does not adjust the critical level based on sample size
#'   (unlike `smeancr()`), instead using a fixed 0.05 level.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{smeancr}}, \code{\link{smean}}, \code{\link{smean2v2}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test if bivariate skipped mean equals (0,0)
#' set.seed(123)
#' m <- matrix(rnorm(100), ncol=2)
#' m <- rbind(m, matrix(c(5,5,-5,-5), ncol=2))  # Add outliers
#' smeancrv2(m)
#'
#' # Use multicore processing
#' smeancrv2(m, MC=TRUE)
#' }
smeancrv2<-function(m,nullv=rep(0,ncol(m)),cop=3,MM=FALSE,SEED=TRUE,
nboot=500,plotit=TRUE,MC=FALSE,xlab="VAR 1",ylab="VAR 2",STAND=TRUE){
#
# m is an n by p matrix
#
# Test hypothesis that multivariate skipped estimators
# are all equal to the null value, which defaults to zero.
# The level of the test is .05.
#
# Eliminate outliers using a projection method
# That is, determine center of data using:
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
# cop=4 MVE
#
# For each point
# consider the line between it and the center
# project all points onto this line, and
# check for outliers using
#
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# Eliminate any outliers and compute means
#  using remaining data.
#
if(SEED)set.seed(2)
m<-elimna(m)
n<-nrow(m)
est=smean(m,MC=MC,cop=cop,STAND=STAND)
crit.level<-.05
if(n<=120)crit.level<-.045
if(n<=80)crit.level<-.04
if(n<=60)crit.level<-.035
if(n<=40)crit.level<-.03
if(n<=30)crit.level<-.025
if(n<=20)crit.level<-.02
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
val<-matrix(NA,ncol=ncol(m),nrow=nboot)
for(j in 1: nboot){
mm<-m[data[j,],]
val[j,]<-smean(mm,MC=MC,cop=cop,STAND=STAND)
}
if(!MC)temp<-pdis(rbind(val,nullv),center=est)
if(MC)temp<-pdisMC(rbind(val,nullv),center=est)
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
if(ncol(m)==2 && plotit){
plot(val[,1],val[,2],xlab=xlab,ylab=ylab)
temp3<-est
points(temp3[1],temp3[2],pch="+")
ic<-round((1-crit.level)*nboot)
if(!MC)temp<-pdis(val,center=est)
if(MC)temp<-pdisMC(val,center=est)
temp.dis<-order(temp)
xx<-val[temp.dis[1:ic],]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(p.value=sig.level)
}


# ============================================================================
# regpca
# ============================================================================
#' Principal Component Analysis for Regression
#'
#' Performs standard principal component analysis (PCA) on a data matrix, with
#' options for correlation or covariance-based analysis, loading display, and
#' scree plot visualization.
#'
#' @param x A numeric matrix or data frame containing the variables for PCA.
#'   Rows with missing values are automatically removed.
#' @param cor Logical; if `TRUE` (default), use the correlation matrix for PCA.
#'   If `FALSE`, use the covariance matrix.
#' @param loadings Logical; if `TRUE` (default), display loadings in the output.
#' @param SCORES Logical; if `TRUE`, return only the principal component scores.
#'   If `FALSE` (default), return the full PCA summary.
#' @param pval Number of principal components to retain (default: all components).
#' @param scree Logical; if `TRUE` (default), display a scree plot showing variance
#'   explained by each component.
#' @param xlab Label for x-axis of scree plot (default: "Principal Component").
#' @param ylab Label for y-axis of scree plot (default: "Proportion of Variance").
#'
#' @return If `SCORES = TRUE`, returns a matrix of principal component scores.
#'   Otherwise, returns the output from `summary.princomp()`, which includes
#'   standard deviations, loadings, and variance proportions.
#'
#' @details
#' This function is a wrapper for `princomp()` that provides convenient options
#' for visualization and output formatting. The scree plot displays:
#' \itemize{
#'   \item Individual variance proportions (solid line, asterisks)
#'   \item Cumulative variance proportions (dashed line, dots)
#' }
#'
#' Missing values are removed using `elimna()` before analysis.
#'
#' @seealso \code{\link[stats]{princomp}}, \code{\link{regIVstr}}
#'
#' @export
#' @examples
#' \dontrun{
#' # PCA on iris data
#' data(iris)
#' x <- iris[, 1:4]
#' regpca(x)  # Full PCA with scree plot
#'
#' # Get only scores
#' scores <- regpca(x, SCORES = TRUE)
#'
#' # Use covariance instead of correlation
#' regpca(x, cor = FALSE)
#' }
regpca<-function(x,cor=TRUE,loadings=TRUE,
SCORES=FALSE,pval=ncol(x),scree=TRUE,xlab="Principal Component",ylab="Proportion of Variance"){
#
# regular PCA, calls princomp
#
x<-elimna(x) # removes any rows having missing values
temp<-princomp(x,cor=cor,scores=TRUE)
if(!SCORES)temp<-summary(temp,loadings=loadings)
if(SCORES){
return(temp$scores)
}
if(scree){
z=temp$sdev
pv=z^2
cs=pv/sum(pv)
cm=cumsum(cs)
plot(rep(c(1:ncol(x)),2),c(cs,cm),type="n",xlab=xlab,ylab=ylab)
points(c(1:ncol(x)),cs,pch="*")
lines(c(1:ncol(x)),cs,lty=1)
points(c(1:ncol(x)),cm,pch=".")
lines(c(1:ncol(x)),cm,lty=2)
}
temp
}


# ============================================================================
# gamindt
# ============================================================================

#' Test for Association Using Generalized Additive Model
#'
#' Tests the hypothesis of no association between predictors and response using
#' a generalized additive model (GAM) with a bootstrap-based permutation test.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values (maximum 4 columns).
#' @param y A numeric vector containing the response variable.
#' @param nboot Number of bootstrap samples for the permutation test
#'   (default: 500).
#'
#' @return A numeric p-value for testing H0: no association between x and y.
#'
#' @details
#' This function tests for association by comparing the observed strength of
#' association (from GAM fit) to its null distribution obtained via bootstrap
#' permutation:
#' \enumerate{
#'   \item Fits a GAM to the data and computes strength of association
#'   \item For each bootstrap sample: randomly pairs x and y observations and
#'     computes strength of association under this permutation
#'   \item Computes p-value as proportion of permuted samples with weaker
#'     association than observed
#' }
#'
#' The strength of association is measured as the square root of explanatory
#' power (variance of fitted values / variance of y).
#'
#' If `xout = TRUE`, outliers in the predictor space are removed before
#' testing using the function specified by `outfun`.
#'
#' @note This test requires the `mgcv` package for GAM fitting.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{gamplot}}, \code{\link{gamplotv2}}, \code{\link[mgcv]{gam}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test for association
#' set.seed(123)
#' x <- seq(-3, 3, length=100)
#' y <- sin(x) + rnorm(100, sd=0.2)
#' gamindt(x, y)  # Should reject H0 (nonlinear association exists)
#'
#' # Test with no association
#' y2 <- rnorm(100)
#' gamindt(x, y2)  # Should not reject H0
#' }
gamindt<-function(x,y,nboot=500,xout=FALSE,outfun=out){
#
# Test the hypothesis of no association based on the fit obtained with
# a generalized additive model
#
m<-elimna(cbind(x,y))
x<-as.matrix(x)
p<-ncol(x)
pp<-p+1
x<-m[,1:p]
y<-m[,pp]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,pp]
}
n=length(y)
data1<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
val=NA
x=as.matrix(x)
for(i in 1:nboot){
val[i]=gamplotv2(x[data1[i,],],y[data2[i,]],plotit=FALSE,pr=FALSE)$Strength.Assoc
}
val=sort(val)
est=gamplotv2(x,y,plotit=FALSE)$Strength.Assoc
p.value=mean((est<val))
p.value
}


# ============================================================================
# gamplotv2
# ============================================================================

#' Plot GAM Regression Surface with Strength of Association
#'
#' Enhanced version of `gamplot()` that fits and plots a GAM regression surface
#' while also computing measures of strength of association and explanatory power,
#' with optional adjustment for independence.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values (maximum 4 columns).
#' @param y A numeric vector containing the response variable.
#' @param sop Logical. If `TRUE`, uses spline smoothing (s()) in the GAM;
#'   if `FALSE`, uses linear terms (default: `FALSE`).
#' @param pyhat Logical. If `TRUE`, returns fitted values; if `FALSE`,
#'   returns strength/explanatory power measures (default: `FALSE`).
#' @param varfun Function to compute variance for explanatory power
#'   (default: `pbvar` for percentage bend variance).
#' @param cor.fun Function to compute correlation as fallback when
#'   explanatory power ≥ 1 (default: `pbcor`).
#' @param ADJ Logical. If `TRUE`, computes adjusted measures accounting for
#'   expected association under independence (default: `FALSE`).
#' @param SCALE Logical. For 3D plots, controls z-axis scaling (default: `FALSE`).
#' @param theta,phi,expand,ticktype Parameters for 3D perspective plots.
#' @param xlab,ylab,zlab Axis labels for plots.
#'
#' @return A list with components:
#'   \item{Strength.Assoc}{Estimated strength of association (square root of
#'     explanatory power).}
#'   \item{Explanatory.power}{Ratio of variance of fitted values to variance
#'     of y (capped at correlation^2 if > 1).}
#'   \item{Strength.Adj}{Adjusted strength of association (if `ADJ = TRUE`),
#'     accounting for expected value under independence.}
#'   \item{Explanatory.Adj}{Adjusted explanatory power (if `ADJ = TRUE`).}
#'
#' @details
#' This function extends `gamplot()` by computing:
#' - **Explanatory power**: var(fitted values) / var(y), measuring proportion
#'   of variance explained
#' - **Strength of association**: sqrt(explanatory power), analogous to
#'   correlation
#'
#' When `ADJ = TRUE`, the function estimates the expected strength under
#' independence via bootstrap permutation and adjusts the observed values:
#' \deqn{Adjusted = (Observed - E[Independent]) / (1 - E[Independent])}
#'
#' The variance and correlation functions (`varfun`, `cor.fun`) can be changed
#' to use different robust estimators.
#'
#' @note Requires `mgcv` package for GAM fitting and `akima` for 3D surface
#'   interpolation.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{gamplot}}, \code{\link{gamindt}}, \code{\link{gamplotINT}}
#'
#' @export
#' @examples
#' \dontrun{
#' # GAM with nonlinear relationship
#' set.seed(123)
#' x <- seq(-3, 3, length=100)
#' y <- sin(x) + rnorm(100, sd=0.2)
#' result <- gamplotv2(x, y, sop=TRUE)
#' result$Strength.Assoc
#'
#' # With adjusted measures
#' result_adj <- gamplotv2(x, y, sop=TRUE, ADJ=TRUE, nboot=100)
#' result_adj$Strength.Adj
#' }
gamplotv2<-function(x,y,sop=FALSE,pyhat=FALSE,eout=FALSE,xout=FALSE,outfun=out,plotit=TRUE,
varfun=pbvar,xlab="X",ylab="",zlab="",theta=50,phi=25,expand=.5,SCALE=FALSE,
cor.fun=pbcor,ADJ=FALSE,nboot=20,pr=TRUE,SEED=TRUE,ticktype="simple"){
#
# Plot regression surface using generalized additive model
#
# sop=F, use lowess
# sop=T, use splines
#
if(ADJ){
if(SEED)set.seed(2)
}
if(pr){
if(!ADJ){
print("To get adjusted estimates of strength of association, use ADJ=T")
print("The strength of association is estimated under independence")
print(" and then rescaled")
}}
library(mgcv)
x<-as.matrix(x)
np<-ncol(x)
np1<-np+1
if(ncol(x)>4)stop("x should have at most four columns of data")
m<-elimna(cbind(x,y))
if(xout && eout)stop("Can't have xout=eout=T")
if(eout){
flag<-outfun(m)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
}
x<-m[,1:np]
x=as.matrix(x)
y<-m[,np1]
if(!sop){
if(ncol(x)==1)fitr<-fitted(gam(y~x[,1]))
if(ncol(x)==2)fitr<-fitted(gam(y~x[,1]+x[,2]))
if(ncol(x)==3)fitr<-fitted(gam(y~x[,1]+x[,2]+x[,3]))
if(ncol(x)==4)fitr<-fitted(gam(y~x[,1]+x[,2]+x[,3]+x[,4]))
}
if(sop){
if(ncol(x)==1)fitr<-fitted(gam(y~s(x[,1])))
if(ncol(x)==2)fitr<-fitted(gam(y~s(x[,1])+s(x[,2])))
if(ncol(x)==3)fitr<-fitted(gam(y~s(x[,1])+s(x[,2])+s(x[,3])))
if(ncol(x)==4)fitr<-fitted(gam(y~s(x[,1])+s(x[,2])+s(x[,3])+s(x[,4])))
}
last<-fitr
if(plotit){
if(ncol(x)==1){
plot(x,fitr,xlab=xlab,ylab=ylab)
}
if(ncol(x)==2){
iout<-c(1:length(fitr))
nm1<-length(fitr)-1
for(i in 1:nm1){
ip1<-i+1
for(k in ip1:length(fitr))if(sum(x[i,]==x[k,])==2)iout[k]<-0
}
fitr<-fitr[iout>=1] # Eliminate duplicate points in the x-y plane
#                 This is necessary when doing three dimensional plots
#                 with the S-PLUS function interp
mkeep<-x[iout>=1,]
fitr<-interp(mkeep[,1],mkeep[,2],fitr)
persp(fitr,theta=theta,phi=phi,expand=expand,xlab="x1",ylab="x2",zlab="",
scale=scale,ticktype=ticktype)
}
}
top=varfun(last)
ep=top/varfun(y)
if(ep>=1)ep=cor.fun(last,y)$cor^2
eta=sqrt(ep)
st.adj=NULL
e.adj=NULL
if(ADJ){
x=as.matrix(x)
val=NA
n=length(y)
data1<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(i in 1:nboot){
temp=gamplotv2.sub(x[data1[i,],],y[data2[i,]],plotit=FALSE)
val[i]=temp$Explanatory.power
}
vindt=median(val)
v2indt=median(sqrt(val))
st.adj=(sqrt(ep)-max(c(0,v2indt)))/(1-max(c(0,v2indt)))
e.adj=(ep-max(c(0,vindt)))/(1-max(c(0,vindt)))
st.adj=max(c(0,st.adj))
e.adj=max(c(0,e.adj))
}
eta=as.matrix(eta)
ep=as.matrix(ep)
dimnames(eta)=NULL
dimnames(ep)=NULL
eta=eta[1]
ep=ep[1]
list(Strength.Assoc=eta,Explanatory.power=ep,
Strength.Adj=st.adj,Explanatory.Adj=e.adj)
}


# ============================================================================
# logreg.P.ci
# ============================================================================

#' Confidence Intervals for Logistic Regression Probabilities
#'
#' Computes pointwise confidence intervals for the probability P(Y=1|X) at each
#' observed x value in a logistic regression model, with optional outlier removal
#' and visualization.
#'
#' @inheritParams common-params
#' @param x A numeric vector or matrix containing predictor values.
#' @param y A binary numeric vector (will be converted to 0/1).
#' @param alpha Significance level for confidence intervals (default: 0.05 for
#'   95% CIs).
#' @param xlab,ylab Axis labels for the plot.
#'
#' @return A list with components:
#'   \item{Strength.Assoc}{Strength of association measured as sd(fitted
#'     probabilities) / sd(y).}
#'   \item{output}{Matrix with columns: X (predictor values), est.p (estimated
#'     probability), ci.low (lower CI bound), ci.up (upper CI bound), sorted by X.}
#'
#' @details
#' This function fits a standard logistic regression model using `glm()` and
#' computes pointwise confidence intervals for P(Y=1|X) based on the standard
#' errors from the linear predictor scale:
#' \enumerate{
#'   \item Fits logistic model: logit(P(Y=1|X)) = b0 + b1*X + ...
#'   \item Obtains standard errors for fitted values on logit scale
#'   \item Constructs normal-theory CIs on logit scale
#'   \item Transforms back to probability scale via inverse logit
#' }
#'
#' If `xout = TRUE`, outliers in the predictor space are removed before fitting
#' using the function specified by `outfun` (default: projection-based outlier
#' detection via `outpro`).
#'
#' If `plotit = TRUE`, creates a plot showing:
#' - Estimated probability curve (solid line)
#' - Confidence bands (dashed lines)
#'
#' @note The output matrix column is named "ci,up" (with comma) due to legacy
#'   naming - this is maintained for backwards compatibility.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{logreg}}, \code{\link{logreg.pred}}, \code{\link{wlogreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Logistic regression with confidence bands
#' set.seed(123)
#' x <- seq(-3, 3, length=100)
#' p <- exp(1 + 2*x) / (1 + exp(1 + 2*x))
#' y <- rbinom(100, 1, p)
#' result <- logreg.P.ci(x, y, plotit = TRUE)
#' head(result$output)
#' }
logreg.P.ci<-function(x,y,alpha=.05,plotit=TRUE,
xlab='X',ylab='P(Y=1|X)',xout=FALSE,outfun=outpro,...){
#
# Assuming the logistic regression model provides an adequate fit,
# compute a confidence interval for P(Y=1|X) for each value stored in x.
#
xx<-elimna(cbind(x,y))
x<-xx[,1]
y<-xx[,2]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1]
y<-m[,2]
}
if(length(unique(y))>2)stop('y should be binary')
# Next convert y to 0 and 1s
n=length(y)
yy=rep(0,n)
y=as.vector(y)
flag=y==max(y)
yy[flag]=1
y=yy

xord=order(x)
x=x[xord]
y=y[xord]
mod1 = glm(y ~ x, family=binomial(link='logit'))
v=predict(mod1,se.fit=TRUE)
top=v$fit+qnorm(1-alpha/2)*v$se.fit
bot=v$fit-qnorm(1-alpha/2)*v$se.fit
p=exp(v$fit)/(1+exp(v$fit))
top=exp(top)/(1+exp(top))
bot=exp(bot)/(1+exp(bot))
est=cbind(x,p,bot,top)
dimnames(est)=list(NULL,c('X','est.p','ci.low','ci,up'))
if(plotit){
plot(c(x,x,x),c(top,bot,p),ylim=c(0,1),type='n',xlab=xlab,ylab=ylab)
lines(x,p)
lines(x,bot,lty=2)
lines(x,top,lty=2)
}
list(Strength.Assoc=sd(p)/sd(y),output=est)
}


# ============================================================================
# Mreglde.sub
# ============================================================================

#' Objective Function for Mreglde Optimization
#'
#' Internal helper function used by \code{\link{Mreglde}} to compute the sum of
#' L1 distances (sum of Euclidean norms of residual vectors) for multivariate
#' regression parameter optimization.
#'
#' @param x Vector containing n, ncx, ncy, followed by flattened X and Y matrices.
#' @param B Flattened parameter vector (intercepts and slope matrix).
#'
#' @return Sum of Euclidean distances of residuals (objective to minimize).
#'
#' @keywords internal
#' @export
Mreglde.sub<-function(x,B){
n=x[1]
ncx=x[2]
ncy=x[3]
nxx=n*ncx
nyy=n*ncy
ncx1=ncx+1
B=matrix(B,nrow=ncx1,ncol=ncy)
iu=nxx+3
xm=matrix(x[4:iu],ncol=ncx)
il=iu+1
ym=matrix(x[il:length(x)],ncol=ncy)
ainit=B[1:ncy]
il=ncy+1
Binit=matrix(B[il:length(B)],nrow=ncx,ncol=ncy)
yhat=matrix(0,nrow=n,ncol=ncy)
for(i in 1:n){
z=as.matrix(xm[i,])
yhat[i,]=t(Binit)%*%z
}
yhat=t(t(yhat)+ainit)
res=ym-yhat
res=sum(sqrt(apply(res^2,1,sum)))
res
}


# ============================================================================
# mgvreg
# ============================================================================

#' Regression with MGV Outlier Detection
#'
#' Performs regression on data after removing outliers identified by the MGV
#' (minimum generalized variance) method. Computes regression coefficients,
#' residuals, and strength of association.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values.
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use after outlier removal (default: `tsreg`).
#' @param cov.fun Covariance function for MGV outlier detection (default: `rmba`).
#'   Determines the robust center used by `outmgv()`.
#' @param se Logical. Currently not used (kept for compatibility).
#' @param varfun Variance function for computing explanatory power (default: `pbvar`).
#' @param corfun Correlation function for computing explanatory power when variance
#'   ratio >= 1 (default: `pbcor`).
#' @param SEED Logical. If `TRUE`, sets seed for reproducibility of `outmgv()`
#'   results (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{coef}{Vector of regression coefficients (intercept and slopes).}
#'   \item{residuals}{Vector of residuals from the regression on cleaned data.}
#'   \item{Strength.Assoc}{Square root of explanatory power (similar to multiple R).}
#'   \item{Explanatory.Power}{Proportion of variance in Y explained by predictors,
#'     analogous to R-squared but using robust variance estimates.}
#'
#' @details
#' This function combines robust outlier detection with regression:
#'
#' \enumerate{
#'   \item Uses `outmgv()` to identify outliers in the (X,Y) space using the
#'     minimum generalized variance method
#'   \item Removes outliers from the data
#'   \item Fits regression using `regfun` (default: Theil-Sen) on cleaned data
#'   \item Computes explanatory power as var(fitted)/var(Y) using robust estimators
#' }
#'
#' The MGV method detects outliers based on robust covariance estimation (specified
#' by `cov.fun`). The default `rmba` uses the median ball algorithm.
#'
#' Explanatory power is computed robustly:
#' \itemize{
#'   \item If var(fitted)/var(Y) < 1: use this ratio directly
#'   \item If ratio >= 1: use squared correlation between fitted and observed
#' }
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{outmgv}}, \code{\link{tsreg}}, \code{\link{rmba}},
#'   \code{\link{mgvfreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Regression with MGV outlier detection
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2 + 3*x[,1] - x[,2] + rnorm(50)
#' # Add outliers
#' y[1:3] <- y[1:3] + 10
#' result <- mgvreg(x, y)
#' result$coef
#' result$Explanatory.Power
#' }
mgvreg<-function(x,y,regfun=tsreg,cov.fun=rmba,se=TRUE,varfun=pbvar,corfun=pbcor,
SEED=TRUE){
#
# Do regression on points not labled outliers
# by the MGV method.
# (This function replaces an older version of mgvreg as of 11/6/06)
#
# SEED=T so that results from outmgv are always duplicated using the same data
#
# In contrast to the old version,
#  when calling outmgv, center of data is determined via
#  the measure of location corresponding to cov.fun, which defaults
#  to the median ball algorithm (MBA)
#
x=as.matrix(x)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
ivec<-outmgv(m,plotit=FALSE,cov.fun=cov.fun,SEED=SEED)$keep
np1<-ncol(x)+1
y=m[ivec,np1]
x=m[ivec,1:ncol(x)]
coef<-regfun(x,y)$coef
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
# gamplotINT
# ============================================================================

#' Plot GAM Regression Surface with Interaction Term
#'
#' Plots a 3D regression surface for two predictors using generalized additive
#' models (GAM) with an interaction term included. This extends `gamplot()` by
#' allowing the predictors to interact.
#'
#' @inheritParams common-params
#' @param x A numeric matrix with exactly 2 columns containing the two predictor
#'   variables.
#' @param y A numeric vector containing the response variable.
#' @param pyhat Logical. If `TRUE`, returns predicted Y values; if `FALSE`,
#'   returns `NULL` (default: `FALSE`).
#' @param plotit Logical. If `TRUE`, creates 3D surface plot (default: `TRUE`).
#' @param theta Azimuthal viewing angle for 3D plot in degrees (default: 50).
#' @param phi Colatitude viewing angle for 3D plot in degrees (default: 25).
#' @param expand Expansion factor for z-axis (default: 0.5).
#' @param SCALE Logical. If `TRUE`, standardizes predictors before fitting
#'   (default: `FALSE`).
#' @param zscale Logical. If `TRUE`, standardizes the response variable
#'   (default: `TRUE`).
#' @param eout Logical. If `TRUE`, removes outliers from joint (X,Y) space
#'   before fitting (default: `FALSE`).
#' @param ticktype Type of tick marks for 3D plot (default: "simple").
#' @param xlab,ylab,zlab Axis labels for the 3D plot.
#'
#' @return If `pyhat = TRUE`, returns a vector of predicted Y values;
#'   otherwise returns `NULL`. Creates a 3D plot as a side effect.
#'
#' @details
#' This function extends `gamplot()` by including an interaction term between
#' the two predictors:
#'
#' \deqn{E(Y|X_1,X_2) = s_1(X_1) + s_2(X_2) + s_{12}(X_1 \times X_2)}
#'
#' where s() represents smooth functions estimated via GAM.
#'
#' The interaction term allows the effect of one predictor to vary with the
#' level of the other predictor, providing more flexibility than the additive
#' model in `gamplot()`.
#'
#' Outlier handling options:
#' \itemize{
#'   \item `xout = TRUE`: Remove outliers from predictor space only
#'   \item `eout = TRUE`: Remove outliers from joint (X,Y) space
#'   \item Cannot set both `xout` and `eout` to TRUE
#' }
#'
#' The 3D plot is created using `persp()` with viewing angles controlled by
#' `theta` (rotation) and `phi` (elevation).
#'
#' Requires the mgcv package for GAM fitting.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{gamplot}}, \code{\link{gamplotv2}}, \code{\link{gamindt}}
#'
#' @export
#' @examples
#' \dontrun{
#' # GAM surface with interaction
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2*x[,1] + 3*x[,2] + x[,1]*x[,2] + rnorm(50)
#' gamplotINT(x, y)
#' }
gamplotINT<-function(x,y,pyhat=FALSE,plotit=TRUE,theta=50,phi=25,expand=.5,xout=FALSE,
SCALE=FALSE,zscale=TRUE,eout=FALSE,outfun=out,ticktype="simple",xlab = "X", ylab = "Y", zlab = "",...){
#
# Plot regression surface, assuming two predictors in
# n by 2 matrix x using gam (generalized additive model)
# Same as gamplot, only a product term is included.
#
if(eout && xout)stop("Not allowed to have eout=xout=T")
x<-as.matrix(x)
if(ncol(x)!=2)stop("x must be an n by 2 matrix")
library(mgcv)
np=ncol(x)
np1=np+1
m<-elimna(cbind(x,y))
x<-m[,1:np]
x<-as.matrix(x)
y<-m[,np1]
if(xout){
flag<-outfun(x,...)$keep
m<-m[flag,]
}
if(eout){
flag<-outfun(m,...)$keep
m<-m[flag,]
}
x1<-m[,1]
x2<-m[,2]
y<-m[,3]
xrem<-m[,1:2]
n<-nrow(x)
fitr<-fitted(gam(y~s(x1)+s(x2)+s(x1,x2)))
allfit<-fitr
if(plotit){
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
}
m<-"Done"
if(pyhat)m<-allfit
m
}


# ============================================================================
# logreg.pred
# ============================================================================

#' Predict Probabilities from Logistic Regression
#'
#' Predicts the probability of success (Y=1) at specified predictor values using
#' a fitted logistic regression model, with options for robust estimation and
#' ridge regression.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing the training predictor values.
#' @param y A binary numeric vector (training response).
#' @param pts A numeric matrix or vector of predictor values at which to predict
#'   probabilities (default: `x`, predicting at training points).
#' @param ROB Logical. If `TRUE`, uses robust weighted logistic regression via
#'   `wlogreg()`; if `FALSE`, uses standard logistic regression (default: `FALSE`).
#' @param ridge Logical. If `TRUE`, uses ridge logistic regression; if `FALSE`,
#'   uses standard or robust logistic regression (default: `FALSE`).
#'
#' @return A numeric vector of predicted probabilities P(Y=1|X) at each point in
#'   `pts`, with values between 0 and 1.
#'
#' @details
#' This function provides a unified interface for predicting probabilities from
#' different types of logistic regression:
#' \itemize{
#'   \item **Standard** (`ROB=FALSE`, `ridge=FALSE`): Uses `logreg()` with
#'     optional outlier detection
#'   \item **Robust** (`ROB=TRUE`, `ridge=FALSE`): Uses `wlogreg()` with
#'     iteratively reweighted least squares
#'   \item **Ridge** (`ridge=TRUE`): Uses `logistic.ridge()` with L2
#'     regularization
#' }
#'
#' Predictions are computed using the inverse logit transformation:
#' \deqn{P(Y=1|X) = \frac{exp(\beta_0 + \beta_1 X_1 + ...)}{1 + exp(\beta_0 + \beta_1 X_1 + ...)}}
#'
#' If `xout = TRUE`, outliers in the training predictors are removed before
#' model fitting using the function specified by `outfun`.
#'
#' @note For single predictor models (p=2 coefficients), `pts` can be a vector.
#'   For multiple predictors, `pts` must be a matrix with the same number of
#'   columns as `x`.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{logreg}}, \code{\link{logreg.P.ci}}, \code{\link{wlogreg}},
#'   \code{\link{logistic.ridge}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Predict probabilities at new points
#' set.seed(123)
#' x <- rnorm(100)
#' p <- exp(1 + 2*x) / (1 + exp(1 + 2*x))
#' y <- rbinom(100, 1, p)
#'
#' # Predict at new x values
#' new_x <- seq(-2, 2, by=0.5)
#' pred_prob <- logreg.pred(x, y, pts=new_x)
#' plot(new_x, pred_prob, type='l')
#' }
logreg.pred<-function(x,y,pts=x,xout=FALSE,outfun=outpro,ROB=FALSE,ridge=FALSE){
#
# logistic regression: estimate the probability of success for points in pts
#  Default is to use pts=x
#
if(!ridge){
if(!ROB)est=logreg(x,y,xout=xout,outfun=outfun)[,1]
else
est=wlogreg(x,y,)$coef
}
if(ridge)est=logistic.ridge(x,y,xout=xout,outfun=outfun,ROB=ROB)$ridge.est
p=length(est)
if(p==2){z=exp(est[1]+est[2]*pts)
pr=z/(1+z)
}
if(p>2){
pr=NA
pts=as.matrix(pts)
if(ncol(pts)==1)pts=t(pts)
n=nrow(pts)
if(!is.matrix(pts))stop('pts should be a matrix')
if(ncol(pts)!=ncol(x))stop('pts should have the same number of col. as x')
for(i in 1:n){
z=exp(est[1]+sum(est[2:p]*pts[i,]))
pr[i]=z/(1+z)
}
}
pr
}


# ============================================================================
# logreg
# ============================================================================

#' Logistic Regression with Outlier Detection
#'
#' Performs logistic regression for binary outcomes, with optional outlier
#' removal from the predictor space and visualization of the fitted model.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values. Can have
#'   multiple columns for multiple predictors.
#' @param y A binary numeric vector (coded as 0/1) or any binary variable that
#'   will be converted to 0/1.
#' @param POLY Logical. If `TRUE` and p=1, assumes columns 2, 3, etc. of x
#'   contain transformations of column 1 (e.g., x^2, x^3) for polynomial
#'   regression (default: `FALSE`).
#' @param duplicate How to handle duplicate points in 2D plots. Options: "error"
#'   (default), "strip", "mean", etc. Passed to `interp()`.
#' @param scale Logical. If `TRUE`, uses scaling in 3D plots (default: `TRUE`).
#' @param expand,theta,phi,ticktype Parameters for 3D perspective plots when p=2.
#' @param xlab,ylab,zlab Axis labels for plots.
#'
#' @return A matrix with coefficient estimates, standard errors, z-statistics,
#'   p-values, and adjusted p-values (using Hochberg's method for multiple
#'   slopes). The first row is the intercept, subsequent rows are slopes.
#'
#' @details
#' This function fits a logistic regression model using `glm()` with optional
#' robust outlier detection in the predictor space before fitting.
#'
#' If `xout = TRUE`, outliers are removed from x before fitting using the
#' function specified by `outfun` (default is projection-based outlier
#' detection via `outpro`).
#'
#' For visualization:
#' - With 1 predictor: plots fitted probabilities vs. x
#' - With 2 predictors: creates 3D perspective plot of fitted surface
#' - Use `POLY=TRUE` for polynomial models with transformed predictors
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{logreg.P.ci}}, \code{\link{logreg.pred}},
#'   \code{\link{wlogreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Simple logistic regression
#' set.seed(123)
#' x <- rnorm(100)
#' p <- exp(1 + 2*x) / (1 + exp(1 + 2*x))
#' y <- rbinom(100, 1, p)
#' logreg(x, y, plotit = TRUE)
#' }
logreg<-function(x,y,xout=FALSE,outfun=outpro,plotit=FALSE,POLY=FALSE,
xlab='X',ylab='Y',zlab='',scale=TRUE ,expand=.5,theta=50,phi=25,
duplicate='error',ticktype='simple',...){
#
# Perform  logistic regression.
# The predictors are assumed to be stored in the n by p matrix x.
# The y values should be 1 or 0.
#
# xout=TRUE will remove outliers from among the x values and then fit
# the regression line.
#  Default:
# One predictor, a mad-median rule is used.
# With more than one, projection method is used.
#
# outfun=out will use MVE method
#
#  plotit=TRUE will plot regression line
#  POLY=T,  will plot regression line assuming predictor
#  is in  col 1 of x and other columns are x (in col 1) raised to some power
#   or some other function of x
#
y=chbin2num(y)
x<-as.matrix(x)
p=ncol(x)
xy=elimna(cbind(x,y))
n=nrow(xy)
x=xy[,1:ncol(x)]
y=xy[,ncol(xy)]
x<-as.matrix(x)

yy=rep(1,n)
vals=sort(unique(y))
if(length(vals)!=2)stop('y should be binary')
flag=y==vals[2]
yy[!flag]=0
y=yy

if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
}
x<-as.matrix(x)
if(p==1 || POLY){
xord=order(x[,1])
x=x[xord,]
y=y[xord]
}
fitit=glm(formula=y~x,family=binomial)
init<-summary(fitit)
if(plotit){
vals=fitted.values(fitit)
if(p==1){
plot(x,y,xlab=xlab,ylab=ylab)
lines(x,vals)
}
if(p==2){
if(!scale)print('With dependence, suggest using scale=TRUE')
fitr=vals
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
}
}
init$coef
p1=p+1
p.adjusted.slopes=c(init$coef[1,1],p.adjust(init$coef[2:p1,4],method='hoch'))
p.adjusted.slopes[1]=NA
a=cbind(init$coef,p.adjusted.slopes)
a
}


# ============================================================================
# mlrregCI
# ============================================================================

#' Bootstrap Confidence Intervals for Multivariate Regression Coefficients
#'
#' Computes percentile bootstrap p-values and confidence intervals for all
#' regression coefficients in a multivariate linear regression using the robust
#' Rousseeuw et al. (2004) estimator.
#'
#' @param x A numeric matrix containing predictor values (independent variables).
#' @param y A numeric matrix containing response values (dependent variables).
#'   Must be multivariate (matrix, not vector).
#' @param nboot Number of bootstrap samples (default: 300).
#' @inheritParams common-params
#' @param op.dis Logical. Reserved for future use (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{estimates}{Matrix of point estimates for all coefficients. Same
#'     structure as output from `mlrreg()$coef`.}
#'   \item{p.values}{Matrix of two-sided p-values for testing H0: coefficient = 0.
#'     Each p-value is computed as 2*min(P(boot > 0), P(boot < 0)).}
#'
#' @details
#' This function performs bootstrap inference for multivariate regression:
#' \enumerate{
#'   \item Fits robust multivariate regression using `mlrreg()` to get point estimates
#'   \item Generates `nboot` bootstrap samples (resampling rows with replacement)
#'   \item Fits `mlrreg()` to each bootstrap sample
#'   \item Computes p-values for each coefficient based on bootstrap distribution
#' }
#'
#' The p-value for each coefficient tests H0: β = 0 using the percentile method:
#' \deqn{p = 2 \times \min(P(\hat{\beta}^* > 0), P(\hat{\beta}^* < 0))}
#'
#' If `MC = TRUE`, uses parallel processing via `mclapply()` for faster bootstrap
#' computation on multicore systems.
#'
#' The coefficient matrix structure matches `mlrreg()`: first row is intercepts
#' (one per response variable), subsequent rows are slopes for each predictor.
#'
#' @references
#' Rousseeuw, P.J., Van Aelst, S., Van Driessen, K., & Agulló, J. (2004).
#' Robust multivariate regression. Technometrics, 46, 293-305.
#'
#' @seealso \code{\link{mlrreg}}, \code{\link{mlrregWtest}}, \code{\link{mulgreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Bootstrap inference for multivariate regression
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- matrix(NA, nrow=50, ncol=2)
#' y[,1] <- 1 + 2*x[,1] - x[,2] + rnorm(50)
#' y[,2] <- -1 + x[,1] + 3*x[,2] + rnorm(50)
#' result <- mlrregCI(x, y, nboot=100)
#' result$estimates
#' result$p.values
#' }
mlrregCI<-function(x,y,nboot=300,MC=FALSE,SEED=TRUE,op.dis=TRUE){
#
#  Based on Rousseeuw et al.
#  multivariate regression estimator
#  compute p-value for each of the parameters using a percentile
#  bootstrap method.
#
if(SEED)set.seed(2)
est=mlrreg(x,y)$coef
pval=est
n=nrow(x)
JK=(ncol(x)+1)*ncol(y)
vals=matrix(0,nrow=nboot,ncol=JK)
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
if(!MC)for(ib in 1:nboot){
vals[ib,]=mlrreg(x[data[ib,],],y[data[ib,],])$coef
}
if(MC){
data=listm(t(data))
vals=mclapply(data,mlrreg.est,x,y,mc.preschedule=TRUE)
vals=t(matl(vals))
}
pv=NULL
for(j in 1:JK){
pv[j]=mean(vals[,j]>0)+.5*mean(vals[,j]==0)
pv[j]=2*min(c(pv[j],1-pv[j]))
}
ic=0
il=1
iu=ncol(x)+1
for(iy in 1:ncol(y)){
pval[,iy]=pv[il:iu]
il=il+ncol(x)+1
iu=iu+ncol(x)+1
}
list(estimates=est,p.values=pval)
}


# ============================================================================
# mlrreg.est
# ============================================================================

#' Bootstrap Helper for mlrregCI
#'
#' Internal helper function used by \code{\link{mlrregCI}} to extract regression
#' coefficients from bootstrap samples for computing confidence intervals.
#'
#' @param data Bootstrap sample indices.
#' @param x Matrix of predictors.
#' @param y Matrix of responses.
#'
#' @return Vector of regression coefficients from \code{mlrreg()} on bootstrap sample.
#'
#' @keywords internal
#' @export
mlrreg.est<-function(data,x,y){
xv=x[data,]
yv=y[data,]
vals=as.vector(mlrreg(xv,yv)$coef)
vals
}


# ============================================================================
# mlrregWtest
# ============================================================================

#' Test All Slopes Equal Zero in Multivariate Regression
#'
#' Tests the global null hypothesis that all slope coefficients equal zero in
#' a multivariate linear regression using a wild bootstrap method with the robust
#' Rousseeuw et al. (2004) estimator.
#'
#' @param x A numeric matrix containing predictor values (independent variables).
#' @param y A numeric matrix containing response values (dependent variables).
#'   Must be multivariate (matrix, not vector).
#' @param nboot Number of bootstrap samples (default: 300).
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{p.value}{P-value for testing H0: all slope coefficients = 0 (intercepts
#'     not tested).}
#'
#' @details
#' This function tests the global hypothesis that X has no predictive power for Y:
#' \deqn{H_0: \beta_1 = \beta_2 = ... = \beta_p = 0}
#' where β₁, β₂, ..., βₚ are the slope coefficients (intercepts not tested).
#'
#' The test uses a wild bootstrap approach:
#' \enumerate{
#'   \item Fits `mlrreg()` to get estimated slopes (excludes intercepts)
#'   \item For each bootstrap sample: resamples Y values (rows) with replacement
#'     while keeping X fixed, fits `mlrreg()`, and extracts slopes
#'   \item Computes Mahalanobis-like distance from slopes to origin (0 vector)
#'   \item P-value = proportion of bootstrap distances ≥ observed distance
#' }
#'
#' This is a multivariate extension of testing whether a regression line has
#' zero slope. The test accounts for the correlation structure among multiple
#' response variables.
#'
#' If `MC = TRUE`, uses parallel processing via `mclapply()` for faster
#' computation on multicore systems.
#'
#' @note This is a one-sided test (large distances indicate rejection). The test
#'   focuses only on slopes; intercepts are not tested.
#'
#' @references
#' Rousseeuw, P.J., Van Aelst, S., Van Driessen, K., & Agulló, J. (2004).
#' Robust multivariate regression. Technometrics, 46, 293-305.
#'
#' @seealso \code{\link{mlrreg}}, \code{\link{mlrregCI}}, \code{\link{mulgreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test if predictors have any effect
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- matrix(rnorm(100), ncol=2)  # No relationship
#' mlrregWtest(x, y, nboot=100)  # Should not reject
#'
#' # With actual relationship
#' y[,1] <- 1 + 2*x[,1] - x[,2] + rnorm(50, sd=0.5)
#' y[,2] <- -1 + x[,1] + 3*x[,2] + rnorm(50, sd=0.5)
#' mlrregWtest(x, y, nboot=100)  # Should reject
#' }
mlrregWtest<-function(x,y,nboot=300,MC=FALSE,SEED=TRUE){
#
#  Test hypothesis that all slopes=0  based on Rousseeuw et al.
#  multivariate regression estimator
#
#  Strategy: a variation of the wild bootstrap method, percentile version.
#
if(SEED)set.seed(2)
estit=mlrreg.subest(y,x)  #YES, y before x
n=nrow(x)
JK=ncol(x)*ncol(y)
vals=matrix(0,nrow=nboot,ncol=JK)
data=list()
for(i in 1:nboot){
bsam=sample(n,replace=TRUE)
data[[i]]=y[bsam,]
}
if(!MC){
vals=lapply(data,mlrreg.subest,x)
}
if(MC){
vals=mclapply(data,mlrreg.subest,x,mc.preschedule=TRUE)
}
vals=t(matl(vals))
nullv=rep(0,JK)
vals=rbind(vals,estit)
cen=rep(0,ncol(vals))
if(MC)dv=pdisMC(vals,center=cen)
if(!MC)dv=pdis(vals,center=cen)
bplus=nboot+1
pv=1-sum(dv[bplus]>=dv[1:nboot])/nboot
list(p.value=pv)
}


# ============================================================================
# mlrreg.subest
# ============================================================================

#' Bootstrap Helper for mlrregWtest
#'
#' Internal helper function used by \code{\link{mlrregWtest}} to extract slope
#' coefficients (excluding intercepts) from bootstrap samples for hypothesis testing.
#'
#' @param data Bootstrap sample (resampled Y values).
#' @param x Matrix of predictors (fixed across bootstrap samples).
#'
#' @return Vector of slope coefficients (intercepts excluded) from \code{mlrreg()}.
#'
#' @keywords internal
#' @export
mlrreg.subest<-function(data,x){
vals=as.vector(mlrreg(x,data)$coef[-1,])
vals
}


# ============================================================================
# regmediate
# ============================================================================

#' Test for Mediation Effect in Regression
#'
#' Tests for mediation by comparing the direct effect (predictor → outcome)
#' with the conditional effect (predictor → outcome | mediator) using bootstrap
#' confidence intervals. Computes CI for the difference in slopes.
#'
#' @inheritParams common-params
#' @param x An n-by-2 numeric matrix where column 1 is the predictor and
#'   column 2 is the mediator variable.
#' @param y A numeric vector containing the outcome variable.
#' @param regfun Regression function to use (default: `tsreg` for Theil-Sen).
#'
#' @return A list with components:
#'   \item{conf.interval}{Bootstrap confidence interval for b11 - b13 (the
#'     mediation effect).}
#'   \item{p.value}{P-value for testing H0: b11 - b13 = 0.}
#'
#' @details
#' Mediation analysis examines whether the effect of a predictor X on outcome Y
#' operates through a mediator M. This function compares two regression models:
#' \itemize{
#'   \item Direct effect: Y = b01 + b11*X + e1
#'   \item Conditional effect: Y = b03 + b13*X + b23*M + e3
#' }
#'
#' If mediation is present, b13 should be smaller than b11 in absolute value,
#' as some of X's effect on Y is now explained by M.
#'
#' The function computes a bootstrap confidence interval for (b11 - b13), which
#' estimates the indirect effect through the mediator. Significant mediation
#' occurs when this CI excludes zero.
#'
#' **Requirements**:
#' - x must have exactly 2 columns: [predictor, mediator]
#' - Use `MC = TRUE` for parallel processing with multiple cores
#' - Use `xout = TRUE` to remove leverage points before analysis
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regmed2}}, \code{\link{tsreg}}, \code{\link{regci}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Mediation analysis
#' set.seed(123)
#' n <- 100
#' x1 <- rnorm(n)          # Predictor
#' m <- 0.5*x1 + rnorm(n)  # Mediator
#' y <- 0.3*x1 + 0.4*m + rnorm(n)  # Outcome (partial mediation)
#' x <- cbind(x1, m)
#' regmediate(x, y)
#' }
regmediate<-function(x,y,regfun=tsreg,nboot=400,alpha=.05,xout=FALSE,outfun=out,MC=FALSE,SEED=TRUE,...){
#
#   In a mediation analysis, two of the linear equations that play a role are
#   y=b_{01} + b_{11}x + e_1
#   y=b_{03} + b_{13}x + b_{23} x_m + e_3
#   where x_m is the mediator variable.
#   An additional assumption is
#   x_m=b_{02} + b_{12}x + \epsilon_2.
#   Goal: Compute a confidence interval for b_{11}-b_{13}
#
#   The default regression method is the Theil-Sen estimator.
#
#   The predictor values are assumed to be in the n-by-2 matrix x, with the
#   mediator variable in column 2.
#   MC=T. A multicore processor will be used.
#   xout=T will remove leverage points using the function indicated by the argument out.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
if(p!=2)stop("Argument x should have two columns")
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
if(MC){
bvec1<-mclapply(data,regbootMC,as.matrix(x[,1]),y,regfun,mc.preschedule=TRUE)
bvec2<-mclapply(data,regbootMC,x,y,regfun,mc.preschedule=TRUE)
}
if(!MC){
bvec1<-lapply(data,regboot,as.matrix(x[,1]),y,regfun)
bvec2<-lapply(data,regboot,x,y,regfun)
}
bvec1=matl(bvec1)
bvec2=matl(bvec2)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
dif=bvec1[2,]-bvec2[2,]
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
sig.level<-NA
temp<-mean(dif<0)
sig.level<-2*(min(temp,1-temp))
bsort<-sort(dif)
regci<-bsort[ilow]
regci[2]<-bsort[ihi]
list(conf.interval=regci,p.value=sig.level)
}


# ============================================================================
# regmed2
# ============================================================================

#' Test Mediation Pathways in Regression
#'
#' Tests the two key regression pathways in mediation analysis: (1) the effect
#' of the predictor on the mediator, and (2) the effect of the mediator on the
#' outcome controlling for the predictor. Uses bootstrap confidence intervals.
#'
#' @inheritParams common-params
#' @param x An n-by-2 numeric matrix where column 1 is the predictor and
#'   column 2 is the mediator variable.
#' @param y A numeric vector containing the outcome variable.
#' @param regfun Regression function to use (default: `tsreg` for Theil-Sen).
#' @param pr Logical. If `TRUE`, prints explanatory text about the output
#'   (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{res1}{Results for predicting mediator from predictor: 95% CI, estimate,
#'     SE, and p-value for slope b12 in mediator ~ predictor.}
#'   \item{res2}{Results for the mediator's effect: 95% CI, estimate, SE, and
#'     p-value for slope b23 in outcome ~ predictor + mediator.}
#'
#' @details
#' Mediation analysis examines three regression equations:
#' \enumerate{
#'   \item **Predictor → Mediator**: \eqn{M = b_{02} + b_{12} X + \epsilon_2}
#'   \item **Direct effect**: \eqn{Y = b_{01} + b_{11} X + e_1}
#'   \item **Conditional effect**: \eqn{Y = b_{03} + b_{13} X + b_{23} M + e_3}
#' }
#'
#' This function tests two critical hypotheses:
#' - **H0: b12 = 0** (predictor doesn't affect mediator)
#' - **H0: b23 = 0** (mediator doesn't affect outcome after controlling for predictor)
#'
#' For mediation to be supported, both hypotheses should typically be rejected,
#' indicating that:
#' 1. The predictor influences the mediator (res1)
#' 2. The mediator influences the outcome even after accounting for the
#'    direct effect of the predictor (res2)
#'
#' **Requirements**:
#' - x must have exactly 2 columns: [predictor, mediator]
#' - Use `MC = TRUE` for parallel processing
#' - Use `xout = TRUE` to remove leverage points
#'
#' @note Unlike `regmediate()` which tests the change in the predictor's
#'   coefficient, this function tests the individual pathway slopes.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regmediate}}, \code{\link{regci}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Mediation analysis
#' set.seed(123)
#' n <- 100
#' x1 <- rnorm(n)          # Predictor
#' m <- 0.7*x1 + rnorm(n)  # Mediator (strong effect from predictor)
#' y <- 0.2*x1 + 0.5*m + rnorm(n)  # Outcome (partial mediation)
#' x <- cbind(x1, m)
#' result <- regmed2(x, y, nboot=200)
#' result$res1  # Test X → M (should be significant)
#' result$res2  # Test M → Y | X (should be significant)
#' }
regmed2<-function(x,y,regfun=tsreg,nboot=400,alpha=.05,xout=FALSE,outfun=out,MC=FALSE,
SEED=TRUE,pr=TRUE,...){
#
#   In a mediation analysis, two of the linear equations that play a role are
#   y=b_{01} + b_{11}x + e_1
#   y=b_{03} + b_{13}x + b_{23} x_m + e_3
#   where x_m is the mediator variable.
#   An additional assumption is
#   x_m=b_{02} + b_{12}x + \epsilon_2.
#   Goal: Test hypotheses  b_{12}=0 and b_{23}=0
#
#   The default regression method is the Theil-Sen estimator.
#
#   The predictor values are assumed to be in the n-by-2 matrix x, with the
#   mediator variable in column 2.
#   MC=T. A multicore processor will be used.
#   xout=T will remove leverage points using the function indicated by the argument out.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
if(p!=2)stop("Argument x should have two columns")
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
if(MC){
temp1=regciMC(x[,1],x[,2],regfun=regfun,nboot=nboot,alpha=alpha,SEED=SEED,pr=FALSE)
temp2=regciMC(x,y,regfun=regfun,nboot=nboot,alpha=alpha,SEED=SEED,pr=FALSE)
}
if(!MC){
temp1=regci(x[,1],x[,2],regfun=regfun,nboot=nboot,alpha=alpha,SEED=SEED,pr=FALSE)
temp2=regci(x,y,regfun=regfun,nboot=nboot,alpha=alpha,SEED=SEED,pr=FALSE)
}
if(pr){
print("Output returned in res1 is for the slope of the regression line")
print("where the goal is to predict the mediator variable given the other")
print("predictor variable stored in column 1 of x.")
print("Output in res2 is for slope of the mediator when both predictors are used.")
}
res1=c(temp1$regci[2,],temp1$p.value[2])
z1=t(as.matrix(res1))
dimnames(z1)=list(NULL,c("ci.low","ci.up",'Estimate','S.E.',"p.value"))
res2=c(temp2$regci[3,],temp2$p.value[3])
z2=t(as.matrix(res2))
dimnames(z2)=list(NULL,c("ci.low","ci.up",'Estimate','S.E.',"p.value"))
list(res1=z1,res2=z2)
}


# ============================================================================
# COVreg
# ============================================================================

#' Regression Estimation via Robust Covariance Matrix
#'
#' Performs regression estimation using robust covariance and location estimators
#' instead of the classical maximum likelihood approach. This provides resistance
#' to outliers in both predictors and response.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values. Can have multiple columns
#'   for multiple predictors.
#' @param y A numeric vector containing the response variable.
#' @param cov.fun Covariance function for robust estimation (default: `MARest`).
#'   Must return a component `$cov` with the covariance matrix.
#' @param loc.fun Location function for robust center estimation (default: `MARest`).
#'   Must return a component `$center` with location estimates.
#'
#' @return A list with components:
#'   \item{coef}{Vector of regression coefficients. First element is intercept,
#'     remaining elements are slopes for each predictor.}
#'   \item{residuals}{Vector of residuals (y - fitted values).}
#'
#' @details
#' Classical regression estimation via maximum likelihood uses the sample
#' covariance matrix, which is sensitive to outliers. This function replaces
#' the classical covariance with a robust estimator:
#'
#' The regression coefficients are derived from the robust covariance matrix
#' AC of (X,Y):
#' \itemize{
#'   \item Slope = solve(AC_XX, AC_XY) where AC_XX is covariance of X and
#'     AC_XY is covariance between X and Y
#'   \item Intercept = center_Y - slopes %*% center_X
#' }
#'
#' Default uses `MARest()` for both covariance and location, but any robust
#' estimator can be substituted (e.g., `rmba`, `covmba`, `covmcd`).
#'
#' If `xout = TRUE`, outliers in the predictor space are removed before
#' estimation using the function specified by `outfun`.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{mlrreg}}, \code{\link{mulgreg}}, \code{\link{tsreg}},
#'   \code{\link{rmba}}, \code{\link{covmba}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Robust regression via covariance matrix
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 2)
#' y <- 2 + 3*x[,1] - x[,2] + rnorm(50)
#' result <- COVreg(x, y)
#' result$coef
#'
#' # With outlier removal
#' COVreg(x, y, xout = TRUE)
#' }
COVreg<-function(x,y,cov.fun=MARest,loc.fun=MARest,xout=FALSE,outfun=out,...){
#
# Regression estimation can be done via the usual maximum likelihood
# covariance matrix. This function uses the same approach
# using a robust covariance matrix instead.
#
# The predictors are assumed to be stored in the n-by-p matrix x.
#
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
x=xy[,1:p]
y=xy[,p1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
AC=cov.fun(cbind(x,y),...)$cov
ma<-AC[1:p,p1]
m<-AC[1:p,1:p]
slope<-solve(m,ma)
mvals<-loc.fun(cbind(x,y))$center
b0<-mvals[p1]-sum(slope%*%mvals[1:p])
res<-y-x%*%slope-b0
list(coef=c(b0,slope),residuals=res)
}


# ============================================================================
# mlrreg
# ============================================================================

#' Robust Multivariate Linear Regression
#'
#' Performs robust multivariate linear regression using the method of Rousseeuw
#' et al. (2004), which combines MCD-based outlier detection with least squares
#' or Theil-Sen regression on the cleaned data.
#'
#' @param x A numeric matrix containing predictor values (independent variables).
#' @param y A numeric matrix containing response values (dependent variables).
#'   Must be multivariate (matrix, not vector).
#' @param cov.fun Covariance function for robust estimation (default: `cov.mcd`).
#'   Must return components `$center` and `$cov`.
#' @param ols.op Logical. If `TRUE`, uses OLS on cleaned data; if `FALSE`, uses
#'   Theil-Sen regression via `mopreg()` (default: `TRUE`).
#' @param mcd.op Logical. If `TRUE`, uses MCD with `quantile.used`; if `FALSE`,
#'   uses `cov.fun` with its default settings (default: `TRUE`).
#' @param quantile.used Quantile for MCD estimation (default: floor(0.75*n),
#'   following Rousseeuw et al.).
#' @param RES Logical. If `TRUE`, returns residuals; if `FALSE`, returns `NULL`
#'   for residuals (default: `FALSE`).
#' @param ... Additional arguments passed to covariance functions.
#'
#' @return A list with components:
#'   \item{coef}{Matrix of regression coefficients. First row is intercepts,
#'     subsequent rows are slopes for each predictor.}
#'   \item{residuals}{Matrix of residuals if `RES = TRUE`, otherwise `NULL`.}
#'
#' @details
#' This function implements the robust multivariate regression method from
#' Rousseeuw et al. (2004):
#' \enumerate{
#'   \item Computes robust location and scatter estimates for (X,Y)
#'   \item Derives initial regression coefficients from robust covariance
#'   \item Computes Mahalanobis distances of residuals
#'   \item Removes outliers (points with distance > chi-square 99th percentile)
#'   \item Re-fits using OLS or Theil-Sen on cleaned data
#' }
#'
#' The `quantile.used` parameter controls the breakdown point. The default
#' (0.75*n) gives a 25% breakdown point as recommended by Rousseeuw et al.
#'
#' @references
#' Rousseeuw, P.J., Van Aelst, S., Van Driessen, K., & Agulló, J. (2004).
#' Robust multivariate regression. Technometrics, 46, 293-305.
#'
#' @seealso \code{\link{mulgreg}}, \code{\link{MULMreg}}, \code{\link{mlrregCI}},
#'   \code{\link{mlrregWtest}}, \code{\link{COVreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Multivariate regression with 2 predictors and 2 responses
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- matrix(NA, nrow=50, ncol=2)
#' y[,1] <- 1 + 2*x[,1] - x[,2] + rnorm(50)
#' y[,2] <- -1 + x[,1] + 3*x[,2] + rnorm(50)
#' result <- mlrreg(x, y)
#' result$coef
#' }
mlrreg<-function(x,y,cov.fun=cov.mcd,ols.op=TRUE,mcd.op=TRUE,
quantile.used=floor(.75*n),RES=FALSE,...){
#
# Do Multivariate regression, using by default the method
#  in Rousseeuw, Van Aelst, Van Driessen Agullo
# Technometrics, 46, 293-305
#
# Note, to use the method recommended by Rousseeuw et al., the argument
# quantile.used=.75*n is used when calling cov.mcd.
#
#  RES=T, the residuals will be returned.
#
# y is assumed to be  multivariate with data stored in a matrix.
#
# an initial fit is found using the measures of scatter and location
# corresponding to cof.fun and mcd.op. If
# mcd.op=T, cov.mcd is used with quanitle.used=.75n
# mcd.op=F, cov.fun is used and defaults to cov.mcd with the
# default value usded by R for the argument quanitle.used
# But any function that returns location and scatter in $center and $cov
# can be used.
#
#  if ols.op=T, OLS is applied after points are removed based on iniital fit
#  if ols.op=F, Theil-Sen is used by calling the function mopreg
#
#  Early version of this function considered estimating
#  explanatory power in terms of the generalized variance
#  of the predicted y values and the observed y values
#  epow.cov determines which robust covariance matrix will be used.
#  This idea has not been explored enough
#  Some choices are:
# cov (the usual generalized variance)
# skipcov
# tbscov
# covout
# covogk
# mgvcov
# mvecov
# mcdcov
#
if(!is.matrix(y))stop("y is not a matrix")
X<-cbind(x,y)
X<-elimna(X)
n<-nrow(X)
qy<-ncol(y)
qx<-ncol(x)
qxp1<-qx+1
tqyqx<-qy+qx
y<-X[,qxp1:tqyqx]
# compute initial estimate of slopes and intercept:
if(!mcd.op)locscat<-cov.fun(X,...)
if(mcd.op)locscat<-cov.mcd(X,quan=quantile.used)
sig<-locscat$cov
mu<-locscat$center
sigxx<-sig[1:qx,1:qx]
sigxy<-sig[1:qx,qxp1:tqyqx]
sigyy<-sig[qxp1:tqyqx,qxp1:tqyqx]
Bhat<-solve(sigxx)%*%sigxy
sige<-sigyy-t(Bhat)%*%sigxx%*%Bhat
sige.inv<-solve(sige)
Ahat<-t(mu[qxp1:tqyqx]-t(Bhat)%*%mu[1:qx])
resL<-matrix(nrow=nrow(X),ncol=qy)
for(i in 1:nrow(X))resL[i,]<-y[i,]-t(Bhat)%*%X[i,1:qx]
for(j in 1:qy)resL[,j]<-resL[,j]-Ahat[j]
drL<-NA
for(i in 1:nrow(X))drL[i]<-t(resL[i,])%*%sige.inv%*%resL[i,]
# In Rousseeuw notation, drL<- is d^2
w<-rep(0,nrow(X))
qdr<-qchisq(.99,qy)
iflag<-(drL<qdr)
w[iflag]<-1
term1<-0
vec<-c(1:nrow(X))
keep<-vec[iflag==1]
X<-X[keep,]
if(ols.op)output<-lsfit(X[,1:qx],X[,qxp1:tqyqx])
if(!ols.op)output<-mopreg(X[,1:qx],X[,qxp1:tqyqx],KEEP=TRUE)
yhat=X[,qxp1:tqyqx]-output$residuals
res=NULL
if(RES)res=output$residuals
#epow=(gvarg(yhat,epow.cov)/gvarg(X[,qxp1:tqyqx],epow.cov))
#list(coef=output$coefficients,residuals=res,E.power=epow,Strength.Assoc=sqrt(epow))
list(coef=output$coefficients,residuals=res)
}


# ============================================================================
# Mreglde
# ============================================================================

#' Multivariate Regression via Least Distance Estimator
#'
#' Performs multivariate regression where parameters are estimated using the
#' least distance estimator (LDE), a robust method that minimizes the sum of
#' Mahalanobis distances between observed and fitted values.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values (independent variables).
#' @param y A numeric matrix containing response values (dependent variables).
#'   Must be multivariate (matrix with multiple columns).
#' @param eout Logical. If `TRUE`, removes outliers from the entire (X,Y) space
#'   before estimation (default: `FALSE`).
#' @param epow.cov Covariance function for computing explanatory power
#'   (default: `mcdcov`). Not currently used in returned output.
#' @param RES Logical. If `TRUE`, returns residuals; if `FALSE`, returns `NULL`
#'   for residuals (default: `FALSE`).
#'
#' @return A list with components:
#'   \item{coef}{Matrix of regression coefficients. First row contains intercepts,
#'     subsequent rows contain slopes for each predictor. Columns correspond to
#'     response variables.}
#'   \item{residuals}{Matrix of residuals (Y - fitted Y) if `RES = TRUE`,
#'     otherwise `NULL`.}
#'
#' @details
#' The least distance estimator (LDE) minimizes:
#' \deqn{\sum_i d^2(y_i, x_i'\beta)}
#' where d() is the Mahalanobis distance based on robust covariance estimation.
#'
#' This provides robustness to outliers in multivariate regression by down-weighting
#' observations with large Mahalanobis distances. The method is particularly useful
#' when both predictors and responses may contain outliers.
#'
#' Predictors are standardized (centered and scaled) before optimization using
#' Nelder-Mead via `nelderv2()`. Initial values come from classical least squares.
#'
#' Options for outlier detection:
#' \itemize{
#'   \item `xout = TRUE`: Remove outliers from predictor space only
#'   \item `eout = TRUE`: Remove outliers from joint (X,Y) space
#'   \item Cannot set both `xout` and `eout` to TRUE
#' }
#'
#' @references
#' Jhun, M., & Choi, I. (2009). Bootstrapping least distance estimator in the
#' multivariate regression model. Computational Statistics & Data Analysis,
#' 53, 4221-4227.
#'
#' @seealso \code{\link{mlrreg}}, \code{\link{mulgreg}}, \code{\link{COVreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Multivariate regression with 2 predictors and 2 responses
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- matrix(NA, nrow=50, ncol=2)
#' y[,1] <- 1 + 2*x[,1] - x[,2] + rnorm(50)
#' y[,2] <- -1 + x[,1] + 3*x[,2] + rnorm(50)
#' result <- Mreglde(x, y, RES = TRUE)
#' result$coef  # Intercepts in row 1, slopes in rows 2-3
#' }
Mreglde<-function(x,y,xout=FALSE,eout=FALSE,outfun=outpro,epow.cov=mcdcov,RES=FALSE,...){
#
# Do multivariate regression where parameters are
# estimated via the least distance estimator.
#  See Jhun and Choi (2009). Comp Stat & Data Analysis, 53, 4221-4227
#
#  RES=T, the residuals will be returned.
#
if(eout){
flag=outfun(cbind(x,y),...)$keep
x=x[flag,]
y=y[flag,]
}
if(xout){
flag=outfun(x,...)$keep
x=x[flag,]
y=y[flag,]
}
npar=(ncol(x)+1)*ncol(y)
xy=elimna(cbind(x,y))
x=xy[,1:ncol(x)]
for(i in 1:ncol(x))x[,i]=(x[,i]-mean(x[,i]))/sqrt(var(x[,i]))
p1=ncol(x)+1
y=xy[,p1:ncol(xy)]
INIT=as.vector(lsfit(x,y)$coef)
xx=c(nrow(x),ncol(x),ncol(y),as.vector(x),as.vector(y))
Bs=nelderv2(xx,npar,Mreglde.sub,START=INIT)
Bs=matrix(Bs,ncol=ncol(y))
dimnames(Bs)<-list(c("INTER",rep("SLOPE",ncol(x))),rep("Y",ncol(Bs)))
yhat=matrix(NA,nrow=nrow(y),ncol=ncol(y))
for(i in 1:nrow(y)){
z=as.matrix(x[i,])
yhat[i,]=t(Bs[2:nrow(Bs),])%*%z
}
yhat=yhat+Bs[1,]
res=NULL
if(RES)res=y-yhat
#epow=gvarg(yhat,epow.cov)/gvarg(y,epow.cov)
list(coef=Bs,residuals=res)
}


# ============================================================================
# mlrreg.Stest
# ============================================================================

#' Test All Slopes Equal Zero in Multivariate Regression
#'
#' Tests the hypothesis that all regression slopes equal zero using Rousseeuw et al.'s
#' robust multivariate regression estimator with bootstrap standard errors and a
#' Hotelling-type test statistic.
#'
#' @param x A numeric matrix containing predictor values (independent variables).
#' @param y A numeric matrix containing response values (dependent variables).
#'   Must be multivariate (matrix with multiple columns).
#' @param nboot Number of bootstrap samples for estimating covariance of slopes
#'   (default: 100).
#' @param SEED Logical. If `TRUE`, sets random seed to 2 for reproducibility
#'   (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{test.stat}{The test statistic (F-distributed under H0).}
#'   \item{p.value}{The p-value for testing H0: all slopes = 0.}
#'   \item{est}{Vector of estimated slope coefficients (column-major order from
#'     coefficient matrix).}
#'
#' @details
#' This function tests the global null hypothesis H0: all slopes = 0 (only intercepts
#' are non-zero) for multivariate regression:
#'
#' \enumerate{
#'   \item Estimates slopes using `mlrreg()` (Rousseeuw et al. robust estimator)
#'   \item Bootstraps the slope estimates to get covariance matrix Sv
#'   \item Computes Hotelling-type statistic: T^2 = est' * Sv^-1 * est / k
#'     where k = number of slope parameters
#'   \item Converts to F-statistic with df1 = k-1, df2 = n-k
#' }
#'
#' The test is robust to outliers because it uses the MCD-based regression
#' estimator and bootstrap inference. With p predictors and q responses, there
#' are k = p*q slope parameters being tested.
#'
#' @references
#' Rousseeuw, P.J., Van Aelst, S., Van Driessen, K., & Agulló, J. (2004).
#' Robust multivariate regression. Technometrics, 46, 293-305.
#'
#' @seealso \code{\link{mlrreg}}, \code{\link{mlrregWtest}}, \code{\link{mlrregCI}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test if predictors have any effect on responses
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- matrix(rnorm(100), ncol=2)  # No relationship
#' mlrreg.Stest(x, y, nboot=200)  # Should not reject H0
#'
#' # Add relationship
#' y[,1] <- 2*x[,1] + rnorm(50)
#' mlrreg.Stest(x, y, nboot=200)  # Should reject H0
#' }
mlrreg.Stest<-function(x,y,nboot=100,SEED=TRUE){
#
#  Test hypothesis that all slopes=0  based on Rousseeuw et al.
#  multivariate regression estimator
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Hotelling type test.
#
if(SEED)set.seed(2)
est=as.vector(mlrreg(x,y)$coef[-1,])
n=nrow(x)
JK=ncol(x)*ncol(y)
vals=matrix(0,nrow=nboot,ncol=JK)
for(i in 1:nboot){
bsam=sample(n,replace=TRUE)
vals[i,]=as.vector(mlrreg(x[bsam,],y[bsam,])$coef[-1,])
}
Sv=cov(vals)
est=as.matrix(est)
k=1/JK
test <- k * crossprod(est, solve(Sv, est))[1, ]
v1=JK-1
v2=n-JK
pval=1-pf(test,v1,v2)
list(test.stat=test,p.value=pval,est=est)
}


# ============================================================================
# qhdplotsm
# ============================================================================

#' Plot Smoothed Quantile Regression Lines
#'
#' Plots one or more smoothed quantile regression lines using the Harrell-Davis
#' estimator via `rplotsm()`. Useful for visualizing conditional quantiles as
#' smooth functions of a predictor.
#'
#' @inheritParams common-params
#' @param x A numeric vector containing the single predictor variable.
#' @param y A numeric vector containing the response variable.
#' @param q Numeric vector of quantiles to plot (default: 0.5 for median).
#'   Can specify multiple quantiles, e.g., `c(0.25, 0.5, 0.75)`.
#' @param xlab,ylab Axis labels for the plot.
#' @param pc Plotting character for data points (default: ".").
#' @param nboot Number of bootstrap samples for smoothing (default: 40).
#' @param fr Span parameter for smoothing, between 0 and 1 (default: 1).
#'   Controls bandwidth of the smoother.
#'
#' @return Invisibly returns NULL. Creates a plot as a side effect.
#'
#' @details
#' This function visualizes conditional quantiles:
#' \enumerate{
#'   \item Plots the scatterplot of (x, y)
#'   \item For each quantile in `q`:
#'     \itemize{
#'       \item Computes quantile regression using Harrell-Davis estimator
#'       \item Smooths the fitted values using running-interval smoother
#'       \item Overlays the smooth curve on the plot
#'     }
#' }
#'
#' Multiple quantiles can be plotted simultaneously to show the conditional
#' distribution of Y given X. Common choices:
#' \itemize{
#'   \item `q = 0.5`: Median regression line
#'   \item `q = c(0.25, 0.5, 0.75)`: Quartile regression lines
#'   \item `q = c(0.1, 0.5, 0.9)`: Shows variability and skewness
#' }
#'
#' The `fr` parameter controls smoothness: smaller values give more local
#' smoothing, larger values give more global smoothing.
#'
#' Only one predictor is allowed. For multiple predictors, use `qhdsm()` or
#' `gamplot()` instead.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{qhdsm}}, \code{\link{rplotsm}}, \code{\link{hd}},
#'   \code{\link{qreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Plot median regression line
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' qhdplotsm(x, y, q = 0.5)
#'
#' # Plot quartile regression lines
#' qhdplotsm(x, y, q = c(0.25, 0.5, 0.75))
#' }
qhdplotsm<-function(x,y,q=.5,xlab="X",ylab="Y",pc=".",
xout=FALSE,outfun=out,nboot=40,fr=1,...){
#
# Plots smooths of quantile regression lines for one or more quantiles
# using rplotsm with Harrell--Davis estimator
#
# q indicates the quantiles to be used.
#
#  EXAMPLE:
#  qhdplotsm(x,y,q=c(.2,.5,.8)) will plot three smooths corresponding to
#  the .2, .5 and .8 quantile regression lines.
#
xy=elimna(cbind(x,y))
x=as.matrix(x)
if(ncol(x)!=1)stop("Only one predictor is allowed")
x=xy[,1]
y=xy[,2]
if(xout){
flag<-outfun(x,...)$keep
x<-x[flag]
y<-y[flag]
}
plot(x,y,xlab=xlab,ylab=ylab,pch=pc)
xord=order(x)
for(j in 1:length(q)){
yhat=rplotsm(x,y,fr=fr,pyhat=TRUE,est=hd,q=q[j],plotit=FALSE,nboot=nboot)$yhat
lines(x[xord],yhat[xord])
}
print("Done")
}


# ============================================================================
# longreg
# ============================================================================

#' Regression Analysis for Longitudinal Data
#'
#' Fits individual regression lines for each subject in longitudinal data,
#' then summarizes the distribution of regression coefficients across subjects
#' using a robust location estimator.
#'
#' @param x A data frame or matrix containing the longitudinal data in long format.
#' @param x.col Numeric. Column number(s) in `x` containing predictor variable(s),
#'   typically time points when measurements were taken.
#' @param y.col Numeric. Column number in `x` containing the outcome variable.
#' @param s.id Numeric. Column number in `x` containing the subject ID variable.
#' @param regfun Regression function to use for individual subject fits
#'   (default: `tsreg` for Theil-Sen regression).
#' @param est Location estimator for summarizing coefficients across subjects
#'   (default: `tmean` for 20% trimmed mean).
#'
#' @return A list with components:
#'   \item{est.S}{Matrix where each row contains regression coefficients
#'     (intercept and slopes) for one subject. Rows correspond to subjects.}
#'   \item{typical.est}{Vector of typical coefficient values across subjects,
#'     computed by applying `est` to each column of `est.S`.}
#'
#' @details
#' This function analyzes growth curves or change over time in longitudinal data:
#'
#' \enumerate{
#'   \item Reshapes data from long to wide format using `long2mat()` and
#'     `longcov2mat()`
#'   \item For each subject, fits regression: Y ~ X using `regfun`
#'   \item Stores coefficients in matrix (one row per subject)
#'   \item Summarizes typical coefficients using `est` (e.g., trimmed mean)
#' }
#'
#' The typical estimates represent:
#' \itemize{
#'   \item Element 1: Typical intercept across subjects
#'   \item Elements 2+: Typical slopes for each predictor
#' }
#'
#' Example use cases:
#' \itemize{
#'   \item Growth curves: Y = size, X = age
#'   \item Longitudinal treatment effects: Y = outcome, X = time
#'   \item Individual differences in change rates
#' }
#'
#' Data must be in long format with one row per observation. For example,
#' the `Orthodont` dataset from the nlme package is in the correct format.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{long2mat}}, \code{\link{longcov2mat}}, \code{\link{tsreg}},
#'   \code{\link{tmean}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Requires nlme package for Orthodont data
#' library(nlme)
#' # Fit growth curves: distance vs age for each subject
#' # Column 2 = age, column 1 = distance, column 4 = Subject ID
#' result <- longreg(Orthodont, x.col=2, y.col=1, s.id=4)
#' result$typical.est  # Typical intercept and slope
#' }
longreg<-function(x,x.col,y.col,s.id,regfun=tsreg,est=tmean){
#
# x is a data frame or matrix
#
# Longitudinal data.
# For each subject, fit a regression line
# using outcome data in col y.col and predictors, usually times
# when measures were taken, in columns indicated by x.col.
# s.id indicates column where subject's id is stored.
#
# Assuming data are stored as for example in the R variable
# Orthodont,
# which can be accessed via the command  library(nlme)
# For this data set, x.col=2 would indicated that the
# participants age at the time of being measured, is used
# to predict the outcome variable.
#
ymat=long2mat(x,s.id,y.col) # matrix, ith row contains outcome y
#                           for the ith subject.
#
xvals=longcov2mat(x,s.id,x.col)# list mode
n=nrow(ymat)
p=length(x.col)+1
outmat=matrix(NA,nrow=n,ncol=p)
for(i in 1:n)outmat[i,]=regfun(as.matrix(xvals[[i]]),ymat[i,])$coef
typval=apply(outmat,2,est)
list(est.S=outmat,typical.est=typval)
}


# ============================================================================
# regYband
# ============================================================================

#' Plot Confidence Band for Regression Predictions
#'
#' Creates a plot with confidence bands or intervals for predicted Y values across
#' the range of X, with options for pointwise or simultaneous coverage.
#'
#' @inheritParams common-params
#' @param x A numeric vector containing predictor values (single predictor only).
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use (default: `tsreg`).
#' @param npts Number of points to use for confidence intervals when `ADJ = FALSE`.
#'   If `NULL`, uses 20 evenly-spaced points between min(x) and max(x).
#' @param ADJ Logical. If `TRUE`, creates simultaneous confidence band with
#'   overall coverage 1-alpha; if `FALSE`, creates pointwise confidence
#'   intervals (default: `TRUE`).
#' @param SCAT Logical. If `TRUE`, plots the scatter plot of data points
#'   (default: `TRUE`).
#' @param nreps Number of replications for computing critical value (default: 1000).
#' @param pch Plotting character for scatter plot (default: ".").
#' @param xlab,ylab Axis labels.
#'
#' @return A matrix with columns: X value, predicted Y, lower CI, upper CI.
#'   The plot is created as a side effect.
#'
#' @details
#' This function creates a visual representation of regression uncertainty:
#' \itemize{
#'   \item With `ADJ = FALSE`: Plots pointwise confidence intervals at `npts`
#'     evenly-spaced x values. Each interval has individual coverage 1-alpha.
#'   \item With `ADJ = TRUE`: Plots a simultaneous confidence band at all
#'     unique x values with overall coverage 1-alpha. Wider than pointwise CIs.
#' }
#'
#' The simultaneous band accounts for multiple comparisons, ensuring that the
#' probability that ALL intervals simultaneously contain the true values is 1-alpha.
#'
#' This function is restricted to a single predictor. For multiple predictors,
#' use `regYci()` directly.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regYci}}, \code{\link{regci}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Confidence band for regression
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' # Simultaneous confidence band
#' regYband(x, y, ADJ = TRUE)
#' # Pointwise confidence intervals
#' regYband(x, y, ADJ = FALSE, npts = 30)
#' }
regYband<-function(x,y,regfun=tsreg,npts=NULL,nboot=100,xout=FALSE,outfun=outpro,SEED=TRUE,
alpha=.05,crit=NULL,xlab="X",ylab="Y",SCAT=TRUE,ADJ=TRUE,pr=TRUE,nreps=1000,
MC=FALSE,pch='.',...){
#
# Plot confidence band for the predicted Y value
# if ADJ=FALSE, plot confidence intervals for
#  npts points between min(x) and max(x)
#  if npts=NULL, then npts=20 is used.
# if ADJ=TRUE, plot confidence band for the predicted Y value for all x values such that
#  the simultaneous probability coverage is .95.
#
#  npts=NULL and ADJ=FALSE, npts will be set equal to 20. That is, computed confidence
#  intervals for 20 point covariate values even space between min(x) and max(x).
#
#
if(!ADJ){
if(is.null(npts))npts=20
if(pr)print('To adjust the confidence band so that the simultaneous probability coverage is .95, set ADJ=TRUE')
}
xy=elimna(cbind(x,y))
x<-as.matrix(x)
p=ncol(x)
if(p!=1)stop("This function assumes a single predictor only")
p1=p+1
vals=NA
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
if(SEED)set.seed(2)
if(!ADJ)pts=seq(min(x),max(x),length.out=npts)
if(ADJ)pts=sort(unique(x))
res=regYci(x,y,pts=pts,regfun=regfun,xout=FALSE,SEED=SEED,alpha=alpha,ADJ=ADJ,nreps=nreps,MC=MC,...)
plot(c(x,pts,pts),c(y,res[,2],res[,3]),xlab=xlab,ylab=ylab,type="n")
abline(regfun(x,y,...)$coef)
if(SCAT)points(x,y,pch=pch)
lines(pts,res[,3],lty=2)
lines(pts,res[,4],lty=2)
res
}


# ============================================================================
# regunstack
# ============================================================================

#' Unstack Data by Groups for Regression Analysis
#'
#' Sorts and reorganizes data into separate groups for one-way ANOVA-style
#' comparisons of regression slopes across multiple groups. Converts a single
#' data frame/matrix into lists of predictor matrices and response vectors
#' organized by group.
#'
#' @param x A matrix or data frame containing all variables including the grouping
#'   variable.
#' @param grp Column index (numeric) indicating which column of `x` contains the
#'   grouping variable.
#' @param xcols Vector of column indices indicating which columns contain the
#'   independent variables (predictors).
#' @param ycol Column index (numeric) indicating which column contains the
#'   dependent variable (response).
#'
#' @return A list with components:
#'   \item{num.grps}{Number of unique groups found.}
#'   \item{x}{List where x[[i]] is the predictor matrix for group i.}
#'   \item{y}{List where y[[i]] is the response vector for group i.}
#'
#' @details
#' This utility function reorganizes data for group-wise regression analysis:
#' \enumerate{
#'   \item Removes rows with missing values
#'   \item Identifies unique group levels in column `grp`
#'   \item Separates data into group-specific predictor matrices and response vectors
#' }
#'
#' Particularly useful when comparing regression relationships (slopes, intercepts)
#' across independent groups using functions like \code{reg1way()}.
#'
#' @note This is primarily a data preparation utility. Missing values are removed
#'   using \code{elimna()} before unstacking.
#'
#' @seealso \code{\link{reg1way}}, \code{\link{elimna}}
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' # Create grouped regression data
#' df <- data.frame(
#'   group = rep(1:3, each=20),
#'   x1 = rnorm(60),
#'   x2 = rnorm(60),
#'   y = rnorm(60)
#' )
#' # Unstack by group
#' result <- regunstack(df, grp=1, xcols=c(2,3), ycol=4)
#' result$num.grps  # Should be 3
#' result$x[[1]]    # Predictors for group 1
#' result$y[[1]]    # Response for group 1
#' }
regunstack<-function(x,grp,xcols,ycol){
#
#  x is assumed to be a matrix or a data frame
#
# sort data in x into group indicated by col grp of x,
# Designed for a one-way ANOVA where goal is to compare slopes
# corresponding to two or more groups.
#
# returns the independent variables in x having list mode
# x[[1]] would be a matrix for group 1, x[[2]] a matrix for group 2, etc
# y[[1]] is the dependent variable for group 1, etc.
#
# xcols indicates the columns of x containing independent variables
# ycol  indicates the column of x containing  dependent variables
#
x=elimna(x)
val=sort(unique(x[,grp]))
xs=list()
ys=list()
for(i in 1:length(val)){
flag=(x[,grp]==val[i])
xs[[i]]=x[flag,xcols]
ys[[i]]=x[flag,ycol]
}
list(num.grps=length(val),x=xs,y=ys)
}


# ============================================================================
# qhdsm
# ============================================================================

#' Quantile Regression Smoother Using Harrell-Davis Estimator
#'
#' Computes and plots quantile regression lines using the Harrell-Davis
#' estimator combined with a running interval smoother and LOESS. This provides
#' a nonparametric estimate of conditional quantiles.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values. Can be
#'   multivariate (but plotting works best with 1-2 predictors).
#' @param y A numeric vector containing the response variable.
#' @param qval A numeric vector specifying the quantiles to estimate (default: 0.5
#'   for median regression). Can specify multiple quantiles, e.g., c(0.25, 0.5, 0.75).
#' @param q Alias for `qval` (for backward compatibility).
#' @param pr Logical. If `TRUE`, prints progress messages.
#' @param fr Span parameter for running interval smoother (default: 0.8 for p=1,
#'   1 for p>1). Controls the fraction of data used in local smoothing.
#' @param LP Logical. If `TRUE`, applies LOESS smoothing after running interval
#'   smoother (default: `FALSE`).
#' @param pyhat Logical. If `TRUE`, returns predicted Y values; if `FALSE`,
#'   returns "DONE" message (default: `FALSE`).
#' @param nmin Minimum number of points required for smoothing (default: 0).
#' @param scale Logical. If `TRUE`, standardizes predictors when p > 1
#'   (default: `TRUE`).
#' @param pr.qhd Logical. If `TRUE`, prints warning about scaling with dependent
#'   predictors (default: `TRUE`).
#' @param theta,phi Viewing angles for 3D plot when p=2 (defaults: 50, 25).
#' @param ticktype Type of tick marks for 3D plot (default: "simple").
#' @param pch Plotting character for scatter plot (default: ".").
#' @param xlab,ylab,zlab Axis labels for plots.
#'
#' @return If `pyhat = TRUE`, returns a matrix of predicted quantile values.
#'   If `pyhat = FALSE`, returns "DONE" message. For p=1, plots the quantile
#'   regression line(s). For p=2, creates a 3D perspective plot.
#'
#' @details
#' This function estimates conditional quantiles using a combination of methods:
#' \enumerate{
#'   \item For each x value, computes the quantile of nearby y values using
#'         Harrell-Davis estimator
#'   \item Applies running interval smoother to reduce variability
#'   \item Optionally applies LOESS smoothing (if LP=TRUE)
#' }
#'
#' For a single predictor (p=1), multiple quantile lines can be plotted
#' simultaneously. For two predictors (p=2), only one quantile can be plotted
#' at a time using a 3D perspective plot.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{qhdsm2g}}, \code{\link{qhdplotsm}}, \code{\link{hd}},
#'   \code{\link{Qreg}}, \code{\link{rplot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Median regression line
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' qhdsm(x, y, qval = 0.5)
#'
#' # Multiple quantile lines (quartiles)
#' qhdsm(x, y, qval = c(0.25, 0.5, 0.75))
#' }
qhdsm<-function(x,y,qval=0.5,q=NULL,pr=FALSE,
xout=FALSE,outfun=outpro,plotit=TRUE,xlab='X',ylab='Y',zlab='Z',pyhat=FALSE,fr=NULL,LP=FALSE,theta=50,phi=25,ticktype='simple',nmin=0,scale=TRUE,pr.qhd=TRUE,pch='.',...){
#
# Compute the quantile regression line for one or more quantiles
# using combination of hd, running interval smoother and LOESS
# That is, determine the qth (qval) quantile of Y given X using the
#
#  plotit=TRUE will plot the lines. WIth p=1 predictor, multiple lines can be plotted.
#  Example: qhdsm(x,y,q=c(.25,.5,.75)) will plot the regression lines for
#   predicting quartiles.
  #
if(!is.null(q))qval=q
x<-as.matrix(x)
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
p<-np-1
if(p>1 & length(q)>1)print('Only first quantile specified can be plotted')
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
if(p==1){
if(is.null(fr))fr=.8
ord=order(x)
x=sort(x)
y=y[ord]
est=matrix(NA,ncol=3,nrow=length(qval))
dimnames(est)=list(NULL,c('q','Inter','Slope'))
#x<-as.matrix(x)
qest=matrix(NA,ncol=length(qval),nrow=length(y))
for(j in 1:length(qval)){
rmd=NA
for(i in 1:length(x))rmd[i]<-hd(y[near(x,x[i],fr)],q=qval[j])
if(LP)rmd=lplot(x,rmd,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
qest[,j]=rmd
}
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab,pch=pch)
for(j in 1:ncol(qest))lines(x,qest[,j])
}
if(!pyhat)qest='DONE'
}
if(p>1){
if(is.null(fr))fr=1
if(p==2){
if(pr.qhd){
if(!scale)print('scale=F is specified. If there is dependence, might want to use scale=TRUE')
}}
qest=rplot(x,y,est=hd,q=qval[1],fr=fr,plotit=plotit,pyhat=pyhat,theta=theta,
phi=phi,scale=scale,SEED=FALSE,varfun=pbvar,xlab=xlab,ylab=ylab,zlab=zlab,
ticktype=ticktype,nmin=nmin,pr=pr)
if(!pyhat)qest='DONE'
if(pyhat)qest=qest$yhat
}
qest
}


# ============================================================================
# qhdsm2g
# ============================================================================
#' Two-Group Quantile Regression Smoother
#'
#' Plots quantile regression smoothers for two independent groups using the
#' Harrell-Davis estimator. Useful for comparing regression curves at a specific
#' quantile (e.g., median regression) across groups.
#'
#' @param x1 Numeric vector of predictor values for group 1.
#' @param y1 Numeric vector of response values for group 1.
#' @param x2 Numeric vector of predictor values for group 2.
#' @param y2 Numeric vector of response values for group 2.
#' @param q Quantile to estimate (default: 0.5 for median). Must be between 0 and 1.
#' @param qval Deprecated; use `q` instead.
#' @param LP Logical; if `TRUE` (default), apply lowess smoothing to the quantile
#'   curves for visualization clarity.
#' @param fr Fraction of data used for local smoothing (default: 0.8). Controls
#'   the bandwidth - larger values produce smoother curves.
#' @param xlab Label for x-axis (default: 'X').
#' @param ylab Label for y-axis (default: 'Y').
#' @inheritParams common-params
#'
#' @return This function is called for its plotting side effect. It creates a
#'   scatter plot with overlaid quantile regression curves for both groups.
#'
#' @details
#' For each group, the function:
#' \enumerate{
#'   \item Computes the qth quantile of Y values near each X value using the
#'     Harrell-Davis estimator via `hd()` and `near()`
#'   \item Optionally smooths the resulting curve using lowess (`lplot()`)
#'   \item Plots the data and fitted curves
#' }
#'
#' Group 1 is plotted with circles and a solid line, group 2 with plus signs
#' and a dashed line.
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' fitting using the function specified by `outfun`.
#'
#' @note This function requires only one predictor variable. For multiple
#' predictors, see other regression methods.
#'
#' @seealso \code{\link{qhdsm}}, \code{\link{qhdsm.pred}}, \code{\link{hd}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare median regression curves for two groups
#' set.seed(123)
#' n <- 100
#' x1 <- runif(n, 0, 10)
#' x2 <- runif(n, 0, 10)
#' y1 <- 2 + 0.5*x1 + rnorm(n)
#' y2 <- 3 + 0.3*x2 + rnorm(n, sd=2)
#' qhdsm2g(x1, y1, x2, y2, q=0.5)
#'
#' # Compare 75th percentile regression
#' qhdsm2g(x1, y1, x2, y2, q=0.75)
#' }
qhdsm2g<-function(x1,y1,x2,y2,q=.5,qval=NULL,LP=TRUE,fr=.8,xlab='X',ylab='Y',xout=FALSE,outfun=outpro,...){
#
# Plot of quantile smoother for two groups using qhdsm
#
# fr controls amount of smoothing
# Missing values are automatically removed.
#
if(!is.null(qval))q=qval
m1<-elimna(cbind(x1,y1))
if(ncol(m1)>3)stop('One covariate only is allowed')
m2<-elimna(cbind(x2,y2))
x1<-m1[,1]
y1<-m1[,2]
x2<-m2[,1]
y2<-m2[,2]
if(xout){
flag<-outfun(m1[,1],plotit=FALSE,...)$keep
m1<-m1[flag,]
x1<-m1[,1]
y1<-m1[,2]
flag<-outfun(m2[,1],plotit=FALSE,...)$keep
m2<-m2[flag,]
x2<-m2[,1]
y2<-m2[,2]
}
flag=order(x1)
x1=x1[flag]
y1=y1[flag]
flag=order(x2)
x2=x2[flag]
y2=y2[flag]
rmd1=NA
rmd2=NA
for(i in 1:length(x1))rmd1[i]<-hd(y1[near(x1,x1[i],fr)],q=q)
for(i in 1:length(x2))rmd2[i]<-hd(y2[near(x2,x2[i],fr)],q=q)
if(LP){
rmd1=lplot(x1,rmd1,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
rmd2=lplot(x2,rmd2,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
plot(c(x1,x2),c(y1,y2),type='n',xlab=xlab,ylab=ylab)
points(x1,y1)
points(x2,y2,pch='+')
lines(x1,rmd1)
lines(x2,rmd2,lty=2)
}


# ============================================================================
# regse
# ============================================================================

#' Bootstrap Standard Errors for Regression Parameters
#'
#' Estimates standard errors and the covariance matrix for regression parameter
#' estimates using bootstrap resampling. Works with any regression estimator
#' (robust or classical).
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values.
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use (default: `tsreg` for Theil-Sen).
#'   Can be any regression function that returns a `coef` component (e.g.,
#'   `tsreg`, `ltsreg`, `opreg`, `qreg`).
#' @param nboot Number of bootstrap samples (default: 200).
#'
#' @return A list with components:
#'   \item{param.estimates}{Vector of regression parameter estimates (intercept
#'     and slopes)}
#'   \item{covar}{Covariance matrix of the parameter estimates. Diagonal elements
#'     are squared standard errors.}
#'   \item{s.e.}{Vector of standard errors for each parameter (square roots of
#'     covariance diagonal)}
#'
#' @details
#' This function uses bootstrap resampling to estimate the sampling variability
#' of regression parameters:
#' \enumerate{
#'   \item Fits the specified regression to the observed data
#'   \item Generates `nboot` bootstrap samples by resampling observations
#'   \item Refits the regression to each bootstrap sample
#'   \item Computes the sample covariance of bootstrap parameter estimates
#' }
#'
#' The bootstrap standard errors are appropriate for making inferences about
#' the regression parameters under general conditions (non-normality,
#' heteroscedasticity).
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' analysis using the function specified by `outfun`.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{tsreg}}, \code{\link{regbtci}}, \code{\link{regci}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Standard errors for Theil-Sen regression
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2 + 3*x[,1] - x[,2] + rnorm(50)
#' regse(x, y)
#'
#' # Standard errors for LTS regression
#' regse(x, y, regfun=ltsreg)
#' }
regse<-function(x,y,xout=FALSE,regfun=tsreg,outfun=outpro,nboot=200,SEED=TRUE,...){
#
#  Estimate the standard errors and
#  covariance matrix associated with the estimates of
#  the regression parameters based on the estimator indicated by the
#  argument
#  regfun:  default is Theil--Sen.
#  So the diagonal elements of the matrix returned by this function
#  are the squared standard errors of the intercept estimator, etc.
#
#  Function returns
# param.estimates: the estimate	of the intercept and slopes
# covar: the covariance matrix	associated with the estimator used
# s.e.:	 the standard errors.
#

if(SEED)set.seed(2)
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
nrem=length(y)
estit=regfun(x,y,xout=xout,...)$coef
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x<-as.matrix(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun,xout=FALSE,...)
#Leverage points already removed.
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
sqe=var(t(bvec))
list(param.estimates=estit,covar=sqe,s.e.=sqrt(diag(sqe)))
}


# ============================================================================
# regGmcp
# ============================================================================

#' Global Multiple Comparisons for Regression Parameters
#'
#' Performs all pairwise comparisons of regression parameters (slopes and/or
#' intercepts) among J independent groups using bootstrap inference and
#' Hochberg's method for familywise error rate (FWE) control.
#'
#' @inheritParams common-params
#' @param x A list of length J, where each element contains a numeric matrix or
#'   vector of predictor values for one group.
#' @param y A list of length J, where each element contains a numeric vector of
#'   response values for one group.
#' @param regfun Regression function to use (default: `tsreg` for Theil-Sen).
#' @param nboot Number of bootstrap samples (default: 100).
#' @param AD Logical; affects bootstrap method (default: `FALSE`).
#' @param STAND Logical; if `TRUE` (default), standardize when detecting outliers.
#' @param alpha Significance level for familywise error rate control (default: 0.05).
#' @param pr Logical; if `TRUE` (default), print warning about leverage points.
#' @param ISO Logical; if `TRUE` (default), compare only slopes (test for parallel
#'   regression lines). If `FALSE`, compare all parameters including intercepts.
#'
#' @return A list with components:
#'   \item{output}{Matrix containing all pairwise comparison results with adjusted
#'     p-values and decision outcomes}
#'   \item{n}{Vector of original sample sizes for each group (after removing missing values)}
#'   \item{n.keep}{Vector of sample sizes after outlier removal (if `xout = TRUE`)}
#'
#' @details
#' This function performs all pairwise global comparisons among J groups:
#' \itemize{
#'   \item For each pair of groups (j, k) with j < k, tests H0: all corresponding
#'     parameters are equal
#'   \item If `ISO = TRUE`: compares only slopes (intercepts ignored), testing for
#'     parallel regression lines
#'   \item If `ISO = FALSE`: compares all parameters (intercepts and slopes)
#'   \item Controls FWE using Hochberg's method across all (J^2-J)/2 pairwise tests
#' }
#'
#' The data structure requires x and y to be lists of length J, where x[[j]] and
#' y[[j]] contain the data for group j.
#'
#' If `xout = TRUE`, leverage points are removed from the predictor space for
#' each group using the function specified by `outfun`.
#'
#' @note For individual parameter comparisons (rather than global tests), use
#' `reg1mcp()` instead.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{reg1way}}, \code{\link{tsreg}}, \code{\link{regIVcom}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare regression lines for 3 groups
#' set.seed(123)
#' n <- 50
#' x1 <- rnorm(n); y1 <- 2 + 3*x1 + rnorm(n)
#' x2 <- rnorm(n); y2 <- 2 + 2*x2 + rnorm(n)  # Different slope
#' x3 <- rnorm(n); y3 <- 3 + 3*x3 + rnorm(n)  # Different intercept
#'
#' # Test for parallel lines (compare slopes only)
#' result <- regGmcp(list(x1, x2, x3), list(y1, y2, y3), ISO=TRUE)
#'
#' # Compare all parameters
#' result2 <- regGmcp(list(x1, x2, x3), list(y1, y2, y3), ISO=FALSE)
#' }
regGmcp<-function(x,y,regfun=tsreg,SEED=TRUE,nboot=100,xout=FALSE,AD=FALSE,
    outfun=outpro,STAND=TRUE,alpha=0.05,pr=TRUE,MC=FALSE,ISO=TRUE,...)
{
#
# If ISO = FALSE:
#  All pairwise comparisons of regression parameters are performed among J independent groups
#  That is, for groups j and k, all j<k, test H_0: all corresponding
#  parameters are equal.
#
#  ISO=TRUE:  compares all slopes, the intercept is ignored. So for a single covariate the
# goal is to test whether the regression lines are parallel.
#
#  For individual parameters, use reg1mcp

#  Perform all pairwise comparisons, where each comparison is based
#  on the global hypothesis that all parameters are equal
#
#  Control FWE via Hochberg's methods for each set of
#  (J^2-J)/2 parameters. That is, control FWE for the intercepts
#  Do the same for the first slope, etc.
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun,
#   which defaults to the projection method.
#
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#   output contains all pairwise comparisons.
#
if(!is.list(x))stop('Argument x should have list mode')
if(!is.list(y))stop('Argument y should have list mode')
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop('Something is wrong.
Number of covariates differs among the groups being compared')
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
if(!xout){
if(pr)print('Might want to consider removing any leverage points')
}
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
outp=matrix(NA,ncol=5,nrow=tot)
x=lapply(x,as.matrix)
xx=list()
yy=list()
iall=0
ivp=c(1,tot)-tot
i=0
for(j in 1:J){
for(k in 1:J){
if(j < k){
i=i+1
xx[[1]]=x[[j]]
xx[[2]]=x[[k]]
yy[[1]]=y[[j]]
yy[[2]]=y[[k]]
if(!ISO){
if(!MC)all=reg1way(xx,yy,regfun=regfun,nboot=nboot,SEED=SEED,AD=TRUE,alpha=alpha,
pr=FALSE,...)
if(MC)all=reg1wayMC(xx,yy,regfun=regfun,nboot=nboot,SEED=SEED,AD=TRUE,alpha=alpha,
pr=FALSE,...)
}
if(ISO){
if(!MC)all=reg1wayISO(xx,yy,regfun=regfun,nboot=nboot,SEED=SEED,AD=TRUE,alpha=alpha,
pr=FALSE,...)
if(MC)all=reg1wayISOMC(xx,yy,regfun=regfun,nboot=nboot,SEED=SEED,AD=TRUE,alpha=alpha,
pr=FALSE,...)
}
if(AD)temp=all$adjusted.p.value
if(!AD)temp=all$p.value
if(is.null(temp))temp=all$p.value
outp[i,1]=j
outp[i,2]=k
outp[i,3]=temp
}}
temp2<-order(0-outp[,3])
icc=c(1:tot)
icc[temp2]=dvec
outp[,4]=icc
}
flag=(outp[,3]<=outp[,4])
outp[,5]=rep(0,tot)
outp[flag,5]=1
dimnames(outp)=list(NULL,c('Group','Group','p.value','p.crit','sig'))
list(n=nv,n.keep=nv.keep,output=outp)
}


# ============================================================================
# regYciCV
# ============================================================================

#' Determine Critical Value for regYci via Simulation
#'
#' Computes a critical value for simultaneous confidence intervals in `regYci()`
#' via simulation. Used to control family-wise error rate when testing at
#' multiple design points simultaneously.
#'
#' @param n Sample size for simulation.
#' @inheritParams common-params
#' @param nboot Number of bootstrap samples for simulation (default: 1000).
#' @param regfun Regression function to use (default: `tsreg`).
#' @param null.value Null hypothesis value for testing (default: 0).
#' @param xout Logical; if TRUE, remove outliers from simulated data before
#'   regression (default: FALSE).
#' @param ... Additional arguments passed to `regfun`.
#'
#' @return Scalar critical value for use in `regYci(..., crit=)`.
#'
#' @details
#' This function simulates data under the null hypothesis to empirically
#' determine the critical value needed for simultaneous inference across
#' multiple design points. The critical value adjusts for the family-wise
#' error rate when making multiple comparisons.
#'
#' The simulation:
#' \itemize{
#'   \item Generates `nboot` datasets with uncorrelated bivariate normal data
#'   \item Fits regression at multiple points for each dataset
#'   \item Computes minimum p-values across all points
#'   \item Returns the alpha quantile of these minimum p-values
#' }
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regYci}}, \code{\link{regYband}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Determine critical value for n=50, alpha=0.05
#' crit <- regYciCV(n=50, alpha=0.05, nboot=500)
#'
#' # Use with regYci for simultaneous inference
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' regYci(x, y, crit=crit)
#' }
regYciCV<-function(n,alpha=.05,nboot=1000,regfun=tsreg,SEED=TRUE,MC=FALSE,null.value=0,xout=FALSE,...){
#
# Determine a critical value for regYci
#
if(SEED)set.seed(2)
mv=NA
chk=0
xy=list()
for (i in 1:nboot)xy[[i]]=rmul(n)
if(!MC)est=lapply(xy,regciCV.sub,regfun=regfun,null.value=null.value,...)
if(MC)est=mclapply(xy,regciCV.sub,regfun=regfun,null.value=null.value,...)
est=as.vector(matl(est))
est=sort(est)
ic=round(alpha*nboot)
crit=est[ic]
crit
}


# ============================================================================
# regYci.sum
# ============================================================================

#' Summarize Confidence Intervals for Regression Predictions
#'
#' Provides a more readable summary of confidence intervals for predicted Y
#' values from `regYci()`. Designed for single predictor models, organizing
#' results by x-values in ascending order.
#'
#' @inheritParams regYci
#' @param pts Points at which to compute confidence intervals (default: uses
#'   all values in `x`).
#'
#' @return A matrix with columns:
#'   \itemize{
#'     \item X: Predictor values (sorted in ascending order)
#'     \item Pred. Y: Predicted Y value at each x
#'     \item Lower.ci: Lower confidence limit
#'     \item Upper.ci: Upper confidence limit
#'     \item p.value: p-value for H0: predicted Y = null.value
#'   }
#'
#' @details
#' This is a wrapper for `regYci()` that:
#' \itemize{
#'   \item Restricts usage to single predictor models (p=1)
#'   \item Sorts output by predictor values for easier interpretation
#'   \item Returns results in matrix format with descriptive column names
#' }
#'
#' All other functionality (bootstrap confidence intervals, outlier removal,
#' adjustment for simultaneous inference) is identical to `regYci()`.
#'
#' @note This function will stop with an error if called with more than one
#'   independent variable.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regYci}}, \code{\link{regYband}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Confidence intervals at each x value
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' regYci.sum(x, y)
#'
#' # At specific points
#' regYci.sum(x, y, pts=seq(-2, 2, by=0.5))
#' }
regYci.sum<-function(x,y,regfun=tsreg,pts=x,nboot=100,xout=FALSE,outfun=out,SEED=TRUE,alpha=.05,crit=NULL,null.value=0,plotPV=FALSE,ADJ=FALSE,MC=FALSE,
scale=FALSE,span=.75,xlab='X',xlab1='X1',xlab2='X2',ylab='p-values',theta=50,phi=25,pch='*',...){
#
# Summarize results from regYci so that results are easier to read.
# single independent variable is assumed.
#
xy=elimna(cbind(x,y))
x<-as.matrix(x)
p=ncol(x)
if(p>1)stop('This function is designed for one independent variable only')
p1=p+1
x<-xy[,1:p]
y<-xy[,p1]
x<-as.matrix(x)
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
x=as.matrix(x)
}
res=regYci(x,y,regfun=regfun,nboot=nboot,null.value=null.value,alpha=alpha,crit=crit,
ADJ=ADJ,MC=MC,plotPV=plotPV,...)
xord=order(pts)
outp=cbind(pts[xord],res[xord,])
dimnames(outp)=list(NULL,c('X','Pred. Y','Lower.ci','Upper.ci','p.value'))
outp
}


# ============================================================================
# regYci
# ============================================================================

#' Confidence Intervals for Predicted Y Values in Regression
#'
#' Computes bootstrap confidence intervals for predicted Y values at specified
#' predictor values, with optional adjustment for simultaneous inference across
#' multiple points (confidence band).
#'
#' @inheritParams common-params
#' @param x A numeric vector or matrix containing predictor values (can be
#'   multivariate).
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use (default: `tsreg` for Theil-Sen).
#'   Can be any regression function that returns coefficients.
#' @param pts Points at which to compute confidence intervals. Default is
#'   `unique(x)` for all unique x values.
#' @param ADJ Logical. If `TRUE`, adjusts critical values so that the simultaneous
#'   probability coverage is 1-alpha (i.e., creates a confidence band). If `FALSE`,
#'   computes pointwise confidence intervals (default: `FALSE`).
#' @param crit Optional user-specified critical value. If `NULL`, computed from
#'   bootstrap (default: `NULL`).
#' @param nreps Number of replications for computing critical value when ADJ=TRUE
#'   (default: 1000). Higher values give more accurate critical values.
#' @param plotPV Logical. If `TRUE`, plots p-values vs. predictor values
#'   (default: `FALSE`).
#' @param SM Logical. If `TRUE`, applies smoothing to p-value plot (default: `FALSE`).
#' @param scale,span,theta,phi,pch Parameters for plotting.
#' @param xlab,xlab1,xlab2,ylab,zlab Axis labels for plots.
#'
#' @return A matrix with columns:
#'   \itemize{
#'     \item Predicted Y value at each point
#'     \item Lower confidence limit
#'     \item Upper confidence limit
#'     \item p-value for testing H0: predicted Y = null.value
#'   }
#'
#' @details
#' This function uses bootstrap methods to construct confidence intervals for
#' the predicted response at specified predictor values:
#' - With `ADJ = FALSE`: Pointwise confidence intervals (each with coverage 1-alpha)
#' - With `ADJ = TRUE`: Simultaneous confidence band (overall coverage 1-alpha)
#'
#' The adjustment for `ADJ=TRUE` has been primarily studied for single predictors.
#' Performance with multiple predictors is less well understood.
#'
#' For single predictors with `regfun=tsreg`, `ols`, or `qreg` and `alpha=0.05`,
#' the adjustment can be computed quickly. Otherwise, `nreps` bootstrap samples
#' are used to estimate the critical value, which can be time-consuming. Use
#' `MC=TRUE` with multiple cores to reduce computation time.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regYband}}, \code{\link{regYci.sum}}, \code{\link{regci}},
#'   \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Confidence intervals for predicted Y values
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' # Pointwise CIs
#' regYci(x, y, pts = c(-1, 0, 1))
#' # Simultaneous confidence band
#' regYci(x, y, pts = seq(-2, 2, by=0.5), ADJ = TRUE)
#' }
regYci<-function(x,y,regfun=tsreg,pts=unique(x),nboot=100,ADJ=FALSE,xout=FALSE,outfun=out,SEED=TRUE,alpha=.05,crit=NULL,null.value=0,plotPV=FALSE,scale=TRUE,span=.75,
xlab='X',xlab1='X1',xlab2='X2',ylab='p-values',zlab='p-values',
theta=50,phi=25,MC=FALSE,nreps=1000,SM=FALSE,pch='*',...){
#
#  Compute confidence interval for the typical value of Y, given X, based on some regression estimator
#  By default,
#  regfun=tsreg meaning that the  Theil--Sen estimator is used.
#
#  ADJ=TRUE,  the critical value is adjusted so that the simultaneous probability coverage is 1-alpha.
#  The adjustment has been studied with one independent variable. It is unknown how well it works with
#  more than one independent variable.
#
# If there is a single independent variable,
#  regfun=tsreg, ols  or qreg, and alpha=.05, an adjustment can be made quickly. Otherwise an
#  adjustment must be computed, which can require relatively high execution time.
#  To reduce execution time, set
#  MC=TRUE, assuming a multi-core processor is available.
#
# nreps: Number of replications used to compute a critical value. Execution time can be high
# MC=TRUE can reduce execution time considerably if a multi-core processor is available.
#

xy=elimna(cbind(x,y))
x<-as.matrix(x)
p=ncol(x)
p1=p+1
vals=NA
x<-xy[,1:p]
y<-xy[,p1]
x<-as.matrix(x)
n=nrow(x)
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
n=nrow(m)
x<-m[,1:p]
y<-m[,p1]
x=as.matrix(x)
}
if(ADJ){
if(n<10)stop('Should have a sample size of at least 10')
if(alpha==.05){
alpha=.01 # assuming tsreg,tsreg_C, tshdreg or qreg are being used.
if(identical(regfun,ols)){
nv=c(10,20,50,100,400)
pval=c(.001,.004, .008, .008, .01)
ipos=sum(nv<=n)
alpha=pval[ipos]
}
if(identical(regfun,tshdreg))alpha=.009
if(identical(regfun,qreg))alpha=.009
crit=qnorm(1-alpha/2)
}
}
if(SEED)set.seed(2)
if(is.null(crit)){
if(!ADJ)crit=qnorm(1-alpha/2)
if(ADJ){
padj=regYciCV(n,nboot=nreps,regfun=regfun,MC=MC,SEED=FALSE,
null.value=0,...)
crit=qnorm(1-padj/2)
}}
sqsd=regYvar(x,y,regfun=regfun,pts=pts,nboot=nboot,SEED=SEED,...)
sd=sqrt(sqsd)
est=regYhat(x,y,regfun=regfun,xr=pts,...)
pv=2*(1-pnorm(abs(est-null.value)/sd))
if(length(pts)==1)est=matrix(c(est,est-crit*sd,est+crit*sd,pv),nrow=1)
if(length(pts)>1)est=cbind(est,est-crit*sd,est+crit*sd,pv)
dimnames(est)=list(NULL,c('Pred. Y','Lower.ci','Upper.ci','p.value'))
if(plotPV){
if(ncol(x)>2)stop('Can plot only with one or two independent variables')
if(ncol(x)==1)plot(pts,pv,xlab=xlab,ylab=ylab,pch=pch)
if(ncol(x)==2){
if(SM)lplot(pts,pv,xlab=xlab1,ylab=xlab2,zlab=zlab,span=span,ticktype='detail',scale=scale,theta=theta,phi=phi)
if(!SM){
library(scatterplot3d)
scatterplot3d(pts[,1],pts[,2],pv,xlab=xlab1,ylab=xlab2,zlab=zlab)
}
}}
if(p==1){
xord=order(pts)
if(length(pts)==1)outp=matrix(c(pts[xord],est[xord,]),nrow=1)
if(length(pts)>1)outp=cbind(pts[xord],est[xord,])
dimnames(outp)=list(NULL,c('X','Pred. Y','Lower.ci','Upper.ci','p.value'))
est=outp
}
est
}


# ============================================================================
# regYciCV2G
# ============================================================================

#' Simulation-Based Critical Value for Two-Group Regression Comparison
#'
#' Determines critical values for `regYci2Gv2()` via simulation, used to control
#' family-wise error rate when comparing regression lines from two independent
#' groups at multiple design points simultaneously.
#'
#' @param n1 Sample size for first group.
#' @param n2 Sample size for second group.
#' @param crit Optional predetermined critical value; if provided, function
#'   computes Type I error rate at this value.
#' @param g Shape parameter for g-and-h distribution (default: 0, normal).
#' @param h Shape parameter for g-and-h distribution (default: 0, normal).
#' @inheritParams regYciCV
#' @param ALL Logical; if TRUE, use all unique x-values as test points; if
#'   FALSE, use `npts` equally-spaced points (default: TRUE).
#' @param pts Optional vector of specific design points for testing (default: NULL).
#' @param npts Number of design points if `ALL=FALSE` (default: 100).
#' @param nmiss Number of missing values to simulate in group 2 (default: 0).
#' @param ... Additional arguments passed to `regfun`.
#'
#' @return A list with components:
#'   \itemize{
#'     \item global.p.value: If `crit` provided, empirical Type I error rate
#'     \item crit.est: Estimated critical value at level `alpha` (Harrell-Davis
#'       quantile estimator)
#'   }
#'
#' @details
#' Simulates data from specified g-and-h distributions (bivariate normal when
#' g=h=0) to determine critical values for simultaneous inference when comparing
#' two regression lines. Accommodates:
#' \itemize{
#'   \item Unequal sample sizes between groups
#'   \item Missing data (via `nmiss` parameter)
#'   \item Non-normal distributions via g-and-h family
#' }
#'
#' The simulation generates paired datasets, computes regression predictions
#' at all/specified design points, and determines the minimum p-value across
#' points for each simulation run.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regYci2Gv2}}, \code{\link{regYciCV}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Critical value for equal groups, n=30 each
#' regYciCV2G(n1=30, n2=30, nboot=500)
#'
#' # Unequal groups with missing data
#' regYciCV2G(n1=40, n2=25, nmiss=5, nboot=500)
#'
#' # Non-normal data (g-and-h with skew and heavy tails)
#' regYciCV2G(n1=50, n2=50, g=0.5, h=0.2, nboot=500)
#' }
regYciCV2G<-function(n1,n2,crit=NULL,g=0,h=0,nboot=1000,regfun=tsreg,ALL=TRUE,
alpha=.05,SEED=TRUE,MC=FALSE,null.value=0,pts=NULL,npts=100,nmiss=0,...){
n=max(n1,n2)
if(nmiss>n)stop('Number of missing values is greater than max(n1,n2)')
if(SEED)set.seed(2)
mv=NA
chk=0
if(n1!=n2)nmiss=max(c(n1,n2))-min(c(n1,n2))
xy=list()
for (i in 1:nboot){
x1=ghdist(n,g=g,h=h)
x2=ghdist(n,g=g,h=h)
if(nmiss>0)x2[1:nmiss]=NA
xx=c(x1,x2)
xx=elimna(xx)
if(is.null(pts)){
if(!ALL)pts=seq(min(xx),max(xx),length.out = npts)
if(ALL)pts=unique(xx)
}
y1=ghdist(n,g=g,h=h)
y2=ghdist(n,g=g,h=h)
xy[[i]]=cbind(x1,y1,x2,y2)
}
if(!MC)est=lapply(xy,regciCV2G.sub,regfun=regfun,null.value=null.value,npts=npts,...)
if(MC)est=mclapply(xy,regciCV2G.sub,regfun=regfun,null.value=null.value,pts=pts,npts=npts,...)
est=as.vector(matl(est))
type1=NULL
if(!is.null(crit))type1=mean(est<=crit)
list(global.p.value=type1,crit.est=hd(est,alpha))
}


# ============================================================================
# regY2G.sub
# ============================================================================

#' Bootstrap Helper for Two-Group Regression Comparisons
#'
#' Internal helper function that extracts the minimum p-value from
#' \code{\link{regYci2Gv2}} for bootstrap-based simultaneous inference.
#'
#' @param xy Matrix with columns: x1, y1, x2, y2 (bootstrap resampled data).
#' @param regfun Regression function to use.
#' @param null.value Null hypothesis value.
#' @param ... Additional arguments passed to \code{regYci2Gv2}.
#'
#' @return Minimum p-value across all covariate values from \code{regYci2Gv2()}.
#'
#' @keywords internal
#' @export
regY2G.sub<-function(xy,regfun,null.value=0,...){
 pv=regYci2Gv2(xy[,1],xy[,2],xy[,3],xy[,4],SEED=FALSE,regfun=regfun,null.value=null.value,...)[,4]
 min(pv)
}


# ============================================================================
# regYci2Gv2
# ============================================================================

#' Two-Group Comparison of Regression Predictions (ANCOVA-Style)
#'
#' Computes confidence intervals for the difference between predicted Y values
#' from two independent groups at specified covariate values, providing an
#' ANCOVA-style analysis using robust regression methods.
#'
#' @param x1,y1 Numeric vectors containing predictor and response for group 1.
#' @param x2,y2 Numeric vectors containing predictor and response for group 2.
#' @param regfun Regression function to use (default: `tsreg` for Theil-Sen).
#' @param pts Covariate values at which to compare groups. If `NULL`, uses
#'   either all unique covariate values (`ALL=TRUE`) or evenly spaced points
#'   (`ALL=FALSE`).
#' @param ALL Logical. If `TRUE`, evaluates at all unique covariate values;
#'   if `FALSE`, uses `npts` evenly spaced points (default: `FALSE`).
#' @param npts Number of evenly spaced evaluation points when `ALL=FALSE`
#'   (default: 25).
#' @param plotit Logical. If `TRUE`, plots data and fitted regression lines
#'   for both groups (default: `TRUE`).
#' @param SCAT Logical. If `TRUE`, includes scatterplot of data points;
#'   if `FALSE`, shows only regression lines (default: `TRUE`).
#' @param pch1,pch2 Plot characters for groups 1 and 2 (defaults: '*', '+').
#' @inheritParams regYci
#' @param ylab2 Y-axis label for main plot (default: 'Y').
#' @param p.crit Critical p-value for simultaneous testing when `ADJ=TRUE`
#'   (default: 0.015).
#' @param nreps Number of replications for computing adjusted critical value
#'   (default: 1000).
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: X, Est.Dif (estimated difference),
#'     Lower.ci, Upper.ci, p.value}
#'   \item{p.crit}{Critical p-value used for testing.}
#'   \item{crit.value}{Critical value used for confidence intervals.}
#'   \item{num.sig}{Number of covariate values where groups differ significantly.}
#'
#' @details
#' This function performs ANCOVA-style analysis by comparing regression functions
#' between two groups:
#' \enumerate{
#'   \item Fits robust regression for each group separately
#'   \item At each covariate value, estimates the difference in predicted Y
#'   \item Computes bootstrap confidence intervals for the difference
#'   \item Optionally adjusts for simultaneous testing across multiple covariate values
#' }
#'
#' With `ADJ=TRUE`, the critical value is adjusted so that the overall family-wise
#' error rate is controlled at level alpha. This adjustment uses a method that
#' typically has better power than Hochberg or Hommel adjustments.
#'
#' For `regfun=tsreg` and `alpha=0.05` with single covariate, a pre-computed
#' adjustment is used (fast). Otherwise, the critical value is computed via
#' simulation (can be slow; use `MC=TRUE` to speed up with multicore processing).
#'
#' @note This function is identical to `anclin()`. Currently restricted to
#'   single covariate (p=1) models.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regYci}}, \code{\link{ancJN}}, \code{\link{regYband}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two groups with ANCOVA-style analysis
#' set.seed(123)
#' x1 <- rnorm(50)
#' y1 <- 2*x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- 2*x2 + 1 + rnorm(50)  # Same slope, different intercept
#'
#' # Pointwise comparisons
#' result <- regYci2Gv2(x1, y1, x2, y2)
#'
#' # Simultaneous confidence band
#' result <- regYci2Gv2(x1, y1, x2, y2, ADJ=TRUE)
#' }
regYci2Gv2<-function(x1,y1,x2,y2,regfun=tsreg,pts=NULL,ALL=FALSE,npts=25,plotit=TRUE,SCAT=TRUE,
pch1='*',pch2='+',
nboot=100,ADJ=FALSE,xout=FALSE,outfun=outpro,SEED=TRUE,p.crit=.015,
alpha=.05,crit=NULL,null.value=0,plotPV=FALSE,scale=TRUE,span=.75,xlab='X',xlab1='X1',xlab2='X2',ylab='p-values',ylab2='Y',theta=50,phi=25,MC=FALSE,nreps=1000,pch='*',...){
#
#  ANCOVA:
#  For two independent groups, compute confidence intervals for difference between
#  the typical value of Y, given X,
#   based on some regression estimator
#  By default,
#  regfun=tsreg meaning that the  Theil--Sen estimator is used.
#
#  The functions anclin and regYci2g are identical.
#
#  In contrast to the function ancJN, this function can deal with a larger number of
#  covariate values and it controls the probability of one or more Type I errors using
#  a method that is better, in terms of power, than using Hochberg or Hommel.
#
#  ADJ=TRUE,  the critical value is adjusted so that the simultaneous
#   probability coverage is 1-alpha.
#   If there is a single covariate,
#  regfun=tsreg or tshdreg, and alpha=.05, an adjustment can be made quickly. Otherwise an
#  adjustment must be computed, which can require relatively high execution time.
#  To reduce execution time, set
#  MC=TRUE, assuming a multi-core processor is available.
#  If n1<20 and n2<100, assuming that n1<n2,
#  an adjusted critical value must be computed even when using the
#  Theil--Sen estimator.
#
#  nboot: the number of bootstrap samples used to estimate standard errors.
#
# nreps: Number of replications used to compute a critical value.
#
#  pts: values for the independent variable where confidence intervals are computed
#  pts=NULL means that 100 points evenly spaced between min(x1,x2) and max(x1,x2) are used.
#  If pts is specified,the function and ADJ=TRUE, the function will  compute an
#  adjusted critical value, which again can result in high execution time.
#
xy=elimna(cbind(x1,y1))
x1<-as.matrix(x1)
p=ncol(x1)
if(p>1)stop('Current version allows one covariate only')
p1=p+1
vals=NA
x1<-xy[,1:p]
y1<-xy[,p1]
x1<-as.matrix(x1)
xy=elimna(cbind(x2,y2))
x2<-as.matrix(x2)
p=ncol(x2)
p1=p+1
vals=NA
x2<-xy[,1:p]
y2<-xy[,p1]
x2<-as.matrix(x2)
n1=length(y1)
n2=length(y2)
n=min(c(n1,n2))
#print(c(n1,n2,n))
if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE)$keep
m<-m[flag,]
n1=nrow(m)
x1<-m[,1:p]
y1<-m[,p1]
x1=as.matrix(x1)
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE)$keep
m<-m[flag,]
n2=nrow(m)
n=min(c(n1,n2))
x2<-m[,1:p]
y2<-m[,p1]
x2=as.matrix(x2)
}
if(is.null(pts)){
xall=unique(c(x1,x2))
if(ALL)pts=xall
if(!ALL)pts=seq(min(xall),max(xall),length.out=npts)
}
if(ADJ){
if(n<10)stop('Should have a sample size of at least 10')
if(alpha==.05){
#if(identical(regfun,tsreg) || identical(regfun,tsreg_C))alpha=p.crit causes an error if WRScpp not installed
#if(identical(regfun,tsreg))alpha=p.crit
alpha=p.crit
crit=qnorm(1-alpha/2)
}
if(!ADJ)p.crit=alpha
if(n<20 & max(c(n1,n2))<100) crit=NULL
if(p>1)crit=NULL
}
if(is.null(crit) & !ADJ)crit=qnorm(1-alpha/2)
if(is.null(crit) & ADJ){
if(SEED)set.seed(2)
print(c(n1,n2))
padj=regYciCV2G(n1,n2,nboot=nreps,regfun=regfun,MC=MC,SEED=FALSE,ALL=ALL,
null.value=null.value,pts=pts,alpha=alpha,...)$crit.est
crit=qnorm(1-padj/2)
p.crit=padj
}
sqsd1=regYvar(x1,y1,regfun=regfun,pts=pts,nboot=nboot,SEED=SEED)
sqsd2=regYvar(x2,y2,regfun=regfun,pts=pts,nboot=nboot,SEED=SEED)
sd=sqrt(sqsd1+sqsd2)
est1=regYhat(x1,y1,regfun=regfun,xr=pts,...)
est2=regYhat(x2,y2,regfun=regfun,xr=pts,...)
pv=2*(1-pnorm(abs(est1-est2-null.value)/sd))
est=cbind(pts,est1-est2,est1-est2-crit*sd,est1-est2+crit*sd,pv)
dimnames(est)=list(NULL,c('X','Est.Dif','Lower.ci','Upper.ci','p.value'))
if(plotit){
plotPV=FALSE
plot(c(x1,x2),c(y1,y2),type="n",xlab=xlab,ylab=ylab2)
reg1=regfun(x1,y1,...)$coef
reg2=regfun(x2,y2,...)$coef
if(SCAT){
points(x1,y1,pch=pch1)
points(x2,y2,pch=pch2)
}
abline(reg1)
abline(reg2,lty=2)
}
if(plotPV){
if(ncol(x1)>2)stop('Can plot only with one or two independent variables')
if(ncol(x1)==1)plot(pts,pv,xlab=xlab,ylab=ylab,pch=pch)
if(ncol(x2)==2)lplot(pts,pv,xlab=xlab1,ylab=xlab2,zlab=ylab,span=span,ticktype='detail',scale=scale,theta=theta,phi=phi)
}
list(output=est,p.crit=p.crit,crit.value=crit,num.sig=sum(est[,5]<=p.crit))
}


# ============================================================================
# MULMreg
# ============================================================================
#' Multivariate Multiple Regression
#'
#' Performs multivariate multiple regression where the response consists of
#' multiple dependent variables. Estimates regression parameters for each
#' response variable separately using a specified regression estimator.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values (must have 2 or more columns).
#' @param y A numeric matrix containing response variables (must have 2 or more
#'   columns, one for each response).
#' @param regfun Regression function to use for each response variable
#'   (default: `MMreg` for MM-regression). Can be any regression function that
#'   returns a `coef` component (e.g., `tsreg`, `ltsreg`, `opreg`).
#'
#' @return A list with component:
#'   \item{coef}{Matrix of regression coefficients with rows corresponding to
#'     parameters (intercept, slopes) and columns corresponding to response
#'     variables. Row names indicate parameter type.}
#'
#' @details
#' This function fits separate regression models for each column of `y`:
#' \itemize{
#'   \item Each response variable is regressed on all predictors in `x`
#'   \item The same regression estimator is used for all responses
#'   \item Results are organized in a coefficient matrix
#' }
#'
#' This is different from truly multivariate regression (which models the
#' covariance structure among responses). Here, each response is modeled
#' independently.
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' analysis using the function specified by `outfun`. The same observations
#' are removed for all response variables.
#'
#' @note Both `x` and `y` must be matrices with at least 2 columns each.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{mlrreg}}, \code{\link{MMreg}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Multivariate regression with 2 predictors and 3 responses
#' set.seed(123)
#' n <- 100
#' x <- matrix(rnorm(n*2), ncol=2)
#' y <- matrix(rnorm(n*3), ncol=3)
#' y[,1] <- 2 + 3*x[,1] - x[,2] + rnorm(n)
#' y[,2] <- -1 + x[,1] + 2*x[,2] + rnorm(n)
#' y[,3] <- 5 - 2*x[,1] + rnorm(n)
#'
#' # Fit using MM-regression
#' result <- MULMreg(x, y)
#' result$coef
#'
#' # Use Theil-Sen instead
#' result2 <- MULMreg(x, y, regfun=tsreg)
#' }
MULMreg<-function(x,y,regfun=MMreg,
xout=FALSE,outfun=outpro,...){
#
# Multivariate regression: simply estimate parameters for
# for each column of Y values based on some multivariate regression
# estimator.
#
#  Use MMreg by default
#
# x and y are assumed to be matrices with two or more columns
#
#
x<-as.matrix(x)
y<-as.matrix(y)
n.keep=nrow(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=FALSE)$keep
x<-x[flag,]
y<-y[flag,]
x<-as.matrix(x)
y<-as.matrix(y)
n.keep=nrow(x)
}
p1=ncol(x)+1
q=ncol(y)
est=matrix(NA,nrow=p1,ncol=q)
dimnames(est)=list(c('Inter',rep('Slope',ncol(x))),NULL)
for(i in 1:q)est[,i]=regfun(x,y[,i],...)$coef
list(coef=est)
}


# ============================================================================
# regIVcom
# ============================================================================

#' Compare Strength of Association for Two Predictor Subsets
#'
#' Compares the strength of association (explanatory power) for two subsets of
#' independent variables in a regression model using bootstrap inference.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values.
#' @param y A numeric vector containing the response variable.
#' @param IV1 Vector of column indices in `x` indicating the first subset of
#'   predictors to compare (default: 1, the first predictor).
#' @param IV2 Vector of column indices in `x` indicating the second subset of
#'   predictors to compare (default: 2, the second predictor).
#' @param regfun Regression function to use (default: `tsreg`). Using `Qreg`
#'   reduces execution time but may reduce power.
#' @param nboot Number of bootstrap samples (default: 200).
#' @param tr Trimming proportion for Winsorized variance (default: 0.2). Controls
#'   the amount of Winsorizing when computing explanatory power.
#'
#' @return A list with components:
#'   \item{n}{Original sample size}
#'   \item{n.keep}{Sample size after outlier removal (if `xout = TRUE`)}
#'   \item{est.1}{Estimated strength measure for predictor subset 1 (Winsorized
#'     variance of fitted values)}
#'   \item{est.2}{Estimated strength measure for predictor subset 2}
#'   \item{e.pow1}{Explanatory power for subset 1 (proportion of Y variance explained)}
#'   \item{e.pow2}{Explanatory power for subset 2}
#'   \item{strength.assoc.1}{Strength of association for subset 1 (sqrt of explanatory power)}
#'   \item{strength.assoc.2}{Strength of association for subset 2}
#'   \item{ratio}{Ratio of strength measures (est.1 / est.2)}
#'   \item{strength.ratio}{Ratio of association strengths (sqrt of ratio)}
#'   \item{p.value}{Bootstrap p-value for testing H0: equal strength}
#'
#' @details
#' This function tests whether two subsets of predictors have equal strength of
#' association with the response:
#' \enumerate{
#'   \item Fits regression using all predictors
#'   \item For each bootstrap sample, computes strength measures for both subsets
#'     based on Winsorized variance of fitted values
#'   \item Compares bootstrap distributions to test H0: equal strength
#' }
#'
#' The subsets IV1 and IV2 must be disjoint (no overlap). For example, `IV1=c(1,2)`
#' and `IV2=3` compares the first two predictors jointly to the third predictor.
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' analysis using the function specified by `outfun`.
#'
#' @note The subsets IV1 and IV2 must not have duplicate values and their combined
#' length cannot exceed the number of columns in `x`.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regIVstr}}, \code{\link{regIVcommcp}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare first predictor to second predictor
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol=3)
#' y <- 2 + 3*x[,1] - 0.5*x[,2] + 0.1*x[,3] + rnorm(100)
#'
#' # Test if first predictor has stronger association than second
#' regIVcom(x, y, IV1=1, IV2=2)
#'
#' # Compare first two predictors jointly to third
#' regIVcom(x, y, IV1=c(1,2), IV2=3)
#' }
regIVcom<-function(x,y,IV1=1,IV2=2,regfun=tsreg,nboot=200,xout=FALSE,outfun=outpro,SEED=TRUE,MC=FALSE,tr=.2,...){
#
# Compare strength of the association for two subsets of independent variables.
# IV1 and IV2 indicate the two sets of independent variables to be compared
# Example: IV1=c(1,2), IV2=3 would compare the first two independent variables to the third.
# Explanatory power is used based on a Winsorized variance.
# tr indicates the amount of Winsorizing
#
#  regfun=Qreg reduces execution time but possibly at the expense of less power.
#
if(sum(duplicated(c(IV1,IV2)))>0)stop('IV1 and IV2 have duplicate values making this method invalid')
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
#if(length(IV1)+length(IV2) != p)stop('ncol(x) should equal the number of variables indicated by IV1 and IV2')
if(length(IV1)+length(IV2) > p)stop('IV!+IV2 should be less than or equal ncol(x)')
if(max(c(IV1,IV2))>p)stop('IV1 or IV2 has a value that exceeds the number of col. in x')
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
nrem=length(y)
if(xout){
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
nkeep=length(y)
#estit=regfun(x,y,xout=xout,...)$coef[2:p1]
nv=length(y)
x<-as.matrix(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(k in 1:2){
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
if(!MC){
bvec<-apply(data,1,regboot,x,y,regfun,xout=FALSE,...)
if(k==1)bvec1=bvec
if(k==2)bvec2=bvec
}
if(MC){
data=listm(t(data))
bvec<-mclapply(data,regbootMC,x,y,regfun,mc.preschedule=TRUE,xout=FALSE,...)
if(k==1)bvec1=matl(bvec)
if(k==2)bvec2=matl(bvec)
data=t(matl(data))
}}
#Leverage points already removed.
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec1=bvec1[2:p1,]  # don't need the intercept
bvec2=bvec2[2:p1,]  # don't need the intercept
v1=NA
v2=NA
for(i in 1:nboot){
v1[i]=regIVcom_sub(bvec1[IV1,i],x[data[i,],IV1],tr=tr)
v2[i]=regIVcom_sub(bvec2[IV2,i],x[data[i,],IV2],tr=tr)
}
pv=bmp(v1,v2)$phat
pv=2*min(c(pv,1-pv))
est=regfun(x,y)$coef[2:p1]
e1=regIVcom_sub(est[IV1],x[,IV1],tr=tr)
e2=regIVcom_sub(est[IV2],x[,IV2],tr=tr)
rat=NA
if(e2>0)rat=e1/e2
ep1=e1/winvar(y,tr=tr)
ep2=e2/winvar(y,tr=tr)
list(n=nrem,n.keep=nkeep,est.1=e1,est.2=e2,e.pow1=ep1,e.pow2=ep2,strength.assoc.1=sqrt(ep1),
strength.assoc.2=sqrt(ep2),
ratio=rat,strength.ratio=sqrt(rat),p.value=pv)
}


# ============================================================================
# regIVcom_sub
# ============================================================================

#' Bootstrap Helper for regIVcom
#'
#' Internal helper function used by \code{\link{regIVcom}} to compute Winsorized
#' variance of fitted values for bootstrap samples when comparing predictor subsets.
#'
#' @param slope Bootstrap sample of slope coefficients.
#' @param x Matrix of predictors.
#' @param tr Trimming proportion for Winsorized variance.
#'
#' @return Winsorized variance of fitted values.
#'
#' @keywords internal
#' @export
regIVcom_sub<-function(slope,x,tr){
yhat=apply(t(slope*t(x)),1,sum)
str=winvar(yhat,tr=tr)
str
}


# ============================================================================
# regIVstr
# ============================================================================

#' Estimate Strength of Association for Each Predictor
#'
#' Estimates the individual strength of association (explanatory power) for each
#' predictor when all predictors are included in the regression model. Provides
#' a measure of each predictor's unique contribution.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values.
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use (default: `tsreg` for Theil-Sen).
#' @param tr Trimming proportion for Winsorized variance (default: 0.2). Controls
#'   the amount of Winsorizing when computing explanatory power.
#'
#' @return A list with components:
#'   \item{explanatory.power}{Vector of explanatory power values for each predictor,
#'     measured as the proportion of Y variance explained by that predictor's
#'     fitted values (based on Winsorized variance)}
#'   \item{explanatory.strength}{Vector of association strengths (square roots of
#'     explanatory power values)}
#'
#' @details
#' This function quantifies how much each predictor contributes to explaining
#' the response when all predictors are in the model:
#' \enumerate{
#'   \item Fits regression using all predictors
#'   \item For each predictor j, computes fitted values using only slope[j] * x[,j]
#'   \item Measures explanatory power as: Winsorized variance of fitted values /
#'     Winsorized variance of Y
#'   \item Returns both explanatory power and its square root (strength)
#' }
#'
#' This differs from univariate regression strength because it uses slopes from
#' the full model (adjusting for other predictors) rather than simple regression
#' slopes.
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' analysis using the function specified by `outfun`.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regIVcom}}, \code{\link{regIVcommcp}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Estimate strength for each predictor
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol=3)
#' y <- 2 + 3*x[,1] - 0.5*x[,2] + 0.1*x[,3] + rnorm(100)
#'
#' result <- regIVstr(x, y)
#' result$explanatory.power    # Proportion of variance explained
#' result$explanatory.strength # Strength of association
#' }
regIVstr<-function(x,y,regfun=tsreg,xout=FALSE,outfun=outpro,tr=.2,...){
#
# Estimate strength of each independent variable
# when all of them are  entered into the model.
#
xy=cbind(x,y)
xy=elimna(xy)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=FALSE)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
p=ncol(x)
p1=p+1
est=regfun(x,y,...)$coef[2:p1]
top=NA
for(j in 1:p)top[j]=regIVcom_sub(est[j],x[,j],tr=tr)
bot=winvar(y,tr=tr)
str=top/bot
list(explanatory.power=str,explanatory.strength=sqrt(str))
}


# ============================================================================
# qhdsm.pred
# ============================================================================
#' Predictions from Quantile Regression Smoother
#'
#' Predicts the qth quantile of Y at specified predictor values using the
#' quantile regression smoother based on the Harrell-Davis estimator. This
#' function provides predictions from the model fitted by `qhdsm()`.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values from the
#'   training data.
#' @param y A numeric vector containing response values from the training data.
#' @param pts A matrix or vector of predictor values at which predictions are
#'   desired (default: `x`, predicting at training points).
#' @param q Quantile to estimate (default: 0.5 for median). Must be between 0 and 1.
#' @param fr Fraction/span parameter controlling the neighborhood size for local
#'   quantile estimation (default: 1). Larger values use more neighbors.
#' @param nmin Minimum number of observations required in a neighborhood for
#'   prediction (default: 1).
#'
#' @return A list with components:
#'   \item{Y.hat}{Vector of predicted Y values (qth quantiles) at each point in `pts`}
#'   \item{nvals}{Vector indicating the number of observations used for each prediction}
#'
#' @details
#' For each point in `pts`, this function:
#' \enumerate{
#'   \item Identifies nearby observations in the training data (x, y) using the
#'     `near()` function with span `fr`
#'   \item Computes the qth quantile of Y values in that neighborhood using the
#'     Harrell-Davis estimator (`hd()`)
#'   \item Returns the predicted quantile and neighborhood size
#' }
#'
#' For single-predictor problems, predictions use `runhat()`. For multiple
#' predictors, `rung3hat()` is used instead.
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' prediction using the function specified by `outfun`.
#'
#' @seealso \code{\link{qhdsm}}, \code{\link{qhdsm2g}}, \code{\link{hd}}, \code{\link{runhat}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Fit and predict median regression
#' set.seed(123)
#' x <- runif(100, 0, 10)
#' y <- 2 + 0.5*x + rnorm(100)
#'
#' # Predict at new points
#' newx <- seq(0, 10, length=20)
#' pred <- qhdsm.pred(x, y, pts=newx, q=0.5)
#' plot(x, y)
#' lines(newx, pred$Y.hat)
#'
#' # Predict 75th percentile
#' pred75 <- qhdsm.pred(x, y, pts=newx, q=0.75)
#' lines(newx, pred75$Y.hat, lty=2)
#' }
qhdsm.pred<-function(x,y,pts=x,q=.5,fr=1,nmin=1,xout=FALSE,outfun=outpro,...){
#
#  Predict the qth quantile of Y based on the values in pts, using the
#  the data in x and y.
#
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(ncol(x)==1){
vals=runhat(x[,1],y,pts=pts,est=hd,q=q,fr=fr,nmin=nmin,...)
nvals=1
for(i in 1:length(pts)){
nvals[i]<-length(y[near(x,pts[i],fr=fr)])
}
}
if(ncol(x)>1){
temp=rung3hat(x,y,pts=pts,est=hd,q=q,fr=fr,...)
vals=temp$rmd
nvals=temp$nval
}
list(Y.hat=vals,nvals=nvals)
}


# ============================================================================
# regbtci
# ============================================================================

#' Bootstrap-t Confidence Intervals for Regression Parameters
#'
#' Computes bootstrap-t confidence intervals and hypothesis tests for regression
#' parameters (intercept and slopes). Works with any regression estimator.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values.
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use (default: `qreg` for quantile regression).
#'   Can be any regression function that returns a `coef` component.
#' @param alpha Significance level for confidence intervals (default: 0.05 for
#'   95% confidence intervals).
#' @param nboot Number of bootstrap samples for estimating standard errors
#'   (default: 300).
#'
#' @return A list with components:
#'   \item{output}{Matrix with rows for each parameter (intercept, slopes) and
#'     columns: `ci.low` (lower CI bound), `ci.up` (upper CI bound), `Estimate`
#'     (point estimate), `S.E.` (standard error), `p-value` (two-sided test of
#'     H0: parameter = 0)}
#'   \item{n}{Original sample size before outlier removal}
#'   \item{n.keep}{Sample size after outlier removal (if `xout = TRUE`)}
#'
#' @details
#' This function uses the bootstrap-t method for inference on regression
#' parameters:
#' \enumerate{
#'   \item Estimates standard errors via bootstrap (`regse()`)
#'   \item Constructs confidence intervals using normal approximation:
#'     estimate ± z(1-α/2) × SE
#'   \item Computes p-values for testing each parameter = 0 using z-statistics
#' }
#'
#' The bootstrap-t approach provides better coverage than percentile bootstrap
#' in many situations, especially with smaller samples.
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' analysis using the function specified by `outfun`.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regse}}, \code{\link{regci}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Bootstrap-t CIs for quantile regression
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2 + 3*x[,1] - x[,2] + rnorm(50)
#' regbtci(x, y)
#'
#' # Use with Theil-Sen regression
#' regbtci(x, y, regfun=tsreg)
#' }
regbtci<-function(x,y,regfun=qreg,alpha=.05,nboot=300,xout=FALSE,outfun=outpro,SEED=TRUE,...){
#
#  Bootstrap-t confidence intervals for regression parameters
#
if(SEED)set.seed(2)
xx<-elimna(cbind(x,y))
np<-ncol(xx)
p<-np-1
y<-xx[,np]
x<-xx[,1:p]
x<-as.matrix(x)
n.orig=length(y)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
}
vlabs='Intercept'
for(j in 2:np)vlabs[j]=paste('Slope',j-1)
regout<-matrix(0,np,5)
dimnames(regout)<-list(vlabs,c('ci.low','ci.up','Estimate','S.E.','p-value'))
val=regse(x,y,regfun=regfun,nboot=nboot,SEED=SEED,...)
tests=val$param.estimates/val$s.e.
pv=2*(1-pnorm(abs(tests)))
est=regfun(x,y,...)
regout[,3]=est$coef
regout[,1]=est$coef-qnorm(1-alpha/2)*val$s.e.
regout[,2]=est$coef+qnorm(1-alpha/2)*val$s.e.
regout[,4]=val$s.e.
regout[,5]=pv
list(output=regout,n=n.orig,n.keep=length(y))
}


# ============================================================================
# regIVcommcp
# ============================================================================

#' Multiple Comparisons of Predictor Strength for All Pairs
#'
#' For each pair of independent variables (predictors) in a multiple regression,
#' compares their relative strength of association when both are included in the
#' model. Performs all pairwise comparisons among predictors.
#'
#' @inheritParams common-params
#' @param x A numeric matrix containing predictor values. Must have at least 2 columns.
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use (default: `tsreg`).
#' @param nboot Number of bootstrap samples for inference (default: 200).
#' @param tr Trimming proportion for computing strength of association (default: 0.2).
#'
#' @return A matrix with one row for each pair of predictors and columns:
#'   \item{IV 1}{Index of first predictor in the pair.}
#'   \item{IV 2}{Index of second predictor in the pair.}
#'   \item{strength.assoc.1}{Strength of association for first predictor.}
#'   \item{strength.assoc.2}{Strength of association for second predictor.}
#'   \item{strength.ratio}{Ratio of strengths (predictor 1 / predictor 2).}
#'   \item{p.value}{P-value for testing H0: equal strength of association.}
#'
#' @details
#' This function performs all pairwise comparisons of predictor importance:
#'
#' For each pair of predictors (i, j):
#' \enumerate{
#'   \item Fits regression with both predictors: Y ~ X_i + X_j
#'   \item Computes strength of association for each predictor separately
#'   \item Tests whether the two predictors have equal strength via bootstrap
#'   \item Returns strength estimates, ratio, and p-value
#' }
#'
#' With p predictors, there are p*(p-1)/2 pairwise comparisons.
#'
#' Strength of association measures how much variance in Y is explained by each
#' predictor when the other is also in the model. This is computed robustly using
#' trimmed variance (controlled by `tr`).
#'
#' This function is useful for:
#' \itemize{
#'   \item Identifying which predictors are most important
#'   \item Detecting redundant predictors (similar strength)
#'   \item Understanding predictor contributions in multiple regression
#' }
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regIVcom}}, \code{\link{regIVstr}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare strength of 3 predictors
#' set.seed(123)
#' x <- matrix(rnorm(150), ncol=3)
#' y <- 2*x[,1] + 0.5*x[,2] + rnorm(50)  # X1 stronger than X2, X3 irrelevant
#' result <- regIVcommcp(x, y, nboot=100)
#' result  # Shows all pairwise comparisons
#' }
regIVcommcp<-function(x,y,regfun = tsreg, nboot = 200,
    xout = FALSE, outfun = outpro, SEED = TRUE, MC = FALSE, tr = 0.2,
    ...){
#
#  For each pair of the independent variables in x, compare strength
#  when both are included in the model.
#
x<-as.matrix(x)
J=ncol(x)
A=(J^2-J)/2
output=matrix(NA,nrow=A,ncol=6)
ic=0
for(i in 1:J){
for(k in 1:J){
if(i<k){
ic=ic+1
res=regIVcom(x[,c(i,k)],y,regfun=regfun,nboot=nboot,xout=xout,
outfun=outfun,SEED=SEED,MC=MC,tr=tr,...)
output[ic,1:2]=c(i,k)
output[ic,3:6]=c(res$strength.assoc.1,res$strength.assoc.2,res$strength.ratio,
res$p.value)
}}}
dimnames(output)=list(NULL,c('IV 1','IV 2','strength.assoc.1','strength.assoc.2',
'strength.ratio','p.value'))
output
}


# ============================================================================
# regstr
# ============================================================================

#' Compute Explanatory Strength of Association
#'
#' Computes a robust measure of explanatory strength of association, analogous
#' to R-squared, using the ratio of variances of fitted to observed values.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values.
#' @param y A numeric vector containing the response variable.
#' @param varfun Variance function for computing strength (default: `winvar`
#'   for Winsorized variance).
#' @param regfun Regression function to use (default: `tsreg`).
#'
#' @return A numeric value between 0 and 1 representing the explanatory strength.
#'   Values close to 1 indicate strong association, values close to 0 indicate
#'   weak association.
#'
#' @details
#' This function computes a robust analog of R-squared:
#'
#' \deqn{Strength = Var(\hat{Y}) / Var(Y)}
#'
#' where Var() is computed using a robust variance estimator (default: Winsorized
#' variance).
#'
#' Unlike classical R-squared which uses least squares and can exceed 1 with
#' robust regression methods, this measure uses robust variance estimation for
#' both fitted values and observed values, providing resistance to outliers.
#'
#' Interpretation:
#' \itemize{
#'   \item 0: No explanatory power (predictors don't explain Y)
#'   \item 1: Perfect explanatory power (predictors fully explain Y)
#'   \item Intermediate values: Proportion of robust variance explained
#' }
#'
#' If `xout = TRUE`, outliers in the predictor space are removed before
#' computing the regression and strength.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regIVstr}}, \code{\link{tsreg}}, \code{\link{winvar}},
#'   \code{\link{mgvreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compute explanatory strength
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2*x[,1] + 3*x[,2] + rnorm(50)
#' regstr(x, y)  # Should be high (strong relationship)
#'
#' # Weak relationship
#' y2 <- rnorm(50)
#' regstr(x, y2)  # Should be low (no relationship)
#' }
regstr<-function(x,y,varfun=winvar,regfun=tsreg,xout=FALSE,outfun=outpro,...){
#
# Exlanatory strength of association; similar to R^2.
#
yhat=regYhat(x,y,regfun=regfun,xout=xout,outfun=outpro,...)
str=varfun(yhat)/varfun(y)
str
}


# ============================================================================
# multireg.prob
# ============================================================================

#' Multinomial Logistic Regression Probabilities
#'
#' Estimates conditional probabilities P(Y=k|X) for multinomial (categorical)
#' outcomes using a multinomial logit model. Optionally plots the predicted
#' probabilities as a function of predictors.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values.
#' @param y A numeric or factor vector containing the categorical response variable
#'   with multiple possible values.
#' @param pts Matrix of predictor values at which to estimate probabilities
#'   (default: `x`, the observed predictor values).
#' @param plotit Logical. If `TRUE`, creates a plot of predicted probabilities
#'   (default: `TRUE`).
#' @param xlab,ylab,zlab Axis labels for plots.
#' @param ticktype Type of tick marks for 3D plot (default: "det").
#' @param vplot Numeric vector specifying which response categories to plot.
#'   If `NULL`, plots the maximum value of Y (default: `NULL`).
#' @param L Logical. For 2 predictors: if `TRUE`, uses `lplot()`; if `FALSE`,
#'   uses `rplot()` (default: `TRUE`).
#' @param scale Logical. For 1 predictor: if `TRUE`, y-axis ranges 0-1; if `FALSE`,
#'   y-axis matches range of probabilities (default: `TRUE`).
#'
#' @return A list with components:
#'   \item{estimates}{Matrix with one row per point in `pts`. First column is
#'     point number, remaining columns are estimated probabilities for each
#'     response category.}
#'   \item{pts}{Matrix of predictor values where probabilities were estimated.}
#'
#' @details
#' This function fits a multinomial logistic regression model:
#'
#' \deqn{P(Y=k|X) = \frac{\exp(\beta_{k0} + \beta_k' X)}{1 + \sum_j \exp(\beta_{j0} + \beta_j' X)}}
#'
#' for k = 1, ..., K-1, where K is the number of response categories. The
#' baseline category (k=0) has P(Y=0|X) = 1 / (1 + sum of other probabilities).
#'
#' The function:
#' \enumerate{
#'   \item Fits multinomial logit model using `multinom()` from nnet package
#'   \item Computes predicted probabilities at each point in `pts`
#'   \item Optionally plots probabilities vs predictors
#' }
#'
#' Plotting behavior:
#' \itemize{
#'   \item 1 predictor: Line plot(s) of probability vs X
#'   \item 2 predictors: 3D surface plot using `lplot()` or `rplot()`
#'   \item >2 predictors: No automatic plot
#' }
#'
#' If `xout = TRUE`, outliers in predictor space are removed before fitting.
#'
#' Requires the nnet package.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{logreg}}, \code{\link{logreg.pred}}, \code{\link[nnet]{multinom}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Multinomial regression with 3 categories
#' set.seed(123)
#' x <- rnorm(100)
#' # Create 3-category outcome
#' y <- cut(2*x + rnorm(100), breaks=3, labels=FALSE)
#' result <- multireg.prob(x, y)
#' result$estimates  # Probabilities for each category
#' }
multireg.prob<-function(x,y,pts=x,xout=FALSE,outfun=outpro,plotit=TRUE,xlab='X',ylab='Prob',zlab='Prob',ticktype='det',vplot=NULL,
L=TRUE,scale=TRUE,...){
#
#
# Returns estimate of P(Y=k|X=pts)
# for all possible values of k and all points stored in pts.
# using a multinomial logit model
#
#  Requires R package nnet
#
# scale =TRUE is the default:
# if  there is only p=1 independent variable, the y-axis of the plot of the regression line will range between 0 and 1.
# This can provide a useful perspective, particularly when there is no association.
#  if scale=TRUE, the y-axis is limited to the range of estimated probabilities.
#
library(nnet)
xy=cbind(x,y)
xy=elimna(xy)
p1=ncol(xy)
p=p1-1
x=xy[,1:p]
y=xy[,p1]
if(p==1){
pts=sort(pts)
}
x=as.matrix(x)
if(xout){
flag=outfun(x,plotit=FALSE,...)$keep
x=x[flag,]
y=y[flag]
}
pts=as.matrix(pts)
npts=nrow(pts)
est=summary(multinom(y~x))$coefficients
x=as.matrix(x)
nv=length(unique(y))
nvm1=nv-1
w=NA
pr=NA
if(is.null(dim(est)))est=matrix(est,nrow=1)
ans=matrix(NA,nrow=npts,ncol=nvm1)
for(k in 1:nrow(pts)){
for(j in 1:nvm1){
w[j]=exp(est[j,1]+sum(est[j,2:p1]*pts[k,]))
}
bot=1+sum(w)
ans[k,]=w/bot
}
v0=1-apply(ans,1,sum)
ptn=c(1:nrow(pts))
res=cbind(ptn,v0,ans)
temp=sort(unique(y))
dimnames(res)=list(NULL,c('pts.no',temp))
if(plotit){
if(is.null(vplot))vplot=max(y)
vplot=vplot+1  # adjustment to match col of res
if(p==1){
nlines=min(ncol(res),6)
nlines=nlines-1
if(scale)plot(c(pts[1:2],rep(pts,length(vplot))),c(0,1,as.vector(res[,vplot])),type='n',xlab=xlab,ylab=ylab)
if(!scale)plot(rep(pts,length(vplot)),as.vector(res[,vplot]),type='n',xlab=xlab,ylab=ylab)
for(k in 1:length(vplot))lines(pts,res[,vplot[k]],lty=k)
}
if(p>1){
if(ylab=='Prob')ylab='Y'
if(p==2){
if(L)lplot(pts,res[,vplot[1]],xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,scale=scale,pr=FALSE)
if(!L)rplot(pts,res[,vplot[1]],xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,scale=scale,pr=FALSE)
}
}
}
list(estimates=res,pts=pts)
}


# ============================================================================
# regIVbinv2_sub
# ============================================================================

#' Bootstrap Helper for Multinomial Logistic Regression
#'
#' Internal helper function that computes standard deviation of predicted
#' probabilities from \code{\link{logreg.pred}} for bootstrap samples.
#'
#' @param data Bootstrap sample indices.
#' @param x Matrix of predictors.
#' @param y Response vector.
#' @param pts Design points at which to predict.
#'
#' @return Standard deviation of predicted probabilities.
#'
#' @keywords internal
#' @export
regIVbinv2_sub<-function(data,x,y,pts){
x=as.matrix(x)
v=logreg.pred(x[data,],y[data],pts=pts)
v=sd(v)
}


# ============================================================================
# reg.hyp.split
# ============================================================================

#' Regression-Based Hypothesis Testing via Hyperplane Split
#'
#' Tests for association between independent variables and response by recursively
#' splitting the design space based on regression hyperplanes, creating a 2×2
#' factorial design for comparing measures of location across regions.
#'
#' @inheritParams common-params
#' @param x Matrix of independent variables (must have p ≥ 2 columns).
#' @param y Numeric vector containing the response variable.
#' @param split.reg Regression function used to define splitting hyperplanes
#'   (default: `Qreg` for median regression). Must return coefficients in `$coef`.
#' @param TR Trimming proportion for trimmed means if `PB=FALSE` (default: 0.2).
#' @param PB Logical; if TRUE, use percentile bootstrap instead of parametric
#'   method (default: FALSE).
#' @param est Measure of location if `PB=TRUE` (default: `tmean`).
#' @param nboot Number of bootstrap samples if `PB=TRUE` (default: 1000).
#' @param pr Logical; if TRUE, print detailed output (default: TRUE).
#' @param method Multiple comparison method: "hoch" (Hochberg) or other options
#'   passed to `lincon()` or `linconpb()` (default: "hoch").
#' @param xout Logical; if TRUE, remove outliers from predictor space before
#'   splitting (default: FALSE).
#' @param ... Additional arguments passed to `split.reg` or comparison functions.
#'
#' @return A list with components:
#'   \itemize{
#'     \item Independent.variables.summary: Summary statistics for each of the
#'       4 regions created by recursive splitting
#'     \item output: Results from `lincon()` or `linconpb()` with pairwise
#'       comparisons of the 4 regions
#'   }
#'
#' @details
#' This function implements a diagnostic technique to detect associations by:
#' \enumerate{
#'   \item Fitting a regression hyperplane (using `split.reg`) to the full data
#'   \item Splitting data into two groups based on residuals (above/below hyperplane)
#'   \item Recursively splitting each group again using the same method
#'   \item Comparing response locations across the resulting 4 regions
#' }
#'
#' If no association exists, the 4 regions should have similar response distributions.
#' Significant differences suggest association between predictors and response.
#'
#' **Split regression choices:**
#' \itemize{
#'   \item `Qreg` (default): Median regression split
#'   \item `Qreg, q=0.25`: Lower quartile regression split
#'   \item `depreg`: Deepest regression estimator split
#'   \item Any function returning `$coef`
#' }
#'
#' @note Requires at least 2 independent variables (p ≥ 2). For binary responses,
#'   use `regbin.hyp.split()` instead.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regbin.hyp.split}}, \code{\link{Qreg}}, \code{\link{lincon}}, \code{\link{linconpb}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Detect association via hyperplane splits
#' set.seed(123)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- 2*x1 + 3*x2 + rnorm(100)
#' reg.hyp.split(cbind(x1, x2), y)
#'
#' # Using bootstrap with median
#' reg.hyp.split(cbind(x1, x2), y, PB=TRUE, est=median, nboot=500)
#'
#' # Split based on lower quartile regression
#' reg.hyp.split(cbind(x1, x2), y, split.reg=Qreg, q=0.25)
#' }
reg.hyp.split<-function(x,y,split.reg=Qreg,TR=.2,alpha = 0.05, PB = FALSE, est = tmean, nboot = 1000, pr = TRUE,
    method = "hoch", xout = FALSE, outfun = outpro, SEED = TRUE, ...){
#
#   Split design space based on the hyperplane associated with the argument
#   split.reg
#   Default is a quantile regression estimate based on the data in x
#   Split the original data then split the results again to get a 2-by-2 ANOVA design
#
# Compare measures of location based on the resulting splits
#
#  Choices for split.reg: any R function that returns coefficients in $coef
#  Ex. split.reg=depreg would use a deepest regression estimator.
#   Could get different split using different quantiles
#   Ex.split.reg=Qreg,q=.25, would split the design space based  .25 quantile hyperplanes.
#   split.reg=mdepreg.coef  would use the deepest regression line estimator.
#
#
p=ncol(x)
xy=elimna(cbind(x,y))
if(xout){
flag<-outfun(xy[,1:p],plotit=FALSE,...)$keep
xy<-xy[flag,]
}
if(identical(est,median))PB=TRUE
if(identical(est,hd))PB=TRUE
if(p<2)stop('Should have two or more independent variables')
pm1=p-1
p1=p+1
hat=reg.pred(xy[,1:pm1],xy[,p],regfun=split.reg,...)
res=xy[,p]-hat
flag=res>0
x1=xy[flag,]
x2=xy[!flag,]
#
hat=reg.pred(x1[,1:pm1],x1[,p],regfun=split.reg,...)
res=x1[,p]-hat
flag=res>0
xy1=x1[flag,]
xy2=x1[!flag,]
#
hat=reg.pred(x2[,1:pm1],x2[,p],regfun=split.reg,...)
res=x2[,p]-hat
flag=res>0
xy3=x2[flag,]
xy4=x2[!flag,]
y=list()
y[[1]]=xy1[,p1]
y[[2]]=xy2[,p1]
y[[3]]=xy3[,p1]
y[[4]]=xy4[,p1]
group=list()
group[[1]]=summary(xy1[,1:p])
group[[2]]=summary(xy2[,1:p])
group[[3]]=summary(xy3[,1:p])
group[[4]]=summary(xy4[,1:p])
if(!PB)a=lincon(y,tr=TR)
if(PB)a=linconpb(y,est=est,nboot=nboot,...)
list(Independent.variables.summary=group,output=a)
}


# ============================================================================
# regbin.hyp.split
# ============================================================================

#' Binary Response Hypothesis Testing via Hyperplane Split
#'
#' Tests for association between independent variables and a binary response by
#' recursively splitting the design space based on regression hyperplanes, then
#' comparing proportions across the resulting regions using methods appropriate
#' for binary data.
#'
#' @inheritParams reg.hyp.split
#' @param y Binary response vector (0/1 or logical).
#' @param method Comparison method for binary data: "SK" (Storer-Kim), "BH"
#'   (Benjamini-Hochberg), or "BY" (Benjamini-Yekutieli) (default: "SK").
#' @param ... Additional arguments passed to `split.reg`.
#'
#' @return A list with components:
#'   \itemize{
#'     \item Independent.variables.summary: Summary statistics for predictors
#'       in each of the 4 regions
#'     \item output: Results from `binmcp()` with pairwise comparisons of
#'       proportions across the 4 regions
#'   }
#'
#' @details
#' Similar to `reg.hyp.split()` but designed specifically for binary responses.
#' The recursive splitting procedure:
#' \enumerate{
#'   \item Fits regression hyperplane to binary response data
#'   \item Splits data based on residuals (above/below hyperplane)
#'   \item Recursively splits each group again
#'   \item Compares success proportions across 4 regions using `binmcp()`
#' }
#'
#' Under the null hypothesis of no association, the 4 regions should have
#' similar success probabilities. Significant differences indicate association.
#'
#' The Storer-Kim (SK) method controls family-wise error rate for binary
#' comparisons and is appropriate when comparing multiple proportions.
#'
#' @note Requires at least 2 independent variables (p ≥ 2). For continuous
#'   responses, use `reg.hyp.split()` instead.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' Storer, B.E., & Kim, C. (1990). Exact properties of some exact test statistics
#' for comparing two binomial proportions. Journal of the American Statistical
#' Association, 85, 146-155.
#'
#' @seealso \code{\link{reg.hyp.split}}, \code{\link{binmcp}}, \code{\link{Qreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test association with binary response
#' set.seed(123)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' # Binary response influenced by predictors
#' p <- plogis(x1 + 2*x2)
#' y <- rbinom(100, 1, p)
#' regbin.hyp.split(cbind(x1, x2), y)
#'
#' # Using Benjamini-Hochberg FDR control
#' regbin.hyp.split(cbind(x1, x2), y, method="BH", nboot=500)
#' }
regbin.hyp.split<-function(x,y,split.reg=Qreg,alpha = 0.05, nboot = 1000,
    method ='SK', xout = FALSE, outfun = outpro, SEED = TRUE, ...){
#
#  y is assumed to be binary
#
#   Split design space based on the hyperplane associated with the argument
#   split.reg
#   Default is a squantile regression estimate based on the data in x
#   Split the original data then split the results again to get a 2-by-2 ANOVA design
#
# Compare binomial distributions using based on the argument
#  method, which defaults to Storer--Kim. To get confidence intervals use
#  method='KMS'
#
#  Choices for split.reg: any R function that returns coefficients in $coef
#  Ex. split.reg=depreg would use a deepest regression estimator.
#   Could get different split using different quantiles
#   Ex.split.reg=Qreg,q=.25, would split the design space based  .25 quantile hyperplanes.
#   split.reg=mdepreg.coef  would use the deepest regression line estimator.
#
#
p=ncol(x)
xy=elimna(cbind(x,y))
if(xout){
flag<-outfun(xy[,1:p],plotit=FALSE,...)$keep
xy<-xy[flag,]
}
if(length(unique(y))>2)stop('y should be binary')
n=length(y)
yy=rep(0,n)
flag=which(y==max(y))
yy[flag]=1
y=yy
if(p<2)stop('Should have two or more independent variables')
pm1=p-1
p1=p+1
hat=reg.pred(xy[,1:pm1],xy[,p],regfun=split.reg,...)
res=xy[,p]-hat
flag=res>0
x1=xy[flag,]
x2=xy[!flag,]
#
hat=reg.pred(x1[,1:pm1],x1[,p],regfun=split.reg,...)
res=x1[,p]-hat
flag=res>0
xy1=x1[flag,]
xy2=x1[!flag,]
#
hat=reg.pred(x2[,1:pm1],x2[,p],regfun=split.reg,...)
res=x2[,p]-hat
flag=res>0
xy3=x2[flag,]
xy4=x2[!flag,]
r=NA
r[1]=sum(xy1[,p1])
r[2]=sum(xy2[,p1])
r[3]=sum(xy3[,p1])
r[4]=sum(xy4[,p1])
n=NA
n[1]=nrow(xy1)
n[2]=nrow(xy2)
n[3]=nrow(xy3)
n[4]=nrow(xy4)
group=list()
group[[1]]=summary(xy1[,1:p])
group[[2]]=summary(xy2[,1:p])
group[[3]]=summary(xy3[,1:p])
group[[4]]=summary(xy4[,1:p])
a=binpair(r,n,method=method,alpha=alpha)
list(Independent.variables.summary=group,output=a)
}


# ============================================================================
# quantregForest
# ============================================================================

#' Quantile Regression Forests
#'
#' Performs robust random forest regression using quantile regression forests,
#' which can estimate conditional quantiles and provide robust predictions
#' that are less sensitive to outliers than standard random forests.
#'
#' @param x A numeric matrix or data frame containing predictor variables.
#'   Cannot contain missing values.
#' @param y A numeric vector containing the response variable. Cannot contain
#'   missing values.
#' @param nthreads Number of threads to use for parallel processing (default: 1).
#' @param keep.inbag Logical. If `TRUE`, keeps in-bag indicators for each tree
#'   (default: `FALSE`). Setting to `TRUE` increases memory usage.
#' @param ... Additional arguments passed to the underlying random forest algorithm.
#'
#' @return A quantile regression forest object containing:
#'   \itemize{
#'     \item Forest structure and parameters
#'     \item Out-of-bag predictions
#'     \item Information for quantile prediction at new points
#'   }
#'   Use with `predict()` to obtain conditional quantile predictions.
#'
#' @details
#' Quantile regression forests extend random forests to estimate conditional
#' quantiles rather than just conditional means:
#' \itemize{
#'   \item Standard random forests estimate E(Y|X) using the mean
#'   \item Quantile regression forests can estimate any conditional quantile Q(tau|X)
#'   \item This provides more complete information about the conditional distribution
#'   \item Particularly useful for heteroscedastic data or when quantile estimates
#'     are of primary interest
#' }
#'
#' The implementation is based on:
#' \itemize{
#'   \item Modified code from GitHub maintained by Loris Michel
#'   \item Uses the methodology from Meinshausen (2006)
#' }
#'
#' Categorical predictors with more than 53 categories are not supported and
#' will trigger an error (limitation inherited from randomForest package).
#'
#' @note Missing values (NA) are not permitted in either predictors or response.
#'   Remove or impute missing values before analysis.
#'
#' @references
#' Meinshausen, N. (2006). Quantile Regression Forests. Journal of Machine Learning
#' Research, 7, 983-999.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regR.Forest}}, \code{\link{KNNreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Quantile regression forest
#' set.seed(123)
#' x <- matrix(rnorm(200), ncol=2)
#' y <- 2*x[,1] - x[,2] + rnorm(100)
#' qrf <- quantregForest(x, y)
#'
#' # Predict median and other quantiles at new points
#' newx <- matrix(rnorm(20), ncol=2)
#' # predict(qrf, newx, what=0.5)  # Median
#' # predict(qrf, newx, what=c(0.1, 0.9))  # 10th and 90th percentiles
#' }
quantregForest <-function(x,y, nthreads = 1, keep.inbag=FALSE, ...){
#
#    This function does robust random Forest regression based on
#    Nicolai Meinshausen (2006) Quantile Regression Forests,
#   Journal of Machine Learning Research, 7, 983-999.
#    The code used here is based on a modification code downloaded from github, which is maintained by
#    Loris Michel <michel@stat.math.ethz.ch>
#
#

x=as.data.frame(x)
  if(is.null(nrow(x)) || is.null(ncol(x)))
    stop(' x contains no data ')
  if( nrow(x) != length(y) )
    stop(' predictor variables and response variable must contain the same number of samples ')

  if (any(is.na(x))) stop('NA not permitted in predictors')
  if (any(is.na(y))) stop('NA not permitted in response')
  ## Check for categorial predictors with too many categories (copied from randomForest package)
   if (is.data.frame(x)) {
        ncat <- sapply(x, function(x) if(is.factor(x) && !is.ordered(x))
                       length(levels(x)) else 1)
      } else {
        ncat <- 1
    }
    maxcat <- max(ncat)
    if (maxcat > 32)
        stop('Can not handle categorical predictors with more than 32 categories.')
  ## Note that crucial parts of the computation
  ## are only invoked by the predict method
  cl <- match.call()
  cl[[1]] <- as.name('quantregForest')
  qrf <- if(nthreads > 1){
    parallelRandomForest(x=x, y=y, nthreads = nthreads,keep.inbag=keep.inbag, ...)
  }else{
    randomForest( x=x,y=y ,keep.inbag=keep.inbag,...)
  }
  nodesX <- attr(predict(qrf,x,nodes=TRUE),'nodes')
  rownames(nodesX) <- NULL
  nnodes <- max(nodesX)
  ntree <- ncol(nodesX)
  n <- nrow(x)
  valuesNodes  <- matrix(nrow=nnodes,ncol=ntree)
  for (tree in 1:ntree){
      shuffledNodes <- nodesX[rank(ind <- sample(1:n,n)),tree]
      useNodes <- sort(unique(as.numeric(shuffledNodes)))
      valuesNodes[useNodes,tree] <- y[ind[match(useNodes,shuffledNodes )]]
  }

  qrf[['call']] <- cl
  qrf[['valuesNodes']] <- valuesNodes
  if(keep.inbag){
  #
    # create a prediction vector with same shape as predictOOBNodes
    predictOOBNodes <- attr(predict(qrf,newdata=x,nodes=TRUE),'nodes')
    rownames(predictOOBNodes) <- NULL
    valuesPredict <- 0*predictOOBNodes
    ntree <- ncol(valuesNodes)
    valuesPredict[qrf$inbag >0] <- NA
    #
    # for each tree and observation sample another observation of the same node
    for (tree in 1:ntree){
      is.oob <- qrf$inbag[,tree] == 0
      n.oob <- sum(is.oob)
      if(n.oob!=0) {
	  y.oob  <- sapply(which(is.oob),
		    function(i) {
			    cur.node <- nodesX[i, tree]
			    y.sampled <- if (length(cur.y <- y[setdiff(which(nodesX[,tree] == cur.node)
			                                               ,i)])!=0) {
			                cur.y[sample(x = 1:length(cur.y), size = 1)]
			                 } else {
			              	   NA
			    	           }
			    return(y.sampled)
		       })
          valuesPredict[is.oob, tree] <- y.oob
      }
    }

    minoob <- min( apply(!is.na(valuesPredict),1,sum))
    if(minoob<10) stop('need to increase number of trees for sufficiently many out-of-bag observations')
    valuesOOB <- t(apply( valuesPredict,1 , function(x) sample( x[!is.na(x)], minoob)))
    qrf[['valuesOOB']] <- valuesOOB
  }
  class(qrf) <- c('quantregForest','randomForest')

  return(qrf)
}


# ============================================================================
# regR.Forest
# ============================================================================
#' Random Forest Regression with Robust Location Estimation
#'
#' Estimates a measure of location (e.g., trimmed mean, median) at predictor
#' values using quantile random forests, with optional LOESS smoothing for
#' visualization and improved predictions.
#'
#' @inheritParams common-params
#' @param x A numeric matrix or vector containing predictor values (training data).
#' @param y A numeric vector containing response values (training data).
#' @param newdata A data frame of predictor values at which predictions are
#'   desired. If `NULL` (default), predictions are made at training data points.
#' @param pts Predictor values for smoothing (default: `x`).
#' @param pyhat Logical; if `TRUE`, return predicted values. If `FALSE` (default),
#'   return `NULL`.
#' @param loc.fun Location function to estimate from the random forest quantile
#'   predictions (default: `tmean` for 20% trimmed mean). Can be any location
#'   estimator like `median`, `mean`, `mom`, etc.
#' @param plotit Logical; if `TRUE` (default), create a plot of the fitted
#'   regression surface or curve.
#' @param span Span parameter for LOESS smoothing (default: 0.75). Controls the
#'   degree of smoothing.
#' @param LP Logical; if `TRUE` (default), apply LOESS smoothing to random forest
#'   predictions. If `FALSE`, return raw random forest predictions.
#' @param pch Plotting character for scatter plot (default: '.').
#' @param ZLIM Logical; controls axis limits for 3D plots.
#' @param scale Logical; if `TRUE` (default), scale axes in plots.
#' @param xlab,ylab,zlab Axis labels for plots.
#' @param ticktype Tick mark type for 3D plots (default: 'simple').
#' @param frame Logical; if `TRUE` (default), draw frame around 3D plots.
#' @param eout Logical; if `TRUE`, remove outliers from fitted values in plots.
#' @param theta,phi Viewing angles for 3D plots (defaults: 50, 25).
#'
#' @return If `pyhat = TRUE`, returns a vector of predicted values (smoothed or
#'   raw depending on `LP`). Otherwise, returns `NULL` (produces plot only).
#'
#' @details
#' This function combines random forest regression with robust location estimation:
#' \enumerate{
#'   \item Fits a quantile random forest to (x, y) using `quantregForest()`
#'   \item Predicts conditional quantiles at `newdata` points
#'   \item Estimates a robust location measure from the quantile predictions
#'   \item Optionally smooths predictions using LOESS for better visualization
#' }
#'
#' For single predictors, creates a 2D scatter plot with fitted curve. For two
#' predictors, creates a 3D surface plot.
#'
#' If `xout = TRUE`, outliers are removed from the predictor space before
#' analysis using the function specified by `outfun`.
#'
#' @note This function requires the `randomForest` package.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{KNNreg}}, \code{\link{qhdsm}}, \code{\link[randomForest]{randomForest}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Random forest regression with 1 predictor
#' set.seed(123)
#' x <- runif(100, 0, 10)
#' y <- sin(x) + rnorm(100, sd=0.3)
#' regR.Forest(x, y, pyhat=TRUE)
#'
#' # With 2 predictors (creates 3D plot)
#' x2 <- matrix(runif(200), ncol=2)
#' y2 <- x2[,1]^2 + x2[,2] + rnorm(100)
#' regR.Forest(x2, y2)
#' }
regR.Forest<-function(x,y,newdata=NULL,pts=x,pyhat=FALSE,loc.fun=tmean,xout=FALSE,plotit=TRUE,outfun=outpro, span = 0.75,LP=TRUE,pch='.',
ZLIM = FALSE, scale = TRUE, xlab = 'X', ylab = 'Y', ticktype='simple',frame=TRUE,eout=FALSE,
    zlab ='', theta = 50, phi = 25,...){

#  Goal: estimate a measure of location for newdata based
#  on the Random Forest method
#   Default, estimate  measure of location for training data x
#
#  loc.fun: a function indicating the measure of location to be estimated.
#  Default is a 20% trimmed mean
#
#  Method, initially use random forest then smooth using LOESS
#  pyhat=TRUE: return the predicted values
#   if LP=FALSE, return the random forest predicted values instead.
#
x<-as.matrix(x)
p=ncol(x)
p1=p+1
if(p==1){
xs=order(x)
x=x[xs]
y=y[xs]
}
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:p]
x<-as.matrix(x)
y<-xx[,p1]
x<-as.data.frame(x)
if(xout){
x<-as.data.frame(x)
flag<-outfun(x,plotit=FALSE)$keep
x<-x[flag,]
y<-y[flag]
x<-as.data.frame(x)
}
if(is.null(newdata))newdata=as.data.frame(x)
library(randomForest)
a=quantregForest(x,y)
res=predict.robust.Forest(a,newdata=newdata,what=loc.fun,...)
if(plotit){
if(p==2)
lplot(x,res,ZLIM=ZLIM,span=span,scale=scale,xlab=xlab,ylab=ylab,zlab,zlab,ticktype=ticktype,frame=frame,theta=theta,
phi=phi,pyhat=pyhat,eout=eout,xout=FALSE,pr=FALSE)
if(p==1){
plot(x[,1],y,xlab=xlab,ylab=ylab,pch=pch)
xs=order(x[,1])
e=lplot.pred(x[xs,1],res[xs],span=span)$yhat
lines(x[,1],e)
}
}
if(LP)res=lplot.pred(x,res,pts=pts)
if(!pyhat)res=NULL
res
}


# ============================================================================
# KNNreg
# ============================================================================

#' K-Nearest Neighbors Regression with Robust Estimation
#'
#' Performs K-nearest neighbors regression using robust location estimators
#' instead of the conventional mean. Predictions are based on a robust estimator
#' of the K nearest neighbors' Y values.
#'
#' @param x A numeric matrix containing predictor values (must have 2 or more columns).
#' @param y A numeric vector containing the response variable.
#' @param pts A matrix of points at which predictions are desired. If `NULL`,
#'   predictions are made at the observed x values (default: `NULL`).
#' @inheritParams common-params
#' @param K Number of nearest neighbors to use (default: 10).
#' @param est Estimator function for combining neighbors' Y values (default: `tmean`).
#'   Can be any location estimator like `median`, `mean`, `mom`, etc.
#' @param cov.fun Covariance function used for computing Mahalanobis distances
#'   to find nearest neighbors (default: `covmcd`).
#'
#' @return A vector of predicted Y values corresponding to rows of `pts` (or `x`
#'   if `pts` is `NULL`).
#'
#' @details
#' K-nearest neighbors regression predicts the response at a point by finding the
#' K closest points in the predictor space and combining their Y values. This
#' function uses several robust enhancements:
#' \itemize{
#'   \item Mahalanobis distance (via robust covariance) instead of Euclidean
#'   \item Robust location estimators (e.g., trimmed mean) instead of mean
#'   \item Optional outlier removal from predictor space (`xout = TRUE`)
#' }
#'
#' The function requires at least 2 predictors. For single-predictor problems,
#' use other smoothing methods like `qhdsm()` or `rplot()`.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{regR.Forest}}, \code{\link{qhdsm}}, \code{\link{rplot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # KNN regression with 2 predictors
#' set.seed(123)
#' x <- matrix(rnorm(200), ncol=2)
#' y <- x[,1]^2 + x[,2] + rnorm(100)
#' # Predict at observed points using 10 nearest neighbors
#' yhat <- KNNreg(x, y, K=10)
#' }
KNNreg<-function(x,y,pts=NULL,K=10,est=tmean,cov.fun=covmcd,
xout=FALSE,outfun=outpro,...){
x<-as.matrix(x)
if(ncol(x)==1)stop('Should have two or more independent variables')
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
if(is.null(pts))pts=x
n=1
if(is.matrix(pts) || is.data.frame(pts))
n=nrow(pts)
if(n==1)pts=matrix(pts,nrow=1)
e=NA
mcov=cov.fun(x)$cov
for(i in 1:n){
e[i]=est(y[nearNN(x,pt=pts[i,],K=K,mcov=mcov,...)])
}
e
}


# ============================================================================
# reg.resid
# ============================================================================

#' Compute Regression Residuals
#'
#' Computes residuals from any regression estimator. Supports robust regression
#' methods and optional outlier removal.
#'
#' @inheritParams common-params
#' @param regfun Regression function to use (default: `tsreg`). Can be any
#'   function that works with `reg.pred()`.
#' @param xout Logical; if TRUE, remove outliers from predictor space before
#'   regression (default: FALSE).
#' @param ... Additional arguments passed to `regfun`.
#'
#' @return A numeric vector of residuals (observed Y - predicted Y).
#'
#' @details
#' This is a utility function that:
#' \itemize{
#'   \item Fits the specified regression model
#'   \item Computes predicted values at observed x-values
#'   \item Returns y - yhat (residuals)
#' }
#'
#' Works with any regression function compatible with `reg.pred()`, including:
#' \itemize{
#'   \item `tsreg` (Theil-Sen)
#'   \item `opreg` (outlier-pruned)
#'   \item `ltsreg` (LTS)
#'   \item `Qreg` (quantile regression)
#' }
#'
#' Primarily used for diagnostics, residual plots, and identifying outliers.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{reg.pred}}, \code{\link{tsreg}}, \code{\link{regcon.out}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Residuals from Theil-Sen regression
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' res <- reg.resid(x, y, regfun=tsreg)
#' plot(x, res)
#' abline(h=0)
#'
#' # Residuals from quantile regression (median)
#' res.q <- reg.resid(x, y, regfun=Qreg, q=0.5)
#' }
reg.resid<-function(x,y,regfun=tsreg,xout=FALSE,outfun=outpro,...){
#
#  Compute residuals using any regression estimator.
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
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
e=reg.pred(x,y,regfun=regfun)
res=as.vector(y-e)
res
}


# ============================================================================
# regIQR
# ============================================================================

#' Interquartile Range for Regression
#'
#' Computes the predicted interquartile range (IQR) of the conditional distribution
#' of Y given X at specified design points using quantile regression.
#'
#' @inheritParams common-params
#' @param xr Points at which to compute IQR (default: observed x values).
#'   Can be a vector (single predictor) or matrix (multiple predictors).
#' @param regfun Regression function to use (default: `Qreg` for quantile
#'   regression).
#' @param ... Additional arguments passed to `regfun`.
#'
#' @return A numeric vector of predicted IQR values at each point in `xr`.
#'
#' @details
#' Computes the difference between the 75th and 25th percentile regression lines:
#' \deqn{IQR(x) = \hat{Y}_{0.75}(x) - \hat{Y}_{0.25}(x)}
#'
#' This provides an estimate of the conditional spread (variability) of Y given
#' X, which can vary across the range of X (heteroscedasticity). Useful for:
#' \itemize{
#'   \item Detecting heteroscedasticity (non-constant variance)
#'   \item Understanding how variability changes with X
#'   \item Robust alternative to conditional variance estimation
#' }
#'
#' Unlike traditional variance estimates, the IQR is robust to outliers and
#' heavy-tailed distributions.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{Qreg}}, \code{\link{qhdsm}}, \code{\link{regYhat}}
#'
#' @export
#' @examples
#' \dontrun{
#' # IQR across range of x (detect heteroscedasticity)
#' set.seed(123)
#' x <- runif(100, 0, 10)
#' y <- 2*x + x*rnorm(100)  # Variance increases with x
#' iqr.vals <- regIQR(x, y)
#' plot(x, iqr.vals, xlab="X", ylab="Conditional IQR")
#'
#' # IQR at specific points
#' regIQR(x, y, xr=c(2, 5, 8))
#' }
regIQR<-function(x,y,xr=x,regfun=Qreg,xout=FALSE,outfun=outpro,...){
#
#
IQR=regYhat(x,y,xr=xr,regfun=regfun,q=.75)-regYhat(x,y,xr=xr,regfun=regfun,q=.25)
IQR
}


# ============================================================================
# qinvreg
# ============================================================================

#' Inverse Quantile Regression
#'
#' Finds the quantile level q such that the quantile regression prediction at
#' a given point equals a specified value. Useful for determining which quantile
#' of Y given X corresponds to a particular target value.
#'
#' @inheritParams common-params
#' @param pt A single design point (scalar for p=1, vector for p>1) at which
#'   to evaluate the quantile regression.
#' @param v Target predicted value to achieve.
#' @param REQMIN Convergence criterion: stop when |prediction - v| < REQMIN
#'   (default: 0.001).
#'
#' @return The quantile level q ∈ (0, 1) such that the q-quantile regression
#'   prediction at `pt` approximately equals `v`.
#'
#' @details
#' This function inverts the quantile regression relationship: instead of
#' predicting Y for a given X and quantile q, it finds the quantile q that
#' produces a target prediction v at a given X.
#'
#' The algorithm performs a binary search over q ∈ (0, 1) to find the quantile
#' level that yields a prediction closest to the target value v.
#'
#' **Applications:**
#' \itemize{
#'   \item Estimating conditional cumulative distribution functions
#'   \item Finding percentile ranks for specific Y values given X
#'   \item Inverse prediction problems in quantile regression
#' }
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{Qreg}}, \code{\link{reg.con.dist}}, \code{\link{qhdsm}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Find quantile level where prediction equals 5
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' # What quantile at x=0 has predicted value 5?
#' q <- qinvreg(x, y, pt=0, v=5)
#' cat("Quantile level:", q, "\n")
#'
#' # Verify: check that Qreg at this quantile predicts v
#' pred <- reg.pred(x, y, xr=0, regfun=Qreg, q=q)
#' cat("Prediction:", pred, "Target:", 5, "\n")
#' }
qinvreg<-function(x,y,pt,v,REQMIN=.001){
#
#  Find q such that for Qreg  Y hat equals v
#
xy=cbind(x,y)
a=nelderv2(xy,1,qinvreg.sub,START=.5,pt=pt,v=v,REQMIN=REQMIN)
#   note: using  optim, even with BFGS method, can result in highly inaccurate values
a
}


# ============================================================================
# qinvreg.sub
# ============================================================================

#' Objective Function for qinvreg Optimization
#'
#' Internal helper function used by \code{\link{qinvreg}} to find the quantile
#' level that produces a target predicted value at a specified design point.
#'
#' @param xy Bootstrap sample (columns: x, y).
#' @param q Quantile level to evaluate.
#' @param pt Design point at which to predict.
#' @param v Target predicted value.
#'
#' @return Absolute difference between predicted value and target (to minimize).
#'
#' @keywords internal
#' @export
qinvreg.sub<-function(xy,q,pt,v){
e=reg.pred(xy[,1],xy[,2],xr=pt,regfun=Qreg,q=q,xout=FALSE)
a=abs(e-v)
a
}


# ============================================================================
# reg.con.dist
# ============================================================================

#' Estimate Conditional Distribution via Quantile Regression
#'
#' Estimates the conditional distribution of Y given X=pts by fitting quantile
#' regressions across a grid of quantile levels. Provides a nonparametric estimate
#' of the conditional CDF at specified design points.
#'
#' @inheritParams common-params
#' @param pts Design point(s) at which to estimate the conditional distribution.
#'   If `NULL`, uses marginal medians of predictors (default: NULL).
#'   For single predictor: scalar. For multiple predictors: vector of length p.
#' @param plotit Logical; if TRUE, plot the estimated conditional density via
#'   kernel density estimation (default: FALSE).
#' @param xlab Label for x-axis if plotting (default: '').
#' @param ylab Label for y-axis if plotting (default: '').
#'
#' @return A numeric vector of length 99 containing predicted Y values
#'   corresponding to quantiles 0.01, 0.02, ..., 0.99 of the conditional
#'   distribution of Y|X=pts.
#'
#' @details
#' This function estimates the conditional distribution by:
#' \enumerate{
#'   \item Fitting quantile regressions at q = 0.01, 0.02, ..., 0.99
#'   \item Evaluating each quantile regression at the specified point(s)
#'   \item Returning the 99 predicted quantiles
#' }
#'
#' The result represents an empirical estimate of the conditional CDF evaluated
#' at percentiles 1% through 99%. This allows visualization and analysis of:
#' \itemize{
#'   \item Conditional distribution shape (skewness, multimodality)
#'   \item Conditional variability
#'   \item Departure from normality given X
#' }
#'
#' If `pts=NULL`, uses the marginal median of each predictor as the evaluation
#' point (i.e., estimates distribution at the "center" of X-space).
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{Qreg}}, \code{\link{qinvreg}}, \code{\link{regIQR}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Estimate conditional distribution at x=0
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' cond.dist <- reg.con.dist(x, y, pts=0, plotit=TRUE)
#'
#' # At median of x
#' cond.dist.med <- reg.con.dist(x, y, pts=NULL, plotit=TRUE)
#'
#' # Examine conditional quantiles
#' summary(cond.dist)
#' }
reg.con.dist<-function(x,y,pts=NULL,plotit=FALSE,xlab='',ylab=''){
#
#  Estimate the conditional distribution of Y given x=pts
#  assuming a linear quantile regression model.
#
#  pts=NULL: the marginal medians of x are used
#
if(is.null(pts)){
x=as.matrix(x)
pts=apply(x,2,median)
}
iv=c(1:99)
v=NA
for(i in 1:99)v[i]=reg.pred(x,y,pts,Qreg,q=i/100)
if(plotit)akerd(v,xlab=xlab,ylab=ylab)
v
}


# ============================================================================
# regcon.out
# ============================================================================

#' Detect Outliers in Regression Allowing Heteroscedasticity
#'
#' Identifies outliers in the response variable given predictors using a method
#' that accommodates heteroscedasticity. Particularly effective at detecting
#' bad leverage points when error variance is non-constant.
#'
#' @inheritParams common-params
#' @param x A numeric vector or matrix containing predictor values.
#' @param y A numeric vector containing the response variable.
#' @param plotit Logical. If `TRUE`, creates a scatterplot highlighting
#'   outliers (default: `TRUE`).
#' @param xlab,ylab Axis labels for the plot.
#'
#' @return A list with components:
#'   \item{n}{Total number of observations.}
#'   \item{n.out}{Number of outliers detected.}
#'   \item{out.id}{Vector of indices identifying which observations are outliers.}
#'   \item{keep}{Vector of indices for non-outlier observations.}
#'
#' @details
#' This function detects outliers using a robust method that handles
#' heteroscedasticity:
#' \enumerate{
#'   \item Estimates conditional quantiles using quantile regression
#'   \item Computes a robust measure of dispersion based on quantile differences
#'   \item Identifies outliers based on their distance from the conditional median
#'     in units of robust conditional standard deviation
#' }
#'
#' Advantages over traditional outlier detection methods:
#' \itemize{
#'   \item Allows for heteroscedasticity (non-constant error variance)
#'   \item Better at detecting bad leverage points under heteroscedasticity
#'   \item Performs similarly to other methods under homoscedasticity
#' }
#'
#' When `plotit=TRUE`, creates a scatterplot with:
#' \itemize{
#'   \item Regular observations shown as points
#'   \item Outliers marked with a different symbol
#' }
#'
#' Missing values are automatically removed before analysis.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{reg.con.dist}}, \code{\link{outpro}}, \code{\link{Qreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Detect outliers with heteroscedastic errors
#' set.seed(123)
#' x <- seq(0, 10, length.out=100)
#' y <- 2*x + rnorm(100, sd=0.5*x)  # Heteroscedastic errors
#' # Add some outliers
#' y[c(10, 50, 90)] <- c(25, -10, 40)
#' result <- regcon.out(x, y)
#' result$out.id  # Indices of outliers
#' }
regcon.out<-function(x,y,plotit=TRUE,xlab='X',ylab='Y'){
#
# Detect outliers among y given x
# in manner that allows heteroscedasticity.
# This improves on other methods for detecting bad leverage points
# when there is heteroscedasticity.. Works about as well as other methods when
# there is homoscedasticity.
#
xx<-cbind(x,y)
if(ncol(xx)!=2)stop('Current version limited to a single independent variable')
xx<-elimna(xx)
x<-xx[,1]
y<-xx[,2]
temp<-NA
xord=order(x)
n=length(x)
vec=keep=c(1:n)
rem=outpro(x)
keepid=rem$keep
iout=rem$out.id
flag=rep(FALSE,n)
for(i in 1:n){
q2=Qreghat(x[keepid],y[keepid],xr=x[i],q=.75)
q1=Qreghat(x[keepid],y[keepid],xr=x[i],q=.25)
iq=q2-q1
top=q2+1.5*iq
bot=q1-1.5*iq
if(y[i]<bot || y[i]>top)flag[i]=TRUE
}
outid <- NULL
if(sum(flag) > 0)outid <- vec[flag] #regression outlier
both=c(iout,outid)
blp=duplicated(both)
if(sum(!blp)>0)
blp=unique(both[blp])
else
 blp=NULL
glp=iout
if(length(blp)>0){
flag=NULL
for(k in 1:length(blp)){
flag=c(flag,which(iout==blp[k]))
}
glp=iout[-flag]
keep=vec[-blp]
}
if(plotit){
plot(x,y,type='n',xlab=xlab,ylab=ylab)
points(x[keep],y[keep],pch='*')
points(x[blp],y[blp],pch='o')
}
list(n=n,n.out=length(iout),res.out.id=outid,keep=keep,good.lev=glp,bad.lev=blp)
}


# ============================================================================
# reghet.blp
# ============================================================================

#' Regression Removing Bad Leverage Points with Heteroscedastic Method
#'
#' Performs robust regression after identifying and removing bad leverage points
#' using a method designed for heteroscedastic errors. Improves resistance to
#' outliers when error variance is non-constant.
#'
#' @inheritParams common-params
#' @param regfun Regression function to use after removing outliers (default:
#'   `tsreg` for Theil-Sen regression).
#' @param HH Logical; if TRUE, use the HH (Hung-Huggins) method for detecting
#'   bad leverage points; if FALSE, use `regcon.out()` (default: TRUE).
#' @param ... Additional arguments passed to `regfun`.
#'
#' @return Result from `regfun` after removing bad leverage points. Typically
#'   a list with regression coefficients and other diagnostic information.
#'
#' @details
#' This function combines two steps:
#' \enumerate{
#'   \item Detect bad leverage points using a heteroscedastic-robust method
#'   \item Fit regression on data with bad leverage points removed
#' }
#'
#' **Bad leverage point detection methods:**
#' \itemize{
#'   \item HH=TRUE: Hung-Huggins method (via `outblp.HH()`)
#'   \item HH=FALSE: Quantile-based method (via `regcon.out()`)
#' }
#'
#' Both methods accommodate heteroscedasticity, making them superior to
#' traditional outlier detection when error variance changes with X.
#'
#' This approach is particularly useful when:
#' \itemize{
#'   \item Error variance increases/decreases with X
#'   \item Suspected bad leverage points exist
#'   \item Standard regression diagnostics fail due to heteroscedasticity
#' }
#'
#' @note Currently limited to single independent variable (p=1).
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{reghet.blp.ci}}, \code{\link{regcon.out}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Regression with heteroscedastic errors and outliers
#' set.seed(123)
#' x <- seq(0, 10, length.out=100)
#' y <- 2*x + rnorm(100, sd=0.5*x)  # Heteroscedastic
#' # Add bad leverage point
#' x[50] <- 12; y[50] <- -5
#' result <- reghet.blp(x, y)
#' result$coef  # Regression coefficients
#' }
reghet.blp<-function(x,y,regfun=tsreg,HH=TRUE,...){

# Eliminate bad leverage points using a heteroscedastic method
# Then estimate the parameters.
#
xx<-cbind(x,y)
xx<-elimna(xx)
if(ncol(xx)!=2)stop('Current version limited to a single independent variable')
x<-xx[,1]
y=xx[,2]
if(HH)id= outblp.HH(x,y,regfun=regfun,plotit=FALSE)$keep
else id=regcon.out(x,y,plotit=FALSE)$keep
e=regfun(x[id],y[id],...)
e
}


# ============================================================================
# reghet.blp.ci
# ============================================================================

#' Confidence Intervals After Removing Bad Leverage Points (Heteroscedastic)
#'
#' Computes bootstrap confidence intervals for regression parameters after
#' identifying and removing bad leverage points using a heteroscedastic-robust
#' method. Provides inference that accounts for both heteroscedasticity and
#' outlying observations.
#'
#' @inheritParams reghet.blp
#' @param nboot Number of bootstrap samples (default: 999).
#' @param BCA Logical; if TRUE, use bias-corrected and accelerated (BCa)
#'   bootstrap; if FALSE, use percentile bootstrap (default: FALSE).
#' @param pr Logical; if TRUE, print informational messages (default: TRUE).
#' @param ... Additional arguments passed to `regfun`.
#'
#' @return
#' If `BCA=FALSE`: Result from `regci()` containing confidence intervals for
#'   both intercept and slope.
#'
#' If `BCA=TRUE`: Result from `reg.bca()` containing confidence interval for
#'   slope only (intercept not provided with BCa method).
#'
#' @details
#' This function provides robust inference by:
#' \enumerate{
#'   \item Detecting bad leverage points using heteroscedastic-robust method
#'   \item Removing identified bad leverage points
#'   \item Computing bootstrap confidence intervals on cleaned data
#' }
#'
#' **Bootstrap methods:**
#' \itemize{
#'   \item Percentile bootstrap (BCA=FALSE): Simpler, provides CIs for all parameters
#'   \item BCa bootstrap (BCA=TRUE): More accurate for small samples, slope only
#' }
#'
#' **When to use BCa:**
#' \itemize{
#'   \item Small sample sizes (n < 50) - more accurate coverage
#'   \item Skewed sampling distributions
#'   \item Only need slope inference (intercept not returned)
#' }
#'
#' The function will print a warning when n < 50 and BCA=FALSE, suggesting
#' BCA=TRUE may be safer.
#'
#' @note Currently limited to single independent variable (p=1).
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' Efron, B., & Tibshirani, R.J. (1993). An Introduction to the Bootstrap.
#' Chapman & Hall/CRC.
#'
#' @seealso \code{\link{reghet.blp}}, \code{\link{regci}}, \code{\link{reg.bca}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Bootstrap CI with heteroscedastic bad leverage point removal
#' set.seed(123)
#' x <- seq(0, 10, length.out=100)
#' y <- 2*x + rnorm(100, sd=0.5*x)
#' x[50] <- 12; y[50] <- -5  # Bad leverage point
#'
#' # Percentile bootstrap (both parameters)
#' ci.perc <- reghet.blp.ci(x, y, nboot=500)
#'
#' # BCa bootstrap (slope only, better for small n)
#' ci.bca <- reghet.blp.ci(x, y, BCA=TRUE, nboot=500)
#' }
reghet.blp.ci<-function(x,y,regfun=tsreg,nboot=999,HH=TRUE,
SEED=TRUE,BCA=FALSE,pr=TRUE,...){

# Eliminate bad leverage points using a heteroscedastic method
# Then compute a confidence interval for the slope
#
#Use bias corrected accelerated bootstrap when BCA=TRUE,
# otherwise use a percentile bootstrap
#
xx<-cbind(x,y)
if(ncol(xx)!=2)stop('Current version limited to a single independent variable')
xx<-elimna(xx)
n=nrow(xx)
if(!BCA & n<50) #print('Might be safer to use BCA=TRUE')
if(BCA & pr)print('Note: when BCA=TRUE, only returns a confidence interval for the slope')
x<-xx[,1]
y=xx[,2]
if(HH)id= outblp.HH(x,y,regfun=regfun,plotit=FALSE)$keep
else id=regcon.out(x,y,plotit=FALSE)$keep
if(BCA)e=reg.bca(x[id],y[id],SEED=SEED,regfun=regfun)
else e=regci(x[id],y[id],SEED=SEED,regfun=regfun,nboot=nboot,...)
e
}


# ============================================================================
# regHH
# ============================================================================

#' Regression After Removing Bad Leverage Points (HH Method)
#'
#' Performs regression after removing bad leverage points identified by the
#' HH (Hadi-Hadi) method via \code{\link{outblp.HH}}. Provides robust regression
#' estimates resistant to high-leverage outliers.
#'
#' @inheritParams common-params
#' @param x A numeric vector containing predictor values (currently limited to
#'   single predictor).
#' @param y A numeric vector containing the response variable.
#' @param regfun Regression function to use after outlier removal (default: `tsreg`
#'   for Theil-Sen). Must return a list with `$coef` component.
#' @param SO Logical. If `TRUE`, returns only slope coefficients (convenient for
#'   bootstrap methods); if `FALSE`, returns full regression object (default: `FALSE`).
#' @param ... Additional arguments passed to `regfun` when `SO=FALSE`.
#'
#' @return A list with component:
#'   \item{coef}{Regression coefficients. If `SO=FALSE`, full output from `regfun`;
#'     if `SO=TRUE`, only slope coefficients.}
#'
#' @details
#' This function combines bad leverage point detection with robust regression:
#' \enumerate{
#'   \item Uses \code{outblp.HH()} to identify high-leverage points that are also
#'     outliers (bad leverage points)
#'   \item Removes identified bad leverage points from the data
#'   \item Fits regression using `regfun` on the cleaned data
#' }
#'
#' Bad leverage points are influential observations that can severely distort
#' regression estimates. The HH method identifies points with:
#' \itemize{
#'   \item High leverage (unusual predictor values)
#'   \item Large residuals (poor fit to the pattern)
#' }
#'
#' The `SO=TRUE` option is useful when embedding this function in bootstrap
#' procedures that need only slope estimates for efficiency.
#'
#' @note Currently restricted to single predictor (p=1) regression. Missing
#'   values are automatically removed.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{outblp.HH}}, \code{\link{reghet.blp}}, \code{\link{tsreg}},
#'   \code{\link{opreg}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Robust regression with HH outlier removal
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2 + 3*x + rnorm(100)
#' # Add bad leverage points
#' x[1:2] <- c(5, -5)
#' y[1:2] <- c(-10, 15)
#'
#' # Compare with and without HH outlier removal
#' fit_all <- tsreg(x, y)
#' fit_HH <- regHH(x, y)
#' fit_all$coef
#' fit_HH$coef
#' }
regHH<-function(x,y,regfun=tsreg,SO=FALSE,...){
#
#
# SO=TRUE, estimate slope only, convenient for some bootstrap methods
#
xy=elimna(cbind(x,y))
if(ncol(xy)!=2)stop('Current version limited to a single independent variable')
id= outblp.HH(xy[,1],xy[,2])$keep
if(!SO)e=regfun(xy[id,1],xy[id,2],...)
else e=regfun(xy[id,1],xy[id,2])$coef
list(coef=e)
}


# ============================================================================
# reg.break
# ============================================================================

#' Estimate Regression Breakpoint (Piecewise Linear Regression)
#'
#' Detects and estimates the breakpoint in a piecewise linear regression
#' relationship where the slope changes abruptly. Uses a robust analog of the
#' Jones and Molitoris (1984) method.
#'
#' @inheritParams common-params
#' @param x A numeric vector containing the predictor variable (currently limited
#'   to single predictor).
#' @param y A numeric vector containing the response variable.
#' @param int Numeric vector specifying the search interval for the breakpoint
#'   (candidate x-values to evaluate). If `NULL`, uses quantiles from `qv` to `1-qv`
#'   divided into `npts` evenly spaced points.
#' @param TEST Logical. Currently not used (reserved for future extension).
#' @param regfun Regression function to use for fitting segments (default: `tsreg`
#'   for Theil-Sen).
#' @param varfun Variance function used to select optimal breakpoint (default: `pbvar`
#'   for percentage bend variance). The breakpoint minimizing pooled variance is selected.
#' @param qv Quantile value for defining default search interval (default: 0.2).
#'   Search interval is [`qv`-quantile, `1-qv`-quantile] of x.
#' @param npts Number of candidate breakpoints to evaluate when `int=NULL` (default: 25).
#' @param ... Additional arguments passed to `outfun` when `xout=TRUE`.
#'
#' @return A list with components:
#'   \item{n}{Sample size after removing missing values (and outliers if `xout=TRUE`).}
#'   \item{break.est}{Estimated breakpoint location (x-value where slope changes).}
#'   \item{coef.est}{2×2 matrix of regression coefficients:
#'     \itemize{
#'       \item Row 1 (Lower): Intercept and slope for x ≤ breakpoint
#'       \item Row 2 (Upper): Intercept and slope for x > breakpoint
#'     }}
#'
#' @details
#' This function implements robust piecewise linear regression:
#' \enumerate{
#'   \item For each candidate breakpoint value:
#'     \itemize{
#'       \item Fits regression to points below breakpoint
#'       \item Adjusts points above breakpoint to connect smoothly
#'       \item Fits regression to adjusted upper points
#'       \item Computes pooled variance of all residuals
#'     }
#'   \item Selects breakpoint minimizing pooled variance
#'   \item Re-fits final model segments at optimal breakpoint
#' }
#'
#' The method is a robust version of Jones & Molitoris (1984), using:
#' \itemize{
#'   \item Robust regression (`regfun`) instead of OLS
#'   \item Robust variance (`varfun`) for model selection
#'   \item Optional outlier removal (`xout=TRUE`)
#' }
#'
#' **Search interval control**:
#' \itemize{
#'   \item Default (`int=NULL`): Search between `qv` and `1-qv` quantiles of x
#'   \item Custom: Provide vector of candidate x-values via `int`
#'   \item Adjust `qv` and `npts` for finer/coarser search
#' }
#'
#' @note Currently restricted to single predictor (p=1). The argument `qv` must
#'   be less than 0.5.
#'
#' @references
#' Jones, B.E. & Molitoris, B.A. (1984). A Statistical Method for Determining
#' the Breakpoint of Two Lines. Analytical Biochemistry, 141, 287-290.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing
#' (5th ed.). Academic Press.
#'
#' @seealso \code{\link{tsreg}}, \code{\link{pbvar}}, \code{\link{outpro}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Piecewise linear data with breakpoint at x=0
#' set.seed(123)
#' x <- seq(-2, 2, length.out=100)
#' y <- ifelse(x < 0, 2 - 3*x, 2 + x) + rnorm(100, sd=0.3)
#'
#' # Estimate breakpoint
#' result <- reg.break(x, y)
#' result$break.est     # Should be near 0
#' result$coef.est      # Slopes for each segment
#'
#' # Visualize
#' plot(x, y)
#' bp <- result$break.est
#' abline(a=result$coef.est[1,1], b=result$coef.est[1,2], col="red")  # Lower
#' segments(bp, result$coef.est[1,1] + result$coef.est[1,2]*bp,
#'          max(x), result$coef.est[2,1] + result$coef.est[2,2]*max(x),
#'          col="blue")  # Upper
#' abline(v=bp, lty=2)  # Breakpoint
#' }
reg.break<-function(x,y,int=NULL,xout=FALSE,TEST=FALSE,regfun=tsreg,outfun=outpro,varfun=pbvar,qv=.2,npts=25,...){
#
# Estimate the break point of a regression line, where the line bends.
# That is, where the  slope suddenly changes.
# Use a. robust analog of the method in
# A Statistical Method for Determining the Breakpoint of Two Lines
#  Jones and Molitoris
# Analytical Biochemistry 287-290 (1984)
#
# int = interval over the IV used to search for the breakpoint.
#
if(qv>.5)stop('argument qv should be less than .5')
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
if(p!=1)stop('Current version limited to a single independent variable')
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
x=as.vector(x)
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
if(is.null(int)){
low=qest(x,qv)
up=qest(x,1-qv)
int=seq(low,up,length.out=npts)
}
nrem=length(y)
if(!is.null(int)){
x0=int
iup=length(x0)

}
v=NA
for(i in 1:iup){
id=x<x0[i]
nid=sum(!id)
e1=regfun(x[id],y[id])
yy=y[!id]-e1$coef[2]*x0[i]-e1$coef[1]
xx=x[!id]-x0[i]
term1=rep(x0[i],nid)
res2=regfun(xx,yy)$residuals 
v[i]=varfun(c(e1$residuals,res2))
}
vor=(order(v))
sel=vor[1]
est=x0[sel]
id=x<=est
elow=regfun(x[id],y[id])$coef
eup=regfun(x[!id],y[!id])$coef
output=matrix(NA,2,2,)
output[1,]=elow
output[2,]=eup
dimnames(output)=list(c('Lower','Upper'),c('Intercept','Slope'))
list(n=length(y),break.est=est,coef.est=output)
}


