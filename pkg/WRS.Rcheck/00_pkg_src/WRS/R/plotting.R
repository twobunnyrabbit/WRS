# ============================================================================
# plotting.R - Visualization and Plotting Functions
# ============================================================================
#
# Comprehensive collection of plotting and visualization functions for
# robust statistical methods.
#
# Categories:
# - Regression plotting (rplot, regplot, etc.)
# - Lowess/loess smoothing plots
# - GAM plotting
# - ANCOVA plotting
# - Group comparison plots
# - Error bar and box plots
# - Functional data plots
# - Interaction plots
# - Logistic/longitudinal regression plots
# - Depth-based plots (bagplot)
# - Distribution/density plots
# - 3D plotting
# - Specialized/utility plots
#
# Total functions: 97 (98 - 1 duplicate linplot excluded)
# ============================================================================


# ============================================================================
# REGRESSION PLOTTING
# ============================================================================

#' Plot Two Regression Lines for Group Comparison
#'
#' Creates a scatter plot with fitted regression lines for two independent
#' groups. Useful for visualizing and comparing regression relationships
#' between groups.
#'
#' @param x1 Numeric vector of predictor values for group 1.
#' @param y1 Numeric vector of response values for group 1.
#' @param x2 Numeric vector of predictor values for group 2.
#' @param y2 Numeric vector of response values for group 2.
#' @param regfun Regression function to use for fitting (default: `tsreg`
#'   for Theil-Sen regression). Can be any regression function that returns
#'   a list with a `coef` component.
#' @inheritParams common-params
#'
#' @details
#' The function removes missing values automatically using `elimna()`.
#' Group 1 regression line is plotted as a solid line, group 2 as a dashed
#' line (lty=2). Points for group 1 use pch1 character, group 2 uses pch2.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect.
#'
#' @seealso \code{\link{tsreg}}, \code{\link{reg2g.p2plot}} for 3D version
#'
#' @export
#' @examples
#' # Compare regression lines for two groups
#' x1 <- rnorm(50)
#' y1 <- 2*x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- 3*x2 + rnorm(50)
#' reg2plot(x1, y1, x2, y2)
reg2plot<-function(x1,y1,x2,y2,regfun=tsreg,xlab="X",ylab="Y",xout=FALSE,outfun=outpro,pch1='.',pch2='+',...){
#
#  For convenience
#  plot two regression lines corresponding to two groups.
#
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
if(xout){
if(identical(outfun,outblp))flag=outblp(x1,y1,plotit=FALSE)$keep
else
flag<-outfun(x1,plotit=FALSE,...)$keep
x1=x1[flag]
y1=y1[flag]
if(identical(outfun,outblp))flag=outblp(x2,y2,plotit=FALSE)$keep
else
flag<-outfun(x2,plotit=FALSE,...)$keep
x2=x2[flag]
y2=y2[flag]
}
plot(c(x1,x2),c(y1,y2),type="n",xlab=xlab,ylab=ylab)
points(x1,y1,pch=pch1)
points(x2,y2,pch=pch2)
abline(regfun(x1,y1,...)$coef)
abline(regfun(x2,y2,...)$coef,lty=2)
}

#' 3D Plot of Two Regression Surfaces for Group Comparison
#'
#' Creates a 3D scatter plot with fitted regression planes for two independent
#' groups when there are two predictor variables. Useful for visualizing and
#' comparing multiple regression relationships between groups.
#'
#' @param x1 Matrix or data frame with 2 columns of predictor values for group 1.
#' @param y1 Numeric vector of response values for group 1.
#' @param x2 Matrix or data frame with 2 columns of predictor values for group 2.
#' @param y2 Numeric vector of response values for group 2.
#' @param regfun Regression function to use for fitting (default: `tsreg`
#'   for Theil-Sen regression). Must return a list with a `coef` component.
#' @param COLOR Logical. If `TRUE`, plots group 1 plane in blue and group 2
#'   in red (default: `TRUE`).
#' @param STAND Logical. If `TRUE`, uses standardized outlier detection
#'   (default: `TRUE`).
#' @param tick.marks Logical. If `TRUE`, adds tick marks to 3D plot (default: `TRUE`).
#' @param type Type of plot. "p" for points, "n" for regression planes only
#'   (default: "p").
#' @inheritParams common-params
#'
#' @details
#' Requires the `scatterplot3d` package. The function removes missing values
#' automatically using `elimna()`. First group regression plane is blue (if COLOR=TRUE),
#' second group is red.
#'
#' @return List with components:
#' \item{coef.group.1}{Regression coefficients for group 1}
#' \item{coef.group.2}{Regression coefficients for group 2}
#'
#' @seealso \code{\link{reg2plot}}, \code{\link{regp2plot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare 3D regression surfaces for two groups
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- 2*x1[,1] + 3*x1[,2] + rnorm(50)
#' x2 <- matrix(rnorm(100), ncol=2)
#' y2 <- 3*x2[,1] + 2*x2[,2] + rnorm(50)
#' reg2g.p2plot(x1, y1, x2, y2)
#' }
reg2g.p2plot<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,xlab="Var 1",ylab="Var 2",zlab="Var 3",regfun=tsreg,COLOR=TRUE,STAND=TRUE,
tick.marks=TRUE,type="p",pr=TRUE,...){
#
# Create a 3D plot of points and plot regression surface for two groups.
#
#  Assumes that the package scatterplot3d has been installed.
#  If not, use the command install.packages("scatterplot3d")
#  assuming you are connected to the web.
#
# The regression method used is specified with the argument
#  regfun.
#
#  type="p", points will be plotted. Use type="n" to get only regression planes plotted
#
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=2)stop("Argument x1 must be stored in a matrix with 2 columns.")
if(ncol(x2)!=2)stop("Argument x2 must be stored in a matrix with 2 columns.")
xy1<-elimna(cbind(x1,y1))
xy2<-elimna(cbind(x2,y2))
if(xout){
if(!STAND)flag1=outfun(xy1[,1:2],plotit=FALSE,...)$keep
if(STAND)flag1=outpro(xy1[,1:2],plotit=FALSE,STAND=TRUE,...)$keep
if(!STAND)flag2=outfun(xy2[,1:2],plotit=FALSE,...)$keep
if(STAND)flag2=outpro(xy2[,1:2],plotit=FALSE,STAND=TRUE,...)$keep
xy1=xy1[flag1,]
xy2=xy2[flag2,]
}
x1=xy1[,1:2]
x2=xy2[,1:2]
y1=xy1[,3]
y2=xy2[,3]
library(scatterplot3d)
temp<-scatterplot3d(rbind(xy1,xy2),xlab=xlab,ylab=ylab,zlab=zlab,tick.marks=tick.marks,type=type)
vals1<-regfun(x1,y1,...)$coef
vals2<-regfun(x2,y2,...)$coef
if(COLOR){
if(pr)print("First group is blue")
temp$plane(vals1,col="blue")
temp$plane(vals2,col="red")
}
if(!COLOR){
temp$plane(vals1)
temp$plane(vals2)
}
list(coef.group.1=vals1,coef.group.2=vals2)
}

#' 3D Plot of Regression Surface
#'
#' Creates a 3D scatter plot with a fitted regression plane for data with
#' two predictor variables. Useful for visualizing multiple regression
#' relationships.
#'
#' @param x Matrix or data frame with 2 columns of predictor values.
#' @param y Numeric vector of response values.
#' @param regfun Regression function to use for fitting (default: `tsreg`
#'   for Theil-Sen regression). Must return a list with a `coef` component.
#' @param COLOR Logical. If `TRUE`, plots regression plane in blue (default: `FALSE`).
#' @param tick.marks Logical. If `TRUE`, adds tick marks to 3D plot (default: `TRUE`).
#' @inheritParams common-params
#'
#' @details
#' Requires the `scatterplot3d` package. The function removes missing values
#' automatically using `elimna()`. The x argument must have exactly 2 columns.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect.
#'
#' @seealso \code{\link{reg2g.p2plot}}, \code{\link{regplot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Plot 3D regression surface
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2*x[,1] + 3*x[,2] + rnorm(50)
#' regp2plot(x, y)
#' }
regp2plot<-function(x,y,xout=FALSE,outfun=out,xlab="Var 1",ylab="Var 2",zlab="Var 3",regfun=tsreg,COLOR=FALSE,tick.marks=TRUE,...){
#
# Create a 3D plot of points and plot regression surface.
# based on the regression estimator indicated by
#  regfun
#
#  Assumes that the package scatterplot3d has been installed.
#  If not, use the command install.packages("scatterplot3d")
#  assuming you are connected to the web.
#
# The regression method used is specified with the argument
#  regfun.
#
#  Package scatterplot3d is required. To install it, use the command
#  install.packages("scatterplot3d")
#  while connected to the web
#
x=as.matrix(x)
if(ncol(x)!=2)stop("Argument x must be stored in a matrix with 2 columns.")
xy<-elimna(cbind(x,y))
if(xout){
flag=outfun(xy[,1:2])$keep
xy=xy[flag,]
}
x=xy[,1:2]
y=xy[,3]
library(scatterplot3d)
temp<-scatterplot3d(xy,xlab=xlab,ylab=ylab,zlab=zlab,tick.marks=tick.marks)
vals<-regfun(x,y,...)$coef
if(COLOR)temp$plane(vals,col="blue")
if(!COLOR)temp$plane(vals)
}

#' General Regression Plot (1D or 2D)
#'
#' Creates a scatter plot with fitted regression line (1 predictor) or
#' 3D regression surface (2 predictors). Automatically adapts to number
#' of predictors.
#'
#' @param x Numeric vector or matrix of predictor values (1 or 2 columns only).
#' @param y Numeric vector of response values.
#' @param regfun Regression function to use for fitting (default: `tsreg`).
#' @param theta Viewing angle (rotation around z-axis) for 3D plots (default: 50).
#' @param phi Viewing angle (colatitude) for 3D plots (default: 25).
#' @param ticktype Tick mark style for 3D plots (default: 'simple').
#' @inheritParams common-params
#'
#' @details
#' For 1 predictor: Creates 2D scatter plot with regression line.
#' For 2 predictors: Creates 3D scatter plot with regression surface using `rplot()`.
#' Maximum 2 predictors allowed.
#'
#' @return Invisibly returns `NULL` (1D) or output from `rplot()` (2D).
#'
#' @seealso \code{\link{regp2plot}}, \code{\link{rplot}}
#'
#' @export
#' @examples
#' # 1D regression plot
#' x <- rnorm(50)
#' y <- 2*x + rnorm(50)
#' regplot(x, y)
regplot<-function(x,y,regfun=tsreg,xlab='X',ylab='Y',zlab='Z',
xout=FALSE,outfun=out,theta=50, phi=25,ticktype='simple',...){
x=as.matrix(x)
xy=elimna(cbind(x,y))
if(ncol(x)>2)stop('One or two predictors only is allowed,')
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
if(xout){
xy=cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag=outfun(x)$keep
x=xy[flag,1:p]
y=xy[flag,p1]
}
if(p==1){
plot(x,y,xlab=xlab,ylab=ylab)
abline(regfun(x,y,...)$coef)
}
if(p==2){
pyhat=regYhat(x,y,regfun=regfun,...)
temp=rplot(x,pyhat,scat=FALSE,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,pr=FALSE)
}
}

#' Regression Interaction Plot
#'
#' Investigates regression interaction by fitting a generalized additive model
#' and plotting residuals. Useful for assessing whether an additive model
#' adequately fits the data or if interaction terms are needed.
#'
#' @param x Matrix with exactly 2 columns of predictor values.
#' @param y Numeric vector of response values.
#' @param adfun Additive model fitting function (default: `adrun`).
#' @param plotfun Plotting function for residuals (default: `lplot`).
#' @param scale Logical. If `TRUE`, scales variables (default: `FALSE`).
#' @inheritParams common-params
#'
#' @details
#' The function first fits an additive model using `adfun`, then plots the
#' residuals using `plotfun`. If the additive model is adequate, residuals
#' should show no systematic pattern. The x argument must have exactly 2 columns.
#'
#' @return Output from the plotting function (usually a plot).
#'
#' @seealso \code{\link{adrun}}, \code{\link{lplot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Check for interaction in regression
#' x <- matrix(rnorm(100), ncol=2)
#' y <- 2*x[,1] + 3*x[,2] + x[,1]*x[,2] + rnorm(50)  # Has interaction
#' riplot(x, y)
#' }
riplot<-function(x,y,adfun=adrun,plotfun=lplot,eout=FALSE,xout=TRUE,scale=FALSE){
#
# Plot used to investigate regression interaction
# (the extent a generalized additive model does not fit data).
# Compute additive fit, plot residuals
# versus x, an n by 2 matrix.
#
if(!is.matrix(x))stop(" x must be a matrix")
if(ncol(x)!=2)stop(" x must have two columns only")
yhat<-adfun(x,y,pyhat=TRUE,eout=eout,xout=xout,plotit=FALSE)
plotfun(x,y-yhat,eout=eout,xout=xout,scale=scale)
}

#' Regression Residual Plot with Running Interval Smoother
#'
#' Fits a smooth using all predictors except one, computes residuals, then
#' plots the smooth of residuals versus the excluded predictor. Useful for
#' examining the contribution of individual predictors after accounting for others.
#'
#' @param x Matrix of predictor values.
#' @param y Numeric vector of response values.
#' @param pv Index of predictor variable to exclude from initial fit and use
#'   for residual plot (default: 1).
#' @param fr Span parameter for smoothing (default: NA, automatically chosen).
#' @param efr Span for error smoothing (default: 0.5).
#' @param theta Viewing angle for 3D plots (default: 50).
#' @param phi Viewing angle for 3D plots (default: 25).
#' @param scale Logical. If `TRUE`, scales variables (default: `TRUE`).
#' @param expand Expansion factor for perspective plots (default: 0.5).
#' @param varfun Function for computing variance (default: `pbvar`).
#' @param STAND Logical. If `TRUE`, standardizes outlier detection (default: `TRUE`).
#' @param nmin Minimum sample size for smoothing (default: 0).
#' @param out Deprecated parameter.
#' @param zscale Logical. If `TRUE`, scales z-axis (default: `FALSE`).
#' @param duplicate How to handle duplicate values (default: 'error').
#' @param ticktype Tick mark style (default: 'simple').
#' @param LP Logical. If `TRUE`, applies lowess smoothing (default: `TRUE`).
#' @inheritParams common-params
#'
#' @details
#' The function first fits a smooth to y using all predictors except the one
#' specified by `pv`, then plots residuals from this fit against the excluded
#' predictor. This helps assess whether the excluded predictor contributes
#' additional information beyond the other predictors.
#'
#' @return Output from `rplot()` applied to residuals.
#'
#' @seealso \code{\link{rplot}}, \code{\link{riplot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Examine residuals after excluding first predictor
#' x <- matrix(rnorm(150), ncol=3)
#' y <- 2*x[,1] + 3*x[,2] + rnorm(50)
#' rplot.res(x, y, pv=1)  # Residuals vs X1 after fitting with X2
#' }
rplot.res<-function(x,y,pv=1,est = tmean, scat = TRUE, fr = NA, plotit = TRUE,
pyhat = FALSE, efr = 0.5, theta = 50, phi = 25, scale = TRUE,
expand = 0.5, SEED = TRUE, varfun = pbvar, outfun = outpro,STAND=TRUE,
nmin = 0, xout = FALSE, out = FALSE, eout = FALSE, xlab='X',
ylab ='Y',zscale=FALSE,zlab=' ', pr=TRUE,duplicate='error',
ticktype='simple',LP=TRUE,...){
#
# Apply rplot excluding the independent variable indicated by the argument
# pv.
# So pv=1 means will exclude the first predictor.
# Fit a smooth using the remaing variables, compute the residuals, then plot
# the smooth using the residuals as the dependent variable and
# the variables indicated by pv as the independent variables.
#
xy=na.omit(cbind(x,y))
p=ncol(x)
p1=p+1
if(xout){
flag=outfun(xy[,1:p],plotit=FALSE,STAND=STAND,...)$keep
xy=xy[flag,]
}
x=xy[,1:p]
y=xy[,p1]
res=y-rplot(x[,1:2],y,est=est,scat=scat,varfun=varfun,expand=expand,nmin=nmin,
pyhat=TRUE,plotit=FALSE,fr=fr,xout=FALSE)$yhat
outp=rplot(x[,pv],res,fr=fr,xout=FALSE,efr=efr,theta=theta,phi=phi,
scale=scale,SEED=SEED,xlab=xlab,ylab=ylab,zlab=zlab,pr=FALSE,
ticktype=ticktype,LP=LP,...)
outp
}

#' Running Interval Smoother with Confidence Band (Trimmed Mean)
#'
#' Plots a running interval smoother based on trimmed means with an approximate
#' simultaneous confidence band. Provides both point estimates and confidence
#' intervals at specified points along the predictor range.
#'
#' @param fr Span parameter controlling amount of smoothing (default: 0.5).
#'   Lower values = less smoothing.
#' @param p.crit Critical p-value for confidence band (default: NA, automatically
#'   estimated). Pre-computed values available for npts=10 or 25 and alpha=0.05.
#' @param scat Logical. If `TRUE`, includes scatter plot of data points (default: `TRUE`).
#' @param npts Number of points where confidence intervals are computed (default: 25).
#'   Pre-computed critical values available for 10 or 25 points.
#' @param low.span Span parameter for lowess smoothing of confidence band (default: 2/3).
#' @param nmin Minimum sample size required in each interval (default: 12).
#' @param LPCI Logical. If `TRUE`, smooths the confidence intervals via lowess
#'   for a smoother appearance (default: `FALSE`).
#' @param MID Logical. If `TRUE`, computes intervals in middle range of data
#'   (default: `TRUE`).
#' @param pch Plotting character for scatter plot (default: '.').
#' @inheritParams common-params
#'
#' @details
#' The function creates a confidence band with simultaneous probability coverage
#' of 1-alpha across `npts` points. For n >= 50, critical values are pre-computed
#' for npts=10 or 25 with alpha=0.05. For other combinations, critical values
#' are estimated via simulation (slower).
#'
#' The running interval smoother uses a trimmed mean (controlled by `tr`) within
#' a moving window (controlled by `fr`). Set LPCI=TRUE for a smoother confidence
#' band appearance.
#'
#' @return List with components:
#' \item{output}{Matrix with columns: pts, y.hat, ci.low, ci.up, n.used}
#' \item{str}{Strength of association measure}
#' \item{n}{Total sample size}
#' \item{n.keep}{Sample size after outlier removal}
#'
#' @seealso \code{\link{rplotCIS}}, \code{\link{rplotCIsmm}}, \code{\link{rplotpbCI}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Running interval smoother with confidence band
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' rplotCI(x, y, fr=0.5)
#'
#' # Smoother confidence band
#' rplotCI(x, y, fr=0.5, LPCI=TRUE)
#' }
rplotCI<-function(x,y,tr=.2,fr=.5,p.crit=NA,plotit=TRUE,scat=TRUE,
SEED=TRUE,pyhat=FALSE,npts=25,xout=FALSE,
xlab='X',ylab='Y',low.span=2/3,nmin=12,pr=TRUE,
outfun=outpro,LPCI=FALSE,MID=TRUE,alpha=.05,pch='.',...){
#
# Confidence interval for running  interval smoother based on a trimmed mean.
# Unlike rplot, includes an approximate  confidence band having simultaneous probability
# coverage equal to 1-alpha.   More precisely, the simultaneous probability
# is for K=npts points
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing
#
# To specify the points where confidence intervals are computed,
# use rplotCIsmm
#
if(pr){
if(!LPCI)print('To get smoother plot, set LPCI=TRUE')
}
m<-cbind(x,y)
if(ncol(m)>2)stop('Only one independent variable can be used')
m<-elimna(m)
x=m[,1]
y=m[,2]
if(xout){
xy=cbind(x,y)
flag=outfun(x,plotit=FALSE)$keep
x=xy[flag,1]
y=xy[flag,2]
}
n.used=NA
n=length(y)
if(n<50)stop('Need at least n=50')
nv=c(50,  60,  70,  80, 100,
150, 200, 300, 400, 500, 600, 800, 1000)
if(npts==25) pv=c(0.004846408,
0.004553274,
0.004236101,
0.004099674,
 0.00353898,   #n=100
 0.003366984,
0.003038767,
 0.003061386,
  0.002793521,
  0.002479689,
 0.002606313,
 0.0026630370,
 0.002836043)
if(npts==10) pv=c(
0.007612451,
0.008383655,
0.006992874,
 0.0068073,
0.005733889,
0.005767139,
0.006130155,
0.005447419,
0.005555951,
0.005228471,
0.005642503,
0.005402162,
0.005569217)
FLAG=FALSE
if(npts==25 || npts==10)FLAG=TRUE
if(alpha!=.05 || !FLAG){
if(is.na(p.crit)){
print('p.crit must be estimated, execution time might be high')
print('Or use the R function rplotCIsmm')
}
p.crit=rplotCITAP.pv(n,tr=tr,fr=fr,alpha=alpha,nmin=nmin,npts=npts,nreps=nreps)
}
rem.n=n
if(n>1000)n=1000
if(is.na(p.crit))p.crit=lplot.pred(1/nv,pv,1/n)$yhat
n=rem.n
xord=order(x)
x=x[xord]
y=y[xord]
infit=rplot(x,y,tr=tr,xout=FALSE,plotit=plotit,LP=LPCI,fr=fr,pr=FALSE,pyhat=TRUE,xlab=xlab,
ylab=ylab)
rmd=infit$pyhat
m<-cbind(x,y)
if(ncol(m)>2)stop('One covariate only is allowed with this function')
m<-elimna(m)
nv=nrow(m)
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
x=m[,1]
y=m[,2]
n.keep=length(y)
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr,nmin=nmin)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
y.hat=NA
civ=matrix(NA,nrow=npts,ncol=2)
for(i in 1:length(pts)){
xx=y[near(x,pt=pts[i],fr)]
doit=trimci(xx,tr=tr,alpha=p.crit,pr=FALSE)
civ[i,]=doit$ci
y.hat[i]=doit$estimate
n.used[i]=doit$n
}
up=civ[,2]
low=civ[,1]
if(plotit){
if(LPCI){
up=lplot(pts,up,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
y.hat=lplot(pts,y.hat,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
low=lplot(pts,low,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
}
lines(pts,up,lty=2)
lines(pts,low,lty=2)
}
if(pyhat){output<-cbind(pts,y.hat,low,up,n.used)
dimnames(output)=list(NULL,c('pts','y.hat','ci.low','ci.up','n.used'))
}
if(!pyhat)output<-'Done'
list(output=output,str=infit$Strength.Assoc,n=nv,n.keep=n.keep)
}

#' Simple Running Interval Smoother with Confidence Band
#'
#' Computes a confidence band for a running interval smoother based on trimmed
#' means using a simple method. Unlike `rplotCI`, this function does not adjust
#' for familywise error rate across multiple points.
#'
#' @param fr Span parameter controlling amount of smoothing (default: 0.8).
#'   Higher values = more smoothing.
#' @param dfmin Minimum degrees of freedom required for confidence intervals
#'   (default: 8). Points with df < dfmin are excluded.
#' @param LP Logical. If `TRUE`, applies lowess smoothing to the confidence
#'   band for smoother appearance (default: `TRUE`).
#' @param pch Plotting character for scatter plot (default: '.').
#' @inheritParams common-params
#'
#' @details
#' This is a simpler alternative to `rplotCI` that computes pointwise confidence
#' intervals without adjusting for simultaneous coverage. The confidence band
#' uses Student's t-distribution based on the degrees of freedom at each point.
#'
#' The function excludes points where degrees of freedom fall below `dfmin` to
#' ensure adequate precision. Set LP=TRUE for a smoother confidence band appearance.
#'
#' Note: `rplotCI` provides better control of familywise error rate (FWE) across
#' multiple points.
#'
#' @return List with components:
#' \item{output}{Matrix with columns: x, y.hat, ci.low, ci.up}
#' \item{str}{Strength of association measure}
#' \item{n}{Total sample size}
#' \item{n.keep}{Sample size after outlier removal}
#'
#' @seealso \code{\link{rplotCI}}, \code{\link{rplotpbCI}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Simple confidence band for running smoother
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' rplotCIS(x, y, fr=0.8)
#' }
rplotCIS<-function(x,y,tr=.2,fr=.8,plotit=TRUE,scat=TRUE,pyhat=FALSE,SEED=TRUE,dfmin=8,
eout=FALSE,xout=FALSE,xlab='x',ylab='y',outfun=out,LP=TRUE,alpha=.05,pch='.',...){
#
# A simple method for computing a confidence band based on
# running  interval smoother and a trimmed mean.
#
#  rplotCI adjusts the band so that FWE=1-alpha
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing
plotit<-as.logical(plotit)
scat<-as.logical(scat)
str=rplot(x,y,tr=tr,xout=xout,plotit=FALSE,LP=LP,fr=fr,pr=FALSE)$Strength.Assoc
m<-cbind(x,y)
if(ncol(m)>2)stop('Only one independent variable can be used')
m<-elimna(m)
nv=nrow(m)
if(eout && xout)stop('Not allowed to have eout=xout=T')
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
x=m[,1]
y=m[,2]
n.keep=length(y)
rmd<-c(1:length(x))
for(i in 1:length(x))rmd[i]<-mean(y[near(x,x[i],fr=fr)],tr=tr)
sedf=runse(x,y,fr=fr,tr=tr,pts=x,SEED=SEED)
df=sedf$df
flag=df>dfmin
se=sedf$se
low=rmd[flag]-qt(1-alpha/2,df[flag])*se[flag]
up=rmd[flag]+qt(1-alpha/2,df[flag])*se[flag]
rmd=rmd[flag]
x=x[flag]
y=y[flag]
if(plotit){
ord=order(x)
x=x[ord]
rmd=rmd[ord]
up=up[ord]
low=low[ord]
if(LP){
rmd=lplot(x,rmd,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
up=lplot(x,up,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
low=lplot(x,low,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
if(scat){
plot(c(x,x),c(y,rmd),xlab=xlab,ylab=ylab,type='n')
lines(x,up,lty=2)
lines(x,low,lty=2)
points(x,y,pch=pch)
}
if(!scat)plot(c(x,x),c(y,rmd),type='n',ylab=ylab,xlab=xlab)
points(x,rmd,type='n')
sx<-sort(x)
xorder<-order(x)
sysm<-rmd[xorder]
lines(sx,sysm)
lines(x,up,lty=2)
lines(x,low,lty=2)
}
if(pyhat){output<-cbind(x,rmd,low,up)
dimnames(output)=list(NULL,c('x','y.hat','ci.low','ci.up'))
}
if(!pyhat)output<-'Done'
list(output=output,str=str,n=nv,n.keep=n.keep)
}

#' Bootstrap Running Interval Smoother with Confidence Band
#'
#' Running interval smoother based on any measure of location with confidence
#' band computed via percentile bootstrap. More flexible than `rplotCI` as it
#' works with any location estimator, not just trimmed means.
#'
#' @param est Estimator function for location (default: `onestep` for one-step
#'   M-estimator). Can be any function like `mean`, `median`, `tmean`, `mom`.
#' @param fr Span parameter controlling amount of smoothing (default: 1).
#'   Higher values = more smoothing.
#' @param LP Logical. If `TRUE`, applies lowess smoothing to the confidence
#'   band for smoother appearance (default: `TRUE`).
#' @inheritParams common-params
#'
#' @details
#' Unlike `rplotCI` which is limited to trimmed means, this function can use
#' any location estimator specified by `est`. Confidence intervals are computed
#' via percentile bootstrap at each point along the curve.
#'
#' The bootstrap approach makes this function more computationally intensive
#' than `rplotCI` or `rplotCIS`, but provides greater flexibility in choice
#' of estimator. Set LP=TRUE for a smoother confidence band appearance.
#'
#' @return List with component:
#' \item{output}{Fitted values at each point (if pyhat=TRUE), otherwise "Done"}
#'
#' @seealso \code{\link{rplotCI}}, \code{\link{rplotCIS}}, \code{\link{onesampb}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Bootstrap confidence band with median
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' rplotpbCI(x, y, est=median, nboot=500)
#'
#' # Using one-step M-estimator
#' rplotpbCI(x, y, est=onestep, fr=0.8)
#' }
rplotpbCI<-function(x,y,est=onestep,fr=1,plotit=TRUE,scat=TRUE,pyhat=FALSE,
xout=FALSE,xlab='x',ylab='y',outfun=out,LP=TRUE,alpha=.05,
nboot=500,SEED=TRUE,...){
#
# running  interval smoother based on any measure of location
# Unlike rplotCI, uses a percentile bootstrap
# method to get a confidence band
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing
plotit<-as.logical(plotit)
scat<-as.logical(scat)
m<-cbind(x,y)
if(ncol(m)>2)stop('Only one independent variable can be used')
m<-elimna(m)
x=m[,1]
y=m[,2]
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
x=m[,1]
y=m[,2]
low=rep(NA,length(y))
up=rep(NA,length(y))
rmd<-NA
for(i in 1:length(x)){
sel=y[near(x,x[i],fr)]
temp=onesampb(sel,est=est,nboot=nboot,alpha=alpha,SEED=SEED,...)
low[i]=temp$ci[1]
up[i]=temp$ci[2]
rmd[i]=temp$estimate
}
all=elimna(cbind(x,low,up,y,rmd))
x=all[,1]
low=all[,2]
up=all[,3]
y=all[,4]
rmd=all[,5]
if(plotit){
ord=order(x)
x=x[ord]
y=y[ord]
rmd=rmd[ord]
up=up[ord]
low=low[ord]
if(LP){
rmd=lplot(x,rmd,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
up=lplot(x,up,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
low=lplot(x,low,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
if(scat){
plot(c(x,x),c(y,rmd),xlab=xlab,ylab=ylab,type='n')
points(x,y)
lines(x,up,lty=2)
lines(x,low,lty=2)
}
if(!scat)plot(c(x,x),c(y,rmd),type='n',ylab=ylab,xlab=xlab)
points(x,rmd,type='n')
sx<-sort(x)
xorder<-order(x)
sysm<-rmd[xorder]
lines(sx,sysm)
lines(x,up,lty=2)
lines(x,low,lty=2)
}
if(pyhat)output<-rmd
if(!pyhat)output<-'Done'
list(output=output)
}

#' Running Interval Smoother with Confidence Band (Median-Based)
#'
#' Running interval smoother based on the Harrell-Davis median estimator with
#' simultaneous confidence band. Provides a robust alternative to `rplotCI` for
#' non-normal data.
#'
#' @param est Estimator function (currently must be `hd` for Harrell-Davis,
#'   default: `hd`).
#' @param fr Span parameter controlling amount of smoothing (default: 0.5).
#'   For strong associations, consider fr=0.2 for less smoothing.
#' @param pts Optional numeric vector of points where CIs are computed (default: NA,
#'   automatically chosen).
#' @param npts Number of points where confidence intervals are computed (default: 25).
#' @param low.span Span parameter for lowess smoothing of confidence band (default: 2/3).
#' @param nmin Minimum sample size required in each interval (default: 16).
#' @param LP Logical. If `TRUE`, smooths the fitted curve via lowess (default: `TRUE`).
#' @param LPCI Logical. If `TRUE`, smooths the confidence intervals via lowess
#'   for smoother appearance (default: `FALSE`).
#' @param MID Logical. If `TRUE`, computes intervals in middle range of data
#'   (default: `TRUE`).
#' @param pch Plotting character for scatter plot (default: '.').
#' @inheritParams common-params
#'
#' @details
#' This function is similar to `rplotCI` but uses the Harrell-Davis median
#' estimator instead of trimmed means, making it more robust for heavy-tailed
#' distributions. The confidence band has simultaneous probability coverage of
#' 1-alpha across `npts` points using Bonferroni adjustment.
#'
#' The current version requires `est=hd`. Future versions may support additional
#' estimators.
#'
#' @return List with components:
#' \item{output}{Matrix with columns: pts, y.hat, ci.low, ci.up}
#' \item{str}{Strength of association measure}
#' \item{n}{Total sample size}
#' \item{n.keep}{Sample size after outlier removal}
#'
#' @seealso \code{\link{rplotCI}}, \code{\link{hd}}, \code{\link{sint}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Median-based running smoother with confidence band
#' x <- rnorm(100)
#' y <- 2*x + rt(100, df=3)  # Heavy-tailed errors
#' rplotCIM(x, y, fr=0.5)
#' }
rplotCIM<-function(x,y,est=hd,fr=.5,p.crit=NA,plotit=TRUE,scat=TRUE,
pyhat=FALSE, pts=NA,npts=25,xout=FALSE,
xlab='X',ylab='Y',low.span=2/3,nmin=16,
outfun=out,LP=TRUE,LPCI=FALSE,MID=TRUE,alpha=.05,pch='.',...){
#
# Confidence interval for running  interval smoother based on a median
# Unlike rplot, includes a confidence band having simultaneous probability
# coverage equal to 1-alpha.
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing. If the association is relatively strong, might want to use fr=.2
#
chk=FALSE
if(identical(est,hd))chk=TRUE
if(!chk)stop('Current version, argument est must be hd')
n=length(y)
if(n<50)stop('Need at least n=50')
xord=order(x)
x=x[xord]
y=y[xord]
infit=rplot(x,y,est=est,xout=xout,plotit=plotit,LP=LP,fr=fr,pr=FALSE,pyhat=TRUE,xlab=xlab,ylab=ylab)
rmd=infit$yhat
m<-cbind(x,y)
if(ncol(m)>2)stop('One covariate only is allowed with this function')
m<-elimna(m)
nv=nrow(m)
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
x=m[,1]
y=m[,2]
n.keep=length(y)
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr,nmin=nmin)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
flag=duplicated(pts)
npts=length(pts)
civ=matrix(NA,nrow=npts,ncol=2)
for(i in 1:length(pts)){
xx=y[near(x,pt=pts[i],fr)]
civ[i,]=sint(xx,alpha=alpha/npts)
}
up=civ[!flag,2]
low=civ[!flag,1]
if(plotit){
if(LPCI){
up=lplot(pts,up,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
low=lplot(pts,low,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
}
pts=pts[!flag]
lines(pts,up,lty=2)
lines(pts,low,lty=2)
}
if(pyhat){output<-cbind(pts,rmd[!flag],low,up)
dimnames(output)=list(NULL,c('pts','y.hat','ci.low','ci.up'))
}
if(!pyhat)output<-'Done'
list(output=output,str=infit$Strength.Assoc,n=nv,n.keep=n.keep)
}

#' Running Interval Smoother with Flexible Confidence Band
#'
#' Computes confidence band for running interval smoother based on trimmed means,
#' using the Studentized maximum modulus distribution. More flexible than `rplotCI`
#' in terms of alpha levels and number of points, though intervals may be slightly
#' wider.
#'
#' @param pts Optional numeric vector of specific points where CIs are computed
#'   (default: NULL, automatically chosen).
#' @param dfmin Minimum degrees of freedom required for CIs (default: 2).
#' @param pch Plotting character for scatter plot (default: '.').
#' @inheritParams rplotCI
#'
#' @details
#' This function offers more flexibility than `rplotCI`: it works with any alpha
#' level and any number of points (not limited to 10 or 25). The confidence band
#' uses the Studentized maximum modulus distribution for simultaneous coverage.
#'
#' Trade-off: `rplotCI` provides tighter confidence intervals when alpha=0.05 and
#' npts=10 or 25, but `rplotCIsmm` offers greater flexibility with minimal
#' computational overhead.
#'
#' @return List with components:
#' \item{output}{Matrix with columns: pts, y.hat, ci.low, ci.up}
#' \item{str}{Strength of association measure}
#' \item{n}{Total sample size}
#' \item{n.keep}{Sample size after outlier removal}
#'
#' @seealso \code{\link{rplotCI}}, \code{\link{qsmm}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Flexible confidence band with custom alpha
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' rplotCIsmm(x, y, alpha=0.01, npts=20)
#' }
rplotCIsmm<-function(x,y,tr=.2,fr=.5,plotit=TRUE,scat=TRUE,pyhat=FALSE,SEED=TRUE,
dfmin=2,pts=NULL,npts=25,nmin=12,
eout=FALSE,xout=FALSE,xlab='x',ylab='y',outfun=out,LP=TRUE,MID=TRUE,alpha=.05,pch='.',...){
#
# Confidence interval for running  interval smoother based on a trimmed mean.
#
#  rplotCI will provide shorter and more accurate confidence intervals but
#  is limited to 10 or 25 points and alpha=.05.
#  This functions returns confidence intervals that are generally a bit wider
#  but it has low execution time if alpha differs from 0.5 or there is interest
#  using something other than 10 or 25 points.
#
# Unlike rplot,a confidence band based on the Studentized maximum modulus dist
# is computed,
# unless alpha is not equal to .05 or the number of confidence intervals
# is greater than npts=28, in which case the distribution of the max of npts
#  random variables is used.
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing
#
xord=order(x)
x=x[xord]
y=y[xord]
if(!is.null(pts))pts=sort(pts)
str=rplot(x,y,tr=tr,xout=xout,plotit=FALSE,LP=LP,fr=fr,pr=FALSE)$Strength.Assoc
m<-cbind(x,y)
if(ncol(m)>2)stop('To get a smooth with more than one covariate, use rplot')
m<-elimna(m)
nv=nrow(m)
if(eout && xout)stop('Not allowed to have eout=xout=T')
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
if(is.null(pts)){
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr,nmin=nmin)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
}
x=m[,1]
y=m[,2]
n.keep=length(y)
if(is.null(pts)){
if(!MID)pts=seq(min(x),max(x),length.out=npts)
vv=idealf(x)
if(MID)pts=seq(vv$ql,vv$qu,length.out=npts)
}
rmd=NA
for(i in 1:length(pts))rmd[i]<-mean(y[near(x,pts[i],fr)],tr=tr)
sedf=runse(x,y,fr=fr,tr=tr,pts=pts,SEED=SEED)
df=sedf$df
flag=df>dfmin
se=sedf$se
ntest=length(df[flag])
mdif=min(df[flag])
crit=NA
dfval=df[flag]
for(it in 1:ntest)crit[it]=qsmm(1-alpha,ntest,dfval[it])
low=rmd[flag]-crit*se[flag]
up=rmd[flag]+crit*se[flag]
ptsall=pts
rmdall=rmd
rmd=rmd[flag]
pts=pts[flag]
if(plotit){
ord=order(x)
x=x[ord]
y=y[ord]
if(LP){
rmd=lplot(pts,rmd,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
up=lplot(pts,up,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
low=lplot(pts,low,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
plot(c(x,pts),c(y,rmd),xlab=xlab,ylab=ylab,type='n')
if(scat)points(x,y,pch=pch)
lines(pts,up,lty=2)
lines(pts,low,lty=2)
lines(pts,rmd)
}
if(pyhat){output<-cbind(pts,rmd,low,up)
dimnames(output)=list(NULL,c('pts','y.hat','ci.low','ci.up'))
}
if(!pyhat)output<-'Done'
list(output=output,str=str,n=nv,n.keep=n.keep)
}

#' Determine Critical P-value for rplotCI via Simulation (TAP Method)
#'
#' @description
#' Internal simulation function to determine the critical p-value threshold for
#' \code{\link{rplotCI}} when using the TAP (Trimmed ANCOVA with Pooling) method.
#' Uses simulation under the null hypothesis of no association.
#'
#' @param n Sample size for simulation.
#' @param nreps Number of simulation replications (default: 2000).
#' @param alpha Significance level (default: 0.05).
#' @param npts Number of evaluation points along the X-axis (default: 25).
#' @inheritParams common-params
#' @param fr Span for running interval smoother (default: 5).
#' @param MC Logical; if `TRUE`, uses parallel processing via \code{mclapply}.
#' @param nmin Minimum number of points in a neighborhood (default: 12).
#' @param SEED Logical; if `TRUE`, sets random seed for reproducibility.
#' @param LP Logical; related to local polynomial settings (default: FALSE).
#'
#' @return Numeric value: the alpha-quantile of the minimum p-values across
#'   simulation replications, representing the critical p-value threshold.
#'
#' @seealso \code{\link{rplotCI}}, \code{\link{rplotCITAP.sub}}
#'
#' @keywords internal
rplotCITAP.pv<-function(n,nreps=2000,alpha=.05,npts=25,tr=.2,fr=5,MC=FALSE,nmin=12,SEED=TRUE,LP=FALSE){
if(SEED)set.seed(2)
pvals=NA
xy=list()
for (i in 1:nreps){
xy[[i]]=rmul(n)
}
if(!MC)pvals=lapply(xy,rplotCITAP.sub,npts=npts,tr=tr,fr=fr,alpha=alpha,nmin=nmin)
if(MC){
pvals=mclapply(xy,rplotCITAP.sub,npts=npts,tr=tr,fr=fr,alpha=alpha,nmin=nmin)
}
pvals=matl(pvals)
pv=hd(pvals,alpha)
pv
}

#' Bootstrap Helper for rplotCITAP.pv Critical P-value Simulation
#'
#' @description
#' Internal helper function for \code{\link{rplotCITAP.pv}}. Computes the minimum
#' p-value across evaluation points when fitting a running interval smoother with
#' confidence band to simulated bivariate normal data.
#'
#' @param xy Matrix with 2 columns: predictor X and response Y (simulated data).
#' @inheritParams common-params
#' @param fr Span for running interval smoother.
#' @param SEED Logical; if `TRUE`, sets random seed.
#' @param nmin Minimum number of points in neighborhood for estimation.
#' @param pts Evaluation points (if NA, determined from ANCOVA output).
#' @param npts Number of evaluation points to use.
#' @param LP Logical; related to local polynomial settings.
#' @param alpha Significance level for trimmed mean confidence intervals.
#' @param xout Logical; if `TRUE`, removes outliers from predictors.
#'
#' @return Numeric value: the minimum p-value across all evaluation points.
#'
#' @keywords internal
rplotCITAP.sub<-function(xy,tr=.2,fr=NA,SEED=TRUE,nmin=12,
pts=NA,npts=25,LP=FALSE,alpha=.05,xout=FALSE,...){
#
# prediction interval running  interval smoother based on a trimmed mean.
# Unlike rplot, includes a confidence band.
#
x=xy[,1]
y=xy[,2]
xord=order(x)
x=x[xord]
y=y[xord]
infit=rplot(x,y,tr=tr,xout=xout,plotit=FALSE,LP=LP,fr=fr,pr=FALSE,pyhat=TRUE,nmin=nmin)
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
rmd=infit$pyhat
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr,nmin=nmin)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
y.hat=NA
pv=NA
civ=matrix(NA,nrow=npts,ncol=2)
for(i in 1:length(pts)){
doit=trimci(y[near(x,pts[i],fr)],tr=tr,alpha=alpha,pr=FALSE)
pv[i]=doit$p.value
}
min(pv)
}

#' Determine Critical P-value for rplotCI via Simulation (Version 2)
#'
#' @description
#' Internal simulation function to determine the critical p-value threshold for
#' \code{\link{rplotCI}} alternative method. Simulates under the null hypothesis
#' of no association and returns the alpha-quantile of minimum p-values.
#'
#' @param n Sample size for simulation.
#' @param nreps Number of simulation replications (default: 4000).
#' @param alpha Significance level (default: 0.05).
#' @inheritParams common-params
#' @param fr Span for running interval smoother (default: 0.5).
#' @param MC Logical; if `TRUE`, uses parallel processing (default: TRUE).
#' @param nmin Minimum number of points in neighborhood (default: 12).
#' @param SEED Logical; if `TRUE`, sets random seed for reproducibility.
#'
#' @return Numeric value: the alpha-quantile of the minimum p-values, used as
#'   a critical threshold for simultaneous inference.
#'
#' @seealso \code{\link{rplotCI}}, \code{\link{rplotCIv2.sub}}
#'
#' @keywords internal
rplotCIv2.pv<-function(n,nreps=4000,alpha=.05,tr=.2,fr=.5,
MC=TRUE,nmin=12,SEED=TRUE){
if(SEED)set.seed(2)
pvals=NA
xy=list()
for (i in 1:nreps){
xy[[i]]=rmul(n)
}
if(!MC)pvals=lapply(xy,rplotCIv2.sub,tr=tr,fr=fr,nmin=nmin)
if(MC){
pvals=mclapply(xy,rplotCIv2.sub,tr=tr,fr=fr,nmin=nmin)
}
pvals=matl(pvals)
pv=hd(pvals,alpha)
pv
}

#' Bootstrap Helper for rplotCIv2.pv Critical P-value Simulation
#'
#' @description
#' Internal helper function for \code{\link{rplotCIv2.pv}}. Computes the minimum
#' p-value across X values where neighborhoods contain at least `nmin` points.
#'
#' @param xy Matrix with 2 columns: predictor X and response Y (simulated data).
#' @param nmin Minimum number of points required in a neighborhood.
#' @param tr Trimming proportion for trimmed mean confidence intervals.
#' @param fr Span for defining neighborhoods around X values.
#'
#' @return Numeric value: the minimum p-value across qualifying evaluation points.
#'
#' @keywords internal
rplotCIv2.sub<-function(xy,nmin,tr,fr){
x=xy[,1]
y=xy[,2]
n=length(y)
nv=NA
for(j in 1:n)nv[j]=sum(near(x,x[j],fr=fr))
pts=x[nv>=nmin]
n.keep=length(pts)
for(j in 1:n.keep)nv[j]=sum(near(x,x[j],fr=fr))
pts=x[nv>=nmin]
rmd=NA
for(i in 1:length(pts))rmd[i]<-trimci(y[near(x,pts[i],fr)],tr=tr,pr=FALSE)$p.value
pv=min(rmd)
pv
}

#' Cross-Validation Prediction Error for Running Interval Smoother
#'
#' Estimates prediction error for a running interval smoother using
#' leave-one-out cross-validation. Useful for assessing model fit quality
#' and comparing different smoothing parameters.
#'
#' @param x Numeric vector or matrix of predictor values.
#' @param y Numeric vector of response values.
#' @param fr Span for running interval smoother. If `NA` (default), uses 0.8
#'   for univariate and 1.0 for multivariate predictors.
#' @param varfun Function to measure variation in predicted Y values
#'   (default: `pbvar` for percentage bend variance).
#' @param est Measure of location used by the smoother (default: `tmean`
#'   for trimmed mean).
#' @inheritParams common-params
#' @param corfun Correlation function (currently unused, for future extension).
#'
#' @details
#' For each observation, the function:
#' 1. Removes that observation
#' 2. Fits a running interval smoother to the remaining data
#' 3. Predicts the removed observation
#' 4. Computes residuals (observed - predicted)
#' 5. Measures variation in residuals using `varfun`
#'
#' Lower values indicate better predictive performance. Can be used with
#' multivariate predictors via `rung3hat()`.
#'
#' @return List with component:
#' \item{VAR.Y.HAT}{Measure of variation in cross-validated residuals.}
#'
#' @seealso \code{\link{rplot}}, \code{\link{runhat}}, \code{\link{rung3hat}}
#'
#' @export
#' @examples
#' # Assess prediction error for different smoothing spans
#' x <- rnorm(100)
#' y <- sin(x) + rnorm(100, sd=0.2)
#' rplotCV(x, y, fr=0.5)
#' rplotCV(x, y, fr=0.8)
rplotCV<-function(x,y,fr=NA,varfun=pbvar,est=tmean,xout=FALSE,outfun=out,eout=FALSE,corfun=pbvar,...){
#
# Estimate prediction error based on
# a running interval smoother in conjunction with
# a leave-one-out cross validation method
#
#  varfun is the measure of variation used on the predicted Y values.
#  est is the measure of location used by the running interval smoother.
#  The estimate is returned in VAR.Y.HAT
#
m=elimna(cbind(x,y))
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m=m[flag,]
}
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=as.matrix(m[,1:p])
y=m[,p1]
vals=NA
if(is.na(fr)){
if(p==1)fr=.8
if(p>1)fr=1
}
if(xout){
keepit<-outfun(x,plotit=FALSE,...)$keep
x<-x[keepit,]
y<-y[keepit]
}
x=as.matrix(x)
for(i in 1:nrow(x)){
if(p==1)vals[i]=runhat(x[-i,],y[-i],fr=fr,est=est,pts=x[i,],...)
if(p>1)vals[i]=rung3hat(x[-i,],y[-i],fr=fr,pts=t(as.matrix(x[i,])))$rmd
}
dif=y-vals
ans=varfun(elimna(dif))
list(VAR.Y.HAT=ans)
}

#' Running Interval Smoother with Strength of Association
#'
#' Applies a running interval smoother and computes the strength of association
#' (explanatory power) between predictors and response. Works with univariate
#' or multivariate predictors and uses bootstrap for robust estimation.
#'
#' @param x Numeric vector or matrix of predictor values.
#' @param y Numeric vector of response values.
#' @param est Measure of location for the smoother (default: `tmean`).
#' @param fr Span for running interval (default: 1).
#' @param plotit Logical; if `TRUE` (default), creates a plot.
#' @param pyhat Logical; if `TRUE`, returns predicted values (default: `FALSE`).
#' @param nboot Number of bootstrap samples for variance estimation (default: 40).
#' @param atr Amount of trimming for bootstrap variance (default: 0).
#' @param nmin Minimum number of points in neighborhood (default: 0).
#' @param pch Plotting character (default: '*').
#' @param outfun Outlier detection function (default: `outpro`).
#' @param eout Logical; if `TRUE`, removes outliers before smoothing.
#' @param xlab Label for x-axis (default: "X").
#' @param ylab Label for y-axis (default: "Y").
#' @param scat Logical; if `TRUE` (default), shows scatter plot.
#' @param SEED Logical; if `TRUE` (default), sets random seed for reproducibility.
#' @param expand Expansion factor for 3D plots (default: 0.5).
#' @param scale Logical; if `TRUE`, uses Mahalanobis distance for multivariate
#'   predictors (default: `FALSE`).
#' @param STAND Logical; if `TRUE` (default), standardizes predictors.
#' @param varfun Function to measure variation (default: `pbvar`).
#' @param pr Logical; if `TRUE` (default), prints warnings about scaling.
#' @param ticktype Tick mark type for 3D plots (default: "simple").
#' @param theta Rotation angle for 3D plots (default: 50).
#' @param phi Colatitude angle for 3D plots (default: 25).
#' @param ... Additional arguments passed to smoother functions.
#'
#' @details
#' The function measures how well the predictors explain the response by:
#' 1. Fitting a running interval smoother to get predicted values
#' 2. Computing variance of predicted values
#' 3. Computing variance of response values
#' 4. Taking ratio: explanatory power = var(predicted) / var(response)
#'
#' Strength of association is the square root of explanatory power.
#' For multivariate predictors (ncol(x) > 1), uses `run3bo()` for 3D plotting.
#'
#' @return List with components:
#' \item{Strength.Assoc}{Square root of explanatory power (correlation-like measure).}
#' \item{Explanatory.Power}{Proportion of variance explained.}
#' \item{yhat}{Predicted values (if `pyhat=TRUE`, otherwise `NULL`).}
#'
#' @seealso \code{\link{runmbo}}, \code{\link{run3bo}}, \code{\link{rplot}}
#'
#' @export
#' @examples
#' # Univariate predictor
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' rplotsm(x, y, fr=0.8)
#'
#' # Multivariate predictors
#' x <- matrix(rnorm(200), ncol=2)
#' y <- x[,1] + x[,2] + rnorm(100, sd=0.5)
#' rplotsm(x, y, scale=TRUE)
rplotsm<-function(x,y,est=tmean,fr=1,plotit=TRUE,pyhat=FALSE,nboot=40,atr=0,nmin=0,pch='*',
outfun=outpro,eout=FALSE,xlab="X",ylab="Y",scat=TRUE,SEED=TRUE,expand=.5,scale=FALSE,STAND=TRUE,
varfun=pbvar,pr=TRUE,ticktype="simple",theta=50,phi=25,...){
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
if(ncol(x)==1){
val<-runmbo(x,y,est=est,scat=scat,fr=fr,plotit=plotit,pyhat=TRUE,STAND=STAND,pch=pch,
xlab=xlab,ylab=ylab,eout=eout,nboot=nboot,outfun=outfun,SEED=SEED,atr=atr,...)
}
if(ncol(x)>1){
if(ncol(x)==2 && !scale){
if(pr){
print("scale=F is specified.")
print("If there is dependence, use scale=T")
}}
if(ncol(x)>2)plotit<-F
val<-run3bo(x,y,est=est,fr=fr,nmin=nmin,plotit=plotit,pyhat=TRUE,phi=phi,
theta=theta,xlab=xlab,ylab=ylab,ticktype=ticktype,STAND=STAND,
SEED=SEED,expand=expand,scale=scale,nboot=nboot,...)
val<-val$output
}
E.power<-varfun(val[!is.na(val)])/varfun(y)
if(!pyhat)val <- NULL
E.power=as.numeric(E.power)
list(Strength.Assoc=sqrt(E.power),Explanatory.Power = E.power, yhat = val)
}

#' Running Interval Smoother for Large Datasets
#'
#' Efficient running interval smoother designed for large sample sizes.
#' Plots the regression surface without a scatter plot by using a random
#' subsample for prediction points.
#'
#' @param x Numeric vector or matrix of predictor values.
#' @param y Numeric vector of response values.
#' @param nsub Size of random subsample for prediction points (default: 5000).
#'   Actual subsample size is `min(n, nsub)` where n is the sample size.
#' @param est Measure of location for the smoother (default: `tmean`).
#' @param fr Span for running interval (default: 1).
#' @param xout Logical; if `TRUE`, removes outliers using projection depth.
#' @param xlab Label for x-axis (default: "X").
#' @param ylab Label for y-axis (default: "Y").
#' @param zlab Label for z-axis (default: "").
#' @param ticktype Tick mark type for 3D plots (default: "simple").
#' @param theta Rotation angle for 3D plots (default: 50).
#' @param phi Colatitude angle for 3D plots (default: 25).
#' @param scale Logical; if `TRUE` (default), scales axes.
#' @param expand Expansion factor for 3D plots (default: 0.5).
#' @param SEED Logical; if `TRUE` (default), sets random seed for reproducibility.
#' @param frame Logical; if `TRUE` (default), draws frame around plot.
#'
#' @details
#' This function is designed for situations with large sample sizes where
#' plotting all points is impractical. It:
#' 1. Randomly selects `nsub` points to serve as prediction locations
#' 2. Uses all data to fit the running interval smoother
#' 3. Predicts at the `nsub` points
#' 4. Plots the smooth surface using `lplot()`
#'
#' For multivariate predictors with outlier removal, uses `outpro.depth()`.
#'
#' @return Invisibly returns the result from `lplot()`.
#'
#' @seealso \code{\link{rplot}}, \code{\link{rplot.pred}}, \code{\link{lplot}}
#'
#' @export
#' @examples
#' # Large dataset example
#' set.seed(123)
#' n <- 10000
#' x <- rnorm(n)
#' y <- sin(x) + rnorm(n, sd=0.3)
#' rplotN(x, y, nsub=1000)
rplotN<-function(x,y,nsub=5000,est=tmean,fr=1,xout=FALSE,xlab='X',ylab='Y',zlab='',ticktype = 'simple',theta = 50, phi = 25, scale = TRUE,
    expand = 0.5, SEED = TRUE,frame=TRUE){
#
 # Running interval smoother, good for large sample sizes or plots of the
 # regression surface without a scatter plot.
 #
 # nsub is size of the random sample of the data used to predict outcome using all of the data
 #  
 if(SEED)set.seed(2)
 x=as.matrix(x)
p=ncol(x)
p1=p+1
 xy=cbind(x,y)
 xy=elimna(xy)
 if(xout){
 flag<-outpro.depth(x,plotit=FALSE)$keep                                
x<-x[flag,]                                                                        
x<-as.matrix(x)                                                                    
y<-y[flag]                                                                         
} 
 n=nrow(xy)
 nsub=min(n,nsub)
 id=sample(n,nsub) 
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
w=rplot.pred(x,y,pts=x[id,],fr=fr)$Y.hat 
a=lplot(x[id,],w,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,frame=frame,phi=phi,theta=theta,scale=scale,pr=FALSE)   
}

#' Running Interval Smoother for Binary Outcomes
#'
#' Applies a running interval smoother specifically designed for binary
#' response variables. Estimates and plots P(Y=1|X) as a function of
#' predictor values. Alternative to logistic regression.
#'
#' @param x Numeric vector or matrix of predictor values.
#' @param y Binary response vector (will be converted to 0/1 if needed).
#' @param est Measure of location for probability estimation (default: `mean`).
#' @param scat Logical; if `TRUE` (default), shows scatter plot.
#' @param fr Span for running interval. If `NULL` (default), uses 0.8 for
#'   univariate and 1.2 for multivariate predictors.
#' @param plotit Logical; if `TRUE` (default), creates a plot.
#' @param pyhat Logical; if `TRUE`, returns predicted probabilities
#'   (default: `FALSE`).
#' @param pts Points at which to compute predictions (default: `x`).
#' @param LP Logical; if `TRUE`, uses lowess/loess smoothing (default: `FALSE`).
#' @param theta Rotation angle for 3D plots (default: 50).
#' @param phi Colatitude angle for 3D plots (default: 25).
#' @param scale Logical; if `TRUE` (default), uses Mahalanobis distance for
#'   multivariate predictors.
#' @param expand Expansion factor for 3D plots (default: 0.5).
#' @param SEED Logical; if `TRUE` (default), sets random seed.
#' @param nmin Minimum number of points in neighborhood (default: 0).
#' @inheritParams common-params
#' @param zlab Label for z-axis (default: "P(Y=1)").
#' @param pr Logical; if `TRUE` (default), prints warnings about scaling.
#' @param duplicate How to handle duplicate predictor values in 3D plots
#'   (default: "error"). Use "strip" if duplicates cause plotting issues.
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @details
#' The function:
#' 1. Converts y to 0/1 coding (max value → 1, min value → 0)
#' 2. Applies running interval smoother with `est=mean` to estimate P(Y=1|X)
#' 3. Plots probability surface
#'
#' For multivariate predictors, warns if binary independent variables are
#' detected (use `rplot.binv2()` instead). Prints warnings about scaling
#' options based on whether independence is expected.
#'
#' @return If `pyhat=TRUE`, returns predictions from `rplot.pred()`.
#'   Otherwise returns "DONE".
#'
#' @seealso \code{\link{rplot}}, \code{\link{rplot.binCI}},
#'   \code{\link{rplot.pred}}
#'
#' @export
#' @examples
#' # Binary logistic-like relationship
#' x <- rnorm(200)
#' p <- plogis(2*x)  # True probabilities
#' y <- rbinom(200, 1, p)
#' rplot.bin(x, y, fr=0.5)
rplot.bin<-function(x,y,est=mean,scat=TRUE,fr=NULL,plotit=TRUE,pyhat=FALSE,pts=x,LP=FALSE,
theta=50,phi=25,scale=TRUE,expand=.5,SEED=TRUE,
nmin=0,xout=FALSE,outfun=outpro,xlab=NULL,ylab=NULL,
zlab='P(Y=1)',pr=TRUE,duplicate='error',...){
#
#  This function applies the running interval smoother, but is designed
#  specifically for situations where y is  binary.
#
# duplicate='error'
# In some situations where duplicate values occur, when plotting with
# two predictors, it is necessary to set duplicate='strip'
#
y=chbin2num(y)
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
x<-as.matrix(x)
if(length(unique(y))!=2)stop('y is not binary')
n=length(y)
Y=rep(0,n)
flag=y==max(y)
Y[flag]=1
y=Y
x<-as.matrix(x)
if(ncol(x)==1){
if(is.null(ylab))ylab='P(Y=1)'
if(is.null(xlab))ylab='X'
if(is.null(fr))fr=.8
a=rplot(x,y,est=mean,xout=xout,outfun=outfun,fr=fr,xlab=xlab,ylab=ylab,pr=FALSE,LP=LP)
}
if(ncol(x)>1){
id=chk4binary(x)
Lid=length(id)
if(Lid>0)stop('Binary independent variables detected. Use rplot.binv2')
if(is.null(xlab))xlab='X1'
if(is.null(ylab))ylab='X2'
if(is.null(fr))fr=1.2
if(ncol(x)==2){
if(scale){
if(pr){print('scale=T is specified.')
print('If there is independence, might want to use scale=F')
a=rplot(x,y,est=mean,xout=xout,outfun=outfun,fr=fr,xlab=xlab,ylab=ylab,zlab=zlab,scale=scale,pr=FALSE)
}}}}
if(!pyhat)val <- 'DONE'
if(pyhat)val=rplot.pred(x,y,pts=pts,est=mean,xout=xout,outfun=outfun,fr=fr)
val
}

#' Confidence Intervals for Binary Outcome Probabilities
#'
#' Computes and plots confidence intervals for P(Y=1|X) at intervals along
#' the predictor range. Provides a robust alternative to logistic regression
#' with interval-based probability estimation.
#'
#' @param x Numeric vector of predictor values.
#' @param y Binary response vector (will be converted to 0/1).
#' @param pts Points defining interval centers (default: `NULL` uses all
#'   unique x values). If specified as vector like `c(-1,0,1,2)`, defines
#'   intervals (-1,0), (0,1), (1,2).
#' @param alpha Significance level for confidence intervals (default: 0.05).
#' @param nmin Minimum number of observations required in interval (default: 5).
#' @inheritParams common-params
#' @param fr Span for defining neighborhoods (default: 0.5).
#' @param tr.plot Logical; if `TRUE`, trims plot to middle 80% of points
#'   (default: `FALSE`).
#' @param method Method for binomial CI: "AC" (Agresti-Coull) or "CP"
#'   (Clopper-Pearson). If `NULL` (default), uses "AC" for n<80, "CP" otherwise.
#' @param plotit Logical; if `TRUE` (default), creates a plot.
#' @param LP Logical; if `TRUE` (default), smooths CI bounds with lowess.
#' @param xlab Label for x-axis (default: "X").
#' @param ylab Label for y-axis (default: "P(Y=1|X)").
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @details
#' For each interval center in `pts`:
#' 1. Finds all x values within `fr` fraction of that point
#' 2. Computes proportion of Y=1 in that neighborhood
#' 3. Computes confidence interval using `binom.conf()`
#' 4. Plots estimate and CI bounds
#'
#' The function automatically selects appropriate binomial CI method based on
#' sample size. Lowess smoothing of bounds (`LP=TRUE`) can improve visualization.
#'
#' @return Matrix with columns:
#' \item{n}{Sample size in interval}
#' \item{low.end}{Lower end of x interval}
#' \item{upper.end}{Upper end of x interval}
#' \item{pts}{Interval center}
#' \item{p.hat}{Estimated P(Y=1|X)}
#' \item{ci.low}{Lower confidence bound}
#' \item{ci.up}{Upper confidence bound}
#'
#' @seealso \code{\link{rplot.bin}}, \code{\link{binom.conf}}
#'
#' @export
#' @examples
#' # Binary outcome with confidence intervals
#' x <- rnorm(200)
#' p <- plogis(2*x)
#' y <- rbinom(200, 1, p)
#' rplot.binCI(x, y, alpha=0.05)
rplot.binCI<-function(x,y,pts=NULL,alpha=.05,nmin=5,xout=FALSE,outfun=outpro,fr=.5,tr.plot=FALSE,
method=NULL,plotit=TRUE,LP=TRUE,xlab='X',ylab='P(Y=1|X)',...){
#
#  An alternative to logistic regression.
#
#  For a collection of intervals among the values in
#  x, compute the probability of success and a confidence based on the corresponding y values
#
#  Default: use the deciles to define the intervals
#
#  Example: pts=c(-1,0,1,2). The intervals would be (-1,0), (0,1), (1,2).
#
y=chbin2num(y)
xx<-elimna(cbind(x,y))
x<-xx[,1]
y<-xx[,2]
if(is.null(pts)){
id=duplicated(x)
pts=x[!id]
}
else plotit=FALSE
if(xout){
m<-cbind(x,y)
if(ncol(m)!=2)stop('Only one  explanatory variable is allowed')
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1]
y<-m[,2]
}
n=length(x)
xor=order(x)
x=x[xor]
y=y[xor]
pts=sort(pts)
npts=length(pts)
#
if(length(unique(y))>2)stop('y should be binary')
#
# Determine which method will be used:
if(is.null(method)){
if(n<80)method='AC'
if(n>=80)method='CP'
}

# Next convert y to 0 and 1s  if not already 0 and 1s
yy=rep(0,n)
yc=as.character(y)
flag=yc==max(yc)
yy[flag]=1
y=yy
#
rmd<-matrix(NA,nrow=npts,ncol=7)
for(i in 1:npts){
isub=near(x,pts[i],fr)
if(sum(isub)>=nmin){
z=y[isub]
v=binom.conf(y=z,method=method,alpha=alpha,pr=FALSE)
rmd[i,1]=v$n
rmd[i,2]=min(x[isub])
rmd[i,3]=max(x[isub])
rmd[i,4]=pts[i]
rmd[i,5]=v$phat
rmd[i,6:7]=v$ci
}}
rs=elimna(rmd)
dimnames(rmd)=list(NULL,c('n','low.end','upper.end','pts','p.hat','ci.low','ci.up'))
if(plotit){
if(tr.plot){
v=quantile(rmd[,4],probs=c(.1,.9),na.rm=TRUE)
flag=(rmd[,4]>=v[1] & rmd[,4]<=v[2])
rmd=rmd[flag,]
}
ys=rs[,5]
plot(rs[,4],rs[,5],ylim=c(0,1),xlab=xlab,ylab=ylab,type='n')
if(LP){
z1=lplot.pred(rmd[,4],rmd[,5],pts=rmd[,2])$yhat
flag=z1>1
z1[flag]=1
flag=z1<0
z1[flag]=0
z2=lplot.pred(rmd[,4],rmd[,6],pts=rmd[,2])$yhat
flag=z2>1
z2[flag]=1
flag=z2<0
z2[flag]=0
z3=lplot.pred(rmd[,4],rmd[,7],pts=rmd[,2])$yhat
flag=z3>1
z3[flag]=1
flag=(z3<0)
z3[flag]=0
lines(rmd[,4],z1)
lines(rmd[,4],z2,lty=2)
lines(rmd[,4],z3,lty=2)
}
if(!LP){
lines(rmd[,4],rmd[,5])
lines(rmd[,4],rmd[,6],lty=2)
lines(rmd[,4],rmd[,7],lty=2)
}}
id=duplicated(rmd[,2:3])
rmd=elimna(rmd[!id,])
rmd
}

#' Compare Regression Estimates vs Running Smoother Estimates
#'
#' Diagnostic plot comparing predicted values from parametric regression
#' against nonparametric running interval smoother. Useful for assessing
#' linearity assumptions.
#'
#' @param x Numeric vector or matrix of predictor values.
#' @param y Numeric vector of response values.
#' @inheritParams common-params
#' @param fr Span for running interval smoother (default: 1).
#' @param est Measure of location for smoother (default: `median`).
#' @param regfun Regression function to use (default: `Qreg` for quantile
#'   regression).
#' @param Qreg.plot Logical; if `TRUE` (default), uses quantile regression
#'   plot. If `FALSE`, uses lowess plot.
#' @param qv Quantiles for quantile regression plot (default: c(0.25, 0.75)).
#' @param SMQ Logical; if `TRUE`, uses smoothed quantile regression
#'   (default: `FALSE`).
#' @param pr Logical; if `TRUE` (default), prints update message.
#' @param xlab Label for x-axis (default: "Reg.Est").
#' @param ylab Label for y-axis (default: "Rplot.Est").
#' @param pch Plotting character (default: '*').
#'
#' @details
#' If the linear model is correct, points should cluster tightly around
#' a line with slope=1 and intercept=0 (shown as dashed reference line).
#' Systematic deviations indicate model misspecification.
#'
#' The function:
#' 1. Computes regression predictions using `regYhat()`
#' 2. Computes smoother predictions using `rplot.pred()`
#' 3. Plots one against the other with reference line y=x
#'
#' @return Invisibly returns `NULL`. Creates diagnostic plot as side effect.
#'
#' @seealso \code{\link{reg.vs.lplot}}, \code{\link{regYhat}},
#'   \code{\link{rplot.pred}}
#'
#' @export
#' @examples
#' # Linear relationship - should follow y=x line
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' reg.vs.rplot(x, y)
#'
#' # Nonlinear relationship - will deviate from y=x
#' y_nl <- x^2 + rnorm(100)
#' reg.vs.rplot(x, y_nl)
reg.vs.rplot<-function(x,y,xout=FALSE,fr=1,est=median,regfun=Qreg,Qreg.plot=TRUE,qv=c(.25,.75),SMQ=FALSE,
pr=TRUE,xlab='Reg.Est',ylab='Rplot.Est',pch='*'){
#
# If the linear model is correct, the plot returned here should be
# tightly clustered around a line having slope=1 and intercept=0, indicated
# by a dashed line.
#
if(pr)print('This function was updated July 2022')
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
x=xy[,1:p]
y=xy[,p1]
e1=regYhat(x,y,xout=xout,regfun=regfun)
e2=rplot.pred(x,y,xout=xout,est=est,fr=fr)$Y.hat
if(Qreg.plot){
if(!SMQ)qplotreg(e1,e2,xlab=xlab,ylab=ylab,pch=pch,qval=qv)
if(SMQ)qhdsm(e1,e2,xlab=xlab,ylab=ylab,pch=pch,qval=qv,LP=TRUE)
}
if(!(Qreg.plot)) lplot(e1,e2,xlab=xlab,ylab=ylab,pc=pch)
abline(0,1,lty=2)
}

#' Compare Regression Estimates vs Lowess Estimates
#'
#' Diagnostic plot comparing predicted values from parametric regression
#' against nonparametric lowess/loess smoother. Useful for assessing
#' linearity assumptions.
#'
#' @param x Numeric vector or matrix of predictor values.
#' @param y Numeric vector of response values.
#' @inheritParams common-params
#' @param Qreg.plot Logical; if `TRUE` (default), uses quantile regression
#'   plot. If `FALSE`, uses lowess plot.
#' @param qv Quantiles for quantile regression plot (default: c(0.25, 0.75)).
#' @param SMQ Logical; if `TRUE`, uses smoothed quantile regression
#'   (default: `FALSE`).
#' @param pch Plotting character (default: '*').
#' @param pr Logical; if `TRUE` (default), prints update message.
#' @param fr Span for running interval in outlier detection (default: 1).
#' @param est Measure of location for smoother (default: `mean`).
#' @param regfun Regression function (default: `tsreg` for Theil-Sen).
#' @param xlab Label for x-axis (default: "Reg.est").
#' @param ylab Label for y-axis (default: "Lplot.est").
#' @param span Span for lowess smoother (default: 0.75).
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @details
#' Similar to `reg.vs.rplot()` but uses lowess/loess smoother instead of
#' running interval smoother. If the linear model is correct, points should
#' cluster tightly around the y=x reference line (shown as dashed).
#'
#' The function:
#' 1. Computes regression predictions using `regYhat()`
#' 2. Computes lowess predictions using `lplot.pred()`
#' 3. Plots one against the other with reference line y=x
#'
#' Systematic deviations from the reference line suggest model misspecification.
#'
#' @return Invisibly returns `NULL`. Creates diagnostic plot as side effect.
#'
#' @seealso \code{\link{reg.vs.rplot}}, \code{\link{regYhat}},
#'   \code{\link{lplot.pred}}
#'
#' @export
#' @examples
#' # Linear relationship
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100)
#' reg.vs.lplot(x, y)
reg.vs.lplot<-function(x,y,xout=FALSE,Qreg.plot=TRUE,qv=c(.25,.75),SMQ=FALSE,pch='*',pr=TRUE,
outfun=outpro,fr=1,est=mean,regfun=tsreg,xlab='Reg.est',ylab='Lplot.est',span=.75,...){
#
#
#
if(pr)print('This function was updated July 2022')
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
x=xy[,1:p]
y=xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
e1=regYhat(x,y,regfun=regfun)
e2=lplot.pred(x,y,,est=est,span=span)$yhat
if(Qreg.plot){
if(!SMQ)qplotreg(e1,e2,xlab=xlab,ylab=ylab,pch=pch,qval=qv)
if(SMQ)qhdsm(e1,e2,xlab=xlab,ylab=ylab,pch=pch,qval=qv,LP=TRUE)
}
if(!(Qreg.plot)) lplot(e1,e2,xlab=xlab,ylab=ylab,pc=pch)
abline(0,1,lty=2)
}

#' 3D Plot of Regression Surface Differences Between Two Groups
#'
#' Fits regression models to two groups with two predictors each, then
#' plots the difference in predicted values as a 3D surface. Useful for
#' visualizing how group differences vary across the predictor space.
#'
#' @param x1 Matrix with 2 columns of predictor values for group 1.
#' @param y1 Numeric vector of response values for group 1.
#' @param x2 Matrix with 2 columns of predictor values for group 2.
#' @param y2 Numeric vector of response values for group 2.
#' @param regfun Regression function (default: `tsreg` for Theil-Sen).
#' @param pts Matrix with 2 columns defining points for prediction.
#'   Default: `x1`.
#' @param xlab Label for x-axis (default: "VAR 1").
#' @param ylab Label for y-axis (default: "VAR 2").
#' @param zlab Label for z-axis (default: "Group 2 minus Group 1").
#' @inheritParams common-params
#' @param ALL Logical; if `TRUE` (default), uses all points from both groups
#'   for prediction (overrides `pts`).
#' @param pts.out Logical; if `TRUE`, removes leverage points from `pts`
#'   using `outfun` (default: `FALSE`).
#' @param SCAT Logical; if `TRUE`, uses `scatterplot3d` package for plotting.
#'   If `FALSE` (default), uses `rplot()`.
#' @param theta Rotation angle for 3D plot (default: 50).
#' @param phi Colatitude angle for 3D plot (default: 25).
#' @param ticktype Tick mark type (default: "simple").
#' @param pr Logical; if `TRUE` (default), prints messages.
#' @param ... Additional arguments passed to `regfun` or `outfun`.
#'
#' @details
#' The function:
#' 1. Fits regression models to both groups
#' 2. Computes predicted values at points in `pts`
#' 3. Computes difference: prediction(group2) - prediction(group1)
#' 4. Creates 3D surface plot of differences
#'
#' Useful for assessing whether group differences are constant across
#' the predictor space (parallel surfaces) or vary (interaction effects).
#'
#' @return Invisibly returns `NULL`. Creates 3D plot as side effect.
#'
#' @seealso \code{\link{reg2plot}}, \code{\link{regYhat}}, \code{\link{rplot}}
#'
#' @export
#' @examples
#' # Two groups with different regression relationships
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- x1[,1] + x1[,2] + rnorm(50)
#' x2 <- matrix(rnorm(100), ncol=2)
#' y2 <- 2*x2[,1] + x2[,2] + rnorm(50)
#' reg2difplot(x1, y1, x2, y2)
reg2difplot<-function(x1,y1,x2,y2,regfun=tsreg,pts=x1,xlab="VAR 1",ylab="VAR 2",zlab="Group 2 minus Group 1",xout=FALSE,outfun=out,ALL=TRUE,pts.out=FALSE,SCAT=FALSE,theta=50,phi=25,ticktype='simple',
pr=TRUE,...){
#
# Fit a regression model to both groups assuming have two predictors.
# Get predicted Y values based on points in pts. By default, use
# pts=x1
#
#  x1 a matrix containing two predictors
#  x2 a matrix containing two predictors
#
#  Compute differences in predicted values and plot the results as a function of the points in pts
#  pts=x1 by default.
#  ALL=T, pts is taken to be all points in x1 and x2.
#
#  pts.out=T will remove leverage points from pts.
#
if(!is.matrix(x1))stop("x1 should be a matrix")
if(!is.matrix(x2))stop("x2 should be a matrix")
if(!is.matrix(pts))stop("pts should be a matrix")
if(ncol(x1)!=2)stop("x1 should be a matrix with two columns")
if(ncol(x2)!=2)stop("x2 should be a matrix with two columns")
if(ncol(pts)!=2)stop("pts should be a matrix with two columns")
if(ALL)pts=rbind(x1,x2)
if(pts.out){
flag=outfun(pts,plotit=FALSE,...)$keep
pts=pts[flag,]
}
e1=regYhat(x1,y1,xout=xout,regfun=regfun,outfun=outfun,xr=pts,...)
e2=regYhat(x2,y2,xout=xout,regfun=regfun,outfun=outfun,xr=pts,...)
if(SCAT){
library(scatterplot3d)
scatterplot3d(cbind(pts,e2-e1),xlab=xlab,ylab=ylab,zlab=zlab)
}
if(!SCAT)rplot(pts,e2-e1,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,pr=FALSE,ticktype=ticktype,prm=FALSE)
}

#' Plot Multiple Quantile Regression Lines
#'
#' Computes and plots quantile regression lines for specified quantiles
#' of the conditional distribution. Allows visualization of how different
#' parts of the response distribution relate to predictors.
#'
#' @param x Numeric vector of predictor values (univariate only).
#' @param y Numeric vector of response values.
#' @param qval Numeric vector of quantiles to plot (default: c(0.2, 0.8)).
#' @param q Alternative specification for `qval` (for backward compatibility).
#' @inheritParams common-params
#' @param plotit Logical; if `TRUE` (default), creates a plot.
#' @param xlab Label for x-axis (default: "X").
#' @param ylab Label for y-axis (default: "Y").
#' @param pch Plotting character (default: '*').
#' @param ... Additional arguments passed to `outfun`.
#'
#' @details
#' For each quantile in `qval`, the function:
#' 1. Fits a quantile regression line using `qreg()`
#' 2. Plots the line on the scatter plot
#'
#' Useful for examining heteroscedasticity and non-constant variance patterns.
#' Only one predictor is allowed.
#'
#' @return Matrix with columns "Inter." (intercept) and "Slope" for each
#'   quantile, with one row per quantile in `qval`.
#'
#' @seealso \code{\link{qregplots}}, \code{\link{qreg}}, \code{\link{Qreg}}
#'
#' @export
#' @examples
#' # Heteroscedastic data
#' x <- rnorm(100)
#' y <- x + (1 + 0.5*x)*rnorm(100)  # variance increases with x
#' qplotreg(x, y, qval=c(0.1, 0.5, 0.9))
qplotreg<-function(x, y,qval=c(.2,.8),q=NULL,plotit=TRUE,xlab="X",ylab="Y",xout=FALSE,
outfun=outpro,pch='*',...){
#
# Compute the quantile regression line for each of the
# quantiles indicated by qval.
# plotit=TRUE, plot the results.
#
if(!is.null(q))qval=q
xy=elimna(cbind(x,y))
if(ncol(xy)>2)stop("Only One Predictor Allowed")
x=xy[,1]
y=xy[,2]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
}
n<-length(qval)
coef<-matrix(NA,ncol=2,nrow=n)
x<-as.matrix(x)
if(ncol(x)>1)stop("This version allows one predictor only.")
if(plotit)plot(x,y,xlab=xlab,ylab=ylab,pch=pch)
for(it in 1:n){
coef[it,]<-qreg(x,y,qval=qval[it],pr=FALSE)$coef
dimnames(coef)=list(NULL,c("Inter.","Slope"))
if(plotit)abline(coef[it,1],coef[it,2])
}
coef
}

#' Plot Quantile Regression Lines Using quantreg Package
#'
#' Computes and plots quantile regression lines using the `quantreg` package.
#' Alternative to `qplotreg()` that uses the Koenker-Bassett approach via
#' `rq()` function.
#'
#' @param x Numeric vector of predictor values (univariate only).
#' @param y Numeric vector of response values.
#' @param qval Numeric vector of quantiles to estimate (default: 0.5 for median).
#' @param q Alternative specification for `qval` (for backward compatibility).
#' @param op Operation code (currently only `op=1` is used).
#' @param pr Logical; if `TRUE`, prints progress messages (default: `FALSE`).
#' @inheritParams common-params
#' @param plotit Logical; if `TRUE`, plots the scatter plot (default: `FALSE`).
#' @param xlab Label for x-axis (default: "X").
#' @param ylab Label for y-axis (default: "Y").
#' @param ... Additional arguments passed to `outfun`.
#'
#' @details
#' Uses the `quantreg` package to fit quantile regression models via the
#' Koenker-Bassett approach. For each quantile in `qval`:
#' 1. Fits quantile regression using `rq()`
#' 2. Plots the regression line
#' 3. Returns coefficient matrix
#'
#' Only one predictor is allowed. This function loads the `quantreg` package.
#'
#' Example: `qregplots(x, y, q=c(.25,.5,.75))` plots regression lines for
#' the quartiles (25th, 50th, 75th percentiles).
#'
#' @return Matrix with columns "q" (quantile), "Inter" (intercept), and
#'   "Slope", with one row per quantile in `qval`.
#'
#' @seealso \code{\link{qplotreg}}, \code{\link{qreg}},
#'   \code{\link[quantreg]{rq}}
#'
#' @export
#' @examples
#' # Plot quartile regression lines
#' x <- rnorm(100)
#' y <- x + rnorm(100)
#' qregplots(x, y, q=c(0.25, 0.5, 0.75), plotit=TRUE)
qregplots<-function(x, y,qval=.5,q=NULL,op=1,pr=FALSE,xout=FALSE,outfun=out,plotit=FALSE,xlab="X",ylab="Y",...){
#
# Compute the quantile regression line for one or more quantiles and plot the results
# That is, determine the qth (qval) quantile of Y given X using the
#  the Koenker-Bassett approach.
#
#  One predictor only is allowed
#
#  v2=T, uses the function rq in the R library quantreg
#  v2=F, uses an older and slower version
#
#  Example: qregplots(x,y,q=c(.25,.5,.75)) will plot the regression lines for
#   predicting quartiles.
#
if(!is.null(q))qval=q
x<-as.matrix(x)
if(ncol(x)!=1)stop("Current version allows only one predictor")
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
p<-np-1
x<-X[,1:p]
x<-as.matrix(x)
y<-X[,np]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
est=matrix(NA,ncol=3,nrow=length(qval))
dimnames(est)=list(NULL,c("q","Inter","Slope"))
library(quantreg)
x<-as.matrix(x)
plot(x,y,xlab=xlab,ylab=ylab)
if(ncol(x)!=1)stop("Current version allows only one predictor")
for(j in 1:length(qval)){
coef=coefficients((rq(y~x,tau=qval[j])))
est[j,1]=qval[j]
est[j,2:3]=coef
abline(coef)
}
list(coef = est)
}

#' Partial Residual Plot for Checking Curvature
#'
#' @description
#' Creates a partial residual plot to check for curvature (nonlinearity) associated
#' with a specific predictor variable in a multiple regression context. Subtracts out
#' the linear predictions from other predictors and plots a smooth against the focal predictor.
#'
#' @param x Matrix of predictor variables (must have 2+ columns).
#' @param y Response variable vector.
#' @param pval Column index of the predictor to check for curvature (default: last column).
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen regression).
#' @param fr Span parameter for the running interval smoother (default: 0.8).
#' @param est Location estimator for running smoother (default: \code{tmean}).
#' @param op Operation mode: 1 for LOWESS (\code{lplot}), other values for running interval
#'   smoother (\code{rungen}) (default: 1).
#' @param xlab Label for X-axis (default: "X").
#' @param ylab Label for Y-axis (default: "Residuals").
#' @param xout Logical; if TRUE, removes outliers based on predictor variables (default: FALSE).
#' @param outfun Outlier detection function (default: \code{out}).
#' @param ... Additional arguments passed to the outlier detection function.
#'
#' @return Invisibly returns NULL. Creates a plot as a side effect.
#'
#' @details
#' A **partial residual plot** is useful for detecting curvature in the relationship
#' between a specific predictor and the response after accounting for other predictors.
#'
#' The function:
#' 1. Fits a regression using all predictors EXCEPT the one indicated by `pval`
#' 2. Computes residuals from this reduced model
#' 3. Plots these residuals against the excluded predictor with a smooth curve
#'
#' If there is curvature (nonlinearity) in the smooth, this suggests that the linear
#' relationship may not be adequate for that predictor.
#'
#' The smoother used depends on the `op` parameter:
#' - `op=1`: Uses LOWESS via \code{lplot}
#' - Otherwise: Uses a running interval smoother via \code{rungen}
#'
#' @seealso \code{\link{lplot}}, \code{\link{rungen}}, \code{\link{tsreg}}
#'
#' @export
#' @examples
#' # Check for curvature in second predictor
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- 2*x1 + 0.5*x2^2 + rnorm(100)
#' X <- cbind(x1, x2)
#' prplot(X, y, pval = 2)  # Should show curvature for x2
prplot<-function(x,y,pval=ncol(x),regfun=tsreg,fr=.8,est=tmean,op=1,
xlab="X",ylab="Residuals",xout=FALSE,outfun=out,...){
#
# Goal: check for curvature associated with predictor
# indicated by pval.
# This is done by creating a partial residual plot.
# That is subtracting out the linear prediction based
# on the other predictors and then
# smooth the result versus the predictor in the column of x indicated by pval
#
x=as.matrix(x)
p=ncol(x)
p1=p+1
temp=elimna(cbind(x,y))
x=temp[,1:p]
y=temp[,p1]
if(xout){
flag<-outfun(x,...)$keep
x<-as.matrix(x)
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(!is.matrix(x))stop("Should have two or more variables stored in a matrix")
flag<-rep(TRUE,ncol(x))
flag[pval]<-FALSE
temp<-regfun(x[,flag],y)$residual
if(op!=1)rungen(x[,!flag],temp,est=est,fr=fr,xlab=xlab,ylab=ylab,...)
if(op==1)lplot(x[,!flag],temp,xlab=xlab,ylab=ylab)
}

#' Forward Response and Residual Plots for Multiple Linear Regression
#'
#' @description
#' Creates diagnostic plots for multiple linear regression: a forward response plot
#' (fitted vs actual values) and a residual plot (fitted vs residuals). Identifies
#' influential observations based on Cook's distance.
#'
#' @param x Matrix of predictor variables.
#' @param Y Response variable vector.
#'
#' @return Invisibly returns NULL. Creates two diagnostic plots as side effects:
#'   1. **Forward Response Plot**: Fitted values vs actual values with y=x reference line.
#'      Influential points (high Cook's distance) are shown as filled squares.
#'   2. **Residual Plot**: Fitted values vs residuals.
#'      Influential points are highlighted.
#'
#' @details
#' This function uses ordinary least squares (OLS) regression via \code{lsfit}.
#' Influential observations are identified using **Cook's distance**, with the threshold
#' set to min(0.5, 2p/n) where p is the number of predictors + 1 and n is sample size.
#'
#' The function creates a 2-panel plot layout:
#' - **Top panel**: Forward response plot shows how well fitted values match observed values.
#'   Points should cluster around the y=x line for a good fit.
#' - **Bottom panel**: Residual plot shows residuals vs fitted values.
#'   Should show random scatter with no pattern if model assumptions hold.
#'
#' Influential points (Cook's distance > threshold) are plotted as filled squares (pch=15).
#' The \code{identify} function is called on both plots, allowing interactive identification
#' of points (click to label, right-click to exit).
#'
#' **Note**: This function uses classical OLS, not robust regression. For robust alternatives,
#' see \code{\link{rplot.res}}.
#'
#' @seealso \code{\link{rplot.res}}, \code{\link{regplot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Multiple regression diagnostic plots
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' Y <- 2 + 3*x1 - 2*x2 + rnorm(100)
#' X <- cbind(x1, x2)
#' MLRplot(X, Y)
#' # Click on points to identify them, right-click to finish
#' }
MLRplot<-function(x, Y)
{
# Forward response plot and residual plot.
# R needs command "library(lqs)" if a robust estimator replaces lsfit.
# Advance the view with the right mouse button.
	x <- as.matrix(x)
	out <- lsfit(x, Y)
	cook <- ls.diag(out)$cooks
	n <- dim(x)[1]
	p <- dim(x)[2] + 1
	tem <- cook > min(0.5, (2 * p)/n)
	bhat <- out$coef
	FIT <- bhat[1] + x %*% bhat[-1]
	par(mfrow = c(2, 1))
	plot(FIT, Y)
	abline(0, 1)
	points(FIT[tem], Y[tem], pch = 15)
	identify(FIT, Y)
	title("Forward Response Plot")
	RES <- Y - FIT
	plot(FIT, RES)
	points(FIT[tem], RES[tem], pch = 15)
	identify(FIT, RES)
	title("Residual Plot")
}


# ============================================================================
# LOWESS/LOESS PLOTS
# ============================================================================

#' Plot Regression Surface Using LOESS (Version 2)
#'
#' @description
#' Creates a regression surface plot using LOESS (locally weighted scatterplot smoothing)
#' with options for outlier removal, strength of association estimation, and 3D visualization
#' for multiple predictors. This version includes robust strength and explanatory power measures.
#'
#' @param x Predictor variable(s). Can be a vector (single predictor) or matrix (multiple predictors).
#'   Maximum of 4 predictors allowed.
#' @param y Response variable (vector).
#' @param span Span parameter for LOESS when using 2+ predictors (default: 0.75).
#'   Controls the degree of smoothing.
#' @param pyhat Logical; if TRUE, returns fitted Y values (default: FALSE).
#' @param eout Logical; if TRUE, eliminates multivariate outliers before fitting (default: FALSE).
#' @param xout Logical; if TRUE, eliminates outliers based on predictor variables only (default: FALSE).
#' @param outfun Outlier detection function (default: \code{out}).
#' @param plotit Logical; if TRUE, creates the plot (default: TRUE).
#' @param expand Expansion factor for 3D plot (default: 0.5).
#' @param low.span Span parameter for LOWESS when using 1 predictor (default: 2/3).
#' @param varfun Variance function for computing explanatory power (default: \code{pbvar}).
#' @param cor.op Logical; if TRUE, uses correlation-based explanatory power (default: FALSE).
#' @param cor.fun Correlation function when cor.op=TRUE (default: \code{pbcor}).
#' @param ADJ Logical; if TRUE, computes adjusted strength/explanatory measures (default: FALSE).
#' @param nboot Number of bootstrap samples for adjusted measures (default: 20).
#' @param scale Logical; if TRUE, scales the 3D plot (default: TRUE).
#' @param xlab Label for X-axis (default: "X").
#' @param ylab Label for Y-axis (default: "Y").
#' @param zlab Label for Z-axis in 3D plots (default: "").
#' @param theta Rotation angle for 3D plot (default: 50).
#' @param phi Colatitude angle for 3D plot (default: 25).
#' @param family Family parameter for LOESS: "gaussian" or "symmetric" (default: "gaussian").
#' @param duplicate How to handle duplicate x values in 3D plots: "error" or "strip" (default: "error").
#' @param pr Logical; if TRUE, prints warnings (default: TRUE).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param ticktype Type of tick marks for 3D plot (default: "simple").
#'
#' @return A list with components:
#'   \item{Strength.Assoc}{Strength of association measure (signed square root of explanatory power).}
#'   \item{Explanatory.power}{Proportion of variance explained by the smoother.}
#'   \item{Strength.Adj}{Adjusted strength of association (if ADJ=TRUE).}
#'   \item{Explanatory.Adj}{Adjusted explanatory power (if ADJ=TRUE).}
#'   \item{yhat.values}{Fitted values (if pyhat=TRUE).}
#'
#' @details
#' This function fits a LOESS regression surface and provides measures of association strength.
#'
#' **For 1 predictor**: Uses LOWESS with the specified low.span parameter.
#'
#' **For 2-4 predictors**: Uses LOESS with interaction terms and creates a 3D surface plot
#' (for 2 predictors only).
#'
#' The **strength of association** is computed as the signed square root of the explanatory
#' power, where the sign indicates the direction of association (increasing vs decreasing).
#'
#' **Explanatory power** is computed as the ratio of variances or as squared correlation,
#' depending on the cor.op parameter.
#'
#' When ADJ=TRUE, adjusted measures account for overfitting by comparing to bootstrap
#' samples where x and y are independent.
#'
#' **Note on duplicates**: When plotting 3D surfaces with 2 predictors, duplicate x values
#' may cause errors. Set duplicate="strip" to remove duplicates.
#'
#' @seealso \code{\link{lplot}}, \code{\link{lplotCI}}, \code{\link{lplotPV}}
#' @export
#' @examples
#' # Single predictor with LOWESS
#' x <- rnorm(100)
#' y <- x + 0.5 * x^2 + rnorm(100, sd = 0.5)
#' lplotv2(x, y, low.span = 0.5)
#'
#' # Two predictors with LOESS 3D surface
#' \dontrun{
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- x1 + x2 + 0.5 * x1 * x2 + rnorm(100)
#' X <- cbind(x1, x2)
#' lplotv2(X, y, span = 0.6)
#' }
lplotv2<-function(x,y,span=.75,pyhat=FALSE,eout=FALSE,xout=FALSE,outfun=out,plotit=TRUE,
expand=.5,low.span=2/3,varfun=pbvar,cor.op=FALSE,cor.fun=pbcor,ADJ=FALSE,nboot=20,
scale=TRUE,xlab="X",ylab="Y",zlab="",theta=50,phi=25,family="gaussian",
duplicate="error",pr=TRUE,SEED=TRUE,ticktype="simple"){
#
# Plot regression surface using LOESS
#
# low.span is the span when lowess is used and there is one predictor
# span is the span when loess is used with two or more predictors
# pyhat=T will return Y hat values
# eout=T will eliminate outliers
# xout=T  will eliminate points where X is an outliers
# family="gaussian"; see the description of the built-in function loess
#
# duplicate="error"
# In some situations where duplicate values occur, when plotting with
# two predictors, it is necessary to set duplicate="strip"
#
st.adj=NULL
e.adj=NULL
if(ADJ){
if(SEED)set.seed(2)
}
si=1
x<-as.matrix(x)
if(!is.matrix(x))stop("x is not a matrix")
d<-ncol(x)
if(d>=2){
if(ncol(x)==2 && !scale){
if(pr){
print("scale=F is specified.")
print("If there is dependence, might use scale=T")
}}
m<-elimna(cbind(x,y))
x<-m[,1:d]
y<-m[,d+1]
if(eout && xout)stop("Can't have both eout and xout = F")
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
}
x<-m[,1:d]
y<-m[,d+1]
if(d==2)fitr<-fitted(loess(y~x[,1]*x[,2],span=span,family=family))
if(d==3)fitr<-fitted(loess(y~x[,1]*x[,2]*x[,3],span=span,family=family))
if(d==4)fitr<-fitted(loess(y~x[,1]*x[,2]*x[,3]*x[,4],span=span,family=family))
if(d>4)stop("Can have at most four predictors")
last<-fitr
if(d==2 && plotit){
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
fitr<-interp(mkeep[,1],mkeep[,2],fitr,duplicate=duplicate)
persp(fitr,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,expand=expand,
scale=scale,ticktype=ticktype)
}}
if(d==1){
m<-elimna(cbind(x,y))
x<-m[,1:d]
y<-m[,d+1]
if(eout && xout)stop("Can't have both eout and xout = F")
if(eout){
flag<-outfun(m)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(x)$keep
m<-m[flag,]
}
x<-m[,1:d]
y<-m[,d+1]
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab)
lines(lowess(x,y,f=low.span))
}
yyy<-lowess(x,y)$y
xxx<-lowess(x,y)$x
if(d==1){
ordx=order(xxx)
yord=yyy[ordx]
flag=NA
for (i in 2:length(yyy))flag[i-1]=sign(yord[i]-yord[i-1])
if(sum(flag)<0)si=-1
}
last<-yyy
chkit<-sum(duplicated(x))
if(chkit>0){
last<-rep(1,length(y))
for(j in 1:length(yyy)){
for(i in 1:length(y)){
if(x[i]==xxx[j])last[i]<-yyy[j]
}}
}
}
E.power<-1
if(!cor.op)E.power<-varfun(last[!is.na(last)])/varfun(y)
if(cor.op || E.power>=1){
if(d==1){
xord<-order(x)
E.power<-cor.fun(last,y[xord])$cor^2
}
if(d>1)E.power<-cor.fun(last,y)$cor^2
}
if(ADJ){
x=as.matrix(x)
val=NA
n=length(y)
data1<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(i in 1:nboot){
temp=lplot.sub(x[data1[i,],],y[data2[i,]],plotit=FALSE,pr=FALSE)
val[i]=temp$Explanatory.power
}
vindt=median(val)
v2indt=median(sqrt(val))
st.adj=(sqrt(E.power)-max(c(0,v2indt)))/(1-max(c(0,v2indt)))
e.adj=(E.power-max(c(0,vindt)))/(1-max(c(0,vindt)))
st.adj=max(c(0,st.adj))
e.adj=max(c(0,e.adj))
}
if(!pyhat)last <- NULL
list(Strength.Assoc=si*sqrt(E.power),Explanatory.power=E.power,
Strength.Adj=st.adj,Explanatory.Adj=e.adj,yhat.values=last)
}

#' Plot Running Interval Smoother for Two Groups
#'
#' @description
#' Creates a scatterplot with running interval smoothers for two independent groups,
#' allowing visual comparison of regression patterns between groups. Uses LOWESS smoothing.
#'
#' @param x1 Predictor values for group 1 (vector).
#' @param y1 Response values for group 1 (vector).
#' @param x2 Predictor values for group 2 (vector).
#' @param y2 Response values for group 2 (vector).
#' @param fr Span for LOWESS smoother (default: 0.8). Controls amount of smoothing.
#' @param est Not used in current implementation (retained for compatibility).
#' @param xlab Label for X-axis (default: "X").
#' @param ylab Label for Y-axis (default: "Y").
#' @param xout Logical; if TRUE, removes outliers based on X values only (default: FALSE).
#' @param eout Logical; if TRUE, removes bivariate outliers (default: FALSE).
#' @param pr Logical; if TRUE, prints suggestion about using xout=TRUE (default: TRUE).
#' @param pch1 Plotting character for group 1 points (default: '*').
#' @param pch2 Plotting character for group 2 points (default: 'o').
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param ... Additional arguments passed to outfun.
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function creates a scatterplot showing data from two groups with separate
#' LOWESS smoothed curves overlaid. The curves allow visual comparison of the
#' regression relationships in the two groups.
#'
#' **Group 1** is plotted with solid line and pch1 symbol.
#' **Group 2** is plotted with dashed line (lty=2) and pch2 symbol.
#'
#' Missing values are automatically removed. Outlier removal (if requested) is
#' performed separately for each group before smoothing.
#'
#' **Note**: The function suggests also examining the plot with xout=TRUE to assess
#' the impact of outliers in the predictor space.
#'
#' @seealso \code{\link{lplot}}, \code{\link{rplot2g}}, \code{\link{reg2plot}}
#' @export
#' @examples
#' # Compare regression patterns in two groups
#' set.seed(123)
#' x1 <- rnorm(50)
#' y1 <- 2 * x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- -1.5 * x2 + rnorm(50)
#' lplot2g(x1, y1, x2, y2)
#'
#' # With outlier removal
#' lplot2g(x1, y1, x2, y2, xout = TRUE)
lplot2g<-function(x1,y1,x2,y2,fr=.8,est=tmean,xlab="X",ylab="Y",xout=FALSE,eout=FALSE,pr=TRUE,
pch1='*',pch2='o',outfun=outpro,...){
#
# Plot of running interval smoother for two groups
#
# fr controls amount of smoothing
# tr is the amount of trimming
#
# Missing values are automatically removed.
#
# sm=T results in using bootstrap bagging when estimating the regression line
# nboot controls number of bootstrap samples
#
if(pr){
if(!xout)print('Suggest also looking at result using xout=TRUE')
}
m1<-elimna(cbind(x1,y1))
if(eout && xout)stop("Can't have both eout and xout = F")
if(eout){
flag<-outfun(m1,plotit=FALSE,...)$keep
m1<-m1[flag,]
}
if(xout){
flag<-outfun(m1[,1],plotit=FALSE,...)$keep
m1<-m1[flag,]
}
x1<-m1[,1]
y1<-m1[,2]
m2<-elimna(cbind(x2,y2))
if(eout){
flag<-outfun(m2,plotit=FALSE,...)$keep
m2<-m2[flag,]
}
if(xout){
flag<-outfun(m2[,1],plotit=FALSE,...)$keep
m2<-m2[flag,]
}
x2<-m2[,1]
y2<-m2[,2]

flag=order(x1)
x1=x1[flag]
y1=y1[flag]
flag=order(x2)
x2=x2[flag]
y2=y2[flag]
temp1<-lplot(x1,y1,pyhat=TRUE,plotit=FALSE,pr=FALSE)$yhat.values
temp2<-lplot(x2,y2,pyhat=TRUE,plotit=FALSE,pr=FALSE)$yhat.values
plot(c(x1,x2),c(y1,y2),type="n",xlab=xlab,ylab=ylab)
points(x1,y1,pch=pch1)
points(x2,y2,pch=pch2)
lines(x1,temp1)
lines(x2,temp2,lty=2)
}

#' LOESS Confidence Band with Simultaneous Coverage
#'
#' @description
#' Creates a LOESS regression plot with simultaneous confidence band that accounts
#' for heteroscedasticity. The confidence intervals are adjusted so that the overall
#' probability coverage is approximately 1-alpha across all points.
#'
#' @param x Predictor variable (vector).
#' @param y Response variable (vector).
#' @param plotit Logical; if TRUE, creates the plot (default: TRUE).
#' @param xlab Label for X-axis (default: 'X').
#' @param ylab Label for Y-axis (default: 'Y').
#' @param p.crit Critical p-value for simultaneous confidence band (default: NULL, auto-determined).
#' @param alpha Familywise Type I error rate (default: 0.05).
#' @param span Span parameter for LOESS (default: NULL, auto-selected based on n).
#' @param CIV Logical; if TRUE, returns confidence interval values instead of plotting (default: FALSE).
#' @param xout Logical; if TRUE, removes outliers based on X values (default: FALSE).
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param pch Plotting character for points (default: '.').
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param nboot Number of bootstrap samples for SE estimation (default: 100).
#' @param pts Points at which to compute confidence intervals (default: NULL, auto-selected).
#' @param npts Number of points for confidence band if pts=NULL (default: 25).
#' @param nreps Number of replications for computing p.crit if needed (default: 2000).
#' @param ... Additional arguments passed to outfun.
#'
#' @return A list with components:
#'   \item{p.crit}{Adjusted critical p-value used for simultaneous coverage.}
#'   \item{Conf.Intervals}{Matrix with columns (X, Y.hat, ci.low, ci.up) if CIV=TRUE, else NULL.}
#'
#' @details
#' This function fits a LOESS curve and constructs a simultaneous confidence band
#' that maintains approximately 1-alpha coverage probability across all points.
#'
#' **Span selection**: If span=NULL, it is chosen based on sample size:
#' - n < 300: span = 2/3
#' - 300 ≤ n < 800: span = 0.5
#' - n ≥ 800: span = 0.3
#'
#' **Critical value**: For alpha=0.05 and n ≤ 2000, p.crit is determined quickly
#' using stored simulation results. Otherwise, it is computed via simulation
#' (which may take some time).
#'
#' The confidence band is computed at npts points ranging from median(x) - 1.5*MAD(x)
#' to median(x) + 1.5*MAD(x), or at user-specified pts.
#'
#' Standard errors are estimated via bootstrap, and the confidence intervals
#' account for heteroscedasticity.
#'
#' @seealso \code{\link{lplot}}, \code{\link{lplotv2}}, \code{\link{rplotCI}}
#' @export
#' @examples
#' # LOESS with simultaneous confidence band
#' set.seed(123)
#' x <- rnorm(100)
#' y <- x + 0.3 * x^2 + rnorm(100, sd = 0.5)
#' lplotCI(x, y, span = 0.6)
#'
#' # Return confidence intervals
#' result <- lplotCI(x, y, CIV = TRUE, plotit = FALSE)
#' head(result$Conf.Intervals)
lplotCI<-function(x,y,plotit=TRUE,xlab='X',ylab='Y',p.crit=NULL,alpha=.05,span=NULL,
CIV=FALSE,xout=FALSE,outfun=outpro, pch='.',SEED=TRUE,nboot=100,pts=NULL,npts=25,nreps=2000,...){
#
#  Confidence band using LOESS
#
#  Method allows heteroscedasticity and adjust the confidence intervals
# so that the simultaneous probabillty coverage is approximately 1-alpha
#
# If CIV=FALSE and plotit=TRUE, creates a plot with the confidence intervals.
# CIV=TRUE, returns the confidence intervals for the points in pts
# pts =NULL, the function picks
# npts points, extending between M-1.5*mad(x) and M+1.5*mad(x)
#
#
# For alpha=0.05, n <=2000 execution time is low. Otherwise
# the adjusted critical value must be computed.
#
#  p.crit=NULL: If alpha=.05, determined quickly, otherwise it is computed.
#
xy=elimna(cbind(x,y))
if(ncol(xy)>2)stop('Current version limited to a single predictor variable')
if(xout){
flag<-outfun(xy[,1],plotit=FALSE,...)$keep
xy<-xy[flag,]
}
n=nrow(xy)
if(is.null(span)){
span=2/3
if(n >=300)span=.5
if(n >=800)span=.3
}
x=xy[,1]
y=xy[,2]
xord=order(x)
y=y[xord]
x=x[xord]
M=median(x)
low=M-1.5*mad(x)
up=M+1.5*mad(x)
if(is.null(pts))pts=seq(low,up,length.out=npts)
if(npts<=5)p.crit=alpha/npts
if(alpha==.05){
if(is.null(p.crit)){
if(n<30)stop('Should have n>=30')
nv=c(30,50,70,100,150, 200,300, 500, 1000, 2000)
pv=c(0.003599898, 0.002661925, 0.002399994, 0.002877103, 0.003000428, 0.003538190,
 0.003872710, 0.004396500, 0.004075000, 0.0045161)

if(npts==25){
if(n<=2000)p.crit=lplot.pred(1/nv,pv,pts=1/n)$yhat
if(n>2000)p.crit=.00452
}}}
if(is.null(p.crit)){
print('p.crit is being computed, this might take some time.')
pts.stand=NULL
if(!is.null(pts))pts.stand=(median(x)-pts)/mad(x)
p.crit=lplotbsepvv3(n,nreps=nreps,npts=npts,pts=pts.stand,alpha=alpha)
}
plx<-predict(loess(y ~ x,span=span), se=TRUE)
se=lplotse(x,y,nboot=nboot,SEED=SEED,pts=pts,span=span)
lfit=lplot.pred(x,y,pts=pts,span=2/3)$yhat
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab,pch=pch)
lines(x,plx$fit)
if(is.null(p.crit))p.crit=alpha
lines(pts,lfit - qt(1-p.crit/2,plx$df)*se, lty=2)
lines(pts,lfit + qt(1-p.crit/2,plx$df)*se, lty=2)
}
ci.low=lfit - qt(1-p.crit/2,plx$df)*se
ci.up=lfit + qt(1-p.crit/2,plx$df)*se
if(!CIV)ci=NULL
if(CIV){
ci=cbind(pts,lfit,ci.low,ci.up)
dimnames(ci)=list(NULL,c('X','Y.hat','ci.low','ci.up'))
}
list(p.crit=p.crit,Conf.Intervals=ci)
}

#' Bootstrap Standard Errors for LOESS Predictions
#'
#' @description
#' Computes bootstrap standard errors for LOESS (locally weighted scatterplot smoothing)
#' predictions at specified points.
#'
#' @param x Predictor variable vector.
#' @param y Response variable vector.
#' @param pts Points at which to estimate standard errors (default: uses x values).
#' @param nboot Number of bootstrap samples (default: 100).
#' @param SEED Logical; if TRUE, sets random seed to 2 for reproducibility (default: TRUE).
#' @param span Span parameter for LOESS (default: 2/3). Controls degree of smoothing.
#'
#' @return Vector of standard errors corresponding to the prediction points in `pts`.
#'
#' @details
#' This is an internal helper function used by other plotting functions to estimate
#' standard errors for LOESS-based smoothers.
#'
#' The function:
#' 1. Orders data by x values
#' 2. Generates bootstrap samples by resampling observations
#' 3. Fits LOESS to each bootstrap sample
#' 4. Computes predictions at specified points
#' 5. Returns the standard deviation of predictions across bootstrap samples
#'
#' @seealso \code{\link{lplotCI}}, \code{\link{lplot.pred}}
#'
#' @keywords internal
#' @examples
#' x <- rnorm(50)
#' y <- x + rnorm(50)
#' se_vals <- lplotse(x, y, nboot = 50)
lplotse<-function(x,y,pts=x,nboot=100,SEED=TRUE,span=2/3){
#
# compute estimae of SE
# return the values corresponding to the order x values
#
xord=order(x)
y=y[xord]
x=x[xord]
if(SEED)set.seed(2)
n=length(y)
ev=matrix(NA,nrow=nboot,ncol=length(pts))
for(i in 1:nboot){
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
ev[i,]=lplot.pred(x[data[i,]],y[data[i,]],pts=pts,span=span)$yhat
}
se=apply(ev,2,sd)
se
}

#' P-Value for LOESS Strength of Association
#'
#' @description
#' Tests for dependence between x and y using the strength of association measure
#' from LOESS regression. Computes a p-value via bootstrap permutation test.
#'
#' @param x Predictor variable(s). Can be a vector (single predictor) or matrix (multiple predictors).
#' @param y Response variable vector.
#' @param span Span parameter for LOESS with 2+ predictors (default: 0.75).
#' @param xout Logical; if TRUE, removes outliers based on predictor variables (default: FALSE).
#' @param pr Logical; if TRUE, prints warnings (default: TRUE).
#' @param outfun Outlier detection function (default: \code{out}).
#' @param nboot Number of bootstrap permutation samples (default: 1000).
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE).
#' @param plotit Logical; if TRUE, creates the LOESS plot (default: TRUE).
#' @param pyhat Logical; if TRUE, returns fitted values (default: FALSE).
#' @param expand Expansion factor for 3D plots (default: 0.5).
#' @param low.span Span for LOWESS with 1 predictor (default: 2/3).
#' @param varfun Variance function for explanatory power (default: \code{pbvar}).
#' @param cor.op Logical; if TRUE, uses correlation for explanatory power (default: FALSE).
#' @param cor.fun Correlation function when cor.op=TRUE (default: \code{pbcor}).
#' @param scale Logical; if TRUE, scales 3D plots (default: FALSE).
#' @param xlab Label for X-axis (default: "X").
#' @param ylab Label for Y-axis (default: "Y").
#' @param zlab Label for Z-axis in 3D plots (default: "").
#' @param theta Rotation angle for 3D plots (default: 50).
#' @param phi Colatitude angle for 3D plots (default: 25).
#' @param family LOESS family: "gaussian" or "symmetric" (default: "gaussian").
#' @param duplicate How to handle duplicate x values: "error" or "strip" (default: "error").
#' @param pc Plotting character for scatter points (default: "*").
#' @param ticktype Type of tick marks for 3D plots (default: "simple").
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @return A list with components:
#'   \item{p.value}{Bootstrap permutation p-value for testing independence.}
#'   \item{Strength.Assoc}{Observed strength of association.}
#'   \item{Explanatory.power}{Observed explanatory power.}
#'   \item{yhat.values}{Fitted values (if pyhat=TRUE).}
#'
#' @details
#' This function tests the null hypothesis of **independence** between x and y.
#'
#' The test statistic is the **strength of association** computed via LOESS:
#' \deqn{Strength = sign \times \sqrt{ExplanatoryPower}}
#'
#' The p-value is computed by:
#' 1. Computing the observed strength of association
#' 2. Generating `nboot` permutation samples where x and y are independent
#'    (achieved by bootstrapping x and y separately)
#' 3. Computing strength for each permutation sample
#' 4. p-value = proportion of permutation strengths >= observed strength
#'
#' A **small p-value** indicates significant evidence of dependence.
#'
#' The function uses \code{\link{lplot}} internally for fitting and visualization.
#'
#' @seealso \code{\link{lplot}}, \code{\link{lplotv2}}, \code{\link{lplotCI}}
#'
#' @export
#' @examples
#' # Test for dependence
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100, sd = 0.5)
#' result <- lplotPV(x, y, nboot = 500)
#' print(result$p.value)  # Should be small (significant dependence)
lplotPV<-function(x,y, span = 0.75, xout = FALSE,pr=TRUE,
    outfun = out,nboot=1000,SEED=TRUE,plotit=TRUE,pyhat = FALSE, expand = 0.5, low.span = 2/3,
    varfun = pbvar, cor.op = FALSE, cor.fun = pbcor, scale = FALSE,
    xlab = "X", ylab = "Y", zlab = "", theta = 50, phi = 25,
    family = "gaussian", duplicate = "error", pc = "*", ticktype = "simple",...){
#
# Compute a p-value based on the Strength of Association estimated via lplot
# If significant, conclude there is dependence.
#
if(SEED)set.seed(2)
x=as.matrix(x)
if(ncol(x)==2 && !scale){
if(pr){
print("scale=F is specified.")
print("If there is dependence, might use scale=T")
}}
vals=NA
nv=ncol(x)
m=elimna(cbind(x,y))
x<-m[,1:nv]
y<-m[,nv+1]
if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:nv]
y<-m[,nv+1]
}
x=as.matrix(x)
est=lplot(x,y,span=span,plotit=plotit,pr=FALSE, pyhat = pyhat,
    outfun = outfun, expand = expand, low.span = low.span,
    varfun = varfun, cor.op =cor.op, cor.fun = cor.fun, scale = scale,
    xlab = xlab, ylab = ylab, zlab =zlab, theta =theta, phi = phi,
    family = family, duplicate = duplicate, pc = pc, ticktype = ticktype,...)
n=nrow(x)
data1<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(i in 1:nboot){
vals[i]=lplot(x[data1[i,],],y[data2[i,]],plotit=FALSE,pr=FALSE)$Strength.Assoc
}
p=mean(est$Strength<vals)
list(p.value=p,Strength.Assoc=est$Strength.Assoc,Explanatory.power=est$Explanatory.power,yhat.values=est$yhat.values)
}

#' Compare Predictor Importance with LOESS (Method 1)
#'
#' @description
#' For two predictors, estimates their relative importance when using LOESS regression.
#' Compares the conditional variability of predictions when one predictor is held fixed.
#'
#' @param x Matrix with 2 columns (two predictor variables).
#' @param y Response variable vector.
#' @param xout Logical; if TRUE, removes multivariate outliers from predictors (default: FALSE).
#' @param pts1 Points at which to evaluate predictor 1 (default: computed from data).
#' @param pts2 Points at which to evaluate predictor 2 (default: computed from data).
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param span Span parameter for LOESS (default: 2/3).
#' @param npts Number of evaluation points if pts1/pts2 not specified (default: 10).
#' @param tr Trimming proportion for Winsorized SD (default: 0.2).
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @return A list with components:
#'   \item{str1}{Strength of predictor 1 at each point in pts2 (conditional SD ratio).}
#'   \item{str2}{Strength of predictor 2 at each point in pts1 (conditional SD ratio).}
#'   \item{p}{Proportion of pairwise comparisons where str1 > str2.}
#'   \item{mean.str1}{Mean strength for predictor 1.}
#'   \item{mean.str2}{Mean strength for predictor 2.}
#'
#' @details
#' This function assesses **predictor importance** by examining conditional variability.
#'
#' For each predictor:
#' 1. Hold the other predictor fixed at several values
#' 2. Compute LOESS predictions varying only the focal predictor
#' 3. Measure the standard deviation of these predictions
#' 4. Normalize by the overall SD of y
#'
#' **Strength** = (SD of conditional predictions) / (SD of y)
#'
#' Higher strength indicates greater importance. The comparison p-value indicates
#' the proportion of comparisons where predictor 1 shows greater strength than predictor 2.
#'
#' If pts1/pts2 are not specified, evaluation points are chosen as:
#' \deqn{Median \pm 1.5 \times MAD}
#'
#' Uses **Winsorized SD** for robust dispersion estimation.
#'
#' @seealso \code{\link{lplotcom2v2}}, \code{\link{lplotcomBCI}}, \code{\link{lplot}}
#'
#' @export
#' @examples
#' # Compare importance of two predictors
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- 2*x1 + 0.5*x2 + rnorm(100)  # x1 more important
#' X <- cbind(x1, x2)
#' result <- lplotcom2(X, y, npts = 5)
#' result$mean.str1  # Strength for x1
#' result$mean.str2  # Strength for x2
#' result$p  # Proportion where x1 > x2
lplotcom2<-function(x,y,xout=FALSE,pts1=NULL,pts2=NULL,outfun=outpro,span=2/3,npts=10,tr=.2,...){
#
# For two independent variables, estimate their relative importance when using LOESS
#
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n.orig=nrow(m)
n.keep=n.orig
d<-ncol(x)
if(d!=2)stop('Current version is for two independent variables only')
if(xout){
flag<-outfun(m[,1:2],plotit=FALSE,...)$keep
m<-m[flag,]
}
n.keep=nrow(m)
M=apply(m,2,median)
SD=apply(m,2,mad)
low=M-1.5*SD
up=M+1.5*SD
if(is.null(pts1))pts1=seq(low[1],up[1],length.out=npts)
if(is.null(pts2))pts2=seq(low[2],up[2],length.out=npts)
e1=NA  #
e2=NA
for(j in 1:length(pts1)){     # Determine strength of x2 given a value stored in pts1.
v2=cbind(rep(pts1[j],n.keep),m[,2])
vals=lplot.pred(m[,1:2],m[,3],v2,span=span)$yhat
vals=elimna(vals)
nv=length(vals)
e2[j]=NA
if(nv>=10)e2[j]=winsd(vals,tr=tr,na.rm=TRUE)/winsd(m[,3],tr=tr,na.rm=TRUE)
}
for(j in 1:length(pts2)){     # Determine strength of x1 given a value stored in pts2.
v1=cbind(m[,1],rep(pts2[j],n.keep))
vals=lplot.pred(m[,1:2],m[,3],v1,span=span)$yhat
vals=elimna(vals)
nv=length(vals)
e1[j]=NA
if(nv>=10)e1[j]=winsd(vals,tr=tr,na.rm=TRUE)/winsd(m[,3],tr=tr,na.rm=TRUE)
 }
p=mean(outer(e1,e2,FUN='-')<0,na.rm=TRUE)
list(str1=e1,str2=e2,p=p,mean.str1=mean(e1),mean.str2=mean(e2))
}

#' Compare Predictor Importance with LOESS (Method 2)
#'
#' @description
#' Alternative version of \code{\link{lplotcom2}} for comparing predictor importance using LOESS.
#' Uses absolute conditional SD instead of ratio to overall SD.
#'
#' @param x Matrix with 2 columns (two predictor variables).
#' @param y Response variable vector.
#' @param xout Logical; if TRUE, removes multivariate outliers from predictors (default: FALSE).
#' @param pts1 Points at which to evaluate predictor 1 (default: computed from data).
#' @param pts2 Points at which to evaluate predictor 2 (default: computed from data).
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param span Span parameter for LOESS (default: 2/3).
#' @param npts Number of evaluation points if pts1/pts2 not specified (default: 10).
#' @param tr Trimming proportion for Winsorized SD (default: 0.2).
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @return A list with components:
#'   \item{str1}{Strength of predictor 1 at each point in pts2 (conditional SD).}
#'   \item{str2}{Strength of predictor 2 at each point in pts1 (conditional SD).}
#'   \item{p}{Proportion of pairwise comparisons where str1 > str2.}
#'   \item{mean.str1}{Mean strength for predictor 1.}
#'   \item{mean.str2}{Mean strength for predictor 2.}
#'
#' @details
#' This function is similar to \code{\link{lplotcom2}} but uses **absolute conditional SD**
#' instead of the ratio to overall SD of y.
#'
#' For each predictor:
#' 1. Hold the other predictor fixed at several values
#' 2. Compute LOESS predictions varying only the focal predictor
#' 3. Measure the Winsorized SD of these predictions (not normalized)
#'
#' **Strength** = Winsorized SD of conditional predictions
#'
#' This version may be more interpretable when comparing the absolute amount of variation
#' explained by each predictor, rather than the proportion of total variation.
#'
#' Higher strength indicates the predictor produces more variability in predictions when
#' the other predictor is held fixed.
#'
#' @seealso \code{\link{lplotcom2}}, \code{\link{lplotcomBCI}}, \code{\link{lplot}}
#'
#' @export
#' @examples
#' # Compare importance using absolute SD
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- 2*x1 + 0.5*x2 + rnorm(100)
#' X <- cbind(x1, x2)
#' result <- lplotcom2v2(X, y, npts = 5)
#' result$mean.str1  # Absolute SD for x1
#' result$mean.str2  # Absolute SD for x2
lplotcom2v2<-function(x,y,xout=FALSE,pts1=NULL,pts2=NULL,outfun=outpro,span=2/3,npts=10,tr=.2,...){
#
# For two independent variables, estimate their relative importance when using LOESS
#
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n.orig=nrow(m)
n.keep=n.orig

d<-ncol(x)
if(d!=2)stop('Current version is for two independent variables only')
if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
} #$
n.keep=nrow(m)
M=apply(m,2,median)
SD=apply(m,2,mad)
low=M-1.5*SD
up=M+1.5*SD
if(is.null(pts1))pts1=seq(low[1],up[1],length.out=npts)
if(is.null(pts2))pts2=seq(low[2],up[2],length.out=npts)
e1=NA
e2=NA
for(j in 1:length(pts1)){     # Determine strength of x2 given a value stored in pts1.
v2=cbind(rep(pts1[j],n.keep),m[,2])
vals=lplot.pred(m[,1:2],m[,3],v2,span=span)$yhat
vals=elimna(vals)
nv=length(vals)
e2[j]=NA
if(nv>=10)e2[j]=winsd(vals,tr=tr,na.rm=TRUE)
}
for(j in 1:length(pts2)){     # Determine strength of x1 given a value stored in pts2.
v1=cbind(m[,1],rep(pts2[j],n.keep))
vals=lplot.pred(m[,1:2],m[,3],v1,span=span)$yhat
vals=elimna(vals)
nv=length(vals)
e1[j]=NA
if(nv>=10)e1[j]=winsd(vals,tr=tr,na.rm=TRUE)
 }
p=mean(outer(e1,e2,FUN='-')<0,na.rm=TRUE)
list(str1=e1,str2=e2,p=p,mean.str1=mean(e1),mean.str2=mean(e2))
}

#' Bootstrap Confidence Interval for Comparing Predictor Importance
#'
#' @description
#' For two predictors, compares their relative importance using LOESS regression
#' and provides bootstrap confidence intervals for the difference in importance.
#'
#' @param x Matrix with 2 columns (two predictor variables).
#' @param y Response variable vector.
#' @param xout Logical; if TRUE, removes multivariate outliers (default: FALSE).
#' @param pts1 Evaluation points for predictor 1 (default: quartiles 0.25, 0.5, 0.75).
#' @param pts2 Evaluation points for predictor 2 (default: quartiles 0.25, 0.5, 0.75).
#' @param p.crit Critical p-value for confidence intervals (default: computed from sample size).
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param span Span parameter for LOESS (default: 2/3).
#' @param npts Number of evaluation points if SEQ=TRUE (default: 10).
#' @param tr Trimming proportion for Winsorized SD (default: 0.2).
#' @param nboot Number of bootstrap samples (default: 500).
#' @param SEED Logical; if TRUE, sets random seed (default: TRUE).
#' @param SEQ Logical; if TRUE, uses equally spaced points instead of quartiles (default: FALSE).
#' @param MAD.OP Logical; if TRUE with SEQ=TRUE, uses MAD for range; otherwise uses quartiles (default: FALSE).
#' @param plotit Logical; if TRUE, creates LOESS plot (default: TRUE).
#' @param ticktype Type of tick marks for 3D plot (default: "simple").
#' @param xlab Label for X1 axis (default: "X1").
#' @param ylab Label for X2 axis (default: "X2").
#' @param zlab Label for Y axis (default: "Y").
#' @param reverse.x1 Logical; if TRUE, reverses order of pts1 (default: FALSE).
#' @param reverse.x2 Logical; if TRUE, reverses order of pts2 (default: FALSE).
#' @param pr Logical; if TRUE, prints warnings (default: FALSE).
#' @param MEDIAN Logical; if TRUE, evaluates only at median of predictors (default: FALSE).
#' @param Q1 Logical; if TRUE, evaluates only at lower quartile (default: FALSE).
#' @param Q2 Logical; if TRUE, evaluates only at upper quartile (default: FALSE).
#' @param alpha Significance level (default: 0.05).
#' @param MC Logical; if TRUE, uses parallel processing via mclapply (default: FALSE).
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @return A list with components:
#'   \item{p.crit}{Critical p-value used for confidence intervals.}
#'   \item{p.value}{P-value for testing whether predictor strengths differ.}
#'   \item{str.x1.given.x2}{Strength of x1 conditional on x2 at each evaluation point.}
#'   \item{str.x2.given.x1}{Strength of x2 conditional on x1 at each evaluation point.}
#'   \item{mean.str1}{Mean strength for predictor 1.}
#'   \item{mean.str2}{Mean strength for predictor 2.}
#'   \item{ci.low}{Lower confidence bound for difference (mean.str1 - mean.str2).}
#'   \item{ci.hi}{Upper confidence bound for difference.}
#'   \item{pts.x1}{Evaluation points used for predictor 1.}
#'   \item{pts.x2}{Evaluation points used for predictor 2.}
#'
#' @details
#' This function extends \code{\link{lplotcom2v2}} by adding bootstrap confidence intervals
#' for the difference in predictor importance.
#'
#' **Evaluation points**:
#' - **Default**: Uses quartiles (0.25, 0.5, 0.75) for both predictors
#' - **MEDIAN=TRUE**: Uses only the median
#' - **Q1=TRUE**: Uses only the lower quartile (0.25)
#' - **Q2=TRUE**: Uses only the upper quartile (0.75)
#' - **SEQ=TRUE**: Uses npts equally spaced points in the interquartile range (or MAD-based range if MAD.OP=TRUE)
#'
#' **Critical p-value** is automatically computed based on sample size if not specified.
#' This adjusts the confidence interval width to maintain proper coverage.
#'
#' The **p-value** tests H0: predictors have equal importance. Small values indicate
#' significant difference in importance.
#'
#' **Requires n >= 50** for reliable results.
#'
#' @seealso \code{\link{lplotcom2}}, \code{\link{lplotcom2v2}}, \code{\link{lplot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two predictors with bootstrap CI
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- 2*x1 + 0.5*x2 + rnorm(100)
#' X <- cbind(x1, x2)
#' result <- lplotcomBCI(X, y, nboot = 200)
#' result$p.value  # Test for difference
#' result$ci.low   # Lower bound for difference
#' result$ci.hi    # Upper bound
#' }
lplotcomBCI<-function(x,y,xout=FALSE,pts1=NULL,pts2=NULL,p.crit=NULL,
outfun=outpro,span=2/3,npts=10,tr=.2,nboot=500,
SEED=TRUE,SEQ=FALSE,MAD.OP=FALSE,plotit=TRUE,ticktype='simple',
xlab='X1',ylab='X2',zlab='Y',reverse.x1=FALSE,reverse.x2=FALSE,pr=FALSE,
MEDIAN=FALSE,Q1=FALSE,Q2=FALSE,alpha=.05,MC=FALSE,...){
#
# For two independent variables, estimate their relative importance when using LOESS
# p.crit is the critical p-value. If not specified, the function returns the approximate 0.05 critical p-value
#
# By default, use the average of the strength of the associations, so essentially a global test based on the quartiles
#  MEDIAN=TRUE, use the median of the independent variables only.
#  Q1=TRUE, use the lower quartile of the independent variables only.
#  Q2=TRUE, use the upper quartile of the independent variables only.
#
#  ADJ.CI=TRUE: Confidence intervals are based on the critical p-value
#  otherwise use alpha
#
if(SEED)set.seed(2)
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n=nrow(m)
x=m[,1:2]
y=m[,3]
if(xout){
flag=outfun(x,plotit=FALSE)$keep
x=x[flag,]
y=y[flag]
n=nrow(x)
}
if(n<50)stop('The sample size must be greater than or equal to 50')
if(MEDIAN){
pts1=median(x[,1])
pts2=median(x[,2])
if(is.null(p.crit)){
if(n<=200)p.crit=regYhat(c(1/50,1/100,1/200),c(.114,.080,.065),1/n)
if(n>200)p.crit=.062
}
}
if(Q1){
pts1=qest(x[,1],.25)
pts2=qest(x[,2],.25)
if(is.null(p.crit)){
if(n<=200)p.crit=regYhat(c(1/50,1/100,1/200),c(.142,.095,.082),1/n)
if(n>200)p.crit=.062
}
}
if(Q2){
pts1=qest(x[,1],.75)
pts2=qest(x[,2],.75)
if(is.null(p.crit)){
if(n<=200)p.crit=regYhat(c(1/50,1/100,1/200),c(.142,.095,.082),1/n)
if(n>200)p.crit=.062
}
}
if(is.null(pts1)){
pts1=qest(x[,1],c(.25,.5,.75))
if(is.null(p.crit)){
if(n<=200)p.crit=regYhat(c(1/50,1/100,1/200),c(.082,.076,.067),1/n)
if(n>200)p.crit=.06
}
if(reverse.x1)pts1=sort(pts1,TRUE)
if(is.null(pts2))pts2=qest(x[,2],c(.25,.5,.75))
if(reverse.x2)pts2=sort(pts2,TRUE)
}
if(SEQ){
if(MAD.OP){
M=apply(m,2,median)
SD=apply(m,2,mad)

low=M-1.5*SD
up=M+1.5*SD
}
else{
low=apply(m,2,qest,.25)
hi=apply(m,2,qest,.75)
}
pts1=seq(low[1],up[1],length.out=npts)
pts2=seq(low[2],up[2],length.out=npts)
}
v1=NA
v2=NA
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
if(!MC){
for(i in 1:nboot){
ib=data[i,]
temp=lplotcom2v2(x[ib,],y[ib],pts1=pts1,pts2=pts2,npts=npts,tr=tr,span=span)
v1[i]=temp$mean.str1
v2[i]=temp$mean.str2
}}
if(MC){
data=listm(t(data))
bvec<-mclapply(data,lplotCIMCv2,x,y,pts1=pts1,pts2=pts2,npts=npts,tr=tr,span=span)
bvec=matl(bvec)  # a 2-by-nboot matrix.
dif=sort(bvec[1,]-bvec[2,])
}
if(!MC)dif=sort(v1-v2)
nbl=length(dif)
#ilow<-round((alpha/2) * nbl)
ilow<-round((p.crit/2) * nbl)
ihi<-nbl - ilow
ilow<-ilow+1
ci.low=dif[ilow]
ci.hi=dif[ihi]
pv=mean(dif<0,na.rm=TRUE)
pv=2*min(pv,1-pv)
est=lplotcom2(x,y,xout=FALSE,pts1=pts1,pts2=pts2,outfun=outfun,span=span,
npts=npts,tr=tr)
if(plotit)lplot(x,y,ticktype=ticktype,xlab=xlab,ylab=ylab,zlab=zlab,pr=pr)
list(p.crit=p.crit,p.value=pv,str.x1.given.x2=est$str1,str.x2.given.x1=est$str2,mean.str1=est$mean.str1,
mean.str2=est$mean.str2,
ci.low=ci.low,ci.hi=ci.hi,pts.x1=pts1,pts.x2=pts2)
}

#' Compare Predictor Importance Using LOESS at Quartile Points (Bootstrap CI)
#'
#' @description
#' Compares the relative importance of two predictors using LOESS smoothing,
#' evaluated at quartile points (25th, 50th, 75th percentiles). Tests all 9
#' combinations of quartile pairs and provides bootstrap confidence intervals.
#'
#' @param x Matrix or data frame with 2 columns containing the predictor variables.
#' @param y Numeric vector of response values.
#' @inheritParams common-params
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param span Span parameter for LOESS (default: 2/3).
#' @param npts Number of points for measuring strength (default: 10).
#' @param nboot Number of bootstrap samples (default: 500).
#' @param SEED Logical; if `TRUE`, sets random seed for reproducibility.
#' @param plotit Logical; if `TRUE`, creates a LOESS plot.
#' @param ticktype Tick type for 3D plot (default: 'simple').
#' @param ADJ.CI Logical; if `TRUE`, adjusts CIs using critical p-value (default: TRUE).
#' @param xlab,ylab,zlab Axis labels for the plot.
#' @param alpha Significance level (default: 0.05).
#' @param MC Logical; if `TRUE`, uses parallel processing.
#'
#' @details
#' This function evaluates predictor importance at all 9 combinations of quartile
#' points (Q1, Q2, Q3) for both predictors. For each combination:
#' - Computes strength of association measures for each predictor
#' - Tests whether the two predictors differ in importance
#' - Provides bootstrap confidence intervals for the difference
#'
#' The critical p-value is adjusted based on sample size when `ADJ.CI = TRUE`.
#'
#' @return List with components:
#' \item{n}{Original sample size.}
#' \item{n.keep}{Sample size after outlier removal.}
#' \item{p.crit}{Critical p-value threshold.}
#' \item{output}{Matrix with 9 rows (one per quartile combination) and columns:
#'   \itemize{
#'     \item \code{pts1}: Evaluation point for predictor 1
#'     \item \code{pts2}: Evaluation point for predictor 2
#'     \item \code{p-value}: Two-sided p-value for difference test
#'     \item \code{str.x1.given.x2}: Strength of predictor 1 given predictor 2
#'     \item \code{str.x2.given.x1}: Strength of predictor 2 given predictor 1
#'     \item \code{ci.low}, \code{ci.hi}: Bootstrap CI for the difference
#'   }}
#'
#' @seealso \code{\link{lplotcom2}}, \code{\link{lplotcom2v2}}, \code{\link{lplotcomBCI}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two predictors at quartile points
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- x1 + 0.5*x2 + rnorm(100)
#' lplotcomBCI9(cbind(x1, x2), y)
#' }
lplotcomBCI9<-function(x,y,xout=FALSE,pr=TRUE,
outfun=outpro,span=2/3,npts=10,tr=.2,nboot=500,
SEED=TRUE,plotit=TRUE,ticktype='simple',ADJ.CI=TRUE,
xlab='X1',ylab='X2',zlab='Y',alpha=.05,MC=FALSE,...){
#
# For two independent variables, estimate their relative importance when using LOESS
# Focus  on the quartiles: none tests based on all possible combinations.
#
p.crit=NA
if(pr){
if(alpha!=.05){
if(pr)print('Critical p-value is taken to be the value of alpha. Unknown how to adjust  when alpha is not .05')
p.crit=alpha
}
if(ADJ.CI)print('Confidence intervals are based on the critical p-value')
}
if(SEED)set.seed(2)
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n=nrow(m)
n.orig=n
x=m[,1:2]
y=m[,3]
if(xout){
flag=outfun(x,plotit=FALSE)$keep
x=x[flag,]
y=y[flag]
n=nrow(x)
}
if(is.na(p.crit)){
if(n<=100)p.crit=regYhat(c(1/50,1/100),c(.042,.025),1/n)
else p.crit=.025
}
output<-matrix(NA,nrow=9,ncol=7)
dimnames(output)=list(NULL,c('pts1','pts2','p-value','str.x1.given.x2','str.x2.given.x1','ci.low','ci.hi'))
pts1=qest(x[,1],c(.25,.5,.75))
pts2=qest(x[,2],c(.25,.5,.75))

v1=matrix(NA,nrow=nboot,ncol=3)
v2=matrix(NA,nrow=nboot,ncol=3)

data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
if(!MC){
for(i in 1:nboot){
ib=data[i,]
for(j in 1:3){
temp=lplotcom2v2(x[ib,],y[ib],pts1=pts1[j],pts2=pts2[j],npts=npts,tr=tr,span=span)
v1[i,j]=temp$mean.str1
v2[i,j]=temp$mean.str2
}}}
if(MC){
data=listm(t(data))
for(j in 1:3){
bvec<-mclapply(data,lplotCIMCv2,x,y,pts1=pts1[j],pts2=pts2[j],npts=npts,tr=tr,span=span)
bvec=matl(bvec)  # a 2-by-nboot matrix.
v1[,j]=bvec[1,]
v2[,j]=bvec[2,]
}
}
pc=matrix(NA,3,3) #rows are for pts1, columns for pts2
ic=0
for(j in 1:3){
for(k in 1:3){
est=lplotcom2(x,y,xout=FALSE,pts1=pts1[j],pts2=pts2[k],outfun=outfun,span=span,
npts=npts,tr=tr)
ic=ic+1
output[ic,1]=pts1[j]
output[ic,2]=pts2[k]
dif=sort(v1[,j]-v2[,k])
nbl=length(dif)
if(ADJ.CI)ilow<-round((p.crit/2) * nbl)
else ilow<-round((alpha/2) * nbl)
ihi<-nbl - ilow
ilow<-ilow+1
ci.low=dif[ilow]
ci.hi=dif[ihi]
pv=mean(dif<0,na.rm=TRUE)
pc[j,k]=2*min(pv,1-pv)
output[ic,3]=pc[j,k]
output[ic,6]=ci.low
output[ic,7]=ci.hi
output[ic,4]=est$mean.str1
output[ic,5]=est$mean.str2
}}
if(plotit)lplot(x,y,ticktype=ticktype,xlab=xlab,ylab=ylab,zlab=zlab)
list(n=n.orig,n.keep=n,p.crit=p.crit,output=output)
}

#' Bootstrap Helper for lplotcom2 (Parallel Processing)
#'
#' @description
#' Internal helper function for parallel bootstrap computation in \code{\link{lplotcomBCI}}.
#' Computes mean strength measures using \code{\link{lplotcom2}}.
#'
#' @param data Bootstrap sample indices.
#' @param x Predictor matrix.
#' @param y Response vector.
#' @param pts1 Evaluation points for predictor 1.
#' @param pts2 Evaluation points for predictor 2.
#' @param npts Number of evaluation points.
#' @param tr Trimming proportion for Winsorized SD.
#' @param span Span parameter for LOESS.
#'
#' @return Vector of length 2 containing mean.str1 and mean.str2.
#'
#' @keywords internal
lplotCIMC<-function(data,x,y,pts1,pts2,npts,tr,span){
temp=lplotcom2(x[data,],y[data],pts1=pts1,pts2=pts2,npts=npts,tr=tr,span=span)
v=c(temp$mean.str1,temp$mean.str2)
}

#' Bootstrap Helper for lplotcom2v2 (Parallel Processing)
#'
#' @description
#' Internal helper function for parallel bootstrap computation in \code{\link{lplotcomBCI}}.
#' Computes mean strength measures using \code{\link{lplotcom2v2}}.
#'
#' @param data Bootstrap sample indices.
#' @param x Predictor matrix.
#' @param y Response vector.
#' @param pts1 Evaluation points for predictor 1.
#' @param pts2 Evaluation points for predictor 2.
#' @param npts Number of evaluation points.
#' @param tr Trimming proportion for Winsorized SD.
#' @param span Span parameter for LOESS.
#'
#' @return Vector of length 2 containing mean.str1 and mean.str2.
#'
#' @keywords internal
lplotCIMCv2<-function(data,x,y,pts1,pts2,npts,tr,span){
temp=lplotcom2v2(x[data,],y[data],pts1=pts1,pts2=pts2,npts=npts,tr=tr,span=span)
v=c(temp$mean.str1,temp$mean.str2)
}

#' Critical P-Value for LOESS Confidence Intervals
#'
#' @description
#' Determines the critical p-value for confidence intervals in \code{\link{lplotCI}}
#' via simulation under the null hypothesis of independence.
#'
#' @param n Sample size.
#' @param nreps Number of simulation replications (default: 2000).
#' @param alpha Significance level (default: 0.05).
#' @param pts Evaluation points (default: computed from simulated data).
#' @param npts Number of evaluation points if pts is NULL (default: 25).
#'
#' @return Critical p-value at the specified alpha level.
#'
#' @details
#' This function simulates the null distribution of minimum p-values across
#' evaluation points when x and y are independent standard normal variables.
#' The critical p-value is the alpha-quantile of this null distribution.
#'
#' @keywords internal
lplotbsepvv3<-function(n,nreps=2000,alpha=0.05,pts=NULL,npts=25){
#
# Determine critical p-value for lplotCI.
#
set.seed(2)
pv=NA
for(i in 1:nreps){
x=rnorm(n)
y=rnorm(n)
xord=order(x)
y=y[xord]
x=x[xord]
M=median(x)
low=M-1.5*mad(x)
up=M+1.5*mad(x)
if(is.null(pts))pts=seq(low,up,length.out=npts)
plx<-predict(loess(y ~ x), se=TRUE)
est=lplot.pred(x,y,pts=pts)$yhat
se=lplotse(x,y,SEED=FALSE,pts=pts)
test=abs(est/se)
pall=2*(1-pt(abs(test),plx$df))
pv[i]=min(elimna(pall))
}
hd(pv,alpha)
}

#' LOESS Smoother Plot for Large Datasets
#'
#' @description
#' Creates a 3D plot of a LOESS regression surface for large datasets by subsampling
#' points for visualization. Good for large sample sizes or plotting regression
#' surfaces without scatterplots.
#'
#' @param x Predictor variable(s). Can be a vector (univariate) or matrix (multivariate).
#' @param y Response variable (vector).
#' @param nsub Size of random subsample used for plotting (default: 5000). The function
#'   uses all data for prediction but plots only a subsample.
#' @param est Measure of location for running interval smoother (default: tmean, 20% trimmed mean).
#' @param fr Span for LOESS smoother (default: 1).
#' @param xout Logical. Remove outliers in predictors using projection method if TRUE (default: FALSE).
#' @param xlab Label for X-axis (default: "X").
#' @param ylab Label for Y-axis (default: "Y").
#' @param zlab Label for Z-axis (default: "").
#' @param ticktype Type of tick marks (default: "simple").
#' @param theta Azimuthal angle for 3D plot viewing (default: 50).
#' @param phi Colatitude angle for 3D plot viewing (default: 25).
#' @param scale Logical. Scale axes to same range if TRUE (default: TRUE).
#' @param pc Point character for plotting (default: " ").
#' @param expand Expansion factor for bounding box (default: 0.5).
#' @param SEED Logical. Set random seed for reproducibility if TRUE (default: TRUE).
#' @param frame Logical. Draw frame around plot if TRUE (default: TRUE).
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function is designed for large datasets where plotting all points would be
#' computationally expensive or visually cluttered. It works by:
#' 1. Randomly subsampling nsub points from the data for visualization
#' 2. Computing LOESS predictions using all data
#' 3. Plotting the subsampled points with predicted values
#'
#' The function uses \code{lplot.pred} to compute predictions and \code{lplot} for
#' 3D visualization. Missing values are automatically removed.
#'
#' **Outlier removal** (xout = TRUE): Uses \code{outpro.depth} projection method
#' to remove outliers before plotting.
#'
#' @seealso \code{\link{lplot}}, \code{\link{lplot.pred}}, \code{\link{rplotN}}
#' @export
#' @examples
#' # Large dataset example
#' set.seed(123)
#' n <- 10000
#' x <- rnorm(n)
#' y <- 2 * x + rnorm(n)
#' lplotN(x, y, nsub = 1000)
lplotN<-function(x,y,nsub=5000,est=tmean,fr=1,xout=FALSE,xlab='X',ylab='Y',zlab='',ticktype = 'simple',theta = 50, phi = 25, scale = TRUE,pc=' ',
    expand = 0.5, SEED = TRUE,frame=TRUE){
#
 # Running interval smoother, good for large sample sizes or plots of the
 # regression surface without a scatter plot.
 #
 # nsub is size of the random sample of the data used to predict outcome using all of the data
 #  
 if(SEED)set.seed(2)
 x=as.matrix(x)
p=ncol(x)
p1=p+1
 xy=cbind(x,y)
 xy=elimna(xy)
  if(xout){
 flag<-outpro.depth(x,plotit=FALSE)$keep                                
x<-x[flag,]                                                                        
x<-as.matrix(x)                                                                    
y<-y[flag]    
xy=cbind(x,y)                                                                     
} 
 n=nrow(xy)
 nsub=min(n,nsub)
 id=sample(n,nsub) 
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
w=lplot.pred(x,y,pts=x[id,],fr=fr)$yhat 
a=lplot(x[id,],w,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,pc=pc,
frame=frame,phi=phi,theta=theta,scale=scale,pr=FALSE)   
}


# ============================================================================
# GAM PLOTS
# ============================================================================

#' Plot Distribution of Linear Contrast
#'
#' @description
#' Plots the bootstrap distribution of a linear contrast c1*X1 + c2*X2 + ... + cJ*XJ
#' from multiple groups, useful for visualizing the distribution of linear combinations
#' of group estimates.
#'
#' @param x Data in list mode or matrix mode. Each element/column represents a group.
#' @param con Contrast coefficients (vector of length J). If not specified or all zeros,
#'   defaults to c(1,1,...,1), plotting the distribution of the sum.
#' @param plotfun Density estimation function for plotting (default: akerd, adaptive kernel estimator).
#' @param nboot Number of bootstrap samples (default: 800).
#' @param plotit Logical. Create plot if TRUE (default: TRUE).
#' @param pyhat Logical. Passed to plotfun for additional output (default: FALSE).
#' @param ... Additional arguments passed to plotfun.
#'
#' @return Output from the density estimation function specified by plotfun.
#'
#' @details
#' This function visualizes the distribution of a linear contrast across multiple groups
#' by bootstrap resampling. It's useful for examining the sampling distribution of
#' linear combinations such as:
#' - Sum of group means: c = (1, 1, ..., 1)
#' - Difference: c = (1, -1)
#' - Custom contrasts: any vector of coefficients
#'
#' **Procedure**:
#' 1. For each of nboot iterations, sample with replacement from each group
#' 2. Compute the linear contrast for each bootstrap sample
#' 3. Plot the resulting distribution using the specified plotfun
#'
#' Missing values are automatically removed from each group.
#'
#' @seealso \code{\link{lin2plot}}, \code{\link{akerd}}, \code{\link{lincon}}
#' @export
#' @examples
#' # Plot distribution of sum of three groups
#' set.seed(123)
#' x1 <- rnorm(30, mean = 0)
#' x2 <- rnorm(30, mean = 1)
#' x3 <- rnorm(30, mean = 2)
#' x <- list(x1, x2, x3)
#' linplot(x)  # Default: sum of all groups
#'
#' # Plot distribution of difference (Group 1 - Group 2)
#' linplot(list(x1, x2), con = c(1, -1))
linplot<-function(x,con=0,plotfun=akerd,nboot=800,plotit=TRUE,pyhat=FALSE,...){
#
#  plot distribtion of the linear contrast
#  c_1X_2+c_2X_2+...
#
#  con contains contrast coefficients. If not specified,
#  con<-c(1,1,...,1)
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
J<-length(x)
tempn<-0
mvec<-NA
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
}
Jm<-J-1
#
# Determine contrast matrix
# If not specified, assume distribution of the sum is to be plotted
#
if(sum(con^2)==0)con<-matrix(1,J,1)
bvec<-matrix(NA,nrow=J,ncol=nboot)
for(j in 1:J){
data<-matrix(sample(x[[j]],size=nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-data
}
bcon<-t(con)%*%bvec #ncon by nboot matrix
bcon<-as.vector(bcon)
dval<-plotfun(bcon,pyhat=pyhat,...)
dval
}

#' Plot Two Distributions Based on Positive and Negative Contrast Coefficients
#'
#' @description
#' Plots two distributions: one for the linear contrast of groups with positive
#' coefficients and one for groups with negative coefficients. Useful for visualizing
#' the two sides of a contrast separately.
#'
#' @param x Data in list mode or matrix mode. Each element/column represents a group.
#' @param con Contrast coefficients (vector of length J, required). Positive coefficients
#'   define one distribution, negative coefficients define the other.
#' @param op Density estimation method for plotting (default: 4, adaptive kernel estimator).
#'   See \code{\link{g2plot}} for options.
#' @param nboot Number of bootstrap samples (default: 800).
#' @param plotit Logical. Create plot if TRUE (default: TRUE).
#' @param pyhat Logical. Passed to plotting function for additional output (default: FALSE).
#'
#' @return Output from \code{g2plot} containing density estimates for both distributions.
#'
#' @details
#' This function separates a contrast into two parts based on the sign of coefficients:
#' - **Distribution 1**: Sum of groups with positive coefficients (ci > 0)
#' - **Distribution 2**: Sum of groups with negative coefficients (ci < 0)
#'
#' This is particularly useful for understanding contrasts where you're comparing
#' one set of groups against another (e.g., treatment vs. control groups).
#'
#' **Example contrast**: c = (1, 1, -1, -1) compares the sum of groups 1 & 2 against
#' the sum of groups 3 & 4. This function plots:
#' - Distribution of X1 + X2
#' - Distribution of X3 + X4
#'
#' The two distributions are plotted using \code{g2plot} with overlaid density estimates.
#' Missing values are automatically removed.
#'
#' **Note**: The con argument is required and must have both positive and negative values.
#'
#' @seealso \code{\link{linplot}}, \code{\link{g2plot}}, \code{\link{lincon}}
#' @export
#' @examples
#' # Compare two sets of groups
#' set.seed(123)
#' x1 <- rnorm(30, mean = 0)
#' x2 <- rnorm(30, mean = 0.5)
#' x3 <- rnorm(30, mean = 1)
#' x4 <- rnorm(30, mean = 1.5)
#' x <- list(x1, x2, x3, x4)
#'
#' # Compare (X1 + X2) vs (X3 + X4)
#' lin2plot(x, con = c(1, 1, -1, -1))
lin2plot<-function(x,con,op=4,nboot=800,plotit=TRUE,pyhat=FALSE){
#
#  plot two distribtions.
#   The first is the distribtion  of the linear contrast
#  c_1X_2+c_2X_2+... c_i>0
#  and the second is the distribution of c_1X_2+c_2X_2+... c_i<0
#
#  con contains contrast coefficients. If not specified,
#  function terminates.
#
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
J<-length(x)
if(J != length(con)){
stop("Number of contrast coefficients must equal the number of groups")
}
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
#
# Determine contrast matrix for positive contrast coefficients
#
flag<-(con<0)
con1<-con
con1[flag]<-0
# Determine contrast matrix for negative contrast coefficients
flag<-(con>0)
con2<-con
con2[flag]<-0
bvec<-matrix(NA,nrow=J,ncol=nboot)
for(j in 1:J){
data<-matrix(sample(x[[j]],size=nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-data
}
bcon1<-t(con1)%*%bvec
bcon2<-t(con2)%*%bvec
bcon1<-as.vector(bcon1)
bcon2<-as.vector(bcon2)
fval<-g2plot(bcon1,bcon2,op=op,rval=15,fr=0.8,aval=0.5,xlab="X",ylab="")
fval
}


# ============================================================================
# GROUP COMPARISON PLOTS
# ============================================================================

#' Plot Density Functions for Two Groups
#'
#' @description
#' Creates overlaid density plots for two independent groups, allowing visual
#' comparison of distributions. Multiple estimation methods are available.
#'
#' @param x1 Data for group 1 (vector).
#' @param x2 Data for group 2 (vector).
#' @param op Option for density estimation method (default: 4):
#'   \itemize{
#'     \item 1: Rosenblatt shifted histogram
#'     \item 2: Kernel density estimate (Epanechnikov kernel)
#'     \item 3: Expected frequency curve
#'     \item 4: Adaptive kernel estimator (recommended)
#'   }
#' @param rval Number of points for density estimation when op=1 (default: 15).
#' @param fr Span parameter for expected frequency curve when op=3 (default: 0.8).
#' @param aval Scaling parameter for adaptive kernel estimator when op=4 (default: 0.5).
#' @param xlab Label for X-axis (default: "X").
#' @param ylab Label for Y-axis (default: "").
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function compares the distributions of two groups by plotting density estimates.
#'
#' **Density estimation methods**:
#' - **op=1**: Rosenblatt shifted histogram using kernel density estimation.
#' - **op=2**: Standard kernel density estimate with Epanechnikov kernel.
#' - **op=3**: Expected frequency curve using running interval smoother (calls \code{rd2plot}).
#'   If graph is ragged, suggests using op=4.
#' - **op=4**: Adaptive kernel estimator (recommended). Uses \code{akerd} for smoother
#'   density estimates that adapt to local data density.
#'
#' **Plot appearance**:
#' - Group 1: solid line (lty=1)
#' - Group 2: dashed line (lty=2)
#'
#' Missing values are automatically removed.
#'
#' @seealso \code{\link{akerd}}, \code{\link{rd2plot}}, \code{\link{g5plot}}
#' @export
#' @examples
#' # Compare two normal distributions
#' set.seed(123)
#' x1 <- rnorm(100, mean = 0, sd = 1)
#' x2 <- rnorm(100, mean = 1, sd = 1.5)
#' g2plot(x1, x2)
#'
#' # Try different estimation methods
#' g2plot(x1, x2, op = 2)  # Standard kernel density
#' g2plot(x1, x2, op = 3)  # Expected frequency curve
g2plot<-function(x1,x2,op=4,rval=15,fr=.8,aval=.5,xlab="X",ylab=""){
#
# plot estimates of the density functions for two groups.
#
# op=1: Use Rosenblatt shifted histogram
#
# op=2:
# Use kernel density estimate
# Using the built-in S+ function density,
#
# op=3: Use expected frequency curve.
#
# op=4: Use adaptive kernel estimator
#
x1<-elimna(x1)
x2<-elimna(x2)
if(op==3){
rd2plot(x1,x2,fr=fr,xlab=xlab,ylab=ylab)
print("Might consider using op=4 if graph is ragged")
}
if(op==2){
tempx<-density(x1,na.rm=TRUE,kernel="epanechnikov")
tempy<-density(x2,na.rm=TRUE,kernel="epanechnikov")
plot(c(tempx$x,tempy$x),c(tempx$y,tempy$y),type="n",xlab=xlab,ylab=ylab)
lines(tempx$x,tempx$y)
lines(tempy$x,tempy$y,lty=2)
}
if(op==1){
        y1 <- sort(x1)
        z1 <- 1
        z2 <- 1
        par(yaxt = "n")
        temp <- floor(0.01 * length(x1))
        if(temp == 0)
                temp <- 5
        ibot <- y1[temp]
        itop <- y1[floor(0.99 * length(x1))]
        xaxis1 <- seq(ibot, itop, length = rval)
        for(i in 1:rval)
                z1[i] <- kerden(x1, 0, xaxis1[i])
        y2 <- sort(x2)
         temp <- floor(0.01 * length(x2))
        if(temp == 0)
                temp <- 5
        ibot <- y2[temp]
        itop <- y2[floor(0.99 * length(x2))]
        xaxis2 <- seq(ibot, itop, length = rval)
        for(i in 1:rval)
                z2[i] <- kerden(x2, 0, xaxis2[i])
plot(c(xaxis1,xaxis2),c(z1,z2), xlab =xlab, ylab =ylab, type = "n")
lines(xaxis1,z1)
lines(xaxis2,z2,lty=2)
}
if(op==4){
x1<-sort(x1)
x2<-sort(x2)
z1<-akerd(x1,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
z2<-akerd(x2,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
plot(c(x1,x2),c(z1,z2), xlab =xlab, ylab =ylab, type = "n")
lines(x1,z1)
lines(x2,z2,lty=2)
}
}

#' Plot Distribution of Pairwise Differences Between Two Groups
#'
#' @description
#' Plots the estimated distribution of all pairwise differences (X - Y) between
#' two independent samples using an adaptive kernel density estimator.
#'
#' @param x Data for group 1 (vector).
#' @param y Data for group 2 (vector).
#' @param xlab Label for X-axis (default: "Difference").
#' @param ylab Label for Y-axis (default: "").
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function computes all n1 × n2 pairwise differences between observations
#' in group 1 (x) and group 2 (y), then plots the distribution using an adaptive
#' kernel density estimator (\code{akerd}).
#'
#' Missing values are automatically removed before computing differences.
#'
#' This plot is useful for visualizing the shift between two distributions and
#' assessing the entire distribution of effect sizes, not just central tendency.
#'
#' @seealso \code{\link{akerd}}, \code{\link{g2plot}}, \code{\link{difQplot}}
#' @export
#' @examples
#' # Plot distribution of differences
#' set.seed(123)
#' x <- rnorm(50, mean = 0, sd = 1)
#' y <- rnorm(50, mean = 1, sd = 1)
#' g2plotdifxy(x, y)
g2plotdifxy<-function(x,y,xlab="Difference",ylab=""){
#
# Plot an estimate of the distribution of X-Y
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
m<-as.vector(outer(x,y,FUN="-"))
akerd(m,xlab=xlab,ylab=ylab)
}

#' Plot Density Functions for Up to Five Groups
#'
#' @description
#' Creates overlaid density plots for up to five independent groups using an
#' adaptive kernel density estimator, allowing visual comparison of distributions.
#'
#' @param x1 Data for group 1 (vector), or a list/matrix/data.frame containing all groups.
#' @param x2 Data for group 2 (vector).
#' @param x3 Data for group 3 (vector, optional, default: NULL).
#' @param x4 Data for group 4 (vector, optional, default: NULL).
#' @param x5 Data for group 5 (vector, optional, default: NULL).
#' @param fr Span parameter for adaptive kernel estimator (default: 0.8).
#' @param aval Scaling parameter for adaptive kernel estimator (default: 0.5).
#' @param xlab Label for X-axis (default: 'X').
#' @param ylab Label for Y-axis (default: '').
#' @param color Vector of colors for the 5 groups (default: rep('black', 5)).
#' @param main Main title for plot (default: NULL).
#' @param sub Subtitle for plot (default: NULL).
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function plots density estimates for up to 5 groups using an adaptive
#' kernel density estimator (\code{akerd}) which adapts to local data density
#' for smoother estimates.
#'
#' **Input formats**:
#' - **Separate vectors**: Provide x1, x2, and optionally x3, x4, x5.
#' - **List/matrix/data.frame**: Provide x1 as a list, matrix, or data.frame
#'   containing all groups (up to 5).
#'
#' **Plot appearance**:
#' - Group 1: solid line (lty=1)
#' - Group 2: dashed line (lty=2)
#' - Group 3: dotted line (lty=3)
#' - Group 4: dotdash line (lty=4)
#' - Group 5: longdash line (lty=5)
#'
#' Each group can have a different color specified via the color parameter.
#'
#' Missing values are automatically removed.
#'
#' @seealso \code{\link{akerd}}, \code{\link{g2plot}}, \code{\link{g5.cen.plot}}
#' @export
#' @examples
#' # Compare three groups
#' set.seed(123)
#' x1 <- rnorm(100, mean = 0, sd = 1)
#' x2 <- rnorm(100, mean = 1, sd = 1.2)
#' x3 <- rnorm(100, mean = -0.5, sd = 0.8)
#' g5plot(x1, x2, x3)
#'
#' # Using a list
#' data_list <- list(x1, x2, x3)
#' g5plot(data_list)
g5plot<-function(x1,x2,x3=NULL,x4=NULL,x5=NULL,fr=.8,aval=.5,xlab='X',ylab='',color=rep('black',5),main=NULL,sub=NULL){
#
# plot estimates of the density functions for up to 5 groups.
# using an adaptive kernel density estimator
#
if(is.matrix(x1)||is.data.frame(x1))x1=listm(x1)
if(is.list(x1)){
x=x1
J=length(x)
ic=0
for(j in 1:J){
ic=ic+1
if(ic==1)x1=x[[1]]
if(ic==2)x2=x[[2]]
if(ic==3)x3=x[[3]]
if(ic==4)x4=x[[4]]
if(ic==5)x5=x[[5]]
}
}
x1<-elimna(x1)
x2<-elimna(x2)
x1<-sort(x1)
x2<-sort(x2)
if(!is.null(x3))x3<-sort(x3)
if(!is.null(x4))x4<-sort(x4)
if(!is.null(x5))x5<-sort(x5)
z3=NULL
z4=NULL
z5=NULL
z1<-akerd(x1,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
z2<-akerd(x2,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
if(!is.null(x3))z3=akerd(x3,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
if(!is.null(x4))z4=akerd(x4,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
if(!is.null(x5))z5=akerd(x5,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
plot(c(x1,x2,x3,x4,x5),c(z1,z2,z3,z4,z5), xlab =xlab, ylab =ylab, type = 'n',main=main,sub=sub)
lines(x1,z1,col=color[1])
lines(x2,z2,lty=2,col=color[2])
if(!is.null(x3))lines(x3,z3,lty=3,col=color[3])
if(!is.null(x4))lines(x4,z4,lty=4,col=color[4])
if(!is.null(x5))lines(x5,z5,lty=5,col=color[5])
}

#' Plot Centered Density Functions for Up to Five Groups
#'
#' @description
#' Creates overlaid density plots for up to five independent groups after centering
#' each group by a specified measure of location (e.g., median or mean). This allows
#' comparison of distributional shapes independent of location shifts.
#'
#' @param x1 Data for group 1 (vector), or a list/matrix/data.frame containing all groups.
#' @param x2 Data for group 2 (vector).
#' @param x3 Data for group 3 (vector, optional, default: NULL).
#' @param x4 Data for group 4 (vector, optional, default: NULL).
#' @param x5 Data for group 5 (vector, optional, default: NULL).
#' @param fr Span parameter for adaptive kernel estimator (default: 0.8).
#' @param aval Scaling parameter for adaptive kernel estimator (default: 0.5).
#' @param xlab Label for X-axis (default: 'X').
#' @param ylab Label for Y-axis (default: '').
#' @param color Vector of colors for the 5 groups (default: rep('black', 5)).
#' @param main Main title for plot (default: NULL).
#' @param sub Subtitle for plot (default: NULL).
#' @param loc.fun Function for computing center (default: \code{median}).
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function is identical to \code{\link{g5plot}} except that each group is
#' centered by subtracting the value computed by loc.fun before plotting.
#'
#' **Centering**: Each group x_i is transformed to x_i - loc.fun(x_i) before
#' density estimation. This removes location differences and allows comparison
#' of distributional shapes (spread, skewness, etc.).
#'
#' Common choices for loc.fun include:
#' - \code{median} (default, robust)
#' - \code{mean}
#' - \code{tmean} (trimmed mean)
#'
#' The plot uses the same line styles and colors as \code{g5plot}.
#'
#' @seealso \code{\link{g5plot}}, \code{\link{g2plot}}, \code{\link{akerd}}
#' @export
#' @examples
#' # Compare shapes after centering by median
#' set.seed(123)
#' x1 <- rnorm(100, mean = 0, sd = 1)
#' x2 <- rnorm(100, mean = 5, sd = 1.5)  # Different location and spread
#' x3 <- rnorm(100, mean = -3, sd = 0.8)
#' g5.cen.plot(x1, x2, x3)
#'
#' # Center by mean instead of median
#' g5.cen.plot(x1, x2, x3, loc.fun = mean)
g5.cen.plot<-function(x1, x2, x3 = NULL, x4 = NULL, x5 = NULL, fr = 0.8,
    aval = 0.5, xlab = 'X', ylab ='', color = rep('black', 5),
    main = NULL, sub = NULL,loc.fun=median){
#
#  Same a g5plot, only center the data based on the
#  measure of location indicated by the argument
#  loc.fun
#
x1=elimna(x1)
x2=elimna(x2)
x1=x1-loc.fun(x1)
x2=x2-loc.fun(x2)
if(!is.null(x3))x3=x3-loc.fun(x3)
if(!is.null(x4))x4=x4-loc.fun(x4)
if(!is.null(x5))x5=x5-loc.fun(x5)
g5plot(x1=x1, x2=x2, x3=x3, x4 = x4, x5 = x5, fr = fr,
    aval = aval, xlab = xlab, ylab =ylab, color = color,
    main = main, sub = sub)
}

#' Group Scatterplot with Jittered Points
#'
#' @description
#' Creates a scatterplot displaying observations from multiple groups with group
#' membership on the X-axis and values on the Y-axis. Useful for visualizing
#' group differences and within-group variation.
#'
#' @param x Data in list or matrix format. If matrix, each column represents a group.
#' @param xlab Label for X-axis (default: "Group").
#' @param ylab Label for Y-axis (default: "").
#' @param xnum Logical; if TRUE, shows numeric group labels on X-axis; if FALSE,
#'   suppresses X-axis labels (default: FALSE).
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function creates a simple scatterplot where the X-axis represents group
#' membership (1, 2, 3, ...) and the Y-axis shows the observed values.
#'
#' If x is a matrix, it is converted to list format using \code{listm}.
#'
#' **Plot appearance**:
#' - Each point represents one observation
#' - X-coordinate indicates group membership
#' - Y-coordinate shows the observed value
#' - By default, X-axis tick labels are suppressed (set xnum=TRUE to show group numbers)
#'
#' This plot is useful for:
#' - Visualizing raw data alongside boxplots or density plots
#' - Identifying outliers within groups
#' - Assessing within-group variability
#' - Comparing group locations and spreads
#'
#' @seealso \code{\link{g2plot}}, \code{\link{g5plot}}, \code{\link{ebarplot}}
#' @export
#' @examples
#' # Create group scatterplot
#' set.seed(123)
#' x1 <- rnorm(30, mean = 0, sd = 1)
#' x2 <- rnorm(30, mean = 1, sd = 1.5)
#' x3 <- rnorm(30, mean = -0.5, sd = 0.8)
#' data_list <- list(x1, x2, x3)
#' gplot(data_list)
#'
#' # Show group numbers on X-axis
#' gplot(data_list, xnum = TRUE)
gplot<-function(x,xlab="Group",ylab="",xnum=FALSE){
if(is.matrix(x))x<-listm(x)
if(!xnum)par(xaxt="n")
mval<-NA
vals<-x[[1]]
gval<-rep(1,length(x[[1]]))
for(j in 2:length(x)){
vals<-c(vals,x[[j]])
gval<-c(gval,rep(j,length(x[[j]])))
}
plot(gval,vals,xlab=xlab,ylab=ylab)
}

#' Plot LOESS Smoothers for Two Groups
#'
#' @description
#' Creates a scatterplot with LOESS (locally weighted scatterplot smoothing) curves
#' for two independent groups, allowing visual comparison of regression trends.
#'
#' @param x1 Predictor variable for group 1 (vector).
#' @param y1 Response variable for group 1 (vector).
#' @param x2 Predictor variable for group 2 (vector).
#' @param y2 Response variable for group 2 (vector).
#' @param f Span parameter for LOESS smoother (default: 2/3). Controls the amount
#'   of smoothing. Larger values create smoother curves.
#' @param SCAT Logical. Include scatterplot if TRUE (default: TRUE). If FALSE,
#'   only the regression lines are plotted.
#' @param xlab Label for X-axis (default: "x").
#' @param ylab Label for Y-axis (default: "y").
#' @param pch Point character for scatterplot (default: '*').
#' @param eout Logical. Currently not used (default: FALSE).
#' @param xout Logical. Currently not used (default: FALSE).
#' @param ... Additional graphical parameters passed to plot.
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function compares regression trends between two groups by overlaying
#' their LOESS smoothers on a single plot. LOESS is a nonparametric method that
#' fits smooth curves through data without assuming a parametric form.
#'
#' **Procedure**:
#' 1. Combine data from both groups into a single scatterplot
#' 2. Fit separate LOESS smoothers to each group using R's \code{lowess} function
#' 3. Overlay both smoothed curves on the plot
#'
#' The span parameter \code{f} controls smoothing:
#' - **Small f** (e.g., 0.3): Less smoothing, follows data more closely
#' - **Large f** (e.g., 0.9): More smoothing, smoother curve
#'
#' Missing values are automatically removed from each group before plotting.
#'
#' **Note**: Uses \code{lowess} (not \code{loess}) for the smoothing algorithm.
#'
#' @seealso \code{\link{lplot2g}}, \code{\link{reg2plot}}, \code{\link{lowess}}
#' @export
#' @examples
#' # Compare regression trends for two groups
#' set.seed(123)
#' n <- 50
#' x1 <- runif(n, 0, 10)
#' y1 <- 2 * x1 + rnorm(n, sd = 2)
#' x2 <- runif(n, 0, 10)
#' y2 <- 1.5 * x2 + 3 + rnorm(n, sd = 2)
#' l2plot(x1, y1, x2, y2)
#'
#' # Just the smoothed lines without scatterplot
#' l2plot(x1, y1, x2, y2, SCAT = FALSE)
l2plot<-function(x1,y1,x2,y2,f=2/3,SCAT=TRUE,xlab="x",ylab="y",pch='*',
eout=FALSE,xout=FALSE,...){
#
# Plot LOESS smoother for two groups
#
# f is the span used by loess
# SCAT=F, scatterplot not created, just the regression lines
# Missing values are automatically removed.
#
m<-elimna(cbind(x1,y1))
x1<-m[,1]
y1<-m[,2]
m<-elimna(cbind(x2,y2))
x2<-m[,1]
y2<-m[,2]
plot(c(x1,x2),c(y1,y2),xlab=xlab,ylab=ylab,pch=pch)
lines(lowess(x1,y1,f=f))
lines(lowess(x2,y2,f=f))
}

#' Summary Plot Panel for Two-Group Comparison
#'
#' @description
#' Creates a 2×2 panel of four diagnostic plots for comparing two independent groups:
#' error bar plot, boxplots, kernel density estimates, and shift function.
#'
#' @param x Data for group 1 (vector), or a matrix/list containing both groups.
#' @param y Data for group 2 (vector, optional if x is a matrix/list, default: NULL).
#' @param xlab Label for X-axis in density plot (default: "X").
#' @param ylab Label for Y-axis in density plot (default: "").
#' @param eblabx Label for X-axis in error bar plot (default: "Groups").
#' @param eblaby Label for Y-axis in error bar plot (default: "").
#' @param nse Number of standard errors for error bars (default: 2).
#'
#' @return None. Creates a 2×2 panel plot.
#'
#' @details
#' This function creates a comprehensive visual summary for two-group comparisons
#' by displaying four complementary plots in a 2×2 layout:
#'
#' **1. Error bar plot** (top-left): Shows means and standard error bars for each group.
#' Uses \code{ebarplot} with nse standard errors.
#'
#' **2. Boxplots** (top-right): Standard boxplots showing medians, quartiles, and outliers.
#'
#' **3. Density estimates** (bottom-left): Overlaid kernel density estimates using
#' \code{g2plot} with adaptive kernel estimator.
#'
#' **4. Shift function** (bottom-right): Shows how quantiles differ between groups
#' using \code{sband} (shift band plot).
#'
#' **Input formats**:
#' - **Two vectors**: Provide x and y separately.
#' - **Matrix**: Provide x as a matrix with 2 columns (y should be NULL).
#' - **List**: Provide x as a list with 2 elements (y should be NULL).
#'
#' After plotting, the graphics parameters are reset to single-plot mode.
#'
#' @seealso \code{\link{ebarplot}}, \code{\link{g2plot}}, \code{\link{sband}}
#' @export
#' @examples
#' # Compare two groups
#' set.seed(123)
#' x <- rnorm(100, mean = 0, sd = 1)
#' y <- rnorm(100, mean = 0.5, sd = 1.2)
#' sumplot2g(x, y)
#'
#' # Using a matrix
#' data_mat <- cbind(x, y)
#' sumplot2g(data_mat)
sumplot2g<-function(x,y=NULL,xlab="X",ylab="",eblabx="Groups",eblaby="",nse=2){
#
# create four plots useful when comparing two groups
# 1. error bars
# 2. boxplots
# 3. kernel density estimates
# 4 shift function
#
if(!is.null(y)){
xy=list()
xy[[1]]=x
xy[[2]]=y
}
if(is.null(y)){
if(is.matrix(x))xy=matl(x)
}
par(mfrow=c(2,2))
par(oma=c(4,0,0,0))
ebarplot(xy,xlab=eblabx,ylab=eblaby,nse=nse)
boxplot(xy)
g2plot(xy[[1]],xy[[2]])
sband(xy[[1]],xy[[2]])
par(mfrow=c(1,1))
}

#' Plot Quantile Sums to Assess Symmetry About Zero
#'
#' @description
#' Creates a plot to assess whether a distribution is symmetric about zero by
#' plotting the sum of q and (1-q) quantiles. Useful for examining the symmetry
#' of difference scores or effect sizes.
#'
#' @param x Data vector, or matrix/data.frame with 2 columns for paired data.
#' @param y Optional second vector for paired data (default: NULL).
#' @param xlab Label for X-axis (default: "Quantile").
#' @param ylab Label for Y-axis (default: "Effect Size").
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function assesses distributional symmetry about zero by plotting the sum
#' of complementary quantiles: Q(q) + Q(1-q) for various values of q.
#'
#' **For paired data** (x and y provided, or x is a 2-column matrix): Computes
#' differences x - y and plots quantile sums.
#'
#' **For single sample** (only x provided): Directly plots quantile sums of x.
#'
#' **Interpretation**:
#' - **Symmetric distribution about zero**: Plot should be approximately a horizontal
#'   line at y = 0 (since Q(q) + Q(1-q) ≈ 0 when symmetric about 0).
#' - **Symmetric about median ≠ 0**: Plot should be approximately horizontal but
#'   not necessarily at y = 0.
#' - **Asymmetric distribution**: Plot will show curvature or trend.
#'
#' The function plots quantile sums for 9 quantiles (0.1 through 0.9) using
#' the Harrell-Davis estimator (\code{hd}), which provides smooth quantile estimates.
#'
#' Missing values are automatically removed.
#'
#' @seealso \code{\link{hd}}, \code{\link{sband}}, \code{\link{g2plotdifxy}}
#' @export
#' @examples
#' # Symmetric distribution about zero
#' set.seed(123)
#' x <- rnorm(100, mean = 0, sd = 1)
#' difQplot(x)
#'
#' # Skewed distribution
#' y <- rexp(100, rate = 1)
#' difQplot(y)
#'
#' # Paired data
#' pre <- rnorm(50, mean = 10, sd = 2)
#' post <- pre + rnorm(50, mean = 1, sd = 1)  # Symmetric shift
#' difQplot(pre, post)
difQplot<-function(x,y=NULL,xlab="Quantile",ylab="Effect Size"){
#
#  Plot that provides perspective on the degree a distribution is symmetric about zero.
#  This function plots the sum of q and 1-q quantiles. If the distributions are symmetric
#  the plot should be approximately a horizontal line. If in addition the median
# of the difference scores is zero, the horizontal line will intercept the y-axis at zero.
#
if(is.null(y))dif=x
if(!is.null(y))dif=x-y
x=elimna(x)
qd=NA
for(i in 1:99)qd[i]=hd(dif,.5-i/200)+hd(dif,.5+i/200)
plot(.5-c(1:99)/200,qd,xlab=xlab,ylab=ylab)
}

#' Split Data by Quantiles and Compare Groups
#'
#' @description
#' Splits the response variable into groups based on quantiles (or specified values)
#' of a predictor variable, then compares and plots the groups using a specified
#' comparison function.
#'
#' @param x Predictor variable (vector) used for splitting into groups.
#' @param y Response variable (vector) to be compared across groups.
#' @param q Quantiles of x to use for splitting (default: c(.25, .5, .75), creating
#'   4 groups based on quartiles). Ignored if vals is specified.
#' @param vals Specific values of x to use for splitting (default: NULL). If provided,
#'   overrides the q argument.
#' @param FUN Function to apply for group comparisons and plotting (default: lincon,
#'   linear contrasts). Can be any function that accepts grouped data.
#' @param ... Additional arguments passed to FUN.
#'
#' @return Output from the comparison function specified by FUN.
#'
#' @details
#' This function implements a common strategy in robust statistics: examining how
#' a response variable behaves across different levels of a predictor by:
#' 1. Splitting the predictor x into groups based on quantiles or specified values
#' 2. Grouping the response y accordingly
#' 3. Comparing the groups using the specified function (FUN)
#'
#' **Default behavior** (vals = NULL):
#' - Uses quantiles specified in q to create groups
#' - q = c(.25, .5, .75) creates 4 groups: (0-25%), (25-50%), (50-75%), (75-100%)
#'
#' **Custom split points** (vals specified):
#' - Groups are formed based on the specified values
#' - Useful for meaningful cutoffs (e.g., clinical thresholds)
#'
#' The grouped data is passed to FUN for comparison and plotting. Common choices:
#' - \code{lincon}: Linear contrasts with multiple comparisons
#' - \code{t1way}: One-way ANOVA with trimmed means
#' - \code{med1way}: One-way ANOVA with medians
#'
#' @seealso \code{\link{split.mat}}, \code{\link{lincon}}, \code{\link{t1way}}
#' @export
#' @examples
#' # Split data into quartiles and compare groups
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 2 * x + rnorm(100)
#' bplot(x, y, FUN = t1way)
#'
#' # Use specific cutoff values
#' bplot(x, y, vals = c(-1, 0, 1), FUN = lincon)
bplot<-function(x,y,q=c(.25,.5,.75),vals=NULL,FUN=lincon,...){
#
#
#  x is a vector
#
#  Split the data in y into groups based on values in x.
#  vals=NULL means quantiles of x will be used, quantiles indicated by argument
#  q
#
#  Next, compare and plot  the
# groups based using the method indicated by the argument
# FUN
#
if(!is.vector(x))stop('Argument x should be a vector')
v=split.mat(x,y,q=q,vals=vals)
a=FUN(v,...)
a
}

 outreg<-function(x,y,regfun=tsreg,outfun=outpro.depth,varfun=pbvar,corfun=pbcor,xout=FALSE,...){
#
# Do regression on points not labled outliers
# based on the  method indicated by outfun
#
#   A more general version of opreg
#
#  Possible alternative choices for outfun:
#  outpro
#  outmgv
#  outmcd
#
#
#  Note: argument xout is not relevant here, but is included to avoid conflicts when using regci.
#
x<-as.matrix(x)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
ivec<-outfun(m,plotit=FALSE)$keep
np1<-ncol(x)+1
coef<-regfun(m[ivec,1:ncol(x)],m[ivec,np1],...)$coef
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
# ERROR BAR/BOX PLOTS
# ============================================================================

#' Error Bar Plot for Group Comparisons
#'
#' @description
#' Creates an error bar plot showing group means (or trimmed means) with error bars
#' representing standard errors. Useful for visualizing group differences with
#' uncertainty estimates.
#'
#' @param x Data in matrix mode (groups as columns) or list mode (groups as elements).
#'   If providing two vectors, see y parameter.
#' @param y Optional second group (vector). If provided, x should also be a vector
#'   (default: NULL).
#' @param nse Number of standard errors for error bar width (default: 2). The error
#'   bars extend ±nse standard errors from the mean.
#' @param liw Lower inner width (passed to plotCI, default: uiw).
#' @param aui Upper confidence limit (default: NULL, computed from data).
#' @param ali Lower confidence limit (default: aui).
#' @param err Direction of error bars: "y" for vertical, "x" for horizontal (default: "y").
#' @param tr Trimming proportion for trimmed mean (default: 0 = regular mean).
#'   Use tr = 0.2 for 20% trimmed mean. Maximum is 0.5 (use ebarplot.med for medians).
#' @param ylim Y-axis limits (default: NULL, automatically determined).
#' @param sfrac Fraction of plot width for error bar caps (default: 0.01).
#' @param gap Gap around points where error bars are not drawn (default: 0).
#' @param add Logical. Add to existing plot if TRUE (default: FALSE).
#' @param col Color for plotting (default: par("col")).
#' @param lwd Line width (default: par("lwd")).
#' @param slty Line type for error bars (default: par("lty")).
#' @param xlab Label for X-axis (default: "Group").
#' @param ylab Label for Y-axis (default: NULL).
#' @param ... Additional graphical parameters passed to plotCI.
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function creates error bar plots for comparing multiple groups. Each group
#' is represented by a point (mean or trimmed mean) with error bars extending
#' ±nse standard errors.
#'
#' **Location measures**:
#' - **tr = 0** (default): Regular arithmetic mean
#' - **0 < tr < 0.5**: Trimmed mean (e.g., tr = 0.2 trims 20% from each tail)
#' - **tr = 0.5**: Median (but use \code{ebarplot.med} instead for proper CIs)
#'
#' **Input formats**:
#' - **Matrix**: Each column represents a group
#' - **List**: Each element represents a group
#' - **Two vectors**: Provide x and y separately
#'
#' **Error bars**: By default, extend ±2 standard errors from the mean. This
#' approximately corresponds to a 95% confidence interval for large samples.
#'
#' Missing values are automatically removed from each group.
#'
#' **Note**: For median-based error bars with distribution-free confidence intervals,
#' use \code{ebarplot.med} instead.
#'
#' @seealso \code{\link{ebarplot.med}}, \code{\link{plotCI}}, \code{\link{tmean}}
#' @export
#' @examples
#' # Error bar plot for three groups
#' set.seed(123)
#' x1 <- rnorm(30, mean = 10, sd = 2)
#' x2 <- rnorm(30, mean = 12, sd = 2)
#' x3 <- rnorm(30, mean = 14, sd = 3)
#' ebarplot(list(x1, x2, x3))
#'
#' # Using trimmed means
#' ebarplot(list(x1, x2, x3), tr = 0.2)
#'
#' # Two-group comparison
#' ebarplot(x1, x2)
ebarplot<-function(x,y=NULL,nse=2, liw = uiw, aui=NULL, ali=aui,
err="y", tr=0,ylim=NULL, sfrac = 0.01, gap=0, add=FALSE,
col=par("col"), lwd=par("lwd"), slty=par("lty"), xlab="Group",
                    ylab=NULL, ...) {
# plots error bars using the data in
# x, which is assumed to be a matrix with J columns (J groups) or
# x has list mode.
# nse indicates how many standard errors to use when plotting.
#
# By default, means are used. To use a trimmed mean, set
# tr to some value between 0 and .5
# So tr=.2 would use a 20% trimmed mean
#
# Missing values are automatically removed.
#
if(tr==.5)stop("For medians, use ebarplot.med")
if(!is.null(y)){
if(is.matrix(x))stop("When y is given, x should not be a matrix")
if(is.list(x))stop("When y is given, x should not be in list mode")
rem=x
x=list()
x[[1]]=rem
x[[2]]=y
}
if(is.matrix(x))x<-listm(x)
mval<-NA
if(!is.list(x) && is.null(y))stop("This function assumes there
 are  two or more groups")
for(j in 1:length(x))mval[j]<-mean(x[[j]],na.rm=TRUE,tr=tr)
se<-NA
#for(j in 1:length(x))se[j]<-sqrt(var(x[[j]],na.rm=TRUE)/length(x[[j]])
for(j in 1:length(x))se[j]<-trimse(x[[j]],na.rm=TRUE,tr=tr)
uiw<-nse*se
plotCI(mval,y=NULL, uiw=uiw, liw = uiw, aui=NULL, ali=aui,
                    err="y", ylim=NULL, sfrac = 0.01, gap=0, add=FALSE,
                    col=par("col"), lwd=par("lwd"), slty=par("lty"), xlab=xlab,
                    ylab=ylab)
}

#' Error Bar Plot for Medians with Distribution-Free Confidence Intervals
#'
#' @description
#' Creates an error bar plot showing group medians with distribution-free confidence
#' intervals. Unlike \code{ebarplot}, this function uses exact nonparametric confidence
#' intervals for medians rather than standard errors.
#'
#' @param x Data in matrix mode (groups as columns) or list mode (groups as elements).
#'   If providing two vectors, see y parameter.
#' @param y Optional second group (vector). If provided, x should also be a vector
#'   (default: NULL).
#' @param alpha Significance level for confidence intervals (default: 0.05 for 95% CIs).
#' @param nse Not used for this function (included for consistency with ebarplot).
#' @param liw Lower inner width (passed to plotCI, default: uiw).
#' @param aui Upper confidence limit (default: NULL, computed from data).
#' @param ali Lower confidence limit (default: aui).
#' @param err Direction of error bars: "y" for vertical, "x" for horizontal (default: "y").
#' @param tr Not used for this function (included for consistency with ebarplot).
#' @param ylim Y-axis limits (default: NULL, automatically determined).
#' @param sfrac Fraction of plot width for error bar caps (default: 0.01).
#' @param gap Gap around points where error bars are not drawn (default: 0).
#' @param add Logical. Add to existing plot if TRUE (default: FALSE).
#' @param col Color for plotting (default: par("col")).
#' @param lwd Line width (default: par("lwd")).
#' @param slty Line type for error bars (default: par("lty")).
#' @param xlab Label for X-axis (default: "Group").
#' @param ylab Label for Y-axis (default: NULL).
#' @param ... Additional graphical parameters passed to plotCI.
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function is specifically designed for median-based comparisons and uses
#' distribution-free (nonparametric) confidence intervals computed via \code{sint}.
#'
#' **Advantages over ebarplot with tr = 0.5**:
#' - Uses exact distribution-free confidence intervals
#' - No normality assumptions required
#' - Valid for small samples and non-normal distributions
#' - Confidence level controlled by alpha parameter
#'
#' **Confidence intervals**: Computed using the sign test approach (\code{sint}),
#' which provides exact coverage probabilities without distributional assumptions.
#'
#' **Input formats**:
#' - **Matrix**: Each column represents a group
#' - **List**: Each element represents a group
#' - **Two vectors**: Provide x and y separately
#'
#' Missing values are automatically removed from each group.
#'
#' **Note**: For means or trimmed means with standard error bars, use \code{ebarplot}.
#'
#' @seealso \code{\link{ebarplot}}, \code{\link{sint}}, \code{\link{plotCI}}
#' @export
#' @examples
#' # Median error bar plot for three groups
#' set.seed(123)
#' x1 <- rexp(30, rate = 1)  # Skewed data
#' x2 <- rexp(30, rate = 0.5)
#' x3 <- rexp(30, rate = 0.8)
#' ebarplot.med(list(x1, x2, x3))
#'
#' # Compare with regular error bars
#' par(mfrow = c(1, 2))
#' ebarplot(list(x1, x2, x3), main = "Mean ± 2SE")
#' ebarplot.med(list(x1, x2, x3), main = "Median with 95% CI")
ebarplot.med<-function(x,y=NULL,alpha=.05,nse=2, liw = uiw, aui=NULL, ali=aui,
err="y", tr=0,ylim=NULL, sfrac = 0.01, gap=0, add=FALSE,
col=par("col"), lwd=par("lwd"), slty=par("lty"), xlab="Group",
                    ylab=NULL, ...) {
# plots error bars using the data in
# x, which is assumed to be a matrix with J columns (J groups) or
# x has list mode.
# nse indicates how many standard errors to use when plotting.
#
# Designed specifically for medians
# Uses distribution-free confidence intervals
#
# Missing values are automatically removed.
#
if(!is.null(y)){
if(is.matrix(x))stop("When y is given, x should not be a matrix")
if(is.list(x))stop("When y is given, x should not be in list mode")
rem=x
x=list()
x[[1]]=rem
x[[2]]=y
}
if(is.matrix(x))x<-listm(x)
mval<-NA
if(!is.list(x) && is.null(y))stop("This function assumes there
 are  two or more groups")
aui=NA
ali=NA
for(j in 1:length(x)){
mval[j]<-median(x[[j]],na.rm=TRUE)
temp=sint(x[[j]],alpha=alpha,pr=FALSE)
ali[j]=temp[1]
aui[j]=temp[2]
}

plotCI(mval,y=NULL, liw = uiw, aui=aui, ali=ali,
                    err="y", ylim=NULL, sfrac = 0.01, gap=0, add=FALSE,
                    col=par("col"), lwd=par("lwd"), slty=par("lty"), xlab=xlab,
                    ylab=ylab)
}

#' Enhanced Boxplot with ggplot2
#'
#' @description
#' Creates a modern, enhanced boxplot using ggplot2 with viridis color scheme and
#' overlaid data points. Saves the plot to a file.
#'
#' @param x Data in matrix or data.frame format. Each column represents a group.
#' @param fileout Output filename for saving the plot (character string).
#'
#' @return None. Saves plot to specified file using \code{ggsave}.
#'
#' @details
#' This function creates publication-quality boxplots with the following features:
#' - **Viridis color scheme**: Color-blind friendly palette
#' - **Overlaid points**: Individual data points shown with jittering
#' - **Modern theme**: Uses \code{hrbrthemes::theme_ipsum()}
#'
#' The plot is automatically saved to the file specified by fileout.
#'
#' **Required packages**: reshape, tidyverse, viridis, ggplot2, hrbrthemes.
#' These packages must be installed before using this function.
#'
#' **Note**: The plot title is hard-coded as "Prediction Errors (D)" and axis
#' labels are blank. This function appears designed for specific use cases and
#' may need customization for general use.
#'
#' @seealso \code{\link{boxplot}}, \code{\link{STRIPchart}}
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' # Create enhanced boxplot
#' x <- matrix(rnorm(300), ncol = 3)
#' colnames(x) <- c("Group1", "Group2", "Group3")
#' box_plot1(x, "boxplot_output.png")
#' }
box_plot1<-function(x,fileout){
library("reshape")
library("tidyverse")
library("viridis")
library("ggplot2")
library("hrbrthemes")

x1<-melt(x)
#tiff(fileout)
#jpeg
plot<-ggplot(x1,aes(x=variable,y=value,fill=variable))+
 geom_boxplot() +
scale_fill_viridis(discrete = TRUE, alpha=0.6) +
geom_jitter(color="black", size=1.5, alpha=0.9) +
theme_ipsum() +
 theme(
   legend.position="none",
   plot.title = element_text(size=11)
  ) +
  ggtitle("             Prediction Errors (D)") +
 xlab("")+
 ylab("")
ggsave(fileout,plot=plot)

}

#' Strip Chart for Matrix Data
#'
#' @description
#' Creates a strip chart (one-dimensional scatterplot) that accepts matrix input.
#' This is a wrapper for \code{stripchart} that converts matrices to list format.
#'
#' @param x Data in matrix format. Each column represents a group. Will be converted
#'   to list mode for plotting.
#' @param method Method for handling overlapping points (default: "overplot"):
#'   \itemize{
#'     \item "overplot": Points plotted on top of each other
#'     \item "jitter": Add random noise to spread points
#'     \item "stack": Stack overlapping points
#'   }
#' @param jitter Amount of jittering when method = "jitter" (default: 0.1).
#' @param offset Stacking offset when method = "stack" (default: 1/3).
#' @param vertical Logical. Plot groups vertically if TRUE, horizontally if FALSE (default: FALSE).
#' @param group.names Labels for groups (optional).
#' @param add Logical. Add to existing plot if TRUE (default: FALSE).
#' @param at Numeric positions for groups (optional).
#' @param xlim X-axis limits (default: NULL, automatically determined).
#' @param ylim Y-axis limits (default: NULL, automatically determined).
#' @param ylab Label for Y-axis (default: NULL).
#' @param xlab Label for X-axis (default: NULL).
#' @param dlab Label for individual data axis (default: "").
#' @param glab Label for group axis (default: "").
#' @param log Specification for log scale: "", "x", "y", or "xy" (default: "").
#' @param pch Point character (default: 0).
#' @param col Color for points (default: par("fg")).
#' @param cex Character expansion factor for points (default: par("cex")).
#' @param axes Logical. Draw axes if TRUE (default: TRUE).
#' @param frame.plot Logical. Draw frame around plot if TRUE (default: axes).
#' @param ... Additional graphical parameters passed to stripchart.
#'
#' @return None. Creates a plot.
#'
#' @details
#' This function extends R's built-in \code{stripchart} to accept matrix input.
#' The standard \code{stripchart} accepts data frames or lists but not matrices,
#' so this wrapper converts matrices to list format before plotting.
#'
#' **Use cases**: Strip charts are particularly useful for:
#' - Small to moderate sample sizes
#' - Showing all individual data points
#' - Comparing distributions across groups
#' - Identifying outliers and patterns
#'
#' **Input format**: Provide data as a matrix where each column represents a group.
#' The function automatically converts to list mode using \code{listm}.
#'
#' @seealso \code{\link{stripchart}}, \code{\link{listm}}, \code{\link{boxplot}}
#' @export
#' @examples
#' # Strip chart for three groups
#' set.seed(123)
#' x <- matrix(rnorm(150), ncol = 3)
#' colnames(x) <- c("Group A", "Group B", "Group C")
#' STRIPchart(x, method = "jitter", vertical = TRUE)
#'
#' # Stack overlapping points
#' STRIPchart(x, method = "stack", vertical = TRUE)
STRIPchart<-function(x,method ='overplot', jitter = 0.1, offset = 1/3,
           vertical = FALSE, group.names, add = FALSE,
           at = NULL, xlim = NULL, ylim = NULL,
           ylab = NULL, xlab = NULL, dlab ='', glab ='',
           log = '', pch = 0, col = par('fg'), cex = par('cex'),
           axes = TRUE, frame.plot = axes, ...){
#
#    Same as stripchart,	only it	accepts	a matrix, unlike stripchart, which
#    allows x to be  a data frame or list mode,	but not	a matrix.
#
if(is.matrix(x))x=listm(x)
stripchart(x,method=method,jitter=jitter,offset = offset,
           vertical = vertical, group.names=group.names, add = add,
           at =at, xlim = xlim, ylim = ylim,
           ylab = ylab, xlab = xlab, dlab = dlab, glab = glab,
           log = log, pch = pch, col = col, cex = cex,
           axes = axes, frame.plot = frame.plot, ...)
}


# ============================================================================
# FUNCTIONAL DATA PLOTS
# ============================================================================

#' Functional Boxplot for Functional Data
#'
#' @description
#' Creates a functional boxplot (Sun & Genton method) for visualizing functional
#' data, identifying outlying curves, and displaying the central 50% region. Works
#' with both functional data objects and matrices.
#'
#' @param fit Functional data. Can be:
#'   \itemize{
#'     \item fd object (from fda package)
#'     \item fdSmooth or fdPar object
#'     \item Matrix (n × p: n curves evaluated at p time points)
#'   }
#' @param x Evaluation points for fd objects (default: NULL, uses 101 equally-spaced points).
#' @param method Depth measure for ordering curves (default: "MBD", Modified Band Depth):
#'   \itemize{
#'     \item "MBD": Modified Band Depth (recommended)
#'     \item "BD2": Band Depth with 2 curves
#'     \item "BD3": Band Depth with 3 curves
#'     \item "Both": Combination of BD2 and MBD
#'   }
#' @param depth Pre-computed depth values (default: NULL, computed from data).
#' @param plot Logical. Create plot if TRUE (default: TRUE).
#' @param prob Probability for central region (default: 0.5 for 50% region). Can be a vector
#'   for multiple regions.
#' @param color Color(s) for central region(s) (default: 6).
#' @param outliercol Color for outlier curves (default: 2, red).
#' @param barcol Color for whiskers and envelope (default: 4, blue).
#' @param fullout Logical. Show all outliers if TRUE (default: FALSE, outliers shown behind envelope).
#' @param factor Whisker extension factor (default: 1.5). Determines outlier threshold.
#' @param xlab Label for X-axis (default: "Time").
#' @param ylab Label for Y-axis (default: "Y").
#' @param ... Additional graphical parameters.
#'
#' @return List with components:
#'   \item{depth}{Depth values for all curves.}
#'   \item{outpoint}{Indices of curves identified as outliers.}
#'
#' @details
#' The functional boxplot extends the traditional boxplot to functional data by:
#' 1. Computing depth for each curve (most central curve has highest depth)
#' 2. Ordering curves by depth
#' 3. Identifying the central 50% region (or other probability)
#' 4. Computing envelope of central region
#' 5. Identifying outliers using whisker rule (factor * envelope height)
#'
#' **Depth measures**:
#' - **MBD** (Modified Band Depth): Recommended. Measures how often a curve is
#'   contained in bands formed by pairs of other curves.
#' - **BD2/BD3**: Standard band depth using 2 or 3 curves
#' - **Both**: Lexicographic combination (BD2 breaks MBD ties)
#'
#' **Outlier detection**: Curves exceeding the envelope ± factor × (envelope range)
#' are flagged as outliers and plotted with different color.
#'
#' **Plot elements**:
#' - Median curve (highest depth): Solid black line
#' - Central region: Shaded polygon
#' - Envelope: Maximum/minimum of non-outlying curves (blue lines)
#' - Whiskers: Vertical bars showing envelope extent (blue bars)
#' - Outliers: Dashed lines (red)
#'
#' Requires the \code{fda} package for functional data objects and depth computations.
#'
#' @references
#' Sun, Y. and Genton, M. G. (2011). Functional Boxplots. Journal of Computational
#' and Graphical Statistics, 20, 316-334.
#'
#' @seealso \code{\link{func.plot}}, \code{\link{Flplot}}, \code{\link{spag.plot}}
#' @export
#' @examples
#' \dontrun{
#' # Functional data example
#' library(fda)
#' # Create some functional data
#' t <- seq(0, 1, length = 51)
#' y <- matrix(0, nrow = 30, ncol = 51)
#' for(i in 1:30) {
#'   y[i,] <- sin(2*pi*t) + rnorm(51, 0, 0.1)
#' }
#' # Add an outlier curve
#' y[1,] <- sin(2*pi*t) + 2
#'
#' fbplot(t(y))
#' }
fbplot<-function(fit,x=NULL,method='MBD',depth=NULL,plot=TRUE,prob=0.5,color=6,outliercol=2,barcol=4,fullout=FALSE, factor=1.5,xlab='Time',ylab='Y',...){

  if(is.fdSmooth(fit) | is.fdPar(fit)){ fit = fit$fd }
	if(is.fd(fit)){
    if(length(x)==0){
      x = seq(fit$basis$rangeval[1],fit$basis$rangeval[2],len=101)
    }
    fit = eval.fd(x,fit)
  }

	tp=dim(fit)[1]
	n=dim(fit)[2]
	if (length(x)==0) {x=1:tp}
  #compute band depth
  if (length(depth)==0){
	if (method=='BD2') {depth=BD2(t(fit))}
	else if (method=='BD3') {depth=BD3(t(fit))}
	else if (method=='MBD') {depth=MBD(t(fit))}
	else if (method=='Both') {depth=round(BD2(t(fit)),4)*10000+MBD(t(fit))}
  }

	dp_s=sort(depth,decreasing=TRUE)
	index=order(depth,decreasing=TRUE)
	if (plot) {
	plot(x,fit[,index[1]],lty=1,lwd=2,col=1,type='l',xlab=xlab,ylab=ylab,...)
	}
	for (pp in 1:length(prob)){
		m=ceiling(n*prob[pp])#at least 50%
		center=fit[,index[1:m]]
		out=fit[,index[(m+1):n]]
		inf=apply(center,1,min)
		sup=apply(center,1,max)

		if (prob[pp]==0.5){ #check outliers
			dist=factor*(sup-inf)
			upper=sup+dist
			lower=inf-dist
			outly=(fit<=lower)+(fit>=upper)
			outcol=colSums(outly)
			remove=(outcol>0)
			#outlier column
			colum=1:n
			outpoint=colum[remove==1]
			out=fit[,remove]
			woout=fit
			good=woout[,(remove==0),drop=FALSE]
			maxcurve=apply(good,1,max)
			mincurve=apply(good,1,min)
			if (sum(outly)>0){
				if (plot) {
				matlines(x,out,lty=2,col=outliercol,type='l',...)
				}
			}
			barval=(x[1]+x[tp])/2
			bar=which(sort(c(x,barval))==barval)[1]
			if (plot) {
			lines(c(x[bar],x[bar]),c(maxcurve[bar],sup[bar]),col=barcol,lwd=2)
		    lines(c(x[bar],x[bar]),c(mincurve[bar],inf[bar]),col=barcol,lwd=2)
			}
		}
		xx=c(x,x[order(x,decreasing=TRUE)])
		supinv=sup[order(x,decreasing=TRUE)]
		yy=c(inf,supinv)
		if (plot) {
		if (prob[pp]==0.5) {polygon(xx,yy,col=color[pp],border=barcol,lwd=2)}
		else {polygon(xx,yy,col=color[pp],border=NA)}
		}
	}
	if (plot) {
	lines(x,fit[,index[1]],lty=1,lwd=2,col=1,type='l')
	lines(x,maxcurve,col=barcol,lwd=2)
	lines(x,mincurve,col=barcol,lwd=2)
	if (fullout) {
		if (sum(outly)>0){
				if (plot) {
				matlines(x,out,lty=2,col=outliercol,type='l',...)
				}
			}
		}
	}
	return(list(depth=depth,outpoint=outpoint))
}

#' Plot Average of Multiple Curves
#'
#' @description
#' Computes and plots the pointwise average (or other measure of central tendency)
#' of multiple functional curves evaluated at the same time points.
#'
#' @param x Matrix (n × p) where rows are curves and columns are time points.
#' @param est Function to compute measure of location at each time point (default: mean).
#'   Can be any function like median, tmean, etc.
#' @param xlab Label for X-axis (default: "Time").
#' @param ylab Label for Y-axis (default: "Y").
#' @param plotit Logical. Create plot if TRUE (default: TRUE).
#'
#' @return Vector of estimated values at each time point.
#'
#' @details
#' This function computes a summary curve from multiple functional observations by
#' applying a location estimator (default: mean) at each time point. Useful for
#' visualizing the average trend across multiple curves.
#'
#' **Procedure**:
#' 1. At each time point (column), compute the specified estimator across all curves (rows)
#' 2. Plot the resulting summary curve
#'
#' **Common estimators**:
#' - **mean**: Arithmetic mean (sensitive to outliers)
#' - **median**: Median (robust to outliers)
#' - **tmean**: Trimmed mean (robust, compromise between mean and median)
#' - **onestep**: One-step M-estimator
#'
#' @seealso \code{\link{FQplot}}, \code{\link{Flplot2g}}, \code{\link{fbplot}}
#' @export
#' @examples
#' # Average of multiple curves
#' set.seed(123)
#' t <- seq(0, 1, length = 50)
#' # Create 20 noisy sine curves
#' x <- matrix(0, nrow = 20, ncol = 50)
#' for(i in 1:20) {
#'   x[i,] <- sin(2*pi*t) + rnorm(50, 0, 0.2)
#' }
#' Flplot(x)
#'
#' # Use median instead of mean
#' Flplot(x, est = median)
Flplot<-function(x,est=mean,xlab='Time',ylab='Y',plotit=TRUE){
#
#  average n curves and plot results
#
es=apply(x,2,est)
if(plotit){
plot(es,xlab=xlab,ylab=ylab,type='n')
lines(es)
}
es
}

#' Plot Median and Quartiles of Multiple Curves
#'
#' @description
#' Computes and plots the pointwise median and quartiles (25th and 75th percentiles)
#' of multiple functional curves, providing a robust summary of central tendency
#' and variability.
#'
#' @param x Matrix (n × p) where rows are curves and columns are time points.
#' @param xlab Label for X-axis (default: "Time").
#' @param ylab Label for Y-axis (default: "Y").
#' @param plotit Logical. Create plot if TRUE (default: TRUE).
#'
#' @return Vector of median values at each time point.
#'
#' @details
#' This function provides a robust summary of functional data by computing the
#' Harrell-Davis estimator of the median and quartiles at each time point. This
#' gives a visualization of the central trend and spread of the curves.
#'
#' **Plot elements**:
#' - **Solid line**: Median (50th percentile) across curves at each time point
#' - **Dashed lines**: Lower (25th) and upper (75th) quartiles
#'
#' **Advantages over Flplot**:
#' - Shows variability (interquartile range) in addition to central tendency
#' - Uses robust Harrell-Davis quantile estimator (\code{hd})
#' - Less sensitive to outlying curves
#'
#' The quartile lines provide information about the spread: wide separation indicates
#' high variability across curves at that time point, while narrow separation indicates
#' curves are similar.
#'
#' @seealso \code{\link{Flplot}}, \code{\link{hd}}, \code{\link{fbplot}}
#' @export
#' @examples
#' # Median and quartiles of multiple curves
#' set.seed(123)
#' t <- seq(0, 1, length = 50)
#' # Create 20 noisy sine curves with varying amplitude
#' x <- matrix(0, nrow = 20, ncol = 50)
#' for(i in 1:20) {
#'   amp <- runif(1, 0.5, 1.5)
#'   x[i,] <- amp * sin(2*pi*t) + rnorm(50, 0, 0.1)
#' }
#' FQplot(x)
FQplot<-function(x,xlab='Time',ylab='Y',plotit=TRUE){
#
# Compute the  median and quartiles of  n curves and plot results
#
es=apply(x,2,hd)
es1=apply(x,2,hd,q=.25)
es2=apply(x,2,hd,q=.75)
if(plotit){
plot(rep(c(1:ncol(x)),3),c(es,es1,es2),xlab=xlab,ylab=ylab,type='n')
lines(es)
lines(es1,lty=2)
lines(es2,lty=2)
}
es
}

#' Plot Average Curves for Two Groups
#'
#' @description
#' Computes and plots pointwise averages (or other location measures) for two
#' independent groups of functional curves, allowing visual comparison of trends.
#'
#' @param x1 Matrix (n1 × p) for group 1, where rows are curves and columns are time points.
#' @param x2 Matrix (n2 × p) for group 2, where rows are curves and columns are time points.
#' @param est Function to compute measure of location at each time point (default: mean).
#' @param xlab Label for X-axis (default: "Time").
#' @param ylab Label for Y-axis (default: "Y").
#' @param plotit Logical. Create plot if TRUE (default: TRUE).
#'
#' @return List with components:
#'   \item{est.1}{Vector of estimated values for group 1 at each time point.}
#'   \item{est.2}{Vector of estimated values for group 2 at each time point.}
#'
#' @details
#' This function compares functional data from two independent groups by computing
#' and overlaying their summary curves. At each time point, the specified estimator
#' is computed separately for each group.
#'
#' **Plot appearance**:
#' - **Solid line**: Group 1 average curve
#' - **Dashed line**: Group 2 average curve
#'
#' **Input requirements**:
#' - Both matrices must have the same number of columns (time points)
#' - Missing values are automatically removed
#'
#' **Common estimators**:
#' - **mean**: Sensitive to outliers but efficient
#' - **median**: Robust to outlying curves
#' - **tmean**: 20% trimmed mean (robust compromise)
#'
#' This is particularly useful for comparing treatment vs. control groups in
#' longitudinal or functional data settings.
#'
#' @seealso \code{\link{Flplot}}, \code{\link{FQplot}}, \code{\link{lplot2g}}
#' @export
#' @examples
#' # Compare two groups of curves
#' set.seed(123)
#' t <- seq(0, 1, length = 50)
#'
#' # Group 1: baseline curves
#' x1 <- matrix(0, nrow = 15, ncol = 50)
#' for(i in 1:15) {
#'   x1[i,] <- sin(2*pi*t) + rnorm(50, 0, 0.2)
#' }
#'
#' # Group 2: shifted curves
#' x2 <- matrix(0, nrow = 15, ncol = 50)
#' for(i in 1:15) {
#'   x2[i,] <- sin(2*pi*t) + 0.5 + rnorm(50, 0, 0.2)
#' }
#'
#' Flplot2g(x1, x2)
Flplot2g<-function(x1,x2,est=mean,xlab='Time',ylab='Y',plotit=TRUE){
#
#  average n curves and plot results
#
x1=elimna(x1)
x2=elimna(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 should have the same number of columns')
x1=elimna(x1)
x2=elimna(x2)
es1=apply(x1,2,est)
es2=apply(x2,2,est)
if(plotit){
plot(rep(1:ncol(x1),2),c(es1,es2),xlab=xlab,ylab=ylab,type='n')
lines(es1)
lines(es2,lty=2)
}
list(est.1=es1,est.2=es2)
}

#' Functional Boxplot (Matrix Input)
#'
#' @description
#' Creates a functional boxplot for matrix data using the Sun & Genton method.
#' This is a wrapper for \code{fbplot} that transposes the matrix to the expected format.
#'
#' @param fit Matrix (n × p) where n = number of curves, p = number of time points.
#'   Rows with missing values are removed.
#' @param x Evaluation points (default: NULL, uses column indices).
#' @param method Depth measure (default: "MBD"). See \code{\link{fbplot}} for options.
#' @param depth Pre-computed depth values (default: NULL).
#' @param plotit Logical. Create plot if TRUE (default: TRUE).
#' @param prob Probability for central region (default: 0.5).
#' @param color Color for central region (default: 6).
#' @param outliercol Color for outliers (default: 2).
#' @param barcol Color for whiskers (default: 4).
#' @param fullout Logical. Show all outliers if TRUE (default: FALSE).
#' @param factor Whisker extension factor (default: 1.5).
#' @param xlim X-axis limits (default: based on data range).
#' @param ylim Y-axis limits (default: based on data range with padding).
#' @param xlab Label for X-axis (default: "Time").
#' @param ylab Label for Y-axis (default: "Y").
#' @param ... Additional graphical parameters.
#'
#' @return List with components:
#'   \item{depth}{Depth values for all curves.}
#'   \item{outpoint}{Indices of curves identified as outliers.}
#'
#' @details
#' This function is a convenience wrapper for \code{fbplot} that accepts matrix
#' input in the natural row-wise format (each row = one curve). The matrix is
#' transposed before passing to \code{fbplot}.
#'
#' **Input format**: Provide data as an n × p matrix where:
#' - **Rows**: Individual curves (subjects)
#' - **Columns**: Time points or evaluation points
#'
#' This is the opposite of \code{fbplot}, which expects columns to be curves.
#'
#' Requires the \code{fda} package. See \code{\link{fbplot}} for full details
#' on the functional boxplot methodology.
#'
#' @seealso \code{\link{fbplot}}, \code{\link{Flplot}}, \code{\link{spag.plot}}
#' @export
#' @examples
#' # Functional boxplot for matrix data
#' set.seed(123)
#' t <- seq(0, 1, length = 50)
#' # Create 20 curves
#' x <- matrix(0, nrow = 20, ncol = 50)
#' for(i in 1:20) {
#'   x[i,] <- sin(2*pi*t) + rnorm(50, 0, 0.1)
#' }
#' # Add outlier
#' x[1,] <- sin(2*pi*t) + 1.5
#'
#' func.plot(x)
func.plot<-function(fit, x = NULL, method ='MBD', depth = NULL, plotit = TRUE,
    prob = 0.5, color = 6, outliercol = 2, barcol = 4, fullout = FALSE,
    factor = 1.5, xlim = c(1, nrow(fit)), ylim = c(min(fit) -
        0.5 * diff(range(fit)), max(fit) + 0.5 * diff(range(fit))),xlab='Time',ylab='Y',
    ...){
#
# functional boxplot for functional data using
# method in Sun and Genton.
#
#
# fit is assumed to be an n-by-p matrix
# n= number of subjects
# p= number points where the function has been evaluated.
#
#  rows with missing values are automatically removed.
#
library(fda)
elimna(fit)
fit=t(fit)
res=fbplot(fit, x = NULL, method = method, depth = depth, plot = plotit,
    prob =prob, color =color, outliercol =outliercol, barcol = barcol,
fullout = fullout, factor = factor, xlim =xlim, ylim = ylim, xlab=xlab,ylab=ylab,...)
res
}

#' Spaghetti Plot for Longitudinal Data
#'
#' Creates a spaghetti plot for longitudinal data stored in a matrix. Each row
#' represents a subject, and each column contains measurements at different time
#' points. Optionally fits linear trends for each subject using a robust regression
#' estimator.
#'
#' @param x An n-by-p matrix where n = number of subjects and p = number of time
#'   points. Each row contains repeated measures for one subject.
#' @param regfun Regression estimator used for linear fit when \code{fit.lin=TRUE}.
#'   Default is \code{\link{tsreg}} (Theil-Sen).
#' @param type Type of plot (see \code{\link[graphics]{plot.default}}): 'l' for
#'   lines (default), 'p' for points, 'b' for both, 'o' for overplotted, 'c' for
#'   lines part of 'b'.
#' @param legend Logical; if TRUE, legend is added. Default is FALSE.
#' @param trace.label Character label for legend. Default uses the name of the
#'   trace factor.
#' @param fixed Logical; passed to \code{\link[graphics]{interaction.plot}}.
#' @param xlab Label for x-axis. Default is 'Time'.
#' @param ylab Label for y-axis. Default is empty string.
#' @param xtick Logical; if TRUE, x-axis tick marks are drawn. Default is FALSE.
#' @param xaxt Character; controls x-axis type (see \code{\link[graphics]{par}}).
#'   Default uses current par setting.
#' @param axes Logical; if TRUE, axes are drawn. Default is TRUE.
#' @param fit.lin Logical; if TRUE, linear fit is plotted for each subject instead
#'   of raw data. Default is FALSE.
#' @param ... Additional plotting parameters.
#'
#' @details
#' The function converts the matrix into long format suitable for
#' \code{\link[graphics]{interaction.plot}}. When \code{fit.lin=TRUE}, it fits
#' a regression line to each subject's data using the specified \code{regfun}
#' (default: Theil-Sen estimator) and plots these fitted values instead of the
#' raw data.
#'
#' @return
#' No return value; creates a plot as a side effect.
#'
#' @seealso
#' \code{\link[graphics]{interaction.plot}}, \code{\link{tsreg}}
#'
#' @examples
#' \dontrun{
#' # Create sample longitudinal data
#' set.seed(123)
#' n <- 20  # 20 subjects
#' p <- 5   # 5 time points
#' x <- matrix(rnorm(n * p, mean = rep(1:p, each = n)), nrow = n)
#'
#' # Basic spaghetti plot
#' spag.plot(x)
#'
#' # With linear fits using Theil-Sen estimator
#' spag.plot(x, fit.lin = TRUE)
#'
#' # With legend
#' spag.plot(x, legend = TRUE, trace.label = "Subject")
#' }
#'
#' @export
spag.plot<-function(x, regfun=tsreg,type = c('l',
    'p', 'b', 'o', 'c'), legend = FALSE, trace.label = deparse(substitute(trace.factor)),
    fixed = FALSE, xlab = 'Time', ylab ='',
    xtick = FALSE, xaxt = par('xaxt'), axes = TRUE, fit.lin=FALSE,...){
#
# Create a spaghetti plot for data stored in a matrix with
# n rows and p columns. The p columns
# contain  measures taken at p times for each subject.
# This function converts x into a form that can be used by interaction.plot
#
#  fit.line=TRUE means that a linear fit is plotted.
#
#  regfun: The linear fit is  based on the regression estimator indicated by
#          regfun. The  default is Theil--Sen estimator
#
#
# type: the type of plot (see plot.default): lines or points or both.
#
x=as.matrix(x)
n=nrow(x)
p=ncol(x)
np=n*p
m=matrix(NA,nrow=np,3)
pvec=c(1:p)
ic=1-p
iu=0
for(i in 1:n){
ic=ic+p
iu=iu+p
m[ic:iu,1]=i  # create Subject id.
m[ic:iu,2]=pvec
m[ic:iu,3]=x[i,]
}
if(!fit.lin)interaction.plot(m[,2],m[,1],m[,3],xlab=xlab,ylab=ylab,legend=legend,
xtick=xtick,xaxt=xaxt,axes=axes)
if(fit.lin){
fit=by(m[,2:3],m[,1],regYval,regfun=regfun)
fit1 <- unlist(fit)
names(fit1) <- NULL
#plotting the linear fit by id
interaction.plot(m[,2],m[,1], fit1,
                  xlab=xlab, ylab=ylab, legend=legend)
}
}


# ============================================================================
# INTERACTION PLOTS
# ============================================================================

#' Interaction Plot for Two-Way Design
#'
#' Creates an interaction plot for a two-way factorial design using location
#' estimators (e.g., means, medians, trimmed means). Useful for visualizing
#' interaction effects between two factors.
#'
#' @param J Number of levels for Factor 1.
#' @param K Number of levels for Factor 2.
#' @param x Data in matrix or list form. If \code{locvec=NULL}, the function
#'   computes location estimates using \code{locfun}.
#' @param locfun Location estimator function to apply to each group. Default is
#'   \code{\link[base]{mean}}. Other options: \code{\link{tmean}},
#'   \code{\link{median}}, etc.
#' @param locvec Optional vector of pre-computed location estimates for all J*K
#'   groups. If provided, \code{locfun} is ignored.
#' @param na.rm Logical; if TRUE, missing values are removed before computing
#'   location estimates. Default is TRUE.
#' @param g1lev Optional vector of custom labels for Factor 1 levels. Default is
#'   1, 2, ..., J.
#' @param g2lev Optional vector of custom labels for Factor 2 levels. Default is
#'   1, 2, ..., K.
#' @param type Type of plot: 'l' for lines (default), 'p' for points, 'b' for
#'   both.
#' @param xlab Label for x-axis (Factor 1). Default is 'Fac 1'.
#' @param ylab Label for y-axis. Default is 'means'.
#' @param trace.label Label for legend (Factor 2). Default is 'Fac 2'.
#' @param ... Additional arguments passed to
#'   \code{\link[graphics]{interaction.plot}}.
#'
#' @details
#' The function creates an interaction plot with Factor 1 on the x-axis and
#' separate lines for each level of Factor 2. If \code{locvec} is not provided,
#' the function computes location estimates for each group using \code{locfun}.
#' The data \code{x} should be organized as a J*K list or matrix where groups
#' are ordered by Factor 1 within Factor 2 (or vice versa depending on how
#' \code{locvec} is structured).
#'
#' @return
#' No return value; creates a plot as a side effect.
#'
#' @seealso
#' \code{\link[graphics]{interaction.plot}}, \code{\link{Qinterplot}},
#' \code{\link{plot.inter}}
#'
#' @examples
#' \dontrun{
#' # Create sample data for 2x3 design
#' set.seed(123)
#' J <- 2  # Factor 1 levels
#' K <- 3  # Factor 2 levels
#' x <- list()
#' for (i in 1:(J*K)) {
#'   x[[i]] <- rnorm(20, mean = i)
#' }
#'
#' # Interaction plot using means
#' interplot(J, K, x)
#'
#' # Using trimmed means
#' interplot(J, K, x, locfun = tmean)
#'
#' # With custom labels
#' interplot(J, K, x, g1lev = c("Low", "High"),
#'           g2lev = c("A", "B", "C"))
#' }
#'
#' @export
interplot<-function(J,K,x,locfun=mean,locvec=NULL,na.rm=TRUE,
g1lev=NULL,g2lev=NULL,type = c("l",
    "p", "b"), xlab = "Fac 1", ylab = "means",trace.label="Fac 2",...){
if(is.null(locvec))locvec=lloc(x,est=locfun,na.rm=na.rm)
if(is.list(locvec))locvec=as.vector(matl(locvec))
if(is.null(g1lev[1])){
g1=c(rep(1,K))
for(j in 2:J)g1=c(g1,rep(j,K))
}
if(!is.null(g1lev)){
g1=c(rep(g1lev[1],K))
for(j in 2:J)g1=c(g1,rep(g1lev[j],K))
}
g1=as.factor(g1)
if(is.null(g2lev[1]))g2=as.factor(rep(c(1:K),J))
if(!is.null(g2lev[1]))g2=as.factor(rep(g2lev,J))
g2=as.factor(g2)
interaction.plot(g1,g2,locvec, xlab = xlab, ylab = ylab,
trace.label=trace.label)
}

#' Quantile-Based Interaction Plot for 2x2 Design
#'
#' Creates an interaction plot for a 2x2 factorial design based on quantiles
#' estimated via the Harrell-Davis estimator. Useful for visualizing interactions
#' at specific quantiles (e.g., median) rather than means.
#'
#' @param x Data in matrix, data frame, or list form containing exactly 4 groups
#'   representing a 2x2 factorial design.
#' @param q Quantile to estimate. Default is 0.5 (median). Must be between 0 and 1.
#'
#' @details
#' This function is specifically designed for 2x2 factorial designs. It computes
#' the specified quantile for each of the four groups using the Harrell-Davis
#' estimator (\code{\link{hd}}) and creates an interaction plot. This is useful
#' when you want to examine interactions at specific quantiles rather than using
#' means or trimmed means.
#'
#' The Harrell-Davis estimator provides a smooth estimate of quantiles that is
#' more efficient than the traditional sample quantile when data are from a
#' continuous distribution.
#'
#' @return
#' No return value; creates a plot as a side effect.
#'
#' @seealso
#' \code{\link{hd}}, \code{\link{interplot}}, \code{\link{plot.inter}}
#'
#' @examples
#' \dontrun{
#' # Create sample data for 2x2 design
#' set.seed(123)
#' x <- list(
#'   rnorm(30, mean = 10),  # Group 1,1
#'   rnorm(30, mean = 12),  # Group 1,2
#'   rnorm(30, mean = 11),  # Group 2,1
#'   rnorm(30, mean = 15)   # Group 2,2
#' )
#'
#' # Median-based interaction plot
#' Qinterplot(x)
#'
#' # First quartile interaction plot
#' Qinterplot(x, q = 0.25)
#'
#' # Third quartile interaction plot
#' Qinterplot(x, q = 0.75)
#' }
#'
#' @export
Qinterplot<-function(x,q=.5){
#
# Plot interactions based on quantiles estimated via the
#  Harrell--Davis estimator
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
if(length(x)!=4)stop('Should have a 2 by 2 design for a total of four groups')
qv=lapply(x,hd,q=q)
qv=as.vector(matl(qv))
interplot(2,2,locvec=qv,xlab='Fac 1',ylab=paste(q,'Quantile'),trace.label='Fac 2')
}

#' Plot Distribution of Interaction Contrasts for 2x2 Design
#'
#' For a 2x2 factorial design, plots the bivariate distribution of two simple
#' effect contrasts: (X1 - X2) versus (X3 - X4). This provides a visual assessment
#' of the interaction pattern by showing the joint distribution of differences
#' across levels of one factor at each level of the other factor.
#'
#' @param x Data in matrix, data frame, or list form containing exactly 4 groups
#'   representing a 2x2 factorial design.
#' @param nreps Number of bootstrap resamples to generate for the distribution.
#'   Default is 500.
#' @param SEED Logical; if TRUE, sets random seed to 2 for reproducibility.
#'   Default is TRUE.
#' @param xlab Label for x-axis. Default is 'DV' (dependent variable).
#' @param ylab Label for y-axis. Default is empty string.
#'
#' @details
#' The function randomly samples one observation from each of the four groups
#' \code{nreps} times and computes L1 = X1 - X2 and L2 = X3 - X4 for each
#' resample. It then plots the bivariate distribution of (L1, L2) using
#' \code{\link{g2plot}}.
#'
#' This visualization helps assess interaction effects by showing:
#' \itemize{
#'   \item Whether the two simple effects tend to differ (interaction present)
#'   \item The joint variability of the two contrasts
#'   \item Potential outliers or unusual patterns
#' }
#'
#' @return
#' No return value; creates a plot as a side effect.
#'
#' @seealso
#' \code{\link{g2plot}}, \code{\link{interplot}}, \code{\link{Qinterplot}}
#'
#' @examples
#' \dontrun{
#' # Create sample data for 2x2 design with interaction
#' set.seed(123)
#' x <- list(
#'   rnorm(30, mean = 10),  # Group 1,1
#'   rnorm(30, mean = 12),  # Group 1,2 (difference = 2)
#'   rnorm(30, mean = 11),  # Group 2,1
#'   rnorm(30, mean = 18)   # Group 2,2 (difference = 7, interaction!)
#' )
#'
#' # Plot interaction contrasts
#' plot.inter(x)
#'
#' # More resamples for smoother distribution
#' plot.inter(x, nreps = 1000)
#' }
#'
#' @export
plot.inter<-function(x,nreps=500,SEED=TRUE,xlab='DV',ylab=''){
#
#  For a 2-by-2 design, plot the distribution of X_1-X_2 and X_3-X_4
#
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
if(J!=4)stop('Length should have four groups')
x=elimna(x)
L1=NA
L2=NA
linv=NA
for(i in 1:nreps){
for(j in 1:J)linv[j]=sample(x[[j]],1)
L1[i]=linv[1]-linv[2]
L2[i]=linv[3]-linv[4]
}
g2plot(L1,L2,xlab=xlab,ylab=ylab)
}

#' Plot Robust Regression Surface with Interaction Term
#'
#' Creates a 3D perspective plot of a regression surface with an interaction term
#' (product of two predictors) using a robust regression estimator. Useful for
#' visualizing nonlinear relationships and interaction effects.
#'
#' @param x An n-by-2 matrix of predictors (must have exactly 2 columns).
#' @param y Vector of response values (same length as rows in \code{x}).
#' @param regfun Robust regression function to use. Default is \code{\link{tsreg}}
#'   (Theil-Sen). Other options: \code{\link{opreg}}, \code{\link{ltsreg}}, etc.
#' @param pyhat Logical; currently not used (legacy parameter).
#' @param eout Logical; currently not used (legacy parameter).
#' @param xout Logical; if TRUE, removes outliers in the predictor space before
#'   fitting. Default is FALSE.
#' @param outfun Outlier detection function to use when \code{xout=TRUE}. Default
#'   is \code{\link{out}}.
#' @param plotit Logical; currently not used (legacy parameter).
#' @param expand Numeric expansion factor for z-axis. See
#'   \code{\link[graphics]{persp}}. Default is 0.5.
#' @param scale Logical; if TRUE (default), the z-axis is scaled. Recommended to
#'   try \code{scale=TRUE} if there is an association.
#' @param xlab Label for first predictor axis. Default is 'X'.
#' @param ylab Label for second predictor axis. Default is 'Y'.
#' @param zlab Label for response (z) axis. Default is empty string.
#' @param theta Viewing angle (azimuthal direction). Default is 50 degrees.
#' @param phi Viewing angle (colatitude). Default is 25 degrees.
#' @param family Currently not used (legacy parameter).
#' @param duplicate How to handle duplicate x-y locations. Default is 'error'.
#'   See \code{\link[akima]{interp}}.
#' @param ticktype Type of tick marks. Default is 'simple'. See
#'   \code{\link[graphics]{persp}}.
#' @param ... Additional arguments passed to \code{outfun}.
#'
#' @details
#' The function fits the regression model: Y = b0 + b1*X1 + b2*X2 + b3*X1*X2
#' using the specified robust regression estimator. The interaction term X1*X2
#' allows the relationship between Y and X1 to vary with X2 (and vice versa).
#'
#' The fitted surface is visualized using a 3D perspective plot via
#' \code{\link[graphics]{persp}} and \code{\link[akima]{interp}} for
#' interpolation. Duplicate points in the x-y plane are removed before plotting.
#'
#' @return
#' No return value; creates a 3D plot as a side effect.
#'
#' @seealso
#' \code{\link{ols.plot.inter}}, \code{\link{tsreg}}, \code{\link{regp2plot}},
#' \code{\link[graphics]{persp}}, \code{\link[akima]{interp}}
#'
#' @examples
#' \dontrun{
#' # Create sample data with interaction
#' set.seed(123)
#' n <- 100
#' x1 <- runif(n, -2, 2)
#' x2 <- runif(n, -2, 2)
#' x <- cbind(x1, x2)
#' y <- 1 + 2*x1 + 3*x2 + 1.5*x1*x2 + rnorm(n)
#'
#' # Plot regression surface with Theil-Sen estimator
#' reg.plot.inter(x, y)
#'
#' # Different viewing angles
#' reg.plot.inter(x, y, theta = 30, phi = 15)
#'
#' # Remove outliers before fitting
#' reg.plot.inter(x, y, xout = TRUE)
#' }
#'
#' @export
reg.plot.inter<-function(x,y, regfun=tsreg,
 pyhat = FALSE, eout = FALSE, xout = FALSE, outfun = out,
    plotit = TRUE, expand = 0.5, scale = TRUE, xlab = "X",
    ylab = "Y", zlab = "", theta = 50, phi = 25, family = "gaussian",
    duplicate = "error",ticktype="simple",...){
#
# Plot regression surface based on the classic interaction model:
#  usual product term
#
#   x is assumed to be a matrix with two columns (two predictors)
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
if(xout){
p=ncol(x)
p1=p+1
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}

if(!scale)print("scale=F. If there is an association, try scale=T")
if(ncol(x)!=2)stop("x should have two columns")
xx=cbind(x,x[,1]*x[,2])
temp=regfun(xx,y)
fitr=y-temp$residuals
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

#' Plot OLS Regression Surface with Interaction Term
#'
#' Creates a 3D perspective plot of an OLS (Ordinary Least Squares) regression
#' surface with an interaction term (product of two predictors). This is the
#' least squares version of \code{\link{reg.plot.inter}}.
#'
#' @param x An n-by-2 matrix of predictors (must have exactly 2 columns).
#' @param y Vector of response values (same length as rows in \code{x}).
#' @param pyhat Logical; currently not used (legacy parameter).
#' @param eout Logical; currently not used (legacy parameter).
#' @param xout Logical; if TRUE, removes outliers in the predictor space before
#'   fitting. Default is FALSE.
#' @param outfun Outlier detection function to use when \code{xout=TRUE}. Default
#'   is \code{\link{out}}.
#' @param plotit Logical; currently not used (legacy parameter).
#' @param expand Numeric expansion factor for z-axis. See
#'   \code{\link[graphics]{persp}}. Default is 0.5.
#' @param scale Logical; if TRUE (default), the z-axis is scaled.
#' @param xlab Label for first predictor axis. Default is 'X'.
#' @param ylab Label for second predictor axis. Default is 'Y'.
#' @param zlab Label for response (z) axis. Default is empty string.
#' @param theta Viewing angle (azimuthal direction). Default is 50 degrees.
#' @param phi Viewing angle (colatitude). Default is 25 degrees.
#' @param family Currently not used (legacy parameter).
#' @param duplicate How to handle duplicate x-y locations. Default is 'error'.
#'   See \code{\link[akima]{interp}}.
#' @param ticktype Type of tick marks. Default is 'simple'. See
#'   \code{\link[graphics]{persp}}.
#' @param ... Additional arguments passed to \code{outfun}.
#'
#' @details
#' The function fits the OLS regression model: Y = b0 + b1*X1 + b2*X2 + b3*X1*X2
#' using ordinary least squares (\code{\link[stats]{lsfit}}). The interaction
#' term X1*X2 allows the relationship between Y and X1 to vary with X2 (and
#' vice versa).
#'
#' The fitted surface is visualized using a 3D perspective plot via
#' \code{\link[graphics]{persp}} and \code{\link[akima]{interp}} for
#' interpolation. Duplicate points in the x-y plane are removed before plotting.
#'
#' For a robust alternative that is less sensitive to outliers, see
#' \code{\link{reg.plot.inter}}.
#'
#' @return
#' No return value; creates a 3D plot as a side effect.
#'
#' @seealso
#' \code{\link{reg.plot.inter}}, \code{\link{regp2plot}},
#' \code{\link[stats]{lsfit}}, \code{\link[graphics]{persp}},
#' \code{\link[akima]{interp}}
#'
#' @examples
#' \dontrun{
#' # Create sample data with interaction
#' set.seed(123)
#' n <- 100
#' x1 <- runif(n, -2, 2)
#' x2 <- runif(n, -2, 2)
#' x <- cbind(x1, x2)
#' y <- 1 + 2*x1 + 3*x2 + 1.5*x1*x2 + rnorm(n)
#'
#' # Plot OLS regression surface
#' ols.plot.inter(x, y)
#'
#' # Different viewing angles
#' ols.plot.inter(x, y, theta = 30, phi = 15)
#'
#' # Remove outliers before fitting
#' ols.plot.inter(x, y, xout = TRUE)
#' }
#'
#' @export
ols.plot.inter<-function(x,y, pyhat = FALSE, eout = FALSE, xout = FALSE, outfun = out,
    plotit = TRUE, expand = 0.5, scale = TRUE, xlab = "X",
    ylab = "Y", zlab = "", theta = 50, phi = 25, family = "gaussian",
    duplicate = "error",ticktype="simple",...){
#
# Plot regression surface based on the classic interaction model:
#  usual product term
#
#   x is assumed to be a matrix with two columns (two predictors)
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
if(ncol(x)!=2)stop("x should have two columns")
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:2]
y<-m[,3]
}
xx=cbind(x,x[,1]*x[,2])
temp=lsfit(xx,y)
fitr=y-temp$residuals
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


# ============================================================================
# LOGISTIC/LONGITUDINAL PLOTS
# ============================================================================

#' Plot Logistic Regression with MLE and/or Robust Estimates
#'
#' @description
#' Plots logistic regression curves for binary outcomes using maximum likelihood
#' estimation (MLE) and/or robust estimation. Supports one or two predictors,
#' creating either 2D line plots or 3D surface plots.
#'
#' @param x Numeric vector or matrix of predictor values. For matrices with
#'   more than one column, only the first 1-2 columns are used.
#' @param y Numeric vector of binary outcomes (0 or 1).
#' @param MLE Logical; if `TRUE`, plots the MLE logistic regression line (solid).
#' @param ROB Logical; if `TRUE`, plots the robust logistic regression line (dashed).
#' @param xlab,ylab,zlab Axis labels. Defaults: NULL (auto-generated), NULL, 'P(Z=1)'.
#' @inheritParams common-params
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param theta,phi Viewing angles for 3D plot (defaults: 50, 25).
#' @param duplicate How to handle duplicate points in 3D (default: "error").
#' @param LP Logical; if `TRUE`, uses LOESS for 3D plotting.
#' @param Lspan Span parameter for LOESS (default: 0.75).
#' @param pyhat Logical; if `TRUE`, returns predicted probabilities.
#' @param LABELS Logical; if `TRUE`, uses data labels for axis names.
#' @param WARN Logical; if `FALSE`, suppresses warnings from robust methods.
#' @param BY Logical; if `TRUE`, uses BY robust estimator; else weighted estimator.
#' @param expand,scale 3D plot parameters for \code{persp()} (defaults: 0.5, TRUE).
#' @param fr Span for LOESS (default: 2).
#' @param ticktype Tick type for 3D plot (default: "simple").
#' @param pr Logical; if `TRUE`, prints messages.
#'
#' @details
#' For **one predictor**:
#' - MLE fit plotted as solid line using \code{\link{logreg}}
#' - Robust fit plotted as dashed line using \code{wlogreg} or \code{BYlogreg}
#'
#' For **two predictors**:
#' - Creates 3D surface plot of predicted probabilities using LOESS
#' - Only works when `LP = TRUE`
#'
#' The function requires the `robustbase` package for robust estimation.
#'
#' @return If `pyhat = TRUE`, returns matrix of predicted probabilities;
#'   otherwise returns "Done".
#'
#' @seealso \code{\link{logreg}}, \code{\link{logreg.pred}}, \code{\link{lplot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Binary outcome data
#' x <- rnorm(100)
#' p <- exp(0.5 + 1.2*x) / (1 + exp(0.5 + 1.2*x))
#' y <- rbinom(100, 1, p)
#'
#' # Plot both MLE and robust fits
#' logreg.plot(x, y, MLE=TRUE, ROB=TRUE)
#' }
logreg.plot<-function(x,y,MLE=TRUE,ROB=FALSE,xlab=NULL,ylab=NULL,zlab='P(Z=1)',xout=FALSE,outfun=outpro,
theta=50,phi=25,duplicate="error",LP=TRUE,Lspan=.75,pyhat=FALSE,LABELS=FALSE,
WARN=FALSE,BY=TRUE,
expand=.5,scale=TRUE,fr=2,ticktype="simple",pr=TRUE,...){
#
# For one predictor, plot logistic regression line
#
#  if x is a matrix with more than one column, plot is  based on data in
#  in column 1.
#
#  MLE=T, will plot usual maximum likelihood estimate using a solid line
#  ROB=T, will plot robust estimate, which is indicated by a
#  dashed line.
#
library(robustbase)
xy=cbind(x,y)
xy=elimna(xy)
p1=ncol(xy)
if(p1>3)stop('Only one or two independent variables can be used')
if(!xout){
if(pr)print('Suggest also looking at result using xout=TRUE')
}
p=p1-1
x=xy[,1:p]
x=as.matrix(x)
y=xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
if(p==1){
if(is.null(ylab))ylab='P(Y=1|X)'
if(is.matrix(x))x=x[,1]
xord=order(x)
xx=x[xord]
yy=y[xord]
est1=logreg(xx,yy)[1:2,1]
if(is.null(xlab))v='X'
if(is.null(ylab))ylab='P(Y=1|X)'
if(LABELS)v=labels(x)[[2]]
if(MLE){
plot(xx,yy,xlab=v[1],ylab=ylab)
phat=logreg.pred(xx,yy,xx)
lines(xx,phat)
}
if(ROB){
if(!WARN)options(warn=-1)
if(!BY)est2=wlogreg(xx,yy)$coef[1:2]
if(BY)est2=BYlogreg(xx,yy)$coef[1:2]
phat2=exp(est2[1]+est2[2]*xx)/(1+exp(est2[1]+est2[2]*xx))
lines(xx,phat2,lty=2)
phat=cbind(xx,phat2)
dimnames(phat)=list(NULL,c(v,'Y.hat'))
if(!WARN)options(warn=0)
}
}
if(p==2){
fitr=logreg.pred(x,y,x)
if(is.null(xlab))v='X'
if(is.null(ylab))v[2]='Y'
if(LABELS)v=labels(x)[[2]]
if(LP)lplot(x,fitr,xlab=v[1],ylab=v[2],zlab=xlab,z=zlab,ticktype=ticktype,theta=theta,phi=phi,pr=FALSE)
phat=cbind(x,fitr)
dimnames(phat)=list(NULL,c(v,'Y.hat'))
}
if(!pyhat)phat<-"Done"
phat
}

#' Plot Individual Regression Lines for Longitudinal Data
#'
#' @description
#' Plots regression lines for longitudinal (repeated measures) data, with one
#' regression line fitted for each subject. Useful for visualizing individual
#' trajectories over time.
#'
#' @param x Data frame or matrix containing longitudinal data in long format.
#' @param x.col Column index indicating the predictor variable (typically time).
#' @param y.col Column index indicating the response variable.
#' @param s.id Column index indicating the subject ID.
#' @param regfun Regression function to use for fitting lines (default: \code{tsreg}
#'   for Theil-Sen regression). Can be any regression function that returns
#'   coefficients.
#' @param scat Logical; if `TRUE`, plots scatter points; if `FALSE`, plots only
#'   regression lines (default: TRUE).
#' @param xlab,ylab Axis labels (defaults: "X", "Y").
#'
#' @details
#' The function expects data in **long format** where each row represents one
#' observation, similar to the `Orthodont` dataset in the `nlme` package.
#'
#' For each subject:
#' 1. Extracts their observations from the data
#' 2. Fits a regression line using `regfun`
#' 3. Plots the fitted line
#'
#' **Note**: Currently only supports a single predictor (typically time or age).
#' The function will stop with an error if multiple predictors are specified.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect showing
#'   individual regression trajectories.
#'
#' @seealso \code{\link{tsreg}}, \code{\link{long2mat}}, \code{\link{longcov2mat}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Example with Orthodont data (requires nlme package)
#' library(nlme)
#' # Plot individual growth curves
#' # Column 2 is age, column 4 is distance, column 1 is subject ID
#' longreg.plot(Orthodont, x.col=2, y.col=4, s.id=1)
#' }
longreg.plot<-function(x,x.col,y.col,s.id,regfun=tsreg,scat=TRUE,xlab="X",
ylab="Y"){
#
# x is a data frame or matrix
#
# Longitudinal data: plot regression lines
#
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
if(p!=2)stop("Plot allows a single covariate only")
outmat=matrix(NA,nrow=n,ncol=p)
datx=NULL
daty=NULL
for(i in 1:n){
outmat[i,]=regfun(as.matrix(xvals[[i]]),ymat[i,])$coef
temp=as.matrix(xvals[[i]])
datx=c(datx,temp)
daty=c(daty,ymat[i,])
}
if(!scat)plot(datx,daty,type="n",xlab=xlab,ylab=ylab)
if(scat)plot(datx,daty,xlab=xlab,ylab=ylab)
for(i in 1:n)abline(outmat[i,1],outmat[i,2])
}


# ============================================================================
# DEPTH/BAGPLOT
# ============================================================================

#' Create a Bagplot for Bivariate Data
#'
#' @description
#' Creates a bagplot, which is a bivariate generalization of the boxplot. Uses
#' data depth to identify the "bag" (analogous to the box in a boxplot),
#' "loop" (fence), and outliers.
#'
#' @param x Matrix or data frame with 2 columns containing bivariate data.
#' @param plotit Logical; if `TRUE`, creates the plot (default: TRUE).
#' @param colorbag Color for the bag region (50\% depth contour). NULL uses default.
#' @param colorloop Color for the loop region (fence). NULL uses default.
#' @param colorchull Color for the convex hull. NULL uses default.
#' @param databag Logical; if `TRUE`, plots points in the bag region.
#' @param dataloop Logical; if `TRUE`, plots points in the loop region.
#' @param plot.fence Logical; if `TRUE`, plots the fence boundary.
#' @param type Character string specifying depth measure (default: 'hdepth'):
#'   \itemize{
#'     \item \code{'hdepth'}: Halfspace depth (Tukey depth)
#'     \item \code{'projdepth'}: Projection depth
#'     \item \code{'sprojdepth'}: Skewness-adjusted projection depth
#'   }
#'
#' @details
#' The bagplot extends the boxplot concept to two dimensions using data depth:
#' - **Bag**: Contains the 50\% deepest points (analogous to the IQR)
#' - **Loop**: Extends from the bag (analogous to the whiskers)
#' - **Outliers**: Points beyond the loop
#'
#' **Requires**: The `mrfDepth` and `ggplot2` packages must be installed.
#'
#' @return A bagplot object from the `mrfDepth` package. Creates a plot as
#'   a side effect.
#'
#' @references
#' Rousseeuw, P.J., Ruts, I., and Tukey, J.W. (1999). The bagplot: A bivariate
#' boxplot. *The American Statistician*, 53, 382-387.
#'
#' @seealso \code{\link{outpro}} for outlier detection using projection methods
#'
#' @export
#' @examples
#' \dontrun{
#' # Bivariate data with outliers
#' x <- cbind(rnorm(100), rnorm(100))
#' x <- rbind(x, c(5, 5))  # Add an outlier
#' Bagplot(x)
#' }
Bagplot<-function(x,plotit=TRUE,colorbag = NULL, colorloop = NULL,
    colorchull = NULL, databag = TRUE, dataloop = TRUE, plot.fence = FALSE,type='hdepth'){
#
#  requires packages mrfDepth and ggplot2
#
#   type = measure of depth: 'hdepth' =  halfspace depth,
#   'projdepth' for projection depth and
#   'sprojdepth'  for skewness-adjusted projection depth.
#
library(mrfDepth)
library(ggplot2)
z=compBagplot(x,type=type)
bagplot(z,colorbag =colorbag, colorloop = colorloop,
    colorchull =colorchull, databag=databag, dataloop =dataloop,
plot.fence = plot.fence)
}

#' Plot Expected Frequency Curves for Two Groups
#'
#' @description
#' Creates smoothed expected frequency curves for two independent groups,
#' useful for comparing distributions. The curves show the relative density
#' at each point based on local neighborhoods.
#'
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param fr Span parameter controlling smoothing (default: 0.8). Larger values
#'   produce smoother curves.
#' @param xlab,ylab Axis labels (defaults: empty strings).
#'
#' @details
#' For each value in a group, the expected frequency is computed as the
#' proportion of observations in a neighborhood defined by \code{fr}. The
#' function then plots smoothed frequency curves for both groups:
#' - Group 1 (x): solid line
#' - Group 2 (y): dashed line
#'
#' Higher curve values indicate regions of higher data density.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect.
#'
#' @seealso \code{\link{g2plot}}, \code{\link{splot}}, \code{\link{near}}
#'
#' @export
#' @examples
#' # Compare distributions of two groups
#' x <- rnorm(100, mean=0)
#' y <- rnorm(100, mean=1)
#' rd2plot(x, y)
rd2plot<-function(x,y,fr=.8,xlab="",ylab=""){
#
# Expected frequency curve
# for two groups.
#
# fr controls amount of smoothing
x<-elimna(x)
y<-elimna(y)
rmdx<-NA
rmdy<-NA
for(i in 1:length(x)){
rmdx[i]<-sum(near(x,x[i],fr))
}
for(i in 1:length(y)){
rmdy[i]<-sum(near(y,y[i],fr))
}
rmdx<-rmdx/length(x)
rmdy<-rmdy/length(y)
plot(c(x,y),c(rmdx,rmdy),type="n",ylab=ylab,xlab=xlab)
sx<-sort(x)
xorder<-order(x)
sysm<-rmdx[xorder]
lines(sx,sysm)
sy<-sort(y)
yorder<-order(y)
sysm<-rmdy[yorder]
lines(sy,sysm,lty=2)
}


# ============================================================================
# DISTRIBUTION/DENSITY PLOTS
# ============================================================================

#' Plot Relative Frequency Distribution for Discrete or Continuous Data
#'
#' @description
#' Creates a frequency plot showing the relative frequency of each unique
#' value in the data. Useful for visualizing discrete distributions or
#' grouped continuous data.
#'
#' @param x Numeric vector of data values.
#' @param op Logical; if `TRUE`, draws lines connecting frequencies (default: TRUE).
#' @param VL Logical; if `TRUE`, draws vertical lines for each value instead
#'   of connected lines (default: FALSE). Only used when `op = TRUE`.
#' @param xlab,ylab Axis labels (defaults: "X", "Rel. Freq.").
#' @param frame.plot Logical; if `TRUE`, draws a frame around the plot (default: TRUE).
#' @param plotit Logical; if `TRUE`, creates the plot (default: TRUE).
#'
#' @details
#' For each unique value in `x`, the function:
#' 1. Computes the frequency (count) and relative frequency (proportion)
#' 2. Plots asterisks at each unique value
#' 3. Optionally connects points with lines or vertical bars
#'
#' **Line types** (when `op = TRUE`):
#' - `VL = FALSE`: Connects frequencies with a line (default)
#' - `VL = TRUE`: Draws vertical lines from 0 to each frequency
#'
#' **No lines** when `op = FALSE` - only asterisks are plotted.
#'
#' @return List with components:
#' \item{obs.values}{Vector of unique observed values in `x`.}
#' \item{n}{Total sample size.}
#' \item{frequencies}{Count of each unique value.}
#' \item{rel.freq}{Relative frequency (proportion) of each unique value.}
#'
#' @seealso \code{\link{splotg5}} for plotting up to 5 groups simultaneously
#'
#' @export
#' @examples
#' # Discrete data
#' x <- sample(1:5, 100, replace=TRUE)
#' splot(x)
#'
#' # With vertical lines
#' splot(x, VL=TRUE)
#'
#' # Return frequency information
#' result <- splot(x, plotit=FALSE)
#' print(result$rel.freq)
splot<-function(x,op=TRUE,VL=FALSE,xlab="X",ylab="Rel. Freq.",frame.plot=TRUE,plotit=TRUE){
#
# Frequency plot
#
# For each unique value in x,
# the relatively frequency is determined and plotted.
#
# op=TRUE a line connecting the relative frequencies is drawn if VL=FALSE.
# VL=TRUE, a vertical line is drawn for each unique value in x;
# the height of the line indicates the relative frequency.
#
# op=FALSE. No lines are drawn
#
# The function returns the sample size as well as the frequencies
# associated with each unique value stored in x.
#
x<-x[!is.na(x)]
temp<-sort(unique(x))
freq<-NA
for(i in 1:length(temp)){
freq[i]<-sum(x==temp[i])
}
rmfreq=freq
nval=sum(freq)
freq<-freq/length(x)
tfreq<-freq
tfreq[1]<-0
tfreq[2]<-max(freq)
if(plotit){
plot(temp,tfreq,xlab=xlab,ylab=ylab,type="n",frame.plot=frame.plot)
points(temp,freq,pch="*")
if(op)
if(!VL)lines(temp,freq)
if(VL){
for(i in 1:length(temp))lines(c(temp[i],temp[i]),c(0,freq[i]))
}}
den=sum(rmfreq)
list(obs.values=temp,n=nval,frequencies=rmfreq,rel.freq=rmfreq/den)
}

#' Plot Relative Frequency Distributions for Up to Five Groups
#'
#' @description
#' Creates overlaid frequency plots for up to five groups, allowing visual
#' comparison of discrete distributions. Each group is plotted with a
#' different line type.
#'
#' @param x1 Numeric vector for group 1 (required).
#' @param x2 Numeric vector for group 2 (optional, default: NULL).
#' @param x3 Numeric vector for group 3 (optional, default: NULL).
#' @param x4 Numeric vector for group 4 (optional, default: NULL).
#' @param x5 Numeric vector for group 5 (optional, default: NULL).
#' @param xlab,ylab Axis labels (defaults: "X", "Rel. Freq.").
#'
#' @details
#' For each group provided, the function computes relative frequencies for
#' each unique value and plots them with different line types:
#' - Group 1: Line type 1 (solid)
#' - Group 2: Line type 2 (dashed)
#' - Group 3: Line type 3 (dotted)
#' - Group 4: Line type 4 (dot-dash)
#' - Group 5: Line type 5 (long-dash)
#'
#' All groups are overlaid on the same plot for easy comparison.
#'
#' @return List with components:
#' \item{freqx1}{Relative frequencies for group 1.}
#' \item{freqx2}{Relative frequencies for group 2 (if provided).}
#' \item{freqx3}{Relative frequencies for group 3 (if provided).}
#' \item{freqx4}{Relative frequencies for group 4 (if provided).}
#' \item{freqx5}{Relative frequencies for group 5 (if provided).}
#'
#' @seealso \code{\link{splot}} for single group frequency plots
#'
#' @export
#' @examples
#' # Compare frequency distributions of three groups
#' x1 <- sample(1:5, 100, replace=TRUE, prob=c(0.1,0.2,0.3,0.2,0.2))
#' x2 <- sample(1:5, 100, replace=TRUE, prob=c(0.3,0.2,0.2,0.2,0.1))
#' x3 <- sample(1:5, 100, replace=TRUE, prob=c(0.2,0.2,0.2,0.2,0.2))
#' splotg5(x1, x2, x3)
splotg5<-function(x1,x2=NULL,x3=NULL,x4=NULL,x5= NULL,xlab="X",ylab="Rel. Freq."){
#
# Frequency plot for up to five variables.
#
#
freqx2=NULL
freqx3=NULL
freqx4=NULL
freqx5=NULL
x1<-x1[!is.na(x1)]
x2<-x2[!is.na(x2)]
x3<-x3[!is.na(x3)]
x4<-x4[!is.na(x4)]
x5<-x5[!is.na(x5)]

xall=c(x1,x2,x3,x4,x5)
xall=xall[!is.na(xall)]
temp=sort(unique(xall))
XL=list(x1,x2,x3,x4,x5)
NN=0
for(j in 1:5)if(!is.null(XL[[j]]))NN=NN+1
freqx1<-NA
for(i in 1:length(temp)){
freqx1[i]<-sum(x1==temp[i])
}
freqx1<-freqx1/length(x1)
if(!is.null(x2)){
freqx2<-NA
for(i in 1:length(temp)){
freqx2[i]<-sum(x2==temp[i])
}
freqx2<-freqx2/length(x2)
}
if(!is.null(x3)){
freqx3<-NA
for(i in 1:length(temp)){
freqx3[i]<-sum(x3==temp[i])
}
freqx3<-freqx3/length(x3)
}
if(!is.null(x4)){
x4<-x4[!is.na(x4)]
freqx4<-NA
for(i in 1:length(temp)){
freqx4[i]<-sum(x4==temp[i])
}
freqx4<-freqx4/length(x4)
}
if(!is.null(x5)){
x5<-x5[!is.na(x5)]
freqx5<-NA
for(i in 1:length(temp)){
freqx5[i]<-sum(x5==temp[i])
}
freqx5<-freqx5/length(x5)
}
X=rep(temp,NN)
pts=c(freqx1,freqx2,freqx3,freqx4,freqx5)
plot(X,pts,type="n",xlab=xlab,ylab=ylab)
points(X,pts)
lines(temp,freqx1)
if(NN>=2)lines(temp,freqx2,lty=2)
if(NN>=3)lines(temp,freqx3,lty=3)
if(NN>=4)lines(temp,freqx4,lty=4)
if(NN>=5)lines(temp,freqx5,lty=5)
}

#' Plot Kernel Density Estimate
#'
#' @description
#' Computes and plots a kernel density estimate across a range of values,
#' providing a smooth estimate of the probability density function.
#'
#' @param x Numeric vector of observations.
#' @param rval Number of evaluation points for the density curve (default: 15).
#' @param xlab,ylab Axis labels (defaults: "X", "Y").
#'
#' @details
#' The function:
#' 1. Removes missing values from `x`
#' 2. Determines evaluation range from 1st to 99th percentiles
#' 3. Computes kernel density at `rval` equally-spaced points using \code{\link{kerden}}
#' 4. Plots the density curve
#'
#' Uses a relatively small number of evaluation points (default 15) for faster
#' computation compared to standard kernel density estimators.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect.
#'
#' @seealso \code{\link{kerden}} for the underlying kernel density computation
#'
#' @export
#' @examples
#' # Plot kernel density estimate
#' x <- rnorm(100)
#' kdplot(x)
#'
#' # Use more evaluation points for smoother curve
#' kdplot(x, rval=50)
kdplot<-function(x,rval=15,xlab="X",ylab="Y"){
#
#   Compute the kernel density estimator for a range of values
#   and plot results.
#
#   x contains vector of observations
#
x<-x[!is.na(x)]  #  Remove any missing values
y<-sort(x)
z<-1
temp<-floor(.01*length(x))
if(temp==0)temp<-5
ibot<-y[temp]
itop<-y[floor(.99*length(x))]
xaxis<-seq(ibot,itop,length=rval)
for(i in 1:rval)z[i]<-kerden(x,0,xaxis[i])
plot(xaxis,z,xlab=xlab,ylab=ylab)
lines(xaxis,z)
}

#' FY Plot with Prediction Intervals
#'
#' @description
#' Creates a fitted-versus-observed (FY) plot for OLS regression with
#' prediction intervals. Points outside the intervals are marked with
#' triangles, helping identify influential observations.
#'
#' @param x Matrix or data frame of predictor variables.
#' @param y Numeric vector of response values.
#' @param alpha Significance level for prediction intervals (default: 0.05).
#'
#' @details
#' The function:
#' 1. Fits OLS regression of `y` on `x`
#' 2. Computes fitted values and residuals
#' 3. Calculates prediction intervals adjusted for:
#'    - Leverage (hat values)
#'    - Small sample correction: `(1 + 15/n) * sqrt(n/(n-p))`
#' 4. Plots fitted vs observed values with:
#'    - Points: actual observations
#'    - Diagonal line: perfect fit
#'    - Triangles: upper and lower prediction limits
#'
#' Observations far from their prediction intervals may be outliers or
#' indicate model misspecification.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect.
#'
#' @seealso \code{\link{rplot.res}} for regression residual plots
#'
#' @export
#' @examples
#' # OLS regression with prediction intervals
#' x <- matrix(rnorm(100), ncol=2)
#' y <- x[,1] + 2*x[,2] + rnorm(50)
#' piplot(x, y)
piplot<-function(x, y, alpha = 0.05)
{
# Makes an FY plot with prediction limits added.
	x <- as.matrix(x)
	p <- dim(x)[2] + 1
	n <- length(y)
	up <- 1:n
	low <- up
	out <- lsfit(x, y)
	tem <- ls.diag(out)
	lev <- tem$hat
	res <- out$residuals
	FIT <- y - res
	Y <- y
	corfac <- (1 + 15/n)*sqrt(n/(n - p))
	val2 <- quantile(res, c(alpha/2, 1 - alpha/2))
	#get lower and upper PI limits for each case
	for(i in 1:n) {
		val <- sqrt(1 + lev[i])
		val3 <- as.single(corfac * val2[1] * val)
		val4 <- as.single(corfac * val2[2] * val)
		up[i] <- FIT[i] + val4
		low[i] <- FIT[i] + val3
	}
	zy <- c(min(low), Y, max(up))
	zx <- c(min(FIT), FIT, max(FIT))
        #change labels so plot labels are good
        ff <- FIT
        yy <- Y
        Y <- zy
        FIT <- zx
	plot(FIT, Y, type = "n")
	points(ff, yy)
	abline(0, 1)
	points(ff, up, pch = 17)
	points(ff, low, pch = 17)
}


# ============================================================================
# 3D PLOTTING
# ============================================================================

#' Create a 3D Scatter or Surface Plot
#'
#' @description
#' Creates a 3D perspective plot (surface plot) using interpolation. Useful
#' for visualizing regression surfaces or any trivariate relationship.
#'
#' @param x Matrix with 2 columns containing the two predictor variables
#'   (X1 and X2).
#' @param y Numeric vector of response values (Z).
#' @param xlab,ylab,zlab Axis labels (defaults: 'X1', 'X2', 'Y').
#' @param theta,phi Viewing angles for the 3D plot (defaults: 50, 25).
#'   `theta` controls azimuthal direction, `phi` controls colatitude.
#' @param duplicate How to handle duplicate (x,y) points (default: "error").
#'   Options: "error", "strip", "mean", etc. (passed to \code{akima::interp}).
#' @param pc Point character (currently unused, default: '*').
#' @param ticktype Tick mark type (default: 'simple').
#' @param expand Expansion factor for Z-axis (default: 0.5).
#'
#' @details
#' The function uses the `akima` package's \code{interp()} function to
#' interpolate irregularly-spaced (X1, X2, Z) data onto a regular grid,
#' then creates a 3D surface plot using \code{persp()}.
#'
#' **Note**: The `scale` variable is referenced but not defined as a parameter.
#' This may cause an error. Typically `scale = TRUE` is desired.
#'
#' @return Invisibly returns the perspective matrix from \code{persp()}.
#'   Creates a 3D plot as a side effect.
#'
#' @seealso \code{\link{regp2plot}}, \code{\link{reg2g.p2plot}} for regression-specific 3D plots
#'
#' @export
#' @examples
#' \dontrun{
#' # Create 3D surface from regression model
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 2*x1 + 3*x2 + rnorm(n, sd=0.5)
#' plot3D(cbind(x1, x2), y)
#' }
plot3D<-function(x,y,xlab='X1',ylab='X2',zlab='Y',theta=50,phi=25,
duplicate='error',pc='*',ticktype='simple',expand=.5){
#
#  A 3D plot: supplied for convenience
#
# Example: plot a regression surface
# x and y generated from regression model with no error term.
#
x=as.matrix(x)
if(ncol(x)!=2)stop('x should have two columns only')
fitr<-interp(x[,1],x[,2],y,duplicate=duplicate)
persp(fitr,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,expand=expand,
scale=scale,ticktype=ticktype)
}


# ============================================================================
# SPECIALIZED/UTILITY PLOTS
# ============================================================================

#' Plot Difference Score Distributions for J-by-2 Between-Within Design
#'
#' @description
#' Plots distributions of difference scores for a J-by-2 mixed design (J
#' independent groups, 2 dependent measures). Each group's difference
#' distribution is overlaid for comparison.
#'
#' @param J Number of independent groups (maximum 5).
#' @param K Number of dependent measures (must be 2).
#' @param x Data in matrix or list format. If matrix, organized as expected
#'   by \code{\link{bwimcp}} (columns alternate between dependent measures
#'   within each group). If list, each element contains data for one condition.
#' @param fr Smoothing parameter for density plots (default: 0.8).
#' @param aval Alpha level for adaptive kernel density (default: 0.5).
#' @param xlab,ylab Axis labels (defaults: 'X', '').
#' @param color Vector of length 5 specifying colors for each group (default:
#'   all black).
#' @param BOX Logical; if `TRUE`, creates boxplots instead of density plots
#'   (default: FALSE).
#'
#' @details
#' For each of the J groups, the function:
#' 1. Extracts data for the two dependent measures
#' 2. Computes difference scores (Measure 1 - Measure 2)
#' 3. Plots all J difference distributions on the same plot
#'
#' Uses \code{\link{g5plot}} for density plots (up to 5 groups) or
#' \code{boxplot} when `BOX = TRUE`.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect.
#'
#' @seealso \code{\link{bwimcp}}, \code{\link{g5plot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 3 groups, 2 time points
#' # Data in matrix format: columns are Time1_Grp1, Time2_Grp1, Time1_Grp2, ...
#' data_mat <- matrix(rnorm(200), ncol=6)
#' bwiJ2plot(J=3, K=2, x=data_mat)
#' }
bwiJ2plot<-function(J,K,x,fr=.8,aval=.5,xlab = 'X', ylab = '',
color = rep('black', 5),BOX=FALSE){
#
# This function is for a J by 2 between by within design
#
# Plot distribution of the difference scores for
# each of the J independpent groups
#
# x: can be a matrix, organized as expected by bwimcp
# or it can have list mode.
#
if(K!=2)stop('Should have only two dependent variables')
if(J>5)stop('Can only have five levels for the independent factor')
      if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}
dif=list()
for(j in 1:5)dif[[j]]=NULL
JK=J*K
m<-matrix(c(1:JK),J,K,byrow=TRUE)
ic<-c(-1,0)
for(j in 1:J){
ic<-ic+2
dif[[j]]=x[[ic[1]]]-x[[ic[2]]]
}
if(!BOX)g5plot(dif[[1]],dif[[2]],dif[[3]],dif[[4]],dif[[5]],fr = fr,
    aval = aval, xlab = xlab, ylab = ylab, color = color)
if(BOX)boxplot(dif)
}

#' Plot Distribution of Linear Contrast for Dependent Groups
#'
#' @description
#' For dependent (repeated measures) data, computes and plots the distribution
#' of a linear contrast across groups. Tests whether the contrast has median
#' zero and optionally tests for symmetry.
#'
#' @param x Matrix or data frame where rows are subjects and columns are
#'   groups/conditions, or a list where each element is a group's data.
#' @param con Numeric vector of contrast coefficients. Must sum to zero.
#' @param xlab,ylab Axis labels (defaults: 'DV', '').
#' @param sym.test Logical; if `TRUE`, tests whether the contrast distribution
#'   is symmetric about zero (default: FALSE).
#'
#' @details
#' For each subject i, computes the contrast: Y_i = sum(c_j * X_ij), where
#' c_j are the contrast coefficients in `con`.
#'
#' The function then:
#' 1. Plots the distribution using adaptive kernel density (\code{\link{akerd}})
#' 2. Tests H0: median = 0 using \code{\link{sintv2}}
#' 3. Computes quantile shift effect size (\code{\link{depQS}})
#' 4. Optionally tests for symmetry using \code{\link{Dqdif}}
#'
#' @return List with components:
#' \item{median}{Estimate of the median contrast.}
#' \item{n}{Sample size.}
#' \item{ci.low, ci.up}{Confidence interval for the median.}
#' \item{p.value}{P-value for testing median = 0.}
#' \item{Q.effect}{Quantile shift effect size.}
#' \item{sym.test}{Symmetry test result (if `sym.test = TRUE`); NULL otherwise.}
#'
#' @seealso \code{\link{linplot}}, \code{\link{sintv2}}, \code{\link{depQS}}, \code{\link{Dqdif}}
#'
#' @export
#' @examples
#' # Three repeated measures
#' x <- matrix(rnorm(150), ncol=3)
#' # Test if measure 1 differs from average of measures 2 and 3
#' dlinplot(x, con=c(1, -0.5, -0.5))
dlinplot<-function(x,con,xlab='DV',ylab='',sym.test=FALSE){
#
# For dependent variables,
# determine distribution of Y_i=sum_j c_jX_j
# and then plot the distribution
#
# The function also tests the hypothesis that Y	has a median of zero.
# sym.test=TRUE: test the hypothesis that Y is symmetric.
#
#  A quantile shift measure of effect size is returned as well.
#
if(is.matrix(con)){
if(ncol(con>1))print('Warning: Argument con should be a vector. Only the first contrast coefficients are used.')
}
con=as.vector(con)
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(is.data.frame(x))x=as.matrix(x)
if(is.list(x))x=matl(x)
x=elimna(x)
n=nrow(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
if(length(con)!=ncol(x))stop('Length of con should equal number of groups')
x=elimna(x)
L=NA
linv=NA
for(i in 1:n){
L[i]=sum(con*x[i,])
}
akerd(L,xlab=xlab,ylab=ylab)
mt=sintv2(L)
sym=NULL
Q=depQS(L)
if(sym.test)sym=Dqdif(L)
list(median=mt$median,n=mt$n,ci.low=mt$ci.low,ci.up=mt$ci.up,
p.value=mt$p.value,Q.effect=Q$Q.effect,sym.test=sym)
}

 dlin.sign<-function(x,con){
 # For dependent variables,
# determine distribution of Y_i=sum_j c_jX_{ij}
# and then do a sign test
con=as.vector(con)
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(is.data.frame(x))x=as.matrix(x)
if(is.list(x))x=matl(x)
x=elimna(x)
n=nrow(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
if(length(con)!=ncol(x))stop('Length of con should equal number of groups')
x=elimna(x)
L=NA
linv=NA
for(i in 1:n){
L[i]=sum(con*x[i,])
}
a=signt(dif=L)
list(Prob_a_value_is_less_than_zerro=a$Prob_x_less_than_y,ci=a$ci,n=a$n,N=a$N,p.value=a$p.value)
 }

#' Diagnostic Plot for Robust PCA
#'
#' @description
#' Creates diagnostic plots for robust principal component analysis (ROBPCA),
#' showing score distances and orthogonal distances to identify outliers.
#' Optionally includes classical PCA diagnostics for comparison.
#'
#' @param robpca.obj Object returned by a robust PCA function (e.g., from
#'   the `rrcov` or `robustbase` package).
#' @param classic Numeric; if 1, also plots classical PCA diagnostics for
#'   comparison (default: 0).
#' @param labod Number of observations with largest orthogonal distances to
#'   label (default: 3).
#' @param labsd Number of observations with largest score distances to label
#'   (default: 3).
#'
#' @details
#' The function creates diagnostic plots showing:
#' - **Score Distance (SD)**: Distance from the robust PCA center in the
#'   principal component space
#' - **Orthogonal Distance (OD)**: Distance from the PCA subspace
#'
#' **Two plot types**:
#' 1. If OD values vary, creates a 2D plot (SD vs OD) with cutoff lines
#' 2. If OD values are constant/tiny, creates an index plot of SD only
#'
#' Outliers appear as points beyond the cutoff lines. The most extreme
#' observations are labeled.
#'
#' If `classic = 1`, creates corresponding plots for classical (non-robust)
#' PCA titled "CPCA".
#'
#' @return Invisibly returns the input `robpca.obj`.
#'
#' @references
#' Hubert, M., Rousseeuw, P.J., and Vanden Branden, K. (2005). ROBPCA: A new
#' approach to robust principal component analysis. *Technometrics*, 47, 64-79.
#'
#' @seealso Functions in `rrcov` or `robustbase` packages for computing ROBPCA
#'
#' @export
#' @examples
#' \dontrun{
#' # Requires rrcov package
#' library(rrcov)
#' data(hbk)
#' res <- PcaHubert(hbk)
#' plot_robpca(res, classic=1)
#' }
plot_robpca<-function(robpca.obj, classic=0, labod=3, labsd=3) {
	diagnosticplot <- !any(robpca.obj$od <= as.vector(1.E-06,mode(robpca.obj$od)))
	if(diagnosticplot == T) {
		xmax <- max(max(robpca.obj$sd), robpca.obj$cutoff$sd)
		ymax <- max(max(robpca.obj$od), robpca.obj$cutoff$od)
		plot(robpca.obj$sd, robpca.obj$od, xlab="Score distance", ylab="Orthogonal distance", xlim=c(0,xmax), ylim=c(0,ymax), type="p")
		abline(v=robpca.obj$cutoff$sd)
		abline(h=robpca.obj$cutoff$od)
		givelabel(robpca.obj, labod, labsd)
	}
	else {
		ymax <- max(max(robpca.obj$sd), robpca.obj$cutoff$sd)
		plot(robpca.obj$sd, xlab="Index", ylab="Score distance", ylim=c(0,ymax), type="p")
		abline(h=robpca.obj$cutoff$sd)
		givelabel(robpca.obj, labod=0, labsd, indexplot=1)
	}
	title("ROBPCA")
	if(classic == 1) {
		diagnosticplot <- !any(robpca.obj$classic$od <= as.vector(1.E-06,mode(robpca.obj$classic$od)))
		if(diagnosticplot == T) {
			xmax <- max(max(robpca.obj$classic$sd), robpca.obj$classic$cutoff$sd)
			ymax <- max(max(robpca.obj$classic$od), robpca.obj$classic$cutoff$od)
			plot(robpca.obj$classic$sd, robpca.obj$classic$od, xlab="Score distance", ylab="Orthogonal distance", xlim=c(0,xmax), ylim=c(0,ymax), type="p")
			abline(v=robpca.obj$classic$cutoff$sd)
			abline(h=robpca.obj$classic$cutoff$od)
			givelabel(robpca.obj$classic, labod, labsd)
		}
		else {
			ymax <- max(max(robpca.obj$classic$sd), robpca.obj$classic$cutoff$sd)
			plot(robpca.obj$classic$sd, xlab="Index", ylab="Score distance", ylim=c(0,ymax), type="p")
			abline(h=robpca.obj$classic$cutoff$sd)
			givelabel(robpca.obj$classic, labod=0, labsd, indexplot=1)
		}
		title("CPCA")
	}
	invisible(robpca.obj)
}

#' Internal Helper for plot_robpca (Label Outliers)
#'
#' @description
#' Internal helper function used by \code{\link{plot_robpca}} to label
#' observations with large score distances or orthogonal distances.
#'
#' @param object ROBPCA object with `od` and `sd` components.
#' @param labod Number of observations with largest orthogonal distances to label.
#' @param labsd Number of observations with largest score distances to label.
#' @param indexplot If 1, uses index plot format; if 0, uses 2D scatter format.
#'
#' @return Invisibly returns the input object. Adds text labels to the current plot.
#'
#' @keywords internal
"givelabel"<-function(object, labod, labsd, indexplot=0) {
	if((labod == 0) && (labsd == 0)) {
		return(invisible(object))
	}
	if(indexplot != 1) {
		order.od <- order(object$od*(-1))
		order.sd <- order(object$sd*(-1))
		if(labod != 0) {
			for(i in 1:labod) {
				lab <- ifelse(is.character(names(object$od)), names(object$od[order.od[i]]), order.od[i])
				text(object$sd[order.od[i]], object$od[order.od[i]]+par("cxy")[2], labels=lab)
			}
		}
		if(labsd != 0) {
			for(i in 1:labsd) {
				lab <- ifelse(is.character(names(object$sd)), names(object$od[order.sd[i]]), order.sd[i])
				text(object$sd[order.sd[i]], object$od[order.sd[i]]+par("cxy")[2], labels=lab)
			}
		}
	}
	else {
		order.sd <- order(object$sd*(-1))
		if(labsd != 0) {
			for(i in 1:labsd) {
				lab <- ifelse(is.character(names(object$sd)), names(object$sd[order.sd[i]]), order.sd[i])
				text(order.sd[i], object$sd[order.sd[i]]+par("cxy")[2], labels=lab)
			}
		}
	}
	return(invisible(object))
}

#' Simple Test Plot (Utility)
#'
#' @description
#' Creates a simple scatter plot of two variables. Primarily for quick
#' visual inspection or testing purposes.
#'
#' @param x Matrix or data frame with at least 2 columns.
#'
#' @details
#' Plots column 1 on the x-axis and column 2 on the y-axis.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect.
#'
#' @keywords internal
#' @examples
#' # Quick scatter plot
#' x <- cbind(rnorm(50), rnorm(50))
#' testplot(x)
testplot<-function(x){
plot(x[,1],x[,2])
}


#' Double Angle Plot (DAP) for Ophthalmology Data
#'
#' @description
#' Creates a specialized "Double Angle Plot" for visualizing ophthalmology
#' astigmatism data. Plots raw data, convex polygons, and centroid on a
#' circular polar coordinate system.
#'
#' @param rawData Data frame with columns "X" and "Y" containing raw observations.
#' @param CRP Data frame with columns "X" and "Y" for the mean convex polygon.
#' @param center Data frame with columns "X" and "Y" for the centroid location.
#' @param DataMean Data frame with columns "X" and "Y" for the dataset convex polygon.
#' @param name Character string for the plot title.
#'
#' @details
#' The Double Angle Plot is used in ophthalmology to visualize astigmatism data
#' in polar coordinates. The function:
#' 1. Centers data around the origin
#' 2. Converts to polar coordinates (angle and radius)
#' 3. Plots concentric circles at radii 1, 2, 3, 4
#' 4. Overlays:
#'    - Raw data points (black dots)
#'    - Mean convex polygon (blue line)
#'    - Dataset convex polygon (purple line)
#'    - Centroid (red square)
#' 5. Labels angles at 0°, 22.5°, 45°, ..., 157.5°
#'
#' **Requires**: The `plotrix` package for drawing circles.
#'
#' @return Invisibly returns `NULL`. Creates a plot as a side effect.
#'
#' @references
#' Used for ophthalmology-specific astigmatism visualization. Related to
#' functions with prefix `oph.` in the WRS package.
#'
#' @export
#' @examples
#' \dontrun{
#' # Requires specific ophthalmology data structure
#' # See oph.* functions for data preparation
#' }
plotDAP<-function(rawData,CRP,center,DataMean,name){
library("plotrix")
#rawDat<-read.table("rawData.txt",sep="\t",header=T)
#realDat<-read.table("liwData.txt",sep="\t",header=F)
colnames(rawData)<-c("X","Y")
colnames(CRP)<-c("X","Y")
colnames(center)<-c("X","Y")
colnames(DataMean)<-c("X","Y")

D1<-data.frame(X1=CRP$X-mean(CRP$X),Y1=CRP$Y-mean(CRP$Y),X=CRP$X,Y=CRP$Y)
tmp1<-D1[D1$Y1>0,]
tmp2<-D1[D1$Y1<0,]
Tmp1<-tmp1[order(tmp1$X1,decreasing=T),]
Tmp2<-tmp2[order(tmp2$X1,decreasing=F),]
Tmp<-rbind(Tmp1,Tmp2)
Tmp$a<-atan2(Tmp$Y1,Tmp$X1)
Tmp<-Tmp[order(Tmp$a,decreasing=T),]
## add an ending point as the starting point, so that the circle is complete
newDat<-rbind(Tmp,Tmp[1,])

D2<-data.frame(X1=DataMean$X-mean(DataMean$X),Y1=DataMean$Y-mean(DataMean$Y),X=DataMean$X,Y=DataMean$Y)
tmp1<-D2[D2$Y1>0,]
tmp2<-D2[D2$Y1<0,]
Tmp1<-tmp1[order(tmp1$X1,decreasing=T),]
Tmp2<-tmp2[order(tmp2$X1,decreasing=F),]
Tmp<-rbind(Tmp1,Tmp2)
Tmp$a<-atan2(Tmp$Y1,Tmp$X1)
Tmp<-Tmp[order(Tmp$a,decreasing=T),]
## add an ending point as the starting point, so that the circle is complete
newDat2<-rbind(Tmp,Tmp[1,])

R<-4
cos45<-cos(pi/4)
#tiff("DoubleAnglePlot.tiff")
plot(newDat$X,newDat$Y,type="p",xlim=c(-5,5),ylim=c(-5,5),col="blue",pch=4,cex=0,frame.plot=F,axes=FALSE,xlab="",ylab="",asp=1,main=paste(name,""))
points(rawData$X,rawData$Y,type="p",col="black",pch=19,cex=0.5)
points(newDat2$X,newDat2$Y,type="p",col="purple",pch=8,cex=0)
points(center$X,center$Y,type="p",col="red",pch=15,cex=1)
lines(newDat$X,newDat$Y,col="blue",lwd=1.3)
lines(newDat2$X,newDat2$Y,col="purple",lwd=2.0)
### lines(D[,1],D[,2],type="p",col="red",pch=0,cex=5)

draw.circle(0,0,1,lty=1,lwd=0.5)
draw.circle(0,0,2,lty=1,lwd=0.5)
draw.circle(0,0,3,lty=1,lwd=0.5)
draw.circle(0,0,4,lty=1,lwd=0.5)
segments(-1*R*cos45,-1*R*cos45,1*R*cos45,R*cos45,lty=1,lwd=0.5)
segments(-1*R*cos45,1*R*cos45,1*R*cos45,-1*R*cos45,lty=1,lwd=0.5)
segments(0,R,0,-1*R,lty=1,lwd=0.5)
segments(R,0,-1*R,0,lty=1,lwd=0.5)

R1<-R+0.6
text(R1,0,paste0("0",intToUtf8(176)))
text(R1*cos45,R1*cos45,paste0("22.5",intToUtf8(176)))
text(0, R1,paste0("45",intToUtf8(176)))
text(-1*R1*cos45,R1*cos45,paste0("67.5",intToUtf8(176)))
text(-1*R1,0,paste0("90",intToUtf8(176)))
text(-1*R1*cos45,-1*R1*cos45,paste0("112.5",intToUtf8(176)))
text(0,-1*R1,paste0("135",intToUtf8(176)))
text(R1*cos45,-1*R1*cos45,paste0("157.5",intToUtf8(176)))
#dev.off()

#legend(x="bottomleft",pch=c(15,4,8),legend=c("Centroid","Mean Convex Polygon", "Dataset Convex Polygon"),col=c("red","blue","purple"))
legend(x="bottomleft",pch=c(15,NA,NA), lty=c(NA,1,1),cex=0.88,legend=c("Centroid","Mean Convex Polygon", "Dataset Convex Polygon"),col=c("red","blue","purple"))
}

