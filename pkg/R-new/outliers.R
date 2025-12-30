# WRS Package - Outlier Detection and Data Depth
# Methods for detecting outliers and computing data depth
#
# This module contains:
#   - Projection-based outlier detection: outpro variants
#   - Classical methods: outbox, outmah, outmve
#   - Modern robust methods: outogk, outDETMCD
#   - Data depth methods: depth variants, bagdepth
#   - Classification and bagging methods
#
# See: Wilcox (2022), Chapter on outlier detection


#' Projection-Based Outlier Detection (Multicore Version)
#'
#' Detects outliers using a modification of the Stahel-Donoho projection method
#' with parallel processing for improved performance. Projects data points onto
#' lines connecting each point to the center, then uses boxplot-based detection
#' on the projected distances.
#'
#' @param m A numeric matrix where rows represent observations and columns
#'   represent variables.
#' @param gval Critical value for outlier detection. Default is
#'   `sqrt(qchisq(0.975, p))` for most centers, or `sqrt(qchisq(0.95, p))` for
#'   cop=1, where p is the number of columns.
#' @param center Optional center point. If NA, computed based on `cop` argument.
#' @inheritParams common-params
#' @param op Logical. If `TRUE`, plots the 0.5 depth contour based on data with
#'   outliers removed. If `FALSE`, plots depth contour without removing outliers.
#' @param MM Logical. If `TRUE`, uses MAD for dispersion. If `FALSE`, uses
#'   interquartile range for outlier detection.
#' @param cop Integer (1-7) specifying method to compute center of data cloud:
#'   \itemize{
#'     \item 1: Deepest point (Tukey median)
#'     \item 2: MCD (Minimum Covariance Determinant) center
#'     \item 3: Coordinate-wise median (default)
#'     \item 4: MVE (Minimum Volume Ellipsoid) center
#'     \item 5: TBS center
#'     \item 6: RMBA (Olive's median ball algorithm)
#'     \item 7: Spatial (L1) median
#'   }
#' @param xlab Label for x-axis in bivariate plots.
#' @param ylab Label for y-axis in bivariate plots.
#' @param STAND Logical. If `TRUE`, standardizes marginal distributions before
#'   checking for outliers (recommended when variables are on different scales).
#' @inheritParams common-params
#' @param ... Additional arguments (not currently used).
#'
#' @return A list with components:
#'   \item{out.id}{Vector of row indices identified as outliers.}
#'   \item{keep}{Vector of row indices for non-outliers.}
#'
#' @details
#' This function is the parallel processing version of `outpro()`, using
#' `mclapply()` for improved performance on multicore systems. For bivariate
#' data with `plotit=TRUE`, creates a scatterplot marking outliers with "o" and
#' non-outliers with "*", and plots the 0.5 depth contour. The Donoho-Gasko
#' (Tukey) median is marked with "+".
#'
#' When using STAND=FALSE with variables on different scales, a warning is
#' printed suggesting STAND=TRUE may be more appropriate.
#'
#' @seealso \code{\link{outpro}}, \code{\link{outproad}}, \code{\link{outproadMC}}
#' @export
outproMC<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
xlab="VAR 1",ylab="VAR 2",STAND=TRUE,tr=.2,q=.5,pr=TRUE,...){
#
# same as function outpro, only it takes advantage of multiple core
# processors
#
# Detect outliers using a modification of the
# Stahel-Donoho  projection method.
#
# Determine center of data cloud, for each point,
# connect it with center, project points onto this line
# and use distances between projected points to detect
# outliers. A boxplot method is used on the
# projected distances.
#
# plotit=T creates a scatterplot when working with
# bivariate data.
#
# op=T
# means the .5 depth contour is plotted
# based on data with outliers removed.
#
# op=F
# means .5 depth contour is plotted without removing outliers.
#
#  MM=F  Use interquatile range when checking for outliers
#  MM=T  uses MAD.
#
#  If value for center is not specified,
#  there are four options for computing the center of the
#  cloud of points when computing projections:
#
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)#  cop=7 uses the spatial (L1) median
#
#  args q and tr having are not used by this function. They are included to deal
#  with situations where smoothers have optional arguments for q and tr
#
#  STAND=T means that marginal distributions are standardized before
#  checking for outliers

#  When using cop=2, 3 or 4, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
#  Donoho-Gasko (Tukey) median is marked with a cross, +.
#
m<-as.matrix(m)
if(pr){
if(!STAND){
if(ncol(m)>1)print('STAND=FALSE. If measures are on different scales, might want to use STAND=TRUE')
}}
if(ncol(m)==1){
dis<-(m-median(m))^2/mad(m)^2
dis<-sqrt(dis)
crit<-sqrt(qchisq(.975,1))
chk<-ifelse(dis>crit,1,0)
vec<-c(1:nrow(m))
outid<-vec[chk==1]
keep<-vec[chk==0]
}
if(ncol(m)>1){
if(STAND)m=standm(m,est=median,scat=mad)
if(is.na(gval) && cop==1)gval<-sqrt(qchisq(.95,ncol(m)))
if(is.na(gval) && cop!=1)gval<-sqrt(qchisq(.975,ncol(m)))
m<-elimna(m) # Remove missing values
if(cop==1 && is.na(center[1])){
if(ncol(m)>2)center<-dmean(m,tr=.5,cop=1)
if(ncol(m)==2){
tempd<-NA
for(i in 1:nrow(m))
tempd[i]<-depth(m[i,1],m[i,2],m)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center<-m[flag,]
if(sum(flag)>1)center<-apply(m[flag,],2,mean)
}}
if(cop==2 && is.na(center[1])){
center<-cov.mcd(m)$center
}
if(cop==4 && is.na(center[1])){
center<-cov.mve(m)$center
}
if(cop==3 && is.na(center[1])){
center<-apply(m,2,median)
}
if(cop==5 && is.na(center[1])){
center<-tbs(m)$center
}
if(cop==6 && is.na(center[1])){
center<-rmba(m)$center
}
if(cop==7 && is.na(center[1])){
center<-spat(m)
}
flag<-rep(0, nrow(m))
outid <- NA
vec <- c(1:nrow(m))
cenmat=matrix(rep(center,nrow(m)),ncol=ncol(m),byrow=TRUE)
Amat=m-cenmat
B=listm(t(Amat))  # so rows are now in B[[1]]...B[[n]]
dis=mclapply(B,outproMC.sub,Amat)
flag=mclapply(dis,outproMC.sub2,MM,gval)
flag=matl(flag)
flag=apply(flag,1,max)
}
if(sum(flag) == 0) outid <- NA
if(sum(flag) > 0)flag<-(flag==1)
outid <- vec[flag]
idv<-c(1:nrow(m))
keep<-idv[!flag]
if(ncol(m)==2){
if(plotit){
plot(m[,1],m[,2],type="n",xlab=xlab,ylab=ylab)
points(m[keep,1],m[keep,2],pch="*")
if(length(outid)>0)points(m[outid,1],m[outid,2],pch="o")
if(op){
tempd<-NA
keep<-keep[!is.na(keep)]
mm<-m[keep,]
for(i in 1:nrow(mm))tempd[i]<-depth(mm[i,1],mm[i,2],mm)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center<-mm[flag,]
if(sum(flag)>1)center<-apply(mm[flag,],2,mean)
m<-mm
}
points(center[1],center[2],pch="+")
x<-m
temp<-fdepth(m,plotit=FALSE)
flag<-(temp>=median(temp))
xx<-x[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}}
list(out.id=outid,keep=keep)
}




#' Helper Function for outproMC - Compute Projection Distances
#'
#' @param B Numeric vector representing a direction for projection.
#' @param Amat Matrix of centered data points.
#'
#' @return Vector of projection distances.
#' @keywords internal
outproMC.sub<-function(B,Amat){
dis<-NA
bot<-sum(B^2)
Bmat=matrix(rep(B,nrow(Amat)),ncol=ncol(Amat),byrow=TRUE)
temp<-apply(Bmat*Amat,1,sum)
temp=matrix(rep(temp,ncol(Amat)),ncol=ncol(Amat))
temp=temp*Bmat/bot
temp=temp^2
dis=apply(temp,1,sum)
dis<-sqrt(dis)
flag=(dis==Inf)
dis[flag]=NA
dis
}

#' Helper Function for outproMC - Identify Outliers from Distances
#'
#' @param dis Vector of projection distances.
#' @param MM Logical. If TRUE, uses MAD; if FALSE, uses IQR.
#' @param gval Critical value multiplier.
#'
#' @return Binary vector indicating outliers (1) and non-outliers (0).
#' @keywords internal
outproMC.sub2<-function(dis,MM,gval){
temp<-idealf(dis)
if(!MM)cu<-median(dis)+gval*(temp$qu-temp$ql)
if(MM)cu<-median(dis)+gval*mad(dis)
outid<-NA
temp2<-(dis> cu)
flag<-rep(0,length(dis))
flag[temp2]<-1
flag
}

#' Projection Outlier Detection with Adjusted Critical Value
#'
#' Adjusts the critical value used by `outpro()` so that the expected proportion
#' of points declared outliers under multivariate normality equals a specified
#' rate. This adjustment is particularly crucial for high-dimensional data (p>9).
#'
#' @param m A numeric matrix where rows are observations and columns are variables.
#' @param center Optional center point. If NA, computed based on `cop`.
#' @inheritParams common-params
#' @param op Logical. If TRUE, plots depth contour with outliers removed.
#' @param MM Logical. If TRUE, uses MAD; if FALSE, uses IQR.
#' @param cop Integer (1-7) specifying center computation method (see outproMC).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param rate Target false positive rate under normality (default: 0.05).
#' @param iter Number of iterations for simulation to adjust critical value
#'   (default: 100).
#' @param ip Number of quantile values to try when adjusting (default: 6).
#' @inheritParams common-params
#' @param STAND Logical. If TRUE, standardizes marginal distributions.
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers detected.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{used.gval}{Adjusted critical value used.}
#'
#' @details
#' The function simulates multivariate normal data and iteratively adjusts the
#' critical value `gval` until the empirical outlier rate matches the target
#' `rate`. This ensures appropriate Type I error control under normality.
#'
#' @seealso \code{\link{outpro}}, \code{\link{outproMC}}, \code{\link{outproadMC}}
#' @export
outproad<-function(m,center=NA,plotit=TRUE,op=TRUE,MM=TRUE,cop=3,
xlab="VAR 1",ylab="VAR 2",rate=.05,iter=100,ip=6,pr=TRUE,SEED=TRUE,STAND=TRUE){
m=elimna(m)
m=as.matrix(m)
n=nrow(m)
if(SEED)set.seed(2)
z=array(rmul(n*iter*ncol(m)),c(iter,n,ncol(m)))
newq=0
gtry=NA
for(itry in 1:ip){
newq=newq+9/10^itry
gtry[itry]=newq
}
gtry=c(.95,.975,gtry[-1])
if(pr)print("Computing adjustment")
for(itry in 1:ip){
val=NA
for(i in 1:iter){
temp=outpro(z[i,,],gval = sqrt(qchisq(gtry[itry],ncol(m))),
center=center,plotit=FALSE,op=op,MM=MM,cop=cop,STAND=STAND)$out.id
val[i]=length(temp)
}
erate=mean(val)/n
if(erate<rate){
newgval=sqrt(qchisq(gtry[itry],ncol(m)))
break
}}
res=outpro(m,gval=newgval,center=center,plotit=TRUE,op=op,MM=MM,
    cop = cop, xlab = "VAR 1", ylab = "VAR 2",STAND=STAND)
 list(n=res$n,n.out=res$n.out,out.id=res$out.id,keep=res$keep,used.gval=newgval)
}

# FOLLOWING CODE IS NO LONGER NEEDED but is retained in case it is desired to use the original version of rdepth

#' Projection Outlier Detection with Adjusted Critical Value (Multicore)
#'
#' Multicore version of `outproad()` that adjusts the critical value for
#' outlier detection to achieve a target false positive rate under normality.
#'
#' @inheritParams outproad
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers detected.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{used.gval}{Adjusted critical value used.}
#'
#' @details
#' This function uses `outproMC()` instead of `outpro()` for parallel processing.
#' Otherwise identical to `outproad()`.
#'
#' @seealso \code{\link{outproad}}, \code{\link{outproMC}}
#' @export
outproadMC<-function(m,center=NA,plotit=TRUE,op=TRUE,MM=TRUE,cop=3,
xlab="VAR 1",ylab="VAR 2",rate=.05,iter=100,ip=6,pr=TRUE,SEED=TRUE){
m=elimna(m)
m=as.matrix(m)
n=nrow(m)
z=array(rmul(n*iter*ncol(m)),c(iter,n,ncol(m)))
newq=0
gtry=NA
for(itry in 1:ip){
newq=newq+9/10^itry
gtry[itry]=newq
}
gtry=c(.95,.975,gtry[-1])
if(pr)print("Computing adjustment")
val=NA
if(SEED)set.seed(2)
for(itry in 1:ip){
for(i in 1:iter){
temp=outproMC(z[i,,],gval = sqrt(qchisq(gtry[itry],ncol(m))),
center=center,plotit=FALSE,op=op,MM=MM,cop=cop)$out.id
val[i]=length(temp)
}
erate=mean(val)/n
if(erate<rate){
newgval=sqrt(qchisq(gtry[itry],ncol(m)))
break
}}
res=outproMC(m,gval=newgval,center=center,plotit=TRUE,op=op,MM=MM,
    cop = cop, xlab = "VAR 1", ylab = "VAR 2")
 list(n=res$n,n.out=res$n.out,out.id=res$out.id,keep=res$keep,used.gval=newgval)
#list(results=res,used.gval=newgval)
}




#' Projection-Based Outlier Detection Using Depth
#'
#' Detects outliers using projection depth. For multivariate data, computes
#' projection depth for each point and identifies outliers based on inverse depth.
#' Designed to handle large sample sizes efficiently.
#'
#' @inheritParams common-params
#' @param ndir Number of random projections to use for computing depth
#'   (default: 1000).
#' @param MM Logical. If TRUE, uses MAD for dispersion in outlier detection.
#' @inheritParams common-params
#' @param xlab Label for x-axis in plots.
#' @param ylab Label for y-axis in plots.
#'
#' @return A list with the same structure as `outpro()`:
#'   \item{out.id}{Indices of detected outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'
#' @details
#' For univariate data, simply calls `outpro()`. For multivariate data, computes
#' projection depth using `prodepth()`, then applies outlier detection to the
#' inverse depths. For bivariate data with `plotit=TRUE`, creates a scatterplot
#' with outliers marked and the median depth contour plotted.
#'
#' @seealso \code{\link{outpro}}, \code{\link{prodepth}}
#' @export
outpro.depth<-function(x,ndir=1000,MM=FALSE,SEED=TRUE,plotit=FALSE,xlab='X',ylab='Y'){
x=elimna(x)
x=as.matrix(x)
if(ncol(x)==1)a=outpro(x)
else{
d=prodepth(x,ndir=ndir,SEED=SEED)
dis=1/d
a=outpro(dis,MM=MM)
}
if(ncol(x)==2){
if(plotit){
id.cen=which(d==max(d))
center=apply(x[id,],2,mean)
plot(x[,1],x[,2],type='n',xlab=xlab,ylab=ylab)
keep=a$keep
points(x[keep,1],x[keep,2],pch=".")
points(center[1],center[2],pch="+")
if(length(a$out.id)>0)points(x[a$out.id,1],x[a$out.id,2])
flag=which(d>=median(d))
xx<-x[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}}
a
}


#' Boxplot-Based Outlier Detection
#'
#' Detects outliers using the boxplot rule with ideal fourths for quartile
#' estimation. Optionally uses Carling's (2000) modified boxplot rule that
#' adjusts for sample size.
#'
#' @inheritParams common-params
#' @param mbox Logical. If TRUE, uses Carling's (2000) modified boxplot rule
#'   with sample-size-adjusted multiplier. If FALSE, uses standard 1.5*IQR rule.
#' @param gval Multiplier for the IQR. If NA, uses 1.5 for standard boxplot or
#'   Carling's formula for modified boxplot.
#' @inheritParams common-params
#' @param STAND Logical. Currently not used in this function.
#'
#' @return A list with components:
#'   \item{out.val}{Values of detected outliers.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers.}
#'   \item{cl}{Lower critical value.}
#'   \item{cu}{Upper critical value.}
#'
#' @details
#' Unlike R's `boxplot()`, this function uses ideal fourths (via `idealf()`)
#' to estimate quartiles. The modified boxplot (mbox=TRUE) uses the multiplier
#' `(17.63*n - 23.64)/(7.74*n - 3.71)` from Carling (2000), which provides
#' better control of the outlier detection rate for finite samples.
#'
#' @references
#' Carling, K. (2000). Resistant outlier rules and the non-Gaussian case.
#' Computational Statistics & Data Analysis, 33, 249-258.
#'
#' @seealso \code{\link{idealf}}, \code{\link{outpro}}
#' @export
outbox<-function(x,mbox=FALSE,gval=NA,plotit=FALSE,STAND=FALSE){
x<-x[!is.na(x)] # Remove missing values
if(plotit)boxplot(x)
n<-length(x)
temp<-idealf(x)
if(mbox){
if(is.na(gval))gval<-(17.63*n-23.64)/(7.74*n-3.71)
cl<-median(x)-gval*(temp$qu-temp$ql)
cu<-median(x)+gval*(temp$qu-temp$ql)
}
if(!mbox){
if(is.na(gval))gval<-1.5
cl<-temp$ql-gval*(temp$qu-temp$ql)
cu<-temp$qu+gval*(temp$qu-temp$ql)
}
flag<-NA
outid<-NA
vec<-c(1:n)
for(i in 1:n){
flag[i]<-(x[i]< cl || x[i]> cu)
}
if(sum(flag)==0)outid<-NULL
if(sum(flag)>0)outid<-vec[flag]
keep<-vec[!flag]
outval<-x[flag]
n.out=sum(length(outid))
list(out.val=outval,out.id=outid,keep=keep,n=n,n.out=n.out,cl=cl,cu=cu)
}


#' 3D Scatterplot with Outlier Highlighting
#'
#' Creates a 3D scatterplot and highlights outliers detected using a specified
#' outlier detection function. Optionally adds a regression plane.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @param xlab Label for x-axis (default: "Var 1").
#' @param ylab Label for y-axis (default: "Var 2").
#' @param zlab Label for z-axis (default: "Var 3").
#' @param reg.plane Logical. If TRUE, adds a regression plane where the first
#'   two columns of `x` are predictors and the third column is the outcome.
#' @param regfun Regression function to use for fitting the plane (default: tsreg).
#' @param COLOR Logical. If TRUE, outliers are shown in red; if FALSE, outliers
#'   are shown with "*" symbol.
#' @param tick.marks Logical. If TRUE, shows tick marks on axes.
#' @param ... Additional arguments passed to `outfun` and `regfun`.
#'
#' @return Invisibly returns NULL. Used for plotting side effects.
#'
#' @details
#' Requires the `scatterplot3d` package. Data must be a matrix or data frame
#' with exactly 3 columns. Outliers are detected using `outfun` and highlighted
#' either in color or with a special symbol.
#'
#' @seealso \code{\link{outpro}}, \code{\link[scatterplot3d]{scatterplot3d}}
#' @export
out3d<-function(x,outfun=outpro,xlab="Var 1",ylab="Var 2",zlab="Var 3",
reg.plane=FALSE,regfun=tsreg,COLOR=FALSE,tick.marks=TRUE,...){
if(!is.matrix(x) && !is.data.frame(x))stop("Data must be stored in a matrix
or data frame with 3 columns.")
if(ncol(x)!=3)stop("Data must be stored in a matrix with 3 columns.")
x=as.matrix(x)
x<-elimna(x)
library(scatterplot3d)
temp<-scatterplot3d(x,xlab=xlab,ylab=ylab,zlab=zlab,tick.marks=tick.marks)
outid<-outfun(x)$out.id
if(!is.na(outid[1])){
if(COLOR){
if(length(outid)==1)temp$points(t(as.matrix(x[outid,])),col="red")
if(length(outid)>1)temp$points(x[outid,],col="red")
}
if(!COLOR){
if(length(outid)==1)temp$points(t(as.matrix(x[outid,])),pch="*")
if(length(outid)>1)temp$points(x[outid,],pch="*")
}
}
if(reg.plane){
vals<-regfun(x[,1:2],x[,3],...)$coef
if(COLOR)temp$plane(vals,col="blue")
if(!COLOR)temp$plane(vals)
}
}


#' Bagplot-Based Outlier Detection
#'
#' Detects outliers using the bagplot method for bivariate data. The bagplot
#' is a bivariate generalization of the boxplot based on halfspace depth.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers detected.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'
#' @details
#' Requires the `mrfDepth` package. Only works with bivariate data (2 columns).
#' Uses `compBagplot()` from mrfDepth to identify outliers as points falling
#' outside the bagplot fence.
#'
#' @seealso \code{\link[mrfDepth]{compBagplot}}
#' @export
outbag<-function(x,plotit=FALSE){
library(mrfDepth)
x=elimna(x)
n=nrow(x)
z=compBagplot(x)
flag=z$datatype[,3]==3
if(sum(flag) == 0) outid <- NA
if(sum(flag) > 0)flag<-(flag==1)
idv<-c(1:n)
outid <- idv[flag]
keep<-idv[!flag]
n.out=length(outid)
list(n=n,n.out=n.out,out.id=outid,keep=keep)
}

 Rdepth<-function(x,y,z=NULL, ndir = NULL){
#
#
# z:
# An m by p+1 matrix containing row wise the hyperplanes for which to compute
# the regression depth. The first column should contain the intercepts.
# If z is not specified, it is set equal to cbind(x,y).
#
# Required: mrfDepth

# For convenience, the arguments correspond to conventions in WRS

x=cbind(x,y)
library(mrfDepth)
a=rdepth(x,z=z,ndir=ndir)
a
}



#' Mahalanobis Distance Outlier Detection
#'
#' Detects outliers using classical Mahalanobis distance. **For demonstration
#' purposes only** - subject to masking effects. Use robust methods like
#' `outpro()` or `outmve()` for actual analysis.
#'
#' @inheritParams common-params
#' @param qval Quantile for chi-squared critical value (default: pnorm(3) ≈ 0.9987,
#'   corresponding to 3 standard deviations in univariate case).
#' @inheritParams common-params
#' @param xlab Label for x-axis in bivariate plots.
#' @param ylab Label for y-axis in bivariate plots.
#'
#' @return A list with components:
#'   \item{out.val}{Values of detected outliers.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{dis}{Mahalanobis distances for all points.}
#'   \item{crit}{Critical value used.}
#'
#' @details
#' Uses classical mean and covariance matrix, making it vulnerable to masking
#' (outliers affecting the center/scatter estimates so they aren't detected).
#' For bivariate data with `plotit=TRUE`, creates a scatterplot with outliers
#' marked with "*".
#'
#' @seealso \code{\link{outmve}}, \code{\link{outpro}}, \code{\link{outDETMCD}}
#' @export
outmah<-function(x,qval=pnorm(3),plotit=TRUE,xlab="VAR 1",ylab="VAR 2"){
x=elimna(x)
x=as.matrix(x)
m=apply(x,2,mean)
v=cov(x)
dis=mahalanobis(x,m,v)
crit<-sqrt(qchisq(qval,ncol(x)))
vec<-c(1:nrow(x))
dis[is.na(dis)]=0
dis<-sqrt(dis)
chk<-ifelse(dis>crit,1,0)
id<-vec[chk==1]
keep<-vec[chk==0]
if(is.matrix(x)){
if(ncol(x)==2 && plotit){
plot(x[,1],x[,2],xlab=xlab,ylab=ylab,type="n")
flag<-rep(TRUE,nrow(x))
flag[id]<-FALSE
points(x[flag,1],x[flag,2])
if(sum(!flag)>0)points(x[!flag,1],x[!flag,2],pch="*")
}}
if(!is.matrix(x))outval<-x[id]
if(is.matrix(x))outval<-x[id,]
list(out.val=outval,out.id=id,keep=keep,dis=dis,crit=crit)
}


#' MVE/MCD-Based Outlier Detection
#'
#' Detects outliers using robust Mahalanobis distance based on either the
#' Minimum Volume Ellipsoid (MVE) or Minimum Covariance Determinant (MCD)
#' estimator of location and scatter.
#'
#' @inheritParams common-params
#' @param mve.flag Logical. If TRUE, uses MVE estimator; if FALSE, uses MCD
#'   estimator (default: TRUE).
#' @inheritParams common-params
#' @inheritParams common-params
#' @param outsym Symbol to use for plotting outliers (default: '*').
#'
#' @return A list with components:
#'   \item{out.id}{Indices of detected outliers.}
#'   \item{keep.id}{Indices of non-outliers.}
#'   \item{dis}{Robust Mahalanobis distances for all points.}
#'   \item{crit}{Critical value used (sqrt of 0.975 quantile of chi-squared).}
#'
#' @details
#' Uses `cov.mve()` or `cov.mcd()` from the MASS package to compute robust
#' center and covariance estimates, then flags points with robust Mahalanobis
#' distance exceeding the critical value. For bivariate data with `plotit=TRUE`,
#' creates a scatterplot with outliers marked.
#'
#' The seed is controlled to ensure reproducibility, then restored.
#'
#' @seealso \code{\link[MASS]{cov.mve}}, \code{\link[MASS]{cov.mcd}},
#'   \code{\link{outDETMCD}}, \code{\link{outpro}}
#' @export
outmve<-function(x,mve.flag=TRUE,plotit=TRUE,SEED=TRUE,outsym='*'){
if(SEED){
oldSeed <- .Random.seed
set.seed(12)
}
if(!is.matrix(x)){
x<-x[!is.na(x)]
dis<-mahalanobis(x,median(x),mad(x)^2)
crit<-sqrt(qchisq(.975,1))
vec<-c(1:length(x))
}
if(is.matrix(x)){
x<-elimna(x) # remove any missing values
if(mve.flag)mve<-cov.mve(x)
if(!mve.flag)mve<-cov.mcd(x)
dis<-mahalanobis(x,mve$center,mve$cov)
crit<-sqrt(qchisq(.975,ncol(x)))
vec<-c(1:nrow(x))
}
dis<-sqrt(dis)
chk<-ifelse(dis>crit,1,0)
id<-vec[chk==1]
keep<-vec[chk==0]
x<-as.matrix(x)
if(plotit && ncol(x)==2){
plot(x[,1],x[,2],xlab="X",ylab="Y",type="n")
flag<-rep(TRUE,nrow(x))
flag[id]<-FALSE
points(x[flag,1],x[flag,2])
if(sum(chk)!=0)points(x[!flag,1],x[!flag,2],pch=outsym)
}
if(SEED) {
    assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
}
list(out.id=id,keep.id=keep,dis=dis,crit=crit)
}


#' MGV-Based Outlier Detection
#'
#' Detects outliers using the MGV (Minimum Generalized Variance) method,
#' which is based on detecting points that inflate the generalized variance.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @param se Logical. If TRUE, standardizes variables using median and MAD
#'   before computing MGV (default: TRUE).
#' @inheritParams common-params
#' @param ndir Number of random projections for computing depth when plotting
#'   (default: 1000).
#' @param cov.fun Covariance function to use for computing center (default: rmba).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @inheritParams common-params
#' @param STAND Included for compatibility when called by other functions (not used).
#' @param ... Additional arguments passed to `outfun`.
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers detected.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'
#' @details
#' Computes MGV values using `mgvar()`, then applies `outfun` to detect outliers.
#' For bivariate data with `plotit=TRUE`, creates a scatterplot with the median
#' depth contour. For p>2, uses modified boxplot rule with adjusted critical value.
#'
#' **Note**: Results can be affected by rounding error in `eigen()` if columns
#' are reordered.
#'
#' @seealso \code{\link{mgvar}}, \code{\link{outbox}}, \code{\link{outmgvad}}
#' @export
outmgv<-function(x,y=NULL,plotit=TRUE,outfun=outbox,se=TRUE,op=1,ndir=1000,
cov.fun=rmba,xlab="X",ylab="Y",SEED=TRUE,STAND=FALSE,...){
if(is.null(y[1]))m<-x
if(!is.null(y[1]))m<-cbind(x,y)
m=elimna(m)
m=as.matrix(m)
nv=nrow(m)
temp<-mgvar(m,se=se,op=op,cov.fun=cov.fun,SEED=SEED)
temp[is.na(temp)]<-0
if(ncol(m)==1){
temp2=outpro(m)
nout=temp2$n.out
keep=temp2$keep
temp2=temp2$out.id
}
if(ncol(m)>1){
if(ncol(m)==2)temp2<-outfun(temp,...)
if(ncol(m)>2){
temp2<-outbox(temp,mbox=TRUE,gval=sqrt(qchisq(.975,ncol(m))))
}
if(plotit && ncol(m)==2){
x<-m[,1]
y<-m[,2]
plot(x,y,type="n",xlab=xlab,ylab=ylab)
points(x[temp2$keep],y[temp2$keep],pch="*")
if(!is.null(temp2$out.id))points(x[temp2$out.id],y[temp2$out.id],pch="o")

d=prodepth(m,ndir=ndir,SEED=SEED)
dis=1/d
id.cen=which(d==max(d))
if(length(id.cen)==1)center=m[id.cen,]
else
center=apply(m[id.cen,],2,mean)
points(center[1],center[2],pch="+")
flag=which(d>=median(d))
xx<-m[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
nout=0
if(!is.na(temp2[1]))nout=length(temp2$out.id)
}
list(n=nv,n.out=nout,out.id=temp2$out.id,keep=temp2$keep)
}


#' MGV Outlier Detection with Adjusted Critical Value
#'
#' Adjusts the critical value used by `outmgv()` to achieve a target false
#' positive rate under multivariate normality. Particularly important for
#' high-dimensional data (p>9).
#'
#' @param m A numeric matrix where rows are observations and columns are variables.
#' @param center Currently not used in the function.
#' @inheritParams common-params
#' @inheritParams common-params
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param rate Target false positive rate under normality (default: 0.05).
#' @param iter Number of simulation iterations for adjusting critical value
#'   (default: 100).
#' @param ip Number of quantile values to try (default: 6).
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers detected.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{used.gval}{Adjusted critical value used.}
#'
#' @details
#' Simulates multivariate normal data and iteratively adjusts the critical value
#' until the empirical outlier rate matches the target `rate`. This provides
#' better Type I error control under normality.
#'
#' @seealso \code{\link{outmgv}}, \code{\link{outmgv.v2}}
#' @export
outmgvad<-function(m,center=NA,plotit=TRUE,op=1,
xlab="VAR 1",ylab="VAR 2",rate=.05,iter=100,ip=6,pr=TRUE){
m=elimna(m)
n=nrow(m)
newgval=sqrt(qchisq(.975,ncol(m)))
z=array(rmul(n*iter*ncol(m)),c(iter,n,ncol(m)))
newq=0
gtry=NA
val=NULL
for(itry in 1:ip){
newq=newq+9/10^itry
gtry[itry]=newq
}
gtry=c(.95,.975,gtry[-1])
if(pr)print("Computing adjustment")
for(itry in 1:ip){
for(i in 1:iter){
temp=outmgv.v2(z[i,,],gval=gval,op=op)$out.id
val[i]=length(temp)
}
erate=mean(val)/n
if(erate<rate){
newgval=sqrt(qchisq(gtry[itry],ncol(m)))
break
}}
res=outmgv(m,gval=newgval,plotit=plotit,op=op, xlab = xlab, ylab = ylab)
#list(results=res,used.gval=newgval)
list(n=res$n,n.out=res$n.out,out.id=res$out.id,keep=res$keep,used.gval=newgval)
}



#' Fast MGV Outlier Detection (Inward Method)
#'
#' Faster version of `outmgv()` using the inward MGV method. Computes generalized
#' variance after removing each point, then applies outlier detection to these
#' deletion diagnostics.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @param se Logical. If TRUE, standardizes variables before analysis.
#' @param ndir Number of random projections for depth computation in plots
#'   (default: 1000).
#' @inheritParams common-params
#' @param ... Additional arguments passed to `outfun`.
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{out.val}{Values of outliers.}
#'   \item{depth.values}{Deletion diagnostic values for all points.}
#'
#' @details
#' Computes generalized variance with each point deleted (via `gvar()`), then
#' identifies outliers as points whose removal substantially changes the
#' generalized variance. Faster than `outmgv()` but based on different diagnostic.
#'
#' @seealso \code{\link{outmgv}}, \code{\link{gvar}}
#' @export
outmgvf<-function(x,y=NA,plotit=TRUE,outfun=outbox,se=TRUE,ndir=1000,SEED=TRUE,...){
if(is.na(y[1]))m<-x
if(!is.na(y[1]))m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing values
if(se){
for(i in 1:ncol(m))m[,i]<-(m[,i]-median(m[,i]))/mad(m[,i])
}
iflag<-rep(TRUE,nrow(m))
dval<-0
for(i in 1:nrow(m)){
dval[i]<-gvar(m[-i,])
}
temp2<-outfun(dval,...)
if(plotit && ncol(m)==2){
flag=which(dval<=median(dval))
x<-m[,1]
y<-m[,2]
plot(x,y,type="n",xlab="X",ylab="Y")
points(x[temp2$keep],y[temp2$keep],pch='*')
d=prodepth(m,ndir=ndir,SEED=SEED)
dis=1/d
id.cen=which(d==max(d))
center=apply(m[id,],2,mean)
points(center[1],center[2],pch="+")
flag=which(d>=median(d))
xx<-m[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
if(!is.null(temp2$out.id))points(x[temp2$out.id],y[temp2$out.id],pch="o")
}
list(n=temp2$n,out.id=temp2$out.id,keep=temp2$keep,out.val=m[temp2$out.id,],depth.values=dval)
}


#' MGV Outlier Detection (Version 2)
#'
#' Alternative version of MGV outlier detection with explicit control over the
#' critical value. Used internally by `outmgvad()`.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @param se Logical. If TRUE, standardizes variables.
#' @inheritParams common-params
#' @param gval Critical value (default: sqrt of 0.975 quantile of chi-squared
#'   with p degrees of freedom).
#' @param cov.fun Covariance function for computing center (default: rmba).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @inheritParams common-params
#' @param ... Additional arguments passed to `outfun`.
#'
#' @return A list with components:
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'
#' @keywords internal
outmgv.v2<-function(x,outfun=outbox,se=TRUE,op=1,
gval=sqrt(qchisq(.975,ncol(x))),
cov.fun=rmba,xlab="X",ylab="Y",SEED=TRUE,...){
m<-x
m=elimna(m)
temp<-mgvar(m,se=se,op=op,cov.fun=cov.fun,SEED=SEED)
temp[is.na(temp)]<-0
temp2<-outbox(temp,mbox=TRUE,gval=sqrt(qchisq(.975,ncol(m))))$out.id
vec<-c(1:nrow(m))
flag<-rep(TRUE,nrow(m))
flag[temp2]<-FALSE
vec<-vec[flag]
vals<-c(1:nrow(m))
keep<-vals[flag]
list(out.id=temp2,keep=keep)
}


#' Mean/SD-Based Outlier Detection
#'
#' Detects univariate outliers using the classical mean and standard deviation
#' approach (z-scores). **For demonstration purposes only** - not robust to outliers.
#'
#' @inheritParams common-params
#' @param crit Critical value for z-scores (default: 2, meaning points more than
#'   2 standard deviations from the mean are flagged).
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers.}
#'   \item{out.value}{Values of outliers.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'
#' @details
#' Only works with univariate data. Uses classical z-scores which are not robust
#' (outliers affect the mean and SD). For robust univariate outlier detection,
#' use `outbox()` or `outpro()`.
#'
#' @seealso \code{\link{outbox}}, \code{\link{outpro}}
#' @export
outms<-function(x,crit=2,plotit=FALSE){
x=elimna(x)
x=as.matrix(x)
if(ncol(x)==1){
z=(x-mean(x))/sd(x)
flag=abs(z)>=crit
out.id=z[flag]
n.out=sum(flag)
nums=c(1:length(x))
keep=nums[!flag]
}
if(ncol(x)>1)stop('Use function	out with outfun=wmean.cov')
list(n=length(x),n.out=n.out,out.value=x[flag],out.id=nums[flag],keep=keep)
}

#' OGK-Based Outlier Detection
#'
#' Detects outliers using the Orthogonalized Gnanadesikan-Kettenring (OGK)
#' robust covariance estimator.
#'
#' @inheritParams common-params
#' @param sigmamu Univariate scale estimator (default: taulc).
#' @param v Bivariate covariance estimator (default: gkcov).
#' @inheritParams common-params
#' @inheritParams common-params
#' @param beta Consistency factor for OGK (default: max(0.95, min(0.99, 1/n+0.94))).
#' @param n.iter Number of iterations for OGK (default: 1).
#' @inheritParams common-params
#' @param ... Additional arguments passed to OGK functions.
#'
#' @return A list with components:
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{distances}{Robust Mahalanobis distances.}
#'
#' @details
#' If `op=TRUE` (recommended), uses robust Mahalanobis distance with the OGK
#' estimator and beta adjusted for approximately 5% outlier rate under normality.
#' If `op=FALSE`, uses the distances from the OGK estimator directly. For
#' bivariate data with `plotit=TRUE`, creates a scatterplot marking outliers.
#'
#' @seealso \code{\link{ogk.pairwise}}, \code{\link{out}}, \code{\link{outDETMCD}}
#' @export
outogk<-function(x,sigmamu=taulc,v=gkcov,op=TRUE,SEED=FALSE,
beta=max(c(.95,min(c(.99,1/nrow(x)+.94)))),n.iter=1,plotit=TRUE,...){
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)
if(!op){
temp<-ogk.pairwise(x,sigmamu=sigmamu,v=v,beta=beta,n.iter=n.iter,...)
vals<-hard.rejection(temp$distances,p=ncol(x),beta=beta,...)
flag<-(vals==1)
vals<-c(1:nrow(x))
outid<-vals[!flag]
keep<-vals[flag]
if(is.matrix(x)){
if(ncol(x)==2 && plotit){
plot(x[,1],x[,2],xlab="X", ylab="Y",type="n")
points(x[flag,1],x[flag,2])
if(sum(!flag)>0)points(x[!flag,1],x[!flag,2],pch="o")
}}}
if(op){
temp<-out(x,cov.fun=ogk,beta=beta,plotit=plotit,SEED=SEED)
outid<-temp$out.id
keep<-temp$keep
}
list(out.id=outid,keep=keep,distances=temp$dis)
}


#' Compute Covariance After Removing Outliers
#'
#' Computes covariance (or covariance matrix) after removing outliers detected
#' by a specified outlier detection function.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#'
#' @return A list with component:
#'   \item{cov}{Covariance (bivariate case) or covariance matrix (multivariate).}
#'
#' @details
#' First applies `outfun` to detect and remove outliers, then computes the
#' covariance matrix of the remaining data using `var()`. For bivariate data,
#' returns the single covariance value.
#'
#' @seealso \code{\link{outogk}}, \code{\link{outpro}}
#' @export
outcov<-function(x,y=NA,outfun=outogk,plotit=FALSE){
if(!is.na(y[1]))x<-cbind(x,y)
keep<-outfun(x,plotit=plotit)$keep
val<-var(x[keep,])
if(ncol(val)==2)val<-val[1,2]
list(cov=val)
}


#' TBS-Based Outlier Detection
#'
#' Detects outliers using the TBS (Tukey Bisquare Scale) robust covariance
#' estimator.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... Additional arguments passed to `out()`.
#'
#' @return A list with components:
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{distances}{Robust Mahalanobis distances.}
#'
#' @details
#' Calls `out()` with `cov.fun=tbs` to compute robust Mahalanobis distances
#' based on the TBS covariance estimator, then identifies outliers.
#'
#' @seealso \code{\link{tbs}}, \code{\link{out}}, \code{\link{outDETMCD}}
#' @export
outtbs<-function(x,SEED=FALSE,plotit=TRUE,xlab="X",ylab="Y",...){
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)
temp<-out(x,cov.fun=tbs,plotit=plotit,SEED=SEED,xlab=xlab,ylab=ylab)
outid<-temp$out.id
keep<-temp$keep
list(out.id=outid,keep=keep,distances=temp$dis)
}


#' Robust Covariance-Based Outlier Detection
#'
#' Detects outliers using robust Mahalanobis distance based on a specified
#' robust covariance estimator. Supports multiple robust methods including
#' DETMCD, MVE, MCD, MBA, and TBS estimators.
#'
#' @inheritParams common-params
#' @param cov.fun Covariance function to use. Options:
#'   \itemize{
#'     \item `DETMCD`: Deterministic MCD (default)
#'     \item `cov.mve`: MVE estimate
#'     \item `cov.mcd`: MCD estimate
#'     \item `covmba2`: MBA (median ball algorithm)
#'     \item `rmba`: Olive's adjustment of MBA
#'     \item `cov.roc`: Rocke's TBS estimator
#'   }
#' @param xlab Label for x-axis in plots.
#' @param ylab Label for y-axis in plots.
#' @param qval Quantile for chi-squared critical value (default: 0.975).
#' @param crit Custom critical value. If NULL, uses sqrt(qchisq(qval, p)).
#' @param KS Logical. If TRUE, preserves the random seed (default: TRUE).
#' @inheritParams common-params
#' @param ... Additional arguments passed to `cov.fun`.
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers detected.}
#'   \item{out.val}{Values of detected outliers.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{dis}{Robust Mahalanobis distances for all points.}
#'   \item{crit}{Critical value used.}
#'
#' @details
#' Computes robust center and covariance using `cov.fun`, then flags points
#' with robust Mahalanobis distance exceeding the critical value. For bivariate
#' data with `plotit=TRUE`, creates a scatterplot with outliers marked with "*".
#'
#' @seealso \code{\link{outmve}}, \code{\link{outogk}}, \code{\link{outtbs}}
#' @export
outDETMCD<-function(x,cov.fun=DETMCD,xlab='X',ylab='Y',qval=.975,
crit=NULL,KS=TRUE,plotit=FALSE,...){
if(is.data.frame(x))x=as.matrix(x)
if(is.list(x))stop('Data cannot be stored in list mode')
nrem=nrow(as.matrix(x))
if(!is.matrix(x)){
dis<-(x-median(x,na.rm=TRUE))^2/mad(x,na.rm=TRUE)^2
if(is.null(crit))crit<-sqrt(qchisq(.975,1))
vec<-c(1:length(x))
}
if(is.matrix(x)){
mve<-cov.fun(elimna(x))
dis<-mahalanobis(x,mve$center,mve$cov)
if(is.null(crit))crit<-sqrt(qchisq(.975,ncol(x)))
vec<-c(1:nrow(x))
}
dis[is.na(dis)]=0
dis<-sqrt(dis)
chk<-ifelse(dis>crit,1,0)
id<-vec[chk==1]
keep<-vec[chk==0]
if(is.matrix(x)){
if(ncol(x)==2 && plotit){
plot(x[,1],x[,2],xlab=xlab,ylab=ylab,type='n')
flag<-rep(TRUE,nrow(x))
flag[id]<-FALSE
points(x[flag,1],x[flag,2])
if(sum(!flag)>0)points(x[!flag,1],x[!flag,2],pch='*')
}}
if(!is.matrix(x))outval<-x[id]
if(is.matrix(x))outval<-x[id,]
n=nrow(as.matrix(x))
n.out=length(id)
list(n=n,n.out=n.out,out.val=outval,out.id=id,keep=keep,dis=dis,crit=crit)
}


#' ICS-Based Outlier Detection
#'
#' Detects outliers using the Invariant Coordinate Selection (ICS) method.
#' Requires the ICSOutlier package.
#'
#' @inheritParams common-params
#' @param n.id Number of outliers expected (must be between 0 and n/2). If
#'   NULL, prints ICS results without identifying specific outliers. Run twice:
#'   first with n.id=NULL to see results, then with n.id set to the number
#'   detected.
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{out.id}{Indices of outliers (only if n.id is specified).}
#'   \item{keep}{Indices of non-outliers (only if n.id is specified).}
#'
#' @details
#' The ICS method finds invariant coordinates via generalized eigenvalue
#' decomposition and detects outliers based on extreme distances in these
#' coordinates. See Archimbaud et al. (CSDA) for details.
#'
#' @references
#' Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2018). ICS for
#' multivariate outlier detection. Computational Statistics & Data Analysis, 128, 200-215.
#'
#' @seealso \code{\link[ICSOutlier]{ics.outlier}}, \code{\link{outpro}}
#' @export
outICS<-function(x,n.id=NULL){
library(ICSOutlier)
x=elimna(x)
n=nrow(x)
v=ics2(x)
print(ics.outlier(v))
if(!is.null(n.id)){
if(n.id>n/2)stop('n.id should be less than n/2')
if(n.id<=0)stop('n.id should be greater than zero')
d=ics.distances(v)
dr=rank(d)
idout=n+1-n.id
id=which(dr>=idout)
j=c(1:n)
}
list(n=n,out.id=id,keep=j[-id])
}


#' Modified Carling Outlier Detection (Skewness-Adjusted)
#'
#' Detects univariate outliers using a modification of Carling's (2000) boxplot
#' rule that adjusts for skewness. Uses asymmetric fences based on distances
#' from median to quartiles.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{out.val}{Values of detected outliers.}
#'   \item{out.id}{Indices of outliers.}
#'   \item{keep}{Indices of non-outliers.}
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of outliers.}
#'   \item{cl}{Lower critical value.}
#'   \item{cu}{Upper critical value.}
#'
#' @details
#' Uses Carling's sample-size-adjusted multiplier but with asymmetric fences
#' to account for skewness: `cl = M - gval*2*(M-QL)` and `cu = M + gval*2*(QU-M)`,
#' where M is the median, QL and QU are the ideal fourths, and gval is Carling's
#' formula.
#'
#' @seealso \code{\link{outbox}}, \code{\link{idealf}}
#' @export
outmc<-function(x,plotit=FALSE){
x=elimna(x)
temp<-idealf(x)
gval<-(17.63*n-23.64)/(7.74*n-3.71)
M=median(x)
cl=M-gval*2*(M-temp$ql)
cu=M+gval*2*(temp$qu-M)
n=length(x)
flag<-NA
outid<-NA
vec<-c(1:n)
for(i in 1:n){
flag[i]<-(x[i]< cl || x[i]> cu)
}
if(sum(flag)==0)outid<-NULL
if(sum(flag)>0)outid<-vec[flag]
keep<-vec[!flag]
outval<-x[flag]
n.out=sum(length(outid))
list(out.val=outval,out.id=outid,keep=keep,n=n,n.out=n.out,cl=cl,cu=cu)
}


#' Outlier Detection with Dummy Variable Column Removal
#'
#' For regression with dummy coding, removes specified columns (typically dummy
#' variables) before applying outlier detection.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @param id Column index to remove before outlier detection (e.g., dummy variables).
#' @inheritParams common-params
#' @param ... Additional arguments passed to `outfun`.
#'
#' @return Output from `outfun` (typically a list with `out.id` and `keep` components).
#'
#' @details
#' Useful when checking for outliers in regression with categorical predictors
#' coded as dummy variables. The dummy variable columns should be excluded from
#' multivariate outlier detection.
#'
#' @seealso \code{\link{outpro}}, \code{\link{out.by.groups}}
#' @export
out.dummy<-function(x,outfun=outpro,id,plotit=FALSE,...){
x=as.matrix(x)
if(ncol(x)==1)stop(' Should have two or more columns')
X=x[,-id]
a=outfun(X,plotit=FALSE)
a
}


#' Group-wise Outlier Detection
#'
#' Divides data into groups and identifies outliers within each group separately,
#' then returns indices of outliers and non-outliers across all groups.
#'
#' @inheritParams common-params
#' @param grp.col Column index indicating group membership.
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @param ... Additional arguments passed to `outfun`.
#'
#' @return A list with components:
#'   \item{out.id}{Row indices of outliers across all groups.}
#'   \item{keep}{Row indices of non-outliers across all groups.}
#'
#' @details
#' Useful for identifying outliers when data contain multiple groups. Applies
#' outlier detection within each group separately, which can be more appropriate
#' than detecting outliers in the pooled data when groups have different
#' distributions.
#'
#' @seealso \code{\link{outpro}}, \code{\link{fac2Mlist}}
#' @export
out.by.groups<-function(x,grp.col,outfun=outpro,pr=TRUE,plotit=FALSE,...){
x=elimna(x)
p=ncol(x)
p1=p+1
pv=c(1:p)
pv=pv[-grp.col]
#pv=c(pv,p1)
n=nrow(x)
ones=c(1:n)
w=cbind(x,ones)
z=fac2Mlist(w,grp.col=grp.col,c(1:p1),pr=FALSE)
MAT=NULL
for(j in 1:length(z)){
m=z[[j]]
a=outfun(m[,pv],plotit=FALSE)
MAT=rbind(MAT,m[a$keep,])
}
keep=MAT[,p1]
ou=ones[-keep]
list(out.id=ou,keep=keep)
}



#' Wrapper for Multiple Outlier Detection Methods
#'
#' Provides a unified interface for calling different outlier detection methods
#' by specifying a method name.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @param regfun Regression function for BLP method (default: tsreg).
#' @inheritParams common-params
#' @param id Column index for DUM method (dummy variable to exclude).
#' @param method Character string specifying the method:
#'   \itemize{
#'     \item `"PRO"`: Projection method (`outpro`)
#'     \item `"PRO.R"`: Fast projection method (`outpro.depth`)
#'     \item `"BLP"`: Bad leverage point detection (`outblp`)
#'     \item `"DUM"`: Projection with dummy variable removal (`out.dummy`)
#'     \item `"MCD"`: MCD-based detection (`outDETMCD`)
#'     \item `"BOX"`: Boxplot method (`outbox`)
#'   }
#'
#' @return Output from the selected method (typically a list with `out.id` and
#'   `keep` components).
#'
#' @details
#' This is a convenience function that dispatches to different outlier detection
#' methods based on the `method` argument. Default method is `"PRO"`.
#'
#' @seealso \code{\link{outpro}}, \code{\link{outpro.depth}}, \code{\link{outblp}},
#'   \code{\link{out.dummy}}, \code{\link{outDETMCD}}, \code{\link{outbox}}
#' @export
out.methods<-function(x,y, regfun = tsreg,plotit=FALSE,id,method=c('PRO','PRO.R','BLP','DUM','MCD','BOX')){
type=match.arg(method)
switch(type,
    PRO=outpro(x,plotit=plotit),         # projection method
    PRO.R=outpro.depth(x),   #projection method   random, lower execution time vs outpro
    BLP=outblp(x,y,regfun=regfun,plotit=FALSE),       # regression method
    DUM=out.dummy(x,y,outfun=outpro.depth,id=id),   #   Detect outliers ignoring  col indicated by argument id
    MCD=outDETMCD(x,plotit=plotit),
    BOX=outbox(x))  # Boxplot method using ideal. fourths
}


#' Bad Leverage Point Detection (Hybrid Method)
#'
#' Detects bad leverage points in regression using a blend of homoscedastic
#' and heteroscedastic methods to control Type I error rates.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @param regfun Regression function to use (default: tsreg).
#' @param omit.col Column indices to omit when checking for bad leverage points
#'   (e.g., dummy variables). If NULL, uses all columns.
#' @inheritParams common-params
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#'
#' @return A list with components:
#'   \item{n}{Sample size.}
#'   \item{n.out}{Number of bad leverage points detected.}
#'   \item{bad.lev}{Indices of bad leverage points.}
#'   \item{keep}{Indices of non-outliers.}
#'
#' @details
#' Combines results from `reglev.gen()` (generalized leverage) and `regcon.out()`
#' (regression outlier detection) to identify bad leverage points that affect
#' regression fits. This hybrid approach helps avoid Type I errors in hypothesis
#' testing.
#'
#' @seealso \code{\link{reglev.gen}}, \code{\link{regcon.out}}
#' @export
outblp.HH<-function(x,y,regfun=tsreg,omit.col=NULL,plotit=TRUE,xlab='X',ylab='Y'){
xy=elimna(cbind(x,y))
n=nrow(xy)
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
if(p>1){
if(!is.null(omit.col))
x=x[,-omit.col]
}
x<-as.matrix(x)
out.id=NULL
temp=reglev.gen(x,y,regfun=regfun,plotit=FALSE)
out.id=temp$bad.lev
temp2=regcon.out(x,y,plotit=FALSE)
vec=keep=c(1:n)
out.id=unique(c(out.id,temp2$bad.lev))
if(length(out.id)>0)keep=vec[-out.id]
n.out=length(out.id)
if(plotit){
plot(x,y,type='n',xlab=xlab,ylab=ylab)
points(x[keep],y[keep],pch='*')
points(x[out.id],y[out.id],pch='o')
}
list(n=n,n.out=n.out,bad.lev=out.id,keep=keep)
}




#' Compute Exact Halfspace Depth for Bivariate Data
#'
#' Computes exact halfspace (Tukey) depth for bivariate data. For each point,
#' calculates the minimum proportion of data in any closed halfspace containing
#' that point.
#'
#' @inheritParams common-params
#' @param pts Points for which to compute depth. If NA, uses `x` (default: NA).
#' @inheritParams common-params
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#'
#' @return Vector of depth values for each point in `pts`.
#'
#' @details
#' Only works with bivariate data (2 columns). Uses the exact algorithm for
#' halfspace depth. With `plotit=TRUE`, creates a scatterplot showing the
#' deepest point (center) marked with "+" and the median depth contour.
#'
#' @seealso \code{\link{fdepthv2}}, \code{\link{prodepth}}
#' @export
depth2<-function(x,pts=NA,plotit=TRUE,xlab="VAR 1",ylab="VAR 2"){
if(ncol(x)!=2)stop("x must be a matrix with 2 columns")
x<-elimna(x)
if(is.na(pts[1]))pts<-x
if(ncol(pts)!=2)stop("Argument pts must be stored as a matrix with 2 columns")
pts<-as.matrix(pts)
ndepth<-NA
for(i in 1:nrow(pts)){
ndepth[i]<-depth(pts[i,1],pts[i,2],x)
}
if(plotit){
m<-x
plot(m,xlab=xlab,ylab=ylab)
flag<-(ndepth==max(ndepth))
if(sum(flag)==1)center<-m[flag,]
if(sum(flag)>1)center<-apply(m[flag,],2,mean)
points(center[1],center[2],pch="+")
temp<-ndepth
flag<-(temp>=median(temp))
xx<-x[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
ndepth
}


#' Compare Two Regression Models Using Depth
#'
#' Compares two regression models by examining the depth of residuals.
#'
#' @param x1,y1 Predictor and response for first model.
#' @param x2,y2 Predictor and response for second model.
#' @inheritParams common-params
#' @param fr Span parameter for running interval smoother (default: 1).
#'
#' @return Depth-based comparison statistic.
#' @keywords internal
depthcom<-function(x1,y1,x2,y2,est=tmean,fr=1){
temp1=depthcomsub(x1,y1,x2,y2,est=est,fr=fr)
temp2=depthcomsub(x2,y2,x1,y1,est=est,fr=fr)
dep=max(c(abs(temp1$dep1-temp1$dep2),abs(temp2$dep1-temp2$dep2)))
dep
}

#' Helper Function for depthcom
#' @keywords internal
depthcomsub<-function(x1,y1,x2,y2,est=tmean,fr=1){
x1=(x1-median(x1))/mad(x1)
x2=(x2-median(x2))/mad(x2)
yh1=runhat(x1,y1,est=tmean,fr=fr)
yh2=runhat(x2,y2,pts=x1,est=tmean,fr=fr)
flag=is.na(yh2)
res1=y1-yh1
res2=y1[!flag]-yh2[!flag]
dep1=resdepth(x1,res1)
dep2=resdepth(x1[!flag],res2)
list(dep1=dep1,dep2=dep2)
}


#' Compare Two Independent Groups Using Depth
#'
#' Compares two independent multivariate groups using data depth methods.
#' Sensitive to differences in both location and scale.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @param MD Logical. If TRUE, uses Mahalanobis depth; if FALSE, uses Tukey
#'   depth (default: FALSE). Automatically set to TRUE if p>2.
#' @inheritParams common-params
#' @param op Logical. If TRUE, prints progress messages.
#' @param fast Logical. If TRUE, uses faster algorithm (`lsqs2.for()`).
#' @inheritParams common-params
#' @param xlab,ylab Axis labels for bivariate plots.
#'
#' @return A list with component:
#'   \item{difci}{Bootstrap confidence interval for depth difference.}
#'
#' @details
#' Uses bootstrap to construct a confidence interval for the difference in
#' depth-based locations. For bivariate data with `plotit=TRUE`, creates a
#' scatterplot with median depth contours for each group.
#'
#' @seealso \code{\link{lsqs2}}, \code{\link{depth2}}
#' @export
depthg2<-function(x,y,alpha=.05,nboot=500,MD=FALSE,plotit=TRUE,op=FALSE,fast=FALSE,SEED=TRUE,
xlab="VAR 1",ylab="VAR 2"){
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
x=elimna(x)
y=elimna(y)
x=as.matrix(x)
y=as.matrix(y)
if(is.matrix(x) && is.matrix(y)){  # YES, code is odd.
nv1<-nrow(x)
nv2<-nrow(y)
if(ncol(x)!=ncol(y))stop("Number of columns of x is not equal to number for y")
if(ncol(x) >2)MD<-T
if(ncol(x)==2 && plotit){
plot(rbind(x,y),type="n",xlab=xlab,ylab=ylab)
points(x,pch="*")
points(y,pch="o")
temp<-NA
for(i in 1:nrow(x)){
temp[i]<-depth(x[i,1],x[i,2],x)
}
flag<-(temp>=median(temp))
xx<-x[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
temp<-NA
for(i in 1:nrow(y)){
temp[i]<-depth(y[i,1],y[i,2],y)
}
flag<-(temp>=median(temp))
xx<-y[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
flag<-(temp>=median(temp))
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,],lty=2)
lines(xx[c(temp[1],temp[length(temp)]),],lty=2)
}
print("Taking bootstrap samples. Please wait.")
data1<-matrix(sample(nv1,size=nv1*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(nv2,size=nv2*nboot,replace=TRUE),nrow=nboot)
qhatd<-NA
dhatb<-NA
for(ib in 1:nboot){
if(op)print(paste("Bootstrap sample ",ib," of ",nboot, "is complete."))
if(!fast)temp<-lsqs2(x[data1[ib,],],y[data2[ib,],],plotit=FALSE,MD=MD)
if(fast)temp<-lsqs2.for(x[data1[ib,],],y[data2[ib,],],plotit=FALSE,MD=MD)
qhatd[ib]<-temp[[1]]-temp[[2]]
}
temp<-sort(qhatd)
lv<-round(alpha*nboot/2)
uv<-nboot-lv
difci<-c(temp[lv+1],temp[uv])
}
#
if(!is.matrix(x) && !is.matrix(y)){
nv1<-length(x)
nv2<-length(y)
print("Taking bootstrap samples. Please wait.")
data1<-matrix(sample(nv1,size=nv1*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(nv2,size=nv2*nboot,replace=TRUE),nrow=nboot)
qhatd<-NA
dhatb<-NA
for(ib in 1:nboot){
if(!fast)temp<-lsqs2(x[data1[ib,]],y[data2[ib,]],plotit=FALSE,MD=MD)
if(fast)temp<-lsqs2.for(x[data1[ib,]],y[data2[ib,]],plotit=FALSE,MD=MD)
qhatd[ib]<-temp[[1]]-temp[[2]]
dhatb[ib]<-(temp[[1]]+temp[[2]])/2
}}
temp<-sort(qhatd)
temp2<-sort(dhatb)
lv<-round(alpha*nboot/2)
uv<-nboot-lv
difci<-c(temp[lv+1],temp[uv])
list(difci=difci)
}

hochberg<-
function(x,x2=NA,cil=NA,con=0,tr=.2,alpha=.05){
#
# A generalization of Hochberg's two-stage method
# method to trimmed mean#
#
# THIS FUNCTION WAS UPDATED FEB., 2024. IT NOW HAS A MORE CONVENIENT AND
# SLIGHTLY MORE ACCURATE METHOD FOR
# COMPUTING THE CRITICAL VALUE; NO NEED TO USE  TABLES AS BEFORE.
#
# x contains first stage data
# x2 contains second stage data
#
# cil is the desired length of the confidence intervals.
# That is, cil is the distance between the upper and lower
# ends of the confidence intervals.
#
x3<-x2
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
J<-length(x)
tempn<-0
svec<-NA
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
svec[j]<-winvar(temp,tr=tr)/(1-2*tr)^2
}
tempt<-floor((1-2*tr)*tempn)
A<-sum(1/(tempt-1))
df<-J/A
if(!is.list(x2) && !is.matrix(x2)){
x2<-list()
for(j in 1:J)x2[[j]]<-NA
}
if(is.na(cil))stop("To proceed, you must specify the maximum length of the confidence intervals.")
#crit<-trange(tempn-1,alpha=alpha,iter=iter,SEED=SEED) #OLD CODE
crit=qtukey(1-alpha,J,df)
#
if(con[1] == 0){
               Jm<-J-1
               ncon <- (J^2 - J)/2
                con <- matrix(0, J, ncon)
                id <- 0
                for(j in 1:Jm) {
                        jp <- j + 1
                        for(k in jp:J) {
                                id <- id + 1
                                con[j, id] <- 1
                                con[k, id] <- 0 - 1
                        }
                }
        }
        ncon <- ncol(con)
avec<-NA
for(i in 1:ncon){
temp<-con[,i]
avec[i]<-sum(temp[temp>0])
}
dvec<-(cil/(2*crit*avec))^2
d<-max(dvec)
n.vec<-NA
for(j in 1:J){
n.vec[j]<-max(tempn[j],floor(svec[j]/d)+1)
print(paste("Need an additional ", n.vec[j]-tempn[j],
" observations for group", j))
}
#
# Do second stage if data are supplied
#
ci.mat=NULL
if(!is.na(x2[1])){
if(is.matrix(x2))x2<-listm(x2)
temp2<-n.vec-tempn
#if(!is.list(x3) && !is.matrix(x3) && sum(temp2)>0)stop("No second stage data supplied; this function is terminating")
if(length(x) != length(x2))warning("Number of groups in first stage data does not match the number in the second stage.")
ci.mat<-NA
if(!is.na(x2[1]) || sum(temp2)==0){
xtil<-NA
nvec2<-NA
for(j in 1:J){
nvec2[j]<-0
temp<-x2[[j]]
if(!is.na(temp[1]))nvec2[j]<-length(x2[[j]])
if(nvec2[j] <n.vec[j]-tempn[j])warning(paste("The required number of observations for group",j," in the second stage is ",n.vec[j]-tempn[j]," but only ",nvec2[j]," are available"))
xtil[j]<-mean(c(x[[j]],x2[[j]]),tr=tr,na.rm=TRUE)
}
ci.mat<-matrix(0,ncol=3,nrow=ncon)
dimnames(ci.mat)<-list(NULL,c("con.num","ci.low","ci.high"))
for(ic in 1:ncon){
ci.mat[ic,1]<-ic
bvec<-con[,ic]*sqrt(svec/(nvec2+tempn))
A<-sum(bvec[bvec>0])
C<-0-sum(bvec[bvec<0])
D<-max(A,C)
ci.mat[ic,2]<-sum(con[,ic]*xtil)-crit*D
ci.mat[ic,3]<-sum(con[,ic]*xtil)+crit*D
}}}
list(ci.mat=ci.mat,con=con)
}


#' Helper Function for Computing Depth Combinations
#' @keywords internal
depths1<-function(m,j){
if(m < j)depths1<-0
else{
if(j==1)depths1<-m
if(j==2)depths1<-(m*(m-1))/2
if(j==3)depths1<-(m*(m-1)*(m-2))/6
}
depths1
}


#' Compute Halfspace Depth (Version 2)
#'
#' Computes an approximation of halfspace depth by projecting onto lines
#' connecting pairs of data points. More accurate than `fdepth()` and handles
#' singular covariance matrices, but slower.
#'
#' @param m Data matrix (n x p).
#' @param pts Points for which to compute depth (default: NA, uses `m`).
#' @inheritParams common-params
#'
#' @return Vector of depth values for each point in `pts` (or `m` if pts=NA).
#'
#' @details
#' Projects data onto lines connecting each pair of distinct points and computes
#' univariate depth of projections. The final depth is the minimum across all
#' projections. Requires O(n^2) space but handles data with singular covariance.
#'
#' For bivariate data with `plotit=TRUE`, creates a scatterplot with the median
#' depth contour. The deepest point (center) is marked with "+".
#'
#' @seealso \code{\link{fdepth}}, \code{\link{depth2}}, \code{\link{prodepth}}
#' @export
fdepthv2<-function(m,pts=NA,plotit=TRUE){
m<-elimna(m) # Remove missing values
if(!is.na(pts[1]))remm<-m
if(!is.matrix(m))dep<-unidepth(m)
if(is.matrix(m)){
nm<-nrow(m)
nt<-nm
nm1<-nm+1
if(!is.na(pts[1])){
if(ncol(m)!=ncol(pts))stop("Number of columns of m is not equal to number of columns for pts")
nt<-nm+nrow(pts)
}}
if(ncol(m)==1)depth<-unidepth(m)
if(ncol(m)>1){
m<-elimna(m) # Remove missing values
nc<-(nrow(m)^2-nrow(m))/2
if(is.na(pts[1]))mdep <- matrix(0,nrow=nc,ncol=nrow(m))
if(!is.na(pts[1])){
mdep <- matrix(0,nrow=nc,ncol=nrow(pts))
}
ic<-0
for (iall in 1:nm){
for (i in 1:nm){
if(iall < i){
ic<-ic+1
B<-m[i,]-m[iall,]
dis<-NA
BB<-B^2
bot<-sum(BB)
if(bot!=0){
if(is.na(pts[1])){
for (j in 1:nrow(m)){
A<-m[j,]-m[iall,]
temp<-sum(A*B)*B/bot
dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
}}
if(!is.na(pts[1])){
m<-rbind(remm,pts)
for (j in 1:nrow(m)){
A<-m[j,]-m[iall,]
temp<-sum(A*B)*B/bot
dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
}}
#
# For ic_th projection, store depths of
# points in mdep[ic,]
#
if(is.na(pts[1]))mdep[ic,]<-unidepth(dis)
if(!is.na(pts[1])){
mdep[ic,]<-unidepth(dis[1:nm],dis[nm1:nrow(m)])
}}
if(bot==0)mdep[ic,]<-rep(0,ncol(mdep))
}}}
dep<-apply(mdep,2,min)
}
if(ncol(m)==2 &&is.na(pts[1])){
flag<-chull(m)
dep[flag]<-min(dep)
}
if(ncol(m)==2){
if(is.na(pts[1]) && plotit){
plot(m)
x<-m
temp<-dep
flag<-(temp>=median(temp))
xx<-x[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}}
dep
}


#' Compute Inward Depth Based on Generalized Variance
#'
#' Computes inward depth for each point based on how much the generalized
#' variance decreases when the point is removed.
#'
#' @param m Data matrix (n x p).
#'
#' @return Vector of inward depth values.
#'
#' @details
#' For each point, removes it from the dataset and computes the generalized
#' variance of the remaining points. Higher values indicate points whose removal
#' increases variance (i.e., central points).
#'
#' @seealso \code{\link{gvar}}, \code{\link{prodepth}}
#' @export
indepth<-function(m){
m<-as.matrix(m)
dep<-NA
n<-nrow(m)
flag<-rep(TRUE,n)
for(i in 1:n){
flag[i]<-FALSE
dep[i]<-gvar(m[flag,])
flag[i]<-TRUE
}
dep
}


#' Compute Projection Depth
#'
#' Computes projection depth as 1/(1+d) where d is the projection distance
#' from `pdis()`.
#'
#' @param m Data matrix (n x p).
#' @param pts Points for which to compute depth (default: `m`).
#' @param MM Logical. If TRUE, uses MAD; if FALSE, uses IQR.
#' @param cop Integer specifying center computation method (see `outproMC`).
#' @param dop Depth option parameter.
#' @param center Optional center point.
#' @inheritParams common-params
#'
#' @return Vector of projection depth values.
#'
#' @seealso \code{\link{pdis}}, \code{\link{prodepth}}, \code{\link{zoudepth}}
#' @export
pdepth<-function(m,pts=m,MM=FALSE,cop=3,dop=1,center=NA, SEED=TRUE){
v=pdis(m,pts=pts,MM=MM,cop=cop,dop=dop,center=center)
v=1/(1+v)
v
}

#' Compute Projection Depth (Fast Approximation)
#'
#' Computes an approximation of projection depth using random projections via
#' the DepthProc package. Much faster than `zoudepth()` but gives approximate
#' results.
#'
#' @inheritParams common-params
#' @param pts Points for which to compute depth (default: `x`).
#' @param ndir Number of random projections to use (default: 1000). Higher
#'   values give more accurate approximations but take longer.
#' @inheritParams common-params
#'
#' @return Vector of projection depth values.
#'
#' @details
#' Uses `depthProjection()` from the DepthProc package. Running the function
#' twice on the same data may give slightly different results unless `SEED=TRUE`
#' (default). This is the recommended depth function for most applications due
#' to its speed.
#'
#' Requires the DepthProc package.
#'
#' @seealso \code{\link{zoudepth}}, \code{\link{pdepth}}, \code{\link[DepthProc]{depthProjection}}
#' @export
prodepth<-function(x,pts=x,ndir=1000,SEED=TRUE){
if(SEED){
oldSeed <- .Random.seed
set.seed(45)
}
library(DepthProc)
res=as.vector(depthProjection(pts,x,ndir=ndir))
if(SEED) {
    assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
}
res
}

#' Compute Regression Depth of a Line
#'
#' Computes the regression depth of a line (specified by intercept and slope)
#' relative to a bivariate dataset. Based on Rousseeuw and Hubert (1996).
#'
#' @param d Vector with two components: d[1] = intercept, d[2] = slope.
#' @inheritParams common-params
#' @inheritParams common-params
#' @param sortx Logical. If FALSE, assumes (x,y) is already sorted by x-coordinates.
#'   Set to TRUE (default) to sort internally.
#'
#' @return Regression depth value (scalar).
#'
#' @details
#' Regression depth measures how deep a regression line is within the bivariate
#' data cloud. It counts the minimum number of points that must be removed to
#' make the line non-fit. Higher values indicate more central lines.
#'
#' @references
#' Rousseeuw, P.J. and Hubert, M. (1996). Regression Depth. Journal of the
#' American Statistical Association, 91, 97-104.
#'
#' @seealso \code{\link{resdepth}}, \code{\link{mregdepth}}
#' @export
rdepth.orig<-function(d, x, y, sortx = TRUE)
{
        if(!is.vector(x) || !is.vector(y)) stop("x and y should be vectors")
        n <- length(x)
        if(n < 2)
                stop("you need at least two observations")
        xy <- cbind(x, y)
        b <- d[1]
        a <- d[2]
        if(sortx)
                xy <- xy[order(xy[, 1], xy[, 2]),  ]
        res <- xy[, 2] - a * xy[, 1] - b
        res[abs(res) < 9.9999999999999995e-08] <- 0
        posres <- res >= 0
        negres <- res <= 0
        lplus <- cumsum(posres)
        rplus <- lplus[n] - lplus
        lmin <- cumsum(negres)
        rmin <- lmin[n] - lmin
        depth <- pmin(lplus + rmin, rplus + lmin)
        min(depth)
}


#' Compute Regression Depth from Residuals
#'
#' Computes regression depth of a fit based on its residuals. Works with any
#' regression or smoothing method (parametric, nonparametric, etc.).
#'
#' @inheritParams common-params
#' @param res Vector of residuals from a regression fit.
#'
#' @return Regression depth value (scalar).
#'
#' @details
#' Based on a modification of Rousseeuw and Hubert (1996). Given residuals from
#' any regression method, computes the depth of the implied fit. Higher values
#' indicate fits that are more central to the data cloud.
#'
#' Missing values in residuals are removed before computation.
#'
#' @references
#' Rousseeuw, P.J. and Hubert, M. (1996). Regression Depth. Technical Report,
#' University of Antwerp.
#'
#' @seealso \code{\link{rdepth.orig}}, \code{\link{mregdepth}}
#' @export
resdepth<-function(x,res)
{
        if(!is.vector(x)) stop("x should be a vector")
        n <- length(x)
        if(n < 2)
                stop("you need at least two observations")
flag=is.na(res)
x=x[!flag]
res[!flag]
xord=order(x)
x=x[xord]
res=res[xord]
        posres <- res >= 0
        negres <- res <= 0
        lplus <- cumsum(posres)
        rplus <- lplus[n] - lplus
        lmin <- cumsum(negres)
        rmin <- lmin[n] - lmin
        depth <- pmin(lplus + rmin, rplus + lmin)
        min(depth)
}

#' Helper Function for Regression Depth from Residuals
#' @keywords internal
resdepth.sub<-function(x,res)
{
        if(!is.vector(x)) stop("x should be vectors")
        n <- length(x)
        if(n < 2)
                stop("you need at least two observations")
flag=is.na(res)
x=x[!flag]
res[!flag]
xord=order(x)
x=x[xord]
res=res[xord]
        posres <- res >= 0
        negres <- res <= 0
        lplus <- cumsum(posres)
        rplus <- lplus[n] - lplus
        lmin <- cumsum(negres)
        rmin <- lmin[n] - lmin
        depth <- pmin(lplus + rmin, rplus + lmin)
        min(depth)
}


#' Compute Zuo's Projection Depth
#'
#' Computes projection depth using the method of Zuo (2003). Uses Nelder-Mead
#' optimization to find the projection with maximum outlyingness.
#'
#' @param m Data matrix (n x p).
#' @param pts Points for which to compute depth (default: `m`).
#' @param zloc Location estimator for projections (default: median).
#' @param zscale Scale estimator for projections (default: mad).
#'
#' @return Vector of (inverse) depth values. Lower values indicate greater depth.
#'
#' @details
#' For each point, finds the projection direction that maximizes the standardized
#' distance from the location. The depth is the inverse of this maximum distance.
#' More computationally intensive than `prodepth()` but provides exact values.
#'
#' @references
#' Zuo, Y. (2003). Projection-based depth functions and associated medians.
#' The Annals of Statistics, 31, 1460-1490.
#'
#' @seealso \code{\link{zoudepth}}, \code{\link{prodepth}}, \code{\link{pdepth}}
#' @export
zdepth<-function(m,pts=m,zloc=median,zscale=mad){
if(!is.matrix(m))stop("argument m should be a matrix")
if(!is.matrix(pts))stop("argument pts should be a matrix")
if(ncol(m)!=ncol(pts))stop("Number of columns for m and pts are not equal")
np<-ncol(m)
val<-NA
for(i in 1:nrow(pts)){
pval<-pts[i,]
START<-rep(1,np)/sqrt(np)
temp<-nelderv2(m,np,FN=zdepth.sub,START=START,zloc=zloc,zscale=zscale,pts=pval)
temp<-temp/sqrt(sum(temp^2))
y<-t(t(m)*temp)
y<-apply(y,1,sum)
ppro<-sum(pval*temp)
val[i]<-abs(ppro-zloc(y))/zscale(y)
}
val
}


#' Helper Function for zdepth - Compute Projection Outlyingness
#' @keywords internal
zdepth.sub<-function(x,theta,zloc=median,zscale=mad,pts=NA){
theta<-theta/sqrt(sum(theta^2))
temp<-t(t(x)*theta)
ppro<-sum(t(t(pts)*theta))
yhat<-apply(temp,1,sum)
val<-0-abs(ppro-zloc(yhat))/zscale(yhat)
val
}

zdist=zdepth


#' Compute Projection Depth Using Zuo's Method
#'
#' Wrapper for `zdepth()` that returns depth values (rather than inverse depth).
#' Uses Nelder-Mead optimization to find projection depth.
#'
#' @inheritParams zdepth
#' @inheritParams common-params
#'
#' @return Vector of projection depth values (higher = more central).
#'
#' @details
#' Computes 1/(1 + zdepth()), converting Zuo's outlyingness measure to a depth
#' measure. More computationally intensive than `prodepth()` but exact.
#'
#' @seealso \code{\link{zdepth}}, \code{\link{prodepth}}
#' @export
zoudepth<-function(x,pts=x, zloc = median, zscale = mad, SEED=TRUE){
res=1/(1+zdepth(x,pts,zloc,zscale))
res
}

#' Compute Univariate Depth
#'
#' Computes halfspace depth for univariate data. The depth of a point is the
#' minimum of the proportions of data on either side of it.
#'
#' @inheritParams common-params
#' @param pts Points for which to compute depth (default: NA, uses `x`).
#'
#' @return Vector of depth values (between 0 and 0.5).
#'
#' @details
#' For univariate data, halfspace depth reduces to min(F(x), 1-F(x)), where F is
#' the empirical CDF. The deepest point (median) has depth 0.5.
#'
#' @seealso \code{\link{depth2}}, \code{\link{prodepth}}
#' @export
unidepth<-function(x,pts=NA){
if(!is.vector(x))stop("x should be a vector")
if(is.na(pts[1]))pts<-x
pup<-apply(outer(pts,x,FUN="<="),1,sum)/length(x)
pdown<-apply(outer(pts,x,FUN="<"),1,sum)/length(x)
pdown<-1-pdown
m<-matrix(c(pup,pdown),nrow=2,byrow=TRUE)
dep<-apply(m,2,min)
dep
}


#' Depth-Based Classification (Discriminant Analysis)
#'
#' Classifies test data into two groups based on data depth. Assigns each test
#' point to the group where it has greater depth.
#'
#' @param train Training data matrix (optional if `x1` and `x2` provided).
#' @param test Test data matrix to classify.
#' @inheritParams common-params
#' @param x1,x2 Training data for groups 1 and 2 (alternative to `train`/`g`).
#' @param depthfun Depth function to use (default: prodepth).
#' @param ... Additional arguments passed to `depthfun`.
#'
#' @return Vector of predicted class labels (1 or 2) for each test observation.
#'
#' @details
#' For each test point, computes its depth relative to each training group.
#' Assigns the point to the group where it has higher depth. Can specify data
#' via `train`/`g` or directly via `x1`/`x2`.
#'
#' @seealso \code{\link{Depth.class.bag}}, \code{\link{prodepth}}
#' @export
discdepth<-function(train=NULL,test=NULL,g,x1=NULL,x2=NULL,depthfun=prodepth,...){
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}
Train=cbind(train,g)
Train=elimna(Train)
p=ncol(train)
p1=p+1
train=Train[,1:p]
g=Train[,p1]
flag=g==min(g)
x1=Train[flag,1:p]
x2=Train[!flag,1:p]
}
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
if(is.null(test))stop('No test data, argument test is NULL')
test=elimna(test)
z=as.matrix(test)
x1=as.matrix(x1)
x2=as.matrix(x2)
z=as.matrix(z)
d1=depthfun(x1,pts=z,...)
d2=depthfun(x2,pts=z,...)
flag=d1>d2
N=nrow(z)
id=rep(2,N)
id[flag]=1
id
}


#' Compute Multiple Regression Depth
#'
#' Approximates regression depth for multiple regression by taking the minimum
#' depth across marginal regressions.
#'
#' @param X Predictor matrix (n x p).
#' @param RES Vector of residuals from a regression fit.
#'
#' @return Approximate regression depth (scalar).
#'
#' @details
#' For each predictor, computes the regression depth of residuals against that
#' predictor. Returns the minimum across all predictors as an approximation of
#' the full multiple regression depth.
#'
#' @seealso \code{\link{resdepth}}
#' @export
mregdepth<-function(X,RES){
X=as.matrix(X)
XRES=elimna(cbind(X,RES))
p=ncol(X)
p1=p+1
vals=NA
for(j in 1:p)vals[j]=resdepth(XRES[,j],XRES[,p1])
mdepthappr=min(vals)
mdepthappr
}



#' Skipped Mean Based on Projection Depth
#'
#' Computes the mean after removing outliers detected using projection depth.
#'
#' @param m Data matrix (n x p).
#'
#' @return Vector of means (length p) after outlier removal.
#'
#' @details
#' Uses `outpro.depth()` to identify outliers via projection depth, then computes
#' the column means of the non-outlier data. A robust location estimator.
#'
#' @seealso \code{\link{outpro.depth}}, \code{\link{tmean}}
#' @export
smean.depth<-function(m){
m=elimna(m)
id=outpro.depth(m)$keep
val=apply(m[id,],2,mean)
val
}


#' Bootstrap Test for Linear Contrasts Using Depth
#'
#' Tests whether specified linear contrasts among J independent groups equal
#' zero, using bootstrap and data depth methods.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @param con Contrast matrix (J x C) where J = number of groups and C = number
#'   of contrasts. If 0, tests all pairwise comparisons (default).
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @param allp Logical. If TRUE and con=0, tests all pairwise contrasts;
#'   if FALSE, tests sequential contrasts (default: TRUE).
#' @param MM Logical. If TRUE, uses MAD in projection distance computation.
#' @inheritParams common-params
#' @param cop Center computation option (see `outproMC`).
#' @inheritParams common-params
#' @inheritParams common-params
#' @param ... Additional arguments passed to `est`.
#'
#' @return A list with components:
#'   \item{p.value}{Bootstrap p-value for testing all contrasts equal zero.}
#'   \item{psihat}{Estimated contrast values.}
#'   \item{con}{Contrast matrix used.}
#'   \item{n}{Sample sizes for each group.}
#'
#' @details
#' Uses bootstrap to test joint hypothesis that all contrasts equal zero.
#' Depth is measured via Mahalanobis distance (op=1), MCD-based Mahalanobis
#' (op=2), or projection distance (op=3, default and recommended).
#'
#' Data can be in list mode or matrix mode (columns = groups).
#'
#' @seealso \code{\link{pdis}}, \code{\link{pdisMC}}
#' @export
pbadepth<-function(x,est=onestep,con=0,alpha=.05,nboot=2000,grp=NA,op=3,allp=TRUE,
MM=FALSE,MC=FALSE,cop=3,SEED=TRUE,na.rm=FALSE,...){
con<-as.matrix(con)
if(is.matrix(x) || is.data.frame(x))x=listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
if(!is.na(grp)){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
x<-xx
}
J<-length(x)
mvec<-NA
nvec=NA
for(j in 1:J){
temp<-x[[j]]
if(na.rm)temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
mvec[j]<-est(temp,...)
nvec[j]=length(temp)
}
Jm<-J-1
d<-ifelse(con==0,(J^2-J)/2,ncol(con))
if(sum(con^2)==0){
if(allp){
con<-matrix(0,J,d)
id<-0
for (j in 1:Jm){
jp<-j+1
for (k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
if(!allp){
con<-matrix(0,J,Jm)
for (j in 1:Jm){
jp<-j+1
con[j,j]<-1
con[jp,j]<-0-1
}}}
bvec<-matrix(NA,nrow=J,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
for(j in 1:J){
#print(paste("Working on group ",j))
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,est,na.rm=na.rm,...) # J by nboot matrix, jth row contains
#                          bootstrapped  estimates for jth group
}
chkna=sum(is.na(bvec))
if(chkna>0){
print("Bootstrap estimates of location could not be computed")
print("This can occur when using an M-estimator")
print("Might try est=tmean")
}
bcon<-t(con)%*%bvec #C by nboot matrix
tvec<-t(con)%*%mvec
tvec<-tvec[,1]
tempcen<-apply(bcon,1,mean)
vecz<-rep(0,ncol(con))
bcon<-t(bcon)
smat<-var(bcon-tempcen+tvec)
temp<-bcon-tempcen+tvec
bcon<-rbind(bcon,vecz)
if(op==1)dv<-mahalanobis(bcon,tvec,smat)
if(op==2){
smat<-cov.mcd(temp)$cov
dv<-mahalanobis(bcon,tvec,smat)
}
if(op==3){
#print("Computing p-value. Might take a while with op=3")
if(!MC)dv<-pdis(bcon,MM=MM,cop=cop)
if(MC)dv<-pdisMC(bcon,MM=MM,cop=cop)
}
bplus<-nboot+1
sig.level<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
list(p.value=sig.level,psihat=tvec,con=con,n=nvec)
}


#' Compare Two Quantile Regression Models Using Depth
#' @keywords internal
Qdepthcom<-function(x1,y1,x2,y2,qval){
temp1=Qdepthcomsub(x1,y1,x2,y2,qval)
temp2=Qdepthcomsub(x2,y2,x1,y1,qval)
dep=max(c(abs(temp1$dep1-temp1$dep2),abs(temp2$dep1-temp2$dep2)))
dep
}

#' Helper Function for Qdepthcom
#' @keywords internal
Qdepthcomsub<-function(x1,y1,x2,y2,qval){
x1=(x1-median(x1))/mad(x1)
x2=(x2-median(x2))/mad(x2)
yh1=qsmcobs(x1,y1,FIT=FALSE,qval=qval,plotit=FALSE)$yhat
temp2=cobs(x2,y2,print.mesg=FALSE,print.warn=FALSE,tau=qval)
yh2=predict(temp2,z=x1)
yh2=yh2[,2]
flag=is.na(yh2)
res1=y1-yh1
res2=y1[!flag]-yh2[!flag]
dep1=resdepth(x1,res1)
dep2=resdepth(x1[!flag],res2)
list(dep1=dep1,dep2=dep2)
}



#' Compare Two Distributions Using SVM and Depth
#'
#' Compares two independent multivariate distributions using Support Vector
#' Machines with leave-one-out cross-validation.
#'
#' @param x1,x2 Data matrices for the two groups.
#' @inheritParams common-params
#' @param depthfun Depth function (default: prodepth). Not currently used in SVM.
#' @inheritParams common-params
#' @param kernel SVM kernel type (default: 'radial'). See `svm()` documentation.
#' @param MISS Logical. If TRUE, returns misclassified vectors.
#' @param TABLE Logical. If TRUE, returns classification confusion table.
#' @param ... Additional arguments (not currently used).
#'
#' @return A list with components:
#'   \item{est.prob}{Estimated probability of correct classification.}
#'   \item{miss.class.vectors}{Matrix of misclassified vectors (if MISS=TRUE).}
#'   \item{TABLE}{Confusion table (if TABLE=TRUE).}
#'
#' @details
#' Uses leave-one-out cross-validation with SVM to estimate classification
#' accuracy. Requires the e1071 package.
#'
#' @references
#' Shao, J. (1993). Linear Model Selection by Cross-Validation. JASA, 88, 486-494.
#'
#' @seealso \code{\link[e1071]{svm}}
#' @export
comdepthsvm<-function(x1,x2,alpha=.05,depthfun=prodepth,
plotit=FALSE,kernel='radial',MISS=FALSE,TABLE=FALSE,...){
library(e1071)
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(1,n1),rep(2,n2))
g1=rep(1,n1)
g2=rep(2,n2)
PCD=NULL
xall=rbind(x1,x2)
n=nrow(xall)
nm1=n-1
rem=NA
for(i in 1:nm1){
svm_model=svm(xall[-i,],as.factor(g[-i]),kernel=kernel)
temp=predict(svm_model,t(as.matrix(xall[i,])))
rem[i]=temp
pick=as.numeric(as.vector(temp[1]))
PCD[i]=pick[1]==g[i]
}
svm_model=svm(xall[-n,],as.factor(g[-n]),kernel=kernel)
temp=predict(svm_model,t(as.matrix(xall[n,])))
rem[n]=temp
pick=as.numeric(as.vector(temp[1]))
pick=as.numeric(as.vector(pick))
PCD[n]=pick[1]==g[n]
MI=NULL
if(MISS){
MI=cbind(xall[!PCD,],g[!PCD])
ir=c(1:n)
idrow=ir[!PCD]
MI=cbind(idrow,MI)
}
tab=NULL
if(TABLE){
tab=table(rem,as.factor(g))
dimnames(tab)=list(c('Pred 1','Pred 2'),c('GRP 1','GRP 2'))
}
list(est.prob=mean(PCD),miss.class.vectors=MI,TABLE=tab)
}



#' Two-by-K ANOVA Using Depth-Based Methods
#'
#' Tests for main effect of Factor A in a 2 x K design using depth of difference
#' vectors. Accounts for the pattern of measures across Factor B levels.
#'
#' @param x1,x2 Data for the two levels of Factor A. Can be matrices (columns =
#'   Factor B levels) or lists.
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @param nmin Minimum sample size (default: 12).
#' @param CR Logical. If TRUE, creates a critical region plot.
#' @param xlab,ylab,zlab Axis labels for plotting (when K=2 or K=3).
#' @inheritParams common-params
#' @param ... Additional arguments passed to `est`.
#'
#' @return A list with components:
#'   \item{p.value}{Bootstrap p-value for testing no main effect of Factor A.}
#'   \item{est1,est2}{Estimated locations for each group.}
#'   \item{dif}{Estimated differences.}
#'   \item{n1,n2}{Sample sizes.}
#'
#' @details
#' For each bootstrap sample, computes the depth of the zero vector relative
#' to the distribution of difference vectors. This tests whether Factor A has
#' an effect while accounting for the multivariate pattern across Factor B levels.
#'
#' @export
aov2depth<-function(x1,x2,est=tmean,nboot=500,SEED=TRUE,nmin=12,CR=FALSE,
xlab=' DIF 1',ylab='DIF 2',zlab='DIF 3',alpha=.05,...){
if(is.matrix(x1)||is.data.frame(x1))x1=listm(x1)
if(is.matrix(x2)||is.data.frame(x2))x2=listm(x2)
J=length(x1)
if(J!=length(x2))stop('x1 and x2 should have same number of groups')
if(SEED)set.seed(2)
for(j in 1:J){
x1[[j]]=na.omit(x1[[j]])
x2[[j]]=na.omit(x2[[j]])
}
n1=mapply(x1,FUN=length)
n2=mapply(x2,FUN=length)
bplus=nboot+1
bvec1=matrix(NA,nrow=nboot,ncol=J)
bvec2=matrix(NA,nrow=nboot,ncol=J)
for(j in 1:J){
data1=matrix(sample(x1[[j]],size=n1[j]*nboot,replace=TRUE),nrow=nboot)
data2=matrix(sample(x2[[j]],size=n2[j]*nboot,replace=TRUE),nrow=nboot)
bvec1[,j]=apply(data1,1,est,...)
bvec2[,j]=apply(data2,1,est,...)
}
difb=bvec1-bvec2
est1=mapply(x1,FUN=est,...)
est2=mapply(x2,FUN=est,...)
dif=est1-est2
m1=var(difb)
nullvec=rep(0,J)
difz=rbind(difb,nullvec)
dis=mahalanobis(difz,dif,m1)
sig=sum(dis[bplus]<=dis)/bplus
if(CR){
dis2<-order(dis[1:nboot])
dis<-sort(dis)
critn<-floor((1-alpha)*nboot)
if(J==2){
plot(difb[,1],difb[,2],xlab=xlab,ylab=ylab)
points(0,0,pch=0)
xx<-difb[dis2[1:critn],]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
if(J==3){
scatterplot3d(difb[dis2[1:critn],],xlab=xlab,ylab=ylab,zlab=zlab,tick.marks=TRUE)
}

}
list(p.value=sig,est1=est1,est2=est2,dif=dif,n1=n1,n2=n2)
}


#' Compute Bagplot-Based Depth
#'
#' Computes depth based on bagdistance from the bagplot method. Requires the
#' mrfDepth package.
#'
#' @inheritParams common-params
#' @param pts Points for which to compute depth (default: NULL, uses `x`).
#' @inheritParams common-params
#'
#' @return Vector of depth values (higher = more central).
#'
#' @details
#' Uses the bagdistance from Rousseeuw's bagplot method. Depth is computed as
#' 1/(bagdistance + 1). The bagplot is a bivariate generalization of the boxplot
#' based on halfspace depth.
#'
#' Requires the mrfDepth package.
#'
#' @seealso \code{\link{outbag}}, \code{\link[mrfDepth]{bagdistance}}
#' @export
bagdepth<-function(x,pts=NULL,SEED=TRUE){
library(mrfDepth)
d=bagdistance(x,pts)$bagdistance
d=1/(d+1)
d
}


#' Estimate Distribution Overlap Using Depth
#'
#' Estimates the extent of overlap between two independent distributions using
#' data depth. Provides a nonparametric measure of effect size.
#'
#' @inheritParams common-params
#' @inheritParams common-params
#' @param fun Depth function to use (default: prodepth).
#' @inheritParams common-params
#' @param xlab,ylab Axis labels for bivariate plots.
#'
#' @return A list with components:
#'   \item{e}{Overall effect size estimate (0 = complete separation, 0.5 = identical distributions).}
#'   \item{e1,e2}{Group-specific effect size estimates.}
#'
#' @details
#' For identical distributions, effect size = 0.5. As distributions separate,
#' effect size approaches 0. Complete separation yields effect size = 0.
#'
#' For bivariate data with `plotit=TRUE`, creates a scatterplot with both groups.
#'
#' @seealso \code{\link{bwdepthMC.ci}}, \code{\link{bwdepth.perm}}, \code{\link{prodepth}}
#' @export
bwdepth<-function(x,y,fun=prodepth,plotit=FALSE,xlab='V1',ylab='V2'){
x=elimna(x)
y=elimna(y)
x=as.matrix(x)
y=as.matrix(y)
n1=nrow(x)
n2=nrow(y)
if(ncol(x)==1){
fun=unidepth
x=as.vector(x)
y=as.vector(y)
}
pdyy=fun(y,y)
pdyx=fun(y,x)
pdxx=fun(x,x)
pdxy=fun(x,y)
v1=NA
v2=NA
ic=0
for(i in 1:n2){
for(j in 1:n1){
ic=ic+1
v1[ic]=pdyy[i]<=pdyx[j]
}}
ic=0
for(j in 1:n1){
for(i in 1:n2){
ic=ic+1
v2[ic]=pdxx[j]<=pdxy[i]
}}
e1=mean(v1)
e2=mean(v2)
e=(n1*e1+n2*e2)/(n1+n2)
x=as.matrix(x)
y=as.matrix(y)
if(plotit){
if(ncol(x)==2){
plot(rbind(x,y),xlab=xlab,ylab=ylab,type='n')
points(x,pch='*')
points(y,pch='o')
}}
list(e=e,e1=e1,e2=e2)
}


#' Bootstrap Confidence Interval for bwdepth Effect Size
#'
#' Computes bootstrap confidence interval for the depth-based effect size from
#' `bwdepth()`.
#'
#' @inheritParams bwdepth
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{n1,n2}{Sample sizes.}
#'   \item{Est}{Point estimate of effect size.}
#'   \item{ci}{Bootstrap confidence interval.}
#'   \item{p.value}{P-value for testing e = 0.5 (identical distributions).}
#'
#' @seealso \code{\link{bwdepth}}, \code{\link{bwdepth.perm}}
#' @export
bwdepthMC.ci<-function(x,y,fun=prodepth,nboot=100,alpha=.05,MC=TRUE,
SEED=TRUE,plotit=FALSE,xlab='V1',ylab='V2'){
if(SEED)set.seed(2)
 crit=qnorm(1-alpha/2)
if(identical(fun,prodepth))MC=FALSE # get odd error otherwise
x=elimna(x)
y=elimna(y)
x=as.matrix(x)
y=as.matrix(y)
n1=nrow(x)
n2=nrow(y)
est=bwdepth(x,y,plotit=plotit,xlab=xlab,ylab=ylab)
id=list()
for(i in 1:nboot)id[[i]]=c(sample(n1,replace=TRUE),sample(n2,replace=TRUE))
if(!MC)BE=lapply(id,bwdepth.sub,x,y,n1,n2,fun=fun)
if(MC)BE=mclapply(id,bwdepth.sub,x,y,n1,n2,fun=fun)
E=matl(BE)
se=sd(E)
c1=est$e-crit*se
c1[2]=est$e+crit*se
test=(est$e-.5)/se
pv=2*(1-pnorm(abs(test)))
list(n1=n1,n2=n2,Est=est$e,ci=c1,p.value=pv)
}


#' Permutation Test for Two Distributions Using Depth
#'
#' Performs permutation test of H0: F = G (identical distributions) using
#' depth-based effect size.
#'
#' @inheritParams bwdepth
#' @param reps Number of permutation resamples (default: 500).
#' @inheritParams common-params
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{Est}{Observed effect size.}
#'   \item{Lower.crit,Upper.crit}{Critical values for the permutation distribution.}
#'
#' @seealso \code{\link{bwdepth}}, \code{\link{bwdepthMC.ci}}
#' @export
bwdepth.perm<-function(x,y,reps=500,
fun=prodepth,alpha=.05,SEED=TRUE){
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
x=as.matrix(x)
y=as.matrix(y)
n1=nrow(x)
np1=n1+1
n2=nrow(y)
n=n1+n2
d=bwdepth(x,y)$e
xy=rbind(x,y)
print(dim(xy))
v=NA
for(i in 1:reps){
ip=sample(n,replace=FALSE)
z=xy[ip,]
v[i]=bwdepth(z[1:n1,],z[np1:n,])$e
}
v=sort(v)
il=round(alpha*reps/2)
iu=reps-il
list(Est=d,Lower.crit=v[il],Upper.crit=v[iu])
}


#' Helper Function for bwdepthMC.ci
#' @keywords internal
bwdepth.sub<-function(id,x,y,n1,n2,fun){
n=n1+n2
np1=n1+1
e=bwdepth(x[id[1:n1],],y[id[np1:n],],fun=fun)$e
e
}


#' Bagged Depth-Based Classification
#'
#' Bootstrap aggregating (bagging) version of depth-based classification.
#' Reduces bias when sample sizes differ between groups.
#'
#' @param train Training data matrix (optional if `x1` and `x2` provided).
#' @param test Test data matrix to classify.
#' @inheritParams common-params
#' @param x1,x2 Training data for groups 1 and 2 (alternative to `train`/`g`).
#' @param depthfun Depth function to use (default: prodepth).
#' @param DIST Logical. If TRUE, uses distance-based classification.
#' @inheritParams common-params
#' @inheritParams common-params
#' @param ... Additional arguments passed to `depthfun`.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Uses bootstrap bagging to reduce classification bias when n1 ≠ n2. For each
#' bootstrap sample, draws equal-sized samples from both groups and classifies
#' test data. Final classification is by majority vote across bootstrap samples.
#'
#' @seealso \code{\link{discdepth}}, \code{\link{KNNbag}}
#' @export
Depth.class.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,DIST=FALSE,nboot=100,SEED=TRUE,...){
if(SEED)set.seed(2)
if(is.null(test))stop('test =NULL, no test data provided')
if(!is.null(train)){
if(is.null(g))stop('Argument g, group labels, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n1,n,replace=TRUE)
id2=sample(n2,n,replace=TRUE)
dvec[i,]=Depth.class(x1=x1[id1,],x2=x2[id2,],test=test,depthfun=depthfun,DIST=DIST,...)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}



#' Bagged Discriminant Depth Classification
#'
#' Bootstrap aggregating version of `discdepth()` for depth-based discrimination.
#'
#' @inheritParams Depth.class.bag
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Similar to `Depth.class.bag()` but uses `discdepth()` for classification.
#' Bagging reduces bias when sample sizes differ.
#'
#' @seealso \code{\link{Depth.class.bag}}, \code{\link{discdepth}}
#' @export
dis.depth.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,nboot=100,SEED=TRUE,...){
if(is.null(test))stop('test =NULL, no test data provided')
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group labels, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n1,n,replace=TRUE)
id2=sample(n2,n,replace=TRUE)

dvec[i,]=discdepth(x1=x1[id1,],x2=x2[id2,],test=test,depthfun=depthfun,...)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}



#' Bagged Projection-Based Classification
#'
#' Bootstrap aggregating version of `pro.class()` for projection-based
#' classification. Handles unbalanced designs via bagging.
#'
#' @inheritParams Depth.class.bag
#' @param kernel Kernel type for SVM (default: 'radial'). Not used in current implementation.
#' @param nboot Number of bootstrap samples (default: 20).
#' @param PR Logical. If TRUE, uses probability-based classification.
#' @inheritParams common-params
#' @param ... Additional arguments.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Bagging version of `pro.class()`. Reduces bias when n1 ≠ n2 by drawing
#' equal-sized bootstrap samples from each group.
#'
#' @seealso \code{\link{pro.class}}, \code{\link{Depth.class.bag}}
#' @export
pro.class.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,kernel='radial',nboot=20,
PR=FALSE,SEED=TRUE,...){
if(SEED)set.seed(2)
if(is.null(train)){
if(is.null(x1) || is.null(x2))stop('train is null and so are x1 and x2')
x1=as.matrix(x1)
x2=as.matrix(x2)
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(1,n1),rep(2,n2))
train=rbind(x1,x2)
}
else{
if(is.null(g))stop('Argument g, group labels, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n,n,replace=TRUE)
id2=sample(n,n,replace=TRUE)
if(!PR)dvec[i,]=pro.class(x1=x1[id1,],x2=x2[id2,],test=test)
if(PR)dvec[i,]=pro.class.probs(x1=x1[id1,],x2=x2[id2,],test=test)$prob.in.second.class
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
if(!PR)idec=chk2>chk1
if(PR){
chk=apply(dvec,2,mean)
idec=chk>.5
}
dec[idec]=2
dec
}




#' Bagged Projection Depth Classification
#'
#' Bootstrap aggregating version of `pro.classPD()` for projection depth-based
#' classification.
#'
#' @inheritParams Depth.class.bag
#' @param rule Classification rule (default: 1).
#' @inheritParams common-params
#' @param ... Additional arguments.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @seealso \code{\link{pro.classPD}}, \code{\link{pro.class.bag}}
#' @export
pro.classPD.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,rule=1,nboot=100,SEED=TRUE,...){
if(is.null(test))stop('Argument test is null, contains  no data')
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group labels, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n1,n,replace=TRUE)
id2=sample(n2,n,replace=TRUE)
dvec[i,]=pro.classPD(x1=x1[id1,],x2=x2[id2,],test=test,rule=rule)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}




# KNNbag
#' Bagged K-Nearest Neighbor Classification Using Depth
#'
#' Bootstrap aggregating version of K-nearest neighbor classification using data
#' depths. Uses bootstrap sampling to address bias when group sizes differ.
#'
#' @inheritParams Depth.class.bag
#' @param depthfun Depth function to use (default: `prodepth`).
#' @inheritParams common-params
#' @param ... Additional arguments.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Performs KNN classification using data depths via `KNNdist()`. Bootstrap
#' bagging is used to reduce bias when group sizes are unequal (n1 != n2).
#' Each bootstrap sample uses balanced group sizes (min(n1,n2)).
#'
#' @seealso \code{\link{KNNdist}}, \code{\link{Depth.class.bag}}
#' @export
KNNbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,nboot=100,SEED=TRUE,...){
#
#  g=class id
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2.
#  The function will then create appropriate labels and store them in g.
#
# KNN classification using data depths.
# KNNdist uses data depths, for the n1!=n2 it can be a bit biased, meaning that
# when there is no association, the probability of a correct classification will be less than .5
# It removes any row vector with missing values
#
#
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}
Train=cbind(train,g)
Train=elimna(Train)
p=ncol(train)
p1=p+1
train=Train[,1:p]
g=Train[,p1]
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x1=elimna(x1)
x2=elimna(x2)
test=elimna(test)
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n1,n,replace=TRUE)
id2=sample(n2,n,replace=TRUE)
dvec[i,]=KNNdist(x1=x1[id1,],x2=x2[id2,],test=test,depthfun=depthfun)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}

#' Bagged Logistic Regression Classification
#'
#' Bootstrap aggregating version of logistic regression classification using
#' smoothing and a decision rule.
#'
#' @inheritParams Depth.class.bag
#' @param sm Logical. If TRUE, uses smoothing in logistic regression (default: TRUE).
#' @param rule Decision threshold (default: 0.5).
#' @inheritParams common-params
#' @param ... Additional arguments.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Performs logistic regression classification via `class.logR()` with bootstrap
#' bagging. Each bootstrap sample uses balanced group sizes (min(n1,n2)) to
#' reduce bias when group sizes differ.
#'
#' @seealso \code{\link{class.logR}}, \code{\link{Depth.class.bag}}
#' @export
LSMbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,
sm=TRUE,rule=.5,nboot=100,SEED=TRUE,...){
#
#  g=class id
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2.
#  The function will then create appropriate labels and store them in g.
#
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n,n,replace=TRUE)
id2=sample(n,n,replace=TRUE)
dvec[i,]=class.logR(x1=x1[id1,],x2=x2[id2,],test=test,sm=sm,rule=rule)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}




#' Bagged Nearest Neighbor Classification
#'
#' Bootstrap aggregating version of nearest neighbor classification. Uses bootstrap
#' sampling to address bias when group sizes differ.
#'
#' @inheritParams Depth.class.bag
#' @inheritParams common-params
#' @param ... Additional arguments.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Performs nearest neighbor classification via `NN.class()` with bootstrap
#' bagging. Each bootstrap sample uses balanced group sizes (min(n1,n2)) to
#' reduce bias when group sizes are unequal (n1 != n2).
#'
#' @seealso \code{\link{NN.class}}, \code{\link{Depth.class.bag}}
#' @export
NNbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,nboot=100,SEED=TRUE,...){
#
#  g=class id
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2.
#  The function will then create appropriate labels and store them in g.
#
# KNN classification using data depths.
# KNNdist uses data depths, for the n1!=n2 it can be a bit biased, meaning that
# when there is no association, the probability of a correct classification will be less than .5
#
#
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group labels, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n1,n,replace=TRUE)
id2=sample(n2,n,replace=TRUE)
dvec[i,]=NN.class(x1=x1[id1,],x2=x2[id2,],test=test)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}


#' Bagged AdaBoost Classification
#'
#' Bootstrap aggregating version of AdaBoost classification. Uses bootstrap
#' sampling to address bias when group sizes differ (n1 != n2).
#'
#' @inheritParams Depth.class.bag
#' @param baselearner Base learner for AdaBoost (default: 'bbs').
#' @inheritParams common-params
#' @param ... Additional arguments.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Performs AdaBoost classification via `class.ada()` with bootstrap bagging.
#' When group sizes are unequal (n1 != n2) and there is no association, the
#' expected probability of correct classification can differ from 0.5. Bootstrap
#' bagging addresses this bias by using balanced samples (min(n1,n2)).
#'
#' @seealso \code{\link{class.ada}}
#' @keywords internal
 class.ada.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,nboot=100,
SEED=TRUE,baselearner='bbs',...){
#
# class.bag: for n1!=n2
# when there is no association, the expected probability of a correct classification can differ from .5
#
#  This function deals with this via bootstrap bagging
#  g=class labels
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2.
#  The function will then create appropriate labels and store them in g.
#
#  nboot: number of bootstrap sample. Using nboot=20,  bias remains with n1=200, n2=100
#
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group labels, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
#test=as.data.frame(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n1,n,replace=TRUE)
id2=sample(n2,n,replace=TRUE)
dvec[i,]=class.ada(x1=x1[id1,],x2=x2[id2,],test=test,baselearner=baselearner)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}


#' Bagged Random Forest Classification
#'
#' Bootstrap aggregating version of random forest classification. Uses bootstrap
#' sampling to address bias when group sizes differ.
#'
#' @inheritParams Depth.class.bag
#' @param depthfun Depth function to use (default: `prodepth`).
#' @param kernel Kernel type for classification (default: 'radial').
#' @inheritParams common-params
#' @param ... Additional arguments.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Performs random forest classification via `class.forest()` with bootstrap
#' bagging. When group sizes are unequal (n1 != n2) and there is no association,
#' standard methods can have biased misclassification probabilities. Bootstrap
#' bagging with balanced samples (min(n1,n2)) addresses this bias.
#'
#' @seealso \code{\link{class.forest}}, \code{\link{Depth.class.bag}}
#' @export
RFbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,kernel='radial',nboot=100,SEED=TRUE,...){
#
#  g=class id
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2.
#  The function will then create appropriate labels and store them in g.
#
# Random forest classification using data depths.
# class., for the n1!=n2 it can be a bit biased, meaning that
# when there is no association, the probability of a correct classification will be less than .5
#
#
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group labels, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n,n,replace=TRUE)
id2=sample(n,n,replace=TRUE)
xs1=as.data.frame(x1[id1,])
xs2=as.data.frame(x2[id2,])
dvec[i,]=class.forest(x1=xs1,x2=xs2,test=test)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}


#' Bagged Support Vector Machine Classification
#'
#' Bootstrap aggregating version of Support Vector Machine (SVM) classification.
#' Addresses misclassification bias when group sizes differ substantially.
#'
#' @inheritParams Depth.class.bag
#' @param depthfun Depth function to use (default: `prodepth`).
#' @param kernel Kernel type for SVM (default: 'radial').
#' @inheritParams common-params
#' @param ... Additional arguments.
#'
#' @return Vector of predicted class labels (1 or 2) for test observations.
#'
#' @details
#' Performs SVM classification via `SVM()` with bootstrap bagging. Unlike standard
#' SVM, this function addresses a critical bias issue: when n1 != n2 and n2/n1 is
#' small, standard SVM has a misclassification probability of approximately
#' n2/(n1+n2) when there is no true association. This bagged version maintains
#' a probability of 0.5 by using balanced bootstrap samples (min(n1,n2)).
#'
#' @seealso \code{\link{SVM}}, \code{\link{Depth.class.bag}}
#' @export
SVMbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,kernel='radial',nboot=100,SEED=TRUE,...){
#
#  g=class id
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2.
#  The function will then create appropriate labels and store them in g.
#
# Support Vector Machine classification method.
# Unlike standard  SVM this function has the following property. Suppose  n1!=n2 and n2/n1 is small. If there is no
# association between the training data and the labels, the probability of a misclassification is .5
# In contrast, using standard SVM, it is approximately n2/(n1+n2)
#
if(is.null(test))stop('Argument test is null, contains  no data')
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
traing=elimna(cbind(train,g))
p=ncol(train)
p1=p+1
train=traing[,1:p]
test=traing[,p1]
if(length(unique(g))>2)stop('Only two groups allowed, g has more than two unique values')
}
x=fac2list(train,g)
x1=x[[1]]
x2=x[[2]]
}
test=as.matrix(test)
n.test=nrow(test)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(ncol(test)!=ncol(x2))stop('test and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
n=min(c(n1,n2))
dvec=matrix(NA,nrow=nboot,ncol=n.test)
for(i in 1:nboot){
id1=sample(n1,n,replace=TRUE)
id2=sample(n2,n,replace=TRUE)
dvec[i,]=SVM(x1=x1[id1,],x2=x2[id2,],test=test,kernel=kernel)
}
dec=rep(1,n.test)
test1=dvec==1
test2=dvec==2
chk1=apply(test1,2,sum)
chk2=apply(test2,2,sum)
idec=chk2>chk1
dec[idec]=2
dec
}


