# ============================================================================
# Effect Size Functions
# ============================================================================
# 
# This module contains general effect size estimation and reporting functions.
# 
# Core functions:
#   - ES.summary, ES.summary.CI - General effect size summaries
#   - qhat, qhatd - Q statistic (nonparametric effect size)
#   - akp.effect - AKP robust effect size
#   - Various factorial/interaction effect sizes
#
# Dependencies: Depends on functions from utils-core, location, bootstrap
# 
# ============================================================================


#' Confidence Intervals for Decile Differences (Dependent Groups)
#'
#' Computes confidence intervals for the difference between deciles of two
#' dependent groups using the Harrell-Davis quantile estimator with bootstrap.
#' The simultaneous probability coverage is 0.95.
#'
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group (must have same length as \code{x}).
#' @param nboot Number of bootstrap samples (default: 200).
#' @param plotit Logical. If \code{TRUE}, creates a plot showing confidence intervals
#'   and point estimates (default: \code{TRUE}).
#' @param plotop Logical. If \code{TRUE}, plots using quantile positions (0.1, 0.2, ...)
#'   instead of actual decile values (default: \code{FALSE}).
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility (default: \code{TRUE}).
#' @param pr Logical. If \code{TRUE}, prints messages including critical value (default: \code{TRUE}).
#' @param xlab Label for x-axis in plot (default: "x (first group)").
#' @param ylab Label for y-axis in plot (default: "Delta").
#'
#' @return A 9×6 matrix where each row corresponds to the i/10 quantile (i=1,...,9).
#'   Columns are:
#'   \itemize{
#'     \item \code{est.1}: Estimated quantile for group 1
#'     \item \code{est.2}: Estimated quantile for group 2
#'     \item \code{est.dif}: Difference (group 1 - group 2)
#'     \item \code{ci.lower}: Lower confidence limit
#'     \item \code{ci.upper}: Upper confidence limit
#'     \item \code{se}: Bootstrap standard error
#'   }
#'
#' @details
#' Uses the Harrell-Davis estimator for quantiles with bootstrap to construct
#' simultaneous confidence intervals across all deciles. The approximate critical
#' value is computed as 37/n^1.4 + 2.75, providing approximate 0.95 family-wise
#' coverage probability.
#'
#' For more control over alpha level or specific quantiles, use \code{qcomdhd} or
#' \code{qdec2ci}.
#'
#' Missing values are automatically removed.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{qcomdhd}}, \code{\link{qdec2ci}}, \code{\link{hd}}
#'
#' @export
#' @examples
#' # Compare depression scores before and after treatment
#' before <- c(30, 35, 28, 41, 33, 29, 38, 32, 36, 31)
#' after <- c(22, 28, 20, 31, 25, 21, 29, 24, 27, 23)
#' result <- shiftdhd(before, after)
#' # Examine decile differences
#' result
shiftdhd<-function(x,y,nboot=200,plotit=TRUE,plotop=FALSE,SEED=TRUE,pr=TRUE,xlab='x (first group)',
ylab='Delta'){
#
if(pr){
print("NOTE:  if the goal is to use an alpha value different from .05")
print("use the function qcomdhd or  qdec2ci")
}
xy=elimna(cbind(x,y))
x=xy[,1]
y=xy[,2]
plotit<-as.logical(plotit)
if(SEED)set.seed(2) # set seed of random number generator so that
#   results can be duplicated.
crit<-37/length(x)^(1.4)+2.75
if(pr)print("The approximate .05 critical value is")
if(pr)print(crit)
m<-matrix(0,9,6)
if(pr)print("Taking Bootstrap Samples. Please wait.")
data<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
xmat<-matrix(x[data],nrow=nboot,ncol=length(x))
ymat<-matrix(y[data],nrow=nboot,ncol=length(x))
for (i in 1:9){
q<-i/10
bvec<-apply(xmat,1,hd,q)-apply(ymat,1,hd,q)
se<-sqrt(var(bvec))
dif<-hd(x,q)-hd(y,q)
m[i,1]=hd(x,q)
m[i,2]=hd(y,q)
m[i,3]<-dif
m[i,4]<-dif-crit*se
m[i,5]<-dif+crit*se
m[i,6]<-se
}
dimnames(m)<-list(NULL,c('est.1','est.2','est.dif','ci.lower','ci.upper','se'))
if(plotit){
if(plotop){
xaxis<-c(1:9)/10
xaxis<-c(xaxis,xaxis)
}
if(!plotop)xaxis<-c(deciles(x),deciles(x))
par(pch="+")
#yaxis<-c(m[,1],m[,2])
yaxis<-c(m[,4],m[,5])
if(!plotop)plot(xaxis,yaxis,ylab=ylab,xlab=xlab)
if(plotop)plot(xaxis,yaxis,ylab="delta",xlab="Deciles")
par(pch="*")
if(!plotop)points(deciles(x),m[,3])
if(plotop)points(c(1:9)/10,m[,3])
}
m
}


#' Bootstrap Helper for Q-Statistic (Dependent Groups)
#'
#' Internal function used by \code{qhat} when working with bootstrap estimates
#' for dependent groups.
#'
#' @param isubx Bootstrap sample indices.
#' @param x Numeric vector for first group.
#' @param y Numeric vector for second group.
#'
#' @return Vector of predicted group classifications.
#'
#' @keywords internal
#' @noRd
qhatds1<-function(isubx,x,y){
#
xx<-x[isubx]
yy<-y[isubx]
group<-disker(xx,yy,x,op=2)$zhat
group
}


#' Q-Statistic Effect Size (Dependent Groups)
#'
#' Estimates Q, a nonparametric measure of effect size based on prediction accuracy,
#' using the .632 bootstrap method for estimating prediction error. This version is
#' for dependent groups.
#'
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group (same length as \code{x}).
#' @param nboot Number of bootstrap samples (default: 50).
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{app.error}: Apparent classification error rate
#'     \item \code{qhat.632}: Q-statistic using .632 estimator (combines apparent
#'       error with bootstrap cross-validation error)
#'   }
#'
#' @details
#' The Q-statistic measures effect size based on how well observations can be
#' classified into their correct group using kernel density estimation. The .632
#' estimator combines the optimistically biased apparent error rate with the
#' pessimistically biased leave-one-out cross-validation estimate:
#'
#' Q.632 = 0.368 × (apparent error) + 0.632 × (cross-validation error)
#'
#' Values near 0 indicate perfect separation (large effect), values near 0.5
#' indicate no separation (no effect), and values near 1 indicate perfect
#' reverse separation.
#'
#' For dependent groups, the function uses paired observations and a kernel
#' discriminant approach.
#'
#' @references
#' Efron, B. & Tibshirani, R.J. (1993). An Introduction to the Bootstrap, pp. 252-254.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{qhat}} for independent groups, \code{\link{disker}} for
#'   kernel discriminant analysis
#'
#' @export
#' @examples
#' # Compare paired observations
#' before <- rnorm(30, mean=10, sd=2)
#' after <- before - rnorm(30, mean=1, sd=1.5)
#' result <- qhatd(before, after)
#' result$qhat.632  # Smaller values = larger effect
qhatd<-function(x,y,nboot=50){
#
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
#    data is an nboot by n matrix containing subscripts for bootstrap sample
bid<-apply(data,1,idb,length(x))
#  bid is a n by nboot matrix. If the jth bootstrap sample from
#  1, ..., n contains the value i, bid[i,j]=0; otherwise bid[i,j]=1
yhat<-apply(data,1,qhatds1,x,y)
bi<-apply(bid,1,sum) # B sub i in notation of Efron and Tibshirani, p. 253
temp<-(bid*yhat)
diff<-apply(temp,1,sum)
temp<-diff/bi
ep0<-sum(temp[!is.na(temp)])/length(y)
aperror<-disker(x,y)$phat  # apparent error
regpre<-.368*aperror+.632*ep0
list(app.error=aperror,qhat.632=regpre)
}


#' Q-Statistic Effect Size (Independent Groups)
#'
#' Estimates Q, a nonparametric measure of effect size based on prediction accuracy,
#' using the .632 bootstrap method for estimating prediction error. This version is
#' for independent groups.
#'
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group (can be different length than \code{x}).
#' @param nboot Number of bootstrap samples (default: 50).
#' @param op Kernel density estimation method (default: 2):
#'   \itemize{
#'     \item \code{1}: Rosenblatt's shifted histogram
#'     \item \code{2}: Adaptive kernel with expected frequency curve initialization
#'   }
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility (default: \code{TRUE}).
#' @param pr.track Logical. If \code{TRUE}, prints iteration progress (default: \code{FALSE}).
#'
#' @return A list with component:
#'   \itemize{
#'     \item \code{qhat.632}: Q-statistic using .632 estimator (weighted average of
#'       apparent error and cross-validation error)
#'   }
#'
#' @details
#' The Q-statistic measures effect size based on how well observations can be
#' classified into their correct group using kernel density estimation. The .632
#' estimator provides a bias-corrected estimate:
#'
#' Q.632 = 0.368 × (apparent error) + 0.632 × (cross-validation error)
#'
#' where apparent error is optimistically biased (too small) and cross-validation
#' error is pessimistically biased (too large).
#'
#' Interpretation:
#' \itemize{
#'   \item Values near 0: Perfect separation (very large effect size)
#'   \item Values near 0.5: No separation (no effect)
#'   \item Values near 1: Perfect reverse separation
#' }
#'
#' The function computes a pooled estimate by classifying observations from both
#' groups and weighting by sample sizes.
#'
#' Missing values are automatically removed.
#'
#' @references
#' Efron, B. & Tibshirani, R.J. (1993). An Introduction to the Bootstrap, pp. 252-254.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{qhatd}} for dependent groups, \code{\link{disker}} for
#'   kernel discriminant analysis
#'
#' @export
#' @examples
#' # Compare two independent groups
#' group1 <- rnorm(50, mean=10, sd=3)
#' group2 <- rnorm(60, mean=12, sd=3)
#' result <- qhat(group1, group2)
#' result$qhat.632  # Smaller values = larger effect
qhat<-function(x,y,nboot=50,op=2,SEED=TRUE,pr.track=FALSE){
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
datax<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
#    datax is an nboot by n matrix containing subscripts for bootstrap sample
#    associated with first group.
datay<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
#    datay is an nboot by m matrix containing subscripts for bootstrap sample
#    associated with second group.
bidx<-apply(datax,1,idb,n=length(x))
#  bidx is a n by nboot matrix. If the jth bootstrap sample from
#  1, ..., n contains the value i, bid[i,j]=0; otherwise bid[i,j]=1
bidy<-apply(datay,1,idb,n=length(y))
temp3<-matrix(0,ncol=length(x),nrow=nboot)
temp5<-matrix(0,ncol=length(y),nrow=nboot)
for(i in 1:nboot){
temp3[i,]<-disker(x[datax[i,]],y[datay[i,]],x,op=op)$zhat
# temp3 contains vector of 0s and 1s, 1 if x[i] is
#  is classified as coming from group 1.
temp5[i,]<-disker(y[datay[i,]],x[datax[i,]],y,op=op)$zhat
if(pr.track)print(paste("Iteration ", i, "of ", nboot," is complete"))
}
temp4<-temp3*t(bidx)
temp4<-apply(temp4,2,sum)/apply(bidx,1,sum)
temp6<-temp5*t(bidy)
temp6<-apply(temp6,2,sum)/apply(bidy,1,sum)
ep0x<-mean(temp4,na.rm=TRUE)  # epsilon hat_x
aperrorx<-disker(x,y,op=op)$phat  # apparent error
regprex<-.368*aperrorx+.632*ep0x
ep0y<-mean(temp6,na.rm=TRUE)
aperrory<-disker(y,x,op=op)$phat  # apparent error
regprey<-.368*aperrory+.632*ep0y
aperror<-(length(x)*aperrorx+length(y)*aperrory)/(length(x)+length(y))
regpre<-(length(x)*regprex+length(y)*regprey)/(length(x)+length(y))
list(qhat.632=regpre)
}


#' AKP Robust Effect Size (Homoscedastic Cohen's d Analog)
#'
#' Computes a robust analog of Cohen's d based on trimmed means and winsorized
#' variances. This is the homoscedastic version recommended by Algina, Keselman,
#' and Penfield (2005).
#'
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group.
#' @param EQVAR Logical. If \code{TRUE}, assumes equal variances and pools the
#'   standard deviation (default: \code{TRUE}). If \code{FALSE}, computes
#'   heteroscedastic versions using each group's SD separately.
#' @param tr Proportion of observations to trim from each tail (default: 0.2 for
#'   20% trimming on each side).
#'
#' @return If \code{EQVAR=TRUE}, a single effect size value. If \code{EQVAR=FALSE},
#'   a vector of length 2 with effect sizes using group 1 SD and group 2 SD as
#'   denominators.
#'
#' @details
#' The AKP effect size is a robust analog of Cohen's d that uses:
#' \itemize{
#'   \item Trimmed means instead of ordinary means (reduces sensitivity to outliers)
#'   \item Winsorized variances instead of ordinary variances
#'   \item A correction factor to make the expected value under normality equal to
#'     the population effect size
#' }
#'
#' The formula is:
#' \deqn{d_{AKP} = c \times \frac{\bar{X}_t - \bar{Y}_t}{s_p}}
#'
#' where \eqn{\bar{X}_t} and \eqn{\bar{Y}_t} are trimmed means, \eqn{s_p} is the
#' pooled winsorized standard deviation, and \eqn{c} is a correction factor.
#'
#' Under normality with equal population means and variances, this effect size has
#' the same expected value as Cohen's d but is more robust to outliers and
#' non-normality.
#'
#' Rough guidelines (similar to Cohen's d):
#' \itemize{
#'   \item Small effect: |d| = 0.2
#'   \item Medium effect: |d| = 0.5
#'   \item Large effect: |d| = 0.8
#' }
#'
#' @references
#' Algina, J., Keselman, H.J., & Penfield, R.D. (2005). An alternative to Cohen's
#' standardized mean difference effect size: A robust parameter and confidence
#' interval in the two independent groups case. Psychological Methods, 10, 317-328.
#'
#' @seealso \code{\link{akp.effect.ci}} for confidence intervals, \code{\link{KMS.ci}}
#'   for heteroscedastic analog
#'
#' @export
#' @examples
#' # Compare two groups
#' group1 <- c(10, 12, 8, 14, 11, 9, 13, 10, 12, 11)
#' group2 <- c(14, 16, 13, 18, 15, 17, 16, 14, 15, 16)
#' akp.effect(group1, group2)
akp.effect<-function(x,y,EQVAR=TRUE,tr=.2){
#
x<-elimna(x)
y<-elimna(y)
n1<-length(x)
n2<-length(y)
s1sq=winvar(x,tr=tr)
s2sq=winvar(y,tr=tr)
spsq<-(n1-1)*s1sq+(n2-1)*s2sq
sp<-sqrt(spsq/(n1+n2-2))
cterm=1
if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
cterm=sqrt(cterm)
if(EQVAR)dval<-cterm*(tmean(x,tr)-tmean(y,tr))/sp
if(!EQVAR){
dval<-cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s1sq)
dval[2]=cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s2sq)
}
dval
}


#' Bootstrap Confidence Interval for AKP Effect Size
#'
#' Computes bootstrap confidence intervals and p-value for the AKP robust effect
#' size (homoscedastic analog of Cohen's d).
#'
#' @inheritParams common-params
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group.
#' @param alpha Significance level for confidence interval (default: 0.05 for 95% CI).
#' @param tr Proportion of trimming (default: 0.2).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility (default: \code{TRUE}).
#' @param null.val Null hypothesis value for testing (default: 0, indicating no effect).
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{akp.effect}: Point estimate of the AKP effect size
#'     \item \code{ci}: Bootstrap confidence interval (percentile method)
#'     \item \code{p.value}: Two-sided p-value testing H0: effect size = null.val
#'   }
#'
#' @details
#' Uses percentile bootstrap to construct confidence intervals for the AKP effect
#' size. The p-value is computed as twice the minimum of the proportion of bootstrap
#' estimates less than or greater than the null value.
#'
#' The AKP effect size is a robust analog of Cohen's d that uses trimmed means and
#' winsorized variances with a correction factor. See \code{\link{akp.effect}} for
#' details on the effect size measure itself.
#'
#' @references
#' Algina, J., Keselman, H.J., & Penfield, R.D. (2005). An alternative to Cohen's
#' standardized mean difference effect size: A robust parameter and confidence
#' interval in the two independent groups case. Psychological Methods, 10, 317-328.
#'
#' @seealso \code{\link{akp.effect}} for point estimates, \code{\link{ES.summary.CI}}
#'   for multiple effect size measures
#'
#' @export
#' @examples
#' # Compare two groups with CI
#' set.seed(123)
#' group1 <- rnorm(30, mean=100, sd=15)
#' group2 <- rnorm(35, mean=108, sd=15)
#' result <- akp.effect.ci(group1, group2, nboot=500)
#' result$akp.effect
#' result$ci
#' result$p.value
akp.effect.ci<-function(x,y,alpha=.05,tr=.2,nboot=1000,SEED=TRUE,null.val=0){
#
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
be.f=NA
for(i in 1:nboot){
X=sample(x,n1,replace=TRUE)
Y=sample(y,n2,replace=TRUE)
be.f[i]=akp.effect(X,Y,tr=tr)
}
L=alpha*nboot/2
U=nboot-L
be.f=sort(be.f)
ci=be.f[L+1]
ci[2]=be.f[U]
est=akp.effect(x,y,tr=tr)
pv=mean(be.f<null.val)+mean(be.f==null.val)
pv=2*min(c(pv,1-pv))
list(akp.effect=est,ci=ci,p.value=pv)
}


#' Q-Statistic Effect Size Using Data Depth
#'
#' Computes the Q-statistic (proportion correctly classified) using data depth
#' methods for classification. This is for multivariate data.
#'
#' @param x1 Matrix or data frame for the first group (rows are observations,
#'   columns are variables).
#' @param x2 Matrix or data frame for the second group (same number of columns as \code{x1}).
#' @param depthfun Depth function to use for classification (default: \code{prodepth}
#'   for projection depth). See \code{\link{discdepth}} for options.
#' @param ... Additional arguments passed to the depth function.
#'
#' @return The apparent classification accuracy (proportion of observations correctly
#'   classified into their true group). Values near 1 indicate perfect separation
#'   (large effect), values near 0.5 indicate no separation (no effect).
#'
#' @details
#' Uses data depth-based discriminant analysis to classify observations. For each
#' observation, compares its depth in group 1 vs group 2 and assigns it to the
#' group where it has greater depth.
#'
#' The Q-statistic is the proportion of observations correctly classified. This is
#' the "apparent" classification accuracy (computed on the training data itself),
#' which is optimistically biased. For bias-corrected estimates, use \code{\link{qhatdepPB}}
#' which applies bootstrap cross-validation.
#'
#' This depth-based approach works for multivariate data and is robust to outliers.
#' It provides a nonparametric alternative to classical discriminant analysis.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{qhatdepPB}} for bootstrap CI, \code{\link{qhat}} for
#'   univariate kernel-based version, \code{\link{discdepth}} for depth discriminant
#'
#' @export
#' @examples
#' # Two bivariate groups
#' set.seed(123)
#' x1 <- cbind(rnorm(30, mean=0), rnorm(30, mean=0))
#' x2 <- cbind(rnorm(35, mean=1), rnorm(35, mean=1))
#' qhatDEP(x1, x2)  # Higher values = larger effect
qhatDEP<-function(x1,x2,depthfun=prodepth,...){
#
x1<-x1[!is.na(x1)]
x2<-x2[!is.na(x2)]
x1=as.matrix(x1)
x2=as.matrix(x2)
tv=c(rep(1,nrow(x1)),rep(2,nrow(x2)))
see=discdepth(x1,x2,z=rbind(x1,x2))
qhat=mean(tv==see)
qhat
}


#' Bootstrap Confidence Interval for Depth-Based Q-Statistic
#'
#' Computes a bootstrap confidence interval for the Q-statistic (classification
#' accuracy) using data depth methods for multivariate data.
#'
#' @param x1 Matrix or data frame for the first group.
#' @param x2 Matrix or data frame for the second group.
#' @param nboot Number of bootstrap samples (default: 500).
#' @param alpha Significance level for confidence interval (default: 0.05 for 95% CI).
#' @param depthfun Depth function to use (default: \code{prodepth} for projection depth).
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility (default: \code{TRUE}).
#' @param ... Additional arguments passed to the depth function.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{estimate}: Point estimate of Q (apparent classification accuracy)
#'     \item \code{ci}: Bootstrap percentile confidence interval
#'   }
#'
#' @details
#' Uses percentile bootstrap to construct a confidence interval for the depth-based
#' Q-statistic. For each bootstrap sample:
#' \enumerate{
#'   \item Resample with replacement from each group independently
#'   \item Compute the Q-statistic (proportion correctly classified)
#' }
#'
#' The confidence interval is formed using the percentiles of the bootstrap distribution.
#'
#' This provides uncertainty quantification for the multivariate effect size measure
#' based on classification accuracy using data depth.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{qhatDEP}} for point estimates, \code{\link{qhat}} for
#'   univariate version, \code{\link{discdepth}}
#'
#' @export
#' @examples
#' # Two bivariate groups with CI
#' set.seed(123)
#' x1 <- cbind(rnorm(30, mean=0, sd=1), rnorm(30, mean=0, sd=1))
#' x2 <- cbind(rnorm(35, mean=0.8, sd=1), rnorm(35, mean=0.8, sd=1))
#' result <- qhatdepPB(x1, x2, nboot=200)
#' result$estimate
#' result$ci
qhatdepPB<-function(x1,x2,nboot=500,alpha=.05,depthfun=prodepth,
SEED=TRUE,...){
#
if(SEED)set.seed(2)
bvec=NA
x1=as.matrix(x1)
x2=as.matrix(x2)
n1=nrow(x1)
n2=nrow(x2)
for(i in 1:nboot){
dat1=sample(n1,n1,replace=TRUE)
dat2=sample(n2,n2,replace=TRUE)
bvec[i]=qhatDEP(x1[dat1,],x2[dat2,],depthfun=depthfun)
}
est=qhatDEP(x1,x2)
bvec=sort(bvec)
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
ci<-bvec[icl]
ci[2]<-bvec[icu]
list(estimate=est,ci=ci)
}


#' Quantile Shift Effect Size for Linear Contrasts
#'
#' Estimates a quantile shift measure of effect size for a linear contrast by estimating
#' the distribution of the linear contrast and computing P(L <= location(L)), where L
#' is the linear contrast.
#'
#' @param x Data in matrix or list mode containing groups.
#' @param con Vector of contrast coefficients (must sum to zero).
#' @param locfun Location estimator function to use (default: \code{\link{tmean}}).
#' @inheritParams ES.summary
#' @param nreps Number of replications for estimating the linear contrast distribution
#'   (default: 200).
#' @inheritParams qhat
#' @param MAIN Logical. Reserved for future use (default: \code{FALSE}).
#' @param INT Logical. Reserved for future use (default: \code{FALSE}).
#' @param POOL Logical. If \code{FALSE}, estimates the distribution of the weighted sum
#'   of observations. If \code{TRUE}, pools data with positive contrast coefficients,
#'   pools data with negative coefficients, then estimates effect size from the difference
#'   distribution (default: \code{FALSE}).
#' @param nmax Maximum number of pairwise differences to compute. If exceeded, uses
#'   resampling approach (default: 10^8).
#' @param ... Additional arguments passed to \code{locfun}.
#'
#' @return A list with component:
#'   \item{Effect.Size}{The quantile shift effect size estimate.}
#'
#' @details
#' This function estimates the effect size as a quantile shift:
#' \deqn{ES = P(L - \theta_L \le \theta_L)}
#' where L is the linear contrast and \eqn{\theta_L} is the location of L.
#'
#' Two approaches are available:
#' \itemize{
#'   \item \strong{POOL = FALSE}: Estimates the distribution of \eqn{\sum c_i X_i} directly
#'   \item \strong{POOL = TRUE}: Pools observations with c_i = 1, pools those with c_i = -1,
#'     then uses the difference distribution
#' }
#'
#' The quantile shift interpretation: values close to 0.5 indicate no effect, values
#' near 0 or 1 indicate strong effects.
#'
#' Contrast coefficients must sum to zero. Missing values are removed.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{LCES}}, \code{\link{linES.sub}}, \code{\link{lin.akp}}
#'
#' @export
#' @examples
#' # Three-group comparison
#' set.seed(123)
#' g1 <- rnorm(20, mean = 0)
#' g2 <- rnorm(20, mean = 0.5)
#' g3 <- rnorm(20, mean = 1)
#' x <- list(g1, g2, g3)
#'
#' # Linear trend contrast
#' con <- c(-1, 0, 1)
#' lin.ES(x, con)
lin.ES<-function(x,con,locfun=tmean,tr=.2,nreps=200,SEED=TRUE,
MAIN=FALSE,INT=FALSE,POOL=FALSE,nmax=10^8,...){
#
#  Estimate a  quantile shift measure of effect size for a linear contrast based
#  on an estimate of the distribution of the linear contrast.
#
#  con = contrast coefficients
#  x is a matrix or has list mode.
#
#  POOL=FALSE: Estimate distribution of sum c_iX_i, c_i are contrast coefficient
#              Then compute a relative shift measure of effect size.
#  POOL=TRUE: Pool data corresponding to c_i=1, do the same for c_i=-1, then
#  use then compute a relative shift measure of effect size.
#
#
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
nv=as.vector(matl(lapply(x,FUN='length')))
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
nv=as.vector(matl(lapply(x,FUN='length')))
if(length(con)!=J)stop('Length of con should equal number of groups')
x=elimna(x)
if(!POOL){
np=prod(nv)
nmin=min(nv)
if(np>nmax){
nmin=min(c(nmin,100))
}
B=list()
M=matrix(NA,nrow=nmin,ncol=J)
for(i in 1:nreps){
for(j in 1:J)M[,j]=sample(x[[j]],nmin)
B[[i]]=M
}
L=lapply(B,linWMWMC.sub,con=con)
ef.size=NA
for(j in 1:length(L))ef.size[j]=linES.sub(L[[j]],locfun=locfun,...)
ef=mean(ef.size)
}
if(POOL){
y=list()
id1=which(con==1)
id2=which(con==-1)
v1=pool.a.list(x[id1])
v2=pool.a.list(x[id2])
if(length(v1)*length(v2)<nmax){
L=outer(pool.a.list(x[id1]),pool.a.list(x[id2]),FUN='-')
ef=linES.sub(L,locfun=locfun,...)
}
if(length(v1)*length(v2)>=nmax){
B=list()
M=matrix(NA,nrow=nmin,ncol=J)
for(i in 1:nreps){
for(j in 1:J)M[,j]=sample(x[[j]],nmin)
B[[i]]=M
}
L=lapply(B,linWMWMC.sub,con=con)
ef.size=NA
for(j in 1:length(L))ef.size[j]=linES.sub(L[[j]],locfun=locfun,...)
ef=mean(ef.size)
}}
list(Effect.Size=ef)
}


#' @keywords internal
linES.sub<-function(L,locfun,...){
est=locfun(L,...)
ef.size=mean(L-est<=est)
ef.size
}


#' Dependent Groups Linear Contrast Effect Sizes
#'
#' For dependent (repeated measures) groups, computes the Algina et al. robust effect
#' size for each linear contrast based on the linear sum of observations.
#'
#' @param x Data matrix where rows are subjects and columns are repeated measures,
#'   or data in list mode.
#' @param con Matrix of contrast coefficients. Each column defines one contrast (must sum to 0).
#'   If \code{NULL}, generates all pairwise contrasts automatically (default: \code{NULL}).
#'
#' @return A list with components:
#'   \item{con}{The contrast matrix used.}
#'   \item{Effect.Size}{Vector of effect size estimates, one per contrast.}
#'
#' @details
#' For repeated measures designs, this function:
#' \enumerate{
#'   \item For each subject i, computes the linear combination: \eqn{S_i = \sum_j c_j X_{ij}}
#'   \item Computes the Algina et al. robust effect size using \code{\link{D.akp.effect}}
#' }
#'
#' If \code{con = NULL}, automatically generates all J(J-1)/2 pairwise contrasts for
#' J repeated measures.
#'
#' The Algina et al. effect size is a robust analog of Cohen's d for dependent groups.
#'
#' Missing values are removed before computation.
#'
#' @references
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). Effect sizes and their intervals:
#' The two-level repeated measures case. Educational and Psychological Measurement, 65(2), 241-258.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{D.akp.effect}}, \code{\link{lin.ES}}, \code{\link{dep.ES.summary}}
#'
#' @export
#' @examples
#' # Repeated measures: 3 time points
#' set.seed(42)
#' time1 <- rnorm(15, mean = 10)
#' time2 <- time1 + rnorm(15, mean = 0.5, sd = 1)
#' time3 <- time1 + rnorm(15, mean = 1, sd = 1)
#' x <- cbind(time1, time2, time3)
#'
#' # All pairwise comparisons
#' rmlinES(x)
#'
#' # Custom contrast: linear trend
#' con <- matrix(c(-1, 0, 1), ncol = 1)
#' rmlinES(x, con)
rmlinES<-function(x, con = NULL){
#
#  Dependent groups:
#  For each linear contrast, compute Algina et al. effect size based on the linear sum
#
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
J=ncol(x)
if(is.null(con)){
C=(J^2-J)/2
con=matrix(0,ncol=C,nrow=J)
ic=0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic=ic+1
con[j,ic]=1
con[k,ic]=-1
}}}}
x=elimna(x)
n=nrow(x)
ES=NA
for (d in 1:ncol(con)){
S=NA
for(i in 1:n)S[i]=sum(con[,d]*x[i,])
ES[d]=D.akp.effect(S,tr=tr)
}
list(con=con,Effect.Size=ES)
}


#' Determine Small/Medium/Large Effect Size Equivalents
#'
#' Calibrates what constitutes "small", "medium", and "large" effect sizes for various
#' effect size measures by computing their values under specified population mean differences.
#'
#' @param REL.M Numeric vector of length 3 containing the mean differences (in standard
#'   deviation units) to be considered small, medium, and large effects (e.g., c(0.2, 0.5, 0.8)
#'   for Cohen's d conventions).
#' @param n Sample size for simulation (default: 10000; capped at 10000 internally).
#' @param reps Number of replications for averaging (default: 10).
#'
#' @return A 6×3 matrix where rows are effect size measures (AKP, EP, QS (median), QStr,
#'   WMW, KMS) and columns are the three magnitude levels (S, M, L). Each cell contains
#'   the calibrated effect size value.
#'
#' @details
#' This function helps determine what values of various robust effect size measures
#' correspond to "small", "medium", and "large" effects by:
#' \enumerate{
#'   \item Simulating from two standard normal distributions
#'   \item Shifting one distribution by the amounts in \code{REL.M}
#'   \item Computing all effect size measures (AKP, EP, QS, QStr, WMW, KMS)
#'   \item Averaging over \code{reps} replications
#' }
#'
#' Typical usage is with \code{REL.M = c(0.2, 0.5, 0.8)} corresponding to Cohen's d
#' conventions for small, medium, and large effects.
#'
#' The resulting calibration can be used with \code{\link{ES.summary}} and related functions.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{ES.summary}}, \code{\link{ES.summary.CI}}
#'
#' @export
#' @examples
#' # Use Cohen's d conventions
#' ES.sum.REL.MAG(REL.M = c(0.2, 0.5, 0.8), reps = 5)
#'
#' # Custom magnitude definitions
#' ES.sum.REL.MAG(REL.M = c(0.3, 0.6, 1.0), reps = 5)
ES.sum.REL.MAG<-function(REL.M,n = 10000,reps=10){
#
#  Determine small medium and large equivalent measures of effect size based on the values in
#  REL.M
#
if(length(REL.M)!=3)stop('Should have three value in REL.M')
if(n>10000)n=10000
x=rnorm(n)
y=rnorm(n)
output=matrix(0,ncol=3,nrow=6)
int=matrix(NA,ncol=3,nrow=6)
dimnames(output)=list(c('AKP','EP','QS (median)','QStr','WMW','KMS'),c('S','M','L'))
for(k in 1:reps){
for(j in 1:3)int[,j]=ES.summary(x,y-REL.M[j],)[,1]
output=output+int
}
output=output/reps
output
}


#' Comprehensive Effect Size Summary (Independent Groups)
#'
#' Computes six different measures of effect size for comparing two independent
#' groups, providing a comprehensive assessment from multiple perspectives.
#'
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group.
#' @param tr Proportion of trimming for AKP and QStr (default: 0.2).
#' @param NULL.V Vector of null values for each effect size measure (default:
#'   c(0, 0, 0.5, 0.5, 0.5, 0) for AKP, EP, QS, QStr, WMW, KMS respectively).
#' @param REL.MAG Optional vector of length 3 specifying AKP values considered
#'   small, medium, and large under normality (e.g., c(0.2, 0.5, 0.8)). If specified,
#'   equivalent values for other effect sizes are computed via simulation.
#' @param REL.M Alias for REL.MAG (for backward compatibility).
#' @param n.est Number of observations to simulate when computing equivalent effect
#'   sizes for REL.MAG (default: 1000000).
#'
#' @return A 6×5 matrix with rows for each effect size measure and columns:
#'   \itemize{
#'     \item \code{Est}: Point estimate of the effect size
#'     \item \code{NULL}: Null hypothesis value (no effect)
#'     \item \code{S}: Small effect benchmark
#'     \item \code{M}: Medium effect benchmark
#'     \item \code{L}: Large effect benchmark
#'   }
#'
#' Row names indicate the six effect size measures:
#' \itemize{
#'   \item \code{AKP}: Homoscedastic robust analog of Cohen's d (trimmed means, winsorized variance)
#'   \item \code{EP}: Explanatory power (proportion of variance explained)
#'   \item \code{QS (median)}: Quantile shift based on median of X-Y distribution
#'   \item \code{QStr}: Quantile shift based on trimmed mean of X-Y distribution
#'   \item \code{WMW}: Wilcoxon-Mann-Whitney type measure, P(X<Y)
#'   \item \code{KMS}: Robust heteroscedastic analog of Cohen's d
#' }
#'
#' @details
#' This function provides a comprehensive effect size assessment by computing six
#' complementary measures:
#'
#' \strong{Standardized difference measures (like Cohen's d):}
#' \itemize{
#'   \item AKP: Assumes equal variances, robust to outliers via trimming/winsorizing
#'   \item KMS: Allows unequal variances, robust via M-estimator and percentage bend variance
#' }
#'
#' \strong{Probabilistic measures:}
#' \itemize{
#'   \item WMW: Probability that a random observation from X is less than one from Y
#'   \item EP: Proportion of variance in combined data explained by group membership
#' }
#'
#' \strong{Quantile shift measures:}
#' \itemize{
#'   \item QS: Based on median, robust to extreme values
#'   \item QStr: Based on trimmed mean, balances robustness and efficiency
#' }
#'
#' Default benchmarks for small/medium/large effects are based on simulation studies
#' assuming normality. Custom benchmarks can be specified via REL.MAG.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{ES.summary.CI}} for confidence intervals, \code{\link{dep.ES.summary}}
#'   for dependent groups, \code{\link{ESfun}} for individual effect size measures
#'
#' @export
#' @examples
#' # Compare two groups
#' group1 <- rnorm(50, mean=100, sd=15)
#' group2 <- rnorm(60, mean=110, sd=15)
#' ES.summary(group1, group2)
#'
#' # With custom benchmarks
#' ES.summary(group1, group2, REL.MAG=c(0.1, 0.3, 0.5))
ES.summary<-function(x,y,tr=.2,NULL.V=c(0,0,.5,.5,.5,0),REL.MAG=NULL, REL.M=NULL,n.est=1000000){
#
if(!is.null(REL.MAG)) REL.M=REL.MAG
SM=c(.2,.14,.55,.55,.55,.1)
MED=c(.5,.34,.64,.64,.64,.25)
LAR=c(.8,.52,.71,.71,.71,.4)
if(is.null(REL.M)){
SM=c(.2,.14,.55,.55,.55,.1)
MED=c(.5,.34,.64,.64,.64,.25)
LAR=c(.8,.52,.71,.71,.71,.4)
}
if(!is.null(REL.M)){
v=ES.sum.REL.MAG(REL.M,n=n.est)
SM=v[,1]
MED=v[,2]
LAR=v[,3]
SM[1]=REL.M[1]
SM[6]=REL.M[1]/2
MED[1]=REL.M[2]
MED[6]=REL.M[2]/2
LAR[1]=REL.M[3]
LAR[6]=REL.M[3]/2
}
a=c('AKP','EP','QS','QStr','WMW','KMS')
output=matrix(NA,ncol=5,nrow=6)
for(j in 1:6){
output[j,1]=ESfun(x,y,method=a[j],tr=tr,pr=FALSE)
}
output[,2:5]=cbind(NULL.V,SM,MED,LAR)
if(output[1,1]<0)output[1,3:5]=-1*output[1,3:5]
for(j in 3:5){
if(!is.na(output[j,1])){
if(output[j,1]<.5) {
dif=output[j,3:5]-.5
output[j,3:5]=.5-dif
}}}
if(output[6,1]<0)output[6,3:5]=-1*output[6,3:5]
dimnames(output)=list(c('AKP','EP','QS (median)','QStr','WMW','KMS'),c('Est','NULL','S','M','L'))
output
}


#' Comprehensive Effect Size Summary with Confidence Intervals
#'
#' Computes six different measures of effect size for comparing two independent
#' groups along with bootstrap confidence intervals and p-values for each measure.
#'
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group.
#' @param tr Proportion of trimming for AKP and QStr (default: 0.2).
#' @param QSfun Function to use for quantile shift (default: \code{median}).
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility (default: \code{TRUE}).
#' @param alpha Significance level for confidence intervals (default: 0.05 for 95% CI).
#' @param nboot Number of bootstrap samples (default: 2000).
#' @param method P-value adjustment method for multiple comparisons (default: "hoch"
#'   for Hochberg's method). See \code{\link[stats]{p.adjust}} for options.
#' @param NULL.V Vector of null values for each effect size measure (default:
#'   c(0, 0, 0.5, 0.5, 0.5, 0) for AKP, EP, QS, QStr, WMW, KMS respectively).
#' @param REL.MAG Optional vector of length 3 specifying AKP values considered
#'   small, medium, and large under normality. If specified, equivalent values
#'   for other effect sizes are computed via simulation.
#' @param REL.M Alias for REL.MAG (for backward compatibility).
#' @param n.est Number of observations to simulate when computing equivalent effect
#'   sizes for REL.MAG (default: 1000000).
#'
#' @return A 6×9 matrix with rows for each effect size measure and columns:
#'   \itemize{
#'     \item \code{Est}: Point estimate of the effect size
#'     \item \code{NULL}: Null hypothesis value (no effect)
#'     \item \code{S}: Small effect benchmark
#'     \item \code{M}: Medium effect benchmark
#'     \item \code{L}: Large effect benchmark
#'     \item \code{ci.low}: Lower confidence limit
#'     \item \code{ci.up}: Upper confidence limit
#'     \item \code{p.value}: Unadjusted p-value
#'     \item \code{adj.p.value}: P-value adjusted for multiple comparisons
#'   }
#'
#' Row names indicate the six effect size measures: AKP, EP, QS (median), QStr, WMW, KMS.
#'
#' @details
#' This function extends \code{\link{ES.summary}} by adding bootstrap confidence
#' intervals and hypothesis tests for each of the six effect size measures:
#'
#' \itemize{
#'   \item \strong{AKP}: Percentile bootstrap CI
#'   \item \strong{EP}: Percentile bootstrap CI, p-value from Yuen's test
#'   \item \strong{QS}: Percentile bootstrap CI and p-value
#'   \item \strong{QStr}: Percentile bootstrap CI and p-value
#'   \item \strong{WMW}: CI and p-value from Cliff's analog (cidv2)
#'   \item \strong{KMS}: Percentile bootstrap CI and p-value
#' }
#'
#' P-values are adjusted for multiple comparisons using the specified method
#' (default: Hochberg's step-up procedure which controls family-wise error rate).
#'
#' This comprehensive approach allows assessment of effect size from multiple
#' perspectives with proper uncertainty quantification.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{ES.summary}} for point estimates only, \code{\link{dep.ES.summary.CI}}
#'   for dependent groups, \code{\link{akp.effect.ci}}, \code{\link{KMS.ci}}
#'
#' @export
#' @examples
#' # Compare two groups with CIs
#' set.seed(123)
#' group1 <- rnorm(40, mean=100, sd=15)
#' group2 <- rnorm(45, mean=108, sd=15)
#' result <- ES.summary.CI(group1, group2, nboot=500)
#' result
#'
#' # Check which effects are significant after adjustment
#' result[result[,"adj.p.value"] < 0.05, ]
ES.summary.CI<-function(x,y,tr=.2,QSfun=median,SEED=TRUE,alpha=.05,nboot=2000,method='hoch',
NULL.V=c(0,0,.5,.5,.5,0),REL.MAG=NULL,REL.M=NULL,n.est=1000000){
#
a=c('AKP','EP','QS','QStr','WMW','KMS')
output=matrix(NA,ncol=9,nrow=6)
dimnames(output)=list(c('AKP','EP','QS (median)','QStr','WMW','KMS'),c('Est','NULL','S','M','L','ci.low','ci.up','p.value','adj.p.value'))
for(j in 1:6){
output[j,1]=ESfun(x,y,QSfun=QSfun,method=a[j],tr=tr,pr=FALSE,SEED=SEED)
}
if(!is.null(REL.MAG))REL.M=REL.MAG
if(is.null(REL.M)){
SM=c(.2,.14,.55,.55,.55,.1)
MED=c(.5,.34,.64,.64,.64,.25)
LAR=c(.8,.52,.71,.71,.71,.4)
}
if(!is.null(REL.M)){
v=ES.sum.REL.MAG(REL.M,n=n.est)
SM=v[,1]
MED=v[,2]
LAR=v[,3]
SM[1]=REL.M[1]
SM[6]=REL.M[1]/2
MED[1]=REL.M[2]
MED[6]=REL.M[2]/2
LAR[1]=REL.M[3]
LAR[6]=REL.M[3]/2
}
output[,2:5]=cbind(NULL.V,SM,MED,LAR)
if(output[1,1]<0)output[1,3:5]=-1*output[1,3:5]
for(j in 3:4){
if(output[j,1]<.5){
dif=output[j,3:5]-.5
output[j,3:5]=.5-dif
}}
if(output[5,1]>.5 & output[5,5]< .5){
output[5,3:5]=1-output[5,3:5]
}
if(output[5,1]<.5 & output[5,5]> .5){
output[5,3:5]=1-output[5,3:5]
}

a=akp.effect.ci(x,y,tr=tr,alpha=alpha,nboot=nboot,SEED=SEED)
output[1,6:7]=a$ci
output[1,8]=a$p.value
a=EPci(x,y,tr=tr,alpha=alpha,SEED=SEED,nboot=nboot)
output[2,6:7]=a$ci
output[2,8]=yuen(x,y,tr=tr)$p.value
a=shiftPBci(x,y,locfun=QSfun,alpha=alpha,nboot=nboot,SEED=SEED)
output[3,6:7]=a$ci
output[3,8]=a$p.value
a=shiftPBci(x,y,locfun=tmean,alpha=alpha,nboot=nboot,SEED=SEED)
output[4,6:7]=a$ci
output[4,8]=a$p.value
a=cidv2(x,y,alpha=alpha)
output[5,6:7]=a$p.ci
output[5,8]=a$p.value
a=KMS.ci(x,y,alpha=alpha,nboot=nboot,SEED=SEED)
output[6,6:7]=a$ci
output[6,8]=a$p.value
if(output[6,1]<0)output[6,3:5]=-1*output[6,3:5]
output[,9]=p.adjust(output[,8],method=method)
output
}


#' Multivariate Effect Size Summary
#'
#' For multivariate data (multiple variables), computes comprehensive effect size
#' summaries for each variable separately when comparing two groups.
#'
#' @param x1 Matrix or data frame for the first group (rows are observations,
#'   columns are variables).
#' @param x2 Matrix or data frame for the second group (must have same number of
#'   columns as \code{x1}).
#'
#' @return A list with p elements (one per variable), where each element contains
#'   the output from \code{\link{ES.summary}} for that variable (a 6×5 matrix with
#'   six effect size measures: AKP, EP, QS, QStr, WMW, KMS).
#'
#' @details
#' This function applies \code{\link{ES.summary}} to each variable (column)
#' independently. For variable j, it computes six effect size measures comparing
#' group 1 vs group 2 on that variable:
#' \itemize{
#'   \item AKP: Robust Cohen's d analog
#'   \item EP: Explanatory power
#'   \item QS (median): Quantile shift based on median
#'   \item QStr: Quantile shift based on trimmed mean
#'   \item WMW: P(X<Y) probability
#'   \item KMS: Heteroscedastic Cohen's d analog
#' }
#'
#' This provides a comprehensive marginal (univariate) effect size assessment for
#' each variable in a multivariate comparison. It does not assess multivariate
#' effects directly.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{ES.summary}} for univariate effect sizes, \code{\link{ES.summary.CI}}
#'   for confidence intervals
#'
#' @export
#' @examples
#' # Two groups, three variables
#' set.seed(123)
#' x1 <- cbind(rnorm(30, 100, 15), rnorm(30, 50, 10), rnorm(30, 25, 5))
#' x2 <- cbind(rnorm(35, 105, 15), rnorm(35, 52, 10), rnorm(35, 25, 5))
#' result <- MUL.ES.sum(x1, x2)
#' result[[1]]  # Effect sizes for variable 1
#' result[[2]]  # Effect sizes for variable 2
#' result[[3]]  # Effect sizes for variable 3
MUL.ES.sum<-function(x1,x2){
#
V=list()
p=ncol(x1)
for(j in 1:p)V[[j]]=ES.summary(x1[,j],x2[,j])
V
}


#' Pairwise Effect Sizes for Independent Groups
#'
#' Computes effect size measures for all pairwise comparisons (or specified contrasts)
#' among J independent groups. For each contrast, pools observations with coefficient
#' +1 and compares them to pooled observations with coefficient -1.
#'
#' @param x Data in matrix, data frame, or list mode. Each column (or list element)
#'   represents one group.
#' @param con Matrix of contrast coefficients with J rows (one per group). Each column
#'   defines one contrast. If \code{NULL} (default), performs all pairwise comparisons.
#' @param fun Function to compute effect sizes (default: \code{ES.summary}). Can also
#'   use \code{ES.summary.CI} for confidence intervals or \code{ESfun} for a single
#'   effect size measure.
#' @param tr Proportion of trimming (default: 0.2). Passed to the effect size function.
#' @param ... Additional arguments passed to \code{fun}.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{con}: Matrix of contrast coefficients used
#'     \item \code{effect.size}: List with one element per contrast, each containing
#'       the output from \code{fun}
#'   }
#'
#' @details
#' For each column of the contrast matrix \code{con}:
#' \enumerate{
#'   \item Groups with coefficient +1 are pooled together
#'   \item Groups with coefficient -1 are pooled together
#'   \item The specified effect size function is applied to compare the two pooled groups
#' }
#'
#' By default, uses \code{fun=ES.summary} which computes six effect size measures
#' (AKP, EP, QS, QStr, WMW, KMS) for each comparison.
#'
#' For individual effect size measures, use \code{fun=ESfun} with the \code{method}
#' argument:
#' \itemize{
#'   \item \code{method='EP'}: Explanatory power
#'   \item \code{method='QS'}: Quantile shift (median)
#'   \item \code{method='QStr'}: Quantile shift (trimmed mean)
#'   \item \code{method='AKP'}: Robust Cohen's d analog
#'   \item \code{method='WMW'}: P(X<Y) probability
#'   \item \code{method='KMS'}: Heteroscedastic Cohen's d analog
#' }
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{ES.summary}}, \code{\link{ES.summary.CI}}, \code{\link{ESfun}},
#'   \code{\link{DEP.PAIR.ES}} for dependent groups, \code{\link{con.all.pairs}}
#'
#' @export
#' @examples
#' # Three independent groups, all pairwise comparisons
#' set.seed(123)
#' group1 <- rnorm(20, mean=100, sd=15)
#' group2 <- rnorm(20, mean=105, sd=15)
#' group3 <- rnorm(20, mean=110, sd=15)
#' x <- cbind(group1, group2, group3)
#'
#' # Get comprehensive effect sizes for all pairs
#' result <- IND.PAIR.ES(x)
#' result$effect.size[[1]]  # Compare groups 1 vs 2
#'
#' # Get only AKP effect size
#' result2 <- IND.PAIR.ES(x, fun=ESfun, method='AKP')
IND.PAIR.ES<-function(x,con=NULL,fun=ES.summary,tr=.2,...){
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
CON=list()
J=length(x)
if(is.null(con))conmat=con.all.pairs(J)
else conmat=con
P=ncol(conmat)
ic=0
for(j in 1:P){
id1=which(conmat[,j]==1)
id2=which(conmat[,j]==-1)
ic=ic+1
if(length(id1)==1)CON[[ic]]=fun(x[[id1]],x[[id2]],tr=tr,...)
else{
z1=pool.a.list(x[id1])
z2=pool.a.list(x[id2])
CON[[ic]]=fun(z1,z2,...)
}}
list(con=conmat,effect.size=CON)
}


#' Row-Column Effect Sizes for Factorial Design
#'
#' For a J×K factorial design, computes effect sizes for all pairwise comparisons
#' among levels of Factor B (columns) within each level of Factor A (rows), and
#' vice versa.
#'
#' @param J Number of levels for Factor A (rows).
#' @param K Number of levels for Factor B (columns).
#' @param x Data in matrix, data frame, or list mode. If matrix/data frame, should
#'   have JK columns arranged as: Factor B levels 1...K for Factor A level 1, then
#'   Factor B levels 1...K for Factor A level 2, etc.
#' @param fun Function to compute effect sizes (default: \code{ES.summary}). See
#'   \code{\link{IND.PAIR.ES}} for options.
#' @param tr Proportion of trimming (default: 0.2). Passed to the effect size function.
#' @param ... Additional arguments passed to \code{fun}.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{A}: List with J elements. Element j contains effect sizes for
#'       all pairwise comparisons among levels of Factor B at level j of Factor A.
#'     \item \code{B}: List with K elements. Element k contains effect sizes for
#'       all pairwise comparisons among levels of Factor A at level k of Factor B.
#'   }
#'
#' @details
#' This function provides a comprehensive effect size analysis for a two-way
#' factorial design by computing:
#' \enumerate{
#'   \item For each level of Factor A: All pairwise effect sizes among Factor B levels
#'   \item For each level of Factor B: All pairwise effect sizes among Factor A levels
#' }
#'
#' This is useful for understanding simple effects and interaction patterns. If
#' Factor B effect sizes differ across levels of Factor A (or vice versa), this
#' suggests an interaction.
#'
#' The \code{fun} argument controls which effect size measures are computed. Default
#' is \code{ES.summary} which provides six measures (AKP, EP, QS, QStr, WMW, KMS)
#' for each comparison.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{IND.PAIR.ES}} for the pairwise comparison engine,
#'   \code{\link{ES.summary}} for effect size measures, \code{\link{inter.ES}}
#'   for interaction effect sizes
#'
#' @export
#' @examples
#' # 2x3 factorial design
#' set.seed(123)
#' J <- 2  # Factor A levels
#' K <- 3  # Factor B levels
#' # Create data: interaction present
#' x <- list()
#' x[[1]] <- rnorm(20, 100, 15)  # A1,B1
#' x[[2]] <- rnorm(20, 105, 15)  # A1,B2
#' x[[3]] <- rnorm(20, 110, 15)  # A1,B3
#' x[[4]] <- rnorm(20, 100, 15)  # A2,B1
#' x[[5]] <- rnorm(20, 100, 15)  # A2,B2
#' x[[6]] <- rnorm(20, 100, 15)  # A2,B3
#'
#' result <- RCES(J, K, x)
#' result$A[[1]]  # Factor B comparisons at Factor A level 1
#' result$B[[1]]  # Factor A comparisons at Factor B level 1
RCES<-function(J,K,x,fun=ES.summary,tr=.2,...){
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
JK=J*K
imat=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
A=list()
B=list()
for(j in 1:J){
z=x[imat[j,]]
A[[j]]=IND.PAIR.ES(z,tr=tr,fun=fun,...)
}
for(k in 1:K){
z=x[imat[,k]]
B[[k]]=IND.PAIR.ES(z,tr=tr,fun=fun,...)
}
list(A=A,B=B)
}


#' Interaction Effect Size for 2×2 Design
#'
#' Computes effect size measures for an interaction in a 2×2 factorial design by
#' comparing the distribution of Factor B differences at level 1 of Factor A vs
#' level 2 of Factor A.
#'
#' @param x Data in matrix, data frame, or list mode with 4 groups arranged as:
#'   A1B1, A1B2, A2B1, A2B2.
#' @param method Effect size method to use (default: "EP" for explanatory power).
#'   Options:
#'   \itemize{
#'     \item \code{'DNT'}: De Neve and Thas method (WMW-type)
#'     \item \code{'PH'}: Patel-Hoel method (Cliff's analog)
#'     \item \code{'EP'}: Explanatory power
#'     \item \code{'QS'}: Quantile shift based on median
#'     \item \code{'QStr'}: Quantile shift based on trimmed mean
#'     \item \code{'AKP'}: Robust Cohen's d analog
#'     \item \code{'WMW'}: Wilcoxon-Mann-Whitney type measure
#'     \item \code{'KMS'}: Heteroscedastic Cohen's d analog
#'   }
#' @param iter Number of iterations for sampling when dataset is large (default: 5).
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility (default: \code{TRUE}).
#' @param tr Proportion of trimming for AKP and QStr methods (default: 0.2).
#' @param pr Logical. If \code{TRUE}, prints progress (default: \code{FALSE}).
#'
#' @return Effect size estimate for the interaction. Interpretation depends on the
#'   method chosen (see Details).
#'
#' @details
#' For a 2×2 design, an interaction means the Factor B effect differs across levels
#' of Factor A. This function quantifies interaction effect size by:
#' \enumerate{
#'   \item At Factor A level 1: Estimate distribution of all pairwise differences
#'     between Factor B levels (X_{A1B1} - X_{A1B2})
#'   \item At Factor A level 2: Estimate distribution of all pairwise differences
#'     between Factor B levels (X_{A2B1} - X_{A2B2})
#'   \item Compute an effect size comparing these two distributions
#' }
#'
#' If the total number of pairwise differences exceeds 1000, the function samples
#' subsets and averages the effect size estimates across iterations.
#'
#' Different methods provide different perspectives:
#' \itemize{
#'   \item \strong{DNT, PH, WMW}: Probabilistic measures (e.g., P(D1 < D2))
#'   \item \strong{EP}: Proportion of variance explained by the interaction
#'   \item \strong{QS, QStr}: Quantile shift measures
#'   \item \strong{AKP, KMS}: Standardized difference measures (like Cohen's d)
#' }
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{interES.2by2}} for multiple effect sizes simultaneously,
#'   \code{\link{ESfun}} for individual effect size computation, \code{\link{RCES}}
#'   for simple effects
#'
#' @export
#' @examples
#' # 2x2 design with interaction
#' set.seed(123)
#' x <- list()
#' x[[1]] <- rnorm(20, 100, 15)  # A1,B1
#' x[[2]] <- rnorm(20, 110, 15)  # A1,B2 (large B effect at A1)
#' x[[3]] <- rnorm(20, 100, 15)  # A2,B1
#' x[[4]] <- rnorm(20, 102, 15)  # A2,B2 (small B effect at A2)
#'
#' # Different effect size measures for the interaction
#' inter.ES(x, method='EP')   # Explanatory power
#' inter.ES(x, method='AKP')  # Robust Cohen's d analog
#' inter.ES(x, method='QS')   # Quantile shift
inter.ES<-function(x,method='EP',iter=5,SEED=TRUE,tr=.2,pr=FALSE){
#
if(is.matrix(x))x=listm(x)
ef=NA
if(length(x)!=4)stop('Limited to a two-by-two design')
x=elimna(x)
FLAG=TRUE
if(method=='DNT'){
ef=WMWinter.est(x,iter=iter,SEED=SEED)
FLAG=FALSE
}
if(method=='PH'){
ef=rimul(2,2,x)$test[,5]
FLAG=FALSE
}
if(FLAG){
nv=as.vector(matl(lapply(x,FUN='length')))
nt=prod(nv)
if(nt>10^3){
if(SEED)set.seed(2)
Nmin1=min(c(nv[1],nv[2],100))
Nmin2=min(c(nv[3],nv[4],100))
for(i in 1:iter){
id1=sample(nv[1],Nmin1)
id2=sample(nv[2],Nmin1)
L1=outer(x[[1]][id1],x[[2]][id2],FUN='-')
id1=sample(nv[3],Nmin2)
id2=sample(nv[4],Nmin2)
L2=outer(x[[3]][id1],x[[4]][id2],FUN='-')
ef[i]=ESfun(L1,L2,method=method,tr=tr,pr=pr,SEED=FALSE)
}}
else{
L1=outer(x[[1]],x[[2]],FUN='-')
L2=outer(x[[3]],x[[4]],FUN='-')
ef=ESfun(L1,L2,method=method,tr=tr,pr=pr)
}
}
ef=mean(ef)
ef
}


#' @keywords internal
inter.TDES.sub<-function(x,method='QS',iter=5,SEED=TRUE,tr=.2,pr=FALSE,switch=FALSE){
#
#
# Measures of effect size for an  interaction in a 2-by-2 design
# For  level 1 of Factor A, estimate the distribution of the
#  the typical difference for levels 1  and 2 Factor B
#  Do the same for level 2 of Factor A, and compute  a measure of
#  effect size based on these two distributions.
#
#  swithch=TRUE, interchange the rows and columns
#
# Choices for the argument method:
# 'DNT',`EP',`QS',`QStr',`AKP',`KMS'
# DNT= De Neve and Thas P(X_1-X_2 < X_3-X_4) so a WMW-type measure
# EP=explanatory power,
#  QS= quantile shift (median,
#  QStr= quantile shift (trimmed mean) ,
#  AKP =trimmed mean version of Cohen's d,
# KMS=heteroscedastic analog of Cohen's d
#
#
if(is.matrix(x))x=listm(x)
if(switch)x=x[1,3,2,4]
ef=NA
if(length(x)!=4)stop('Limited to a two-by-two design')
x=elimna(x)
FLAG=TRUE
if(method=='DNT'){
ef=WMWinter.est(x,iter=iter,SEED=SEED)
FLAG=FALSE
}
if(FLAG){
nv=as.vector(matl(lapply(x,FUN='length')))
nt=prod(nv)
if(nt>10^3){
if(SEED)set.seed(2)
Nmin1=min(c(nv[1],nv[2],100))
Nmin2=min(c(nv[3],nv[4],100))
for(i in 1:iter){
id1=sample(nv[1],Nmin1)
id2=sample(nv[2],Nmin1)
L1=outer(x[[1]][id1],x[[2]][id2],FUN='-')
id1=sample(nv[3],Nmin2)
id2=sample(nv[4],Nmin2)
L2=outer(x[[3]][id1],x[[4]][id2],FUN='-')
ef[i]=ESfun(L1,L2,method=method,tr=tr,pr=pr,SEED=FALSE)
}}
else{
L1=outer(x[[1]],x[[2]],FUN='-')
L2=outer(x[[3]],x[[4]],FUN='-')
ef=ESfun(L1,L2,method=method,tr=tr,pr=pr)
}
}
ef=mean(ef)
ef
}


#' Interaction Effect Sizes for 2×2 Designs
#'
#' Computes a collection of effect size measures for interactions in a 2×2 factorial
#' design. Effect sizes are computed for each row separately, then differenced to
#' assess the interaction.
#'
#' @param x Data in matrix, data frame, or list mode. For a 2×2 design, should contain
#'   4 groups in row-major order: (1,1), (1,2), (2,1), (2,2).
#' @inheritParams ES.summary
#' @param SW Logical. If \code{TRUE}, transposes the design (reverses rows and columns)
#'   (default: \code{FALSE}).
#'
#' @return A 6×4 matrix with rows corresponding to different effect size measures:
#'   \itemize{
#'     \item AKP: Homoscedastic robust analog of Cohen's d
#'     \item EP: Explanatory power
#'     \item QS (median): Quantile shift based on median of X-Y distribution
#'     \item QStr: Quantile shift based on trimmed mean of X-Y distribution
#'     \item KMS: Robust heteroscedastic analog of Cohen's d
#'     \item PH (WMW): Patel-Hoel using Cliff's analog of Wilcoxon-Mann-Whitney
#'   }
#'   Columns are:
#'   \itemize{
#'     \item NULL: Null value for each measure
#'     \item Est 1: Effect size for first row (compare columns within row 1)
#'     \item Est 2: Effect size for second row (compare columns within row 2)
#'     \item Diff: Interaction estimate (Est 1 - Est 2)
#'   }
#'
#' @details
#' For a 2×2 design, this function:
#' \enumerate{
#'   \item Computes the effect size comparing columns within the first row
#'   \item Computes the effect size comparing columns within the second row
#'   \item Takes the difference to estimate the interaction
#' }
#'
#' The interaction is assessed as the difference between the effect of Factor B at
#' level 1 of Factor A versus the effect of Factor B at level 2 of Factor A.
#'
#' Missing values are removed before computation.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{inter.ES}}, \code{\link{interJK.ESmul}}, \code{\link{ESfun}}
#'
#' @export
#' @examples
#' # 2x2 design: Factor A (treatment), Factor B (time)
#' set.seed(42)
#' # Group 1,1: Control at Time 1
#' g11 <- rnorm(20, mean = 10, sd = 2)
#' # Group 1,2: Control at Time 2
#' g12 <- rnorm(20, mean = 11, sd = 2)
#' # Group 2,1: Treatment at Time 1
#' g21 <- rnorm(20, mean = 10, sd = 2)
#' # Group 2,2: Treatment at Time 2 (larger increase)
#' g22 <- rnorm(20, mean = 14, sd = 2)
#'
#' x <- list(g11, g12, g21, g22)
#' interES.2by2(x)
interES.2by2<-function(x,tr=.2,SW=FALSE){
#
# Estimate a collection of effect sizes
# for the first row of a 2-by-2 design
# do the same for the second row
# return estimates of the differences
#
#  AKP: Homoscedastic robust analog of Cohen's d
#  EP:  Explanatory power
#  QS:  Quantile shift based on the median of the distribution of X-Y,
#  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
#  KMS:  Robust heteroscedastic analog  of Cohen's d
#  PH: Patel--Hoel, uses Cliff'a analog of Wilcoxon--Mann--Whitney
#
#  switch=TRUE: reverses rows and columns

if(is.matrix(x)  || is.data.frame(x))x=listm(x)
if(SW)x=x[c(1,3,2,4)]
J=length(x)
if(J!=4)stop('Should have four groups; designed for a 2-by-2 ANOVA only')
a=c('AKP','EP','QS','QStr','KMS','WMW')
output=matrix(NA,ncol=4,nrow=6)
output[,1]=c(0.0,0.0,0.5,0.5,0.0,0.5)
for(j in 1:6){
output[j,2]=ESfun(x[[1]],x[[2]],method=a[j],tr=tr,pr=FALSE)
output[j,3]=ESfun(x[[3]],x[[4]],method=a[j],tr=tr,pr=FALSE)
output[j,4]=output[j,2]-output[j,3]
}
dimnames(output)=list(c('AKP','EP','QS (median)','QStr',
'KMS','PH'),c('NULL','Est 1','Est 2','Diff'))
output
}


#' Interaction Effect Sizes for J×K Designs
#'
#' Computes effect size measures for all tetrad interactions in a J×K factorial design.
#' For each 2×2 sub-design (tetrad), computes the interaction effect size using
#' the specified method.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in matrix or list mode. Expected to contain J*K groups in row-major order.
#' @param method Effect size method to use. See \code{\link{ESfun}} for options. Common choices:
#'   \code{'QS'} (quantile shift, median-based), \code{'AKP'}, \code{'EP'}, \code{'KMS'}, \code{'WMW'}
#'   (default: \code{'QS'}).
#' @inheritParams ES.summary
#' @inheritParams qhat
#'
#' @return A list with component:
#'   \item{EFFECT.est}{A matrix with one row per tetrad interaction. Columns are:
#'     \itemize{
#'       \item Factor A, Factor A: The two levels of Factor A being compared
#'       \item Factor B, Factor B: The two levels of Factor B being compared
#'       \item Effect Size 1: Effect size for first level of Factor A
#'       \item Effect Size 2: Effect size for second level of Factor A
#'       \item Diff: Interaction effect size (difference between the two)
#'     }
#'   }
#'
#' @details
#' For a J×K design, this function examines all possible 2×2 sub-designs (tetrads).
#' For each tetrad defined by rows j, j' and columns k, k':
#' \enumerate{
#'   \item Compute the effect of Factor B (columns k vs k') at Factor A level j
#'   \item Compute the effect of Factor B (columns k vs k') at Factor A level j'
#'   \item Take the difference to estimate the interaction
#' }
#'
#' The number of tetrad interactions examined is C(J,2) × C(K,2) = [J(J-1)/2] × [K(K-1)/2].
#'
#' Missing values are automatically removed.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{interES.2by2}}, \code{\link{inter.ES}}, \code{\link{ESfun}}
#'
#' @export
#' @examples
#' # 3x2 design
#' set.seed(42)
#' x <- list(
#'   rnorm(15, 10), rnorm(15, 11),  # Row 1
#'   rnorm(15, 10), rnorm(15, 13),  # Row 2 (larger effect)
#'   rnorm(15, 10), rnorm(15, 12)   # Row 3
#' )
#' result <- interJK.ESmul(J = 3, K = 2, x, method = 'QS')
#' result$EFFECT.est
interJK.ESmul<-function(J,K,x,method='QS',tr=.2,SEED=TRUE){
#
#  Compute measures of effect size for interactions associated with
#  in J-by-K design.
#  This is done for all relevant tetrad cells using interES.2by2
#  Missing values are automatically removed.
#
#  Methods, see the R function ESfun
#  Defaults to quantile shfit
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
CCJ<-(J^2-J)/2
CCK<-(K^2-K)/2
CC<-CCJ*CCK
JK=J*K
test<-matrix(NA,CC,7)
x=elimna(x)
mat=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
dimnames(test)<-list(NULL,c("Factor A","Factor A","Factor B","Factor B","Effect Size  1","Effect Size  2","Diff"))
jcom<-0
for (j in 1:J){
for (jj in 1:J){
if (j < jj){
for (k in 1:K){
for (kk in 1:K){
if (k < kk){
jcom<-jcom+1
test[jcom,1]<-j
test[jcom,2]<-jj
test[jcom,3]<-k
test[jcom,4]<-kk
id1=mat[j,k]
id2=mat[j,kk]
a=ESfun(x[[id1]],x[[id2]],method=method,tr=tr,pr=FALSE,SEED=SEED)
id1=mat[jj,k]
id2=mat[jj,kk]
b=ESfun(x[[id1]],x[[id2]],method=method,tr=tr,pr=FALSE,SEED=SEED)
test[jcom,5:7]<-c(a,b,a-b)
}}}}}}
list(EFFECT.est=test)
}


#' Linear Contrast Effect Sizes (Multiple Measures)
#'
#' For each linear contrast (column of \code{con}), computes four different measures
#' of effect size: quantile shift (median-based), quantile shift (trimmed mean-based),
#' AKP robust Cohen's d, and a sign-test analog.
#'
#' @param x Data in matrix, data frame, or list mode containing groups.
#' @param con Matrix of contrast coefficients. Each column defines one contrast (must sum to 0).
#' @param nreps Number of replications for estimating the distribution of linear contrasts
#'   (default: 200).
#' @inheritParams lin.ES
#'
#' @return A list with components:
#'   \item{EST}{A 4×(ncol(con)+1) matrix. Rows are effect size measures (QS, Qstr, AKP, SIGN).
#'     First column contains null values; remaining columns contain effect size estimates
#'     for each contrast.}
#'   \item{con}{The contrast matrix used.}
#'
#' @details
#' This function computes four different effect size measures for each linear contrast:
#' \enumerate{
#'   \item \strong{QS}: Quantile shift based on the median of the linear contrast distribution
#'   \item \strong{Qstr}: Quantile shift based on the trimmed mean of the linear contrast distribution
#'   \item \strong{AKP}: Robust generalization of Cohen's d using \code{\link{lin.akp}}
#'   \item \strong{SIGN}: Analog of the sign test via \code{\link{linsign}}
#' }
#'
#' The quantile shift measures estimate P(L <= location(L)), where L is the linear contrast.
#' Null values represent the expected value under no effect (typically 0.5 for quantile
#' shifts, 0.0 for AKP).
#'
#' Contrasts must sum to zero. Missing values are removed before computation.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{lin.ES}}, \code{\link{lin.akp}}, \code{\link{linsign}}
#'
#' @export
#' @examples
#' # Three-group comparison
#' set.seed(123)
#' g1 <- rnorm(20, mean = 0)
#' g2 <- rnorm(20, mean = 0.5)
#' g3 <- rnorm(20, mean = 1)
#' x <- list(g1, g2, g3)
#'
#' # Pairwise contrasts
#' con <- matrix(c(1, -1, 0,
#'                  1, 0, -1,
#'                  0, 1, -1), ncol = 3)
#' LCES(x, con)
LCES<-function(x,con,nreps=200,tr=.2,SEED=TRUE){
#
# For each column of con,  compute four measures of effect size:
# quantile shift based on median
# quantile shift based on a trimmed mean
# AKP generalization of Cohen's d
# SIGN: analog of the sign test.
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
x=elimna(x)
con=as.matrix(con)
d=ncol(con)
mat=matrix(NA,nrow=4,ncol=d)
LAB=NULL
for(i in 1:d){
LAB[i]=paste('Con',i)
mat[1,i]=lin.ES(x,as.vector(con[,i]),locfun=median,nreps=nreps)$Effect.Size
mat[2,i]=lin.ES(x,as.vector(con[,i]),locfun=mean,nreps=nreps,tr=tr)$Effect.Size
mat[3,i]=lin.akp(x,con[,i],locfun=mean,nreps=nreps,tr=tr)$Effect.Size
mat[4,i]=linsign(x,con[,i],nreps=nreps)
}
mat=cbind(c(0.5,0.5,0.0,0.5),mat)
LAB=c('NULL',LAB)
dimnames(mat)=list(c('QS','Qstr','AKP','SIGN'),LAB)
list(EST=mat,con=con)
}


#' Quantile Estimation (Alternative to Harrell-Davis)
#'
#' Estimates the qth quantile using an alternative method that can offer advantages
#' over the Harrell-Davis estimator when comparing extreme quantiles with heavy-tailed
#' distributions.
#'
#' @param x Numeric vector.
#' @param q Quantile to estimate (default: 0.5 for median).
#'
#' @return The estimated quantile value (scalar).
#'
#' @details
#' This estimator uses a binomial-weighted approach that can provide better performance
#' than the Harrell-Davis estimator in certain situations, particularly when:
#' \itemize{
#'   \item Comparing extreme quantiles (e.g., q < 0.1 or q > 0.9)
#'   \item Distributions have heavy tails
#' }
#'
#' The method weights sorted observations using binomial probabilities, with special
#' handling for small sample sizes (n <= 2).
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{hd}} for the Harrell-Davis quantile estimator
#'
#' @export
#' @examples
#' # Estimate median
#' x <- c(1, 3, 5, 7, 9, 100)  # outlier present
#' qno.est(x, q = 0.5)
#'
#' # Compare with Harrell-Davis for extreme quantile
#' qno.est(x, q = 0.9)
#' hd(x, q = 0.9)
qno.est<-function(x,q=.5){
#
#  Estimate of the qth quantile
#  In some situations, offers a distinct advantage over the Harrell-Davis estimator when
# comparing extreme quantiles and distributions have heavy tails.
#
n<-length(x)
x<-sort(x)
s<-numeric()
ifelse(n>2, {for(g in 1:(n-2)){
s[g]<-x[g+1]*(dbinom(g, size=n, prob=q)*(1-q)+dbinom(g+1, size=n, prob=q)*q)
}
sum(s,na.rm = TRUE)
t1<-(2*dbinom(0, size=n, prob=q)*q+dbinom(1, size=n, prob=q)*q)*x[1]
t2<-(2*(1-q)*dbinom(n, size=n, prob=q)+dbinom(n-1, size=n, prob=q)*(1-q))*x[n]
t3<-dbinom(0, size=n, prob=q)*(2-3*q)*x[2]-dbinom(0, size=n, prob=q)*(1-q)*x[3]-
dbinom(n, size=n, prob=q)*q*x[n-2]+dbinom(n, size=n, prob=q)*(3*q-1)*x[n-1]
quan<-sum(s,na.rm = T)+t1+t2+t3},
ifelse(n==2,{quan <- (1-q)*x[1]+q*x[2]},quan<-x))
quan
}


#' Between-Within Design: Factor A Effect Sizes
#'
#' For a between-within (mixed) factorial design, computes effect sizes for all
#' pairwise comparisons of Factor A (between-subjects) at each level of Factor B
#' (within-subjects).
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (within-subjects factor).
#' @param x Data in matrix or list mode. Expected to contain J*K groups in row-major order.
#'   \code{x[[1]]} is level (1,1), \code{x[[2]]} is level (1,2), ..., \code{x[[K]]} is level (1,K),
#'   \code{x[[K+1]]} is level (2,1), etc.
#' @inheritParams ES.summary
#' @param pr Logical. If \code{TRUE}, prints informational messages (default: \code{TRUE}).
#' @param fun Effect size summary function to use (default: \code{\link{ES.summary}}).
#'   Use \code{\link{ES.summary.CI}} for confidence intervals.
#' @param ... Additional arguments passed to \code{fun} (e.g., \code{REL.MAG} for custom
#'   small/medium/large calibration).
#'
#' @return A list with K components (one per level of Factor B):
#'   \item{B[[k]]}{Effect size results for all pairwise comparisons of Factor A at
#'     level k of Factor B, as returned by \code{\link{IND.PAIR.ES}}.}
#'
#' @details
#' For each level of the within-subjects Factor B, this function computes pairwise
#' effect sizes comparing all levels of the between-subjects Factor A.
#'
#' The structure is:
#' \itemize{
#'   \item \code{B[[1]]}: Effect sizes for Factor A comparisons at Factor B level 1
#'   \item \code{B[[2]]}: Effect sizes for Factor A comparisons at Factor B level 2
#'   \item ...and so on for all K levels
#' }
#'
#' Effect sizes can be customized via the \code{fun} argument and \code{REL.MAG} parameter
#' (e.g., \code{REL.MAG = c(0.1, 0.3, 0.5)} for custom small/medium/large calibration).
#'
#' The function assumes \code{fun} has a \code{tr} argument for trimming.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{bw.es.B}}, \code{\link{bw.es.I}}, \code{\link{IND.PAIR.ES}}, \code{\link{ES.summary}}
#'
#' @export
#' @examples
#' # 2x3 mixed design: 2 between-subjects groups, 3 time points
#' set.seed(42)
#' x <- list(
#'   rnorm(15, 10), rnorm(15, 11), rnorm(15, 12),  # Group 1, times 1-3
#'   rnorm(15, 10), rnorm(15, 12), rnorm(15, 15)   # Group 2, times 1-3 (larger increase)
#' )
#' result <- bw.es.A(J = 2, K = 3, x)
#' # B[[1]]: Between-group effect at time 1
#' # B[[2]]: Between-group effect at time 2
#' # B[[3]]: Between-group effect at time 3
bw.es.A<-function(J,K,x,tr=.2,pr=TRUE,fun=ES.summary,...){
#
#  Between-by-within design.
#
# Assumed that fun function has an argument tr
#
#Using REL.M, can change default values for small, medium and large
# Example REL.M=c(.1,.3,.5)
#
# For each level of Factor B,  compute  effect sizes
# for all pairs of  levels of  Factor A .
#
#  The R variable x is assumed to contain the raw
#  x stored in list mode. x[[1]] contains the x
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the x for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the x for level 1,K
#  x[[K+1]] is the x for level 2,1, x[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested, but a subset of the x can be analyzed
#  using grp
#
#
if(pr){
if(!identical(fun,ES.summary.CI))print('To get confidence intervals, set the argument fun=ES.summary.CI')
}
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK=J*K
mat=matrix(c(1:JK),J,K,byrow=TRUE)
B=list()
for(k in 1:K){
B[[k]]=IND.PAIR.ES(x[mat[,k]],fun=fun,tr=tr,...)
if(k==1){
if(pr){
print('B[[1]] contains pairwise measures of effect size for all levels of Factor A')
print(' and level 1 of Factor B')
print(' B[[2]]   contains pairwise measures of effect size for all levels of Factor A')
print('and level 2 of Factor B')
}}}
list(B=B)
}


#' Effect Size Summary for Dependent Groups
#'
#' Computes four different measures of effect size for comparing two dependent
#' groups (repeated measures or matched pairs) based on difference scores.
#'
#' @param x Numeric vector. If \code{y} is provided, this is the first group/time point.
#'   If \code{y=NULL}, this should contain the difference scores.
#' @param y Numeric vector for the second group/time point, or \code{NULL} if \code{x}
#'   contains difference scores (default: \code{NULL}).
#' @param tr Proportion of trimming for AKP (default: 0.2).
#' @param alpha Significance level (default: 0.05) - currently not used in this version.
#' @param REL.MAG Optional vector of length 3 specifying AKP values considered small,
#'   medium, and large under normality (default: c(0.1, 0.3, 0.5)). Equivalent values
#'   for other effect sizes are computed via simulation.
#' @param SEED Logical. If \code{TRUE}, sets random seed for simulations (default: \code{TRUE}).
#' @param nboot Number of bootstrap samples (default: 2000) - currently not used in this version.
#'
#' @return A 4×5 matrix with rows for each effect size measure and columns:
#'   \itemize{
#'     \item \code{NULL}: Null hypothesis value (no effect)
#'     \item \code{Est}: Point estimate of the effect size
#'     \item \code{S}: Small effect benchmark
#'     \item \code{M}: Medium effect benchmark
#'     \item \code{L}: Large effect benchmark
#'   }
#'
#' Row names indicate the four effect size measures:
#' \itemize{
#'   \item \code{AKP}: Robust standardized difference (analog of Cohen's d for paired data)
#'   \item \code{QS (median)}: Quantile shift based on median of difference scores
#'   \item \code{QStr}: Quantile shift based on trimmed mean of difference scores
#'   \item \code{SIGN}: P(X<Y), probability that first observation is less than second
#' }
#'
#' @details
#' For dependent groups, effect sizes are based on the distribution of difference
#' scores (X - Y). The four measures assess different aspects:
#'
#' \itemize{
#'   \item \strong{AKP}: Standardized mean difference using robust estimators (trimmed
#'     mean and winsorized variance). Analogous to Cohen's d but more robust.
#'   \item \strong{QS}: Proportion of difference distribution below the median. Values
#'     < 0.5 indicate X tends to be larger than Y.
#'   \item \strong{QStr}: Similar to QS but based on trimmed mean, balancing robustness
#'     and efficiency.
#'   \item \strong{SIGN}: Probability that a randomly selected pair has X < Y. Related
#'     to the sign test.
#' }
#'
#' Default benchmarks assume normality and are based on simulation. Custom benchmarks
#' can be specified via REL.MAG.
#'
#' For confidence intervals and p-values, use \code{\link{dep.ES.summary.CI}}.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{dep.ES.summary.CI}} for confidence intervals, \code{\link{ES.summary}}
#'   for independent groups, \code{\link{D.akp.effect}}, \code{\link{depQS}}
#'
#' @export
#' @examples
#' # Pre-post comparison
#' before <- c(85, 90, 88, 92, 87, 89, 91, 86, 90, 88)
#' after <- c(88, 93, 90, 95, 91, 92, 94, 89, 93, 91)
#' dep.ES.summary(before, after)
#'
#' # Using difference scores directly
#' diff_scores <- before - after
#' dep.ES.summary(diff_scores)
dep.ES.summary<-function(x,y=NULL,tr=.2, alpha=.05, REL.MAG=NULL,SEED=TRUE,nboot=2000){
#

ecom=c(0.10,  0.54,  0.54,  0.46, 0.30,  0.62,
   0.62,  0.38, 0.50,  0.69,  0.69,  0.31)
REL.EF=matrix(NA,4,3)

if(!is.null(y))x=x-y
x=elimna(x)
n=length(x)
output=matrix(NA,ncol=5,nrow=4)
dimnames(output)=list(c('AKP','QS (median)','QStr','SIGN'),c('NULL','Est','S','M','L'))
output[1,1:2]=c(0,D.akp.effect(x,tr=tr))
output[2,1:2]=c(0.5,depQS(x)$Q.effect)
output[3,1:2]=c(0.5,depQS(x,locfun=mean,tr=tr)$Q.effect)
output[4,1:2]=c(0.5,mean(x[x!=0]<0))
if(is.null(REL.MAG)){
REL.MAG=c(.1,.3,.5)
REL.EF=matrix(ecom,4,3)
}
else{
if(length(REL.MAG)!=3)stop('REL.MAG should have three values')
if(SEED)set.seed(2)
z=rnorm(1000000)
v=dep.ES.summary.sub(z+REL.MAG[1],tr=tr)
REL.EF[,1]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[2],tr=tr)
REL.EF[,2]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[3],tr=tr)
REL.EF[,3]=v[,2]
}
if(output[1,2]<0)REL.EF[1,]=0-REL.EF[1,]
if(output[2,2]<0.5)REL.EF[2,]=.5-(REL.EF[2,]-.5)
if(output[3,2]<0.5)REL.EF[3,]=.5-(REL.EF[3,]-.5)
if(output[4,2]>0.5)REL.EF[4,]=.5-(REL.EF[4,]-.5)
output[,3:5]=REL.EF
output
}


#' @keywords internal
dep.ES.summary.sub<-function(x,y=NULL,tr=.2){
#
#
#  Used  to determine equivalent effect size based on specified standard deviations
#

if(!is.null(y))x=x-y
output=matrix(NA,ncol=2,nrow=4)
dimnames(output)=list(c('AKP','QS (median)','QStr','SIGN'),c('NULL','Est'))
output[1,1:2]=c(0,D.akp.effect(x))
output[2,1:2]=c(0.5,depQS(x)$Q.effect)
output[3,1:2]=c(0.5,depQS(x,locfun=mean,tr=tr)$Q.effect)
output[4,1:2]=c(0.5,mean(x[x!=0]<0))
output
}


#' Effect Size Summary with Confidence Intervals (Dependent Groups)
#'
#' Computes four different measures of effect size for comparing two dependent
#' groups along with bootstrap confidence intervals and p-values for each measure.
#'
#' @param x Numeric vector. If \code{y} is provided, this is the first group/time point.
#'   If \code{y=NULL}, this should contain the difference scores.
#' @param y Numeric vector for the second group/time point, or \code{NULL} if \code{x}
#'   contains difference scores (default: \code{NULL}).
#' @param tr Proportion of trimming for AKP (default: 0.2).
#' @param alpha Significance level for confidence intervals (default: 0.05 for 95% CI).
#' @param REL.MAG Optional vector of length 3 specifying AKP values considered small,
#'   medium, and large under normality (default: c(0.1, 0.3, 0.5)). Equivalent values
#'   for other effect sizes are computed via simulation.
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility (default: \code{TRUE}).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param AUTO Logical. Passed to binomial confidence interval computation for SIGN
#'   effect size (default: \code{FALSE}).
#'
#' @return A 4×8 matrix with rows for each effect size measure and columns:
#'   \itemize{
#'     \item \code{NULL}: Null hypothesis value (no effect)
#'     \item \code{Est}: Point estimate of the effect size
#'     \item \code{S}: Small effect benchmark
#'     \item \code{M}: Medium effect benchmark
#'     \item \code{L}: Large effect benchmark
#'     \item \code{ci.low}: Lower confidence limit
#'     \item \code{ci.up}: Upper confidence limit
#'     \item \code{p.value}: P-value for testing no effect
#'   }
#'
#' Row names indicate the four effect size measures: AKP, QS (median), QStr, SIGN.
#'
#' @details
#' This function extends \code{\link{dep.ES.summary}} by adding bootstrap confidence
#' intervals and hypothesis tests for each of the four effect size measures based on
#' difference scores:
#'
#' \itemize{
#'   \item \strong{AKP}: Percentile bootstrap CI and p-value from D.akp.effect.ci
#'   \item \strong{QS (median)}: Percentile bootstrap CI and p-value from depQSci
#'   \item \strong{QStr}: Percentile bootstrap CI and p-value based on trimmed mean
#'   \item \strong{SIGN}: Exact binomial CI and p-value (P(X<Y) based on signs)
#' }
#'
#' The benchmarks for small/medium/large effects are computed by simulating from
#' normal distributions with specified standardized mean differences (REL.MAG) and
#' determining equivalent values for the other effect size measures.
#'
#' This provides a comprehensive assessment of paired data effect sizes with proper
#' uncertainty quantification.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{dep.ES.summary}} for point estimates only, \code{\link{ES.summary.CI}}
#'   for independent groups, \code{\link{D.akp.effect.ci}}, \code{\link{depQSci}}
#'
#' @export
#' @examples
#' # Pre-post comparison with CIs
#' set.seed(123)
#' before <- rnorm(25, mean=100, sd=15)
#' after <- before + rnorm(25, mean=5, sd=10)
#' result <- dep.ES.summary.CI(before, after, nboot=500)
#' result
#'
#' # Check significance
#' result[result[,"p.value"] < 0.05, ]
dep.ES.summary.CI<-function(x,y=NULL,tr=.2, alpha=.05, REL.MAG=NULL,SEED=TRUE,nboot=1000,AUTO=FALSE){
#

ecom=c(0.10,  0.54,  0.54,  0.46, 0.30,  0.62,
   0.62,  0.38, 0.50,  0.69,  0.69,  0.31)
REL.EF=matrix(NA,4,3)

if(!is.null(y))x=x-y
x=elimna(x)
n=length(x)
output=matrix(NA,ncol=8,nrow=4)
dimnames(output)=list(c('AKP','QS (median)','QStr','SIGN'),c('NULL','Est','S','M','L','ci.low','ci.up','p.value'))
output[1,1:2]=c(0,D.akp.effect(x))
output[2,1:2]=c(0.5,depQS(x)$Q.effect)
output[3,1:2]=c(0.5,depQS(x,locfun=mean,tr=tr)$Q.effect)
output[4,1:2]=c(0.5,mean(x[x!=0]<0))
if(is.null(REL.MAG)){
REL.MAG=c(.1,.3,.5)
REL.EF=matrix(ecom,4,3)
}
else{
if(length(REL.MAG)!=3)stop('REL.MAG should have three values')
if(SEED)set.seed(2)
z=rnorm(1000000)
v=dep.ES.summary.sub(z+REL.MAG[1],tr=tr)
REL.EF[,1]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[2],tr=tr)
REL.EF[,2]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[3],tr=tr)
REL.EF[,3]=v[,2]
}
if(output[1,2]<0)REL.EF[1,]=0-REL.EF[1,]
if(output[2,2]<0.5)REL.EF[2,]=.5-(REL.EF[2,]-.5)
if(output[3,2]<0.5)REL.EF[3,]=.5-(REL.EF[3,]-.5)
if(output[4,2]>0.5)REL.EF[4,]=.5-(REL.EF[4,]-.5)
output[,3:5]=REL.EF
a=D.akp.effect.ci(x,alpha=alpha,SEED=SEED,tr=tr,nboot=nboot)
output[1,6:7]=a$ci
output[1,8]=a$p.value
#output[1,6:7]=D.akp.effect.ci(x,alpha=alpha,SEED=SEED,tr=tr,nboot=nboot)$ci
a=depQSci(x,alpha=alpha,SEED=SEED,nboot=nboot)
output[2,6:7]=a$ci
output[2,8]=a$p.value
a=depQSci(x,locfun=tmean,alpha=alpha, SEED=SEED,tr=tr,nboot=nboot)
output[3,6:7]=a$ci
output[3,8]=a$p.value
Z=sum(x<0)
nm=length(x[x!=0])
a=binom.conf.pv(Z,nm,alpha=alpha,AUTO=AUTO,pr=FALSE)
output[4,6:7]=a$ci
output[4,8]=a$p.value
output
}


#' Between-Within Design: Factor B Effect Sizes
#'
#' For a between-within (mixed) factorial design, computes effect sizes for Factor B
#' (within-subjects) comparisons at each level of Factor A (between-subjects), with options
#' to pool data or compute confidence intervals.
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (within-subjects factor).
#' @param x Data in matrix or list mode. Expected to contain J*K groups in row-major order.
#'   \code{x[[1]]} is level (1,1), \code{x[[2]]} is level (1,2), ..., \code{x[[K]]} is level (1,K),
#'   \code{x[[K+1]]} is level (2,1), etc.
#' @inheritParams ES.summary
#' @param POOL Logical. If \code{TRUE}, pools data over levels of Factor A before computing
#'   effect sizes for Factor B (default: \code{FALSE}).
#' @param OPT Logical. Reserved for future use (default: \code{FALSE}).
#' @param CI Logical. If \code{TRUE}, computes confidence intervals for effect sizes
#'   (default: \code{FALSE}).
#' @inheritParams qhat
#' @param REL.MAG Optional vector of length 3 for custom small/medium/large effect size
#'   calibration (e.g., \code{c(0.1, 0.3, 0.5)}) (default: \code{NULL} uses standard values).
#' @param pr Logical. If \code{TRUE}, prints informational messages (default: \code{TRUE}).
#'
#' @return A list with component:
#'   \item{A}{If \code{POOL = FALSE}, a list of length J containing effect size results for
#'     Factor B comparisons at each level of Factor A (via \code{\link{DEP.PAIR.ES}}).
#'     If \code{POOL = TRUE}, a single set of results for pooled data.}
#'
#' @details
#' This function examines the effect of the within-subjects Factor B:
#' \itemize{
#'   \item \strong{POOL = FALSE}: For each level j of Factor A, computes all pairwise
#'     comparisons among the K levels of Factor B
#'   \item \strong{POOL = TRUE}: Pools all subjects across Factor A levels, then computes
#'     pairwise comparisons for Factor B
#' }
#'
#' Use \code{CI = TRUE} to obtain bootstrap confidence intervals for effect sizes.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{bw.es.A}}, \code{\link{bw.es.I}}, \code{\link{DEP.PAIR.ES}}
#'
#' @export
#' @examples
#' # 2x3 mixed design: 2 between-subjects groups, 3 time points
#' set.seed(42)
#' n <- 15
#' # Create correlated repeated measures
#' baseline <- rnorm(n)
#' x <- list(
#'   baseline + rnorm(n, 0, 0.5),         # Group 1, time 1
#'   baseline + rnorm(n, 0.5, 0.5),       # Group 1, time 2
#'   baseline + rnorm(n, 1, 0.5),         # Group 1, time 3
#'   baseline + rnorm(n, 0, 0.5),         # Group 2, time 1
#'   baseline + rnorm(n, 1, 0.5),         # Group 2, time 2
#'   baseline + rnorm(n, 2, 0.5)          # Group 2, time 3
#' )
#' result <- bw.es.B(J = 2, K = 3, x)
#' # A[[1]]: Time effects for Group 1
#' # A[[2]]: Time effects for Group 2
bw.es.B<-function(J,K,x,tr=.2,POOL=FALSE,OPT=FALSE,CI=FALSE,SEED=TRUE,REL.MAG=NULL,pr=TRUE){
#
#  Between-by-within design.
#
# For each level of Factor A,  compute  effect sizes
# for all j<m levels of  Factor B .
#
# If POOL=TRUE, pool the data over levels of Factor B
#
# If the data are stored in a matrix or data frame, it is converted to list mode.
#   Once in list,
# x[[1]] contains the x
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the x for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the x for level 1,K
#  x[[K+1]] is the x for level 2,1, x[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested, but a subset of the x can be analyzed
#  using grp
#
#
if(pr){
if(!CI)print('To get confidence intervals, set the argument CI=TRUE')
}
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK=J*K
mat=matrix(c(1:JK),J,K,byrow=TRUE)
if(POOL){
Y=list()
for(k in 1:K){
Y[[k]]=pool.a.list(x[mat[,k]])
}
A=DEP.PAIR.ES(Y,tr=tr,SEED=SEED,REL.MAG=REL.MAG,CI=CI)
}
if(!POOL){
if(OPT){
REL.EF=matrix(NA,nrow=4,ncol=3)
if(!is.null(REL.MAG))
if(length(REL.MAG)!=3)stop('REL.MAG should have three values')
if(SEED)set.seed(2)
z=rnorm(1000000)
v=dep.ES.summary.sub(z+REL.MAG[1],tr=tr)
REL.EF[,1]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[2],tr=tr)
REL.EF[,2]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[3],tr=tr)
REL.EF[,3]=v[,2]
}
A=list()
for(j in 1:J){
A[[j]]=DEP.PAIR.ES(x[mat[j,]],tr=tr,SEED=SEED,REL.MAG=REL.MAG,CI=CI)
}}
list(A=A)
}


#' Between-Within Design: Interaction Effect Sizes
#'
#' For a between-within (mixed) factorial design, computes interaction effect sizes
#' by comparing difference scores (within-subjects effects) between levels of the
#' between-subjects factor.
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (within-subjects factor).
#' @param x Data in matrix or list mode. Expected to contain J*K groups in row-major order.
#'   \code{x[[1]]} is level (1,1), \code{x[[2]]} is level (1,2), ..., \code{x[[K]]} is level (1,K),
#'   \code{x[[K+1]]} is level (2,1), etc.
#' @inheritParams ES.summary
#' @param OPT Logical. Reserved for future use (default: \code{FALSE}).
#' @inheritParams qhat
#' @param CI Logical. If \code{TRUE}, computes confidence intervals for interaction effect sizes
#'   (default: \code{FALSE}).
#' @param alpha Significance level for confidence intervals when \code{CI = TRUE} (default: 0.05).
#' @param REL.MAG Optional vector of length 3 for custom small/medium/large effect size
#'   calibration (default: \code{NULL}).
#' @param pr Logical. If \code{TRUE}, prints informational messages (default: \code{TRUE}).
#'
#' @return A list with components:
#'   \item{Interaction.ES}{List of effect size results for each interaction contrast,
#'     as returned by \code{\link{ES.summary}} or \code{\link{ES.summary.CI}}.}
#'   \item{con}{The contrast matrix defining the interactions (from \code{\link{con2way}}).}
#'
#' @details
#' This function assesses interactions by:
#' \enumerate{
#'   \item Using \code{\link{con2way}} to generate all interaction contrasts
#'   \item For each interaction contrast, computing difference scores for Factor B within
#'     different levels of Factor A
#'   \item Comparing these difference scores between levels of Factor A using
#'     \code{\link{ES.summary}} or \code{\link{ES.summary.CI}}
#' }
#'
#' Each interaction compares whether the effect of Factor B differs across levels of
#' Factor A.
#'
#' Use \code{CI = TRUE} to obtain bootstrap confidence intervals.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{bw.es.A}}, \code{\link{bw.es.B}}, \code{\link{con2way}}, \code{\link{ES.summary}}
#'
#' @export
#' @examples
#' # 2x2 mixed design with interaction
#' set.seed(42)
#' n <- 20
#' # Group 1: small time effect
#' g1_t1 <- rnorm(n, 10, 2)
#' g1_t2 <- g1_t1 + rnorm(n, 0.5, 1)
#' # Group 2: large time effect (interaction)
#' g2_t1 <- rnorm(n, 10, 2)
#' g2_t2 <- g2_t1 + rnorm(n, 2, 1)
#'
#' x <- list(g1_t1, g1_t2, g2_t1, g2_t2)
#' result <- bw.es.I(J = 2, K = 2, x)
#' result$Interaction.ES
bw.es.I<-function(J,K,x,tr=.2,OPT=FALSE,SEED=TRUE,CI=FALSE, alpha=.05, REL.MAG=NULL,pr=TRUE){
#
#  Between-by-within design.
#
# Effect size based on an interaction: compare difference scores
#
# If POOL=TRUE, pool the data over levels of Factor B
#
# If the data are stored in a matrix or data frame, it is converted to list mode.
#   Once in list,
# x[[1]] contains the x
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the x for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the x for level 1,K
#  x[[K+1]] is the x for level 2,1, x[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested, but a subset of the x can be analyzed
#  using grp
#
#
if(pr){
if(!CI)print('To get confidence intervals, set the argument CI=TRUE')
}
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK=J*K
mat=matrix(c(1:JK),J,K,byrow=TRUE)
Int=list()
mat=con2way(J,K)$conAB
NC=ncol(mat)
for(j in 1:NC){
id=which(mat[,j]!=0)
dif1=x[[id[1]]]-x[[id[2]]]
dif2=x[[id[3]]]-x[[id[4]]]
if(!CI)Int[[j]]=ES.summary(dif1,dif2,tr=tr)
if(CI)Int[[j]]=ES.summary.CI(dif1,dif2,tr=tr,alpha=alpha)
}
list(Interaction.ES=Int,con=mat)
}


#' Within-Within Design: Effect Sizes for Factor B
#'
#' For a within-within (fully repeated measures) factorial design, computes effect sizes
#' for all pairwise comparisons of Factor B at each level of Factor A.
#'
#' @param J Number of levels for Factor A (first within-subjects factor).
#' @param K Number of levels for Factor B (second within-subjects factor).
#' @param x Data in matrix or list mode. Expected to contain J*K repeated measures in
#'   row-major order. \code{x[[1]]} is level (1,1), \code{x[[2]]} is level (1,2), etc.
#' @inheritParams ES.summary
#' @param CI Logical. If \code{TRUE}, computes confidence intervals for effect sizes
#'   (default: \code{FALSE}).
#' @inheritParams qhat
#' @param REL.MAG Optional vector of length 3 for custom small/medium/large effect size
#'   calibration (default: \code{NULL}).
#'
#' @return A list of length J, where each element contains effect size results for all
#'   pairwise comparisons of Factor B at that level of Factor A (via \code{\link{DEP.PAIR.ES}}).
#'
#' @details
#' For a fully repeated measures J×K design, this function computes effect sizes for
#' Factor B (within-subjects) at each level of Factor A (within-subjects).
#'
#' The structure is:
#' \itemize{
#'   \item Element 1: All pairwise comparisons among K levels of Factor B at Factor A level 1
#'   \item Element 2: All pairwise comparisons among K levels of Factor B at Factor A level 2
#'   \item ...and so on for all J levels
#' }
#'
#' This is similar to \code{\link{bw.es.B}} but for a fully within-subjects design where
#' all participants experience all J*K conditions.
#'
#' Use \code{CI = TRUE} to obtain bootstrap confidence intervals.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{bw.es.B}}, \code{\link{DEP.PAIR.ES}}, \code{\link{wwlin.es}}
#'
#' @export
#' @examples
#' # 2x3 within-within design: 2 tasks, 3 time points
#' # All subjects complete all 6 conditions
#' set.seed(42)
#' n <- 15
#' baseline <- rnorm(n)
#' x <- list(
#'   baseline + rnorm(n, 0, 0.5),    # Task 1, time 1
#'   baseline + rnorm(n, 0.5, 0.5),  # Task 1, time 2
#'   baseline + rnorm(n, 1, 0.5),    # Task 1, time 3
#'   baseline + rnorm(n, 0, 0.5),    # Task 2, time 1
#'   baseline + rnorm(n, 1, 0.5),    # Task 2, time 2
#'   baseline + rnorm(n, 2, 0.5)     # Task 2, time 3
#' )
#' result <- ww.es(J = 2, K = 3, x)
#' # [[1]]: Time comparisons for Task 1
#' # [[2]]: Time comparisons for Task 2
ww.es<-function(J,K,x,tr=.2,CI=FALSE,SEED=TRUE,REL.MAG=NULL){
#
#  within-by-within design.
#
# For each level of Factor A,  compute  effect sizes
# for all j<m levels of  Factor B .
#
#
# If the data are stored in a matrix or data frame, it is converted to list mode.
#   Once in list,
# x[[1]] contains the x
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the x for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the x for level 1,K
#  x[[K+1]] is the x for level 2,1, x[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested
#
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK=J*K
mat=matrix(c(1:JK),J,K,byrow=TRUE)
if(!is.null(REL.MAG)){
REL.EF=matrix(NA,nrow=4,ncol=3)
if(length(REL.MAG)!=3)stop('REL.MAG should have three values')
if(SEED)set.seed(2)
z=rnorm(1000000)
v=dep.ES.summary.sub(z+REL.MAG[1],tr=tr)
REL.EF[,1]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[2],tr=tr)
REL.EF[,2]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[3],tr=tr)
REL.EF[,3]=v[,2]
}
A=list()
B=list()
AB=list()
for(j in 1:J)A[[j]]=DEP.PAIR.ES(x[mat[j,]],tr=.2,SEED=SEED,REL.MAG=REL.MAG,CI=CI)
for(k in 1:K)B[[k]]=DEP.PAIR.ES(x[mat[,k]],tr=.2,SEED=SEED,REL.MAG=REL.MAG,CI=CI)
MAT=con2way(J,K)$conAB
for(k in 1:ncol(MAT)){
id=which(MAT[,k]!=0)
dif1=x[[id[1]]]-x[[id[2]]]
dif2=x[[id[3]]]-x[[id[4]]]
if(!CI)AB[[k]]=ES.summary(dif1,dif2)
if(CI)AB[[k]]=ES.summary.CI(dif1,dif2)
}
list(Factor.A=A,Factor.B=B,INT=AB,conA=con.all.pairs(J),conB=con.all.pairs(K),AB=MAT)
}


#' Dependent Groups Linear Contrast Effect Sizes with Confidence Intervals
#'
#' For dependent (repeated measures) groups, computes four measures of effect size with
#' bootstrap confidence intervals for each linear contrast of the J variables.
#'
#' @param x Data matrix where rows are subjects and columns are repeated measures,
#'   or data in list mode.
#' @param con Matrix of contrast coefficients. Each column defines one contrast (must sum to 0).
#'   If \code{NULL}, generates all pairwise contrasts automatically (default: \code{NULL}).
#' @inheritParams ES.summary
#' @param REL.MAG Optional vector of length 3 specifying AKP values considered small,
#'   medium, and large (e.g., c(0.1, 0.3, 0.5)). Corresponding values for other effect
#'   sizes are computed via simulation (default: \code{NULL}).
#' @inheritParams qhat
#' @inheritParams ES.summary.CI
#'
#' @return A list with components:
#'   \item{con}{The contrast matrix used.}
#'   \item{output}{List with one element per contrast column, each containing output from
#'     \code{\link{dep.ES.summary.CI}} for that linear contrast.}
#'
#' @details
#' This function generalizes \code{\link{dep.ES.summary.CI}} to handle linear contrasts
#' of dependent variables. For each contrast in \code{con}:
#' \enumerate{
#'   \item Computes the linear combination for each subject: \eqn{L_i = \sum_j c_j X_{ij}}
#'   \item Computes four effect size measures with bootstrap CIs using \code{\link{dep.ES.summary.CI}}:
#'     \itemize{
#'       \item \strong{AKP}: Robust standardized difference (analog of Cohen's d)
#'       \item \strong{QS (median)}: Quantile shift based on median of L
#'       \item \strong{QStr}: Quantile shift based on trimmed mean of L
#'       \item \strong{SIGN}: P(L < 0), probability the linear contrast is negative
#'     }
#' }
#'
#' If \code{con = NULL}, automatically generates all J(J-1)/2 pairwise contrasts.
#'
#' Special case: If \code{x} has 2 columns and \code{con = c(1, -1)}, produces the same
#' results as \code{dep.ES.summary.CI(x[,1], x[,2])}.
#'
#' Under no effect, the distribution of the linear contrast is symmetric about zero.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{dep.ES.summary.CI}}, \code{\link{rmlinES}}, \code{\link{DEP.PAIR.ES}}
#'
#' @export
#' @examples
#' # Three time points
#' set.seed(42)
#' time1 <- rnorm(20, 10)
#' time2 <- time1 + rnorm(20, 0.5, 1)
#' time3 <- time1 + rnorm(20, 1, 1)
#' x <- cbind(time1, time2, time3)
#'
#' # Linear trend contrast
#' con <- matrix(c(-1, 0, 1), ncol = 1)
#' result <- deplin.ES.summary.CI(x, con, nboot = 500)
#' result$output[[1]]
deplin.ES.summary.CI<-function(x,con=NULL,tr=.2,REL.MAG=NULL,SEED=TRUE,nboot=1000){
#
#  For J dependent variables,
#  compute four measures of effect size based on a linear contrast of the J variables specified by the argument
#  con
#
#  Generalizes dep.ES.summary.CI
#  Example:
#  If x is a matrix with two columns and con=c(1,-1), get the same results as dep.ES.summary.CI
#
#  By default, do all pairwise comparisons
#
#  Measures of effect size:
#
#  AKP: robust standardized difference similar to  Cohen's d
#  QS:  Quantile shift based on the median of the distribution of difference scores,
#  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
#  SIGN:  P(X<Y), probability that for a random pair, the first is less than the second.
#
#  If there is  no effect, distribution of the  linear contrast is  symmetric about zero.
#  Under normality, suppose a shift of .2, .5 and .8 standard deviation of the difference score
#  is considered small, medium and large. The corresponding values for QS and SIGN are printed.
#
#
#  REL.MAG: suppose it is decided that AKP values .1, .3 and .5 are viewed as small, medium and large under normality Then
# REL.MAG=c(.1,.3,.5)  the default,
# means that the function will compute the corresponding values for the measures of effect size used here.
#
CON=list()
if(is.list(x))x=matl(x)
x=elimna(x)
n=length(x)
J=ncol(x)
if(is.null(con))con=con.all.pairs(J)
for(k in 1:ncol(con)){
X=con[1,k]*x[,1]
for(j in 2:J)X=X+con[j,k]*x[,j]
output=dep.ES.summary.CI(X,tr=tr,REL.MAG=REL.MAG,SEED=SEED,nboot=nboot)
CON[[k]]=output
}
list(con=con, output=CON)
}


#' Pairwise Effect Sizes for Dependent Groups
#'
#' Computes effect size measures for all pairwise comparisons among J dependent
#' groups (repeated measures). For each pair, computes effect sizes based on the
#' difference scores.
#'
#' @param x Data in matrix, data frame, or list mode. Each column (or list element)
#'   represents one repeated measure/time point. Rows are matched observations.
#' @param tr Proportion of trimming for AKP effect size (default: 0.2).
#' @param CI Logical. If \code{TRUE}, computes bootstrap confidence intervals and
#'   p-values for each effect size (default: \code{FALSE}).
#' @param REL.MAG Optional vector of length 3 specifying AKP values considered small,
#'   medium, and large under normality (default: c(0.1, 0.3, 0.5)). Equivalent values
#'   for other effect sizes are computed via simulation.
#' @param SEED Logical. If \code{TRUE}, sets random seed for reproducibility (default: \code{TRUE}).
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{groups.compared}: Matrix with 2 columns showing which groups were
#'       compared in each row
#'     \item \code{Effect.size.estimates}: List with one element per pairwise comparison,
#'       each containing output from \code{dep.ES.summary} (if \code{CI=FALSE}) or
#'       \code{dep.ES.summary.CI} (if \code{CI=TRUE})
#'   }
#'
#' @details
#' For J dependent groups, computes all J(J-1)/2 pairwise comparisons. For each pair
#' of groups (i,j), computes effect sizes based on the distribution of difference
#' scores X_i - X_j using four measures:
#' \itemize{
#'   \item \strong{AKP}: Robust standardized difference (analog of Cohen's d)
#'   \item \strong{QS (median)}: Quantile shift based on median of differences
#'   \item \strong{QStr}: Quantile shift based on trimmed mean of differences
#'   \item \strong{SIGN}: P(X_i < X_j), probability first is less than second
#' }
#'
#' When \code{CI=FALSE}, provides point estimates with small/medium/large benchmarks.
#' When \code{CI=TRUE}, adds bootstrap confidence intervals and p-values.
#'
#' This function is useful for repeated measures designs where you want to assess
#' the magnitude of change between all time points or conditions.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{dep.ES.summary}}, \code{\link{dep.ES.summary.CI}},
#'   \code{\link{IND.PAIR.ES}} for independent groups
#'
#' @export
#' @examples
#' # Four time points (repeated measures)
#' set.seed(123)
#' n <- 20
#' time1 <- rnorm(n, mean=100, sd=15)
#' time2 <- time1 + rnorm(n, mean=3, sd=10)
#' time3 <- time2 + rnorm(n, mean=2, sd=10)
#' time4 <- time3 + rnorm(n, mean=1, sd=10)
#' x <- cbind(time1, time2, time3, time4)
#'
#' # Get effect sizes for all pairwise comparisons
#' result <- DEP.PAIR.ES(x)
#' result$groups.compared  # Shows which groups were compared
#' result$Effect.size.estimates[[1]]  # Time 1 vs Time 2
#'
#' # With confidence intervals
#' result_ci <- DEP.PAIR.ES(x, CI=TRUE)
DEP.PAIR.ES<-function(x,tr=.2,CI=FALSE,
REL.MAG=NULL,
SEED=TRUE){
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
J=length(x)
conmat=con.all.pairs(J)
ES=list()
Ncon=ncol(conmat)
ic=0

grp.com=matrix(NA,Ncon,2)
for(j in 1:Ncon){
id1=which(conmat[,j]==1)
id2=which(conmat[,j]==-1)
grp.com[j,]=c(id1,id2)
if(!CI)ES[[j]]=dep.ES.summary(x[[id1]],y=x[[id2]],tr=tr,REL.MAG=REL.MAG,SEED=SEED)
if(CI)ES[[j]]=dep.ES.summary.CI(x[[id1]],y=x[[id2]],tr=tr,REL.MAG=REL.MAG,SEED=SEED)
}
list(groups.compared=grp.com,Effect.size.estimates=ES)
}


#' Within-Within Design: Linear Contrast Effect Sizes
#'
#' For a within-within (fully repeated measures) factorial design, computes effect sizes
#' for linear contrasts with bootstrap confidence intervals.
#'
#' @param J Number of levels for Factor A (first within-subjects factor).
#' @param K Number of levels for Factor B (second within-subjects factor).
#' @param x Data in matrix or list mode. Expected to contain J*K repeated measures in
#'   row-major order.
#' @inheritParams ES.summary
#' @param REL.MAG Optional vector of length 3 for custom small/medium/large effect size
#'   calibration (default: \code{NULL}).
#' @inheritParams qhat
#' @inheritParams ES.summary.CI
#'
#' @return Effect size results with bootstrap confidence intervals for linear contrasts
#'   of the J*K repeated measures, as computed by \code{\link{deplin.ES.summary.CI}}.
#'
#' @details
#' This function computes effect sizes based on linear sums of the random variables in
#' a fully within-subjects J×K factorial design. It applies \code{\link{deplin.ES.summary.CI}}
#' to compute robust effect sizes with bootstrap confidence intervals.
#'
#' The effect sizes are based on linear contrasts of the J*K conditions, accounting for
#' the dependence structure in repeated measures data.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{wwwlin.es}}, \code{\link{deplin.ES.summary.CI}}, \code{\link{ww.es}}
#'
#' @export
#' @examples
#' # 2x2 within-within design
#' set.seed(42)
#' n <- 15
#' baseline <- rnorm(n)
#' x <- list(
#'   baseline + rnorm(n, 0, 0.5),    # Condition (1,1)
#'   baseline + rnorm(n, 0.5, 0.5),  # Condition (1,2)
#'   baseline + rnorm(n, 0, 0.5),    # Condition (2,1)
#'   baseline + rnorm(n, 1, 0.5)     # Condition (2,2)
#' )
#' result <- wwlin.es(J = 2, K = 2, x, nboot = 500)
wwlin.es<-function(J,K,x,tr = 0.2, REL.MAG = NULL, SEED = TRUE, nboot = 1000){
#
#  # For within-by-within
#
# Effect sizes based on linear sum of the random variables.
# Simplest case, compute effect sizes on x-y, difference scores
#
con=con2way(J,K)
A=deplin.ES.summary.CI(x,con=con$conA,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
B=deplin.ES.summary.CI(x,con=con$conB,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
AB=deplin.ES.summary.CI(x,con=con$conAB,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
list(Factor.A=A,Factor.B=B,Interactions=AB)
}


#' Within-Within-Within Design: Linear Contrast Effect Sizes
#'
#' For a within-within-within (fully repeated measures) three-way factorial design,
#' computes effect sizes for linear contrasts with bootstrap confidence intervals for
#' all main effects, two-way, and three-way interactions.
#'
#' @param J Number of levels for Factor A (first within-subjects factor).
#' @param K Number of levels for Factor B (second within-subjects factor).
#' @param L Number of levels for Factor C (third within-subjects factor).
#' @param x Data in matrix or list mode. Expected to contain J*K*L repeated measures.
#' @inheritParams ES.summary
#' @param REL.MAG Optional vector of length 3 for custom small/medium/large effect size
#'   calibration (default: \code{NULL}).
#' @inheritParams qhat
#' @inheritParams ES.summary.CI
#'
#' @return A list with components:
#'   \item{Factor.A, Factor.B, Factor.C}{Main effect sizes for each factor}
#'   \item{Inter.AB, Inte.AC, Inter.BC}{Two-way interaction effect sizes}
#'   \item{Inter.ABC}{Three-way interaction effect sizes}
#'
#' @details
#' This function handles three-way fully repeated measures designs by computing effect
#' sizes for:
#' \itemize{
#'   \item Three main effects (A, B, C)
#'   \item Three two-way interactions (A×B, A×C, B×C)
#'   \item One three-way interaction (A×B×C)
#' }
#'
#' Uses \code{\link{con3way}} to generate appropriate contrast matrices, then applies
#' \code{\link{deplin.ES.summary.CI}} for each set of contrasts.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{wwlin.es}}, \code{\link{deplin.ES.summary.CI}}, \code{\link{con3way}}
#'
#' @export
#' @examples
#' # 2x2x2 within design
#' set.seed(42)
#' n <- 12
#' baseline <- rnorm(n)
#' x <- lapply(1:8, function(i) baseline + rnorm(n, mean = i/4, sd = 0.5))
#' result <- wwwlin.es(J = 2, K = 2, L = 2, x, nboot = 500)
#' # Access main effects
#' result$Factor.A
#' # Access interactions
#' result$Inter.AB
wwwlin.es<-function(J,K,L,x,tr = 0.2, REL.MAG = NULL, SEED = TRUE, nboot = 1000){
#
# For within-by-within-by-within
#
# Effect sizes based on linear sum of the random variables.
# Simplest case, compute effect sizes based on x-y, difference scores
#
con=con3way(J,K,L)
A=deplin.ES.summary.CI(x,con=con$conA,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
B=deplin.ES.summary.CI(x,con=con$conB,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
C=deplin.ES.summary.CI(x,con=con$conC,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
AB=deplin.ES.summary.CI(x,con=con$conAB,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
AC=deplin.ES.summary.CI(x,con=con$conAC,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
BC=deplin.ES.summary.CI(x,con=con$conBC,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
ABC=deplin.ES.summary.CI(x,con=con$conABC,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
list(Factor.A=A,Factor.B=B,Factor.C=C,Inter.AB=AB,Inte.AC=AC,Inter.BC=BC,Inter.ABC=ABC)
}


#' Between-Within-Within Design: Factor A Effect Sizes
#'
#' For a between-within-within (mixed) three-way factorial design, computes effect sizes
#' for all pairwise comparisons of Factor A (between-subjects) at each combination of
#' Factors B and C (within-subjects).
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (first within-subjects factor).
#' @param L Number of levels for Factor C (second within-subjects factor).
#' @param x Data in matrix or list mode. Expected to contain J*K*L groups.
#' @param fun Effect size function to use (default: \code{\link{KMS.ci}} for KMS effect size
#'   with 20\% trimming).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param ... Additional arguments passed to \code{fun}.
#'
#' @return A matrix with one row per pairwise comparison of Factor A levels at each
#'   B×C combination. Columns include Factor A levels being compared, Factor B level,
#'   Factor C level, sample sizes, effect size estimate, and confidence interval limits.
#'
#' @details
#' This function is designed for three-way mixed designs with:
#' \itemize{
#'   \item Factor A: Between-subjects (J independent groups)
#'   \item Factor B: Within-subjects (K repeated measures)
#'   \item Factor C: Within-subjects (L repeated measures)
#' }
#'
#' For each of the K×L combinations of within-subjects factors, computes all pairwise
#' effect sizes comparing the J levels of the between-subjects factor.
#'
#' The default uses KMS robust effect size (\code{\link{KMS.ci}}) which is a robust
#' analog of Cohen's d.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{bwwA.es.sub}}, \code{\link{KMS.ci}}, \code{\link{bw.es.A}}
#'
#' @export
#' @examples
#' # 2x2x2 B-W-W design
#' set.seed(42)
#' n <- 15
#' # 2 between groups, 2x2 within factors
#' x <- list(
#'   rnorm(n, 10), rnorm(n, 11), rnorm(n, 12), rnorm(n, 13),  # Group 1
#'   rnorm(n, 10), rnorm(n, 12), rnorm(n, 14), rnorm(n, 16)   # Group 2
#' )
#' result <- bwwA.es(J = 2, K = 2, L = 2, x, nboot = 500)
#' result
bwwA.es<-function(J,K,L,x,fun=KMS.ci,nboot=1000,...){
#
#  For every two levels  of Factor A, compute  KMS shift effect size  assuming 20% trim
#  and do this for each
#  level of Factors B and C.
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JKL=J*K*L
KL=K*L
id=matrix(c(1:JKL),ncol=KL,byrow=TRUE)
con=con.all.pairs(J)
A=list()
N=(J^2-J)/2
for(i in 1:N){
w=which(con[,i]!=0)
isel=id[w,]
A[[i]]=bwwA.es.sub(2,K,L,x[isel],fun=fun,nboot=nboot,...)
}
list(A=A,con=con)
}


#' @keywords internal
bwwA.es.sub<-function(J,K,L,x,fun=QSci,nboot=1000,...){
#
#  Effect sizes for the between factor, computed
#  for each level of the within factors
#
if(J!=2)stop('Must have J=2')
JKL=J*K*L
KL=K*L
ic=0
LOW=1
UP=1+KL
id=matrix(c(1:JKL),ncol=KL,byrow=TRUE)
ES=matrix(NA,nrow=KL,ncol=7)
for(k in 1:K){
for(l in 1:L){
ic=ic+1
ES[ic,1]=k
ES[ic,2]=l
isel=id[,ic]
d=IND.PAIR.ES(x[id[,ic]],fun=fun,...)$effect.size
temp=c(d[[1]]$n1,d[[1]]$n2,d[[1]]$effect.size,d[[1]][4]$ci[1],d[[1]][4]$ci[2])
ES[ic,3:7]=temp
}}
dimnames(ES)=list(NULL,c('B.Level','C.Level','n1','n2','Effect.Size','ci.low','ci.up'))
ES
}


#' Identify the Cell with Largest (or Smallest) Probability
#'
#' For a multinomial distribution based on observed cell frequencies, determines whether
#' a decision can be made about which cell has the highest (or lowest) probability,
#' controlling for familywise error rate.
#'
#' @param x Numeric vector containing the cell frequencies.
#' @param alpha Nominal familywise error rate (default: 0.05).
#' @param LARGEST Logical. If \code{TRUE}, identifies the cell with the largest probability.
#'   If \code{FALSE}, identifies the cell with the smallest probability (default: \code{TRUE}).
#' @param method Method for computing confidence intervals. Options: \code{'AC'} (Agresti-Coull),
#'   or other methods supported by \code{\link{cell.com}} (default: \code{'AC'}).
#' @param p.crit Optional vector of critical p-values for comparisons. If \code{NULL},
#'   computed via \code{\link{best.cell.crit}} (default: \code{NULL}).
#' @param AUTO Logical. Passed to \code{\link{cell.com}} for automatic method selection
#'   (default: \code{FALSE}).
#' @param iter Number of iterations for computing critical p-values via simulation if
#'   \code{p.crit = NULL} (default: 2000).
#' @inheritParams qhat
#' @param pr Logical. If \code{TRUE}, prints informational message about critical p-values
#'   (default: \code{TRUE}).
#'
#' @return An S4 object of class \code{'BIN'} with slots:
#'   \item{Cell.with.largest.estimate (or smallest)}{The cell number with the largest (or smallest)
#'     estimated probability.}
#'   \item{Larger.than (or smaller.than)}{Vector of cell numbers for which the identified cell
#'     has been determined to have significantly larger (or smaller) probability.}
#'   \item{n}{Total sample size.}
#'   \item{output}{A matrix with one row per comparison. Columns:
#'     \itemize{
#'       \item Largest.Est (or Smallest.Est): Estimated probability for the identified cell
#'       \item CELL: Cell number being compared
#'       \item Est: Estimated probability for comparison cell
#'       \item Dif: Difference in probabilities
#'       \item ci.low, ci.up: Confidence interval for the difference
#'       \item p.value: Individual p-value
#'       \item p.crit: Critical p-value (adjusted for multiple comparisons)
#'     }
#'   }
#'
#' @details
#' This function addresses the question: "Can we determine which cell has the highest
#' (or lowest) probability while controlling the familywise error rate?"
#'
#' The procedure:
#' \enumerate{
#'   \item Identifies the cell with the largest (or smallest) observed frequency
#'   \item Compares this cell to all others using confidence intervals
#'   \item Uses critical p-values (via \code{\link{best.cell.crit}}) to control FWER
#'   \item Returns which cells are significantly different
#' }
#'
#' Missing values in \code{x} are removed before computation.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{best.cell.crit}}, \code{\link{cell.com}}
#'
#' @export
#' @examples
#' # Multinomial with 5 cells
#' # Cell frequencies suggest cell 3 has highest probability
#' x <- c(15, 18, 45, 20, 12)
#' result <- BEST.cell(x)
#' result
BEST.cell<-function(x,alpha=.05,LARGEST=TRUE,method='AC',p.crit=NULL,AUTO=FALSE,iter=2000,SEED=TRUE,pr=TRUE){
#
#  For a multinomial distribution, can a decision be made about
#  about which cell has the highest probability?
#
#  PV if specified, is a N by iter matrix of p-values that can be computed via best.cell.crit
#  N=number of cells
# x  Assumed to contain the  cell frequencies
#
if(pr)print('Confidence intervals are based on the critical p-values')
if(SEED)set.seed(2)
x=elimna(x)
n=sum(x)
NCELL=length(x)
NCm1=NCELL-1
xor=order(x,decreasing = LARGEST)
IND.pv=NA
ic=0
CI=matrix(NA,nrow=NCm1,ncol=2)
for(j in 2:NCELL){
ic=ic+1
IND.pv[ic]=cell.com.pv(x,xor[1],xor[j])
}
if(is.null(p.crit))p.crit=best.cell.crit(n,NCELL,LARGEST=LARGEST,iter=iter,AUTO=FALSE,SEED=SEED)
output=matrix(NA,nrow=NCm1,8)
output[,1]=rep(x[xor[1]]/n,NCm1)
output[,2]=xor[2:NCELL]
output[,3]=x[xor[2:NCELL]]/n
output[,4]=output[,1]-output[,3]
ic=0
for(j in 2:NCELL){
ic=ic+1
CI=cell.com(x,xor[1],xor[j],AUTO=AUTO,method=method,alpha=p.crit[ic])
output[ic,5:6]=CI$ci
}
output[,7]=IND.pv
output[,8]=p.crit
dimnames(output)=list(NULL,c('Largest.Est','CELL','Est','Dif','ci.low','ci.up','p.value','p.crit'))
flag=IND.pv<=p.crit
id=output[flag,2]
setClass('BIN',slots=c('Cell.with.largest.estimate','Larger.than','n','output'))
put=new('BIN',Cell.with.largest.estimate=xor[1],Larger.than=id,n=n,output=output)

if(!LARGEST){
dimnames(output)=list(NULL,c('Smallest.Est','CELL','Est','Dif','ci.low','ci.up','p.value','p.crit'))
setClass('BIN',slots=c('Cell.with.smallest.estimate','smaller.than','n','output'))
put=new('BIN',Cell.with.smallest.estimate=xor[1],smaller.than=id,n=n,output=output)
}
put
}


#' KMS Robust Effect Size Using M-Estimator
#'
#' Computes a robust effect size using the KMS (Kulinskaya-Morgenthaler-Staudte) method
#' based on an M-estimator of location and percentage bend variance.
#'
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group.
#'
#' @return A list with components:
#'   \item{effect.size}{The KMS effect size estimate}
#'   \item{Cohen.d.equiv}{Equivalent Cohen's d value (= 2 * effect.size)}
#'
#' @details
#' This function implements a robust generalization of Cohen's d using:
#' \itemize{
#'   \item One-step M-estimator (\code{\link{onestep}}) for location
#'   \item Percentage bend variance (\code{\link{pbvar}}) for scale
#' }
#'
#' The relationship to traditional Cohen's d values is approximately:
#' \itemize{
#'   \item Cohen d = 0.2 (small) corresponds to KMS \eqn{\approx} 0.1
#'   \item Cohen d = 0.5 (medium) corresponds to KMS \eqn{\approx} 0.25
#'   \item Cohen d = 0.8 (large) corresponds to KMS \eqn{\approx} 0.4
#' }
#'
#' Missing values are removed before computation.
#'
#' @references
#' Kulinskaya, E., Morgenthaler, S. & Staudte, R. (2008). Meta Analysis: A guide to
#' calibrating and combining statistical evidence. John Wiley & Sons. p. 177.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @seealso \code{\link{onestep}}, \code{\link{pbvar}}, \code{\link{akp.effect}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 0)
#' y <- rnorm(30, mean = 0.5)
#' KMS.ES.M(x, y)
#'
#' # With outliers - robust to contamination
#' x2 <- c(x, 10, 15)
#' y2 <- c(y, -10)
#' KMS.ES.M(x2, y2)
KMS.ES.M<-function(x,y){
#
# Computes the robust effect size using a simple generalization of the method in
# Kulinskaya, E., Morgenthaler, S. & Staudte, R. (2008).
# Meta Analysis: A guide to calibrating and combining statistical evidence  p. 177
# based	on an M-estimator and percentage bend variance
#Cohen d=.2, .5 .8 correspond to .1, .25 and .4') (KMS p. 180)

x<-elimna(x)
y<-elimna(y)
n1<-length(x)
n2<-length(y)
N=n1+n2
q=n1/N
s1sq=pbvar(x)
s2sq=pbvar(y)
t1=onestep(x)
t2=onestep(y)
top=q*s1sq+(1-q)*s2sq
bot=q*(1-q)
sigsq=top/bot  #  Quantity in brackets KMS p. 176 eq 21.1
varrho=s2sq/s1sq
d1=(t1-t2)/sqrt(sigsq)
list(effect.size=d1,Cohen.d.equiv=2*d1)
}


#' @keywords internal
ES.summary.sub<-function(x,n1,n2){
id1=c(1:n1)
n1p=n1+1
N=n1+n2
id2=c(n1p:N)
a=ES.summary.CI(x[id1],x[id2],SEED=F)[,8]
}

