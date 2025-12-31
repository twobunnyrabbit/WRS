# ==============================================================================
# Two-Sample Comparison Functions
# ==============================================================================
#
# This module contains functions for comparing two groups (independent or paired).
#
# Main functions:
#   - yuen(), yuend() - Trimmed mean comparisons (independent/dependent)
#   - pb2gen() - Percentile bootstrap for two groups
#   - wmw*() - Wilcoxon-Mann-Whitney and related tests
#   - cid*() - Cliff's analog and related methods
#   - qcomhd*() - Quantile comparisons using Harrell-Davis
#   - trimpb*() - Bootstrap trimmed mean comparisons
#   - two*() - Various two-sample comparison methods
#
# Extracted: 2025-12-30
# Functions: 103
# ==============================================================================

#' Compare Multiple Quantiles Between Two Groups (Multicore Version)
#'
#' Compares quantiles between two independent groups using the Harrell-Davis
#' estimator with percentile bootstrap. Uses Hochberg's method to control
#' the familywise error rate across multiple quantile comparisons.
#'
#' @inheritParams common-params
#' @param est Estimator function for quantiles (default: `hd`, the Harrell-Davis estimator).
#' @param q Vector of quantiles to compare (default: c(.1, .25, .5, .75, .9)).
#' @param xlab Label for x-axis in plot (default: "Group 1").
#' @param ylab Label for y-axis in plot (default: "Est.1-Est.2").
#' @param ADJ.CI Logical. If `TRUE`, adjusts confidence intervals based on the
#'   significance level used by the corresponding test statistic (default: `TRUE`).
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item `q`: The quantile being compared
#'     \item `n1`, `n2`: Sample sizes for groups 1 and 2
#'     \item `est.1`, `est.2`: Quantile estimates for each group
#'     \item `est.1_minus_est.2`: Difference between group estimates
#'     \item `ci.low`, `ci.up`: Confidence interval bounds for the difference
#'     \item `p_crit`: Critical p-value based on Hochberg's method
#'     \item `p-value`: P-value for the comparison
#'     \item `signif`: "YES" if significant, "NO" otherwise
#'   }
#'
#' @details
#' This function uses the percentile bootstrap method via `pb2genMC` to compare
#' quantiles. Tied values are allowed. When comparing lower or upper quartiles,
#' both power and Type I error rates compare favorably to other methods.
#'
#' The function uses Hochberg's sequentially rejective procedure to control the
#' familywise error rate across all quantile comparisons. If `ADJ.CI = TRUE`,
#' confidence intervals are adjusted to maintain consistency with the test results.
#'
#' A plot is created showing the point estimates and confidence intervals for
#' each quantile difference.
#'
#' @note This function prints a message suggesting that `qcomhd` with `MC=TRUE`
#'   can also be used for the same analysis.
#'
#' @seealso \code{\link{qcomhd}}, \code{\link{pb2genMC}}, \code{\link{hd}}
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(50, mean = 0, sd = 1)
#' y <- rnorm(50, mean = 0.5, sd = 1.2)
#'
#' # Compare multiple quantiles
#' result <- qcomhdMC(x, y, nboot = 1000)
#' print(result)
#'
#' # Compare specific quantiles only
#' qcomhdMC(x, y, q = c(.25, .5, .75), nboot = 1000)
#' }
qcomhdMC<-function(x,y,est=hd,q=c(.1,.25,.5,.75,.9),nboot=4000,plotit=TRUE,SEED=TRUE,xlab="Group 1",ylab="Est.1-Est.2",alpha=.05,ADJ.CI=TRUE){
#
# Compare quantiles using pb2gen
# via hd estimator. Tied values are allowed.
#
#  ADJ.CI=TRUE means that the confidence intervals are adjusted based on the level used by the corresponding
#  test statistic. If a test is performed with at the .05/3 level, for example, the confidence returned has
#  1-.05/3 probability coverage.
#
# When comparing lower or upper quartiles, both power and the probability of Type I error
# compare well to other methods that have been derived.
# q: can be used to specify the quantiles to be compared
# q defaults to comparing the .1,.25,.5,.75, and .9 quantiles
#
#   Function returns p-values and critical p-values based on Hochberg's method.
#
if(SEED)set.seed(2)
print('Can also use the function qcomhd with the argument MC=TRUE')
pv=NULL
output=matrix(NA,nrow=length(q),ncol=10)
dimnames(output)<-list(NULL,c("q","n1","n2","est.1","est.2","est.1_minus_est.2","ci.low","ci.up","p_crit","p-value"))
for(i in 1:length(q)){
output[i,1]=q[i]
output[i,2]=length(elimna(x))
output[i,3]=length(elimna(y))
output[i,4]=hd(x,q=q[i])
output[i,5]=hd(y,q=q[i])
output[i,6]=output[i,4]-output[i,5]
temp=qcom.sub(x,y,nboot=nboot,q=q[i],SEED=FALSE,alpha=alpha)
output[i,7]=temp$ci[1]
output[i,8]=temp$ci[2]
output[i,10]=temp$p.value
}
temp=order(output[,10],decreasing=TRUE)
zvec=alpha/c(1:length(q))
output[temp,9]=zvec
if(ADJ.CI){
for(i in 1:length(q)){
temp=pb2genMC(x,y,nboot=nboot,est=est,q=q[i],SEED=FALSE,alpha=output[i,9],pr=FALSE)
output[i,7]=temp$ci[1]
output[i,8]=temp$ci[2]
output[i,10]=temp$p.value
}
temp=order(output[,10],decreasing=TRUE)
}
output <- data.frame(output)
output$signif=rep("YES",nrow(output))
for(i in 1:nrow(output)){
if(output[temp[i],10]>output[temp[i],9])output$signif[temp[i]]="NO"
#if(output[temp[i],10]<=output[temp[i],9])break
}
if(plotit){
xax=rep(output[,4],3)
yax=c(output[,6],output[,7],output[,8])
plot(xax,yax,xlab=xlab,ylab=ylab,type="n")
points(output[,4],output[,6],pch="*")
lines(output[,4],output[,6])
points(output[,4],output[,7],pch="+")
points(output[,4],output[,8],pch="+")
}
output
}

#' Multivariate Two-Sample Behrens-Fisher Test Using Bootstrap
#'
#' Compares trimmed means between two independent groups for multiple measures
#' (multivariate Behrens-Fisher problem). Uses the percentile-t bootstrap
#' method with simultaneous confidence intervals controlling familywise error rate.
#'
#' @inheritParams common-params
#'
#' @details
#' For each of two independent groups with p measures per subject, this function
#' compares the trimmed means of the first measure, the second measure, and so on,
#' resulting in p comparisons between the two groups.
#'
#' The percentile-t bootstrap method is used to compute simultaneous confidence
#' intervals with family-wise error control. By default, 20% trimming is used
#' with 599 bootstrap samples.
#'
#' Data can be provided as n-by-p matrices or in list mode. Missing values are
#' not allowed.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `psihat`: Matrix with columns: `con.num` (contrast number),
#'       `psihat` (difference in trimmed means), `ci.lower`, `ci.upper`
#'     \item `teststat`: Matrix with test statistics and standard errors
#'     \item `critical.value`: Bootstrap critical value for simultaneous inference
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' # Create multivariate data for two groups
#' set.seed(123)
#' g1 <- matrix(rnorm(50 * 3), ncol = 3)  # 50 subjects, 3 measures
#' g2 <- matrix(rnorm(60 * 3, mean = 0.3), ncol = 3)
#'
#' # Compare trimmed means for all measures
#' result <- twomanbt(g1, g2)
#' print(result$psihat)
#' }
twomanbt<-function(x,y,tr=.2,alpha=.05,nboot=599){
#
#   Two-sample Behrens-Fisher problem.
#
#   For each of two independent groups,
#   have p measures for each subject. The goal is to compare the
#   trimmed means of the first measure, the trimmed means for the second
#   and so on.   So there are a total of p comparisons between the two
#   groups, one for each measure.
#
#   The percentile t bootstrap method is used to
#   compute a .95 confidence interval.
#
#   By default, 20% trimming is used with B=599 bootstrap samples.
#
#   x contains the data for the first group; it
#    can be an n by J matrix or it can have list mode.
#   y contains the data for the second group.
#
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(!is.list(y) && !is.matrix(y))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
# put the data in an n by p matrix
matx<-matrix(0,length(x[[1]]),length(x))
for (j in 1:length(x))matx[,j]<-x[[j]]
}
if(is.list(y)){
# put the data in an n by p matrix
maty<-matrix(0,length(y[[1]]),length(y))
for (j in 1:length(y))maty[,j]<-y[[j]]
}
if(is.matrix(x)){
matx<-x
}
if(is.matrix(y)){
maty<-y
}
if(ncol(matx)!=ncol(maty))stop("The number of variables for group one is not equal to the number for group 2")
if(sum(is.na(mat)>=1))stop("Missing values are not allowed.")
J<-ncol(mat)
connum<-ncol(matx)
bvec<-matrix(0,connum,nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
xcen<-matrix(0,nrow(matx),ncol(matx))
ycen<-matrix(0,nrow(maty),ncol(maty))
for (j in 1:connum)xcen[,j]<-matx[,j]-mean(matx[,j],tr) #Center data
for (j in 1:connum)ycen[,j]<-maty[,j]-mean(maty[,j],tr) #Center data
print("Taking bootstrap samples. Please wait.")
bootx<-sample(nrow(matx),size=nrow(matx)*nboot,replace=TRUE)
booty<-sample(nrow(maty),size=nrow(maty)*nboot,replace=TRUE)
matval<-matrix(0,nrow=nboot,ncol=connum)
for (j in 1:connum){
datax<-matrix(xcen[bootx,j],ncol=nrow(matx))
datay<-matrix(ycen[booty,j],ncol=nrow(maty))
paste("Working on variable", j)
top<- apply(datax, 1., mean, tr) - apply(datay, 1., mean, tr)
botx <- apply(datax, 1., trimse, tr)
boty <- apply(datay, 1., trimse, tr)
matval[,j]<-abs(top)/sqrt(botx^2. + boty^2.)
}
bvec<-apply(matval,1,max)
icrit<-round((1-alpha)*nboot)
bvec<-sort(bvec)
crit<-bvec[icrit]
psihat<-matrix(0,ncol=4,nrow=connum)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol=3,nrow=connum)
dimnames(test)<-list(NULL,c("con.num","test","se"))
for(j in 1:ncol(matx)){
temp<-yuen(matx[,j],maty[,j],tr=tr)
test[j,1]<-j
test[j,2]<-abs(temp$test)
test[j,3]<-temp$se
psihat[j,1]<-j
psihat[j,2]<-mean(matx[,j],tr)-mean(maty[,j])
psihat[j,3]<-mean(matx[,j],tr)-mean(maty[,j])-crit*temp$se
psihat[j,4]<-mean(matx[,j],tr)-mean(maty[,j])+crit*temp$se
}
list(psihat=psihat,teststat=test,critical.value=crit)
}

#' Confidence Interval for Difference Between Two Binomial Proportions
#'
#' Computes a confidence interval for the difference between two independent
#' binomial proportions (p1 - p2) using Beal's method.
#'
#' @param r1 Number of successes in group 1 (default: sum of non-NA values in `x`).
#' @param n1 Sample size for group 1 (default: length of non-NA values in `x`).
#' @param r2 Number of successes in group 2 (default: sum of non-NA values in `y`).
#' @param n2 Sample size for group 2 (default: length of non-NA values in `y`).
#' @inheritParams common-params
#'
#' @details
#' If vectors `x` and `y` are provided, the function treats them as binary data
#' (0 = failure, 1 = success) and computes `r1`, `n1`, `r2`, and `n2` automatically.
#' Otherwise, these values must be provided explicitly.
#'
#' Beal's method provides accurate confidence intervals for the difference between
#' two binomial proportions without relying on large-sample approximations.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `ci`: Two-element vector with lower and upper confidence bounds
#'     \item `p1`: Estimated proportion for group 1 (r1/n1)
#'     \item `p2`: Estimated proportion for group 2 (r2/n2)
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' # Using summary statistics
#' twobici(r1 = 15, n1 = 50, r2 = 25, n2 = 60)
#'
#' # Using binary data vectors
#' x <- rbinom(50, 1, 0.3)
#' y <- rbinom(60, 1, 0.4)
#' twobici(x = x, y = y)
#' }
twobici<-function(r1=sum(elimna(x)),n1=length(elimna(x)),r2=sum(elimna(y)),n2=length(elimna(y)),
x=NA,y=NA,alpha=.05){
#
# Compute confidence interval for p1-p2,
# the difference between probabilities of
# success for a two binomials using Beal's method.
#
# r is number of successes
# n is sample size
# if x contains data, r1 is taken to be the
# number of 1s in x and n1 is length(x)
#
if(length(r1)>1)stop("r1 must be a single number, not a vector")
if(length(n1)>1)stop("n1 must be a single number, not a vector")
if(length(r2)>1)stop("r2 must be a single number, not a vector")
if(!is.na(sum(r1)) || !is.na(sum(n1)) || !is.na(sum(r2)) || !is.na(sum(n2))){
if(r1<0 || n1<0)stop("Both r1 and n1 must be greater than 0")
if(r1 > n1)stop("r1 can't be greater than n1")
if(r2<0 || n2<0)stop("Both r2 and n2 must be greater than 0")
if(r2 > n2)stop("r2 can't be greater than n2")
}
if(!is.na(sum(x))){
r1<-sum(x)
n1<-length(x)
}
if(!is.na(sum(y))){
r2<-sum(y)
n2<-length(y)
}
a<-(r1/n1)+(r2/n2)
b<-(r1/n1)-(r2/n2)
u<-.25*((1/n1)+(1/n2))
v<-.25*((1/n1)-(1/n2))
V<-u*((2-a)*a-b^2)+2*v*(1-a)*b
crit<-qchisq(1-alpha/2,1)
A<-sqrt(crit*(V+crit*u^2*(2-a)*a+crit*v^2*(1-a)^2))
B<-(b+crit*v*(1-a))/(1+crit*u)
ci<-NA
ci[1]<-B-A/(1+crit*u)
ci[2]<-B+A/(1+crit*u)
p1<-r1/n1
p2<-r2/n2
list(ci=ci,p1=p1,p2=p2)
}

#' Percentile Bootstrap CI for Difference in Trimmed Means
#'
#' Computes a percentile bootstrap confidence interval for the difference
#' between two trimmed means. Independent groups are assumed.
#'
#' @inheritParams common-params
#' @param WIN Logical. If `TRUE`, winsorizes data before bootstrapping (default: `FALSE`).
#' @param win Amount of Winsorizing to apply when `WIN = TRUE` (default: 0.1).
#' @param op Plot type parameter for `g2plot` when `plotit = TRUE` (default: 4).
#'
#' @details
#' This function computes a bootstrap confidence interval for the difference
#' between trimmed means of two independent groups. Missing values are automatically removed.
#'
#' When `WIN = TRUE`, data are winsorized before bootstrapping. The amount of trimming
#' should be at least 0.2 when winsorizing. Sample sizes less than 15 may result in
#' poor control of Type I error when winsorizing.
#'
#' If `plotit = TRUE`, creates a plot comparing the bootstrap distributions using `g2plot`.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `Est1`: Trimmed mean for group 1
#'     \item `Est2`: Trimmed mean for group 2
#'     \item `p.value`: P-value for the test
#'     \item `ci`: Two-element vector with confidence interval bounds
#'     \item `est.dif`: Estimated difference (Est1 - Est2)
#'   }
#'
#' @seealso \code{\link{trimpb}}, \code{\link{pb2gen}}
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(60, mean = 0.5)
#'
#' # Basic usage
#' trimpb2(x, y)
#'
#' # With winsorizing
#' trimpb2(x, y, WIN = TRUE, win = 0.1)
#'
#' # With plotting
#' trimpb2(x, y, plotit = TRUE)
#' }
trimpb2<-function(x,y,tr=.2,alpha=.05,nboot=2000,WIN=FALSE,win=.1,plotit=FALSE,op=4,
SEED=TRUE){
#
#   Compute a 1-alpha confidence interval for
#   the difference between two 20% trimmed means.
#   Independent groups are assumed.
#
#   The default number of bootstrap samples is nboot=2000
#
#   tr is the amount of trimming
#
#   win is the amount of Winsorizing before bootstrapping
#   when WIN=T.
#
#   Missing values are automatically removed.
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
if(WIN){
if(win>tr)stop("Cannot Winsorize more than you trim")
if(tr < .2){print("When Winsorizing, the amount of trimming")
print("should be at least .2")
}
if(min(c(length(x),length(y))) < 15){
print ("Warning: Winsorizing with sample sizes less than 15")
print("can result in poor control over the probability of a Type I error")
}
x<-winval(x,win)
y<-winval(y,win)
}
xx<-list()
xx[[1]]<-x
xx[[2]]<-y
e1=mean(xx[[1]],tr=tr)
e2=mean(xx[[2]],tr=tr)
#est.dif<-tmean(xx[[1]],tr=tr)-tmean(xx[[2]],tr=tr)
est.dif=e1-e2
crit<-alpha/2
temp<-round(crit*nboot)
icl<-temp+1
icu<-nboot-temp
bvec<-matrix(NA,nrow=2,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:2){
data<-matrix(sample(xx[[j]],size=length(xx[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,mean,tr) # Bootstrapped trimmed means for jth group
}
top<-bvec[1,]-bvec[2,]
test<-sum(top<0)/nboot+.5*sum(top==0)/nboot
if(test > .5)test<-1-test
top<-sort(top)
ci<-NA
ci[1]<-top[icl]
ci[2]<-top[icu]
if(plotit)g2plot(bvec[1,],bvec[2,],op=op)
list(Est1=e1,Est2=e2,p.value=2*test,ci=ci,est.dif=est.dif)
}

#' Bootstrap CI for Difference Between Two Regression Slopes
#'
#' Computes a 95% confidence interval for the difference between two OLS
#' regression slopes corresponding to two independent groups.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#'
#' @details
#' Uses an adjusted percentile bootstrap method that gives good results when
#' the error term is heteroscedastic. The function uses 599 bootstrap samples.
#'
#' WARNING: If the number of bootstrap samples is altered, it is unknown how to
#' adjust the confidence interval when n1 + n2 < 250. The confidence interval
#' bounds are automatically adjusted based on the combined sample size.
#'
#' Missing values are automatically removed. Only one covariate is allowed.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `b1`: OLS slope estimate for group 1
#'     \item `b2`: OLS slope estimate for group 2
#'     \item `ci`: 95% confidence interval for the difference (b1 - b2)
#'   }
#'
#' @seealso \code{\link{tworegwb}}, \code{\link{reg2ci}}
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' x1 <- rnorm(50)
#' y1 <- 2 + 3 * x1 + rnorm(50)
#' x2 <- rnorm(60)
#' y2 <- 2 + 2.5 * x2 + rnorm(60)
#'
#' result <- twolsreg(x1, y1, x2, y2)
#' print(result)
#' }
twolsreg<-function(x1,y1,x2,y2){
#
#   Compute a .95 confidence interval for
#   the difference between two regression slopes,
#   estimated via least squares and
#    corresponding to two independent groups.
#
#   This function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
#
#   WARNING: If the number of boostrap samples is altered, it is
#   unknown how to adjust the confidence interval when n1+n2 < 250.
#
nboot<-599  #Number of bootstrap samples
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples; please wait")
xy=elimna(cbind(x1,y1))
if(ncol(xy)>2)stop("This function only allows one covariate")
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data1,1,twolsregsub,x1,y1) # A 1 by nboot matrix.
data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
bvec2<-apply(data2,1,twolsregsub,x2,y2) # A 1 by nboot matrix.
bvec<-bvec1-bvec2
ilow<-15
ihi<-584
if(length(y1)+length(y2) < 250){
ilow<-14
ihi<-585
}
if(length(y1)+length(y2) < 180){
ilow<-11
ihi<-588
}
if(length(y1)+length(y2) < 80){
ilow<-8
ihi<-592
}
if(length(y1)+length(y2) < 40){
ilow<-7
ihi<-593
}
bsort<-sort(bvec)
b1<-lsfit(x1,y1)$coef[2]
b2<-lsfit(x2,y2)$coef[2]
ci<-c(bsort[ilow],bsort[ihi])
list(b1=b1,b2=b2,ci=ci)
}

#' Mann-Whitney U Test with Adjusted P-Value
#'
#' Performs the Mann-Whitney (Wilcoxon rank-sum) test and returns both the
#' usual p-value and an adjusted p-value using the Hodges-Ramsey-Wechsler method.
#'
#' @inheritParams common-params
#'
#' @details
#' This function computes the Mann-Whitney U test statistic and returns two p-values:
#' the usual asymptotic p-value and an adjusted p-value using the Hodges, Ramsey and
#' Wechsler (1990) method that provides better finite-sample performance.
#'
#' Missing values are automatically removed before analysis.
#'
#' The function also returns p.hat, the estimated probability that a randomly
#' selected value from y exceeds a randomly selected value from x.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.value`: Standard asymptotic p-value for Mann-Whitney test
#'     \item `adj.p.value`: Adjusted p-value using Hodges-Ramsey-Wechsler method
#'     \item `p.hat`: Estimate of P(Y > X), computed as U/(n*m)
#'   }
#'
#' @references Hodges, J. L., Ramsey, P. H., & Wechsler, S. (1990). Improved
#'   significance probabilities of the Wilcoxon test. Journal of Educational
#'   Statistics, 15(3), 249-265.
#'
#' @seealso \code{\link{wmwloc}}, \code{\link{wmwpb}}, \code{\link{wmw.RZR}}
#'
#' @export
#' @examples
#' \dontrun{
#' x <- rnorm(30)
#' y <- rnorm(40, mean = 0.5)
#'
#' result <- wmw(x, y)
#' print(result)
#' }
wmw<-function(x,y){
#
# Do Mann-Whitney test
# Return the usual p-value followed by adjusted
# p-value using Hodges, Ramsey and Wechsler (1990) method
# (See Wilcox, 2003, p. 559.)
#
x=elimna(x)
y=elimna(y)
m<-length(x)
n<-length(y)
com<-rank(c(x,y))
xp1<-length(x)+1
x<-com[1:length(x)]
y<-com[xp1:length(com)]
u<-sum(y)-n*(n+1)/2
sigsq<-m*n*(n+m+1)/12
yv<-(u+.5-m*n/2)/sqrt(sigsq)
kv<-20*m*n*(m+n+1)/(m^2+n^2+n*m+m+n)
S<-yv^2
T1<-S-3
T2<-(155*S^2-416*S-195)/42
cv<-1+T1/kv+T2/kv^2
sighrw<-2*(1-pnorm(abs(cv*yv)))
z<-(u-(.5*m*n))/sqrt(sigsq)
sig<-2*(1-pnorm(abs(z)))
list(p.value=sig,adj.p.value=sighrw,p.hat=u/(n*m))
}

#' Bootstrap CI for Difference Between Two Pearson Correlations
#'
#' Computes a 95% confidence interval for the difference between two Pearson
#' correlations corresponding to two independent groups using adjusted percentile
#' bootstrap.
#'
#' @param x1 First variable for group 1.
#' @param y1 Second variable for group 1.
#' @param x2 First variable for group 2.
#' @param y2 Second variable for group 2.
#' @inheritParams common-params
#'
#' @details
#' Uses an adjusted percentile bootstrap method (599 samples) that gives good
#' results when the error term is heteroscedastic. The confidence interval bounds
#' are automatically adjusted based on combined sample size.
#'
#' WARNING: If the number of bootstrap samples is altered, adjustment of CI
#' bounds for n1 + n2 < 250 is unknown.
#'
#' Missing values are automatically removed.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `r1`: Pearson correlation for group 1
#'     \item `r2`: Pearson correlation for group 2
#'     \item `ci`: 95% confidence interval for r1 - r2
#'   }
#'
#' @seealso \code{\link{tworhobt}}, \code{\link{twohc4cor}}
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' x1 <- rnorm(50)
#' y1 <- 0.5 * x1 + rnorm(50)
#' x2 <- rnorm(60)
#' y2 <- 0.3 * x2 + rnorm(60)
#'
#' result <- twopcor(x1, y1, x2, y2)
#' print(result)
#' }
twopcor<-function(x1,y1,x2,y2,SEED=TRUE){
#
#   Compute a .95 confidence interval for
#   the difference between two Pearson
#   correlations corresponding to two independent
#   goups.
#
#   This function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
#
#   WARNING: If the number of bootstrap samples is altered, it is
#   unknown how to adjust the confidence interval when n1+n2 < 250.
#
nboot<-599  #Number of bootstrap samples
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
X<-elimna(cbind(x1,y1))
x1<-X[,1]
y1<-X[,2]
X<-elimna(cbind(x2,y2))
x2<-X[,1]
y2<-X[,2]
data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data1,1,pcorbsub,x1,y1) # A 1 by nboot matrix.
data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
bvec2<-apply(data2,1,pcorbsub,x2,y2) # A 1 by nboot matrix.
bvec<-bvec1-bvec2
ilow<-15
ihi<-584
if(length(y1)+length(y2) < 250){
ilow<-14
ihi<-585
}
if(length(y1)+length(y2) < 180){
ilow<-11
ihi<-588
}
if(length(y1)+length(y2) < 80){
ilow<-8
ihi<-592
}
if(length(y1)+length(y2) < 40){
ilow<-7
ihi<-593
}
bsort<-sort(bvec)
r1<-cor(x1,y1)
r2<-cor(x2,y2)
ci<-c(bsort[ilow],bsort[ihi])
list(r1=r1,r2=r2,ci=ci)
}

#' Bootstrap-t Test for Comparing Two Independent Correlations
#'
#' Compares two independent correlations using a bootstrap-t method with the
#' HC4 heteroscedasticity-consistent estimator.
#'
#' @param X1 First variable for group 1.
#' @param Y1 Second variable for group 1.
#' @param X2 First variable for group 2.
#' @param Y2 Second variable for group 2.
#' @inheritParams common-params
#'
#' @details
#' Uses a bootstrap-t approach with the HC4 estimator to compare correlations
#' between two independent groups. This method provides better control of Type I
#' error rates than methods based on Fisher's z-transformation, especially with
#' non-normal data.
#'
#' The function standardizes the variables before computing correlations and
#' uses the HC4 method for robust standard error estimation.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `test`: Test statistic
#'     \item `crit.val`: Bootstrap critical values
#'     \item `p.value`: P-value for the test
#'   }
#'
#' @seealso \code{\link{twopcor}}, \code{\link{twohc4cor}}
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' X1 <- rnorm(50)
#' Y1 <- 0.5 * X1 + rnorm(50)
#' X2 <- rnorm(60)
#' Y2 <- 0.3 * X2 + rnorm(60)
#'
#' result <- tworhobt(X1, Y1, X2, Y2)
#' print(result)
#' }
tworhobt<-function(X1,Y1,X2,Y2,alpha=.05,nboot=499,SEED=TRUE){
#
# compare two independent correlations using a bootstrap-t method in conjunction with the HC4 estimator
#
if(SEED)set.seed(2)
r1=cor(X1,Y1)
r2=cor(X2,Y2)
n1=length(X1)
n2=length(X2)
v=NA
Nboot=nboot+1
for(i in 1:Nboot){
if(i<=nboot){
id1=sample(n1,n1,replace=TRUE)
id2=sample(n2,n2,replace=TRUE)
}
if(i==Nboot){
id1=c(1:n1)
id2=c(1:n2)
}
x1=X1[id1]
y1=Y1[id1]
x2=X2[id2]
y2=Y2[id2]
x1=(x1-mean(x1))/sd(x1)
y1=(y1-mean(y1))/sd(y1)
x2=(x2-mean(x2))/sd(x2)
y2=(y2-mean(y2))/sd(y2)
temp1=olshc4(x1,y1)
temp2=olshc4(x2,y2)
if(i<=nboot)v[i]=(temp1$ci[2,2]-r1-temp2$ci[2,2]+r2)/sqrt(temp1$ci[2,6]^2+temp2$ci[2,6]^2)
if(i==Nboot)v[i]=(temp1$ci[2,2]-temp2$ci[2,2])/sqrt(temp1$ci[2,6]^2+temp2$ci[2,6]^2)
}
ibot<-round(alpha*nboot/2)
itop<-nboot-ibot+1
ibot=ibot+1    #adjusted so that p-value and confidence interval give consistent results.
vs=sort(v[1:nboot])
crit=c(vs[ibot],vs[itop])
test=v[Nboot]
if(test<0)G=mean(test>v[1:nboot])
if(test>=0)G=mean(test<v[1:nboot])
pv=2*G
if(pv>1)pv=1
if(pv<0)pv=0
list(test=test,crit.val=crit,p.value=pv)
}

#' Test for Equality of Two Binomial Proportions
#'
#' Tests the hypothesis that two independent binomial distributions have equal
#' probability of success using the Storer-Kim method.
#'
#' @param r1 Number of successes in group 1.
#' @param n1 Sample size for group 1.
#' @param r2 Number of successes in group 2.
#' @param n2 Sample size for group 2.
#' @inheritParams common-params
#'
#' @details
#' The Storer-Kim method provides an exact test for comparing two binomial
#' proportions without relying on large-sample approximations. It has better
#' control of Type I error rates than asymptotic methods, especially with
#' small sample sizes.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.value`: P-value for the test
#'     \item `p1`: Estimated proportion for group 1 (r1/n1)
#'     \item `p2`: Estimated proportion for group 2 (r2/n2)
#'     \item `est.dif`: Estimated difference (p1 - p2)
#'   }
#'
#' @seealso \code{\link{twobici}}, \code{\link{twobicipv}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test if success rates differ
#' twobinom(r1 = 15, n1 = 50, r2 = 25, n2 = 60)
#'
#' # Using binary data vectors
#' x <- rbinom(50, 1, 0.3)
#' y <- rbinom(60, 1, 0.4)
#' twobinom(x = x, y = y)
#' }
twobinom<-function(r1=sum(elimna(x)),n1=length(elimna(x)),r2=sum(elimna(y)),n2=length(elimna(y)),x=NA,y=NA,alpha=.05){
#
# Test the hypothesis that two independent binomials have equal
# probability of success using the Storer--Kim method.
#
# r1=number of successes in group 1
# n1=number of observations in group 1
#
n1p<-n1+1
n2p<-n2+1
n1m<-n1-1
n2m<-n2-1
chk<-abs(r1/n1-r2/n2)
x<-c(0:n1)/n1
y<-c(0:n2)/n2
phat<-(r1+r2)/(n1+n2)
m1<-outer(x,y,"-")
m2<-matrix(1,n1p,n2p)
flag<-(abs(m1)>=chk)
m3<-m2*flag
b1<-1
b2<-1
xv<-c(1:n1)
yv<-c(1:n2)
xv1<-n1-xv+1
yv1<-n2-yv+1
dis1<-c(1,pbeta(phat,xv,xv1))
dis2<-c(1,pbeta(phat,yv,yv1))
pd1<-NA
pd2<-NA
for(i in 1:n1)pd1[i]<-dis1[i]-dis1[i+1]
for(i in 1:n2)pd2[i]<-dis2[i]-dis2[i+1]
pd1[n1p]<-phat^n1
pd2[n2p]<-phat^n2
m4<-outer(pd1,pd2,"*")
test<-sum(m3*m4)
list(p.value=test,p1=r1/n1,p2=r2/n2,est.dif=r1/n1-r2/n2)
}

#' Multiple Comparisons via Percentile Bootstrap (Rom's Method)
#'
#' Computes confidence intervals for linear contrasts of trimmed means across
#' multiple independent groups using the percentile bootstrap method. Uses
#' Rom's (1990) method to control the familywise error rate.
#'
#' @inheritParams common-params
#' @param x Data in list mode (each element is a group) or matrix (columns are groups).
#' @param crit Critical value for Rom's method (default: `NA`, computed automatically for tr=0.2).
#' @param con Matrix of contrast coefficients (rows = groups, columns = contrasts).
#'   If 0, performs all pairwise comparisons (default: 0).
#' @param tr Trimming proportion (default: 0.2).
#' @param alpha Significance level (default: 0.05).
#' @param nboot Number of bootstrap samples (default: 2000).
#' @param grp Vector specifying which groups to analyze (default: `NA`, uses all groups).
#' @param WIN Logical. If `TRUE`, uses Winsorized means instead of trimmed means (default: `FALSE`).
#' @param win Winsorizing proportion when `WIN=TRUE` (default: 0.1).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `num.sig`: Number of significant contrasts
#'     \item `con`: Matrix of contrast coefficients used
#'     \item `psihat`: Vector of contrast estimates
#'     \item `test`: Test statistics
#'     \item `crit.val`: Critical value(s) used
#'     \item `p.value`: P-values for each contrast
#'     \item Various other components depending on contrast type
#'   }
#'
#' @details
#' This function performs multiple comparisons of trimmed means using Rom's (1990)
#' sequentially rejective procedure, which controls the familywise Type I error rate.
#' The percentile bootstrap is used to obtain p-values and confidence intervals.
#'
#' By default, all pairwise comparisons are performed. Custom linear contrasts can
#' be specified via the `con` matrix, where each column represents one contrast.
#'
#' When `WIN=TRUE`, Winsorized means are used instead of trimmed means. In this
#' case, the trimming proportion should typically be at least 0.2.
#'
#' @references
#' Rom, D. M. (1990). A sequentially rejective test procedure based on a modified
#' Bonferroni inequality. Biometrika, 77, 663-665.
#'
#' @seealso \code{\link{lincon}}, \code{\link{pb2gen}}, \code{\link{tmean}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare four groups
#' set.seed(123)
#' x <- list(
#'   rnorm(30, mean = 0),
#'   rnorm(30, mean = 0.5),
#'   rnorm(30, mean = 0.3),
#'   rnorm(30, mean = 0.8)
#' )
#'
#' # All pairwise comparisons
#' result <- mcppb20(x, nboot = 1000)
#' print(result)
#' }
mcppb20<-function(x,crit=NA,con=0,tr=.2,alpha=.05,nboot=2000,grp=NA,WIN=FALSE,
win=.1){
#
#   Compute a 1-alpha confidence interval for a set of d linear contrasts
#   involving trimmed means using the percentile bootstrap method.
#   Independent groups are assumed.
#
#   The data are assumed to be stored in x in list mode.  Thus,
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J, say.
#
#   By default, all pairwise comparisons are performed, but contrasts
#   can be specified with the argument con.
#   The columns of con indicate the contrast coefficients.
#   Con should have J rows, J=number of groups.
#   For example, con[,1]=c(1,1,-1,-1,0,0) and con[,2]=c(,1,-1,0,0,1,-1)
#   will test two contrasts: (1) the sum of the first two trimmed means is
#   equal to the sum of the second two, and (2) the difference between
#   the first two is equal to the difference between the trimmed means of
#   groups 5 and 6.
#
#   The default number of bootstrap samples is nboot=2000
#
#
con<-as.matrix(con)
if(is.matrix(x)){
xx<-list()
for(i in 1:ncol(x)){
xx[[i]]<-x[,i]
}
x<-xx
}
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
if(!is.na(sum(grp))){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[1]]]
x<-xx
}
J<-length(x)
tempn<-0
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
}
Jm<-J-1
d<-ifelse(sum(con^2)==0,(J^2-J)/2,ncol(con))
if(is.na(crit) && tr != .2){
print("A critical value must be specified when")
stop("the amount of trimming differs from .2")
}
if(WIN){
if(tr < .2){
print("Warning: When Winsorizing, the amount")
print("of trimming should be at least .2")
}
if(win > tr)stop("Amount of Winsorizing must <= amount of trimming")
if(min(tempn) < 15){
print("Warning: Winsorizing with sample sizes")
print("less than 15 can result in poor control")
print("over the probability of a Type I error")
}
for (j in 1:J){
x[[j]]<-winval(x[[j]],win)
}
}
if(is.na(crit)){
if(d==1)crit<-alpha/2
if(d==2 && alpha==.05 && nboot==1000)crit<-.014
if(d==2 && alpha==.05 && nboot==2000)crit<-.014
if(d==3 && alpha==.05 && nboot==1000)crit<-.009
if(d==3 && alpha==.05 && nboot==2000)crit<-.0085
if(d==3 && alpha==.025 && nboot==1000)crit<-.004
if(d==3 && alpha==.025 && nboot==2000)crit<-.004
if(d==3 && alpha==.01 && nboot==1000)crit<-.001
if(d==3 && alpha==.01 && nboot==2000)crit<-.001
if(d==4 && alpha==.05 && nboot==2000)crit<-.007
if(d==5 && alpha==.05 && nboot==2000)crit<-.006
if(d==6 && alpha==.05 && nboot==1000)crit<-.004
if(d==6 && alpha==.05 && nboot==2000)crit<-.0045
if(d==6 && alpha==.025 && nboot==1000)crit<-.002
if(d==6 && alpha==.025 && nboot==2000)crit<-.0015
if(d==6 && alpha==.01 && nboot==2000)crit<-.0005
if(d==10 && alpha==.05 && nboot<=2000)crit<-.002
if(d==10 && alpha==.05 && nboot==3000)crit<-.0023
if(d==10 && alpha==.025 && nboot<=2000)crit<-.0005
if(d==10 && alpha==.025 && nboot==3000)crit<-.001
if(d==15 && alpha==.05 && nboot==2000)crit<-.0016
if(d==15 && alpha==.025 && nboot==2000)crit<-.0005
if(d==15 && alpha==.05 && nboot==5000)crit<-.0026
if(d==15 && alpha==.025 && nboot==5000)crit<-.0006
}
if(is.na(crit) && alpha==.05)crit<-0.0268660714*(1/d)-0.0003321429
if(is.na(crit))crit<-alpha/(2*d)
if(d> 10 && nboot <5000){
print("Warning: Suggest using nboot=5000")
print("when the number of contrasts exceeds 10.")
}
icl<-round(crit*nboot)+1
icu<-round((1-crit)*nboot)
if(sum(con^2)==0){
con<-matrix(0,J,d)
id<-0
for (j in 1:Jm){
jp<-j+1
for (k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
psihat<-matrix(0,ncol(con),6)
dimnames(psihat)<-list(NULL,c("con.num","psihat","se","ci.lower",
"ci.upper","p-value"))
if(nrow(con)!=length(x)){
print("The number of groups does not match")
stop("the number of contrast coefficients.")
}
bvec<-matrix(NA,nrow=J,ncol=nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
for(j in 1:J){
paste("Working on group ",j)
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,mean,tr) # Bootstrapped trimmed means for jth group
}
test<-NA
for (d in 1:ncol(con)){
top<-0
for (i in 1:J){
top<-top+con[i,d]*bvec[i,]
}
test[d]<-(sum(top>0)+.5*sum(top==0))/nboot
test[d]<-min(test[d],1-test[d])
top<-sort(top)
psihat[d,4]<-top[icl]
psihat[d,5]<-top[icu]
}
for (d in 1:ncol(con)){
psihat[d,1]<-d
testit<-lincon(x,con[,d],tr,pr=FALSE)
psihat[d,6]<-2*test[d]
psihat[d,2]<-testit$psihat[1,2]
psihat[d,3]<-testit$test[1,4]
}
list(psihat=psihat,crit.p.value=2*crit,con=con)
}

#' Compare Variances of Two Dependent Groups
#'
#' Tests whether two dependent (paired) groups have equal variances using
#' a bootstrap approach based on the correlation between U = X - Y and V = X + Y.
#'
#' @param x Numeric vector for group 1 (or first measurement).
#' @param y Numeric vector for group 2 (or second measurement). Must have
#'   the same length as `x` since observations are paired.
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `n`: Sample size (number of complete pairs)
#'     \item `ci`: Bootstrap confidence interval for the correlation between U and V
#'   }
#'
#' @details
#' This function tests for equality of variances in paired data using a clever
#' transformation:
#'
#' 1. Computes U = X - Y (difference scores)
#' 2. Computes V = X + Y (sum scores)
#' 3. Tests whether cor(U, V) = 0
#'
#' Under the null hypothesis of equal variances (Var(X) = Var(Y)), the
#' correlation between U and V equals zero. If the variances differ, this
#' correlation will be non-zero.
#'
#' The function uses `pcorb` to compute a bootstrap confidence interval for
#' the correlation coefficient with 599 bootstrap samples. If the confidence
#' interval excludes zero, we reject the hypothesis of equal variances.
#'
#' **Advantages**:
#' \itemize{
#'   \item Works with dependent/paired data
#'   \item Distribution-free (no normality assumption)
#'   \item Robust to outliers
#' }
#'
#' @note
#' Missing values are automatically removed using pairwise deletion. The number
#' of bootstrap samples is fixed at 599.
#'
#' @seealso \code{\link{comvar2}} for independent groups,
#'   \code{\link{pcorb}} for the underlying bootstrap correlation test
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate paired data with unequal variances
#' set.seed(123)
#' x <- rnorm(50, mean = 0, sd = 1)
#' y <- rnorm(50, mean = 0, sd = 2)  # Larger variance
#'
#' # Test for equal variances
#' result <- comvar2d(x, y)
#' print(result)
#'
#' # If CI excludes 0, variances differ
#' if (result$ci[1] > 0 || result$ci[2] < 0) {
#'   cat("Variances are significantly different\n")
#' }
#' }
comvar2d<-function(x,y,SEED=TRUE){
#
#  Compare the variances of two dependent groups.
#
nboot<-599
m<-cbind(x,y)
m<-elimna(m) # Remove missing values
U<-m[,1]-m[,2]
V<-m[,1]+m[,2]
ci<-pcorb(U,V,SEED=SEED)$ci
list(n=nrow(m),ci=ci)
}

#' Multivariate WMW Test Using Projection onto Line Connecting Centers
#'
#' Compares two multivariate groups by projecting data onto the line connecting
#' the group centers, then estimating P(X < Y) using the projected distances.
#' This provides a distribution-free multivariate two-sample test.
#'
#' @param m1 Matrix or data frame for group 1 (must have 2+ columns).
#' @param m2 Matrix or data frame for group 2 (must have 2+ columns).
#' @param plotit Logical. If `TRUE`, creates a plot of projected data (default: `TRUE`).
#' @param cop Method for computing the center of each group (default: 3):
#'   \itemize{
#'     \item 1 = Donoho-Gasko (Tukey) median
#'     \item 2 = MCD (Minimum Covariance Determinant) center
#'     \item 3 = Median of marginal distributions
#'     \item 4 = Modified M-estimator (smean)
#'   }
#' @param alpha Significance level for confidence interval (default: 0.05).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param pop Plot type when `plotit = TRUE` (default: 4):
#'   \itemize{
#'     \item 1 = Two dotplots of projected distances
#'     \item 2 = Boxplots
#'     \item 3 = Expected frequency curve
#'     \item 4 = Adaptive kernel density
#'   }
#' @param fr Fraction parameter for span in plots (default: 0.8).
#' @param pr Logical. If `TRUE`, prints detailed output (default: `FALSE`).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param tr Trimming proportion for center computation when `cop = 3` (default: 0.5).
#' @param NC Logical. If `FALSE`, critical values not computed (default: `TRUE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.value`: P-value for testing H0: P(X < Y) = 0.5
#'     \item `ci`: Bootstrap confidence interval for P(X < Y)
#'     \item `p.hat`: Estimated probability that group 1 < group 2
#'     \item `n1`, `n2`: Sample sizes
#'     \item `proj1`, `proj2`: Projected distances for each group
#'   }
#'
#' @details
#' This function reduces the multivariate two-sample problem to a univariate one
#' by projecting both groups onto the line connecting their centers. The projection
#' approach is:
#'
#' 1. Compute the center of each group using the method specified by `cop`
#' 2. Project all observations onto the line connecting these two centers
#' 3. Use the WMW test on the projected distances
#' 4. Estimate P(projected distance from group 1 < projected distance from group 2)
#'
#' The choice of center (`cop`) affects robustness:
#' \itemize{
#'   \item **cop = 1**: Donoho-Gasko median (most robust, marked with +)
#'   \item **cop = 2**: MCD center (robust, good for elliptical data)
#'   \item **cop = 3**: Componentwise median (simple, fairly robust)
#'   \item **cop = 4**: Modified M-estimator (balance of efficiency and robustness)
#' }
#'
#' Bootstrap confidence intervals are computed for the probability estimate.
#'
#' @note
#' Data must have at least 2 columns. For univariate data, use `wmw.ref.dif`
#' or related functions. Missing values are automatically removed.
#'
#' @seealso \code{\link{mwmw}} for related multivariate WMW approach,
#'   \code{\link{wmw.ref.dif}} for univariate WMW,
#'   \code{\link{smean}} for modified M-estimator
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate bivariate data
#' set.seed(123)
#' m1 <- cbind(rnorm(50), rnorm(50))
#' m2 <- cbind(rnorm(50, mean = 0.5), rnorm(50, mean = 0.3))
#'
#' # Compare groups with default settings (median of marginals)
#' result <- mulwmw(m1, m2)
#'
#' # Use Donoho-Gasko median
#' mulwmw(m1, m2, cop = 1)
#'
#' # Use MCD center with boxplots
#' mulwmw(m1, m2, cop = 2, pop = 2)
#'
#' # Increase bootstrap samples
#' mulwmw(m1, m2, nboot = 2000)
#' }
mulwmw<-function(m1,m2,plotit=TRUE,cop=3,alpha=.05,nboot=1000,pop=4,fr=.8,pr=FALSE,SEED=TRUE,tr=.5,NC=TRUE){
#
#
# Determine center correpsonding to two
# independent groups, project all  points onto line
# connecting the centers,
# then based on the projected distances,
# estimate p=probability that a randomly sampled
# point from group 1 is less than a point from group 2
# based on the projected distances.
#
# plotit=TRUE creates a plot of the projected data
# pop=1 plot two dotplots based on projected distances
# pop=2 boxplots
# pop=3 expected frequency curve.
# pop=4 adaptive kernel density
#
#  There are three options for computing the center of the
#  cloud of points when computing projections:
#  cop=1 uses Donoho-Gasko median
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#
#  When using cop=2 or 3, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
#  NC=F: critical values not computed
#
#  Donoho-Gasko (Tukey) median is marked with a cross, +.
#
if(is.null(dim(m1))||dim(m1)[2]<2){print("Data are assumed to be stored in")
print(" a matrix or data frame having two or more columns.")
stop(" For univariate data, use the function outbox or out")
}
m1<-elimna(m1) # Remove missing values
m2<-elimna(m2)
n1=nrow(m1)
n2=nrow(m2)
if(cop==1){
if(ncol(m1)>2){
center1<-dmean(m1,tr=.5)
center2<-dmean(m2,tr=.5)
}
if(ncol(m1)==2){
tempd<-NA
for(i in 1:nrow(m1))
tempd[i]<-depth(m1[i,1],m1[i,2],m1)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center1<-m1[flag,]
if(sum(flag)>1)center1<-apply(m1[flag,],2,mean)
for(i in 1:nrow(m2))
tempd[i]<-depth(m2[i,1],m2[i,2],m2)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center2<-m2[flag,]
if(sum(flag)>1)center2<-apply(m2[flag,],2,mean)
}}
if(cop==2){
center1<-cov.mcd(m1)$center
center2<-cov.mcd(m2)$center
}
if(cop==3){
center1<-apply(m1,2,mean,tr=tr)
center2<-apply(m2,2,mean,tr=tr)
}
if(cop==4){
center1<-smean(m1)
center2<-smean(m2)
}
center<-(center1+center2)/2
B<-center1-center2
if(sum(center1^2)<sum(center2^2))B<-(0-1)*B
BB<-B^2
bot<-sum(BB)
disx<-NA
disy<-NA
if(bot!=0){
for (j in 1:nrow(m1)){
AX<-m1[j,]-center
tempx<-sum(AX*B)*B/bot
disx[j]<-sign(sum(AX*B))*sqrt(sum(tempx^2))
}
for (j in 1:nrow(m2)){
AY<-m2[j,]-center
tempy<-sum(AY*B)*B/bot
disy[j]<-sign(sum(AY*B))*sqrt(sum(tempy^2))
}
}
if(plotit){
if(pop==1){
par(yaxt="n")
xv<-rep(2,length(disx))
yv<-rep(1,length(disy))
plot(c(disx,disy),c(xv,yv),type="n",xlab="",ylab="")
xv<-rep(1.6,length(disx))
yv<-rep(1.4,length(disy))
points(disx,xv)
points(disy,yv)
}
if(pop==2)boxplot(disx,disy)
if(pop==3)rd2plot(disx,disy,fr=fr)
if(pop==4)g2plot(disx,disy,fr=fr)
}
m<-outer(disx,disy,FUN="-")
m<-sign(m)
phat<-(1-mean(m))/2
if(bot==0)phat<-.5
if(pr)print("Computing critical values")
m1<-t(t(m1)-center1)
m2<-t(t(m2)-center2)
v1=c(NA,NA)
if(NC)v1<-mulwmwcrit(m1,m2,cop=cop,alpha=alpha,iter=nboot,pr=pr,SEED=SEED)
list(phat=phat,lower.crit=v1[1],upper.crit=v1[2],n1=n1,n2=n2)
}

#' Determine Critical Values for Multivariate WMW Test (Internal Helper)
#'
#' Computes bootstrap critical values for the multivariate Wilcoxon-Mann-Whitney
#' test implemented in `mulwmw`. This is an internal helper function.
#'
#' @param mm1 Matrix of data for group 1 (rows are observations, columns are variables).
#' @param mm2 Matrix of data for group 2 (rows are observations, columns are variables).
#' @param plotit Logical. If `TRUE`, creates a plot (default: `TRUE`).
#' @param cop Integer specifying centering method (1=halfspace depth, 2=MCD, 3=trimmed mean, 4=marginal median; default: 3).
#' @param iter Number of bootstrap iterations (default: 1000).
#' @param alpha Significance level (default: 0.05).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param pr Logical. If `TRUE`, prints progress messages (default: `FALSE`).
#'
#' @return A two-element vector with lower and upper critical values.
#'
#' @keywords internal
#' @seealso \code{\link{mulwmw}}
mulwmwcrit<-function(mm1,mm2,plotit=TRUE,cop=3,iter=1000,alpha=.05,SEED=TRUE,pr=FALSE){
#
#
# Determine critical value for the function mulwmw
#
if(!is.matrix(mm1))stop("Data are assumed to be stored in a matrix having two or more columns. For univariate data, use the function outbox or out")
#if(is.na(SEED))set.seed(2)
#if(!is.na(SEED))set.seed(SEED)
if(SEED)set.seed(2)
val<-NA
n1<-nrow(mm1)
n2<-nrow(mm2)
for(it in 1:iter){
ivec1<-sample(c(1:n1),replace=TRUE)
ivec2<-sample(c(1:n2),replace=TRUE)
m1<-mm1[ivec1,]
m2<-mm2[ivec2,]
if(cop==1){
if(ncol(m1)>2){
center1<-dmean(m1,tr=.5)
center2<-dmean(m2,tr=.5)
}
if(ncol(m1)==2){
tempd<-NA
for(i in 1:nrow(m1))
tempd[i]<-depth(m1[i,1],m1[i,2],m1)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center1<-m1[flag,]
if(sum(flag)>1)center1<-apply(m1[flag,],2,mean)
for(i in 1:nrow(m2))
tempd[i]<-depth(m2[i,1],m2[i,2],m2)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center2<-m2[flag,]
if(sum(flag)>1)center2<-apply(m2[flag,],2,mean)
}}
if(cop==2){
center1<-cov.mcd(m1)$center
center2<-cov.mcd(m2)$center
}
if(cop==3){
center1<-apply(m1,2,median)
center2<-apply(m2,2,median)
}
center<-(center1+center2)/2
B<-center1-center2
if(sum(center1^2)>sum(center2^2))B<-(0-1)*B
BB<-B^2
bot<-sum(BB)
disx<-NA
disy<-NA
if(bot!=0){
for (j in 1:nrow(m1)){
AX<-m1[j,]-center
tempx<-sum(AX*B)*B/bot
disx[j]<-sign(sum(AX*B))*sqrt(sum(tempx^2))
}
for (j in 1:nrow(m2)){
AY<-m2[j,]-center
tempy<-sum(AY*B)*B/bot
disy[j]<-sign(sum(AY*B))*sqrt(sum(tempy^2))
}}
m<-outer(disx,disy,FUN="-")
m<-sign(m)
val[it]<-(1-mean(m))/2
if(bot==0)val[it]<-.5
if(pr)print(paste("Iteration ",it," of ",iter," complete"))
}
val<-sort(val)
low<-round(alpha*iter/2)+1
up<-iter-low
crit<-NA
crit[1]<-val[low]
crit[2]<-val[up]
crit
}

#' Multivariate Analog of Wilcoxon-Mann-Whitney Test
#'
#' Computes a multivariate effect size measure analogous to the WMW test
#' using halfspace depth. Provides a distribution-free way to compare
#' two multivariate groups.
#'
#' @param m1 Matrix or data frame for group 1 (must have 2+ columns).
#' @param m2 Matrix or data frame for group 2 (must have 2+ columns).
#' @param cop Method for computing the center when approximating halfspace depth
#'   (default: 5):
#'   \itemize{
#'     \item 1 = Halfspace median
#'     \item 2 = MCD (Minimum Covariance Determinant)
#'     \item 3 = Marginal medians
#'     \item 4 = MVE (Minimum Volume Ellipsoid)
#'     \item 5 = Skipped mean
#'   }
#' @param pr Logical. If `TRUE`, prints critical values and guidance (default: `TRUE`).
#' @param plotit Logical. If `TRUE`, creates a plot (default: `TRUE`).
#' @param pop Plot type when `plotit = TRUE` (default: 1):
#'   \itemize{
#'     \item 1 = Scatterplot
#'     \item 2 = Expected frequency curve
#'     \item 3 = Adaptive kernel density
#'   }
#' @param fr Fraction parameter for plotting (default: 0.8).
#' @param op Operation method (default: 1). Controls depth calculation approach.
#' @param dop Depth operation method for halfspace depth approximation (default: 1):
#'   \itemize{
#'     \item 1 = Method A1 approximation
#'     \item 2 = Method A2 approximation
#'   }
#'
#' @return A list with components:
#'   \itemize{
#'     \item `phat`: Effect size measure (relative depth of zero vector)
#'     \item `center`: Center of the difference distribution
#'     \item `crit.val`: Critical values for alpha = 0.1, 0.05, 0.025, 0.01
#'   }
#'
#' @details
#' This function extends the WMW test to multivariate data using halfspace depth.
#' The approach:
#'
#' 1. Computes all pairwise differences between groups: D_ij = X_i - Y_j
#' 2. Calculates the halfspace depth of the zero vector relative to these differences
#' 3. The effect size `phat` is the relative depth of zero
#'
#' **Interpretation**:
#' \itemize{
#'   \item `phat` near 1 indicates the groups are similar (zero is central to differences)
#'   \item `phat` near 0 indicates the groups differ substantially
#'   \item Reject H0 (groups equal) if `phat` <= critical value
#' }
#'
#' Critical values are provided for alpha levels 0.1, 0.05, 0.025, and 0.01,
#' based on the asymptotic distribution. The test statistic approximately
#' follows a normal distribution for large samples.
#'
#' When plotting (`plotit = TRUE`), the function displays the distribution of
#' pairwise differences, with the center marked by "o" and the zero vector by "+".
#'
#' @note
#' For univariate data (single column), use `cid` or `bmp` instead.
#' Data must have matching numbers of columns in both groups.
#'
#' @seealso \code{\link{mulwmw}} for projection-based multivariate WMW,
#'   \code{\link{cid}} for univariate Cliff's analog,
#'   \code{\link{dmean}} for Donoho-Gasko median
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate bivariate data
#' set.seed(123)
#' m1 <- cbind(rnorm(40), rnorm(40))
#' m2 <- cbind(rnorm(40, mean = 0.3), rnorm(40, mean = 0.2))
#'
#' # Multivariate WMW test
#' result <- mwmw(m1, m2)
#' print(result)
#'
#' # Use MCD center
#' mwmw(m1, m2, cop = 2)
#'
#' # With adaptive kernel density plot
#' mwmw(m1, m2, pop = 3)
#' }
mwmw<-function(m1,m2,cop=5,pr=TRUE,plotit=TRUE,pop=1,fr=.8,op=1,dop=1){
#
# Compute measure of effect size, p,
# a multivariate analog of Wilcoxon-Mann-Whitney p
#
# When plotting:
# pop=1 Use scatterplot
# pop=2 Use expected frequency curve.
# pop=3 Use adaptive kernel density
#
# dop=1, use method A1 approximation of halfspace depth
# dop=2, use method A2 approximation of halfspace depth
#
# cop determines how center of data is determined when
# approximating halfspace depth
# cop=1, Halfspace median
# cop=2, MCD
# cop=3, marginal medians
# cop=4, MVE
# cop=5, skipped mean
#
if(is.null(dim(m1)))stop("m1 is not a matrix or data frame")
if(is.null(dim(m2)))stop("m2 is not a matrix or data frame")
if(ncol(m1)!=ncol(m2))stop("number of columns for m1 and m2 are not equal")
if(ncol(m1)==1)stop("Use R function cid or bmp")
nn<-min(c(nrow(m1),nrow(m2)))
mdif<-matrix(as.vector(outer(m1[,1],m2[,1],"-")),ncol=1)
for(j in 2:ncol(m1)){
mdif<-cbind(mdif,matrix(as.vector(outer(m1[,j],m2[,j],"-")),ncol=1))
}
if(op==1){
if(ncol(m1)==2)temp2<-depth2(rbind(mdif,c(rep(0,ncol(m1)))))
#if(ncol(m1)==3)temp2<-depth3(rbind(mdif,c(rep(0,ncol(m1)))))
if(ncol(m1)>2){
if(cop==1)center<-dmean(mdif,tr=.5,dop=dop)
if(cop==2)center<-cov.mcd(mdif)$center
if(cop==3)center<-apply(mdif,2,median)
if(cop==4)center<-cov.mve(mdif)$center
if(cop==5)center<-smean(mdif)
temp2<-fdepth(rbind(mdif,c(rep(0,ncol(m1)))))
}}
if(op==2){
temp2<-pdis(rbind(mdif,c(rep(0,ncol(m1)))))
temp2<-1/(temp2+1)
}
center<-dmean(mdif,tr=.5,dop=dop)
phat<-temp2[nrow(mdif)+1]/max(temp2)
# phat is relative depth of zero vector
# Determine critical value
crit<-NA
alpha<-c(.1,.05,.025,.01)
crit[1]<-1-1.6338/sqrt(nn)
crit[2]<-1-1.8556/sqrt(nn)
crit[3]<-1-2.0215/sqrt(nn)
crit[4]<-1-2.1668/sqrt(nn)
if(pr){
print("For alpha=.1,.05,.025,.01, the correspoding critical values are")
print(crit)
print("Reject if phat is less than or equal to the critical value")
}
if(plotit && ncol(m1)==2){
if(pop==2)rdplot(mdif,fr=fr)
if(pop==1){
plot(mdif[,1],mdif[,2],xlab="VAR 1",ylab="VAR 2",type="n")
points(mdif[,1],mdif[,2],pch=".")
points(center[1],center[2],pch="o")
points(0,0,pch="+")
}
if(pop==3)akerdmul(mdif,fr=fr)
}
list(phat=phat,center=center,crit.val=crit)
}

#' Compare Two Dependent Pearson Correlations (Non-Overlapping)
#'
#' @description
#' Computes a confidence interval for the difference between two dependent Pearson
#' correlations in the non-overlapping case, where the correlations involve completely
#' different pairs of variables.
#'
#' @inheritParams common-params
#' @param x Matrix or data frame with 2 columns.
#' @param y Matrix or data frame with 2 columns.
#' @param HC4 Logical; if TRUE, uses HC4 heteroscedasticity correction. Required for
#'   alpha not equal to 0.05.
#'
#' @return A list with components:
#'   \item{est.1}{Correlation between x[,1] and x[,2]}
#'   \item{est.2}{Correlation between y[,1] and y[,2]}
#'   \item{ci.lower}{Lower confidence limit for the difference}
#'   \item{ci.upper}{Upper confidence limit for the difference}
#'
#' @details
#' This function compares the correlation between x[,1] and x[,2] to the correlation
#' between y[,1] and y[,2]. This is the "non-overlapping" case because the two
#' correlations involve completely different pairs of variables (though the data
#' may be from the same subjects).
#'
#' The method deals with heteroscedasticity and non-normality. For alpha != 0.05,
#' you must use HC4=TRUE.
#'
#' @references
#' Wilcox, R.R. (2009). Comparing Pearson correlations: Dealing with heteroscedasticity
#' and non-normality. Communications in Statistics - Simulation and Computation, 38,
#' 2220-2234.
#'
#' @seealso \code{\link{TWOpov}}, \code{\link{twoDNOV}}, \code{\link{tworhobt}}
#'
#' @export
#' @examples
#' set.seed(123)
#' n <- 50
#' x <- cbind(rnorm(n), rnorm(n))
#' y <- cbind(rnorm(n), rnorm(n))
#'
#' result <- TWOpNOV(x, y)
#' result$ci.lower
#' result$ci.upper
TWOpNOV<-function(x,y,HC4=FALSE,alpha=.05){
#
#   Compute a .95 confidence interval
#   for the difference between two dependent Pearson correlations,
#   non-overlapping case.
#
#    Both x and y are assumed to be matrices with two columns.
#   The function compares the correlation between x[,1] and x[,2]
#   to the correlation between y[,1] and y[,2].
#
#  For simulation results, see Wilcox (2009).
#  COMPARING PEARSON CORRELATIONS: DEALING WITH
#  HETEROSCEDASTICITY AND NON-NORMALITY, Communications in Statistics--Simulations
#   and Computations, 38, 2220-2234.
#
#
if(!HC4 && alpha!=.05)stop('For alpha not equal to .05, must use HC4=TRUE')
#if(!is.matrix(x))stop("x should be a matrix")
#if(!is.matrix(y))stop("y should be a matrix")
if(ncol(x)!=2)stop("x should be a matrix or data a frame with 2 columns")
if(ncol(y)!=2)stop("y should be a matrix or a data frame with 2 columns")
xy=elimna(cbind(x,y))
x1=xy[,1]
x2=xy[,2]
y1=xy[,3]
y2=xy[,4]
r12=cor(x1,x2)
r13=cor(x1,y1)
r14=cor(x1,y2)
r23=cor(x2,y1)
r24=cor(x2,y2)
r34=cor(y1,y2)
term1=.5*r12*r34*(r13^2+r14^2+r23^2+r24^2)
term2=r12*r13*r14+r12*r23*r24+r13*r23*r34+r14*r24*r34
corhat=(term1+r13*r24+r14*r23-term2)/((1-r12^2)*(1-r34^2))
if(!HC4)temp=pcorbv4(x1,x2,SEED=FALSE)
if(HC4)temp=pcorhc4(x1,x2,alpha=alpha)
ci12=temp$ci[1]
ci12[2]=temp$ci[2]
if(!HC4)temp=pcorbv4(y1,y2,SEED=FALSE)
if(HC4)temp=pcorhc4(y1,y2,alpha=alpha)
ci34=temp$ci[1]
ci34[2]=temp$ci[2]
terml=2*corhat*(r12-ci12[1])*(ci34[2]-r34)
termu=2*corhat*(ci12[2]-r12)*(r34-ci34[1])
L=r12-r34-sqrt((r12-ci12[1])^2+(ci34[2]-r34)^2-terml)
U=r12-r34+sqrt((r12-ci12[2])^2+(ci34[1]-r34)^2-termu)
if(ZCI){
if(is.na(L) || is.na(U))L=U=0
}
list(est.1=r12,est.2=r34,ci.lower=L,ci.upper=U)
}

#' Compare Two Dependent Pearson Correlations (Overlapping)
#'
#' @description
#' Computes a confidence interval for the difference between two dependent Pearson
#' correlations in the overlapping case, where both correlations share a common variable.
#'
#' @inheritParams common-params
#' @param x Matrix or data frame with 2 columns.
#' @param y Numeric vector.
#' @param CN Logical; if TRUE, uses Chen-Nadarajah method.
#' @param BOOT Logical; if TRUE, uses bootstrap method.
#' @param nboot Number of bootstrap samples (default 499).
#' @param ZCI Logical; internal parameter for handling NA values in CI.
#'
#' @return A list with components:
#'   \item{est.1}{Correlation between x[,1] and y}
#'   \item{est.2}{Correlation between x[,2] and y}
#'   \item{ci.lower}{Lower confidence limit for the difference}
#'   \item{ci.upper}{Upper confidence limit for the difference}
#'   \item{p.value}{P-value for testing equality (if BOOT=TRUE)}
#'
#' @details
#' This function compares the correlation of x[,1] with y to the correlation of x[,2]
#' with y. This is the "overlapping" case because both correlations share the common
#' variable y.
#'
#' The default method uses a bootstrap approach (BOOT=TRUE). Alternatively, CN=TRUE
#' uses the Chen-Nadarajah method. The method deals with heteroscedasticity and
#' non-normality.
#'
#' @references
#' Wilcox, R.R. (2009). Comparing Pearson correlations: Dealing with heteroscedasticity
#' and non-normality. Communications in Statistics - Simulation and Computation, 38,
#' 2220-2234.
#'
#' @seealso \code{\link{TWOpNOV}}, \code{\link{twoDcorR}}, \code{\link{tworhobt}}
#'
#' @export
#' @examples
#' set.seed(123)
#' n <- 50
#' x <- cbind(rnorm(n), rnorm(n))
#' y <- rnorm(n)
#'
#' result <- TWOpov(x, y)
#' result$ci.lower
#' result$ci.upper
TWOpov<-function(x,y,alpha=.05,CN=FALSE,BOOT=TRUE, nboot=499,SEED=TRUE,ZCI=FALSE){
#
# Comparing two dependent correlations: Overlapping case
#
# x is assumed to be a matrix with 2 columns
#
#  Compare correlation of x[,1] with y to x[,2] with y
#
#  returns a confidence stored in
#  ci
#
if(ncol(x)!=2)stop('x should be a matrix with two columns')
x1y=elimna(cbind(x[,1],y))
x2y=elimna(cbind(x[,2],y))
xx=elimna(x)
r12=cor(x1y[,1],x1y[,2])
r13=cor(x2y[,1],x2y[,2])
r23=cor(xx[,1],xx[,2])
if(!BOOT){
ci12=pcorhc4(x1y[,1],x1y[,2],alpha=alpha,CN=CN)$ci
ci13=pcorhc4(x2y[,1],x2y[,2],alpha=alpha,CN=CN)$ci
}
if(BOOT){
ci12=rhohc4bt(x1y[,1],x1y[,2],alpha=alpha,SEED=SEED,nboot=nboot)$ci
ci13=rhohc4bt(x2y[,1],x2y[,2],alpha=alpha,SEED=SEED,nboot=nboot)$ci
}
corhat=((r23-.5*r12*r13)*(1-r12^2-r13^2-r23^2)+r23^3)/((1-r12^2)*(1-r13^2))
term1=2*corhat*(r12-ci12[1])*(ci13[2]-r13)
term2=2*corhat*(r12-ci12[2])*(ci13[1]-r13)
L=r12-r13-sqrt((r12-ci12[1])^2+(ci13[2]-r13)^2-term1)
U=r12-r13+sqrt((r12-ci12[2])^2+(ci13[1]-r13)^2-term2)
if(ZCI){
if(is.na(L) || is.na(U))L=U=0
}
list(est.rho1=r12,est.rho2=r13,dif=r12-r13,ci=c(L,U))
}

#' Percentile Bootstrap for Trimmed Mean Comparison
#'
#' Enhanced version of trimmed mean comparison with bootstrap,  supporting both
#' one-sample and two-sample designs (independent or dependent groups).
#'
#' @inheritParams common-params
#' @param WIN Logical. If `TRUE`, uses Winsorizing (default: `FALSE`).
#' @param win Amount of Winsorizing when `WIN = TRUE` (default: 0.1).
#' @param pop Plot type: 1=expected frequency, 2=kernel density, 3=boxplot, 4=stem-and-leaf (default: 1).
#' @param null.value Null hypothesis value for testing (default: 0).
#' @param xlab Label for x-axis (default: "X").
#' @param fr Fraction parameter (default: NA).
#'
#' @details
#' This function compares trimmed means using percentile bootstrap. It can handle:
#' \itemize{
#'   \item One-sample test when only `x` is provided
#'   \item Independent two-sample test when `x` and `y` are separate vectors
#'   \item Dependent (paired) test when `y = NULL` and `x` is a two-column matrix
#' }
#'
#' Missing values are automatically removed.
#'
#' When `plotit = TRUE`, creates a plot of the bootstrap distribution with type
#' specified by `pop`.
#'
#' @return A list with components including estimates, confidence interval, and p-value.
#'
#' @seealso \code{\link{trimpb2}}, \code{\link{pb2gen}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Two independent groups
#' x <- rnorm(50)
#' y <- rnorm(60, mean = 0.5)
#' trimpb(x, y)
#'
#' # Dependent groups (paired)
#' before <- rnorm(30)
#' after <- before + rnorm(30, mean = 0.3)
#' trimpb(cbind(before, after))
#' }
trimpb<-function(x,y=NULL,tr=.2,alpha=.05,nboot=2000,WIN=FALSE,win=.1,
plotit=FALSE,pop=1,null.value=0,pr=TRUE,xlab="X",fr=NA,SEED=TRUE){
#
#   Compute a 1-alpha confidence interval for
#   a trimmed mean.
#
#   The default number of bootstrap samples is nboot=2000
#
#   win is the amount of Winsorizing before bootstrapping
#   when WIN=T.
#
#   Missing values are automatically removed.
#
#  nv is null value. That test hypothesis trimmed mean equals nv
#
#  plotit=TRUE gives a plot of the bootstrap values
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate
#  pop=3 boxplot
#  pop=4 stem-and-leaf
#  pop=5 histogram
#  pop=6 adaptive kernel density estimate.
#
#  fr controls the amount of smoothing when plotting the bootstrap values
#  via the function rdplot. fr=NA means the function will use fr=.8
#  (When plotting bivariate data, rdplot uses fr=.6 by default.)
#
#  If y is not null, the function uses x-y; so can be used for two dependent variables.
#
if(pr){
print("The p-value returned by this function is based on the")
print("null value specified by the argument null.value, which defaults to 0")
}
if(!is.null(y))x=x-y
x<-x[!is.na(x)]
if(WIN){
if(win > tr)stop("The amount of Winsorizing must be <= to the amount of trimming")
x<-winval(x,win)
}
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
bvec<-NA
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,mean,tr) # Bootstrapped trimmed means
bvec<-sort(bvec)
#p.value<-sum(bvec<null.value)/nboot
p.value<-mean(bvec<null.value)+.5*mean(bvec==null.value)
p.value<-2*min(p.value,1-p.value)
ci<-NA
ci[1]<-bvec[icl]
ci[2]<-bvec[icu]
if(plotit){
if(pop==1)rdplot(as.vector(bvec),fr=fr,xlab=xlab)
if(pop==2)kdplot(as.vector(bvec),rval=rval)
if(pop==3)boxplot(as.vector(bvec))
if(pop==4)stem(as.vector(bvec))
if(pop==5)hist(as.vector(bvec))
if(pop==6)akerd(as.vector(bvec),xlab=xlab)
}
list(estimate=mean(x,tr=tr),ci=ci,p.value=p.value)
}

#' Yuen's Test for Dependent Groups with Missing Values
#'
#' Compares trimmed means of two dependent (paired) variables allowing for
#' missing values. Uses marginal trimmed means computed from all available data.
#'
#' @inheritParams common-params
#'
#' @details
#' This function extends Yuen's test for dependent groups to handle missing data.
#' Unlike standard paired tests that delete cases with any missing values, this
#' function uses all available data by computing marginal trimmed means.
#'
#' Pairs with both values missing are deleted. For pairs with one value missing,
#' the observed value contributes to its marginal trimmed mean. A weighted
#' combination of these estimates is used.
#'
#' If `y` is not supplied, `x` is assumed to be a two-column matrix.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `estimate`: Estimated difference in trimmed means
#'     \item `test`: Test statistic
#'     \item `se`: Standard error of the difference
#'   }
#'
#' @seealso \code{\link{yuend}}, \code{\link{yuendv2}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Data with missing values
#' x <- c(10, 12, NA, 15, 18, 20, NA, 25)
#' y <- c(8, NA, 14, 13, 16, 19, 22, 23)
#'
#' result <- yuendna(x, y)
#' print(result)
#' }
yuendna<-function(x,y=NULL,tr=.2,alpha=.05){
#
#  Compare the trimmed means of two dependent random variables
#  using the data in x and y.
#  The default amount of trimming is 20%
#
# If y is not supplied, this function assumes x is a matrix with 2 columns.
#
#  pairs of observations, for which one value is missing, are NOT deleted.
#  Marginal trimmed means are compared
#  using all available data.
#
if(is.null(y)){
if(!is.matrix(x))stop("y is null and x is not a matrix")
y=x[,2]
x=x[,1]
}
if(length(x)!=length(y))stop("The number of observations must be equal")
m<-cbind(x,y)
# first eliminate any rows with both values missing.
flag=(apply(is.na(m),1,sum)==2)
m=m[!flag,]
x<-m[,1]
y<-m[,2]
flagx=is.na(y) # Indicates observed x values for which y is missing
flagy=is.na(x) # Indicates the y values for which x is missing
m<-elimna(m)   # m has data where both values are available--no missing values
n=nrow(m)
n1=sum(flagx)  # number of x values for which y is missing
n2=sum(flagy)
h=n-2*floor(tr*n)
h1=n1-2*floor(tr*n1)
h2=n2-2*floor(tr*n2)
xbarn=mean(x,tr=tr,na.rm=TRUE)
xbarn1=0
if(h1>0)xbarn1=mean(x[flagx],tr=tr)
ybarn=mean(y[!flagy],tr=tr,na.rm=TRUE)
ybarn1=0
if(h2>0)ybarn1=mean(y[flagy],tr=tr)
lam1=h/(h+h1)
lam2=h/(h+h2)
est=lam1*xbarn-lam2*ybarn+(1-lam1)*xbarn1-(1-lam2)*ybarn1
sex=trimse(elimna(x),tr=tr)
sey=trimse(elimna(y),tr=tr)
q1<-(n-1)*winvar(m[,1],tr)
q2<-(n-1)*winvar(m[,2],tr)
q3<-(n-1)*wincor(m[,1],m[,2],tr)$cov
sen=sqrt((lam1^2*q1+lam2^2*q2-2*lam1*lam2*q3)/(h*(h-1)))
SE=sqrt(sen^2+(1-lam1)^2*sex^2+(1-lam2)^2*sey^2)
test=est/SE
list(estimate=est,test=test,se=SE)
}

#' Yuen's Test for Independent Groups with Visualization
#'
#' Performs Yuen's test for comparing trimmed means between two independent
#' groups with optional visualization and various display options.
#'
#' @inheritParams common-params
#' @param plotfun Plotting function to use (default: `splot`).
#' @param op Logical. If `TRUE`, prints output (default: `TRUE`).
#' @param VL Logical. If `TRUE`, adds vertical lines to plot (default: `TRUE`).
#' @param cor.op Logical. If `TRUE`, displays correlation information (default: `FALSE`).
#' @param loc.fun Location function for plotting (default: `median`).
#' @param xlab Label for x-axis (default: "Groups").
#' @param ylab Label for y-axis (default: "").
#' @param PB Logical. If `TRUE`, adds percentile bootstrap confidence interval (default: `FALSE`).
#'
#' @details
#' This is an enhanced version of Yuen's test that provides visualization options
#' and additional diagnostic information. The test compares trimmed means of two
#' independent groups using Welch's approach for unequal variances.
#'
#' When `plotit = TRUE`, creates a plot comparing the two groups. When `PB = TRUE`,
#' also computes a percentile bootstrap confidence interval for comparison.
#'
#' Missing values are automatically removed.
#'
#' @return A list with Yuen test results including estimates, test statistic,
#'   p-value, confidence interval, and degrees of freedom.
#'
#' @seealso \code{\link{yuen}}, \code{\link{yuenbt}}
#'
#' @export
#' @examples
#' \dontrun{
#' x <- rnorm(50)
#' y <- rnorm(60, mean = 0.5)
#'
#' # Basic test with plot
#' yuenv2(x, y, plotit = TRUE)
#'
#' # With bootstrap CI
#' yuenv2(x, y, plotit = TRUE, PB = TRUE)
#' }
yuenv2<-function(x,y=NULL,tr=.2,alpha=.05,plotit=FALSE,plotfun=splot,op=TRUE, VL=TRUE,cor.op=FALSE, loc.fun=median,
xlab="Groups",ylab="",PB=FALSE,nboot=100, SEED=TRUE){
#
#  Perform Yuen's test for trimmed means on the data in x and y.
#  The default amount of trimming is 20%
#  Missing values (values stored as NA) are automatically removed.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuen$ci.
#  The significance level is returned in yuen$p.value
#
#  For an omnibus test with more than two independent groups,
#  use t1way.
#
#   Unlike the function yuen, a robust heteroscedastic measure
#   of effect size is returned.
#  PB=FALSE means that a Winsorized variation of prediction error is used to measure effect size.
#  PB=TRUE:  A percentage bend measure of variation is used instead.
#
if(tr==.5)stop("Use medpb to compare medians.")
if(tr>.5)stop("Can't have tr>.5")
if(is.null(y)){
if(is.matrix(x) || is.data.frame(x)){
y=x[,2]
x=x[,1]
}
if(is.list(x)){
y=x[[2]]
x=x[[1]]
}
}
x<-x[!is.na(x)]  # Remove any missing values in x
y<-y[!is.na(y)]  # Remove any missing values in y
n1=length(x)
n2=length(y)
h1<-length(x)-2*floor(tr*length(x))
h2<-length(y)-2*floor(tr*length(y))
q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
crit<-qt(1-alpha/2,df)
m1=mean(x,tr)
m2=mean(y,tr)
mbar=(m1+m2)/2
dif=m1-m2
low<-dif-crit*sqrt(q1+q2)
up<-dif+crit*sqrt(q1+q2)
test<-abs(dif/sqrt(q1+q2))
yuen<-2*(1-pt(test,df))
xx=c(rep(1,length(x)),rep(2,length(y)))
if(h1==h2){
pts=c(x,y)
top=var(c(m1,m2))
#
if(!PB){
if(tr==0)cterm=1
if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
bot=winvar(pts,tr=tr)/cterm
e.pow=top/bot
if(!is.na(e.pow)){
if(e.pow>1){
x0=c(rep(1,length(x)),rep(2,length(y)))
y0=c(x,y)
e.pow=wincor(x0,y0,tr=tr)$cor^2
}
}
}
#
if(PB){
bot=pbvar(pts)
e.pow=top/bot
}
#
}
if(n1!=n2){
N=min(c(n1,n2))
vals=0
if(SEED)set.seed(2)
for(i in 1:nboot)vals[i]=yuen.effect(sample(x,N),sample(y,N),tr=tr)$Var.Explained
e.pow=loc.fun(vals)
}
if(plotit){
plot(xx,pts,xlab=xlab,ylab=ylab)
if(op)
points(c(1,2),c(m1,m2))
if(VL)lines(c(1,2),c(m1,m2))
}
list(ci=c(low,up),n1=n1,n2=n2,
p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,
crit=crit,df=df,Var.Explained=e.pow,Effect.Size=sqrt(e.pow))
}

#' Explanatory Measure of Effect Size for Yuen's Test with CI
#'
#' Computes an explanatory measure of effect size for Yuen's test along with
#' a bootstrap confidence interval.
#'
#' @inheritParams common-params
#'
#' @details
#' The effect size measure indicates the degree to which distributions differ
#' in a way that is directly interpretable. Bootstrap methods are used to
#' compute the confidence interval.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `estimate`: Effect size estimate
#'     \item `ci`: Bootstrap confidence interval for the effect size
#'   }
#'
#' @seealso \code{\link{yuen.effect}}, \code{\link{akp.effect}}
#'
#' @export
#' @examples
#' \dontrun{
#' x <- rnorm(50)
#' y <- rnorm(60, mean = 0.5)
#' yuen.effect.ci(x, y)
#' }
yuen.effect.ci<-function(x,y,SEED=TRUE,nboot=400,tr=.2,alpha=.05){
#
# Compute a 1-alpha  confidence interval
# for a robust, heteroscedastic  measure of effect size
#  The absolute value of the measure of effect size is used.
#
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
x=elimna(x)
y=elimna(y)
bvec=0
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(x)*nboot,replace=TRUE),nrow=nboot)
for(i in 1:nboot){
bvec[i]=yuenv2(datax[i,],datay[i,],tr=tr,SEED=FALSE)$Effect.Size
}
bvec<-sort(abs(bvec))
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
ci<-NA
ci[1]<-bvec[icl]
pchk=yuen(x,y,tr=tr)$p.value
if(pchk>alpha)ci[1]=0
ci[2]<-bvec[icu]
if(ci[1]<0)ci[1]=0
es=abs(yuenv2(x,y,tr=tr)$Effect.Size)
list(CI=ci,Effect.Size=es)
}

#' Explanatory Measure of Effect Size for Yuen's Test
#'
#' Computes an explanatory measure of effect size when comparing trimmed means.
#'
#' @inheritParams common-params
#' @inheritParams yuenv2
#'
#' @details
#' This effect size measure is an explanatory measure that indicates the
#' probability-based interpretation of the difference between trimmed means.
#' Same as `yuen` but additionally computes explanatory power and related
#' effect size measures.
#'
#' @return A list with effect size estimate and related information.
#'
#' @seealso \code{\link{yuen.effect.ci}}, \code{\link{akp.effect}}
#'
#' @export
#' @examples
#' \dontrun{
#' x <- rnorm(50)
#' y <- rnorm(60, mean = 0.5)
#' yuen.effect(x, y)
#' }
yuen.effect<-function(x,y,tr=.2,alpha=.05,plotit=FALSE,
plotfun=splot,op=TRUE,VL=TRUE,cor.op=FALSE,
xlab="Groups",ylab="",PB=FALSE){
#
#  Same as yuen, only it computes explanatory power and the related
# measure of effect size. Only use this with n1=n2. Called by yuenv2
# which allows n1!=n2.
#
#
#  Perform Yuen's test for trimmed means on the data in x and y.
#  The default amount of trimming is 20%
#  Missing values (values stored as NA) are automatically removed.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuen$ci.
#  The p-valueis returned in yuen$p.value
#
#  For an omnibus test with more than two independent groups,
#  use t1way.
#  This function uses winvar from chapter 2.
#
if(tr==.5)stop("Use medpb to compare medians.")
if(tr>.5)stop("Can't have tr>.5")
x<-x[!is.na(x)]  # Remove any missing values in x
y<-y[!is.na(y)]  # Remove any missing values in y
h1<-length(x)-2*floor(tr*length(x))
h2<-length(y)-2*floor(tr*length(y))
q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
crit<-qt(1-alpha/2,df)
m1=mean(x,tr)
m2=mean(y,tr)
mbar=(m1+m2)/2
dif=m1-m2
low<-dif-crit*sqrt(q1+q2)
up<-dif+crit*sqrt(q1+q2)
test<-abs(dif/sqrt(q1+q2))
yuen<-2*(1-pt(test,df))
xx=c(rep(1,length(x)),rep(2,length(y)))
pts=c(x,y)
top=var(c(m1,m2))
#
if(!PB){
if(tr==0)cterm=1
if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
bot=winvar(pts,tr=tr)/cterm
}
if(PB)bot=pbvar(pts)/1.06
#
e.pow=top/bot
if(e.pow>1){
x0=c(rep(1,length(x)),rep(2,length(y)))
y0=c(x,y)
e.pow=wincor(x0,y0,tr=tr)$cor^2
}
if(plotit){
plot(xx,pts,xlab=xlab,ylab=ylab)
if(op)
points(c(1,2),c(m1,m2))
if(VL)lines(c(1,2),c(m1,m2))
}
list(ci=c(low,up),p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,
crit=crit,df=df,Var.Explained=e.pow,Effect.Size=sqrt(e.pow))
}

#' WMW Location Estimates for All Pairs of Groups
#'
#' Computes WMW-based location estimates for all pairwise group comparisons.
#'
#' @inheritParams common-params
#'
#' @details
#' This function computes pairwise WMW location differences for multiple groups.
#' Primarily an internal helper function.
#'
#' @keywords internal
wmwloc2<-function(x,est=median){
#
# Compute loc2dif for all pairs of groups
#
if(is.matrix(x))x=listm(x)
locvec=NULL
ic=0
J=length(x)
for(j in 1:J){
for(k in 1:J){
if (j<k){
ic=ic+1
locvec[ic]=loc2dif(x[[j]],x[[k]],est=est)
}}}
locvec
}

#' Multiple Group Comparisons Using Cliff's Method
#'
#' Performs all pairwise comparisons using a variation of Cliff's method based
#' on the median of X-Y. Controls familywise error rate using Hochberg's method.
#'
#' @inheritParams common-params
#' @param g Column number containing group/factor variable (default: `NULL`).
#' @param dp Column number containing dependent variable (default: `NULL`).
#'
#' @details
#' Uses p = P(X < Y) as an effect size measure and tests whether p = 0.5.
#' Performs all pairwise comparisons with familywise error rate controlled
#' via Hochberg's method.
#'
#' Data can be provided as:
#' \itemize{
#'   \item Matrix where columns are groups (when `g = NULL`)
#'   \item List where each element is a group (when `g = NULL`)
#'   \item Data frame/matrix with group indicator in column `g` and data in column `dp`
#' }
#'
#' @return Matrix of pairwise comparison results with adjusted p-values.
#'
#' @seealso \code{\link{cid}}, \code{\link{cidmul}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Matrix input
#' x <- matrix(rnorm(150), ncol = 3)
#' cidM(x)
#'
#' # With grouping variable
#' data <- data.frame(group = rep(1:3, each = 50), value = rnorm(150))
#' cidM(data, g = 1, dp = 2)
#' }
cidM<-function(x,nboot=1000,alpha=.05,MC=FALSE,SEED=TRUE,g=NULL,dp=NULL){
#
# Variation of Cliff method based on median of X-Y
# i.e., use p=P(X<Y) as effect size.
# test p=.5
# All pairwise comparisons performed.
# FWE controlled via Hochberg method.
# x can be a matrix (columns are groups) or have list mode
#
#   g=NULL, x is assumed to be a matrix or have list mode
#   if g is specifed, it is assumed that column g of x is
#   a factor variable and that the dependent variable of interest is in column
#   dp of x, which can be a matrix or data frame.
#
if(!is.null(g)){
if(is.null(dp))stop("Specify a value for dp, the column containing the data")
x=fac2list(x[,dp],x[,g])
}
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x=listm(x)
chk=tlist(x)
if(chk!=0)print("Warning: tied values detected. Suggest using cidmulv2")
J=length(x)
L=(J^2-J)/2
CC=L
pvec=NA
boot=list()
MAT=matrix(NA,nrow=nboot,ncol=L)
for(i in 1:nboot){
jcom=0
for (j in 1:J){
boot[[j]]=sample(x[[j]],size=length(x[[j]]),replace=TRUE)
}
MAT[i,]=wmwloc2(boot)
}
#
pvec=NA
test<-matrix(NA,CC,8)
dimnames(test)<-list(NULL,c("Group","Group","p-value","p.crit",
"P(X<Y)","P(X=Y)","P(X>Y)","p.hat"))
dvec<-alpha/c(1:CC)
for(j in 1:J){
for(k in 1:J){
if(j<k){
jcom=jcom+1
p.value=mean(MAT[,jcom]>0)+.5*mean(MAT[,jcom]==0)
pvec[jcom]=2*min(c(p.value,1-p.value))
if(is.na(pvec[jcom]))pvec=1
test[jcom,1]<-j
test[jcom,2]<-k
test[jcom,3]<-pvec[jcom]
test[jcom,5:7]<-cid(x[[j]],x[[k]])$summary.dvals
test[jcom,8]<-test[jcom,5]+.5*test[jcom,6]
}}}
temp2<-order(0-test[,3])
test[temp2,4]=dvec
list(test=test)
}

#' Multiple Comparisons Using Cliff's Method
#'
#' @description
#' Performs all pairwise comparisons among J independent groups using Cliff's method
#' for computing P(X<Y). Unlike some alternatives, this method allows tied values.
#' Family-wise error rate is controlled using the Studentized maximum modulus distribution.
#'
#' @inheritParams common-params
#' @param x Data in list mode, matrix, or data frame. If matrix/data frame, each column
#'   represents a group. Length(x) or ncol(x) corresponds to the total number of groups J.
#' @param g Optional grouping variable (factor) when data is in long format.
#' @param dp Column index containing the data when using grouping variable \code{g}.
#' @param pr Logical; if TRUE, prints a message suggesting \code{cidmulv2} for better power.
#'
#' @return A list with components:
#'   \item{n}{Vector of sample sizes for each group}
#'   \item{test}{Matrix with columns: Group, Group, d (dominance measure), ci.lower,
#'     ci.upper, p.hat (P(X<Y)), p-value}
#'
#' @details
#' For each pair of groups (j,k) where j < k, the function computes Cliff's dominance
#' measure d = P(X>Y) - P(X<Y) and estimates p.hat = P(X<Y). The familywise Type I
#' error probability is controlled using a critical value from the Studentized maximum
#' modulus distribution.
#'
#' The default alpha is 0.05; any other value results in using alpha = 0.01.
#'
#' Note: \code{cidmulv2} may provide better power and is recommended.
#'
#' @seealso \code{\link{cidmulv2}}, \code{\link{cid}}, \code{\link{cidv2}}, \code{\link{wmwaov}}
#'
#' @export
#' @examples
#' # Three groups
#' set.seed(123)
#' g1 <- rnorm(20, mean = 5)
#' g2 <- rnorm(20, mean = 5.5)
#' g3 <- rnorm(20, mean = 6)
#' x <- list(g1, g2, g3)
#'
#' result <- cidmul(x)
#' result$test
cidmul<-function(x,alpha=.05,g=NULL,dp=NULL,pr=TRUE){
#
#  Perform Cliff's method for all pairs of J independent groups.
#  Unlike the function meemul, ties are allowed.
#  The familywise type I error probability is controlled by using
#  a critical value from the Studentized maximum modulus distribution.
#
#  The data are assumed to be stored in $x$ in list mode.
#  Length(x) is assumed to correspond to the total number of groups, J.
#  It is assumed all groups are independent.
#
#  Missing values are automatically removed.
#
#  The default value for alpha is .05. Any other value results in using
#  alpha=.01.
#
if(pr)print('cidmulv2 might provide better power')
if(!is.null(g)){
if(is.null(dp))stop("Specify a value for dp, the column containing the data")
x=fac2list(x[,dp],x[,g])
}
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
J<-length(x)
CC<-(J^2-J)/2
test<-matrix(NA,CC,7)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
}
dimnames(test)<-list(NULL,c("Group","Group","d","ci.lower","ci.upper",
"p.hat","p-value"))
jcom<-0
crit<-smmcrit(200,CC)
if(alpha!=.05)crit<-smmcrit01(200,CC)
alpha<-1-pnorm(crit)
n=matl(lapply(x,length))
for (j in 1:J){
for (k in 1:J){
if (j < k){
temp<-cid(x[[j]],x[[k]],alpha,plotit=FALSE)
temp2<-cidv2(x[[j]],x[[k]],alpha,plotit=FALSE)
jcom<-jcom+1
test[jcom,1]<-j
test[jcom,2]<-k
test[jcom,3]<-temp$d
test[jcom,4]<-temp$cl
test[jcom,5]<-temp$cu
test[jcom,6]<-temp$phat
test[jcom,7]<-temp2$p.value
}}}
list(n=n,test=test)
}

#' Compare Two Independent Groups Using Medians
#'
#' @description
#' Compares medians of two independent groups using a percentile bootstrap method that
#' performs well when there are tied values. This is a robust alternative to parametric
#' tests for comparing central tendency.
#'
#' @inheritParams common-params
#' @param x Either a numeric vector (when \code{y} is provided), or a matrix/data frame
#'   with two columns, or a list with two elements.
#' @param y Optional numeric vector for the second group. If NULL, \code{x} must contain
#'   both groups.
#' @param nboot Number of bootstrap samples (default 2000).
#'
#' @return A list with components:
#'   \item{n1}{Sample size of group 1}
#'   \item{n2}{Sample size of group 2}
#'   \item{p.value}{P-value for testing equality of medians}
#'   \item{ci}{Confidence interval for the difference in medians}
#'   \item{est1}{Median of group 1}
#'   \item{est2}{Median of group 2}
#'   \item{est.dif}{Difference in medians (est1 - est2)}
#'
#' @details
#' The function uses a percentile bootstrap approach which is particularly robust when
#' dealing with tied values, a common issue when using medians. Missing values are
#' automatically removed.
#'
#' @seealso \code{\link{pb2gen}}, \code{\link{wmw}}, \code{\link{cid}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 5)
#' y <- rnorm(25, mean = 5.5)
#'
#' # Compare medians
#' result <- medpb2(x, y)
#' result$p.value
#' result$ci
#'
#' # Data in matrix format
#' dat <- cbind(x[1:25], y)
#' medpb2(dat)
medpb2<-function(x,y=NULL,alpha=.05,nboot=2000,SEED=TRUE){
#
#   Compare 2 independent groups using medians.
#
#   A percentile bootstrap method is used, which performs well when
#   there are tied values.
#
#   The data are assumed to be stored in x and y. If y=NULL, x is assumed to have two columns.
#
#   Missing values are automatically removed.
#
if(is.null(y)){
if(is.matrix(x) || is.data.frame(x)){
y=x[,2]
x=x[,1]
}
if(is.list(x)){
y=x[[2]]
x=x[[1]]
}
}
x=elimna(x)
y=elimna(y)
xx<-list()
xx[[1]]<-x
xx[[2]]<-y
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
est1=median(xx[[1]])
est2=median(xx[[2]])
est.dif<-median(xx[[1]])-median(xx[[2]])
crit<-alpha/2
temp<-round(crit*nboot)
icl<-temp+1
icu<-nboot-temp
bvec<-matrix(NA,nrow=2,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:2){
data<-matrix(sample(xx[[j]],size=length(xx[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,median) # Bootstrapped medians for jth group
}
top<-bvec[1,]-bvec[2,]
test<-sum(top<0)/nboot+.5*sum(top==0)/nboot
if(test > .5)test<-1-test
top<-sort(top)
ci<-NA
ci[1]<-top[icl]
ci[2]<-top[icu]
list(n1=length(x),n2=length(y),p.value=2*test,ci=ci,est1=est1,est2=est2,
est.dif=est.dif)
}

#' Cliff's Analog of the Wilcoxon-Mann-Whitney Test
#'
#' @description
#' Computes a confidence interval for P(X<Y) for two independent groups using Cliff's (1996)
#' method. This is a robust alternative to the WMW test that handles tied values and provides
#' estimates of stochastic superiority. Also reports a confidence interval for the dominance
#' measure P(X>Y) - P(X<Y).
#'
#' @inheritParams common-params
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param plotit Logical; if TRUE, creates a plot of the difference scores D = X - Y.
#' @param pop Integer (0-6) specifying plot type when plotit=TRUE: 0=adaptive kernel density,
#'   1=expected frequency curve, 2=kernel density (Rosenblatt), 3=boxplot, 4=stem-and-leaf,
#'   5=histogram, 6=kernel density estimate.
#' @param fr Argument passed to rdplot when pop=1 (default 0.8).
#' @param rval Argument passed to kdplot when pop=2 (default 15).
#' @param xlab Label for x-axis in plot.
#' @param ylab Label for y-axis in plot.
#'
#' @return A list with components:
#'   \item{n1}{Sample size of group 1}
#'   \item{n2}{Sample size of group 2}
#'   \item{cl}{Lower confidence limit for P(X>Y) - P(X<Y)}
#'   \item{cu}{Upper confidence limit for P(X>Y) - P(X<Y)}
#'   \item{d}{Estimate of P(X>Y) - P(X<Y)}
#'   \item{sqse.d}{Squared standard error of d}
#'   \item{phat}{Estimate of P(X<Y)}
#'   \item{summary.dvals}{Matrix with P(X<Y), P(X=Y), P(X>Y)}
#'   \item{ci.p}{Confidence interval for P(X<Y)}
#'
#' @details
#' The method stems from Cliff (1996, p. 140, eq 5.12) and allows tied values. The function
#' computes all pairwise differences D = X - Y and estimates P(X<Y), which equals 0.5 when
#' the distributions are identical.
#'
#' When the distribution of D is symmetric about zero, the test for symmetry can be performed
#' and provides insight into how the tails differ. Use \code{\link{cbmhd}} to compare lower
#' and upper quantiles of D directly, or \code{\link{qwmwhd}} to apply the method across
#' a range of quantiles.
#'
#' For large sample sizes (product of sample sizes > 1,000,000), use the bmp function instead.
#'
#' @references
#' Cliff, N. (1996). Ordinal methods for behavioral data analysis. Psychology Press.
#'
#' @seealso \code{\link{cidv2}}, \code{\link{cidM}}, \code{\link{cidmul}}, \code{\link{wmw}},
#'   \code{\link{cbmhd}}, \code{\link{qwmwhd}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 5)
#' y <- rnorm(25, mean = 5.5)
#'
#' # Basic usage
#' result <- cid(x, y)
#' result$phat  # P(X < Y)
#' result$ci.p  # CI for P(X < Y)
#'
#' # With plot
#' cid(x, y, plotit = TRUE, pop = 5)  # histogram of differences
cid<-function(x,y,alpha=.05,plotit=FALSE,pop=0,fr=.8,rval=15,xlab="",ylab=""){
#
# For two independent groups,
#  compute a confidence interval for  P(X<Y). The method stems from
#  Cliff, 1996, p. 140, eq 5.12. Tied values are allowed.
#
#  To compare the lower and upper quantiles of the distribution of D=X-Y,
#  use cbmhd.
#
#  This function also  reports a 1-alpha confidence interval for
#  P(X>Y)-P(X<Y)
#
#  Let xi_q be the qth quantile of the distribution of D, q<.5. To test xi_q +xi_1-q=0, use cbmhd
#  If nothing is going on, D is symmetric about zero, so this function tests for symmetry and provides
# some sense of how the tails of the distribution D differ.
#  qwmwhd applies the method using a range of q values
#
#
#  plotit=TRUE creates a plot of the difference scores.
#  pop=0 adaptive kernel density estimate
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate (Rosenblatt's shifted histogram)
#  pop=3 boxplot
#  pop=4 stem-and-leaf
#  pop=5 histogram
#  pop=6  kernel density estimate
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
if(length(x)*length(y)>10^6)stop('Use bmp with a large sample size. If using rimul, use ribmp instead')
m<-outer(x,y,FUN="-")
msave<-m
m<-sign(m)
d<-mean(m)
phat<-(1-d)/2
flag=TRUE
if(phat==0 || phat==1)flag=FALSE
q0<-sum(msave==0)/length(msave)
qxly<-sum(msave<0)/length(msave)
qxgy<-sum(msave>0)/length(msave)
c.sum<-matrix(c(qxly,q0,qxgy),nrow=1,ncol=3)
dimnames(c.sum)<-list(NULL,c("P(X<Y)","P(X=Y)","P(X>Y)"))
if(flag){
sigdih<-sum((m-d)^2)/(length(x)*length(y)-1)
di<-NA
for (i in 1:length(x))di[i]<-sum(x[i]>y)/length(y)-sum(x[i]<y)/length(y)
dh<-NA
for (i in 1:length(y))dh[i]<-sum(y[i]>x)/length(x)-sum(y[i]<x)/length(x)
sdi<-var(di)
sdh<-var(dh)
sh<-((length(y)-1)*sdi+(length(x)-1)*sdh+sigdih)/(length(x)*length(y))
zv<-qnorm(alpha/2)
cu<-(d-d^3-zv*sqrt(sh)*sqrt((1-d^2)^2+zv^2*sh))/(1-d^2+zv^2*sh)
cl<-(d-d^3+zv*sqrt(sh)*sqrt((1-d^2)^2+zv^2*sh))/(1-d^2+zv^2*sh)
}
if(!flag){
sh=NULL
nm=max(c(length(x),length(y)))
if(phat==1)bci=binomci(nm,nm,alpha=alpha)
if(phat==0)bci=binomci(0,nm,alpha=alpha)
}
if(plotit){
if(pop==1 || pop==0){
if(length(x)*length(y)>2500){
print("Product of sample sizes exceeds 2500.")
print("Execution time might be high when using pop=0 or 1")
print("If this is case, might consider changing the argument pop")
}}
if(pop==0)akerd(as.vector(msave),xlab=xlab,ylab=ylab)
if(pop==1)rdplot(as.vector(msave),fr=fr,xlab=xlab,ylab=ylab)
if(pop==2)kdplot(as.vector(msave),rval=rval,xlab=xlab,ylab=ylab)
if(pop==3)boxplot(as.vector(msave))
if(pop==4)stem(as.vector(msave))
if(pop==5)hist(as.vector(msave),xlab=xlab)
if(pop==6)skerd(as.vector(msave))
}
if(flag)pci=c((1-cu)/2,(1-cl)/2)
if(!flag){
pci=bci$ci
cl=1-2*pci[2]
cu=1-2*pci[1]
}
list(n1=length(x),n2=length(y),cl=cl,cu=cu,d=d,sqse.d=sh,phat=phat,summary.dvals=c.sum,ci.p=pci)
}

#' Cliff's Analog of WMW Test with P-Value
#'
#' @description
#' Computes a p-value for Cliff's analog of the Wilcoxon-Mann-Whitney test. This is a variation
#' of \code{\link{cid}} that focuses on hypothesis testing rather than just confidence intervals.
#'
#' @inheritParams common-params
#' @inheritParams cid
#'
#' @return A list with components:
#'   \item{n1}{Sample size of group 1}
#'   \item{n2}{Sample size of group 2}
#'   \item{p.value}{P-value for testing P(X>Y) - P(X<Y) = 0}
#'   \item{phat}{Estimate of P(X<Y)}
#'   \item{ci}{Confidence interval for P(X<Y)}
#'   \item{d.ci}{Confidence interval for P(X>Y) - P(X<Y)}
#'   \item{summary.dvals}{Matrix with P(X<Y), P(X=Y), P(X>Y)}
#'
#' @details
#' This function is similar to \code{\link{cid}} but determines the p-value by searching
#' for the alpha level at which the confidence interval for d = P(X>Y) - P(X<Y) includes
#' zero. The null hypothesis is that d = 0, which is equivalent to the distributions
#' being identical.
#'
#' To compare lower and upper quantiles of D = X - Y, use \code{\link{cbmhd}}.
#'
#' @seealso \code{\link{cid}}, \code{\link{wmw}}, \code{\link{cbmhd}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 5)
#' y <- rnorm(25, mean = 5.3)
#'
#' result <- cidv2(x, y)
#' result$p.value
#' result$phat
cidv2<-function(x,y,alpha=.05,plotit=FALSE,pop=0,fr=.8,rval=15,xlab='',ylab=''){
#
#   p-value for Cliff's analog of WMW test
#
#  To compare the lower and upper quantiles of the distribution of D=X-Y,
#  use cbmhd.
#
if(length(x)*length(y)>10^6)stop('Use bmp with a large sample size.')
nullval<-0
ci<-cid(x,y,alpha=alpha,plotit=plotit,pop=pop,fr=fr,rval=rval)
FLAG=TRUE
if(ci$phat==0 || ci$phat==1)FLAG=FALSE
if(FLAG){
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-cid(x,y,alpha=alph[i],plotit=FALSE)
if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
}
p.value<-irem/100
if(p.value<=.01){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-cid(x,y,alpha=alph[i],plotit=FALSE,xlab=xlab,ylab=ylab)
if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
}}
if(p.value<=.001){
alph<-seq(.0001,.001,.0001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-cid(x,y,alpha=alph[i],plotit=FALSE)
if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
}}
phat<-(1-ci$d)/2
pci=c((1-ci$cu)/2,(1-ci$cl)/2)
d.ci=c(ci$cl,ci$cu)
dval=cid(x,y)$summary.dvals
}
if(!FLAG){
D=bmp(x,y)
p.value=D$p.value
d.ci=NA
pci=D$ci.p
phat=D$phat
dval=ci$summary.dvals
}
list(n1=length(elimna(x)),n2=length(elimna(y)),d.hat=ci$d,d.ci=d.ci,p.value=p.value,p.hat=phat,p.ci=pci,summary.dvals=dval)
}

#' Two-Way ANOVA Multiple Comparisons with Trimmed Means (Bootstrap)
#'
#' @description
#' Performs multiple comparisons for a J-by-K factorial design using
#' trimmed means and bootstrap methods. Tests main effects and interactions.
#'
#' @inheritParams common-params
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in list mode or matrix format. If list, `x[[1]]` is level (1,1),
#'   `x[[2]]` is level (1,2), etc. Matrix format will be converted to list.
#' @param grp Subset of groups to analyze (default: all groups `1:p`).
#' @param p Total number of groups J*K (default: J*K).
#' @param tr Trimming proportion (default: 0.2 for 20% trimming).
#' @param nboot Number of bootstrap samples (default: NA, uses function default).
#' @param alpha Significance level (default: 0.05).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param pr Logical. If `TRUE`, prints progress messages (default: `TRUE`).
#' @param bhop Logical. If `TRUE`, uses Benjamini-Hochberg procedure (default: `FALSE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `Factor.A`: Multiple comparison results for Factor A main effects
#'     \item `Factor.B`: Multiple comparison results for Factor B main effects
#'     \item `Factor.AB`: Multiple comparison results for interaction effects
#'     \item `bhop`: Echo of the bhop parameter
#'     \item `SEED`: Echo of SEED setting
#'   }
#'
#' @details
#' This function performs pairwise comparisons for all effects in a two-way design
#' using trimmed means and bootstrap-t methods. The data organization follows the
#' convention where groups are ordered first by Factor A, then by Factor B.
#'
#' @seealso \code{\link{pbtrmcp}}, \code{\link{con2way}}, \code{\link{t2way}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Create 2x3 factorial design data
#' set.seed(123)
#' data <- lapply(1:6, function(i) rnorm(20, mean = i))
#'
#' # Perform two-way ANOVA with multiple comparisons
#' result <- pb2trmcp(J = 2, K = 3, x = data, nboot = 500)
#'
#' # View Factor A comparisons
#' result$Factor.A
#' }
pb2trmcp<-function(J,K,x,grp=c(1:p),p=J*K,tr=.2,nboot=NA,alpha=.05,SEED=TRUE,pr=TRUE,
bhop=FALSE){
#
#  Perform a J by K anova using trimmed means with
#  for two independent groups using a bootstrap-t method
#
#  tr=.2 is default trimming
#
#
#  The R variable data is assumed to contain the raw
#  data stored in list mode. data[[1]] contains the data
#  for the first level of both factors: level 1,1.
#  data[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  data[[K]] is the data for level 1,K
#  data[[K+1]] is the data for level 2,1, data[2K] is level 2,K, etc.
#
#  It is assumed that data has length JK, the total number of
#  groups being tested, but a subset of the data can be analyzed
#  using grp
#
if(SEED)set.seed(2)
if(is.list(x))x<-elimna(matl(x))
if(is.matrix(x))x<-elimna(x)
data<-x
if(is.matrix(data))data<-listm(data)
if(!is.list(data))stop("Data are not stored in list mode or a matrix")
if(p!=length(data)){
print("The total number of groups, based on the specified levels, is")
print(p)
print("The number of groups stored in x is")
print(length(data))
print("Warning: These two values are not equal")
}
if(p!=length(grp))stop("Apparently a subset of the groups was specified that does not match the total number of groups indicated by the values for J and K.")
temp=con2way(J,K)
conA<-temp$conA
conB<-temp$conB
conAB<-temp$conAB
if(pr)print("Taking bootstrap samples")
Factor.A<-pbtrmcp(x,con=conA,tr=tr,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
Factor.B<-pbtrmcp(x,con=conB,tr=tr,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
Factor.AB<-pbtrmcp(x,con=conAB,tr=tr,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.AB=Factor.AB,bhop=bhop,SEED=FALSE)
}

#' WMW-Based Global Test for J Independent Groups
#'
#' @description
#' Extension of the Wilcoxon-Mann-Whitney (WMW) test to J groups using a bootstrap approach.
#' Tests the global null hypothesis that all pairwise probabilities p_jk = P(X_j < X_k) equal
#' 0.5 for all j < k. This provides a robust omnibus test for comparing multiple groups.
#'
#' @inheritParams common-params
#' @param x Either a matrix, data frame, or list where each column/element represents a group.
#' @param est The estimator to use (default is \code{median}). Can also use \code{hd} for
#'   Harrell-Davis estimator when dealing with tied values.
#' @param nboot Number of bootstrap samples (default 500).
#' @param MC Logical; if TRUE, uses parallel processing via mclapply.
#' @param MM Logical; if TRUE, uses a different distance measure.
#'
#' @return The p-value for the global test.
#'
#' @details
#' The function computes P(X<Y) as an effect size measure for all pairs of groups and
#' performs a bootstrap-based global test. When tied values are detected and the default
#' median estimator is used, a warning suggests using \code{est=hd} or the function
#' \code{cidmulv2} instead.
#'
#' The test uses a projection-based distance measure to assess whether the pairwise
#' probabilities differ from the null hypothesis values.
#'
#' @seealso \code{\link{wmw}}, \code{\link{cidmul}}, \code{\link{cidmulv2}}, \code{\link{cid}}
#'
#' @export
#' @examples
#' # Three groups
#' set.seed(123)
#' g1 <- rnorm(20, mean = 5)
#' g2 <- rnorm(20, mean = 5.3)
#' g3 <- rnorm(20, mean = 5.6)
#' x <- list(g1, g2, g3)
#'
#' # Test for differences
#' wmwaov(x)
#'
#' # With Harrell-Davis estimator
#' wmwaov(x, est = hd)
wmwaov<-function(x,est=median,nboot=500,MC=FALSE,SEED=TRUE,MM=FALSE){
#
# Extension of WMW to J groups
# i.e., use p=P(X<Y) as effect size.
# test p_{jk}=.5 all j<k
#
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x=listm(x)
chk=tlist(x)
if(chk!=0){
if(identical(est,median))print("Warning: tied values detected. Suggest using est=hd or the function cidmulv2")
}
J=length(x)
L=(J^2-J)/2
ic=0
pvec=NA
boot=list()
MAT=matrix(NA,nrow=nboot,ncol=L)
for(i in 1:nboot){
for (j in 1:J){
boot[[j]]=sample(x[[j]],size=length(x[[j]]),replace=TRUE)
}
MAT[i,]=wmwloc2(boot,est=est)
}
zero=rep(0,L)
bconB=rbind(MAT,zero)
if(MC)dv=pdisMC(bconB,MM=MM)
if(!MC)dv=pdis(bconB,MM=MM)
bplus<-nboot+1
p.value<-1-sum(dv[bplus]>dv[1:nboot])/nboot-.5*sum(dv[bplus]==dv[1:nboot])/nboot
p.value
}

#' Compare Two Independent Discrete Distributions
#'
#' @description
#' Compares two independent variables in terms of their probability functions by testing
#' P(X=x) = P(Y=x) for each observed value. This method is particularly useful for highly
#' discrete data. Multiple testing is controlled using Hochberg's method by default.
#'
#' @inheritParams common-params
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param KMS Logical; if TRUE, uses the Kulinskaya, Morgenthaler and Staudte (2010)
#'   variance-stabilizing method for comparing binomial proportions and provides confidence
#'   intervals.
#' @param plotit Logical; if TRUE, creates a plot comparing relative frequencies.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis (default "Rel. Freq.").
#' @param method Multiple testing correction method passed to p.adjust (default 'hoch'
#'   for Hochberg's method).
#' @param pr Logical; if TRUE, prints results.
#'
#' @return A data frame with columns:
#'   \item{Value}{The observed value}
#'   \item{p1.est}{Proportion in group 1}
#'   \item{p2.est}{Proportion in group 2}
#'   \item{p1-p2}{Difference in proportions}
#'   \item{ci.low}{Lower confidence limit (if KMS=TRUE)}
#'   \item{ci.up}{Upper confidence limit (if KMS=TRUE)}
#'   \item{p.value}{P-value for testing P(X=value) = P(Y=value)}
#'   \item{p.adj}{Adjusted p-value using specified method}
#'
#' @details
#' For each unique value that occurs in either group, the function tests whether the
#' probability of that value is the same in both groups. When KMS=FALSE (default), uses
#' the Storer and Kim method. When KMS=TRUE, uses a variance-stabilizing transformation
#' appropriate for binomial data and provides confidence intervals.
#'
#' The plot displays relative frequencies for both groups with error bars.
#'
#' @references
#' Kulinskaya, E., Morgenthaler, S. and Staudte, R. (2010). Variance Stabilizing the
#' Difference of two Binomial Proportions. American Statistician, 64, 350-356.
#' DOI:10.1198/tast.2010.09096
#'
#' @seealso \code{\link{twobinom}}, \code{\link{twobici}}
#'
#' @export
#' @examples
#' # Discrete data example
#' set.seed(123)
#' x <- sample(1:5, 50, replace = TRUE, prob = c(0.3, 0.2, 0.2, 0.2, 0.1))
#' y <- sample(1:5, 50, replace = TRUE, prob = c(0.1, 0.2, 0.3, 0.2, 0.2))
#'
#' result <- binband(x, y)
#'
#' # Using KMS method with confidence intervals
#' binband(x, y, KMS = TRUE)
binband<-function(x,y,KMS=FALSE,alpha=.05, plotit=TRUE,xlab="X",    #ADJ.P=FALSE, old code deleted.
ylab="Rel. Freq.", method='hoch',pr=TRUE){
#
#  Comparing two independent variables in terms of their probability function.
#  For each value that occurs, say x, test P(X=x)=P(Y=x)
#  So this method is useful when dealing with highly discrete data.
#
#  If KMS=TRUE, use Kulinskaya, Morgenthaler and Staudte (2010)
#   method for comparing binomials
# Kulinskaya, E., Morgenthaler, S. and Staudte, R. (2010).
# Variance Stabilizing the Difference of two Binomial
#  Proportions. American Statistician, 64,
#  350--356 DOI:10.1198/tast.2010.09096

#  Otherwise use Storer and Kim.
#
#  method='hoch': p-values are adjusted via Hochberg's method
#
if(!KMS){
if(pr)print('To get confidence intervals, set KMS=TRUE')
}
x=elimna(x)
y=elimna(y)
vals=sort(unique(c(x,y)))
ncon=length(vals)
n1=length(x)
n2=length(y)
p.values=NA
adj=1
cv=1
if(!KMS){
output=matrix(NA,ncol=6,nrow=length(vals))
dimnames(output)=list(NULL,c("Value","p1.est","p2.est","p1-p2","p.value","p.adj"))
}
if(KMS){
output=matrix(NA,ncol=8,nrow=length(vals))
dimnames(output)=list(NULL,c("Value","p1.est","p2.est","p1-p2","ci.low","ci.up","p.value",
"p.adj"))
}
for(i in 1:length(vals)){
x1=sum(x==vals[i])
y1=sum(y==vals[i])
if(!KMS){
output[i,5]=twobinom(x1,n1,y1,n2)$p.value
output[i,2]=x1/n1
output[i,3]=y1/n2
output[i,1]=vals[i]
output[i,4]=output[i,2]-output[i,3]
}
if(KMS){
temp=bi2KMSv2(x1,n1,y1,n2)
output[i,1]=vals[i]
output[i,5]=temp$ci[1]
output[i,6]=temp$ci[2]
output[i,2]=x1/n1
output[i,3]=y1/n2
output[i,4]=output[i,2]-output[i,3]
output[i,7]=temp$p.value
}}
# Determine adjusted  critical p-value using Hochberg method
ncon=length(vals)
dvec=alpha/c(1:ncon)

if(KMS){
output[,8]=p.adjust(output[,7],method=method)
}
if(!KMS){
output[,6]=p.adjust(output[,5],method=method)
}
if(plotit)splotg2(x,y, xlab=xlab, ylab=ylab)
output
}

#' Test Equality of Regression Slopes for Two Independent Groups (Bootstrap)
#'
#' @description
#' Tests whether two independent groups have equal regression slopes using
#' a bootstrap method. Supports outlier detection and robust estimation.
#'
#' @inheritParams common-params
#' @param x1 Predictor values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Predictor values for group 2.
#' @param y2 Response values for group 2.
#' @param nboot Number of bootstrap samples (default: 599).
#' @param RAD Logical. If `TRUE`, uses robust analog of distance (default: `FALSE`).
#' @param alpha Significance level (default: 0.05).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param xout Logical. If `TRUE`, removes outliers using `outfun` (default: `FALSE`).
#' @param outfun Outlier detection function (default: \code{out}).
#'
#' @return The p-value for testing equality of slopes.
#'
#' @details
#' This function tests H0: the two groups have equal regression slopes using
#' an OLS-based bootstrap test. The test is performed by comparing the interaction
#' term in a combined regression model. Outliers in the predictor can optionally
#' be removed before analysis.
#'
#' @seealso \code{\link{olswbtest}}, \code{\link{reg2ci}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate data with different slopes
#' set.seed(123)
#' x1 <- rnorm(30)
#' y1 <- 2 * x1 + rnorm(30)
#' x2 <- rnorm(30)
#' y2 <- 0.5 * x2 + rnorm(30)
#'
#' # Test for equal slopes
#' tworegwb(x1, y1, x2, y2, nboot = 500)
#' }
tworegwb<-function(x1,y1,x2,y2,nboot=599,RAD=FALSE,alpha=.05,SEED=TRUE,xout=FALSE,
outfun=out){
#
# Simple regression (one predictor)
# Test H_0: two independent groups have equal slopes.
#
xy=elimna(cbind(x1,y1))
if(ncol(xy)>2)stop("This function only allows one covariate")
if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE)$keep
m<-m[flag,]
x1<-m[,1]
y1<-m[,2]
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE)$keep
m<-m[flag,]
x2<-m[,1]
y2<-m[,2]
}
x=c(x1,x2)
y=c(y1,y2)
g=c(rep(0,length(x1)),rep(1,length(x2)))
xgy=elimna(cbind(x,g,x*g,y))
xg=xgy[,1:3]
y=xgy[,4]
res=olswbtest(xg,y,nboot=nboot,SEED=SEED,RAD=RAD,alpha=alpha)
res[3,6]
}

#' Compare Two Discrete Distributions (Stouffer-Kolmogorov Approach)
#'
#' Compares two independent discrete variables by testing whether their
#' probability mass functions are equal: P(X=x) = P(Y=x) for all x.
#' Uses a bootstrap approach based on multinomial resampling.
#'
#' @inheritParams common-params
#' @param x Vector of discrete observations for group 1.
#' @param y Vector of discrete observations for group 2.
#' @param alpha Significance level (default: 0.05).
#' @param nboot Number of bootstrap samples (default: 500).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `test`: Test statistic (sum of squared differences in proportions)
#'     \item `p.value`: Bootstrap p-value
#'   }
#'
#' @details
#' This function provides a global test for comparing discrete distributions.
#' It tests the null hypothesis that the two groups have identical probability
#' mass functions across all observed values.
#'
#' The test statistic is the sum of squared differences between the empirical
#' probability mass functions. Bootstrap resampling under the null hypothesis
#' (using pooled proportions) provides the reference distribution.
#'
#' Note: This method appears to have no advantage over the chi-square test
#' implemented in `disc2com`. For testing the hypothesis separately at each
#' value, use `binband` instead.
#'
#' @note Requires the `mc2d` package for multinomial sampling.
#'
#' @seealso \code{\link{disc2com}}, \code{\link{binband}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two discrete distributions
#' set.seed(123)
#' x <- sample(1:5, 100, replace = TRUE, prob = c(0.2, 0.2, 0.3, 0.2, 0.1))
#' y <- sample(1:5, 100, replace = TRUE, prob = c(0.15, 0.25, 0.3, 0.2, 0.1))
#'
#' # Test for equality of distributions
#' result <- disc2comSK(x, y, nboot = 1000)
#' print(result)
#' }
disc2comSK<-function(x,y,alpha=.05,nboot=500,SEED=TRUE){
#
#  Comparing two independent variables in terms of their probability function.
#  A global test of P(X=x)=P(Y=x) for all x.
#  Appears  to have no advantage over a chi-square test done by the R function disc2com
#
#  The R function binband tests this hypothesis for each x.
#
library(mc2d)
x=elimna(x)
y=elimna(y)
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
vals=sort(unique(c(x,y)))
n1=length(x)
n2=length(y)
K=length(vals)
C1=NULL
C2=NULL
HT=NULL
for(i in 1:K){
C1[i]=sum(x==vals[i])
C2[i]=sum(y==vals[i])
HT[i]=(C1[i]+C2[i])/(n1+n2)
}
p1hat=C1/n1
p2hat=C2/n2
test=sum((p1hat-p2hat)^2)
tv=NULL
TB=NA
VP=NA
for(ib in 1:nboot){
xx=rmultinomial(n1,1,HT)
yy=rmultinomial(n2,1,HT)
B1=NA
B2=NA
BP=NA
for(i in 1:K){
B1[i]=sum(xx[,i])
B2[i]=sum(yy[,i])
}
B1hat=B1/n1
B2hat=B2/n2
TB[ib]=sum((B1hat-B2hat)^2)
}
pv=1-mean(test>TB)-.5*mean(test==TB)
list(test=test,p.value=pv)
}

#' Estimate Location for Distribution of Pairwise Differences
#'
#' @description
#' Estimates a measure of location for the distribution of all pairwise differences
#' X - Y. This is a helper function used by various WMW-based procedures.
#'
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param na.rm Logical; if TRUE, removes missing values (default TRUE).
#' @param est Estimator function to apply to the pairwise differences (default \code{median}).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return The estimated location measure for the distribution of X - Y.
#'
#' @details
#' The function computes all pairwise differences between x and y values (creating
#' an outer product), then applies the specified estimator. This is useful for
#' WMW-type analyses that don't assume independence.
#'
#' @seealso \code{\link{wmwpb}}, \code{\link{wmwloc2}}, \code{\link{cid}}
#'
#' @export
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(2, 3, 4)
#'
#' # Median of all pairwise differences
#' wmwloc(x, y)
#'
#' # Using mean instead
#' wmwloc(x, y, est = mean)
wmwloc<-function(x,y,na.rm=TRUE,est=median,...){
#
# Estimate the median of the distribution of x-y
#
if(na.rm){
x<-x[!is.na(x)]
y<-y[!is.na(y)]
}
m<-outer(x,y,FUN="-")
est=est(m,na.rm=TRUE,...)
est
}

#' Test Symmetry of Distribution Using Multiple Quantiles
#'
#' @description
#' Provides perspective on whether the distribution of D = X - Y is symmetric about zero
#' by plotting and testing the sum of q and (1-q) quantiles for multiple values of q.
#' If the distribution is symmetric, the plot should be approximately a horizontal line.
#'
#' @inheritParams common-params
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param q Numeric vector of quantiles to test (default seq(0.05, 0.40, 0.05)).
#'   All values must be less than 0.5.
#' @param xlab Label for x-axis (default "Quantile").
#' @param ylab Label for y-axis (default "Sum of q and 1-q Quantiles").
#' @param plotit Logical; if TRUE, creates a plot showing estimates and confidence intervals.
#' @param nboot Number of bootstrap samples (default 1000).
#'
#' @return A list with components:
#'   \item{n}{Vector of sample sizes c(n1, n2)}
#'   \item{output}{Data frame with columns: quantile, Est.1, Est.2, SUM, ci.low, ci.up,
#'     p_crit (critical p-value), p-value, and signif (YES/NO)}
#'
#' @details
#' For each quantile q in the argument \code{q}, the function computes a confidence
#' interval for the sum of the q-th and (1-q)-th quantiles of D = X - Y. If the
#' distribution is symmetric about zero, this sum should be zero for all q.
#'
#' Family-wise error rate (FWE) is controlled using Hochberg's method, which determines
#' critical p-values based on the specified alpha level. The plot shows estimates (marked
#' with *) and confidence interval limits (marked with +).
#'
#' This function internally calls \code{\link{cbmhd}} for each quantile.
#'
#' @seealso \code{\link{cbmhd}}, \code{\link{cid}}, \code{\link{cidv2}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 0.5)
#'
#' # Test symmetry across multiple quantiles
#' result <- qwmwhd(x, y)
#'
#' # Custom quantiles
#' qwmwhd(x, y, q = seq(0.1, 0.4, 0.1))
qwmwhd<-function(x,y,q=seq(5,40,5)/100,xlab="Quantile",ylab="Sum of q and 1-q Quantiles",plotit=TRUE,alpha=.05,nboot=1000,SEED=TRUE){
#
#  Plot that provides perspective on the degree a distribution is symmetric about zero.
#  This function plots the sum of q and 1-q quantiles of the distribution of D=X-Y, X and Y independent.
#  A 1-alpha confidence interval for the sum is indicated by a +
#  If the distribution is symmetric
#  the plot should be approximately a horizontal line.
#
#  FWE is controlled via Hochberg's method, which was used to determine critical
#  p-values based on the argument
#  alpha.
#
#  Can alter the quantiles compared via the argument
#  q
#  q must be less than .5
#
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
output=matrix(NA,ncol=8,nrow=length(q))
dimnames(output)=list(NULL,c("quantile","Est.1","Est.2","SUM","ci.low","ci.up","p_crit","p-value"))
for(i in 1:length(q)){
test=cbmhd(x,y,q=q[i],plotit=FALSE,nboot=nboot,SEED=SEED)
output[i,1]=q[i]
output[i,2]=test$Est1
output[i,3]=test$Est2
output[i,4]=test$sum
output[i,8]=test$p.value
output[i,5]=test$ci[1]
output[i,6]=test$ci[2]
}
temp=order(output[,8],decreasing=TRUE)
zvec=alpha/c(1:length(q))
output[temp,7]=zvec
output <- data.frame(output)
output$signif=rep("YES",nrow(output))
for(i in 1:nrow(output)){
if(output[temp[i],8]>output[temp[i],7])output$signif[temp[i]]="NO"
if(output[temp[i],8]<=output[temp[i],7])break
}
if(plotit){
plot(rep(q,3),c(output[,4],output[,5],output[,6]),type="n",xlab=xlab,ylab=ylab)
points(q,output[,6],pch="+")
points(q,output[,5],pch="+")
points(q,output[,4],pch="*")
}
list(n=c(n1,n2),output=output)
}

#' Yuen's Test for Dependent Groups with Effect Size
#'
#' @description
#' Performs Yuen's test for comparing trimmed means of two dependent groups and returns
#' an effect size measure similar to the one used by \code{\link{yuenv2}}. This is an
#' extension of \code{yuend} that adds effect size estimation.
#'
#' @inheritParams common-params
#' @param x Either a matrix/data frame with two columns (for paired data), or the first
#'   vector when \code{y} is specified.
#' @param y Optional second vector for paired data. If NULL, \code{x} must have two columns.
#' @param null.value Null value for the hypothesis test (default 0).
#' @param pr Logical; if TRUE, prints informational messages about effect size interpretation.
#'
#' @return A list with components:
#'   \item{test}{Test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p.value}{P-value}
#'   \item{conf.int}{Confidence interval for the difference in trimmed means}
#'   \item{estimate}{Difference in trimmed means}
#'   \item{effsize}{Effect size: (est.dif - null.value) / sd, where sd is rescaled
#'     Winsorized variance}
#'
#' @details
#' The effect size is computed as (est.dif - null.value) / sd, where sd is a Winsorized
#' variance rescaled to estimate the standard deviation under normality. This effect size
#' is similar to the one used by \code{\link{yuenv2}}.
#'
#' To get an effect size based on the difference scores directly, use \code{\link{trimciv2}}
#' instead.
#'
#' @seealso \code{\link{yuend}}, \code{\link{yuenv2}}, \code{\link{trimciv2}}, \code{\link{yuen}}
#'
#' @export
#' @examples
#' # Paired data in two columns
#' set.seed(123)
#' pre <- rnorm(30, mean = 100, sd = 15)
#' post <- pre + rnorm(30, mean = 5, sd = 10)
#' dat <- cbind(pre, post)
#'
#' # Test with effect size
#' result <- yuendv2(dat, tr = 0.2)
#' result$p.value
#' result$effsize
#'
#' # Or with separate vectors
#' yuendv2(pre, post, tr = 0.2)
yuendv2<-function(x, y, tr = 0.2, alpha = 0.05,null.value=0,pr=TRUE){
#
#  Same as yuend, only it also returns a measure of
#  effect size similar to the one used by yuenv2.
# To get a measure of effect size based on the difference scores,
#  use the function trimci or trimcipb
#  (est.dif - null.value)/sd
#  For trimmed means, sd is a Winsorized variance
#  rescaled so that it estimates the standard deviation under normality
#
if(pr)print('This version returns an effect size similar to what is used by yuenv2')
if(pr)print('To get a measure of effect size based on the difference scores, use trimciv2')
if(tr<0)stop('tr must be between 0 and .5')
if(tr>.5)stop('tr must be between 0 and .5')
res=yuend(x=x,y=y,tr=tr,alpha=alpha)
#
#if(tr==0)term=1
#if(tr>0)term=sqrt(area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr)
#epow=(res$dif-null.value)*term/sqrt(winvar(x-y,tr=tr,na.rm=TRUE))
epow=yuenv2(x,y,tr=tr)$Effect.Size
list(ci=res$ci,p.value=res$p.value,est1=res$est1,est2=res$est2,dif=res$dif,se=res$se,
teststat=res$teststat,n=res$n,df=res$df,Effect.Size=epow)
}

#' ANCOVA Using WMW Method for Two Independent Groups
#'
#' @description
#' Compares two independent groups using an ANCOVA-type method in conjunction with
#' Cliff's improvement on the Wilcoxon-Mann-Whitney test. Makes no parametric assumptions
#' about the form of regression lines; uses a running interval smoother instead.
#'
#' @inheritParams common-params
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#' @param fr1 Span for the running interval smoother for group 1 (default 1).
#' @param fr2 Span for the running interval smoother for group 2 (default 1).
#' @param sm Logical; deprecated parameter (was for bootstrap bagging).
#' @param est Estimator function for location (default \code{tmean}).
#' @param plotit Logical; if TRUE, creates plots of the regression lines.
#' @param pts Optional vector of covariate values at which to compare groups.
#' @param xout Logical; if TRUE, removes outliers in the covariate.
#' @param outfun Outlier detection function (default \code{out}).
#' @param LP Logical; if TRUE, uses running interval smoother followed by LOESS.
#' @param ... Additional arguments passed to the outlier detection function.
#'
#' @return A matrix with columns:
#'   \item{X}{Covariate value}
#'   \item{n1}{Sample size for group 1 at this X}
#'   \item{n2}{Sample size for group 2 at this X}
#'   \item{p.hat}{Estimated P(Y1 < Y2) at this X}
#'   \item{ci.low}{Lower confidence limit}
#'   \item{ci.hi}{Upper confidence limit}
#'   \item{p.value}{P-value for testing P(Y1 < Y2) = 0.5}
#'   \item{p.crit}{Critical p-value after adjustment}
#'
#' @details
#' This function compares two groups while adjusting for a covariate using a robust
#' nonparametric approach. At specified covariate values (either automatically chosen
#' or user-specified via \code{pts}), the function estimates P(Y1 < Y2) using Cliff's
#' method within neighborhoods defined by the running interval smoother.
#'
#' The function allows only one covariate. For multiple covariates, see related functions.
#'
#' @seealso \code{\link{cid}}, \code{\link{cidv2}}, \code{\link{wmw}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Example with simulated data
#' set.seed(123)
#' n <- 50
#' x1 <- rnorm(n)
#' y1 <- x1 + rnorm(n)
#' x2 <- rnorm(n)
#' y2 <- x2 + 0.5 + rnorm(n)
#'
#' result <- ancovaWMW(x1, y1, x2, y2)
#' }
ancovaWMW<-function(x1,y1,x2,y2,fr1=1,fr2=1,alpha=.05,sm=FALSE,est=tmean,
plotit=TRUE,pts=NA,xout=FALSE,outfun=out,LP=TRUE,...){
#
# Compare two independent  groups using the ancova method in conjunction
# with Cliff's improvement on the Wilcoxon-Mann-Whitney test.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#
#  OLD version: sm=TRUE will use bootstrap bagging when plotting the regression lines
#  The plot is based on measure of location indicated by the argument
#  est. Default is the Harrell-Davis estimate of the median.  Not working, took this out.
#
#   LP=TRUE: use running interval smoother followed by LOESS
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
dv.sum=NULL
if(is.na(pts[1])){
npt<-5
CC=5
isub<-c(1:5)  # Initialize isub
test<-c(1:5)
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
n1<-1
n2<-1
vecn<-1
for(i in 1:length(x1))n1[i]<-length(y1[near(x1,x1[i],fr1)])
for(i in 1:length(x1))n2[i]<-length(y2[near(x2,x1[i],fr2)])
for(i in 1:length(x1))vecn[i]<-min(n1[i],n2[i])
sub<-c(1:length(x1))
isub[1]<-min(sub[vecn>=12])
isub[5]<-max(sub[vecn>=12])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
mat<-matrix(NA,5,8)
dimnames(mat)<-list(NULL,c('X','n1','n2','p.hat','ci.low','ci.hi','p.value','p.crit'))
for (i in 1:5){
g1<-y1[near(x1,x1[isub[i]],fr1)]
g2<-y2[near(x2,x1[isub[i]],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test<-cidv2(g1,g2,alpha=alpha)
dv.sum=rbind(dv.sum,test$summary.dvals)
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
mat[i,4]<-test$p.hat
mat[i,5]<-test$p.ci[1]
mat[i,6]<-test$p.ci[2]
mat[i,7]<-test$p.value
}}
if(!is.na(pts[1])){
CC=length(pts)
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),8)
dimnames(mat)<-list(NULL,c('X','n1','n2','p.hat','ci.low','ci.hi','p.value','p.crit'))
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test=cidv2(g1,g2,alpha=alpha)
dv.sum=rbind(dv.sum,test$summary.dvals)
mat[i,1]<-pts[i]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
if(length(g1)<=5)print(paste('Warning, there are',length(g1),' points corresponding to the design point X=',pts[i]))
if(length(g2)<=5)print(paste('Warning, there are',length(g2),' points corresponding to the design point X=',pts[i]))
mat[i,4]<-test$p.hat
mat[i,5]<-test$p.ci[1]
mat[i,6]<-test$p.ci[2]
mat[i,7]<-test$p.value
}}
dvec<-alpha/c(1:CC)
temp2<-order(0-mat[,6])
mat[temp2,8]=dvec
if(plotit){
runmean2g(x1,y1,x2,y2,fr=fr1,est=est,sm=sm,xout=FALSE,LP=LP,...)
}
list(output=mat,summary=dv.sum)
}

#' Multiple Comparisons for K Independent Tests (Linear Contrasts/Trimmed Means)
#'
#' Performs a step-down multiple comparisons procedure for K independent tests
#' comparing trimmed means. Uses the Fisher method combined with Hochberg's
#' adjustment. Requires that the tests be independent.
#'
#' @param x Data in matrix or list format (optional if `x1` and `x2` provided).
#'   If matrix: 2K columns (K pairs). If list: length 2K.
#' @param x1 Data for group 1 (matrix with K columns or list of length K).
#' @param x2 Data for group 2 (matrix with K columns or list of length K).
#' @param tr Trimming proportion (default: 0.2).
#' @param alpha Significance level (default: 0.05).
#' @param pr Logical. If `TRUE`, prints results (default: `TRUE`).
#' @param opt Integer option for data format (default: 1).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `pvalues`: Vector of p-values from Yuen's test for each variable
#'     \item `adj.pvalues`: Adjusted p-values using Hochberg's method
#'     \item `num.sig`: Number of significant tests
#'     \item `test.stats`: Test statistics for each comparison
#'   }
#'
#' @details
#' This function performs Yuen's test (trimmed mean comparison) for each of K
#' independent variables, then combines the p-values using the Fisher method
#' and adjusts for multiple comparisons using Hochberg's procedure.
#'
#' The tests MUST be independent for this method to be valid. This is a
#' specialized version of `twoKgen` specifically for trimmed mean comparisons.
#'
#' @seealso \code{\link{twoKgen}}, \code{\link{yuen}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two groups on multiple independent variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 5)
#' x2 <- matrix(rnorm(120, mean = 0.3), ncol = 5)
#'
#' # Apply Yuen's test to each variable
#' result <- twoKlin(x1 = x1, x2 = x2, tr = 0.2)
#' print(result)
#' }
twoKlin<-function(x=NULL,x1=NULL,x2=NULL,tr=.2,alpha=.05,pr=TRUE,opt=1){
#
#  A step-down MCP based on K independent tests.
#  It is essential that the tests are independent.
#
#  Use Fisher method based on p-values coupled with Hochberg
#
# Data are assumed to be stored in two R variables, x1 and x2 or in one
#  R variable, x
#
# If stored in x1 and x2, they are assumed to be matrices with K columns
# or to have list mode, both having length K.
#
# If the data are stored in x,
# x is assumed to have 2K columns if a matrix or length 2K if it has list mode.
#
# If data are stored in x1 and x2, for each column, compute a p-value.
# That is, perform a test based on the data in column 1 of x1 and x2,
# followed by a test using the data in column 2 of x1 and x2, etc.
#
# If data are stored in x, the first test is based
# on the data in columns 1 and K+1,
# the second test is based on columns 2 and K+2, etc.
#
#  opt=1  Fisher's method
#  opt=2  Chen-Nadarajah method
#  opt=3  Max method
#
if(is.null(x[1])){
if(is.matrix(x1))x=cbind(x1,x2)
if(is.list(x1))x=c(x1,x2)
}
if(is.matrix(x))x=listm(x)
crit=NA
n1=NA
n2=NA
if(is.matrix(x) || is.data.frame(x))K2=ncol(x)
if(is.list(x))K2=length(x)
K=floor(K2/2)
if(2*K!=K2)stop('Total number of groups, K2, should be an even number')
ic=0
ic2=K
pv=NULL
for(i in 1:K){
ic=ic+1
ic2=ic2+1
testit=yuen(x[[ic]],x[[ic2]],tr=tr,alpha=alpha)
n1[ic]=testit$n1
n2[ic]=testit$n2
pv[ic]=testit$p.value
}
pick=NULL
v=order(pv)
ic=0
for(i in K:1){
K2=2*K
flag=TRUE
if(opt==1){
i2=i*2
if(i==K)res=(0-2)*sum(log(pv))  # Fisher test statistic
if(i<K)res=(0-2)*sum(log(pv[-pick]))  # Fisher test statistic
pvF=1-pchisq(res,i2)   #Fisher p-value based on all tests.
}
if(opt==2){
if(i==K)res=sum(qnorm(pv/2)^2)  # C-N test
if(i<K)res=sum(qnorm(pv[-pick]/2)^2)
pvF=1-pchisq(res,i)
}
if(opt==3){
if(i==K)res=max(pv)
if(i<K)res=max(pv[-pick])
pvF=pbeta(res,i,1)
}
if(pvF>alpha)flag=TRUE
if(pvF<=alpha/(K+1-i)){
ic=ic+1
pick=c(pick,v[ic])
flag=FALSE
if(pv[v[ic]]>alpha)flag=TRUE
}
if(flag)break
}
Decision=rep('Not Sig',length(pv))
if(!is.null(pick))Decision[pick]='Reject'
nsig=sum(length(pick))
list(n1=n1,n2=n2,p.values=pv,
Decisions=as.matrix(Decision),num.sig=nsig)
}

#' P-Value for Comparing Two Independent Binomial Proportions (Beal's Method)
#'
#' @description
#' Computes a p-value for comparing two independent binomial proportions
#' using Beal's method based on confidence interval inversion.
#'
#' @inheritParams common-params
#' @param r1 Number of successes in group 1 (default: `sum(x)` if x provided).
#' @param n1 Sample size for group 1 (default: `length(x)` if x provided).
#' @param r2 Number of successes in group 2 (default: `sum(y)` if y provided).
#' @param n2 Sample size for group 2 (default: `length(y)` if y provided).
#' @param x Binary data vector for group 1 (optional).
#' @param y Binary data vector for group 2 (optional).
#' @param alpha Significance level (default: 0.05).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.value`: P-value for testing equality of proportions
#'     \item `ci`: Confidence interval for difference in proportions
#'     \item `p1`: Proportion for group 1
#'     \item `p2`: Proportion for group 2
#'   }
#'
#' @details
#' This function uses Beal's method for comparing two independent binomial
#' proportions by inverting confidence intervals. The p-value is found by
#' determining the smallest alpha level at which the confidence interval
#' excludes zero.
#'
#' @seealso \code{\link{twobici}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two proportions
#' x <- rbinom(30, 1, 0.3)
#' y <- rbinom(30, 1, 0.6)
#' twobicipv(x = x, y = y)
#' }
twobicipv<-function(r1=sum(x),n1=length(x),r2=sum(y),n2=length(y),x=NA,y=NA,alpha=.05){
#
# Compute a p-value based on Beal's method for comparing two independent
# binomials.
#
alph=seq(.001,.999,.001)
for(i in 1:length(alph)){
pv=alph[i]
chk=twobici(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=y,alpha=alph[i])$ci #$
if(chk[1]>0 && chk[2]>0)break
if(chk[1]<0 && chk[2]<0)break
}
reg=twobici(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=y,alpha=alpha)
list(p.value=pv,ci=reg$ci,p1=reg$p1,p2=reg$p2)
}

#' Internal Helper for twoDcorR (Bootstrap Subsample)
#'
#' @keywords internal
twoDcorR_sub<-function(data,x,y,corfun=wincor,...){
#
# Used by TwoDcorR
#
rv=corfun(x[data,1],y[data],...)$cor
rv[2]=corfun(x[data,2],y[data],...)$cor
rv
}

#' Compare Two Robust Dependent Correlations (Overlapping)
#'
#' @description
#' Compares two robust dependent correlations in the overlapping case using bootstrap.
#' By default, uses Winsorized correlation but other robust correlation measures can
#' be specified.
#'
#' @inheritParams common-params
#' @param x Matrix or data frame with 2 columns.
#' @param y Numeric vector shared by both correlations.
#' @param corfun Correlation function to use (default \code{wincor} for Winsorized correlation).
#' @param nboot Number of bootstrap samples (default 500).
#' @param MC Logical; if TRUE, uses parallel processing via mclapply.
#' @param outfun Outlier detection function (default \code{outpro}).
#' @param ... Additional arguments passed to the correlation function.
#'
#' @return A list with components:
#'   \item{est.rho1}{Correlation between x[,1] and y}
#'   \item{est.rho2}{Correlation between x[,2] and y}
#'   \item{ci}{Bootstrap confidence interval for the difference}
#'   \item{p.value}{P-value for testing equality of correlations}
#'
#' @details
#' This function compares the correlation of x[,1] with y to the correlation of x[,2]
#' with y using a bootstrap approach. This is the "overlapping" case because both
#' correlations share the common variable y.
#'
#' The function uses robust correlation measures (Winsorized correlation by default)
#' which are less sensitive to outliers than Pearson correlation. The bootstrap
#' provides valid inference without normality assumptions.
#'
#' @seealso \code{\link{twoDNOV}}, \code{\link{TWOpov}}, \code{\link{tworhobt}},
#'   \code{\link{wincor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' n <- 50
#' x <- cbind(rnorm(n), rnorm(n))
#' y <- rnorm(n)
#'
#' # Using Winsorized correlation
#' result <- twoDcorR(x, y)
#' result$ci
#' result$p.value
twoDcorR<-function(x,y,corfun=wincor,alpha=.05,nboot=500,SEED=TRUE,MC=FALSE,outfun=outpro,...){
#
# Comparing two robust dependent correlations: Overlapping case
# Winsorized correlation is used by default.
#
# x is assumed to be a matrix with 2 columns
#
#  Compare correlation of x[,1] with y to x[,2] with y
#
#  The confidence interval is returned in ci
#  The estimates of the correlations are returned in est.rho1 and est.rho2
#
if(nrow(x)!=length(y))stop('x and y have different sample sizes; should be equal')
if(ncol(x)!=2)stop('Argument x should have two columns')
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1:2]
y=m1[,3]
est<-cor2xy(x,y,corfun=corfun,...)$cor
r12=est[1]
r13=est[2]
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
bvec<-mclapply(data,twoDcorR_sub,x,y,corfun,...)
}
if(!MC)bvec<-lapply(data,twoDcorR_sub,x,y,corfun,...)
mat=matrix(NA,nrow=nboot,ncol=2)
for(i in 1:nboot)mat[i,]=bvec[[i]]
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
bsort<-sort(mat[,1]-mat[,2])
ci12<-1
ci12[1]<-bsort[ilow]
ci12[2]<-bsort[ihi]
pv=mean(bsort<0)+.5*mean(bsort==0)
pv=2*min(c(pv,1-pv))
list(est.rho1=r12,est.rho2=r13,ci=ci12,p.value=pv)
}

#' Compare Two Robust Dependent Correlations (Non-Overlapping)
#'
#' @description
#' Compares two robust dependent correlations in the non-overlapping case using bootstrap.
#' By default, uses Winsorized correlation but other robust correlation measures can
#' be specified.
#'
#' @inheritParams common-params
#' @param x Matrix or data frame with 2 columns.
#' @param y Matrix or data frame with 2 columns.
#' @param corfun Correlation function to use (default \code{wincor} for Winsorized correlation).
#' @param nboot Number of bootstrap samples (default 500).
#' @param MC Logical; if TRUE, uses parallel processing via mclapply.
#'
#' @return A list with components:
#'   \item{est.rho1}{Correlation between x[,1] and x[,2]}
#'   \item{est.rho2}{Correlation between y[,1] and y[,2]}
#'   \item{est.dif}{Difference between the two correlations}
#'   \item{ci}{Bootstrap confidence interval for the difference}
#'   \item{p.value}{P-value for testing equality of correlations}
#'
#' @details
#' This function compares the correlation between x[,1] and x[,2] to the correlation
#' between y[,1] and y[,2]. This is the "non-overlapping" case because the two
#' correlations involve completely different pairs of variables.
#'
#' The function uses robust correlation measures (Winsorized correlation by default)
#' which are less sensitive to outliers than Pearson correlation. The bootstrap
#' provides valid inference without normality assumptions.
#'
#' @seealso \code{\link{twoDcorR}}, \code{\link{TWOpNOV}}, \code{\link{tworhobt}},
#'   \code{\link{wincor}}
#'
#' @export
#' @examples
#' set.seed(123)
#' n <- 50
#' x <- cbind(rnorm(n), rnorm(n))
#' y <- cbind(rnorm(n), rnorm(n))
#'
#' # Using Winsorized correlation
#' result <- twoDNOV(x, y)
#' result$ci
#' result$p.value
twoDNOV<-function(x,y,corfun=wincor,alpha=.05,nboot=500,SEED=TRUE,MC=FALSE){
#
# Comparing two robust dependent correlations: Non-overlapping case
# Winsorized correlation is used by default.
#
# Both x and y are assumed to be a matrix with 2 columns
#
#  Compare correlation of x[,1] x[,2] to the correlation between
#   y[,1] and  y[,2]
#
if(nrow(x)!=nrow(y))stop('x and y have different sample sizes; should be equal')
m1=cbind(x,y)
if(ncol(m1)!=4)stop('Both x and y should have two columns')
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1:2]
y=m1[,3:4]
r12=corfun(x[,1],x[,2])$cor
r13=corfun(y[,1],y[,2])$cor
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(nrow(y),size=nrow(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
bvec1<-mclapply(data,corbsub,x[,1],x[,2],corfun)
bvec2<-mclapply(data,corbsub,y[,1],y[,2],corfun)
}
if(!MC){
bvec1<-lapply(data,corbsub,x[,1],x[,2],corfun)
bvec2<-lapply(data,corbsub,y[,1],y[,2],corfun)
}
mat1=matl(bvec1)
mat2=matl(bvec2)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
bsort<-sort(mat1-mat2)
ci12<-bsort[ilow]
ci12[2]<-bsort[ihi]
ci12
pv=mean(bsort<0)
pv=2*min(c(pv,1-pv))
list(est.rho1=r12,est.rho2=r13,est.dif=r12-r13,ci=ci12,p.value=pv)
}

#' Bootstrap Confidence Interval for Location of Pairwise Differences
#'
#' @description
#' Computes a bootstrap confidence interval for a measure of location associated with
#' the distribution of all pairwise differences X - Y. This method is useful for both
#' independent and dependent data and uses the WMW-type measure of location.
#'
#' @inheritParams common-params
#' @param x Either a numeric vector (when \code{y} is provided), or a matrix/data frame
#'   with two columns.
#' @param y Optional numeric vector. If NULL, \code{x} must be a matrix/data frame.
#' @param est Estimator function to use (default \code{median}). Applied to the distribution
#'   of all pairwise differences.
#' @param nboot Number of bootstrap samples (default 2000).
#' @param pr Logical; if TRUE, prints messages.
#' @param na.rm Logical; if TRUE, removes missing values.
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list with components:
#'   \item{estimate}{Estimated location of X - Y differences}
#'   \item{ci}{Bootstrap confidence interval}
#'   \item{p.value}{P-value for testing that the location equals zero}
#'
#' @details
#' The function computes all pairwise differences between values in x and y, then
#' applies the specified estimator (default median) to this distribution. The bootstrap
#' procedure resamples from x and y separately, so the method works for both independent
#' and dependent data.
#'
#' For a non-bootstrap confidence interval, see \code{loc2dif.ci}.
#'
#' @seealso \code{\link{wmwloc}}, \code{\link{cid}}, \code{\link{wmw}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 5)
#' y <- rnorm(30, mean = 5.3)
#'
#' # Default median estimator
#' result <- wmwpb(x, y)
#' result$ci
#' result$p.value
#'
#' # Using Harrell-Davis estimator
#' wmwpb(x, y, est = hd, q = 0.5)
wmwpb<-function(x,y=NULL,est=median,alpha=.05,nboot=2000,SEED=TRUE,pr=TRUE,
na.rm=TRUE,...){
#
#   Compute a bootstrap confidence interval for a
#   measure of location associated with
#   the distribution of x-y,
#   est indicates which measure of location will be used
#   x and y are possibly dependent
#
#   loc2dif.ci  computes a non-bootstrap confidence interval
#
if(is.null(y[1])){
if(!is.matrix(x) & !is.data.frame(x))stop('With y missing, x should be a matrix')
y=x[,2]
x=x[,1]
}
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data1<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-NA
for(i in 1:nboot)bvec[i]<-wmwloc(x[data1[i,]],y[data2[i,]],est=est,na.rm=na.rm,...)
bvec<-sort(bvec)
low<-round((alpha/2)*nboot)+1
up<-nboot-low
temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
sig.level<-2*(min(temp,1-temp))
estdiff=wmwloc(x,y,est=est,na.rm=na.rm,...)
list(estimate=estdiff,ci=c(bvec[low],bvec[up]),p.value=sig.level)
}

#' Compare Two Dependent Correlations with P-Value (Overlapping Case)
#'
#' @description
#' Compares two dependent Pearson correlations in the overlapping case
#' (when both correlations share a common variable). Returns both confidence
#' interval and p-value.
#'
#' @inheritParams common-params
#' @param x Matrix with 2 columns containing the two variables to correlate with y.
#' @param y Vector to correlate with both columns of x.
#' @param alpha Significance level (default: 0.05).
#' @param CN Logical. If `TRUE`, uses Zou's confidence interval method (default: `FALSE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.value`: P-value for testing equality of correlations
#'     \item `est.rho1`: Correlation between x[,1] and y
#'     \item `est.rho2`: Correlation between x[,2] and y
#'     \item `ci`: Confidence interval for difference in correlations
#'   }
#'
#' @details
#' This function extends \code{TWOpov} by computing a p-value via confidence
#' interval inversion. Tests whether cor(x[,1], y) = cor(x[,2], y) in the
#' overlapping case where y is shared.
#'
#' @seealso \code{\link{TWOpov}}, \code{\link{TWOpNOVPV}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two correlations sharing a common variable
#' set.seed(123)
#' x <- matrix(rnorm(60), ncol = 2)
#' y <- 0.5 * x[,1] + rnorm(30)
#' TWOpovPV(x, y)
#' }
TWOpovPV<-function(x,y,alpha=.05,CN=FALSE){
#
# Comparing two dependent correlations: Overlapping case
#
# x is assumed to be a matrix with 2 columns
#
#  Compare correlation of x[,1] with y to x[,2] with y
#  returns a confidence interval stored in
#  ci
#
# This function is exactly like TWOpov, only it returns a p-value as well.
#
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-TWOpov(x,y,alpha=alph[i],CN=CN,ZCI=TRUE)$ci
if(sign(chkit[1]*chkit[2])==1)break
}
p.value<-irem/100
if(p.value<=.1){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
irem=i
p.value<-alph[i]
chkit<-TWOpov(x,y,alpha=alph[i],CN=CN,ZCI=TRUE)$ci
if(sign(chkit[1]*chkit[2])==1)break
}}
if(p.value<=.1){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-TWOpov(x,y,alpha=alph[i],CN=CN,ZCI=TRUE)$ci
if(sign(chkit[1]*chkit[2])==1)break
}}
res=TWOpov(x,y,alpha=alpha,CN=CN)
list(p.value=p.value,est.rho1=res$est.rho1,est.rho2=res$est.rho2,ci=res$ci)
}

#' Compare Two Dependent Correlations with P-Value (Non-Overlapping Case)
#'
#' @description
#' Compares two dependent Pearson correlations in the non-overlapping case
#' using HC4 heteroscedasticity correction. Returns both confidence interval
#' and p-value.
#'
#' @inheritParams common-params
#' @param x Matrix with 2 columns for first correlation.
#' @param y Matrix with 2 columns for second correlation.
#' @param HC4 Logical. If `TRUE`, uses HC4 heteroscedasticity correction (default: `TRUE`).
#' @param alpha Significance level (default: 0.05).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.value`: P-value for testing equality of correlations
#'     \item `est.rho1`: Correlation between x[,1] and x[,2]
#'     \item `est.rho2`: Correlation between y[,1] and y[,2]
#'     \item `ci.lower`: Lower bound of confidence interval
#'     \item `ci.upper`: Upper bound of confidence interval
#'   }
#'
#' @details
#' Tests whether cor(x[,1], x[,2]) = cor(y[,1], y[,2]) for dependent groups
#' in the non-overlapping case (no variables shared between the two correlations).
#' Uses HC4 method for heteroscedasticity-robust inference.
#'
#' Reference: Wilcox (2009). Comparing Pearson Correlations: Dealing with
#' Heteroscedasticity and Non-Normality. Communications in Statistics--Simulations
#' and Computations, 38, 2220-2234.
#'
#' @seealso \code{\link{TWOpNOV}}, \code{\link{TWOpovPV}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two correlations with no shared variables
#' set.seed(123)
#' x <- matrix(rnorm(60), ncol = 2)
#' y <- matrix(rnorm(60), ncol = 2)
#' TWOpNOVPV(x, y)
#' }
TWOpNOVPV<-function(x,y,HC4=TRUE,alpha=.05){
#
# Comparing two dependent correlations: Non-overlapping case
#
#   Compute a .95 confidence interval
#   for the difference between two dependent Pearson correlations,
#   non-overlapping case.
#
#    Both x and y are assumed to be matrices with two columns.
#   The function compares the correlation between x[,1] and x[,2]
#   to the correlation between y[,1] and y[,2].
#
#  For simulation results, see Wilcox (2009).
#  COMPARING PEARSON CORRELATIONS: DEALING WITH
#  HETEROSCEDASTICITY AND NON-NORMALITY, Communications in Statistics--Simulations
#   and Computations, 38, 2220-2234.
#
# This function is exactly like TWOpNOV, only it returns a p-value as well.
#
#  Note: To get a p-value, HC4=TRUE must be used.
#
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}
p.value<-irem/100
if(p.value<=.1){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
p.value<-alph[i]
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}}}
p.value<-irem/100
if(p.value<=.1){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}}
if(p.value<=.001){
alph<-seq(.0001,.001,.0001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}}
if(p.value<=.001){
alph<-seq(.0001,.001,.0001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}}
res=TWOpNOV(x,y,alpha=alpha,HC4=TRUE)
ci=c(res$ci.lower,res$ci.upper)
list(p.value=p.value,est.1=res$est.1,est.2=res$est.2,ci=ci) #ci.lower=res$ci.lower,ci.upper=res$ci.upper)
}

#' Compare Two Independent Pearson Correlations (HC4 Method)
#'
#' @description
#' Compares two independent Pearson correlations using the HC4
#' heteroscedasticity-robust method.
#'
#' @inheritParams common-params
#' @param x1 First variable for group 1.
#' @param y1 Second variable for group 1.
#' @param x2 First variable for group 2.
#' @param y2 Second variable for group 2.
#' @param alpha Significance level (default: 0.05).
#'
#' @return P-value for testing equality of correlations.
#'
#' @details
#' Tests whether the Pearson correlation in group 1 equals that in group 2
#' using HC4 standard errors for heteroscedasticity-robust inference.
#' Variables are standardized before comparison.
#'
#' @seealso \code{\link{olshc4}}, \code{\link{twopcor}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare correlations from two independent groups
#' set.seed(123)
#' x1 <- rnorm(30); y1 <- 0.5 * x1 + rnorm(30)
#' x2 <- rnorm(40); y2 <- 0.2 * x2 + rnorm(40)
#' twohc4cor(x1, y1, x2, y2)
#' }
twohc4cor<-function(x1,y1,x2,y2,alpha=.05){
#
#   Compare two independent Pearson correlations using the HC4 method
#
#
X<-elimna(cbind(x1,y1))
x1<-X[,1]
y1<-X[,2]
X<-elimna(cbind(x2,y2))
x2<-X[,1]
y2<-X[,2]
x1=(x1-mean(x1))/sd(x1)
y1=(y1-mean(y1))/sd(y1)
x2=(x2-mean(x2))/sd(x2)
y2=(y2-mean(y2))/sd(y2)
temp1=olshc4(x1,y1)
temp2=olshc4(x2,y2)
test=(temp1$ci[2,2]-temp2$ci[2,2])/sqrt(temp1$ci[2,6]^2+temp2$ci[2,6]^2)
df=length(x1)+length(x2)-4
pv=2*(1-pt(abs(test),df))
pv
}

#' Functional Data Analysis: Compare Two Groups at Multiple Time Points
#'
#' @description
#' Compares two groups using trimmed means with percentile bootstrap for
#' functional data (multiple measurements over time or space). Tests are
#' performed at specified time points with FWE control.
#'
#' @inheritParams common-params
#' @param x1 Matrix (n1-by-p) for group 1, where p is number of time points.
#' @param x2 Matrix (n2-by-p) for group 2, same number of columns as x1.
#' @param tr Trimming proportion (default: 0.2 for 20% trimming).
#' @param pts Vector of time point indices for comparisons (default: NULL for automatic selection).
#' @param npts Number of evenly-spaced points if pts=NULL (default: 25).
#' @param plotit Logical. If `TRUE`, plots results (default: `TRUE`).
#' @param alpha Significance level (default: 0.05).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param nboot Number of bootstrap samples (default: 2000).
#' @param xlab X-axis label (default: 'T').
#' @param ylab Y-axis label (default: 'Est.dif').
#' @param FBP Logical. If `TRUE`, uses family-wise error control (default: `TRUE`).
#' @param method Multiple comparison adjustment method (default: 'hochberg').
#' @param COLOR Logical. If `TRUE`, uses color in plot (default: `TRUE`).
#'
#' @return A list with test results at each time point including estimates,
#'   confidence intervals, and adjusted p-values.
#'
#' @details
#' Designed for functional data where many measurements are taken over time.
#' Uses Yuen's test with percentile bootstrap at multiple time points.
#' Controls family-wise error rate across time points.
#'
#' @seealso \code{\link{yuenbt}}, \code{\link{pb2gen}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Functional data: 30 subjects, 50 time points
#' set.seed(123)
#' x1 <- matrix(rnorm(30 * 50), nrow = 30)
#' x2 <- matrix(rnorm(30 * 50, mean = 0.5), nrow = 30)
#' funyuenpb(x1, x2, npts = 10, nboot = 500)
#' }
funyuenpb<-function(x1,x2,tr=.2,pts=NULL,npts=25,plotit=TRUE,alpha=.05,
SEED=TRUE,
nboot=2000,xlab='T',ylab='Est.dif',FBP=TRUE,method='hochberg',COLOR=TRUE){
#
#  x1 and x2 are n-by-p matrices,
#  Designed for functional data.
#  For example, p measures taken over time where  p is typically large
#
#  Goal: at speficied times, compare the two groups.
#  pts: Can specify time points where comparisons are to be made
#  if pts=NULL, pick
#  npts points evenly space between min and max time points
#
p=ncol(x1)
pm1=p-1
if(p!=ncol(x2))stop('ncol(x1) is not equal to ncol(x2)')
n1=nrow(x1)
n2=nrow(x2)
if(SEED)set.seed(2)
if(is.null(pts)){
np=round(p/npts)
if(np==0)np=1
pts=seq(2,pm1,np)
notpts=-1*length(pts)
pts=pts[c(-1,notpts)]
}
npts=length(pts)
xsub1=x1[,pts]
xsub2=x2[,pts]
res=NA
dif=NA
bvals=matrix(nrow=nboot,ncol=npts)
for(j in 1:nboot){
data1=sample(n1,size=n1,replace=TRUE)
data2<-sample(n2,size=n2,replace=TRUE)
bvals[j,]=apply(xsub1[data1,],2,tmean,tr=tr)-apply(xsub2[data2,],2,tmean,tr=tr)
}
bsort=apply(bvals,2,sort)
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
op=matrix(NA,nrow=length(pts),ncol=7)
dimnames(op)=list(NULL,c('est 1','est 2','dif','p.value',
'adjust.p.value','ci.low','ci.hi'))
op[,1]=apply(xsub1,2,tmean,tr=tr)
op[,2]=apply(xsub2,2,tmean,tr=tr)
op[,3]=op[,1]-op[,2]
bsort=apply(bvals,2,sort)
bs=bvals<0
pv=apply(bs,2,mean)
pv2=rbind(pv,1-pv)
pv2=apply(pv2,2,min)
op[,4]=2*pv2
#flag0=op[,4]==0
#op[flag0,4]=.004
op[,5]=p.adjust(op[,4],method=method)
op[,6]=bsort[icl,]
op[,7]=bsort[icu,]
if(plotit){
if(!FBP){
xlow=c(1:nrow(op))
xax=rep(c(1:nrow(op)),3)
rplot(xlow,op[,3],xlab=xlab,ylab=ylab,scat=FALSE)
plot(xax,as.vector(op[,c(3,6,7)]),type='n',xlab=xlab,ylab=ylab)
lines(xlow,op[,3])
lines(xlow,op[,6],lty=2)
lines(xlow,op[,7],lty=2)
}
if(FBP){
par(mfrow=c(1,2))
if(COLOR)FBplot(x1)
if(!COLOR)func.out(x1)
lines(medcurve(x2))
if(COLOR)FBplot(x2)
if(!COLOR)func.out(x2)
lines(medcurve(x1))
par(mfrow=c(1,1))
}}
op=cbind(pts,op)
op
}

#' Linear Contrast with Wilcoxon-Mann-Whitney Method
#'
#' Determines the distribution of a linear contrast Y_i = sum_j c_j*X_j and
#' estimates P(Y < 0) and a measure of location based on the specified function.
#'
#' @param x Data matrix, data frame, or list. Each column (for matrix/data frame)
#'   or element (for list) represents a different group.
#' @param con Vector of contrast coefficients. Must sum to zero.
#' @param locfun Function for computing the measure of location (default: `median`).
#'   Can be any function that takes a vector and returns a scalar (e.g., `mean`, `tmean`).
#' @param nreps Number of replications for determining the distribution (default: 100).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p`: Estimated probability that the linear contrast is less than zero
#'     \item `center`: Estimated measure of location of the linear contrast
#'   }
#'
#' @details
#' This function creates a linear contrast across J groups and uses a resampling
#' approach to estimate the distribution of the contrast. For each replication,
#' it samples `nmin` observations from each group (where `nmin` is the minimum
#' group size) and computes the linear contrast.
#'
#' The contrast coefficients must sum to zero (a requirement for linear contrasts).
#' The function estimates:
#' \itemize{
#'   \item The probability that the contrast is negative
#'   \item A measure of central tendency (default: median) of the contrast distribution
#' }
#'
#' This is useful for testing specific hypotheses about group differences using
#' a distribution-free approach based on the WMW method.
#'
#' @note Missing values are automatically removed from all groups.
#'
#' @seealso \code{\link{linWMWpb}} for bootstrap confidence intervals,
#'   \code{\link{wmw.ref.dif}} for WMW-based comparisons
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare three groups with contrast (-1, 0.5, 0.5)
#' set.seed(123)
#' x1 <- rnorm(30, mean = 0)
#' x2 <- rnorm(30, mean = 0.5)
#' x3 <- rnorm(30, mean = 0.5)
#' x <- list(x1, x2, x3)
#'
#' # Test if group 1 differs from average of groups 2 and 3
#' result <- linWMW(x, con = c(-1, 0.5, 0.5))
#' print(result)
#'
#' # Use mean instead of median as location measure
#' linWMW(x, con = c(-1, 0.5, 0.5), locfun = mean)
#' }
linWMW<-function(x,con,locfun=median,nreps=100,SEED=TRUE){
#
# Determine distribution of Y_i=sum_j c_jX_j
# Then estimate P(Y<0) and measure of location
# based on
# locfun, which defaults to the median.
#
con=as.vector(con)
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
if(length(con)!=J)stop('Length of con should equal number of groups')
x=elimna(x)
nv=as.vector(matl(lapply(x,FUN='length')))
nmin=min(nv)
est=NA
p=NA
B=list()
M=matrix(NA,nrow=nmin,ncol=J)
for(i in 1:nreps){
for(j in 1:J)M[,j]=sample(x[[j]],nmin)
B[[i]]=M
}
L=lapply(B,linWMWMC.sub,con=con)
est=lapply(L,locfun)
p=lapply(L,linWMWMC.sub2)
est=as.vector(matl(est))
p=as.vector(matl(p))
list(p=mean(p),center=mean(est))
}

#' Bootstrap Confidence Interval for Interaction Test Using WMW
#'
#' Computes a percentile bootstrap confidence interval for the interaction
#' effect in a 2x2 design using the WMW approach. Extends `interWMW` with
#' bootstrap inference.
#'
#' @param x Matrix with four columns or list with four elements, representing
#'   the four groups in a 2x2 design.
#' @param nreps Number of replications for WMW estimation (default: 100).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param nboot Number of bootstrap samples (default: 500).
#' @param alpha Significance level for confidence interval (default: 0.05).
#' @param nmax Maximum number of comparisons to prevent memory issues (default: 10^8).
#' @param MC Logical. If `TRUE`, uses multicore processing via `mclapply` (default: `TRUE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.est`: Point estimate of P(X1 - X2 < X3 - X4)
#'     \item `ci`: Bootstrap confidence interval for the probability
#'     \item `p.value`: P-value for testing H0: P = 0.5 (no interaction)
#'     \item `row.results`: Matrix with estimates for each row combination
#'   }
#'
#' @details
#' This function extends `interWMW` by adding bootstrap confidence intervals.
#' For each bootstrap sample, it resamples with replacement from each of the
#' four groups and computes the interaction probability.
#'
#' The confidence interval uses the percentile method from the bootstrap
#' distribution. The p-value tests whether the probability differs from 0.5,
#' which would indicate an interaction effect.
#'
#' **Interpretation**:
#' \itemize{
#'   \item CI including 0.5: No evidence of interaction
#'   \item CI excluding 0.5: Evidence of interaction
#'   \item Direction of effect shown by whether P is above or below 0.5
#' }
#'
#' When `MC = TRUE`, bootstrap samples are processed in parallel using multiple
#' cores, which can significantly speed up computation.
#'
#' @note
#' Data must have exactly four groups. Missing values are automatically removed.
#' The function requires the `parallel` package for multicore functionality.
#'
#' @seealso \code{\link{interWMW}} for point estimation without bootstrap,
#'   \code{\link{interWMWAP}} for adjusted p-values,
#'   \code{\link{linWMWpb}} for linear contrasts
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x2 factorial design with interaction
#' set.seed(123)
#' x1 <- rnorm(30, mean = 0)
#' x2 <- rnorm(30, mean = 0.5)
#' x3 <- rnorm(30, mean = 0.3)
#' x4 <- rnorm(30, mean = 1.5)  # Larger effect
#' x <- list(x1, x2, x3, x4)
#'
#' # Bootstrap test for interaction
#' result <- interWMWpb(x, nboot = 1000)
#' print(result)
#'
#' # Use more bootstrap samples with multicore
#' interWMWpb(x, nboot = 2000, MC = TRUE)
#'
#' # Single-core processing
#' interWMWpb(x, nboot = 1000, MC = FALSE)
#' }
interWMWpb<-function(x,nreps=100,SEED=TRUE,nboot=500,alpha=.05,nmax=10^8,MC=TRUE){
#
#
#
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
J<-length(x)
if(J!=4)stop('Number of groups should be four')
nv=lapply(x,length)
y=list()
pv=NA
N=max(pool.a.list(nv))
mat=matrix(NA,nrow=N,ncol=4)
for(i in 1:nboot){
for(j in 1:4)mat[1:nv[[j]],j]=sample(x[[j]],nv[[j]],replace=TRUE)
y[[i]]=mat
}
if(!MC)pv=lapply(y,interWMWpb.lsub)
if(MC)pv=mclapply(y,interWMWpb.lsub)
pv=pool.a.list(pv)
est=interWMW(x,nreps=nreps,SEED=SEED,nmax=nmax)
pv=sort(pv)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=pv[ilow]
ci[2]=pv[ihi]
pval=mean(pv<.5)+.5*mean(pv==.5)
pval=2*min(c(pval,1-pval))
list(p.est=est$p.est,ci=ci,p.value=pval,row.results=est$results.4.rows)
}

#' Internal Helper for interWMWpb (Bootstrap Subsample)
#'
#' @keywords internal
interWMWpb.lsub<-function(x,nreps=nreps){
v=interWMW(x,nreps=nreps,SEED=FALSE)$p.est
v
}

#' Bootstrap Confidence Interval for Linear Contrast with WMW Method
#'
#' Computes a percentile bootstrap confidence interval for the probability
#' that a linear contrast is less than zero, using the WMW approach.
#'
#' @param x Data matrix, data frame, or list. Each column (for matrix/data frame)
#'   or element (for list) represents a different group.
#' @param con Vector of contrast coefficients. Must sum to zero and have the
#'   same length as the number of groups.
#' @param nreps Number of replications for WMW estimation (default: 100).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param nboot Number of bootstrap samples (default: 500).
#' @param alpha Significance level for confidence interval (default: 0.05).
#' @param MC Logical. If `TRUE`, uses multicore processing via `mclapply` (default: `FALSE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.est`: Point estimate of P(contrast < 0)
#'     \item `ci`: Two-element vector with confidence interval bounds
#'     \item `p.value`: P-value for testing H0: P(contrast < 0) = 0.5
#'   }
#'
#' @details
#' This function extends `linWMW` by adding bootstrap confidence intervals.
#' For each bootstrap sample, it resamples with replacement from each group,
#' computes the linear contrast using `linWMW`, and estimates the probability
#' that the contrast is negative.
#'
#' The confidence interval is constructed using the percentile method from
#' the bootstrap distribution. The p-value tests whether the probability
#' differs from 0.5 (which would indicate no difference).
#'
#' When `MC = TRUE`, the bootstrap samples are processed in parallel using
#' multiple cores, which can significantly speed up computation for large
#' datasets or many bootstrap samples.
#'
#' @note Missing values are automatically removed from all groups before analysis.
#'
#' @seealso \code{\link{linWMW}} for the basic estimation without bootstrap,
#'   \code{\link{wmw.ref.dif}} for two-group WMW comparisons
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare three groups with contrast (-1, 0.5, 0.5)
#' set.seed(123)
#' x1 <- rnorm(30, mean = 0)
#' x2 <- rnorm(30, mean = 0.3)
#' x3 <- rnorm(30, mean = 0.3)
#' x <- list(x1, x2, x3)
#'
#' # Get bootstrap CI for P(contrast < 0)
#' result <- linWMWpb(x, con = c(-1, 0.5, 0.5), nboot = 1000)
#' print(result)
#'
#' # Use multicore processing for faster computation
#' linWMWpb(x, con = c(-1, 0.5, 0.5), nboot = 2000, MC = TRUE)
#' }
linWMWpb<-function(x,con,nreps=100,SEED=TRUE,nboot=500,alpha=.05,MC=FALSE){
#
# Compute a confidence interval for the probability that a linear contrast
# is less than zero.
#
con=as.vector(con)
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
J<-length(x)
if(J!=length(con))stop('Number of groups should be equal to the number of rows in con')
nv=lapply(x,length)
N=max(pool.a.list(nv))
mat=matrix(NA,nrow=N,ncol=J)
y=list()
pv=NA
est=linWMW(x,con=con,nreps=nreps,SEED=SEED)$p
for(i in 1:nboot){
#for(j in 1:J)y[[j]]=sample(x[[j]],nv[[j]],replace=TRUE)
for(j in 1:J)mat[1:nv[[j]],j]=sample(x[[j]],nv[[j]],replace=TRUE)
y[[i]]=mat
}
if(!MC)pv=lapply(y,linWMWpb.lsub,con=con,nreps=nreps,SEED=SEED)
if(MC){
pv=mclapply(y,linWMWpb.lsub,con=con,nreps=nreps,SEED=SEED)
}
pv=pool.a.list(pv)
pv=sort(pv)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=pv[ilow]
ci[2]=pv[ihi]
pval=mean(pv<.5)+.5*mean(pv==.5)
pval=2*min(c(pval,1-pval))
list(p.est=est,ci=ci,p.value=pval)
}

#' Internal Helper for linWMWpb (Bootstrap Subsample)
#'
#' @keywords internal
linWMWpb.lsub<-function(x,nreps=nreps,con=con,SEED=SEED){
v=linWMW(x,nreps=nreps,con=con,SEED=SEED)$p
v
}

#' Compute Delta for WMW/Cliff's Method Given Target Probability
#'
#' Determines the shift value delta such that P(X - delta < Y) = q.
#' This provides an estimate of the qth quantile of the sampling distribution
#' for the probability estimator used in Cliff's analog and WMW methods.
#'
#' @param x Vector of observations for group 1.
#' @param y Vector of observations for group 2.
#' @param q Target probability (quantile level).
#'
#' @return The shift value delta satisfying P(X - delta < Y) = q.
#'
#' @details
#' This function uses numerical optimization (via `nelder`) to find the shift
#' value that achieves a specified target probability. This is useful for
#' understanding the sampling distribution of probability estimates in
#' nonparametric methods.
#'
#' @keywords internal
#' @seealso \code{\link{cid}}, \code{\link{WMW2med.sub}}, \code{\link{nelder}}
WMW2med<-function(x,y,q){
#
# If P(X<Y)=q, determine the value delta such that
#   P(X-delta-Y<=0.0)=q
#  So an estimate of delta yields an estimate of the qth quantile
#  of the sampling distribution of  the estimator p used by Cliff's (and the WMW) method.
#
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
n=max(n1,n2)
X=matrix(NA,nrow=n+1,ncol=2)
X[1:n1,1]=x
X[1:n2,2]=y
X[n+1,]=q
v=nelder(X,1,FN=WMW2med.sub,START=0)
v
}

#' Objective Function for WMW2med Optimization (Internal)
#'
#' Internal helper function for `WMW2med` that computes the absolute difference
#' between the current probability estimate and the target quantile.
#'
#' @param X Matrix containing data and target quantile.
#' @param delta Shift parameter being optimized.
#'
#' @return Absolute difference between estimated and target probability.
#'
#' @keywords internal
#' @seealso \code{\link{WMW2med}}
WMW2med.sub<-function(X,delta){
n=nrow(X)
n1m=n-1
pv=cid(X[1:n1m,1]-delta,X[1:n1m,2])$phat
dif=abs(pv-X[n,1])
dif
}

#' Multiple Comparisons Procedure for K Independent Tests (General Function)
#'
#' Performs a step-down multiple comparisons procedure for K independent tests
#' using the Fisher method combined with Hochberg's adjustment. Requires that
#' the tests be independent.
#'
#' @param x Data in matrix or list format (optional if `x1` and `x2` provided).
#' @param x1 Data for group 1 (matrix or list, columns/elements are variables).
#' @param x2 Data for group 2 (matrix or list, columns/elements are variables).
#' @param func Test function to apply to each variable pair (default: `cidv2`).
#' @param alpha Significance level (default: 0.05).
#' @param pr Logical. If `TRUE`, prints results (default: `TRUE`).
#' @param opt Integer option for data format (default: 1).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `pvalues`: Vector of p-values for each test
#'     \item `adj.pvalues`: Adjusted p-values using Hochberg's method
#'     \item `num.sig`: Number of significant tests
#'     \item Additional components from the test function
#'   }
#'
#' @details
#' This function applies a specified test function to each of K variables
#' independently, then combines the p-values using the Fisher method and
#' adjusts for multiple comparisons using Hochberg's procedure.
#'
#' The tests MUST be independent for this method to be valid. For dependent
#' tests, other adjustment methods should be used.
#'
#' @seealso \code{\link{twoKlin}}, \code{\link{cidv2}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two groups on multiple independent variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 5)
#' x2 <- matrix(rnorm(120, mean = 0.3), ncol = 5)
#'
#' # Apply Cliff's test to each variable
#' result <- twoKgen(x1 = x1, x2 = x2, func = cidv2)
#' print(result)
#' }
twoKgen<-function(x=NULL,x1=NULL,x2=NULL,func=cidv2,alpha=.05,pr=TRUE,opt=1){
#
#  A step-down MCP based on K independent tests.
#  It is essential that the tests are independent.
#
#  Use Fisher method based on p-values coupled with Hochberg
#
# Data are assumed to be stored in two R variables, x1 and x2 or in one
#  R variable, x
#
# If stored in x1 and x2, they are assumed to be matrices with K columns
# or to have list mode, both having length K.
#
# If the data are stored in x,
# x is assumed to have 2K columns if a matrix or length 2K if it has list mode.
#
# If data are stored in x1 and x2, for each column, compute a p-value.
# That is, perform a test based on the data in column 1 of x1 and x2,
# followed by a test using the data in column 2 of x1 and x2, etc.
#
# If data are stored in x, the first test is based
# on the data in columns 1 and K+1,
# the second test is based on columns 2 and K+2, etc.
#
#  opt=1  Fisher's method
#  opt=2  Chen-Nadarajah method
#  opt=3  Max method
#
if(is.null(x[1])){
if(is.matrix(x1))x=cbind(x1,x2)
if(is.list(x1))x=c(x1,x2)
}
if(is.matrix(x))x=listm(x)
crit=NA
n1=NA
n2=NA
if(is.matrix(x) || is.data.frame(x))K2=ncol(x)
if(is.list(x))K2=length(x)
K=floor(K2/2)
if(2*K!=K2)stop('Total number of groups, K2, should be an even number')
ic=0
ic2=K
pv=NULL
for(i in 1:K){
ic=ic+1
ic2=ic2+1
testit=func(x[[ic]],x[[ic2]],alpha=alpha)
n1[ic]=testit$n1
n2[ic]=testit$n2
pv[ic]=testit$p.value
}
pick=NULL
v=order(pv)
ic=0
for(i in K:1){
K2=2*K
flag=TRUE
if(opt==1){
i2=i*2
if(i==K)res=(0-2)*sum(log(pv))  # Fisher test statistic
if(i<K)res=(0-2)*sum(log(pv[-pick]))  # Fisher test statistic
pvF=1-pchisq(res,i2)   #Fisher p-value based on all tests.
}
if(opt==2){
if(i==K)res=sum(qnorm(pv/2)^2)  # C-N test
if(i<K)res=sum(qnorm(pv[-pick]/2)^2)
pvF=1-pchisq(res,i)
}
if(opt==3){
if(i==K)res=max(pv)
if(i<K)res=max(pv[-pick])
pvF=pbeta(res,i,1)
}
if(pvF>alpha)flag=TRUE
if(pvF<=alpha/(K+1-i)){
ic=ic+1
pick=c(pick,v[ic])
flag=FALSE
if(pv[v[ic]]>alpha)flag=TRUE
}
if(flag)break
}
Decision=rep('Not Sig',length(pv))
if(!is.null(pick))Decision[pick]='Reject'
nsig=sum(length(pick))
list(n1=n1,n2=n2,p.values=pv,
Decisions=as.matrix(Decision),num.sig=nsig)
}

#' Interaction Test Using WMW Approach for 2x2 Designs
#'
#' Tests for interaction in a 2x2 ANOVA design using a Wilcoxon-Mann-Whitney
#' approach. Estimates P(X1 - X2 < X3 - X4) where groups 1-2 form one factor
#' level and groups 3-4 form another.
#'
#' @param x Matrix with four columns or list with four elements, representing
#'   the four groups in a 2x2 design.
#' @param locfun Function for computing location measure (default: `median`).
#' @param nreps Number of replications for distribution estimation (default: 200).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param nmax Maximum number of comparisons to prevent memory issues (default: 10^8).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p.est`: Estimated probability P(X1 - X2 < X3 - X4)
#'     \item `results.4.rows`: Matrix with estimates for each row combination
#'   }
#'
#' @details
#' This function provides a robust, distribution-free test for interaction in
#' a 2x2 factorial design. The approach:
#'
#' 1. Computes differences within the first pair: D1 = X1 - X2
#' 2. Computes differences within the second pair: D2 = X3 - X4
#' 3. Estimates P(D1 < D2) using a resampling approach
#'
#' **Interpretation**:
#' \itemize{
#'   \item P ≈ 0.5 suggests no interaction (similar effects across factor levels)
#'   \item P far from 0.5 suggests an interaction is present
#'   \item The direction indicates which factor level has a larger effect
#' }
#'
#' The function handles heteroscedasticity (unequal variances) and does not
#' assume normality. It uses a non-parametric estimation of the distributions
#' of the difference scores.
#'
#' For each replication, the function:
#' \itemize{
#'   \item Samples with replacement from each group to match the minimum group size
#'   \item Computes the differences D1 and D2
#'   \item Estimates the probability that D1 < D2
#' }
#'
#' The final estimate is the average across all replications.
#'
#' @note
#' The `nmax` parameter prevents memory overflow when group sizes are very large.
#' If the total number of pairwise comparisons would exceed `nmax`, a warning
#' is issued. Missing values are automatically removed.
#'
#' @seealso \code{\link{interWMWpb}} for bootstrap confidence intervals,
#'   \code{\link{interWMWAP}} for adjusted p-values,
#'   \code{\link{wmw.ref.dif}} for two-group WMW
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x2 factorial design data
#' set.seed(123)
#' # Group 1 (Factor A level 1, Factor B level 1)
#' x1 <- rnorm(30, mean = 0)
#' # Group 2 (Factor A level 1, Factor B level 2)
#' x2 <- rnorm(30, mean = 0.5)
#' # Group 3 (Factor A level 2, Factor B level 1)
#' x3 <- rnorm(30, mean = 0.3)
#' # Group 4 (Factor A level 2, Factor B level 2)
#' x4 <- rnorm(30, mean = 1.2)  # Interaction: larger effect at level 2
#'
#' x <- list(x1, x2, x3, x4)
#'
#' # Test for interaction
#' result <- interWMW(x)
#' print(result)
#'
#' # Use mean instead of median
#' interWMW(x, locfun = mean)
#' }
interWMW<-function(x,locfun=median,nreps=200,SEED=TRUE,nmax=10^8){
#
#  Goal: estimate P(X_1-X_2 < X_3-X_4).
#
# That is, dealing with an interaction in a 2-by-2 ANOVA design based on
# a Wilcoxon--Mann--Whitney approach but allow heteroscedasticity.
#
#  Strategy: estimate the distribution of X_1-X_2, non-parametrically  do the same
#  for X_3-X_4, then estimate P(X_1-X_2< X_3-X_4)
#
#  x should be a matrix with four columns or have list mode with length=4
#
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
J<-length(x)
if(J!=4)stop('Number of groups should be four')
#nx=pool.a.list(lapply(x,FUN='length'))
LL=list()
LL[[1]]=outer(x[[1]],x[[2]],FUN='-')
LL[[2]]=outer(x[[3]],x[[4]],FUN='-')
nv=c(length(LL[[1]]),length(LL[[2]]))
ntot=nv[1]*nv[2]
if(ntot<=nmax)p=bmp(LL[[1]],LL[[2]])$phat
else{
nmin=min(nv)
est=NA
p=NA
pest=NA
B=list()
M=matrix(NA,nrow=nmin,ncol=2)
for(i in 1:nreps){
for(j in 1:2)M[,j]=sample(LL[[j]],nmin)
B[[i]]=M
pest[i]=mean(M[,1]<M[,2])
}
L=lapply(B,linWMWMC.sub,con=c(1,-1))
est=lapply(L,locfun)
p=lapply(L,linWMWMC.sub2)
est=as.vector(matl(est))
p=as.vector(matl(p))
}
#
# NOTE:
# When computing a confidence interval
#L1=outer(x[[1]],x[[2]],FUN='-')
#L2=outer(x[[3]],x[[4]],FUN='-')
#bm=bmp(L1,L2)
# does not work due to the dependence among the values in L1 as well as L2
#
row=matrix(NA,nrow=2,ncol=6)
dimnames(row)=list(c('Row 1','Row 2'),c('n1','n2','p-value','p.hat','ci.low','ci.up'))
p1=cidv2(x[[1]],x[[2]])
p2=cidv2(x[[3]],x[[4]])
temp=matl(p1)
row[1,]=c(temp[1,c(1,2,5,6,7)],temp[2,7])
temp=matl(p2)
row[2,]=c(temp[1,c(1,2,5,6,7)],temp[2,7])
list(p.est=mean(p),results.4.rows=row)
}

#' WMW Interaction Effect Size for 2-by-2 Design
#'
#' @description
#' Estimates a Wilcoxon-Mann-Whitney type interaction effect size for a
#' 2-by-2 factorial design: P(Z < Z*) where Z = X1-X2 and Z* = X3-X4.
#'
#' @inheritParams common-params
#' @param x Data in matrix or list format with 4 groups.
#' @param iter Number of iterations for large sample approximation (default: 10).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#'
#' @return Estimated effect size P(Z < Z*).
#'
#' @details
#' For computational efficiency with large samples (n1*n2*n3*n4 > 1000),
#' uses a sampling approximation. Otherwise computes exact effect size.
#' The effect size measures interaction on the probability scale.
#'
#' @seealso \code{\link{WMWinterci}}, \code{\link{interWMWAP}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x2 design data
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 1), rnorm(20, 1), rnorm(20, 3))
#' WMWinter.est(x)
#' }
WMWinter.est<-function(x,iter=10,SEED=TRUE){
#
# For a 2-by-2 design,
# Estimate P(Z<Z*)
# Z = X_1-X_2
# Z* =X_3-X_4
#
#  Compute a WMW type measure of effect	size.
#
ef=NA
if(is.matrix(x))x=listm(x)
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
ef[i]=pxly(L1,L2,iter=iter,SEED=SEED)
}}
if(nt<=10^3){
L1=outer(x[[1]],x[[2]],FUN='-')
L2=outer(x[[3]],x[[4]],FUN='-')
ef=pxly(L1,L2,iter=iter,SEED=SEED)
}
ef=mean(ef)
ef
}

#' WMW Interaction Test for 2-by-2 Design (Adaptive Percentile Bootstrap)
#'
#' @description
#' Tests for interaction in a 2-by-2 factorial design using a WMW-type measure
#' with adaptive percentile bootstrap. Better than \code{WMWinterci} for
#' small unequal sample sizes with heteroscedasticity.
#'
#' @inheritParams common-params
#' @param x Data in matrix or list format with 4 groups.
#' @param nreps Number of replications for adaptive method (default: 100).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param nboot Number of bootstrap samples (default: 500).
#' @param alpha Significance level (default: 0.05).
#' @param nmax Maximum sample size product (default: 10^8).
#' @param MC Logical. If `TRUE`, uses parallel processing (default: `TRUE`).
#'
#' @return Test results including effect size estimate and confidence interval.
#'
#' @details
#' Tests interaction using P(X1-X2 < X3-X4). The adaptive percentile bootstrap
#' provides better Type I error control than \code{WMWinterci} when sample
#' sizes are unequal and heteroscedasticity is present.
#'
#' @seealso \code{\link{WMWinterci}}, \code{\link{WMWinter.est}}, \code{\link{linWMWpb}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x2 factorial design
#' set.seed(123)
#' data <- list(rnorm(15), rnorm(20, 1), rnorm(25, 1), rnorm(18, 3))
#' interWMWAP(data, nboot = 500)
#' }
interWMWAP<-function(x,nreps=100,SEED=TRUE,nboot=500,alpha=.05,nmax=10^8,MC=TRUE){
#
#  Interaction in a 2-by-2 design using P(X_1-X_2<X_3-X_4)
#  Compared to WMWinterci, this method is  better at dealing with
#  small unequal sample sizes when there is heteroscedasticity.
#
#
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
J<-length(x)
if(J!=4)stop('Number of groups should be four')
TV=linWMWpb(x,con=c(1,-1,-1,1),nreps=nreps,SEED=SEED,nboot=nboot,
alpha = alpha, MC = MC)
TV
}

#' WMW Interaction Test for 2-by-2 Design with Confidence Interval
#'
#' @description
#' Tests for interaction in a 2-by-2 factorial design using a Mann-Whitney
#' type effect measure: P(Z < Z*) where Z = X1-X2 and Z* = X3-X4.
#'
#' @inheritParams common-params
#' @param x Data in matrix or list format with 4 groups.
#' @param alpha Significance level (default: 0.05).
#' @param SW Logical. If `TRUE`, uses Studentized Wilcoxon method (default: `FALSE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `estimate`: Estimated P(Z < Z*)
#'     \item `ci`: Confidence interval
#'     \item `p.value`: P-value for testing no interaction
#'     \item `test.stat`: Test statistic
#'   }
#'
#' @details
#' Implements the Mann-Whitney type interaction effect from De Neve & Thas (2016).
#' For small unequal sample sizes with heteroscedasticity, consider using
#' \code{interWMWAP} instead.
#'
#' Reference: De Neve & Thas (2016). A Mann-Whitney type effect measure
#' of interaction for factorial designs. Communications in Statistics - Theory
#' and Methods, DOI: 10.1080/03610926.2016.1263739
#'
#' @seealso \code{\link{interWMWAP}}, \code{\link{WMWinter.est}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x2 factorial design
#' set.seed(123)
#' data <- list(rnorm(20), rnorm(20, 1), rnorm(20, 1), rnorm(20, 3))
#' WMWinterci(data)
#' }
WMWinterci<-function(x,alpha=0.05,SW=FALSE){
#
#
# Jan De Neve & Olivier Thas (2016): A Mann--Whitney type effect measure
# of interaction for factorial designs, Communications in Statistics - Theory and Methods, DOI:
#10.1080/03610926.2016.1263739
#
# For a 2-by-2 design,
# Estimate P(Z<Z^*)
# Z = X_1-X-2
# Z^* =X_3-X_4
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
if(SW)x=x[c(1,3,2,4)]
J<-length(x)
if(J!=4)stop('Number of groups should be four')
bhat=WMWinter.est(x)
gb=qnorm(bhat)
nx=pool.a.list(lapply(x,FUN='length'))
N=prod(nx)
I1=0
I2=0
I3=0
I4=0
for(i in 1:length(x[[1]]))I1=I1+Ifun(x[[1]][1],x[[2]],x[[3]],x[[4]],bhat)^2
for(i in 1:length(x[[2]]))I2=I2+Ifun(x[[1]],x[[2]][i],x[[3]],x[[4]],bhat)^2
for(i in 1:length(x[[3]]))I3=I3+Ifun(x[[1]],x[[2]],x[[3]][i],x[[4]],bhat)^2
for(i in 1:length(x[[4]]))I4=I4+Ifun(x[[1]],x[[2]],x[[3]],x[[4]][i],bhat)^2
deriv=1/dnorm(qnorm(bhat))
sighat=(deriv/N)^2*(I1+I2+I3+I4)
crit=qnorm(1-alpha/2)
ci=gb-crit*sqrt(sighat)
ci[2]=gb+crit*sqrt(sighat)
ci=pnorm(ci)
test=(gb-0.)/sqrt(sighat)
pv=2*(1-pnorm(abs(test)))
list(n=nx,p.hat=bhat,p.value=pv,ci=ci)
}

#' Yuen's Test with Quantile Shift Effect Size
#'
#' Performs Yuen's test for comparing trimmed means between two independent
#' groups, with an additional robust quantile shift (QS) measure of effect size.
#'
#' @param x Vector of data for group 1, or a matrix/data frame with two columns,
#'   or a list with two elements.
#' @param y Vector of data for group 2. Not needed if `x` is a matrix/data frame/list.
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming from each tail).
#' @param alpha Significance level for confidence interval (default: 0.05).
#' @param plotit Logical. If `TRUE`, creates a plot comparing the two groups (default: `FALSE`).
#' @param op Logical. If `TRUE`, includes output (default: `TRUE`).
#' @param cor.op Logical. For correlation-related operations (default: `FALSE`).
#' @param loc.fun Function for computing location measure (default: `median`).
#' @param pr Logical. If `TRUE`, prints interpretation guidance (default: `TRUE`).
#' @param xlab Label for x-axis when plotting (default: "X").
#' @param ylab Label for y-axis when plotting (default: " ").
#'
#' @return A list with components:
#'   \itemize{
#'     \item `ci`: Confidence interval for the difference between trimmed means
#'     \item `n1`, `n2`: Sample sizes for groups 1 and 2
#'     \item `p.value`: P-value for Yuen's test
#'     \item `dif`: Difference between trimmed means (group 1 - group 2)
#'     \item `se`: Standard error of the difference
#'     \item `teststat`: Test statistic value
#'     \item `crit`: Critical value from t-distribution
#'     \item `df`: Degrees of freedom
#'     \item `Q.Effect`: Quantile shift effect size measure
#'   }
#'
#' @details
#' This function extends the basic Yuen test (see `yuen`) by adding a robust
#' quantile shift (QS) effect size measure. The quantile shift provides an
#' alternative to Cohen's d that is more robust to outliers and non-normality.
#'
#' The function prints an interpretation guide relating Cohen's d values to
#' approximate Q.Effect values:
#' \itemize{
#'   \item Cohen's d = 0 corresponds to Q.Effect ≈ 0.5 (no effect)
#'   \item Cohen's d = 0.2 corresponds to Q.Effect ≈ 0.55 (small effect)
#'   \item Cohen's d = 0.5 corresponds to Q.Effect ≈ 0.65 (medium effect)
#'   \item Cohen's d = 0.8 corresponds to Q.Effect ≈ 0.70 (large effect)
#' }
#'
#' These correspondences assume normality and homoscedasticity. The QS effect
#' size is often more interpretable and robust than standardized mean differences.
#'
#' @note
#' Use `medpb` to compare medians (when `tr = 0.5`). Values of `tr > 0.5` are
#' not allowed. Missing values are automatically removed.
#'
#' @seealso \code{\link{yuen}} for basic Yuen test,
#'   \code{\link{yuend}} for dependent groups,
#'   \code{\link{yuend.QS.ci}} for bootstrap CI with QS effect size,
#'   \code{\link{shiftQS}} for computing quantile shift
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(40, mean = 0, sd = 1)
#' y <- rnorm(50, mean = 0.5, sd = 1.2)
#'
#' # Yuen's test with quantile shift effect size
#' result <- yuenQS(x, y)
#' print(result)
#'
#' # With plot
#' yuenQS(x, y, plotit = TRUE)
#'
#' # Less trimming
#' yuenQS(x, y, tr = 0.1)
#'
#' # Using matrix input
#' data <- cbind(x[1:40], y[1:40])
#' yuenQS(data)
#' }
yuenQS<-function(x,y=NULL,tr=.2,alpha=.05, plotit=FALSE,op=TRUE,
cor.op=FALSE,loc.fun=median,pr=TRUE,xlab='X',ylab=' ' ){
#
#  Perform Yuen's test for trimmed means on the data in x and y.
#  The default amount of trimming is 20%
#  Missing values (values stored as NA) are automatically removed.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuen$ci.
#  The significance level is returned in yuen$p.value
#
#   Unlike the function yuen, a robust quantile shift measure
#   of effect size is returned.
#
if(pr){
print('Note: Under normality and homoscedasticity, Cohen d= 0, .2, .5, .8')
print('correspond approximately to Q.Effect = 0.5, 0.55, 0.65 and 0.70, respectively')
}
if(tr==.5)stop('Use medpb to compare medians.')
if(tr>.5)stop('cannot have tr>.5')
if(is.null(y)){
if(is.matrix(x) || is.data.frame(x)){
y=x[,2]
x=x[,1]
}
if(is.list(x)){
y=x[[2]]
x=x[[1]]
}
}
x<-x[!is.na(x)]  # Remove any missing values in x
y<-y[!is.na(y)]  # Remove any missing values in y
n1=length(x)
n2=length(y)
h1<-length(x)-2*floor(tr*length(x))
h2<-length(y)-2*floor(tr*length(y))
q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
crit<-qt(1-alpha/2,df)
m1=mean(x,tr)
m2=mean(y,tr)
mbar=(m1+m2)/2
dif=m1-m2
low<-dif-crit*sqrt(q1+q2)
up<-dif+crit*sqrt(q1+q2)
test<-abs(dif/sqrt(q1+q2))
yuen<-2*(1-pt(test,df))
e.pow=shiftQS(x,y,tmean,tr=tr)$Q.Effect
if(plotit){
g2plot(x,y,xlab=xlab,ylab=ylab)
}
list(ci=c(low,up),n1=n1,n2=n2,
p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,
crit=crit,df=df,Q.Effect=e.pow)
}

#' ANCOVA Using Wilcoxon-Mann-Whitney Method with Design Points
#'
#' Performs ANCOVA comparing two groups using the Wilcoxon-Mann-Whitney analog.
#' Estimates P(Y1 < Y2 | X) at multiple covariate values with simultaneous
#' confidence bands. Uses Method S for selecting design points based on minimum
#' sample size.
#'
#' @param x1 Vector of covariate values for group 1.
#' @param y1 Vector of outcome values for group 1.
#' @param x2 Vector of covariate values for group 2.
#' @param y2 Vector of outcome values for group 2.
#' @param fr1 Span parameter for determining neighborhoods in group 1 (default: 1).
#' @param fr2 Span parameter for determining neighborhoods in group 2 (default: 1).
#' @param nmin Minimum number of observations needed in each neighborhood (default: 8).
#' @param Ycrit Logical. If `TRUE`, uses Y-based critical values (default: `FALSE`).
#' @param alpha Significance level for simultaneous confidence band (default: 0.05).
#' @param plotit Logical. If `TRUE`, creates a plot with confidence band (default: `TRUE`).
#' @param pts Vector of covariate points for evaluation. If `NA`, automatically determined
#'   (default: `NA`).
#' @param span Span parameter for loess smoothing when plotting (default: 2/3).
#' @param sm Logical. If `TRUE`, smooths the plot using lowess (default: `TRUE`).
#' @param BOTH Logical. If `TRUE`, considers both groups when determining range
#'   (default: `TRUE`).
#' @param pr Logical. If `TRUE`, prints additional information (default: `TRUE`).
#' @param xout Logical. If `TRUE`, removes outliers from covariates (default: `FALSE`).
#' @param outfun Function for detecting outliers (default: `out`).
#' @param MC Logical. If `TRUE`, uses multicore processing (default: `FALSE`).
#' @param npts Number of covariate points to evaluate (default: 25).
#' @param p.crit Critical p-value for simultaneous testing. If `NULL`, determined
#'   automatically (default: `NULL`).
#' @param EST Logical. If `TRUE`, estimates critical values (default: `FALSE`).
#' @param SCAT Logical. If `TRUE`, includes scatter plot (default: `TRUE`).
#' @param xlab Label for x-axis (default: "X").
#' @param ylab Label for y-axis (default: "P.hat").
#' @param pc Plotting character (default: ".").
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `output`: Matrix with columns for covariate values, estimates, and confidence bounds
#'     \item `p.crit`: Critical p-value used for simultaneous inference
#'     \item Other components from `ancovaWMW`
#'   }
#'
#' @details
#' This function extends ANCOVA to use a robust WMW-based approach. At each
#' covariate value, it estimates P(Y1 < Y2) - the probability that an observation
#' from group 1 is less than an observation from group 2, given the covariate value.
#'
#' **Method S** (used by this function) chooses covariate points based on having
#' at least `nmin` observations in the neighborhood. This differs from **Method Q**
#' (see `ancdetwmwQ`) which selects points between quantiles.
#'
#' The function creates a simultaneous confidence band across all covariate values,
#' controlling the familywise error rate at level `alpha`. Critical p-values are
#' either specified via `p.crit` or determined automatically based on simulation
#' results for different sample sizes.
#'
#' When `plotit = TRUE`, the function creates a plot showing:
#' \itemize{
#'   \item Point estimates of P(Y1 < Y2) across covariate values
#'   \item Simultaneous confidence band
#'   \item Optional smoothing via lowess
#' }
#'
#' @note
#' This function only allows one covariate. For multiple covariates, use other
#' ANCOVA methods. Missing values are automatically removed.
#'
#' @seealso \code{\link{ancdetwmwQ}} for Method Q (quantile-based points),
#'   \code{\link{ancovaWMW}} for the underlying WMW ANCOVA computation,
#'   \code{\link{ancdet}} for ANCOVA with trimmed means
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 50
#' x1 <- rnorm(n)
#' y1 <- 2 + 0.5*x1 + rnorm(n)
#' x2 <- rnorm(n)
#' y2 <- 3 + 0.5*x2 + rnorm(n)
#'
#' # ANCOVA with WMW method
#' result <- ancdetwmw(x1, y1, x2, y2)
#'
#' # With fewer design points
#' ancdetwmw(x1, y1, x2, y2, npts = 15)
#'
#' # Remove outliers in covariates
#' ancdetwmw(x1, y1, x2, y2, xout = TRUE)
#' }
ancdetwmw<-function(x1,y1,x2,y2,fr1=1,fr2=1,nmin=8,Ycrit=FALSE,
alpha=.05,plotit=TRUE,pts=NA,span=2/3,sm=TRUE,BOTH=TRUE,
pr=TRUE,xout=FALSE,outfun=out,MC=FALSE,
npts=25,p.crit=NULL,EST=FALSE,
SCAT=TRUE,xlab='X',ylab='P.hat',pc='.',...){
#
#  Like the function ancdet, only use analog of Wilcoxon--Mann--Whitney
#  plot=TRUE:  plot  estimates P.hat plus a
# confidence band having simultaneous probability coverage 1-alpha
#
#   Method S: choose covariate points based on nmin.
#  In contrast method Q, performed by ancdetwmwQ, picks points between the q and 1-q quantiles.
#
#
#  span = the span when using loess to plot the regression line.
#
# npts = number of  covariate values to be used
#
# sm=TRUE will smooth the plot using lowess
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
xor1=order(x1)
xor2=order(x2)
x1=x1[xor1]
x2=x2[xor2]
y1=y1[xor1]
y2=y2[xor2]
n1<-1
n2<-1
vecn<-1
for(i in 1:length(x1))n1[i]<-length(y1[near(x1,x1[i],fr1)])
for(i in 1:length(x1))n2[i]<-length(y2[near(x2,x1[i],fr2)])
for(i in 1:length(x1))vecn[i]<-min(n1[i],n2[i])
sub<-c(1:length(x1))
bot<-min(sub[vecn>=nmin])
itop<-max(sub[vecn>=nmin])
xbot=x1[bot]
xup=x1[itop]

if(BOTH){
vecn=1
n1=1
n2=1
for(i in 1:length(x2))n1[i]<-length(y1[near(x1,x1[i],fr1)])
for(i in 1:length(x2))n2[i]<-length(y2[near(x2,x1[i],fr2)])
for(i in 1:length(x2))vecn[i]<-min(n1[i],n2[i])
sub<-c(1:length(x2))
bot<-max(sub[vecn>=nmin])
itop<-min(sub[vecn>=nmin])
xbot[2]=x2[itop]  #CORRECT, need to switch
xup[2]=x2[bot]
}
xbot=max(xbot)
xup=min(xup)
pts=seq(xbot,xup,length.out=npts)
if(alpha!=.05)EST=TRUE
if(is.null(p.crit)){
nv=c(30,  50,  60,  70,  80, 100, 150, 200,
300, 400, 500, 600, 800)
if(Ycrit)pv=c(0.00824497,0.00581, 0.005435089, 0.004763079,
0.00416832, 0.004406774, 0.00388228,0.003812836,0.003812836,0.003453055, 0.003625061,
.003372966, 0.003350022)
if(!Ycrit)
pv=c(0.008566, # 30
0.0083847, # 50
0.006758,  # 60
0.006871,   # 70
0.006157,  # 80
0.006629, #100
0.006629, #  150
0.004681, # 200
0.004537,  # 300
0.004952, # 400
 0.004294,  # 500
 0.004288,  # 600
 0.004148)
n1= length(y1)
 n2=length(y2)
p.crit=(lplot.pred(1/nv,pv,1/n1)$yhat+lplot.pred(1/nv,pv,1/n2)$yhat)/2
p.crit=(alpha/.05)*p.crit # Crude approximation when alpha != .05, tends to be conservative.
}
temp=ancovaWMW(x1,y1,x2,y2,pts=pts,fr1=fr1,fr2=fr2,alpha=p.crit,plotit=plotit)
res=temp$output[,1:7]
if(plotit){
x=res[,1]
y=res[,4]
minx=min(x)
maxx=max(x)
plot(c(minx,maxx,x),c(0,1,y),xlab=xlab,ylab=ylab,type='n')
points(x,y,pch=pc)
if(!sm){lines(res[,1],res[,5],lty=2)
lines(res[,1],res[,6],lty=2)
lines(res[,1],res[,4])
}
if(sm){
plin=lplot.pred(res[,1],res[,4],span=span)$yhat
lines(res[,1],plin)
low.line=lplot.pred(res[,1],res[,5],span=span)$yhat
lines(res[,1],low.line,lty=2)
up.line=lplot.pred(res[,1],res[,6],span=span)$yhat
lines(res[,1],up.line,lty=2)
}

}
sig=rep(0,nrow(res))
sig[res[,7]<=p.crit]=1
sig=as.matrix(sig,ncol=1)
dimnames(sig)=list(NULL,'Sig.Dif')
res=cbind(res,sig)
list(p.crit=p.crit,output=res,summary=temp$summary,num.sig=sum(sig))
}

#' ANCOVA Using WMW at Design Points with Simultaneous Confidence Band
#'
#' @description
#' Performs ANCOVA using a Wilcoxon-Mann-Whitney analog at multiple design points
#' with a simultaneous confidence band. Like \code{ancdet} but uses WMW instead
#' of means.
#'
#' @inheritParams common-params
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#' @param fr1 Span for group 1 (default: 1).
#' @param fr2 Span for group 2 (default: 1).
#' @param nmin Minimum sample size at each design point (default: 8).
#' @param q Quantile for selecting covariate range (default: 0.05).
#' @param alpha Significance level (default: 0.05).
#' @param plotit Logical. If `TRUE`, creates plot (default: `TRUE`).
#' @param pts Design points (default: NA for automatic selection).
#' @param span Span for lowess smoothing (default: 2/3).
#' @param sm Logical. If `TRUE`, applies smoothing (default: `TRUE`).
#' @param xout Logical. If `TRUE`, removes outliers (default: `FALSE`).
#' @param outfun Outlier detection function (default: \code{out}).
#' @param MC Logical. If `TRUE`, uses parallel processing (default: `FALSE`).
#' @param npts Number of design points if pts not specified (default: 25).
#' @param p.crit Critical p-value for confidence band (default: NULL for automatic).
#' @param SCAT Logical. If `TRUE`, shows scatter plot (default: `TRUE`).
#' @param xlab X-axis label (default: 'X').
#' @param ylab Y-axis label (default: 'P.hat').
#' @param pc Plot character (default: '.').
#' @param ... Additional arguments passed to outlier function.
#'
#' @return A list with WMW estimates and confidence bands at design points.
#'
#' @details
#' Estimates P(Y1 < Y2 | X = x) at multiple design points with simultaneous
#' confidence band. Design points chosen between q and 1-q quantiles of
#' the covariate, requiring at least nmin observations at each point.
#'
#' @seealso \code{\link{ancovaWMW}}, \code{\link{ancdet}}
#'
#' @export
#' @examples
#' \dontrun{
#' # ANCOVA with WMW
#' set.seed(123)
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- x2 + 0.5 + rnorm(50)
#' ancdetwmwQ(x1, y1, x2, y2, npts = 15)
#' }
ancdetwmwQ<-function(x1,y1,x2,y2,fr1=1,fr2=1,nmin=8,q=.05,
alpha=.05,plotit=TRUE,pts=NA,span=2/3,sm=TRUE, xout=FALSE,outfun=out,MC=FALSE,
npts=25,p.crit=NULL,
SCAT=TRUE,xlab='X',ylab='P.hat',pc='.',...){
#
#  Like the function ancdet, only use analog of Wilcoxon--Mann--Whitney
#  plot=TRUE:  plot  estimates P.hat plus a
# confidence band having simultaneous probability coverage 1-alpha
#
#  span = the span when using loess to plot the regression line.
#
# npts = number of  covariate values to be used
#
# sm=TRUE will smooth the plot using lowess
#
#  Covariate points are chosen that lie between the q and 1-q quantiles
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
xor1=order(x1)
xor2=order(x2)
x1=x1[xor1]
x2=x2[xor2]
y1=y1[xor1]
y2=y2[xor2]
n1<-1
n2<-1
vecn<-1

xbot=max(qest(x1,q),qest(x2,q))
xup=min(qest(x1,1-q),qest(x2,1-q))

pts=seq(xbot,xup,length.out=npts)

nchk1=0
for(i in 1:length(pts))nchk1[i]=length(y1[near(x1,pts[i],fr1)])
nchk2=0
for(i in 1:length(pts))nchk2[i]=length(y2[near(x2,pts[i],fr2)])
flag1=nchk1>=nmin
flag2=nchk2>=nmin
flag=as.logical(flag1*flag2)
pts=pts[flag]
if(is.null(p.crit)){
nv=c(30,  50,  60,  70,  80, 100, 150, 200,
300, 400, 500, 600, 800)
pv=c(0.00824497,0.00581, 0.005435089, 0.004763079,
0.00416832, 0.004406774, 0.00388228,0.003812836,0.003812836,0.003453055, 0.003625061,
.003372966, 0.003350022)
n1= length(y1)
 n2=length(y2)
p.crit=(lplot.pred(1/nv,pv,1/n1)$yhat+lplot.pred(1/nv,pv,1/n2)$yhat)/2
p.crit=(alpha/.05)*p.crit # Crude approximation when alpha != .05, tends to be conservative.
}
temp=ancovaWMW(x1,y1,x2,y2,pts=pts,fr1=fr1,fr2=fr2,alpha=p.crit,plotit=plotit)
res=temp$output
if(plotit){
x=res[,1]
y=res[,4]
minx=min(x)
maxx=max(x)
plot(c(minx,maxx,x),c(0,1,y),xlab=xlab,ylab=ylab,type='n')
points(x,y,pch=pc)
if(!sm){lines(res[,1],res[,5],lty=2)
lines(res[,1],res[,6],lty=2)
lines(res[,1],res[,4])
}
if(sm){
plin=lplot.pred(res[,1],res[,4],span=span)$yhat
lines(res[,1],plin)
low.line=lplot.pred(res[,1],res[,5],span=span)$yhat
lines(res[,1],low.line,lty=2)
up.line=lplot.pred(res[,1],res[,6],span=span)$yhat
lines(res[,1],up.line,lty=2)
}

}
sig=rep(0,nrow(res))
sig[res[,7]<=p.crit]=1
sig=as.matrix(sig,ncol=1)
dimnames(sig)=list(NULL,'Sig.Dif')
res=cbind(res,sig)
list(p.crit=p.crit,output=res,summary=temp$summary,num.sig=sum(sig),p.crit=p.crit)
}

#' Multiple Cliff's Analog Tests for Multivariate Data
#'
#' For two independent p-variate random variables, performs Cliff's analog
#' of the Wilcoxon-Mann-Whitney test for each variable separately.
#'
#' @param x1 Matrix or data frame for group 1 (rows are observations, columns are variables).
#' @param x2 Matrix or data frame for group 2 (rows are observations, columns are variables).
#' @param alpha Significance level (default: 0.05).
#' @param BMP Logical. If `TRUE`, uses `bmp` function; if `FALSE`, uses `cidv2` (default: `FALSE`).
#'
#' @return A matrix with one row per variable and columns:
#'   \itemize{
#'     \item `n1`, `n2`: Sample sizes for groups 1 and 2
#'     \item `p-value`: P-value for the comparison
#'     \item `p.hat`: Estimate of P(X < Y)
#'     \item `ci.low`, `ci.up`: Confidence interval bounds
#'     \item `P(X<Y)`, `P(X=Y)`, `P(X>Y)`: Summary of pairwise comparisons
#'   }
#'
#' @details
#' This function applies Cliff's analog test to each variable in a multivariate
#' dataset. It provides a robust alternative to multivariate t-tests when comparing
#' distributions rather than just means.
#'
#' The function does NOT adjust for multiple comparisons; users should apply
#' appropriate corrections (e.g., Bonferroni, FDR) if needed.
#'
#' @seealso \code{\link{cidv2}}, \code{\link{bmp}}, \code{\link{cid}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two groups on multiple variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 5)
#' x2 <- matrix(rnorm(120, mean = 0.3), ncol = 5)
#'
#' # Apply Cliff's test to each variable
#' result <- cidMULT(x1, x2)
#' print(result)
#' }
cidMULT<-function(x1,x2,alpha=.05,BMP=FALSE){
#
#
#  Deals with  two independent p-variate random variables
# for each variable, compare the two groups with Cliff's analog of the
# Wilcoxon--Mann--Whitney test.
#
p=ncol(x1)
V=list()
res=matrix(NA,nrow=p,ncol=9)
for(j in 1:p){
if(!BMP)a=cidv2(x1[,j],x2[,j])
if(BMP)a=bmp(x1[,j],x2[,j])
res[j,1]=a$n1
res[j,2]=a$n2
res[j,3]=a$p.value
if(!BMP)res[j,4]=a$p.hat
if(BMP)res[j,4]=a$phat
if(!BMP)res[j,5:6]=a$p.ci
#if(BMP)res[j,5:6]=a$ci.p
res[j,7:9]=a$summary.dvals
}
L=NULL
for(j in 1:p)L[j]=paste('Var',j)
dimnames(res)=list(L,c('n1','n2','p-value','p.hat','ci.low','ci.up','P(X<Y)','P(X=Y)','P(X>Y)'))
res
}

#' Pool Data for Two-Way Design (Internal Helper)
#'
#' For a two-way design, pools data for each level of Factor A across levels
#' of Factor B, and vice versa. Internal helper for two-way effect size calculations.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in matrix or list format (J*K groups).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `A`: List of length J, each element containing pooled data across Factor B
#'     \item `B`: List of length K, each element containing pooled data across Factor A
#'   }
#'
#' @keywords internal
#' @seealso \code{\link{twowayESM}}, \code{\link{pool.a.list}}
twoway.pool<-function(J,K,x){
#
#  For a two-way design,for each level of Factor A, pool over B
#  Do the same for Factor B
#
#  Return results in A and B having list mode,
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
JK=J*K
imat=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
A=list()
B=list()
for(j in 1:J)A[[j]]=pool.a.list(x[imat[j,]])
for(k in 1:K)B[[k]]=pool.a.list(x[imat[,k]])
list(A=A,B=B)
}

#' Two-Way Effect Size Measures
#'
#' For a two-way design, computes effect size measures for each level of Factor A
#' (pooling over Factor B) and for each level of Factor B (pooling over Factor A).
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in matrix or list format (J*K groups).
#' @param fun Function to compute effect sizes (default: `ES.summary`).
#' @param ... Additional arguments passed to `fun`.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `FactorA`: Effect size results for Factor A (pooled over B)
#'     \item `FactorB`: Effect size results for Factor B (pooled over A)
#'   }
#'
#' @details
#' This function pools data appropriately for a two-way design and computes
#' effect size measures using the specified function (typically `ES.summary`
#' or related functions from the effect size module).
#'
#' @seealso \code{\link{twoway.pool}}, \code{\link{ES.summary}}, \code{\link{IND.PAIR.ES}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Two-way design: 2 x 3
#' set.seed(123)
#' J <- 2; K <- 3
#' x <- matrix(rnorm(120), ncol = J*K)
#'
#' # Compute effect sizes for main effects
#' result <- twowayESM(J, K, x)
#' print(result)
#' }
twowayESM<-function(J,K,x,fun=ES.summary,...){
#
# For each level of Factor A, pool data over levels of Factor B, compute
# measures of effect size.
#
# Do the same for Factor B
#
#  argument fun, see the function IND.PAIR.ES
#
a=twoway.pool(J,K,x)
A=IND.PAIR.ES(a$A,fun=fun,...)
B=IND.PAIR.ES(a$B,fun=fun,...)
list(Factor.A=A,Factor.B=B)
}

#' Compare Variation in Tails for Dependent Groups After Centering
#'
#' @description
#' Compares the marginal distributions of two dependent groups by examining
#' variation in the tails after centering the data using a location measure.
#'
#' @inheritParams common-params
#' @param x Observations for group 1.
#' @param y Observations for group 2.
#' @param loc.fun Location function for centering (default: \code{median}).
#' @param CI Logical. If `TRUE`, computes confidence intervals (default: `FALSE`).
#' @param plotit Logical. If `TRUE`, creates plot (default: `TRUE`).
#' @param xlab X-axis label (default: 'First Group').
#' @param ylab Y-axis label (default: 'Est.1 - Est.2').
#' @param ylabQCOM Y-axis label for QCOM (default: 'Est.2 - Est.1').
#' @param sm Logical. If `TRUE`, applies smoothing (default: `TRUE`).
#' @param QCOM Logical. If `TRUE`, uses quantile comparison method (default: `TRUE`).
#' @param q Quantiles to compare (default: c(0.1, 0.25, 0.75, 0.9)).
#' @param MC Logical. If `TRUE`, uses parallel processing (default: `FALSE`).
#' @param nboot Number of bootstrap samples (default: 2000).
#' @param PR Logical. If `TRUE`, prints interpretation guide (default: `TRUE`).
#' @param ... Additional arguments passed to location function.
#'
#' @return Results comparing tail variation between groups.
#'
#' @details
#' Centers both groups using the specified location function, then compares
#' the distributions focusing on tail behavior. Useful for detecting differences
#' in variability patterns after removing location differences.
#'
#' @seealso \code{\link{Dqcomhd}}, \code{\link{lband}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare tail variation in paired data
#' set.seed(123)
#' x <- rt(50, df = 5)
#' y <- rt(50, df = 3) # Heavier tails
#' RMcomvar.locdis(x, y)
#' }
RMcomvar.locdis<-function(x,y,
loc.fun=median,CI=FALSE,plotit=TRUE,xlab='First Group',
ylab='Est.1 - Est.2',ylabQCOM='Est.2 - Est.1',sm=TRUE,QCOM=TRUE,q=c(.1,.25,.75,.9),MC=FALSE,nboot=2000,PR=TRUE,...){
#
#  Compare the marginal distributions of two dependent groups in terms of the
#  variation in the tails using all of the quantiles
#  after centering the data.
#
#  CI=FALSE, suppresses confidence intervals
#
if(!QCOM){
if(PR){
print('Interpretation: when using QCOM=F:  If values in  q.sig.greater are less than .5')
print('this indicates more variation in the lower tail for group 1')
print('Interpretation: If values in  q.sig.greater are greater than .5')
print('This indicates more variation in the lower tail for group 2')

print('Interpretation: If values in  q.sig.less are less than .5')
print('this indicates more variation in the upper tail for group 2')
print('Interpretation: If values in  q.sig.less are greater than .5')
print('This indicates more variation in the upper tail for group 1')
}
}
x=elimna(x)
y=elimna(y)
mx=loc.fun(x,...)
my=loc.fun(y,...)
X=x-mx
Y=y-my
if(!QCOM){
a=lband(X,Y,plotit=plotit,xlab=xlab,ylab=ylabQCOM,sm=sm,CI=CI)
if(!CI)a$m=NULL
}
else{
a=Dqcomhd(X,Y,q=q,nboot=nboot,plotit=plotit,xlab=xlab,ylab=ylab)
}
a
}

#' Confidence Interval for Quantile Shift Effect Size (Dependent Groups)
#'
#' Computes a bootstrap confidence interval for the quantile shift measure of
#' effect size for dependent (paired) groups. The quantile shift is a robust
#' measure of how much one distribution is shifted relative to another.
#'
#' @inheritParams common-params
#' @param x Either a vector of observations for group 1, or a two-column matrix/data frame
#'   where the first column is group 1 and the second is group 2.
#' @param y Vector of observations for group 2 (optional if `x` is a matrix).
#' @param tr Trimming proportion (default: 0.2 for 20% trimming).
#' @param alpha Significance level for confidence interval (default: 0.05 for 95% CI).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param ... Additional arguments passed to `tmean`.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `Q.effect`: The quantile shift effect size estimate
#'     \item `ci`: Two-element vector with lower and upper confidence interval bounds
#'   }
#'
#' @details
#' This function uses the percentile bootstrap method to construct a confidence
#' interval for the quantile shift measure. The quantile shift (Q) is computed
#' using the `shiftes` function with trimmed means as the location estimator.
#'
#' The bootstrap resamples pairs of observations (preserving the dependency
#' structure) and computes the quantile shift for each bootstrap sample. The
#' confidence interval is then obtained from the percentiles of the bootstrap
#' distribution.
#'
#' @seealso \code{\link{shiftes}}, \code{\link{yuend}}, \code{\link{tmean}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate paired data with a shift
#' set.seed(123)
#' x <- rnorm(30)
#' y <- x + rnorm(30, mean = 0.5)
#'
#' # Compute quantile shift effect size with CI
#' result <- yuend.QS.ci(x, y, nboot = 1000)
#' print(result)
#' }
yuend.QS.ci<-function(x,y=NULL,tr=.2,alpha=.05,nboot=1000,SEED=TRUE,...){
#
# confidence interval for the quantile shift measure of effect size.
#
if(is.null(y)){
if(ncol(x)!=2)stop('If y is null, x should have	two columns')
}
if(SEED)set.seed(2)
xy=cbind(x,y)
xy=elimna(xy)
n1=nrow(xy)
v=NA
for(i in 1:nboot){
id=sample(n1,replace=TRUE)
v[i]=shiftes(xy[id,1],xy[id,2],locfun=tmean,tr=tr)$Q.Effect
}
v=sort(v)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=v[ilow]
ci[2]=v[ihi]
est=shiftes(xy[,1],xy[,2],locfun=tmean,tr=tr)$Q.Effect
list(Q.effect=est,ci=ci)
}

#' Multiple Comparisons of Variances for Independent Groups
#'
#' Performs all pairwise comparisons of variances for J independent groups
#' using a slight extension of the HC4 version of the Morgan-Pitman test.
#' Controls the familywise error rate using a specified p-value adjustment method.
#'
#' @param x Either a matrix/data frame (each column is a group) or a list of vectors.
#' @param method Method for p-value adjustment (default: "hoch" for Hochberg's method).
#'   Other options include "holm", "bonferroni", "BH", etc. See \code{\link[stats]{p.adjust}}.
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#'
#' @return A matrix with one row per comparison and columns:
#'   \itemize{
#'     \item `Var`: Index of first variable
#'     \item `Var`: Index of second variable
#'     \item `SD 1`: Standard deviation of first group
#'     \item `SD 2`: Standard deviation of second group
#'     \item `Dif`: Difference in standard deviations (SD1 - SD2)
#'     \item `p.value`: Unadjusted p-value
#'     \item `Adj.p.value`: Adjusted p-value using specified method
#'   }
#'
#' @details
#' This function compares the variances of J independent variables using all
#' pairwise comparisons. The Morgan-Pitman test (via `varcom.IND.MP`) is used
#' for each comparison, which is a robust heteroscedastic test for comparing
#' two variances.
#'
#' The function performs C = J(J-1)/2 pairwise comparisons and adjusts the
#' p-values to control the familywise error rate.
#'
#' @seealso \code{\link{varcom.IND.MP}}, \code{\link{comvar2}}, \code{\link[stats]{p.adjust}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare variances of three groups
#' set.seed(123)
#' x1 <- rnorm(30, sd = 1)
#' x2 <- rnorm(30, sd = 1.5)
#' x3 <- rnorm(30, sd = 2)
#' x <- cbind(x1, x2, x3)
#'
#' # Perform all pairwise variance comparisons
#' result <- comvar.mcp(x)
#' print(result)
#' }
comvar.mcp<-function(x,method='hoch',SEED=TRUE){
#
#  Compare the variances of J indepenent variables.
#  Perform all pairwise comparisons using
#   a slight extension of  HC4 vesion of the Morgan-Pitman test
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
J=length(x)
CC=(J^2-J)/2
output<-matrix(0,CC,7)
dimnames(output)<-list(NULL,c('Var','Var','SD 1','SD 2','Dif','p.value','Adj.p.value'))
ic=0
for (j in 1:J){
for (k in 1:J){
if (j < k){
ic=ic+1
a=varcom.IND.MP(x[[j]],x[[k]],SEED=SEED)
a=pool.a.list(a)
output[ic,1]=j
output[ic,2]=k
output[ic,3:4]=a[1:2]
output[ic,3:4]=sqrt(output[ic,3:4])
output[ic,5]=output[ic,3]-output[ic,4]
output[ic,6]=a[3]
}}}
output[,7]=p.adjust(output[,6],method=method)
output
}

#' Compare Two Dependent Groups via Percentile Bootstrap
#'
#' Compares measures of location for two dependent (paired) groups using a
#' percentile bootstrap method. By default, uses trimmed means as the measure
#' of location.
#'
#' @inheritParams common-params
#' @param x Vector of observations for group 1, or a two-column matrix/data frame.
#' @param y Vector of observations for group 2 (optional if `x` is a matrix).
#' @param alpha Significance level (default: 0.05).
#' @param est Estimator function for location (default: `tmean`).
#' @param plotit Logical. If `TRUE`, creates a plot (default: `FALSE`).
#' @param dif Logical. If `TRUE`, analyzes difference scores; if `FALSE`, compares
#'   marginal trimmed means (default: `TRUE`).
#' @param nboot Number of bootstrap samples (default: `NA`, which uses 1000).
#' @param xlab Label for x-axis (default: "Group 1").
#' @param ylab Label for y-axis (default: "Group 2").
#' @param pr Logical. If `TRUE`, prints whether difference scores were used (default: `TRUE`).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A matrix with one row containing:
#'   \itemize{
#'     \item If `dif=FALSE`: `Est.1`, `Est.2`, `Est.dif`, `p.value`, `ci.lower`, `ci.upper`
#'     \item If `dif=TRUE`: `Est.typical.dif`, `p.value`, `ci.lower`, `ci.upper`
#'   }
#'
#' @details
#' This is a wrapper function that calls `rmmcppb` for convenience. It provides
#' a simpler interface for the common case of comparing two dependent groups.
#'
#' When `dif=TRUE`, the function analyzes the distribution of difference scores
#' (X - Y). When `dif=FALSE`, it compares the marginal distributions of X and Y.
#'
#' @seealso \code{\link{rmmcppb}}, \code{\link{yuend}}, \code{\link{tmean}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate paired data
#' set.seed(123)
#' pre <- rnorm(30, mean = 100)
#' post <- pre + rnorm(30, mean = 5, sd = 10)
#'
#' # Compare using difference scores
#' two.dep.pb(pre, post, nboot = 1000)
#'
#' # Compare marginal means
#' two.dep.pb(pre, post, dif = FALSE, nboot = 1000)
#' }
two.dep.pb<-function(x,y=NULL,alpha=.05,est=tmean,plotit=FALSE,dif=TRUE,
nboot=NA,xlab='Group 1',ylab='Group 2',pr=TRUE,SEED=TRUE,...){
#
#  Two dependent groups
# Compare measures of location  via a percentile bootstrap.
#  Trimmed mean	used by	default.
#
#  Calls rmmcppb, provided for convenience
#  nboot, number of bootstrap samples defaults to 1000
#
if(pr){
if(dif)print('dif=TRUE, difference scores were used')
if(!dif)print('dif=FALSE, marginal trimmed means  were used')
}
if(is.null(y)){
if(ncol(x)!=2)stop('y is null  so x should have two columns')
}
if(!is.null(y)){
xy=cbind(x,y)
xy=elimna(xy)
x=xy[,1]
y=xy[,2]
}
e=apply(cbind(x,y),2,est,...)
a=rmmcppb(x,y,est=est,nboot=nboot,alpha=alpha,SR=FALSE,SEED=SEED,
plotit=plotit,dif=dif,BA=FALSE,pr=FALSE,...)$output
if(!dif){
output=matrix(c(e[1],e[2],a[1,2],a[1,3],a[1,5],a[1,6]),nrow=1)
dimnames(output)=list(NULL,c('Est.1','Est.2','Est.dif','p.value','ci.lower','ci.upper'))
}
if(dif){
output=matrix(c(a[1,2],a[1,3],a[1,5],a[1,6]),nrow=1)
dimnames(output)=list(NULL,c('Est.typical.dif','p.value','ci.lower','ci.upper'))
}
output
}

#' Estimate P(X < Y) for WMW Test (Internal Helper)
#'
#' Internal helper function to estimate P(X < Y) for two independent random
#' variables. Used in bootstrap and simulation contexts.
#'
#' @param m Matrix with two columns (group 1 and group 2).
#' @param M Alternative matrix (default: same as `m`).
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#'
#' @return Scalar estimate of P(X < Y).
#'
#' @details
#' This is an internal helper function used primarily for bootstrap resampling
#' in WMW-based methods. It samples from the provided matrices and computes
#' the probability estimate using `bmp`.
#'
#' @keywords internal
#' @seealso \code{\link{bmp}}, \code{\link{wmw}}
wmw.est.only<-function(m,M=m,n1,n2){
#
# For two independent random variables,
# estimate P(X<Y)
#
# m = matrix, ncol=2
#
if(sum(is.na(m[,1])!=n1))x1=sample(M[,1])
else
x1=sample(m[,1])
if(sum(is.na(m[,2])!=n1))x2=sample(M[,2])
else
x2=sample(m[,2])
e=bmp(x1,x2)$phat
e
}

#' Reiczigel et al. Improved Wilcoxon-Mann-Whitney Test
#'
#' Performs the Reiczigel et al. (2005) improvement of the Wilcoxon-Mann-Whitney
#' test using a bootstrap approach to better control Type I error rates.
#'
#' @inheritParams common-params
#' @param x Vector of observations for group 1.
#' @param y Vector of observations for group 2.
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#'
#' @return A list with component:
#'   \itemize{
#'     \item `p.value`: P-value from the improved WMW test
#'   }
#'
#' @details
#' This function implements the improved WMW test proposed by Reiczigel et al. (2005).
#' The method uses a bootstrap resampling approach after centering the data based
#' on rank-based location estimates.
#'
#' The test procedure:
#' 1. Ranks all observations from both groups
#' 2. Computes a location shift estimate using trimmed means on ranks (tr=0)
#' 3. Centers the data by subtracting the shift estimate
#' 4. Performs bootstrap resampling to obtain the null distribution
#' 5. Computes a two-sided p-value
#'
#' This approach can provide better Type I error control than the classical WMW test
#' in certain situations.
#'
#' @references
#' Reiczigel, J., Zakarias, I., & Rozsa, L. (2005). A bootstrap test of stochastic
#' equality of two populations. The American Statistician, 59(2), 156-161.
#'
#' @seealso \code{\link{wmw}}, \code{\link{yuen}}
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(30)
#' y <- rnorm(35, mean = 0.3)
#'
#' # Perform improved WMW test
#' result <- wmw.RZR(x, y, nboot = 1000)
#' print(result)
#' }
wmw.RZR<-function(x,y,nboot=1000,SEED=TRUE){
#
#   Perform the Reiczigel et al. (2005) improvement of of the
# Wilcoxon--Mann--Whitney  test
#
if(SEED)set.seed(2)
val=0
n1=length(x)
n2=length(y)
xy=rank(c(x,y))
N=n1+n2
n1p1=n1+1
a=yuen(xy[1:n1],xy[n1p1:N],tr=0)$teststat
#print(a)
LOC=loc2dif(x,y)
x=x-a
y=y-a
#print(yuen(x,x2,tr=0))
bval=0
for(i in 1:nboot){
z1=sample(x,n1,replace=TRUE)
z2=sample(y,n2,replace=TRUE)
XY=rank(c(z1,z2))
bval[i]=yuen(XY[1:n1],XY[n1p1:N],tr=0)$teststat
}
#print(bval[1:10])
pv1=mean(a>bval)
pv2=mean(a<bval)

pv=2*min(pv1,pv2)
list(p.value=pv)
}

#' Compare Variances for Astigmatism Data
#'
#' Compares variances between two groups for astigmatism data (or any paired
#' bivariate data). Specifically designed to compare variances of the first
#' variable between groups and the second variable between groups.
#'
#' @param m1 Matrix or data frame with two columns (two variables for group 1).
#' @param m2 Matrix or data frame with two columns (two variables for group 2).
#' @param method Method for p-value adjustment (default: "holm").
#'   See \code{\link[stats]{p.adjust}} for options.
#'
#' @return A matrix with two rows (one per variable comparison) and columns:
#'   \itemize{
#'     \item `VAR 1`: Variance estimate for group 1
#'     \item `VAR 2`: Variance estimate for group 2
#'     \item `Ratio`: Ratio of variances (VAR 1 / VAR 2)
#'     \item `p.value`: Unadjusted p-value
#'     \item `p.adjusted`: Adjusted p-value
#'   }
#'
#' @details
#' This function is specifically designed for astigmatism data where each
#' subject/group has two measurements. It compares:
#' - Row 1: Variances of first variable (m1[,1] vs m2[,1])
#' - Row 2: Variances of second variable (m1[,2] vs m2[,2])
#'
#' The Morgan-Pitman test (`varcom.IND.MP`) is used for each comparison,
#' and p-values are adjusted for multiple testing.
#'
#' @seealso \code{\link{varcom.IND.MP}}, \code{\link{comvar.mcp}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Astigmatism data: two measurements per eye
#' set.seed(123)
#' m1 <- matrix(rnorm(60, sd = 1), ncol = 2)   # Group 1
#' m2 <- matrix(rnorm(60, sd = 1.5), ncol = 2) # Group 2
#'
#' # Compare variances for both measurements
#' result <- comvar2.astig(m1, m2)
#' print(result)
#' }
comvar2.astig<-function(m1,m2,method='holm'){
#
#  This function is designed specifically for dealing with astigmatism.
#
#  m= matrix or data frame, four columns
#
# compare variances of m[,1] vs m[,3] as well as m[,2] vs m[,4]
#
output=matrix(NA,2,5)
a=varcom.IND.MP(m1[,1], m2[,1])
output[1,1]=a$est1
output[1,2]=a$est2
output[1,3]=a$est1/a$est2
output[1,4]=a$p.value
a=varcom.IND.MP(m1[,2], m2[,2])
output[2,1]=a$est1
output[2,2]=a$est2
output[2,3]=a$est1/a$est2
output[2,4]=a$p.value
output[,5]=p.adjust(output[,4],method=method)
dimnames(output)=list(NULL,c('VAR 1','VAR 2','Ratio','p.value','p.adjusted') )
output
}

#' WMW-Based ANCOVA with Bootstrap Standard Errors
#'
#' Performs ANCOVA (analysis of covariance) based on a heteroscedastic analog
#' of the Wilcoxon-Mann-Whitney test. Estimates P(Y1 < Y2) at specified
#' covariate values with bootstrap standard errors.
#'
#' @inheritParams common-params
#' @param x1 Matrix or vector of covariate(s) for group 1.
#' @param y1 Vector of response variable for group 1.
#' @param x2 Matrix or vector of covariate(s) for group 2.
#' @param y2 Vector of response variable for group 2.
#' @param pts Matrix of covariate values at which to compare groups. If `NULL`,
#'   uses values from `ancova.KMS`.
#' @param nboot Number of bootstrap samples (default: 100).
#' @param SEED Logical. If `TRUE`, sets random seed for reproducibility (default: `TRUE`).
#' @param MC Logical. If `TRUE`, uses parallel processing via `mclapply` (default: `FALSE`).
#' @param null.value Null hypothesis value for P(Y1 < Y2) (default: 0.5).
#' @param xout Logical. If `TRUE`, removes outliers in covariates (default: `FALSE`).
#' @param outfun Function for outlier detection (default: `outpro`).
#' @param alpha Significance level (default: 0.05).
#' @param ... Additional arguments passed to outlier detection function.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `Est`: Vector of P(Y1 < Y2) estimates at each covariate value
#'     \item `SE`: Bootstrap standard errors
#'     \item `test.stat`: Test statistics (standardized differences)
#'     \item `conf.int`: Matrix of confidence intervals
#'     \item `p.value`: P-values for each comparison
#'   }
#'
#' @details
#' This function extends the Wilcoxon-Mann-Whitney approach to ANCOVA by
#' estimating P(Y1 < Y2) conditional on specific covariate values. Bootstrap
#' resampling is used to obtain standard errors and confidence intervals.
#'
#' The method allows for heteroscedasticity and does not assume linearity or
#' normality. It provides a robust alternative to classical ANCOVA.
#'
#' @seealso \code{\link{wmw.anc}}, \code{\link{ancova.KMS}}, \code{\link{outpro}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate data with covariate
#' set.seed(123)
#' n1 <- 40; n2 <- 40
#' x1 <- rnorm(n1)
#' y1 <- 2*x1 + rnorm(n1)
#' x2 <- rnorm(n2)
#' y2 <- 2*x2 + rnorm(n2, mean = 1)
#'
#' # Compare at specific covariate values
#' pts <- matrix(c(-1, 0, 1), ncol = 1)
#' result <- wmw.ancbse(x1, y1, x2, y2, pts = pts, nboot = 200)
#' print(result)
#' }
wmw.ancbse<-function(x1,y1,x2,y2,pts,nboot=100,SEED=TRUE,MC=FALSE,null.value=.5,xout=FALSE,
outfun=outpro,alpha=.05,...){
#
#  ANCOVA based on a heteroscedastic analog of the
#  Wilcoxon--Mann-Whitney test
#
#  pts indicates the covariance values  for which the groups will be compared.

x1<-as.matrix(x1)
p1<-ncol(x1)+1
p<-ncol(x1)
xy<-cbind(x1,y1)
xy<-elimna(xy)
x1<-xy[,1:p]
y1<-xy[,p1]

x2<-as.matrix(x2)
if(p!=ncol(x2))stop('Number of col. for x1 is not equal to the number of col. for x2')
xy<-cbind(x2,y2)
xy<-elimna(xy)
x2<-xy[,1:p]
y2<-xy[,p1]

if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE,...)$keep
m<-m[flag,]
x1<-m[,1:p]
y1<-m[,p1]
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE,...)$keep
m<-m[flag,]
x2<-m[,1:p]
y2<-m[,p1]
}

if(is.null(pts))pts=ancova.KMS(x1,y1,x2,y2,plotit=FALSE)[,1]
if(SEED)set.seed(2)
e=wmw.anc(x1,y1,x2,y2,pts=pts)
npt=length(pts)
ci=matrix(NA,npt,2)
LAB=NULL
for(j in 1:npt)LAB[j]=paste('pts',j)
dimnames(ci)=list(LAB,c('ci.low','ci.up'))
E=matrix(NA,nboot,npt)
n1=length(y1)
n2=length(y2)
bl=list()
for(k in 1:nboot){
id1=sample(n1,replace=TRUE)
id2=sample(n2,replace=TRUE)
bl[[k]]=list(x1[id1],y1[id1],x2[id2],y2[id2])
}
if(!MC)tv=lapply(bl,wmw.ancbse.sub,pts)
else{
tv=mclapply(bl,wmw.ancbse.sub,pts)
}
E=t(matl(tv))
se=apply(E,2,sd)
ci[,1]<-e-qnorm(1-alpha/2)*se
ci[,2]<-e+qnorm(1-alpha/2)*se
test<-(e-null.value)/se
sig<-2*(1-pnorm(abs(test)))
list(Est=e,SE=se,test.stat=test,conf.int=ci,p.value=sig)
}

#' Bootstrap Helper for WMW ANCOVA (Internal)
#'
#' Internal helper function for `wmw.ancbse` that computes WMW ANCOVA estimates
#' for bootstrap samples.
#'
#' @param m List containing bootstrapped data: x1, y1, x2, y2.
#' @param pts Matrix of covariate values.
#'
#' @return Vector of P(Y1 < Y2) estimates at each covariate value.
#'
#' @keywords internal
#' @seealso \code{\link{wmw.ancbse}}, \code{\link{wmw.anc}}
wmw.ancbse.sub<-function(m,pts){
v=wmw.anc(m[[1]],m[[2]],m[[3]],m[[4]],pts=pts)
v
}

#' Estimate P(Y1 < Y2) Conditional on Covariate Values
#'
#' Estimates P(Y1 < Y2) given that the covariate X equals specific values.
#' Used as the core computation for WMW-based ANCOVA methods.
#'
#' @inheritParams common-params
#' @param x1 Matrix or vector of covariate(s) for group 1.
#' @param y1 Vector of response variable for group 1.
#' @param x2 Matrix or vector of covariate(s) for group 2.
#' @param y2 Vector of response variable for group 2.
#' @param pts Matrix of covariate values at which to estimate P(Y1 < Y2).
#' @param xout Logical. If `TRUE`, removes outliers in covariates (default: `FALSE`).
#' @param outfun Function for outlier detection (default: `outpro`).
#'
#' @return Vector of P(Y1 < Y2) estimates, one for each row of `pts`.
#'
#' @details
#' For each specified covariate value in `pts`, this function estimates the
#' probability that a randomly selected observation from group 1 is less than
#' a randomly selected observation from group 2, conditional on the covariate
#' value.
#'
#' The estimation uses regression-based conditional distributions via
#' `reg.con.dist` to approximate the conditional distributions at each
#' covariate value, then computes the probability empirically.
#'
#' @keywords internal
#' @seealso \code{\link{wmw.ancbse}}, \code{\link{reg.con.dist}}, \code{\link{outpro}}
wmw.anc<-function(x1,y1,x2,y2,pts,xout=FALSE,outfun=outpro){
#
#  Estimate P(y1<y2)  given that the covariate X is equal to the values in pts.
#
x1<-as.matrix(x1)
p1<-ncol(x1)+1
p<-ncol(x1)
xy<-cbind(x1,y1)
xy<-elimna(xy)
x1<-xy[,1:p]
y1<-xy[,p1]

pts=as.matrix(pts)
if(p>1){
if(p!=ncol(pts))stop('pts should be a matrix with',paste(p),'columns')
}

x2<-as.matrix(x2)
if(p!=ncol(x2))stop('Number of col. for x1 is not equal to the number of col. for x2')
xy<-cbind(x2,y2)
xy<-elimna(xy)
x2<-xy[,1:p]
y2<-xy[,p1]

if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE,...)$keep
m<-m[flag,]
x1<-m[,1:p]
y1<-m[,p1]
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE,...)$keep
m<-m[flag,]
x2<-m[,1:p]
y2<-m[,p1]
}

e=NA
for(i in 1:nrow(pts)){
d1=reg.con.dist(x1,y1,pts=pts[i,])
d2=reg.con.dist(x2,y2,pts=pts[i,])
p=NA
for(j in 1:99)p[j]=mean(d1[j]<=d2)
e[i]=mean(p)
}
e
}

#' Compare Single Quantile Between Two Groups (Internal Helper)
#'
#' Internal helper function for `qcomhd` that compares a single quantile
#' between two independent groups using the Harrell-Davis estimator and
#' percentile bootstrap.
#'
#' @param x Vector of observations for group 1.
#' @param y Vector of observations for group 2.
#' @param q Quantile to compare (single value between 0 and 1).
#' @param alpha Significance level (default: 0.05).
#' @param nboot Number of bootstrap samples (default: 2000).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param MC Logical. If `TRUE`, uses parallel processing (default: `TRUE`).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `est.1`: Harrell-Davis estimate for group 1
#'     \item `est.2`: Harrell-Davis estimate for group 2
#'     \item `ci`: Confidence interval for the difference
#'     \item `p.value`: P-value for the test
#'     \item `sq.se`: Squared standard error (variance)
#'     \item `n1`, `n2`: Sample sizes
#'   }
#'
#' @keywords internal
#' @seealso \code{\link{qcomhd}}, \code{\link{hd}}
qcomhd.sub<-function(x,y,q,alpha=.05,nboot=2000,SEED=TRUE,MC=TRUE){
#
x<-x[!is.na(x)] # Remove any missing values in x
y<-y[!is.na(y)] # Remove any missing values in y
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
datax=listm(t(datax))
datay=listm(t(datay))
if(MC){
bvecx<-mclapply(datax,hd,q,mc.preschedule=TRUE)
bvecy<-mclapply(datay,hd,q,mc.preschedule=TRUE)
}
else{
bvecx<-lapply(datax,hd,q)
bvecy<-lapply(datay,hd,q)
}
bvecx=as.vector(matl(bvecx))
bvecy=as.vector(matl(bvecy))
bvec<-sort(bvecx-bvecy)
low<-round((alpha/2)*nboot)+1
up<-nboot-low
temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
sig.level<-2*(min(temp,1-temp))
se<-var(bvec)
list(est.1=hd(x,q),est.2=hd(y,q),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}

#' Compare Quantiles Between Two Independent Groups
#'
#' @description
#' Compares quantiles between two independent groups using the Harrell-Davis estimator
#' and percentile bootstrap. This provides a robust method for comparing multiple quantiles
#' simultaneously with control for multiple testing via Hochberg's method.
#'
#' @inheritParams common-params
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param est Estimator function to use (default \code{hd} for Harrell-Davis). Can also
#'   use other quantile estimators.
#' @param q Numeric vector of quantiles to compare (default c(0.1, 0.25, 0.5, 0.75, 0.9)).
#' @param nboot Number of bootstrap samples (default 4000).
#' @param plotit Logical; if TRUE, creates a plot showing estimates and confidence intervals.
#' @param xlab Label for x-axis in plot.
#' @param ylab Label for y-axis in plot.
#' @param ADJ.CI Logical; if TRUE (default), confidence intervals are adjusted to match the
#'   alpha level used by the corresponding test statistic (e.g., if testing at 0.05/3 level,
#'   CI has 1-0.05/3 coverage).
#' @param MC Logical; if TRUE, uses parallel processing via mclapply.
#'
#' @return A matrix with columns:
#'   \item{q}{The quantile being compared}
#'   \item{n1}{Sample size of group 1}
#'   \item{n2}{Sample size of group 2}
#'   \item{est.1}{Quantile estimate for group 1}
#'   \item{est.2}{Quantile estimate for group 2}
#'   \item{est.1_minus_est.2}{Difference in quantile estimates}
#'   \item{ci.low}{Lower confidence limit}
#'   \item{ci.up}{Upper confidence limit}
#'   \item{p-value}{P-value for testing equality}
#'   \item{adj.p.value}{Adjusted p-value using Hochberg's method}
#'
#' @details
#' The function uses the Harrell-Davis estimator which handles tied values well. For
#' comparing lower or upper quartiles, both power and Type I error control compare
#' favorably to other methods.
#'
#' Family-wise error rate is controlled using Hochberg's method across all quantile
#' comparisons. The plot shows differences (marked with *) connected by lines, with
#' confidence interval limits marked with +.
#'
#' @seealso \code{\link{qcomhdMC}}, \code{\link{pb2gen}}, \code{\link{pb2genMC}}, \code{\link{hd}}
#'
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 0.5)
#'
#' # Compare default quantiles
#' result <- qcomhd(x, y)
#'
#' # Compare specific quantiles
#' qcomhd(x, y, q = c(0.25, 0.5, 0.75), plotit = TRUE)
qcomhd<-function(x,y,est=hd,q=c(.1,.25,.5,.75,.9),nboot=4000,plotit=TRUE,SEED=TRUE,xlab='Group 1',ylab='Est.1-Est.2',alpha=.05,ADJ.CI=TRUE,MC=FALSE){
#
# Compare quantiles using pb2gen using trimmed version of the Harrell-Davis estimator
# Tied values are allowed.
#
#  ADJ.CI=TRUE means that the confidence intervals are adjusted based on the level used by the corresponding
#  test statistic. If a test is performed with at the .05/3 level, for example, the confidence returned has
#  1-.05/3 probability coverage.
#
# When comparing lower or upper quartiles, both power and the probability of Type I error
# compare well to other methods that have been derived.
# q: can be used to specify the quantiles to be compared
# q defaults to comparing the .1,.25,.5,.75, and .9 quantiles
#
#   Function returns p-values and critical p-values based on Hochberg's method.
#

if(SEED)set.seed(2)
pv=NULL
output=matrix(NA,nrow=length(q),ncol=10)
dimnames(output)<-list(NULL,c('q','n1','n2','est.1','est.2','est.1_minus_est.2','ci.low','ci.up','p-value','adj.p.value'))
for(i in 1:length(q)){
output[i,1]=q[i]
output[i,2]=length(elimna(x))
output[i,3]=length(elimna(y))
output[i,4]=hd(x,q=q[i])
output[i,5]=hd(y,q=q[i])
output[i,6]=output[i,4]-output[i,5]
temp=qcomhd.sub(x,y,nboot=nboot,q=q[i],SEED=FALSE,alpha=alpha,MC=MC)
output[i,7]=temp$ci[1]
output[i,8]=temp$ci[2]
output[i,9]=temp$p.value
}
temp=order(output[,9],decreasing=TRUE)
zvec=alpha/c(1:length(q))
zvec[temp]=zvec
if(ADJ.CI){
for(i in 1:length(q)){
if(!MC)temp=pb2gen(x,y,nboot=nboot,est=est,q=q[i],SEED=FALSE,alpha=zvec[i],pr=FALSE)
else
temp=pb2genMC(x,y,nboot=nboot,est=est,q=q[i],SEED=FALSE,alpha=zvec[i],pr=FALSE)
output[i,7]=temp$ci[1]
output[i,8]=temp$ci[2]
output[i,9]=temp$p.value
}
temp=order(output[,10],decreasing=TRUE)
}
output[,10]=p.adjust(output[,9],method='hoch')
if(plotit){
xax=rep(output[,4],3)
yax=c(output[,6],output[,7],output[,8])
plot(xax,yax,xlab=xlab,ylab=ylab,type='n')
points(output[,4],output[,6],pch='*')
lines(output[,4],output[,6])
points(output[,4],output[,7],pch='+')
points(output[,4],output[,8],pch='+')
}
output
}

#' Plot WMW ANCOVA: P(Y1 < Y2 | X) vs X
#'
#' @description
#' Plots WMW-based conditional probability estimates P(Y1 < Y2 | X=x)
#' as a function of the covariate X.
#'
#' @inheritParams common-params
#' @param x1 Covariate for group 1.
#' @param y1 Response for group 1.
#' @param x2 Covariate for group 2.
#' @param y2 Response for group 2.
#' @param q1 Quantile range for group 1 covariate (default: c(0.1, 0.9)).
#' @param q2 Quantile range for group 2 covariate (default: c(0.1, 0.9)).
#' @param npts Number of design points (default: 20).
#' @param pts Optional vector of design points (default: NULL for automatic).
#' @param xout Logical. If `TRUE`, removes outliers (default: `FALSE`).
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param xlab X-axis label (default: 'X').
#' @param ylab Y-axis label (default: 'P(Y1<Y2)').
#' @param ... Additional arguments passed to outlier function.
#'
#' @return A list with the covariate range used.
#'
#' @details
#' Single covariate only. Estimates the conditional WMW effect size at
#' multiple covariate values and creates a plot. Design points are chosen
#' within the overlapping quantile ranges of both groups.
#'
#' @seealso \code{\link{wmw.anc}}, \code{\link{wmw.ancp2}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Plot conditional WMW effect
#' set.seed(123)
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- x2 + 0.5 + rnorm(50)
#' wmw.anc.plot(x1, y1, x2, y2, npts = 15)
#' }
wmw.anc.plot<-function(x1,y1,x2,y2,q1=c(.1,.9),q2=c(.1,.9),npts=20,
pts=NULL,xout=FALSE,outfun=outpro,xlab='X',ylab='P(Y1<Y2)',...){
#
#    Plot estimates of P(x)= P(Y1 <Y2|x)
#    Current version, single covariate assumed
#
x1<-as.matrix(x1)
p1<-ncol(x1)+1
if(p1!=2)stop('Single covariate only is allowed')
p<-ncol(x1)
xy<-cbind(x1,y1)
xy<-elimna(xy)
x1<-xy[,1:p]
y1<-xy[,p1]
if(xout){
m<-cbind(x1,y1)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x1,plotit=FALSE,...)$keep
m<-m[flag,]
x1<-m[,1:p]
y1<-m[,p1]
m<-cbind(x2,y2)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x2,plotit=FALSE,...)$keep
m<-m[flag,]
x2<-m[,1:p]
y2<-m[,p1]
}
if(is.null(pts)){
v1=qest(x1,q1)
v2=qest(x2,q2)
bot=max(v1[1],v2[1])
top=min(v1[2],v2[2])
pts=seq(bot,top,length.out=npts)
}
v=NA
for(i in 1:length(pts))v[i]=wmw.anc(x1,y1,x2,y2,pts[i])
plot(pts,v,xlab=xlab,ylab=ylab,ylim=c(0,1))
list(Range=c(pts[1],pts[length(pts)]))
}

#' WMW-Based ANCOVA with Bootstrap for Multiple Covariates
#'
#' @description
#' Performs ANCOVA using a heteroscedastic Wilcoxon-Mann-Whitney analog
#' with bootstrap standard errors. Supports multiple covariates.
#'
#' @inheritParams common-params
#' @param x1 Covariate matrix for group 1.
#' @param y1 Response vector for group 1.
#' @param x2 Covariate matrix for group 2.
#' @param y2 Response vector for group 2.
#' @param pts Design points matrix (default: NULL for automatic selection).
#' @param nboot Number of bootstrap samples (default: 100).
#' @param alpha Significance level (default: 0.05).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param MC Logical. If `TRUE`, uses parallel processing (default: `FALSE`).
#' @param npts Number of design points if pts=NULL (default: 30).
#' @param profun Projection depth function (default: \code{prodepth}).
#' @param BOTH Logical. If `TRUE`, selects points from pooled data (default: `TRUE`).
#' @param plotit Logical. If `TRUE`, creates plot (default: `TRUE`).
#' @param xlab X-axis label (default: 'X1').
#' @param ylab Y-axis label (default: 'X2').
#' @param method Multiple testing adjustment method (default: 'hoch').
#' @param xout Logical. If `TRUE`, removes outliers (default: `FALSE`).
#' @param outfun Outlier detection function (default: \code{outpro}).
#'
#' @return Matrix with WMW estimates, standard errors, and p-values at design points.
#'
#' @details
#' Assumes a linear model is reasonable. Uses projection depth to select
#' design points and bootstrap to estimate standard errors. Tests are
#' adjusted for multiple comparisons.
#'
#' @seealso \code{\link{wmw.ancp2}}, \code{\link{wmw.ancbsep2.sub}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Two-covariate ANCOVA
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 2)
#' y1 <- rowSums(x1) + rnorm(50)
#' x2 <- matrix(rnorm(100), ncol = 2)
#' y2 <- rowSums(x2) + 0.5 + rnorm(50)
#' wmw.ancbsep2(x1, y1, x2, y2, nboot = 200)
#' }
wmw.ancbsep2<-function(x1,y1,x2,y2,pts=NULL,nboot=100,alpha=.05,SEED=TRUE,MC=FALSE,
npts=30,profun=prodepth,BOTH=TRUE,plotit=TRUE,xlab='X1',ylab='X2',method='hoch',
xout=FALSE,outfun=outpro){
#
#  ANCOVA based on a heteroscedastic analog of the
#  Wilcoxon--Mann-Whitney test
#. Assume a linear model is reasonable
#
if(SEED)set.seed(2)
x1=as.matrix(x1)
p=ncol(x1)
p1=p+1
m1=elimna(cbind(x1,y1))
x1=m1[,1:p]
y1=m1[,p1]
x2=as.matrix(x2)
p=ncol(x2)
p1=p+1
m2=elimna(cbind(x2,y2))
x2=m2[,1:p]
y2=m2[,p1]
if(xout){
flag<-outfun(x1,plotit=FALSE)$keep
m<-m[flag,]
x1<-m[,1:p]
y1<-m[,p1]
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE)$keep
m<-m[flag,]
x2<-m[,1:p]
y2<-m[,p1]
}
n1=length(y1)
n2=length(y2)
X=rbind(x1,x2)
if(is.null(pts)){
if(!BOTH){
d=profun(x1)
ior=order(d)
npts=min(min(n1),npts)
id=seq(1,min(n1),length.out=npts)
id=floor(id)
pts=x1[ior[id],]
}
if(BOTH){
n=n1+n2
X=rbind(x1,x2)
d=profun(X)
ior=order(d)
npts=min(min(n),npts)
id=seq(1,n,length.out=npts)
id=floor(id)
pts=X[ior[id],]
}
}
pts<-as.matrix(pts)
e=wmw.ancp2(x1,y1,x2,y2,pts=pts)
npt=nrow(pts)
ci=matrix(NA,npt,2)
E=matrix(NA,nboot,npt)
n1=length(y1)
n2=length(y2)
bl=list()
for(k in 1:nboot){
id1=sample(n1,replace=TRUE)
id2=sample(n2,replace=TRUE)
bl[[k]]=list(x1[id1,],y1[id1],x2[id2,],y2[id2])
}
if(!MC)tv=lapply(bl,wmw.ancbsep2.sub,pts)
else{
tv=mclapply(bl,wmw.ancbsep2.sub,pts)
}
E=t(matl(tv))
se=apply(E,2,sd)
test=(e-.5)/se
se=apply(E,2,sd)
ci[,1]<-e-qnorm(1-alpha/2)*se
ci[,2]<-e+qnorm(1-alpha/2)*se
sig<-2*(1-pnorm(abs(test)))
padj=p.adjust(sig,method=method)
output=cbind(pts[,1],pts[,2],e,se,test,ci[,1],ci[,2],sig,padj)
dimnames(output)=list(NULL,c('Pt1','Pt2','Est','SE','test','ci.low','ci.up','p-value','p.adj'))
if(plotit){
if(!BOTH)plot(x1[,1],x1[,2],xlab=xlab,ylab=ylab,pch='.')
if(BOTH)plot(X[,1],X[,2],xlab=xlab,ylab=ylab,pch='.')
points(pts[,1],pts[,2],pch='*')
flag=output[,8]<=alpha
if(sum(flag)>0)points(pts[flag,1],pts[flag,2],pch='o')
}
output
}

#' Internal Helper for wmw.ancbsep2 Bootstrap
#'
#' @description
#' Internal function used by \code{wmw.ancbsep2} to compute WMW estimates
#' from bootstrap samples.
#'
#' @param m List with bootstrap sample data.
#' @param pts Design points matrix.
#'
#' @return Vector of WMW estimates at design points.
#'
#' @keywords internal
wmw.ancbsep2.sub<-function(m,pts){
v=wmw.ancp2(m[[1]],m[[2]],m[[3]],m[[4]],pts=pts)
v
}

#' 3D Surface Plot of WMW ANCOVA Effect Sizes
#'
#' @description
#' Creates a 3D surface plot of WMW-based effect sizes across two covariates.
#'
#' @inheritParams common-params
#' @param x1 Covariate matrix for group 1 (2 columns).
#' @param y1 Response vector for group 1.
#' @param x2 Covariate matrix for group 2 (2 columns).
#' @param y2 Response vector for group 2.
#' @param pts Design points matrix (default: NULL for automatic).
#' @param xlab X-axis label (default: 'X1').
#' @param ylab Y-axis label (default: 'X2').
#' @param zlab Z-axis label (default: 'Effect Size').
#' @param REV Logical. If `TRUE`, reverses perspective (default: `FALSE`).
#' @param xout Logical. If `TRUE`, removes outliers (default: `FALSE`).
#' @param outfun Outlier detection function (default: \code{outpro}).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param theta Azimuthal viewing angle (default: 50).
#' @param phi Colatitude viewing angle (default: 25).
#'
#' @return Invisibly returns the effect size matrix.
#'
#' @details
#' Requires exactly 2 covariates. Creates a 3D perspective plot showing
#' how the WMW effect size varies across the covariate space.
#'
#' @seealso \code{\link{wmw.ancp2}}, \code{\link{wmw.ancbsep2}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 3D plot of WMW effect
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 2)
#' y1 <- rowSums(x1) + rnorm(50)
#' x2 <- matrix(rnorm(100), ncol = 2)
#' y2 <- rowSums(x2) + rnorm(50)
#' ancovap2.wmw.plot(x1, y1, x2, y2)
#' }
ancovap2.wmw.plot<-function(x1,y1,x2,y2,pts=NULL,xlab='X1',ylab='X2',zlab='Effect Size',REV=FALSE,
xout=FALSE,outfun=outpro,SEED=TRUE, theta = 50, phi = 25){
#
#
# Two covariates, plot the Wilcoxon--Mann--Whitney  measure of effect size
# using a smoother
#
#
xy=elimna(cbind(x1,y1))
if(ncol(xy)!=3)stop('Only two covariates can  be used')
x1=xy[,1:2]
y1=xy[,3]
xy=elimna(cbind(x2,y2))
x2=xy[,1:2]
y2=xy[,3]
if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE)$keep
m<-m[flag,]
x1<-m[,1:2]
y1<-m[,3]
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE)$keep
m<-m[flag,]
x2<-m[,1:2]
y2<-m[,3]
}
if(is.null(pts))pts=rbind(x1,x2)
N=nrow(pts)
e=wmw.ancp2(x1,y1,x2,y2,pts=pts)
if(N<25){
library(scatterplot3d)
scatterplot3d(pts[,1],pts[,2],e,xlab=xlab,ylab=ylab,zlab=zlab)
}
if(N>=25){
if(!REV)f=lplot(pts,e,xlab=ylab,ylab=xlab,zlab=zlab,ticktype='det',pr=FALSE,theta=theta,phi=phi)
else f=lplot(pts[,c(2,1)],e,xlab=xlab,ylab=ylab,zlab=zlab,ticktype='det',pr=FALSE,theta=theta,phi=phi)
}
list(Number_of_points_used_is=N)
}

#' Conditional WMW Effect Size at Design Points
#'
#' @description
#' Estimates the conditional WMW effect size at specified design points
#' for two independent groups with covariates.
#'
#' @inheritParams common-params
#' @param x1 Covariate matrix/vector for group 1.
#' @param y1 Response vector for group 1.
#' @param x2 Covariate matrix/vector for group 2.
#' @param y2 Response vector for group 2.
#' @param pts Design points matrix (default: NULL for 3 automatic points).
#' @param xout Logical. If `TRUE`, removes outliers (default: `FALSE`).
#' @param outfun Outlier detection function (default: \code{outpro}).
#'
#' @return Vector of WMW effect size estimates at each design point.
#'
#' @details
#' For regression lines from two independent groups, estimates the conditional
#' WMW effect size P(Y1 < Y2 | X = x) at each point in pts. If pts is NULL,
#' automatically selects three design points based on the data.
#'
#' @seealso \code{\link{wmw.ancbsep2}}, \code{\link{wmw.anc}}, \code{\link{reg.con.dist}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Estimate conditional WMW effect
#' set.seed(123)
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- x2 + 0.5 + rnorm(50)
#' pts <- matrix(c(-1, 0, 1), ncol = 1)
#' wmw.ancp2(x1, y1, x2, y2, pts = pts)
#' }
wmw.ancp2<-function(x1,y1,x2,y2,pts=NULL,xout=FALSE,outfun=outpro){
#
#  For the regression lines corresponding to  two independent groups
#  estimate the conditional  WMW effect size for each point in pts
#
# pts=NULL: three  points used that are determined based on the data
#
x1=as.matrix(x1)
p=ncol(x1)
p1=p+1
m1=elimna(cbind(x1,y1))
x1=m1[,1:p]
y1=m1[,p1]
x2=as.matrix(x2)
p=ncol(x2)
p1=p+1
m2=elimna(cbind(x2,y2))
x2=m2[,1:p]
y2=m2[,p1]
if(xout){
m<-cbind(x1,y1)
if(identical(outfun,reglev))flag=outfun(x1,y1,plotit=FALSE)$keep
else
flag<-outfun(x1,plotit=FALSE)$keep
m<-m[flag,]
x1<-m[,1:p]
y1<-m[,p1]
m<-cbind(x2,y2)
if(identical(outfun,reglev))flag=outfun(x2,y2,plotit=FALSE)$keep
else
flag<-outfun(x2,plotit=FALSE)$keep
m<-m[flag,]
x2<-m[,1:p]
y2<-m[,p1]
}
if(is.null(pts[1])){
x1<-as.matrix(x1)
pts<-ancdes(x1)
pts=unique(pts)
}
pts<-as.matrix(pts)
if(is.null(pts)){
x1<-as.matrix(x1)
pts<-ancdes(x1)
pts=unique(pts)
}
e=NA
PV=NA
n1=length(y1)
n2=length(y2)
for(i in 1:nrow(pts)){
d1=reg.con.dist(x1,y1,pts=pts[i,])
d2=reg.con.dist(x2,y2,pts=pts[i,])
p=NA
for(j in 1:99)p[j]=mean(d1[j]<=d2)
e[i]=mean(p)
}
e
}

#' Test Symmetry of Difference Distribution Using Reference Points
#'
#' @description
#' Tests whether the distribution of X-Y is symmetric about zero using
#' reference points. Can test either at quantiles or specified values.
#'
#' @inheritParams common-params
#' @param x Observations from group 1.
#' @param y Observations from group 2.
#' @param q Quantile for reference point if pts=NULL (default: 0.25).
#' @param pts Reference point(s) (default: NULL for quantile-based).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param alpha Significance level (default: 0.05).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param plotit Logical. If `TRUE`, plots distribution (default: `FALSE`).
#' @param xlab X-axis label (default: 'Difference').
#' @param ylab Y-axis label (default: 'Density').
#' @param estfun Quantile estimator function (default: \code{hdmq}).
#' @param plotfun Plotting function (default: \code{kerSORT}).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `L`: Estimated P(X-Y < -pts)
#'     \item `U`: Estimated P(X-Y > pts)
#'     \item `Est.dif`: Difference U - L
#'     \item `ci`: Bootstrap confidence interval
#'     \item `p.value`: P-value for symmetry test
#'   }
#'
#' @details
#' If X and Y have identical distributions, X-Y is symmetric about zero
#' and the sum of the q and 1-q quantiles equals zero. Tests this
#' using bootstrap. Should not be used with tied values.
#'
#' @seealso \code{\link{wmw.ref.mul}}, \code{\link{wmw.det}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test symmetry of difference distribution
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 0.3)
#' wmw.ref.dif(x, y, q = 0.25, nboot = 500)
#' }
wmw.ref.dif<-function(x,y,q=.25,pts=NULL,nboot=1000,alpha=.05,SEED=TRUE,
plotit=FALSE,xlab='Difference',ylab='Density',estfun=hdmq,
plotfun=kerSORT){
#
# If pts is specified, the goal is to make inferences about
# P(x-y< -pts)-P(x-y > pts)
#  using a percentile bootstrap method
#
# If pts is not specified, and make inferences
#  about the 1-q and q quantiles. If X and Y
# have identical distributions,  D=X-Y is symmetric about zero and the sum of the
# 1-q and qth quantiles is zero. Should not  be used when there are tied values
#
# If QC=FALSE and pts=NULL,
#  take pts to be estimate of the q quantile of D.
#
# if pts is not NULL, QC=FALSE is used
#
# Output:
# L=P(x-y< -pts)
# U = P(x-y > pts)
# Est.dif=U-L
#
QC=TRUE
if(!is.null(pts))QC=FALSE
if(is.null(pts)){
if(sum(q<.5)!=length(q))stop('All q values should be <=.5')
}
if(SEED)set.seed(2)
d=NA
x<-x[!is.na(x)]
y<-y[!is.na(y)]
n1=length(x)
n2=length(y)
if(!QC){
if(!is.null(pts)){
e=wmw.det(x,y,refp=pts,plotit=plotit,xlab=xlab,ylab=ylab,plotfun=plotfun)
est=e$dif
L=e$L
U=e$U
}
else{
pts=qest(outer(x,y,FUN='-'),q=q)
e=wmw.det(x,y,refp=pts,plotit=plotit,xlab=xlab,ylab=ylab,plotfun=plotfun)
est=e$dif
L=e$L
U=e$U
}}

if(QC){
d=outer(x,y,FUN='-')
d=as.vector(d)
qv=estfun(d,q=c(q,1-q))
if(plotit)plotfun(d,xlab=xlab,ylab=ylab)
est=qv[1]+qv[2]
L=qv[1]
U=qv[2]
}
for(i in 1:nboot){
id1=sample(n1,replace=TRUE)
id2=sample(n2,replace=TRUE)
if(!QC)d[i]=wmw.det(x[id1],y[id2],refp=pts,plotit=FALSE)$dif
else{
qv=estfun(outer(x[id1],y[id2],FUN='-'),q=c(q,1-q))
d[i]=qv[1]+qv[2]
}}
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
dif=sort(d)
ci=dif[icl]
ci[2]=dif[icu]
pv=mean(dif<0)+.5*mean(dif==0)
pv<-2*min(pv,1-pv)
list(L=L,U=U,Est.dif=est,ci=ci,p.value=pv)
}

#' Test Symmetry at Multiple Reference Points
#'
#' @description
#' Tests symmetry of the X-Y distribution at multiple reference points using
#' bootstrap. Performs multiple tests with p-value adjustment.
#'
#' @inheritParams common-params
#' @param x Observations from group 1.
#' @param y Observations from group 2.
#' @param refp Vector of reference points (default: NULL for quantile-based).
#' @param pts Alias for refp (default: NULL).
#' @param q Quantiles for reference points if refp=NULL (default: seq(0.6, 0.9, 0.1)).
#' @param center Logical. If `TRUE`, centers differences at median (default: `FALSE`).
#' @param estfun Quantile estimator (default: \code{hdmq} for Harrell-Davis).
#' @param alpha Significance level (default: 0.05).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param method Multiple testing adjustment method (default: 'BH').
#' @param plotit Logical. If `TRUE`, plots distribution (default: `FALSE`).
#' @param xlab X-axis label (default: 'Difference').
#' @param ylab Y-axis label (default: 'Density').
#' @param plotfun Plotting function (default: \code{kerSORT}).
#'
#' @return Matrix with columns for reference points, tail probabilities,
#'   differences, confidence intervals, p-values, and adjusted p-values.
#'
#' @details
#' Extension of \code{wmw.ref.dif} to multiple reference points with
#' family-wise error control via p-value adjustment.
#'
#' @seealso \code{\link{wmw.ref.dif}}, \code{\link{wmw.QC.mul}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Multiple reference points test
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 0.3)
#' wmw.ref.mul(x, y, nboot = 500)
#' }
wmw.ref.mul<-function(x,y,refp=NULL,pts=NULL,q=seq(.6,.9,.1),  center=FALSE, estfun=hdmq, alpha=.05,nboot=1000,SEED=TRUE,method='BH',plotit=FALSE,
xlab='Difference',ylab='Density',
plotfun=kerSORT){
#
#
# For multiple reference points, refp,
# make inferences about P(x-y< -refp) vs P(x-y > refp)
# refp can be constants chosen by the user. If not specified,
# refp are taken to be the .6(.1).9 estimated quantiles of the distribution of X-Y
#
# pts can be used to indicate specified reference points, refp
#
# To use the Harrell-Davis estimator, set estfun=hdmq
#
if(SEED)set.seed(2)
if(!is.null(pts))refp=pts
x=elimna(x)
y=elimna(y)
if(is.null(refp)){
m=outer(x,y,FUN='-')
m=as.vector(m)
morig=m

if(center)m=m-median(m)
refp=estfun(m,q)
}
np=length(refp)
output<-matrix(NA,np,8)
dimnames(output)=list(NULL,c('Pts','P(x-y<-Pts)' ,'P(x-y>Pts)','Dif','ci.low','ci.up','p.value','p.adj'))
for(i in 1:np){
e=wmw.ref.dif(x,y,pts=refp[i],alpha=alpha,nboot=nboot,SEED=FALSE)
output[i,1:7]=c(refp[i],e$L,e$U,e$Est.dif,e$ci[1],e$ci[2],e$p.value)
}
output[,8]=p.adjust(output[,7],method=method)
if(plotit)plotfun(as.vector(morig),xlab=xlab,ylab=ylab)
output
}

#' Test Symmetry at Multiple Quantile-Based Reference Points
#'
#' @description
#' Tests symmetry of X-Y distribution at multiple quantile-based reference
#' points. Similar to \code{wmw.ref.mul} but specifically for quantile-based testing.
#'
#' @inheritParams common-params
#' @param x Observations from group 1.
#' @param y Observations from group 2.
#' @param q Vector of quantiles less than 0.5 (default: seq(0.1, 0.4, 0.1)).
#' @param estfun Quantile estimator (default: \code{hdmq} for Harrell-Davis).
#' @param alpha Significance level (default: 0.05).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param method Multiple testing adjustment method (default: 'BH').
#' @param plotit Logical. If `TRUE`, plots distribution (default: `FALSE`).
#' @param xlab X-axis label (default: 'Difference').
#' @param ylab Y-axis label (default: 'Density').
#' @param plotfun Plotting function (default: \code{kerSORT}).
#'
#' @return Matrix with test results at each quantile including adjusted p-values.
#'
#' @details
#' For each quantile q < 0.5, tests whether the q and 1-q quantiles of X-Y
#' sum to zero (as they should under symmetry). Uses bootstrap with
#' multiple comparison adjustment.
#'
#' @seealso \code{\link{wmw.ref.mul}}, \code{\link{wmw.ref.dif}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test at multiple quantiles
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 0.2)
#' wmw.QC.mul(x, y, q = c(0.1, 0.25), nboot = 500)
#' }
wmw.QC.mul<-function(x,y,q=seq(.1,.4,.1), estfun=hdmq, alpha=.05,nboot=1000,SEED=TRUE,method='BH',plotit=FALSE,
xlab='Difference',ylab='Density',
plotfun=kerSORT){
#
#
# For multiple reference quantiles, q>.5,
# make inferences about P(x-y< -refp) vs P(x-y > refp)
#. where refp is the q<.5 quantile
# refp can be constants chosen by the user. If not specified,
# refp are taken to be the .6(.1).9 estimated quantiles of the distribution of X-Y
#
# To use the Harrell-Davis estimator, set estfun=hdmq which is the default
#
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
np=length(q)
output<-matrix(NA,np,8)
dimnames(output)=list(NULL,c('q','q.quant' ,'1-q.quant','Sum','ci.low','ci.up','p.value','p.adj'))
for(i in 1:np){
e=wmw.QC(x,y,q=q[i],alpha=alpha,nboot=nboot,SEED=FALSE)
output[i,1:7]=c(q[i],e$L,e$U,e$Est.dif,e$ci[1],e$ci[2],e$p.value)
}
output[,8]=p.adjust(output[,7],method=method)
if(plotit)plotfun(as.vector(morig),xlab=xlab,ylab=ylab)
output
}

#' Kolmogorov-Smirnov Effect Size for 2x2 Interaction
#'
#' Computes Kolmogorov-Smirnov (KMS) effect sizes for a 2x2 factorial design
#' to assess the interaction effect. Compares effect sizes at two levels of
#' Factor A using specified differences for Factor B.
#'
#' @param x Data in matrix, data frame, or list format containing 4 groups.
#' @param DIF1 Specified difference for Factor B at level 1 of Factor A.
#' @param DIF2 Specified difference for Factor B at level 2 of Factor A.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `KMS.Effect1`: KMS effect size at level 1 of Factor A
#'     \item `KMS.Effect2`: KMS effect size at level 2 of Factor A
#'     \item `Difference`: Difference between the two effect sizes
#'   }
#'
#' @details
#' This function quantifies interaction effects in a 2x2 design using the
#' Kolmogorov-Smirnov effect size measure. For each level of Factor A, it
#' computes the KMS effect size between the two levels of Factor B given
#' the specified difference values.
#'
#' The difference between the two effect sizes provides a measure of the
#' interaction. By common interpretation, KMS effect sizes of .1, .25, and .4
#' are considered small, medium, and large. A difference of .15 or more between
#' the two effect sizes might be viewed as practically important.
#'
#' @seealso \code{\link{kms.effect}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x2 design data
#' set.seed(123)
#' x <- list(
#'   rnorm(30, mean = 0),      # A1B1
#'   rnorm(30, mean = 0.5),    # A1B2
#'   rnorm(30, mean = 0.2),    # A2B1
#'   rnorm(30, mean = 1.0)     # A2B2
#' )
#'
#' # Assess interaction using KMS effect sizes
#' result <- twoway.inter.2.delta(x, DIF1 = 0.5, DIF2 = 0.8)
#' print(result)
#' }
twoway.inter.2.delta<-function(x,DIF1,DIF2){
#
# 2-by-2 interaction
# For each level of Factor A, specify differences  DIF1 and DIF2  for the two levels of Factor B
# determine what the  KMS effect size  for the two leve
#
# The difference between the estimates, what is important depends on the situaiton.
# By a common point of view, .1, .25 and .4 are small, medium and large KMS effect sizes.
# So difference of .15 between the two measures of effect size  might be viewed as being important.
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
x=elimna(x)
if(length(x)!= 4)stop('Should have four groups')
e1=kms.effect(x[[1]],x[[2]],DIF=DIF1)$effect.size
e2=kms.effect(x[[3]],x[[4]],DIF=DIF2)$effect.size
dif=e1-e2
list(KMS.Effect1=e1,KMS.Effect2=e2,Difference=dif)
}

#' Test Symmetry with TOST (Two One-Sided Tests) Extension
#'
#' @description
#' Tests symmetry of X-Y distribution with TOST equivalence testing extension.
#' Can test both for difference from zero and equivalence to zero.
#'
#' @inheritParams common-params
#' @param x Observations from group 1.
#' @param y Observations from group 2.
#' @param q Quantile for reference point if pts=NULL (default: 0.25).
#' @param pts Reference point (default: NULL for quantile-based).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param alpha Significance level (default: 0.05).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param plotit Logical. If `TRUE`, creates plot (default: `TRUE`).
#' @param xlab X-axis label (default: 'Difference').
#' @param ylab Y-axis label (default: 'Density').
#' @param estfun Quantile estimator (default: \code{hdmq}).
#' @param plotfun Plotting function (default: NULL function).
#' @param eqbound Equivalence bound for TOST (default: NULL for no TOST).
#'
#' @return A list with symmetry test results and optionally TOST results:
#'   \itemize{
#'     \item `L`, `U`, `Est.dif`: As in \code{wmw.ref.dif}
#'     \item `ci90`, `ci95`: 90% and 95% confidence intervals
#'     \item `p.value`: P-value for difference test
#'     \item `TOST`: TOST equivalence test results (if eqbound specified)
#'   }
#'
#' @details
#' Extends \code{wmw.ref.dif} with TOST procedure for testing equivalence.
#' Provides both 90% CIs (for TOST) and 95% CIs (for regular inference).
#'
#' @seealso \code{\link{wmw.ref.dif}}, \code{\link{wmw.ref.mul.TOST}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Test with TOST
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 0.05)
#' wmw.ref.dif.TOST(x, y, eqbound = 0.1, nboot = 500)
#' }
wmw.ref.dif.TOST <- function(
    x, y, q = .25, pts = NULL, nboot = 1000, alpha = .05, SEED = TRUE,
    plotit = T, xlab = 'Difference', ylab = 'Density', estfun = hdmq,
    plotfun = function(...)NULL, eqbound = NULL) {
  # TOST extension plus ggplot visual
  
  QC = TRUE
  if (!is.null(pts)) QC = FALSE
  if (is.null(pts)) {
    if (sum(q < .5) != length(q)) stop('All q values should be <=.5')
  }
  if (SEED) set.seed(2)
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n1 = length(x); n2 = length(y)
  
  if (!QC) {
    if (!is.null(pts)) {
      e = wmw.det(x, y, refp = pts, plotit = plotit, xlab = xlab, ylab = ylab, plotfun = plotfun)
      est = e$dif; L = e$L; U = e$U
    } else {
      pts = qest(outer(x, y, FUN = '-'), q = q)
      e = wmw.det(x, y, refp = pts, plotit = plotit, xlab = xlab, ylab = ylab, plotfun = plotfun)
      est = e$dif; L = e$L; U = e$U
    }
  } else {
    d = outer(x, y, FUN = '-')
    d = as.vector(d)
    qv = estfun(d, q = c(q, 1 - q))
    if (plotit) plotfun(d, xlab = xlab, ylab = ylab)
    est = qv[1] + qv[4]; L = qv[5]; U = qv[4]
  }
  
  boot.dif = numeric(nboot)
  for (i in 1:nboot) {
    id1 = sample(n1, replace = TRUE)
    id2 = sample(n2, replace = TRUE)
    if (!QC) {
      boot.dif[i] = wmw.det(x[id1], y[id2], refp = pts, plotit = FALSE)$dif
    } else {
      qv = estfun(outer(x[id1], y[id2], FUN = '-'), q = c(q, 1 - q))
      boot.dif[i] = qv[1] + qv[4]
    }
  }
  
  crit <- alpha# CAUTION! I remove alpha/2 because for TOST we actually need 90% CIs
  icl90 <- round(crit * nboot) + 1
  icu90 <- nboot - icl90
  dif.sorted <- sort(boot.dif)
  ci90 <- c(dif.sorted[icl90], dif.sorted[icu90])
  crit2<-alpha/2#I calculate also 95% CIs for other potential uses of difference and this Ci
  icl95 <- round(crit2 * nboot) + 1
  icu95 <- nboot - icl95
  ci95<-c(dif.sorted[icl95], dif.sorted[icu95])
  pv <- mean(dif.sorted < 0) + .5 * mean(dif.sorted == 0)
  pv <- 2 * min(pv, 1 - pv)
  
  # TOST section
  tost.res <- NULL
  if (!is.null(eqbound)) {

    pval.lower <- mean(boot.dif <= -eqbound) + 0.5 * mean(boot.dif == -eqbound)
    pval.upper <- mean(boot.dif >= eqbound) + 0.5 * mean(boot.dif == eqbound)
    tost.reject <- (pval.lower < alpha) & (pval.upper < alpha)# CAUTION! I use alpha because is unilateral
    tost.res <- list(
      eqbound = eqbound,
      p.value.lower = pval.lower,
      p.value.upper = pval.upper,
      equivalence = tost.reject
    )
  }
  
  # TOST-friendly report
  TOST.details <- NULL
  if (!is.null(eqbound)) {
    TOST.details <- list(
      Est.dif = est,
      ci90=ci90,
      eqbound = eqbound,
      lower_tail = list(
        H0 = paste0("uc0u916  u8804  ", -eqbound),
        rejected = ifelse(tost.res$p.value.lower < alpha, "Yes", "No"),
        p.value = tost.res$p.value.lower
      ),
      upper_tail = list(
        H0 = paste0("uc0u916  u8805  ", eqbound),
        rejected = ifelse(tost.res$p.value.upper < alpha, "Yes", "No"),
        p.value = tost.res$p.value.upper
      ),
      equivalence = ifelse(tost.res$equivalence, "YES", "NO")
    )
  }
  
  # Plot
  plot_obj <- NULL
  if (!is.null(eqbound)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) 
      stop("Package 'ggplot2' is required for plotting.")
    df <- data.frame(
      group = "Estimate",
      estimate = est,
      lower = ci90[1],
      upper = ci90[2]
    )
    equivalence_df <- data.frame(xmin = -eqbound, xmax = eqbound, ymin = -Inf, ymax = Inf)
    
    xlimits=max(abs(ci90))
    xlimits=c(0-xlimits, xlimits+eqbound)
    
    require(ggplot2)
    plot_obj <- ggplot(df, aes(x = estimate, y = 1)) +
      geom_rect(data = equivalence_df, inherit.aes = FALSE,aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax),
                fill = "gray80", alpha = 0.4) +
      geom_segment(data=NULL,aes(x = lower, xend = upper,y=1,yend=1), linewidth = 0.5) +
      geom_point(size = 3, fill = "white",shape=21) +
      labs(y = NULL, x = "Difference",
           title = "Estimate (point & CI) with Equivalence Zone") +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
      scale_x_continuous(limits=c(xlimits[1],xlimits[2]))+
      scale_y_continuous(breaks=NULL, labels=NULL)+
      theme_minimal() 
  }
  
  return(list(
    L = L, U = U, Est.dif = est, ci95=ci95, p.value = pv,
    TOST = TOST.details,
    plot = plot_obj
  ))
}

#' Test Symmetry at Multiple Points with TOST Extension
#'
#' @description
#' Tests symmetry at multiple reference points with TOST equivalence testing.
#' Extension of \code{wmw.ref.mul} with equivalence bounds.
#'
#' @inheritParams common-params
#' @param x Observations from group 1.
#' @param y Observations from group 2.
#' @param refp Vector of reference points (default: NULL for quantile-based).
#' @param pts Alias for refp (default: NULL).
#' @param q Quantiles if refp=NULL (default: seq(0.6, 0.9, 0.1)).
#' @param center Logical. If `TRUE`, centers at median (default: `FALSE`).
#' @param estfun Quantile estimator (default: \code{hdmq}).
#' @param alpha Significance level (default: 0.05).
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param SEED Logical. If `TRUE`, sets random seed (default: `TRUE`).
#' @param method Multiple comparison adjustment (default: 'BH').
#' @param plotit Logical. If `TRUE`, plots distribution (default: `FALSE`).
#' @param xlab X-axis label (default: 'Difference').
#' @param ylab Y-axis label (default: 'Density').
#' @param plotfun Plotting function (default: \code{kerSORT}).
#' @param eqbound Equivalence bound for TOST (default: 0.05).
#'
#' @return Data frame with test results at each reference point including:
#'   90% CIs, 95% CIs, p-values, adjusted p-values, and TOST equivalence results.
#'
#' @details
#' Performs \code{wmw.ref.dif.TOST} at multiple reference points with
#' family-wise error control. Provides comprehensive output including
#' both difference testing and equivalence testing results.
#'
#' @seealso \code{\link{wmw.ref.dif.TOST}}, \code{\link{wmw.ref.mul}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Multiple points with TOST
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 0.05)
#' wmw.ref.mul.TOST(x, y, eqbound = 0.1, nboot = 500)
#' }
wmw.ref.mul.TOST <- function(x, y, refp = NULL, pts = NULL, q = seq(.6, .9, .1), center = FALSE,estfun = hdmq,
                             alpha = .05, nboot = 1000, SEED = TRUE, method = 'BH',plotit = FALSE, xlab = 'Difference', ylab = 'Density', plotfun = kerSORT,
                             eqbound = 0.05) {
  if (SEED) set.seed(2)
  
  if (!is.null(pts)) refp <- pts
  
  x <- elimna(x)
  y <- elimna(y)
  
  if (is.null(refp)) {
    m <- outer(x, y, FUN = '-')
    m <- as.vector(m)
    if (center) m <- m - median(m)
    refp <- estfun(m, q)
  }
  
  np <- length(refp)
  output <- data.frame(
    Pts = numeric(np),
    P_less = numeric(np),
    P_greater = numeric(np),
    Dif = numeric(np),
    ci90.low = numeric(np),
    ci90.up = numeric(np),
    ci95.low = numeric(np),
    ci95.up = numeric(np),
    p.value = numeric(np),
    p.adj = numeric(np),
    TOST_equivalent = character(np),
    TOST_lower_pval=numeric(np),
    TOST_upper_pval=numeric(np),
    TOST_pval=numeric(np))
  colnames(output) <- c('Pts', 'P(x-y<-Pts)', 'P(x-y>Pts)', 'Dif', 'ci90.low', 'ci90.up','ci95.low', 'ci95.up',
                        'p.value', paste0('p.adj.', method), 'TOST.equivalent','TOST_lower_pval','TOST_upper_pval','TOST_pval')
  
  for (i in 1:np) {
    res <- wmw.ref.dif.TOST(
      x, y, q = q[1], pts = refp[i], nboot = nboot, alpha = alpha, SEED = FALSE,
      plotit = FALSE, xlab = xlab, ylab = ylab, estfun = estfun, plotfun = plotfun,
      eqbound = eqbound)
    
    output[i, 1] <- refp[i]
    output[i, 2] <- res$L
    output[i, 3] <- res$U
    output[i, 4] <- res$Est.dif
    output[i, 5] <- res$TOST$ci90[1]
    output[i, 6] <- res$TOST$ci90[2]
    output[i, 7] <- res$ci95[1]
    output[i, 8] <- res$ci95[2]
    output[i, 9] <- res$p.value
    output[i, 10] <- res$p.value
    output[i, 11] <- res$TOST$equivalence
    output[i, 12] <-  res$TOST$lower_tail$p.value
    output[i, 13] <-  res$TOST$upper_tail$p.value
    output[i, 14] <-  max(c(res$TOST$lower_tail$p.value, res$TOST$upper_tail$p.value))
  }
  output[, 10] <- p.adjust(output[, 9],method=method)
  return(as.data.frame(output))
}


# ==== pb2genMC ====
#' Percentile Bootstrap for Two-Group Comparisons (Multicore Version)
#'
#' Computes a percentile bootstrap confidence interval for the difference
#' between any two parameters for independent groups. Uses multicore processing
#' via `mclapply` for improved performance. This is the parallel version of `pb2gen`.
#'
#' @inheritParams common-params
#' @param pr Logical. If `TRUE`, prints progress message (default: `FALSE`).
#' @param ... Additional arguments passed to the estimator function `est`.
#'
#' @details
#' This function is identical to `pb2gen` except that it uses `mclapply` from
#' the `parallel` package to process bootstrap samples in parallel across
#' multiple CPU cores. This can provide substantial speedups, especially with
#' large datasets or many bootstrap samples.
#'
#' The function is a general-purpose two-sample bootstrap that can compare any
#' parameter between independent groups. By default, it compares M-estimators
#' of location (`est = onestep`).
#'
#' Common choices for `est`:
#' \itemize{
#'   \item `onestep` - One-step M-estimator (default)
#'   \item `mean` - Arithmetic mean
#'   \item `median` - Median
#'   \item `tmean` - Trimmed mean (use with `trim = 0.2`)
#'   \item `mom` - Modified one-step M-estimator
#'   \item `hd` - Harrell-Davis estimator (use with `q = 0.5` for median)
#' }
#'
#' Missing values are automatically removed. The function constructs a
#' percentile bootstrap confidence interval and computes a p-value for
#' testing whether the difference equals zero.
#'
#' **Performance Note**: Multicore processing is most beneficial when:
#' \itemize{
#'   \item `nboot` is large (e.g., >= 2000)
#'   \item The estimator function is computationally intensive
#'   \item Multiple CPU cores are available
#' }
#'
#' On single-core systems, `pb2gen` may be faster due to reduced overhead.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `est.1`: Estimate for group 1
#'     \item `est.2`: Estimate for group 2
#'     \item `ci`: Two-element vector with confidence interval bounds
#'     \item `p.value`: P-value for testing H0: difference = 0
#'     \item `sq.se`: Squared standard error of the bootstrap distribution
#'     \item `n1`, `n2`: Sample sizes for groups 1 and 2
#'   }
#'
#' @note
#' The `parallel` package must be loaded for multicore functionality.
#' Results may vary slightly across runs even with `SEED = TRUE` due to
#' the random number generation in parallel processing.
#'
#' @seealso \code{\link{pb2gen}} for single-core version,
#'   \code{\link{trimpb2}} for trimmed means specifically,
#'   \code{\link{yuenbt}} for bootstrap-t method
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(60, mean = 0.5)
#'
#' # Compare M-estimators (default) using multicore
#' pb2genMC(x, y, nboot = 4000)
#'
#' # Compare means
#' pb2genMC(x, y, est = mean, nboot = 4000)
#'
#' # Compare medians
#' pb2genMC(x, y, est = median, nboot = 4000)
#'
#' # Compare trimmed means
#' pb2genMC(x, y, est = mean, trim = 0.2, nboot = 4000)
#'
#' # Compare specific quantiles using Harrell-Davis
#' pb2genMC(x, y, est = hd, q = 0.75, nboot = 4000)  # Third quartile
#' }
 pb2genMC<-function(x,y,alpha=.05,nboot=2000,est=onestep,SEED=TRUE,pr=FALSE,...){
#
#   Compute a bootstrap confidence interval for the
#   the difference between any two parameters corresponding to
#   independent groups.
#   By default, M-estimators are compared.
#   Setting est=mean, for example, will result in a percentile
#   bootstrap confidence interval for the difference between means.
#   Setting est=onestep will compare M-estimators of location.
#   The default number of bootstrap samples is nboot=2000
#
x<-x[!is.na(x)] # Remove any missing values in x
y<-y[!is.na(y)] # Remove any missing values in y
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
#
datax=t(datax)
datay=t(datay)
datax=listm(datax)
datay=listm(datay)
bvecx<-mclapply(datax,est,mc.preschedule=TRUE,...)
bvecy<-mclapply(datay,est,mc.preschedule=TRUE,...)
bvec=sort(matl(bvecx)-matl(bvecy))
low<-round((alpha/2)*nboot)+1
up<-nboot-low
temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
sig.level<-2*(min(temp,1-temp))
se<-var(bvec)
list(est.1=est(x,...),est.2=est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}

# ==== pb2gen ====
#' Percentile Bootstrap for Two-Group Comparisons (General Purpose)
#'
#' Computes a percentile bootstrap confidence interval for the difference
#' between any two parameters for independent groups. This is a general-purpose
#' two-sample bootstrap function.
#'
#' @inheritParams common-params
#' @param ... Additional arguments passed to the estimator function `est`.
#'
#' @details
#' This is a flexible function that can compare any parameter between two
#' independent groups using the percentile bootstrap method. By default, it
#' compares M-estimators of location (`est = onestep`).
#'
#' Common choices for `est`:
#' \itemize{
#'   \item `onestep` - One-step M-estimator (default)
#'   \item `mean` - Arithmetic mean
#'   \item `median` - Median
#'   \item `tmean` - Trimmed mean (use with `trim = 0.2`, for example)
#'   \item `mom` - Modified one-step M-estimator
#'   \item `hd` - Harrell-Davis estimator (use with `q = 0.5` for median)
#' }
#'
#' Missing values are automatically removed. The function uses `nboot` bootstrap
#' samples (default: 2000) to construct the confidence interval and compute the
#' p-value.
#'
#' @return A list with components:
#'   \itemize{
#'     \item `est.1`: Estimate for group 1
#'     \item `est.2`: Estimate for group 2
#'     \item `est.dif`: Estimated difference (est.1 - est.2)
#'     \item `ci`: Two-element vector with confidence interval bounds
#'     \item `p.value`: P-value for testing H0: difference = 0
#'     \item `sq.se`: Squared standard error of the bootstrap distribution
#'     \item `n1`, `n2`: Sample sizes for groups 1 and 2
#'   }
#'
#' @seealso \code{\link{pb2genMC}}, \code{\link{trimpb2}}, \code{\link{yuenbt}}
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(50)
#' y <- rnorm(60, mean = 0.5)
#'
#' # Compare M-estimators (default)
#' pb2gen(x, y)
#'
#' # Compare means
#' pb2gen(x, y, est = mean)
#'
#' # Compare medians
#' pb2gen(x, y, est = median)
#'
#' # Compare trimmed means
#' pb2gen(x, y, est = mean, trim = 0.2)
#'
#' # Compare specific quantiles using Harrell-Davis
#' pb2gen(x, y, est = hd, q = 0.25)  # First quartile
#' }
pb2gen<-function(x,y,alpha=.05,nboot=2000,est=onestep,SEED=TRUE,pr=FALSE,...){
#
#   Compute a bootstrap confidence interval for the
#   the difference between any two parameters corresponding to
#   independent groups.
#   By default, M-estimators are compared.
#   Setting est=mean, for example, will result in a percentile
#   bootstrap confidence interval for the difference between means.
#   Setting est=onestep will compare M-estimators of location.
#   The default number of bootstrap samples is nboot=2000
#
x<-x[!is.na(x)] # Remove any missing values in x
y<-y[!is.na(y)] # Remove any missing values in y
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvecx<-apply(datax,1,est,...)
bvecy<-apply(datay,1,est,...)
bvec<-sort(bvecx-bvecy)
low<-round((alpha/2)*nboot)+1
up<-nboot-low
temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
sig.level<-2*(min(temp,1-temp))
se<-var(bvec)
list(est.1=est(x,...),est.2=est(y,...),est.dif=est(x,...)-est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}
