# ANCOVA Functions
# Analysis of Covariance methods for robust statistics
# Extracted from WRS package v0.45

# This module contains:
# - Core ANCOVA functions (independent groups)
# - Dependent/repeated measures ANCOVA
# - Bootstrap and permutation methods
# - Effect size computation
# - Global tests and multiple comparisons
# - Regression-based ANCOVA
# - Plotting and visualization
# - Detection and diagnostics

# anc2COV.CV
# difqci.mul (dependency - for tailsci.mul alias)

#' Difference Quantile Confidence Intervals for Multiple Reference Points
#'
#' Computes confidence intervals for P(X-Y < refp) at multiple reference points.
#' Uses the Harrell-Davis estimator by default to estimate quantiles of the
#' distribution of X-Y differences.
#'
#' @inheritParams common-params
#' @param refp Reference points for the probability estimates. If NULL, uses quantile
#'   estimates based on the q parameter.
#' @param pts Alternative specification for reference points (synonymous with refp).
#' @param q Quantile values to use for automatic reference point selection
#'   (default: seq(.1,.4,.1)).
#' @param estfun Estimator function for quantiles (default: hdmq, Harrell-Davis).
#' @param center Logical. If TRUE, centers the differences by subtracting the median.
#' @param method Method for p-value adjustment (default: 'BH' for Benjamini-Hochberg).
#' @param xlab Label for x-axis when plotit=TRUE.
#' @param ylab Label for y-axis when plotit=TRUE.
#' @param plotfun Plotting function to use (default: kerSORT).
#'
#' @return Matrix with columns for quantiles, reference points, probability estimates,
#'   and confidence interval bounds.
#'
#' @export
#' @examples
#' \dontrun{
#' x <- rnorm(50, mean=5)
#' y <- rnorm(50, mean=3)
#' difqci.mul(x, y, nboot=500)
#' }
difqci.mul<-function(x,y,refp=NULL,pts=NULL,q=seq(.1,.4,.1), estfun=hdmq, center=FALSE, alpha=.05,nboot=1000,SEED=TRUE,method='BH',plotit=FALSE,
xlab='Difference',ylab='Density',
plotfun=kerSORT){
#
#
# For multiple reference points, refp,
# make inferences about P(x-y< refp)
# refp can be constants chosen by the user. If not specified,
# refp are taken to be the .6(.1).9 estimated quantiles of the distribution of X-Y
#
# pts can be used to indicate specified reference points, refp
#
# If refp=NULL, reference points are based on estimates of the q quantiles
# Default is q=.1(.1).4
#
# Her,  use the Harrell-Davis estimator by default,  estfun=hdmq
#
if(SEED)set.seed(2)
if(!is.null(pts))refp=pts
x=elimna(x)
y=elimna(y)
FLAG=TRUE
m=outer(x,y,FUN='-')
m=as.vector(m)
morig=m

if(center)m=m-median(m)
if(is.null(refp)){
FLAG=FALSE
refp=estfun(m,q)
}
np=length(refp)
if(!FLAG){
output<-matrix(NA,np,5)
dimnames(output)=list(NULL,c('q','Pts','P(x-y<Pts)' ,'ci.low','ci.up'))
output[,1]=q
for(i in 1:np){
e=difqci(x,y,pts=refp[i],alpha=alpha,nboot=nboot,SEED=FALSE)
output[i,2:5]=c(refp[i],e$Est,e$ci[1],e$ci[2])
}}
if(FLAG){
output<-matrix(NA,np,4)
dimnames(output)=list(NULL,c('Pts','P(x-y<Pts)' ,'ci.low','ci.up'))
for(i in 1:np){
e=difqci(x,y,pts=refp[i],alpha=alpha,nboot=nboot,SEED=FALSE)
output[i,1:4]=c(refp[i],e$Est,e$ci[1],e$ci[2])
}}
if(plotit)plotfun(as.vector(morig),xlab=xlab,ylab=ylab)
output
}

tailci.mul=difqci.mul


#' Effect Size Function - Multiple Methods
#'
#' Computes effect sizes for two-group comparisons using various robust methods.
#' Provides a unified interface to multiple effect size measures including explanatory
#' power, quantile shift, and robust analogues of Cohen's d.
#'
#' @inheritParams common-params
#' @param QSfun Location function for quantile shift method (default: median).
#' @param method Character string specifying the effect size method:
#'   \itemize{
#'     \item 'EP': Explanatory power based on trimmed means
#'     \item 'QS': Quantile shift based on medians
#'     \item 'QStr': Quantile shift based on trimmed means
#'     \item 'AKP': Robust analog of Cohen's d (Algina-Keselman-Penfield)
#'     \item 'WMW': P(X<Y) - Wilcoxon-Mann-Whitney type measure
#'     \item 'KMS': Heteroscedastic analog of Cohen's d
#'   }
#'
#' @return Numeric value representing the effect size.
#'
#' @keywords internal
ESfun<-function(x,y,QSfun=median,method=c('EP','QS','QStr','AKP','WMW','KMS'),tr=.2,pr=TRUE,SEED=TRUE){
type=match.arg(method)
switch(type,
    EP=yuenv2(x,y,tr=tr,SEED=SEED)$Effect.Size,         #Explanatory power
    QS=shiftQS(x,y,locfun=QSfun)$Q.Effect,   #Quantile shift based on the medians
    QStr=yuenQS(x,y,tr=tr,pr=pr)$Q.Effect,       #Based on trimmed means
    AKP=akp.effect(x,y,tr=tr),                      #Robust analog of Cohen's d
    WMW=pxly(x,y,SEED=SEED),          #  P(X<Y)
    KMS=kms.effect(x,y,tr=tr)$effect.size)  #  Heteroscedastic analog of Cohen's d
}

#' Normal Distribution Variance Helper
#'
#' Internal helper function that computes x^2 * dnorm(x).
#'
#' @param x Numeric vector.
#' @return Numeric vector of transformed values.
#' @keywords internal
dnormvar<-function(x){
x^2*dnorm(x)
}


#' Critical Values for ANCOVA with Two Covariates
#'
#' Determines the null distribution of the test statistic used by ancov2COV.
#' Uses simulation to generate critical values for two-group ANCOVA with
#' two covariates.
#'
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param iter Number of iterations for simulation (default: 1000).
#' @inheritParams common-params
#' @param FRAC Fraction parameter for depth-based methods (default: 0.5).
#' @param TPM Logical. If TRUE, uses truncated product method.
#' @param tau Threshold for TPM method (default: 0.05).
#'
#' @return List with critical value and vector of simulated test statistics.
#'
#' @export
#' @examples
#' \dontrun{
#' anc2COV.CV(30, 30, iter=500)
#' }
anc2COV.CV<-function(n1,n2,iter=1000,alpha=.05,FRAC=.5,SEED=TRUE,MC=FALSE,
TPM=FALSE,tau=.05){
#
# Determine null distribution of the test statistic
# used by ancov2COV.
#
if(SEED)set.seed(2)
n=max(c(n1,n2))
data=list()
for(i in 1:iter){
x=rmul(n,p=6)
N=n-min(c(n1,n2))
if(N>0)x[1:N,1:3]=NA
data[[i]]=x
}
if(MC){
res=mclapply(data,anc2COV.sub,FRAC=FRAC,TPM=TPM,tau=tau)
}
if(!MC)res=lapply(data,anc2COV.sub,FRAC=FRAC,TPM=TPM,tau=tau)
M=as.vector(list2mat(res))
M=sort(M)
ic=round(alpha*iter)
crit=M[ic]
list(crit.val=crit,M=M)
}


#' Helper Function for anc2COV.CV
#'
#' Internal subroutine for computing test statistics in anc2COV.CV simulations.
#'
#' @param data Simulated data matrix.
#' @param FRAC Fraction parameter for depth-based methods.
#' @param TPM Logical. If TRUE, uses truncated product method.
#' @param tau Threshold for TPM method.
#'
#' @return Test statistic value.
#' @keywords internal
anc2COV.sub<-function(data,FRAC=.5,TPM=FALSE,tau=.05){
val=ancovampG(data[,1:2],data[,3],data[,4:5],data[,6],DH=TRUE,SEED=FALSE,test=yuen,FRAC=FRAC,cov.fun=covl)$results[,3]
val=elimna(val)
if(!TPM)val=mean(val,na.rm=TRUE)
if(TPM)val=sum(log(val[val<=tau]))
val
}

#' ANCOVA for Binary Outcomes - Two Groups
#'
#' Compares probability of success for two independent groups across covariate values.
#' Uses a running-interval smoother coupled with binomial comparison methods.
#' Provides familywise error rate control through p-value adjustment.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Binary outcome (0/1) for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Binary outcome (0/1) for group 2.
#' @inheritParams common-params
#' @param FWE.alpha Familywise error rate for significance testing (default: 0.05).
#' @param pts Design points for comparison. If NA, uses npts equally spaced points.
#' @param fr1 Span for running interval smoother for group 1 (default: 0.8).
#' @param fr2 Span for running interval smoother for group 2 (default: 0.8).
#' @param npts Number of design points if pts=NA (default: 5).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param nmin Minimum sample size required in each interval (default: 12).
#' @param method P-value adjustment method (default: 'BH' for Benjamini-Hochberg).
#' @param bin.method Binomial comparison method: 'KMS', 'ECP', 'ZHZ', or 'SK'.
#' @param reps Number of replications for critical value estimation (default: 5000).
#'
#' @return List with components:
#'   \item{n1}{Sample size for group 1}
#'   \item{n2}{Sample size for group 2}
#'   \item{FWE.crit}{Critical value for familywise error rate control}
#'   \item{output}{Matrix with results at each design point}
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(100)
#' y1 <- rbinom(100, 1, plogis(x1))
#' x2 <- rnorm(100)
#' y2 <- rbinom(100, 1, plogis(x2 - 0.5))
#' anc.2gbin(x1, y1, x2, y2, npts=5)
#' }
anc.2gbin<-function(x1,y1,x2,y2,alpha=.05,FWE.alpha=.05,
pts=NA,fr1=.8,fr2=.8,npts=5,xlab='X',ylab='Est. Dif',xout=FALSE,
outfun=out,nmin=12,plotit=TRUE,method='BH',bin.method='KMS',SEED=TRUE,reps=5000,pr=TRUE){
#
# Compare probability of success for two independent  groups given a value for some covariate. 
# A running-interval smoother is used coupled with the technique indicated by the argument  
#  bin.method. See R function binom2g
# 
# method: how to adjust the p-value to control FWE, see R function p.adjust
# default is Benjamini--Hochberg.  method='hoch' is Hochberg
#
# pts=NULL, use npts covariate values, default is 10, spread over the values in x1.
#     
# 
#
xy=elimna(cbind(x1,y1))
if(pr){
if(nrow(xy)<36)print('With a sample size less 36, computational errors might occur')
}
x1=xy[,1]
y1=xy[,2]

xy=elimna(cbind(x2,y2))
if(pr){
if(nrow(xy)<36)print('With a sample size less 36, computational errors might occur')
}

x2=xy[,1]
y2=xy[,2]

isub=0
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
N1=NA
N2=NA
xorder=order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
if(is.na(pts[1])){
n1<-1
n2<-1
N1=NA
N2=NA
vecn<-1
isub=0
for(i in 1:length(x1))n1[i]<-length(y1[near(x1,x1[i],fr1)])
for(i in 1:length(x1))n2[i]<-length(y2[near(x2,x1[i],fr2)])
for(i in 1:length(x1))vecn[i]<-min(n1[i],n2[i])
sub<-c(1:length(x1))
isub[1]<-min(sub[vecn>=nmin])
isub[2]<-max(sub[vecn>=nmin])
bot=x1[isub[1]]
top=x1[isub[2]]
pts=seq(bot,top,length.out=npts)
}
if(bin.method!='SK')output=matrix(NA,nrow=length(pts),ncol=8)
else output=matrix(NA,nrow=length(pts),ncol=6)

for(i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
N1[i]=length(g1)
N2[i]=length(g2)
v=binom2g(x=g1,y=g2,method=bin.method,alpha=alpha,SEED=SEED) #correct when using ECP
#                     because only using this function for two groups.
if(bin.method=='KMS' || bin.method=='ECP' ){
output[i,1:7]=c(pts[i],v$p1,v$p2,v$est.dif,v$ci[1],v$ci[2],v$p.value)
}
if(bin.method=='ZHZ')output[i,1:7]=c(pts[i],v$p.hat[1],v$p.hat[2],v$CI[1],v$CI[2],v$CI[3],v$CI[4])
if(bin.method=='ECP')output[i,1:7]=c(pts[i],v$output[1,3],v$output[1,4],v$output[1,5],
v$output[1,6],
v$output[1,7],v$output[1,8])
if(bin.method=='SK')output[i,1:5]=c(pts[i],v$p1,v$p2,v$est.dif,v$p.value)
}
if(bin.method!='SK')
dimnames(output)=list(NULL,c('pts','p1','p2','est.dif','ci.lower','ci.upper','p.value','p.adjusted'))
else
dimnames(output)=list(NULL,c('pts','p1','p2','est.dif','p.value','p.adjusted'))
if(plotit){
plot(c(pts,pts,pts),c(output[,4],output[,5],output[,6]),type='n',xlab=xlab,ylab=ylab,ylim=c(-1,1))
points(pts,output[,4])
if(bin.method!='SK'){
points(pts,output[,5],pch='+')
points(pts,output[,6],pch='+')
}
}
pcritFWE=NULL
if(bin.method=='ECP'){
p=(output[,2]+output[,3])/2
pcritFWE=binKMS.mpts.crit(p,N1,N2,iter=reps,SEED=SEED,FWE=FWE.alpha)
}
if(bin.method!='SK')output[,8]=p.adjust(output[,7],method=method)
else
output[,6]=p.adjust(output[,5],method=method)
list(n1=length(y1),n2=length(y2),FWE.crit=pcritFWE,output=output)
}






#' ANCOVA Bootstrap Bagging with Median
#'
#' Compares two independent groups using ANCOVA based on bootstrap bagging with
#' running interval smoother. Uses median as the estimator for the dependent variable.
#' No parametric assumptions about regression lines.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Dependent variable for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Dependent variable for group 2.
#' @param fr1 Span for running interval smoother for group 1 (default: 1).
#' @param fr2 Span for running interval smoother for group 2 (default: 1).
#' @inheritParams common-params
#' @param pts Design points for comparison. If NA, uses 5 empirically chosen points.
#'
#' @return List with:
#'   \item{output}{Matrix with design points, sample sizes, differences, CIs, and p-values}
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50)
#' y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- x2 + 0.5 + rnorm(50)
#' ancbbmed(x1, y1, x2, y2, nboot=500)
#' }
ancbbmed<-function(x1,y1,x2,y2,fr1=1,fr2=1,nboot=100,pts=NA,plotit=TRUE,
SEED=TRUE,alpha=.05){
#
# Compare two independent  groups using an ancova method
# based in part on a bootstrap bagging estimate of the dependent variable.
# No assumption is made about the form of the regression
# lines--a running interval smoother is used.
# Confidence intervals are computed using a percentile bootstrap
# method. Comparisons are made at five empirically chosen design points.
#
#  Assume data are in x1 y1 x2 and y2
#
if(SEED)set.seed(2)
if(is.na(pts[1])){
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
mat<-matrix(NA,5,7)
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","ci.low","ci.hi","p.value"))
gv1<-vector("list")
for (i in 1:5){
j<-i+5
temp1<-y1[near(x1,x1[isub[i]],fr1)]
temp2<-y2[near(x2,x1[isub[i]],fr2)]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(temp1)
mat[i,3]<-length(temp2)
mat[,4]<-runmbo(x1,y1,pts=x1[isub],pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=tmean)-
runmbo(x2,y2,pts=x1[isub],pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=median)
gv1[[i]]<-temp1
gv1[[j]]<-temp2
}
I1<-diag(5)
I2<-0-I1
con<-rbind(I1,I2)
estmat1<-matrix(nrow=nboot,ncol=length(isub))
estmat2<-matrix(nrow=nboot,ncol=length(isub))
data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
#
for(ib in 1:nboot){
estmat1[ib,]=runmbo(x1[data1[ib,]],y1[data1[ib,]],pts=x1[isub],
pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=median)
estmat2[ib,]=runmbo(x2[data2[ib,]],y2[data2[ib,]],pts=x1[isub],
pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=median)
}
dif<-(estmat1<estmat2)
dif0<-(estmat1==estmat2)
pvals=apply(dif,2,mean,na.rm=TRUE)+.5*apply(dif0,2,mean,na.rm=TRUE)
tmat<-rbind(pvals,1-pvals)
pvals=2*apply(tmat,2,min)
mat[,7]<-pvals
for(ij in 1:length(isub)){
dif<-estmat1[,ij]-estmat2[,ij]
dif<-elimna(dif)
nbad<-length(dif)
lo<-round(nbad*alpha/2)
hi<-nbad-lo
dif<-sort(dif)
mat[ij,5]<-dif[lo]
mat[ij,6]<-dif[hi]
}
}
if(!is.na(pts[1])){
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
if(n1[i]<=5)print(paste("Warning, there are",n1[i]," points corresponding to the design point X=",pts[i]))
if(n2[i]<=5)print(paste("Warning, there are",n2[i]," points corresponding to the design point X=",pts[i]))
}
mat<-matrix(NA,length(pts),7)
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","ci.low","ci.hi",
"p.value"))
gv<-vector("list",2*length(pts))
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
j<-i+length(pts)
gv[[i]]<-g1
gv[[j]]<-g2
}
I1<-diag(length(pts))
I2<-0-I1
con<-rbind(I1,I2)
isub=c(1:length(pts))
estmat1<-matrix(nrow=nboot,ncol=length(isub))
estmat2<-matrix(nrow=nboot,ncol=length(isub))
data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
est1=runmbo(x1,y1,pts=pts,pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=median)
est2=runmbo(x2,y2,pts=pts,pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=median)
mat[,4]<-est1-est2
for(ib in 1:nboot){
estmat1[ib,]=runmbo(x1[data1[ib,]],y1[data1[ib,]],pts=pts,pyhat=TRUE,plotit=FALSE,
SEED=FALSE,est=median)
estmat2[ib,]=runmbo(x2[data2[ib,]],y2[data2[ib,]],pts=pts,pyhat=TRUE,plotit=FALSE,
SEED=FALSE,est=median)
}
dif<-(estmat1<estmat2)
dif0<-(estmat1==estmat2)
pvals=apply(dif,2,mean,na.rm=TRUE)+.5*apply(dif0,2,mean,na.rm=TRUE)
tmat<-rbind(pvals,1-pvals)
pvals=2*apply(tmat,2,min)
#
mat[,1]<-pts
mat[,2]<-n1
mat[,3]<-n2
mat[,7]<-pvals
for(ij in 1:length(pts)){
dif<-sort(estmat1[,ij]-estmat2[,ij])
dif<-elimna(dif)
nbad<-length(dif)
lo<-round(nbad*alpha/2)
hi<-nbad-lo
mat[ij,5]<-dif[lo]
mat[ij,6]<-dif[hi]
}
}
if(plotit)
runmean2g(x1,y1,x2,y2,fr=fr1,est=median,sm=TRUE)
list(output=mat)
}







#' ANCOVA Bootstrap Bagging with Percentile Bootstrap
#'
#' Compares two independent groups using ANCOVA with bootstrap bagging and
#' running interval smoother. Provides more general estimator options than ancbbmed.
#' Computes percentile bootstrap confidence intervals.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Dependent variable for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Dependent variable for group 2.
#' @param fr1 Span for running interval smoother for group 1 (default: 1).
#' @inheritParams common-params
#' @param fr2 Span for running interval smoother for group 2 (default: 1).
#' @param pts Design points for comparison. If NA, uses 5 empirically chosen points.
#' @param SCAT Logical. If FALSE, suppresses scatterplot (default: TRUE).
#' @param pch1 Plotting character for group 1 (default: '+').
#' @param pch2 Plotting character for group 2 (default: 'o').
#' @param RNA Logical. If TRUE, removes NA values when computing bagged estimates.
#' @param sm Logical. If FALSE, doesn't use bootstrap bagging when plotting.
#' @param LP Logical. If TRUE, keeps leverage points when plotting.
#'
#' @return List with:
#'   \item{output}{Matrix with design points, sample sizes, differences, CIs, p-values, and critical p-values}
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50)
#' y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- x2 + 0.5 + rnorm(50)
#' ancbbpb(x1, y1, x2, y2, nboot=500)
#' }
ancbbpb<-function(x1,y1,x2,y2,fr1=1,est=tmean,fr2=1,nboot=200,pts=NA,plotit=TRUE,SCAT=TRUE,
pch1='+',pch2='o',
SEED=TRUE,alpha=.05,RNA=TRUE,sm=FALSE,LP=TRUE,xout=FALSE,outfun=outpro,...){
#
# Compare two independent  groups using an ancova method.
# A running-interval smooth is used to estimate the regression lines and is
# based in part on bootstrap bagging.
#
#  This function is limited to two groups and one covariate.
#
# No assumption is made about the parametric form of the regression
# lines.
# Confidence intervals are computed using a percentile bootstrap
# method. Comparisons are made at five empirically chosen design points when
# pts=NA. To compare groups at specified x values, use pts.
# Example: pts=c(60,70,80) will compare groups at the three design points
# 60, 70 and 80.
#
#   xout=F, when plotting, keep leverage points
#   sm=F, when plotting, do not use bootstrap bagging
#
#  Assume data are in x1 y1 x2 and y2
#
# fr1 and fr2 are the spans used by the smooth.
#
#  SCAT=FALSE will suppress the scatterplot when plotting the regression lines.
#
# RNA=F, when computing bagged estimate, NA values are not removed
#  resulting in no estimate of Y at the specified design point,
#  RNA=T, missing values are removed and the remaining values are used.
#
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
#
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
if(SEED)set.seed(2)
flag=TRUE
if(is.na(pts[1])){
flag=FALSE
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
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","ci.low","ci.hi","p.value","p.crit"))
gv1<-vector("list")
for (i in 1:5){
j<-i+5
temp1<-y1[near(x1,x1[isub[i]],fr1)]
temp2<-y2[near(x2,x1[isub[i]],fr2)]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(temp1)
mat[i,3]<-length(temp2)
mat[,4]<-runmbo(x1,y1,pts=x1[isub],pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=est,RNA=RNA)-
runmbo(x2,y2,pts=x1[isub],pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=est,RNA=RNA)
gv1[[i]]<-temp1
gv1[[j]]<-temp2
}
I1<-diag(5)
I2<-0-I1
con<-rbind(I1,I2)
estmat1<-matrix(nrow=nboot,ncol=length(isub))
estmat2<-matrix(nrow=nboot,ncol=length(isub))
data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
#
for(ib in 1:nboot){
estmat1[ib,]=runmbo(x1[data1[ib,]],y1[data1[ib,]],pts=x1[isub],
pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=est,...)
estmat2[ib,]=runmbo(x2[data2[ib,]],y2[data2[ib,]],pts=x1[isub],
pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=est,...)
}
dif<-(estmat1<estmat2)
dif0<-(estmat1==estmat2)
pvals=apply(dif,2,mean,na.rm=TRUE)+.5*apply(dif0,2,mean,na.rm=TRUE)
tmat<-rbind(pvals,1-pvals)
pvals=2*apply(tmat,2,min)
mat[,7]<-pvals
for(ij in 1:length(isub)){
dif<-estmat1[,ij]-estmat2[,ij]
dif<-elimna(dif)
nbad<-length(dif)
lo<-round(nbad*alpha/2)
hi<-nbad-lo
dif<-sort(dif)
mat[ij,5]<-dif[lo]
mat[ij,6]<-dif[hi]
}
}
if(!is.na(pts[1])){
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
if(n1[i]<=5)print(paste("Warning, there are",n1[i]," points corresponding to the design point X=",pts[i]))
if(n2[i]<=5)print(paste("Warning, there are",n2[i]," points corresponding to the design point X=",pts[i]))
}
mat<-matrix(NA,length(pts),7)
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","ci.low","ci.hi",
"p.value"))
gv<-vector("list",2*length(pts))
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
j<-i+length(pts)
gv[[i]]<-g1
gv[[j]]<-g2
}
I1<-diag(length(pts))
I2<-0-I1
con<-rbind(I1,I2)
isub=c(1:length(pts))
estmat1<-matrix(nrow=nboot,ncol=length(isub))
estmat2<-matrix(nrow=nboot,ncol=length(isub))
data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
est1=runmbo(x1,y1,pts=pts,pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=est,...)
est2=runmbo(x2,y2,pts=pts,pyhat=TRUE,plotit=FALSE,SEED=FALSE,est=est,...)
mat[,4]<-est1-est2
for(ib in 1:nboot){
estmat1[ib,]=runmbo(x1[data1[ib,]],y1[data1[ib,]],pts=pts,pyhat=TRUE,plotit=FALSE,
SEED=FALSE,est=est,...)
estmat2[ib,]=runmbo(x2[data2[ib,]],y2[data2[ib,]],pts=pts,pyhat=TRUE,plotit=FALSE,
SEED=FALSE,est=est,...)
}
dif<-(estmat1<estmat2)
dif0<-(estmat1==estmat2)
pvals=apply(dif,2,mean,na.rm=TRUE)+.5*apply(dif0,2,mean,na.rm=TRUE)
tmat<-rbind(pvals,1-pvals)
pvals=2*apply(tmat,2,min)
#
mat[,1]<-pts
mat[,2]<-n1
mat[,3]<-n2
mat[,7]<-pvals
for(ij in 1:length(pts)){
dif<-sort(estmat1[,ij]-estmat2[,ij])
dif<-elimna(dif)
nbad<-length(dif)
lo<-round(nbad*alpha/2)
hi<-nbad-lo
mat[ij,5]<-dif[lo]
mat[ij,6]<-dif[hi]
}
}
temp2<-order(0-pvals)
zvec=alpha/c(1:length(pvals))
if(flag)mat[temp2,7]=zvec
if(!flag)mat[temp2,8]=zvec
if(plotit)
runmean2g(x1,y1,x2,y2,fr=fr1,est=est,sm=sm,LP=LP,xout=FALSE,SCAT=SCAT,pch1=pch1,pch2=pch2,...)
#                outliers already removed if argument xout=T
list(output=mat)
}


#' Select Best Group - Multiple Group Comparisons
#'
#' For J independent groups, identifies the group with the largest trimmed mean
#' and makes an inference about whether it is significantly larger than all other groups.
#' Controls familywise error rate through adjusted critical p-values.
#'
#' @param x Matrix or list of data where each column (or list element) is a group.
#' @param p.crit Critical p-values for comparisons. If NULL, determined via simulation
#'   to control familywise error rate at level alpha.
#' @inheritParams common-params
#' @param iter Number of iterations for simulation to determine critical values (default: 5000).
#' @param NEG Logical. If TRUE, negates all values (to find smallest instead of largest).
#'
#' @return S4 object of class 'BIN' with slots:
#'   \item{Group.with.largest.estimate}{Index of group with largest estimate}
#'   \item{Larger.than}{Groups that the best group is significantly larger than}
#'   \item{n}{Sample sizes}
#'   \item{output}{Matrix with detailed comparison results}
#'
#' @note Confidence intervals DO NOT necessarily have simultaneous probability
#'   coverage 1-alpha using the adjusted level.
#'
#' @export
#' @examples
#' \dontrun{
#' x <- list(rnorm(30), rnorm(30, mean=0.5), rnorm(30, mean=1))
#' result <- anc.best(x, alpha=0.05)
#' }
anc.best<-function(x,p.crit=NULL,alpha=.05,tr=.2,iter=5000,SEED=TRUE,NEG=FALSE){
#
#
#  For J independent groups
#  Identify group with largest trimmed mean
#  Make a decision if every  p.value<=p.crit
#
#  p.crit:  If NULL, critical p-values are determined so that that FWE is alpha
#  This is done using a simulation to determine the null distribution based on
#  iter=5000 replications.
#
#
#  Returns:
#   Best='No Decision' if not significant
#   Best= the group with largest measure if a decision can be made.
#
#   Confidence intervals DO NOT necessarily  have simultaneous probability coverage 1-alpha
#  using the adjusted level.
#
x=elimna(x)
if(is.matrix(x))x=listm(x)
J=length(x)
if(J<3)stop('Should have 3 or more groups')
if(NEG)for(j in 1:J)x[[j]]=0-x[[j]]
Jm1=J-1
est=lapply(x,tmean,tr=tr)
n=lapply(x,length)
est=matl(est)
n=as.vector(matl(n))
R=order(est,decreasing = TRUE)
pvec=NA
if(is.null(p.crit)){
v=anc.best.crit(J=J,n=n,iter=iter,alpha=alpha,SEED=SEED)
p.crit=v$fin.crit
pvdist=v$pvdist
}
output<-matrix(NA,Jm1,8)
dimnames(output)=list(NULL,c('Est.Best','Grp','Est','Dif','ci.low','ci.up','p.value','p.crit'))
for(i in 2:J){
im1=i-1
a=yuen(x[[R[1]]],x[[R[[i]]]],alpha=p.crit[im1])
pvec[im1]=a$p.value
output[im1,]=c(a$est.1,R[[i]],a$est.2,a$dif,a$ci[1],a$ci[2],a$p.value,p.crit[im1])
}
Best='No Decisions'
flag=sum(output[,7]<=output[,8])
id=output[,7]<=output[,8]
if(sum(id>0))Best=output[id,2]
if(flag==Jm1)Best='All'
#setClass('BIN',slots=c('Group.with.largest.estimate','Select.Best.p.value','Larger.than','n','output'))
setClass('BIN',slots=c('Group.with.largest.estimate','Larger.than','n','output'))
#put=new('BIN',Group.with.largest.estimate=R[[1]],Select.Best.p.value=dpv,Larger.than=Best,n=n,output=output)
put=new('BIN',Group.with.largest.estimate=R[[1]],Larger.than=Best,n=n,output=output)
put
}


#' Determine Critical P-Values for anc.best
#'
#' Simulates the null distribution to determine critical p-values that control
#' the familywise error rate for the anc.best function.
#'
#' @param J Number of groups.
#' @param n Sample size (single value for equal n, or vector of length J for unequal n).
#' @inheritParams common-params
#' @param iter Number of simulation iterations (default: 5000).
#'
#' @return List with:
#'   \item{fin.crit}{Vector of critical p-values}
#'   \item{pvdist}{Matrix of simulated p-value distributions}
#'
#' @export
#' @examples
#' \dontrun{
#' anc.best.crit(J=4, n=30, alpha=0.05)
#' }
anc.best.crit<-function(J,n=30,alpha=.05,tr=.2,iter=5000,SEED=TRUE){
#
#  Determine critical p-values for anc.best
#
if(SEED)set.seed(2)
Jm1=J-1
rem=matrix(NA,iter,Jm1)
for(k in 1:iter){
if(length(n)==1){
x=rmul(n,p=J)
x=listm(x)
}
else{
x=list()
if(length(n)!=J)stop('J is not equal to the length  of n')
for(j in 1:J)x[[j]]=rnorm(n[j])
}
rem[k,]=anc.best.ex(x,tr=tr)
}
#
init=apply(rem,2,qest,alpha)
z=optim(0,anc.best.fun,init=init,iter=iter,rem=rem,Jm1=Jm1,alpha=alpha,method='Brent',lower=0,upper=1)
fin.crit=z$par*init
list(fin.crit=fin.crit,pvdist=rem)
}


#' Determine Critical P-Values for anc.best - Detailed Version
#'
#' Similar to anc.best.crit but computes critical p-values for a range of alpha levels
#' (from 0.001 to 0.99). Useful for sensitivity analysis and understanding the relationship
#' between nominal alpha and critical p-values.
#'
#' @param J Number of groups.
#' @param n Sample size (single value for equal n, or vector of length J for unequal n).
#' @inheritParams common-params
#' @param iter Number of simulation iterations (default: 5000).
#'
#' @return Matrix where each row corresponds to an alpha level and columns contain
#'   critical p-values for each of the J-1 comparisons.
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' # Compute critical values for range of alpha levels
#' crit_matrix <- anc.best.crit.det(J=4, n=30, iter=1000)
#' }
anc.best.crit.det<-function(J,n,alpha=.05,tr=.2,iter=5000,SEED=TRUE){
#
#  Determine critical p-values for anc.best
#
if(SEED)set.seed(2)
Jm1=J-1
rem=matrix(NA,iter,Jm1)
for(k in 1:iter){
if(length(n)==1){
x=rmul(n,p=J)
x=listm(x)
}
else{
x=list()
if(length(n)!=J)stop('J is not equal to the length  of n')
for(j in 1:J)x[[j]]=rnorm(n[j])
}
rem[k,]=anc.best.ex(x,tr=tr)
}
aval=c(seq(.001,.1,.001),seq(.011,.99,.01))
na=length(aval)
fin.crit=matrix(NA,na,Jm1)
for(i in 1:na){
init=apply(rem,2,qest,aval[i])
z=optim(0,anc.best.fun,init=init,iter=iter,rem=rem,Jm1=Jm1,alpha=aval[i],method='Brent',lower=0,upper=1)
fin.crit[i,]=z$par*init
}
fin.crit
}




#' Extract P-Values for Best Group Comparisons
#'
#' Helper function for anc.best.crit. Computes p-values comparing the group with
#' the largest trimmed mean to all other groups.
#'
#' @inheritParams common-params
#' @param x Matrix or list of data where each column (or list element) is a group.
#'
#' @return Vector of p-values from comparing best group to each other group.
#'
#' @keywords internal
anc.best.ex<-function(x,tr=.2){
#
#  Used by anc.best.crit
pvec=NA
x=elimna(x)
if(is.matrix(x))x=listm(x)
J=length(x)
est=lapply(x,tmean,tr=tr)
est=matl(est)
R=order(est,decreasing = TRUE)
ic=0
for(j in 2:J){
ic=ic+1
pvec[ic]=yuen(x[[R[1]]],x[[R[[j]]]],tr=tr)$p.value
}
pvec
}


#' Select Best Group with Adjusted P-Values
#'
#' Similar to anc.best but uses adjusted p-values (via p.adjust) instead of
#' critical p-values determined by simulation. Identifies the group with the largest
#' trimmed mean and provides inference with familywise error control.
#'
#' @param x Matrix or list of data where each column (or list element) is a group.
#' @param rem Matrix of simulated p-value distributions. If NULL, computed via simulation.
#' @inheritParams common-params
#' @param iter Number of iterations for simulation if rem is NULL (default: 5000).
#'
#' @return S4 object of class 'BIN' with slots:
#'   \item{Group.with.largest.estimate}{Index of group with largest estimate}
#'   \item{Larger.than}{Groups that the best group is significantly larger than}
#'   \item{n}{Sample sizes}
#'   \item{output}{Matrix with detailed comparison results including adjusted p-values}
#'
#' @note Confidence intervals have simultaneous probability coverage 1-alpha.
#'
#' @export
#' @examples
#' \dontrun{
#' x <- list(rnorm(30), rnorm(30, mean=0.5), rnorm(30, mean=1))
#' result <- anc.bestH(x, alpha=0.05)
#' }
anc.bestH<-function(x,rem=NULL,alpha=.05,tr=.2,iter=5000,SEED=TRUE){
#
#
#  Identify group with largest trimmed mean
#  Make a decision if every  p.value<=p.crit
#
#  p.crit:  If NULL, critical p-values are determined so that that FWE is alpha
#  This is done using a simulation to determine the null distribution based on
#  iter=5000 replications.
#
#
#  Returns:
#   Best='No Decision' if not significant
#   Best= the group with largest measure if a decision can be made.
#
#   Confidence intervals having simultaneous probability coverage 1-alpha
#  using the adjusted level.
#
x=elimna(x)
if(is.matrix(x))x=listm(x)
J=length(x)
if(J<3)stop('Should have 3 or more groups')
Jm1=J-1
est=lapply(x,tmean,tr=tr)
n=lapply(x,length)
est=matl(est)
n=as.vector(matl(n))
R=order(est,decreasing = TRUE)
pvec=NA
if(is.null(rem)){
pvdist=anc.bestH.crit(J=J,n=n,iter=iter,alpha=alpha,SEED=SEED)
}
output<-matrix(NA,Jm1,8)
dimnames(output)=list(NULL,c('Est.Best','Grp','Est','Dif','ci.low','ci.up','p.value','p.adj'))
for(i in 2:J){
im1=i-1
a=yuen(x[[R[1]]],x[[R[[i]]]],alpha=qest(pvdist[,im1],alpha))
pvec[im1]=mean(pvdist<=a$p.value)
output[im1,1:7]=c(a$est.1,R[[i]],a$est.2,a$dif,a$ci[1],a$ci[2],a$p.value)
}
output[,8]=p.adjust(output[,7])
Best='No Decisions'
id=output[,8]<=alpha
if(sum(id>0))Best=output[id,2]
if(sum(id)==Jm1)Best='All'
setClass('BIN',slots=c('Group.with.largest.estimate','Larger.than','n','output'))
put=new('BIN',Group.with.largest.estimate=R[[1]],Larger.than=Best,n=n,output=output)
put
}


#' Select Best Group - Bootstrap Percentile Method
#'
#' Bootstrap percentile version of anc.best. Identifies the group with the largest
#' location estimate (default: trimmed mean) and makes inference using bootstrap
#' percentile confidence intervals with familywise error rate control.
#'
#' @param x Matrix or list of data where each column (or list element) is a group.
#' @param loc.fun Location estimator function (default: tmean for trimmed mean).
#' @param p.crit Critical p-values for comparisons. If NULL, determined via simulation.
#' @inheritParams common-params
#' @param iter Number of iterations for determining critical values (default: 5000).
#' @param ... Additional arguments passed to loc.fun.
#'
#' @return List with:
#'   \item{Group.with.largest.est}{Index of group with largest estimate}
#'   \item{Best}{Index of best group if significant decision can be made, 'No Decision' otherwise}
#'   \item{n}{Sample sizes}
#'   \item{output}{Matrix with detailed comparison results}
#'
#' @export
#' @examples
#' \dontrun{
#' x <- list(rnorm(30), rnorm(30, mean=0.5), rnorm(30, mean=1))
#' result <- anc.bestpb(x, nboot=2000, alpha=0.05)
#' }
anc.bestpb<-function(x,loc.fun=tmean,nboot=3000,p.crit=NULL,alpha=.05,iter=5000,SEED=TRUE,...){
#
#
#  Identify group with largest trimmed mean
#  Make a decision if every  p.value<=p.crit
#
#  p.crit:  If NULL, critical p-values are determined so that that FWE is alpha
#  This is done using a simulation to determine the null distribution based on
#  iter=5000 replications.
#
#
#  Returns:
#   Best='No Decision' if not significant
#   Best= the group with largest measure if a decision can be made.
#
#   Confidence intervals having simultaneous probability coverage 1-alpha
#  using the adjusted level.
#
x=elimna(x)
if(is.matrix(x))x=listm(x)
J=length(x)
if(J<3)stop('Should have 3 or more groups')
Jm1=J-1
est=lapply(x,loc.fun,...)
n=lapply(x,length)
est=matl(est)
n=as.vector(matl(n))
R=order(est,decreasing = TRUE)
pvec=NA
if(is.null(p.crit))p.crit=anc.best.crit(J=J,n=n,iter=iter,alpha=alpha,SEED=SEED)$fin.crit
output<-matrix(NA,Jm1,8)
dimnames(output)=list(NULL,c('Est.Best','Grp','Est','Dif','ci.low','ci.up','p.value','p.crit'))
for(i in 2:J){
im1=i-1
a=pb2gen(x[[R[1]]],x[[R[[i]]]],alpha=p.crit[im1],nboot=nboot,est=loc.fun,SEED=SEED,...)
pvec[im1]=a$p.value
output[im1,]=c(a$est.1,R[[i]],a$est.2,a$est.dif,a$ci[1],a$ci[2],a$p.value,p.crit[im1])
}
Best='No Decision'
flag=sum(pvec<=p.crit)
if(flag==Jm1)Best=R[[1]]
list(Group.with.largest.est=R[[1]],Best=Best,n=n,output=output)
}

#' Select Best Group with Overall P-Value - Bootstrap Version
#'
#' Similar to anc.bestpb but also computes an overall p-value for the decision about
#' which group has the largest measure of location. Uses bootstrap percentile method.
#'
#' @param x Matrix or list of data where each column (or list element) is a group.
#' @param loc.fun Location estimator function (default: tmean for trimmed mean).
#' @inheritParams common-params
#' @param iter Number of iterations for determining critical values (default: 5000).
#' @param ... Additional arguments passed to loc.fun.
#'
#' @return List with:
#'   \item{Group.with.largest.est}{Index of group with largest estimate}
#'   \item{Best}{Index of best group if significant, 'No Decision' otherwise}
#'   \item{overall.pvalue}{Overall p-value for the selection decision}
#'   \item{n}{Sample sizes}
#'   \item{output}{Matrix with detailed comparison results}
#'
#' @note The overall p-value is computed by comparing the observed decision pattern
#'   to its null distribution over a range of alpha levels.
#'
#' @export
#' @examples
#' \dontrun{
#' x <- list(rnorm(30), rnorm(30, mean=0.5), rnorm(30, mean=1))
#' result <- anc.bestpb.PV(x, nboot=2000)
#' result$overall.pvalue
#' }
 anc.bestpb.PV<-function(x,loc.fun=tmean,nboot=2000,alpha=.05,iter=5000,SEED=TRUE,...){
#
#
#  Identify group with largest trimmed mean
#  Make a decision if every  p.value<=p.crit
#
#  p.crit:  If NULL, critical p-values are determined so that that FWE is alpha
#  This is done using a simulation to determine the null distribution based on
#  iter=5000 replications.
#
#
#  Returns:
#  a p-value related to making a decision about which group has the largest measure of location.
#   Best='No Decision' if not significant
#   Best= the group with largest measure if a decision can be made.
#
#   Confidence intervals having simultaneous probability coverage 1-alpha
#  using the adjusted level.
#
x=elimna(x)
if(is.matrix(x))x=listm(x)
J=length(x)
if(J<3)stop('Should have 3 or more groups')
Jm1=J-1
est=lapply(x,loc.fun,...)
n=lapply(x,length)
est=matl(est)
n=as.vector(matl(n))
R=order(est,decreasing = TRUE)
pvec=NA
aval=c(seq(.001,.1,.001),seq(.11,.99,.01))
id=which(aval==alpha)
if(length(id)==0)stop('alpha be one of the values .001(.001).1 or 11(.01).99')
v=anc.best.crit.det(J=J,n=n,iter=iter,alpha=alpha,SEED=SEED)
p.crit=v[id,]
if(is.null(p.crit))p.crit=anc.best.crit(J=J,n=n,iter=iter,alpha=alpha,SEED=SEED)
output<-matrix(NA,Jm1,8)
dimnames(output)=list(NULL,c('Est.Best','Grp','Est','Dif','ci.low','ci.up','p.value','p.crit'))
for(i in 2:J){
im1=i-1
a=pb2gen(x[[R[1]]],x[[R[[i]]]],alpha=p.crit[im1],nboot=nboot,est=loc.fun,SEED=SEED,...)
pvec[im1]=a$p.value
output[im1,]=c(a$est.1,R[[i]],a$est.2,a$est.dif,a$ci[1],a$ci[2],a$p.value,p.crit[im1])
}
# Determine p-value for overall decision
na=length(aval)
for(i in 1:na){
chk=sum(output[,7]<=v[i,])
pv=aval[i]
if(chk==Jm1)break
}
Best='No Decisions'
flag=sum(output[,7]<=output[,8])
id=output[,7]<=output[,8]
if(sum(id>0))Best=output[id,2]
if(flag==Jm1)Best='All'
setClass('BIN',slots=c('Group.with.largest.estimate','Select.Best.p.value','Larger.than','n','output'))
put=new('BIN',Group.with.largest.estimate=R[[1]],Select.Best.p.value=pv,Larger.than=Best,n=n,output=output)
put
}



#' Compute Crossover Region for Johnson-Neyman Technique
#'
#' Internal helper function for Johnson-Neyman analysis. Determines the crossover
#' region (range of covariate values) where group differences in regression lines
#' are not statistically significant.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#'
#' @return Vector of length 2 containing the lower and upper bounds of the crossover region.
#'
#' @keywords internal
ancCR<-function(x1,y1,x2,y2){
v=optim(0,JNH_sub1,x1=x1,y1=y1,x2=x2,y2=y2,method='BFGS')$par
v[2]=optim(0,JNH_sub2,x1=x1,y1=y1,x2=x2,y2=y2,method='BFGS')$par
a=min(v)
v=c(a,max(v))
}



#' ANCOVA Effect Sizes for Dependent Groups with Multiple Covariates
#'
#' For dependent (paired/repeated measures) groups with two or more covariates,
#' estimates effect sizes at specified design points. Uses a running interval
#' smoother approach to identify nearest neighbors and compute effect sizes
#' based on their response values.
#'
#' @param x1 Matrix of covariate values (must have 2+ columns).
#' @param y1 Response values for condition 1.
#' @param y2 Response values for condition 2 (paired with y1).
#' @param fr1 Span for condition 1 running interval smoother (default: 1.5).
#' @param fr2 Span for condition 2 running interval smoother (default: 1.5).
#' @inheritParams common-params
#' @param pts Matrix of design points where effect sizes should be computed.
#'   If NULL, uses significant points from ancdetM4.
#' @param xout Logical. If TRUE, remove outliers from covariates before analysis.
#' @param cov.fun Covariance function for determining nearest neighbors (default: skip.cov).
#' @param ... Additional arguments passed to outfun or cov.fun.
#'
#' @return List with:
#'   \item{pts}{Matrix of design points analyzed}
#'   \item{ES}{Matrix of effect sizes (EP, QS, QStr, AKP) at each point}
#'   \item{SML.class}{Classification of effect sizes (small/medium/large)}
#'
#' @note Minimum sample size for each group is 50.
#'
#' @export
#' @examples
#' \dontrun{
#' # Paired data with 2 covariates
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- rnorm(50)
#' y2 <- y1 + rnorm(50, mean=0.5)
#' result <- ancDEP.MULC.ES(x1, y1, y2, fr1=1.2, fr2=1.2)
#' }
ancDEP.MULC.ES<-function(x1,y1,y2,fr1=1.5,fr2=1.5,tr=.2,pts=NULL,xout=FALSE,outfun=outpro,cov.fun=skip.cov,...){
#
#
#  Dependent groups
#  For two or more covariates, estimate effect sizes for
#  a collection of points.
#
#  That is, for each point of interest, determine
#  a cloud of points close to it and based on the
#  corresponding y values, compute measures of effect size
#
#  If pts=NULL
#  the significant points returned by
#  ancdetM4 are used
#
x2=x1
p=ncol(x1)
if(p<2)stop('This function is for two or more covariates')
p1=p+1
#if(ncol(x2)!=p)stop('x1 and x2 do not have the same number of columns')
xy=elimna(cbind(x1,y1))
x1=xy[,1:p]
y1=xy[,p1]
xy=elimna(cbind(x2,y2))
x2=xy[,1:p]
y2=xy[,p1]
if(min(length(y1),length(y2))<50)stop('The minimum sample size must be greater than or equal to  50')
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag,]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag,]
y2<-y2[flag]
}
guide1=NA
guide2=NA
x1=as.matrix(x1)
x2=as.matrix(x2)
a=Dancovamp(x1=x1,y1=y1,x2=x2,y2=y2,fr1=fr1,fr2=fr2,tr=tr,pts=pts,...)
if(is.null(pts))pts=a$pts
M=NA
SML.class=NA
if(!is.null(pts)){
if(is.vector(pts))pts=matrix(pts,nrow=1)
nr=nrow(pts)
M=matrix(NA,nr,4)
SML.class=matrix(0,nr,4)
m1<-cov.fun(x1)
m2<-cov.fun(x2)
for(i in 1:nr){
id1=near3d(x1,pts[i,],fr1,m1)
id2=near3d(x2,pts[i,],fr2,m2)
if(sum(id1)<10 ||sum(id2)<10)print(paste('For point',j,'not enough nearest neighbors'))
if(sum(id1)>=10 & sum(id2)>10){
ES=dep.ES.summary(y1[id1],y2[id2])
M[i,]=ES[,2]
}}
dimnames(M)=list(NULL,names(ES[,1]))
dum1=rnorm(50)
dum2=rnorm(50)
guide1=dep.ES.summary(dum1,dum2+3)[,-2]
guide2=dep.ES.summary(dum1,dum2-3)[,-2]
}
if(!is.na(M[1])){
for(i in 1:nr){
flag=M[i,] <= guide1[,2] || M[i,]>= guide2[,2]
SML.class[i,flag]=1
flag=M[i,] <= guide1[,3] || M[i,]>= guide2[,3]
SML.class[i,flag]=2
flag=M[i,] <= guide1[,4] || M[i,]>= guide2[,4]
SML.class[i,flag]=3
}}
leg1='0=At most small, 1=between S and M, 2=between M and L, 3=greater than L'
list(ES.REL.MAG.G1.less.than.G2=guide1, ES.REL.MAG.G1.greater.than.G2=guide2,Est=M,legend.4.SML.class=leg1,SML.class=SML.class)
}


#' Choose Design Points for ANCOVA Analysis
#'
#' Selects design points for ANCOVA analysis based on data depth. Can either return
#' the deepest point plus points on the 0.5 depth contour, or return a specified
#' fraction of the deepest points.
#'
#' @param x Matrix or data frame of covariate values (n by p matrix).
#' @param depfun Depth function to use (default: fdepth for functional data depth).
#' @param DH Logical. If FALSE, returns deepest point and 0.5 depth contour points.
#'   If TRUE, returns FRAC deepest points.
#' @param FRAC Fraction between 0 and 1. When DH=TRUE, returns the 1-FRAC deepest points.
#'   For example, FRAC=0.2 returns the deepest 80% of points.
#' @param ... Additional arguments passed to depfun.
#'
#' @return Matrix of selected design points with same number of columns as x.
#'
#' @note Useful for selecting covariate values where group comparisons should be made.
#'
#' @export
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(100), ncol=2)
#' # Get deepest point and 0.5 contour
#' pts1 <- ancdes(x, DH=FALSE)
#' # Get deepest 80% of points
#' pts2 <- ancdes(x, DH=TRUE, FRAC=0.2)
#' }
ancdes<-function(x,depfun=fdepth,DH=FALSE,FRAC=.5,...){
#
#  Choose points for design of an ANCOVA
#  x is the n by p matrix m.
#
#  FRAC some value between 0 and 1.
#
#    FRAC is  the fraction of the least deep points that will not be returned when
#   DH=TRUE
#  That is, return 1-FRAC deepest points.
#   For example, FRAC=.2 means that  the deepest 80% of the
#   data will be returned.
#
#   DH=F, return deepest point and those points on the
#   .5 depth contour
#
if(is.data.frame(x))x=as.matrix(x)
if(!is.matrix(x))stop("x must be a matrix or a data frame")
temp<-depfun(x,plotit=FALSE,...)
temp2<-order(temp)
if(!DH){
val<-matrix(x[temp2[length(temp)],],ncol=ncol(x))
nmid<-round(length(temp)/2)
id2<-(temp[temp2[nmid]]==temp)
val2<-matrix(x[id2,],ncol=ncol(x))
if(!is.matrix(val2))val2<-t(as.matrix(val2))
val<-rbind(val,val2)
}
if(DH){
bot=round(length(temp)*FRAC)
val=matrix(x[temp2[bot:length(temp)],],ncol=ncol(x))
}
val=elimna(val)
val
}



#' Detailed ANCOVA Analysis with Multiple Design Points
#'
#' Enhanced version of ancova that provides detailed analysis at multiple design points
#' along the covariate range. Computes group differences and confidence bands with
#' familywise error rate control.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#' @param fr1 Span for group 1 running interval smoother (default: 1).
#' @param fr2 Span for group 2 running interval smoother (default: 1).
#' @inheritParams common-params
#' @param method Effect size measure: "EP" (explanatory power), "QS" (quantile shift),
#'   "AKP" (trimmed mean analog of Cohen's d), or "WMW" (P(X<Y)).
#' @param plot.dif Logical. If TRUE, plots the difference with confidence band.
#' @param pts Design points for analysis. If NA, determined automatically.
#' @param sm Logical. Use smoothing for estimates.
#' @param pr Logical. Print progress information.
#' @param xout Logical. Remove outliers from covariates before analysis.
#' @param npts Number of covariate values to use (default: 25).
#' @param p.crit Critical p-value for familywise error control. If NULL, determined by simulation.
#' @param nreps Number of simulation replications for determining p.crit (default: 5000).
#' @param EST Logical. If TRUE, estimate critical value via simulation rather than using stored values.
#' @param SCAT Logical. Show scatterplot of data.
#' @param xlab,ylab Axis labels for plots.
#' @param pch1,pch2 Plot characters for groups 1 and 2.
#' @param ... Additional arguments.
#'
#' @return List containing analysis results with adjusted confidence intervals
#'   having simultaneous probability coverage 1-alpha.
#'
#' @note One covariate only is allowed with this function.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- 0.5*x2 + rnorm(50)
#' result <- ancdet(x1, y1, x2, y2, method="EP", npts=20)
#' }
ancdet<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,
alpha=.05,method="EP",plotit=TRUE,plot.dif=FALSE,pts=NA,sm=FALSE,
pr=TRUE,xout=FALSE,outfun=out,MC=FALSE,
npts=25,p.crit=NULL,nreps=5000,SEED=TRUE,EST=FALSE,
SCAT=TRUE,xlab='X',ylab='Y',pch1='*',pch2='+',...){
#
#  Like the function ancova, but a more detailed analysis
#  plot.dif=TRUE:  plot difference in the estimates plus a
# confidence band having simultaneous probability coverate 1-alpha
#
# npts = number of  covariate values to be used
#
#  Argument method indicates which measure of effect size will be used
#  EP: explanatory measure of effect size
#  QS: quantile shift measure of effect size
#  AKP:  trimmed mean Winsorized variance analog of Cohen's d
#  WMW:  P(X<Y)
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(is.null(p.crit))set.seed(2)
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
if(xout){
flag<-outfun(x1)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2)$keep
x2<-x2[flag]
y2<-y2[flag]
}
res1=ancova(x1,y1,x2,y2,pr=FALSE,plotit=FALSE,method=method)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)

if(alpha!=.05)EST=TRUE
if(is.null(p.crit)){
if(!EST){
nv=c(30,  50,  60,  70,  80, 100,
150, 200, 300, 400, 500, 600, 800)
pv=c(0.00824497,0.00581, 0.005435089, 0.004763079,
0.00416832, 0.004406774, 0.00388228,0.003812836,0.003812836,0.003453055, 0.003625061,
.003372966, 0.003350022)
n1= length(y1)
 n2=length(y2)
p.crit=(lplot.pred(1/nv,pv,1/n1)$yhat+lplot.pred(1/nv,pv,1/n2)$yhat)/2
}
if(EST)p.crit=ancdet.pv(length(y1),length(y2),nreps=nreps,tr=tr,npts=npts,MC=MC)
}

if(plot.dif)plotit=FALSE
critv=qnorm(1-p.crit/2)
res=ancova(x1,y1,x2,y2,fr1=fr1,fr2=fr2,tr=tr,alpha=alpha,pr=FALSE,plotit=plotit,pts=pts,SCAT=SCAT)$output
res[,7]=res[,4]-critv*res[,6] # adjust confidence interval based on adjusted p-value
res[,8]=res[,4]+critv*res[,6] # adjust confidence interval based on adjusted p-value
if(plot.dif){
yhat=plot(c(res[,1],res[,1],res[,1]),c(res[,4],res[,7],res[,8]),type='n',xlab=xlab,ylab=ylab)
z1=lplot(res[,1],res[,4],plotit=FALSE,pyhat=TRUE)$yhat
z2=lplot(res[,1],res[,7],plotit=FALSE,pyhat=TRUE)$yhat
z3=lplot(res[,1],res[,8],plotit=FALSE,pyhat=TRUE)$yhat
lines(res[,1],z1)
lines(res[,1],z2,lty=2)
lines(res[,1],z3,lty=2)
}
sig=rep(0,nrow(res))
sig[res[,9]<=p.crit]=1
sig=as.matrix(sig,ncol=1)
dimnames(sig)=list(NULL,'Sig.Dif')
res=cbind(res,sig)
list(p.crit=p.crit,output=res[,-10],num.sig=sum(sig),p.crit=p.crit)
}



#' ANCOVA Detection Method for Two Covariates
#'
#' Performs two-group ANCOVA comparisons using two covariates with familywise
#' error rate control. Uses method MC3 from Wilcox (2017) with an improvement
#' on Hochberg's method. Compares groups at the deepest half of covariate values.
#'
#' @param x1 Matrix of covariate values for group 1 (n1 x 2).
#' @param y1 Dependent variable for group 1.
#' @param x2 Matrix of covariate values for group 2 (n2 x 2).
#' @param y2 Dependent variable for group 2.
#' @param fr1 Span for running interval smoother for group 1 (default: 1).
#' @param fr2 Span for running interval smoother for group 2 (default: 1).
#' @inheritParams common-params
#' @param test Test function to use: yuen (default), qcomhd, or qcomhdMC.
#' @param q Quantile for qcomhd test (default: 0.5).
#' @param op Logical. If TRUE, plots estimated differences for covariate points.
#' @param pts Design points (currently not fully implemented).
#' @param sm Logical. Smoothing parameter.
#' @param FRAC Fraction for depth-based methods (default: 0.5 for deepest half).
#' @param p.crit Critical p-value for significance. If NULL, computed via simulation.
#' @param nreps Number of replications for p.crit computation (default: 2000).
#' @param FAST Logical. If TRUE, uses pre-computed critical values when possible.
#' @param ticktype Tick type for lplot ('simple' or 'detail').
#' @param xlab Label for first covariate axis.
#' @param ylab Label for second covariate axis.
#' @param zlab Label for z-axis when op=TRUE.
#' @param pch1 Plot character for non-significant points (default: '*').
#' @param pch2 Plot character for significant points (default: '+').
#'
#' @return List with components:
#'   \item{num.sig}{Number of significant points found}
#'   \item{p.crit}{Critical p-value used}
#'   \item{points.used}{Matrix of covariate points and test results}
#'   \item{sig.points}{Matrix of significant covariate points}
#'
#' @references
#' Wilcox, R.R. (2017). Introduction to Robust Estimation and Hypothesis Testing (4th ed.).
#'
#' @note Minimum sample size of 50 per group required.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- x1[,1] + x1[,2] + rnorm(50)
#' x2 <- matrix(rnorm(100), ncol=2)
#' y2 <- x2[,1] + x2[,2] + 0.5 + rnorm(50)
#' ancdet2C(x1, y1, x2, y2)
#' }
ancdet2C<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,test=yuen,q=.5,
alpha=.05,plotit=TRUE,op=FALSE,pts=NA,sm=FALSE,FRAC=.5,
pr=TRUE,xout=FALSE,outfun=outpro,MC=FALSE,
p.crit=NULL,nreps=2000,SEED=TRUE,FAST=TRUE,ticktype='detail',
xlab='X1',ylab='X2',zlab='Y',pch1='*',pch2='+',...){
#
#  method MC3 in Wilcox (2017, Intro to Robust Estimation and Hypothesis Testing, 4th ed.)
#  Multiple comparisons using an improvement on Hochberg to control FWE
#
#  Like ancdet, only two covariate values can be used.
#  Like method MC2,
#  use the deepest half of the covariate values.
#
#  politit=TRUE. Plot covariate points. Significant points are indicated by
#  pch='+'
#
#  test can have one of three values: yuen (default), qcomhd or qcomhdMC
#
if(ncol(as.matrix(x1))!=2)stop('Two covariates only can be used')
if(is.null(p.crit))set.seed(2)
xy=elimna(cbind(x1,y1))
x1=xy[,1:2]
y1=xy[,3]
xy=elimna(cbind(x2,y2))
x2=xy[,1:2]
y2=xy[,3]
if(min(length(y1),length(y2))<50)stop('The minimum sample size must be greater than or equal to  50')
if(xout){
flag<-outfun(x1,plotit=FALSE)$keep
x1<-x1[flag,]
y1<-y1[flag]
flag<-outfun(x2,plotit=FALSE)$keep
x2<-x2[flag,]
y2<-y2[flag]
}
if(FAST){
if(FRAC==.5){
if(is.null(p.crit)){
if(alpha==.05){
nv=c(50,  55,  60,  70,  80, 100, 200, 300, 400, 500, 600,800)
pv=c(0.004585405, 0.003199894, 0.002820089, 0.002594342, 0.002481210, 0.001861313,
  0.001419821, 0.001423000, 0.001313700, 0.001351900, 0.001075, 0.00095859)
 n1= length(y1)
 n2=length(y2)
p.crit=(lplot.pred(1/nv,pv,1/n1)$yhat+lplot.pred(1/nv,pv,1/n2)$yhat)/2
  # Using K=n1 points, i.e., K=n1 tests are performed
if(max(n1,n2)>max(nv)){
p.crit=min(pv)
print('Warning: p.crit has not been computed exactly for sample sizes greater than 800')
if(n1>800)p.crit1=regYhat(1/pv[8:12,1],pv[8:12,2],1/n1)
if(n1<=800)p.crit1=lplot.pred(1/nv,pv,1/n1)$yhat
if(n2>800)p.crit1=regYhat(1/pv[8:12,1],pv[8:12,2],1/n2)
if(n2<=800)p.crit1=lplot.pred(1/nv,pv,1/n2)$yhat
p.crit=(p.crit1+p.crit2)/2
}
}}}}
res1=ancov2COV(x1,y1,x2,y2,DETAILS=TRUE,pr=FALSE,FRAC=FRAC,tr=tr,test=test,q=q,MC=MC)
if(is.null(p.crit))p.crit=ancdet2C.pv(length(y1),length(y2),MC=MC,nreps=nreps,
SEED=SEED)
LL=length(ncol(res1$all.results))
if(LL==1)num.sig=sum(res1$all.results[,3]<=p.crit)
if(LL==0)num.sig=NA
sig.points=NA
if(LL==1){
flag=res1$all.results[,3]<=p.crit
sig.points=res1$all.points.used[flag,1:2]
}
if(plotit){
if(!op){
if(pr)print('To plot the estimated differences for the covariate points used, set op=TRUE')
if(LL==0)plot(res1$all.points.used[,1],res1$all.points.used[,2],xlab=xlab,ylab=ylab)
if(LL==1){
plot(res1$all.points.used[,1],res1$all.points.used[,2],type='n',xlab=xlab,ylab=ylab)
points(res1$all.points.used[!flag,1],res1$all.points.used[!flag,2],pch=pch1)
points(res1$all.points.used[flag,1],res1$all.points.used[flag,2],pch=pch2)
}}
if(op)
lplot(res1$all.points.used[,1:2],res1$all.points.used[,3],xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype)
}
list(num.sig=num.sig,p.crit=p.crit,points.used=cbind(res1$all.points.used[,1:3],res1$all.results),sig.points=sig.points)
}



#' Critical P-Value for ancdet2C
#'
#' Determines the critical p-value for ancdet2C via simulation under the null
#' hypothesis. Uses the Harrell-Davis estimator to estimate the appropriate
#' quantile of the null distribution.
#'
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param nreps Number of simulation replications (default: 2000).
#' @inheritParams common-params
#' @param FRAC Fraction for depth-based methods (default: 0.5).
#'
#' @return Critical p-value (numeric).
#'
#' @keywords internal
ancdet2C.pv<-function(n1,n2,nreps=2000,alpha=.05,FRAC=.5,tr=.2,MC=FALSE,SEED=TRUE){
pvals=NA
xy=list()
n=max(c(n1,n2))
nmiss=n-min(c(n1,n2))
for (i in 1:nreps){
xy[[i]]=rmul(n,p=6)
xy[[i]][1:nmiss,1:3]=NA
}
if(!MC)pvals=lapply(xy,ancdet2C.sub,tr=tr,FRAC=FRAC)
if(MC){
pvals=mclapply(xy,ancdet2C.sub,tr=tr,FRAC=FRAC)
}
pvals=matl(pvals)
pv=hd(pvals,alpha)
pv
}


#' Internal Subroutine for ancdet2C.pv
#'
#' Helper function for ancdet2C.pv simulation. Computes test statistic
#' from simulated data.
#'
#' @param xy Simulated data matrix (6 columns: 2 covariates + y for each group).
#' @inheritParams common-params
#' @param FRAC Fraction for depth-based methods.
#'
#' @return Minimum p-value from ancov2COV test.
#'
#' @keywords internal
ancdet2C.sub<-function(xy,tr=.2,FRAC=.5){
#
xy1=elimna(xy[,1:3])
xy2=elimna(xy[,4:6])
x1=xy1[,1:2]
y1=xy1[,3]
x2=xy2[,1:2]
y2=xy2[,3]
res1=ancov2COV(x1,y1,x2,y2,pr=FALSE,FRAC=FRAC)$min.p.value
res1
}


#' ANCOVA Detection Method for Multiple Covariates
#'
#' Extension of ancdet2C that handles more than two covariates. Performs
#' two-group ANCOVA comparisons with familywise error rate control for
#' multidimensional covariate spaces.
#'
#' @param x1 Matrix of covariate values for group 1 (n1 x p).
#' @param y1 Dependent variable for group 1.
#' @param x2 Matrix of covariate values for group 2 (n2 x p).
#' @param y2 Dependent variable for group 2.
#' @param fr1 Span for running interval smoother for group 1 (default: 1).
#' @param fr2 Span for running interval smoother for group 2 (default: 1).
#' @inheritParams common-params
#' @param pts Design points matrix (if NA, uses unique x1 values).
#' @param BOTH Logical. If TRUE, uses unique values from both groups as design points.
#'
#' @return List with components:
#'   \item{num.sig}{Number of significant points}
#'   \item{p.crit}{Critical value for adjusted p-values}
#'   \item{output}{Matrix of test results}
#'   \item{sel.pts}{Selected design points used}
#'   \item{sig.points}{Significant design points}
#'
#' @note
#' - Minimum sample size of 50 per group required
#' - Requires at least 12 nearest neighbors at each design point
#' - Uses Hochberg adjustment for n.sel <= 25, FWE control otherwise
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- matrix(rnorm(150), ncol=3)
#' y1 <- rowSums(x1) + rnorm(50)
#' x2 <- matrix(rnorm(150), ncol=3)
#' y2 <- rowSums(x2) + 0.5 + rnorm(50)
#' ancdetM4(x1, y1, x2, y2)
#' }
ancdetM4<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,
alpha=.05, pts=NA,
pr=TRUE,xout=FALSE,outfun=outpro,MC=FALSE,BOTH=FALSE,...){
#
#  Like ancdet2C, only  more than two covariate values can be used.
#
p=ncol(x1)
p1=p+1
if(ncol(x2)!=p)stop('x1 and x2 do not have the same number of columns')
xy=elimna(cbind(x1,y1))
x1=xy[,1:p]
y1=xy[,p1]
xy=elimna(cbind(x2,y2))
x2=xy[,1:p]
y2=xy[,p1]
if(min(length(y1),length(y2))<50)stop('The minimum sample size must be greater than or equal to  50')
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag,]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag,]
y2<-y2[flag]
}
x1=as.matrix(x1)
x2=as.matrix(x2)
m1<-covmcd(x1)
m2<-covmcd(x2)
flag=NA
if(!BOTH)pts=unique(x1)
if(BOTH)pts=unique(rbind(x1,x2))
output=matrix(NA,nrow=nrow(pts),ncol=8)
dimnames(output)=list(NULL,c('n1','n2','est.1','est.2','ci.low','ci.up','p.value','dif'))
ic=0
set.pts=NULL
for(i in 1:nrow(pts)){
nval1<-length(y1[near3d(x1,pts[i,],fr1,m1)])
nval2<-length(y2[near3d(x2,pts[i,],fr2,m2)])
flag.chk=as.logical((nval1>11)*(nval2>11))
flag[i]=flag.chk
if(flag.chk){
Y1=y1[near3d(x1,pts[i,],fr=fr1,m=m1)]
Y2=y2[near3d(x2,pts[i,],fr=fr2,m=m2)]
temp=yuen(Y1,Y2,tr=tr)
temp=pool.a.list(temp[1:7])
ic=ic+1
output[ic,]=temp
}}
sel.pts=NULL
sig.points=NULL
if(ic>0){
n.sel=sum(flag)
output=output[1:n.sel,]
sel.pts=pts[flag,]
}

if(sum(flag)==0){
print('Could not find any point with 12 or more nearest neighbors')
output=matrix(NA,nrow=1,ncol=8)
n.sel=0
num.sig=0
}
id.sig=NULL
padj=NULL
p.crit=NULL
if(n.sel>0){
sel.id=c(1:n.sel)
if(n.sel<=25){
if(n.sel==1)padj=output[7]
else
padj=p.adjust(output[,7],method='hoch')
flag=padj<=alpha
if(sum(flag)==1)sig.points=sel.pts[flag]
if(sum(flag)>1)sig.points=sel.pts[flag,]
num.sig=sum(flag)
id.sig=sel.id[flag]
}
if(n.sel>25){
if(n.sel<=100)p.crit=0.0806452604/n.sel-0.0002461736
if(n.sel>100)p.crit=6.586286e-02/n.sel+4.137143e-05
flag=output[,7]<=p.crit
if(sum(flag)>0)sig.points=sel.pts[flag,]
}
num.sig=sum(flag)
id.sig=sel.id[flag]
}
list(selected.points=sel.pts,output=output,significant.points=sig.points,num.sig=num.sig,id.sig=id.sig)
}




#' Critical P-Value for ancdet
#'
#' Determines the critical p-value for ancdet via simulation under the null
#' hypothesis. Uses the Harrell-Davis estimator to find the appropriate
#' quantile of the null distribution.
#'
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param nreps Number of simulation replications (default: 2000).
#' @inheritParams common-params
#' @param npts Number of design points to use (default: 25).
#'
#' @return Critical p-value (numeric).
#'
#' @keywords internal
ancdet.pv<-function(n1,n2,nreps=2000,alpha=.05,npts=25,tr=.2,MC=FALSE,SEED=TRUE){
if(SEED)set.seed(2)
pvals=NA
xy=list()
n=max(c(n1,n2))
nmiss=n-min(c(n1,n2))
for (i in 1:nreps){
xy[[i]]=rmul(n,p=4)
xy[[i]][1:nmiss,1:2]=NA
}
if(!MC)pvals=lapply(xy,ancdet.sub,npts=npts,tr=tr)
if(MC){
pvals=mclapply(xy,ancdet.sub,npts=npts,tr=tr)
}
pvals=matl(pvals)
pv=hd(pvals,alpha)
pv
}


#' Internal Subroutine for ancdet.pv
#'
#' Helper function for ancdet.pv simulation. Performs ANCOVA at multiple
#' design points and returns the minimum p-value.
#'
#' @param xy Simulated data matrix (4 columns: x,y for each group).
#' @inheritParams common-params
#' @param plot.dif Logical. If TRUE, plots differences.
#' @param pts Design points (if NA, computed automatically).
#' @param sm Logical. Smoothing parameter.
#' @param LP Logical. Line plot parameter.
#' @param npts Number of design points (default: 25).
#' @param p.crit Critical p-value.
#' @param nreps Number of replications.
#' @param SCAT Logical. Scatter plot parameter.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param pch1 Plot character for group 1.
#' @param pch2 Plot character for group 2.
#'
#' @return Minimum p-value across all design points.
#'
#' @keywords internal
ancdet.sub<-function(xy,tr=.2,
alpha=.05,plotit=FALSE,plot.dif=FALSE,pts=NA,sm=FALSE,
pr=TRUE,xout=FALSE,outfun=out,LP=TRUE,
npts=25,p.crit=NULL,nreps=2000,
SCAT=TRUE,xlab='X',ylab='Y',pch1='*',pch2='+',...){
#
#  Like ancova, but a more detailed analysis based on using
# npts covariate values
#
xy1=elimna(xy[,1:2])
xy2=elimna(xy[,3:4])
x1=xy1[,1]
y1=xy1[,2]
x2=xy2[,1]
y2=xy2[,2]
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
res1=ancova(x1,y1,x2,y2,pr=FALSE,plotit=FALSE)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
res=ancova(x1,y1,x2,y2,tr=tr,alpha=alpha,plotit=FALSE,
pr=FALSE,pts=pts,skip.crit=TRUE)$output
res.out=min(res[,9])
res.out
}


#' Plot Difference Between Two Groups with Confidence Band - ANCOVA
#'
#' Compares two independent groups using ANCOVA with running interval smoother.
#' Creates a plot showing the difference between groups with a confidence band
#' having simultaneous probability coverage 1-alpha.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#' @param fr1 Span for group 1 running interval smoother (default: 1).
#' @param fr2 Span for group 2 running interval smoother (default: 1).
#' @inheritParams common-params
#' @param pr Logical. Print progress information.
#' @param xout Logical. Remove outliers from covariates before analysis.
#' @param LP Logical. Use lowess/loess smoothing for the difference plot.
#' @param nmin Minimum number of values close to a design point (default: 8).
#' @param scat Logical. Show scatterplot of data.
#' @param xlab,ylab Axis labels for the plot.
#' @param report Logical. If TRUE, return detailed output instead of just plotting.
#' @param ... Additional arguments passed to outfun.
#'
#' @return If report=TRUE, returns list with detailed results. Otherwise creates plot.
#'
#' @note No parametric assumption is made about the form of regression lines.
#'   One covariate only is allowed with this function.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- 0.5*x2 + rnorm(50, mean=1)
#' ancdifplot(x1, y1, x2, y2, alpha=0.05)
#' }
ancdifplot<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,alpha=.05,pr=TRUE,xout=FALSE,outfun=out,LP=TRUE,
nmin=8,scat=TRUE,xlab='X',ylab='Y',report=FALSE,...){
#
# Compare two independent  groups using the ancova method
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#
#  nmin indicates minimun number of values close to a point
#
#  Similar to ancova, only compute a confidence band for the difference and plot it.
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
flago<-outfun(x1,...)$keep
x1<-x1[flago]
y1<-y1[flago]
flag<-outfun(x2,...)$keep
x2<-x2[flago]
y2<-y2[flago]
}
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
flag=vecn>=nmin
ptsum=sum(flag)
est=NA
low=NA
up=NA
ic=0
xp1=NA
xp2=NA
pv=NA
for (i in 1:length(x1)){
if(flag[i]){
g1<-y1[near(x1,x1[i],fr1)]
g2<-y2[near(x2,x2[i],fr2)]
test<-yuen(g1,g2,tr=tr)
ic=ic+1
xp1[ic]=x1[i]
xp2[ic]=x2[i]
est[ic]=test$dif
low[ic]=test$ci[1]
up[ic]=test$ci[2]
pv[ic]=test$p.value
}}
#print(length(pv))
#print(length(xp1))
if(LP){
xy=elimna(cbind(xp1,est,low,up,pv))
est=lplot(xy[,1],xy[,2],plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
up=lplot(xy[,1],xy[,4],plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
low=lplot(xy[,1],xy[,3],plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
if(!report)output='DONE'
plot(c(x1,x2),c(y1,y2),xlab=xlab,ylab=ylab,type='n')
if(!LP){
lines(xp1,up,lty=2)
lines(xp1,low,lty=2)
lines(xp1,est)
if(scat)points(c(x1,x2),c(y1,y2))
if(report){
output=cbind(xp1,est,low,up,pv)
dimnames(output)=list(NULL,c(xlab,'est.dif','lower.ci','upper.ci','p-value'))
}}
if(LP){
lines(xy[,1],up,lty=2)
lines(xy[,1],low,lty=2)
lines(xy[,1],est)
 if(scat)points(c(x1,x2),c(y1,y2))
if(report){
output=cbind(xy[,1],est,low,up,xy[,5])
dimnames(output)=list(NULL,c(xlab,'est.dif','lower.ci','upper.ci','p-value'))
}
}
output
}

#' Compute Effect Sizes for ANCOVA
#'
#' Computes effect sizes at specified design points for comparing two groups
#' while accounting for one or more covariates. Supports multiple effect size
#' measures and can handle both single and multiple covariates.
#'
#' @param x1 Covariate value(s) for group 1 (vector or matrix).
#' @param y1 Response values for group 1.
#' @param x2 Covariate value(s) for group 2 (vector or matrix).
#' @param y2 Response values for group 2.
#' @param method Effect size measure: "EP" (explanatory power), "QS" (quantile shift),
#'   "QStr" (quantile shift with trimmed means), "AKP" (robust Cohen's d),
#'   "WMW" (P(X<Y)), or "KMS" (heteroscedastic analog of Cohen's d).
#' @param fr1 Span for group 1 running interval smoother (default: 1).
#' @param fr2 Span for group 2 running interval smoother (default: 1).
#' @inheritParams common-params
#' @param pts Design points for analysis. If NA, determined automatically based on data depth
#'   (for multiple covariates) or evenly spaced points (for single covariate).
#' @param npts Number of design points to use for single covariate (default: 25).
#' @param nmin Minimum number of values close to a design point (default: 12).
#' @param FRAC Fraction for selecting deep points when pts=NA and multiple covariates (default: 0.2).
#' @inheritParams common-params
#' @param xlab,ylab,zlab Axis labels for plots.
#' @param ticktype Tick mark type for 3D plots (default: 'det').
#'
#' @return Matrix with columns for design points ('pts') and effect sizes ('ES').
#'   For single covariate: 2-column matrix. Creates plot if plotit=TRUE.
#'
#' @note For multiple covariates, design points are chosen based on data depth if not specified.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- 0.5*x2 + rnorm(50, mean=0.5)
#' result <- ancES(x1, y1, x2, y2, method="EP", npts=20)
#' }
ancES<-function(x1,y1,x2,y2,method='EP', fr1=1,fr2=1,tr=.2,pts=NA,SEED=TRUE,npts=NULL,nmin=12,
FRAC=.2,plotit=TRUE,xlab='X1',ylab='X2',zlab='ES',ticktype='det'){
#
#  Compute effect size. If
# plotit=TRUE, plot the results.
#
# pts are the covariate points to be used.
# If more than one covariate, by default,
# design points are chosen based on depth of points in x1 if pts=NA
#
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
if(p!=ncol(x2))stop('x1 and x2 have different number of columns')
if(p==1){
if(is.null(npts))npts=25
if(is.na(pts[1])){
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
LOW<-min(sub[vecn>=nmin])
UP<-max(sub[vecn>=nmin])
pts=seq(x1[LOW],x1[UP],length.out=npts)
}
ES=NA
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
ES[i]=ESfun(g1,g2,method=method,pr=FALSE)
}
res=cbind(pts,ES)
dimnames(res)=list(NULL,c('pts','ES'))
if(plotit)lplot(res[,1],res[,2],xlab=xlab,ylab='ES')
}
if(p>1){
if(SEED)set.seed(2) # now cov.mve always returns same result
ES=NA
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
#
if(is.na(pts[1])){
x1<-as.matrix(x1)
pts<-ancdes(x1,FRAC=FRAC,DH=TRUE)
}
pts<-as.matrix(pts)
n1<-1
n2<-1
vecn<-1
mval1<-cov.mve(x1)
mval2<-cov.mve(x2)
for(i in 1:nrow(pts)){
n1[i]<-length(y1[near3d(x1,pts[i,],fr1,mval1)])
n2[i]<-length(y2[near3d(x2,pts[i,],fr2,mval2)])
}
flag<-rep(TRUE,nrow(pts))
for(i in 1:nrow(pts))if(n1[i]<10 || n2[i]<10)flag[i]<-F
flag=as.logical(flag)
pts<-pts[flag,]
if(sum(flag)==1)pts<-t(as.matrix(pts))
if(sum(flag)==0)stop('No comparable design points found, might increase span.')
for (i in 1:nrow(pts)){
g1<-y1[near3d(x1,pts[i,],fr1,mval1)]
g2<-y2[near3d(x2,pts[i,],fr2,mval2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
ES[i]=ESfun(g1,g2,method=method,tr=tr,pr=FALSE)
}
if(p==2){
if(plotit) lplot(pts,ES,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype)
}
res=cbind(pts,ES)
vlabs=NA
for(j in 1:ncol(pts))vlabs[j]=paste('X',j)
dimnames(res)=list(NULL,c(vlabs,'ES'))
}
res
}



#' @title ANCOVA Effect Size as Function of Covariates (Alternative Version)
#'
#' @description
#' Compute effect size as a function of the covariates for two independent groups.
#' This is an alternative implementation with additional plotting options.
#'
#' @param x1 Covariate values for group 1 (vector or matrix).
#' @param y1 Outcome values for group 1.
#' @param x2 Covariate values for group 2 (vector or matrix).
#' @param y2 Outcome values for group 2.
#' @param method Effect size method: 'EP' (explanatory power), 'AKP' (trimmed-Winsorized Cohen's d),
#'   'QS' (quantile shift), or 'WMW' (P(Y1<Y2|X)). Default is 'EP'.
#' @param fr1 Span for group 1 in running interval smoother. Default is 1.
#' @param fr2 Span for group 2 in running interval smoother. Default is 1.
#' @param tr Trim proportion for 'AKP' method. Default is 0.2.
#' @param pts Covariate points at which to compute effect size. If NA (default), points are chosen automatically.
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#' @param npts Number of design points for single covariate. Default is NULL (uses 25).
#' @param nmin Minimum sample size in each group at each design point. Default is 12.
#' @param pch Plotting character for points. Default is '*'.
#' @param FRAC Fraction for selecting design points with multiple covariates. Default is 0.2.
#' @param plotit Logical. If TRUE (default), plot the effect size function.
#' @param xlab Label for x-axis. Default is 'X1'.
#' @param ylab Label for y-axis. Default is 'X2'.
#' @param zlab Label for z-axis (effect size). Default is 'ES'.
#' @param ticktype Type of tick marks for 3D plots. Default is 'det'.
#'
#' @return A matrix with covariate values and corresponding effect sizes.
#'
#' @export
#' @keywords ancova effect-size
# ancES
ancES<-function(x1,y1,x2,y2,method='EP', fr1=1,fr2=1,tr=.2,pts=NA,SEED=TRUE,npts=NULL,nmin=12,pch='*',
FRAC=.2,plotit=TRUE,xlab='X1',ylab='X2',zlab='ES',ticktype='det'){
#
# ANCOVA
#  Compute effect size as a function of the covariates. If
# plotit=TRUE, plot the results.
#
# pts are the covariate points to be used.
# If more than one covariate, by default,
# design points are chosen based on depth of points in x1 if pts=NA. The
#  range of values when there is a single covariate follows the approach used by the function ancova.
#
#  method:  four choices can be used.
#
#   'EP': explanatory power.
#  'AKP': trimmed-Winsorized analog of Cohen's d
#   'QS': quantile shift
#   'WMW': P(Y1<Y2|X),
#
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
if(p!=ncol(x2))stop('x1 and x2 have different number of columns')
if(p==1){
if(is.null(npts))npts=25
if(is.na(pts[1])){
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
LOW<-min(sub[vecn>=nmin])
UP<-max(sub[vecn>=nmin])
pts=seq(x1[LOW],x1[UP],length.out=npts)
}
ES=NA
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
ES[i]=ESfun(g1,g2,method=method,pr=FALSE,tr=tr)
}
res=cbind(pts,ES)
dimnames(res)=list(NULL,c('pts','ES'))
if(plotit){
vp=0
if(method=='AKP' || method=='QS')vp=-1
plot(c(res[1,1],res[1,1],res[,1]),c(1,vp,res[,2]),xlab=xlab,ylab='ES',type='n')
v=lplot(res[,1],res[,2],xlab=xlab,ylab='ES',plotit=FALSE,pyhat=TRUE)$yhat.values
points(res[,1],res[,2],pch=pch)
lines(res[,1],v)
}}
if(p>1){
if(SEED)set.seed(2) # now cov.mve always returns same result
ES=NA
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
#
if(is.na(pts[1])){
x1<-as.matrix(x1)
pts<-ancdes(x1,FRAC=FRAC,DH=TRUE)
}
pts<-as.matrix(pts)
n1<-1
n2<-1
vecn<-1
mval1<-cov.mve(x1)
mval2<-cov.mve(x2)
for(i in 1:nrow(pts)){
n1[i]<-length(y1[near3d(x1,pts[i,],fr1,mval1)])
n2[i]<-length(y2[near3d(x2,pts[i,],fr2,mval2)])
}
flag<-rep(TRUE,nrow(pts))
for(i in 1:nrow(pts))if(n1[i]<10 || n2[i]<10)flag[i]<-F
flag=as.logical(flag)
pts<-pts[flag,]
if(sum(flag)==1)pts<-t(as.matrix(pts))
if(sum(flag)==0)stop('No comparable design points found, might increase span.')
for (i in 1:nrow(pts)){
g1<-y1[near3d(x1,pts[i,],fr1,mval1)]
g2<-y2[near3d(x2,pts[i,],fr2,mval2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
ES[i]=ESfun(g1,g2,method=method,tr=tr,pr=FALSE)
}
if(p==2){
if(plotit) lplot(pts,ES,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype)
}
res=cbind(pts,ES)
vlabs=NA
for(j in 1:ncol(pts))vlabs[j]=paste('X',j)
dimnames(res)=list(NULL,c(vlabs,'ES'))
}
res
}


#' ANCOVA Effect Size with Confidence Band
#'
#' Computes effect size for two independent groups across covariate values with
#' simultaneous confidence bands. Uses running interval smoother with bootstrap
#' inference and familywise error rate control.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#' @param fr1 Span for running interval smoother for group 1 (default: 1).
#' @param fr2 Span for running interval smoother for group 2 (default: 1).
#' @param method Effect size method: 'WMW' (default), 'AKP', or 'QS'.
#'   \itemize{
#'     \item 'WMW': P(Y1<Y2|X) - probability of superiority
#'     \item 'AKP': Robust analog of Cohen's d (trimmed-Winsorized)
#'     \item 'QS': Quantile shift based on medians
#'   }
#' @inheritParams common-params
#' @param FAST Logical. If TRUE, uses fast computation mode.
#' @param npts Number of design points (default: 25).
#' @param nmin Minimum sample size at each design point (default: 12).
#' @param SCAT Logical. If TRUE, includes scatter plot.
#' @param pc Plot character for points (default: '.').
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#'
#' @return List with components:
#'   \item{output}{Matrix with design points, effect sizes, and confidence intervals}
#'   \item{crit.val}{Critical value for simultaneous inference}
#'
#' @note
#' - Only one covariate is allowed
#' - Confidence intervals have simultaneous probability coverage ~0.95 for npts=25
#' - FWE controlled using method similar to TAP (Wilcox, 2017)
#' - For npts <= 15, uses Bonferroni adjustment instead
#' - Bootstrap samples obtained by resampling (x,y) pairs
#'
#' @references
#' Wilcox, R.R. (2017). Introduction to Robust Estimation and Hypothesis Testing (4th ed.).
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50)
#' y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- x2 + 0.5 + rnorm(50)
#' ancESband(x1, y1, x2, y2, method='WMW', nboot=500)
#' }
ancESband<-function(x1=NULL,y1=NULL,x2=NULL,y2=NULL,fr1=1,fr2=1,method='WMW',
pr=TRUE,FAST=TRUE,alpha=.05,plotit=TRUE,xlab='X',ylab='ES',npts=25,
xout=FALSE,outfun=out,nboot=500,SEED=TRUE,
nmin=12,SCAT=TRUE,pc='.',...){
#
# Compute a measure of effect size for
#  two independent  groups when there is a single covariate.
#
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Confidence intervals are computed so that the simultaneous probability
#  coverage is approximately .95 when npts=25 covariate points are used.
#
#  Three methods can be used:
#
#  'AKP': trimmed-Winsorized analog of Cohen's d
#   'QS': quantile shift
#   'WMW': P(Y1<Y2|X), values for X are stored in pts computed as follows:
#   The lowest and highest points among the covariate value that will be used are determined as in
#  (Wilcox, 2017, Intro to Robust Estimation and Hypothesis Testing, 4th ed., section 12.2).
#  Then npts points equally spaced between the lowest and highest are used.
#
# One covariate only is allowed.
#
#  Bootstrap samples are obtained by resampling
#  from c(x1,y1) and c(x2,y2) rather than conditioning on the x value as done by
#  the function ancova.
#
#
#  Effect sizess are estimated  at each covariate value stored in pts
#  npts: number of  covariate points to be used, which defaults to 25.
#
#  Family-wise error  (FWE) rate is controlled in a manner similar to method TAP (Wilcox, 2017).
#  assuming that FWE=.05, which is based on npts 25 points. If npts<=15 use a Bonferroni method instead.
#
#   x1 y1 are measures for group 1
#   x2 y2 are measures for group 2
#
if(SEED)set.seed(2)
if(method=='QS'){
if(pr){
print('Note: Q.effect is a probabilistic measure of effect size: a shift of the median.')
print(' Cohen d= 0.0 0.2, 0.5 and 0.8 correspond to Q.effect= 0.5, 0.556,  0.638 and  0.714,  respectively.')
print('  Cohen d = -0.2, -0.5 and -0.8; Q.effect=0.45,  0.35 and 0.30, respectively.')
}}

if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(method=='EP')stop('Using method EP not recommended at this time')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
xy1=elimna(cbind(x1,y1))
xy2=elimna=cbind(x2,y2)
x1=xy1[,1]
y1=xy1[,2]
x2=xy2[,1]
y2=xy2[,2]
n1.in=nrow(xy1)
n2.in=nrow(xy2)

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
n1n=length(y1)
n2n=length(y2)
nv=c(30,  50,  60,  70,  80, 100,
150, 200, 300, 400, 500, 600, 800)
pv=c(0.00824497,0.00581, 0.005435089, 0.004763079,
0.00416832, 0.004406774, 0.00388228,0.003812836,0.003812836,0.003453055, 0.003625061,
.003372966, 0.003350022)
p.crit=(lplot.pred(1/nv,pv,1/n1n)$yhat+lplot.pred(1/nv,pv,1/n2n)$yhat)/2
if(alpha!=.05){
p.crit=p.crit*alpha/.05 # A crude adjustment
}
if(npts<=15)p.crit=alpha/npts
qmin=nboot*p.crit
bmin=ceiling(1/p.crit)
if(qmin<1){
stop(paste('nboot must be at least ',bmin))
}
EST=ancES(x1,y1,x2,y2,plotit=FALSE,npts=npts,method=method)
pts=EST[,1]
MAT=matrix(NA,nrow=nboot,ncol=length(pts))
for(i in 1:nboot){
id1=sample(n1n,n1n,replace=TRUE)
id2=sample(n2n,n2n,replace=TRUE)
B=ancES(x1[id1],y1[id1],x2[id2],y2[id2],plotit=FALSE,method=method,pts=pts,npts=npts,SEED=FALSE)
MAT[i,]=B[,2]
}

flag1=MAT<.5
flag2=MAT==.5
pv1=apply(flag1,2,mean,na.rm=TRUE)
pv2=apply(flag2,2,mean,na.rm=TRUE)
pv=pv1+.5*pv2
one.m.pv=1-pv
pv=2*apply(rbind(pv,one.m.pv),2,min)
ci.low=NA
ci.up=NA
qlow=p.crit/2
qhi=1-p.crit/2
ci.low=NA
ci.up=NA
for(i in 1:length(pts)){
if(!is.na(pts[i]))ci.low[i]=qest(MAT[,i],qlow)
if(!is.na(pts[i]))ci.up[i]=qest(MAT[,i],qhi)
}
pvm=matrix(NA,nrow=length(pts),ncol=5)
pvm[,1]=pts
pvm[,2]=EST[,2]
pvm[,3]=pv
pvm[,4]=ci.low
pvm[,5]=ci.up
num.sig=sum(pv<p.crit)
dimnames(pvm)=list(NULL,c('pts','Est.QS','p.values','ci.low','ci.up'))
if(plotit){
plot(c(pts,pts,pts),c(pvm[,c(2,4,5)]),type='n',xlab=xlab,ylab=ylab)
lines(pts,pvm[,2])
lines(pts,pvm[,4],lty=2)
lines(pts,pvm[,5],lty=2)
}

list(output=pvm,n=c(n1.in,n2.in),p.crit=p.crit,num.sig=num.sig)
}

#' ANCOVA Effect Size Summary
#'
#' Computes multiple effect size measures at design points where regression lines
#' are compared. Uses ES.summary.CI to provide comprehensive effect size information
#' at each comparison point.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#' @param fr1 Span for running interval smoother for group 1 (default: 1).
#' @param fr2 Span for running interval smoother for group 2 (default: 1).
#' @inheritParams common-params
#' @param pts Design points. If NA, uses 5 empirically determined points.
#' @param nmin Minimum sample size at each design point (default: 12).
#' @param NULL.V Null values for effect size measures (default: c(0,0,0.5,0.5,0.5,0)).
#' @param REL.M Relative magnitude parameter for ES.summary.CI (default: NULL).
#' @param n.est Number of estimates for simulation (default: 1e6).
#'
#' @return List with components:
#'   \item{n1}{Sample sizes at each design point for group 1}
#'   \item{n2}{Sample sizes at each design point for group 2}
#'   \item{pts}{Design points used}
#'   \item{ES.4.Each.pt}{List of effect size results for each design point}
#'
#' @note Only one covariate is allowed.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50)
#' y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50)
#' y2 <- x2 + 0.5 + rnorm(50)
#' anc.ES.sum(x1, y1, x2, y2, nboot=500)
#' }
anc.ES.sum<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,alpha=.05,pts=NA,SEED=TRUE,nboot=1000,
pr=TRUE,xout=FALSE,outfun=out, nmin=12,NULL.V = c(0, 0, 0.5, 0.5, 0.5, 0), REL.M = NULL, n.est = 1e+06,...){
#
#
#  For each point where the regression lines are compared,
#  compute several measures of effect size  via the R function ESsummary.CI
#
#  Results for the ith point are returned in ES.4.Each.pt[[i]]
#
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
A=list()
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
FLAG=TRUE
if(is.na(pts[1])){
FLAG=FALSE
npt<-5
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
isub[1]<-min(sub[vecn>=nmin])
isub[5]<-max(sub[vecn>=nmin])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
pts=NA
n1=NA
n2=NA
for (i in 1:5){
g1<-y1[near(x1,x1[isub[i]],fr1)]
g2<-y2[near(x2,x1[isub[i]],fr2)]
pts[i]=x1[isub[i]]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
n1[i]=length(g1)
n2[i]=length(g2)
A[[i]]=ES.summary.CI(g1,g2,tr=tr,SEED=SEED,alpha=alpha,nboot=nboot,NULL.V=NULL.V, REL.M =REL.M,n.est=n.est)
}}
if(FLAG){
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
A[[i]]=ES.summary.CI(g1,g2,tr=tr,SEED=SEED,alpha=alpha,nboot=nboot,NULL.V=NULL.V, REL.M =REL.M,n.est=n.est)
}}
list(n1=n1,n2=n2,pts=pts,ES.4.Each.pt=A)
}


#' Global Test for ANCOVA - Equal Location Across Covariate Values
#'
#' Performs a global test that measures of location are equal between two groups
#' across all chosen covariate values. Uses bootstrap percentile method with
#' familywise error rate control.
#'
#' @param x1 Covariate values for group 1 (vector or matrix).
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2 (vector or matrix).
#' @param y2 Response values for group 2.
#' @param xout Logical. If TRUE, remove outliers from covariates before analysis.
#' @inheritParams common-params
#' @param p.crit Critical p-value for global test. If NULL, determined by simulation.
#' @param nreps Number of simulation replications for determining p.crit (default: 500).
#' @param pr Logical. Print progress information.
#' @param CR Logical. If TRUE and number of points is 2 or 3, plot confidence region.
#' @param nmin Minimum number of values close to a design point (default: 12).
#' @param pts Design points for analysis. If NULL, function selects five points automatically.
#' @param fr1 Span for group 1 running interval smoother (default: 1).
#' @param fr2 Span for group 2 running interval smoother (default: 1).
#' @param SCAT Logical. Show scatterplot of data.
#' @param pch1,pch2 Plot characters for groups 1 and 2.
#' @param xlab,ylab Axis labels for plots.
#' @param LP Logical. Use lowess/loess smoothing for plots.
#' @param cpp Logical. Use C++ function for critical p-value (requires WRScpp package).
#' @param ... Additional arguments passed to outfun or est.
#'
#' @return List with:
#'   \item{test}{Global test statistic}
#'   \item{p.value}{P-value for global test}
#'   \item{p.crit}{Critical p-value used}
#'   \item{reject}{Logical indicating whether null hypothesis is rejected}
#'   \item{pts}{Design points used in analysis}
#'
#' @note Works well with alpha=0.05. Performance for alpha < 0.05 is uncertain.
#'   Reject null hypothesis if p.value < p.crit.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- 0.5*x2 + rnorm(50, mean=0.5)
#' result <- ancGLOB(x1, y1, x2, y2, nreps=1000)
#' result$reject  # TRUE if groups differ
#' }
ancGLOB<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,est=tmean,p.crit=NULL,nreps=500,
alpha=.05,pr=TRUE,nboot=500,SEED=TRUE,MC=FALSE,CR=FALSE, nmin=12,pts=NULL,fr1=1,fr2=1,plotit=TRUE,SCAT=TRUE,pch1='+',pch2='o',
xlab='X',ylab='Y',LP=TRUE,cpp=FALSE,...){
#
#  Like the function ancova, only performs a global test that the measures of location
#  are equal among all the covariate values that are chosen.
#
#  pts = NULL, the function picks five covariate values.
#  iter=500 means that when the critical p-value is determined, simulations with 500
#  replications are used to determine the critical p-value.
#
#  Reject if the p-value is less than the critical p-value.
#  Works well with alpha=.05. Uncertain about alpha <.05.
#
#  cpp=TRUE, a C++ function is used to determine the critical p-value
#  assuming the library WRScpp has been installed.  This is done as follows:
#  install.packages('devtools')
#  library("devtools")
#  install_github( "WRScpp", "mrxiaohe")
#
#  CR=TRUE: If number of points is two or three, plot 1-alpha confidence region
#
#
if(CR)plotit=FALSE # Can't plot both regression lines and confidence region
if(SEED)set.seed(2)
iter=nreps
pts.flag=is.null(pts)
if(!is.null(pts))cpp=FALSE
x1<-as.matrix(x1)
p1<-ncol(x1)+1
p<-ncol(x1)
if(p>1)stop('Current version is for one independent variable only')
xy<-cbind(x1,y1)
xy<-elimna(xy)
x1<-xy[,1:p]
y1<-xy[,p1]
xy<-cbind(x2,y2)
xy<-elimna(xy)
x2<-xy[,1:p]
y2<-xy[,p1]
if(xout){
m<-cbind(x1,y1)
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
N1=length(y1)
N2=length(y2)
if(is.null(pts[1])){
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
isub[1]<-min(sub[vecn>=nmin])
isub[5]<-max(sub[vecn>=nmin])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
pts=x1[isub]
g1=list()
g2=list()
for (i in 1:5){
g1[[i]]<-y1[near(x1,x1[isub[i]],fr1)]
g2[[i]]<-y2[near(x2,x1[isub[i]],fr2)]
}}
if(!is.null(pts[1])){
if(length(pts)<2)stop('Should have at least two points (use the R function ancova)')
g1=list()
g2=list()
for (i in 1:length(pts)){
g1[[i]]<-y1[near(x1,pts[i],fr1)]
g2[[i]]<-y2[near(x2,pts[i],fr2)]
}
}
p.alpha=NULL
if(is.null(p.crit)){
if(pts.flag){
if(cpp){
library(WRScpp)
ve=ancGLOB_pv_C(N1,N2,est=est,iter=iter,fr1=fr1,fr2=fr2,nboot=nboot,SEED=SEED,...)
v=hd(ve,q=alpha)
}
else{
 v=ancGLOB_pv(N1,N2,est=est,iter=iter,fr1=fr1,fr2=fr2,nboot=nboot,
PRM=FALSE,SEED=SEED,alpha=alpha,xlab=xlab,ylab=ylab,...)$p.crit
}
}
if(!pts.flag)v=ancGLOB_pv_pts(x1,x2,pts=pts,nmin=nmin,iter=iter,est=est,fr1=fr1,fr2=fr2,
nboot=nboot,SEED=SEED,alpha=alpha,MC=MC)$p.crit
}
if(!is.null(p.crit))v=p.crit
res=aov2depth(g1,g2,est=est,SEED=SEED,CR=CR,alpha=v,...)
if(pr)print('Reject if p.test is less than p.crit')
if(plotit)runmean2g(x1,y1,x2,y2,fr=fr1,est=est,xout=FALSE,LP=LP,xlab=xlab,ylab=ylab,
SCAT=SCAT,pch1=pch1,pch2=pch2,...)
list(p.test=res$p.value,p.crit=v,est1=res$est1,est2=res$est2,dif=res$dif,pts=pts,n1=res$n1,n2=res$n2)
}



#' Determine Critical P-Value for ancGLOB
#'
#' Computes the critical p-value for ancGLOB via simulation under the null
#' hypothesis (normal distribution). Uses the Harrell-Davis estimator on
#' simulated p-values.
#'
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param est Estimator function (default: tmean).
#' @param fr1 Span for group 1 (default: 0.8).
#' @param fr2 Span for group 2 (default: 0.8).
#' @inheritParams common-params
#' @param iter Number of simulation iterations (default: 1000).
#' @param nmin Minimum sample size (default: 12).
#' @param PRM Logical. If TRUE, returns all computed p-values.
#' @param pts Design points (default: NULL for automatic selection).
#'
#' @return List with:
#'   \item{p.crit}{Critical p-value for specified alpha}
#'   \item{prm}{All p-values (if PRM=TRUE), otherwise NULL}
#'   \item{ef.iter}{Effective number of iterations}
#'
#' @keywords internal
ancGLOB_pv<-function(n1,n2,est=tmean,fr1=.8,fr2=.8,nboot=500,SEED=TRUE,iter=1000,nmin=12,MC=TRUE,alpha=.05,PRM=FALSE,pts=NULL,...){
#
#  Determine critical p-value when using the function ancGLOB
#  Strategy: generage data from a normal distribution, NULL true
#  compute p-value, repeat
#  iter times (iter=100 is default)
#  (a larger choice for iter is recommended. To reduce execution time use ancGLOB_pv_C
#
# returns:
# p.crit, the critical p-value for the specified alpha value
# if PRM=T, all p-values that were computed.
# ef.iter, the actual number of iterations, which might differ from iter
# due to sample sizes where it makes no sense to compute a p-value
# based on the  generated data.
#
if(SEED)set.seed(45)
bvec=list()
np1=min(c(n1,n2))+1
nmax=max(c(n1,n2))
for(i in 1:iter){
bvec[[i]]=rmul(nmax,p=4)
if(n1!=n2)bvec[[i]][np1:nmax,1:2]=NA
}
if(MC){
prm=mclapply(bvec,ancGLOB_sub2,fr1=fr1,fr2=fr2,est=est,SEED=SEED,nboot=nboot,pts=pts,...)
}
if(!MC)prm=lapply(bvec,ancGLOB_sub2,fr1=fr1,fr2=fr2,est=est,SEED=SEED,nboot=nboot,pts=pts,...)
prm=elimna(as.vector(matl(prm)))
ef.iter=length(prm)
p.crit=hd(prm,alpha)
prm=sort(elimna(prm))
if(!PRM)prm=NULL
list(p.crit=p.crit,prm=prm,ef.iter=ef.iter)
}


#' Determine Critical P-Value for ancGLOB with Specified Design Points
#'
#' Computes critical p-value for ancGLOB when design points are pre-specified.
#' Like ancGLOB_pv but uses actual covariate values from x1 and x2 with
#' specified design points.
#'
#' @param x1 Covariate values for group 1.
#' @param x2 Covariate values for group 2.
#' @param est Estimator function (default: tmean).
#' @param fr1 Span for group 1 (default: 1).
#' @param fr2 Span for group 2 (default: 1).
#' @inheritParams common-params
#' @param iter Number of simulation iterations (default: 1000).
#' @param nmin Minimum sample size (default: 12).
#' @param PRM Logical. If TRUE, returns all computed p-values.
#' @param pts Design points (required). Example: c(1,4,6) compares at X=1, 4, and 6.
#'
#' @return List with:
#'   \item{p.crit}{Critical p-value for specified alpha}
#'   \item{prm}{All p-values (if PRM=TRUE), otherwise NULL}
#'   \item{ef.iter}{Effective number of iterations}
#'
#' @keywords internal
ancGLOB_pv_pts<-function(x1,x2,est=tmean,fr1=1,fr2=1,nboot=500,SEED=TRUE,iter=1000,nmin=12,MC=TRUE,alpha=.05,PRM=FALSE,pts=NULL,...){
#
#  Determine critical p-value when using the function ancGLOB and pts is specified.
#  Strategy: generage data from a normal distribution, NULL true
#  compute p-value, repeat
#  iter times (iter=1000 is default)
#
#  pts is used to indicate the covariate values where comparisons are to be made.
#  Example: pts=c(1,4,6) will compare regression lines at X=1, 4 and 6
# if pts is not specified, the function terminates with an error.
#
#
# returns:
# p.crit, the critical p-value for the specified alpha value
# if PRM=T, all p-values that were computed.
# ef.iter, the actual number of interations, which might differ from iter
# due to sample sizes where it makes no sense to compute a p-value
# based on the  generated data.
#
# Like ancGLOB_pv, only pts is specified and use data in x1 and x2
#
if(is.null(pts[1]))stop('pts is null, use ancGLOB_pv')
x1=elimna(x1)
x2=elimna(x2)
n1=length(x1)
n2=length(x2)

if(SEED)set.seed(45)
bvec=list()
np1=min(c(n1,n2))+1
nmax=max(c(n1,n2))
for(i in 1:iter){
bvec[[i]]=rmul(nmax,p=4)
if(n1!=n2)bvec[[i]][np1:nmax,1:2]=NA
bvec[[i]][1:n1,1]=x1
bvec[[i]][1:n2,3]=x2
}
prm=NA
if(MC){
prm=mclapply(bvec,ancGLOB_sub4,fr1=fr1,fr2=fr2,est=est,SEED=SEED,nboot=nboot,pts=pts,...)
}
#if(!MC)prm=lapply(bvec,ancGLOB_sub4,fr1=fr1,fr2=fr2,est=est,SEED=SEED,nboot=nboot,pts=pts,...)
if(!MC){
for(ij in 1:length(bvec)){
bv=as.matrix(bvec[[ij]])
prm[ij]=ancGLOB_sub4(bv,fr1=fr1,fr2=fr2,est=est,SEED=SEED,nboot=nboot,pts=pts,nmin=nmin,...)
}}
prm=elimna(as.vector(matl(prm)))
ef.iter=length(prm)
p.crit=hd(prm,alpha)
prm=sort(elimna(prm))
if(!PRM)prm=NULL
list(p.crit=p.crit,prm=prm,ef.iter=ef.iter)
}


#' Internal Subroutine for ancGLOB_pv
#'
#' Helper function that extracts p-value from ancGLOB_sub3 for simulation.
#'
#' @param bvec Bootstrap vector data.
#' @param fr1,fr2 Span parameters.
#' @param est Estimator function.
#' @param SEED Logical for random seed.
#' @param nmin Minimum sample size.
#' @param nboot Number of bootstrap samples.
#' @param pts Design points.
#' @param ... Additional arguments.
#'
#' @return P-value from ancGLOB_sub3.
#'
#' @keywords internal
ancGLOB_sub2<-function(bvec,fr1=fr1,fr2=fr2,est=est,SEED=SEED,nmin=12,nboot=nboot,pts=pts,...){
p=ancGLOB_sub3(bvec[,1],bvec[,2],bvec[,3],bvec[,4],est=est,SEED=SEED,fr1=fr1,fr2=fr2,nboot=nboot,
plotit=FALSE,nmin=12,pts=pts,...)$p.value
p
}



#' @title Global ANCOVA Test Helper Function (Version 3)
#'
#' @description
#' Internal helper function for ancGLOB that performs global ANCOVA test at selected design points.
#' Computes test statistics and p-values using bootstrap methods.
#'
#' @inheritParams ancGLOB
#' @param pcrit Deprecated parameter (use p.crit instead).
#' @param p.crit Critical p-value for determining number of design points. Default is NULL.
#' @param iter Number of iterations for determining critical value. Default is 100.
#' @param LP Logical. If TRUE (default), use LOWESS for plotting.
#' @param ... Additional arguments passed to estimation function.
#'
#' @return A list with components:
#'   \item{p.value}{Global test p-value}
#'   \item{est1}{Estimates for group 1 at each design point}
#'   \item{est2}{Estimates for group 2 at each design point}
#'   \item{dif}{Differences between group estimates}
#'   \item{pts}{Design points used}
#'   \item{n1}{Sample sizes for group 1}
#'   \item{n2}{Sample sizes for group 2}
#'
#' @keywords internal
# ancGLOB_sub3
ancGLOB_sub3<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,est=tmean,pcrit=NULL,p.crit=NULL,iter=100,
nboot=500,SEED=TRUE,MC=FALSE,nmin=12,pts=NULL,fr1=1,fr2=1,plotit=TRUE,xlab='X',ylab='Y',LP=TRUE,...){
#
#
if(SEED)set.seed(2)
x1<-as.matrix(x1)
p1<-ncol(x1)+1
p<-ncol(x1)
if(p>1)stop('Current version is for one independent variable only')
xy<-cbind(x1,y1)
xy<-elimna(xy)
x1<-xy[,1:p]
y1<-xy[,p1]
xy<-cbind(x2,y2)
xy<-elimna(xy)
x2<-xy[,1:p]
y2<-xy[,p1]
if(xout){
m<-cbind(x1,y1)
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
N1=length(y1)
N2=length(y2)
if(is.null(pts[1])){
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
isub[1]<-min(sub[vecn>=nmin])
isub[5]<-max(sub[vecn>=nmin])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
pts=x1[isub]
g1=list()
g2=list()
for (i in 1:5){
g1[[i]]<-y1[near(x1,x1[isub[i]],fr1)]
g2[[i]]<-y2[near(x2,x1[isub[i]],fr2)]
}}
if(!is.null(pts[1])){
if(length(pts)<2)stop('Should have at least two points (With one point, use the R function ancova)')
g1=list()
g2=list()
for (i in 1:length(pts)){
g1[[i]]<-y1[near(x1,pts[i],fr1)]
g2[[i]]<-y2[near(x2,pts[i],fr2)]
}
}
n1=lapply(g1,length)
nv=(min(as.vector(matl(n1))))
res=aov2depth(g1,g2,est=est,SEED=SEED,nboot=nboot,...)
if(plotit)runmean2g(x1,y1,x2,y2,nboot=nboot,fr=fr1,est=est,xout=xout,LP=LP,...)
list(p.value=res$p.value,est1=res$est1,est2=res$est2,dif=res$dif,pts=pts,n1=res$n1,n2=res$n2)
}









#' @title Global ANCOVA Test Helper Function (Version 4)
#'
#' @description
#' Internal helper function for ancGLOB. Similar to ancGLOB_sub3 but with slight variations
#' in implementation for bootstrap resampling procedures.
#'
#' @inheritParams ancGLOB_sub3
#'
#' @return A list with test results including p-value, estimates, and design points.
#'
#' @keywords internal
# ancGLOB_sub4
ancGLOB_sub4<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,est=tmean,pcrit=NULL,p.crit=NULL,iter=100,
nboot=500,SEED=TRUE,MC=FALSE,nmin=12,pts=NULL,fr1=1,fr2=1,plotit=TRUE,xlab='X',ylab='Y',LP=TRUE,...){
#
#
if(SEED)set.seed(2)
x1<-as.matrix(x1)
p1<-ncol(x1)+1
p<-ncol(x1)
if(p>1)stop('Current version is for one independent variable only')
xy<-cbind(x1,y1)
xy<-elimna(xy)
x1<-xy[,1:p]
y1<-xy[,p1]
xy<-cbind(x2,y2)
xy<-elimna(xy)
x2<-xy[,1:p]
y2<-xy[,p1]
if(xout){
m<-cbind(x1,y1)
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
N1=length(y1)
N2=length(y2)
if(is.null(pts[1])){
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
isub[1]<-min(sub[vecn>=nmin])
isub[5]<-max(sub[vecn>=nmin])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
pts=x1[isub]
g1=list()
g2=list()
for (i in 1:5){
g1[[i]]<-y1[near(x1,x1[isub[i]],fr1)]
g2[[i]]<-y2[near(x2,x1[isub[i]],fr2)]
}}
if(!is.null(pts[1])){
if(length(pts)<2)stop('Should have at least two points (use the R function ancova)')
g1=list()
g2=list()
for (i in 1:length(pts)){
g1[[i]]<-y1[near(x1,pts[i],fr1)]
g2[[i]]<-y2[near(x2,pts[i],fr2)]
}}
n1=lapply(g1,length)
nv=(min(as.vector(matl(n1))))
res=aov2depth(g1,g2,est=est,SEED=SEED,nboot=nboot,...)
if(plotit)runmean2g(x1,y1,x2,y2,nboot=nboot,fr=fr1,est=est,xout=xout,LP=LP,...)
list(p.value=res$p.value,est1=res$est1,est2=res$est2,dif=res$dif,pts=pts,n1=res$n1,n2=res$n2)
}









#' @title Global ANCOVA Bootstrap Wrapper (Version 4b)
#'
#' @description
#' Internal helper function that wraps ancGLOB_sub5 for bootstrap resampling.
#' Extracts p-value from bootstrap samples.
#'
#' @param bvec Bootstrap sample matrix with columns for x1, y1, x2, y2.
#' @param fr1 Span for group 1.
#' @param fr2 Span for group 2.
#' @param est Estimator function.
#' @param SEED Logical for random seed.
#' @param nmin Minimum sample size.
#' @param nboot Number of bootstrap samples.
#' @param pts Design points.
#' @param ... Additional arguments.
#'
#' @return P-value from ancGLOB_sub5.
#'
#' @keywords internal
# ancGLOB_sub4
ancGLOB_sub4<-function(bvec,fr1=fr1,fr2=fr2,est=est,SEED=SEED,nmin=12,nboot=nboot,pts=pts,...){
p=ancGLOB_sub5(bvec[,1],bvec[,2],bvec[,3],bvec[,4],est=est,SEED=SEED,fr1=fr1,fr2=fr2,nboot=nboot,nmin=12,pts=pts,...)
p
}




#' @title Global ANCOVA Test Helper Function (Version 5)
#'
#' @description
#' Internal helper function for ancGLOB that requires pre-specified design points.
#' Performs global test at the given points without automatic point selection.
#'
#' @inheritParams ancGLOB_sub3
#'
#' @return P-value from the global test.
#'
#' @keywords internal
# ancGLOB_sub5
ancGLOB_sub5<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,est=tmean,pcrit=NULL,p.crit=NULL,iter=100,
nboot=500,SEED=TRUE,MC=FALSE,nmin=12,pts=NULL,fr1=1,fr2=1,xlab='X',ylab='Y',LP=TRUE,...){
#
#
if(is.null(pts))stop('pts should be specified')
if(SEED)set.seed(2)
x1<-as.matrix(x1)
p1<-ncol(x1)+1
p<-ncol(x1)
if(p>1)stop('Current version is for one independent variable only')
xy<-cbind(x1,y1)
xy<-elimna(xy)
x1<-xy[,1:p]
y1<-xy[,p1]
xy<-cbind(x2,y2)
xy<-elimna(xy)
x2<-xy[,1:p]
y2<-xy[,p1]
if(xout){
m<-cbind(x1,y1)
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
N1=length(y1)
N2=length(y2)
if(length(pts)<2)stop('Should have at least two points (With one point, use the R function ancova)')
g1=list()
g2=list()
for (i in 1:length(pts)){
g1[[i]]<-y1[near(x1,pts[i],fr1)]
g2[[i]]<-y2[near(x2,pts[i],fr2)]
}
n1=lapply(g1,length)
nv=(min(as.vector(matl(n1))))
res=aov2depth(g1,g2,est=est,SEED=SEED,nboot=nboot,nmin=nmin,...)$p.value
res
}



#' @title Global Test of Regression Parameter Equality (Parametric Bootstrap)
#'
#' @description
#' Test the hypothesis that all regression parameters are equal for two independent groups
#' using bootstrap estimation of standard errors and a Johansen-type test statistic.
#'
#' @param x1 Covariate values for group 1 (vector or matrix).
#' @param y1 Outcome values for group 1.
#' @param x2 Covariate values for group 2 (vector or matrix).
#' @param y2 Outcome values for group 2.
#' @param regfun Regression function to use. Default is tsreg (Theil-Sen estimator).
#' @param nboot Number of bootstrap samples. Default is 100.
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#' @param xout Logical. If TRUE, remove outliers from covariates. Default is FALSE.
#' @param eout Logical. If TRUE, remove outliers from (X,Y) pairs. Default is FALSE.
#' @param outfun Function to detect outliers. Default is outpro.
#' @param STAND Logical. If TRUE (default), standardize test statistics.
#' @param plotit Logical. If TRUE (default), plot the data with fitted regression lines.
#' @param xlab Label for x-axis. Default is "X".
#' @param ylab Label for y-axis. Default is "Y".
#' @param ISO Logical. If TRUE, ignore intercept and test only slope parameters. Default is FALSE.
#' @param ... Additional arguments passed to regression function.
#'
#' @return Output from reg1way or reg1wayISO containing test statistics and p-values.
#'
#' @export
#' @keywords ancova regression
# ancGpar
ancGpar<-function(x1,y1,x2,y2,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,eout=FALSE,outfun=outpro,STAND=TRUE,plotit=TRUE,xlab="X",ylab="Y",ISO=FALSE,...){
#
#  Test hypothesis that for two independent groups, all regression parameters are equal
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen type test statistic.
#
x1=as.matrix(x1)
p=ncol(x1)
p1=p+1
xy=elimna(cbind(x1,y1))
x1=xy[,1:p]
y1=xy[,p1]
x2=as.matrix(x2)
p=ncol(x2)
p1=p+1
xy=elimna(cbind(x2,y2))
x2=xy[,1:p]
y2=xy[,p1]
if(plotit){
xx1=x1
yy1=y1
xx2=x2
yy2=y2
if(ncol(as.matrix(x1))==1){
if(eout){
flag=outfun(cbind(x1,y1),plotit=FALSE,...)$keep
xx1=x1[flag]
yy1=y1[flag]
flag=outfun(cbind(x2,y2),plotit=FALSE,...)$keep
xx2=x2[flag]
yy2=y2[flag]
}
if(xout){
flag=outfun(xx1,plotit=FALSE,...)$keep
xx1=x1[flag]
yy1=y1[flag]
flag=outfun(xx2,plotit=FALSE,...)$keep
xx2=x2[flag]
yy2=y2[flag]
}
plot(c(xx1,xx2),c(yy1,yy2),type="n",xlab=xlab,ylab=ylab)
points(xx1,yy1)
points(xx2,yy2,pch="+")
abline(regfun(xx1,yy1,...)$coef)
abline(regfun(xx2,yy2,...)$coef,lty=2)
}}
x=list()
y=list()
x[[1]]=x1
x[[2]]=x2
y[[1]]=y1
y[[2]]=y2
if(!ISO)output=reg1way(x,y,regfun=regfun,nboot=nboot,xout=xout,outfun=outfun,SEED=SEED,STAND=STAND,...)
if(ISO)output=reg1wayISO(x,y,regfun=regfun,nboot=nboot,xout=xout,outfun=outfun,SEED=SEED,STAND=STAND,...)
output
}


#' @title Global Test of Regression Parameter Equality (Parallel Processing)
#'
#' @description
#' Test the hypothesis that all regression parameters are equal for two independent groups
#' using bootstrap estimation with parallel processing via mclapply.
#'
#' @inheritParams ancGpar
#'
#' @return Output from reg1wayMC or reg1wayISOMC containing test statistics and p-values.
#'
#' @export
#' @keywords ancova regression multicore
# ancGparMC
ancGparMC<-function(x1,y1,x2,y2,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,eout=FALSE,outfun=outpro,
STAND=TRUE,plotit=TRUE,xlab="X",ylab="Y",ISO=FALSE,...){
#
#  Test hypothesis that for two independent groups, all regression parameters are equal
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen type test statistic.
#
#  ISO=TRUE, ignore intercept, test only the slope parameters.
#
x1=as.matrix(x1)
p=ncol(x1)
p1=p+1
xy=elimna(cbind(x1,y1))
x1=xy[,1:p]
y1=xy[,p1]
x2=as.matrix(x2)
p=ncol(x2)
p1=p+1
xy=elimna(cbind(x2,y2))
x2=xy[,1:p]
y2=xy[,p1]
if(plotit){
xx1=x1
yy1=y1
xx2=x2
yy2=y2
if(ncol(as.matrix(x1))==1){
if(eout){
flag=outfun(cbind(x1,y1),plotit=FALSE,...)$keep
xx1=x1[flag]
yy1=y1[flag]
flag=outfun(cbind(x2,y2),plotit=FALSE,...)$keep
xx2=x2[flag]
yy2=y2[flag]
}
if(xout){
flag=outfun(xx1,plotit=FALSE,...)$keep
xx1=x1[flag]
yy1=y1[flag]
flag=outfun(xx2,plotit=FALSE,...)$keep
xx2=x2[flag]
yy2=y2[flag]
}
plot(c(xx1,xx2),c(yy1,yy2),type="n",xlab=xlab,ylab=ylab)
points(xx1,yy1)
points(xx2,yy2,pch="+")
abline(regfun(xx1,yy1,...)$coef)
abline(regfun(xx2,yy2,...)$coef,lty=2)
}}
x=list()
y=list()
x[[1]]=x1
x[[2]]=x2
y[[1]]=y1
y[[2]]=y2
if(!ISO)output=reg1wayMC(x,y,regfun=regfun,nboot=nboot,xout=xout,outfun=outfun,
SEED=SEED,STAND=STAND,...)
if(ISO)output=reg1wayISOMC(x,y,regfun=regfun,nboot=nboot,xout=xout,outfun=outfun,
SEED=SEED,STAND=STAND,...)
output
}




#' @title Grid-Based ANCOVA for Two Independent Groups
#'
#' @description
#' Split data into grid cells based on two covariates and compare trimmed means
#' between two independent groups within each cell. Provides multiple comparison adjustments.
#'
#' @param x1 Matrix of covariates for group 1 (must have 2 columns).
#' @param y1 Outcome values for group 1.
#' @param x2 Matrix of covariates for group 2 (must have 2 columns).
#' @param y2 Outcome values for group 2.
#' @param alpha Alpha level for confidence intervals. Default is 0.05.
#' @param Qsplit1 Quantile(s) for splitting first covariate. Default is 0.5 (median).
#' @param Qsplit2 Quantile(s) for splitting second covariate. Default is 0.5 (median).
#' @param SV1 Alternative to Qsplit1: specific values for splitting first covariate. Default is NULL.
#' @param SV2 Alternative to Qsplit2: specific values for splitting second covariate. Default is NULL.
#' @param tr Trim proportion for trimmed mean. Default is 0.2.
#' @param PB Logical. If TRUE, use percentile bootstrap for p-values. Default is FALSE.
#' @param est Estimator function. Default is tmean (trimmed mean).
#' @param nboot Number of bootstrap samples if PB=TRUE. Default is 1000.
#' @param CI Logical. If TRUE, compute confidence intervals for effect sizes. Default is FALSE.
#' @param xout Logical. If TRUE, remove outliers from covariates. Default is FALSE.
#' @param outfun Function to detect outliers. Default is outpro.
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#' @param ... Additional arguments.
#'
#' @return A list with components:
#'   \item{GRID.INTERVALS}{Matrix defining the grid cell boundaries}
#'   \item{GRID}{Matrix with test results for each cell (sample sizes, estimates, differences, p-values)}
#'   \item{Effect.Sizes}{List of effect size summaries for each cell}
#'
#' @export
#' @keywords ancova grid
#'
#' @examples
#' \dontrun{
#'   # Split data at medians of two covariates
#'   result <- anc.grid(x1, y1, x2, y2)
#'
#'   # Split at quartiles for first covariate, median for second
#'   result <- anc.grid(x1, y1, x2, y2, Qsplit1=c(.25,.5,.75))
#' }
# anc.grid
anc.grid<-function(x1,y1,x2,y2, alpha=.05,
#IV=c(1,2),
Qsplit1=.5,Qsplit2=.5, SV1=NULL,SV2=NULL,
tr=.2,PB=FALSE,est=tmean,nboot=1000,CI=FALSE,
xout=FALSE,outfun=outpro,SEED=TRUE,...){
#
# Two independent groups.
# Split on two independent variables based on data in x1. Compare the corresponding regions
#
#
# Qsplit: split the independent variable based on the
#         quantiles indicated by Qsplit
#  Example
#   Qsplit1=c(.25,.5,.75)
#   Qsplit2=.5
#   would split based on the quartiles for the first independent variable and the median
#   for the second independent variable
#
#  Alternatively, the data can be split based in values stored in the arguments
#  SV1 and SV2.
#

#  Then test the hypothesis of equal measures of location
#  IV[1]: indicates the column containing the first independent variable to use.
#  IV[2]:  indicates the column containing the second independent variable to use.
#
#  if(length(unique(y)>2))stop('y should be binary')
p=ncol(x1)
if(p==1)stop('There should be two or more independent variables')
p1=p+1
if(p!=ncol(x2))stop('x2 and x1 do not have the same of variables, ncol(x1)!=ncol(x2)')

if(ncol(x1) != 2 || ncol(x2) !=2)stop('Should have two covariates')

xy1<-elimna(cbind(x1,y1))
x1<-xy1[,1:p]
y1<-xy1[,p1]

xy2<-elimna(cbind(x2,y2))
x2<-xy2[,1:p]
y2<-xy2[,p1]
ES=list()
if(xout){
flag<-outfun(x1,plotit=FALSE)$keep
x1<-x1[flag,]
y1<-y1[flag]
flag<-outfun(x2,plotit=FALSE)$keep
x2<-x2[flag,]
y2<-y2[flag]
}
J=length(Qsplit1)+1
K=length(Qsplit2)+1
if(!is.null(SV1))J=length(SV1)+1
if(!is.null(SV2))K=length(SV2)+1

JK=J*K
MAT=matrix(1:JK,J,K,byrow=TRUE)
z=list()
group=list()
N.int=J
N.int2=K

NG=N.int*N.int2
GRID=matrix(NA,NG,9)
GI=matrix(NA,NG,4)  # grid intervals
L1=NULL
L2=NULL
qv=quantile(x1[,1],Qsplit1)
if(!is.null(SV1))qv=SV1
qv=c(min(x1[,1]),qv,max(x1[,1]))
qv2=quantile(x2[,2],Qsplit2)
if(!is.null(SV2))qv2=SV2
qv2=c(min(x2[,2]),qv2,max(x2[,2]))
ic=0
for(j in 1:N.int){
j1=j+1
xsub1.1=binmat(xy1,1,qv[j],qv[j1])  # split, group 1
xsub1.2=binmat(xy2,1,qv[j],qv[j1])  #split, group 2
for(k in 1:N.int2){
k1=k+1
xsub2.1=binmat(xsub1.1,2,qv2[k],qv2[k1])
xsub2.2=binmat(xsub1.2,2,qv2[k],qv2[k1])
ic=ic+1
if(length(xsub2.1[,3])<=7  || length(xsub2.2[,3])<=7)print('Not enough data in one or more  grids')
GI[ic,]=c(qv[j],qv[j1],qv2[k],qv2[k1])

if(length(xsub2.1[,3])>7  || length(xsub2.2[,3])>7){
a=yuen(xsub2.1[,3],xsub2.2[,3],tr=tr,alpha=alpha)
a=pool.a.list(a)
a=a[c(1:4,8,5:7)]
if(PB){
pbv=trimpb2(xsub2.1[,3],xsub2.2[,3],tr=tr,alpha=alpha,nboot=nboot)
pbv=pool.a.list(pbv)
a[6:8]=pbv[c(2,3,1)]
}
GRID[ic,1:8]=a[1:8]
if(!CI)ES[[ic]]=ES.summary(xsub2.1,xsub2.2,tr=tr)
if(CI)ES[[ic]]=ES.summary.CI(xsub2.1,xsub2.2,tr=tr)
}}
}
dimnames(GI)=list(NULL,c('Int.1.low','Int.1.up','Int.2.low','Int.2.up'))
GRID[,9]=p.adjust(GRID[,8],method='hoch')
dimnames(GRID)=list(NULL,c('n1','n2','est.1','est.2','dif','ci.low','ci.up','p.value','adj.p.value'))
list(GRID.INTERVALS=GI,GRID=GRID, Effect.Sizes=ES)
}


#' @title Grid-Based ANCOVA for Binary Outcomes
#'
#' @description
#' Split data into grid cells based on two covariates and compare binary outcomes
#' between two independent groups using various methods (KMS, SK, or ECP).
#'
#' @param x1 Matrix of covariates for group 1 (must have 2 columns).
#' @param y1 Binary outcome values for group 1.
#' @param x2 Matrix of covariates for group 2 (must have 2 columns).
#' @param y2 Binary outcome values for group 2.
#' @param alpha Alpha level for tests. Default is 0.05.
#' @param method Method for comparing binary outcomes: 'KMS', 'SK', or 'ECP'. Default is 'KMS'.
#' @param Qsplit1 Quantile(s) for splitting first covariate. Default is 0.5 (median).
#' @param Qsplit2 Quantile(s) for splitting second covariate. Default is 0.5 (median).
#' @param SV1 Alternative to Qsplit1: specific values for splitting first covariate. Default is NULL.
#' @param SV2 Alternative to Qsplit2: specific values for splitting second covariate. Default is NULL.
#' @param xout Logical. If TRUE, remove outliers from covariates. Default is FALSE.
#' @param outfun Function to detect outliers. Default is outpro.
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#' @param ... Additional arguments.
#'
#' @return A list with grid intervals and test results for each cell.
#'
#' @export
#' @keywords ancova binary grid
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'
#' @examples
#' \dontrun{
#'   # Compare binary outcomes across grid cells
#'   result <- anc.grid.bin(x1, y1, x2, y2, method='KMS')
#' }
# anc.grid.bin
anc.grid.bin<-function(x1,y1,x2,y2, alpha=.05,method='KMS',
Qsplit1=.5,Qsplit2=.5, SV1=NULL,SV2=NULL,
xout=FALSE,outfun=outpro,SEED=TRUE,...){
#
# Two independent groups.
# Split on two independent variables based on data in x1. Compare the corresponding regions
#
#
# Qsplit: split the independent variable based on the
#         quantiles indicated by Qsplit
#  Example
#   Qsplit1=c(.25,.5,.75)
#   Qsplit2=.5
#   would split based on the quartiles for the first independent variable and the median
#   for the second independent variable
#
#  The argument method can be 'KMS, 'SK' or 'ECP'
#  See the 5th edition of Wilcox, Intro to Robust Estimation and Hypothesis Testing
#  details.
#
#  Alternatively, the data can be split based in values stored in the arguments
#  SV1 and SV2.
#
if(identical(method,'ZHZ'))stop('Argument method should be KMS, SK or ECP')
if(length(unique(y1))>2)stop('y1 should be binary')
if(length(unique(y2))>2)stop('y2 should be binary')
p=ncol(x1)
if(p==1)stop('There should be two or more independent variables')
p1=p+1
if(p!=ncol(x2))stop('x2 and x1 do not have the same of variables, ncol(x1)!=ncol(x2)')
if(ncol(x1) != 2 || ncol(x2) !=2)stop('Should have two covariates')
xy1<-elimna(cbind(x1,y1))
x1<-xy1[,1:p]
y1<-xy1[,p1]
xy2<-elimna(cbind(x2,y2))
x2<-xy2[,1:p]
y2<-xy2[,p1]
if(xout){
flag<-outfun(x1,plotit=FALSE)$keep
x1<-x1[flag,]
y1<-y1[flag]
flag<-outfun(x2,plotit=FALSE)$keep
x2<-x2[flag,]
y2<-y2[flag]
}
J=length(Qsplit1)+1
K=length(Qsplit2)+1
if(!is.null(SV1))J=length(SV1)+1
if(!is.null(SV2))K=length(SV2)+1
JK=J*K
MAT=matrix(1:JK,J,K,byrow=TRUE)
z=list()
group=list()
N.int=J
N.int2=K
NG=N.int*N.int2
GRID=matrix(NA,NG,9)
GI=matrix(NA,NG,4)  # grid intervals
L1=NULL
L2=NULL
qv=quantile(x1[,1],Qsplit1)
if(!is.null(SV1))qv=SV1
qv=c(min(x1[,1]),qv,max(x1[,1]))
qv2=quantile(x2[,2],Qsplit2)
if(!is.null(SV2))qv2=SV2
qv2=c(min(x2[,2]),qv2,max(x2[,2]))
ic=0
for(j in 1:N.int){
j1=j+1
xsub1.1=binmat(xy1,1,qv[j],qv[j1])  # split, group 1
xsub1.2=binmat(xy2,1,qv[j],qv[j1])  #split, group 2
for(k in 1:N.int2){
k1=k+1
xsub2.1=binmat(xsub1.1,2,qv2[k],qv2[k1])
xsub2.2=binmat(xsub1.2,2,qv2[k],qv2[k1])
ic=ic+1
if(length(xsub2.1[,3])<=7  || length(xsub2.2[,3])<=7)print('Not enough data in one or more  grids')
GI[ic,]=c(qv[j],qv[j1],qv2[k],qv2[k1])

if(length(xsub2.1[,3])>7  || length(xsub2.2[,3])>7){
a=binom2g(sum(xsub2.1[,3]),length(xsub2.1[,3]),
sum(xsub2.2[,3]),length(xsub2.2[,3]), method=method,alpha=alpha)
if(identical(method,'KMS')){
a=pool.a.list(a)
#print(a)
a=c(length(xsub2.1[,3]),length(xsub2.2[,3]),a[c(3:5,1:2,6)])
}
if(identical(method,'SK')){
a=c(length(xsub2.1[,3]),length(xsub2.2[,3]),a$p1,a$p2,a$p1-a$p2,NA,NA,a$p.value)
}
if(identical(method,'ECP')){
a=c(length(xsub2.1[,3]),length(xsub2.2[,3]),a$output[1,3:8])
}
GRID[ic,1:8]=a[1:8]
}}
}
dimnames(GI)=list(NULL,c('Int.1.low','Int.1.up','Int.2.low','Int.2.up'))
GRID[,9]=p.adjust(GRID[,8],method='hoch')
dimnames(GRID)=list(NULL,c('n1','n2','est.1','est.2','dif','ci.low','ci.up','p.value','adj.p.value'))
list(GRID.INTERVALS=GI,GRID=GRID)
}


# ancJN
#' Johnson-Neyman ANCOVA Method
#'
#' Compares two independent groups using the Johnson-Neyman technique, which identifies
#' regions of the covariate where group differences are statistically significant. Uses
#' robust regression (default: Theil-Sen) and provides simultaneous confidence intervals.
#'
#' @param x1 Covariate values for group 1 (must be single covariate).
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2 (must be single covariate).
#' @param y2 Response values for group 2.
#' @param pts Design points where regression lines are compared. If NULL, determined automatically.
#' @param Npts Number of design points to use (default: 5).
#' @param Dpts Logical. If FALSE, points are uniformly spaced between min and max values.
#'   If TRUE, points chosen as in ancova function.
#' @param regfun Regression function to use (default: tsreg for Theil-Sen).
#' @param fr1 Span for group 1 running interval smoother (default: 1).
#' @param fr2 Span for group 2 running interval smoother (default: 1).
#' @param SCAT Logical. Show scatterplot of data.
#' @param pch1,pch2 Plot characters for groups 1 and 2.
#' @param xlab,ylab Axis labels for plots.
#' @param xout Logical. Remove outliers from covariates before analysis.
#' @inheritParams common-params
#' @param ... Additional arguments passed to regfun or outfun.
#'
#' @return Matrix with columns: X (design points), Est1, Est2 (estimates), DIF (difference),
#'   TEST (test statistic), se (standard error), ci.low, ci.hi, p.value, p.adjust.
#'
#' @note Only one covariate is allowed. For multiple covariates, use ancJNmp.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- 2*x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- x2 + rnorm(50, mean=1)
#' result <- ancJN(x1, y1, x2, y2, Npts=7)
#' # Check where groups differ significantly
#' result[result[,'p.value'] < 0.05, ]
#' }
ancJN<-function(x1,y1,x2,y2,pts=NULL,Npts=5,Dpts=FALSE,regfun=tsreg,fr1=1,fr2=1,SCAT=TRUE,pch1='*',pch2='+',
alpha=.05,plotit=TRUE,xout=FALSE,outfun=out,nboot=100,SEED=TRUE,xlab='X',ylab='Y',...){
#
# Compare the regression lines of two independent groups at specified design points
# using a robust regression estimator.
# By default, use the Theil--Sen estimator
#
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  Dpts=FALSE: Five covariate points are chosen uniformly space between the smallest and largest
#                         values observed.
# Dpts=TRUE:  Five covariate points are chosen in the same manner as done by the function ancova
#
# Npts: number of points used
#
if(identical(outfun,boxplot))stop('Use outfun=outbox')
if(SEED)set.seed(2)
FLAG=pts
xy=elimna(cbind(x1,y1))
if(ncol(xy)>2)stop('Only one covariate is allowed. Use ancJNmp')
x1=xy[,1]
y1=xy[,2]
nv1=length(y1)
xy=elimna(cbind(x2,y2))
if(ncol(xy)>2)stop('Only one covariate is allowed. Use ancJNmp')
x2=xy[,1]
y2=xy[,2]
nv2=length(y2)
if(xout){
m<-cbind(x1,y1)
p1=ncol(m)
p=p1-1
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
if(is.null(pts[1])){
xall=unique(c(x1,x2))
pts=seq(min(xall),max(xall),length.out=Npts)
if(Dpts){
npt<-5
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
pts=x1[isub]
}
npts=length(pts)
mat<-matrix(NA,npts,10)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value','adj.p.values'))
mat[,1]=pts
sqsd1=regYvar(x1,y1,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
sqsd2=regYvar(x2,y2,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
est1=regYhat(x1,y1,xr=pts,regfun=regfun,xout=FALSE,outfun=outfun,...)
est2=regYhat(x2,y2,xr=pts,regfun=regfun,xout=FALSE,outfun=outfun,...)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd1+sqsd2)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,9]=pv
crit<-smmcrit(Inf,5)
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
if(!is.null(FLAG)){
n1=NA
n2=NA
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),10)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value','p.adjust'))
mat[,1]<-pts
sqsd1=regYvar(x1,y1,regfun=regfun,pts=pts,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
sqsd2=regYvar(x2,y2,regfun=regfun,pts=pts,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
est1=regYhat(x1,y1,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
est2=regYhat(x2,y2,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd1+sqsd2)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,9]=pv
crit<-smmcrit(Inf,length(pts))
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
reg1=regfun(x1,y1,...)$coef
reg2=regfun(x2,y2,...)$coef
if(plotit){
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
plot(c(x1,x2),c(y1,y2),type='n',xlab=xlab,ylab=ylab)
if(SCAT){
points(x1,y1,pch=pch1)
points(x2,y2,pch=pch2)
}
abline(reg1)
abline(reg2,lty=2)
}
mat[,10]=p.adjust(mat[,9],method='hoch')
list(n=c(nv1,nv2),intercept.slope.group1=reg1,intercept.slope.group2=reg2,output=mat)
}


#' @title Johnson-Neyman ANCOVA for Linear Contrasts
#'
#' @description
#' Perform Johnson-Neyman ANCOVA for linear contrasts among J independent groups
#' using robust regression (default: Theil-Sen estimator). Compares regression lines
#' at specified design points.
#'
#' @param x Covariate values in list or matrix form (J groups).
#' @param y Outcome values in list or matrix form (J groups).
#' @param pts Design points where regression lines are compared. If NULL, points are chosen automatically.
#' @param con Contrast matrix for linear contrasts. If NULL, uses default contrasts.
#' @param regfun Regression function to use. Default is tsreg (Theil-Sen estimator).
#' @param nmin Minimum sample size required at each design point. Default is 12.
#' @param npts Number of equally-spaced design points if pts=NULL. Default is 5.
#' @param alpha Alpha level for confidence intervals. Default is 0.05.
#' @param xout Logical. If TRUE, remove outliers from covariates. Default is FALSE.
#' @param outfun Function to detect outliers. Default is out.
#' @param nboot Number of bootstrap samples for variance estimation. Default is 100.
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#' @param pr Logical. If TRUE (default), print results.
#' @param ... Additional arguments passed to regression function.
#'
#' @return A list with test results for linear contrasts at each design point.
#'
#' @details
#' For adjusted p-values controlling FWE, critical values from simulation studies
#' can be obtained from files p6n50.csv and p6n100.csv at the Wilcox lab website.
#'
#' @export
#' @keywords ancova contrasts
#'
#' @examples
#' \dontrun{
#'   # Three groups
#'   x <- list(rnorm(30), rnorm(30), rnorm(30))
#'   y <- list(2*x[[1]]+rnorm(30), x[[2]]+rnorm(30), x[[3]]+2+rnorm(30))
#'   result <- ancJN.LC(x, y, npts=7)
#' }
# ancJN.LC
ancJN.LC<-function(x,y,pts=NULL,con=NULL,regfun=tsreg,nmin=12,npts=5,
alpha=.05,xout=FALSE,outfun=out,nboot=100,SEED=TRUE,pr=TRUE,...){
#
# ANCOVA: Linear contrasts
#  J independent groups
# using a robust regression estimator.
# By default, use the Theil--Sen estimator
#
#  Assume data are in
#   x and y: list mode with length J or matrices with J columns
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  npts=5 points are used that are  equally spaced
#  npts=25 would use 25 points equally spaced.
#
#  To get adjusted p-values that control FWE, set
#    p6n50=p5n50 and p6n100=p6n100
#     where the R variables p6n5  and p6n100 contain the data in the files
#    p6n5.csv  and p6n100.csv, which are stored at
#    https://dornsife.usc.edu/cf/labs/wilcox/wilcox-faculty-display.cfm
#     in the directory labeled datasets
#     These  adjusted p-values are based on an estimate of the null distribution of the p-values using 10,000 replications.
#
#
if(identical(outfun,boxplot))stop('Use outfun=outbox')
if(SEED)set.seed(2)
FLAG=pts
if(is.matrix(x) || is.data.frame(x))x=listm(x)
if(is.matrix(y) || is.data.frame(y))y=listm(y)
J=length(x)
if(is.null(con))con=con.all.pairs(J)
con=as.matrix(con)

for(j in 1:J){
xy=cbind(x[[j]],y[[j]])
xy=elimna(xy)
x[[j]]=xy[,1]
y[[j]]=xy[,2]
}

LEV=seq(.001,.1,.001)
#  Critical p-values:
PC=c( 0.0004715554, 0.0010600032, 0.0014646097, 0.0018389954, 0.0022100770, 0.0026392236
 , 0.0030845986, 0.0034799136, 0.0039374557, 0.0045480591, 0.0052206792, 0.0057883125
 , 0.0062902648, 0.0068322005, 0.0074001036, 0.0079687149, 0.0085694544, 0.0091944370
 , 0.0097899444, 0.0103380231, 0.0108787002, 0.0114575017, 0.0120748898, 0.0127022411
 , 0.0133151567, 0.0138929664, 0.0144234632, 0.0149174499, 0.0153969778, 0.0158737390
 , 0.0163486640, 0.0168226766, 0.0172990400, 0.0177788317, 0.0182590407, 0.0187350584
 , 0.0192063600, 0.0196827391, 0.0201828045, 0.0207226980, 0.0213060840, 0.0219224603
 , 0.0225496919, 0.0231594943, 0.0237279940, 0.0242453799, 0.0247167117, 0.0251561867
 , 0.0255817130, 0.0260113587, 0.0264602401, 0.0269376359, 0.0274445303, 0.0279729461
 , 0.0285099972, 0.0290459058, 0.0295790691, 0.0301128895, 0.0306482248, 0.0311800075
 , 0.0317011819, 0.0322091285, 0.0327081943, 0.0332067908, 0.0337124366, 0.0342287456
 , 0.0347555836, 0.0352908030, 0.0358314032, 0.0363734449, 0.0369117409, 0.0374406811
 , 0.0379564126, 0.0384591067, 0.0389535674, 0.0394475575, 0.0399489437, 0.0404633647
 , 0.0409932052, 0.0415375872, 0.0420929727, 0.0426543011, 0.0432164300, 0.0437752831
 , 0.0443282924, 0.0448743464, 0.0454136998, 0.0459479133, 0.0464794759, 0.0470108970
 , 0.0475436440, 0.0480776649, 0.0486119396, 0.0491457340, 0.0496797020, 0.0502161233
 , 0.0507581656, 0.0513086116, 0.0518687047, 0.0524377104)
 PC100=c( 0.0002966929, 0.0007169931, 0.0011232960, 0.0014665833, 0.0018499989, 0.0021811159, 0.0025574707, 0.0030353289, 0.0035106136,
0.0039372118, 0.0043468911, 0.0047755738, 0.0052219494, 0.0056692252, 0.0061033036, 0.0065103109, 0.0068936246, 0.0072767310
, 0.0076766388, 0.0080950473, 0.0085363254, 0.0090079848, 0.0095103514, 0.0100413399, 0.0105923616, 0.0111373687, 0.0116501421
, 0.0121316345, 0.0126067868, 0.0130940807, 0.0135865147, 0.0140673848, 0.0145340038, 0.0149959496, 0.0154605814, 0.0159310004
, 0.0164110067, 0.0169026157, 0.0173999561, 0.0178905940, 0.0183633096, 0.0188126308, 0.0192375443, 0.0196398972, 0.0200251178
, 0.0204024531, 0.0207816755, 0.0211673500, 0.0215562764, 0.0219417240, 0.0223201432, 0.0226935235, 0.0230670390, 0.0234465536
, 0.0238377401, 0.0242449082, 0.0246692644, 0.0251088323, 0.0255604045, 0.0260209857, 0.0264872622, 0.0269547078, 0.0274182559
, 0.0278743138, 0.0283226660, 0.0287673658, 0.0292163846, 0.0296798737, 0.0301671058, 0.0306828917, 0.0312249485, 0.0317836349
, 0.0323445314, 0.0328929974, 0.0334187668, 0.0339185027, 0.0343952967, 0.0348558545, 0.0353073007, 0.0357552402, 0.0362033565
, 0.0366537740, 0.0371074184, 0.0375642517, 0.0380236755, 0.0384852821, 0.0389497084, 0.0394190068, 0.0398960179, 0.0403828416
, 0.0408792748, 0.0413823211, 0.0418872318, 0.0423894946, 0.0428866160, 0.0433788246, 0.0438685719, 0.0443592827, 0.0448539520
, 0.0453540386)

LV10=seq(.11,.99,.01)

PC2=c( 0.05847911,0.06353127,0.06846620,0.07310904,0.07881193,0.08499317
 ,0.09026885,0.09574866,0.10124244,0.10716984,0.11236755,0.11770113
 ,0.12277641,0.12816436,0.13376686,0.13841793,0.14410387,0.14968273
,0.15593615,0.16189668,0.16835575,0.17408581,0.17990687,0.18615084
,0.19255775,0.19789818,0.20461883,0.21144604,0.21805108,0.22337111
,0.22987129,0.23652924,0.24281857,0.24931816,0.25527023,0.26170124
,0.26891556,0.27539636,0.28170631,0.28851563,0.29594635,0.30318192
,0.31045830,0.31845106,0.32544425,0.33439997,0.34149133,0.34720162
,0.35473714,0.36186855,0.36966857,0.37744663,0.38495268,0.39335818
,0.40139728,0.40976955,0.41707974,0.42741844,0.43614133,0.44525646
,0.45212461,0.46081343,0.46949317,0.47873498,0.48961349,0.50015426
,0.50931802,0.52037963,0.53148887,0.54267461,0.55344746,0.56735092
,0.57914832,0.58912604,0.60177045,0.61639495,0.62953494,0.64515539
,0.66071262,0.67259936,0.68883014,0.70680224,0.72531254,0.74408131
,0.76444699,0.79049224,0.81775009,0.84817210,0.88710556)

PC100v2=c( 0.05001203,0.05433510,0.05923540,0.06449529,0.06983196,0.07506101
 ,0.07980333,0.08542380,0.09160087,0.09662086,0.10166432,0.10676042
,0.11225467,0.11892794,0.12375516,0.12925937,0.13452963,0.14002073
,0.14507355,0.15091410,0.15606790,0.16249155,0.16797190,0.17373478
,0.17882609,0.18493900,0.19193563,0.19735772,0.20348841,0.20935194
,0.21635380,0.22225273,0.22983262,0.23667483,0.24413437,0.25020140
,0.25702517,0.26403520,0.27014313,0.27804532,0.28442682,0.29088905
,0.29821809,0.30537715,0.31289211,0.32086820,0.32832301,0.33561776
,0.34266364,0.34882764,0.35668900,0.36457176,0.37183425,0.37923028
,0.38652374,0.39483644,0.40189416,0.41095361,0.42066139,0.43029008
,0.43909925,0.44914388,0.45986646,0.47000245,0.48024114,0.48927833
,0.50263045,0.51392305,0.52498817,0.53697186,0.54793443,0.55984674
,0.57110087,0.58246374,0.59547592,0.61054637,0.62515843,0.63817892
,0.65448463,0.67194322,0.69024711,0.70584982,0.72322165,0.74211690
,0.76534343,0.78913944,0.81819273,0.84990380,0.89251684)

LV=c(LEV,LV10)
PC50=c(PC,PC2)
PC100.all=c(PC100,PC100v2)
n=lapply(y,length)
n=as.vector(matl(n))
nmin=min(n)
if(nmin<=75)cp4=lplot.pred(LV,PC50,alpha)$yhat
else
cp4=lplot.pred(LV,PC100.all,alpha)$yhat
crit=qnorm(1-cp4/2)
if(xout){
for(j in 1:J){
flag=outfun(x[[j]],plotit=FALSE,...)$keep
m<-cbind(x[[j]],y[[j]])
p1=ncol(m)
p=p1-1
m<-m[flag,]
x[[j]]<-m[,1:p]
y[[j]]<-m[,p1]
}}
if(!is.null(pts))npts=length(pts)
if(is.null(pts[1])){
xall=lapply(x,unique)
L=lapply(xall,min)
U=lapply(xall,max)
L=matl(L)
U=matl(U)
L=max(L)
U=min(U)
if(L>=U)stop('The range of covariate values is not sufficiently similar among the groups')
pts=seq(max(L),min(U),length.out=npts)
}
NT=ncol(con)
CON=list()
mat<-matrix(NA,npts,8)
dimnames(mat)<-list(NULL,c('X','Est','TEST','se','ci.low','ci.hi','p.value','Adj.p.value'))
mat[,1]=pts
sqsd=list()
est=list()
for(j in 1:J){
sqsd[[j]]=regYvar(x[[j]],y[[j]],pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
est[[j]]=regYhat(x[[j]],y[[j]],xr=pts,regfun=regfun,xout=FALSE,outfun=outfun,...)
}
est=matl(est)
sqsd=matl(sqsd)
for(K in 1:NT){
for(k in 1:npts){
EST=0
SD=0
EST=sum(con[,K]*est[k,])
SD=sum(con[,K]^2*sqsd[k,])
sd=sqrt(SD)
mat[k,4]=sd
tests=EST/sd
mat[k,3]=tests
pv=2*(1-pnorm(abs(tests)))
mat[k,7]=pv
mat[k,5]=EST-crit*sd
mat[k,6]=EST+crit*sd
mat[k,2]=EST
mat[k,8]=NA
# Compute a p-value
if(nmin<=75){
flag=mat[k,7]>=PC50
ID=which(flag==TRUE)
ic=max(ID,1)
mat[k,8]=LV[ic]
}
else{
flag=mat[k,7]>=PC100.all
ID=which(flag==TRUE)
ic=max(ID,1)
mat[k,8]=LV[ic]
}
}
CON[[K]]=mat
}
pts=as.matrix(pts,ncol=1)
g.est=cbind(pts,est)
LAB='X'
J1=J+1
for(j in 2:J1)LAB[j]=paste('GRP',j-1)
dimnames(g.est)=list(NULL,LAB)
list(n=n,crit.p.value=cp4,CON=CON,con=con,GRP.est=g.est)
}


#' @title Johnson-Neyman ANCOVA for Multiple Covariates (Deprecated)
#'
#' @description
#' **DEPRECATED: This function has been replaced by ancJNPVAL.**
#'
#' Compare two independent groups using a generalization of Johnson-Neyman ANCOVA
#' that allows multiple covariates.
#'
#' @param x1 Matrix of covariates for group 1.
#' @param y1 Outcome values for group 1.
#' @param x2 Matrix of covariates for group 2.
#' @param y2 Outcome values for group 2.
#' @param regfun Regression function. Default is qreg (quantile regression).
#' @param p.crit Critical p-value for testing. Default is NULL (computed automatically).
#' @param DEEP Logical. If TRUE, use deepest half of covariate points. Default is FALSE.
#' @param WARN Logical. Display warnings. Default is FALSE.
#' @param plotit Logical. If TRUE, plot results for 2 covariates. Default is TRUE.
#' @param xlab,ylab Axis labels for plots.
#' @param null.value Null hypothesis value for difference. Default is 0.
#' @param FRAC Fraction parameter for depth calculations. Default is 0.5.
#' @param cov1 Logical. If TRUE, use covariate points from x1. Default is FALSE.
#' @param SMM Logical. Use Studentized maximum modulus distribution. Default is TRUE.
#' @param ALL Logical. If TRUE, return all results; if FALSE, only significant results. Default is TRUE.
#' @param pr Logical. Print messages. Default is TRUE.
#' @param alpha Alpha level. Default is 0.05.
#' @param nreps Number of replications for computing critical value. Default is 1000.
#' @param MC Logical. Use parallel processing. Default is FALSE.
#' @param pts Design points. If NULL, chosen automatically.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param nboot Number of bootstrap samples. Default is 100.
#' @param xout Logical. Remove outliers from covariates. Default is FALSE.
#' @param outfun Outlier detection function. Default is outpro.
#' @param ... Additional arguments.
#'
#' @return This function stops with an error message directing users to ancJNPVAL.
#'
#' @export
#' @keywords ancova deprecated
# ancJNmp
ancJNmp<-function(x1,y1,x2,y2,regfun=qreg,p.crit=NULL,DEEP=FALSE,WARN=FALSE,
plotit=TRUE,xlab='X1',ylab='X2',null.value=0,FRAC=.5,cov1=FALSE,SMM=TRUE,ALL=TRUE,pr=TRUE,
alpha=.05,nreps=1000, MC=FALSE, pts=NULL,SEED=TRUE,nboot=100,xout=FALSE,outfun=outpro,...){
#
# Compare two independent  groups using a generalization of the ancts function that
#  allows more than one covariate.
#
# DEEP=FALSE: If pts=NULL, design points are chosen to be deepest  point in
# x1  plus points on the .5 depth contour.
#
# DEEP=TRUE, choose deepest half of c(x1,x2) and use critical p-value indicated by
# p.crit, the critical p-value,which defaults to .015 when alpha=.05.
# If alpha!=.05, p.crit must be computed, which can require high execution time.
# MC=TRUE will reduce execution time considerably.
#
#  cov1=TRUE: the covariates that are used are taken to be the points in x1. If
#
#  plotit=TRUE: if p=2 covariates, plot covariate points with non-significant points indicated by * and
#  significant points by +

# (This function replaces anctsmp, which does not have an option for using the deepest half of covariate points.)
#
if(SEED)set.seed(2) # now cov.mve always returns same result

stop('This function has been replaced by an improved ancJNPVAL')
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have a different number of columns')
p=ncol(x1)
p1=p+1
m1=elimna(cbind(x1,y1))
x1=m1[,1:p]
y1=m1[,p1]
m2=elimna(cbind(x2,y2))
x2=m2[,1:p]
y2=m2[,p1]
#
if(xout){
m<-cbind(x1,y1)
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
nv=c(n1,n2)
if(cov1){
pts=unique(x1)
if(alpha==.05)p.crit=0.00676
}
if(DEEP)pts=NULL
if(!is.null(pts[1])){
p.crit=NULL
DEEP=FALSE
}
if(is.null(pts[1])){
if(!DEEP){
x1<-as.matrix(x1)
pts<-ancdes(unique(rbind(x1,x2)))
p.crit=NULL
}
if(DEEP){
pts=ancov2COV(x1,y1,x2,y2,DETAILS=TRUE,cr=.27,pr=FALSE,FRAC=FRAC)$all.points.used[,1:2]
}}
pts<-as.matrix(pts)
ntests=nrow(pts)
mat<-matrix(NA,ntests,8)
dimnames(mat)<-list(NULL,c('Est 1', 'Est 2','DIF','TEST','se','ci.low','ci.hi','p.value'))
sqsd1=regYvar(x1,y1,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
sqsd2=regYvar(x2,y2,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
#  xout=F because leverage points have already been removed.
est1=regYhat(x1,y1,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
est2=regYhat(x2,y2,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
mat[,1]=est1
mat[,2]=est2
est=est1-est2
mat[,3]=est
sd=sqrt(sqsd1+sqsd2)
mat[,5]=sd
tests=(est1-est2)/sd
mat[,4]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,8]=pv
crit=NULL
if(!cov1){
if(!DEEP){
if(ntests==1)crit=qnorm(1-alpha/2)
if(length(pts)>1){
if(ntests<=28){
if(alpha==.05)crit<-smmcrit(Inf,ntests)
if(alpha==.01)crit<-smmcrit01(Inf,ntests)
}
if(ntests>28 || is.null(crit))crit=smmvalv2(dfvec=rep(Inf,nrow(pts)),alpha=alpha)
}}}
if(cov1){
if(!DEEP){
if(alpha==.05)p.crit=0.00676
if(alpha!=.05)p.crit=ancJNmpcp(n1,n2,alpha=alpha,regfun=regfun,nreps=nreps,MC=MC,cov1=cov1)$pc.est
crit=qnorm(1-p.crit/2)
}}
if(DEEP){
if(p==2){
p.crit=.012
if(alpha!=.05)p.crit=ancJNmpcp(n1,n2,alpha=alpha,regfun=regfun,nreps=nreps,MC=MC,cov1=cov1)$pc.est
crit=qnorm(1-p.crit/2)
}
if(p>2){
if(length(pts)>1){
if(SMM){
if(ntests<=28){
if(alpha==.05)crit<-smmcrit(Inf,ntests)
if(alpha==.01)crit<-smmcrit01(Inf,ntests)
}
if(ntests>28 || is.null(crit))crit=smmvalv2(dfvec=rep(Inf,nrow(pts)),alpha=alpha)
}
if(!SMM){
p.crit=ancJNmpcp(n1,n2,alpha=alpha,regfun=regfun,nreps=nreps,MC=MC,cov1=FALSE)
crit=qnorm(1-p.crit/2)
}}}}
mat[,6]=est-crit*sd
mat[,7]=est+crit*sd
flag=rep(FALSE,nrow(mat))
flag.chk1=as.logical(mat[,6]>null.value)
flag.chk2=(mat[,7]<null.value)
flag.chk=(flag.chk1+flag.chk2>0)
num.sig=sum(flag.chk)
if(p==2){
if(plotit){
plot(pts[,1],pts[,2],xlab=xlab,ylab=ylab,type='n')
flag[flag.chk]=TRUE
points(pts[!flag,1],pts[!flag,2],pch='*')
points(pts[flag,1],pts[flag,2],pch='+')  #significant points
}}
output.sig=NULL
if(p==2){
if(num.sig>0){
output.sig=matrix(NA,nrow=num.sig,ncol=8)
output.sig[,1]=pts[flag,1]
output.sig[,2]=pts[flag,2]
output.sig[,3]=mat[flag,1]
output.sig[,4]=mat[flag,2]
output.sig[,5]=mat[flag,3]
output.sig[,6]=mat[flag,6]
output.sig[,7]=mat[flag,7]
output.sig[,8]=mat[flag,8]
dimnames(output.sig)<-list(NULL,c('COV 1','COV 2','Est 1', 'Est 2','DIF','ci.low','ci.hi','p.value'))
if(!ALL){
mat=NULL
pts=NULL
}
if(pr){
if(ALL)print('To get  only the results for all covariate points where this is a significant result, set ALL=FALSE')
}
}}
list(n=nv,num.sig=num.sig,p.crit=p.crit,points=pts,output.sig=output.sig,output=mat)
}


#' @title Compute Critical P-Value for ancJNmp
#'
#' @description
#' Internal function to compute the critical p-value for ancJNmp via simulation.
#' Used to control familywise error rate in multiple comparisons.
#'
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param regfun Regression function. Default is qreg (quantile regression).
#' @param CPP Logical. Use C++ implementation (requires WRScpp package). Default is FALSE.
#' @param nreps Number of simulation replications. Default is 1000.
#' @param alpha Desired familywise error rate. Default is 0.05.
#' @param MC Logical. Use parallel processing via mclapply. Default is FALSE.
#' @param SEED Logical. Set random seed for reproducibility. Default is TRUE.
#' @param cov1 Logical. Use covariate points from first group. Default is FALSE.
#'
#' @return A list with component pc.est (the estimated critical p-value).
#'
#' @keywords internal
# ancJNmpcp
ancJNmpcp<-function(n1,n2,regfun=qreg,CPP=FALSE,nreps=1000,alpha=.05,MC=FALSE,
SEED=TRUE,cov1=FALSE){
if(CPP)library(WRScpp)
if(SEED)set.seed(2)
x=list()
n=max(c(n1,n2))
nmiss=n-min(c(n1,n2))
for(i in 1:nreps){
x[[i]]=rmul(n,p=6)
if(n1<n2){
if(nmiss>0)x[[i]][1:nmiss,1:3]=NA
}
if(n1>n2){
if(nmiss>0)x[[i]][1:nmiss,4:6]=NA
}
}

if(!MC)vals=lapply(x,ancJNmpcp.sub,regfun=regfun,cov1=cov1)
if(MC)vals=mclapply(x,ancJNmpcp.sub,regfun=regfun,cov1=cov1)
vals=as.vector(matl(vals))
pc.est=hd(vals,alpha)
list(pc.est=pc.est)
}


#' @title Helper Function for ancJNmpcp Simulation
#'
#' @description
#' Internal helper function that processes a single simulation replicate for ancJNmpcp.
#' Extracts the minimum p-value from ancJNmp analysis.
#'
#' @param x Simulated data matrix with 6 columns (3 for each group).
#' @param regfun Regression function.
#' @param cov1 Logical. Use covariate points from first group.
#'
#' @return Minimum p-value from the analysis.
#'
#' @keywords internal
# ancJNmpcp.sub
ancJNmpcp.sub<-function(x,regfun=qreg,cov1=FALSE){
pts=NULL
z=elimna(x[,1:3])
z2=elimna(x[,4:6])
if(cov1)pts=z[,1:2]
res1=ancJNmp(z[,1:2],z[,3],z2[,1:2],z2[,3],SEED=TRUE,plotit=FALSE,
regfun=regfun,pts=pts)$output
v=min(res1[,8])
}



#' @title Johnson-Neyman ANCOVA with Multiple Covariates (Improved Version)
#'
#' @description
#' Compare two independent groups using Johnson-Neyman ANCOVA allowing multiple covariates.
#' This is an improved replacement for ancJNmp. Controls familywise error rate across
#' multiple design points.
#'
#' @param x1 Matrix of covariates for group 1 (2-6 columns).
#' @param y1 Outcome values for group 1.
#' @param x2 Matrix of covariates for group 2 (must match x1 in columns).
#' @param y2 Outcome values for group 2.
#' @param regfun Regression function to use. Default is MMreg (MM-regression).
#' @param p.crit Critical p-value for controlling familywise error. If NULL (default), computed automatically.
#' @param DEEP Logical. If TRUE (default), use deepest half of combined covariate points.
#'   If FALSE, use depth-based design point selection.
#' @param plotit Logical. If TRUE (default), plot results for 2 covariates.
#' @param xlab,ylab Axis labels for plots.
#' @param null.value Null hypothesis value for difference. Default is 0.
#' @param WARNS Logical. If FALSE (default), suppress warnings during computation.
#' @param alpha Familywise error rate. Default is 0.05. Also supports 0.01.
#' @param pts Design points matrix (p columns). If NULL, chosen automatically.
#' @param SEED Logical. Set random seed for reproducibility. Default is TRUE.
#' @param nboot Number of bootstrap samples for variance estimation. Default is 100.
#' @param xout Logical. Remove outliers from covariates. Default is FALSE.
#' @param outfun Outlier detection function. Default is outpro.
#' @param ... Additional arguments passed to regression function.
#'
#' @return A list with components:
#'   \item{n1, n2}{Sample sizes for groups 1 and 2}
#'   \item{num.sig}{Number of significant comparisons}
#'   \item{p.crit}{Critical p-value used}
#'   \item{points}{Design points used}
#'   \item{output}{Matrix with estimates, tests, and p-values for each point}
#'   \item{significant.points}{Covariate points with significant differences}
#'
#' @details
#' For DEEP=FALSE, design points are chosen as the deepest point in combined data
#' plus points on the 0.5 depth contour. For DEEP=TRUE, uses deepest half of points.
#'
#' When alpha=0.05 or 0.01 and number of covariates ≤6, p.crit is determined automatically
#' to control familywise error. For >6 covariates, adjustment method is unknown.
#'
#' @export
#' @keywords ancova multivariate
#'
#' @examples
#' \dontrun{
#'   # Two covariates
#'   x1 <- matrix(rnorm(100), 50, 2)
#'   y1 <- x1[,1] + x1[,2] + rnorm(50)
#'   x2 <- matrix(rnorm(100), 50, 2)
#'   y2 <- x2[,1] + x2[,2] + 1 + rnorm(50)
#'   result <- ancJNPVAL(x1, y1, x2, y2)
#' }
# ancJNPVAL
ancJNPVAL<-function(x1,y1,x2,y2,regfun=MMreg,p.crit=NULL,DEEP=TRUE,
plotit=TRUE,xlab='X1',ylab='X2',null.value=0,WARNS=FALSE,
alpha=.05, pts=NULL,SEED=TRUE,nboot=100,xout=FALSE,outfun=outpro,...){
#
# Compare two independent  groups using a generalization of the ancts function that
#  allows more than one covariate.
#
# Design points can be specified via the argument
# pts: a matrix with p=ncol(x1) columns.
#
# DEEP=FALSE: If pts=NULL, design points are chosen to be deepest  point in
# rbind(x1,x2)  plus points on the .5 depth contour.
#
# DEEP=TRUE, choose deepest half of the points in rbind(x1,x2) and use critical p-value indicated by
# p.crit.
#
#  alpha=.05, refers to the desired probability of one or more Type I errors. If
#  p.crit=NULL,
#  when alpha=.05 or .01 and number of covariates is <=6, p.crit is
#  determined quickly by this function. That is, the familywise error will be approximately alpha.
#
# If number of covariates is > 6, unknown how to adjust p.crit to control familywise error.
#
#  plotit=TRUE: if p=2 covariates, plot covariate points with
#  non-significant points indicated by * and  significant points by +

# (This function replaces anctsmp, which does not have an option for
#  using the deepest half of the covariate points.)
#
if(SEED)set.seed(2)
if(!is.null(pts[1]))DEEP=FALSE
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have a different number of columns')
p=ncol(x1)
if(p==1)stop('Should have at least two covariates')
if(p>6)stop('Current version is limited to six covariates or less')
p1=p+1
m1=elimna(cbind(x1,y1))
x1=m1[,1:p]
y1=m1[,p1]
m2=elimna(cbind(x2,y2))
x2=m2[,1:p]
y2=m2[,p1]
#
if(xout){
m<-cbind(x1,y1)
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
if(is.null(pts[1])){
if(!DEEP){
x1<-as.matrix(x1)
pts<-ancdes(unique(rbind(x1,x2)))
p.crit=NULL
}
if(DEEP){
xall=unique(rbind(x1,x2))
pd=pdepth(xall)
id.keep=which(pd>median(pd))
pts=xall[id.keep,]
pts=unique(pts)
}}
pts<-as.matrix(pts)
ntests=nrow(pts)
mat<-matrix(NA,ntests,8)
dimnames(mat)<-list(NULL,c('Est 1', 'Est 2','DIF','TEST','se','ci.low','ci.hi','p.value'))
if(!WARNS)options(warn=-1)
sqsd1=regYvar(x1,y1,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
sqsd2=regYvar(x2,y2,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
#  xout=F because leverage points have already been removed.
est1=regYhat(x1,y1,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
est2=regYhat(x2,y2,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
if(!WARNS)options(warn=0) #turn warnings back on
mat[,1]=est1
mat[,2]=est2
est=est1-est2
mat[,3]=est
sd=sqrt(sqsd1+sqsd2)
mat[,5]=sd
tests=(est1-est2)/sd
mat[,4]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,8]=pv
crit=NULL
if(ntests==1)crit=qnorm(1-alpha/2)
if(nrow(pts)>1){
if(ntests<=25){
if(alpha==.05)crit<-smmcrit(Inf,ntests)
if(alpha==.01)crit<-smmcrit01(Inf,ntests)
}}
if(ntests>25){
pvals.05=c(NA,.00615847,0.002856423,.00196,0.001960793,0.001120947)
pvals.01=c(NA,0.001006744,0.000237099,0.0003169569,0.0002031497,9.442465e-05)
if(alpha==.05){
crit=qnorm(1-pvals.05[p]/2)
p.crit=pvals.05[p]
}
if(alpha==.01){
crit=qnorm(1-pvals.01[p]/2)
p.crit=pvals.01[p]
}
}
mat[,6]=est-crit*sd
mat[,7]=est+crit*sd
flag=rep(FALSE,nrow(mat))
flag.chk1=as.logical(mat[,6]>null.value)
flag.chk2=(mat[,7]<null.value)
flag.chk=(flag.chk1+flag.chk2>0)
#if(!is.null(p.crit))num.sig=sum(mat[,8]<=p.crit)
num.sig=sum(mat[,8]<=p.crit)
if(p==2){
if(plotit){
plot(pts[,1],pts[,2],xlab=xlab,ylab=ylab,type='n')
flag[flag.chk]=TRUE
points(pts[!flag,1],pts[!flag,2],pch='*')
points(pts[flag,1],pts[flag,2],pch='+')
}}
sig.points=NULL
if(!is.null(p.crit)){
if(num.sig>0){
pick=which(mat[,8]<=p.crit)
sig.points=pts[pick,]
}}
list(n1=n1,n2=n2,num.sig=num.sig,p.crit=p.crit,points=pts,output=mat, significant.points=sig.points)
}



#' @title ANCOVA with Adjusted Confidence Intervals for Single Covariate
#'
#' @description
#' For two independent groups, compute confidence intervals for the difference between
#' typical values of Y given X, based on robust regression. Controls familywise error rate
#' across multiple covariate values using improved methods (better power than Hochberg/Hommel).
#'
#' @param x1 Covariate values for group 1 (single covariate only).
#' @param y1 Outcome values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Outcome values for group 2.
#' @param regfun Regression function. Default is tsreg (Theil-Sen estimator).
#' @param pts Covariate values where CIs are computed. If NULL, uses npts evenly-spaced points.
#' @param ALL Logical. If TRUE, return all results; if FALSE, only significant results. Default is FALSE.
#' @param npts Number of design points if pts=NULL. Default is 25.
#' @param plotit Logical. Plot results. Default is TRUE.
#' @param SCAT Logical. Show scatterplot. Default is TRUE.
#' @param pch1,pch2 Plotting characters for groups.
#' @param nboot Number of bootstrap samples for SE estimation. Default is 100.
#' @param ADJ Logical. Adjust critical value for simultaneous coverage. Default is TRUE.
#' @param xout Logical. Remove outliers from covariates. Default is FALSE.
#' @param outfun Outlier detection function. Default is outpro.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param p.crit Critical p-value. Default is 0.015.
#' @param alpha Familywise error rate. Default is 0.05.
#' @param crit Critical value. If NULL, computed automatically.
#' @param null.value Null hypothesis value. Default is 0.
#' @param plotPV Logical. Plot p-values. Default is FALSE.
#' @param scale,span Parameters for LOWESS smoothing of p-values.
#' @param xlab,ylab,ylab2 Axis labels.
#' @param MC Logical. Use parallel processing for critical value computation. Default is FALSE.
#' @param nreps Number of replications for computing critical value. Default is 1000.
#' @param pch Plotting character. Default is '*'.
#' @param ... Additional arguments passed to regression function.
#'
#' @return A list with confidence intervals, test results, and plots.
#'
#' @details
#' This function is identical to regYci2g. It handles larger numbers of covariate values
#' than ancJN. For alpha=0.05, adjustment is quick. For other alpha values or small samples,
#' critical value computation can be slow; use MC=TRUE to speed up.
#'
#' @export
#' @keywords ancova regression
#'
#' @examples
#' \dontrun{
#'   x1 <- rnorm(50); y1 <- 2*x1 + rnorm(50)
#'   x2 <- rnorm(50); y2 <- x2 + 1 + rnorm(50)
#'   result <- anclin(x1, y1, x2, y2, npts=10)
#' }
# anclin
anclin<-function(x1,y1,x2,y2,regfun=tsreg,pts=NULL,ALL=FALSE,npts=25,plotit=TRUE,SCAT=TRUE,
pch1='*',pch2='+',
nboot=100,ADJ=TRUE,xout=FALSE,outfun=outpro,SEED=TRUE,p.crit=.015,
alpha=.05,crit=NULL,null.value=0,plotPV=FALSE,scale=TRUE,span=.75,xlab='X',ylab='p-values',ylab2='Y',MC=FALSE,nreps=1000,pch='*',...){
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
#   A single covariate is assumed.
#   If alpha=.05, an adjustment can be made quickly. Otherwise an
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
if(p>1)stop('Current version allows one covariate only')
p1=p+1
vals=NA
x2<-xy[,1:p]
y2<-xy[,p1]
x2<-as.matrix(x2)
n1=length(y1)
n2=length(y2)
n=min(c(n1,n2))
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
chk.test=abs(est1-est2-null.value)/sd
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
plot(pts,pv,xlab=xlab,ylab=ylab,pch=pch)
}
list(output=est,p.crit=p.crit,crit.value=crit,num.sig=sum(est[,5]<=p.crit))
}


#' @title ANCOVA Quantile Shift Effect Size
#'
#' @description
#' Compute quantile shift (QS) effect sizes for two independent groups (control vs experimental).
#' The QS effect size is the conditional quantile of the control group corresponding to the
#' median of Y for the experimental group, given X.
#'
#' @param x1 Covariate values for control group (single covariate).
#' @param y1 Outcome values for control group.
#' @param x2 Covariate values for experimental group.
#' @param y2 Outcome values for experimental group.
#' @param pts Covariate points where effect sizes are computed. If NULL, chosen automatically.
#' @param xout Logical. Remove outliers from covariates. Default is FALSE.
#' @param ALL Logical. If TRUE, use all unique covariate values. Default is FALSE (use 3 points).
#' @param npts Number of points if pts=NULL and ALL=TRUE. Default is 10.
#' @param outfun Outlier detection function. Default is outpro.
#' @param REQMIN Minimum required value for qinvreg. Default is 0.001.
#' @param ... Additional arguments.
#'
#' @return Matrix with columns: Pts (covariate values), Y.hat4ExpGrp (predicted median for experimental group),
#'   QS.Effect.Size (quantile shift effect size).
#'
#' @export
#' @keywords ancova effect-size
# anclin.QS
anclin.QS<-function(x1,y1,x2,y2,pts=NULL,xout=FALSE,ALL=FALSE,npts=10,outfun=outpro,REQMIN=.001,...){
#
#  x1, y1 is the control group
#  x2  y2 is the experimental group
#
#  For Exp group, estimate the median of Y given the x values stored in
#  pts
#  pts=NULL: If ALL=TRUE, 20 points are chosen  by this function
#         otherwise three points are used.
#
#  The QS effect size is the conditional quantile of the control group corresponding
#  to the median of Y, given x, for the experimental group.
#
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
if(p>1)stop('Current version allows one covariate only')
p1=p+1
vals=NA
x2<-xy[,1:p]
y2<-xy[,p1]
x2<-as.matrix(x2)
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
n1=length(y1)
n2=length(y2)
n=min(c(n1,n2))
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
if(!is.null(pts))npts=length(pts)
if(is.null(pts)){
xall=unique(c(x1,x2))
if(ALL)pts=xall
if(!ALL){
L1=qest(x1,.2)
L2=qest(x2,.2)
U1=qest(x1,.8)
U2=qest(x2,.8)
L=max(L1,L2)
U=min(U1,U2)
if(ALL)pts=seq(L,U,length.out=npts)
else{pts=c(L,(L+U)/2,U)
npts=3
}
}}
e=reg.pred(x2,y2,xr=pts,regfun=Qreg,q=.5,xout=FALSE)
qs=NA
for(i in 1:npts){
qs[i]=qinvreg(x1,y1,pts[i],e[i],REQMIN=REQMIN)
}
M=cbind(pts,e,qs)
dimnames(M)=list(NULL,c('Pts','Y.hat4ExpGrp','QS.Effect.Size'))
M
}


#' @title Bootstrap Confidence Interval for Quantile Shift Effect Size
#'
#' @description
#' Compute bootstrap confidence intervals for the conditional quantile shift (QS)
#' effect size in ANCOVA for two independent groups.
#'
#' @param x1 Covariate values for control group (single covariate).
#' @param y1 Outcome values for control group.
#' @param x2 Covariate values for experimental group.
#' @param y2 Outcome values for experimental group.
#' @param alpha Alpha level for confidence intervals. Default is 0.05.
#' @param pts Covariate points where CIs are computed. If NULL, chosen automatically.
#' @param xout Logical. Remove outliers from covariates. Default is FALSE.
#' @param ALL Logical. If TRUE, use many points. Default is FALSE (use 3 points for speed).
#' @param npts Number of points if pts=NULL and ALL=TRUE. Default is 10.
#' @param outfun Outlier detection function. Default is outpro.
#' @param nboot Number of bootstrap samples. Default is 200.
#' @param MC Logical. Use parallel processing via mclapply. Default is TRUE.
#' @param REQMIN Minimum required value for qinvreg. Default is 0.01.
#' @param SEED Logical. Set random seed for reproducibility. Default is TRUE.
#' @param ... Additional arguments.
#'
#' @return Matrix with confidence intervals for QS effect sizes at each covariate point.
#'
#' @export
#' @keywords ancova effect-size bootstrap
# anclin.QS.CIpb
anclin.QS.CIpb<-function(x1,y1,x2,y2,alpha=.05,pts=NULL,xout=FALSE,ALL=FALSE,npts=10,outfun=outpro,nboot=200,
MC=TRUE,REQMIN=.01,SEED=TRUE,...){
#
#  ANCOVA
#
#   Compute a confidence interval for the conditional quantile shift
#   measure of effect size
#  x1, y1 is the control group
#  x2  y2 is the experimental group
#
#  for Exp group, estimate the median of Y given x for values stored in
#  pts
#  pts=NULL: If ALL=TRUE, 20 points are chosen  by this function
#         otherwise three points are used to reduce execution time.
#
#  The QS effect size is the conditional quantile of the control group corresponding
#  to the median of Y for the experimental group.
#
#
if(SEED)set.seed(2)
xy=elimna(cbind(x1,y1))
x1=as.matrix(x1)
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
if(p>1)stop('Current version allows one covariate only')
p1=p+1
vals=NA
x2<-xy[,1:p]
y2<-xy[,p1]
x2<-as.matrix(x2)
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
n1=length(y1)
n2=length(y2)
n=min(c(n1,n2))
if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE,...)$keep
m<-m[flag,]
n1=nrow(m)
x1<-m[,1:p]
y1<-m[,p1]
x1=as.matrix(x1)
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE,...)$keep
m<-m[flag,]
n2=nrow(m)
n=min(c(n1,n2))
x2<-m[,1:p]
y2<-m[,p1]
x2=as.matrix(x2)
}
if(!is.null(pts))npts=length(pts)
if(is.null(pts)){
xall=unique(c(x1,x2))
if(ALL)pts=xall
if(!ALL){
L1=qest(x1,.2)
L2=qest(x2,.2)
U1=qest(x1,.8)
U2=qest(x2,.8)
L=max(L1,L2)
U=min(U1,U2)
if(ALL)pts=seq(L,U,length.out=npts)
else{pts=c(L,(L+U)/2,U)
npts=3
}
}}
v=NA
m=list()
for(i in 1:nboot){
id1=sample(n1,replace=TRUE)
id2=sample(n2,replace=TRUE)
m[[i]]=list(x1[id1],y1[id1],x2[id2],y2[id2])
}
if(!MC)v=lapply(m,anclinQS.sub,pts=pts,npts=npts,...)
if(MC){
v=mclapply(m,anclinQS.sub,pts=pts,npts=npts,...)
}
v=matl(v)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
nv=nrow(v)
CI=matrix(NA,nv,2)
pv=NA
for(i in 1:nv){
pv[i]=mean(v[i,]<0.5)
pv[i]=2*min(pv[i],1-pv[i])
sv=sort(v[i,])
CI[i,1]=sv[ilow]
CI[i,2]=sv[ihi]
}
output=matrix(NA,nrow=nv,ncol=4)
output[,1]=pts
output[,2]=pv
output[,3:4]=CI
e=anclin.QS(x1,y1,x2,y2,pts=pts)
e=as.matrix(e[,2:3])
if(nv==1)e=t(e)
output=cbind(output,e)
dimnames(output)=list(NULL,c('X','p.value','ci.low','ci.hi','Median.ExpGrp','QS.effect'))
output
}

#' @title Plot Quantile Shift Effect Sizes
#'
#' @description
#' Plot quantile shift (QS) effect sizes as a function of covariate values for two independent groups.
#' The QS effect size is the conditional quantile of the control group corresponding to the median
#' of Y for the experimental group, given X.
#'
#' @param x1 Covariate values for control group (single covariate).
#' @param y1 Outcome values for control group.
#' @param x2 Covariate values for experimental group.
#' @param y2 Outcome values for experimental group.
#' @param pts Covariate points where effect sizes are plotted. If NULL, chosen automatically.
#' @param q Quantile for determining covariate range. Default is 0.1.
#' @param xout Logical. Remove outliers from covariates. Default is FALSE.
#' @param ALL Logical. If TRUE, use many points. Default is TRUE.
#' @param npts Number of points if pts=NULL. Default is 10.
#' @param line Logical. Connect points with lines. Default is TRUE.
#' @param xlab,ylab Axis labels.
#' @param outfun Outlier detection function. Default is outpro.
#' @param REQMIN Minimum required value for qinvreg. Default is 0.001.
#' @param ... Additional arguments.
#'
#' @return Matrix with columns: Pts, Y.hat4ExpGrp, QS.Effect.Size (returned invisibly).
#'
#' @export
#' @keywords ancova plot effect-size
# anclinQS.plot
anclinQS.plot<-function(x1,y1,x2,y2,pts=NULL,q=0.1,xout=FALSE,ALL=TRUE,npts=10,line=TRUE,
xlab='X',ylab='QS.Effect',outfun=outpro,REQMIN=.001,...){
#
#  x1, y1 is the control group
#  x2  y2 is the experimental group
#
#  For Exp group, estimate the median of Y given the x values stored in
#  pts
#  pts=NULL: If ALL=TRUE, 20 points are chosen  by this function
#         otherwise three points are used.
#
#  The QS effect size is the conditional quantile of the control group corresponding
#  to the median of Y, given x, for the experimental group.
#   The function plots estimates of the QS effect size for the points in pts
#
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
if(p>1)stop('Current version allows one covariate only')
p1=p+1
vals=NA
x2<-xy[,1:p]
y2<-xy[,p1]
x2<-as.matrix(x2)
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
n1=length(y1)
n2=length(y2)
n=min(c(n1,n2))
if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE,...)$keep
m<-m[flag,]
n1=nrow(m)
x1<-m[,1:p]
y1<-m[,p1]
x1=as.matrix(x1)
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE,...)$keep
m<-m[flag,]
n2=nrow(m)
n=min(c(n1,n2))
x2<-m[,1:p]
y2<-m[,p1]
x2=as.matrix(x2)
}
if(!is.null(pts))npts=length(pts)
if(is.null(pts)){
if(q<=0 || q>=1)stop('Argument q must be greater than 0 and less than 1')
qu=1-q
L1=qest(x1,q)
L2=qest(x2,q)
U1=qest(x1,qu)
U2=qest(x2,qu)
L=max(L1,L2)
U=min(U1,U2)
if(ALL)pts=seq(L,U,length.out=npts)
else{pts=c(L,(L+U)/2,U)
npts=3
}}
e=reg.pred(x2,y2,xr=pts,regfun=Qreg,q=.5,xout=FALSE)
qs=NA
for(i in 1:npts){
qs[i]=qinvreg(x1,y1,pts[i],e[i],REQMIN=REQMIN)
}
M=cbind(pts,e,qs)
if(line){
plot(pts,qs,xlab=xlab,ylab=ylab,ylim=c(0,1),type='n')
lines(pts,qs)
}
else
plot(pts,qs,xlab=xlab,ylab=ylab,ylim=c(0,1))
dimnames(M)=list(NULL,c('Pts','Y.hat4ExpGrp','QS.Effect.Size'))
M
}


#' @title Helper for anclinQS.plot Bootstrap
#'
#' @description
#' Internal helper function for bootstrap computation in anclinQS.plot.
#'
#' @param m List with bootstrap sample data.
#' @param pts Design points.
#' @param npts Number of points.
#' @param ... Additional arguments.
#'
#' @return Vector of QS effect sizes.
#'
#' @keywords internal
# anclinQS.sub
anclinQS.sub<-function(m,pts,npts=npts,...){
v=anclin.QS(m[[1]],m[[2]],m[[3]],m[[4]],pts=pts,npts=npts,...)[,3]
v
}


#' @title ANCOVA for Logistic Regression (Binary Outcomes)
#'
#' @description
#' Wrapper function for ANCOVA with binary (logistic) outcomes using ancovampG with twobinom test.
#'
#' @inheritParams ancovampG
#'
#' @return List with components: points (design points), results (test results), num.sig (number of significant tests).
#'
#' @export
#' @keywords ancova logistic
# anclog
anclog<-function(x1,y1,x2,y2,fr1=1,fr2=1,
alpha=.05, pts=NULL,SEED=TRUE,DH=FALSE,FRAC=.5,cov.fun=skip.cov,
pr=FALSE,q=.5,plotit=FALSE,pv=FALSE,theta=50,xlab=' ',ylab=' ',SCAT=FALSE,zlab=' ',...){
res=ancovampG(x1=x1,y1=y1,x2=x2,y2=y2,fr1=fr1,fr2=fr2, tr=.2,
alpha=alpha, pts=pts,SEED=SEED,test=twobinom,DH=DH,FRAC=FRAC,cov.fun=cov.fun,
pr=pr,q=q,plotit=plotit,pv=pv,theta=theta,xlab=xlab,ylab=ylab,SCAT=SCAT,zlab=zlab,...)
list(points=res$points,results=res$results,num.sig=res$num.sig)
}


#' @title ANCOVA Effect Sizes for Multiple Covariates
#'
#' @description
#' Estimate effect sizes at specified covariate points for two or more covariates.
#' For each point, identifies nearby observations and computes effect size measures.
#'
#' @param x1 Matrix of covariates for group 1 (2 or more columns).
#' @param y1 Outcome values for group 1.
#' @param x2 Matrix of covariates for group 2.
#' @param y2 Outcome values for group 2.
#' @param fr1 Span for group 1 running interval smoother. Default is 1.
#' @param fr2 Span for group 2 running interval smoother. Default is 1.
#' @param tr Trim proportion. Default is 0.2.
#' @param pts Matrix of covariate points where effect sizes are computed. If NULL, uses significant points from ancdetM4.
#' @param xout Logical. Remove outliers from covariates. Default is FALSE.
#' @param outfun Outlier detection function. Default is outpro.
#' @param ... Additional arguments.
#'
#' @return List with effect size estimates at each design point.
#'
#' @details
#' Minimum sample size must be at least 50 for each group.
#'
#' @export
#' @keywords ancova effect-size multivariate
# ancM.COV.ES
ancM.COV.ES<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,pts=NULL,xout=FALSE,outfun=outpro,...){
#
#
#  For two or more covariates, estimate effect sizes for
#  a collection of points.
#
#  That is, for each point of interest, determine
#  a cloud of points close to it and based on the
#  corresponding y values, compute measures of effect size
#
#  If pts=NULL
#  the significant points returned by
#  ancdetM4 are used
#
p=ncol(x1)
if(p<2)stop('This function is for two or more covariates')
p1=p+1
if(ncol(x2)!=p)stop('x1 and x2 do not have the same number of columns')
xy=elimna(cbind(x1,y1))
x1=xy[,1:p]
y1=xy[,p1]
xy=elimna(cbind(x2,y2))
x2=xy[,1:p]
y2=xy[,p1]
if(min(length(y1),length(y2))<50)stop('The minimum sample size must be greater than or equal to  50')
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag,]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag,]
y2<-y2[flag]
}
guide1=NA
guide2=NA
x1=as.matrix(x1)
x2=as.matrix(x2)
a=ancdetM4(x1=x1,y1=y1,x2=x2,y2=y2,fr1=fr1,fr2=fr2,tr=tr,pts=pts,...)
if(is.null(pts))pts=a$significant.points
M=NA
SML.class=NA
if(!is.null(pts)){
if(is.vector(pts))pts=matrix(pts,nrow=1)
nr=nrow(pts)
M=matrix(NA,nr,6)
SML.class=matrix(0,nr,6)
m1<-covmcd(x1)
m2<-covmcd(x2)
for(i in 1:nr){
id1=near3d(x1,pts[i,],fr1,m1)
id2=near3d(x2,pts[i,],fr2,m2)
if(sum(id1)<10 ||sum(id2)<10)print(paste('For point',j,'not enough nearest neighbors'))
if(sum(id1)>=10 & sum(id2)>10){
ES=ES.summary(y1[id1],y2[id2])
M[i,]=ES[,1]
}}
dimnames(M)=list(NULL,names(ES[,1]))
dum1=rnorm(50)
dum2=rnorm(50)
guide1=ES.summary(dum1,dum2+3)[,-1]
guide2=ES.summary(dum1,dum2-3)[,-1]
}
if(!is.na(M[1])){
for(i in 1:nr){
flag=M[i,] <= guide1[,2] || M[i,]>= guide2[,2]
SML.class[i,flag]=1
flag=M[i,] <= guide1[,3] || M[i,]>= guide2[,3]
SML.class[i,flag]=2
flag=M[i,] <= guide1[,4] || M[i,]>= guide2[,4]
SML.class[i,flag]=3
}}
leg1='0=At most small, 1=between S and M, 2=between M and L, 3=greater than L'
list(ES.REL.MAG.G1.less.than.G2=guide1, ES.REL.MAG.G1.greater.than.G2=guide2,Est=M,legend.4.SML.class=leg1,SML.class=SML.class)
}


#' @title ANCOVA for Multiple Groups with Multiple Covariates
#'
#' @description
#' ANCOVA for J groups (J >= 2) with two or more covariates, using trimmed means or medians.
#' Provides omnibus tests or pairwise multiple comparisons with various estimation methods.
#'
#' @param x Covariate data. Either a matrix with Jp columns (p covariates per group) or a list of
#'   J matrices each with p columns.
#' @param y Outcome data. Either a matrix with J columns or a list of length J.
#' @param pool Logical. If TRUE (default), pool data when determining design points and scatter.
#'   If FALSE, use data from group jcen.
#' @param jcen Group index to use as center when pool=FALSE. Default is 1.
#' @param fr Span for running interval smoother. Default is 1.
#' @param depfun Function to compute depth of points. Default is fdepth (projection depth).
#' @param nmin Minimum sample size required for each group at design points. Default is 8.
#' @param op Operation mode: 1=omnibus test with trimmed means, 2=omnibus test with medians,
#'   3=multiple comparisons with trimmed means (recommended), 4=multiple comparisons with medians.
#'   Default is 3.
#' @param tr Trim proportion for trimmed means. Default is 0.2.
#' @param pts Matrix of design points (p columns). If NULL, chosen automatically.
#' @param SEED Logical. Set random seed for reproducibility. Default is TRUE.
#' @param pr Logical. Print informational messages. Default is TRUE.
#' @param cop Parameter for ancdes function. Default is 3.
#' @param con Contrast specification. Default is 0.
#' @param nboot Number of bootstrap samples for op=3 or op=4. Default is NA.
#' @param alpha Alpha level for confidence intervals. Default is 0.05.
#' @param bhop Logical. Use Benjamini-Hochberg procedure. Default is FALSE.
#'
#' @return A list with test results at each design point and a matrix of sample sizes.
#'   For op=1 or 2, returns matrix with columns: point, test.stat, p-value.
#'   For op=3 or 4, returns list with results for each design point.
#'
#' @details
#' Requires at least two covariates (for single covariate, use ancmg1).
#' Design points with sample size < nmin are eliminated.
#' op=2 not recommended for discrete data with ties; use op=4 instead.
#'
#' @export
#' @keywords ancova multivariate
#'
#' @examples
#' \dontrun{
#'   # Three groups, two covariates
#'   x1 <- cbind(rnorm(30), rnorm(30))
#'   x2 <- cbind(rnorm(30), rnorm(30))
#'   x3 <- cbind(rnorm(30), rnorm(30))
#'   x <- list(x1, x2, x3)
#'   y <- list(rnorm(30), rnorm(30, 0.5), rnorm(30, 1))
#'   result <- ancmg(x, y, op=3)
#' }
# ancmg
ancmg<-function(x,y,pool=TRUE,jcen=1,fr=1,depfun=fdepth,nmin=8,op=3,tr=.2,pts=NULL,
SEED=TRUE,pr=TRUE,cop=3,con=0,nboot=NA,alpha=.05,bhop=FALSE){
#
# ANCOVA
# for two or more groups based on trimmed means or medians
# Two or more  covariates is assumed.
#
# op=1 use omnibus test for trimmed means, with trimming given by tr
# op=2 use omnibus test for medians.
#       (Not recommended when there are tied values, use op=4)
# op=3 multiple comparisons using trimming and percentile bootstrap.
#     This method seems best for general use.
# op=4 multiple comparisons using medians and percentile bootstrap
#
# y is matrix with J columns, so have J groups.
# or y can have list mode with length J
#
# x is a matrix with Jp columns, so first p columns
# correspond to the p covariates in the first group, etc.
# Or,
# x can have list mode with length J and each component
# being a matrix with p columns.
# So if covariates for group 1 are in the matrix m1
# x[[1]]<-m1 will store them in x, x having list mode
#
# nmin is the minimum sample size allowed for any group
# when testing hypotheses.
# If a design point results in a sample size <nmin,
# that point is eliminated.
#
#  pool=T means pool the data when determining the center of the
#  design points and the measure of scatter when applying the smooth
#
#  pool=F, does not pool but rather use the data from group
#  jcen to determine the center and the measure of scatter
#
#  pts, a matrix that can be used to specify design points to be used
#       number of columns should equal number of covariates.
#
#  depfun determines how depth of a point is determined,
#  default is projection depth
#
#  The output includes a matrix of sample sizes. The ith row
#  corresponds to the ith point used to compare groups.
#  The jth column indicates the number of points (the sample size)
#  that was found for the jth group. That is, how many points
#  in the jth group were found that are close to the design point
#  under consideration.
#
#  The command
#   see=ancmg
#  results in see[[1]] containing the results for the first point, see[[2]] the results for the 2nd point etc.
#
output<-NULL
if(SEED)set.seed(2) # set the seed so that MVE always gives same result
if(pr){
if(op==1)print("Trimmed means are to be compared. For medians, use op=2")
if(op==2)print("Medians are to be compared. For trimmed means, use op=1")
if(op==3)print("trimmed means are compared. For medians, use op=4")
if(op==4)print("medians are compared. For trimmed means, use op=3")
}
nval<-NA
if(is.matrix(y))J<-ncol(y)
if(is.list(y))J<-length(y)
if(is.list(x))pval<-ncol(x[[1]])
if(J==1)stop("Only have one group stored in y")
if(pval==1)stop("For one covariate only, using ancgm1")
if(is.matrix(x)){
if(ncol(x)%%J!=0)stop("Number of columns of x should be a multiple of ncol(y)")
}
if(is.matrix(x)){
pval<-ncol(x)/J
if(pval==1)stop("For one covariate only, using ancgm1")
temp<-seq(1,ncol(x),pval)
js<-temp[jcen]
jcenp<-js+pval-1
if(jcenp > ncol(x))stop("jcen has an invalid value")
xcen<-x[,js:jcenp]
}
if(is.list(x))xcen<-x[[jcen]]
if(pool){
if(is.matrix(x))xval<-stackit(x,pval)
if(is.list(x))xval<-stacklist(x)
mval<-cov.mve(xval)
if(is.null(pts))pts<-ancdes(xval,depfun=depfun,cop=cop)
}
if(!pool){
if(is.null(pts))pts<-ancdes(xcen,depfun=depfun,cop=cop)
mval<-cov.mve(xcen)
}
npts=1
if(is.matrix(pts))npts=nrow(pts)
nval<-matrix(NA,ncol=J,nrow=npts)
icl<-0-pval+1
icu<-0
for(j in 1:J){
icl<-icl+pval
icu<-icu+pval
for(i in 1:nrow(pts)){
if(is.matrix(x)  && is.matrix(y)){
nval[i,j]<-length(y[near3d(x[,icl:icu],pts[i,],fr,mval),j])
}
if(is.matrix(x)  && is.list(y)){
tempy<-y[[j]]
nval[i,j]<-length(tempy[near3d(x[,icl:icu],pts[i,],fr,mval)])
}
if(is.list(x)  && is.matrix(y)){
xm<-as.matrix(x[[j]])
nval[i,j]<-length(y[near3d(xm,pts[i,],fr,mval),j])
}
if(is.list(x)  && is.list(y)){
tempy<-y[[j]]
xm<-as.matrix(x[[j]])
nval[i,j]<-length(tempy[near3d(xm,pts[i,],fr,mval)])
}
#
}}
flag<-rep(TRUE,nrow(pts))
for(i in 1:npts){
if(min(nval[i,])<nmin)flag[i]<-FALSE
}
nflag<-FALSE
if(sum(flag)==0){
print("Warning: No design points found with large enough sample size")
nflag<-TRUE
}
flag=as.logical(flag)
if(!nflag){
pts<-pts[flag,] # eliminate points for which the sample size is too small
nval<-nval[flag,]
if(!is.matrix(pts))pts<-t(as.matrix(pts))
output<-matrix(NA,nrow=nrow(pts),ncol=3)
dimnames(output)<-list(NULL,c("point","test.stat","p-value"))
if(op==3 || op==4)output<-list()
}
for(i in 1:nrow(pts)){
if(op==1 || op==2)output[i,1]<-i
icl<-0-pval+1
icu<-0
yval<-list()
for(j in 1:J){
icl<-icl+pval
icu<-icu+pval
if(is.matrix(x)  && is.matrix(y)){
yval[[j]]<-y[near3d(x[,icl:icu],pts[i,],fr,mval),j]
}
if(is.matrix(x)  && is.list(y)){
tempy<-y[[j]]
yval[[j]]<-tempy[near3d(x[,icl:icu],pts[i,],fr,mval)]
}
if(is.list(x)  && is.matrix(y)){
yval[[j]]<-y[near3d(x[[j]],pts[i,],fr,mval),j]
}
if(is.list(x)  && is.list(y)){
tempy<-y[[j]]
yval[[j]]<-tempy[near3d(x[[j]],pts[i,],fr,mval)]
}
#
}
if(op==1)temp<-t1way(yval,tr=tr)
if(op==2)
{
print("WARNING: NOT RECOMMENDED FOR DISCRETE DATA WITH TIES")
print("RECOMMENDATION: Set the argument op=4")
temp<-med1way(yval,SEED=SEED,pr=FALSE)
}
if(op==1 || op==2){
conout=NULL
output[i,3]<-temp$p.value
output[i,2]<-temp$TEST
}
if(op==3){
output[[i]]=linconpb(yval,alpha=alpha,SEED=SEED,con=con,bhop=bhop,tr=tr,nboot=nboot)
}
if(op==4){output[[i]]<-medpb(yval,alpha=alpha,SEED=SEED,con=con,bhop=bhop,nboot=nboot)
}
}
if(nflag)output<-NULL
if(op==1 || op==2)tempout=output
if(op==3 || op==4){
tempout=list()
conout=list()
for(j in 1:length(output))tempout[[j]]=output[[j]]$output
for(j in 1:length(output))conout[[j]]=output[[j]]$con
}
#list(points.chosen=pts,sample.sizes=nval,output=tempout,contrast.coef=conout)
list(points.chosen=pts,sample.sizes=nval,point=tempout,contrast.coef=conout[[1]])
}


#' @title ANCOVA for Multiple Groups with Single Covariate
#'
#' @description
#' ANCOVA for J groups (J >= 2) with a single covariate, using trimmed means or medians.
#' For two or more covariates, use ancmg instead. For two groups, also consider ancova and ancGpar.
#'
#' @param x Covariate data. Matrix with J columns (one per group) or list of length J.
#' @param y Outcome data. Matrix with J columns or list of length J.
#' @param pool Logical. If TRUE (default), pool data when determining design points.
#' @param jcen Group index to use as center when pool=FALSE. Default is 1.
#' @param fr Span for running interval smoother. Default is 1.
#' @param depfun Function to compute depth (not used for single covariate). Default is fdepth.
#' @param nmin Minimum sample size required at design points. Default is 8.
#' @param op Operation mode: 1=omnibus trimmed means, 2=omnibus medians,
#'   3=multiple comparisons trimmed means (recommended), 4=multiple comparisons medians. Default is 3.
#' @param tr Trim proportion. Default is 0.2.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param pr Logical. Print messages. Default is TRUE.
#' @param pts Design points. If NA (default), three points are chosen automatically.
#' @param con Contrast specification. Default is 0.
#' @param nboot Number of bootstrap samples. Default is NA.
#' @param alpha Alpha level. Default is 0.05.
#' @param bhop Logical. Use Benjamini-Hochberg. Default is FALSE.
#'
#' @return List with components: points.chosen, sample.sizes, point (test results), contrast.coef.
#'
#' @details
#' Single covariate version of ancmg. op=2 not recommended for tied values; use op=4.
#'
#' @export
#' @keywords ancova
#'
#' @examples
#' \dontrun{
#'   x <- cbind(rnorm(30), rnorm(30), rnorm(30))
#'   y <- cbind(rnorm(30), rnorm(30, 0.5), rnorm(30, 1))
#'   result <- ancmg1(x, y)
#' }
# ancmg1
ancmg1<-function(x,y,pool=TRUE,jcen=1,fr=1,depfun=fdepth,nmin=8,op=3,tr=.2,
SEED=TRUE,pr=TRUE,pts=NA,con=0,nboot=NA,alpha=.05,bhop=FALSE){
#
# ANCOVA
# for two or more groups based on trimmed means or medians
# Single  covariate is assumed.
#
# FOR TWO OR MORE COVARIATES, USE ANCMG
#
# for two groups and one covariate, also consider ancova and ancGpar
#
# op=1 use omnibus test for trimmed means, with trimming given by tr
# op=2 use omnibus test for medians.
#       (Not recommended when there are tied values, use op=4)
# op=3 multiple comparisons using trimming and percentile bootstrap.
#     This method seems best for general use.
# op=4 multiple comparisons using medians and percentile bootstrap
#
# y is matrix with J columns, so have J groups.
# or y can have list mode with length J
#
# x is a matrix with Jp columns, so first p columns
# correspond to the p covariates in the first group, etc.
#
# Or,
# x can have list mode with length J
#
# nmin is the minimum sample size allowed for any group
# when testing hypotheses.
# If a design point results in a sample size <nmin,
# that point is eliminated.
#
#  pool=T means pool the data when determining the center of the
#  design points and the measure of scatter when applying the smooth
#
#  pool=F, does not pool but rather use the data from group
#  jcen to determine center and the measure of scatter
#
#  The output includes a matrix of sample sizes. The ith row
#  corresponds to the ith point used to compare groups.
#  The jth column indicates the number of points (the sample size)
#  that was found for the jth group. That is, how many points
#  in the jth group were found that are close to the design point
#  under consideration.
#
if(SEED)set.seed(2) # set the seed so that MVE always gives same result
if(pr){
if(op==1)print("Trimmed means are to be compared. For medians, use op=2")
if(op==2)print("Medians are to be compared. For trimmed means, use op=1")
if(op==3)print("20% trimmed means are compared. For medians, use op=4")
if(op==4)print("medians are compared. For 20% trimmed means, use op=3")
}
temp=elimna2g(x,y)
x=temp$x
y=temp$y
output<-NULL
conout=NULL
nval<-NA
if(is.data.frame(x))x=as.matrix(x)
if(is.data.frame(y))y=as.matrix(y)
if(is.matrix(y))J<-ncol(y)
if(is.matrix(x))Jx=ncol(x)
if(is.list(y))J<-length(y)
if(is.list(x))Jx<-length(x)
if(J!=Jx)stop("Only one covariate allowed. Number of covariates not equal to number of groups")
if(J==1)stop("Only have one group stored in y")
if(is.matrix(x)){
if(ncol(x)%%J!=0)stop("Number of columns of x should be a multiple of ncol(y)")
}
if(is.matrix(x)){
xcen<-x[,jcen]
}
if(is.list(x))xcen<-x[[jcen]]
if(is.na(pts[1])){
if(pool){
if(is.matrix(x))xval<-stackit(x,1)
if(is.list(x))xval<-stacklist(x)
temp<-idealf(xval)
pts<-temp$ql
pts[2]<-median(xval)
pts[3]<-temp$qu
}
if(!pool){
temp<-idealf(xcen)
pts<-temp$ql
pts[2]<-median(xval)
pts[3]<-temp$qu
}}
nval<-matrix(NA,ncol=J,nrow=length(pts))
for(j in 1:J){
for(i in 1:length(pts)){
if(is.matrix(x)  && is.matrix(y)){
nval[i,j]<-length(y[near(x[,j],pts[i],fr=fr),j])
}
if(is.matrix(x)  && is.list(y)){
tempy<-y[[j]]
nval[i,j]<-length(tempy[near(x[,j],pts[i],fr=fr)])
}
if(is.list(x)  && is.matrix(y)){
xm<-as.matrix(x[[j]])
xm=as.vector(xm)
nval[i,j]<-length(y[near(xm,pts[i],fr=fr),j])
}
if(is.list(x)  && is.list(y)){
tempy<-y[[j]]
xm<-as.matrix(x[[j]])
xm=as.vector(xm)
nval[i,j]<-length(tempy[near(xm,pts[i],fr=fr)])
}
#
}}
flag<-rep(TRUE,length(pts))
for(i in 1:length(pts)){
if(min(nval[i,])<nmin)flag[i]<-FALSE
}
nflag<-FALSE
if(sum(flag)==0){
print("Warning: No design points found with large enough sample size")
nflag<-TRUE
}
if(!nflag){
pts<-pts[flag] # eliminate points for which the sample size is too small
nval<-nval[flag,]
if(!is.matrix(pts))pts<-t(as.matrix(pts))
output<-matrix(NA,nrow=length(pts),ncol=3)
dimnames(output)<-list(NULL,c("point","test.stat","p-value"))
if(op==3 || op==4)output<-list()
}
for(i in 1:length(pts)){
if(op==1 || op==2)output[i,1]<-i
yval<-list()
for(j in 1:J){
if(is.matrix(x)  && is.matrix(y)){
yval[[j]]<-y[near(x[,j],pts[i],fr=fr),j]
}
if(is.matrix(x)  && is.list(y)){
tempy<-y[[j]]
yval[[j]]<-tempy[near(x[,j],pts[i],fr=fr)]
}
if(is.list(x)  && is.matrix(y)){
yval[[j]]<-y[near3d(x[[j]],pts[i],fr=fr),j]
}
if(is.list(x)  && is.list(y)){
tempy<-y[[j]]
yval[[j]]<-tempy[near(x[[j]],pts[i],fr=fr)]
}
#
}
if(op==1)temp<-t1way(yval,tr=tr)
if(op==2)temp<-med1way(yval,SEED=SEED,pr=FALSE)

if(op==1 || op==2){
output[i,2]<-temp$TEST
output[i,3]<-temp$p.value
}
if(op==3){
output[[i]]<-linconpb(yval,alpha=alpha,SEED=SEED,con=con,bhop=bhop,est=tmean,nboot=nboot)
}
if(op==4){output[[i]]<-medpb(yval,alpha=alpha,SEED=SEED,con=con,bhop=bhop,nboot=nboot)
}
}
if(op==1 || op==2)tempout=output
if(nflag)output<-NULL
if(op==3 || op==4){
conout=list()
tempout=list()
for(j in 1:length(output))tempout[[j]]=output[[j]]$output
for(j in 1:length(output))conout[[j]]=output[[j]]$con
}
list(points.chosen=pts,sample.sizes=nval,point=tempout,contrast.coef=conout[[1]])
}


#' @title Power Estimation for ancmg1
#'
#' @description
#' Estimate statistical power for ancmg1 through simulation, assuming standard normal
#' distributions with first group having shifted mean.
#'
#' @param n Vector of sample sizes (length indicates number of groups).
#' @param del Mean shift for first group. Default is 0.2.
#' @param alpha Significance level. Default is 0.05.
#' @param iter Number of simulation iterations. Default is 100.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param ADJ Logical. Use adjusted p-values. Default is FALSE.
#'
#' @return Estimated power (proportion of rejections).
#'
#' @export
#' @keywords ancova power
#'
#' @examples
#' \dontrun{
#'   power <- ancmg1.power(n=c(30,30,30), del=0.3, iter=500)
#' }
# ancmg1.power
ancmg1.power<-function(n,del=.2,alpha=.05,iter=100,SEED=TRUE,ADJ=FALSE){
#
# n sample sizes, length of n indicates number of groups
#. Estimate power with no data
# Simulate assuming standard normal distributions but first group has a mean del
#
#
J=length(n)
x=list()
y=list()
chk=0
if(SEED)set.seed(2)
for(i in 1:iter){
for(j in 1:J){
x[[j]]=rnorm(n[j])
y[[j]]=rnorm(n[j])
}
y[[1]]=y[[1]]+del
a=ancmg1(x,y,pr=FALSE)
pv=NA
K=length(a$points)
for(k in 1:K){
if(!ADJ)pv[k]=min(a$point[[k]][,3])
else
pv[k]=min(a$point[[k]][,7])
}
if(min(pv)<=alpha)chk=chk+1
}
chk/iter
}




#' @title ANCOVA Multiple Comparisons (Deprecated)
#'
#' @description
#' **DEPRECATED: This function has been eliminated. Use ancmppb instead.**
#'
#' @param x1,y1,x2,y2,fr1,fr2,alpha,pts,est,nboot,bhop,SEED,... Parameters (see ancmppb).
#'
#' @return Error message directing to ancmppb.
#'
#' @export
#' @keywords ancova deprecated
# ancmpbpb
ancmpbpb<-function(x1,y1,x2,y2,fr1=1,fr2=1,alpha=.05,pts=NA,est=tmean,nboot=NA,
bhop=FALSE,SEED=TRUE,...){
print("This function has been eliminated. Please use ancmppb instead.")
}



#' @title ANCOVA Multiple Comparisons with Percentile Bootstrap
#'
#' @description
#' Compare two independent groups using ANCOVA with multiple covariates and percentile bootstrap.
#' Uses running interval smoother without parametric assumptions about regression form.
#'
#' @param x1 Matrix of covariates for group 1.
#' @param y1 Outcome values for group 1.
#' @param x2 Matrix of covariates for group 2.
#' @param y2 Outcome values for group 2.
#' @param fr1 Span for group 1 smoother. Default is 1.
#' @param fr2 Span for group 2 smoother. Default is 1.
#' @param alpha Significance level. Default is 0.05.
#' @param pts Design points matrix. If NA (default), chosen based on depth in x1.
#' @param est Estimator function. Default is tmean (trimmed mean).
#' @param nboot Number of bootstrap samples. Default is NA.
#' @param bhop Logical. Use Benjamini-Hochberg procedure. Default is TRUE.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param cov.fun Covariance function for finding nearest points. Default is skip.
#'   Options include cov.mve, cov.mcd, etc.
#' @param cop Parameter for ancdes. Default is NULL.
#' @param COV.both Logical. Use covariance from both groups. Default is FALSE.
#' @param pr Logical. Print messages. Default is TRUE.
#' @param ... Additional arguments passed to cov.fun.
#'
#' @return List with test results, design points, and comparison summaries.
#'
#' @export
#' @keywords ancova bootstrap multivariate
#'
#' @examples
#' \dontrun{
#'   x1 <- matrix(rnorm(100), 50, 2)
#'   y1 <- rnorm(50)
#'   x2 <- matrix(rnorm(100), 50, 2)
#'   y2 <- rnorm(50, 0.5)
#'   result <- ancmppb(x1, y1, x2, y2)
#' }
# ancmppb
ancmppb<-function(x1,y1,x2,y2,fr1=1,fr2=1,alpha=.05,pts=NA,est=tmean,nboot=NA,
bhop=TRUE,SEED=TRUE,cov.fun=skip,cop=NULL,COV.both=FALSE,pr=TRUE,...){
#
# Compare two independent  groups using the ancova method
# with multiple covariates.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
# Design points are chosen based on depth of points in x1 if pts=NA
#  Assume data are in x1 y1 x2 and y2
#
#  cov.fun determines the location and
#  scatter matrix used to find closest points to
#  a design point. It is used by ancdes.
#
#  Choices for cov.fun include
#  cov.mve
#  cov.mcd
#  rmba
#  skip
#  tbs
#
#if(pr)print("For the old version of this function, use ancmpbpb")
x1=as.matrix(x1)
y1=as.matrix(y1)
if(ncol(x1)==1)stop("Use a function designed for one covariate only")
x2=as.matrix(x2)
y2=as.matrix(y2)
if(ncol(x1)!=ncol(x2))
stop("Number of covariates must be the same for each group")
xy=elimna(cbind(x1,y1))
p=ncol(x1)
p1=p+1
x1=xy[,1:p]
y1=xy[,p1]
xy=elimna(cbind(x2,y2))
x2=xy[,1:p]
y2=xy[,p1]
x1=as.matrix(x1)
x2=as.matrix(x2)
mval1=cov.fun(x1)
mval2=cov.fun(x2)
if(is.na(pts[1])){
x1<-as.matrix(x1)
if(!COV.both){
if(!is.null(cop))pts<-ancdes(x1,cop=cop)
if(is.null(cop))pts=ancdes(x1,center=mval1$center)
}
if(COV.both){
if(!is.null(cop))pts<-ancdes(rbind(x1,x2),cop=cop)
if(is.null(cop))pts=ancdes(rbind(x1,x2),center=mval1$center)
}
}
pts<-as.matrix(pts)
if(nrow(pts)>=29){
print("WARNING: More than 28 design points")
print("Only first 28 are used.")
pts<-pts[1:28,]
}
n1<-1
n2<-1
vecn<-1
for(i in 1:nrow(pts)){
n1[i]<-length(y1[near3d(x1,pts[i,],fr1,mval1)])
n2[i]<-length(y2[near3d(x2,pts[i,],fr2,mval2)])
}
flag<-rep(TRUE,nrow(pts))
for(i in 1:nrow(pts))if(n1[i]<10 || n2[i]<10)flag[i]<-FALSE
flag=as.logical(flag)
pts<-pts[flag,]
if(sum(flag)==1)pts<-t(as.matrix(pts))
if(sum(flag)==0)stop("No comparable design points found, might increase span.")
mat<-matrix(NA,nrow(pts),7)
dimnames(mat)<-list(NULL,c("n1","n2","DIF","TEST","se","ci.low","ci.hi"))
g1<-list()
ip<-nrow(pts)
ncom<-0
nc2<-ip
con<-matrix(0,nrow=2*ip,ncol=nrow(pts))
for (i in 1:nrow(pts)){
ip<-ip+1
ncom<-ncom+1
nc2<-nc2+1
con[ncom,i]<-1
con[nc2,i]<-0-1
temp<-y1[near3d(x1,pts[i,],fr1,mval1)]
g1[[i]]<-temp[!is.na(temp)]
temp<-y2[near3d(x2,pts[i,],fr2,mval2)]
g1[[ip]]<-temp[!is.na(temp)]
}
flag.est=FALSE
if(identical(est,onestep))flag.est=TRUE
if(identical(est,mom))flag.est=TRUE
if(flag.est)mat<-pbmcp(g1,alpha=alpha,nboot=nboot,est=est,con=con,bhop=bhop,SEED=SEED,...)
if(!flag.est)mat<-linconpb(g1,alpha=alpha,nboot=nboot,est=est,con=con,bhop=bhop,SEED=SEED,...)
list(points=pts,output=mat)
}



#' @title Plot Quantile Shift Effect Size (No Control Group)
#'
#' @description
#' Plot quantile shift (QS) measure of effect size without designating a control group.
#' Compares two groups using quantile-based methods.
#'
#' @param x1,y1 Covariate and outcome for group 1 (single covariate).
#' @param x2,y2 Covariate and outcome for group 2.
#' @param pts Design points. If NULL, chosen automatically.
#' @param q Lower quantile for determining points range. Default is 0.1.
#' @param xout Logical. Remove outliers. Default is FALSE.
#' @param ALL Logical. If TRUE, use many points. Default is TRUE.
#' @param npts Number of points if pts=NULL. Default is 10.
#' @param line Logical. Connect points with lines. Default is TRUE.
#' @param xlab,ylab Axis labels.
#' @param outfun Outlier detection function. Default is outpro.
#' @param ... Additional arguments.
#'
#' @return Matrix with columns: Pts, QS.Effect.Size.
#'
#' @export
#' @keywords ancova effect-size plot
# ancNCE.QS.plot
ancNCE.QS.plot<-function(x1,y1,x2,y2,pts=NULL,q=0.1,xout=FALSE,ALL=TRUE,npts=10,line=TRUE,
xlab='X',ylab='QS.Effect',outfun=outpro,...){
#
# Plot quantile shift measure of effect size
# No control group
#
#    q = lower quantile used to determine the points used,
#
#  pts=NULL: If ALL=TRUE, 20 points are chosen  by this function
#         otherwise three points are used.
#
#
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
if(p>1)stop('Current version allows one covariate only')
p1=p+1
vals=NA
x2<-xy[,1:p]
y2<-xy[,p1]
x2<-as.matrix(x2)
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
n1=length(y1)
n2=length(y2)
n=min(c(n1,n2))
if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE,...)$keep
m<-m[flag,]
n1=nrow(m)
x1<-m[,1:p]
y1<-m[,p1]
x1=as.matrix(x1)
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE,...)$keep
m<-m[flag,]
n2=nrow(m)
n=min(c(n1,n2))
x2<-m[,1:p]
y2<-m[,p1]
x2=as.matrix(x2)
}
if(!is.null(pts))npts=length(pts)
if(is.null(pts)){
if(q<=0 || q>=1)stop('Argument q must be greater than 0 and less than 1')
qu=1-q
L1=qest(x1,q)
L2=qest(x2,q)
U1=qest(x1,qu)
U2=qest(x2,qu)
L=max(L1,L2)
U=min(U1,U2)
if(ALL)pts=seq(L,U,length.out=npts)
else{pts=c(L,(L+U)/2,U)
npts=3
}}
qs=QSanc(x1,y1,x2,y2,pts=pts)
qs=as.vector(matl(qs))
M=cbind(pts,qs)
if(line){
plot(pts,qs,xlab=xlab,ylab=ylab,ylim=c(0,1),type='n')
lines(pts,qs)
}
else
plot(pts,qs,xlab=xlab,ylab=ylab,ylim=c(0,1))
dimnames(M)=list(NULL,c('Pts','QS.Effect.Size'))
M
}


# ancom
#' Omnibus ANCOVA Test
#'
#' Performs an omnibus ANCOVA test comparing two independent groups. Can test either
#' that the regression surface is a horizontal plane through the origin, or that two
#' regression surfaces are parallel. Uses robust depth-based approach.
#'
#' @param x1 Covariate values for group 1 (vector or matrix).
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2 (vector or matrix).
#' @param y2 Response values for group 2.
#' @param dchk Logical. If TRUE, eliminates design space points with halfspace depth of zero.
#' @inheritParams common-params
#' @param plotfun Plotting function to use (default: rplot). Use rplotsm for smoothed plots via bagging.
#' @param PARTEST Logical. If FALSE, tests that regression surface is horizontal plane through origin.
#'   If TRUE, tests that two regression surfaces are parallel.
#' @param tr Trimming proportion. tr=0 recommended for general use; tr>0 may result in poor
#'   Type I error control (default: 0).
#' @param ... Additional arguments passed to plotfun.
#'
#' @return List with:
#'   \item{dstat}{Test statistic}
#'   \item{critd}{Critical value}
#'
#' @note tr=0 is strongly recommended. Higher trimming may compromise Type I error control.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- matrix(rnorm(100), ncol=2)
#' y1 <- rowSums(x1) + rnorm(50)
#' x2 <- matrix(rnorm(100), ncol=2)
#' y2 <- rowSums(x2) + rnorm(50, mean=0.5)
#' # Test for horizontal plane
#' result1 <- ancom(x1, y1, x2, y2, PARTEST=FALSE)
#' # Test for parallel surfaces
#' result2 <- ancom(x1, y1, x2, y2, PARTEST=TRUE)
#' }
ancom<-function(x1,y1,x2,y2,dchk=FALSE,plotit=TRUE,plotfun=rplot,nboot=500,
alpha=.05,SEED=TRUE,PARTEST=FALSE,tr=0,...){
#
# Omnibus ANCOVA
# tr=0 is recommended for general use. tr>0 might result in
#   poor control over the probability of a Type I error.
#  PARTEST=T will test the hypothesis of parallel regression lines.
#
# Setting plotfun=rplotsm will smooth the plots via bagging
#
# dchk=T, points in design space with a halfspace of zero are eliminated
#
# PARTEST=F tests hypothesis that regression surface is a horizontal
# plane through the origin
# PARTEST=T tests the hypothesis that the two regression surfaces
# are parallel.
#
flag1<-rep(TRUE,length(y1))
flag2<-rep(TRUE,length(y2))
if(dchk){
dep1<-fdepth(x2,x1) # depth of points in x1 relative to x2
dep2<-fdepth(x1,x2)
flag1<-(dep1>0)
flag2<-(dep2>0)
}
n1<-sum(flag1)
n2<-sum(flag2)
n<-n1+n2
y<-c(n2*y1[flag1]/n,0-n1*y2[flag2]/n)
x1<-as.matrix(x1)
x1<-x1[flag1,]
x2<-as.matrix(x2)
x2<-x2[flag2,]
x1<-as.matrix(x1)
x2<-as.matrix(x2)
x<-rbind(x1,x2)
if(plotit){
if(ncol(x)<=2)plotfun(x,y,...)
}
if(PARTEST)output<-indt(x,y,nboot=nboot,SEED=SEED)
if(!PARTEST)output<-indt0(x,y,nboot=nboot,alpha=alpha,SEED=SEED)
list(dstat=output$dstat,critd=output$critd)
}

#' @title ANCOVA with Two Covariates and Depth-Based Design Points
#'
#' @description
#' ANCOVA for two independent groups with two covariates, using no parametric assumptions
#' about the regression surface. Uses deeply nested design points and combines p-values
#' across multiple locations.
#'
#' @param x1 Matrix of two covariates for group 1.
#' @param y1 Outcome values for group 1.
#' @param x2 Matrix of two covariates for group 2.
#' @param y2 Outcome values for group 2.
#' @param tr Trim proportion. Default is 0.2.
#' @param test Test function to use. Default is yuen (trimmed means test).
#' @param cr Critical value. If NULL, computed via simulation.
#' @param pr Logical. Print results. Default is TRUE.
#' @param DETAILS Logical. Return detailed output. Default is FALSE.
#' @param cp.value Logical. Compute critical p-value. Default is FALSE.
#' @param plotit Logical. Create plots. Default is FALSE.
#' @param xlab,ylab,zlab Axis labels for plots.
#' @param span Span for LOESS smoothing. Default is 0.75.
#' @param PV Logical. If TRUE, plot p-values; if FALSE, plot differences. Default is TRUE.
#' @param FRAC Fraction of deepest points to use. Default is 0.5.
#' @param MC Logical. Use multicore processing. Default is FALSE.
#' @param q Quantile for median-based tests. Default is 0.5.
#' @param iter Number of iterations for computing critical value. Default is 1000.
#' @param alpha Significance level. Default is 0.05.
#' @param TPM Logical. Use truncated product method (Zaykin et al., 2002). Default is FALSE.
#' @param tau Parameter for truncated product method. Default is 0.05.
#' @param est Estimator function. Default is tmean.
#' @param fr Span for running interval smoother. Default is 1.
#' @param ... Additional arguments passed to test function.
#'
#' @return List with global test statistic, critical value, p-values, and design points.
#'
#' @details
#' Global test statistic is the average of p-values from tests at each design point.
#' With TPM=TRUE, uses truncated product method instead of averaging.
#'
#' @references
#' Zaykin, D.V., Zhivotovsky, L.A., Westfall, P.H., & Weir, B.S. (2002).
#' Truncated product method for combining p-values. Genetic Epidemiology, 22, 170-185.
#'
#' @export
#' @keywords ancova multivariate
#'
#' @examples
#' \dontrun{
#'   x1 <- matrix(rnorm(100), 50, 2)
#'   y1 <- rowSums(x1) + rnorm(50)
#'   x2 <- matrix(rnorm(100), 50, 2)
#'   y2 <- rowSums(x2) + rnorm(50, 0.5)
#'   result <- ancov2COV(x1, y1, x2, y2)
#' }
# ancov2COV
ancov2COV<-function(x1,y1,x2,y2,tr=.2,test=yuen,cr=NULL,pr=TRUE,DETAILS=FALSE,cp.value=FALSE,
plotit=FALSE,xlab='X',ylab='Y',zlab=NULL,span=.75,PV=TRUE,FRAC=.5,MC=FALSE,q=.5,
iter=1000,alpha=.05,TPM=FALSE,tau=.05,est=tmean,fr=1,...){
#
# ANCOVA two covariates, no parametric assumption about the regression surface.
# Use all design points nested deeply within the cloud of data.
# Global test statistic is the average of the p-values based on
# Yuen's test performed at each of the deepest
# points where comparisons can be made.
#
# TPM=TRUE: replace average of the p-values with the test statistic
#  studied by
#  Zaykin, D. V., Zhivotovsky, L. A., Westfall, P. H., & Weir, B.S. (2002).
#   Truncated product method for combining p-values.
#    Genetic Epidemiology 22, 170--185.
#
#  x1 and x2 assumed to be a matrix or data frame with two columns
#
#   if plotit=TRUE then if
#   PV=TRUE create a plot of the p.values as a function of the two covariates
#           using LOESS.
#   if PV=FALSE, plot the difference between the dependent variable as a function of
#         the covariates
#
#   By default, Yuen's test is used, but other tests can be used via the argument
#   test
#
#   pr=TRUE: warning messages will be printed
#
#   DETAILS=TRUE: all p.values are reported for all covariate points used.
#
#   span: the span used by LOESS

#   fr: span used by rung3hat when estimating the difference between
#       predicted Y for group 1 minus predicted Y for group 2.
#
#   FRAC is the fraction of least deep covariate points that are ignored
#
#   MC=TRUE: use a multicore processor to compute a critical value and global p.value
#
#   com.p.value=TRUE: compute p.value based on the global hypothesis of no differences.
#
#   iter=1000: number of iterations used to compute a critical value or global p.value
#
#      Function returns:
#   test.stat: the test statistic. there are two allowed choices: yuen or qcomhd
#   crit.p.val:  the critical value, reject if test.stat<=crit.p.val
#   min.p.val.point:  the values of the covariate that had the smallest p-value
#   min.p.value:  the minimum p-value among all p-values that were computed.
#
com.p.value=cp.value
if(pr)print('Reject if test.stat is less than or equal to crit.value')
if(FRAC<=0 || FRAC >=1)stop('FRAC should be a value between 0 and 1.')
if(ncol(x1)!=2)stop('Should have two covariates')
xy1=elimna(cbind(x1,y1))
x1=xy1[,1:2]
y1=xy1[,3]
xy2=elimna(cbind(x2,y2))
x2=xy2[,1:2]
y2=xy2[,3]
n=min(c(nrow(x1),nrow(x2)))
if(n<50){
if(pr)print('Warning: sample size is less than 50; critical value unknown')
}
if(is.null(cr)){
if(n>=50 & n<=80)cr=as.vector(regYhat(c(50,75),c(.23,.264),xr=n))
if(n>80)cr=.27
}
flag0=is.null(cr)
flag1=FRAC!=.5
flag3=flag0+flag1+com.p.value+TPM
if(flag3>0){
comp.pv=anc2COV.CV(nrow(x1),nrow(x2),iter=iter,MC=MC,TPM=TPM,tau=tau)
MV=sort(comp.pv$M)
ic=round(alpha*iter)
cr=MV[ic]
}
DONE=FALSE
if(identical(test,qcomhd)){
val=ancovampG(x1,y1,x2,y2,DH=TRUE,SEED=TRUE,test=qcomhd,q=q,FRAC=FRAC)
DONE=TRUE
}
if(identical(test,qcomhdMC)){
val=ancovampG(x1,y1,x2,y2,DH=TRUE,SEED=TRUE,test=qcomhdMC,q=q,FRAC=FRAC)
DONE=TRUE
}
if(!DONE){
if(identical(test,yuen))val=ancovampG(x1,y1,x2,y2,DH=TRUE,SEED=TRUE,test=yuen,tr=tr,FRAC=FRAC)
#if(!identical(test,yuen))val=ancovampG(x1,y1,x2,y2,DH=TRUE,SEED=TRUE,test=test,FRAC=FRAC,...)
}
est.dif=rung3hat(x1,y1,fr=fr,pts=val$points,est=est)$rmd-rung3hat(x2,y2,pts=val$points,fr=fr,est=est)$rmd
pavg=mean(val$results[,3])
if(TPM){
vals=val$results[,3]
vals=elimna(vals)
pavg=sum(log(vals[vals<=tau]))
}
mpv=which(val$results[,3]==min(val$results[,3]))
points=val$points
results=val$results
rem.res=results[mpv,3]
rem.points=points[mpv,]
points=cbind(points,est.dif)
dimnames(points)=list(NULL,c('COV 1','COV 2','EST.DIF'))
if(plotit){
if(is.null(zlab)){
if(PV)zlab='P-Value'
if(!PV)zlab='Est.Dif'
}
if(PV)lplot(points[,1:2],results[,3],xlab=xlab,ylab=ylab,zlab=zlab,tick='det',span=span)
if(!PV)lplot(points[,1:2],est.dif,xlab=xlab,ylab=ylab,zlab=zlab,tick='det',span=span)
}
nk=nrow(points)
if(!DETAILS){
points=NULL
results=NULL
}
pval=NULL
if(com.p.value || TPM)pval=1-mean(pavg<=comp.pv$M)
list(num.points.used=nk,test.stat=pavg,crit.value=cr,GLOBAL.p.value=pval,min.p.val.point=rem.points,min.p.value=rem.res,all.points.used=points,all.results=results[,1:3])
}


# ancova
#' Dependent Samples ANCOVA Effect Size
#'
#' Compares two dependent (paired/repeated measures) regression lines using a robust
#' measure of effect size that accounts for the correlation between paired observations.
#' Uses bootstrap estimation of the covariance matrix to compute a marginal measure
#' of dispersion.
#'
#' @param x1 A numeric vector of covariate values at time 1 (or condition 1).
#' @param y1 A numeric vector of outcome values at time 1.
#' @param x2 A numeric vector of covariate values at time 2 (or condition 2),
#'   paired with x1.
#' @param y2 A numeric vector of outcome values at time 2, paired with y1.
#' @param pts A numeric vector specifying design points where effect sizes are
#'   computed. If `NULL`, three points are determined based on data (see `QM` parameter).
#' @inheritParams common-params
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#' @param OM Outlier detection method (not currently implemented in shown code).
#'   Intended choices: 'PRO' (projection), 'BLP' (regression-based), 'MCD', 'BOX'.
#' @param id Column index to ignore when detecting outliers (used with 'DUM' method).
#' @param QM Logical. If `TRUE` (default) and `pts = NULL`, selects three points based
#'   on quantiles. If `FALSE`, uses default points from \code{\link{ancova}}.
#' @param Nboot Number of bootstrap samples for estimating covariance (default: 100).
#' @param xlab Label for x-axis in plot (default: 'X').
#' @param ylab Label for y-axis in plot (default: 'Y').
#' @param pch Plotting character (default: '.').
#'
#' @return A matrix with columns:
#' \describe{
#'   \item{pts}{The design points where effect sizes were computed.}
#'   \item{Effect.size}{The robust effect size at each design point, computed as
#'     \eqn{\sqrt{2} \times (\hat{y}_1 - \hat{y}_2) / \sqrt{n \times SQE}}, where
#'     SQE accounts for variances and covariance of the paired estimates.}
#' }
#'
#' @details
#' This function computes effect sizes for dependent (paired) regression lines.
#' Unlike independent samples ANCOVA, this accounts for the correlation between
#' paired observations through bootstrap estimation of the covariance matrix.
#'
#' The effect size at each design point is:
#' \deqn{ES = \sqrt{2} \times \frac{\hat{y}_1 - \hat{y}_2}{\sqrt{n(Var_1 + Var_2 - 2Cov)}}}
#'
#' where variances and covariance are estimated via bootstrap.
#'
#' If `plotit = TRUE`, displays both regression lines using \code{\link{reg2plot}}.
#'
#' @seealso \code{\link{ancovad.ESci}} for confidence intervals,
#'   \code{\link{Dancova}} for dependent ANCOVA hypothesis testing
#'
#' @export
#' @examples
#' \dontrun{
#' # Paired data example
#' set.seed(123)
#' n <- 30
#' x1 <- rnorm(n)
#' y1 <- 2 + 3*x1 + rnorm(n)
#' x2 <- x1 + rnorm(n, 0, 0.5)  # Correlated covariates
#' y2 <- 2.5 + 3*x2 + rnorm(n)  # Paired outcomes
#'
#' result <- ancovad.ES(x1, y1, x2, y2, Nboot = 500)
#' print(result)
#' }
ancovad.ES<-function(x1,y1,x2,y2,pts=NULL,plotit=TRUE,regfun=tsreg,
SEED=TRUE,OM='PRO',id=NULL,
QM=TRUE,xout=FALSE,Nboot=100,outfun=outpro,xlab='X',ylab='Y',pch='.'){
#
#  Compare two dependent regression lines via a robust measure of effect sizde
#
# pts=NULL: three  points used that are determined based on the data
#
# OM indicates which method is used to remove leverage points. Choices:
 #   PRO=outpro(x,plotit=plotit),         # projection method
 #   PRO.R=outpro.depth(x),   #projection method   random, lower execution time vs outpro
 #   BLP=outblp(x,y,regfun=regfun,plotit=FALSE),       # regression method
 #   DUM=out.dummy(x,y,outfun=outpro.depth,id=id),   #   Detect outliers ignoring  col indicated by argument id
 #   MCD=outDETMCD(x,plotit=plotit),   Uses a robust analog of Mahalanobis distance
   # BOX, Boxplot method using ideal fourths
#
if(SEED)set.seed(2)
xy=elimna(cbind(x1,y1,x2,y2))
n=nrow(xy)
if(ncol(xy)!=4)stop('Only one covariate can  be used')
x1=xy[,1]
y1=xy[,2]
x2=xy[,3]
y2=xy[,4]
if(xout){
flag=c(1:n)
flag1=out.methods(x1,y1,regfun=regfun,plotit=FALSE,id=id)$out.id
flag2=out.methods(x2,y2,regfun=regfun,plotit=FALSE,id=id)$out.id
flag.out=unique(c(flag1,flag2))
if(length(flag.out)>0)flag=flag[-flag.out]
x1<-x1[flag]
y1<-y1[flag]
x2<-x2[flag]
y2<-y2[flag]
}
if(is.null(pts)){
if(QM){
q1=qest(x1,q)
q2=qest(x2,q)
pts=max(q1,q2)
u1=qest(x1,1-q)
u2=qest(x2,1-q)
up=min(u1,u2)
pts[2]=mean(c(pts[1],up))
pts[3]=up
}
if(!QM)pts=ancova(x1,y1,x2,y2,pr=FALSE,plotit=FALSE,SEED=FALSE)$output[,1]
}

# Use a bootstrap estimate of the covariance matrix to get an estimate of
# a marginal measure  of dispersion.
n=length(y1)
if(n!=length(y2))stop('Time one sample size not equal to time two sample size')
npts=length(pts)
e1=matrix(NA,Nboot,npts)
e2=matrix(NA,Nboot,npts)
for(i in 1:Nboot){
flag=sample(n,replace=TRUE)
e1[i,]=regYhat(x1[flag],y1[flag],xr=pts,regfun=regfun)
e2[i,]=regYhat(x2[flag],y2[flag],xr=pts,regfun=regfun)
}
CV=NA
sq1=apply(e1,2,var)
sq2=apply(e2,2,var)
for(j in 1:npts)CV[j]=cov(e1[,j],e2[,j])
SQE=sq1+sq2-2*CV
bot=n*SQE
es=(regYhat(x1,y1,xr=pts,regfun=regfun)-regYhat(x2,y2,xr=pts,regfun=regfun))/sqrt(bot)
es=sqrt(2)*es
mat=cbind(pts,es)
if(plotit)reg2plot(x1,y1,x2,y2,regfun=regfun,xlab=xlab,ylab=ylab)
dimnames(mat)=list(NULL,c('pts','Effect.size'))
mat
}


#' Dependent Samples ANCOVA Effect Size with Confidence Intervals
#'
#' Computes robust effect sizes for two dependent (paired/repeated measures) regression
#' lines with bootstrap confidence intervals and hypothesis tests. Tests whether the
#' effect size significantly differs from zero at specified design points.
#'
#' @param x1 A numeric vector of covariate values at time 1 (or condition 1).
#' @param y1 A numeric vector of outcome values at time 1.
#' @param x2 A numeric vector of covariate values at time 2 (or condition 2),
#'   paired with x1.
#' @param y2 A numeric vector of outcome values at time 2, paired with y1.
#' @param pts A numeric vector specifying design points where effect sizes are
#'   computed. If `NULL`, points are selected automatically (see `QM` parameter).
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#' @inheritParams common-params
#' @param QM Logical. If `FALSE` (default) and `pts = NULL`, uses default points
#'   from \code{\link{ancovad.ES}}. If `TRUE`, selects three points based on
#'   quantiles determined by `ql`.
#' @param ql Quantile level for determining range of design points when `QM = TRUE`
#'   (default: 0.2, uses 0.2 and 0.8 quantiles).
#' @param xlab Label for x-axis in plot (default: 'Pts').
#' @param ylab Label for y-axis in plot (default: 'Y').
#' @param method Method for adjusting p-values for multiple comparisons (default: 'hoch'
#'   for Hochberg). See \code{\link{p.adjust}} for options.
#'
#' @return A matrix with columns:
#' \describe{
#'   \item{pts}{The design points where effect sizes were computed.}
#'   \item{Est.}{Effect size estimate at each point.}
#'   \item{SE}{Bootstrap standard error of the effect size.}
#'   \item{Test.Stat}{Test statistic (Est./SE).}
#'   \item{ci.low}{Lower confidence interval bound.}
#'   \item{ci.up}{Upper confidence interval bound.}
#'   \item{p-value}{Unadjusted p-value for testing H0: effect size = 0.}
#'   \item{p.adjusted}{P-value adjusted for multiple comparisons.}
#' }
#'
#' @details
#' This function extends \code{\link{ancovad.ES}} by adding inferential statistics.
#' For each design point, it computes:
#' \enumerate{
#'   \item Effect size using \code{\link{ancovad.ES}}
#'   \item Bootstrap standard error via \code{ancovad.ES.SEpb}
#'   \item Confidence intervals assuming normality
#'   \item P-values for testing whether effect size differs from zero
#'   \item Adjusted p-values for multiple comparisons
#' }
#'
#' The heteroscedastic effect size accounts for the conditional dispersion of y
#' given x, rather than using an overall measure of location.
#'
#' If `plotit = TRUE`, displays effect size estimates with confidence bands.
#'
#' @seealso \code{\link{ancovad.ES}}, \code{\link{Dancova}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Paired data with effect size inference
#' set.seed(123)
#' n <- 30
#' x1 <- rnorm(n)
#' y1 <- 2 + 3*x1 + rnorm(n)
#' x2 <- x1 + rnorm(n, 0, 0.5)
#' y2 <- 2.5 + 3*x2 + rnorm(n)
#'
#' result <- ancovad.ESci(x1, y1, x2, y2, nboot = 500)
#' print(result)
#' }
ancovad.ESci<-function(x1,y1,x2,y2,pts=NULL,regfun=tsreg,alpha=.05,nboot=100,SEED=TRUE,
QM=FALSE,ql=.2,
xout=FALSE,outfun=outpro,xlab='Pts',ylab='Y',method='hoch',plotit=TRUE){
#
#  Two dependent groups.
#
#  For each specified value for x, compute a heteroscedastic measure of effect
#  So if  x=2, rather can compare the groups using some specified measure of location, use a
#  measure of effect size that takes into account the conditional measure of dispersion of y given x.
#
# if pts=NULL and
# QM=FALSE:  pick covariate points based on default points used by the R funtcion ancovad
#
# ql determines quantiles of x that form the range of points
#  pts can be used to specify the points x, if NULL, the function picks three values
#  The function tests the hypothesis that the measure of effect is zero, no effect.
#
#  iter=100: number of replications used to estimate the standard error.
#
#
xy=elimna(cbind(x1,y1,x2,y2))
n=nrow(xy)
if(ncol(xy)!=4)stop('Only one covariate can  be used')
x1=xy[,1]
y1=xy[,2]
x2=xy[,3]
y2=xy[,4]
if(xout){
flag=c(1:n)
flag1=out.methods(x1,y1,regfun=regfun,plotit=FALSE,id=id)$out.id
flag2=out.methods(x2,y2,regfun=regfun,plotit=FALSE,id=id)$out.id
flag.out=unique(c(flag1,flag2))
if(length(flag.out)>0)flag=flag[-flag.out]
x1<-x1[flag]
y1<-y1[flag]
x2<-x2[flag]
y2<-y2[flag]
}
if(is.null(pts)){
if(QM){
qu=1-ql
q1=qest(x1,ql)
q2=qest(x2,ql)
pts=max(q1,q2)
u1=qest(x1,qu)
u2=qest(x2,qu)
up=min(u1,u2)
pts[2]=mean(c(pts[1],up))
pts[3]=up
}
pts=unique(pts)
if(!QM)pts=ancovad.ES(x1,y1,x2,y2,regfun=regfun,plotit=FALSE,SEED=FALSE)[,1]
}
npts=length(pts)
RES=matrix(NA,nrow=npts,ncol=8)
SE=ancovad.ES.SEpb(x1,y1,x2,y2,regfun=regfun,nboot=nboot,pts=pts,SEED=FALSE)
a=ancovad.ES(x1,y1,x2,y2,regfun=regfun,pts=pts,SEED=FALSE,plotit=plotit,xlab=xlab,ylab=ylab)
RES[,1]=a[,1]
RES[,2]=a[,2]
RES[,3]=SE
RES[,5]=RES[,2]-qnorm(1-alpha/2)*RES[,3]
RES[,6]=RES[,2]+qnorm(1-alpha/2)*RES[,3]
test=RES[,2]/RES[,3]
pv=2*(1-pnorm(abs(test)))
RES[,7]=pv
RES[,4]=test
dimnames(RES)=list(NULL,c('pts','Est.','SE','Test.Stat','ci.low','ci.up','p-value','p.adjusted'))
RES[,8]=p.adjust(RES[,7],method=method)
if(plotit){
xa=c(pts,pts,pts)
ya=c(RES[,2],RES[,5],RES[,6])
plot(xa,ya,xlab=xlab,ylab='ES',type='n')
lines(pts,RES[,5],lty=2)
lines(pts,RES[,2])
lines(pts,RES[,6],lty=2)
}
RES
}




# ancovad.ES.SEpb
ancovad.ES.SEpb<-function(x1,y1,x2,y2,nboot=100,regfun=tsreg,pts=0,SEED=TRUE){
#
#  Estimate standard error
#
n1=length(x1)
npts=length(pts)
if(SEED)set.seed(2)
v=matrix(NA,nrow=nboot,ncol=npts)
for(i in 1:nboot){
id1=sample(n1,replace=TRUE)
X1=x1[id1]
Y1=y1[id1]
X2=x2[id1]
Y2=y2[id1]
v[i,]=ancovad.ES(X1,Y1,X2,Y2,regfun=regfun,pts=pts,plotit=FALSE,SEED=FALSE)[,2]
}
se=apply(v,2,sd,na.rm=TRUE)
se
}



# ancova.ES
#' @title ANCOVA Effect Size (KMS Analog) with Automatic Design Points
#'
#' @description
#' Compute KMS analog effect size for comparing two independent regression lines.
#' If design points not specified, automatically selects points within the overlapping
#' quantile range of both groups' covariate distributions.
#'
#' @param x1 Covariate values for group 1
#' @param y1 Response values for group 1
#' @param x2 Covariate values for group 2
#' @param y2 Response values for group 2
#' @param pts Optional vector of design points. If NULL, uses npts equally-spaced
#'   points in the overlapping range determined by ql quantiles
#' @param plotit Logical; if TRUE, plot effect sizes (default: TRUE)
#' @param xout Logical; if TRUE, remove outliers before analysis (default: FALSE)
#' @param outfun Function for detecting outliers (default: outpro)
#' @param xlab Label for x-axis (default: 'X')
#' @param ylab Label for y-axis (default: 'KMS.effect')
#' @param npts Number of design points to use when pts=NULL (default: 10)
#' @param ql Quantile level for determining overlapping range (default: .1)
#' @param line Logical; if TRUE, connect points with line in plot (default: TRUE)
#' @param ylim Y-axis limits for plot (default: NULL for automatic)
#' @param ... Additional arguments (currently unused)
#'
#' @return Matrix with columns:
#'   \item{pts}{Design points}
#'   \item{Effect.size}{KMS effect size at each point}
#'
#' @note
#' - Only one covariate supported
#' - Design points chosen in [max(q_ql^1, q_ql^2), min(q_{1-ql}^1, q_{1-ql}^2)]
#'   where q_p^j is the p-th quantile of group j
#'
#' @export
#' @examples
#' # Generate data with group differences
#' set.seed(111)
#' x1 <- runif(50, 0, 10)
#' y1 <- 2 + 0.5*x1 + rnorm(50, 0, 1.5)
#' x2 <- runif(50, 0, 10)
#' y2 <- 4 + 0.5*x2 + rnorm(50, 0, 2)
#'
#' # Compute effect sizes at 10 design points
#' es <- ancova.ES(x1, y1, x2, y2, plotit = FALSE)
#' print(es)
#'
#' # Use custom design points
#' es2 <- ancova.ES(x1, y1, x2, y2, pts = seq(2, 8, by = 1), plotit = FALSE)
#'
ancova.ES<-function(x1,y1,x2,y2,pts=NULL,plotit=TRUE,xout=FALSE,outfun=outpro,xlab='X',ylab='KMS.effect',
npts=10,ql=.1,line=TRUE,ylim=NULL,...){
#
#  Comparing two independent regression lines.
#
#  Compute an analog of the KMS measure of effect size
#  for the points indicated by pts.
# pts=NULL: three  points used that are determined based on the data
#
#. plotit=TRUE, plot measures of effect size
#
xy=elimna(cbind(x1,y1))
if(ncol(xy)!=2)stop('Only one covariate can  be used')
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
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
#if(is.null(pts))pts=ancova(x1,y1,x2,y2,pr=FALSE,plotit=FALSE,SEED=FALSE)$output[,1]
if(is.null(pts)){
if(ql<=0 || ql>=1)stop('Argument q must be greater than 0 and less than 1')
qu=1-ql
L1=qest(x1,ql)
L2=qest(x2,ql)
U1=qest(x1,qu)
U2=qest(x2,qu)
L=max(L1,L2)
U=min(U1,U2)
pts=seq(L,U,length.out=npts)
}
s1sq=regIQRsd(x1,y1,pts=pts)
s2sq=regIQRsd(x2,y2,pts=pts)
e1=regYhat(x1,y1,xr=pts,regfun=Qreg)
e2=regYhat(x2,y2,xr=pts,regfun=Qreg)
v1=s1sq^2
v2=s2sq^2
n1=length(y1)
n2=length(y2)
N=n1+n2
q=n1/N
top=(1-q)*v1+q*v2
bot=q*(1-q)
sigsq=top/bot  #  Quantity in brackets KMS p. 176 eq 21.1
es=(e1-e2)/sqrt(sigsq)
mat=cbind(pts,es)
#if(plotit)reg2plot(x1,y1,x2,y2,regfun=Qreg,xlab=xlab,ylab=ylab)
if(plotit)anclinKMS.plot(x1,y1,x2,y2,pts=pts,line=line,xlab=xlab,ylab=ylab,ylim=ylim)
dimnames(mat)=list(NULL,c('pts','Effect.size'))
mat
}

 anclinKMS.plot<-function(x1,y1,x2,y2,pts=NULL,q=0.1,xout=FALSE,ALL=TRUE,npts=10,line=TRUE,
xlab='X',ylab='KMS.Effect',outfun=outpro,ylim=NULL,...){
#
# Plot KMS measure of effect size
#  pts=NULL: If ALL=TRUE, 10 points are chosen  by this function
#         otherwise three points are used.
#
#  The KMS effect size is a heteroscedastic robust analog of Cohen's d
#
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
if(p>1)stop('Current version allows one covariate only')
p1=p+1
vals=NA
x2<-xy[,1:p]
y2<-xy[,p1]
x2<-as.matrix(x2)
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
n1=length(y1)
n2=length(y2)
n=min(c(n1,n2))
if(xout){
m<-cbind(x1,y1)
flag<-outfun(x1,plotit=FALSE,...)$keep
m<-m[flag,]
n1=nrow(m)
x1<-m[,1:p]
y1<-m[,p1]
x1=as.matrix(x1)
m<-cbind(x2,y2)
flag<-outfun(x2,plotit=FALSE,...)$keep
m<-m[flag,]
n2=nrow(m)
n=min(c(n1,n2))
x2<-m[,1:p]
y2<-m[,p1]
x2=as.matrix(x2)
}
if(!is.null(pts))npts=length(pts)
if(is.null(pts)){
if(q<=0 || q>=1)stop('Argument q must be greater than 0 and less than 1')
qu=1-q
L1=qest(x1,q)
L2=qest(x2,q)
U1=qest(x1,qu)
U2=qest(x2,qu)
L=max(L1,L2)
U=min(U1,U2)
if(ALL)pts=seq(L,U,length.out=npts)
else{pts=c(L,(L+U)/2,U)
npts=3
}}

#e=reg.pred(x2,y2,xr=pts,regfun=Qreg,q=.5,xout=FALSE)
qs=ancova.ES(x1,y1,x2,y2,pts=pts,plotit=FALSE)[,2]
M=cbind(pts,qs)
if(is.null(ylim))
ylim=min(-.8,min(qs))
ylim[2]=max(.8,max(qs))
if(line){
plot(pts,qs,xlab=xlab,ylab=ylab,ylim=ylim,type='n')
lines(pts,qs)
}
else
plot(pts,qs,xlab=xlab,ylab=ylab,ylim=ylim)
dimnames(M)=list(NULL,c('Pts','KMS.Effect.Size'))
M
}


# ancova.ES.SEpb
ancova.ES.SEpb<-function(x1,y1,x2,y2,nboot=100,pts=0,SEED=TRUE){
#
#  Estimate standard error assuming normaliy
#
v=NA
n1=length(x1)
n2=length(x2)
npts=length(pts)
if(SEED)set.seed(2)
v=matrix(NA,nrow=nboot,ncol=npts)
for(i in 1:nboot){
id1=sample(n1,replace=TRUE)
id2=sample(n2,replace=TRUE)
X1=x1[id1]
Y1=y1[id1]
X2=x2[id2]
Y2=y2[id2]
v[i,]=ancova.KMS(X1,Y1,X2,Y2,pts=pts,plotit=FALSE)[,2]
}
se=apply(v,2,sd)
se
}


# ancovaG
#' @title Generalized ANCOVA with Running Interval Smoother
#'
#' @description
#' Compare two independent groups using a generalized ANCOVA approach that allows
#' any hypothesis testing method. Uses running interval smoothers without making
#' parametric assumptions about the regression line forms. Groups are compared at
#' specified design points along the covariate range.
#'
#' @param x1 Covariate values for group 1
#' @param y1 Response values for group 1
#' @param x2 Covariate values for group 2
#' @param y2 Response values for group 2
#' @param fr1 Span parameter for group 1 smoother (default: 1)
#' @param fr2 Span parameter for group 2 smoother (default: 1)
#' @inheritParams common-params
#' @param plotit Logical; if TRUE, plot the running means (default: FALSE)
#' @param pts Optional vector of design points where groups are compared. If NULL,
#'   uses 5 equally-spaced points where both groups have at least 12 observations
#' @param sm Logical; if TRUE, use bootstrap bagging for smoother curves (default: FALSE)
#' @param pr Logical; if TRUE, print results (default: TRUE)
#' @param xout Logical; if TRUE, remove outliers before analysis (default: FALSE)
#' @param outfun Function for detecting outliers (default: out)
#' @param test Hypothesis testing function to use for group comparisons (default: medpb2)
#' @param ... Additional arguments passed to the test function
#'
#' @return A list with two components:
#'   \item{[[1]]}{Matrix with columns X (design point), n1 (sample size group 1),
#'     n2 (sample size group 2)}
#'   \item{[[2]]}{List of test results from the specified test function at each design point}
#'
#' @export
#' @examples
#' # Generate example data with covariate effect
#' set.seed(123)
#' x1 <- runif(50, 0, 10)
#' y1 <- 2 + 0.5*x1 + rnorm(50)
#' x2 <- runif(50, 0, 10)
#' y2 <- 3 + 0.5*x2 + rnorm(50)  # Different intercept
#'
#' # Compare groups using default median-based test
#' result <- ancovaG(x1, y1, x2, y2)
#'
#' # Use custom design points and trimmed mean test
#' result2 <- ancovaG(x1, y1, x2, y2, pts = c(2, 5, 8), test = yuen)
#'
ancovaG<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,alpha=.05,plotit=FALSE,pts=NULL,sm=FALSE,
pr=TRUE,xout=FALSE,outfun=out,test=medpb2,...){
#
# This function generalizes the R function ancova so that any hypothesis testing method
# can be used to compare groups at specified design points.
#
# Compare two independent  groups using the ancova method coupled with method
# indicated by the argument test.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#
#  sm=T will create smooths using bootstrap bagging.
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
output=list()
if(is.null(pts[1])){
mat<-matrix(NA,5,3)
dimnames(mat)<-list(NULL,c("X","n1","n2"))
npt<-5
isub<-c(1:5)  # Initialize isub
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
for (i in 1:5){
mat[i,1]=x1[isub[i]]
g1<-y1[near(x1,x1[isub[i]],fr1)]
g2<-y2[near(x2,x1[isub[i]],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
mat[i,2]=length(g1)
mat[i,3]=length(g2)
output[[i]]<-test(g1,g2,...)
}}
if(!is.null(pts[1])){
mat<-matrix(NA,length(pts),3)
dimnames(mat)<-list(NULL,c("X","n1","n2"))
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
for (i in 1:length(pts)){
mat[i,1]=pts[i]
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
mat[i,2]=length(g1)
mat[i,3]=length(g2)
output[[i]]<-test(g1,g2,...)
}}
if(plotit)
runmean2g(x1,y1,x2,y2,fr=fr1,est=mean,tr=tr,sm=sm,xout=xout,outfun=outfun,...)
list(mat,output)
}


#' @title ANCOVA Effect Size Using KMS Analog
#'
#' @description
#' Compare two independent regression lines using an analog of the KMS
#' (Kulinskaya, Morgenthaler, & Staudte) heteroscedastic measure of effect size.
#' This provides a robust analog of Cohen's d that accounts for conditional variance
#' of y given x at specified design points.
#'
#' @param x1 Covariate values for group 1
#' @param y1 Response values for group 1
#' @param x2 Covariate values for group 2
#' @param y2 Response values for group 2
#' @param pts Optional vector of design points for comparison. If NULL, uses three
#'   points determined from the data via ancova()
#' @param plotit Logical; if TRUE, plot the regression lines (default: TRUE)
#' @param xout Logical; if TRUE, remove outliers before analysis (default: FALSE)
#' @param outfun Function for detecting outliers (default: outpro)
#' @param xlab Label for x-axis (default: 'X')
#' @param ylab Label for y-axis (default: 'Y')
#'
#' @return Matrix with columns:
#'   \item{pts}{Design points where effect size is computed}
#'   \item{Effect.size}{KMS effect size analog at each point}
#'
#' @note Only one covariate is supported. The function uses quantile regression
#'   (via Qreg) for estimating regression values and regIQRsd for standard deviation
#'   estimation based on IQR.
#'
#' @export
#' @examples
#' # Generate data with different regression lines
#' set.seed(123)
#' x1 <- runif(40, 0, 10)
#' y1 <- 2 + 0.5*x1 + rnorm(40, 0, 2)
#' x2 <- runif(40, 0, 10)
#' y2 <- 4 + 0.5*x2 + rnorm(40, 0, 2)  # Different intercept
#'
#' # Compute KMS effect size at default points
#' es <- ancova.KMS(x1, y1, x2, y2, plotit = FALSE)
#' print(es)
#'
#' # Specify custom design points
#' es2 <- ancova.KMS(x1, y1, x2, y2, pts = c(2, 5, 8), plotit = FALSE)
#'
ancova.KMS<-function(x1,y1,x2,y2,pts=NULL,plotit=TRUE,xout=FALSE,outfun=outpro,xlab='X',ylab='Y'){
#
#  Comparing two independent regression lines.
#  using an analog of the KMS measure of effect size
#  for the points indicated by pts.
# pts=NULL: three  points used that are determined based on the data
#
xy=elimna(cbind(x1,y1))
if(ncol(xy)!=2)stop('Only one covariate can  be used')
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
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
if(is.null(pts))pts=ancova(x1,y1,x2,y2,pr=FALSE,plotit=FALSE,SEED=FALSE)$output[,1]
s1sq=regIQRsd(x1,y1,pts=pts)
s2sq=regIQRsd(x2,y2,pts=pts)
e1=regYhat(x1,y1,xr=pts,regfun=Qreg)
e2=regYhat(x2,y2,xr=pts,regfun=Qreg)
v1=s1sq^2
v2=s2sq^2
n1=length(y1)
n2=length(y2)
N=n1+n2
q=n1/N
top=(1-q)*v1+q*v2
bot=q*(1-q)
sigsq=top/bot  #  Quantity in brackets KMS p. 176 eq 21.1
es=(e1-e2)/sqrt(sigsq)
mat=cbind(pts,es)
if(plotit)reg2plot(x1,y1,x2,y2,regfun=Qreg,xlab=xlab,ylab=ylab)
dimnames(mat)=list(NULL,c('pts','Effect.size'))
mat
}


#' @title ANCOVA KMS Effect Size with Bootstrap Confidence Intervals
#'
#' @description
#' Compute heteroscedastic KMS effect size measure (robust analog of Cohen's d) for
#' comparing two independent regression lines at specified design points. Provides
#' bootstrap-based confidence intervals and hypothesis tests for zero effect size.
#'
#' @param x1 Covariate values for group 1
#' @param y1 Response values for group 1
#' @param x2 Covariate values for group 2
#' @param y2 Response values for group 2
#' @param pts Optional vector of design points. If NULL, uses three points determined
#'   either from ancova() (default) or from quantiles if QM=TRUE
#' @inheritParams common-params
#' @param nboot Number of bootstrap samples for standard error estimation (default: 100)
#' @param SEED Logical; if TRUE, set random seed for reproducibility (default: TRUE)
#' @param QM Logical; if TRUE, determine design points using quantiles of x rather than
#'   ancova() (default: FALSE)
#' @param ql Quantile level for determining point range when QM=TRUE (default: .2)
#' @param xout Logical; if TRUE, remove outliers before analysis (default: FALSE)
#' @param outfun Function for detecting outliers (default: outpro)
#' @param xlab Label for x-axis (default: 'Pts')
#' @param ylab Label for y-axis (default: 'Y')
#' @param method P-value adjustment method for multiple comparisons (default: 'hoch')
#' @param plotit Logical; if TRUE, plot effect sizes and confidence bands (default: TRUE)
#'
#' @return Matrix with columns:
#'   \item{pts}{Design points}
#'   \item{Est.}{KMS effect size estimate}
#'   \item{Test.Stat}{Test statistic (Est./SE)}
#'   \item{ci.low}{Lower confidence limit}
#'   \item{ci.up}{Upper confidence limit}
#'   \item{p-value}{Unadjusted p-value for testing zero effect}
#'   \item{p.adjusted}{P-value adjusted for multiple comparisons}
#'
#' @note
#' - Only one covariate is supported
#' - Standard errors are multiplied by 0.89 when min(n1,n2) < 100 to reduce bias
#'
#' @export
#' @examples
#' # Generate data with group differences
#' set.seed(456)
#' x1 <- runif(50, 0, 10)
#' y1 <- 2 + 0.5*x1 + rnorm(50, 0, 2)
#' x2 <- runif(50, 0, 10)
#' y2 <- 4 + 0.5*x2 + rnorm(50, 0, 2.5)
#'
#' # Compute effect sizes with CIs
#' result <- ancova.KMSci(x1, y1, x2, y2, nboot = 50, plotit = FALSE)
#' print(result)
#'
ancova.KMSci<-function(x1,y1,x2,y2,pts=NULL,alpha=.05,nboot=100,SEED=TRUE,QM=FALSE,ql=.2,
xout=FALSE,outfun=outpro,xlab='Pts',ylab='Y',method='hoch',plotit=TRUE){
#
#  Two independent groups.
#
#  For each specified value for x, compute a heteroscedastic measure of effect size,
# robust analog of Cohen's d
#  So if  x=2, rather can compare the groups using some specified measure of location, use a
#  measure of effect size that takes into account the conditional measure of dispersion of y given x.
#
# ql determines quantiles of x that form the range of points
#  pts can be used to specify the points x, if NULL, the function picks three values
#  The function tests the hypothesis that the measure of effect is zero, no effect.
#
#  iter=100: number of replications used to estimate the standard error.
#
#
xy=elimna(cbind(x1,y1))
if(ncol(xy)!=2)stop('Only one covariate can  be used')
x1=xy[,1]
y1=xy[,2]
n1=nrow(xy)
xy=elimna(cbind(x2,y2))
n2=nrow(xy)
x2=xy[,1]
y2=xy[,2]
if(xout){
m<-cbind(x1,y1)
if(identical(outfun,reglev))flag=outfun(x1,y1,plotit=FALSE)$keep
else
flag<-outfun(x1,plotit=FALSE)$keep
m<-m[flag,]
x1<-m[,1]
y1<-m[,2]
m<-cbind(x2,y2)
if(identical(outfun,reglev))flag=outfun(x2,y2,plotit=FALSE)$keep
else
flag<-outfun(x2,plotit=FALSE)$keep
m<-m[flag,]
x2<-m[,1]
y2<-m[,2]
}
n1=length(y1)
n2=length(y2)
if(is.null(pts)){
if(QM){
qu=1-ql
q1=qest(x1,ql)
q2=qest(x2,ql)
pts=max(q1,q2)
u1=qest(x1,qu)
u2=qest(x2,qu)
up=min(u1,u2)
pts[2]=mean(c(pts[1],up))
pts[3]=up
}
if(!QM)pts=ancova(x1,y1,x2,y2,pr=FALSE,plotit=FALSE,SEED=FALSE)$output[,1]
}
pts=unique(pts)
npts=length(pts)
RES=matrix(NA,nrow=npts,ncol=7)
SE=ancova.ES.SEpb(x1,y1,x2,y2,nboot=nboot,pts=pts,SEED=SEED)
if(min(n1,n2)<100)SE=.89*SE #Empirically determined to reduce bias
RES[,1:2]=ancova.ES(x1,y1,x2,y2,pts=pts,plotit=plotit,xlab=xlab,ylab=ylab)
RES[,4]=RES[,2]-qnorm(1-alpha/2)*SE
RES[,5]=RES[,2]+qnorm(1-alpha/2)*SE
test=RES[,2]/SE
pv=2*(1-pnorm(abs(test)))
RES[,6]=pv
RES[,3]=test
dimnames(RES)=list(NULL,c('pts','Est.','Test.Stat','ci.low','ci.up','p-value','p.adjusted'))
RES[,7]=p.adjust(RES[,6],method=method)
if(plotit){
xa=c(pts,pts,pts)
ya=c(RES[,2],RES[,4],RES[,5])
plot(xa,ya,xlab=xlab,ylab='KMS',type='n')
lines(pts,RES[,4],lty=2)
lines(pts,RES[,2])
lines(pts,RES[,5],lty=2)
}
RES
}




#' @title Plot ANCOVA KMS Effect Size as Function of Covariate
#'
#' @description
#' Plot the robust KMS measure of effect size as a function of covariate values.
#' If no design points are specified, uses all unique covariate values within the
#' overlapping range of both groups.
#'
#' @param x1 Covariate values for group 1
#' @param y1 Response values for group 1
#' @param x2 Covariate values for group 2
#' @param y2 Response values for group 2
#' @param pts Optional vector of design points. If NULL, uses unique values from
#'   x1 and x2 within their overlapping range
#' @param xlab Label for x-axis (default: 'X')
#' @param ylab Label for y-axis (default: 'Effect Size')
#' @param xout Logical; if TRUE, remove outliers before analysis (default: FALSE)
#' @param outfun Function for detecting outliers (default: outpro)
#' @param pch Point character for plotting if line=FALSE (default: 'x')
#' @param line Logical; if TRUE, connect points with lines (default: TRUE)
#'
#' @return NULL (creates plot as side effect)
#'
#' @note Only one covariate is supported.
#'
#' @export
#' @examples
#' # Generate data with varying effect size
#' set.seed(789)
#' x1 <- runif(60, 0, 10)
#' y1 <- 2 + 0.3*x1 + rnorm(60, 0, 1.5)
#' x2 <- runif(60, 0, 10)
#' y2 <- 3 + 0.4*x1 + rnorm(60, 0, 2)
#'
#' # Plot effect size across covariate range
#' ancova.KMS.plot(x1, y1, x2, y2)
#'
#' # Use specific design points
#' ancova.KMS.plot(x1, y1, x2, y2, pts = seq(1, 9, by = 1))
#'
ancova.KMS.plot<-function(x1,y1,x2,y2,pts=NULL,xlab='X',ylab='Effect Size',xout=FALSE,outfun=outpro,pch='x',line=TRUE){
#
#
# Plot the robust KMS measure of effect size for the covariate values in pts
#
#  pts=NULL, use the uniques values in x1 and x2
#
xy=elimna(cbind(x1,y1))
if(ncol(xy)!=2)stop('Only one covariate can  be used')
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
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
if(is.null(pts)){
id1=duplicated(x1)
id2=duplicated(x2)
X1=sort(x1[!id1])
X2=sort(x2[!id2])
n1=length(X1)
n2=length(X2)
low=max(X1[1],X2[1])
up=min(X1[n1],X2[n2])
X12=sort(c(X1,X2))
flag=(X12>=low & X12<=up)
pts=X12[flag]
}
e=ancova.ES(x1,y1,x2,y2,pts=pts,plotit=FALSE)
plot(e[,1],e[,2],xlab=xlab,ylab=ylab,type='n')
if(line)lines(e[,1],e[,2])
else
points(e[,1],e[,2],pch=pch)
}


#' @title Median-Based ANCOVA for Two Independent Groups
#'
#' @description
#' Compare two independent groups using median-based ANCOVA with running interval
#' smoothers. No parametric assumptions are made about regression line forms.
#' Confidence intervals are adjusted to control familywise error rate across
#' design points (via Studentized Maximum Modulus distribution).
#'
#' @param x1 Covariate values for group 1
#' @param y1 Response values for group 1
#' @param x2 Covariate values for group 2
#' @param y2 Response values for group 2
#' @param fr1 Span parameter for group 1 smoother (default: 1)
#' @param fr2 Span parameter for group 2 smoother (default: 1)
#' @inheritParams common-params
#' @param plotit Logical; if TRUE, plot running medians (default: TRUE)
#' @param pts Optional vector of design points. If NA, uses 5 equally-spaced points
#'   where both groups have at least 12 observations. Maximum 28 points allowed.
#' @param sm Logical; if TRUE, use bootstrap bagging for smoother curves (default: FALSE)
#' @param pr Logical; if TRUE, print warning about p-values (default: TRUE)
#'
#' @return List with components:
#'   \item{output}{Matrix with columns: X (design point), n1, n2, DIF (median
#'     difference), TEST (test statistic), se (standard error), ci.low, ci.hi
#'     (adjusted confidence limits), p.value}
#'   \item{crit}{Critical value from Studentized Maximum Modulus distribution}
#'
#' @note
#' - Confidence intervals control familywise Type I error rate
#' - P-values are NOT adjusted for multiple comparisons
#' - Warns if fewer than 6 observations at any design point
#'
#' @export
#' @examples
#' # Generate data with different medians
#' set.seed(234)
#' x1 <- runif(50, 0, 10)
#' y1 <- 5 + 0.3*x1 + rt(50, df=3)  # Heavy tails
#' x2 <- runif(50, 0, 10)
#' y2 <- 7 + 0.3*x2 + rt(50, df=3)  # Different intercept
#'
#' # Compare at default 5 points
#' result <- ancovam(x1, y1, x2, y2, plotit = FALSE)
#' print(result$output)
#'
#' # Custom design points
#' result2 <- ancovam(x1, y1, x2, y2, pts = c(2, 5, 8), plotit = FALSE)
#'
ancovam<-function(x1,y1,x2,y2,fr1=1,fr2=1,alpha=.05,plotit=TRUE,pts=NA,sm=FALSE,
pr=TRUE){
#
# Compare two independent  groups using an ancova method
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
# This function is designed specifically for
# MEDIANS
#
#  Assume data are in x1 y1 x2 and y2
#
if(pr){
print("NOTE: Confidence intervals are adjusted to control the probability")
print("of at least one Type I error.")
print("But p-values are not")
}
if(is.na(pts[1])){
npt<-5
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
mat<-matrix(NA,5,9)
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","se","ci.low","ci.hi","p.value"))
critv<-NA
critv=qsmm(1-alpha,5,500)
for (i in 1:5){
g1<-y1[near(x1,x1[isub[i]],fr1)]
g2<-y2[near(x2,x1[isub[i]],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test<-msmed(g1,g2)
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
mat[i,4]<-median(g1)-median(g2)
mat[i,5]<-test$test[3]
mat[i,6]<-test$test[5]
cilow<-mat[i,4]-critv*mat[i,6]
cihi<-mat[i,4]+critv*mat[i,6]
mat[i,7]<-cilow
mat[i,8]<-cihi
mat[i,9]<-test$test[6]
}}
if(!is.na(pts[1])){
if(length(pts)>=29)stop("At most 28 points can be compared")
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),9)
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","se","ci.low","ci.hi",
"p.value"))
critv<-NA
if(length(pts)>=2){
#if(alpha==.05)critv<-smmcrit(500,length(pts))
#if(alpha==.01)critv<-smmcrit01(500,length(pts))
#if(is.na(critv))critv<-smmval(rep(999,length(pts)),alpha=alpha)
critv=qsmm(1-alpha,length(pts),500)
}
if(length(pts)==1)critv<-qnorm(1-alpha/2)
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test<-msmed(g1,g2)
mat[i,1]<-pts[i]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
if(length(g1)<=5)print(paste("Warning, there are",length(g1)," points corresponding to the design point X=",pts[i]))
if(length(g2)<=5)print(paste("Warning, there are",length(g2)," points corresponding to the design point X=",pts[i]))
mat[i,4]<-median(g1)-median(g2)
mat[i,5]<-test$test[3]
mat[i,6]<-test$test[5]
cilow<-mat[i,4]-critv*mat[i,6]
cihi<-mat[i,4]+critv*mat[i,6]
mat[i,7]<-cilow
mat[i,8]<-cihi
mat[i,9]<-test$test[6]
}}
if(plotit)
runmean2g(x1,y1,x2,y2,fr=fr1,est=median,sm=sm)
list(output=mat,crit=critv)
}



#' @title ANCOVA for Multiple Covariates with Trimmed Means
#'
#' @description
#' Compare two independent groups using ANCOVA with multiple covariates and running
#' interval smoothers. Uses Yuen's trimmed mean test at design points chosen based
#' on data depth. No parametric assumptions about regression surface forms.
#' Confidence intervals control familywise error rate.
#'
#' @param x1 Matrix of covariate values for group 1 (two or more columns)
#' @param y1 Response values for group 1
#' @param x2 Matrix of covariate values for group 2 (two or more columns)
#' @param y2 Response values for group 2
#' @param fr1 Span parameter for group 1 smoother (default: 1)
#' @param fr2 Span parameter for group 2 smoother (default: 1)
#' @inheritParams common-params
#' @param pts Optional matrix of design points (one row per point). If NA, uses
#'   points with deepest depth in x1 (via ancdes)
#' @param SEED Logical; if TRUE, set seed for reproducible cov.mve results (default: TRUE)
#' @param plotit Logical; if TRUE and exactly 2 covariates, plot design points
#'   marking significant ones (default: FALSE)
#' @param FWE Logical; if TRUE, use familywise error rate for determining
#'   significance; if FALSE, use individual p-values (default: TRUE)
#' @param xlab Label for first covariate in plot (default: 'V1')
#' @param ylab Label for second covariate in plot (default: 'V2')
#'
#' @return List with components:
#'   \item{points}{Matrix of design points used for comparisons}
#'   \item{output}{Matrix with columns: n1, n2, DIF (trimmed mean difference),
#'     TEST (test statistic), se, ci.low, ci.hi (adjusted CI), p.value, adj.p.value
#'     (FWE-adjusted p-value)}
#'   \item{crit}{Critical value from Studentized Maximum Modulus distribution}
#'   \item{sig.pts}{Matrix of significant design points (based on FWE setting)}
#'
#' @note
#' - Design points with fewer than 10 observations in either group are excluded
#' - Warnings issued if fewer than 6 observations at any design point
#' - Uses MVE (Minimum Volume Ellipsoid) for robust covariance estimation
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate data with two covariates
#' set.seed(567)
#' n1 <- 60; n2 <- 60
#' x1 <- matrix(rnorm(n1*2), ncol=2)
#' y1 <- 2 + 0.5*x1[,1] + 0.3*x1[,2] + rnorm(n1)
#' x2 <- matrix(rnorm(n2*2), ncol=2)
#' y2 <- 3 + 0.5*x2[,1] + 0.3*x2[,2] + rnorm(n2)
#'
#' # Compare at data-driven design points
#' result <- ancovamp(x1, y1, x2, y2, plotit = TRUE)
#' print(result$output)
#' }
#'
ancovamp<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,alpha=.05,pts=NA,SEED=TRUE,plotit=FALSE,FWE=TRUE,
xlab='V1',ylab='V2'){
#
# Compare two independent  groups using the ancova method when there are multiple covariates.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
# Design points are chosen based on depth of points in x1 if pts=NA
#  Assume data are in x1 y1 x2 and y2
#  x1 and x2 should be matrices with two or more columns.
#
# FWE=TRUE, controls FWE rate when indicating significant points.
#
if(SEED)set.seed(2) # now cov.mve always returns same result
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
#
#
#
if(is.na(pts[1])){
x1<-as.matrix(x1)
pts<-ancdes(x1)
pts=unique(pts)
}
pts<-as.matrix(pts)
n1<-1
n2<-1
vecn<-1
mval1<-cov.mve(x1)
mval2<-cov.mve(x2)
for(i in 1:nrow(pts)){
n1[i]<-length(y1[near3d(x1,pts[i,],fr1,mval1)])
n2[i]<-length(y2[near3d(x2,pts[i,],fr2,mval2)])
}
flag<-rep(TRUE,nrow(pts))
for(i in 1:nrow(pts))if(n1[i]<10 || n2[i]<10)flag[i]<-FALSE
flag=as.logical(flag)
pts<-pts[flag,]
if(sum(flag)==1)pts<-t(as.matrix(pts))
if(sum(flag)==0)stop('No comparable design points found, might increase span.')
mat<-matrix(NA,nrow(pts),9)
dimnames(mat)<-list(NULL,c('n1','n2','DIF','TEST','se','ci.low','ci.hi','p.value',
'adj.p.value'))
for (i in 1:nrow(pts)){
g1<-y1[near3d(x1,pts[i,],fr1,mval1)]
g2<-y2[near3d(x2,pts[i,],fr2,mval2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test<-yuen(g1,g2,tr=tr)
mat[i,1]<-length(g1)
mat[i,2]<-length(g2)
if(length(g1)<=5)print(paste('Warning, there are',length(g1),' points corresponding to the design point X=',pts[i,]))
if(length(g2)<=5)print(paste('Warning, there are',length(g2),' points corresponding to the design point X=',pts[i,]))
mat[i,3]<-test$dif
mat[i,4]<-test$teststat
mat[i,5]<-test$se
mat[i,8]<-test$p.value
mat[i,9]=1-psmm(abs(test$teststat),nrow(pts),test$df)
if(nrow(pts)>=2)critv<- qsmm(1-alpha,nrow(pts),test$df)           #smmcrit(test$df,nrow(pts))
if(nrow(pts)==1)critv<-qt(.975,test$df)
cilow<-test$dif-critv*test$se
cihi<-test$dif+critv*test$se
mat[i,6]<-cilow
mat[i,7]<-cihi
}
sig.pts=NA
id=NULL
if(sum(mat[,9]<=alpha)>0){
if(FWE)id=which(mat[,9]<=alpha)
if(!FWE)id=which(mat[,8]<=alpha)
sig.pts=pts[id,]
}
if(plotit){
plot(pts[,1],pts[,2],xlab=xlab,ylab=ylab,type='n')
if(FWE)id2=which(mat[,9]>alpha)
if(!FWE)id2=which(mat[,8]>alpha)
points(pts[id2,1],pts[id2,2],pch='*')
if(length(id)>0)points(pts[id,1],pts[id,2],pch='+')
}
sig.pts=unique(sig.pts)
list(points=pts,output=mat,crit=critv,sig.pts=sig.pts)
}

rplot2g<-runmean2g


#' @title Generalized ANCOVA for Multiple Covariates with Custom Tests
#'
#' @description
#' Generalize ancovamp to allow any hypothesis testing method for comparing two
#' independent groups with multiple covariates. Uses running interval smoothers
#' without parametric assumptions. Design points chosen based on data depth.
#' Controls familywise error rate using Hochberg's method.
#'
#' @param x1 Matrix of covariate values for group 1 (can have multiple columns)
#' @param y1 Response values for group 1
#' @param x2 Matrix of covariate values for group 2
#' @param y2 Response values for group 2
#' @param fr1 Span parameter for group 1 smoother (default: 1)
#' @param fr2 Span parameter for group 2 smoother (default: 1)
#' @inheritParams common-params
#' @param pts Optional matrix of design points. If NULL, uses points with deepest
#'   depth in x1 (controlled by DH and FRAC)
#' @param SEED Logical; if TRUE, set seed for reproducibility with bootstrap/MVE (default: TRUE)
#' @param test Hypothesis testing function to use: yuen, qcomhd, qcomhdMC, or
#'   binom2g (default: yuen)
#' @param DH Logical; if TRUE, use only deepest (1-FRAC) design points (default: FALSE)
#' @param FRAC Fraction of deepest points to use when DH=TRUE (default: .5)
#' @param cov.fun Covariance function for depth computation (default: skip.cov)
#' @param ZLIM Logical; if TRUE, set z-axis limits in plot (default: TRUE)
#' @param pr Logical; if TRUE, print progress (default: FALSE)
#' @param q Quantile for qcomhd tests (default: .5 for median)
#' @param plotit Logical; if TRUE and 2 covariates, create 3D plot of p-values (default: FALSE)
#' @param LP Logical; if TRUE, use lplot for 3D visualization; if FALSE, use scatterplot3d (default: FALSE)
#' @param theta Rotation angle for 3D plot (default: 50)
#' @param xlab Label for first covariate (default: ' X1')
#' @param ylab Label for second covariate (default: 'X2 ')
#' @param SCAT Logical; not currently used (default: FALSE)
#' @param zlab Label for z-axis in 3D plot (default: 'p.value ')
#' @param ticktype Type of tick marks for 3D plot (default: 'detail')
#' @param ... Additional arguments passed to test function
#'
#' @return List with components:
#'   \item{points}{Matrix of design points used}
#'   \item{results}{Matrix with columns: n1, n2, p.value, crit.p.value (Hochberg),
#'     Sig (1 if significant by Hochberg, 0 otherwise), est.dif (estimated difference)}
#'   \item{num.sig}{Number of significant comparisons}
#'
#' @note
#' - Design points with fewer than 10 observations in either group are excluded
#' - P-values adjusted using Hochberg's sequentially rejective method
#' - Requires scatterplot3d package if plotit=TRUE and LP=FALSE
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate data with two covariates
#' set.seed(890)
#' n1 <- 70; n2 <- 70
#' x1 <- matrix(rnorm(n1*2), ncol=2)
#' y1 <- 2 + 0.5*x1[,1] + 0.3*x1[,2] + rnorm(n1)
#' x2 <- matrix(rnorm(n2*2), ncol=2)
#' y2 <- 3 + 0.5*x2[,1] + 0.3*x2[,2] + rnorm(n2)
#'
#' # Compare using trimmed means (default)
#' result <- ancovampG(x1, y1, x2, y2)
#' print(result$results)
#'
#' # Compare medians using quantile comparison
#' result2 <- ancovampG(x1, y1, x2, y2, test = qcomhd, q = 0.5)
#' }
#'
ancovampG<-function(x1,y1,x2,y2,fr1=1,fr2=1, tr=.2,
alpha=.05, pts=NULL,SEED=TRUE,test=yuen,DH=FALSE,FRAC=.5,cov.fun=skip.cov,ZLIM=TRUE,
pr=FALSE,q=.5,plotit=FALSE,LP=FALSE,theta=50,xlab=' X1',ylab='X2 ',SCAT=FALSE,zlab='p.value ',ticktype='detail',...){
#
#  ANCOVA:
#
#  This function generalizes the R function ancovamp
#  so that any hypothesis testing method
#  can be used to compare groups at specified design points.
#
# No parametric assumption is made about the form of
# the regression surface--a running interval smoother is used.
# Design points are chosen based on depth of points in x1 if pts=NULL
#  Assume data are in x1 y1 x2 and y2, can have more than one covariate
#
#  test: argument test determines the method that will be used to compare groups.
#       two choices: yuen, qcomhd qcomhdMC
#        Example: test=qcomhd would compare medians using a percentile bootstrap
#   q: controls the quantile used by qcomhd.
#
#  pts: a matrix of design points at which groups are compared
#
#  DH=TRUE, groups compared at the deepest (1-FRAC) design points.
#  if DH=TRUE, there are two covariates and plot=TRUE, plot a smooth with dependent variable=p.values if pv=TRUE
#  or the estimated difference in the measures of location if pv=FALSE
#  If SCAT=TRUE, instead create a scatterplot of the points used in pts, the covariate values
#  and mark the significant ones with *
#
#  theta can be use to rotate the plot.
#
#  SEED=TRUE sets the seed for the random number generator
#       so that same result is always returned when
#        using a bootstrap method or when using cov.mve or cov.mcd
#
#   cov.fun: returns covariance matrix in $cov (e.g.
#   skipcov does not return it in $cov, but skip.cov does. So cov.mve could be used)
#
#   Returns:
#   designs points where comparisons were made.
#   n's used, p-values
#   crit.p.value: critical p-value based on Hochberg's method for controlling FWE
#   sig=1 if a signficant result based on Hochberg; 0 otherwise
#
t.sel=0
if(identical(test,yuen))t.sel=1
if(identical(test,qcomhd))t.sel=2
if(identical(test,qcomhdMC))t.sel=2
if(identical(test,binom2g))t.sel=3
if(t.sel==0)stop('Argument test should be either yuen, qcomhd, qcomhd or binom2g')
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
#
#
#
if(is.null(pts[1])){
x1<-as.matrix(x1)
pts<-ancdes(x1,DH=DH,FRAC=FRAC)
pts=unique(pts)
}
pts<-as.matrix(pts)
n1<-1
n2<-1
vecn<-1
mval1<-cov.fun(x1)
mval2<-cov.fun(x2)
for(i in 1:nrow(pts)){
n1[i]<-length(y1[near3d(x1,pts[i,],fr1,mval1)])
n2[i]<-length(y2[near3d(x2,pts[i,],fr2,mval2)])
}
flag<-rep(TRUE,nrow(pts))
for(i in 1:nrow(pts))if(n1[i]<10 || n2[i]<10)flag[i]<-F
flag=as.logical(flag)
pts<-pts[flag,]
if(sum(flag)==1)pts<-t(as.matrix(pts))
dd=NULL
if(sum(flag)==0){
print('No comparable design points found, might increase span.')
pts=NULL
mat=NULL
dd=NULL
}
if(sum(flag)>0){
mat<-matrix(NA,nrow(pts),6)
mat[,5]=0
dimnames(mat)<-list(NULL,c('n1','n2','p.value','crit.p.value','Sig','est.dif'))
output=list()
for (i in 1:nrow(pts)){
g1<-y1[near3d(x1,pts[i,],fr1,mval1)]
g2<-y2[near3d(x2,pts[i,],fr2,mval2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
temp=NULL
if(identical(test,qcomhd))temp=qcomhd(g1,g2,q=q,plotit=FALSE)
if(identical(test,qcomhdMC))temp=qcomhdMC(g1,g2,q=q,plotit=FALSE)
if(identical(test,yuen))temp=yuen(g1,g2,tr=tr)
if(identical(test,binom2g))temp=binom2g(x=g1,y=g2)
if(is.null(temp$p.value))print('Argument test should be yuen, or qcomhd or qcomhdMC')
mat[i,3]=temp$p.value
output[[i]]=temp
mat[i,1]<-length(g1)
mat[i,2]<-length(g2)
if(t.sel==1)mat[i,6]=mean(g1,tr=tr)-mean(g2,tr=tr)
if(t.sel==2)mat[i,6]=hd(g1,q=q)-hd(g2,q=q)
if(t.sel==3)mat[i,6]=mean(g1)-mean(g2)
if(length(g1)<=5)print(paste('Warning, there are',length(g1),' points corresponding to the design point X=',pts[i,]))
if(length(g2)<=5)print(paste('Warning, there are',length(g2),' points corresponding to the design point X=',pts[i,]))
}
npt=nrow(pts)
dvec=alpha/c(1:npt)
temp2<-order(0-mat[,3])
sigvec<-(mat[temp2,3]>=dvec)
dd=0
if(sum(sigvec)<npt)dd<-npt-sum(sigvec) #number that are sig.
mat[temp2,4]=dvec
flag=mat[,3]<=mat[,4]
if(sum(flag)>0)mat[flag,5]=1
}
if(plotit){
if(!LP){
library(scatterplot3d)
scatterplot3d(pts[,1],pts[,2],mat[,3],xlab=xlab, ylab=ylab,zlab='p.value',zlim=c(0,1))
}
if(LP)lplot(pts,mat[,3],xlab=xlab, ylab=ylab,zlab='p.value',theta=theta,ZLIM=ZLIM,ticktype=ticktype)
}
list(points=pts,results=mat,num.sig=dd)
}


#' ANCOVA Effect Size (KMS Analog) for Two Covariates
#'
#' Compares two independent regression lines based on an analog of the KMS
#' (Kraemer-Morgan-Schork) measure of effect size for data with two covariates.
#' The effect size is computed at specified covariate value combinations.
#'
#' @param x1 A matrix with two columns containing the covariate values for group 1.
#' @param y1 A numeric vector containing the outcome variable for group 1.
#' @param x2 A matrix with two columns containing the covariate values for group 2.
#' @param y2 A numeric vector containing the outcome variable for group 2.
#' @param pts A matrix with two columns specifying the covariate value combinations
#'   at which to compute effect sizes. If `NULL`, points are automatically selected
#'   based on data depth.
#' @param BOTH Logical. If `TRUE` (default), combines x1 and x2 when selecting
#'   design points. If `FALSE`, uses only x1.
#' @param npts Number of design points to use when `pts = NULL` (default: 20).
#' @param profun Function for computing depth of points. Default is `prodepth`
#'   (uses product of marginal depths). Alternative: `pdepth` (deterministic method,
#'   may be slow for large n).
#' @inheritParams common-params
#'
#' @return A matrix with columns for each covariate (X1, X2) and the KMS effect
#'   size. Each row corresponds to one covariate value combination specified in
#'   `pts` or automatically selected.
#'
#' @details
#' The function implements a heteroscedastic effect size measure for comparing two
#' regression lines with two covariates. The KMS effect size at a given covariate
#' combination is computed as:
#'
#' \deqn{ES = \frac{\hat{y}_1 - \hat{y}_2}{\sqrt{\sigma^2_{pooled}}}}
#'
#' where \eqn{\hat{y}_1} and \eqn{\hat{y}_2} are predicted values from quantile
#' regression, and \eqn{\sigma^2_{pooled}} is a weighted pooled variance estimate.
#'
#' When `pts = NULL`, the function automatically selects design points based on
#' depth, ranging from the deepest to least deep points in the covariate space.
#'
#' @seealso \code{\link{ancovap2.KMSci}} for confidence intervals,
#'   \code{\link{ancovap2.KMS.plot}} for visualization
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate example data with two covariates
#' set.seed(123)
#' n1 <- 50; n2 <- 50
#' x1 <- matrix(rnorm(n1 * 2), ncol = 2)
#' y1 <- x1[,1] + x1[,2] + rnorm(n1)
#' x2 <- matrix(rnorm(n2 * 2), ncol = 2)
#' y2 <- x2[,1] + x2[,2] + 0.5 + rnorm(n2)
#'
#' # Compute KMS effect sizes
#' result <- ancovap2.KMS(x1, y1, x2, y2)
#' print(result)
#' }
ancovap2.KMS<-function(x1,y1,x2,y2,pts=NULL,BOTH=TRUE,npts=20,profun=prodepth,
xout=FALSE,outfun=outpro){
#
#  Comparing two independent regression lines.
#  based on an analog of the KMS measure of effect size
#  for the points indicated by pts.
# pts=NULL: three  points used that are determined based on the data
#
# profun=prdepth, random projections are used to measure the depth of a point
#              =pdepth    would use a deterministic method, might have high execution time
#                               if n is large.
#
# BOTH=TRUE: combine x1 and x2 when picking points, otherwise use x1
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
if(is.null(pts)){
if(!BOTH){
d=profun(x1)
ior=order(d)
id=seq(1,min(n1),length.out=npts)
id=floor(id)
pts=x1[ior[id],]
}
if(BOTH){
X=rbind(x1,x2)
d=profun(X)
ior=order(d)
id=seq(1,min(n1),length.out=npts)
id=floor(id)
pts=X[ior[id],]
}
}
pts<-as.matrix(pts)
s1sq=regIQRsd(x1,y1,pts=pts)
s2sq=regIQRsd(x2,y2,pts=pts)
e1=regYhat(x1,y1,xr=pts,regfun=Qreg)
e2=regYhat(x2,y2,xr=pts,regfun=Qreg)
v1=s1sq^2
v2=s2sq^2
n1=length(y1)
n2=length(y2)
N=n1+n2
q=n1/N
top=(1-q)*v1+q*v2
bot=q*(1-q)
sigsq=top/bot  #  Quantity in brackets KMS p. 176 eq 21.1
es=(e1-e2)/sqrt(sigsq)
mat=cbind(pts,es)
lab=NA
for(i in 1:p)lab[i]=paste('X',i)
dimnames(mat)=list(NULL,c(lab,'KMS'))
mat
}


#' ANCOVA KMS Effect Size with Bootstrap Confidence Intervals (Two Covariates)
#'
#' Computes KMS effect sizes for two independent groups with two covariates,
#' along with bootstrap confidence intervals and hypothesis tests. Tests whether
#' the effect size is significantly different from zero at specified covariate
#' value combinations.
#'
#' @param x1 A matrix with two columns containing the covariate values for group 1.
#' @param y1 A numeric vector containing the outcome variable for group 1.
#' @param x2 A matrix with two columns containing the covariate values for group 2.
#' @param y2 A numeric vector containing the outcome variable for group 2.
#' @param pts A matrix with two columns specifying the covariate value combinations.
#'   If `NULL`, points are automatically selected (see `SIMPLE` parameter).
#' @inheritParams common-params
#' @param npts Number of design points when `pts = NULL` and `SIMPLE = FALSE`
#'   (default: 20). Points are evenly spaced from deepest to least deep.
#' @param SIMPLE Logical. When `pts = NULL`, if `TRUE`, uses quartiles of the
#'   marginal distributions of group 1 to determine covariate points. If `FALSE`
#'   (default), uses depth-based selection via `profun`.
#' @param PLOT.ADJ Logical. If `TRUE`, highlights points using adjusted p-values.
#'   If `FALSE` (default), uses unadjusted p-values for highlighting.
#' @param xlab Label for first covariate axis in plot (default: 'X1').
#' @param ylab Label for second covariate axis in plot (default: 'X2').
#' @param BOTH Logical. If `TRUE` (default) and `SIMPLE = FALSE`, combines x1 and
#'   x2 when selecting design points. If `FALSE`, uses only x1.
#' @param profun Function for computing depth of points when `SIMPLE = FALSE`.
#'   Default is `prodepth` (product of marginal depths). Alternative: `pdepth`
#'   (random projections).
#' @param method Method for adjusting p-values for multiple comparisons. Options
#'   include 'hoch' (Hochberg, default), 'BH' (Benjamini-Hochberg), 'bonferroni',
#'   etc. See \code{\link{p.adjust}}.
#'
#' @return A list with components:
#' \describe{
#'   \item{pts}{Matrix of covariate value combinations where effect sizes were computed.}
#'   \item{output}{Matrix with columns: Est. (effect size estimate), Test.Stat
#'     (test statistic), ci.low (lower CI bound), ci.up (upper CI bound),
#'     p-value (unadjusted), p.adjusted (adjusted for multiple comparisons).}
#' }
#'
#' @details
#' For each covariate value combination in `pts`, the function:
#' \enumerate{
#'   \item Computes the KMS effect size using \code{\link{ancovap2.KMS}}
#'   \item Estimates standard errors via percentile bootstrap
#'   \item Applies a sample-size adjustment factor to improve accuracy
#'   \item Computes confidence intervals and p-values
#'   \item Adjusts p-values for multiple comparisons
#' }
#'
#' The plot displays all tested points with asterisks, and circles significant
#' points (those with p ≤ 0.05).
#'
#' @seealso \code{\link{ancovap2.KMS}}, \code{\link{ancovap2.KMS.plot}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n1 <- 50; n2 <- 50
#' x1 <- matrix(rnorm(n1 * 2), ncol = 2)
#' y1 <- x1[,1] + x1[,2] + rnorm(n1)
#' x2 <- matrix(rnorm(n2 * 2), ncol = 2)
#' y2 <- x2[,1] + x2[,2] + 0.5 + rnorm(n2)
#'
#' # Compute effect sizes with CIs
#' result <- ancovap2.KMSci(x1, y1, x2, y2, nboot = 500)
#' print(result$output)
#' }
ancovap2.KMSci<-function(x1,y1,x2,y2,pts=NULL,alpha=.05,nboot=100,SEED=TRUE,npts=20,
SIMPLE=FALSE,PLOT.ADJ=FALSE,
plotit=TRUE,xlab='X1',ylab='X2',BOTH=TRUE,profun=prodepth,
xout=FALSE,outfun=outpro,method='hoch'){
#
#  Two independent groups, have two covariates.
#
#  For each specified value for x, stored in pts, compute a heteroscedastic measure of effect
#
# if pts=NULL
# SIMPLE=TRUE: use the quartiles of the marginal distributions of group 1
#  to determine the covariate points used,
# SIMPLE=FALSE
#  points are chosen based on the depths of the points, which is computed
#  by the R function indicated by the argument profun.
# The default  is profun=prodepth.
# To use a  random collection of projections, set
# profun=pdepth.depth
#
# npts=20  When SIMPLE=FALSE, means 20 points are selected evenly spaced
#  between the deepest point and the
#  least deep point.
#
#  The function tests the hypothesis that the measure of effect is zero, no effect.
#
#  iter=100: number of replications used to estimate the standard error.
#
# BOTH=TRUE and SIMPLE=FALSE:
# combine x1 and x2 when picking points, otherwise use x1
#
FLAG=FALSE
if(!is.null(pts))FLAG=TRUE
p=ncol(x1)
if(p!=2)stop('Current version is limited to two covariates')
p1=p+1
xy=elimna(cbind(x1,y1))
x1=xy[,1:p]
y1=xy[,p1]
n1=nrow(xy)
xy=elimna(cbind(x2,y2))
n2=nrow(xy)
x2=xy[,1:p]
y2=xy[,p1]
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
n1=length(y1)
n2=length(y2)
n=n1+n2
if(is.null(pts)&SIMPLE) pts=cbind(qest(x1[,1],c(.25,.5,.75)),qest(x1[,2],c(.25,.5,.75)))

if(is.null(pts)){
if(!BOTH){
d=profun(x1)
ior=order(d)
id=seq(1,n1,length.out=npts)
id=floor(id)
pts=x1[ior[id],]
}
if(BOTH){
X=rbind(x1,x2)
d=profun(X)
ior=order(d)
id=seq(1,n,length.out=npts)
id=floor(id)
pts=X[ior[id],]
}
}
adj=matrix(c(20,     0.673367,
30,     0.8048804,
40,    0.8452348,
50,      0.8702816,
75,        0.8975298,
100,     0.9231938,
125,      0.9363285,
150,      0.940,
175,   0.9438881,
200,   0.9492541,
250,  0.9546365,
300,  0.9527324),byrow=TRUE,ncol=2)
nmid=(n1+n2)/2
if(max(n1,n2)>300)b.adj=.975
else
b.adj=lplot.pred(1/adj[,1],adj[,2],1/nmid)$yhat
npts=nrow(pts)
RES=matrix(NA,nrow=npts,ncol=6)
SE=ancovap2.KMS.SEpb(x1,y1,x2,y2,nboot=nboot,pts=pts,SEED=SEED)
SE=b.adj*SE
RES[,1]=ancovap2.KMS(x1,y1,x2,y2,pts=pts)[,p1]
RES[,3]=RES[,1]-qnorm(1-alpha/2)*SE
RES[,4]=RES[,1]+qnorm(1-alpha/2)*SE
test=RES[,1]/SE
pv=2*(1-pnorm(abs(test)))
RES[,5]=pv
RES[,2]=test
dimnames(RES)=list(NULL,c('Est.','Test.Stat','ci.low','ci.up','p-value','p.adjusted'))
RES[,6]=p.adjust(RES[,5],method=method)
ip=which(RES[,5]<=.05)
sig.output=NULL
sig.points=NULL
if(length(ip)>0){
sig.output=RES[ip,]
sig.points=pts[ip,]
}
if(FLAG)sig.output=RES
if(plotit){
plot(x1[,1],x1[,2],xlab=xlab,ylab=ylab,pch='.') #type='n')
if(PLOT.ADJ)ip=which(RES[,6]<=.05)
else ip=which(RES[,5]<=.05)
points(pts[,1],pts[,2],pch='*')
if(length(ip)>0)points(pts[ip,],pch='o')
}
RES
list(pts=pts,output=RES)
}




#' Plot ANCOVA KMS Effect Size as Function of Two Covariates
#'
#' Creates a 3D visualization of the KMS effect size as a function of two
#' covariates when comparing two independent regression lines. Automatically
#' removes leverage points for robust estimation.
#'
#' @param x1 A matrix with two columns containing the covariate values for group 1.
#' @param y1 A numeric vector containing the outcome variable for group 1.
#' @param x2 A matrix with two columns containing the covariate values for group 2.
#' @param y2 A numeric vector containing the outcome variable for group 2.
#' @param pts A matrix with two columns specifying covariate value combinations.
#'   If `NULL` (default), uses the deepest 90% of the pooled data from x1 and x2.
#' @param xlab Label for first covariate axis (default: 'X1').
#' @param ylab Label for second covariate axis (default: 'X2').
#' @param zlab Label for vertical axis showing effect size (default: 'Effect Size').
#' @inheritParams common-params
#' @param theta Azimuthal viewing angle for 3D plot in degrees (default: 50).
#'   Controls horizontal rotation of the plot.
#' @param phi Colatitude viewing angle for 3D plot in degrees (default: 25).
#'   Controls vertical tilt of the plot.
#' @param REV Logical. If `FALSE` (default), uses covariates in original order.
#'   If `TRUE`, reverses the order of covariates. This can affect the appearance
#'   of the LOESS surface fitted to the effect sizes.
#'
#' @return A list with component:
#' \describe{
#'   \item{Number_of_points_used_is}{The number of covariate value combinations
#'     used in computing and plotting the effect sizes.}
#' }
#'
#' @details
#' The function creates a 3D plot showing how the KMS effect size varies across
#' the two-dimensional covariate space. The plotting method depends on the number
#' of points:
#' \itemize{
#'   \item If N < 25: Uses \code{scatterplot3d} for a simple 3D scatter plot
#'   \item If N ≥ 25: Uses \code{lplot} to fit a LOESS surface and create a
#'     smoother visualization
#' }
#'
#' Leverage points are automatically detected and removed using the function
#' specified by `outfun` (default: \code{outpro}) to ensure robust effect size
#' estimates.
#'
#' The `REV` parameter can be useful when the orientation of the LOESS surface
#' doesn't clearly show the relationship - reversing the covariate order may
#' provide a clearer view.
#'
#' @seealso \code{\link{ancovap2.KMS}}, \code{\link{ancovap2.KMSci}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n1 <- 50; n2 <- 50
#' x1 <- matrix(rnorm(n1 * 2), ncol = 2)
#' y1 <- x1[,1] + x1[,2] + rnorm(n1)
#' x2 <- matrix(rnorm(n2 * 2), ncol = 2)
#' y2 <- x2[,1] + x2[,2] + 0.5 + rnorm(n2)
#'
#' # Create 3D plot of effect sizes
#' ancovap2.KMS.plot(x1, y1, x2, y2)
#' }
ancovap2.KMS.plot<-function(x1,y1,x2,y2,pts=NULL,xlab='X1',ylab='X2',zlab='Effect Size',
xout=FALSE,outfun=outpro,SEED=TRUE, theta = 50, phi = 25,REV=FALSE){
#
#
# Two covariates, plot the KMS measure of effect size
# The function automatically removes leverage points.
#
#  The function computes the KMS measure of effect size for the points in
#  pts
#  and plots the results. if
#. pts=NULL, the function picks the deepest 90% of the pooled
#   data in x1 and x2
#
#  REV=FALSE: The plot created by LOESS is impacted by which independent
#  variable is first in the matrix
#. pts
# To switch which is first, set REV=TRUE
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
e=ancovap2.KMS(x1,y1,x2,y2,pts=pts)[,3]
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


#Wrap-upp functions

#Spherical Equivalent Prediction Error Dataset
#Choose this file (SEQ_PE):  source(file.choose())
#for dependent dataset:  SEQ_PE(dependent=T)
#for independent dataset:  SEQ_PE(dependent=F)


# ancovap2.KMS.SEpb
ancovap2.KMS.SEpb<-function(x1,y1,x2,y2,nboot=100,pts=NULL,SEED=TRUE){
#
#  Estimate standard error
#
if(is.null(pts))stop('No points were specified')
n1=nrow(x1)
n2=nrow(x2)
p=ncol(x1)+1
npts=nrow(pts)
if(SEED)set.seed(2)
v=matrix(NA,nrow=nboot,ncol=nrow(pts))
for(i in 1:nboot){
id1=sample(n1,replace=TRUE)
id2=sample(n2,replace=TRUE)
X1=x1[id1,]
Y1=y1[id1]
X2=x2[id2,]
Y2=y2[id2]
v[i,]=ancovap2.KMS(X1,Y1,X2,Y2,pts=pts)[,p]
}
se=apply(v,2,sd)
se
}


# ancovap2.wmw.plot
#' @title Bootstrap ANCOVA with Unconditional Resampling (Version 2)
#'
#' @description
#' Compare two independent groups using ANCOVA with bootstrap resampling from (x,y)
#' pairs rather than conditioning on x. This unconditional resampling approach tends
#' to have more power than standard ancova(). Uses running interval smoothers with
#' no parametric assumptions. Supports only one covariate.
#'
#' @param x1 Covariate values for group 1
#' @param y1 Response values for group 1
#' @param x2 Covariate values for group 2
#' @param y2 Response values for group 2
#' @param fr1 Span parameter for group 1 smoother (default: 1)
#' @param fr2 Span parameter for group 2 smoother (default: 1)
#' @param p.crit Critical p-value for significance. If NULL, determined via simulation
#'   or fast approximation (default: NULL)
#' @param padj Logical; if TRUE, use adjusted p-values via Hochberg/Hommel; if FALSE,
#'   use critical p-value approach for better power (default: TRUE)
#' @param pr Logical; if TRUE, print messages (default: TRUE)
#' @param method P-value adjustment method: 'hochberg' or 'hommel' (default: 'hochberg')
#' @param FAST Logical; if TRUE and alpha=.05, use fast approximation for p.crit;
#'   if FALSE, compute via simulation (default: TRUE)
#' @param est Estimator function to use (default: tmean for trimmed mean)
#' @inheritParams common-params
#' @param plotit Logical; if TRUE, plot running smoothers (default: TRUE)
#' @param xlab Label for x-axis (default: 'X')
#' @param ylab Label for y-axis (default: 'Y')
#' @param qpts Logical; if TRUE, use quantiles of x1 as design points; if FALSE,
#'   use ancova() default points (default: FALSE)
#' @param qvals Quantile values for design points when qpts=TRUE (default: c(.25,.5,.75))
#' @param sm Logical; if TRUE, use bootstrap bagging for smoother curves (default: FALSE)
#' @param xout Logical; if TRUE, remove outliers in x before analysis (default: FALSE)
#' @param eout Logical; if TRUE, eliminate all outliers when plotting (default: FALSE)
#' @param outfun Function for detecting outliers (default: out)
#' @param LP Logical; if TRUE, smooth running interval smoother again using lplot (default: TRUE)
#' @param nboot Number of bootstrap samples (default: 500)
#' @param SEED Logical; if TRUE, set seed for reproducibility (default: TRUE)
#' @param nreps Number of replications for determining critical p-value (default: 2000)
#' @param MC Logical; if TRUE, use parallel processing via mclapply (default: FALSE)
#' @param nmin Minimum sample size at design points (default: 12)
#' @param q Quantile for Harrell-Davis estimator when est=hd (default: .5 for median)
#' @param SCAT Logical; if TRUE, add scatterplot (default: TRUE)
#' @param pch1 Point character for group 1 in plot (default: '*')
#' @param pch2 Point character for group 2 in plot (default: '+')
#' @param ... Additional arguments passed to estimator function
#'
#' @return List with components:
#'   \item{pts}{Design points where comparisons made}
#'   \item{p.value}{P-values at each design point}
#'   \item{p.crit}{Critical p-value (if padj=FALSE)}
#'   \item{num.sig}{Number of significant comparisons}
#'   \item{n1}{Sample size group 1}
#'   \item{n2}{Sample size group 2}
#'
#' @note
#' - Only one covariate supported
#' - FAST=TRUE uses approximation valid for alpha=.05 and n <= 800
#' - For faster execution with large samples, set MC=TRUE (requires multicore)
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate data with covariate relationship
#' set.seed(345)
#' x1 <- runif(60, 0, 10)
#' y1 <- 2 + 0.5*x1 + rnorm(60, 0, 1.5)
#' x2 <- runif(60, 0, 10)
#' y2 <- 3 + 0.5*x2 + rnorm(60, 0, 1.5)
#'
#' # Compare using bootstrap ANCOVA
#' result <- ancovaV2(x1, y1, x2, y2, nboot = 200, plotit = FALSE)
#' print(result$p.value)
#' }
#'
ancovaV2<-function(x1=NULL,y1=NULL,x2=NULL,y2=NULL,fr1=1,fr2=1,p.crit=NULL,padj=TRUE, pr=TRUE,
method='hochberg',FAST=TRUE,
est=tmean,alpha=.05,plotit=TRUE,xlab='X',ylab='Y',qpts=FALSE,qvals=c(.25,.5,.75),sm=FALSE,
xout=FALSE,eout=FALSE,outfun=out,LP=TRUE,nboot=500,SEED=TRUE,nreps=2000,MC=FALSE,
nmin=12,q=.5,SCAT=TRUE,pch1='*',pch2='+',...){
#
# Compare two independent  groups using the ancova method.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
# One covariate only is allowed.
#
#  Like ancova, only bootstrap samples are obtained by resampling
#  from c(x1,y1) and c(x2,y2) rather than conditioning on the x value as done by ancova.
#   This function tends to have more power than ancova.
#
#
#  padj=TRUE, p-values are adjusted using the method indicated by the argument
#  method. By default, Hochberg's method is used Setting method='hommel' will use Hommel's method.
# By default, adjusted p-values are not reported; by default the critical p-value, p_c, is used to get better power.
#
#   x1 y1 are measures for group 1
#   x2 y2 are measures for group 2
#
#  LP=T, when plotting, running interval smoother is smoothed again using lplot.
#  sm=T will create smooths using bootstrap bagging.
#
#  qvals: Determines which quantiles of the covariate x1 will be used when
#   qpts=TRUE
#  qpts=FALSE means covariate chosen as done by the function ancova
#
#  q=.5 means when est=hd (Harrell-Davis estimator), median is estimated.
#
#  eout=TRUE will eliminate all outliers when plotting.
#
#  nreps: indicates number of replications used to determine a critical value and p-value.
#
if(SEED)set.seed(2)
if(pr){
if(!FAST){
if(!MC & is.null(p.crit))print('For faster execution time, set MC=TRUE, assuming a multi-core processor is available')
}}
#if(padj)p.crit=0
if(padj)p.crit=NULL
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
xy1=elimna(cbind(x1,y1))
xy2=elimna=cbind(x2,y2)
n1=nrow(xy1)
n2=nrow(xy2)
if(plotit){
ef=identical(est,hd)
if(!ef)runmean2g(xy1[,1],xy1[,2],xy2[,1],xy2[,2],fr=fr1,est=est,sm=sm,xout=xout,LP=LP,eout=eout,
xlab=xlab,ylab=ylab,SCAT=SCAT,pch1=pch1,pch2=pch2,...)
if(ef)runmean2g(xy1[,1],xy1[,2],xy2[,1],xy2[,2],fr=fr1,est=hd,sm=sm,xout=xout,LP=LP,q=q,eout=eout,
xlab=xlab,ylab=ylab,SCAT=SCAT,pch1=pch1,pch2=pch2,...)
}
if(is.null(p.crit)){
if(FAST){
if(alpha==.05){
nm=max(c(n1,n2))
if(nm<=800){
nv=c(50,60,80,100,200,300,500,800)
if(qpts){
pv=c(.02709,.0283,.0306,.02842,.02779,.02410,.02683,.01868,.02122)
p.crit=lplot.pred(1/nv,pv,1/n1)$yhat
}
if(!qpts){
pv=c(.020831,.017812,.015796,.014773,.012589,.015664,.011803,.012479)
p.crit=lplot.pred(1/nv,pv,1/n1)$yhat
}}
}}}
if(is.null(p.crit)){
p.crit=ancovaV2.pv(n1,n2,nreps=nreps,MC=MC,qpts=qpts,est=est,qvals=qvals,SEED=SEED,
alpha=alpha,nboot=nboot)$p.crit
}
pts=NULL
if(qpts)for(i in 1:length(qvals))pts=c(pts,qest(xy1[,1],qvals[i]))
if(!qpts)pts=ancova(x1,y1,x2,y2,pr=FALSE,plotit=FALSE)$output[,1]
if(SEED)set.seed(2)
ef=identical(est,hd)
est1=NA
est2=NA
J=length(pts)
est1=matrix(NA,nrow=nboot,ncol=J)
est2=matrix(NA,nrow=nboot,ncol=J)
#
data1=matrix(sample(n1,size=n1*nboot,replace=TRUE),ncol=nboot,nrow=n1)
data2=matrix(sample(n2,size=n2*nboot,replace=TRUE),ncol=nboot,nrow=n2)
if(!MC){
if(!ef){
est1=apply(data1,2,DancGLOB_sub,xy=xy1[,1:2],pts=pts,est=est,fr=fr1,nmin=nmin,...)
est2=apply(data2,2,DancGLOB_sub,xy=xy2[,1:2],pts=pts,est=est,fr=fr2,nmin=nmin,...)
}
if(ef){
est1=apply(data1,2,DancGLOB_sub,xy=xy1[,1:2],pts=pts,est=hd,fr=fr1,nmin=nmin,q=q,...)
est2=apply(data2,2,DancGLOB_sub,xy=xy2[,1:2],pts=pts,est=hd,fr=fr2,nmin=nmin,q=q,...)
}
est1=t(as.matrix(est1))
est2=t(as.matrix(est2))
}
if(MC){
data1=listm(data1)
data2=listm(data2)
if(!ef){
est1=mclapply(data1,DancGLOB_sub,xy=xy1[,1:2],pts=pts,est=est,fr=fr1,nmin=nmin,...)
est2=mclapply(data2,DancGLOB_sub,xy=xy2[,1:2],pts=pts,est=est,fr=fr2,nmin=nmin,...)
}
if(ef){
est1=mclapply(data1,DancGLOB_sub,xy=xy1[,1:2],pts=pts,est=hd,fr=fr1,nmin=nmin,q=q,...)
est2=mclapply(data2,DancGLOB_sub,xy=xy2[,1:2],pts=pts,est=hd,fr=fr2,nmin=nmin,q=q,...)
}
est1=t(matl(est1))
est2=t(matl(est2))
}
pv=NULL
if(J==1){
est1=t(as.matrix(est1))
est2=t(as.matrix(est2))
}
for(j in 1:J){
pv[j]=mean(est1[,j]<est2[,j],na.rm=TRUE)+.5*mean(est1[,j]==est2[,j],na.rm=TRUE)
pv[j]=2*min(c(pv[j],1-pv[j]))
}
pvadj=rep(NA,length(pts))
if(padj)pvadj=p.adjust(pv,method=method)
pvm=cbind(pts,pv,pvadj)
dimnames(pvm)=list(NULL,c('X','p.values','p.adjusted'))
list(output=pvm,n=c(n1,n2),p.crit=p.crit)
}

# ancovaV2.pv
ancovaV2.pv<-function(n1,n2,nreps=2000,MC=FALSE,qpts=FALSE,qvals = c(0.25, 0.5, 0.75),
nboot=500,SEED=TRUE,est=tmean,alpha=.05){
iter=nreps
if(SEED)set.seed(45)
xy=list()
for(i in 1:iter){
xy[[i]]=list()
xy[[i]][[1]]=rnorm(n1)
xy[[i]][[2]]=rnorm(n1)
xy[[i]][[3]]=rnorm(n2)
xy[[i]][[4]]=rnorm(n2)
}
if(!MC)pv=lapply(xy,ancovaV2pv.sub,qpts=qpts,qvals=qvals,nboot=nboot,MC=FALSE,est=est)
if(MC){
pv=mclapply(xy,ancovaV2pv.sub,qpts=qpts,qvals=qvals,nboot=nboot,MC=FALSE,est=est)
}
pv=as.vector(matl(pv))
p=hd(pv,q=alpha)
list(p.crit=p)
}


# ancovaV2pv.sub
ancovaV2pv.sub<-function(xy,qpts=FALSE,qvals = c(0.25, 0.5, 0.75),nboot=500,MC=TRUE,
est=tmean){
res=ancovaV2(xy[[1]],xy[[2]],xy[[3]],xy[[4]],est=est,plotit=FALSE,p.crit=.03,SEED=TRUE,qpts=qpts,
nboot=nboot,MC=MC)
rm=min(res$output[,2])
rm
}



#' Compare Regression Lines at Design Points - Parametric ANCOVA
#'
#' Compares regression lines of two independent groups at specified design points
#' using robust regression. Handles both single and multiple covariates. Default
#' uses Theil-Sen estimator.
#'
#' @param x1 Covariate value(s) for group 1 (vector or matrix).
#' @param y1 Response values for group 1.
#' @param x2 Covariate value(s) for group 2 (vector or matrix).
#' @param y2 Response values for group 2.
#' @param pts Design points where regression lines are compared. For p>1 predictors,
#'   should be a matrix with p columns. If NULL, determined automatically.
#' @param regfun Regression function to use (default: tsreg for Theil-Sen).
#' @param fr1 Span for group 1 running interval smoother (default: 1).
#' @param fr2 Span for group 2 running interval smoother (default: 1).
#' @inheritParams common-params
#' @param xout Logical. Remove outliers from covariates before analysis.
#' @param xlab,ylab Axis labels for plots.
#' @param ... Additional arguments passed to ancts (single covariate) or anctsmp (multiple covariates).
#'
#' @return Results from ancts (single covariate) or anctsmp (multiple covariates).
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- 2*x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- x2 + rnorm(50, mean=1)
#' result <- ancpar(x1, y1, x2, y2, regfun=tsreg)
#' }
#' Parametric ANCOVA Wrapper (Single or Multiple Covariates)
#'
#' Wrapper function that automatically calls the appropriate robust ANCOVA function
#' based on the number of covariates. Calls \code{\link{ancts}} for one covariate
#' or \code{\link{anctsmp}} for multiple covariates.
#'
#' @param x1 A numeric vector or matrix of covariate values for group 1.
#' @param y1 A numeric vector of outcome values for group 1.
#' @param x2 A numeric vector or matrix of covariate values for group 2.
#' @param y2 A numeric vector of outcome values for group 2.
#' @param pts Design points for comparisons. For single covariate: a vector.
#'   For multiple covariates: a matrix with one column per covariate.
#'   If `NULL`, points are chosen automatically.
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#' @param fr1 Span for group 1 running interval (default: 1).
#' @param fr2 Span for group 2 running interval (default: 1).
#' @inheritParams common-params
#' @param xlab,ylab Axis labels for plots.
#' @param ... Additional arguments passed to \code{\link{ancts}} or \code{\link{anctsmp}}.
#'
#' @return The output from either \code{\link{ancts}} (for 1 covariate) or
#'   \code{\link{anctsmp}} (for >1 covariates). See those functions for details.
#'
#' @details
#' This is a convenience wrapper that automatically selects the appropriate
#' robust ANCOVA function:
#' \itemize{
#'   \item If x1 and x2 are vectors (ncol = 1): calls \code{\link{ancts}}
#'   \item If x1 and x2 are matrices (ncol > 1): calls \code{\link{anctsmp}}
#' }
#'
#' Both functions use robust regression (default: Theil-Sen estimator) to
#' compare regression lines/surfaces at specified design points.
#'
#' @seealso \code{\link{ancts}}, \code{\link{anctsmp}}
#'
#' @export
#' @examples
#' \dontrun{
#' # One covariate
#' x1 <- rnorm(30); y1 <- 2 + 3*x1 + rnorm(30)
#' x2 <- rnorm(30); y2 <- 1 + 3*x2 + rnorm(30)
#' result1 <- ancpar(x1, y1, x2, y2)
#'
#' # Two covariates
#' x1_mat <- matrix(rnorm(60), ncol = 2)
#' y1 <- 2 + 3*x1_mat[,1] + x1_mat[,2] + rnorm(30)
#' x2_mat <- matrix(rnorm(60), ncol = 2)
#' y2 <- 1 + 3*x2_mat[,1] + x2_mat[,2] + rnorm(30)
#' result2 <- ancpar(x1_mat, y1, x2_mat, y2)
#' }
ancpar<-function(x1,y1,x2,y2,pts=NULL,regfun=tsreg,fr1=1,fr2=1,alpha=.05,plotit=TRUE,xout=FALSE,outfun=out,nboot=100,SEED=TRUE,xlab="X",ylab="Y",...){
#
# Compare the regression lines of two independent groups at specified design points.
# By default, use the Theil--Sen estimator
#
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#  For p>1 predictors, pts should be a matrix with p columns
#
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop("x1 and x2 have different number of columns")
if(ncol(x1)==1)output=ancts(x1,y1,x2,y2,pts=pts,regfun=regfun,fr1=fr1,fr2=fr2,alpha=alpha,
plotit=plotit,xout=xout,outfun=outfun,nboot=nboot,SEED=SEED,xlab=xlab,ylab=ylab,...)
if(ncol(x1)>1)output=anctsmp(x1,y1,x2,y2,regfun=regfun,alpha=alpha,pts=pts,SEED=SEED,xout=xout,outfun=outfun,nboot=nboot,...)
output
}


 ols.coef<-function(x,y,xout=FALSE){
 # In some cases, want the OLS estimate returned in $coef
 res=ols(x,y,xout=xout)$coef[,1]
 list(coef=res)
 }



#' ANCOVA with Percentile Bootstrap and Running Interval Smooth
#'
#' Compares two independent groups using robust ANCOVA with percentile bootstrap
#' combined with running interval smoother. Design points can be specified or
#' determined automatically.
#'
#' @param x1 Covariate values for group 1.
#' @param y1 Response values for group 1.
#' @param x2 Covariate values for group 2.
#' @param y2 Response values for group 2.
#' @param est Estimator function (default: hd for Harrell-Davis).
#' @param pts Design points for comparisons. If NA, five points chosen automatically.
#' @param fr1 Span for group 1 running interval smoother (default: 1).
#' @param fr2 Span for group 2 running interval smoother (default: 1).
#' @inheritParams common-params
#' @param nmin Minimum number of values close to a design point (default: 12).
#' @param xout Logical. Remove outliers from covariates before analysis.
#' @param LP Logical. Use lowess/loess smoothing for plots (default: TRUE).
#'   Set to FALSE when using onestep or mom estimators.
#' @param xlab,ylab Axis labels for plots.
#' @param pch1,pch2 Plot characters for groups 1 and 2.
#' @param ... Additional arguments passed to outfun or est.
#'
#' @return List with bootstrap results and confidence intervals at each design point.
#'
#' @note LP is automatically set to FALSE when using onestep or mom estimators.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- 0.5*x2 + rnorm(50, mean=0.5)
#' result <- ancpb(x1, y1, x2, y2, est=hd, nboot=2000)
#' }
ancpb<-function(x1,y1,x2,y2,est=hd,pts=NA,fr1=1,fr2=1,nboot=NA,nmin=12,alpha=.05,xout=FALSE,outfun=outpro,plotit=TRUE,LP=TRUE,xlab='X',ylab='Y',pch1='*',pch2='+',...){
#
# Compare two independent  groups using an ancova method
# with a percentile bootstrap combined with a running interval
# smooth.
#
#  Assume data are in x1 y1 x2 and y2
#  Comparisons are made at the design points contained in the vector
#  pts
#
flag.est=FALSE
if(identical(est,onestep))flag.est=TRUE
if(flag.est)LP=FALSE   # Get an error when using onestep in conjunction with LP=T
if(identical(est,mom))flag.est=TRUE
xy1=elimna(cbind(x1,y1))
x1=xy1[,1]
y1=xy1[,2]
xy2=elimna(cbind(x2,y2))
x2=xy2[,1]
y2=xy2[,2]
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
npt<-5
gv1<-vector("list")
if(is.na(pts[1])){
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
isub[1]<-min(sub[vecn>=nmin])
isub[5]<-max(sub[vecn>=nmin])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
mat<-matrix(NA,5,3)
dimnames(mat)<-list(NULL,c("X","n1","n2"))
for (i in 1:5){
j<-i+5
temp1<-y1[near(x1,x1[isub[i]],fr1)]
temp2<-y2[near(x2,x1[isub[i]],fr2)]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(temp1)
mat[i,3]<-length(temp2)
gv1[[i]]<-temp1
gv1[[j]]<-temp2
}
I1<-diag(npt)
I2<-0-I1
con<-rbind(I1,I2)
if(flag.est)test<-pbmcp(gv1,alpha=alpha,nboot=nboot,est=est,con=con,...)
if(!flag.est)test<-linconpb(gv1,alpha=alpha,nboot=nboot,est=est,con=con,...)
}
#
if(!is.na(pts[1])){
npt<-length(pts)
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),3)
dimnames(mat)<-list(NULL,c("X","n1","n2"))
gv<-vector("list",2*length(pts))
for (i in 1:length(pts)){
j<-i+npt
temp1<-y1[near(x1,pts[i],fr1)]
temp2<-y2[near(x2,pts[i],fr2)]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
mat[i,1]<-pts[i]
if(length(temp1)<=5)paste("Warning, there are",length(temp1)," points corresponding to the design point X=",pts[i])
if(length(temp2)<=5)paste("Warning, there are",length(temp2)," points corresponding to the design point X=",pts[i])
mat[i,2]<-length(temp1)
mat[i,3]<-length(temp2)
gv1[[i]]<-temp1
gv1[[j]]<-temp2
}
I1<-diag(npt)
I2<-0-I1
con<-rbind(I1,I2)
if(flag.est)test<-pbmcp(gv1,alpha=alpha,nboot=nboot,est=est,con=con,...)
if(!flag.est)test<-linconpb(gv1,alpha=alpha,nboot=nboot,est=est,con=con,...)
}
if(plotit){
runmean2g(x1,y1,x2,y2,fr=fr1,est=est,LP=LP,xlab=xlab,ylab=ylab,pch1=pch1,pch2=pch2,...)
}
list(mat=mat,output=test$output,con=test$con,num.sig=test$num.sig)
}


#' Plot Effect Size as Function of Covariate
#'
#' Creates an effect size curve showing how the effect size varies across
#' covariate values. Uses running interval smoothing and optionally adds
#' confidence intervals.
#'
#' @param x1,y1 Covariate and outcome for group 1.
#' @param x2,y2 Covariate and outcome for group 2.
#' @param fr1,fr2 Spans for nearest neighbor selection (default: 1).
#' @inheritParams common-params
#' @param pts Design points where effect sizes are computed (default: x1 values).
#' @param method Effect size measure: 'QS' (quantile shift, default), 'EP'
#'   (explanatory power), 'AKP' (robust Cohen's d), 'WMW' (P(X<Y)), 'KMS'.
#' @param CI Logical. If `TRUE`, computes confidence intervals and p-values.
#' @param pts.only Logical. If `TRUE` (default), plots only the estimates.
#'   If `FALSE`, adds LOESS smoother.
#' @param low.span Span parameter for LOESS smoother (default: 2/3).
#' @param nmin Minimum number of neighbors required at each design point (default: 12).
#' @param xlab,ylab Axis labels.
#' @param pch Plotting character (default: '*').
#' @param ... Additional arguments.
#'
#' @return If `CI = FALSE`, returns 'Done'. If `CI = TRUE`, returns matrix
#'   with columns: X, Est, ci.low, ci.hi, p.value.
#'
#' @details
#' For each design point where both groups have ≥ nmin observations nearby,
#' computes an effect size using the specified method. The default 'QS'
#' (quantile shift) method provides a robust measure of effect size.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- x2 + 0.5 + rnorm(50)
#' anc.plot.es(x1, y1, x2, y2, method = 'QS', CI = TRUE)
#' }
anc.plot.es<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,alpha=.05,plotit=TRUE,pts=x1,method='QS',CI=FALSE,
pr=TRUE,xout=FALSE,outfun=out,xlab='X',ylab='Effect.Size',pch='*',pts.only=TRUE,low.span=2/3,
nmin=12,...){

#  Plot effect size curve. Done for each point in x1 for which the number of nearest neighbors for
# both x1 and x2 is > nmin
#  nmim default =12
#
#  pts.only=TRUE: plot the estimates
#  pts.only=FALSE: add a smoother to the points using LOESS
#  low.span control the span
#
#  fr1 and fr2 are the spans when looking for the nearest neighbors
#   see function near
#
if(pr){
print('Effect size is based on the argument method, default is quantile shift measure of effect size')
print('Other options: EP, explanatory power; AKP, robust analog of Cohen d; WMW, P(X<Y)')
print('KMS, robust analog of Cohen d')
}

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

xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
e=NA
x=NA
ic=0
n1=NA
ci.low=NA
ci.hi=NA
pv=NA
n=length(pts)
for(i in 1:n){
ysub1=y1[near(x1,pts[i],fr1)]
ysub2=y2[near(x2,pts[i],fr2)]
n1=length(ysub1)
n2=length(ysub2)
if(n1>=nmin & n2>=nmin){
ic=ic+1
e[ic]=ESfun(ysub1,ysub2,method=method)
x[ic]=pts[i]
if(CI){
temp=ESfun.CI(ysub1,ysub2,method=method)
if(identical(method,'WMW')){
ci.low[ic]=temp$p.ci[1]
ci.hi[ic]=temp$p.ci[2]
pv[ic]=temp$p.value
}
if(!identical(method,'WMW')){
ci.low[i]=temp$ci[1]
ci.hi[i]=temp$ci[2]
pv[ic]=temp$p.value
}
}}}
if(plotit){
if(pts.only)plot(x,e,pch=pch,xlab=xlab,ylab=ylab)
else
lplot(x,e,pr=FALSE,xlab=xlab,ylab=ylab,low.span=low.span)
M='Done'
}
if(CI){
M=cbind(x,e,ci.low,ci.hi,pv)
dimnames(M)=list(NULL,c('X','Est','ci.low','ci.hi','p.value'))
}
M
}


#' ANCOVA Using Depths of Smooths
#'
#' Compares two independent group regression lines using the depths of smooths
#' rather than traditional hypothesis tests. Uses bootstrap to determine
#' significance based on smooth depth.
#'
#' @param x1,y1 Covariate and outcome for group 1.
#' @param x2,y2 Covariate and outcome for group 2.
#' @param crit.mat Optional matrix of bootstrap critical values. If `NULL`,
#'   computes via bootstrap.
#' @inheritParams common-params
#' @param REP.CRIT Logical. If `TRUE`, returns bootstrap critical values matrix.
#' @param LP Logical. Use LOESS for plot (default: `TRUE`).
#' @param est Estimator function (default: \code{tmean} for trimmed mean).
#' @param fr Span for running interval smoother. If `NULL`, chosen automatically
#'   based on sample size.
#' @param sm Logical. If `TRUE`, uses bootstrap bagging for smooths.
#' @param xlab,ylab Axis labels.
#' @param ... Additional arguments.
#'
#' @return List with:
#' \describe{
#'   \item{p.value}{Bootstrap p-value for testing equal regression lines.}
#'   \item{crit.mat}{Bootstrap critical values (if `REP.CRIT = TRUE`).}
#'   \item{test.depth}{Observed depth statistic.}
#' }
#'
#' @details
#' Uses running interval smoother to create nonparametric regression estimates.
#' Tests equality of regression lines by comparing depth of observed smooth
#' difference to bootstrap distribution under null hypothesis.
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- x2 + 0.5 + rnorm(50)
#' result <- ancsm(x1, y1, x2, y2, nboot = 500)
#' }
ancsm<-function(x1,y1,x2,y2,crit.mat=NULL,nboot=200,SEED=TRUE,REP.CRIT=FALSE,LP=TRUE,
est=tmean,fr=NULL,plotit=TRUE,sm=FALSE,xout=FALSE,outfun=out,xlab="X",ylab="Y",...){
#
# Compare two nonparametric
# regression lines corresponding to two independent groups
#  using the depths of smooths.
# One covariate only is allowed.
#
# A running interval smoother is used.
#
#  sm=T will create smooths using bootstrap bagging.
#
if(ncol(as.matrix(x1))>1)stop("One covariate only is allowed")
if(xout){
flag1=outfun(x1,...)$keep
flag2=outfun(x2,...)$keep
x1=x1[flag1]
y1=y1[flag1]
x2=x2[flag2]
y2=y2[flag2]
}
xy=elimna(cbind(x1,y1))
x1=xy[,1]
xord=order(x1)
x1=x1[xord]
y1=xy[xord,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
xord=order(x2)
x2=x2[xord]
y2=xy[xord,2]
n1=length(y1)
n2=length(y2)
if(is.null(fr)){
fr=1
if(min(n1,n2)>150)fr=.2
if(max(n1,n2)<35)fr=.5
}
if(SEED)set.seed(2)
if(is.null(crit.mat[1])){
crit.val=NA
yall=c(y1,y2)
xall=c(x1,x2)
nn=n1+n2
il=n1+1
for(i in 1:nboot){
data=sample(nn,nn,TRUE)
yy1=yall[data[1:n1]]
yy2=yall[data[il:nn]]
xx1=xall[data[1:n1]]
xx2=xall[data[il:nn]]
crit.mat[i]=depthcom(xx1,yy1,xx2,yy2,est=est,fr=fr)
}}
if(plotit)runmean2g(x1,y1,x2,y2,fr=fr,est=est,sm=sm,xlab=xlab,ylab=ylab,LP=LP,...)
dep=depthcom(x1,y1,x2,y2,est=est,fr=fr)
n=min(n1,n2)
pv=1-mean(crit.mat<dep)
if(!REP.CRIT)crit.mat=NULL
list(p.value=pv,crit.mat=crit.mat,test.depth=dep)
}



#' General ANCOVA with Running Interval Smoother
#'
#' Compares two independent groups using ANCOVA with trimmed means and running
#' interval smoother. Makes no assumptions about regression line form.
#'
#' @param x1,y1 Covariate and outcome for group 1.
#' @param x2,y2 Covariate and outcome for group 2.
#' @param pts Design points for comparisons (required). Maximum 28 points allowed.
#' @param fr1,fr2 Spans for running interval selection (default: 1).
#' @inheritParams common-params
#'
#' @return List with:
#' \describe{
#'   \item{output}{Matrix with X, n1, n2, DIF, TEST, se, ci.low, ci.hi for each point.}
#'   \item{crit}{Critical value used for confidence intervals (Rom's method).}
#' }
#'
#' @details
#' Uses Yuen's test for trimmed means at each design point. No parametric
#' regression assumption - uses running interval smoother to get local estimates.
#' Critical values adjusted for multiple comparisons via Rom's method (max 28 points).
#'
#' @export
#' @examples
#' \dontrun{
#' x1 <- rnorm(50); y1 <- x1 + rnorm(50)
#' x2 <- rnorm(50); y2 <- x2 + 0.5 + rnorm(50)
#' pts <- seq(min(c(x1,x2)), max(c(x1,x2)), length.out = 5)
#' result <- anctgen(x1, y1, x2, y2, pts)
#' }
anctgen<-function(x1,y1,x2,y2,pts,fr1=1,fr2=1,tr=.2){
#
# Compare two independent  groups using the ancova method
# in chapter 9. No assumption is made about the form of the regression
# lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#  Comparisons are made at the design points contained in the vector
#  pts
#
#  Comparisons can be made using at most 28 design points, otherwise
#  a critical value for controlling the experimentwise type I error cannot
#  be computed.
#
if(length(pts)>=29)stop("At most 28 points can be compared")
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),8)
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","se","ci.low","ci.hi"))
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test<-yuen(g1,g2,tr=tr)
mat[i,1]<-pts[i]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
mat[i,4]<-test$dif
mat[i,5]<-test$teststat
mat[i,6]<-test$se
if(length(pts)>=2)critv<-smmcrit(test$df,length(pts))
if(length(pts)==1)critv<-qt(.975,test$df)
cilow<-test$dif-critv*test$se
cihi<-test$dif+critv*test$se
mat[i,7]<-cilow
mat[i,8]<-cihi
}
list(output=mat,crit=critv)
}


# ancts
#' Theil-Sen ANCOVA for Two Independent Groups
#'
#' Compares regression lines of two independent groups at specified design points
#' using a robust regression estimator. By default, uses the Theil-Sen estimator,
#' which is highly resistant to outliers.
#'
#' @param x1 A numeric vector of covariate values for group 1.
#' @param y1 A numeric vector of outcome values for group 1.
#' @param x2 A numeric vector of covariate values for group 2.
#' @param y2 A numeric vector of outcome values for group 2.
#' @param pts A numeric vector specifying design points where regression lines
#'   are compared. If `NULL`, points are chosen automatically (see `Dpts` and `Npts`).
#' @param Dpts Logical. If `FALSE` (default), `Npts` covariate points are chosen
#'   uniformly spaced between minimum and maximum observed values. If `TRUE`,
#'   points are chosen using the same strategy as \code{\link{ancova}}.
#' @param Npts Number of design points to use when `pts = NULL` (default: 5).
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#'   Other options include \code{tshdreg} (Theil-Sen with Harrell-Davis),
#'   \code{opreg} (outlier-pruned), etc.
#' @param fr1 Span for selecting observations near design points in group 1 (default: 1).
#' @param fr2 Span for selecting observations near design points in group 2 (default: 1).
#' @param SCAT Logical. If `TRUE` (default), displays scatter plot points.
#' @param pch1 Plotting character for group 1 points (default: '*').
#' @param pch2 Plotting character for group 2 points (default: '+').
#' @inheritParams common-params
#' @param xlab Label for x-axis (default: "X").
#' @param ylab Label for y-axis (default: "Y").
#' @param ... Additional arguments passed to the regression function.
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Vector of sample sizes for the two groups.}
#'   \item{intercept.slope.group1}{Regression coefficients (intercept, slope) for group 1.}
#'   \item{intercept.slope.group2}{Regression coefficients (intercept, slope) for group 2.}
#'   \item{output}{Matrix with columns: X (design point), Est1 (group 1 estimate),
#'     Est2 (group 2 estimate), DIF (difference), TEST (test statistic), se (standard error),
#'     ci.low (lower CI), ci.hi (upper CI), p.value, p.adjust (adjusted p-value).}
#' }
#'
#' @details
#' The function compares two regression lines using robust regression at each design point.
#' The Theil-Sen estimator is robust to outliers and does not require normality assumptions.
#'
#' Standard errors are estimated via bootstrap. Critical values are adjusted for multiple
#' comparisons using Rom's method when testing at multiple design points.
#'
#' If `plotit = TRUE`, displays both groups' data and fitted regression lines (group 1: solid,
#' group 2: dashed).
#'
#' @seealso \code{\link{anctspb}} for percentile bootstrap version,
#'   \code{\link{anctsmp}} for multiple covariates, \code{\link{anctsmcp}} for
#'   multiple groups
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two groups with Theil-Sen regression
#' set.seed(123)
#' x1 <- rnorm(30)
#' y1 <- 2 + 3*x1 + rnorm(30)
#' x2 <- rnorm(30)
#' y2 <- 1 + 3*x2 + rnorm(30)
#'
#' result <- ancts(x1, y1, x2, y2)
#' print(result$output)
#' }
ancts<-function(x1,y1,x2,y2,pts=NULL,Dpts=FALSE,Npts=5,regfun=tsreg,fr1=1,fr2=1,SCAT=TRUE,pch1='*',pch2='+',
alpha=.05,plotit=TRUE,xout=FALSE,outfun=out,nboot=100,SEED=TRUE,xlab="X",ylab="Y",...){
#
# Compare the regression lines of two independent groups at specified design points
# using a robust regression estimator.
# By default, use the Theil--Sen estimator
#
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  Dpts=FALSE: Five covariate points are chosen uniformly space between the smallest and largest
#                         values observed.
# Dpts=TRUE:  Five covariate points are chosen in the same manner as done by the function ancova
#
if(identical(outfun,boxplot))stop('Use outfun=outbox')
if(SEED)set.seed(2)
FLAG=pts
xy=elimna(cbind(x1,y1))
if(ncol(xy)>2)stop("Only one covariate is allowed. Use ancJNmp")
x1=xy[,1]
y1=xy[,2]
nv1=length(y1)
xy=elimna(cbind(x2,y2))
if(ncol(xy)>2)stop("Only one covariate is allowed. Use ancJNmp")
x2=xy[,1]
y2=xy[,2]
nv2=length(y2)
if(xout){
m<-cbind(x1,y1)
p1=ncol(m)
p=p1-1
if(identical(outfun,outblp))flag=outblp(x1,y1,plotit=FALSE)$keep
else
flag<-outfun(x1,plotit=FALSE)$keep
m<-m[flag,]
x1<-m[,1:p]
y1<-m[,p1]
m<-cbind(x2,y2)
if(identical(outfun,outblp))flag=outblp(x2,y2,plotit=FALSE)$keep
else
flag<-outfun(x2,plotit=FALSE)$keep
m<-m[flag,]
x2<-m[,1:p]
y2<-m[,p1]
}
if(is.null(pts[1])){
xall=unique(c(x1,x2))
pts=seq(min(xall),max(xall),length.out=Npts)
if(Dpts){
npt<-5
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
pts=x1[isub]
}
mat<-matrix(NA,5,10)
dimnames(mat)<-list(NULL,c("X","Est1","Est2","DIF","TEST","se","ci.low","ci.hi","p.value",'adj.p.values'))
mat[,1]=pts
sqsd1=regYvar(x1,y1,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
sqsd2=regYvar(x2,y2,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
est1=regYhat(x1,y1,xr=pts,regfun=regfun,xout=FALSE,outfun=outfun,...)
est2=regYhat(x2,y2,xr=pts,regfun=regfun,xout=FALSE,outfun=outfun,...)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd1+sqsd2)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,9]=pv
crit<-smmcrit(Inf,5)
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
if(!is.null(FLAG)){
n1=NA
n2=NA
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),10)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value','p.adjust'))
mat[,1]<-pts
sqsd1=regYvar(x1,y1,regfun=regfun,pts=pts,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
sqsd2=regYvar(x2,y2,regfun=regfun,pts=pts,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
est1=regYhat(x1,y1,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
est2=regYhat(x2,y2,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd1+sqsd2)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,9]=pv
crit<-smmcrit(Inf,length(pts))
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
reg1=regfun(x1,y1,...)$coef
reg2=regfun(x2,y2,...)$coef
if(plotit){
if(xout){
if(identical(outfun,outblp))flag=outblp(x1,y1,plotit=FALSE)$keep
else
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
if(identical(outfun,outblp))flag=outblp(x2,y2,plotit=FALSE)$keep
else
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
plot(c(x1,x2),c(y1,y2),type="n",xlab=xlab,ylab=ylab)
if(SCAT){
points(x1,y1,pch=pch1)
points(x2,y2,pch=pch2)
}
abline(reg1)
abline(reg2,lty=2)
}
mat[,10]=p.adjust(mat[,9],method='hoch')
list(n=c(nv1,nv2),intercept.slope.group1=reg1,intercept.slope.group2=reg2,output=mat)
}


#' Theil-Sen ANCOVA Multiple Comparisons for J Groups
#'
#' Performs pairwise comparisons of regression lines for J independent groups
#' using robust regression (default: Theil-Sen estimator). All pairwise comparisons
#' are conducted via calls to \code{\link{ancts}}.
#'
#' @param x A list where `x[[j]]` contains the covariate values for group j.
#'   Can be a matrix if multiple covariates, but `ancts` restricts to single covariate.
#' @param y A list where `y[[j]]` contains the outcome values for group j.
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#' @inheritParams common-params
#' @param pts A numeric vector specifying design points where regression lines
#'   are compared. If `NULL`, points are chosen automatically by `ancts`.
#' @param fr1 Span for selecting observations near design points in first group of each pair.
#' @param fr2 Span for selecting observations near design points in second group of each pair.
#' @param ... Additional arguments passed to \code{\link{ancts}}.
#'
#' @return The function prints results for all pairwise comparisons but does not return
#'   a value. Each comparison shows the output from \code{\link{ancts}} including
#'   effect sizes, test statistics, confidence intervals, and p-values at each design point.
#'
#' @details
#' For J groups, performs J(J-1)/2 pairwise comparisons. Each comparison tests whether
#' the regression lines differ at specified design points using robust regression.
#'
#' Results are printed to the console for each pair, showing which groups are being
#' compared (e.g., "Group 1 Group 2") followed by the detailed output from `ancts`.
#'
#' @seealso \code{\link{ancts}}, \code{\link{anctsmp}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare three groups
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20), rnorm(20))
#' y <- list(2 + 3*x[[1]] + rnorm(20),
#'           1 + 3*x[[2]] + rnorm(20),
#'           2 + 4*x[[3]] + rnorm(20))
#'
#' anctsmcp(x, y, nboot = 500)
#' }
anctsmcp<-function(x,y,regfun=tsreg,nboot=599,alpha=0.05,pts=NULL,
SEED=TRUE,xout=FALSE,outfun=out,fr1=1,fr2=1,...){
#
# Like reg2ci only x1 etc have list mode containing
# data for J>1 groups. For all pairs of groups are compared via a
# call to ancova.
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
res=ancts(x[[j]],y[[j]],x[[k]],y[[k]],regfun=regfun,pts=pts,
nboot=nboot,alpha=alpha,fr1=fr1,fr2=fr2,
plotit=FALSE,xout=xout,outfun=outfun,WARN=FALSE,...)
print(paste('Group', j,'Group', k))
print(res)
}}}
}


#' Theil-Sen ANCOVA for Multiple Covariates
#'
#' Compares regression lines of two independent groups using robust regression
#' (default: Theil-Sen) when there are multiple covariates. Generalizes
#' \code{\link{ancts}} to handle more than one covariate.
#'
#' @param x1 A matrix containing covariate values for group 1 (one covariate per column).
#' @param y1 A numeric vector of outcome values for group 1.
#' @param x2 A matrix containing covariate values for group 2 (one covariate per column).
#' @param y2 A numeric vector of outcome values for group 2.
#' @param regfun Regression function to use (default: \code{tsreg} for Theil-Sen).
#' @inheritParams common-params
#' @param pts A matrix specifying design points (covariate combinations) where
#'   regression surfaces are compared. If `NULL`, points are chosen based on depth
#'   using \code{\link{ancdes}}.
#' @param ... Additional arguments passed to the regression function.
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Vector of sample sizes for the two groups.}
#'   \item{points}{Matrix of design points where comparisons were made.}
#'   \item{output}{Matrix with columns: Est 1 (group 1 estimate), Est 2
#'     (group 2 estimate), DIF (difference), TEST (test statistic), se (standard error),
#'     ci.low (lower CI), ci.hi (upper CI), p.value.}
#' }
#'
#' @details
#' This function extends \code{\link{ancts}} to multiple covariates. Design points
#' are selected based on data depth when `pts = NULL`, ensuring comparisons are made
#' at representative covariate combinations.
#'
#' Standard errors are estimated via bootstrap. Critical values are adjusted for
#' multiple comparisons using either Rom's method (for ≤28 tests) or Simes' method
#' (for >28 tests).
#'
#' The Theil-Sen estimator provides robust regression that is resistant to outliers
#' and does not require normality assumptions.
#'
#' @seealso \code{\link{ancts}} for single covariate, \code{\link{anctsmcp}}
#'   for multiple groups
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two groups with two covariates
#' set.seed(123)
#' x1 <- matrix(rnorm(60), ncol = 2)
#' y1 <- 2 + 3*x1[,1] + 2*x1[,2] + rnorm(30)
#' x2 <- matrix(rnorm(60), ncol = 2)
#' y2 <- 1 + 3*x2[,1] + 2*x2[,2] + rnorm(30)
#'
#' result <- anctsmp(x1, y1, x2, y2, nboot = 500)
#' print(result$output)
#' }
anctsmp<-function(x1,y1,x2,y2,regfun=tsreg,
alpha=.05,pts=NULL,SEED=TRUE,nboot=100,xout=FALSE,outfun=out,...){
#
# Compare two independent  groups using a generalization of the ancts function that
#  allows more than one covariate.
#
# Design points are chosen based on depth of points in x1 if pts=NULL
#  Assume data are in x1 y1 x2 and y2
#
if(SEED)set.seed(2) # now cov.mve always returns same result
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
#
if(xout){
m<-cbind(x1,y1)
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
nv=c(length(y1),length(y2))
if(is.null(pts[1])){
x1<-as.matrix(x1)
pts<-ancdes(x1)
}
pts<-as.matrix(pts)
ntests=nrow(pts)
mat<-matrix(NA,ntests,8)
dimnames(mat)<-list(NULL,c("Est 1", "Est 2","DIF","TEST","se","ci.low","ci.hi","p.value"))
sqsd1=regYvar(x1,y1,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
sqsd2=regYvar(x2,y2,pts=pts,regfun=regfun,nboot=nboot,SEED=FALSE,xout=FALSE,outfun=outfun,...)
#  xout=F because Leverage points have already been removed.
est1=regYhat(x1,y1,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
est2=regYhat(x2,y2,regfun=regfun,xr=pts,xout=FALSE,outfun=outfun,...)
mat[,1]=est1
mat[,2]=est2
est=est1-est2
mat[,3]=est
sd=sqrt(sqsd1+sqsd2)
mat[,5]=sd
tests=(est1-est2)/sd
mat[,4]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,8]=pv
crit=NULL
if(ntests==1)crit=qnorm(1-alpha/2)
if(length(pts)>1){
if(ntests<=28){
if(alpha==.05)crit<-smmcrit(Inf,ntests)
if(alpha==.01)crit<-smmcrit01(Inf,ntests)
}
if(ntests>28)crit=smmvalv2(dfvec=rep(Inf,nrow(pts)),alpha=alpha)
if(is.null(crit))crit=smmvalv2(dfvec=rep(Inf,nrow(pts)),alpha=alpha)
}
mat[,6]=est-crit*sd
mat[,7]=est+crit*sd
list(n=nv,points=pts,output=mat)
}


#' Theil-Sen ANCOVA with Percentile Bootstrap
#'
#' Compares regression lines of two independent groups at specified design points
#' using robust regression with percentile bootstrap confidence intervals and p-values.
#' Similar to \code{\link{ancts}} but uses percentile bootstrap, which can be more
#' accurate when there are tied values in the outcome variable.
#'
#' @param x1 A numeric vector of covariate values for group 1.
#' @param y1 A numeric vector of outcome values for group 1.
#' @param x2 A numeric vector of covariate values for group 2.
#' @param y2 A numeric vector of outcome values for group 2.
#' @param pts A numeric vector specifying design points where regression lines
#'   are compared. If `NULL`, five points are chosen using the same strategy as
#'   \code{\link{ancova}}, ensuring adequate sample sizes near each point.
#' @param regfun Regression function to use (default: \code{tshdreg} for Theil-Sen
#'   with Harrell-Davis estimator). Other options include \code{tsreg}, \code{opreg}, etc.
#' @param fr1 Span for selecting observations near design points in group 1 (default: 1).
#' @param fr2 Span for selecting observations near design points in group 2 (default: 1).
#' @inheritParams common-params
#' @param xlab Label for x-axis (default: 'X').
#' @param ylab Label for y-axis (default: 'Y').
#' @param ... Additional arguments passed to the regression function.
#'
#' @return A list with component:
#' \describe{
#'   \item{output}{Matrix with columns: X (design point), Est1 (group 1 estimate),
#'     Est2 (group 2 estimate), DIF (difference), ci.low (lower CI bound),
#'     ci.hi (upper CI bound), p.value.}
#' }
#'
#' @details
#' This function is similar to \code{\link{ancts}} but uses a percentile bootstrap
#' method for inference instead of estimating standard errors. The percentile bootstrap
#' can provide more accurate inferences when:
#' \itemize{
#'   \item There are tied values in the outcome variable
#'   \item The distribution of differences is skewed
#'   \item Standard error estimation is problematic
#' }
#'
#' The p-value is computed using the bootstrap distribution: P(difference = 0)
#' is estimated by the proportion of bootstrap samples where group 1 < group 2,
#' plus half the proportion where group 1 = group 2.
#'
#' If `plotit = TRUE`, displays both groups' data and fitted regression lines
#' (group 1: circles, solid line; group 2: plus signs, dashed line).
#'
#' @seealso \code{\link{ancts}} for standard version, \code{\link{anctsmp}}
#'   for multiple covariates
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare two groups with percentile bootstrap
#' set.seed(123)
#' x1 <- rnorm(30)
#' y1 <- round(2 + 3*x1 + rnorm(30))  # Rounded to create ties
#' x2 <- rnorm(30)
#' y2 <- round(1 + 3*x2 + rnorm(30))
#'
#' result <- anctspb(x1, y1, x2, y2, nboot = 1000)
#' print(result$output)
#' }
anctspb<-function(x1,y1,x2,y2,pts=NULL,regfun=tshdreg,fr1=1,fr2=1,alpha=.05,plotit=TRUE,xout=FALSE,outfun=outpro,nboot=500,SEED=TRUE,xlab='X',ylab='Y',...){
#
# Compare the regression lines of two independent groups
# at specified design points using a robust regression estimator.
#
#  Like ancts but uses
#   a percentile bootstrap method is used.
# This might help when there are tied values among the dependent variable.
#
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
if(SEED)set.seed(2)
FLAG=pts
xy=elimna(cbind(x1,y1))
if(ncol(xy)>2)stop('Only one covariate is allowed')
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
if(ncol(xy)>2)stop('Only one covariate is allowed')
x2=xy[,1]
y2=xy[,2]
if(is.null(pts[1])){
npt<-5
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
mat<-matrix(NA,length(pts),9)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value'))
pts=x1[isub]
}
mat<-matrix(NA,length(pts),7)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','ci.low','ci.hi','p.value'))
mat[,1]<-pts
bvec1=matrix(NA,nrow=nboot,ncol=length(pts))
bvec2=matrix(NA,nrow=nboot,ncol=length(pts))
x1=as.matrix(x1)
x2=as.matrix(x2)
p1=ncol(x1)+1
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
for(ib in 1:nboot){
bvec1[ib,]=regYsub(x1[data[ib,],],y1[data[ib,]],pts,p1=p1,regfun=regfun,...)
}
data<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
for(ib in 1:nboot){
bvec2[ib,]=regYsub(x2[data[ib,],],y2[data[ib,]],pts,p1=p1,regfun=regfun,...)
}
dif=bvec1<bvec2
L=apply(dif,2,mean)
E=bvec1==bvec2
T=apply(E,2,mean)
pvec=L+.5*T
pop=1-pvec
pb=cbind(pvec,pop)
pv=2*apply(pb,1,min)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ciL=NA
ciU=NA
difb=bvec1-bvec2
for(i in 1:length(pts)){
bs=sort(difb[,i])
ciL[i]=bs[ilow]
ciU[i]=bs[ihi]
}
est1=regYhat(x1,y1,xr=pts,xout=xout,outfun=outfun,regfun=regfun,...)
est2=regYhat(x2,y2,xr=pts,xout=xout,outfun=outfun,regfun=regfun,...)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
mat[,5]=ciL
mat[,6]=ciU
mat[,7]=pv
if(!is.null(FLAG[1])){
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
est1=regYhat(x1,y1,regfun=regfun,xr=pts,xout=xout,outfun=outfun,...)
est2=regYhat(x2,y2,regfun=regfun,xr=pts,xout=xout,outfun=outfun,...)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
mat[,7]=pv
}
if(plotit){
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
plot(c(x1,x2),c(y1,y2),type='n',xlab=xlab,ylab=ylab)
points(x1,y1,pch='o')
points(x2,y2,pch='+')
abline(regfun(x1,y1)$coef)
abline(regfun(x2,y2)$coef,lty=2)
}
list(output=mat)
}

# CLASSanc
CLASSanc<-function(x1,y1,x2,y2,xout=FALSE,outfun=out,...){
#
#  Perform classic ANCOVA
#  x1 and x2 contain covariate
#  NOT RECOMMENDED. ONLY INCLUDED IN CASE
#  YOU WANT TO COMPARE THE RESULTS WITH A ROBUST METHOD
#
#  ONE COVARIATE ONLY
#
x1=as.matrix(x1)
if(ncol(x1)!=1)stop("This function allows one covariate only")
if(xout){
flag=outfun(x1,plotit=FALSE)$keep
x1=x1[flag]
y1=y1[flag]
flag=outfun(x2,plotit=FALSE)$keep
x2=x2[flag]
y2=y2[flag]
}
x=c(x1,x2)
y=c(y1,y2)
g=c(rep(1,length(y1)),rep(2,length(y2)))
model=lm(y~as.factor(g)*x)
res1=summary.aov(model)
model=lm(y~as.factor(g)+x)
res2=summary.aov(model)
list(slope.test=res1,ancova=res2)
}


# DancCR
DancCR<-function(x1,y1,x2,y2){
v=optim(0,Dancols_sub1,x1=x1,y1=y1,x2=x2,y2=y2,method='BFGS')$par
v[2]=optim(0,Dancols_sub2,x1=x1,y1=y1,x2=x2,y2=y2,method='BFGS')$par
a=min(v)
v=c(a,max(v))
}


# Dancdet
Dancdet<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,DIF=TRUE,
alpha=.05,plotit=TRUE,plot.dif=FALSE,pts=NA,sm=FALSE,
pr=TRUE,xout=FALSE,outfun=out,MC=FALSE,
npts=25,p.crit=NULL,nreps=2000,SEED=TRUE,
SCAT=TRUE,xlab='X',ylab='Y',pch1='*',pch2='+',...){
#
#  ANCOVA for dependent groups.
#
#  Like Dancova, but a more detailed analysis
#  plot.dif=TRUE:  plot difference in the estimates plus a
# confidence band having simultaneous probability coverate 1-alpha
#
# npts = number of  covariate values to be used
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x1)!=length(x2))stop('x1 and y2 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
if(length(y1)!=length(y2))stop('y1 and y2 have different lengths')
xy=elimna(cbind(x1,y1,x2,y2))
x1=xy[,1]
y1=xy[,2]
x2=xy[,3]
y2=xy[,4]
#if(is.null(p.crit))set.seed(2)
if(alpha!=.05)p.crit=ancdet.pv(length(y1),length(y2),nreps=nreps,tr=tr,npts=npts,MC=MC)
else{
if(n>800)p.crit=0.00335002
if(n <= 800){
nv=c(30,  50,  60,  70,  80, 100,
150, 200, 300, 400, 500, 600, 800)
pv=c(0.00824497,0.00581, 0.005435089, 0.004763079,
0.00416832, 0.004406774, 0.00388228,0.003812836,0.003812836,0.003453055, 0.003625061,
.003372966, 0.003350022)
n1= length(y1)
 n2=length(y2)
p.crit=(lplot.pred(1/nv,pv,1/n1)$yhat+lplot.pred(1/nv,pv,1/n2)$yhat)/2
}}

res1=ancova(x1,y1,x2,y2,pr=FALSE,plotit=FALSE)$output # Get lowest and covariate
#   values where comparisons can be made.
pts=seq(res1[1,1],res1[5,1],length.out=npts)
#if(is.null(p.crit))p.crit=ancdet.pv(length(y1),length(y2),nreps=nreps,tr=tr,npts=npts,MC=MC)
if(plot.dif)plotit=FALSE
res=Dancova(x1,y1,x2,y2,fr1=fr1,fr2=fr2,tr=tr,alpha=p.crit,
DIF=DIF,plotit=plotit,pts=pts,SCAT=SCAT)$output
if(plot.dif){
yhat=plot(c(res[,1],res[,1],res[,1]),c(res[,3],res[,6],res[,7]),type='n',xlab=xlab,ylab=ylab)
z1=lplot(res[,1],res[,3],plotit=FALSE,pyhat=T)$yhat
z2=lplot(res[,1],res[,6],plotit=FALSE,pyhat=T)$yhat
z3=lplot(res[,1],res[,7],plotit=FALSE,pyhat=T)$yhat
lines(res[,1],z1)
lines(res[,1],z2,lty=2)
lines(res[,1],z3,lty=2)
}
sig=rep(0,nrow(res))
sig[res[,8]<=p.crit]=1
sig=as.matrix(sig)
res=res[,-9]
dimnames(sig)=list(NULL,'Sig.Dif')
res=cbind(res,sig)
list(p.crit=p.crit,output=res,num.sig=sum(sig),p.crit=p.crit)
}


# DancGLOB_pv
DancGLOB_pv<-function(n,est=tmean,fr1=1,fr2=1,nboot=500,
SEED=TRUE,iter=1000,nmin=12,MC=TRUE,alpha=.05,PRM=FALSE,qvals=c(.25,.5,.75),
cpp=FALSE,...){
#
#  Determine critical p-value when using the function DancovaV2  (formerly DancGLOB).
#  Strategy: generate data from a normal distribution, NULL true
#  compute minimun p-value, repeat
#  iter times (iter=1000 is default)
#
# qvals indicates the quantiles associated with the first group that will be used
# as the covariate values where the regression lines are to be compared
#
# returns:
# p.crit, the critical p-value for the specified alpha value
# if PRM=T, all p-values that were computed.
# ef.iter, the actual number of iterations, which might differ from iter
# due to sample sizes where it makes no sense to compute a p-value
# based on the  generated data.
#
if(!cpp){
if(SEED)set.seed(45)
bvec=list()
for(i in 1:iter)bvec[[i]]=rmul(n,p=4)
prm=NA
pv=lapply(bvec,DancGLOBv2,est=est,fr1=fr1,fr2=fr2,nboot=nboot,SEED=FALSE,nmin=nmin,MC=MC,qvals=qvals,plotit=FALSE,...)
pv=as.vector(matl(pv))
p=hd(pv,q=alpha)
}
if(cpp){
library(WRScpp)
p=DancGLOB_pv_C(n=n,est=est,fr1=fr1,fr2=fr2,SEED=SEED,nboot=nboot,qvals=qvals,
nmin=nmin,alpha=alpha,PRM=PRM,...)$p.crit
}
list(p.crit=p)
}



# DancGLOB_sub
DancGLOB_sub<-function(data,xy=xy,pts=pts,est=est,fr=fr,nmin=nmin,...){
x1=xy[data,1]
y1=xy[data,2]
xye=elimna(cbind(x1,y1))
est1=runhat(xye[,1],xye[,2],pts=pts,est=est,fr=fr,nmin=nmin,...)
est1
}



# DancGLOB_sub
DancGLOB_sub<-function(data,xy=NULL,pts=pts,est=est,fr=fr,nmin=nmin,...){
x1=xy[data,1]
y1=xy[data,2]
xye=elimna(cbind(x1,y1))
est1=runhat(xye[,1],xye[,2],pts=pts,est=est,fr=fr,nmin=nmin,...)
est1
}



# DancGLOBv2
DancGLOBv2<-function(xy=NULL,x1=NULL,y1=NULL,x2=NULL,y2=NULL,fr1=1,fr2=1,
est=tmean,alpha=.05,plotit=TRUE,pts=NULL,qvals=c(.25,.5,.75),sm=FALSE,
xout=FALSE,outfun=out,DIF=FALSE,LP=TRUE,nboot=500,SEED=TRUE,
nmin=12,MC=FALSE,...){
#
# Compare two dependent  groups using the ancova method.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
# One covariate only is allowed.
#
#  Assume data are in xy having four columns: x1, y1, x2 and y2. Or can have the
#  data stored in four separate variables:
#   x1 y1 x2 and y2
#  either matrices or in list mode.
#
#   x1 y1 are measures at time 1
#   x2 y2 are measures at time 2
#
#  LP=T, when plotting, running interval smoother is smoothed again using lplot.
#  sm=T will create smooths using bootstrap bagging.
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  DIF=FALSE: marginal trimmed means are compared
#  DIF=TRUE: Trimmed means of difference scores are used.
#
if(!is.null(x1[1])){
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x1)!=length(x2))stop('x1 and y2 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
if(length(y1)!=length(y2))stop('y1 and y2 have different lengths')
xy=elimna(cbind(x1,y1,x2,y2))
}
if(is.null(pts)){
for(i in 1:length(qvals))pts=c(pts,qest(xy[,1],qvals[i]))
}
if(SEED)set.seed(2)
n=nrow(xy)
est1=NA
est2=NA
J=length(pts)
est1=matrix(NA,nrow=nboot,ncol=J)
est2=matrix(NA,nrow=nboot,ncol=J)

data=matrix(sample(n,size=n*nboot,replace=TRUE),ncol=nboot,nrow=n)
if(!MC){
est1=apply(data,2,DancGLOB_sub,xy=xy[,1:2],pts=pts,est=est,fr=fr1,nmin=nmin,...)
est2=apply(data,2,DancGLOB_sub,xy=xy[,3:4],pts=pts,est=est,fr=fr2,nmin=nmin,...)
est1=t(as.matrix(est1))
est2=t(as.matrix(est2))
}

if(MC){
data=listm(data)
est1=mclapply(data,DancGLOB_sub,xy=xy[,1:2],pts=pts,est=est,fr=fr1,nmin=nmin,...)
est2=mclapply(data,DancGLOB_sub,xy=xy[,3:4],pts=pts,est=est,fr=fr2,nmin=nmin,...)
est1=t(matl(est1))
est2=t(matl(est2))
}

e1=runhat(xy[,1],xy[,2],pts=pts,est=est,fr=fr1,...)
e2=runhat(xy[,3],xy[,4],pts=pts,est=est,fr=fr2,...)
dif=e1-e2

pv=NA
for(j in 1:J){
pv[j]=mean(est1[,j]<est2[,j],na.rm=TRUE)+.5*mean(est1[,j]==est2[,j],na.rm=TRUE)
pv[j]=2*min(c(pv[j],1-pv[j]))
}
if(plotit){
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
runmean2g(xy[,1],xy[,2],xy[,3],xy[,4],fr=fr1,est=tmean,sm=sm,xout=FALSE,...)
}
pv=min(pv)
pv
}


# Danc.grid
Danc.grid<-function(x,y1,y2, alpha=.05,DIF=TRUE,METHOD='TR',AUTO=TRUE,PVSD=FALSE,
Qsplit1=.5,Qsplit2=.5, SV1=NULL,SV2=NULL,
tr=.2,PB=FALSE,est=tmean,nboot=1000,
xout=FALSE,outfun=outpro,SEED=TRUE,...){
#
# Two covariates. That is, x is a   matrix with two columns
#
# Two dependent groups.
# Split on that covariate variables based on data in x.
#  This is done for each covariate. Compare the corresponding the groups for each region.
#
#
# Qsplit: split on the two covariate variables based on the
#         quantiles indicated by Qsplit
#  Example
#   Qsplit1=c(.25,.5,.75)
#   Qsplit2=.5
#   would split based on the quartiles for the first independent variable and the median
#   for the second independent variable
#
#  Alternatively, the data can be split based in values stored in the arguments
#  SV1 and SV2.
#
p=ncol(x)
if(p==1)stop('There should be two  independent variables')
p1=p+1
xy1<-elimna(cbind(x,y1,y2))
x<-xy1[,1:2]
y1<-xy1[,3]
y2<-xy1[,4]
ES=list()
if(xout){
flag<-outfun(2,plotit=FALSE,...)$keep
x<-x[flag,1:2]
y1<-y1[flag]
y2<-y2[flag]
}

IZ=matrix(c(1,2,3,6,
1,2,3,4,
1,3,4,5,
1,2,3,4,
1,2,3,6),byrow=TRUE,nrow=5)
LAB=c('TR','TRPB','MED','AD','SIGN')


J=length(Qsplit1)+1
K=length(Qsplit2)+1
if(!is.null(SV1))J=length(SV1)+1
if(!is.null(SV2))K=length(SV2)+1
JK=J*K
MAT=matrix(1:JK,J,K,byrow=TRUE)
z=list()
group=list()
N.int=J
N.int2=K
NG=N.int*N.int2
if(!DIF)GRID=matrix(NA,NG,8)
if(DIF)GRID=matrix(NA,NG,6)
GI=matrix(NA,NG,4)  # grid intervals
L1=NULL
L2=NULL
qv=quantile(x[,1],Qsplit1)
if(!is.null(SV1))qv=SV1
qv=c(min(x[,1]),qv,max(x[,1]))
qv2=quantile(x[,2],Qsplit2)
if(!is.null(SV2))qv2=SV2
qv2=c(min(x[,2]),qv2,max(x[,2]))
ic=0
for(j in 1:N.int){
j1=j+1
xsub1.1=binmat(xy1,1,qv[j],qv[j1])  # split, covariate 1
for(k in 1:N.int2){
k1=k+1
xsub2.1=binmat(xsub1.1,2,qv2[k],qv2[k1])   # now split, using covariate 2
ic=ic+1
if(length(xsub2.1[,3])<=7)print('Not enough data in one or more  grids')
GI[ic,]=c(qv[j],qv[j1],qv2[k],qv2[k1])
if(!DIF){
if(length(xsub2.1[,3])>7){
a=yuend(xsub2.1[,3],xsub2.1[,4],tr=tr,alpha=alpha)
GRID[ic,1]=a$n
GRID[ic,2]=a$est1
GRID[ic,3]=a$est2
GRID[ic,4]=a$dif
GRID[ic,5]=a$ci[1]
GRID[ic,6]=a$ci[2]
GRID[ic,7]=a$p.value
if(PB){
pbv=trimpb2(xsub2.1[,3],xsub2.1[,4],tr=tr,alpha=alpha,nboot=nboot)
GRID[ic,5]=pbv$ci[1]
GRID[ic,6]=pbv$ci[2]
GRID[ic,7]=pbv$p.value
}
}}
if(DIF){
if(length(xsub2.1[,3])>7){
d=dep.dif.fun(xsub2.1[,3],xsub2.1[,4],tr=tr,alpha=alpha,method=METHOD,AUTO=AUTO,PVSD=PVSD)
J=which(LAB==METHOD)
d=pool.a.list(d)
GRID[ic,1]=length(xsub2.1[,3])
GRID[ic,2:5]=d[IZ[J,]]
}}
ES[[ic]]=dep.ES.summary.CI(xsub2.1[,3],xsub2.1[,4],tr=tr)
}}
dimnames(GI)=list(NULL,c('Int.1.low','Int.1.up','Int.2.low','Int.2.up'))
if(DIF){
dimnames(GRID)=list(NULL,c('n','Est','ci.low','ci.up','p.value','adj.p.value'))
GRID[,6]=p.adjust(GRID[,5],method='hoch')
}
if(!DIF){
dimnames(GRID)=list(NULL,c('n','est.1','est.2','DIF','ci.low','ci.up','p.value','adj.p.value'))
GRID[,8]=p.adjust(GRID[,7],method='hoch')
}
list(GRID.INTERVALS=GI,GRID=GRID, Effect.Sizes=ES)
}


# Dancols
Dancols<-function(x1,y1,x2,y2,pts=NULL,fr1=1,fr2=1,alpha=.05,plotit=TRUE,xout=FALSE,outfun=out,nboot=100,SEED=TRUE,xlab='X',ylab='Y',CR=FALSE,...){
#
# Compare the OLS regression lines of two dependent (within) groups
# at specified design points
#
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#  If not specified, points are chosen for you.
#
#  CR=TRUE: determine interval outside of which the lines cross.
#  (Analog of Johnson--Neyman method)
#
#   OUTPUT:
#  cross.interval indicates interval outside of which the lines have crossed.
#  output cr.quant.grp1 indicates that quantiles of group 1 corresponding to
#  to the end of the intervals returned in cross.interval
#
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(SEED)set.seed(2)
FLAG=pts
X=elimna(cbind(x1,y1,x2,y2))
if(ncol(X)>4)stop('Only one covariate is allowed')
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
flag1=outfun(x1,SEED=SEED,...)$out.id
flag2=outfun(x2,SEED=SEED,...)$out.id
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
n.keep=length(y1)
if(is.null(pts[1])){
npt<-5
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
pts=x1[isub]
pts=unique(pts)
npt=nrow(as.matrix(pts))
mat<-matrix(NA,npt,9)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value'))
mat[,1]=pts
sqsd=difregYvar(x1,y1,x2,y2,regfun=lsfit,pts=pts,nboot=nboot,SEED=SEED)
est1=regYhat(x1,y1,xr=pts,regfun=lsfit) #Note: if xout=T, leverage points already removed
est2=regYhat(x2,y2,xr=pts,regfun=lsfit)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
df=length(y1)-1
pv=2*(1-pt(abs(tests),df))
mat[,9]=pv
crit<-smmcrit(df,5)
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
if(!is.null(FLAG)){
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
pts=unique(pts)
mat<-matrix(NA,length(pts),9)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value'))
mat[,1]<-pts
sqsd=difregYvar(x1,y1,x2,y2,regfun=lsfit,pts=pts,nboot=nboot,SEED=SEED)
est1=regYhat(x1,y1,xr=pts,regfun=lsfit,,...)
est2=regYhat(x2,y2,xr=pts,regfun=lsfit,,...)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
df=length(y1)-1
pv=2*(1-pt(abs(tests),df))
mat[,9]=pv
crit<-smmcrit(df,length(pts))
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
if(plotit){
plot(c(x1,x2),c(y1,y2),type='n',xlab=xlab,ylab=ylab)
points(x1,y1,pch='o')
points(x2,y2,pch='+')
abline(lsfit(x1,y1)$coef)
abline(lsfit(x2,y2)$coef,lty=2)
}
int=NULL
crq=NULL
crq2=NULL
if(CR){
if(ncol(as.matrix(x1))>1)stop('CR=T only allowed with one covariate')
int=DancCR(x1,y1,x2,y2)
crq=mean(x1<=int[1])
crq[2]=mean(x1<=int[2])
crq2=mean(x2<=int[1])
crq2[2]=mean(x2<=int[2])
}

list(n=n,n.keep=n.keep,output=mat,cross.interval=int,cr.quant.grp1=crq,
cr.quant.grp2=crq2)
}

# Dancols_sub
Dancols_sub<-function(x1,y1,x2,y2,pts=NULL,fr1=1,fr2=1,alpha=.05,plotit=FALSE,xout=FALSE,outfun=out,nboot=100,SEED=TRUE,xlab='X',ylab='Y',...){
#
# Compare the OLS regression lines of two dependent (within) groups
# at specified design points
#
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#  If not specified, points are chosen for you.
#
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(SEED)set.seed(2)
FLAG=pts
X=elimna(cbind(x1,y1,x2,y2))
if(ncol(X)>4)stop('Only one covariate is allowed')
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
flag1=outfun(x1,SEED=SEED,...)$out.id
flag2=outfun(x2,SEED=SEED,...)$out.id
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
n.keep=length(y1)
if(is.null(pts[1])){
npt<-5
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
mat<-matrix(NA,5,9)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value'))
pts=x1[isub]
mat[,1]=pts
sqsd=difregYvar(x1,y1,x2,y2,regfun=lsfit,pts=pts,nboot=nboot,SEED=SEED)
est1=regYhat(x1,y1,xr=pts,regfun=lsfit) #Note: if xout=T, leverage points already removed
est2=regYhat(x2,y2,xr=pts,regfun=lsfit)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
df=length(y1)-1
pv=2*(1-pt(abs(tests),df))
mat[,9]=pv
crit<-smmcrit(df,5)
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
if(!is.null(FLAG)){
n1=1
n2=1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),9)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value'))
mat[,1]<-pts
sqsd=difregYvar(x1,y1,x2,y2,regfun=lsfit,pts=pts,nboot=nboot,SEED=SEED)
est1=regYhat(x1,y1,xr=pts,regfun=lsfit,,...)
est2=regYhat(x2,y2,xr=pts,regfun=lsfit,,...)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
df=length(y1)-1
pv=2*(1-pt(abs(tests),df))
mat[,9]=pv
crit<-smmcrit(df,length(pts))
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
if(plotit){
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
plot(c(x1,x2),c(y1,y2),type='n',xlab=xlab,ylab=ylab)
points(x1,y1,pch='o')
points(x2,y2,pch='+')
abline(lsfit(x1,y1)$coef)
abline(lsfit(x2,y2)$coef,lty=2)
}
list(n=n,n.keep=n.keep,output=mat)
}

# Dancols_sub1
Dancols_sub1<-function(pts,x1,y1,x2,y2){
#
#
ci=abs(Dancols_sub(x1,y1,x2,y2,pts=pts)$output[1,7])
ci
}

# Dancols_sub2
Dancols_sub2<-function(pts,x1,y1,x2,y2){
#
#
ci=abs(Dancols_sub(x1,y1,x2,y2,pts=pts)$output[1,8])
ci
}

# Dancova
#' ANCOVA for Dependent (Paired/Repeated Measures) Groups
#'
#' Compares two dependent (paired or repeated measures) groups using a robust ANCOVA
#' approach with running interval smoother. No parametric assumptions about regression
#' line form are required.
#'
#' @param x1 Covariate values for condition 1 (single covariate only).
#' @param y1 Response values for condition 1.
#' @param x2 Covariate values for condition 2 (default: same as x1).
#' @param y2 Response values for condition 2 (paired with y1).
#' @param fr1 Span for condition 1 running interval smoother (default: 1).
#' @param fr2 Span for condition 2 running interval smoother (default: 1).
#' @inheritParams common-params
#' @param pts Design points where regression lines are compared. If NA, five points chosen automatically.
#' @param sm Logical. If TRUE, create smooths using bootstrap bagging.
#' @param xout Logical. Remove outliers from covariates before analysis.
#' @param DIF Logical. If FALSE, marginal trimmed means are compared.
#'   If TRUE, trimmed means of difference scores are used.
#' @param LP Logical. Use lowess/loess smoothing for plots.
#' @param xlab,ylab Axis labels for plots.
#' @param ... Additional arguments passed to outfun.
#'
#' @return Matrix with columns: X (design points), n (sample size), DIF (difference),
#'   TEST (test statistic), se (standard error), ci.low, ci.hi, p.value, p.adjust.
#'   Also creates plot if plotit=TRUE.
#'
#' @note Data must be paired (same length for all vectors). One covariate only.
#'
#' @export
#' @examples
#' \dontrun{
#' # Paired data (e.g., pre-post with covariate)
#' x <- rnorm(50)  # Same covariate for both conditions
#' y1 <- x + rnorm(50)  # Condition 1
#' y2 <- x + rnorm(50, mean=0.5)  # Condition 2
#' result <- Dancova(x, y1, x, y2)
#' }
Dancova<-function(x1,y1,x2=x1,y2,fr1=1,fr2=1,tr=.2,alpha=.05,plotit=TRUE,pts=NA,
sm=FALSE,xout=FALSE,outfun=out,DIF=FALSE,LP=FALSE,xlab='X',ylab='Y',...){
#
# Compare two dependent  groups using a method  similar to the one used by the R function ancova
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#
#  sm=T will create smooths using bootstrap bagging.
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  DIF=FALSE: marginal trimmed means are compared
#  DIF=TRUE: Trimmed means of difference scores are used.
#
if(!is.null(x2)){
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
if(length(x1)!=length(x2))stop('x1 and y2 have different lengths')

if(length(y1)!=length(y2))stop('y1 and y2 have different lengths')
xy=elimna(cbind(x1,y1,x2,y2))
x1=xy[,1]
y1=xy[,2]
x2=xy[,3]
y2=xy[,4]
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
if(!is.na(pts[1]))mat=Dancovapts(x1,y1,x2,y2,fr1=fr1,fr2=fr2,tr=tr,alpha=alpha,
plotit=FALSE,pts=pts,sm=sm,xout=xout,outfun=outfun,DIF=DIF,...)
if(is.na(pts[1])){
npt<-5
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
n=length(y1)
ivals=c(1:n)
for(i in 1:length(x1))n1[i]<-length(y1[near(x1,x1[i],fr1)])
for(i in 1:length(x1))n2[i]<-length(y2[near(x2,x1[i],fr2)])
for(i in 1:length(x1))vecn[i]<-min(n1[i],n2[i])
sub<-c(1:length(x1))
isub[1]<-min(sub[vecn>=12])
isub[5]<-max(sub[vecn>=12])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
mat<-matrix(NA,5,9)
dimnames(mat)<-list(NULL,c('X','n','DIF','TEST','se','ci.low','ci.hi','p.value','p.adjust'))
for (i in 1:5){
t1=near(x1,x1[isub[i]],fr1)
t2=near(x2,x1[isub[i]],fr2)
iv1=ivals[t1]
iv2=ivals[t2]
pick=unique(c(iv1,iv2))
mat[i,2]<-length(y1[pick])
if(!DIF)test<-yuend(y1[pick],y2[pick],tr=tr)
if(DIF)test<-trimci(y1[pick]-y2[pick],tr=tr,pr=FALSE)
mat[i,1]<-x1[isub[i]]
if(!DIF){
mat[i,4]<-test$teststat
mat[i,3]<-test$dif
}
if(DIF){
mat[i,4]<-test$test.stat
mat[i,3]<-test$estimate
}
mat[i,5]<-test$se
mat[i,6]<-test$ci[1]
mat[i,7]<-test$ci[2]
mat[i,8]<-test$p.value
}
temp2<-order(0-mat[,8])
bot=c(1:nrow(mat))
dvec=sort(alpha/bot,decreasing=TRUE)
#mat[temp2,9]=dvec
mat[,9]=p.adjust(mat[,8],method='hoch')
}
if(plotit){
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
runmean2g(x1,y1,x2,y2,fr=fr1,est=tmean,sm=sm,xout=FALSE,LP=LP,xlab=xlab,ylab=ylab,...)
}}
list(output=mat)
}




#' @title Effect Size Summary for Dependent Groups ANCOVA
#'
#' @description
#' Compute multiple measures of effect size for dependent groups ANCOVA based on
#' difference scores at specified design points. Uses running interval smoother without
#' parametric assumptions about regression form.
#'
#' @inheritParams common_ancova_dep_params
#' @param REL.MAG Magnitude threshold for relative effect sizes. Default is NULL.
#' @param nboot Number of bootstrap samples for confidence intervals. Default is 1000.
#' @param ... Additional arguments passed to dep.ES.summary.CI.
#'
#' @details
#' For each design point (covariate value), computes comprehensive effect size measures
#' based on difference scores (y1 - y2). If design points are not specified, automatically
#' selects 5 points where both regression lines can be reliably estimated.
#'
#' @return List with components:
#' \describe{
#'   \item{n}{Sample sizes at each design point}
#'   \item{pts}{Design points used for comparison}
#'   \item{ES.4.Each.pt}{List of effect size summaries for each point from dep.ES.summary.CI}
#' }
#'
#' @seealso \code{\link{Dancova}}, \code{\link{dep.ES.summary}}
#'
#' @export
#' @keywords ancova dependent effect-size
#'
#' @examples
#' \dontrun{
#'   # Paired data with covariate
#'   x1 <- rnorm(30)
#'   y1 <- x1 + rnorm(30)
#'   x2 <- x1  # Same covariate (dependent design)
#'   y2 <- x1 + 0.5 + rnorm(30)
#'   result <- Dancova.ES.sum(x1, y1, x2, y2)
#'   # Shows effect sizes at 5 design points
#' }
# Dancova.ES.sum
Dancova.ES.sum<-function(x1,y1,x2=x1,y2,fr1=1,fr2=1,tr=.2,alpha=.05,pts=NA,xout=FALSE,outfun=out,
REL.MAG=NULL, SEED=TRUE,nboot=1000,...){
#
#  Compute measures of effect size based on difference scores.
#  This is done for each covariate value where the regression lines are compared as indicated the the argument
#  pts
#   This is done via the R function dep.ES.summary.CI
#
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#
A=list()
N=NA
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x1)!=length(x2))stop('x1 and y2 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
if(length(y1)!=length(y2))stop('y1 and y2 have different lengths')
xy=elimna(cbind(x1,y1,x2,y2))
x1=xy[,1]
y1=xy[,2]
x2=xy[,3]
y2=xy[,4]
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
FLAG=TRUE
n=length(y1)
ivals=c(1:n)
if(is.na(pts[1])){
FLAG=FALSE
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
for (i in 1:5){
t1=near(x1,x1[isub[i]],fr1)
t2=near(x2,x1[isub[i]],fr2)
iv1=ivals[t1]
iv2=ivals[t2]
pick=unique(c(iv1,iv2))
N[i]=length(pick)
pts[i]=x1[isub[i]]
A[[i]]=dep.ES.summary.CI(y1[pick],y2[pick], tr=tr, alpha=alpha, REL.MAG=REL.MAG, SEED=SEED,nboot=nboot)
}}
if(FLAG){
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
npts=length(pts)
for(i in 1:npts){
t1=near(x1,pts[i],fr1)
t2=near(x2,pts[i],fr2)
iv1=ivals[t1]
iv2=ivals[t2]
pick=unique(c(iv1,iv2))
N[i]=length(pick)
A[[i]]=dep.ES.summary.CI(y1[pick],y2[pick], tr=tr, alpha=alpha, REL.MAG=REL.MAG, SEED=SEED,nboot=nboot)
}}
list(n=N,pts=pts,ES.4.Each.pt=A)
}


#' @title Dependent Groups ANCOVA with Multiple Covariates and Design Points
#'
#' @description
#' Compare two dependent groups using nonparametric ANCOVA with multiple covariates.
#' Tests differences at multiple design points using running interval smoother.
#'
#' @param x1 Matrix of covariates for condition 1.
#' @param y1 Outcome values for condition 1.
#' @param x2 Matrix of covariates for condition 2. Default is NULL (uses x1).
#' @param y2 Outcome values for condition 2.
#' @inheritParams common_ancova_params
#' @param tr Amount of trimming (0 to 0.5). Default is 0.2.
#' @param alpha Significance level. Default is 0.05.
#' @param pts Design points matrix. If NULL, chosen based on depth of x1.
#' @param SEED Logical. Set random seed for reproducibility. Default is TRUE.
#' @param DIF Logical. If TRUE, analyze difference scores; if FALSE, compare marginal estimates. Default is TRUE.
#' @param cov.fun Covariance function for finding nearest points. Default is skipcov.
#'   Options: skipcov, tbscov, covogk, mgvcov, mvecov, mcdcov, wincov.
#' @param ... Additional arguments passed to cov.fun.
#'
#' @details
#' Extension of Dancova to multiple covariates. At each design point, compares groups
#' using points nearest to the design point in multivariate covariate space.
#' If DIF=TRUE, uses trimmed mean of difference scores; if DIF=FALSE, compares
#' marginal trimmed means with Yuen's test.
#'
#' @return List with components:
#' \describe{
#'   \item{pts}{Design points matrix}
#'   \item{output}{Matrix with test results for each design point (n, estimates, test statistics, CI, p-values)}
#' }
#'
#' @seealso \code{\link{Dancova}}, \code{\link{ancovamp}}
#'
#' @export
#' @keywords ancova dependent multivariate
#'
#' @examples
#' \dontrun{
#'   x1 <- matrix(rnorm(100), 50, 2)
#'   y1 <- rnorm(50)
#'   y2 <- y1 + 0.5 + rnorm(50)  # Paired outcomes
#'   result <- Dancovamp(x1, y1, x1, y2)  # x2=x1 for dependent design
#' }
# Dancovamp
Dancovamp<-function(x1,y1,x2=NULL,y2,fr1=1,fr2=1,tr=0.2,alpha=0.05, pts=NULL,SEED=TRUE,DIF=TRUE,cov.fun=skipcov,...){
#
# Compare two dependent  groups using a nonparametric ANCOVA method.
# Multiple covariates are allowed.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
# Design points are chosen based on depth of points in x1 if pts=NULL
#  Assume data are in x1 y1 x2 and y2
#
#  Choices for cov.fun include
# skipcov
# tbscov
# covogk
# mgvcov
# mvecov
# mcdcov
# wincov
#
if(is.null(x2))x2=x1
flag=identical(cov.fun,cov.mve)
if(flag)if(SEED)set.seed(2) # now cov.mve always returns same result
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 should have same number of columns')
if(ncol(x1)==1)stop('For one covariate, use Dancova')
if(nrow(x1)!=nrow(x2))stop('x1 and x2 should have same number of rows')
if(length(y1)!=length(y2))stop('y1 and y2 should have same length')
p=ncol(x1)
p1=p+1
m1=elimna(cbind(x1,y1,x2,y2))
x1=m1[,1:p]
y1=m1[,p1]
p2=p1+1
p3=p2+p-1
p4=p3+1
x2=m1[,p2:p3]
y2=m1[,p4]
if(is.null(pts[1])){
x1<-as.matrix(x1)
x2<-as.matrix(x2)
pts<-ancdes(x1)
}
pts<-as.matrix(pts)
flag<-rep(TRUE,nrow(pts))
if(!DIF){
mat<-matrix(NA,nrow(pts),10)
dimnames(mat)<-list(NULL,c('n','est1','est2','DIF','TEST','se','ci.low','ci.hi','p.value','p.adj'))
}
if(DIF){
mat<-matrix(NA,nrow(pts),8)
dimnames(mat)<-list(NULL,c('n','DIF','TEST','se','ci.low','ci.hi','p.value','p.adj'))
}
n<-1
vecn<-1
mval1<-cov.funl(cov.fun(x1,...))
mval2<-cov.funl(cov.fun(x2,...))
for(i in 1:nrow(pts)){
t1=near3d(x1,pts[i,],fr1,mval1)
t2=near3d(x2,pts[i,],fr2,mval2)
pick=as.logical(t1*t2)
n[i]<-length(y1[pick])
if(n[i]<5)flag[i]<-FALSE
if(n[i]>=5){
if(!DIF){
test<-yuend(y1[pick],y2[pick],tr=tr,alpha=alpha)
mat[i,2]=test$est1
mat[i,3]=test$est2
mat[i,4]=test$dif
mat[i,5]=test$teststat
mat[i,6]=test$se
mat[i,7]=test$ci[1]
mat[i,8]=test$ci[2]
mat[i,9]=test$p.value
}
if(DIF){
test<-trimci(y1[pick]-y2[pick],tr=tr,pr=FALSE,alpha=alpha)
mat[i,2]=test$estimate
mat[i,3]=test$test.stat
mat[i,4]=test$se
mat[i,5]=test$ci[1]
mat[i,6]=test$ci[2]
mat[i,7]=test$p.value
}
}
mat[i,1]<-n[i]
}
if(!DIF)mat[,10]=p.adjust(mat[,9],method='hoch')
if(DIF)mat[,8]=p.adjust(mat[,7],method='hoch')
if(sum(flag)==0)print('No comparable design points found, might increase span.')
list(pts=pts,output=mat)
}



#' @title Dependent Groups ANCOVA with Percentile Bootstrap
#'
#' @description
#' Compare two dependent groups using nonparametric ANCOVA with percentile bootstrap.
#' Uses running interval smoother without parametric assumptions.
#'
#' @inheritParams common_ancova_dep_params
#' @param est Estimator function. Default is hd (Harrell-Davis).
#' @param alpha Significance level. Default is 0.05.
#' @param nboot Number of bootstrap samples. Default is 500.
#' @param pr Logical. Print progress messages. Default is TRUE.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param plotit Logical. Create comparison plot. Default is TRUE.
#' @param sm Logical. Use smooth (bootstrap bagging). Default is FALSE.
#' @param xout Logical. Remove outliers before analysis. Default is FALSE.
#' @param outfun Outlier detection function. Default is out.
#' @param DIF Logical. If TRUE, analyze difference scores; if FALSE, compare marginal estimates. Default is FALSE.
#' @param na.rm Logical. Remove NA values. Default is TRUE.
#' @param ... Additional arguments.
#'
#' @return List with test results and design points.
#'
#' @seealso \code{\link{Dancova}}, \code{\link{ancovapb}}
#'
#' @export
#' @keywords ancova dependent bootstrap
#'
#' @examples
#' \dontrun{
#'   x <- rnorm(30)
#'   y1 <- x + rnorm(30)
#'   y2 <- x + 0.5 + rnorm(30)
#'   result <- Dancovapb(x, y1, x, y2)
#' }
# Dancovapb
Dancovapb<-function(x1,y1,x2=x1,y2,fr1=1,fr2=1,est=hd,alpha=.05,nboot=500,pr=TRUE,SEED=TRUE,
plotit=TRUE, pts=NA,sm=FALSE,xout=FALSE,outfun=out,DIF=FALSE,na.rm=TRUE,...){
#
# Compare two dependent  groups using the ancova method
# (a method similar to the one used by the  R function ancova).
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
# percentile bootstrap method is used.
#
# est indicates estimator to be used; Harrell-Davis median estimator is default.
#
#  Assume data are in x1 y1 x2 and y2
#
#  sm=T will create smooths using bootstrap bagging.
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#  pts=NA means five points will be picked empirically.
#
#
if(SEED)set.seed(2)
if(DIF & !na.rm){
if(pr)stop('With na.rm=TRUE, must have DIF=FALSE')
}
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x1)!=length(x2))stop('x1 and y2 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
if(length(y1)!=length(y2))stop('y1 and y2 have different lengths')
xy=elimna(cbind(x1,y1,x2,y2))
x1=xy[,1]
y1=xy[,2]
x2=xy[,3]
y2=xy[,4]
if(is.na(pts[1])){
npt<-5
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
mat<-matrix(NA,5,7)
dimnames(mat)<-list(NULL,c('X','n','DIF','ci.low','ci.hi','p.value','p.adjust'))
for (i in 1:5){
t1=near(x1,x1[isub[i]],fr1)
t2=near(x2,x1[isub[i]],fr2)
pick=as.logical(t1*t2)
if(!na.rm)test=rmmismcp(y1[pick],y2[pick],est=est,,nboot=nboot,alpha=alpha,pr=FALSE,
plotit=FALSE,SEED=FALSE,...)
if(na.rm){
test=rmmcppb(y1[pick],y2[pick],est=est,dif=DIF,nboot=nboot,plotit=FALSE,alpha=alpha,
pr=FALSE,SEED=SEED,...)
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(y1[pick])
mat[i,3]<-test$output[,2]
mat[i,4]<-test$output[,5]
mat[i,5]<-test$output[,6]
mat[i,6]<-test$output[,3]
}
if(!na.rm){
test=rmmismcp(y1[pick],y2[pick],est=est,nboot=nboot,alpha=alpha,pr=FALSE,
plotit=FALSE,SEED=SEED,...)
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(y1[pick])
mat[i,3]<-test$output[,3]
mat[i,4]<-test$output[,6]
mat[i,5]<-test$output[,7]
mat[i,6]<-test$output[,4]
}
}
temp2<-order(0-mat[,6])
bot=c(1:nrow(mat))
dvec=sort(alpha/bot,decreasing=TRUE)
#mat[temp2,7]=dvec
mat[,7]=p.adjust(mat[,6],method='hoch')
}
if(!is.na(pts[1])){
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
# First check sample size
#
flage=rep(TRUE,length(pts))
for (i in 1:length(pts)){
t1<-near(x1,pts[i],fr1)
t2<-near(x2,pts[i],fr2)
pick=as.logical(t1*t2)
if(sum(pick)<=5){print(paste('Warning: there are',sum(pick),' points corresponding to the design point X=',pts[i]))
flage[i]=FALSE
}}
pts=pts[flage]
mat<-matrix(NA,length(pts),7)
dimnames(mat)<-list(NULL,c('X','n','DIF','ci.low','ci.hi',
'p.value','p.crit'))
for (i in 1:length(pts)){
t1<-near(x1,pts[i],fr1)
t2<-near(x2,pts[i],fr2)
pick=as.logical(t1*t2)
#print(y1[pick])
test=rmmcppb(y1[pick],y2[pick],est=est,dif=DIF,plotit=FALSE,alpha=alpha,pr=FALSE,SEED=FALSE,...)
mat[i,3]<-test$output[,2]
mat[i,1]<-pts[i]
mat[i,2]<-length(y1[pick])
mat[i,4]<-test$output[,5]
mat[i,5]<-test$output[,6]
mat[i,6]<-test$output[,3]
}
#temp2<-order(0-mat[,6])
mat[,7]=p.adjust(mat[,6],method='hoch')
bot=c(1:nrow(mat))
dvec=sort(alpha/bot,decreasing=TRUE)
#mat[temp2,7]=dvec
}
if(plotit){
runmean2g(x1,y1,x2,y2,fr=fr1,est=est,sm=sm,xout=xout,outfun=outfun,,...)
}
list(output=mat)
}

#' @title Dependent Groups ANCOVA at User-Specified Design Points
#'
#' @description
#' Compare two dependent groups using nonparametric ANCOVA at user-specified design points.
#' Uses running interval smoother without parametric regression assumptions.
#'
#' @inheritParams common_ancova_dep_params
#' @inheritParams common_ancova_params
#' @param pts Design points vector where comparisons should be made.
#' @param sm Logical. Use smoothed estimates via bootstrap bagging. Default is FALSE.
#' @param DIF Logical. If TRUE, use trimmed means of difference scores; if FALSE, compare marginal trimmed means. Default is FALSE.
#' @param LP Logical. Use linear programming for control. Default is TRUE.
#' @param ... Additional arguments.
#'
#' @details
#' User specifies exact design points where regression lines are compared.
#' At each point, identifies observations from both conditions nearest to that covariate value.
#' For those matched observations, either compares marginal trimmed means (DIF=FALSE)
#' or analyzes trimmed mean of difference scores (DIF=TRUE).
#'
#' @return Matrix with test results for each design point.
#'
#' @seealso \code{\link{Dancova}}, \code{\link{ancova}}
#'
#' @export
#' @keywords ancova dependent
#'
#' @examples
#' \dontrun{
#'   x <- rnorm(30)
#'   y1 <- x + rnorm(30)
#'   y2 <- x + 0.5 + rnorm(30)
#'   result <- Dancovapts(x, y1, x, y2, pts = c(-1, 0, 1))
#' }
# Dancovapts
Dancovapts<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,alpha=.05,plotit=TRUE,pts=NA,sm=FALSE,xout=FALSE,outfun=out,DIF=FALSE,LP=TRUE,...){
#
# Compare two dependent  groups using a method  similar to the one used by the R function ancova
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#
#  sm=T will create smooths using bootstrap bagging.
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  DIF=FALSE: marginal trimmed means are compared
#  DIF=TRUE: Trimmed means of difference scores are used.
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x1)!=length(x2))stop('x1 and y2 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
if(length(y1)!=length(y2))stop('y1 and y2 have different lengths')
xy=elimna(cbind(x1,y1,x2,y2))
x1=xy[,1]
y1=xy[,2]
x2=xy[,3]
y2=xy[,4]
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
xorder<-order(x2)
y2<-y2[xorder]
x2<-x2[xorder]
n=length(y1)
npts=length(pts)
mat<-matrix(NA,nrow=npts,nco=9)
dimnames(mat)<-list(NULL,c('X','n','DIF','TEST','se','ci.low','ci.hi',
'p.value','p.crit'))
for (i in 1:npts){
t1=near(x1,pts[i],fr1)
t2=near(x2,pts[i],fr2)
ivals=c(1:n)
iv1=ivals[t1]
iv2=ivals[t2]
pick=unique(c(iv1,iv2))
mat[i,2]<-length(y1[pick])
if(!DIF)test<-yuend(y1[pick],y2[pick],tr=tr,alpha=alpha)
if(DIF)test<-trimci(y1[pick]-y2[pick],tr=tr,pr=FALSE,alpha=alpha)
mat[i,1]<-pts[i]
if(!DIF){
mat[i,4]<-test$teststat
mat[i,3]<-test$dif
}
if(DIF){
mat[i,4]<-test$test.stat
mat[i,3]<-test$estimate
}
mat[i,5]<-test$se
mat[i,6]<-test$ci[1]
mat[i,7]<-test$ci[2]
mat[i,8]<-test$p.value
}
temp2<-order(0-mat[,8])
bot=c(1:nrow(mat))
dvec=sort(alpha/bot,decreasing=TRUE)
mat[temp2,9]=dvec
mat
}

#' @title Dependent Groups ANCOVA Version 2 with Improved Bootstrap
#'
#' @description
#' Compare two dependent groups using nonparametric ANCOVA with improved bootstrap method.
#' Similar to Dancova but with enhanced power through different bootstrap resampling strategy.
#'
#' @param x1,y1,x2,y2 Covariate and outcome vectors for each condition.
#' @param xy Alternatively, data matrix with columns [x1, y1, x2, y2].
#' @inheritParams common_ancova_dep_params
#' @param est Estimator function. Default is tmean (trimmed mean).
#' @param alpha Significance level. Default is 0.05.
#' @param plotit Logical. Create diagnostic plot. Default is TRUE.
#' @param xlab,ylab Axis labels for plot.
#' @param qvals Quantiles for design point selection. Default is c(0.25, 0.5, 0.75).
#' @param sm Logical. Use smoothing. Default is FALSE.
#' @param DIF Logical. Analyze difference scores if TRUE. Default is FALSE.
#' @param LP Logical. Use linear programming. Default is TRUE.
#' @param method P-value adjustment method. Default is 'hochberg'.
#' @param nboot Number of bootstrap samples. Default is 500.
#' @param SEED,nreps,MC,cpp Bootstrap control parameters.
#' @param SCAT,pch1,pch2 Plotting parameters.
#' @param nmin Minimum sample size at design points. Default is 12.
#' @param q Quantile for certain estimators. Default is 0.5.
#' @param ... Additional arguments.
#'
#' @details
#' Like Dancova but with improved power. Bootstrap samples are obtained by resampling
#' from the full data (x1, y1, x2, y2) rather than conditioning on x values.
#' This unconditional bootstrap often provides better power.
#'
#' @return List with test results and design points.
#'
#' @seealso \code{\link{Dancova}}
#'
#' @export
#' @keywords ancova dependent bootstrap
#'
#' @examples
#' \dontrun{
#'   x <- rnorm(30)
#'   y1 <- x + rnorm(30)
#'   y2 <- x + 0.5 + rnorm(30)
#'   result <- DancovaV2(x, y1, x, y2)
#' }
# DancovaV2
DancovaV2<-function(x1=NULL,y1=NULL,x2=NULL,y2=NULL,xy=NULL,fr1=1,fr2=1,
est=tmean,alpha=.05,plotit=TRUE,xlab='X',ylab='Y',qvals=c(.25,.5,.75),sm=FALSE,
xout=FALSE,eout=FALSE,outfun=out,DIF=FALSE,LP=TRUE,method='hochberg',
nboot=500,SEED=TRUE,nreps=2000,MC=TRUE,cpp=FALSE,
SCAT=TRUE,pch1='*',pch2='+',
nmin=12,q=.5,...){
#
# Compare two dependent  groups using the ancova method.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Like Dancova, only bootstrap samples are obtained by resampling
#  from c(x1,y1,x2,y2) rather than conditioning on the x value as done by Dancova.
#   This function tends to have more power than Dancova.
#
# One covariate only is allowed.
#
# method='hochberg
# By default,  family wise error rate is controlled by Hochberg's methoe

# To get critical p-value, need the following commands to get access to the software.
# library(`devtools')
# install_github( `WRScpp', `mrxiaohe')

#  Assume data are in xy having four columns: x1, y1, x2 and y2.
#
#  Or can have the
#  data stored in four separate variables:
#   x1 y1 x2 and y2
#
#   x1 y1 are measures at time 1
#   x2 y2 are measures at time 2
#
#  LP=T, when plotting, running interval smoother is smoothed again using lplot.
#  sm=T will create smooths using bootstrap bagging.
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  q=.5 means when est=hd (Harrell-Davis estimator), median is estimated.
#
#  eout=TRUE will eliminate all outliers when plotting.
#
if(SEED)set.seed(2)
iter=nreps
if(!is.null(x1[1])){
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop('x1 and y1 have different lengths')
if(length(x1)!=length(x2))stop('x1 and y2 have different lengths')
if(length(x2)!=length(y2))stop('x2 and y2 have different lengths')
if(length(y1)!=length(y2))stop('y1 and y2 have different lengths')
xy=cbind(x1,y1,x2,y2)
}
n=nrow(elimna(xy))
if(plotit){
ef=identical(est,hd)
if(!ef)runmean2g(xy[,1],xy[,2],xy[,3],xy[,4],fr=fr1,est=est,sm=sm,xout=xout,LP=LP,eout=eout,
xlab=xlab,ylab=ylab,SCAT=SCAT,pch1=pch1,pch2=pch2,...)
if(ef)runmean2g(xy[,1],xy[,2],xy[,3],xy[,4],fr=fr1,est=hd,sm=sm,xout=xout,LP=LP,q=q,eout=eout,
xlab=xlab,ylab=ylab,SCAT=SCAT,pch1=pch1,pch2=pch2,...)
}

#eliminate this code and use Hochberg instead
#if(is.null(p.crit)){
#if(cpp)library(WRScpp)
#p.crit=DancGLOB_pv(n,fr1=fr1,fr2=fr2,nboot=nboot,est=est,SEED=SEED,iter=iter,
#nmin=nmin,MC=MC,alpha=alpha,qvals=qvals,cpp=cpp)$p.crit
#
#}

pts=NULL
#if(is.null(pts)){
for(i in 1:length(qvals))pts=c(pts,qest(xy[,1],qvals[i]))
#}
if(SEED)set.seed(2)
ef=identical(est,hd)
n=nrow(xy)
est1=NA
est2=NA
J=length(pts)
est1=matrix(NA,nrow=nboot,ncol=J)
est2=matrix(NA,nrow=nboot,ncol=J)
#
data=matrix(sample(n,size=n*nboot,replace=TRUE),ncol=nboot,nrow=n)
if(!MC){
if(!ef){
est1=apply(data,2,DancGLOB_sub,xy=xy[,1:2],pts=pts,est=est,fr=fr1,nmin=nmin,...)
est2=apply(data,2,DancGLOB_sub,xy=xy[,3:4],pts=pts,est=est,fr=fr2,nmin=nmin,...)
}
if(ef){
est1=apply(data,2,DancGLOB_sub,xy=xy[,1:2],pts=pts,est=hd,fr=fr1,nmin=nmin,q=q,...)
est2=apply(data,2,DancGLOB_sub,xy=xy[,3:4],pts=pts,est=hd,fr=fr2,nmin=nmin,q=q,...)
}
est1=t(as.matrix(est1))
est2=t(as.matrix(est2))
}
if(MC){
data=listm(data)
if(!ef){
est1=mclapply(data,DancGLOB_sub,xy=xy[,1:2],pts=pts,est=est,fr=fr1,nmin=nmin,...)
est2=mclapply(data,DancGLOB_sub,xy=xy[,3:4],pts=pts,est=est,fr=fr2,nmin=nmin,...)
}
if(ef){
est1=mclapply(data,DancGLOB_sub,xy=xy[,1:2],pts=pts,est=hd,fr=fr1,nmin=nmin,q=q,...)
est2=mclapply(data,DancGLOB_sub,xy=xy[,3:4],pts=pts,est=hd,fr=fr2,nmin=nmin,q=q,...)
}
est1=t(matl(est1))
est2=t(matl(est2))
}
pv=NA
for(j in 1:J){
pv[j]=mean(est1[,j]<est2[,j],na.rm=TRUE)+.5*mean(est1[,j]==est2[,j],na.rm=TRUE)
pv[j]=2*min(c(pv[j],1-pv[j]))
}
pad=p.adjust(pv,method=method)
pvm=cbind(pts,pv,pad)
dimnames(pvm)=list(NULL,c('X','p.values','p.adjusted'))
list(output=pvm,n=n)
}

#' @title Dependent Groups ANCOVA Using Theil-Sen Regression
#'
#' @description
#' Compare regression lines of two dependent groups using robust regression (Theil-Sen by default).
#' A robust Johnson-Neyman method for dependent groups.
#'
#' @inheritParams common_ancova_dep_params
#' @param pts Design points for comparison. If NULL, automatically selected. Default is NULL.
#' @param regfun Regression function. Default is tsreg (modified Theil-Sen).
#'   Can also use tshdreg for tied values.
#' @param alpha Significance level. Default is 0.05.
#' @param plotit Logical. Create comparison plot. Default is TRUE.
#' @param xout Logical. Remove outliers. Default is FALSE.
#' @param outfun Outlier detection function. Default is out.
#' @param nboot Number of bootstrap samples. Default is 100.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param xlab,ylab Axis labels for plot.
#' @param pr Logical. Print warning messages. Default is TRUE.
#' @param ... Additional arguments passed to regfun.
#'
#' @details
#' For OLS regression, use Dancols instead. This function uses robust regression
#' to compare dependent groups at specified design points. The default Theil-Sen
#' estimator is resistant to outliers. For data with tied values, consider regfun=tshdreg.
#'
#' @return List with test results and design points.
#'
#' @seealso \code{\link{Dancols}}, \code{\link{tsreg}}, \code{\link{tshdreg}}
#'
#' @export
#' @keywords ancova dependent robust regression
#'
#' @examples
#' \dontrun{
#'   x <- rnorm(30)
#'   y1 <- 2*x + rnorm(30)
#'   y2 <- 2*x + 0.5 + rnorm(30)  # Parallel shift
#'   result <- Dancts(x, y1, x, y2)
#' }
# Dancts
Dancts<-function(x1,y1,x2,y2,pts=NULL,regfun=tsreg,fr1=1,fr2=1,alpha=.05,plotit=TRUE,xout=FALSE,
outfun=out,nboot=100,SEED=TRUE,xlab='X',ylab='Y',pr=TRUE,...){
#
# Compare the regression lines of two dependent  groups using
# the robust regression indicated by the argument
# regfun. Default is modified Theil--Sen estimator
#
#  Comparisons are done at specified design points
#  This is a robust Johnson-Neyman method for dependent groups.
#
#  For OLS, use Dancols
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(SEED)set.seed(2)
FLAG=pts
X=elimna(cbind(x1,y1,x2,y2))
if(ncol(X)>4)stop('Only one covariate is allowed')
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
if(xout){
if(identical(outfun,outblp)){
flag1=outblp(x1,y1,plotit=FALSE)$bad.lev
flag2=outblp(x2,y2,plotit=FALSE)$bad.lev
}
else{
flag1=outfun(x1)$out.id
flag2=outfun(x2)$out.id
}
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
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
if(is.null(pts[1])){
npt<-5
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
mat<-matrix(NA,5,9)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value'))
pts=x1[isub]
mat[,1]=pts
sqsd=difregYvar(x1,y1,x2,y2,regfun=regfun,pts=pts,nboot=nboot,SEED=SEED)
est1=regYhat(x1,y1,xr=pts,regfun=regfun) #Note: if xout=T, leverage points already removed
est2=regYhat(x2,y2,xr=pts,regfun=regfun)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,9]=pv
crit<-smmcrit(Inf,5)
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
if(!is.null(FLAG)){
n1=1
n2=1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),9)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','TEST','se','ci.low','ci.hi','p.value'))
mat[,1]<-pts
sqsd=difregYvar(x1,y1,x2,y2,regfun=regfun,pts=pts,nboot=nboot,SEED=SEED)
est1=regYhat(x1,y1,xr=pts,regfun=regfun)
est2=regYhat(x2,y2,xr=pts,regfun=regfun)
mat[,2]=est1
mat[,3]=est2
est=est1-est2
mat[,4]=est
sd=sqrt(sqsd)
mat[,6]=sd
tests=(est1-est2)/sd
mat[,5]=tests
pv=2*(1-pnorm(abs(tests)))
mat[,9]=pv
crit<-smmcrit(Inf,length(pts))
mat[,7]=est-crit*sd
mat[,8]=est+crit*sd
}
if(plotit){
#if(xout){       #Leverage points already removed if xout=TRUE
#flag<-outfun(x1,...)$keep
#x1<-x1[flag]
#y1<-y1[flag]
#flag<-outfun(x2,...)$keep
#x2<-x2[flag]
#y2<-y2[flag]
#}
plot(c(x1,x2),c(y1,y2),type='n',xlab=xlab,ylab=ylab)
points(x1,y1,pch='o')
points(x2,y2,pch='+')
abline(regfun(x1,y1)$coef)
abline(regfun(x2,y2)$coef,lty=2)
}
list(output=mat)
}



#' @title Dependent Groups ANCOVA Using Theil-Sen with Percentile Bootstrap
#'
#' @description
#' Compare regression lines of two dependent groups using robust regression with
#' percentile bootstrap confidence intervals. A robust Johnson-Neyman method for dependent groups.
#'
#' @inheritParams common_ancova_dep_params
#' @param pts Design points for comparison. If NULL, automatically selected.
#' @param regfun Regression function. Default is tsreg (modified Theil-Sen).
#' @param alpha Significance level. Default is 0.05.
#' @param plotit Logical. Create comparison plot. Default is TRUE.
#' @param xout Logical. Remove outliers based on x values. Default is FALSE.
#' @param SCAT Logical. Show scatter plot. Default is TRUE.
#' @param outfun Outlier detection function. Default is outpro.
#' @param BLO Logical. Use biweight estimator for location/scale. Default is FALSE.
#' @param nboot Number of bootstrap samples. Default is 500.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param xlab,ylab Axis labels for plot.
#' @param pr Logical. Print warning messages. Default is TRUE.
#' @param eout Logical. Remove error outliers. Default is FALSE.
#' @param ... Additional arguments passed to regfun.
#'
#' @details
#' Similar to Dancts but uses percentile bootstrap for confidence intervals and hypothesis tests.
#' Provides more accurate inference when sampling distributions are skewed.
#' For OLS regression, use Dancols instead.
#'
#' @return List with test results and design points.
#'
#' @seealso \code{\link{Dancts}}, \code{\link{DanctspbMC}}, \code{\link{tsreg}}
#'
#' @export
#' @keywords ancova dependent robust bootstrap
#'
#' @examples
#' \dontrun{
#'   x <- rnorm(30)
#'   y1 <- 2*x + rnorm(30)
#'   y2 <- 2*x + 0.5 + rnorm(30)
#'   result <- Danctspb(x, y1, x, y2, nboot = 1000)
#' }
# Danctspb
Danctspb<-function(x1,y1,x2,y2,pts=NULL,regfun=tsreg,fr1=1,fr2=1,alpha=.05,plotit=TRUE,xout=FALSE,SCAT=TRUE,
outfun=outpro,BLO=FALSE,nboot=500,SEED=TRUE,xlab='X',ylab='Y',pr=TRUE,eout=FALSE,...){
#
# Compare the regression lines of two dependent  groups at specified design points using
# the robust regression estimator indicated by the argument
# regfun. Default is modified Theil--Sen estimator
#
#  Comparisons are done at specified design points
#  This is a robust Johnson-Neyman method for dependent groups.
#
#  For OLS, use Dancols
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
# Uses bootstrap samples based on resamples of the points followed by a regression fit.
# In contrast, Dancts uses bootstrap estimate of the se of Yhat followed by a pivotal test
# statistic.
#
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(SEED)set.seed(2)
FLAG=pts
X=elimna(cbind(x1,y1,x2,y2))
if(ncol(X)>4)stop('Only one covariate is allowed')
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
flagF=identical(regfun,tsreg)
if(identical(regfun,tshdreg))flagF=FALSE
if(flagF){
if(pr){
if(sum(duplicated(y1)>0))print('Duplicate values detected; tshdreg might have more power than tsreg')
pr=FALSE
}
if(pr){
if(sum(duplicated(y2)>0))print('Duplicate values detected; tshdreg might have more power than tsreg')
}}
if(is.null(pts[1])){
npt<-5
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
pts=x1[isub]
}
for(i in 1:length(pts)){
n1<-1
n2<-1
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),7)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','ci.low','ci.hi','p.value'))
mat[,1]=pts
n=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=n*nboot,replace=TRUE),nrow=nboot)
est1=apply(data,1,Danctspb.sub,x1,y1,xr=pts,regfun=regfun,xout=FALSE,...)
est2=apply(data,1,Danctspb.sub,x2,y2,xr=pts,regfun=regfun,xout=FALSE,...)
mat[,2]=regYhat(x1,y1,xr=pts,regfun=regfun,...)
mat[,3]=regYhat(x2,y2,xr=pts,regfun=regfun,...)
est=est1-est2
if(!is.matrix(est))est=matrix(est,nrow=1)
mat[,4]=mat[,2]-mat[,3]
pv1=apply(est<0,1,mean,na.rm=TRUE)
pv2=apply(est==0,1,mean,na.rm=TRUE)
pv=pv1+.5*pv2
pv1m=1-pv
pv=2*apply(cbind(pv,pv1m),1,min)
mat[,7]=pv
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
estsort=apply(est,1,sort)
mat[,5]=estsort[ilow,]
mat[,6]=estsort[ihi,]
if(plotit){
if(eout && xout)stop('Cannot have both eout and xout = F')
if(eout){
flag<-outfun(cbind(x1,y1),plotit=FALSE,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(cbind(x2,y2),plotit=FALSE,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
plot(c(x1,x2),c(y1,y2),type='n',xlab=xlab,ylab=ylab)
if(SCAT)points(x1,y1,pch='o')
if(SCAT)points(x2,y2,pch='+')
abline(regfun(x1,y1)$coef)
abline(regfun(x2,y2)$coef,lty=2)
}
list(output=mat)
}


#' @title Dependent Groups ANCOVA Using Theil-Sen with Bootstrap (Multicore)
#'
#' @description
#' Compare regression lines of two dependent groups using robust regression with
#' percentile bootstrap. Multicore-enabled version using parallel processing.
#'
#' @inheritParams common_ancova_dep_params
#' @param pts Design points for comparison. If NULL, automatically selected.
#' @param regfun Regression function. Default is tshdreg (Theil-Sen for tied values).
#' @param alpha Significance level. Default is 0.05.
#' @param SCAT Logical. Show scatter plot. Default is TRUE.
#' @param plotit Logical. Create comparison plot. Default is TRUE.
#' @param xout Logical. Remove outliers based on x values. Default is FALSE.
#' @param outfun Outlier detection function. Default is outpro.
#' @param nboot Number of bootstrap samples. Default is 500.
#' @param SEED Logical. Set random seed. Default is TRUE.
#' @param xlab,ylab Axis labels for plot.
#' @param WARN Logical. Show warnings. Default is FALSE.
#' @param pr Logical. Print progress messages. Default is TRUE.
#' @param eout Logical. Remove error outliers. Default is FALSE.
#' @param ... Additional arguments passed to regfun.
#'
#' @details
#' Similar to Danctspb but uses parallel processing for faster bootstrap computation.
#' Bootstrap samples are obtained from the full data rather than conditioning on x values,
#' similar to approach used in regci. For OLS regression, use Dancols instead.
#'
#' @return List with test results and design points.
#'
#' @seealso \code{\link{Danctspb}}, \code{\link{Dancts}}, \code{\link{tshdreg}}
#'
#' @export
#' @keywords ancova dependent robust bootstrap parallel
#'
#' @examples
#' \dontrun{
#'   x <- rnorm(30)
#'   y1 <- 2*x + rnorm(30)
#'   y2 <- 2*x + 0.5 + rnorm(30)
#'   result <- DanctspbMC(x, y1, x, y2, nboot = 1000)
#' }
# DanctspbMC
DanctspbMC<-function(x1,y1,x2,y2,pts=NULL,regfun=tshdreg,fr1=1,fr2=1,alpha=.05,SCAT=TRUE,
plotit=TRUE,xout=FALSE,outfun=outpro,nboot=500,SEED=TRUE,xlab='X',ylab='Y',WARN=FALSE,pr=TRUE,eout=FALSE,...){
#
# Compare the regression lines of two dependent  groups at specified design points using
# the robust regression estimator indicated by the argument
# regfun. Default is modified Theil--Sen estimator
#
#  Similar to Dancts, which uses a bootstrap estimate of se of Y hat
#  Here, do bootstrap based on bootstrap samples from the data
#  as done for example by regci
#
#  Comparisons are done at specified design points
#  This is a robust Johnson-Neyman method for dependent groups.
#
#  For OLS, use Dancols
#  Assume data are in x1 y1 x2 and y2
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(SEED)set.seed(2)
FLAG=pts
X=elimna(cbind(x1,y1,x2,y2))
if(ncol(X)>4)stop('Only one covariate is allowed')
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
flagF=identical(regfun,tsreg)
if(identical(regfun,tshdreg))flagF=FALSE
if(flagF){
if(pr){
if(sum(duplicated(y1)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
pr=FALSE
}
if(pr){
if(sum(duplicated(y2)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
}}
if(is.null(pts[1])){
npt<-5
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
pts=x1[isub]
}
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),7)
dimnames(mat)<-list(NULL,c('X','Est1','Est2','DIF','ci.low','ci.hi','p.value'))
mat[,1]=pts
n=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=n*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
est1=mclapply(data,Danctspb.sub,x1,y1,xr=pts,regfun=regfun,xout=FALSE,...)
est2=mclapply(data,Danctspb.sub,x2,y2,xr=pts,regfun=regfun,xout=FALSE,...)
est1=matl(est1)
est2=matl(est2)
mat[,2]=regYhat(x1,y1,xr=pts,regfun=regfun,...)
mat[,3]=regYhat(x2,y2,xr=pts,regfun=regfun,...)
est=est1-est2
if(!is.matrix(est))est=matrix(est,nrow=1)
mat[,4]=mat[,2]-mat[,3]
pv1=apply(est<0,1,mean,na.rm=TRUE)
pv2=apply(est==0,1,mean,na.rm=TRUE)
pv=pv1+.5*pv2
pv1m=1-pv
pv=2*apply(cbind(pv,pv1m),1,min)
mat[,7]=pv
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
estsort=apply(est,1,sort)
mat[,5]=estsort[ilow,]
mat[,6]=estsort[ihi,]
if(plotit){
if(eout && xout)stop('Cannot have both eout and xout = F')
if(eout){
flag<-outfun(cbind(x1,y1),plotit=FALSE,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(cbind(x2,y2),plotit=FALSE,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
plot(c(x1,x2),c(y1,y2),type='n',xlab=xlab,ylab=ylab)
if(SCAT)points(x1,y1,pch='o')
if(SCAT)points(x2,y2,pch='+')
abline(regfun(x1,y1)$coef)
abline(regfun(x2,y2)$coef,lty=2)
}
list(output=mat)
}


# Danctspb.sub
DEPanc<-function(x1,y1,y2,fr1=1,tr=.2,alpha=.05,plotit=TRUE,DISDIF=FALSE,DIF=TRUE,
pts=NULL,sm=FALSE,xout=FALSE,outfun=out,nboot=500){
#
# Compare two dependent  groups using a covariate
#
# x1 is the covariate and
# y1 and y2 are the two measures. For instance time 1 and time 2.
#
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  fr1 is span for running interval smoother
#
#  sm=T will create smooths using bootstrap bagging.
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  If DISDIF=T: 1. compare groups using median of distribution of D=Y1-Y2
#               2. if na.rm=T, case wise deletion is used, otherwise all of the data are used.
#
#   Also see the R function DEPancB, which includes alternative methods for handling missing values
#
m=cbind(x1,y1,y2)
flag=is.na(x1)
m=m[!flag,]
if(is.null(pts[1])){
npt<-5
isub<-c(1:5)  # Initialize isub
test<-c(1:5)
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
y2<-y2[xorder]
vecn<-1
for(i in 1:length(x1))vecn[i]<-length(y1[near(x1,x1[i],fr1)])
sub<-c(1:length(x1))
isub[1]<-min(sub[vecn>=12])
isub[5]<-max(sub[vecn>=12])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
}
if(!is.null(pts[1]))isub=c(1:length(pts))
mat<-matrix(NA,length(isub),8)
dimnames(mat)<-list(NULL,c("X","n","DIF","TEST","se","ci.low","ci.hi",
"p.value"))
for (i in 1:length(isub)){
if(is.null(pts)){
ch=near(x1,x1[isub[i]],fr1)
mat[i,1]=x1[isub[i]]
}
if(!is.null(pts)){
ch=near(x1,pts[i],fr1)
mat[i,1]=pts[i]
}
mat[i,2]=sum(ch)
if(!DISDIF){
if(!DIF){
test<-yuend(m[ch,2],m[ch,3],tr=tr)
mat[i,3]=mean(m[ch,2],tr=tr)-mean(m[ch,3],tr=tr)
mat[i,4]<-test$teststat
mat[i,5]<-test$se
mat[i,6]<-test$ci[1]
mat[i,7]<-test$ci[2]
mat[i,8]<-test$siglevel
}
if(DIF){
test=trimci(m[ch,2]-m[ch,3],tr=tr,pr=FALSE)
mat[i,3]=mean(m[ch,2]-m[ch,3],tr=tr)
mat[i,4]<-test$test.stat
mat[i,5]<-test$se
mat[i,6]<-test$ci[1]
mat[i,7]<-test$ci[2]
mat[i,8]<-test$p.value
}}
if(DISDIF){
test=l2drmci(m[ch,2:3],pr=FALSE,nboot=nboot,na.rm=na.rm)
mat[i,3]<-loc2dif(m[ch,2],m[ch,3],na.rm=na.rm)
mat[i,4]<-NA
mat[i,5]<-NA
mat[i,6]<-test$ci[1]
mat[i,7]<-test$ci[2]
mat[i,8]<-test$p.value
}}
if(plotit)
runmean2g(x1,y1,x1,y2,fr=fr1,est=mean,tr=tr,sm=sm,xout=xout,outfun=outfun)
list(output=mat)
}



# DEPancpb
DEPancpb<-function(x1,y1,y2,fr1=1,est=tmean,alpha=.05,plotit=TRUE,DISDIF=FALSE,DIF=TRUE,TLS=FALSE,SEED=TRUE,
pts=NULL,sm=FALSE,xout=FALSE,outfun=out,nboot=500,pr=FALSE,na.rm=TRUE,xlab="Group 1", ylab="Group 2",...){
#
# Compare two dependent  groups using a covariate
#
#  same as DEPanc, only use bootstrap methods in all cases.
#
# x1 is the covariate and
# y1 and y2 are the two measures. For instance time 1 and time 2.
#
#   case wise deletion of missing values used by default.
#   To use all of the data not missing, set DIF=F and na.rm=F
#   For the special case where the goal is to compare means, also set TLS=T
#   (But this can produce an error if too many missing values)
#
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  TLS=F, use percentile bootstrap when DIF=FALSE;
#  otherwise (TLS=TRUE) use Lin-Stivers method for means
#  fr1 is span for running interval smoother
#
#  sm=T will create smooths using bootstrap bagging.
#
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
m=cbind(x1,y1,y2)
flag=is.na(x1)
if(na.rm)m=elimna(m)
if(!na.rm){
m=m[!flag,]
}
x1=m[,1]
y1=m[,2]
y2=m[,3]
if(is.null(pts[1])){
npt<-5
isub<-c(1:5)  # Initialize isub
test<-c(1:5)
xorder<-order(x1)
y1<-y1[xorder]
x1<-x1[xorder]
y2<-y2[xorder]
vecn<-1
for(i in 1:length(x1))vecn[i]<-length(y1[near(x1,x1[i],fr1)])
sub<-c(1:length(x1))
isub[1]<-min(sub[vecn>=12])
isub[5]<-max(sub[vecn>=12])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
}
if(!is.null(pts[1]))isub=c(1:length(pts))
mat<-matrix(NA,length(isub),6)
dimnames(mat)<-list(NULL,c("X","n","DIF","ci.low","ci.hi",
"p.value"))
for (i in 1:length(isub)){
if(is.null(pts)){
ch=near(x1,x1[isub[i]],fr1)
mat[i,1]=x1[isub[i]]
}
if(!is.null(pts)){
ch=near(x1,pts[i],fr1)
mat[i,1]=pts[i]
}
mat[i,2]=sum(ch)
if(!DISDIF){
if(!DIF){
if(!TLS){
test=rmmismcp(m[ch,2],m[ch,3],alpha=alpha,SEED=SEED,est=est,plotit = FALSE,
    grp = grp, nboot = 500, xlab = xlab, ylab = ylab, pr = pr, ...)
mat[i,3]=est(m[ch,2],na.rm=TRUE)-est(m[ch,3],na.rm=TRUE)
mat[i,4]<-test$output[1,6]
mat[i,5]<-test$output[1,7]
mat[i,6]<-test$output[1,4]
}
if(TLS){
test=rm2miss(m[ch,2],m[ch,3], nboot = nboot, alpha = alpha, SEED = SEED)
mat[i,3]=mean(m[ch,2],na.rm=TRUE)-mean(m[ch,3],na.rm=TRUE)
mat[i,4]<-test$ci[1]
mat[i,5]<-test$ci[2]
mat[i,6]<-test$p.value
}}
if(DIF){
test=onesampb(m[ch,2]-m[ch,3],est=est,nboot=nboot,alpha=alpha,SEED=SEED,...)
mat[i,3]=est(m[ch,2]-m[ch,3],na.rm=TRUE,...)
mat[i,4]<-test$ci[1]
mat[i,5]<-test$ci[2]
mat[i,6]<-test$p.value
}}
if(DISDIF){
test=l2drmci(m[ch,2:3],pr=FALSE,nboot=nboot,na.rm=na.rm)
mat[i,3]<-loc2dif(m[ch,2],m[ch,3],na.rm=na.rm)
mat[i,4]<-test$ci[1]
mat[i,5]<-test$ci[2]
mat[i,6]<-test$p.value
}}
if(plotit)
runmean2g(x1,y1,x1,y2,fr=fr1,est=est,sm=sm,xout=xout,outfun=outfun,...)
list(output=mat)
}



# oancpb
oancpb<-function(x1,y1,x2,y2,est=tmean,tr=.2,pts=NA,fr1=1,fr2=1,nboot=600,
alpha=.05,plotit=TRUE,SEED=TRUE,PRO=FALSE,...){
#
# Compare two independent  groups using an ancova method
# with a percentile bootstrap combined with a running interval
# smooth.
#
#  CURRENTLY SEEMS THAT THE R FUNCTION ancGLOB is better.
#
#  This function performs an omnibus test using data corresponding
#  to K design points  specified by the argument pts. If
#  pts=NA, K=5 points are chosen for you (see Introduction to Robust
#  Estimation and Hypothesis Testing.)
#  Null hypothesis is that conditional distribution of Y, given X for first
#  group,  minus the conditional distribution of Y, given X for second
#  group is equal to zero.
#  The strategy is to choose K specific X values
#  and then test the hypothesis that all K differences are zero.
#
#  If you want to choose specific X values, Use the argument
#  pts
#  Example: pts=c(1,3,5) will use X=1, 3 and 5.
#
#  For multiple comparisons using these J points, use ancpb
#
#  Assume data are in x1 y1 x2 and y2
#
# PRO=F, means Mahalanobis distance is used.
# PRO=T, projection distance is used.
#
#  fr1 and fr2 are the spans used to fit a smooth to the data.
#
stop('USE ancGLOB')
#
#
gv1<-vector("list")
if(is.na(pts[1])){
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
for (i in 1:5){
j<-i+5
temp1<-y1[near(x1,x1[isub[i]],fr1)]
temp2<-y2[near(x2,x1[isub[i]],fr2)]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
gv1[[i]]<-temp1
gv1[[j]]<-temp2
}
#
loc<-NA
if(SEED)set.seed(2)
bvec<-matrix(NA,nrow=nboot,ncol=5)
for(j in 1:5){
k<-j+5
loc[j]<-est(gv1[[j]])-est(gv1[[k]])
xx<-matrix(sample(gv1[[j]],size=length(gv1[[j]])*nboot,replace=TRUE),
nrow=nboot)
yy<-matrix(sample(gv1[[k]],size=length(gv1[[k]])*nboot,replace=TRUE),
nrow=nboot)
bvec[,j]<-apply(xx,1,FUN=est,...)-apply(yy,1,FUN=est,...)
}
nullv<-rep(0,5)
if(!PRO){
mvec<-apply(bvec,2,FUN=mean)
m1<-var(t(t(bvec)-mvec+loc))
temp<-mahalanobis(rbind(bvec,nullv),loc,m1)
}
if(PRO){
temp<-pdis(rbind(bvec,nullv))
}
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
}
if(!is.na(pts[1])){
npt<-length(pts)
n1<-1
n2<-1
vecn<-1
mat<-matrix(NA,nrow=2*length(pts),ncol=3)
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
gv<-vector("list",2*length(pts))
for (i in 1:length(pts)){
j<-i+npt
temp1<-y1[near(x1,pts[i],fr1)]
temp2<-y2[near(x2,pts[i],fr2)]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
mat[i,1]<-pts[i]
if(length(temp1)<=10)print(paste("Warning, there are",length(temp1)," points corresponding to the design point X=",pts[i]))
if(length(temp2)<=10)print(paste("Warning, there are",length(temp2)," points corresponding to the design point X=",pts[i]))
mat[i,2]<-length(temp1)
mat[i,3]<-length(temp2)
gv1[[i]]<-temp1
gv1[[j]]<-temp2
}
loc<-NA
if(SEED)set.seed(2)
bvec<-matrix(NA,nrow=nboot,ncol=npt)
for(j in 1:npt){
k<-j+npt
loc[j]<-est(gv1[[j]])-est(gv1[[k]])
xx<-matrix(sample(gv1[[j]],size=length(gv1[[j]])*nboot,replace=TRUE),
nrow=nboot)
yy<-matrix(sample(gv1[[k]],size=length(gv1[[k]])*nboot,replace=TRUE),
nrow=nboot)
bvec[,j]<-apply(xx,1,FUN=est,...)-apply(yy,1,FUN=est,...)
}
nullv<-rep(0,npt)
if(!PRO){
mvec<-apply(bvec,2,FUN=mean)
m1<-var(t(t(bvec)-mvec+loc))
temp<-mahalanobis(rbind(bvec,nullv),loc,m1)
}
if(PRO)temp<-pdis(rbind(bvec,nullv))
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
}
if(plotit)runmean2g(x1,y1,x2,y2,fr=fr1,est=est,...)
list(p.value=sig.level)
}


# Qancsm
Qancsm<-function(x1,y1,x2,y2,crit.mat=NULL,nboot=200,SEED=TRUE,REP.CRIT=FALSE,
qval=.5,q=NULL,xlab="X",ylab="Y",plotit=TRUE,pr=TRUE,xout=FALSE,outfun=out,...){
#
# Compare two nonparametric
# regression lines corresponding to two independent groups
#  using the depths of smooths.
#
# NULL hypothesis: regression lines are identical in terms of the median
# of Y, given$X, for all X
# The method is based on comparing the depth of the fitted regression lines
# and is essentially a slight variation of the method in Wilcox
# (in press) Journal of Data Science.
#
# One covariate only is allowed.
#
if(ncol(as.matrix(x1))>1)stop("One covariate only is allowed")
if(!is.null(q))qval=q
if(xout){
flag1=outfun(x1)$keep
flag2=outfun(x2)$keep
x1=x1[flag1]
y1=y1[flag1]
x2=x2[flag2]
y2=y2[flag2]
}
if(SEED)set.seed(2)
xy=elimna(cbind(x1,y1))
x1=xy[,1]
xord=order(x1)
x1=x1[xord]
y1=xy[xord,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
xord=order(x2)
x2=x2[xord]
y2=xy[xord,2]
n1=length(y1)
n2=length(y2)
if(is.null(crit.mat[1])){
if(pr)print("Determining critical value. This might take a while")
crit.val=NA
yall=c(y1,y2)
xall=c(x1,x2)
nn=n1+n2
il=n1+1
for(i in 1:nboot){
data=sample(nn,nn,T)
yy1=yall[data[1:n1]]
yy2=yall[data[il:nn]]
xx1=xall[data[1:n1]]
xx2=xall[data[il:nn]]
crit.mat[i]=Qdepthcom(xx1,yy1,xx2,yy2,qval=qval)
}}
dep=Qdepthcom(x1,y1,x2,y2,qval=qval)
pv=1-mean(crit.mat<dep)
if(!REP.CRIT)crit.mat=NULL
if(plotit){
plot(c(x1,x2),c(y1,y2),type="n",xlab=xlab,ylab=ylab)
temp1=cobs(x1,y1,print.mesg=FALSE,print.warn=FALSE,tau=qval)
temp2=cobs(x2,y2,print.mesg=FALSE,print.warn=FALSE,tau=qval)
points(x1,y1)
points(x2,y2,pch="+")
lines(x1,temp1$fitted)
lines(x2,temp2$fitted,lty=2)
}
list(p.value=pv,crit.mat=crit.mat,test.depth=dep)
}



# QSanc
QSanc<-function(x1,y1,x2,y2,pts,xout=FALSE,outfun=outpro){
#
#  Estimate quantile shift measure of effect size given that the covariate X is equal to the values in pts.
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
e[i]=shiftQS(d1,d2)
}
e
}

# QS.ancbse
QS.ancbse<-function(x1,y1,x2,y2,pts,nboot=100,SEED=TRUE,MC=FALSE,null.value=.5,
xout=FALSE,outfun=outpro,alpha=.05,...){
#
#  ANCOVA based on quantile shift measure of effect size.
#
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
if(SEED)set.seed(2)
e=QSanc(x1,y1,x2,y2,pts=pts)
e=as.vector(matl(e))
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
if(!MC)temp=lapply(bl,QS.ancbse.sub,pts)
else{
temp=mclapply(bl,QS.ancbse.sub,pts)
}
tv=list()
for(j in 1:nboot)tv[[j]]=as.vector(matl(temp[[j]][1:npt]))
E=matl(tv)
E=as.matrix(E)
se=apply(E,1,sd)
ci[,1]<-e-qnorm(1-alpha/2)*se
ci[,2]<-e+qnorm(1-alpha/2)*se
test<-(e-null.value)/se
sig<-2*(1-pnorm(abs(test)))
list(Est=e,SE=se,test.stat=test,conf.int=ci,p.value=sig)
}


# QS.ancbse.sub
QS.ancbse.sub<-function(m,pts){
v=QSanc(m[[1]],m[[2]],m[[3]],m[[4]],pts=pts)
v
}


# QS.ancp2pb
rmanc.best<-function(x,alpha=.05,tr=.2,iter=5000,SEED=TRUE){
#
#
#  Identify group with largest trimmed mean
#  Make a decision if every  p.value<=p.crit
#
#  p.crit is determined via
#  a simulation to determine the null distribution based on
#  iter=5000 replications.
#
#
#  Returns:
#   Best='No Decision' if not significant
#   Best= the group with largest measure if a decision can be made.
#
#   Confidence intervals having simultaneous probability coverage 1-alpha
#  using the adjusted level.
#
x=elimna(x)
flag=TRUE
if(is.list(x))stop('x should be a matrix or a data frame')
J=ncol(x)
if(J<3)stop('Should have 3 or more groups')
Jm1=J-1
est=apply(x,2,tmean,tr=tr)
n=nrow(x)
est=matl(est)
R=order(est,decreasing = TRUE)
pvec=NA
p.crit=rmanc.best.crit(x,iter=iter,alpha=alpha,SEED=SEED)
output<-matrix(NA,Jm1,8)
dimnames(output)=list(NULL,c('Est.Best','Grp','Est','Dif','ci.low','ci.up','p.value','p.crit'))
for(i in 2:J){
im1=i-1
a=yuend(x[,R[1]],x[,R[i]],alpha=p.crit[im1])
pvec[im1]=a$p.value
output[im1,]=c(a$est1,R[[i]],a$est2,a$dif,a$ci[1],a$ci[2],a$p.value,p.crit[im1])
}
Best='No Decisions'
flag=sum(output[,7]<=output[,8])
id=output[,7]<=output[,8]
if(sum(id>0))Best=output[id,2]
if(flag==Jm1)Best='All'
setClass('BIN',slots=c('Group.with.largest.estimate','Larger.than','n','output'))
put=new('BIN',Group.with.largest.estimate=R[[1]],Larger.than=Best,n=n,output=output)
put
}



# rmanc.best.crit
rmanc.best.crit<-function(x,alpha=.05,tr=.2,iter=5000,SEED=TRUE,...){
#
#  Determine critical p-values for rmanc.best
#
if(SEED)set.seed(2)
J=ncol(x)
n=nrow(x)
Jm1=J-1
rem=matrix(NA,iter,Jm1)
A=winall(x,tr=tr)$cov
for(k in 1:iter){
xs=mvrnorm(n,mu=rep(0,J),Sigma=A)
rem[k,]=rmanc.best.ex(xs,tr=tr)
}
init=apply(rem,2,qest,alpha)
z=optim(0,anc.best.fun,init=init,iter=iter,rem=rem,Jm1=Jm1,alpha=alpha,method='Brent',lower=0,upper=1)
fin.crit=z$par*init
fin.crit
}


# rmanc.best.crit.det
rmanc.best.crit.det<-function(x,tr=.2,iter=5000,SEED=TRUE,...){
#
#  Determine critical p-values for rmanc.best
#
if(SEED)set.seed(2)
J=ncol(x)
n=nrow(x)
Jm1=J-1
rem=matrix(NA,iter,Jm1)
A=winall(x,tr=tr)$cov
for(k in 1:iter){
xs=mvrnorm(n,mu=rep(0,J),Sigma=A)
rem[k,]=rmanc.best.ex(xs,tr=tr)
}

aval=c(seq(.001,.1,.001),seq(.011,.99,.01))
na=length(aval)
fin.crit=matrix(NA,na,Jm1)
for(i in 1:na){
init=apply(rem,2,qest,aval[i])
z=optim(0,anc.best.fun,init=init,iter=iter,rem=rem,Jm1=Jm1,alpha=aval[i],method='Brent',lower=0,upper=1)
fin.crit[i,]=z$par*init
}
fin.crit
}


# rmanc.best.DO
rmanc.best.DO<-function(x,tr=.2,...){
#
#  Determine whether it is reasonable to
#  decide which group has largest measure of location
#
#
if(is.list(x))x=matl(x)
x=elimna(x)
x<-listm(x)
J=length(x)
e=lapply(x,tmean,tr)
e=pool.a.list(e)
id=which(e==max(e))
id=id[1]
e=lapply(x,tmean,tr)
e=pool.a.list(e)
id=which(e==max(e))
CON=conCON(J,id)$conCON
a=rmmcp(x,con=CON,dif=FALSE,tr=tr)
pv=max(a$test[,3])
list(Best.Group=id,Est.=e,p.value=pv)
}






# rmanc.best.ex
rmanc.best.ex<-function(x,tr=.2){
#
#
pvec=NA
x=elimna(x)
if(is.matrix(x))x=listm(x)
J=length(x)
est=lapply(x,tmean,tr=tr)
est=matl(est)
R=order(est,decreasing = TRUE)
ic=0
for(j in 2:J){
ic=ic+1
pvec[ic]=yuend(x[[R[1]]],x[[R[[j]]]],tr=tr)$p.value
}
pvec
}


# rmanc.bestPB
rmanc.bestPB<-function(x,alpha=.05,est=tmean,iter=5000,SEED=TRUE,nboot=2000,PB=FALSE,...){
#
#
#  For J dependent groups,
#  identify the group with largest trimmed mean
#  Make a decision if every  p.value<=p.crit
#
#  p.crit is determined via
#  a simulation to determine the null distribution based on
#  iter=5000 replications.
#
#   PB=FALSE: Determine critical values via rmanc.best.crit.  Faster execution time but can differ somewhat
#   from values based on PB method
#
#
#  Returns:
#   Best='No Decision' if not significant
#   Best= the group with largest measure of location if a decision can be made.
#
#   Confidence intervals having simultaneous probability coverage 1-alpha
#  using the adjusted level.
#
if(is.list(x))x=matl(x)
x=elimna(x)
flag=TRUE
J=ncol(x)
if(J<3)stop('Should have 3 or more groups')
Jm1=J-1
est=apply(x,2,tmean,tr=tr)
n=nrow(x)
est=matl(est)
R=order(est,decreasing = TRUE)
pvec=NA
if(!PB)p.crit=rmanc.best.crit(x,iter=iter,alpha=alpha,SEED=SEED)
if(PB)p.crit=rmanc.best.critPB(x,iter=iter,alpha=alpha,SEED=SEED)
output<-matrix(NA,Jm1,8)
dimnames(output)=list(NULL,c('Est.Best','Grp','Est','Dif','ci.low','ci.up','p.value','p.crit'))
for(i in 2:J){
im1=i-1
a=yuend(x[,R[1]],x[,R[i]],alpha=p.crit[im1],tr=tr)
pvec[im1]=a$p.value
output[im1,]=c(a$est1,R[[i]],a$est2,a$dif,a$ci[1],a$ci[2],a$p.value,p.crit[im1])
}
Best='No Decisions'
flag=sum(output[,7]<=output[,8])
id=output[,7]<=output[,8]
if(sum(id>0))Best=output[id,2]
if(flag==Jm1)Best='All'
setClass('BIN',slots=c('Group.with.largest.estimate','Larger.than','n','output'))
put=new('BIN',Group.with.largest.estimate=R[[1]],Larger.than=Best,n=n,output=output)
put
}


# rmanc.best.PV
rmanc.best.PV<-function(x,alpha=.05,tr=.2,iter=5000,SEED=TRUE){
#
#
#  For J dependent groups,
#  identify the group with largest trimmed mean
#  Make a decision if every  p.value<=p.crit
#
#  p.crit is determined via
#  a simulation to determine the null distribution based on
#  iter=5000 replications.
#
#
#  Returns:
#   Best='No Decision' if not significant
#   Best= the group with largest measure of location if a decision can be made.
#
#   Confidence intervals having simultaneous probability coverage 1-alpha
#  using the adjusted level.
#
x=elimna(x)
flag=TRUE
if(is.list(x))stop('x should be a matrix or a data frame')
J=ncol(x)
if(J<3)stop('Should have 3 or more groups')
Jm1=J-1
est=apply(x,2,tmean,tr=tr)
n=nrow(x)
est=matl(est)
R=order(est,decreasing = TRUE)
pvec=NA

aval=c(seq(.001,.1,.001),seq(.11,.99,.01))
id=which(aval==alpha)
if(length(id)==0)stop('alpha be one one values .001(.001).1 or 11(.01).99')
v=rmanc.best.crit.det(x,iter=iter,alpha=alpha,tr=tr,SEED=SEED)
p.crit=v[id,]
output<-matrix(NA,Jm1,8)
dimnames(output)=list(NULL,c('Est.Best','Grp','Est','Dif','ci.low','ci.up','p.value','p.crit'))
for(i in 2:J){
im1=i-1
a=yuend(x[,R[1]],x[,R[i]],alpha=p.crit[im1],tr=tr)
pvec[im1]=a$p.value
output[im1,]=c(a$est1,R[[i]],a$est2,a$dif,a$ci[1],a$ci[2],a$p.value,p.crit[im1])
}

# Determine p-value for overall decision
na=length(aval)
for(i in 1:na){
chk=sum(output[,7]<=v[i,])
pv=aval[i]
if(chk==Jm1)break
}
Best='No Decisions'
flag=sum(output[,7]<=output[,8])
id=output[,7]<=output[,8]
if(sum(id>0))Best=output[id,2]
if(flag==Jm1)Best='All'
setClass('BIN',slots=c('Group.with.largest.estimate','Select.Best.p.value','Larger.than','n','output'))
put=new('BIN',Group.with.largest.estimate=R[[1]],Select.Best.p.value=pv,Larger.than=Best,n=n,output=output)
put
}


# smeancrv2
