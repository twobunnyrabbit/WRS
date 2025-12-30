# ==============================================================================
# Correlation and Association Functions
# ==============================================================================
#
# This module contains functions for computing correlations and associations:
#   - Pearson correlation: pbcor(), pcor(), wincor()
#   - Spearman correlation: scor(), scorci(), spear()
#   - Kendall's tau: tau(), tauci()
#   - Percentage bend correlation: pbcor(), pcorhc4()
#   - Skipped correlations: scorci(), scorall()
#   - Multiple correlation: mscor(), mscorci()
#   - Regression-based: correg(), scorreg()
#
# Extracted: 2025-12-30
# Functions: 96
# ==============================================================================

wincor<-function(x,y=NULL,tr=.2){
#   Compute the Winsorized correlation between x and y.
#
#   tr is the amount of Winsorization
#   This function also returns the Winsorized covariance
#
#    Pairwise deletion of missing values is performed.
#
#   x is a vector, or it can be a matrix with two columns when y=NULL
#

if(!is.null(y[1])){
m=cbind(x,y)
}
else m=x
m<-elimna(m)
nval=nrow(m)
if(ncol(m)==2){
a=wincor.sub(m[,1],m[,2],tr=tr)
wcor=a$cor
wcov=a$cov
sig=a$p.value
}
if(ncol(m)>2){
#if(is.data.frame(m))m=as.matrix(m)
if(!is.matrix(m))stop("The data must be stored in a n by p matrix")
wcor<-matrix(1,ncol(m),ncol(m))
wcov<-matrix(0,ncol(m),ncol(m))
siglevel<-matrix(NA,ncol(m),ncol(m))
for (i in 1:ncol(m)){
ip<-i
for (j in ip:ncol(m)){
val<-wincor.sub(m[,i],m[,j],tr)
wcor[i,j]<-val$cor
wcor[j,i]<-wcor[i,j]
if(i==j)wcor[i,j]<-1
wcov[i,j]<-val$cov
wcov[j,i]<-wcov[i,j]
if(i!=j){
siglevel[i,j]<-val$p.value
siglevel[j,i]<-siglevel[i,j]
}
}}
sig=siglevel
}
list(n=nval,cor=wcor,cov=wcov,p.value=sig)
}

wincor.sub<-function(x,y,tr=tr){
sig<-NA
g<-floor(tr*length(x))
xvec<-winval(x,tr)
yvec<-winval(y,tr)
wcor<-cor(xvec,yvec)
wcov<-var(xvec,yvec)
if(sum(x==y)!=length(x)){
test<-wcor*sqrt((length(x)-2)/(1.-wcor^2))
sig<-2*(1-pt(abs(test),length(x)-2*g-2))
}
list(cor=wcor,cov=wcov,p.value=sig)
}

tau<-function(x,y=NULL,alpha=.05){
#
#   Compute Kendall's tau plus a 1-alpha confidence interval
#   using the method recommended by Long and Cliff (1997).
#
#  tau-b provides an adjustment for ties. However, no adjustment for tied
#  values is made here due to arguments made by N. Cliff (1996) Ordinal Methods for
#  Behavioral Data Analysis. Mahwah, NJ: Erlbaum. pp. 36-37.
#
#   y=NULL, assume x is a matrix with two columns
#
if(is.null(y))m=elimna(x)
if(!is.null(y))m=elimna(cbind(x,y)) # casewise deletion of missing values.
x=m[,1]
y=m[,2]
xdif<-outer(x,x,FUN="-")
ydif<-outer(y,y,FUN="-")
tv<-sign(xdif)*sign(ydif)
n<-length(x)
df=n-2
dbar<-apply(tv,1,sum)/(n-1)
tau<-sum(tv)/(n*(n-1))
A<-sum((dbar-tau)^2)/(n-1)
B<-(n*(n-1)*(-1)*tau^2+sum(tv^2))/(n^2-n-1)
C<-(4*(n-2)*A+2*B)/(n*(n-1))
#crit<-qnorm(alpha/2)
crit=qt(alpha/2,df)
cilow<-tau+crit*sqrt(C)
cihi<-tau-crit*sqrt(C)
test<-tau/sqrt((2*(2*n+5))/(9*n*(n-1)))
siglevel<-2*(1-pt(abs(test),df=df))
list(cor=tau,ci=c(cilow,cihi),p.value=siglevel)
}

tauall<-function(m){
#
#    Compute Kendall's tau for the
#    data in the n-by-p matrix m.
#
#    This function also returns the two-sided significance level
#    for all pairs of variables, plus a test of zero correlations
#    among all pairs. (See chapter 6 for details.)
#
if(!is.matrix(m))stop("Data must be stored in an n by p matrix")
taum<-matrix(0,ncol(m),ncol(m))
siglevel<-matrix(NA,ncol(m),ncol(m))
for (i in 1:ncol(m)){
ip1<-i
for (j in ip1:ncol(m)){
if(i<j){
pbc<-tau(m[,i],m[,j])
taum[i,j]<-pbc$cor
taum[j,i]<-pbc$cor
siglevel[i,j]<-pbc$p.value
siglevel[j,i]<-siglevel[i,j]
}
}
}
list(taum=taum,p.value=siglevel)
}

pbcor<-function(x,y,beta=.2){
#   Compute the percentage bend correlation between x and y.
#
#   beta is the bending constant for omega sub N.
#
if(length(x)!=length(y))stop("The vectors do not have equal lengths")
m1=cbind(x,y)
m1<-elimna(m1)
nval=nrow(m1)
x<-m1[,1]
y<-m1[,2]
#  Have eliminated missing values
temp<-sort(abs(x-median(x)))
omhatx<-temp[floor((1-beta)*length(x))]
temp<-sort(abs(y-median(y)))
omhaty<-temp[floor((1-beta)*length(y))]
a<-(x-pbos(x,beta))/omhatx
b<-(y-pbos(y,beta))/omhaty
a<-ifelse(a<=-1,-1,a)
a<-ifelse(a>=1,1,a)
b<-ifelse(b<=-1,-1,b)
b<-ifelse(b>=1,1,b)
pbcor<-sum(a*b)/sqrt(sum(a^2)*sum(b^2))
test<-pbcor*sqrt((length(x) - 2)/(1 - pbcor^2))
sig<-2*(1 - pt(abs(test),length(x)-2))
list(cor=pbcor,test=test,p.value=sig,n=nval)
}

corb<-function(x,y,corfun=pbcor,nboot=599,alpha=.05,plotit=FALSE,xlab='X',ylab='Y',SEED=TRUE,...){
#
#   Compute a 1-alpha confidence interval for a correlation.
#   The default correlation is the percentage bend.
#
#   The function corfun is any R function that returns a
#   correlation coefficient in corfun$cor. The functions pbcor and
#   wincor follow this convention.
#
#   When using Pearson's correlation, and when n<250, use
#   lsfitci instead.
#
#   The default number of bootstrap samples is nboot=599
#
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1]
y<-m1[,2]
est<-corfun(x,y,...)$cor
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,corbsub,x,y,corfun,...) # A 1 by nboot matrix.
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
bsort<-sort(bvec)
corci<-1
corci[1]<-bsort[ilow]
corci[2]<-bsort[ihi]
phat <- sum(bvec < 0)/nboot
sig <- 2 * min(phat, 1 - phat)
if(plotit)outpro(cbind(x,y),xlab=xlab,ylab=ylab,plotit=TRUE)
list(cor.ci=corci,p.value=sig,cor.est=est)
}

runcor<-function(x,y,z,fr=1,corflag=FALSE,corfun=pbcor,plotit=TRUE,rhat=FALSE){
#
# Estimate how the correlation between  x and y varies with  z
#
# running correlation using interval method
#
# fr controls amount of smoothing
#
# corfun is the  correlation to be used. It is assumed that
# corfun is an R function that returns a correlation coefficient
# in corfun$cor
#
# To use Pearsons correlation, set corflag=T
#
temp<-cbind(x,y,z) # Eliminate any rows with missing values
temp<-elimna(temp)
x<-temp[,1]
y<-temp[,2]
z<-temp[,3]
plotit<-as.logical(plotit)
rmd<-NA
if(!corflag){
for(i in 1:length(x)){
flag<-near(z,z[i],fr)
if(sum(flag)>2)rmd[i]<-corfun(x[flag],y[flag])$cor
}}
if(corflag){
for(i in 1:length(x)){
flag<-near(z,z[i],fr)
if(sum(flag)>2)rmd[i]<-cor(x[flag],y[flag])
}}
if(plotit){
plot(c(max(z),min(z),z),c(1,-1,rmd),xlab="Modifier",ylab="Correlation",type="n")
sz<-sort(z)
zorder<-order(z)
sysm<-rmd[zorder]
lines(sz,sysm)
}
if(!rhat)rmd<-"Done"
rmd
}

pcorb<-function(x,y,SEED=TRUE){
#   Compute a .95 confidence interval for Pearson's correlation coefficient.
#
#   This function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
#
nboot<-599  #Number of bootstrap samples
xy<-elimna(cbind(x,y))
x<-xy[,1]
y<-xy[,2]
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples; please wait")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,pcorbsub,x,y) # A 1 by nboot matrix.
ilow<-15
ihi<-584
if(length(y) < 250){
ilow<-14
ihi<-585
}
if(length(y) < 180){
ilow<-11
ihi<-588
}
if(length(y) < 80){
ilow<-8
ihi<-592
}
if(length(y) < 40){
ilow<-7
ihi<-593
}
bsort<-sort(bvec)
r<-cor(x,y)
ci<-c(bsort[ilow],bsort[ihi])
list(r=r,ci=ci)
}

correg.sub<-function(X,theta,corfun=tau){
np<-ncol(X)
p<-np-1
x<-X[,1:p]
y<-X[,np]
temp<-t(t(x)*theta)
yhat<-apply(temp,1,sum)
yhat<-yhat
res<-y-yhat
val<-sum(abs(taureg(x,res,corfun=corfun)$cor))
val
}

correg<-function(x,y,corfun=tau,loc.fun=median){
#
# A generalization of the Theil-Sen estimator
# Rather than use Kendall's tau, can use an alternative
# correlation via the argument corfun.
# loc.fun determines how the intercept is computed;
#
# The Nelder-Mead method is used rather than
# Gauss-Seidel.
#
#
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
N<-np-1
temp<-tsreg(x,y)$coef
START<-temp[2:np]
temp<-nelderv2(X,N,FN=correg.sub,START=START,corfun=corfun)
x <- as.matrix(x)
alpha <- loc.fun(y - x %*% temp)
coef <- c(alpha,temp)
res <- y - x %*% temp - alpha
list(coef = coef, residuals = res)
}

tauloc<-function(x,cval=4.5){
#
# Compute the tau measure of location as described in
# Yohai and Zamar (JASA, 83, 406-413).
#
x<-elimna(x)
s<-qnorm(.75)*mad(x)
y<-(x-median(x))/s
W<-(1-(y/cval)^2)^2
flag<-(abs(W)>cval)
W[flag]<-0
val<-sum(W*x)/sum(W)
val
}

tauvar<-function(x,cval=3){
#
# Compute the tau measure of scale as described in
# Yohai and Zamar (JASA, 1988, 83, 406-413).
# The computational method is described in Maronna and Zamar
# (Technometrics, 2002, 44, 307-317)
#  see p. 310
#
x<-elimna(x)
s<-qnorm(.75)*mad(x)
y<-(x-tauloc(x))/s
cvec<-rep(cval,length(x))
W<-apply(cbind(y^2,cvec^2),1,FUN="min")
val<-s^2*sum(W)/length(x)
val
}

taulc<-function(x,mu.too=FALSE){
#
val<-tauvar(x)
if(mu.too){
val[2]<-val
val[1]<-tauloc(x)
}
val
}

erho.bt <- function(p,c1,M)
#   expectation of rho(d) under chi-squared p
    return(chi.int(p,2,M)/2
        +(M^2/2+c1*(5*c1+16*M)/30)*chi.int2(p,0,M+c1)
        +(M^2/2-M^2*(M^4-5*M^2*c1^2+15*c1^4)/(30*c1^4))*(
chi.int(p,0,M+c1)-chi.int(p,0,M))
        +(1/2+M^4/(2*c1^4)-M^2/c1^2)*(chi.int(p,2,M+c1)-chi.int(p,2,M))
        +(4*M/(3*c1^2)-4*M^3/(3*c1^4))*(chi.int(p,3,M+c1)-chi.int(p,3,M))
        +(3*M^2/(2*c1^4)-1/(2*c1^2))*(chi.int(p,4,M+c1)-chi.int(p,4,M))
        -(4*M/(5*c1^4))*(chi.int(p,5,M+c1)-chi.int(p,5,M))
        +(1/(6*c1^4))*(chi.int(p,6,M+c1)-chi.int(p,6,M)))
chi.int <- function(p,a,c1)
#   partial expectation d in (0,c1) of d^a under chi-squared p
  return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*pchisq(c1^2,p+a) )

erho.bt.lim <- function(p,c1)
#   expectation of rho(d) under chi-squared p
  return(chi.int(p,2,c1)+c1^2*chi.int2(p,0,c1))
erho.bt.lim.p <- function(p,c1)
#   derivative of erho.bt.lim wrt c1
  return(chi.int.p(p,2,c1)+c1^2*chi.int2.p(p,0,c1)+2*c1*chi.int2(p,0,c1))


rejpt.bt.lim <- function(p,r){
#   find p-value of translated biweight limit c
#   that gives a specified breakdown
    c1 <- 2*p
    iter <- 1
    crit <- 100
    eps <- 1e-5
    while ((crit > eps)&(iter<100))
    {
        c1.old <- c1
        fc <- erho.bt.lim(p,c1) - c1^2*r
        fcp <- erho.bt.lim.p(p,c1) - 2*c1*r
        c1 <- c1 - fc/fcp
        if (c1 < 0)  c1 <- c1.old/2
        crit <- abs(fc)
        iter <- iter+1
    }
    return(c(c1,pchisq(c1^2,p),log10(1-pchisq(c1^2,p))))
}

erho.bt.lim.p <- function(p,c1)
#   derivative of erho.bt.lim wrt c1
  return(chi.int.p(p,2,c1)+c1^2*chi.int2.p(p,0,c1)+2*c1*chi.int2(p,0,c1))


rejpt.bt.lim <- function(p,r){
#   find p-value of translated biweight limit c
#   that gives a specified breakdown
    c1 <- 2*p
    iter <- 1
    crit <- 100
    eps <- 1e-5
    while ((crit > eps)&(iter<100))
    {
        c1.old <- c1
        fc <- erho.bt.lim(p,c1) - c1^2*r
        fcp <- erho.bt.lim.p(p,c1) - 2*c1*r
        c1 <- c1 - fc/fcp
        if (c1 < 0)  c1 <- c1.old/2
        crit <- abs(fc)
        iter <- iter+1
    }
    return(c(c1,pchisq(c1^2,p),log10(1-pchisq(c1^2,p))))
}

ecor<-function(x,y,pcor=FALSE,regfun=tsreg,corfun=pbcor,outkeep=FALSE,outfun=outmgvf){
#
# Estimate the explanatory correlation between x and y
#
# It is assumed that x is a vector or a matrix having one column only
xx<-elimna(cbind(x,y)) # Remove rows with missing values
x<-xx[,1]
y<-xx[,2]
x<-as.matrix(x)
if(ncol(x) > 1)stop("x must be a vector or matrix with one column")
flag<-rep(TRUE,nrow(x))
if(!outkeep){
temp<-outfun(cbind(x,y))$out.id
flag[temp]<-FALSE
}
coef<-regfun(x,y)$coef
ip<-ncol(x)+1
yhat<-x %*% coef[2:ip] + coef[1]
if(pcor)epow2<-cor(yhat[flag],y[flag])^2
if(!pcor)epow2<-corfun(yhat[flag],y[flag])$cor^2
ecor<-sqrt(epow2)*sign(coef[2])
ecor
}

ocor<-function(x,y,corfun=pbcor,outfun=outmgvf,pcor=FALSE,plotit=FALSE){
#
#  Compute a correlation when outliers are ignored.
#
xx<-elimna(cbind(x,y)) # Remove rows with missing values
x<-xx[,1]
y<-xx[,2]
flag<-rep(TRUE,length(x))
temp<-outfun(cbind(x,y),plotit=plotit)$out.id
flag[temp]<-FALSE
if(pcor)ocor<-cor(x[flag],y[flag])
if(!pcor)ocor<-corfun(x[flag],y[flag])$cor
list(cor=ocor)
}

cori<-function(x,y,z,pt=median(z),fr=.8,est=onestep,corfun=pbcor,testit=FALSE,
nboot=599,sm=FALSE,xlab="X",ylab="Y",...){
#
# Split the data according to whether z is < or > pt, then
# use runmean2g to plot a smooth of the regression
# lines corresponding to these two groups.
#
# If testit=T, the hypothesis of equal correlations is tested using the
#  the R function twocor
#
m<-cbind(x,y,z)
m<-elimna(m)
x<-m[,1]
y<-m[,2]
z<-m[,3]
flag<-(z<pt)
runmean2g(x[flag],y[flag],x[!flag],y[!flag],fr=fr,est=est,sm=sm,
xlab=xlab,ylab=ylab,...)
output<-"Done"
if(testit){
output<-twocor(x[flag],y[flag],x[!flag],y[!flag],corfun=corfun,nboot=nboot,plotit=FALSE)
}
output
}

pcor<-function(x,y=NA){
if(!is.na(y[1]))temp<-wincor(x,y,tr=0)
if(is.na(y[1]))temp<-winall(x,tr=0)
list(cor=temp$cor,p.value=temp$p.value)
}

mscor<-function(m,corfun=spear,cop=3,MM=FALSE,gval=NA,ap=TRUE,pw=TRUE,STAND=TRUE,
outfun=outpro,alpha=.05){
#
# m is an n by p matrix
#
# Compute a skipped correlation matrix
#
#  corfun indicates the correlation to be used
#  corfun=pcor uses Pearson's correlation
#  corfun=spear uses Spearman's correlation
#
#  When calling outpro,
#  STAND=T means marginals are first standardized.
# This function returns the p by p matrix of correlations
#
# Method: Eliminate outliers using a projection technique.
# That is, compute Donoho-Gasko median, for each point
# consider the line between it and the median,
# project all points onto this line, and
# check for outliers using a boxplot rule.
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# cop determines how center of the scatterplot is
# estimated; see the function outpro.
# cop=l Donoho-Gasko halfspace median
# cop=2 MCD measure of location
# cop=3 marginal medians
# cop=4 MVE measure of location
#
# Eliminate any outliers and compute
# correlations using remaining data.
#
# gval is critical value for determining whether a point
# is an outlier. It is determined automatically if not specified,
# assuming that Spearman's correlation is used. Critical
# values when using some other correlation have not been
# determined.
#
# Hypothesis of zero correlations tested with FWE=.05
#
# AGRUMENTS:
# MM; see function outpro
# ap=T all pairwise comparisons are tested
# ap=F first variable is tested versus all others
# (for a total of p-1 tests).
# pw=T, print message about high execution time
# pw=F, suppress the message.
#
if(alpha!=.05)stop('For alpha other than .05, use mscorpb or mscorpbMC')
m<-elimna(m)
p<-ncol(m)
pm<-p-1
n<-nrow(m)
if(p<2)stop("Something wrong; number of variables is < 2")
if(pw && cop==1){
print("If execution time is too high,")
print("use cop=2 or 4 rather than 1")
}
if(ap){
inter<-c(2.374,2.780,3.030,3.208,3.372,3.502,3.722,3.825,3.943)
slope<-c(5.333,8.8,25.67,32.83,51.53,75.02,111.34,123.16,126.72)
expo<-c(-1,-1,-1.2,-1.2,-1.3,-1.4,-1.5,-1.5,-1.5)
if(p>10){
qvec<-NA
for(i in 1:9)qvec[i]<-inter[i]+slope[i]*n^expo[i]
pval<-c(2:10)
temp<-lsfit(pval,qvec)$coef
}
}
if(!ap){
inter<-c(2.374,2.54,2.666,2.92,2.999,3.097,3.414,3.286,3.258)
slope<-c(5.333,8.811,14.89,20.59,51.01,52.15,58.498,64.934,59.127)
expo<-c(-1,-1,-1.2,-1.2,-1.5,-1.5,-1.5,-1.5,-1.5)
if(p>10){
qvec<-NA
for(i in 1:9)qvec[i]<-inter[i]+slope[i]*n^expo[i]
pval<-c(1:9)
temp<-lsfit(pval,qvec)$coef
}
}
if(p<=10)crit<-inter[pm]+slope[pm]*n^expo[pm]
if(p>10)crit<-temp[2]*p+temp[1]
if(cop!=1 && is.na(gval))gval<-sqrt(qchisq(.975,ncol(m)))
temp<-outfun(m,plotit=FALSE,MM=MM,gval=gval,cop=cop,STAND=STAND)$keep
mcor<-corfun(m[temp,])$cor
test<-abs(mcor*sqrt((nrow(m)-2)/(1-mcor^2)))
diag(test) <- NA
if(!ap){
test<-as.matrix(test[1,])
}
list(cor=mcor,crit.val=crit,test.stat=test)
}

rhom<-function(x,y,op=1,op2=FALSE,tr=.2,plotit=TRUE,xlab="NA",ylab="NA",zlab="ABS(res)",
est=median,sm=FALSE,SEED=TRUE,xout=FALSE,outfun=outpro,...){
# For regression model, Y=m(X)+s(X)e,
# where s(X) models heteroscedasticity, and e has median 0,
# test hypothesis s(X)=1 for any X
#
# For p>1, method tests for each p whether residuals and x_j
# have a horizontal regression line.
#
# op2=F, tests for homogeneity using running interval smoother
# op2=T, test of independence based on Y-M(Y), M(Y) some measure
#        of location given by argument est.
#  In general, op2=T should NOT be used when the goal is to test
#  the hypothesis of a homoscedastic error term.
#
# op=1 test using regression method (function regci)
# op=2 test using Winsorized correlation
#      tr is amount of Winsorizing. A heteroscedastic bootstrap method is used. wincor is not asymptotically correct.
# op=3 test using a wild boostrap method
#
x<-as.matrix(x)
p<-ncol(x)
pp<-p+1
xy<-elimna(cbind(x,y))
x<-xy[,1:p]
y<-xy[,pp]
x<-as.matrix(x)
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,pp]
x<-as.matrix(x)
}
output<-NA
if(ncol(x)==1){
if(!op2)res<-y-runhat(x[,1],y,est=est,pts=x)
if(op2)res<-y-est(y)
if(op==1)output<-regci(x,abs(res),SEED=SEED,pr=FALSE)$regci[2,5]
if(op==2)output<-wincorci(x,abs(res),tr=tr,SEED=SEED)$p.value
if(op==3)output<-indt(x,abs(res),SEED=SEED)$p.value.d
}
if(ncol(x)>1){
pv<-ncol(x)+1
if(!op2)res<-y-rung3hat(x,y,est=est,pts=x)$rmd
if(op2)res<-y-est(y)
if(op==1)output<-regci(x,abs(res),pr=FALSE)$regci[2:pv,5]
if(op==2)output<-winall(cbind(x,abs(res)),tr=tr)$p.values[1:ncol(x),pv]
if(op==3)output<-indt(x,abs(res),SEED=SEED)$p.value.d
}
if(plotit){
if(ncol(x)==1){
if(xlab=='NA')xlab="X"
if(ylab=='NA')ylab="ABS(res)"
if(!sm)rungen(x,abs(res),est=est,xlab=xlab,ylab=ylab)
if(sm)runmbo(x,abs(res),est=est,xlab=xlab,ylab=ylab)
}
if(ncol(x)==2){
if(xlab=='NA')xlab="X1"
if(ylab=='NA')ylab="X2"
if(sm)rung3d(x,abs(res),est=est,xlab=xlab,ylab=ylab,zlab=zlab)
if(!sm)run3bo(x,abs(res),est=est,xlab=xlab,ylab=ylab,zlab=zlab)
}}
list(p.value=output)
}

tbscor<-function(x,y=NA){
#
# Compute a correlation coefficient using the TBS measure of scatter
#
if(!is.na(y[1]))x<-cbind(x,y)
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)
n<-nrow(x)
p<-ncol(x)
temp<-tbs(x)$cov
val<-matrix(NA,p,p)
for(j in 1:p){
for(k in 1:p){
val[j,k]<-temp[k,j]/sqrt(temp[k,k]*temp[j,j])
}}
test<-abs(val*sqrt((n-2)/(1-val^2)))
if(p==2){
val<-val[1,2]
p.value<-c("Greater than .1")
crit<-20.20/n+1.89
if(test>=crit)p.value<-c("Less than .1")
crit<-30.41/n+2.21
if(test>=crit)p.value<-c("Less than .05")
crit<-39.72/n+2.5
if(test>=crit)p.value<-c("Less than .025")
crit<-58.55/n+2.80
if(test>=crit)p.value<-c("Less than .01")
}
list(cor=val,test.stat=test,p.value=p.value)
}

erho.bt <- function(p,c1,M)
#   expectation of rho(d) under chi-squared p
    return(chi.int(p,2,M)/2
        +(M^2/2+c1*(5*c1+16*M)/30)*chi.int2(p,0,M+c1)
        +(M^2/2-M^2*(M^4-5*M^2*c1^2+15*c1^4)/(30*c1^4))*(
chi.int(p,0,M+c1)-chi.int(p,0,M))
        +(1/2+M^4/(2*c1^4)-M^2/c1^2)*(chi.int(p,2,M+c1)-chi.int(p,2,M))
        +(4*M/(3*c1^2)-4*M^3/(3*c1^4))*(chi.int(p,3,M+c1)-chi.int(p,3,M))
        +(3*M^2/(2*c1^4)-1/(2*c1^2))*(chi.int(p,4,M+c1)-chi.int(p,4,M))
        -(4*M/(5*c1^4))*(chi.int(p,5,M+c1)-chi.int(p,5,M))
        +(1/(6*c1^4))*(chi.int(p,6,M+c1)-chi.int(p,6,M)))
chi.int <- function(p,a,c1)
#   partial expectation d in (0,c1) of d^a under chi-squared p
  return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*pchisq(c1^2,p+a) )

erho.bt.lim <- function(p,c1)
#   expectation of rho(d) under chi-squared p
  return(chi.int(p,2,c1)+c1^2*chi.int2(p,0,c1))
erho.bt.lim.p <- function(p,c1)
#   derivative of erho.bt.lim wrt c1
  return(chi.int.p(p,2,c1)+c1^2*chi.int2.p(p,0,c1)+2*c1*chi.int2(p,0,c1))


rejpt.bt.lim <- function(p,r){
#   find p-value of translated biweight limit c
#   that gives a specified breakdown
    c1 <- 2*p
    iter <- 1
    crit <- 100
    eps <- 1e-5
    while ((crit > eps)&(iter<100))
    {
        c1.old <- c1
        fc <- erho.bt.lim(p,c1) - c1^2*r
        fcp <- erho.bt.lim.p(p,c1) - 2*c1*r
        c1 <- c1 - fc/fcp
        if (c1 < 0)  c1 <- c1.old/2
        crit <- abs(fc)
        iter <- iter+1
    }
    return(c(c1,pchisq(c1^2,p),log10(1-pchisq(c1^2,p))))
}

erho.bt.lim.p <- function(p,c1)
#   derivative of erho.bt.lim wrt c1
  return(chi.int.p(p,2,c1)+c1^2*chi.int2.p(p,0,c1)+2*c1*chi.int2(p,0,c1))


rejpt.bt.lim <- function(p,r){
#   find p-value of translated biweight limit c
#   that gives a specified breakdown
    c1 <- 2*p
    iter <- 1
    crit <- 100
    eps <- 1e-5
    while ((crit > eps)&(iter<100))
    {
        c1.old <- c1
        fc <- erho.bt.lim(p,c1) - c1^2*r
        fcp <- erho.bt.lim.p(p,c1) - 2*c1*r
        c1 <- c1 - fc/fcp
        if (c1 < 0)  c1 <- c1.old/2
        crit <- abs(fc)
        iter <- iter+1
    }
    return(c(c1,pchisq(c1^2,p),log10(1-pchisq(c1^2,p))))
}

erho.bt <- function(p,c1,M)
#   expectation of rho(d) under chi-squared p
    return(chi.int(p,2,M)/2
        +(M^2/2+c1*(5*c1+16*M)/30)*chi.int2(p,0,M+c1)
        +(M^2/2-M^2*(M^4-5*M^2*c1^2+15*c1^4)/(30*c1^4))*(
chi.int(p,0,M+c1)-chi.int(p,0,M))
        +(1/2+M^4/(2*c1^4)-M^2/c1^2)*(chi.int(p,2,M+c1)-chi.int(p,2,M))
        +(4*M/(3*c1^2)-4*M^3/(3*c1^4))*(chi.int(p,3,M+c1)-chi.int(p,3,M))
        +(3*M^2/(2*c1^4)-1/(2*c1^2))*(chi.int(p,4,M+c1)-chi.int(p,4,M))
        -(4*M/(5*c1^4))*(chi.int(p,5,M+c1)-chi.int(p,5,M))
        +(1/(6*c1^4))*(chi.int(p,6,M+c1)-chi.int(p,6,M)))
chi.int <- function(p,a,c1)
#   partial expectation d in (0,c1) of d^a under chi-squared p
  return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*pchisq(c1^2,p+a) )

erho.bt.lim <- function(p,c1)
#   expectation of rho(d) under chi-squared p
  return(chi.int(p,2,c1)+c1^2*chi.int2(p,0,c1))
erho.bt.lim.p <- function(p,c1)
#   derivative of erho.bt.lim wrt c1
  return(chi.int.p(p,2,c1)+c1^2*chi.int2.p(p,0,c1)+2*c1*chi.int2(p,0,c1))


rejpt.bt.lim <- function(p,r){
#   find p-value of translated biweight limit c
#   that gives a specified breakdown
    c1 <- 2*p
    iter <- 1
    crit <- 100
    eps <- 1e-5
    while ((crit > eps)&(iter<100))
    {
        c1.old <- c1
        fc <- erho.bt.lim(p,c1) - c1^2*r
        fcp <- erho.bt.lim.p(p,c1) - 2*c1*r
        c1 <- c1 - fc/fcp
        if (c1 < 0)  c1 <- c1.old/2
        crit <- abs(fc)
        iter <- iter+1
    }
    return(c(c1,pchisq(c1^2,p),log10(1-pchisq(c1^2,p))))
}

erho.bt.lim.p <- function(p,c1)
#   derivative of erho.bt.lim wrt c1
  return(chi.int.p(p,2,c1)+c1^2*chi.int2.p(p,0,c1)+2*c1*chi.int2(p,0,c1))


rejpt.bt.lim <- function(p,r){
#   find p-value of translated biweight limit c
#   that gives a specified breakdown
    c1 <- 2*p
    iter <- 1
    crit <- 100
    eps <- 1e-5
    while ((crit > eps)&(iter<100))
    {
        c1.old <- c1
        fc <- erho.bt.lim(p,c1) - c1^2*r
        fcp <- erho.bt.lim.p(p,c1) - 2*c1*r
        c1 <- c1 - fc/fcp
        if (c1 < 0)  c1 <- c1.old/2
        crit <- abs(fc)
        iter <- iter+1
    }
    return(c(c1,pchisq(c1^2,p),log10(1-pchisq(c1^2,p))))
}

scor<-function(x,y=NULL,corfun=pcor,gval=NA,plotit=FALSE,op=TRUE,MM=FALSE,cop=3,xlab='VAR 1',
ylab='VAR 2',STAND=TRUE,pr=TRUE,SEED=TRUE,MC=FALSE,RAN=FALSE){
#
# Compute a skipped correlation coefficient.
#
# Eliminate outliers using a projection method
# That is, compute Donoho-Gasko median, for each point
# consider the line between it and the median,
# project all points onto this line, and
# check for outliers using a boxplot rule.
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# For information about the argument cop, see the function
# outpro.
#
# Eliminate any outliers and compute correlation using
# remaining data.
#
#  MC=TRUE, the multicore version of outpro is used
#
# corfun=pcor means Pearson's correlation is used.
# corfun=spear means Spearman's correlation is used.
# corfun=tau means Kendall tau is used.
#
#.  RAN=TRUE uses random projections instead, which results in faster execution time
#
if(SEED){
oldSeed <- .Random.seed
set.seed(12) # So when using MVE or MCD, get consistent results
}
if(identical(corfun,wincor))corfun=winall
if(is.null(y[1]))m<-x
if(!is.null(y[1]))m<-cbind(x,y)
m<-elimna(m)
if(!RAN){
if(!MC)temp<-outpro(m,gval=gval,plotit=plotit,op=op,cop=cop,MM=MM,
xlab=xlab,ylab=ylab,STAND=STAND,pr=pr)$keep
if(MC)temp<-outproMC(m,gval=gval,plotit=plotit,op=op,cop=cop,MM=MM,
xlab=xlab,ylab=ylab,STAND=STAND,pr=pr)$keep
}
if(RAN)temp=outpro.depth(m,MM=MM,plotit=plotit)$keep
tcor<-corfun(m[temp,])$cor
if(!is.null(dim((tcor))))tcor<-tcor[1,2]
test<-abs(tcor*sqrt((nrow(m)-2)/(1-tcor**2)))
if(ncol(m)!=2)diag(test)<-NA
crit<-6.947/nrow(m)+2.3197
if(SEED) {
    assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
}
list(cor=tcor,test.stat=test,crit.05=crit)
}

ogkcor<-function(x,y=NA,n.iter=1,sigmamu=taulc,v=gkcov,beta=.9,...){
#
# Compute robust (weighted) correlation matrix in Maronna and Zamar
# (2002, Technometrics, eq. 7).
#
# n.iter number of iterations. 1 seems to be best
# sigmamu computes a robust measure of location and scale for
#  data stored in a single vector.
#  v robust correlation coefficient
#  estloc, a robust measure of location
#
if(!is.na(y[1]))x<-cbind(x,y)
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)
n<-nrow(x)
p<-ncol(x)
val<-matrix(NA,p,p)
temp<-ogk(x,sigmamu=sigmamu,v=v,n.iter=n.iter,beta=beta,...)$cov
J=(p^2-p)/2
p.values=matrix(NA,nrow=J,ncol=3)
info=matrix(NA,nrow=J,ncol=4)
dimnames(p.values)<-list(NULL,c("VAR","VAR","p.value"))
dimnames(info)<-list(NULL,c("VAR","VAR","COR","Test.Stat"))
ic=0
for(j in 1:p){
for(k in 1:p){
val[j,k]<-temp[j,k]/sqrt(temp[k,k]*temp[j,j])
test.stat<-abs(val*sqrt((n-2)/(1-val^2)))
if(j<k){
test=test.stat[j,k]
ic=ic+1
p.values[ic,1]=j
p.values[ic,2]=k
info[ic,1]=j
info[ic,2]=k
info[ic,3]=val[j,k]
info[ic,4]=test
p.value=c("Greater than .1")
crit<-4.8/n+2.72
if(test>=crit)p.value<-c("Less than .1")
crit<-15.49/n+2.68
if(test>=crit)p.value<-c("Less than .05")
crit<-14.22/n+3.26
if(test>=crit)p.value<-c("Less than .025")
crit<-24.83/n+3.74
if(test>=crit)p.value<-c("Less than .01")
p.values[ic,3]=p.value
}}}
list(cor=val,test.results=info,p.values=p.values)
}

pcorhc4sub<-function(x,y,CN=FALSE){
#
#   Compute a .95 confidence interval for Pearson's correlation coefficient.
#   using the HC4 method
#
# CN=T degrees of freedom are infinite, as done by Cribari-Neto (2004)
# CN=F degrees of freedom are n-p
#
xy<-elimna(cbind(x,y))
x<-xy[,1]
y<-xy[,2]
z1=(x-mean(x))/sqrt(var(x))
z2=(y-mean(y))/sqrt(var(y))
ans=olshc4sub(z1,z2,CN=CN)
ci=ans$ci[2,3:4]
ci
}

pcorhc4<-function(x,y,alpha=.05,CN=FALSE,HC3=FALSE){
#
#   Compute a .95 confidence interval for Pearson's correlation coefficient.
#   using the HC4 method
#
# CN=F, degrees of freedom are n-p; seems better for general use.
# CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
#
print('Can return meaningless confidence interval due to outliers')
xy<-elimna(cbind(x,y))
x<-xy[,1]
y<-xy[,2]
z1=(x-mean(x))/sqrt(var(x))
z2=(y-mean(y))/sqrt(var(y))
ans=olshc4(z1,z2,alpha=alpha,CN=CN,HC3=HC3)
list(r=ans$ci[2,2],ci=ans$ci[2,3:4],p.value=ans$ci[2,5],test.stat=ans$test.stat)
}

smcorcom<-function(x1,y1,x2,y2,nboot=200,pts=NA,plotit=TRUE,
SEED=TRUE,varfun=pbvar,xout=TRUE,outfun=out,...){
#
# Compare strength of association of pairs of variables associated with
# two independent  group.
# The strength of the association is based on Cleveland's LOWESS
# smoother coupled with a  robust analog of explanatory power.
#
# The method generalizes the goal of compared the usual
# coefficient of determination associated with two independent groups.
#
#  Assume data are in x1 y1 x2 and y2
#
# Reject at the .05 level if the reported p-value is less than or
# equal to p.crit, which is returned by the function.
#
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
}
m<-elimna(cbind(x2,y2))
x2<-m[,1]
y2<-m[,2]
if(xout){
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
if(SEED)set.seed(2)
estmat1=NA
estmat2=NA
data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
#
for(ib in 1:nboot){
estmat1[ib]=lplot(x1[data1[ib,]],y1[data1[ib,]],plotit=FALSE,
varfun=varfun)$Explanatory.power
estmat2[ib]=lplot(x2[data2[ib,]],y2[data2[ib,]],
varfun=varfun,plotit=FALSE)$Explanatory.power
}
dif<-(estmat1<estmat2)
dif0<-(estmat1==estmat2)
p.value=mean(dif)+.5*mean(dif0)
p.value=2*min(c(p.value,1-p.value))
n1=length(y1)
n2=length(y2)
p1=.05
p2=.05
temp1=tsreg(c(100,200),c(.08,.05))$coef
temp2=tsreg(c(50,100),c(.21,.08))$coef
temp3=tsreg(c(30,50),c(.3,.21))$coef
if(n1<200)p1=temp1[1]+temp1[2]*n1
if(n1<100)p1=temp2[1]+temp2[2]*n1
if(n1<50)p1=temp3[1]+temp3[2]*n1
if(n1<30)p1=.3
if(n2<200)p2=temp1[1]+temp1[2]*n2
if(n2<100)p2=temp2[1]+temp2[2]*n2
if(n2<50)p2=temp3[1]+temp3[2]*n2
if(n2<30)p2=.3
pcrit=(n2*p1+n1*p2)/(n1+n2)
names(pcrit)=NULL
if(plotit)lplot2g(x1,y1,x2,y2)
list(p.value=p.value,pcrit.05=pcrit)
}

pcorbv4<-function(x,y,SEED=TRUE){
#   Compute a .95 confidence interval for Pearson's correlation coefficient.
#
#   This function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
#
#   WARNING: If the number of boostrap samples is altered, it is
#   unknown how to adjust the confidence interval when n < 250.
#   (An obvious guess seems to work well, but no formal investigations
#    have been performed.)
#
nboot<-599  #Number of bootstrap samples
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
xy<-elimna(cbind(x,y))
x<-xy[,1]
y<-xy[,2]
#print("Taking bootstrap samples; please wait")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,pcorbsub,x,y) # A 1 by nboot matrix.
ilow<-13
ihi<-586
if(length(y) < 250){
ilow<-12
ihi<-587
}
if(length(y) < 180){
ilow<-11
ihi<-588
}
if(length(y) < 80){
ilow<-8
ihi<-592
}
if(length(y) < 40){
ilow<-7
ihi<-593
}
bsort<-sort(bvec)
r<-cor(x,y)
ci<-c(bsort[ilow],bsort[ihi])
list(r=r,ci=ci)
}

corbMC<-function(x,y,corfun=pbcor,nboot=599,alpha=.05,SEED=TRUE,...){
#
#   Compute a .95 confidence interval for a correlation.
#   The default correlation is the percentage bend.
#
#   The function corfun is any R function that returns a
#   correlation coefficient in corfun$cor. The functions pbcor and
#   wincor follow this convention.
#
#   When using Pearson's correlation, and when n<250, use
#   lsfitci instead.
#
#   The default number of bootstrap samples is nboot=599
#
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1]
y<-m1[,2]
est<-corfun(x,y,...)$cor
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec<-mclapply(data,corbsub,x,y,corfun,...)
bvec=matl(bvec)   # A 1 by nboot matrix.
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
bsort<-sort(bvec)
corci<-1
corci[1]<-bsort[ilow]
corci[2]<-bsort[ihi]
phat <- sum(bvec < 0)/nboot
sig <- 2 * min(phat, 1 - phat)
list(cor.ci=corci,p.value=sig,cor.est=est)
}

scorci<-function(x,y,nboot=1000,alpha=.05,V2=TRUE,SEED=TRUE,plotit=TRUE,STAND=TRUE,
corfun=pcor,pr=TRUE,cop=3,RAN=FALSE,...){
#
#   Compute a 1-alpha confidence interval for the skipped correlation.
#   alpha=0.05 is the default.
#   By default, Pearson's correlation is computed after outliers are removed via the R function outdoor
#   corfun=spear, for example would replace Pearson's correlation with Spearman's correlation.
#
#   The default number of bootstrap samples is nboot=1000
#
if(pr){
print('As of Sept. 4, 2019, an improved version of this function is used when n<120. To use the old version, set V2=FALSE')
}
if(ncol(as.matrix(x))!=1)stop('x should be a single vector')
if(!V2){
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1]
y<-m1[,2]
est<-scor(x,y,plotit=plotit,STAND=STAND,corfun=corfun,SEED=SEED,cop=cop,pr=FALSE,RAN=RAN,...)$cor
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec<-lapply(data,scorsubMC,x,y,STAND=STAND,corfun=corfun,cop=cop,RAN=RAN,...)
bvec=matl(bvec)   # A 1 by nboot matrix.
bvec=as.vector(bvec)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
bsort<-sort(bvec)
corci<-1
corci[1]<-bsort[ilow]
corci[2]<-bsort[ihi]
phat <- sum(bvec < 0)/nboot
sig <- 2 * min(phat, 1 - phat)
}
if(V2){
a=scorregciH(x,y,nboot=nboot,alpha=alpha,pr=FALSE,SEED=SEED,STOP=FALSE)
est=a$Estimates[1]
sig=a$Estimates[2]
corci=a$confidence.int[2:3]
chk=sign(corci[1]*corci[2])
if(chk>0){
if(sig>alpha)sig=a$Estimates[2]
if(sig>alpha)sig=.95*alpha
}}
list(cor.ci=corci,p.value=sig,cor.est=est)
}

scorsubMC<-function(isub,x,y,pr=FALSE,STAND=TRUE,corfun=corfun,cop=cop,CPP=FALSE,RAN=FALSE,...){
isub=as.vector(isub)
if(!CPP)corbsub<-scor(x[isub],y[isub],plotit=FALSE,pr=FALSE,STAND=STAND,corfun=corfun,cop=cop,
SEED=FALSE,RAN=RAN,...)$cor
if(CPP)stop('Need to use RStudio with WRScpp installed and use the file WRSC++')
corbsub
}

qcorp1<-function(x,y,qest=hd,q=.5,xout=FALSE,outfun=outpro,plotit=FALSE,...){
#
# Compute a measure of the strength of the association
# based on the quantile regression lines
#
X=cbind(x,y)
X=elimna(X)
x<-as.matrix(x)
p=ncol(x)
x=X[,1:p]
p1=p+1
y=X[,p1]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
X=cbind(x,y)
}
est=qreg(x,y,q=q)$coef
top=qreg.sub(X,est,qval=q)
null=qest(y,q)
v=c(null,rep(0,p))
bot=qreg.sub(X,v,qval=q)
ce=sqrt(1-top/bot)
if(p==1)ce=sign(est[2])*ce
list(cor=ce)
}

scorciMC<-function(x,y,nboot=1000,alpha=.05,V2=TRUE,SEED=TRUE,plotit=TRUE,STAND=TRUE,corfun=pcor,pr=TRUE,cop=3,...){
#
#   Compute a 1-alpha confidence interval for the skipped correlation.
#   alpha=0.05 is the default.
#   By default, Pearson's correlation is computed after outliers are removed via the R function outdoor
#   corfun=spear, for example would replace Pearson's correlation with Spearman's correlation.
#
#   The default number of bootstrap samples is nboot=1000
#
#   This function uses the R package parallel
#
if(pr){
print('As of Sept. 4, 2019, an improved version of this function is used when n<120. To use the old version, set V2=FALSE')
}
if(ncol(as.matrix(x))!=1)stop('x should be a single vector')
if(!V2){
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1]
y<-m1[,2]
est<-scor(x,y,plotit=plotit,STAND=STAND,corfun=corfun,SEED=SEED,cop=cop,pr=FALSE,...)$cor
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec<-mclapply(data,scorsubMC,x,y,STAND=STAND,corfun=corfun,cop=cop,...)
bvec=matl(bvec)   # A 1 by nboot matrix.
bvec=as.vector(bvec)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
bsort<-sort(bvec)
corci<-1
corci[1]<-bsort[ilow]
corci[2]<-bsort[ihi]
phat <- sum(bvec < 0)/nboot
sig <- 2 * min(phat, 1 - phat)
}
if(V2){
a=scorregciH(x,y,nboot=nboot,alpha=alpha,pr=FALSE,SEED=SEED,STOP=FALSE)
est=a$Estimates[1]
sig=a$Estimates[2]
corci=a$confidence.int[2:3]
chk=sign(corci[1]*corci[2])
if(chk>0){
if(sig>alpha)sig=a$Estimates[2]
if(sig>alpha)sig=.95*alpha
}}
list(cor.ci=corci,p.value=sig,cor.est=est)
}

tauci<-function(x,y,nboot=1000,alpha=.05,SEED=TRUE,MC=FALSE){
if(!MC)res=corb(x,y,corfun=tau,nboot=nboot,alpha=alpha,SEED=SEED)
if(MC)res=corbMC(x,y,corfun=tau,nboot=nboot,alpha=alpha,SEED=SEED)
res
}

tscor<-function(x,y,xout = FALSE, outfun = out, varfun = winvar,
WARN = TRUE, HD = FALSE, ...){
#
#  Correlation coefficient (explanatory measure of association)
#  based on the Theil--Sen estimator
#
#   To get a p.value, use the R function corb
#
temp=tsreg(x,y,varfun=varfun,xout=xout,outfun=outfun,HD=HD)
val=sign(temp$coef[2])*temp$Strength.Assoc
list(cor=val)
}

tscorci<-function(x,y,nboot=599,alpha=.05,SEED=TRUE,MC=FALSE){
if(!MC)res=corb(x,y,corfun=tscor,nboot=nboot,alpha=alpha,SEED=SEED)
if(MC)res=corbMC(x,y,corfun=tscor,nboot=nboot,alpha=alpha,SEED=SEED)
res
}

wincorci<-function(x,y,nboot=1000,alpha=.05,SEED=TRUE,MC=FALSE,tr=0.2){
if(!MC)res=corb(x,y,corfun=wincor,nboot=nboot,alpha=alpha,SEED=SEED,tr=tr)
if(MC)res=corbMC(x,y,corfun=wincor,nboot=nboot,alpha=alpha,SEED=SEED,tr=tr)
res
}

qcor<-function(x,y,q=.5,qfun=qest,xout=FALSE,outfun=outpro){
#
# Compute quantile correlation as in Li, Li and Tsai, JASA 2015
#
if(xout){
flag<-outfun(x,plotit=plotit)$keep
x<-x[flag]
y<-y[flag]
}
dif=y-qfun(x,q)
flag=dif<0
psi=q-flag
qcov=mean(psi*(x-mean(x)))
qc=qcov/sqrt((q-q^2)*var(x))
list(cor=qc,cov=qcov)
}

corCOMmcp<-function(x,y,corfun=wincor,alpha=.05,nboot=500,SEED=TRUE,MC=FALSE,xout=FALSE,outfun=outpro,method='hommel',...){
#
# Comparing robust dependent correlations: Overlapping case
# That is, have two or more independent variables, compare
# cor(y,x_j) to cor(y,x_k) for each j<k
# Winsorized correlation is used by default.
# Hommel's method is used to control FWE.
#
# x is assumed to be a matrix
#
#
if(nrow(x)!=length(y))stop('x and y have different sample sizes; should be equal')
p=ncol(x)
p1=p+1
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1:p]
y=m1[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
est=NA
for(j in 1:p)est[j]=corfun(x[,j],y)$cor
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
bvec<-mclapply(data,corCOMmcp_sub,x,y,corfun,...)
}
if(!MC)bvec<-lapply(data,corCOMmcp_sub,x,y,corfun,...)
K=(p^2-p)/2
output=matrix(NA,nrow=K,ncol=8)
dimnames(output)=list(NULL,c('IV','IV','Est.1','Est.2','ci.low','ci.hi','p.value','adj.p.value'))
mat=matrix(NA,nrow=nboot,ncol=p)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
for(i in 1:nboot)mat[i,]=bvec[[i]]
ic=0
for(j in 1:p){
for(k in 1:p){
if(j<k){
ic=ic+1
output[ic,1]=j
output[ic,2]=k
output[ic,3]=est[j]
output[ic,4]=est[k]
bsort<-sort(mat[,j]-mat[,k])
output[ic,5]=bsort[ilow]
output[ic,6]=bsort[ihi]
pv=mean(bsort<0)+.5*mean(bsort==0)
output[ic,7]=2*min(c(pv,1-pv))
}}}
output[,8]=p.adjust(output[,7],method=method)
list(results=output)
}

regcor<-function(x,y,regfun=winreg,corfun=wincor,xout=FALSE,outfun=outpro,
plotit=FALSE,xlab='Y',ylab='Y.hat',...){
#
# Estimate strength of association based on some robust regression
# estimator
#
yhat=regYhat(x,y,regfun=regfun,xout=xout,outfun=outfun,,...)
est=corfun(yhat,y)
if(plotit)plot(y,yhat,xlab=xlab,ylab=ylab)
list(cor=est$cor,cor.sq=est$cor^2)
}

mscorpb<-function(x,corfun=pcor,nboot=500,alpha=0.05,SEED=TRUE,outfun=outpro,pr=TRUE){
#
# For p-variate data, test the hypothesis that the
# skipped correlation is zero for all pairs of variables.
# The probability of one or more Type I errors is indicated by the
# argument
# alpha
#
if(pr)print('Here, the p-value is the smallest alpha value for which one or more hypotheses are rejected')
if(SEED)set.seed(2)
x=elimna(x)
n=nrow(x)
tval=NA
for(i in 1:nboot){
y=apply(x,2,sample,replace=TRUE)
v=mscor(y,corfun=corfun,outfun=outfun)
tval[i]=max(abs(elimna(as.vector(v$test.stat))))
}
crit=hd(tval,q=1-alpha)
test=mscor(x,outfun=outfun)
res=test$test.stat
flag=upper.tri(res)
num.sig=sum(abs(res[flag])>=crit)
mtest=max(abs(res[flag]))
hdPV=optimize(hdpv,interval=c(.001,.999),dat=tval,obs=mtest)
if(is.na(num.sig))num.sig=0
list(n=n,cor=test$cor,test.stats=res,crit.val=crit,
num.sig=num.sig,p.value=1-hdPV$minimum)
}

mscorpbMC<-function(x,corfun=pcor,nboot=500,alpha=0.05,SEED=TRUE,WARN=FALSE,
outfun=outpro,pr=TRUE){
#
# For p-variate data, test the hypothesis that the
# skipped correlation is zero for all pairs of variables.
# The probability of one or more Type I errors is indicated by the
# argument
# alpha
#
if(pr)print('Here, the p-value is the smallest alpha value for which one or more hypotheses are rejected')
if(SEED)set.seed(2)
x=elimna(x)
n=nrow(x)
tval=NA
y=list()
for(i in 1:nboot)y[[i]]=apply(x,2,sample,replace=TRUE)
v=mclapply(y,mscor,corfun=corfun,outfun=outpro)
for(i in 1:nboot)tval[i]=max(abs(elimna(as.vector(v[[i]]$test.stat))))

crit=hd(tval,q=1-alpha)
test=mscor(x)
res=test$test.stat
flag=upper.tri(res)
num.sig=sum(abs(res[flag])>=crit)
mtest=max(abs(res[flag]))
if(!WARN)options(warn=-1)
hdPV=optimize(hdpv,interval=c(.001,.999),dat=tval,obs=mtest)
if(!WARN)options(warn=0)
if(is.na(num.sig))num.sig=0
list(n=n,cor=test$cor,test.stats=res,crit.val=crit,
num.sig=num.sig,p.value=1-hdPV$minimum)
}

mscorci.cr<-function(n,p,iter=500,corfun=pcor,alpha=c(.05,.025,.01),TV=FALSE,SEED=TRUE){
#
# Determine critical p-values for the function mscorci
# Returns the estimate of the distribution of the null minimum p-value
#  plus the critical p-values corresponding to the levels indicated by
#  alpha.
#
if(SEED)set.seed(65)
x=list()
for(i in 1:iter)x[[i]]=rmul(n,p=p)
tval=mclapply(x,mscorci.cr.sub,corfun=corfun,nboot=iter)
tval=list2vec(tval)
crit.p=NA
for(j in 1:length(alpha))crit.p[j]=hd(tval,alpha[j])
if(!TV)tval=NULL
list(crit.p.values=crit.p,tval=tval)
}

mscorci.cr.sub<-function(x,corfun,nboot=500){
v=mscorci(x,SEED=FALSE,corfun=corfun,nboot=nboot,crit.pv=1)$p.values
mp=min(as.vector(v),na.rm=T)
mp
}

scorv2<-function(x,y=NULL,corfun=pcor,gval=NA,plotit=FALSE,op=TRUE,cop=3,xlab="VAR 1",
ylab="VAR 2",STAND=TRUE,pr=TRUE,SEED=TRUE,MC=FALSE){
#
# Compute a skipped correlation coefficient.
#
# Eliminate outliers using a projection method
# That is, compute Donoho-Gasko median, for each point
# consider the line between it and the median,
# project all points onto this line, and
# check for outliers using a boxplot rule.
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# For information about the argument cop, see the function
# outpro.
#
# Eliminate any outliers and compute correlation using
# remaining data.
#
#   Nearly the same as scor, but does not reset the SEED, which corrects problems with other functions
#
#  MC=TRUE, the multicore version of outpro is used
#
# corfun=pcor means Pearson's correlation is used.
# corfun=spear means Spearman's correlation is used.
# corfun=tau means Kendall tau is used.
if(SEED){
set.seed(12) # So when using MVE or MCD, get consistent results
}
if(is.null(y[1]))m<-x
if(!is.null(y[1]))m<-cbind(x,y)
m<-elimna(m)
if(!MC)temp<-outpro(m,gval=gval,plotit=plotit,op=op,cop=cop,
xlab=xlab,ylab=ylab,STAND=STAND,pr=pr)$keep
if(MC)temp<-outproMC(m,gval=gval,plotit=plotit,op=op,cop=cop,
xlab=xlab,ylab=ylab,STAND=STAND,pr=pr)$keep
tcor<-corfun(m[temp,])$cor
if(!is.null(dim((m))))tcor<-tcor[1,2]
test<-abs(tcor*sqrt((nrow(m)-2)/(1-tcor**2)))
if(ncol(m)!=2)diag(test)<-NA
crit<-6.947/nrow(m)+2.3197
list(cor=tcor,test.stat=test,crit.05=crit)
}

scorreg<-function(x,y,corfun=spear,cop=3,MM=FALSE,gval=NA,
outfun=outpro,alpha=.05,MC=NULL,SEED=TRUE,ALL=TRUE,...){
#
# x is an n by p matrix
#
# Compute a skipped correlation matrix between y and each variable in x.
#
#  corfun indicates the correlation to be used
#  corfun=pcor uses Pearson's correlation
#  corfun=spear uses Spearman's correlation
#
#  ALL=TRUE: eliminate all outliers among cbind(x,y)
#  ALL=FALSE: skipped correlation is computed for each x[,j] and y. So outliers are eliminated only
#  for these  two variables and this done for j=1,...p, p=number of predictors.
#
# This function returns the p by p matrix of correlations
#
# Method: Eliminate outliers using a projection technique.
# That is, compute Donoho-Gasko median, for each point
# consider the line between it and the median,
# project all points onto this line, and
# check for outliers using a boxplot rule.
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# cop determines how center of the scatterplot is
# estimated; see the function outpro.
# cop=l Donoho-Gasko halfspace median
# cop=2 MCD measure of location
# cop=3 marginal medians
# cop=4 MVE measure of location
#
# gval is critical value for determining whether a point
# is an outlier. It is determined automatically if not specified,
# assuming that Spearman's correlation is used. Critical
# values when using some other correlation have not been
# determined.
#
m<-elimna(cbind(x,y))
m=as.matrix(m)
p1<-ncol(m)
p=p1-1
n<-nrow(m)
e=NA
if(!ALL){
if(is.null(MC)){
if(n>=200)MC=TRUE
else MC=FALSE
}
e=NA
for(j in 1:p)e[j]=scorv2(m[,j],m[,p1],MC=MC,corfun=corfun,SEED=SEED)$cor.value
}
if(ALL){
if(n<500)
flag=outpro(m,cop=cop,plotit=FALSE)$keep
else  flag=outpro.depth(m,plotit=FALSE,SEED=SEED)$keep
xy=m[flag,]
for(j in 1:p)e[j]=corfun(xy[,j],xy[,p1],...)$cor
}
list(cor=e)
}

mscorci<-function(x,y=NULL,nboot=1000,alpha=c(.05,.025,.01),SEED=TRUE,
STAND=TRUE,corfun=pcor,outfun=outpro, crit.pv=NULL,
pvals=NULL,hoch=FALSE,iter=500,pval.SEED=TRUE,pr=TRUE){
#
#  For p-variate data, test the hypothesis of a zero skipped correlation for each pair of variables in a manner
#  that controls the probability of one or more Type I errors.
#
#   The function also returns confidence intervals for each of the skipped correlations when hoch=FALSE.
#   alpha=0.05 is the default.
#   By default, Pearson's correlation is computed after outliers are removed via the R function indicated by
#   outfun, which defaults to a projection-type method.
#   corfun=spear, for example would replace Pearson's correlation with Spearman's correlation.
#
#  alpha=c(.05,.025,.01) is the default, meaning that when determining critical p-values, this is done for
#                                     for alpha .05, 0.25 and .01. So can use different alpha values if desired.
#  For other purposes the family wise error (FWE) rate is taken to be
#   alpha[1]=.05 by default. So setting the argument alpha=.01, FWE is taken to be .01 and a critical p-value is
#   computed for the.01 level only.
#
#   The default number of bootstrap samples is
#   nboot=500
#
#  hoch=TRUE is the default in order to reduce execution time.
#   If  n>=60, this might suffice when testing at the 0.05 level. But power might be increased by using
#   hoch=FALSE at the expense of higher execution time.
#
#   If alpha is less than .05, say .025 or .01,  hoch=FALSE is recommended.
#
#   Note: confidence intervals are reported only when hoch=FALSE.
#
#   pvals can be used to supply a vector of p-values estimating the distribution of the minimum p-value among the tests that are
#   are performed when all hypotheses are true.
#
#  iter=500: number of replications used to estimate the distribution of the minimum p-value.
# Or use the argument crit.pv as indicated below.
# Note: in the journal article dealing with this method, iter=1000 was used.

#  By default
#   pvals=NULL, the functions computes these values if the p-values suggest that there might be
#   significant results and hoch=FALSE; this can result in high execution time.
#   The pvals are computed via the R function
#   mscorci.cr(n,p,iter=500,corfun=pcor,alpha=alpha,SEED=TRUE).
#
#  Critical p-values are a function of n and p. Once known, can supply them via the argument
#  crit.pv as follows:
#
#   pv=mscorci.cr(n,p)$crit.p.values
#   mscorci(x,crit.pv=pv)
#
#
#
if(pr){
if(!hoch){print('To reduce execution time, critical p-values are not computed when the observed p.values are too large to')
print('reject at the 0.05 level. To compute them any way, use the R function mscorci.cr')
}
if(hoch){
print('Hochberg adjusted p-values are used.')
print('This is reasonable when n>120 and alpha=.05. Otherwise suggest using hoch=FALSE')
print('With hoch=TRUE, unadjusted 1-alphaj[1] confidence intervals are reported')
}}
if(SEED)set.seed(2)
if(!is.null(y))x=cbind(x,y)
x<-elimna(x)  # Eliminate rows with missing values
nval=nrow(x)
p=ncol(x)
J=(p^2-p)/2
est<-mscor(x,STAND=STAND,corfun=corfun,outfun=outfun)$cor
flag=upper.tri(est)
est=est[flag]
data<-matrix(sample(nval,size=nval*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec<-lapply(data,scorci.sub,x,STAND=STAND,corfun=corfun,outfun=outfun)
bvec=matl(bvec)   # A J by nboot matrix.

phat=0
sig=0
for(j in 1:J){
phat[j]=sum(bvec[j,] < 0)/nboot
sig[j] <- 2 * min(phat[j], 1 - phat[j])
}
#
# Compute critical p-values if any of the p-values are sufficiently small.
#
FLAG=FALSE
if(p==2 && sig[1]<=.15){
FLAG=TRUE
if(hoch){
if(pr)print('If the p.value is <=.15, suggest using hoch=FALSE')
}}
if(p>2){
if(min(sig)<=alpha[1])FLAG=TRUE
}
if(FLAG){
if(is.null(crit.pv)){
if(!hoch){
if(pr)print('Computing critical p-values. Execution time might require several minutes')
temp=mscorci.cr(nval,p,iter=iter,corfun=corfun,alpha=alpha,SEED=pval.SEED,TV=TRUE) #returns tval in case want to adjust p-values.
#                             Need to add code to do this. (See mscorpbMC for how this might be done.)
crit.pv=temp$crit.p.values
}}}
ci.mat=matrix(NA,nrow=J,ncol=4)
dimnames(ci.mat)=list(NULL,c('Var i','Var j','ci.low','ci.up'))
for(j in 1:J)bvec[j,]<-sort(bvec[j,])
if(J==1)bvec=as.matrix(bvec)
ic=0
if(is.null(crit.pv))crit.pv=alpha[1]
for(j in 1:p){
for(k in 1:p){
if(j<k){
ic=ic+1
ci.mat[ic,1]=j
ci.mat[ic,2]=k
ihi<-floor((1-crit.pv[1]/2)*nboot+.5)
ilow<-floor((crit.pv[1]/2)*nboot+.5)
ci.mat[ic,3]=bvec[ic,ilow]
ci.mat[ic,4]=bvec[ic,ihi]
}}}

p.mat=matrix(NA,nrow=p,ncol=p)
p.mat[flag]=sig
adj.p=NULL
if(hoch)adj.p=p.adjust(sig,method='hochberg')
list(p.values=p.mat,cor.est=est,confidence.int=ci.mat,critical.p.values=crit.pv,hoch.adjusted.p.values=adj.p)
}

scorci.sub<-function(data,x,STAND=STAND,corfun=corfun,outfun=outfun){
est<-mscor(x[data,],STAND=STAND,corfun=corfun)$cor
flag=upper.tri(est)
est=est[flag]
est
}

scorregciH<-function(x,y,nboot=500,alpha=.05,SEED=TRUE,
corfun=pcor,outfun=outpro, crit.pv=NULL,ALL=TRUE,MC=TRUE,
pvals=NULL,iter=500,pval.SEED=TRUE,pr=TRUE,STOP=TRUE){
#
#  For explanatory variables, test the hypothesis of a zero skipped correlation between y  and  each
#   explanatory variable in a manner
#  that controls the probability of one or more Type I errors.
#
#   Use Hochberg adjusted critical p-values based on adjusted p-values when when there is p=1 predictor.
#  This controls, approximately, the probability of one or more Type I errors.
#
#   By default, Pearson's correlation is computed after outliers are removed via the R function indicated by
#   outfun, which defaults to a projection-type method.
#   corfun=spear, for example would replace Pearson's correlation with Spearman's correlation.
#
#  alpha=.05 is the default, meaning that for each correlation, .95 confidence intervals are computed.
# That is, the simultaneous probability is not .05.
#
#   The default number of bootstrap samples is
#   nboot=500
#
#
if(pr){
print('Hochberg adjusted p-values are used that are designed so that FWE is approximately  alpha.')
print('Each confidence interval has, approximately, 1-alpha probability coverage')
print('So it is possible that Hochberg fails to rejects when confidence intervals do not')
}
if(SEED)set.seed(2)
xy=elimna(cbind(x,y))
x=as.matrix(x)
p=ncol(x)
if(p==1){
if(STOP)stop('With a single independent variable, use scorci')
}
p1=p+1
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
n=nrow(x)
est<-scorreg(x,y,corfun=corfun,outfun=outfun,SEED=FALSE,ALL=ALL)$cor
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(!MC)bvec<-lapply(data,scorreg.sub,xy,corfun=corfun,outfun=outfun,ALL=ALL)
if(MC){
bvec<-mclapply(data,scorreg.sub,xy,corfun=corfun,outfun=outfun,ALL=ALL)
}
bvec=matl(bvec)   # A p by nboot matrix.
phat=0
sig=0
sigadj=0
for(j in 1:p){
phat[j]=sum(bvec[j,] < 0)/nboot
sig[j]=2 * min(phat[j], 1 - phat[j])
if(n<=40){vv=p.crit.n30(alpha[1],sig[j])
# Next, adjust the p-values
crit.p=vv$crit.p.value
sigadj[j]=vv$adj.p.value
}
if(n>40){
if(n<=70){
vv=p.crit.n60(alpha[1],sig[j])
sigadj[j]=vv$adj.p.value
crit.p=vv$crit.p.value
}}
if(n>70){
if(n<=100){
vv=p.crit.n80(alpha[1],sig[j])
sigadj[j]=vv$adj.p.value
crit.p=vv$crit.p.value
}
}
if(n>100){
if(n<=120)
{
vv=p.crit.n100(alpha[1],sig[j])
crit.p=vv$crit.p.value
sigadj[j]=vv$adj.p.value
}}
if(n>120){  # no adjustment
sigadj[j]=sig[j]  #i.e., no adjustment
crit.p=alpha[1]
}}
hadj=p.adjust(sig,method='hoch')
ci.mat=matrix(NA,nrow=p,ncol=3)
dimnames(ci.mat)=list(NULL,c('Var','ci.low','ci.up'))
for(j in 1:p)bvec[j,]<-sort(bvec[j,])
if(p==1)bvec=as.matrix(bvec)
ic=0
if(is.null(crit.pv))crit.pv=alpha[1]
for(j in 1:p){
ic=ic+1
ci.mat[ic,1]=j
ihi<-floor((1-crit.p[1]/2)*nboot+.5)
ilow<-floor((crit.p[1]/2)*nboot+.5)
ci.mat[ic,2]=bvec[ic,ilow]
ci.mat[ic,3]=bvec[ic,ihi]
}
p.mat=matrix(NA,nrow=p,ncol=4)
p.mat[,1]=est
p.mat[,2]=sig
p.mat[,3]=sigadj
adj.p=p.adjust(sigadj,method='hochberg')
p.mat[,4]=adj.p
dimnames(p.mat)=list(NULL,c('Est.','p-value','adjusted p.value','Hoch.adjusted.p.value'))
list(Estimates=p.mat,confidence.int=ci.mat)
}

scorreg.sub<-function(data,xy,corfun=corfun,outfun=outfun,ALL=ALL){
p1=ncol(xy)
p=p1-1
est<-scorreg(xy[data,1:p],xy[data,p1],corfun=corfun,SEED=FALSE,ALL=ALL)$cor
est
}

mscorciH<-function(x,nboot=1000,alpha=.05,SEED=TRUE,method='hoch',
corfun=pcor,outfun=outpro, crit.pv=NULL,ALL=TRUE,MC=TRUE,pr=TRUE){
#
#  Test the hypothesis of a zero skipped correlation for each pair of variables in
#  x, an n-by-p matrix.
#
#   Use Hochberg adjusted critical p-values based on adjusted p-values to control the probability of one or more Type I errors.
#
#   The function also returns 1-alpha confidence intervals for each of the skipped correlations
#   alpha=0.05 is the default.
#   By default, Pearson's correlation is computed after outliers are removed via the R function indicated by
#   outfun, which defaults to a projection-type method.
#   corfun=spear, for example would replace Pearson's correlation with Spearman's correlation.
#
#   The default number of bootstrap samples is
#   nboot=500
#
#
if(pr){
print('Each confidence interval has, approximately, 1-alpha probability coverage')
}
if(SEED)set.seed(2)
xy=elimna(x)
x=as.matrix(x)
p=ncol(x)
p1=p+1
J=(p^2-p)/2
x=as.matrix(x)
n=nrow(x)
est<-mscor(x,corfun=corfun,outfun=outfun)$cor
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(!MC)bvec<-lapply(data,scorci.sub,x,corfun=corfun,outfun=outfun,STAND=TRUE)
if(MC){
bvec<-mclapply(data,scorci.sub,x,corfun=corfun,outfun=outfun,STAND=TRUE)
}
bvec=matl(bvec)   # A J by nboot matrix. J=(p^2-p)/2, p=number of IV variables.

#
# Compute critical p-values when p=1 and then use Hochberg adjustment
#

phat=0
sig=matrix(NA,p,p)
sigadj=matrix(NA,p,p)
ic=0
for(j in 1:p){
for(k in 1:p){
if(j<k){
ic=ic+1
phat=sum(bvec[ic,] < 0)/nboot
sig[j,k]=2 * min(phat, 1 - phat)
if(n<=40){
vv=p.crit.n30(alpha[1],sig[j,k])
crit.p=vv$crit.p.value
sigadj[j,k]=vv$adj.p.value  # Adjust each p-value when sample size is small.
}
if(n>40){
if(n<=70){
vv=p.crit.n60(alpha[1],sig[j,k])
sigadj[j,k]=vv$adj.p.value
crit.p=vv$crit.p.value
}
}
if(n>70){
if(n<=100){
vv=p.crit.n80(alpha[1],sig[j,k])
sigadj[j,k]=vv$adj.p.value
crit.p=vv$crit.p.value
}
}
if(n>100){
if(n<=120)
{
vv=p.crit.n100(alpha[1],sig[j,k])
crit.p=vv$crit.p.value
sigadj[j,k]=vv$adj.p.value
}}
if(n>120){  # no adjustment
sigad[j,k]=sig[j,k]  #i.e., no adjustment
crit.p=alpha
}}}}
ci.mat=matrix(NA,nrow=J,ncol=7)
dimnames(ci.mat)=list(NULL,c('Var i','Var j','Est','ci.low','ci.up','P-value','FWE Adjusted p-value' ))
crit.pv=crit.p
for(j in 1:J)bvec[j,]<-sort(bvec[j,])
if(J==1)bvec=as.matrix(bvec)
ic=0
if(is.null(crit.pv))crit.pv=alpha[1]
for(j in 1:p){
for(k in 1:p){
if(j<k){
ic=ic+1
ci.mat[ic,1]=j
ci.mat[ic,2]=k
#ihi<-floor((1-alpha/2)*nboot+.5)
ihi<-floor((1-crit.pv[1]/2)*nboot+.5)
#ilow<-floor((alpha/2)*nboot+.5)
ilow<-floor((crit.pv[1]/2)*nboot+.5)
ci.mat[ic,3]=est[j,k]
ci.mat[ic,4]=bvec[ic,ilow]
ci.mat[ic,5]=bvec[ic,ihi]
ci.mat[ic,6]=sigadj[j,k]
}}}
flag=upper.tri(sigadj)
p.mat=matrix(NA,nrow=p,ncol=p)
p.mat[flag]=p.adjust(sigadj[flag],method=method)
ic=0
for(j in 1:p){
for(k in 1:p){
if(j<k){
ic=ic+1
ci.mat[ic,7]=p.mat[j,k]
}}}
list(output=ci.mat,critical.p.value=crit.pv)
}

scorregci<-function(x,y,nboot=500,alpha=c(.05,.025,.01),SEED=TRUE,
corfun=pcor,outfun=outpro, crit.pv=NULL,ALL=TRUE,MC=TRUE,
pvals=NULL,hoch=TRUE,iter=500,pval.SEED=TRUE,pr=TRUE,...){
#
#  For explanatory variables, test the hypothesis of a zero skipped correlation between y  and  each
#   explanatory variable in a manner
#  that controls the probability of one or more Type I errors.
#
#   The function also returns confidence intervals for each of the skipped correlations when hoch=FALSE.
#   alpha=0.05 is the default.
#   By default, Pearson's correlation is computed after outliers are removed via the R function indicated by
#   outfun, which defaults to a projection-type method.
#   corfun=spear, for example would replace Pearson's correlation with Spearman's correlation.
#
#  alpha=c(.05,.025,.01) is the default, meaning that when determining critical p-values, this be done for
#                                     for alpha .05, 0.25 and .01. So can use different alpha values if desired.
#  For other purposes the family wise error (FWE) rate is taken to be
#   alpha[1]=.05 by default. So setting the argument alpha=.01, FWE is taken to be .01 and a critical p-value is
#   computed for the.01 level only.
#
#   The default number of bootstrap samples is
#   nboot=500
#
#  hoch=TRUE is the default in order to reduce execution time.
#   If  n>=60, this might suffice when testing at the 0.05 level. But power might be increased by using
#   hoch=FALSE at the expense of higher execution time.
#
#   If alpha is less than .05, say .025 or .01,  hoch=FALSE is recommended.
#
#   Note: confidence intervals are reported only when hoch=FALSE.
#
#   pvals can be used to supply a vector of p-values estimating the distribution of the minimum p-value among the tests that are
#   are performed when all hypotheses are true.
#
#  iter=500: number of replications used to estimate the distribution of the minimum p-value.
# Or use the argument crit.pv as indicated below.
# Note: in the journal article dealing with this method, iter=1000 was used.

#  By default
#   pvals=NULL, the functions computes these values if the p-values suggest that there might be
#   significant results and hoch=FALSE; this can result in high execution time.
#   The pvals are computed via the R function
#   mscorci.cr(n,p,iter=500,corfun=pcor,alpha=alpha,SEED=TRUE).
#
#  Critical p-values are a function of n and p. Once known, can supply them via the argument
#  crit.pv as follows:
#
#   pv=scorregci.cr(n,p)$crit.p.values
#   scorregci(x,crit.pv=pv)
#
#  When hoch=TRUE, unadjusted confidence intervals are returned.
#
#
#
if(pr){
if(!hoch){print('To reduce execution time, critical p-values are not computed when the observed p.values are too large to')
print('reject at the 0.05 level. To compute them any way, use the R function scorregci.cr')
}
if(hoch){
print('Hochberg adjusted p-values are used.')
print('This is reasonable when n>=60 and alpha=.05. Otherwise suggest using hoch=FALSE')
print('To get confidence intervals, use hoch=FALSE')
}}
if(SEED)set.seed(2)
xy=elimna(cbind(x,y))
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
n=nrow(x)
est<-scorreg(x,y,corfun=corfun,outfun=outfun,SEED=FALSE,ALL=ALL,...)$cor
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(!MC)bvec<-lapply(data,scorreg.sub,xy,corfun=corfun,outfun=outfun,ALL=ALL,...)
if(MC){
bvec<-mclapply(data,scorreg.sub,xy,corfun=corfun,outfun=outfun,ALL=ALL,...)
}
bvec=matl(bvec)   # A p by nboot matrix.

phat=0
sig=0
for(j in 1:p){
phat[j]=sum(bvec[j,] < 0)/nboot
sig[j] <- 2 * min(phat[j], 1 - phat[j])
}
#
# Compute critical p-values if any of the p-values are sufficiently small.
#
FLAG=FALSE
if(p==2 && sig[1]<=.15){
FLAG=TRUE
if(hoch){
if(pr)print('If the p.value is <=.15, suggest using hoch=FALSE')
}}
if(p>2){
if(min(sig)<=alpha[1])FLAG=TRUE
}
if(FLAG){
if(is.null(crit.pv)){
if(!hoch){
if(pr)print('Computing critical p-values. Execution time might require several minutes')
temp=scorregci.cr(nval,p,iter=iter,corfun=corfun,alpha=alpha,SEED=pval.SEED,TV=TRUE,ALL=ALL) #returns tval in case want to adjust p-values.
#                             Need to add code to do this. (See mscorpbMC for how this might be done.)
crit.pv=temp$crit.p.values
}}}
ci.mat=matrix(NA,nrow=p,ncol=3)
dimnames(ci.mat)=list(NULL,c('Var','ci.low','ci.up'))
for(j in 1:p)bvec[j,]<-sort(bvec[j,])
if(p==1)bvec=as.matrix(bvec)
ic=0
if(is.null(crit.pv))crit.pv=alpha[1]
for(j in 1:p){
ic=ic+1
ci.mat[ic,1]=j
ihi<-floor((1-crit.pv[1]/2)*nboot+.5)
ilow<-floor((crit.pv[1]/2)*nboot+.5)
ci.mat[ic,2]=bvec[ic,ilow]
ci.mat[ic,3]=bvec[ic,ihi]
}

p.mat=matrix(NA,nrow=p,ncol=3)
p.mat[,1]=est
p.mat[,2]=sig
adj.p=NULL
if(hoch){
adj.p=p.adjust(sig,method='hochberg')
p.mat[,3]=adj.p
}
dimnames(p.mat)=list(NULL,c('Est.','p-value','adjusted p.value'))
list(Estimates=p.mat,confidence.int=ci.mat,critical.p.values=crit.pv)
}

scorreg.sub<-function(data,xy,corfun=corfun,outfun=outfun,ALL=ALL,...){
p1=ncol(xy)
p=p1-1
est<-scorreg(xy[data,1:p],xy[data,p1],corfun=corfun,SEED=FALSE,ALL=ALL,...)$cor
est
}

scorreg.cr<-function(n,p,iter=500,nboot=500,corfun=pcor,alpha=c(.05,.025,.01),TV=FALSE,ALL=TRUE,SEED=TRUE,outfun=outpro){
#
# Determine critical p-values for the function scorregci
# Returns the estimate of the distribution of the null minimum p-value
#  plus the critical p-values corresponding to the levels indicated by
#  alpha.
#
#  p = number or predictors
#
#  Function assumes that a multicore processor is used and that the R package parallel has been installed.
#
if(SEED)set.seed(65)
x=list()
p1=p+1
for(i in 1:iter){
x[[i]]=rmul(n,p=p1)
}
tval=mclapply(x,scorreg.cr.sub,p=p,corfun=corfun,nboot=nboot,ALL=ALL)
tval=list2vec(tval)
crit.p=NA
for(j in 1:length(alpha))crit.p[j]=hd(tval,alpha[j])
if(!TV)tval=NULL
list(crit.p.values=crit.p,tval=tval)
}

scorreg.cr.sub<-function(x,corfun,p=p,nboot=500,ALL=ALL,outfun=outfun){
p1=p+1
v=scorregci(x[,1:p],x[,p1],SEED=FALSE,corfun=corfun,nboot=nboot,crit.pv=1,pr=FALSE,hoch=TRUE,ALL=ALL,outfun=outfun)$Estimates[,2]
mp=min(as.vector(v),na.rm=T)
mp
}

scorall<-function(x,outfun=outpro,corfun=pcor,RAN=FALSE,...){
#
# Eliminate outliers and compute a correlation based on the
# remaining data.
#
x=elimna(x)
if(!RAN)flag=outpro(x)$keep
if(RAN)flag=outpro.depth(x)$keep
est=corfun(x[flag,],...)$cor
est
}

corxy<-function(x,y,corfun=wincor,alpha=.05,nboot=599,...){
#
# x is an n by p matrix
#
# Compute a correlation matrix between y and each variable in x.
#
#  corfun indicates the correlation to be used
#  corfun=pbcor uses Percentage bend  correlation
#  corfun=spear uses Spearman's correlation
#  corfun=tau	uses Kendall's tau.
#
#
m<-elimna(cbind(x,y))
m=as.matrix(m)
p1<-ncol(m)
p=p1-1
n<-nrow(m)
e=matrix(NA,nrow=p,ncol=4)
for(j in 1:p){
v=corb(x[,j],y,corfun=corfun,alpha=alpha,nboot=nboot,...)
e[j,1:2]=v$cor.ci
e[j,3]=v$p.value
e[j,4]=v$cor.est
}
dimnames(e)=list(NULL,c('ci.low','ci.hi','p-value','est'))
e
}

perm.rho<-function(x,y,alpha=.05,nboot=1000,SEED=TRUE){
#
# Do a  permutation test based on Pearson's correlation
# Diciccio--Romano version of a permuation test (JASA, 2017, 112, 1211-1220)
#

# The default number of permutations is nboot=1000
#
if(SEED)set.seed(2)
xx<-cbind(x,y)
xx=elimna(xx)
x=xx[,1]
y=xx[,2]
n=length(x)
tval<-perm.rho.sub(x,y)
vec<-c(1:length(xx))
v1<-length(x)+1
difb<-NA
tv<-NA
for(i in 1:nboot){
id=sample(n,n)
tv[i]<-perm.rho.sub(x,y[id])
}
tv<-sort(tv)
icl<-floor((alpha/2)*nboot+.5)
icu<-floor((1-alpha/2)*nboot+.5)
reject<-0
if(tval>=tv[icu] || tval <=tv[icl])reject<-1
list(teststat=tval,lower.crit=tv[icl],upper.crit=tv[icu],reject=reject)
}

rhohc4bt<-function(X1,Y1,alpha=.05,nboot=2999,SEED=TRUE){
#
# A p-value for rhohc4bt
# Pearson's correlations using a bootstrap-t method in conjunction with the HC4 estimator
#
if(SEED)set.seed(2)
print('Can return meaningless confidence interval due to outliers')
XY=cbind(X1,Y1)
XY=elimna(XY)
X1=XY[,1]
Y1=XY[,2]
r1=cor(X1,Y1)
n1=length(X1)
v=NA
CI=NULL
alph=seq(0.01,0.99,.01)
for(i in 1:nboot){
id1=sample(n1,n1,replace=TRUE)
x1=X1[id1]
y1=Y1[id1]
x1=(x1-mean(x1))/sd(x1)
y1=(y1-mean(y1))/sd(y1)
temp1=olshc4(x1,y1)
v[i]=(temp1$ci[2,2]-r1)/temp1$ci[2,6]
}
vs=sort(v)
na=length(alph)
temp1=olshc4(X1,Y1)
for(i in 1:na){
ibot<-round(alph[i]*nboot/2)
itop<-nboot-ibot+1
ibot=ibot+1
low=r1-abs(vs[itop])*temp1$ci[2,6]
up=r1+abs(vs[ibot])*temp1$ci[2,6]
if(sign(low*up)==1)break
}
# Determine CI
bot<-round(alpha[i]*nboot/2)
itop<-nboot-ibot+1
ibot=ibot+1
low=r1-abs(vs[itop])*temp1$ci[2,6]
up=r1+abs(vs[ibot])*temp1$ci[2,6]
low=max(-1,low)
up=min(up,1)
list(est=r1,p.value=alph[i],ci=c(low,up))
}

corCOM.DVvsIV<-function(x,y,com.p.dist=FALSE,corfun=wincor,iter=200,PV=NULL,pr=TRUE,neg.col=NULL,LARGEST=TRUE,
alpha=.05,nboot=500,SEED=TRUE,MC=FALSE,xout=FALSE,outfun=outpro,FWE.method='hoch',...){
#
# Regresiion:
# Consider the IV with the largest correlation estimate with the DV
# Is it reasonable to decide that it has the highest population
# correlation?
#
# That is, have two or more independent variables, compare
# cor(y,x_I) to cor(y,x_k) for all k!=I, where
# cor(i,x_I) is the highest correlation
# Winsorized correlation is used by default.
# Hochberg's method is used to control FWE.
#
# x is assumed to be a matrix or data frame
#
#  Possible alternative choices for corfun include:
#  spear
#  tau
#  pbcor
#  bicor
#  scor
#  mve.cor
#  mcd.cor
#
#
if(nrow(x)!=length(y))stop('x and y have different sample sizes; should be equal')
p=ncol(x)
if(p<=2)stop('Should have 3 or more independent variables. With two, use TWOpov')
pm1=p-1
p1=p+1
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
nval=nrow(x)
x=neg.colM(x,neg.col)
L=c(seq(.001,.1,.001),seq(.11,.99,.01))

if(!is.null(PV)){
rem=matrix(NA,length(L),pm1)
for(k in 1:pm1){
for(i in 1:length(L))rem[i,k]=hd(PV[,k],q=L[i])
}}
if(p>6)com.p.dist=TRUE
if(com.p.dist){
if(is.null(PV)){
if(pr)print('Computing the null distribution can take several minutes')
PV=corCOM.DVvsIV.crit(p,nval,iter=iter,MC=MC,...)
rem=matrix(NA,length(L),pm1)
for(k in 1:pm1){
for(i in 1:length(L))rem[i,k]=hd(PV[,k],q=L[i])
}}}

P3=matrix(c( 0.002175049 ,1.335108e-05
 ,0.005872039 ,7.055924e-05
 ,0.010619227 ,2.497530e-04
 ,0.016286140 ,5.345789e-04
 ,0.022398371 ,7.976283e-04
 ,0.028262446 ,1.277063e-03
 ,0.033345906 ,1.627195e-03
 ,0.037456504 ,1.939087e-03
 ,0.040694218 ,2.234302e-03
 ,0.043307931 ,2.546080e-03
 ,0.045568092 ,2.903338e-03
 ,0.047696567 ,3.323658e-03
 ,0.049845675 ,3.813582e-03
 ,0.052103950 ,4.371974e-03
 ,0.054511258 ,4.993196e-03
 ,0.057074761 ,5.668892e-03
 ,0.059782746 ,6.388620e-03
 ,0.062615440 ,7.140126e-03
 ,0.065552496 ,7.909937e-03
 ,0.068577121 ,7.684452e-03
 ,0.071677317 ,9.451369e-03
 ,0.074845101 ,1.020104e-02
 ,0.078074759 ,1.092742e-02
 ,0.081360990 ,1.162839e-02
 ,0.084697501 ,1.230540e-02
 ,0.088076254 ,1.296274e-02
 ,0.091487319 ,1.360651e-02
 ,0.094919200 ,1.424359e-02
 ,0.098359413 ,1.488081e-02
 ,0.101795143 ,1.552446e-02
 ,0.105213802 ,1.617997e-02
 ,0.108603375 ,1.685188e-02
 ,0.111952481 ,1.754394e-02
 ,0.115250189 ,1.825919e-02
 ,0.118485658 ,1.900005e-02
 ,0.121647762 ,1.976830e-02
 ,0.124724850 ,2.056506e-02
 ,0.127704758 ,2.139068e-02
 ,0.130575140 ,2.224462e-02
 ,0.133324080 ,2.312538e-02
 ,0.135940914 ,2.403049e-02
 ,0.138417122 ,2.495656e-02
 ,0.140747141 ,2.589933e-02
 ,0.142928988 ,2.685388e-02
 ,0.144964592 ,2.781478e-02
 ,0.146859829 ,2.877632e-02
 ,0.148624253 ,2.973273e-02
 ,0.150270597 ,3.067838e-02
 ,0.151814117 ,3.160798e-02
 ,0.153271855 ,3.251674e-02
 ,0.154661905 ,3.340049e-02
 ,0.156002725 ,3.425580e-02
 ,0.157312540 ,3.508002e-02
 ,0.158608833 ,3.587127e-02
 ,0.159907941 ,3.662845e-02
 ,0.161224730 ,3.735119e-02
 ,0.162572340 ,3.803977e-02
 ,0.163961992 ,3.869502e-02
 ,0.165402841 ,3.931829e-02
 ,0.166901885 ,3.991132e-02
 ,0.168463921 ,4.047615e-02
 ,0.170091558 ,4.101511e-02
 ,0.171785292 ,4.153070e-02
 ,0.173543642 ,4.202554e-02
 ,0.175363353 ,4.250234e-02
 ,0.177239642 ,4.296380e-02
 ,0.179166497 ,4.341259e-02
 ,0.181137006 ,4.385123e-02
 ,0.183143686 ,4.428208e-02
 ,0.185178820 ,4.470724e-02
 ,0.187234759 ,4.512852e-02
 ,0.189304201 ,4.554735e-02
 ,0.191380404 ,4.596483e-02
 ,0.193457363 ,4.638162e-02
 ,0.195529910 ,4.679804e-02
 ,0.197593772 ,4.721406e-02
 ,0.199645563 ,4.762936e-02
 ,0.201682735 ,4.804340e-02
 ,0.203703491 ,4.845551e-02
 ,0.205706665 ,4.886496e-02
 ,0.207691581 ,4.927106e-02
 ,0.209657914 ,4.967325e-02
 ,0.211605535 ,5.007115e-02
 ,0.213534385 ,5.046464e-02
 ,0.215444346 ,5.085386e-02
 ,0.217335153 ,5.123927e-02
 ,0.219206331 ,5.162162e-02
 ,0.221057158 ,5.200196e-02
 ,0.222886667 ,5.238156e-02
 ,0.224693687 ,5.276190e-02
 ,0.226476899 ,5.314459e-02
 ,0.228234929 ,5.353130e-02
 ,0.229966452 ,5.392373e-02
 ,0.231670301 ,5.432350e-02
 ,0.233345582 ,5.473211e-02
 ,0.234991770 ,5.515090e-02
 ,0.236608793 ,5.558098e-02
 ,0.238197083 ,5.602325e-02
 ,0.239757602 ,5.647831e-02
 ,0.241291828 ,5.694653e-02
 ,0.255725281 ,6.227133e-02
 ,0.269643155 ,6.812649e-02
 ,0.282572599 ,7.379864e-02
 ,0.294230675 ,7.937555e-02
 ,0.305279532 ,7.560157e-02
 ,0.316231246 ,9.290848e-02
 ,0.327118515 ,1.006091e-01
 ,0.337725949 ,1.077792e-01
 ,0.348507973 ,1.143072e-01
 ,0.360275335 ,1.203813e-01
 ,0.372732081 ,1.259302e-01
 ,0.384865494 ,1.309854e-01
 ,0.396170468 ,1.359841e-01
 ,0.406525702 ,1.414512e-01
 ,0.416179917 ,1.474082e-01
 ,0.425830260 ,1.532244e-01
 ,0.435731315 ,1.584113e-01
 ,0.445574831 ,1.634671e-01
 ,0.455403556 ,1.693949e-01
 ,0.465473322 ,1.764807e-01
 ,0.475493015 ,1.840579e-01
 ,0.485104441 ,1.914634e-01
 ,0.494622601 ,1.987244e-01
 ,0.504479896 ,2.062144e-01
 ,0.514530980 ,2.140417e-01
 ,0.524379809 ,2.218487e-01
 ,0.533811067 ,2.291090e-01
 ,0.542964071 ,2.357491e-01
 ,0.552415670 ,2.422511e-01
 ,0.562702410 ,2.489342e-01
 ,0.573662536 ,2.556239e-01
 ,0.584499138 ,2.621385e-01
 ,0.594424353 ,2.684997e-01
 ,0.603091449 ,2.746540e-01
 ,0.610580104 ,2.804293e-01
 ,0.617291412 ,2.857687e-01
 ,0.623819392 ,2.908378e-01
 ,0.630625766 ,2.959208e-01
 ,0.637780329 ,3.012891e-01
 ,0.645128388 ,3.071430e-01
 ,0.652676811 ,3.135295e-01
 ,0.660706577 ,3.202344e-01
 ,0.669513580 ,3.269345e-01
 ,0.679074179 ,3.335610e-01
 ,0.688924859 ,3.403140e-01
 ,0.698359030 ,3.473033e-01
 ,0.706803246 ,3.544525e-01
 ,0.714138033 ,3.617598e-01
 ,0.720773425 ,3.694019e-01
 ,0.727376023 ,3.775318e-01
 ,0.734427437 ,3.860345e-01
 ,0.741910594 ,3.945174e-01
 ,0.749301181 ,4.026368e-01
 ,0.756038316 ,4.104358e-01
 ,0.762108143 ,4.183084e-01
 ,0.767991314 ,4.266920e-01
 ,0.774125688 ,4.358112e-01
 ,0.780642168 ,4.456474e-01
 ,0.787654735 ,4.560432e-01
 ,0.795560943 ,4.667688e-01
 ,0.804494568 ,4.776577e-01
 ,0.813689728 ,4.888006e-01
 ,0.822246533 ,5.003803e-01
 ,0.830074100 ,5.121099e-01
 ,0.837601422 ,5.231564e-01
 ,0.845257402 ,5.330505e-01
 ,0.853115769 ,5.423039e-01
 ,0.860852876 ,5.518459e-01
 ,0.868387995 ,5.622222e-01
 ,0.875848963 ,5.737404e-01
 ,0.882870692 ,5.869893e-01
 ,0.889442420 ,6.017184e-01
 ,0.896636569 ,6.161664e-01
 ,0.904626820 ,6.294051e-01
 ,0.911772218 ,6.420279e-01
 ,0.917549883 ,6.554093e-01
 ,0.923211670 ,6.714817e-01
 ,0.928746246 ,6.890702e-01
 ,0.934160521 ,7.060012e-01
 ,0.941177134 ,7.215806e-01
 ,0.949363922 ,7.357055e-01
 ,0.956847590 ,7.494285e-01
 ,0.963599180 ,7.621246e-01
 ,0.969055149 ,7.763082e-01
 ,0.973136108 ,7.936805e-01
 ,0.976868136 ,7.142080e-01
 ,0.981258840 ,7.460053e-01
 ,0.988264909 ,7.786117e-01
 ,0.994433844 ,9.192672e-01),
 byrow=TRUE,ncol=2)


P4=matrix(c(0.02825196 ,0.002694346 ,2.999028e-06
 ,0.03650906 ,0.004146146 ,4.263801e-05
 ,0.04206787 ,0.005668266 ,1.964571e-04
 ,0.04755122 ,0.007287032 ,4.847254e-04
 ,0.05291511 ,0.009334544 ,8.277385e-04
 ,0.05813244 ,0.011977936 ,1.156019e-03
 ,0.06348677 ,0.015026277 ,1.458581e-03
 ,0.06908051 ,0.018151721 ,1.739547e-03
 ,0.07468576 ,0.021131256 ,1.999549e-03
 ,0.08003211 ,0.023887028 ,2.244613e-03
 ,0.08506408 ,0.026419261 ,2.483839e-03
 ,0.08992250 ,0.028745172 ,2.721838e-03
 ,0.09476611 ,0.030872849 ,2.960090e-03
 ,0.09964373 ,0.032802489 ,3.201924e-03
 ,0.10449749 ,0.034540459 ,3.452190e-03
 ,0.10923896 ,0.036109254 ,3.712771e-03
 ,0.11381116 ,0.037544213 ,3.980112e-03
 ,0.11820149 ,0.038882089 ,4.247472e-03
 ,0.12242309 ,0.040152826 ,4.509743e-03
 ,0.12649443 ,0.041380182 ,4.767249e-03
 ,0.13043191 ,0.042587785 ,5.026383e-03
 ,0.13425237 ,0.043803447 ,5.297105e-03
 ,0.13797722 ,0.045057120 ,5.588809e-03
 ,0.14163298 ,0.046373077 ,5.906670e-03
 ,0.14524737 ,0.047760788 ,6.249983e-03
 ,0.14884258 ,0.049209682 ,6.612920e-03
 ,0.15242802 ,0.050690819 ,6.986926e-03
 ,0.15599509 ,0.052164992 ,7.363452e-03
 ,0.15951646 ,0.053593885 ,7.735951e-03
 ,0.16295105 ,0.054949857 ,8.100655e-03
 ,0.16625418 ,0.056221112 ,8.456323e-03
 ,0.16938955 ,0.057411310 ,8.803438e-03
 ,0.17233932 ,0.058534904 ,9.143341e-03
 ,0.17510887 ,0.059610687 ,9.477586e-03
 ,0.17772492 ,0.060655906 ,9.807601e-03
 ,0.18022813 ,0.061682484 ,1.013458e-02
 ,0.18266319 ,0.062695748 ,1.045946e-02
 ,0.18506945 ,0.063695282 ,1.078300e-02
 ,0.18747477 ,0.064677064 ,1.110569e-02
 ,0.18989354 ,0.065635990 ,1.142782e-02
 ,0.19232847 ,0.066568052 ,1.174942e-02
 ,0.19477450 ,0.067471733 ,1.207037e-02
 ,0.19722304 ,0.068348479 ,1.239045e-02
 ,0.19966526 ,0.069202350 ,1.270946e-02
 ,0.20209341 ,0.070039138 ,1.302729e-02
 ,0.20450072 ,0.070865253 ,1.334392e-02
 ,0.20688005 ,0.071686698 ,1.365930e-02
 ,0.20922248 ,0.072508287 ,1.397321e-02
 ,0.21151643 ,0.073333231 ,1.428515e-02
 ,0.21374770 ,0.074163062 ,1.459419e-02
 ,0.21590048 ,0.074997845 ,1.489902e-02
 ,0.21795919 ,0.075836565 ,1.519811e-02
 ,0.21991052 ,0.076677610 ,1.548999e-02
 ,0.22174535 ,0.077519237 ,1.577353e-02
 ,0.22346014 ,0.078359986 ,1.604830e-02
 ,0.22505752 ,0.079198980 ,1.631478e-02
 ,0.22654605 ,0.080036104 ,1.657442e-02
 ,0.22793930 ,0.080872053 ,1.682959e-02
 ,0.22925445 ,0.081708279 ,1.708334e-02
 ,0.23051059 ,0.082546845 ,1.733898e-02
 ,0.23172712 ,0.083390212 ,1.759979e-02
 ,0.23292230 ,0.084240996 ,1.786853e-02
 ,0.23411229 ,0.085101710 ,1.814722e-02
 ,0.23531044 ,0.085974525 ,1.843693e-02
 ,0.23652714 ,0.086861072 ,1.873778e-02
 ,0.23776992 ,0.087762307 ,1.904896e-02
 ,0.23904386 ,0.088678461 ,1.936895e-02
 ,0.24035207 ,0.089609084 ,1.969576e-02
 ,0.24169630 ,0.090553177 ,2.002710e-02
 ,0.24307742 ,0.091509412 ,2.036069e-02
 ,0.24449595 ,0.092476402 ,2.069439e-02
 ,0.24595237 ,0.093452991 ,2.102637e-02
 ,0.24744741 ,0.094438523 ,2.135516e-02
 ,0.24898216 ,0.095433051 ,2.167968e-02
 ,0.25055813 ,0.096437453 ,2.199923e-02
 ,0.25217708 ,0.097453412 ,2.231345e-02
 ,0.25384088 ,0.098483278 ,2.262224e-02
 ,0.25555107 ,0.099529803 ,2.292570e-02
 ,0.25730850 ,0.100595780 ,2.322410e-02
 ,0.25911286 ,0.101683637 ,2.351787e-02
 ,0.26096228 ,0.102795033 ,2.380753e-02
 ,0.26285307 ,0.103930520 ,2.409377e-02
 ,0.26477965 ,0.105089312 ,2.437741e-02
 ,0.26673470 ,0.106269188 ,2.465942e-02
 ,0.26870964 ,0.107466555 ,2.494090e-02
 ,0.27069514 ,0.108676633 ,2.522304e-02
 ,0.27268195 ,0.109893752 ,2.550708e-02
 ,0.27466156 ,0.111111714 ,2.579422e-02
 ,0.27662692 ,0.112324166 ,2.608553e-02
 ,0.27857284 ,0.113524958 ,2.638190e-02
 ,0.28049622 ,0.114708449 ,2.668392e-02
 ,0.28239603 ,0.115869742 ,2.699191e-02
 ,0.28427298 ,0.117004848 ,2.730585e-02
 ,0.28612913 ,0.118110784 ,2.762538e-02
 ,0.28796727 ,0.119185602 ,2.794990e-02
 ,0.28979037 ,0.120228378 ,2.827855e-02
 ,0.29160108 ,0.121239153 ,2.861033e-02
 ,0.29340136 ,0.122218858 ,2.894414e-02
 ,0.29519222 ,0.123169211 ,2.927891e-02
 ,0.29697369 ,0.124092605 ,2.961360e-02
 ,0.31397659 ,0.132621459 ,3.287664e-02
 ,0.32905317 ,0.141687415 ,3.635344e-02
 ,0.34142264 ,0.150239874 ,4.011146e-02
 ,0.35320262 ,0.157769982 ,4.336456e-02
 ,0.36543162 ,0.165856637 ,4.597872e-02
 ,0.37719191 ,0.174946656 ,4.927422e-02
 ,0.38887463 ,0.183829329 ,5.322584e-02
 ,0.40152779 ,0.191918219 ,5.706844e-02
 ,0.41555978 ,0.200422337 ,6.072351e-02
 ,0.42908873 ,0.209907793 ,6.413420e-02
 ,0.44052818 ,0.219253257 ,6.745447e-02
 ,0.45132279 ,0.227775172 ,7.077917e-02
 ,0.46314082 ,0.234996714 ,7.417266e-02
 ,0.47488247 ,0.241586411 ,7.752946e-02
 ,0.48469259 ,0.248330403 ,8.082463e-02
 ,0.49326509 ,0.255625615 ,8.409529e-02
 ,0.50168814 ,0.263258846 ,8.746989e-02
 ,0.51028363 ,0.271055706 ,9.086717e-02
 ,0.51912560 ,0.278856970 ,9.420999e-02
 ,0.52769780 ,0.286320658 ,9.786763e-02
 ,0.53603821 ,0.293266278 ,1.020140e-01
 ,0.54478766 ,0.299467332 ,1.062560e-01
 ,0.55366271 ,0.305690735 ,1.102368e-01
 ,0.56272309 ,0.312768360 ,1.140490e-01
 ,0.57233129 ,0.319963119 ,1.178623e-01
 ,0.58200310 ,0.326558078 ,1.219016e-01
 ,0.59099024 ,0.332901385 ,1.263491e-01
 ,0.59937345 ,0.339303862 ,1.310123e-01
 ,0.60758641 ,0.345851059 ,1.357467e-01
 ,0.61598033 ,0.352190085 ,1.405864e-01
 ,0.62528643 ,0.357913761 ,1.455560e-01
 ,0.63543095 ,0.363356622 ,1.504331e-01
 ,0.64471050 ,0.369132056 ,1.551387e-01
 ,0.65228925 ,0.375308548 ,1.599456e-01
 ,0.65934374 ,0.381636254 ,1.650710e-01
 ,0.66706291 ,0.388370558 ,1.704709e-01
 ,0.67515257 ,0.395534711 ,1.758367e-01
 ,0.68284073 ,0.402421344 ,1.812059e-01
 ,0.69038336 ,0.408623429 ,1.868699e-01
 ,0.69845919 ,0.414460275 ,1.926662e-01
 ,0.70715142 ,0.420523246 ,1.981833e-01
 ,0.71574865 ,0.427206110 ,2.033596e-01
 ,0.72320765 ,0.434498031 ,2.084776e-01
 ,0.72954754 ,0.442149909 ,2.137674e-01
 ,0.73571186 ,0.449959300 ,2.191326e-01
 ,0.74226506 ,0.457421738 ,2.245519e-01
 ,0.74905392 ,0.464237364 ,2.303037e-01
 ,0.75612564 ,0.470902809 ,2.367213e-01
 ,0.76388865 ,0.477994733 ,2.435881e-01
 ,0.77159609 ,0.485500173 ,2.500961e-01
 ,0.77797709 ,0.493064413 ,2.561486e-01
 ,0.78328352 ,0.500852021 ,2.622563e-01
 ,0.78864834 ,0.509132560 ,2.684670e-01
 ,0.79436056 ,0.517227408 ,2.749325e-01
 ,0.80015524 ,0.524694907 ,2.817636e-01
 ,0.80624435 ,0.532053355 ,2.884733e-01
 ,0.81293829 ,0.540118580 ,2.947507e-01
 ,0.82015085 ,0.548921459 ,3.010033e-01
 ,0.82763569 ,0.557637151 ,3.078365e-01
 ,0.83505891 ,0.566310996 ,3.155196e-01
 ,0.84198778 ,0.575253522 ,3.239425e-01
 ,0.84848454 ,0.583976152 ,3.326209e-01
 ,0.85494638 ,0.592377395 ,3.411542e-01
 ,0.86153456 ,0.600715905 ,3.495319e-01
 ,0.86786548 ,0.608554390 ,3.581296e-01
 ,0.87340081 ,0.615887245 ,3.669205e-01
 ,0.87845635 ,0.623939370 ,3.759445e-01
 ,0.88379941 ,0.632832612 ,3.856236e-01
 ,0.88951674 ,0.641696896 ,3.957400e-01
 ,0.89459349 ,0.651106256 ,4.056838e-01
 ,0.89882091 ,0.662461064 ,4.146894e-01
 ,0.90377294 ,0.675025139 ,4.226494e-01
 ,0.90970927 ,0.687246375 ,4.314995e-01
 ,0.91531400 ,0.699233434 ,4.438556e-01
 ,0.92062598 ,0.710363456 ,4.575252e-01
 ,0.92625531 ,0.720885287 ,4.697992e-01
 ,0.93147783 ,0.730627162 ,4.799491e-01
 ,0.93661779 ,0.740878409 ,4.909327e-01
 ,0.94256155 ,0.753110750 ,5.054689e-01
 ,0.94879244 ,0.766049888 ,5.218112e-01
 ,0.95441461 ,0.779461673 ,5.382110e-01
 ,0.95922582 ,0.795415271 ,5.580688e-01
 ,0.96461303 ,0.812122645 ,5.830378e-01
 ,0.96954495 ,0.827981003 ,6.131022e-01
 ,0.97446302 ,0.845818141 ,6.404610e-01
 ,0.97912580 ,0.862982643 ,6.619244e-01
 ,0.98539520 ,0.887825781 ,6.905761e-01
 ,0.98986795 ,0.916261737 ,7.236242e-01
 ,0.99487474 ,0.960335712 ,7.689094e-01),
 byrow=TRUE,ncol=3)

P5=matrix(c( 0.05652412 ,0.01305194 ,0.002390178 ,3.378326e-12
 ,0.06574073 ,0.01487696 ,0.003342965 ,7.244389e-11
 ,0.07151676 ,0.01720008 ,0.004648512 ,7.773423e-10
 ,0.07511946 ,0.01981078 ,0.006063858 ,5.567153e-09
 ,0.07828081 ,0.02244535 ,0.007485693 ,2.994702e-08
 ,0.08180130 ,0.02494366 ,0.008893920 ,1.291037e-07
 ,0.08583995 ,0.02730528 ,0.010278443 ,4.647955e-07
 ,0.09025786 ,0.02960906 ,0.011622679 ,1.437865e-06
 ,0.09481793 ,0.03192766 ,0.012912280 ,3.903397e-06
 ,0.09930536 ,0.03429573 ,0.014138866 ,9.451086e-06
 ,0.10358876 ,0.03671695 ,0.015295959 ,2.067608e-05
 ,0.10762929 ,0.03918039 ,0.016375245 ,4.130899e-05
 ,0.11145588 ,0.04167054 ,0.017367679 ,7.605707e-05
 ,0.11512861 ,0.04416935 ,0.018267991 ,1.300654e-04
 ,0.11870678 ,0.04665486 ,0.019078858 ,2.080362e-04
 ,0.12222956 ,0.04910107 ,0.019812158 ,3.132019e-04
 ,0.12571005 ,0.05148055 ,0.020487010 ,4.464557e-04
 ,0.12913936 ,0.05376876 ,0.021125903 ,6.059430e-04
 ,0.13249584 ,0.05594791 ,0.021750598 ,7.872902e-04
 ,0.13575514 ,0.05800867 ,0.022379035 ,9.844397e-04
 ,0.13889807 ,0.05994945 ,0.023023715 ,1.190860e-03
 ,0.14191516 ,0.06177395 ,0.023691420 ,1.400786e-03
 ,0.14480785 ,0.06348830 ,0.024383868 ,1.610156e-03
 ,0.14758751 ,0.06509889 ,0.025098853 ,1.817036e-03
 ,0.15027304 ,0.06661130 ,0.025831529 ,2.021495e-03
 ,0.15288821 ,0.06803034 ,0.026575628 ,2.225055e-03
 ,0.15545896 ,0.06936082 ,0.027324476 ,2.429942e-03
 ,0.15801110 ,0.07060840 ,0.028071770 ,2.638347e-03
 ,0.16056820 ,0.07178026 ,0.028812074 ,2.851885e-03
 ,0.16315001 ,0.07288527 ,0.029541074 ,3.071323e-03
 ,0.16577118 ,0.07393380 ,0.030255603 ,3.296558e-03
 ,0.16844049 ,0.07493712 ,0.030953535 ,3.526795e-03
 ,0.17116069 ,0.07590682 ,0.031633586 ,3.760793e-03
 ,0.17392879 ,0.07685417 ,0.032295108 ,3.997135e-03
 ,0.17673683 ,0.07778972 ,0.032937919 ,4.234427e-03
 ,0.17957306 ,0.07872305 ,0.033562189 ,4.471420e-03
 ,0.18242312 ,0.07966272 ,0.034168393 ,4.707065e-03
 ,0.18527142 ,0.08061630 ,0.034757295 ,4.940502e-03
 ,0.18810221 ,0.08159041 ,0.035329958 ,5.171027e-03
 ,0.19090055 ,0.08259081 ,0.035887740 ,5.398048e-03
 ,0.19365295 ,0.08362233 ,0.036432265 ,5.621047e-03
 ,0.19634778 ,0.08468879 ,0.036965369 ,5.839553e-03
 ,0.19897548 ,0.08579287 ,0.037489019 ,6.053142e-03
 ,0.20152856 ,0.08693588 ,0.038005212 ,6.261435e-03
 ,0.20400154 ,0.08811768 ,0.038515895 ,6.464128e-03
 ,0.20639083 ,0.08933656 ,0.039022878 ,6.661029e-03
 ,0.20869457 ,0.09058924 ,0.039527788 ,6.852104e-03
 ,0.21091258 ,0.09187106 ,0.040032038 ,7.037528e-03
 ,0.21304622 ,0.09317615 ,0.040536826 ,7.217739e-03
 ,0.21509838 ,0.09449787 ,0.041043146 ,7.393470e-03
 ,0.21707346 ,0.09582915 ,0.041551812 ,7.565765e-03
 ,0.21897727 ,0.09716299 ,0.042063485 ,7.735961e-03
 ,0.22081697 ,0.09849292 ,0.042578697 ,7.905643e-03
 ,0.22260093 ,0.09981337 ,0.043097865 ,8.076563e-03
 ,0.22433851 ,0.10112002 ,0.043621312 ,8.250539e-03
 ,0.22603978 ,0.10240999 ,0.044149265 ,8.429335e-03
 ,0.22771524 ,0.10368199 ,0.044681863 ,8.614546e-03
 ,0.22937544 ,0.10493623 ,0.045219159 ,8.807490e-03
 ,0.23103064 ,0.10617433 ,0.045761127 ,9.009122e-03
 ,0.23269046 ,0.10739901 ,0.046307667 ,9.219978e-03
 ,0.23436362 ,0.10861381 ,0.046858620 ,9.440156e-03
 ,0.23605766 ,0.10982269 ,0.047413783 ,9.669322e-03
 ,0.23777881 ,0.11102960 ,0.047972925 ,9.906756e-03
 ,0.23953181 ,0.11223815 ,0.048535807 ,1.015141e-02
 ,0.24131994 ,0.11345120 ,0.049102198 ,1.040201e-02
 ,0.24314491 ,0.11467064 ,0.049671889 ,1.065710e-02
 ,0.24500695 ,0.11589718 ,0.050244706 ,1.091518e-02
 ,0.24690483 ,0.11713029 ,0.050820511 ,1.117476e-02
 ,0.24883595 ,0.11836818 ,0.051399206 ,1.143446e-02
 ,0.25079643 ,0.11960797 ,0.051980725 ,1.169304e-02
 ,0.25278126 ,0.12084583 ,0.052565028 ,1.194945e-02
 ,0.25478439 ,0.12207730 ,0.053152083 ,1.220290e-02
 ,0.25679899 ,0.12329752 ,0.053741852 ,1.245283e-02
 ,0.25881759 ,0.12450159 ,0.054334269 ,1.269891e-02
 ,0.26083237 ,0.12568485 ,0.054929219 ,1.294106e-02
 ,0.26283537 ,0.12684313 ,0.055526514 ,1.317939e-02
 ,0.26481879 ,0.12797294 ,0.056125873 ,1.341420e-02
 ,0.26677527 ,0.12907169 ,0.056726901 ,1.364589e-02
 ,0.26869809 ,0.13013771 ,0.057329079 ,1.387498e-02
 ,0.27058146 ,0.13117030 ,0.057931750 ,1.410201e-02
 ,0.27242063 ,0.13216971 ,0.058534126 ,1.432755e-02
 ,0.27421212 ,0.13313704 ,0.059135301 ,1.455210e-02
 ,0.27595371 ,0.13407414 ,0.059734269 ,1.477611e-02
 ,0.27764451 ,0.13498346 ,0.060329960 ,1.499992e-02
 ,0.27928490 ,0.13586793 ,0.060921287 ,1.522377e-02
 ,0.28087643 ,0.13673079 ,0.061507187 ,1.544778e-02
 ,0.28242170 ,0.13757549 ,0.062086680 ,1.567195e-02
 ,0.28392419 ,0.13840552 ,0.062658920 ,1.589617e-02
 ,0.28538803 ,0.13922434 ,0.063223239 ,1.612023e-02
 ,0.28681780 ,0.14003533 ,0.063779188 ,1.634387e-02
 ,0.28821831 ,0.14084167 ,0.064326561 ,1.656674e-02
 ,0.28959443 ,0.14164635 ,0.064865403 ,1.678851e-02
 ,0.29095083 ,0.14245215 ,0.065396011 ,1.700879e-02
 ,0.29229189 ,0.14326158 ,0.065918911 ,1.722723e-02
 ,0.29362151 ,0.14407696 ,0.066434823 ,1.744351e-02
 ,0.29494305 ,0.14490033 ,0.066944618 ,1.765732e-02
 ,0.29625927 ,0.14573355 ,0.067449265 ,1.786841e-02
 ,0.29757232 ,0.14657825 ,0.067949777 ,1.807657e-02
 ,0.29888374 ,0.14743586 ,0.068447155 ,1.828162e-02
 ,0.30019451 ,0.14830758 ,0.068942336 ,1.848344e-02
 ,0.31328362 ,0.15791008 ,0.073880466 ,2.031038e-02
 ,0.32658914 ,0.16858633 ,0.078745227 ,2.182686e-02
 ,0.34081183 ,0.17902263 ,0.083679715 ,2.336938e-02
 ,0.35529268 ,0.18923993 ,0.089192608 ,2.528373e-02
 ,0.36855379 ,0.19911099 ,0.095098319 ,2.734008e-02
 ,0.37994080 ,0.20799604 ,0.100761289 ,2.928966e-02
 ,0.39062100 ,0.21671444 ,0.105774473 ,3.127971e-02
 ,0.40257143 ,0.22609178 ,0.110115264 ,3.337375e-02
 ,0.41621275 ,0.23569341 ,0.113958744 ,3.539360e-02
 ,0.43046818 ,0.24490377 ,0.117643306 ,3.736688e-02
 ,0.44486669 ,0.25330316 ,0.121679236 ,3.972222e-02
 ,0.45923530 ,0.26068382 ,0.126584071 ,4.280986e-02
 ,0.47296930 ,0.26745795 ,0.132471529 ,4.640810e-02
 ,0.48590094 ,0.27415019 ,0.138753880 ,5.001684e-02
 ,0.49834493 ,0.28093903 ,0.144554010 ,5.345159e-02
 ,0.51026131 ,0.28775596 ,0.149452848 ,5.679297e-02
 ,0.52119211 ,0.29432490 ,0.153686512 ,5.996537e-02
 ,0.53082589 ,0.30059787 ,0.157685666 ,6.275756e-02
 ,0.53943126 ,0.30692404 ,0.161762480 ,6.518147e-02
 ,0.54774724 ,0.31353121 ,0.166120746 ,6.753031e-02
 ,0.55642024 ,0.32037096 ,0.170885289 ,7.007491e-02
 ,0.56550181 ,0.32736781 ,0.176132244 ,7.279915e-02
 ,0.57464883 ,0.33444881 ,0.181843360 ,7.550823e-02
 ,0.58368440 ,0.34149272 ,0.187714538 ,7.816054e-02
 ,0.59259791 ,0.34836154 ,0.193239542 ,8.092069e-02
 ,0.60135765 ,0.35491577 ,0.198212454 ,8.389308e-02
 ,0.60997532 ,0.36111908 ,0.202943209 ,8.697786e-02
 ,0.61843676 ,0.36715577 ,0.207920518 ,9.000649e-02
 ,0.62674513 ,0.37332556 ,0.213468203 ,9.287236e-02
 ,0.63517287 ,0.37979485 ,0.219617217 ,9.553586e-02
 ,0.64411712 ,0.38648760 ,0.226087437 ,9.803810e-02
 ,0.65355691 ,0.39327005 ,0.232464891 ,1.005292e-01
 ,0.66286128 ,0.40022900 ,0.238568656 ,1.032200e-01
 ,0.67142296 ,0.40767245 ,0.244631428 ,1.062392e-01
 ,0.67945227 ,0.41574779 ,0.250960847 ,1.095391e-01
 ,0.68784816 ,0.42420043 ,0.257504402 ,1.129935e-01
 ,0.69717228 ,0.43269300 ,0.263984690 ,1.165367e-01
 ,0.70693776 ,0.44108963 ,0.270350361 ,1.201702e-01
 ,0.71609195 ,0.44926200 ,0.276782149 ,1.239319e-01
 ,0.72398625 ,0.45697627 ,0.283281861 ,1.278724e-01
 ,0.73066318 ,0.46420271 ,0.289549848 ,1.319571e-01
 ,0.73650853 ,0.47133730 ,0.295322586 ,1.360494e-01
 ,0.74191863 ,0.47888254 ,0.300607528 ,1.400640e-01
 ,0.74721366 ,0.48699345 ,0.305518411 ,1.440467e-01
 ,0.75262380 ,0.49549336 ,0.310111010 ,1.480611e-01
 ,0.75824105 ,0.50420009 ,0.314535705 ,1.520995e-01
 ,0.76409670 ,0.51305053 ,0.319195032 ,1.561658e-01
 ,0.77030423 ,0.52194190 ,0.324545735 ,1.603832e-01
 ,0.77696430 ,0.53068348 ,0.330761451 ,1.648740e-01
 ,0.78392814 ,0.53921911 ,0.337687815 ,1.695305e-01
 ,0.79087142 ,0.54764584 ,0.345040336 ,1.740919e-01
 ,0.79764992 ,0.55599003 ,0.352535977 ,1.784667e-01
 ,0.80440374 ,0.56418827 ,0.359947893 ,1.828362e-01
 ,0.81125920 ,0.57217204 ,0.367126284 ,1.874504e-01
 ,0.81807396 ,0.57997214 ,0.374007327 ,1.923635e-01
 ,0.82455298 ,0.58778776 ,0.380706775 ,1.973552e-01
 ,0.83054464 ,0.59568049 ,0.387507714 ,2.021553e-01
 ,0.83619452 ,0.60335949 ,0.394648281 ,2.067057e-01
 ,0.84182788 ,0.61067252 ,0.402210773 ,2.111209e-01
 ,0.84769382 ,0.61810882 ,0.410130422 ,2.155355e-01
 ,0.85379684 ,0.62646991 ,0.418239217 ,2.201354e-01
 ,0.85982991 ,0.63596326 ,0.426436850 ,2.251539e-01
 ,0.86534864 ,0.64569215 ,0.434743378 ,2.307076e-01
 ,0.87020343 ,0.65434970 ,0.443196556 ,2.368692e-01
 ,0.87472363 ,0.66158858 ,0.451892254 ,2.437247e-01
 ,0.87942450 ,0.66828595 ,0.460851245 ,2.509675e-01
 ,0.88437025 ,0.67525025 ,0.469770141 ,2.580093e-01
 ,0.88904758 ,0.68246978 ,0.478325822 ,2.648763e-01
 ,0.89322701 ,0.68982898 ,0.486536728 ,2.723872e-01
 ,0.89745058 ,0.69753791 ,0.494776780 ,2.812096e-01
 ,0.90239053 ,0.70572617 ,0.503530501 ,2.908907e-01
 ,0.90800180 ,0.71430596 ,0.513032547 ,3.003753e-01
 ,0.91346775 ,0.72311282 ,0.523668882 ,3.096709e-01
 ,0.91813323 ,0.73249598 ,0.535903771 ,3.196147e-01
 ,0.92221702 ,0.74332982 ,0.549315193 ,3.304979e-01
 ,0.92622555 ,0.75519029 ,0.562491098 ,3.426984e-01
 ,0.93049876 ,0.76676088 ,0.574311153 ,3.573481e-01
 ,0.93559016 ,0.77792982 ,0.584778988 ,3.742742e-01
 ,0.94158878 ,0.78972489 ,0.594741057 ,3.918491e-01
 ,0.94715113 ,0.80271250 ,0.605197756 ,4.107232e-01
 ,0.95168452 ,0.81675031 ,0.616562519 ,4.307354e-01
 ,0.95643627 ,0.83145104 ,0.629586884 ,4.491496e-01
 ,0.96150476 ,0.84592600 ,0.643095195 ,4.683523e-01
 ,0.96619142 ,0.86139307 ,0.658098224 ,4.899430e-01
 ,0.97232678 ,0.87780728 ,0.682404452 ,5.109344e-01
 ,0.97800785 ,0.89418487 ,0.715182394 ,5.334721e-01
 ,0.98168187 ,0.90711236 ,0.750757190 ,5.668841e-01
 ,0.98661473 ,0.92175945 ,0.795670468 ,6.221022e-01
 ,0.99161624 ,0.94482117 ,0.844332168 ,6.873383e-01),
 byrow=TRUE,ncol=4)

 P6=matrix(c( 0.02845841 ,0.009644671 ,0.0009664242 ,0.0007502456 ,2.594514e-11
 ,0.03862415 ,0.017910204 ,0.0031407471 ,0.0015589418 ,5.022928e-10
 ,0.05220254 ,0.025693640 ,0.0057938158 ,0.0021160814 ,4.868380e-09
 ,0.06621145 ,0.032227944 ,0.0082558224 ,0.0025765248 ,3.151134e-08
 ,0.07867811 ,0.037545052 ,0.0103200867 ,0.0030739790 ,1.532944e-07
 ,0.08886734 ,0.041725975 ,0.0120483568 ,0.0036182088 ,5.980879e-07
 ,0.09687788 ,0.044882223 ,0.0135631342 ,0.0041650257 ,1.950299e-06
 ,0.10316830 ,0.047217776 ,0.0149788321 ,0.0046727336 ,5.469913e-06
 ,0.10824058 ,0.048995043 ,0.0163941682 ,0.0051213088 ,1.347713e-05
 ,0.11250832 ,0.050462880 ,0.0178877558 ,0.0055125816 ,2.965286e-05
 ,0.11627805 ,0.051812087 ,0.0195084853 ,0.0058643506 ,5.903392e-05
 ,0.11976844 ,0.053169834 ,0.0212684847 ,0.0062028177 ,1.075073e-04
 ,0.12313219 ,0.054614572 ,0.0231443453 ,0.0065546842 ,1.807643e-04
 ,0.12647380 ,0.056192980 ,0.0250870063 ,0.0069403944 ,2.829186e-04
 ,0.12986454 ,0.057930892 ,0.0270368653 ,0.0073701618 ,4.152024e-04
 ,0.13335383 ,0.059838004 ,0.0289391057 ,0.0078435357 ,5.752174e-04
 ,0.13697580 ,0.061909087 ,0.0307548517 ,0.0083519589 ,7.570695e-04
 ,0.14075114 ,0.064124538 ,0.0324658987 ,0.0088828750 ,9.524158e-04
 ,0.14468601 ,0.066452103 ,0.0340732641 ,0.0094238009 ,1.152134e-03
 ,0.14877041 ,0.068850570 ,0.0355916148 ,0.0099652366 ,1.348124e-03
 ,0.15297783 ,0.071275211 ,0.0370422456 ,0.0105019844 ,1.534731e-03
 ,0.15726699 ,0.073684014 ,0.0384468508 ,0.0110330581 ,1.709452e-03
 ,0.16158577 ,0.076043401 ,0.0398233432 ,0.0115607019 ,1.872830e-03
 ,0.16587656 ,0.078332254 ,0.0411839642 ,0.0120890846 ,2.027675e-03
 ,0.17008218 ,0.080543571 ,0.0425352335 ,0.0126230926 ,2.177920e-03
 ,0.17415147 ,0.082683715 ,0.0438790091 ,0.0131674283 ,2.327430e-03
 ,0.17804364 ,0.084769760 ,0.0452139704 ,0.0137260428 ,2.479054e-03
 ,0.18173091 ,0.086825700 ,0.0465370496 ,0.0143018329 ,2.634074e-03
 ,0.18519926 ,0.088878323 ,0.0478445734 ,0.0148965065 ,2.792108e-03
 ,0.18844742 ,0.090953392 ,0.0491330541 ,0.0155105477 ,2.951390e-03
 ,0.19148459 ,0.093072518 ,0.0503996631 ,0.0161432460 ,3.109303e-03
 ,0.19432745 ,0.095250927 ,0.0516424566 ,0.0167927877 ,3.263029e-03
 ,0.19699700 ,0.097496165 ,0.0528604171 ,0.0174564124 ,3.410139e-03
 ,0.19951573 ,0.099807706 ,0.0540533669 ,0.0181306383 ,3.549062e-03
 ,0.20190550 ,0.102177381 ,0.0552217983 ,0.0188115440 ,3.679327e-03
 ,0.20418614 ,0.104590523 ,0.0563666642 ,0.0194950810 ,3.801605e-03
 ,0.20637484 ,0.107027661 ,0.0574891666 ,0.0201773871 ,3.917550e-03
 ,0.20848607 ,0.109466576 ,0.0585905770 ,0.0208550641 ,4.029514e-03
 ,0.21053195 ,0.111884485 ,0.0596721117 ,0.0215253931 ,4.140185e-03
 ,0.21252275 ,0.114260111 ,0.0607348702 ,0.0221864686 ,4.252232e-03
 ,0.21446746 ,0.116575429 ,0.0617798326 ,0.0228372461 ,4.367992e-03
 ,0.21637417 ,0.118816908 ,0.0628079009 ,0.0234775060 ,4.489244e-03
 ,0.21825036 ,0.120976170 ,0.0638199623 ,0.0241077482 ,4.617087e-03
 ,0.22010296 ,0.123050046 ,0.0648169547 ,0.0247290363 ,4.751923e-03
 ,0.22193833 ,0.125040101 ,0.0657999179 ,0.0253428107 ,4.893526e-03
 ,0.22376206 ,0.126951747 ,0.0667700240 ,0.0259506919 ,5.041168e-03
 ,0.22557886 ,0.128793104 ,0.0677285875 ,0.0265542902 ,5.193786e-03
 ,0.22739239 ,0.130573768 ,0.0686770614 ,0.0271550377 ,5.350146e-03
 ,0.22920525 ,0.132303638 ,0.0696170311 ,0.0277540524 ,5.508989e-03
 ,0.23101900 ,0.133991919 ,0.0705502132 ,0.0283520430 ,5.669140e-03
 ,0.23283434 ,0.135646375 ,0.0714784671 ,0.0289492606 ,5.829573e-03
 ,0.23465133 ,0.137272882 ,0.0724038164 ,0.0295454968 ,5.989428e-03
 ,0.23646965 ,0.138875254 ,0.0733284760 ,0.0301401298 ,6.148001e-03
 ,0.23828893 ,0.140455337 ,0.0742548720 ,0.0307322114 ,6.304704e-03
 ,0.24010890 ,0.142013293 ,0.0751856416 ,0.0313205878 ,6.459029e-03
 ,0.24192967 ,0.143548021 ,0.0761236012 ,0.0319040403 ,6.610507e-03
 ,0.24375172 ,0.145057640 ,0.0770716723 ,0.0324814318 ,6.758695e-03
 ,0.24557590 ,0.146539976 ,0.0780327640 ,0.0330518431 ,6.903176e-03
 ,0.24740333 ,0.147992990 ,0.0790096154 ,0.0336146829 ,7.043585e-03
 ,0.24923515 ,0.149415131 ,0.0800046099 ,0.0341697603 ,7.179648e-03
 ,0.25107232 ,0.150805565 ,0.0810195786 ,0.0347173115 ,7.311226e-03
 ,0.25291533 ,0.152164307 ,0.0820556131 ,0.0352579785 ,7.438362e-03
 ,0.25476396 ,0.153492230 ,0.0831129092 ,0.0357927443 ,7.561312e-03
 ,0.25661711 ,0.154791002 ,0.0841906620 ,0.0363228330 ,7.680559e-03
 ,0.25847267 ,0.156062941 ,0.0852870246 ,0.0368495882 ,7.796806e-03
 ,0.26032752 ,0.157310848 ,0.0863991407 ,0.0373743441 ,7.910946e-03
 ,0.26217762 ,0.158537805 ,0.0875232506 ,0.0378983036 ,8.024013e-03
 ,0.26401815 ,0.159746996 ,0.0886548636 ,0.0384224367 ,8.137120e-03
 ,0.26584379 ,0.160941546 ,0.0897889847 ,0.0389474062 ,8.251387e-03
 ,0.26764897 ,0.162124387 ,0.0909203766 ,0.0394735277 ,8.367871e-03
 ,0.26942825 ,0.163298178 ,0.0920438371 ,0.0400007620 ,8.487505e-03
 ,0.27117659 ,0.164465254 ,0.0931544716 ,0.0405287383 ,8.611045e-03
 ,0.27288969 ,0.165627620 ,0.0942479401 ,0.0410567998 ,8.739032e-03
 ,0.27456420 ,0.166786971 ,0.0953206650 ,0.0415840663 ,8.871781e-03
 ,0.27619794 ,0.167944735 ,0.0963699857 ,0.0421095038 ,9.009371e-03
 ,0.27779004 ,0.169102127 ,0.0973942565 ,0.0426319942 ,9.151667e-03
 ,0.27934099 ,0.170260209 ,0.0983928834 ,0.0431503996 ,9.298339e-03
 ,0.28085262 ,0.171419936 ,0.0993663043 ,0.0436636169 ,9.448902e-03
 ,0.28232802 ,0.172582199 ,0.1003159175 ,0.0441706205 ,9.602751e-03
 ,0.28377144 ,0.173747838 ,0.1012439683 ,0.0446704936 ,9.759205e-03
 ,0.28518812 ,0.174917647 ,0.1021534027 ,0.0451624485 ,9.917545e-03
 ,0.28658401 ,0.176092356 ,0.1030476997 ,0.0456458386 ,1.007705e-02
 ,0.28796561 ,0.177272586 ,0.1039306913 ,0.0461201648 ,1.023703e-02
 ,0.28933969 ,0.178458806 ,0.1048063812 ,0.0465850770 ,1.039685e-02
 ,0.29071301 ,0.179651274 ,0.1056787690 ,0.0470403750 ,1.055595e-02
 ,0.29209211 ,0.180849979 ,0.1065516870 ,0.0474860075 ,1.071386e-02
 ,0.29348304 ,0.182054592 ,0.1074286553 ,0.0479220717 ,1.087020e-02
 ,0.29489121 ,0.183264423 ,0.1083127575 ,0.0483488124 ,1.102471e-02
 ,0.29632114 ,0.184478407 ,0.1092065408 ,0.0487666200 ,1.117722e-02
 ,0.29777639 ,0.185695102 ,0.1101119397 ,0.0491760269 ,1.132768e-02
 ,0.29925947 ,0.186912720 ,0.1110302259 ,0.0495777007 ,1.147610e-02
 ,0.30077178 ,0.188129170 ,0.1119619808 ,0.0499724340 ,1.162261e-02
 ,0.30231363 ,0.189342139 ,0.1129070920 ,0.0503611291 ,1.176737e-02
 ,0.30388432 ,0.190549177 ,0.1138647708 ,0.0507447791 ,1.191065e-02
 ,0.30548226 ,0.191747804 ,0.1148335893 ,0.0511244432 ,1.205270e-02
 ,0.30710505 ,0.192935620 ,0.1158115351 ,0.0515012188 ,1.219386e-02
 ,0.30874972 ,0.194110415 ,0.1167960815 ,0.0518762105 ,1.233443e-02
 ,0.31041287 ,0.195270270 ,0.1177842705 ,0.0522504975 ,1.247476e-02
 ,0.31209087 ,0.196413653 ,0.1187728072 ,0.0526251016 ,1.261515e-02
 ,0.31378008 ,0.197539481 ,0.1197581619 ,0.0530009554 ,1.275592e-02
 ,0.33073053 ,0.207947746 ,0.1287029996 ,0.0569307733 ,1.424281e-02
 ,0.34766587 ,0.218180516 ,0.1351175921 ,0.0609818273 ,1.602065e-02
 ,0.36497276 ,0.229372730 ,0.1404381947 ,0.0646360301 ,1.803053e-02
 ,0.38214131 ,0.240774517 ,0.1460402987 ,0.0678929574 ,1.984669e-02
 ,0.39892707 ,0.251527352 ,0.1517699541 ,0.0710908941 ,2.140172e-02
 ,0.41477035 ,0.261312042 ,0.1577197965 ,0.0744655082 ,2.290802e-02
 ,0.42985193 ,0.270169908 ,0.1649510495 ,0.0779547358 ,2.442240e-02
 ,0.44459125 ,0.278214867 ,0.1733065465 ,0.0816344914 ,2.596410e-02
 ,0.45830629 ,0.285626732 ,0.1811614844 ,0.0857820745 ,2.747529e-02
 ,0.47040113 ,0.293106158 ,0.1877528247 ,0.0903640622 ,2.897806e-02
 ,0.48149178 ,0.301009058 ,0.1935701784 ,0.0951257701 ,3.068337e-02
 ,0.49230962 ,0.308762383 ,0.1992769670 ,0.0998239781 ,3.268333e-02
 ,0.50293366 ,0.315833032 ,0.2049784969 ,0.1044290674 ,3.485126e-02
 ,0.51329043 ,0.322327249 ,0.2104162141 ,0.1090755934 ,3.704227e-02
 ,0.52322279 ,0.328643258 ,0.2155878925 ,0.1137388083 ,3.919411e-02
 ,0.53250776 ,0.334876631 ,0.2207245049 ,0.1182751586 ,4.126276e-02
 ,0.54136758 ,0.341122307 ,0.2258390954 ,0.1227003370 ,4.321786e-02
 ,0.55057692 ,0.347936973 ,0.2307672090 ,0.1272211470 ,4.510598e-02
 ,0.56083373 ,0.355714012 ,0.2355401547 ,0.1320061535 ,4.697462e-02
 ,0.57188476 ,0.364046026 ,0.2403678821 ,0.1369508295 ,4.879706e-02
 ,0.58258251 ,0.372072572 ,0.2454048413 ,0.1418771565 ,5.057703e-02
 ,0.59206406 ,0.379279106 ,0.2507385738 ,0.1468283555 ,5.241157e-02
 ,0.60043277 ,0.385862014 ,0.2565083544 ,0.1518653755 ,5.442717e-02
 ,0.60824065 ,0.392470283 ,0.2627890921 ,0.1568425686 ,5.672558e-02
 ,0.61587328 ,0.399665003 ,0.2693240828 ,0.1615409971 ,5.932340e-02
 ,0.62350556 ,0.407442902 ,0.2757405639 ,0.1659131190 ,6.211202e-02
 ,0.63130400 ,0.415294312 ,0.2820549478 ,0.1701469869 ,6.493040e-02
 ,0.63950683 ,0.422769887 ,0.2885097791 ,0.1744470541 ,6.766549e-02
 ,0.64822576 ,0.429800809 ,0.2949604623 ,0.1788535306 ,7.024569e-02
 ,0.65718853 ,0.436524548 ,0.3009940240 ,0.1832922406 ,7.260984e-02
 ,0.66586026 ,0.443072284 ,0.3066417800 ,0.1877498418 ,7.477539e-02
 ,0.67385499 ,0.449570002 ,0.3124210302 ,0.1923869531 ,7.686614e-02
 ,0.68116860 ,0.456177204 ,0.3186629632 ,0.1973744564 ,7.901179e-02
 ,0.68817020 ,0.463016405 ,0.3252135392 ,0.2026604329 ,8.127321e-02
 ,0.69545456 ,0.470153870 ,0.3317877025 ,0.2080443155 ,8.365480e-02
 ,0.70339322 ,0.477628656 ,0.3382869934 ,0.2133829471 ,8.614100e-02
 ,0.71168363 ,0.485354548 ,0.3446879623 ,0.2185816498 ,8.874283e-02
 ,0.71966391 ,0.493111522 ,0.3508452160 ,0.2235543145 ,9.151748e-02
 ,0.72706877 ,0.500759585 ,0.3565974060 ,0.2283233375 ,9.452512e-02
 ,0.73411394 ,0.508436756 ,0.3620578719 ,0.2330030834 ,9.779298e-02
 ,0.74096543 ,0.516596717 ,0.3676326358 ,0.2376654675 ,1.013167e-01
 ,0.74754062 ,0.525572730 ,0.3736020261 ,0.2423531847 ,1.050525e-01
 ,0.75389357 ,0.534868176 ,0.3798150608 ,0.2471555672 ,1.089039e-01
 ,0.76044477 ,0.543419973 ,0.3859717103 ,0.2520806896 ,1.127538e-01
 ,0.76755285 ,0.550743471 ,0.3920256074 ,0.2568981812 ,1.165591e-01
 ,0.77509023 ,0.557269510 ,0.3981237467 ,0.2613608036 ,1.203644e-01
 ,0.78258820 ,0.563602135 ,0.4043785159 ,0.2656200877 ,1.241878e-01
 ,0.78956970 ,0.569952952 ,0.4108934218 ,0.2702271522 ,1.280316e-01
 ,0.79578357 ,0.576313169 ,0.4177625748 ,0.2756524098 ,1.320053e-01
 ,0.80125054 ,0.582816465 ,0.4248741636 ,0.2818778472 ,1.362441e-01
 ,0.80604156 ,0.589718190 ,0.4319593756 ,0.2884428281 ,1.407053e-01
 ,0.81021918 ,0.597012063 ,0.4388877066 ,0.2948281054 ,1.452044e-01
 ,0.81409451 ,0.604298241 ,0.4456550600 ,0.3008228462 ,1.496867e-01
 ,0.81820723 ,0.611181771 ,0.4521718937 ,0.3065488400 ,1.543527e-01
 ,0.82291940 ,0.617566520 ,0.4583418105 ,0.3121847541 ,1.593959e-01
 ,0.82816550 ,0.623606817 ,0.4642293076 ,0.3178061250 ,1.646657e-01
 ,0.83358635 ,0.629645757 ,0.4699982684 ,0.3234987538 ,1.698169e-01
 ,0.83888900 ,0.636078808 ,0.4758026388 ,0.3294433637 ,1.747666e-01
 ,0.84407903 ,0.643093954 ,0.4817498902 ,0.3359109428 ,1.797028e-01
 ,0.84924411 ,0.650551157 ,0.4878352970 ,0.3432352168 ,1.847478e-01
 ,0.85426493 ,0.658159827 ,0.4939669306 ,0.3514380763 ,1.899971e-01
 ,0.85905395 ,0.665825808 ,0.5002414181 ,0.3599104495 ,1.956613e-01
 ,0.86383966 ,0.673703674 ,0.5069596693 ,0.3680164470 ,2.017933e-01
 ,0.86878464 ,0.681723132 ,0.5141410970 ,0.3757733259 ,2.081422e-01
 ,0.87369919 ,0.689369371 ,0.5213869817 ,0.3836108960 ,2.144345e-01
 ,0.87848296 ,0.696307162 ,0.5284253727 ,0.3917481718 ,2.205671e-01
 ,0.88333342 ,0.703025539 ,0.5354226066 ,0.3997804393 ,2.267723e-01
 ,0.88852444 ,0.710392130 ,0.5425937450 ,0.4071416647 ,2.334748e-01
 ,0.89422273 ,0.718655581 ,0.5500700208 ,0.4139631116 ,2.405566e-01
 ,0.90022856 ,0.727393600 ,0.5580678239 ,0.4208102626 ,2.474752e-01
 ,0.90634711 ,0.736233079 ,0.5663716281 ,0.4277100246 ,2.541958e-01
 ,0.91250898 ,0.745181126 ,0.5745840700 ,0.4341555467 ,2.613699e-01
 ,0.91826744 ,0.754393531 ,0.5831569010 ,0.4404809807 ,2.695945e-01
 ,0.92342137 ,0.763619264 ,0.5926394880 ,0.4480172195 ,2.787520e-01
 ,0.92839321 ,0.772228663 ,0.6031756664 ,0.4576132755 ,2.886374e-01
 ,0.93344816 ,0.780039421 ,0.6154126922 ,0.4690720632 ,2.991463e-01
 ,0.93844307 ,0.787892760 ,0.6293505344 ,0.4817779858 ,3.100554e-01
 ,0.94305951 ,0.796994339 ,0.6430075903 ,0.4942618672 ,3.217248e-01
 ,0.94732447 ,0.807084776 ,0.6556595336 ,0.5058265444 ,3.343702e-01
 ,0.95204731 ,0.817744490 ,0.6691307167 ,0.5191457332 ,3.476188e-01
 ,0.95754260 ,0.828993245 ,0.6843916089 ,0.5349257669 ,3.612037e-01
 ,0.96262769 ,0.840489010 ,0.7010625113 ,0.5506951842 ,3.762375e-01
 ,0.96744559 ,0.851185162 ,0.7179714379 ,0.5654664242 ,3.939733e-01
 ,0.97198841 ,0.862041313 ,0.7355142033 ,0.5808629729 ,4.145818e-01
 ,0.97581215 ,0.873028712 ,0.7534527187 ,0.6022223081 ,4.363161e-01
 ,0.98074122 ,0.884033783 ,0.7720247523 ,0.6336556646 ,4.631840e-01
 ,0.98656778 ,0.899552193 ,0.7975923713 ,0.6670769076 ,4.902773e-01
 ,0.99129231 ,0.924710622 ,0.8318770069 ,0.7077418476 ,5.161952e-01
 ,0.99566295 ,0.958284024 ,0.8782468473 ,0.7720592633 ,5.667889e-01),
 byrow=TRUE,ncol=5)

if(is.null(PV)){
if(!com.p.dist){
if(p==3)rem=P3
if(p==4)rem=P4
if(p==5)rem=P5
if(p==6)rem=P6
} }

est=NA
for(j in 1:p)est[j]=corfun(x[,j],y,...)$cor
id=which(est==max(est))
R=order(est,decreasing=LARGEST)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
bvec<-mclapply(data,corCOMmcp_sub,x,y,corfun,...)
}
if(!MC)bvec<-lapply(data,corCOMmcp_sub,x,y,corfun,...)
output=matrix(NA,nrow=pm1,ncol=8)
if(LARGEST)dimnames(output)=list(NULL,c('IV','Est.', 'Largest.Est', 'Dif','ci.low','ci.hi','p.value','adj.p.value'))
if(!LARGEST)dimnames(output)=list(NULL,c('IV','Est.', 'Smallest.Est', 'Dif','ci.low','ci.hi','p.value','adj.p.value'))
mat=matrix(NA,nrow=nboot,ncol=p)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
for(i in 1:nboot)mat[i,]=bvec[[i]]
for(j in 2:p){
k=j-1
output[k,1]=R[j]
output[k,3]=est[R[1]]
output[k,2]=est[R[j]]
bsort<-sort(mat[,R[1]]-mat[,R[j]])
output[k,4]=est[R[1]]-est[R[j]]
output[k,5]=bsort[ilow]
output[k,6]=bsort[ihi]
pv=mean(bsort<0)+.5*mean(bsort==0)
output[k,7]=2*min(c(pv,1-pv))
flag=output[k,7]>=rem[,k]
ID=which(flag==TRUE)
ic=max(ID,1)
output[k,7]=L[ic]
}
Best='No Decision'
CH=R[1]
names(CH)='IV.w.Largest.Est'
if(!LARGEST)names(CH)='IV.w.Smallest.Est'
output[,8]=p.adjust(output[,7],method=FWE.method)
if(sum(output[,8]<=alpha)==pm1)Best='Decide'
list(CH,Conclusion=Best,results=output)
}

corREGorder<-function(x,y,com.p.dist=FALSE, corfun=wincor,iter=1000,PV=NULL,pr=TRUE,
alpha=.05,nboot=500,SEED=TRUE,MC=FALSE,xout=FALSE,outfun=outpro,method='hoch',...){
#
# Regresion:
#
# Have two or more independent variables, compare
# cor(y,x_I) to cor(y,x_k) for all k!=I, where
# cor(i,x_I) is the highest
# Winsorized correlation is used by default.
# Hochberg's method is used to control FWE.
#
# x is assumed to be a matrix or data frame
#
#
if(nrow(x)!=length(y))stop('x and y have different sample sizes; should be equal')
p=ncol(x)
pm1=p-1
p1=p+1
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
x=m1[,1:pm1]
y=m1[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
nval=nrow(x)

if(pr){
if(!com.p.dist){
if(is.null(PV)){
if(p>6  || nval>350)print('Might need to use com.p.dist=TRUE')
}}}

if(is.null(rem))rem=corREGorder.crit(p,nval,iter=iter)
x<-m1[,1:p]
y=m1[,p1]

L=c(seq(.001,.1,.001),seq(.11,.99,.01))
if(!is.null(PV)){
rem=matrix(NA,length(L),pm1)
for(k in 1:pm1){
for(i in 1:length(L))rem[i,k]=hd(PV[,k],q=L[i])
}}

if(p>6)com.p.dist=TRUE
if(com.p.dist){
if(is.null(PV)){
if(pr)print('Computing the null distribution can take several minutes')
PV=corREGorder.crit(p,nval,iter=iter,MC=MC)
rem=matrix(NA,length(L),pm1)
for(k in 1:pm1){
for(i in 1:length(L))rem[i,k]=hd(PV[,k],q=L[i])
}}}
P3=matrix(c(0.002175049 ,0.01338189
 ,0.005872039 ,0.02090430
 ,0.010619227 ,0.02738247
 ,0.016286140 ,0.03219214
 ,0.022398371 ,0.03579440
 ,0.028262446 ,0.03877950
 ,0.033345906 ,0.04156497
 ,0.037456504 ,0.04436593
 ,0.040694218 ,0.04726051
 ,0.043307931 ,0.05025906
 ,0.045568092 ,0.05334328
 ,0.047696567 ,0.05648043
 ,0.049845675 ,0.05962933
 ,0.052103950 ,0.06274849
 ,0.054511258 ,0.06580645
 ,0.057074761 ,0.06878943
 ,0.059782746 ,0.07170227
 ,0.062615440 ,0.07456228
 ,0.065552496 ,0.07738946
 ,0.068577121 ,0.08019727
 ,0.071677317 ,0.08298752
 ,0.074845101 ,0.08575039
 ,0.078074759 ,0.08846894
 ,0.081360990 ,0.09112571
 ,0.084697501 ,0.09370929
 ,0.088076254 ,0.09621874
 ,0.091487319 ,0.09866517
 ,0.094919200 ,0.10107020
 ,0.098359413 ,0.10346231
 ,0.101795143 ,0.10587191
 ,0.105213802 ,0.10832639
 ,0.108603375 ,0.11084604
 ,0.111952481 ,0.11344138
 ,0.115250189 ,0.11611213
 ,0.118485658 ,0.11884783
 ,0.121647762 ,0.12162960
 ,0.124724850 ,0.12443290
 ,0.127704758 ,0.12723046
 ,0.130575140 ,0.12999536
 ,0.133324080 ,0.13270348
 ,0.135940914 ,0.13533548
 ,0.138417122 ,0.13787790
 ,0.140747141 ,0.14032355
 ,0.142928988 ,0.14267126
 ,0.144964592 ,0.14492507
 ,0.146859829 ,0.14709314
 ,0.148624253 ,0.14918647
 ,0.150270597 ,0.15121765
 ,0.151814117 ,0.15319972
 ,0.153271855 ,0.15514520
 ,0.154661905 ,0.15706537
 ,0.156002725 ,0.15896974
 ,0.157312540 ,0.16086577
 ,0.158608833 ,0.16275881
 ,0.159907941 ,0.16465206
 ,0.161224730 ,0.16654682
 ,0.162572340 ,0.16844268
 ,0.163961992 ,0.17033786
 ,0.165402841 ,0.17222951
 ,0.166901885 ,0.17411410
 ,0.168463921 ,0.17598769
 ,0.170091558 ,0.17784632
 ,0.171785292 ,0.17968623
 ,0.173543642 ,0.18150411
 ,0.175363353 ,0.18329733
 ,0.177239642 ,0.18506405
 ,0.179166497 ,0.18680333
 ,0.181137006 ,0.18851520
 ,0.183143686 ,0.19020068
 ,0.185178820 ,0.19186171
 ,0.187234759 ,0.19350113
 ,0.189304201 ,0.19512259
 ,0.191380404 ,0.19673036
 ,0.193457363 ,0.19832922
 ,0.195529910 ,0.19992426
 ,0.197593772 ,0.20152065
 ,0.199645563 ,0.20312347
 ,0.201682735 ,0.20473749
 ,0.203703491 ,0.20636692
 ,0.205706665 ,0.20801528
 ,0.207691581 ,0.20968521
 ,0.209657914 ,0.21137841
 ,0.211605535 ,0.21309554
 ,0.213534385 ,0.21483623
 ,0.215444346 ,0.21659918
 ,0.217335153 ,0.21838219
 ,0.219206331 ,0.22018238
 ,0.221057158 ,0.22199635
 ,0.222886667 ,0.22382040
 ,0.224693687 ,0.22565074
 ,0.226476899 ,0.22748367
 ,0.228234929 ,0.22931582
 ,0.229966452 ,0.23114424
 ,0.231670301 ,0.23296658
 ,0.233345582 ,0.23478109
 ,0.234991770 ,0.23658668
 ,0.236608793 ,0.23838287
 ,0.238197083 ,0.24016976
 ,0.239757602 ,0.24194791
 ,0.241291828 ,0.24371824
 ,0.255725281 ,0.26121652
 ,0.269643155 ,0.27818483
 ,0.282572599 ,0.29361598
 ,0.294230675 ,0.30801197
 ,0.305279532 ,0.32217337
 ,0.316231246 ,0.33525837
 ,0.327118515 ,0.34679924
 ,0.337725949 ,0.35781816
 ,0.348507973 ,0.36999019
 ,0.360275335 ,0.38380923
 ,0.372732081 ,0.39753304
 ,0.384865494 ,0.40899196
 ,0.396170468 ,0.41839542
 ,0.406525702 ,0.42826562
 ,0.416179917 ,0.44045085
 ,0.425830260 ,0.45430132
 ,0.435731315 ,0.46791777
 ,0.445574831 ,0.47988279
 ,0.455403556 ,0.48962935
 ,0.465473322 ,0.49746681
 ,0.475493015 ,0.50444076
 ,0.485104441 ,0.51166945
 ,0.494622601 ,0.51972051
 ,0.504479896 ,0.52847494
 ,0.514530980 ,0.53743526
 ,0.524379809 ,0.54615592
 ,0.533811067 ,0.55455128
 ,0.542964071 ,0.56289375
 ,0.552415670 ,0.57149093
 ,0.562702410 ,0.58046156
 ,0.573662536 ,0.58975098
 ,0.584499138 ,0.59931329
 ,0.594424353 ,0.60929957
 ,0.603091449 ,0.61973844
 ,0.610580104 ,0.63005679
 ,0.617291412 ,0.63940637
 ,0.623819392 ,0.64740125
 ,0.630625766 ,0.65425937
 ,0.637780329 ,0.66050475
 ,0.645128388 ,0.66675110
 ,0.652676811 ,0.67358269
 ,0.660706577 ,0.68141986
 ,0.669513580 ,0.69030010
 ,0.679074179 ,0.69978580
 ,0.688924859 ,0.70926800
 ,0.698359030 ,0.71838801
 ,0.706803246 ,0.72709029
 ,0.714138033 ,0.73536725
 ,0.720773425 ,0.74316348
 ,0.727376023 ,0.75049022
 ,0.734427437 ,0.75739163
 ,0.741910594 ,0.76379473
 ,0.749301181 ,0.76965778
 ,0.756038316 ,0.77526923
 ,0.762108143 ,0.78110525
 ,0.767991314 ,0.78741712
 ,0.774125688 ,0.79407849
 ,0.780642168 ,0.80072386
 ,0.787654735 ,0.80714732
 ,0.795560943 ,0.81355663
 ,0.804494568 ,0.82019692
 ,0.813689728 ,0.82693733
 ,0.822246533 ,0.83343489
 ,0.830074100 ,0.83948569
 ,0.837601422 ,0.84515904
 ,0.845257402 ,0.85070450
 ,0.853115769 ,0.85637797
 ,0.860852876 ,0.86248746
 ,0.868387995 ,0.86942987
 ,0.875848963 ,0.87706919
 ,0.882870692 ,0.88443942
 ,0.889442420 ,0.89118298
 ,0.896636569 ,0.89779842
 ,0.904626820 ,0.90423655
 ,0.911772218 ,0.91069863
 ,0.917549883 ,0.91754899
 ,0.923211670 ,0.92378977
 ,0.928746246 ,0.92901361
 ,0.934160521 ,0.93466485
 ,0.941177134 ,0.94142063
 ,0.949363922 ,0.94828957
 ,0.956847590 ,0.95420934
 ,0.963599180 ,0.95946403
 ,0.969055149 ,0.96451675
 ,0.973136108 ,0.97019334
 ,0.976868136 ,0.97565816
 ,0.981258840 ,0.98193948
 ,0.988264909 ,0.98843665
 ,0.994433844 ,0.99456404),
 byrow=TRUE,ncol=2)

P4=matrix(c( 0.01551012 ,0.04057726 ,0.04378631
 ,0.02786269 ,0.06546848 ,0.04944117
 ,0.03553717 ,0.08711494 ,0.05573588
 ,0.04023774 ,0.10496114 ,0.06252382
 ,0.04411887 ,0.12025416 ,0.06944572
 ,0.04807201 ,0.13340708 ,0.07602629
 ,0.05230897 ,0.14435931 ,0.08201922
 ,0.05685652 ,0.15314426 ,0.08742154
 ,0.06171756 ,0.16004283 ,0.09233439
 ,0.06689524 ,0.16549047 ,0.09686627
 ,0.07238520 ,0.16994001 ,0.10110915
 ,0.07816280 ,0.17377046 ,0.10514421
 ,0.08417136 ,0.17725430 ,0.10904320
 ,0.09031706 ,0.18056527 ,0.11285974
 ,0.09647542 ,0.18380356 ,0.11661935
 ,0.10250880 ,0.18702173 ,0.12031681
 ,0.10828903 ,0.19024325 ,0.12392259
 ,0.11371720 ,0.19347300 ,0.12739517
 ,0.11873498 ,0.19670257 ,0.13069417
 ,0.12332546 ,0.19991412 ,0.13379041
 ,0.12750572 ,0.20308477 ,0.13667102
 ,0.13131499 ,0.20619197 ,0.13933985
 ,0.13480237 ,0.20921901 ,0.14181443
 ,0.13801702 ,0.21215938 ,0.14412125
 ,0.14100194 ,0.21501886 ,0.14629078
 ,0.14379136 ,0.21781520 ,0.14835323
 ,0.14641091 ,0.22057551 ,0.15033563
 ,0.14887943 ,0.22333212 ,0.15226034
 ,0.15121153 ,0.22611802 ,0.15414467
 ,0.15342004 ,0.22896253 ,0.15600151
 ,0.15551793 ,0.23188795 ,0.15784045
 ,0.15751958 ,0.23490765 ,0.15966900
 ,0.15944129 ,0.23802542 ,0.16149379
 ,0.16130115 ,0.24123624 ,0.16332148
 ,0.16311848 ,0.24452789 ,0.16515915
 ,0.16491299 ,0.24788319 ,0.16701445
 ,0.16670375 ,0.25128243 ,0.16889527
 ,0.16850811 ,0.25470562 ,0.17080920
 ,0.17034067 ,0.25813429 ,0.17276281
 ,0.17221242 ,0.26155271 ,0.17476099
 ,0.17413006 ,0.26494845 ,0.17680637
 ,0.17609563 ,0.26831228 ,0.17889891
 ,0.17810649 ,0.27163767 ,0.18103585
 ,0.18015567 ,0.27491998 ,0.18321179
 ,0.18223251 ,0.27815551 ,0.18541909
 ,0.18432367 ,0.28134075 ,0.18764844
 ,0.18641426 ,0.28447176 ,0.18988943
 ,0.18848908 ,0.28754393 ,0.19213127
 ,0.19053376 ,0.29055199 ,0.19436336
 ,0.19253579 ,0.29349030 ,0.19657586
 ,0.19448531 ,0.29635330 ,0.19876015
 ,0.19637555 ,0.29913601 ,0.20090917
 ,0.19820304 ,0.30183457 ,0.20301766
 ,0.19996751 ,0.30444661 ,0.20508235
 ,0.20167155 ,0.30697148 ,0.20710203
 ,0.20332016 ,0.30941038 ,0.20907756
 ,0.20492013 ,0.31176625 ,0.21101180
 ,0.20647945 ,0.31404361 ,0.21290945
 ,0.20800667 ,0.31624825 ,0.21477685
 ,0.20951041 ,0.31838693 ,0.21662164
 ,0.21099886 ,0.32046697 ,0.21845237
 ,0.21247950 ,0.32249597 ,0.22027815
 ,0.21395877 ,0.32448147 ,0.22210809
 ,0.21544200 ,0.32643073 ,0.22395091
 ,0.21693331 ,0.32835049 ,0.22581448
 ,0.21843563 ,0.33024683 ,0.22770540
 ,0.21995078 ,0.33212513 ,0.22962873
 ,0.22147965 ,0.33398991 ,0.23158777
 ,0.22302229 ,0.33584491 ,0.23358389
 ,0.22457818 ,0.33769307 ,0.23561658
 ,0.22614641 ,0.33953655 ,0.23768352
 ,0.22772593 ,0.34137682 ,0.23978072
 ,0.22931573 ,0.34321472 ,0.24190279
 ,0.23091509 ,0.34505052 ,0.24404323
 ,0.23252371 ,0.34688402 ,0.24619468
 ,0.23414188 ,0.34871459 ,0.24834929
 ,0.23577053 ,0.35054128 ,0.25049895
 ,0.23741133 ,0.35236284 ,0.25263561
 ,0.23906664 ,0.35417779 ,0.25475146
 ,0.24073945 ,0.35598448 ,0.25683915
 ,0.24243332 ,0.35778113 ,0.25889193
 ,0.24415220 ,0.35956584 ,0.26090380
 ,0.24590024 ,0.36133667 ,0.26286960
 ,0.24768166 ,0.36309169 ,0.26478508
 ,0.24950046 ,0.36482898 ,0.26664698
 ,0.25136029 ,0.36654674 ,0.26845307
 ,0.25326420 ,0.36824330 ,0.27020217
 ,0.25521450 ,0.36991721 ,0.27189420
 ,0.25721259 ,0.37156729 ,0.27353012
 ,0.25925888 ,0.37319266 ,0.27511196
 ,0.26135269 ,0.37479283 ,0.27664276
 ,0.26349226 ,0.37636767 ,0.27812647
 ,0.26567472 ,0.37791746 ,0.27956791
 ,0.26789620 ,0.37944287 ,0.28097264
 ,0.27015188 ,0.38094492 ,0.28234681
 ,0.27243611 ,0.38242495 ,0.28369706
 ,0.27474260 ,0.38388453 ,0.28503033
 ,0.27706456 ,0.38532538 ,0.28635373
 ,0.27939486 ,0.38674928 ,0.28767441
 ,0.28172621 ,0.38815796 ,0.28899935
 ,0.30368556 ,0.40162893 ,0.30368088
 ,0.32094974 ,0.41407715 ,0.32178330
 ,0.33413706 ,0.42686681 ,0.33940421
 ,0.34519638 ,0.44106918 ,0.35258534
 ,0.35553655 ,0.45434941 ,0.36200235
 ,0.36622977 ,0.46608461 ,0.37040153
 ,0.37786945 ,0.47769452 ,0.37925941
 ,0.39037354 ,0.48941876 ,0.38854815
 ,0.40339828 ,0.50079642 ,0.39825958
 ,0.41579925 ,0.51188178 ,0.40920109
 ,0.42617477 ,0.52311541 ,0.42158463
 ,0.43455327 ,0.53487858 ,0.43453761
 ,0.44208791 ,0.54714403 ,0.44723309
 ,0.45001905 ,0.55923271 ,0.45932100
 ,0.45896481 ,0.57027401 ,0.47072048
 ,0.46860843 ,0.57999513 ,0.48143650
 ,0.47834450 ,0.58893700 ,0.49173131
 ,0.48781034 ,0.59778542 ,0.50198123
 ,0.49695736 ,0.60684063 ,0.51232433
 ,0.50614965 ,0.61608565 ,0.52269385
 ,0.51570653 ,0.62525226 ,0.53300494
 ,0.52544301 ,0.63382683 ,0.54336275
 ,0.53507224 ,0.64135887 ,0.55402012
 ,0.54466956 ,0.64795069 ,0.56477936
 ,0.55450045 ,0.65423352 ,0.57492627
 ,0.56470767 ,0.66082598 ,0.58412311
 ,0.57520589 ,0.66798477 ,0.59274171
 ,0.58563795 ,0.67557212 ,0.60127760
 ,0.59553109 ,0.68316811 ,0.60997760
 ,0.60465763 ,0.69037519 ,0.61886863
 ,0.61309294 ,0.69714665 ,0.62782311
 ,0.62096747 ,0.70381712 ,0.63662138
 ,0.62843358 ,0.71066866 ,0.64498294
 ,0.63586964 ,0.71749020 ,0.65265906
 ,0.64375353 ,0.72384456 ,0.65957668
 ,0.65212793 ,0.72967232 ,0.66583193
 ,0.66049812 ,0.73529566 ,0.67159225
 ,0.66842274 ,0.74091198 ,0.67713840
 ,0.67595398 ,0.74641608 ,0.68293746
 ,0.68342673 ,0.75170365 ,0.68940237
 ,0.69097598 ,0.75681883 ,0.69648797
 ,0.69842249 ,0.76180033 ,0.70370713
 ,0.70563560 ,0.76670054 ,0.71064719
 ,0.71281012 ,0.77173041 ,0.71735349
 ,0.72033581 ,0.77705000 ,0.72412267
 ,0.72848382 ,0.78241378 ,0.73111595
 ,0.73708729 ,0.78733284 ,0.73832883
 ,0.74551054 ,0.79160726 ,0.74567870
 ,0.75319837 ,0.79552364 ,0.75291314
 ,0.76015349 ,0.79957172 ,0.75966877
 ,0.76670769 ,0.80408855 ,0.76584461
 ,0.77313827 ,0.80924036 ,0.77174026
 ,0.77962403 ,0.81516519 ,0.77765648
 ,0.78632238 ,0.82177046 ,0.78360390
 ,0.79330110 ,0.82852505 ,0.78963534
 ,0.80032049 ,0.83486410 ,0.79614367
 ,0.80690242 ,0.84071720 ,0.80341723
 ,0.81285318 ,0.84643390 ,0.81116026
 ,0.81843087 ,0.85232962 ,0.81884205
 ,0.82389821 ,0.85834921 ,0.82628132
 ,0.82925515 ,0.86412609 ,0.83357420
 ,0.83452095 ,0.86947172 ,0.84068096
 ,0.84009069 ,0.87468947 ,0.84740433
 ,0.84641252 ,0.88012125 ,0.85369747
 ,0.85337749 ,0.88550274 ,0.85984420
 ,0.86056675 ,0.89031731 ,0.86629374
 ,0.86782542 ,0.89466471 ,0.87323085
 ,0.87514888 ,0.89925142 ,0.88037505
 ,0.88234515 ,0.90465670 ,0.88731009
 ,0.88906611 ,0.91058956 ,0.89382801
 ,0.89489509 ,0.91600949 ,0.89989887
 ,0.89964387 ,0.92037164 ,0.90556334
 ,0.90392518 ,0.92425694 ,0.91107573
 ,0.90874993 ,0.92837270 ,0.91667095
 ,0.91433102 ,0.93276449 ,0.92229467
 ,0.92009330 ,0.93716660 ,0.92780387
 ,0.92614330 ,0.94129808 ,0.93326004
 ,0.93274464 ,0.94523346 ,0.93882150
 ,0.93904729 ,0.94970628 ,0.94451555
 ,0.94498605 ,0.95507385 ,0.95042484
 ,0.95172002 ,0.96051456 ,0.95623595
 ,0.95904715 ,0.96546866 ,0.96109440
 ,0.96528014 ,0.96983221 ,0.96583206
 ,0.96981308 ,0.97446466 ,0.97124046
 ,0.97390026 ,0.97945170 ,0.97556419
 ,0.97783109 ,0.98385517 ,0.97915772
 ,0.98192499 ,0.98887894 ,0.98323212
 ,0.98701317 ,0.99271102 ,0.98660775
 ,0.99226134 ,0.99482149 ,0.99079872),
 byrow=TRUE,ncol=3)

 P5=matrix(c( 0.05652412 ,0.1064700 ,0.05210169 ,0.005493583
 ,0.06574073 ,0.1193620 ,0.06365731 ,0.013235067
 ,0.07151676 ,0.1298226 ,0.07635689 ,0.024181097
 ,0.07511946 ,0.1389007 ,0.09074149 ,0.036214033
 ,0.07828081 ,0.1482306 ,0.10635395 ,0.047401883
 ,0.08180130 ,0.1583439 ,0.12166792 ,0.056974651
 ,0.08583995 ,0.1691003 ,0.13539360 ,0.065203466
 ,0.09025786 ,0.1802110 ,0.14713767 ,0.072752928
 ,0.09481793 ,0.1914175 ,0.15726850 ,0.080170849
 ,0.09930536 ,0.2024707 ,0.16646093 ,0.087701224
 ,0.10358876 ,0.2130961 ,0.17530679 ,0.095324733
 ,0.10762929 ,0.2230221 ,0.18413231 ,0.102886098
 ,0.11145588 ,0.2320476 ,0.19299485 ,0.110212699
 ,0.11512861 ,0.2400992 ,0.20177360 ,0.117183521
 ,0.11870678 ,0.2472403 ,0.21027870 ,0.123747085
 ,0.12222956 ,0.2536340 ,0.21833522 ,0.129906857
 ,0.12571005 ,0.2594800 ,0.22582830 ,0.135695457
 ,0.12913936 ,0.2649560 ,0.23271445 ,0.141152226
 ,0.13249584 ,0.2701803 ,0.23901157 ,0.146309897
 ,0.13575514 ,0.2752025 ,0.24477977 ,0.151189960
 ,0.13889807 ,0.2800167 ,0.25010143 ,0.155803801
 ,0.14191516 ,0.2845862 ,0.25506478 ,0.160156707
 ,0.14480785 ,0.2888690 ,0.25975221 ,0.164252701
 ,0.14758751 ,0.2928368 ,0.26423308 ,0.168098998
 ,0.15027304 ,0.2964849 ,0.26856022 ,0.171709323
 ,0.15288821 ,0.2998329 ,0.27276929 ,0.175105614
 ,0.15545896 ,0.3029190 ,0.27688008 ,0.178317905
 ,0.15801110 ,0.3057918 ,0.28089924 ,0.181382526
 ,0.16056820 ,0.3085024 ,0.28482365 ,0.184339066
 ,0.16315001 ,0.3110976 ,0.28864393 ,0.187226779
 ,0.16577118 ,0.3136161 ,0.29234767 ,0.190081080
 ,0.16844049 ,0.3160868 ,0.29592209 ,0.192930735
 ,0.17116069 ,0.3185293 ,0.29935598 ,0.195796047
 ,0.17392879 ,0.3209550 ,0.30264093 ,0.198688172
 ,0.17673683 ,0.3233687 ,0.30577201 ,0.201609481
 ,0.17957306 ,0.3257704 ,0.30874798 ,0.204554768
 ,0.18242312 ,0.3281569 ,0.31157121 ,0.207513022
 ,0.18527142 ,0.3305232 ,0.31424738 ,0.210469505
 ,0.18810221 ,0.3328636 ,0.31678506 ,0.213407838
 ,0.19090055 ,0.3351728 ,0.31919516 ,0.216311897
 ,0.19365295 ,0.3374465 ,0.32149031 ,0.219167333
 ,0.19634778 ,0.3396823 ,0.32368416 ,0.221962623
 ,0.19897548 ,0.3418799 ,0.32579080 ,0.224689627
 ,0.20152856 ,0.3440415 ,0.32782411 ,0.227343686
 ,0.20400154 ,0.3461714 ,0.32979739 ,0.229923353
 ,0.20639083 ,0.3482759 ,0.33172307 ,0.232429882
 ,0.20869457 ,0.3503624 ,0.33361256 ,0.234866606
 ,0.21091258 ,0.3524387 ,0.33547633 ,0.237238310
 ,0.21304622 ,0.3545126 ,0.33732401 ,0.239550684
 ,0.21509838 ,0.3565905 ,0.33916450 ,0.241809893
 ,0.21707346 ,0.3586772 ,0.34100612 ,0.244022277
 ,0.21897727 ,0.3607754 ,0.34285661 ,0.246194145
 ,0.22081697 ,0.3628855 ,0.34472307 ,0.248331635
 ,0.22260093 ,0.3650055 ,0.34661180 ,0.250440593
 ,0.22433851 ,0.3671314 ,0.34852808 ,0.252526433
 ,0.22603978 ,0.3692577 ,0.35047594 ,0.254593979
 ,0.22771524 ,0.3713776 ,0.35245788 ,0.256647268
 ,0.22937544 ,0.3734840 ,0.35447483 ,0.258689342
 ,0.23103064 ,0.3755696 ,0.35652605 ,0.260722057
 ,0.23269046 ,0.3776277 ,0.35860928 ,0.262745921
 ,0.23436362 ,0.3796525 ,0.36072093 ,0.264760011
 ,0.23605766 ,0.3816393 ,0.36285653 ,0.266761969
 ,0.23777881 ,0.3835848 ,0.36501105 ,0.268748094
 ,0.23953181 ,0.3854871 ,0.36717945 ,0.270713533
 ,0.24131994 ,0.3873455 ,0.36935709 ,0.272652546
 ,0.24314491 ,0.3891604 ,0.37154010 ,0.274558844
 ,0.24500695 ,0.3909334 ,0.37372567 ,0.276425957
 ,0.24690483 ,0.3926667 ,0.37591220 ,0.278247625
 ,0.24883595 ,0.3943630 ,0.37809926 ,0.280018169
 ,0.25079643 ,0.3960256 ,0.38028752 ,0.281732830
 ,0.25278126 ,0.3976578 ,0.38247846 ,0.283388055
 ,0.25478439 ,0.3992633 ,0.38467408 ,0.284981717
 ,0.25679899 ,0.4008454 ,0.38687646 ,0.286513250
 ,0.25881759 ,0.4024077 ,0.38908746 ,0.287983708
 ,0.26083237 ,0.4039536 ,0.39130828 ,0.289395740
 ,0.26283537 ,0.4054864 ,0.39353919 ,0.290753487
 ,0.26481879 ,0.4070095 ,0.39577932 ,0.292062415
 ,0.26677527 ,0.4085258 ,0.39802648 ,0.293329080
 ,0.26869809 ,0.4100386 ,0.40027722 ,0.294560863
 ,0.27058146 ,0.4115507 ,0.40252683 ,0.295765673
 ,0.27242063 ,0.4130649 ,0.40476958 ,0.296951643
 ,0.27421212 ,0.4145841 ,0.40699893 ,0.298126831
 ,0.27595371 ,0.4161107 ,0.40920786 ,0.299298952
 ,0.27764451 ,0.4176472 ,0.41138921 ,0.300475144
 ,0.27928490 ,0.4191957 ,0.41353603 ,0.301661783
 ,0.28087643 ,0.4207583 ,0.41564194 ,0.302864352
 ,0.28242170 ,0.4223369 ,0.41770144 ,0.304087372
 ,0.28392419 ,0.4239330 ,0.41971011 ,0.305334383
 ,0.28538803 ,0.4255480 ,0.42166488 ,0.306607981
 ,0.28681780 ,0.4271829 ,0.42356406 ,0.307909888
 ,0.28821831 ,0.4288386 ,0.42540740 ,0.309241059
 ,0.28959443 ,0.4305155 ,0.42719605 ,0.310601803
 ,0.29095083 ,0.4322137 ,0.42893242 ,0.311991913
 ,0.29229189 ,0.4339327 ,0.43061997 ,0.313410782
 ,0.29362151 ,0.4356718 ,0.43226306 ,0.314857513
 ,0.29494305 ,0.4374295 ,0.43386662 ,0.316331000
 ,0.29625927 ,0.4392040 ,0.43543592 ,0.317829987
 ,0.29757232 ,0.4409928 ,0.43697631 ,0.319353103
 ,0.29888374 ,0.4427930 ,0.43849299 ,0.320898866
 ,0.30019451 ,0.4446010 ,0.43999077 ,0.322465682
 ,0.31328362 ,0.4621878 ,0.45451950 ,0.338826479
 ,0.32658914 ,0.4768889 ,0.46817448 ,0.354761127
 ,0.34081183 ,0.4889837 ,0.48092037 ,0.369582715
 ,0.35529268 ,0.4997010 ,0.49462451 ,0.384419904
 ,0.36855379 ,0.5096083 ,0.51023274 ,0.398755180
 ,0.37994080 ,0.5184617 ,0.52647200 ,0.411411797
 ,0.39062100 ,0.5263287 ,0.54163677 ,0.423038765
 ,0.40257143 ,0.5344950 ,0.55461543 ,0.434627278
 ,0.41621275 ,0.5440195 ,0.56505570 ,0.446414188
 ,0.43046818 ,0.5547192 ,0.57379461 ,0.458794883
 ,0.44486669 ,0.5653656 ,0.58229160 ,0.471543195
 ,0.45923530 ,0.5750124 ,0.59101203 ,0.483555129
 ,0.47296930 ,0.5839699 ,0.59924570 ,0.494652841
 ,0.48590094 ,0.5928607 ,0.60697069 ,0.505540897
 ,0.49834493 ,0.6017113 ,0.61531580 ,0.516236508
 ,0.51026131 ,0.6105778 ,0.62505014 ,0.526149184
 ,0.52119211 ,0.6200780 ,0.63576320 ,0.535122517
 ,0.53082589 ,0.6304166 ,0.64643772 ,0.543660340
 ,0.53943126 ,0.6407196 ,0.65638365 ,0.552285332
 ,0.54774724 ,0.6500929 ,0.66564380 ,0.561135243
 ,0.55642024 ,0.6586570 ,0.67450792 ,0.570308867
 ,0.56550181 ,0.6670994 ,0.68302769 ,0.580053861
 ,0.57464883 ,0.6756273 ,0.69143503 ,0.590349969
 ,0.58368440 ,0.6838335 ,0.70025072 ,0.600590044
 ,0.59259791 ,0.6914073 ,0.70958149 ,0.610032910
 ,0.60135765 ,0.6984377 ,0.71882369 ,0.618579939
 ,0.60997532 ,0.7051144 ,0.72721163 ,0.626728815
 ,0.61843676 ,0.7115369 ,0.73447399 ,0.634788345
 ,0.62674513 ,0.7178305 ,0.74089008 ,0.642664491
 ,0.63517287 ,0.7242564 ,0.74693265 ,0.650337943
 ,0.64411712 ,0.7310644 ,0.75290477 ,0.658048604
 ,0.65355691 ,0.7382701 ,0.75884913 ,0.665979304
 ,0.66286128 ,0.7456401 ,0.76482766 ,0.674038165
 ,0.67142296 ,0.7528643 ,0.77103158 ,0.682044900
 ,0.67945227 ,0.7597670 ,0.77745500 ,0.689918781
 ,0.68784816 ,0.7663382 ,0.78372370 ,0.697608584
 ,0.69717228 ,0.7725694 ,0.78940484 ,0.705026668
 ,0.70693776 ,0.7783837 ,0.79438682 ,0.712091065
 ,0.71609195 ,0.7837138 ,0.79895071 ,0.718801927
 ,0.72398625 ,0.7886294 ,0.80353613 ,0.725335927
 ,0.73066318 ,0.7934308 ,0.80838726 ,0.731990069
 ,0.73650853 ,0.7984516 ,0.81339086 ,0.738959992
 ,0.74191863 ,0.8036594 ,0.81827024 ,0.746236537
 ,0.74721366 ,0.8086236 ,0.82286981 ,0.753653660
 ,0.75262380 ,0.8129806 ,0.82718266 ,0.760952429
 ,0.75824105 ,0.8167987 ,0.83122205 ,0.767851565
 ,0.76409670 ,0.8204615 ,0.83502651 ,0.774128504
 ,0.77030423 ,0.8243454 ,0.83871807 ,0.779748097
 ,0.77696430 ,0.8286333 ,0.84243305 ,0.784919804
 ,0.78392814 ,0.8333291 ,0.84622101 ,0.789918249
 ,0.79087142 ,0.8383467 ,0.85001054 ,0.794848310
 ,0.79764992 ,0.8435563 ,0.85367430 ,0.799645577
 ,0.80440374 ,0.8488329 ,0.85720937 ,0.804300457
 ,0.81125920 ,0.8540607 ,0.86077075 ,0.808937266
 ,0.81807396 ,0.8590979 ,0.86444685 ,0.813668806
 ,0.82455298 ,0.8638895 ,0.86820083 ,0.818594382
 ,0.83054464 ,0.8685076 ,0.87202674 ,0.823874450
 ,0.83619452 ,0.8729580 ,0.87586813 ,0.829667963
 ,0.84182788 ,0.8772068 ,0.87951882 ,0.836051379
 ,0.84769382 ,0.8813966 ,0.88292402 ,0.842860562
 ,0.85379684 ,0.8857317 ,0.88639513 ,0.849605729
 ,0.85982991 ,0.8902404 ,0.89031782 ,0.855810573
 ,0.86534864 ,0.8948524 ,0.89481091 ,0.861500169
 ,0.87020343 ,0.8994987 ,0.89970471 ,0.867139594
 ,0.87472363 ,0.9041210 ,0.90468146 ,0.872950260
 ,0.87942450 ,0.9088721 ,0.90940470 ,0.878638514
 ,0.88437025 ,0.9139336 ,0.91367910 ,0.883816792
 ,0.88904758 ,0.9189854 ,0.91755701 ,0.888324527
 ,0.89322701 ,0.9235781 ,0.92120549 ,0.892386541
 ,0.89745058 ,0.9278084 ,0.92477201 ,0.896563267
 ,0.90239053 ,0.9320707 ,0.92830567 ,0.901221708
 ,0.90800180 ,0.9364567 ,0.93180350 ,0.906121919
 ,0.91346775 ,0.9405897 ,0.93545855 ,0.911051876
 ,0.91813323 ,0.9441004 ,0.93949383 ,0.916500400
 ,0.92221702 ,0.9470722 ,0.94376981 ,0.922684763
 ,0.92622555 ,0.9500208 ,0.94792822 ,0.928942208
 ,0.93049876 ,0.9535137 ,0.95168820 ,0.934711316
 ,0.93559016 ,0.9573979 ,0.95514668 ,0.940041161
 ,0.94158878 ,0.9611006 ,0.95847721 ,0.945997043
 ,0.94715113 ,0.9646607 ,0.96144331 ,0.952681937
 ,0.95168452 ,0.9683622 ,0.96448503 ,0.957991575
 ,0.95643627 ,0.9719767 ,0.96793903 ,0.962249417
 ,0.96150476 ,0.9751701 ,0.97130901 ,0.967311290
 ,0.96619142 ,0.9782121 ,0.97505959 ,0.972704298
 ,0.97232678 ,0.9821594 ,0.97956025 ,0.977392196
 ,0.97800785 ,0.9866658 ,0.98474153 ,0.982039939
 ,0.98168187 ,0.9907911 ,0.98956488 ,0.987315734
 ,0.98661473 ,0.9937451 ,0.99413910 ,0.991894481
 ,0.99161624 ,0.9972316 ,0.99730242 ,0.996798896),
 byrow=TRUE,ncol=4)

 P6=matrix(c(0.02845841 ,0.1323308 ,0.1003219 ,0.1384688 ,0.05477254
 ,0.03862415 ,0.1504455 ,0.1516471 ,0.1577896 ,0.05960211
 ,0.05220254 ,0.1708109 ,0.1948779 ,0.1744245 ,0.06519134
 ,0.06621145 ,0.1900433 ,0.2222013 ,0.1878881 ,0.07081226
 ,0.07867811 ,0.2067988 ,0.2375538 ,0.1993997 ,0.07609186
 ,0.08886734 ,0.2210344 ,0.2471487 ,0.2097566 ,0.08088784
 ,0.09687788 ,0.2332937 ,0.2552052 ,0.2193390 ,0.08524504
 ,0.10316830 ,0.2441869 ,0.2636724 ,0.2284077 ,0.08932865
 ,0.10824058 ,0.2541494 ,0.2731186 ,0.2372217 ,0.09334220
 ,0.11250832 ,0.2634030 ,0.2834439 ,0.2459903 ,0.09745755
 ,0.11627805 ,0.2720146 ,0.2942621 ,0.2547909 ,0.10177306
 ,0.11976844 ,0.2799750 ,0.3051103 ,0.2635490 ,0.10630380
 ,0.12313219 ,0.2872606 ,0.3155798 ,0.2720913 ,0.11099790
 ,0.12647380 ,0.2938675 ,0.3253870 ,0.2802295 ,0.11576708
 ,0.12986454 ,0.2998218 ,0.3343885 ,0.2878280 ,0.12051851
 ,0.13335383 ,0.3051771 ,0.3425552 ,0.2948321 ,0.12517872
 ,0.13697580 ,0.3100060 ,0.3499315 ,0.3012592 ,0.12970519
 ,0.14075114 ,0.3143922 ,0.3565967 ,0.3071707 ,0.13408620
 ,0.14468601 ,0.3184238 ,0.3626399 ,0.3126418 ,0.13833266
 ,0.14877041 ,0.3221871 ,0.3681465 ,0.3177403 ,0.14246686
 ,0.15297783 ,0.3257612 ,0.3731941 ,0.3225181 ,0.14651218
 ,0.15726699 ,0.3292136 ,0.3778517 ,0.3270116 ,0.15048655
 ,0.16158577 ,0.3325961 ,0.3821812 ,0.3312470 ,0.15440006
 ,0.16587656 ,0.3359439 ,0.3862375 ,0.3352459 ,0.15825605
 ,0.17008218 ,0.3392769 ,0.3900693 ,0.3390287 ,0.16205392
 ,0.17415147 ,0.3426022 ,0.3937200 ,0.3426162 ,0.16579194
 ,0.17804364 ,0.3459191 ,0.3972286 ,0.3460302 ,0.16946914
 ,0.18173091 ,0.3492229 ,0.4006301 ,0.3492928 ,0.17308558
 ,0.18519926 ,0.3525091 ,0.4039573 ,0.3524272 ,0.17664150
 ,0.18844742 ,0.3557754 ,0.4072408 ,0.3554574 ,0.18013586
 ,0.19148459 ,0.3590236 ,0.4105092 ,0.3584099 ,0.18356508
 ,0.19432745 ,0.3622586 ,0.4137885 ,0.3613127 ,0.18692243
 ,0.19699700 ,0.3654879 ,0.4171007 ,0.3641954 ,0.19019849
 ,0.19951573 ,0.3687203 ,0.4204629 ,0.3670873 ,0.19338240
 ,0.20190550 ,0.3719641 ,0.4238853 ,0.3700154 ,0.19646371
 ,0.20418614 ,0.3752258 ,0.4273710 ,0.3730020 ,0.19943431
 ,0.20637484 ,0.3785096 ,0.4309154 ,0.3760622 ,0.20229001
 ,0.20848607 ,0.3818167 ,0.4345067 ,0.3792026 ,0.20503150
 ,0.21053195 ,0.3851457 ,0.4381272 ,0.3824203 ,0.20766451
 ,0.21252275 ,0.3884927 ,0.4417550 ,0.3857039 ,0.21019930
 ,0.21446746 ,0.3918524 ,0.4453661 ,0.3890342 ,0.21264952
 ,0.21637417 ,0.3952189 ,0.4489369 ,0.3923868 ,0.21503090
 ,0.21825036 ,0.3985867 ,0.4524457 ,0.3957345 ,0.21735972
 ,0.22010296 ,0.4019509 ,0.4558750 ,0.3990498 ,0.21965153
 ,0.22193833 ,0.4053081 ,0.4592123 ,0.4023072 ,0.22192010
 ,0.22376206 ,0.4086562 ,0.4624507 ,0.4054849 ,0.22417673
 ,0.22557886 ,0.4119939 ,0.4655884 ,0.4085662 ,0.22642994
 ,0.22739239 ,0.4153207 ,0.4686284 ,0.4115399 ,0.22868536
 ,0.22920525 ,0.4186356 ,0.4715768 ,0.4144002 ,0.23094602
 ,0.23101900 ,0.4219366 ,0.4744418 ,0.4171462 ,0.23321268
 ,0.23283434 ,0.4252203 ,0.4772317 ,0.4197810 ,0.23548425
 ,0.23465133 ,0.4284812 ,0.4799538 ,0.4223106 ,0.23775824
 ,0.23646965 ,0.4317114 ,0.4826136 ,0.4247430 ,0.24003124
 ,0.23828893 ,0.4349012 ,0.4852138 ,0.4270871 ,0.24229916
 ,0.24010890 ,0.4380394 ,0.4877545 ,0.4293518 ,0.24455754
 ,0.24192967 ,0.4411139 ,0.4902331 ,0.4315456 ,0.24680175
 ,0.24375172 ,0.4441127 ,0.4926456 ,0.4336759 ,0.24902706
 ,0.24557590 ,0.4470248 ,0.4949865 ,0.4357493 ,0.25122873
 ,0.24740333 ,0.4498409 ,0.4972506 ,0.4377710 ,0.25340210
 ,0.24923515 ,0.4525542 ,0.4994334 ,0.4397457 ,0.25554261
 ,0.25107232 ,0.4551608 ,0.5015319 ,0.4416770 ,0.25764589
 ,0.25291533 ,0.4576600 ,0.5035450 ,0.4435684 ,0.25970787
 ,0.25476396 ,0.4600538 ,0.5054743 ,0.4454232 ,0.26172490
 ,0.25661711 ,0.4623474 ,0.5073235 ,0.4472447 ,0.26369385
 ,0.25847267 ,0.4645479 ,0.5090988 ,0.4490365 ,0.26561225
 ,0.26032752 ,0.4666643 ,0.5108083 ,0.4508025 ,0.26747840
 ,0.26217762 ,0.4687065 ,0.5124617 ,0.4525472 ,0.26929144
 ,0.26401815 ,0.4706847 ,0.5140696 ,0.4542750 ,0.27105140
 ,0.26584379 ,0.4726089 ,0.5156436 ,0.4559907 ,0.27275919
 ,0.26764897 ,0.4744884 ,0.5171951 ,0.4576988 ,0.27441663
 ,0.26942825 ,0.4763314 ,0.5187352 ,0.4594038 ,0.27602636
 ,0.27117659 ,0.4781447 ,0.5202743 ,0.4611098 ,0.27759181
 ,0.27288969 ,0.4799335 ,0.5218217 ,0.4628199 ,0.27911708
 ,0.27456420 ,0.4817019 ,0.5233854 ,0.4645368 ,0.28060691
 ,0.27619794 ,0.4834523 ,0.5249719 ,0.4662620 ,0.28206660
 ,0.27779004 ,0.4851863 ,0.5265858 ,0.4679962 ,0.28350186
 ,0.27934099 ,0.4869043 ,0.5282303 ,0.4697390 ,0.28491882
 ,0.28085262 ,0.4886062 ,0.5299067 ,0.4714888 ,0.28632382
 ,0.28232802 ,0.4902914 ,0.5316147 ,0.4732435 ,0.28772339
 ,0.28377144 ,0.4919591 ,0.5333528 ,0.4749997 ,0.28912404
 ,0.28518812 ,0.4936083 ,0.5351181 ,0.4767537 ,0.29053216
 ,0.28658401 ,0.4952382 ,0.5369065 ,0.4785012 ,0.29195381
 ,0.28796561 ,0.4968483 ,0.5387135 ,0.4802377 ,0.29339456
 ,0.28933969 ,0.4984381 ,0.5405341 ,0.4819586 ,0.29485931
 ,0.29071301 ,0.5000073 ,0.5423631 ,0.4836596 ,0.29635210
 ,0.29209211 ,0.5015561 ,0.5441954 ,0.4853366 ,0.29787602
 ,0.29348304 ,0.5030846 ,0.5460263 ,0.4869864 ,0.29943305
 ,0.29489121 ,0.5045932 ,0.5478516 ,0.4886062 ,0.30102400
 ,0.29632114 ,0.5060823 ,0.5496679 ,0.4901943 ,0.30264856
 ,0.29777639 ,0.5075524 ,0.5514723 ,0.4917498 ,0.30430529
 ,0.29925947 ,0.5090041 ,0.5532629 ,0.4932728 ,0.30599175
 ,0.30077178 ,0.5104381 ,0.5550384 ,0.4947644 ,0.30770464
 ,0.30231363 ,0.5118548 ,0.5567981 ,0.4962267 ,0.30943996
 ,0.30388432 ,0.5132549 ,0.5585419 ,0.4976622 ,0.31119330
 ,0.30548226 ,0.5146391 ,0.5602700 ,0.4990746 ,0.31295994
 ,0.30710505 ,0.5160078 ,0.5619830 ,0.5004678 ,0.31473518
 ,0.30874972 ,0.5173619 ,0.5636814 ,0.5018461 ,0.31651448
 ,0.31041287 ,0.5187020 ,0.5653660 ,0.5032140 ,0.31829362
 ,0.31209087 ,0.5200289 ,0.5670370 ,0.5045763 ,0.32006892
 ,0.31378008 ,0.5213431 ,0.5686948 ,0.5059374 ,0.32183726
 ,0.33073053 ,0.5339107 ,0.5844004 ,0.5201736 ,0.33893240
 ,0.34766587 ,0.5455699 ,0.5977995 ,0.5358098 ,0.35618985
 ,0.36497276 ,0.5565563 ,0.6088557 ,0.5512613 ,0.37599206
 ,0.38214131 ,0.5670711 ,0.6182221 ,0.5653348 ,0.39630223
 ,0.39892707 ,0.5774423 ,0.6269864 ,0.5781859 ,0.41294824
 ,0.41477035 ,0.5887255 ,0.6365914 ,0.5904167 ,0.42604441
 ,0.42985193 ,0.6012620 ,0.6472327 ,0.6023352 ,0.43772949
 ,0.44459125 ,0.6137931 ,0.6578100 ,0.6136908 ,0.44944166
 ,0.45830629 ,0.6247817 ,0.6680975 ,0.6238971 ,0.46148163
 ,0.47040113 ,0.6338022 ,0.6784202 ,0.6325018 ,0.47293114
 ,0.48149178 ,0.6415214 ,0.6879435 ,0.6396385 ,0.48317936
 ,0.49230962 ,0.6486243 ,0.6955894 ,0.6462754 ,0.49290670
 ,0.50293366 ,0.6556637 ,0.7014553 ,0.6530889 ,0.50315292
 ,0.51329043 ,0.6629905 ,0.7067582 ,0.6598945 ,0.51404269
 ,0.52322279 ,0.6704942 ,0.7127512 ,0.6664369 ,0.52474488
 ,0.53250776 ,0.6779803 ,0.7196826 ,0.6730249 ,0.53501309
 ,0.54136758 ,0.6852235 ,0.7269392 ,0.6804293 ,0.54550010
 ,0.55057692 ,0.6920138 ,0.7339733 ,0.6891016 ,0.55652673
 ,0.56083373 ,0.6985598 ,0.7407423 ,0.6983890 ,0.56770747
 ,0.57188476 ,0.7051878 ,0.7473342 ,0.7069275 ,0.57851707
 ,0.58258251 ,0.7118141 ,0.7535933 ,0.7139242 ,0.58875425
 ,0.59206406 ,0.7181377 ,0.7592423 ,0.7197399 ,0.59840912
 ,0.60043277 ,0.7241198 ,0.7641102 ,0.7253388 ,0.60737553
 ,0.60824065 ,0.7302008 ,0.7682440 ,0.7313809 ,0.61559579
 ,0.61587328 ,0.7369912 ,0.7720497 ,0.7377497 ,0.62335187
 ,0.62350556 ,0.7444895 ,0.7761923 ,0.7438837 ,0.63120528
 ,0.63130400 ,0.7518477 ,0.7811232 ,0.7494397 ,0.63946844
 ,0.63950683 ,0.7582114 ,0.7868131 ,0.7545253 ,0.64791259
 ,0.64822576 ,0.7634978 ,0.7929559 ,0.7594398 ,0.65618729
 ,0.65718853 ,0.7682302 ,0.7991861 ,0.7643607 ,0.66415285
 ,0.66586026 ,0.7729213 ,0.8051205 ,0.7692685 ,0.67170795
 ,0.67385499 ,0.7777031 ,0.8105010 ,0.7741140 ,0.67868318
 ,0.68116860 ,0.7824831 ,0.8154188 ,0.7789630 ,0.68495291
 ,0.68817020 ,0.7872791 ,0.8202914 ,0.7838794 ,0.69056808
 ,0.69545456 ,0.7921821 ,0.8255185 ,0.7888123 ,0.69590049
 ,0.70339322 ,0.7971795 ,0.8311156 ,0.7937833 ,0.70157680
 ,0.71168363 ,0.8021996 ,0.8366622 ,0.7989904 ,0.70809578
 ,0.71966391 ,0.8071535 ,0.8415997 ,0.8045273 ,0.71550505
 ,0.72706877 ,0.8119058 ,0.8456326 ,0.8101695 ,0.72339718
 ,0.73411394 ,0.8163800 ,0.8489277 ,0.8155690 ,0.73113224
 ,0.74096543 ,0.8207198 ,0.8519149 ,0.8205280 ,0.73821297
 ,0.74754062 ,0.8252428 ,0.8548854 ,0.8251003 ,0.74459805
 ,0.75389357 ,0.8301610 ,0.8578935 ,0.8295502 ,0.75057740
 ,0.76044477 ,0.8353000 ,0.8609909 ,0.8341226 ,0.75637074
 ,0.76755285 ,0.8401622 ,0.8643055 ,0.8388245 ,0.76200276
 ,0.77509023 ,0.8444077 ,0.8678815 ,0.8434774 ,0.76755578
 ,0.78258820 ,0.8482099 ,0.8716303 ,0.8478628 ,0.77338458
 ,0.78956970 ,0.8520035 ,0.8754454 ,0.8518595 ,0.77990334
 ,0.79578357 ,0.8559626 ,0.8792329 ,0.8555642 ,0.78708237
 ,0.80125054 ,0.8599094 ,0.8828409 ,0.8592369 ,0.79436547
 ,0.80604156 ,0.8636647 ,0.8861239 ,0.8631591 ,0.80125859
 ,0.81021918 ,0.8672377 ,0.8890983 ,0.8674684 ,0.80774304
 ,0.81409451 ,0.8707223 ,0.8919124 ,0.8719953 ,0.81404504
 ,0.81820723 ,0.8742131 ,0.8946943 ,0.8764030 ,0.82035059
 ,0.82291940 ,0.8778013 ,0.8975358 ,0.8804936 ,0.82673907
 ,0.82816550 ,0.8815756 ,0.9005363 ,0.8842174 ,0.83309877
 ,0.83358635 ,0.8855560 ,0.9037186 ,0.8875872 ,0.83911146
 ,0.83888900 ,0.8896385 ,0.9070040 ,0.8907332 ,0.84454564
 ,0.84407903 ,0.8937433 ,0.9103602 ,0.8937858 ,0.84960856
 ,0.84924411 ,0.8979929 ,0.9137611 ,0.8967297 ,0.85480217
 ,0.85426493 ,0.9025422 ,0.9170020 ,0.8996128 ,0.86038748
 ,0.85905395 ,0.9072085 ,0.9198988 ,0.9027462 ,0.86618304
 ,0.86383966 ,0.9115513 ,0.9226146 ,0.9064442 ,0.87189147
 ,0.86878464 ,0.9153553 ,0.9254626 ,0.9106255 ,0.87743613
 ,0.87369919 ,0.9187912 ,0.9285133 ,0.9149145 ,0.88292991
 ,0.87848296 ,0.9221554 ,0.9315980 ,0.9190752 ,0.88839609
 ,0.88333342 ,0.9255519 ,0.9345754 ,0.9230635 ,0.89375641
 ,0.88852444 ,0.9289291 ,0.9374606 ,0.9268703 ,0.89911378
 ,0.89422273 ,0.9325025 ,0.9403301 ,0.9305157 ,0.90464514
 ,0.90022856 ,0.9367256 ,0.9431860 ,0.9339928 ,0.91036426
 ,0.90634711 ,0.9416530 ,0.9460100 ,0.9373097 ,0.91620446
 ,0.91250898 ,0.9465501 ,0.9489720 ,0.9406272 ,0.92204347
 ,0.91826744 ,0.9504197 ,0.9522643 ,0.9440552 ,0.92778824
 ,0.92342137 ,0.9530621 ,0.9556061 ,0.9473789 ,0.93323589
 ,0.92839321 ,0.9551899 ,0.9585219 ,0.9505219 ,0.93796702
 ,0.93344816 ,0.9574912 ,0.9611978 ,0.9538301 ,0.94208767
 ,0.93844307 ,0.9600012 ,0.9642909 ,0.9574326 ,0.94613939
 ,0.94305951 ,0.9625902 ,0.9679124 ,0.9608894 ,0.95006010
 ,0.94732447 ,0.9652668 ,0.9714955 ,0.9636651 ,0.95357045
 ,0.95204731 ,0.9681020 ,0.9746738 ,0.9659957 ,0.95718333
 ,0.95754260 ,0.9712236 ,0.9775883 ,0.9688384 ,0.96124457
 ,0.96262769 ,0.9744829 ,0.9802432 ,0.9724819 ,0.96510323
 ,0.96744559 ,0.9773397 ,0.9826635 ,0.9763755 ,0.96890250
 ,0.97198841 ,0.9802640 ,0.9855485 ,0.9794530 ,0.97335289
 ,0.97581215 ,0.9841504 ,0.9886425 ,0.9820909 ,0.97807428
 ,0.98074122 ,0.9873238 ,0.9909320 ,0.9853537 ,0.98318562
 ,0.98656778 ,0.9900431 ,0.9927363 ,0.9882823 ,0.98878283
 ,0.99129231 ,0.9932818 ,0.9951413 ,0.9923500 ,0.99385907
 ,0.99566295 ,0.9970155 ,0.9968195 ,0.9963426 ,0.99728643),
 byrow=TRUE,ncol=5)

 if(is.null(PV)){
 if(!com.p.dist){
if(p==3)rem=P3
if(p==4)rem=P4
if(p==5)rem=P5
if(p==6)rem=P6
}}
est=NA
for(j in 1:p)est[j]=corfun(x[,j],y)$cor
#id=which(est==max(est))
R=order(est,decreasing=TRUE)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
bvec<-mclapply(data,corCOMmcp_sub,x,y,corfun,...)
}
if(!MC)bvec<-lapply(data,corCOMmcp_sub,x,y,corfun,...)
output=matrix(NA,nrow=pm1,ncol=9)
dimnames(output)=list(NULL,c('IV.1','IV.2','Est.1',
'Est.2','Dif','ci.low','ci.hi','p.value','adj.p.value'))
mat=matrix(NA,nrow=nboot,ncol=p)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
for(i in 1:nboot)mat[i,]=bvec[[i]]
for(j in 2:p){
k=j-1
output[k,1]=R[k]
output[k,2]=R[j]
output[k,3]=est[R[k]]
output[k,4]=est[R[j]]
bsort<-sort(mat[,R[k]]-mat[,R[j]])
output[k,5]=est[R[k]]-est[R[j]]
output[k,6]=bsort[ilow]
output[k,7]=bsort[ihi]
pv=mean(bsort<0)+.5*mean(bsort==0)
output[k,8]=2*min(c(pv,1-pv))
flag=output[k,8]>=rem[,k]
ID=which(flag==TRUE)
ic=max(ID,1)
output[k,8]=L[ic]
}
output[,9]=p.adjust(output[,8],method=method)
list(results=output)
}

corREGorder.crit<-function(p,n,corfun=wincor,iter=1000,nboot=1000,SEED=TRUE,MC=FALSE,pr=TRUE,...){
#
# Estimate null distribution of	the p-values for corREGorde
#
#
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
if(pr)print('Execution time might take several minutes')
pm1=p-1
rem=matrix(NA,iter,pm1)
p1=p+1
for(I in 1:iter){
x=rmul(n,p)
y=rnorm(n)
est=NA
for(j in 1:p)est[j]=corfun(x[,j],y)$cor
R=order(est,decreasing=TRUE)
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
bvec<-mclapply(data,corCOMmcp_sub,x,y,corfun,...)
}
if(!MC)bvec<-lapply(data,corCOMmcp_sub,x,y,corfun,...)
output=matrix(NA,nrow=pm1,ncol=8)
dimnames(output)=list(NULL,c('IV','Larges.Est','Est.2','Dif','ci.low','ci.hi','p.value','adj.p.value'))
mat=matrix(NA,nrow=nboot,ncol=p)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
for(i in 1:nboot)mat[i,]=bvec[[i]]
for(j in 2:p){
k=j-1
bsort<-sort(mat[,R[k]]-mat[,R[j]])
pv=mean(bsort<0)+.5*mean(bsort==0)
rem[I,k]=2*min(c(pv,1-pv))
}
}
rem
}

corCOM.DVvsIV.crit<-function(p,n,corfun=wincor,iter=1000,nboot=500,SEED=TRUE,MC=FALSE,...){
#
# Null p-value distribution for	corCOM.DVvsIV
#
#  p: number of	independent variables
#  n: sample size
#
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
pm1=p-1
rem=matrix(NA,iter,pm1)
p1=p+1
for(I in 1:iter){
x=rmul(n,p)
y=rnorm(n)
est=NA
for(j in 1:p)est[j]=corfun(x[,j],y)$cor
id=which(est==max(est))
R=order(est,decreasing=TRUE)
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
bvec<-mclapply(data,corCOMmcp_sub,x,y,corfun,...)
}
if(!MC)bvec<-lapply(data,corCOMmcp_sub,x,y,corfun,...)
mat=matrix(NA,nrow=nboot,ncol=p)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
for(i in 1:nboot)mat[i,]=bvec[[i]]
for(j in 2:p){
k=j-1
bsort<-sort(mat[,R[1]]-mat[,R[j]])
pv=mean(bsort<0)+.5*mean(bsort==0)
rem[I,k]=2*min(c(pv,1-pv))
}
}
rem
}

bicorM<-function(x){
a=bicovM(x)
a=cov2cor(a)
a
}

bicor<-function(x,y){
a=bicovM(cbind(x,y))
a=cov2cor(a)
list(cor=a[1,2])
}

corregci<-function(x,y,corfun=wincor,nboot=599,alpha=.05,SEED=TRUE,pr=TRUE,...){
#
#  Deals with correlations between some dependent variable y and p independent variables, x
#  Compute confidence intervals for correlation coefficients and p-values when testing the
#  hypothesis of a zero correlation,
#  Also reported are adjusted p-values based on Hochberg's method.
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=599
#
#   corfun can be any R function that returns a correlation having the form
#   the vector corfun$cor. Examples are pbcor and bicor, spear and tau.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
nrem=length(y)
estit=NA
for(j in 1:p)estit[j]=corfun(x[,j],y,...)$cor
nv=length(y)
x<-as.matrix(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,corregci.sub,x,y,corfun=corfun)
#Leverage points already removed.
# bvec is a p by nboot matrix. The first row
#                     contains the bootstrap correlations for the first IV,, the second row
#                     contains the bootstrap values for first predictor, etc.
regci<-matrix(0,p,5)
vlabs=NA
for(j in 1:p)vlabs[j]=paste("ind.var",j)
i#vlabs[1:p]=labels(x)[[2]]
dimnames(regci)<-list(vlabs,c("ci.low","ci.up","Estimate","p-value",'Adj.p.value'))
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
se<-NA
pvec<-NA
for(i in 1:p){
bsort<-sort(bvec[i,])
pvec[i]<-sum(bvec[i,]<0)/nboot  #+.5*sum(bvec[i,]==0)/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
regci[i,1]<-bsort[ilow]
regci[i,2]<-bsort[ihi]
}
regci[,3]=estit
pvec<-2*pvec
regci[,4]=pvec
regci[,5]=p.adjust(pvec,method='hoch')
num.sig=sum(regci[,5]<=alpha)
list(output=regci,n=nrem,num.sig=num.sig)
}

mcd.cor<-function(x,y){
xy=cbind(x,y)
a=MCDCOR(x=xy)[1,2]
list(cor=a)
}

corCOM.PMDPCD<-function(n,p,rho=0,delta=.3,corfun=wincor,LARGEST=TRUE,alpha=.05,
x=NULL,y=NULL,iter=500,pr=TRUE,SEED=TRUE,MC=TRUE,FUN=mean,...){
#
#  Given n and p, the number of explanatory variables,
#  determine the probability of making a decision about which independent variable has the
#  the largest correlation with the dependent variable
#  in the context of an indifference zone
#
#
#  All independent variables have a same  correlation with dependent variable indicated by the argument
#  rho,
#  except the first independent variable, which has correlation
#  rho+delta.
#  Default is rho=0 and delta= .3
#  (Cohen's suggestion for a  small medium and large  correlation
#  are .1, .3 and .5)
# This is designed for situations where the goal is to make a decision about which IV has the largest correlation.
#
#  If LARGEST=TRUE;  the function default to rho=0 for the first IV and delta for the remaining IVs
#
#  Possible alternative choices for corfun include:
#  spear
#  tau
#  pbcor
#  bicor
#  scor
#  mve.cor
#  mcd.cor
#
if(pr)print('Execution time can be high')
use.cor=FALSE
if(!is.null(x) & !is.null(y)){
xy=cbind(y,x)
R=COR.ROB(xy,method=COR.method)
rho=FUN(R[upper.tri(R)])
n=nrow(x)
p=ncol(x)
use.cor=TRUE
}

if(rho+delta>1)stop('rho+delta is greater than 1')
if(rho+delta<0-1)stop('rho+delta is less than -1')
if(SEED)set.seed(2)
if(delta<0-1 || delta>1)stop('rho + delta should be between -1 and 1')
PMD=0
PCD=0
p1=p+1
if(LARGEST){
COV=matrix(rho,p1,p1)
COV[1,2]=COV[2,1]=rho+delta
}
if(!LARGEST){
COV=matrix(rho+delta,p1,p1)
COV[1,2]=COV[2,1]=rho
}
diag(COV)=1
if(use.cor)COV=R  #Over rule using an indifference zone; use estimate of the correlation matrix
x=list()
for(i in 1:iter)x[[i]]=rmulnorm(n,p1,COV)
if(!MC)a=lapply(x,corCOM.PMDPCD.sub,corfun=corfun,LARGEST=LARGEST,...)
if(MC)a=mclapply(x,corCOM.PMDPCD.sub,corfun=corfun,LARGEST=LARGEST,...)
for(i in 1:iter){
if(a[[i]]$Conclusion=='Decide'){
PMD=PMD+1
if(a[[i]][1]==1)PCD=PCD+1
}}
PCD.CI=binom.conf(PCD,PMD,alpha=alpha,pr=FALSE)$ci
PMD.CI=binom.conf(PMD,iter,alpha=alpha,pr=FALSE)$ci
PCD=PCD/max(1,PMD)
PMD=PMD/iter
list(PMD=PMD,PMD.CI=PMD.CI, PCD=PCD, PCD.CI=PCD.CI)
}

corCOM.PMDPCD.sub<-function(x,corfun,LARGEST=LARGEST,...){
p1=ncol(x)
a=corCOM.DVvsIV(x[,2:p1],x[,1],corfun=corfun,SEED=FALSE,LARGEST=LARGEST,...)
a
}

corREG.best<-function(x,y,corfun=wincor,alpha=.05,nboot=500, neg.col=NULL,LARGEST=TRUE, SEED=TRUE,MC=FALSE,xout=FALSE,outfun=outpro,...){
#
# Can a decision be made about which IV
# has the strongest correlation with the DV
# Winsorized correlation is used by default.
#
# x is assumed to be a matrix
#
#
if(nrow(x)!=length(y))stop('x and y have different sample sizes; should be equal')
p=ncol(x)
p1=p+1
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1:p]
y=m1[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x=neg.colM(x,neg.col)
est=NA
for(j in 1:p)est[j]=corfun(x[,j],y)$cor
if(LARGEST)ID=which(est==max(est))
if(!LARGEST)ID=which(est==min(est))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
bvec<-mclapply(data,corCOMmcp_sub,x,y,corfun)
}
IB=NA
if(!MC)bvec<-lapply(data,corCOMmcp_sub,x,y,corfun)
for(i in 1:nboot){
if(LARGEST)IB[i]=which(bvec[[i]]==max(bvec[[i]]))
if(!LARGEST)IB[i]=which(bvec[[i]]==min(bvec[[i]]))
}
PC=mean(IB==ID)
PC=2*min(PC,1-PC)
list(Est.=est,p.value=PC)
}

PcorREG.best.DO<-function(x,y,neg.col=NULL,
LARGEST=TRUE,xout=FALSE,outfun=outpro,...){
#
# Can a decision be made about which IV
# has the strongest Pearson correlation with the DV
#
# x is assumed to be a matrix or data frame
#
#
if(nrow(x)!=length(y))stop('x and y have different sample sizes; should be equal')
p=ncol(x)
p1=p+1
pm1=p-1
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1:p]
y=m1[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x=neg.colM(x,neg.col)
est=NA
for(j in 1:p)est[j]=cor(x[,j],y)
if(LARGEST)ID=which(est==max(est))
if(!LARGEST)ID=which(est==min(est))
a=matrix(NA,nrow=pm1,ncol=7)
dimnames(a)=list(NULL,c('Best.IV','IV','Est.best','Est','dif','ci.low','ci.up'))
ic=0
for(j in 1:p){
if(j!=ID){
ic=ic+1
b=TWOpov(x[,c(ID,j)],y)
a[ic,]=c(ID,j,b$est.rho1,b$est.rho2,b$dif,b$ci[1],b$ci[2])
}}
chk=sign(a[,6]*a[,7])
D='No Decision'
if(sum(chk)==pm1)D=paste('Decide IV',ID,' is best')
list(output=a, Result=D)
}

corREG.DO<-function(x,y,corfun=wincor,alpha=.05,nboot=500,SEED=TRUE,neg.col=NULL,
LARGEST=TRUE,xout=FALSE,outfun=outpro,...){
#
# Can a decision be made about which IV
# has the strongest correlation with the DV
# Winsorized correlation is used by default.
#
#  Uses the max  p-value method
# x is assumed to be a matrix
#
#
if(nrow(x)!=length(y))stop('x and y have different sample sizes; should be equal')
p=ncol(x)
p1=p+1
pm1=p-1
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1:p]
y=m1[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x=neg.colM(x,neg.col)
est=NA
for(j in 1:p)est[j]=corfun(x[,j],y)$cor

if(LARGEST)ID=which(est==max(est))
if(!LARGEST)ID=which(est==min(est))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#

a=corCOMmcp(x,y,SEED=SEED,...)$results
flag=a[,1:2]==ID
pv=NA
ic=0
FLAG=rep(FALSE,nrow(a))
for(j in 1:nrow(a)){
if(flag[j,1] || flag[j,2]) FLAG[j]=TRUE
}
a=a[FLAG,1:7]
pv=max(a[,7])
dc='NO'
if(pv <=alpha)dc='YES'
list(ID=ID,Est.=est,results=a,p.value=pv,Decide=dc)
}

cor.skip.com<-function(x,y,corfun=wincor,outfun=outpro,alpha=.05,nboot=1000,SEED=TRUE,...){
#
#  Regression,  two explanatory variables
#  Compare the skipped correlation with the dependent variable.
#  That is, deal with the overlapping case.
#
#  Limited to 0.95 confidence interval for the difference between the
#  measures of association.
#
#  Method: Bias corrected accelerated bootstrap
#
library(bcaboot)
if(SEED)set.seed(2)
m=cbind(x,y)
m=elimna(m)
id=outfun(m[,c(1,3)],plotit=FALSE)$keep
e1=corfun(m[id,1],m[id,3])$cor
id=outfun(m[,c(1,3)],plotit=FALSE)$keep
e2=corfun(m[id,2],m[id,3])$cor
dif=e1-e2
a=bcajack(m,nboot,skip.o.lap,corfun=corfun,outfun=outfun,alpha=alpha/2,verbose=FALSE,...)
ci=c(a$lims[1,1],a$lims[3,1])
list(Est1=e1,Est2=e2,dif=dif,ci=ci)
}

corskip.comPV<-function(x,y,corfun=wincor,outfun=outpro,alpha=.05,nboot=1000,SEED=TRUE,...){
#
#  Regression,  two explanatory variables
#  Compare the skipped correlation with the dependent variable.
#  That is, deal with the overlapping case.
#
#  Returns a p-value plus confidence interval
#
#  Method: Bias corrected accelerated bootstrap
#
library(bcaboot)
if(SEED)set.seed(2)
ALPHA=c(seq(.001,.1,.001),seq(.011,.05,.01),seq(.11,.99,.01))
al=length(ALPHA)
il=0
m=cbind(x,y)
m=elimna(m)
if(ncol(m)!=3)stop('Must have exactly two independent variables')
est1=corfun(m[,1],m[,3],...)$cor
est2=corfun(m[,2],m[,3],...)$cor
dif=est1-est2
a=bcajack(m,nboot,skip.o.lap,corfun=corfun,outfun=outfun,alpha=ALPHA/2,verbose=FALSE,...)
iu=nrow(a$lims)+1
for(j in 1:al){
il=il+1
iu=iu-1
pv=ALPHA[j]
ci=c(a$lims[il,1],a$lims[iu,1])
if(a$lims[il,1]>0. || a$lims[iu,1]<0.)break
}
A=cor.skip.com(x,y,corfun=corfun,outfun=outpro,alpha=.05,nboot=nboot,SEED=SEED,...)
ci=c(a$lims[1,1],a$lims[3,1])
list(n=nrow(m),Est1=est1,Est2=est2,difference=dif,ci.low=A[1],ci.upper=A[2],p.value=pv)
}

rmdif.scores<-function(x){
#
# Compute all pairwise difference scores
#

  if(!is.matrix(x) & !is.data.frame(x))stop('x should be matrix or data frame')
  x=elimna(x)
  n=nrow(x)
  J=ncol(x)
  ALL=(J^2-J)/2
  M=matrix(NA,nrow=n,ncol=ALL)
  ic=0
  for(j in 1:J){
    for(k in 1:J){
      if(j<k){
        ic=ic+1
        M[,ic]=x[,j]-x[,k]
      }}}
  M
}

smeancr.cord.oph<-function(m,nullv=rep(0,ncol(m)),cop=3,MM=FALSE,SEED=TRUE,
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
list(p.value=sig.level,boot.vals=val,center=est)
}

smeancr.cord<-function(m,nullv=rep(0,ncol(m)),cop=3,MM=FALSE,SEED=TRUE,
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

meancr.cord.oph<-function(m,nullv=rep(0,ncol(m)),cop=3,MM=FALSE,SEED=TRUE,tr=0,
nboot=500,plotit=TRUE,MC=FALSE,xlab="VAR 1",ylab="VAR 2",STAND=TRUE){
#
# m is an n by p matrix
#
# Test hypothesis that the means
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
#est=smean(m,MC=MC,cop=cop,STAND=STAND)
est=apply(m,2,mean,tr=tr)
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
#val[j,]<-smean(mm,MC=MC,cop=cop,STAND=STAND)
val[j,]=apply(mm,2,mean)
}
if(!MC)temp<-pdis(rbind(val,nullv),center=est)
if(MC)temp<-pdisMC(rbind(val,nullv),center=est)
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
list(p.value=sig.level,boot.vals=val,center=est)
}

part.cor<-function(x,y,z,corfun=wincor,regfun=MMreg,plotit=FALSE,xout=FALSE,GEN=TRUE,BOOT=TRUE,SEED=TRUE,nboot=599,
XOUT.blp=TRUE,plot.out=FALSE,
outfun=outpro,plotfun=plot,xlab='Res 1',ylab='Res 2',...){
#
# Robust partial correlation.
# Uses the correlation between the residuals of  x with z and y with  z
#
# Default is a Winsorized  correlation between x and y controlling for z
# XOUT.blp=TRUE means that if any bad leverage points are detected, they are removed
#
#
# If XOUT.blp=FALSE
# and
# xout=TRUE remove leverage points. If
# GEN =TRUE, remove only bad leverage when dealing with the association between
#. x and z as well as y and z. In contrast,
#  XOUT.blp=TRUE  permanently removes bad leverage points associated the
#  regression line for x and y, where x is the independent variable.
#
# if z contains a dummy variable, can ignore the corresponding col when removing outliers
#.  Example
# part.cor(x,y,z,GEN=FALSE,outfun=out.dummy,id=2,xout=TRUE)
#
# Examples:
# part.cor(x,y,z,regfun=MMreg,corfun=wincor)
# part.cor(x,y,z,regfun=MMreg,corfun=cor.test,method='kendall')
# part.cor(x,y,z,regfun=MMreg,corfun=cor.test,method='spear')
# part.cor(x,y,z,regfun=MMreg,corfun=scor) #skpped correlation,
#
#
xyz=elimna(cbind(x,y,z))
p3=ncol(xyz)
p1=p3-1
x=xyz[,1]
y=xyz[,2]
z=xyz[,3:p3]
z=as.matrix(z)
if(XOUT.blp){
id=outblp(x,y)$keep
x=x[id]
y=y[id]
z=z[id,]
xout=FALSE
}
if(xout){
if(GEN){
e1=reg.reglev(z,x,regfun=regfun)$coef
e2=reg.reglev(z,y,regfun=regfun)$coef
}
else{
e1=regfun(z,x,xout=xout,outfun=outfun,...)$coef
e2=regfun(z,y,xout=xout,outfun=outfun,...)$coef
}
}
if(!xout){
e1=regfun(z,x)$coef
e2=regfun(z,y)$coef
}
z=as.matrix(z)
res1=x-z%*%e1[2:p1]-e1[1]
res2=y-z%*%e2[2:p1]-e2[1]
if(plotit){
if(plot.out){
id=outpro(res1)$keep
res1=res1[id]
res2=res2[id]
}
if(identical(plotfun,plot))plot(res1,res2,xlab=xlab,ylab=ylab)
else plotfun(res1,res2,xlab=xlab,ylab=ylab,pr=FALSE)
}
if(BOOT)est=corb(res1,res2,corfun,SEED=SEED,nboot=nboot)
else
est=corfun(res1,res2)
est
}

corblp.EP<-function(x,y,regfun=tsreg,varfun=pbvar,plotit=FALSE,...){
#
# Correlation based on a robust regression estimator with bad
# leverage points removes
#
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
id=reglev.gen(x,y,regfun=regfun,plotit=plotit)$keep
e=reg.pred(x[id,],y[id],regfun=regfun)
top=varfun(e)
bot=varfun(y)
rsq=top/bot
rsq=min(rsq,1)
rest=NA
if(p==1){
est=regfun(x[id],y[id])$coef[2]
rest=sign(est)*sqrt(rsq)
}
list(cor=rest,Rsq=rsq)
}

corblp.ci<-function(x,y,regfun=tsreg,varfun=pbvar,nboot=100,alpha=.05,outfun=outpro.depth,SEED=TRUE,
plotit=FALSE,...){
#
# Correlation, basically a robust version of explanatory power,
#  based on a robust regression estimator with bad
# leverage points removes
#
if(SEED)set.seed(2)
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
if(p!=1)stop('Only a single independent variable is allowed')
x=xy[,1]
y=xy[,2]
id=reglev.gen(x,y,regfun=regfun,plotit=plotit,outfun=outpro.depth)$keep
X=x[id]
Y=y[id]
v=NA
n=length(Y)
bot=varfun(Y)
for(i in 1:nboot){
id=sample(n,replace=TRUE)
e=reg.pred(X[id],Y[id],regfun=regfun)
top=varfun(e)
rsq=top/bot
rsq=min(rsq,1)
est=regfun(X[id],Y[id])$coef[2]
rest=sign(est)*sqrt(rsq)
v[i]=rest
}
se=sd(v)
est=corblp.EP(x,y,regfun=regfun,varfun=varfun)
test=est$cor/se
sig<-2*(1-pnorm(abs(test)))
crit=qnorm(1-alpha/2)
ci=est$cor-crit*se
ci=max(ci,-1)
ci[2]=est$cor+crit*se
ci[2]=min(ci[2],1)
list(cor=est$cor,test=test,p.value=sig,ci=ci)
}

corblp.bca.C<-function(x,y,regfun=tsreg,varfun=pbvar,nboot=1000,alpha=.05,outfun=outpro.depth,SEED=TRUE,
plotit=FALSE,...){
#
# Correlation based on a robust regression estimator with bad
# leverage points removes

library(bcaboot)
if(SEED)set.seed(2)
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
if(p!=1)stop('Only a single independent variable is allowed')
x=xy[,1]
y=xy[,2]
n=length(y)
est=corblp(x,y,regfun=regfun,varfun=varfun)$cor
a=bcajack2(xy,1000,corblp.sub,alpha=alpha/2,regfun=regfun,varfun=varfun)
ci=c(a$lims[1,1],a$lims[3,1])
list(cor=est,ci=ci)
}

corblppb<-function(x,y,regfun=tsreg,varfun=pbvar,nboot=1000,alpha=.05,outfun=outpro.depth,SEED=TRUE,
plotit=FALSE,...){
#
# Correlation based on a robust regression estimator with bad
# leverage points removes
#
if(SEED)set.seed(2)
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
if(p!=1)stop('Only a single independent variable is allowed')
x=xy[,1]
y=xy[,2]
v=NA
n=length(y)
for(i in 1:nboot){
id=sample(n,replace=TRUE)
v[i]=corblp(x[id],y[id],regfun=regfun,varfun=varfun)$cor
}
P=mean(v<0)
pv=2*min(P,1-P)
sv=sort(v)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci<-sv[ilow]
ci[2]<-sv[ihi]
est=corblp(x,y,regfun=regfun,varfun=varfun)
list(cor.est=est$cor,p.value=pv,ci=ci)
}

qcorp1.ci<-function(x,y,q=.5,alpha=.05,nboot=599,SEED=TRUE, xout=TRUE,
method='PRO',regfun=Qreg){
#
# Confidence interval for a quantile regression measure of association
#
#
if(SEED)set.seed(2)
 xy=elimna(cbind(x,y))
p1=ncol(xy)
if(p1>2)stop('Current version is for  a single independent variable only')
#, use the R function  corblp.ci')
x=xy[,1]
y=xy[,2]
 if(xout){
x<-as.matrix(x)
flag<-out.methods(x,y,plotit=FALSE,method=method,regfun=regfun)$keep
x<-x[flag,]
y<-y[flag]
n.keep=length(y)
}
p=p1-1
x=xy[,1:p]
x=as.matrix(x)
y=xy[,p1]
 v=NA
 n=nrow(xy)
 if(n<=40)xout=FALSE # to avoid computational issues when n is small.
 for(i in 1:nboot){
 id=sample(n,replace=TRUE)
v[i]=qcorp1(xy[id,1:p],xy[id,p1],q=q,xout=xout)$cor
 }
 mv=mean(v)
 v=sort(v)
 ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
 ci=v[ilow]
 ci[2]=v[ihi]
 e=qcorp1(x,y,q=q,xout=xout)$cor
 pv=mean(v<0)+.5*mean(v==0)
 pv=2*min(pv,1-pv)
 list(Est.=e,ci=ci,p.value=pv)
 }

qcor.ci<-function(x,y,q=.5,alpha=.05,nboot=599,SEED=TRUE, xout=TRUE,
method='PRO',regfun=Qreg){
#
# Confidence interval for a quantile regression measure of association derived by Li et al.
#
#
if(SEED)set.seed(2)
 xy=elimna(cbind(x,y))
p1=ncol(xy)
if(p1>2)stop('Current version is for  a single independent variable only')
#, use the R function  corblp.ci')
x=xy[,1]
y=xy[,2]
 if(xout){
x<-as.matrix(x)
flag<-out.methods(x,y,plotit=FALSE,method=method,regfun=regfun)$keep
x<-x[flag,]
y<-y[flag]
n.keep=length(y)
}
p=p1-1
x=xy[,1:p]
x=as.matrix(x)
y=xy[,p1]
 v=NA
 n=nrow(xy)
 if(n<=40)xout=FALSE # to avoid computational issues when n is small.
 for(i in 1:nboot){
 id=sample(n,replace=TRUE)
v[i]=qcor(xy[id,1:p],xy[id,p1],q=q,xout=xout)$cor
 }
 mv=mean(v)
 v=sort(v)
 ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
 ci=v[ilow]
 ci[2]=v[ihi]
  pv=mean(v<0)+.5*mean(v==0)
 pv=2*min(pv,1-pv)
 e=qcor(x,y,q=q)$cor    # Already did or did not remove leverage points
 list(Est.=e,ci=ci,p.value=pv)
 }

qcor.ep.ci<-function(x,y,q=.5,alpha=.05,nboot=599,SEED=TRUE, xout=TRUE,
method='PRO',regfun=Qreg){
#
# Confidence interval for a quantile regression measure of association
#
#
if(SEED)set.seed(2)
 xy=elimna(cbind(x,y))
p1=ncol(xy)
if(p1>2)stop('Current version is for  a single independent variable only')
#, use the R function  corblp.ci')
 if(xout){
x<-as.matrix(x)
flag<-out.methods(x,y,plotit=FALSE,method=method,regfun=regfun)$keep
x<-x[flag,]
y<-y[flag]
n.keep=length(y)
}
p=p1-1
x=xy[,1:p]
x=as.matrix(x)
y=xy[,p1]
 v=NA
 n=nrow(xy)
 if(n<=40)xout=FALSE # to avoid computational issues when n is small.
 for(i in 1:nboot){
 id=sample(n,replace=TRUE)
v[i]=qcor.ep(xy[id,1:p],xy[id,p1],q=q,xout=xout)$cor
 }
 mv=mean(v)
 v=sort(v)
 ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
 ci=v[ilow]
 ci[2]=v[ihi]
 e=qcor.ep(x,y,q=q,xout=xout)$cor
  pv=mean(v<0)+.5*mean(v==0)
 pv=2*min(pv,1-pv)
 list(Est.=e,ci=ci,p.value=pv)
 }

qcor.R<-function(x,y,q=c(.25,.5,.75),alpha=.05,nboot=1000,SEED=TRUE, xout=TRUE,
method='PRO',regfun=Qreg){
#
# Correlation based on the quantile regression estimator.
# This version is based on the ratio of the  loss function for the full model versus the null case
# See Wilcox (2022, section 11.9).
#
nq=length(q)
res=matrix(NA,nq,5)
for(j in 1:nq){
a=qcorp1.ci(x,y,q=q[j],alpha=alpha,nboot=nboot,SEED=SEED,xout=xout)
res[j,]=c(q[j],a$Est.,a$ci[1],a$ci[2],a$p.value)
}
dimnames(res)=list(NULL,c('q','Est','ci.low','ci.up','p-value'))
res
}

qcor.EP<-function(x,y,q=c(.25,.5,.75),alpha=.05,nboot=1000,SEED=TRUE, xout=TRUE,
method='PRO',regfun=Qreg){
#
# Correlation based on the quantile regression estimator.
# This version is based on explanatory power.
# See Wilcox (2022, section 11.9).
#
nq=length(q)
res=matrix(NA,nq,5)
for(j in 1:nq){
a=qcor.ep.ci(x,y,q=q[j],alpha=alpha,nboot=nboot,SEED=SEED,xout=xout)
res[j,]=c(q[j],a$Est.,a$ci[1],a$ci[2],a$p.value)
}
dimnames(res)=list(NULL,c('q','Est','ci.low','ci.up','p-value'))
res
}

qcor.ep<-function(x,y,qest=hd,q=.5,xout=FALSE,method='PRO',regfun=MMreg,
plotit=FALSE,...){
#
# Compute a measure of the strength of the association in terms of explanatory power and
# based on the quantile regression line
#
X=cbind(x,y)
X=elimna(X)
x<-as.matrix(x)
p=ncol(x)
x=X[,1:p]
p1=p+1
y=X[,p1]
if(xout){
x<-as.matrix(x)
flag<-out.methods(x,y,plotit=plotit,method=method,regfun=regfun)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
X=cbind(x,y)
}
est=qreg(x,y,q=q)$coef
pred=reg.pred(x,y,regfun=Qreg,q=q)
top=pbvar(pred)
bot=pbvar(y)
EP=top/bot
if(p==1)ce=sign(est[2])*sqrt(EP)
list(cor=ce,Explanatory.power=EP)
}

