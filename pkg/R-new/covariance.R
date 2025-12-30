# covariance.R
#
# Covariance and Scatter Matrix Estimation Functions
# Part of WRS package refactoring
#
# This module contains functions for:
# - Robust covariance estimation (OGK, MVE, MCD, S-estimators, MBA)
# - Winsorized and trimmed covariance matrices
# - Skipped covariance (outlier removal)
# - Bootstrap covariance estimation
# - Distance covariance and ROC-based covariance
# - Mixed design covariance (between-within)
# - Median-based covariance methods
# - Utility functions for covariance matrices
#
# Extracted: 2025-12-30
# Number of functions: 43


bicov<-function(x,y){
#
# compute biweight midcovariance of x and y
#
mx<-median(x)
my<-median(y)
ux<-abs((x-mx)/(9*qnorm(.75)*mad(x)))
uy<-abs((y-my)/(9*qnorm(.75)*mad(y)))
aval<-ifelse(ux<=1,1,0)
bval<-ifelse(uy<=1,1,0)
top<-sum(aval*(x-mx)*(1-ux^2)^2*bval*(y-my)*(1-uy^2)^2)
top<-length(x)*top
botx<-sum(aval*(1-ux^2)*(1-5*ux^2))
boty<-sum(bval*(1-uy^2)*(1-5*uy^2))
bi<-top/(botx*boty)
bi
}

mscov<-function(m,STAND=TRUE){
#
# m is an n by p matrix
#
# Compute a skipped covariance matrix
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
# Eliminate any outliers and compute covariances
#  using remaining data.
#
m<-elimna(m)
temp<-outpro(m,plotit=FALSE,STAND=STAND)$keep
mcor<-var(m[temp,])
mcor
}

cov2med<-function(x,y=NA,q=.5){
#
# Estimate the covariance between two dependent
# order statistics
# By default, q=.5 meaning that an estimate of
# of covariance is made when a single order statistic
# is used to estimate the median.
# y=NA, function returns squared standard error.
#
if(is.na(y[1]))val<-qse(x,q=q,op=3)^2
if(!is.na(y[1])){
if(sum((x-y)^2)==0)val<-qse(x,q=q,op=3)^2
if(sum((x-y)^2)>0){
n<-length(x)
m<-floor(q*n+.5)
yord<-sort(y)
flag<-(y<=yord[m])
xord<-sort(x)
xq<-xord[m]
yord<-sort(y)
yq<-yord[m]
flag1<-(x<=xq)
flag2<-(y<=yq)
A<-mean(flag1*flag2)
flag1<-(x<=xq)
flag2<-(y>yq)
B<-mean(flag1*flag2)
flag1<-(x>xq)
flag2<-(y<=yq)
C1<-mean(flag1*flag2)
flag1<-(x>xq)
flag2<-(y>yq)
D1<-mean(flag1*flag2)
fx<-akerd(x,pts=xq,plotit=FALSE,pyhat=TRUE)
fy<-akerd(y,pts=yq,plotit=FALSE,pyhat=TRUE)
v1<-(q-1)^2*A
v2<-(q-1)*q*B
v3<-(q-1)*q*C1
v4<-q*q*D1
val<-((v1+v2+v3+v4)/(fx*fy))/n
}}
val
}

covmmed<-function(x,p=length(x),grp=c(1:p),q=.5){
#
#  Estimate the covariance matrix for the sample medians
#  based on a SINGLE order statistic, using
#  the data in the R variable x.
# (x[[1]] contains the data for group 1, x[[2]] the data for group 2, etc.)
#  The function returns a p by p matrix of covariances, the diagonal
#  elements being equal to the squared standard error of the sample
#  trimmed means, where p is the number of groups to be included.
#  By default, all the groups in x are used, but a subset of
#  the groups can be used via grp.  For example, if
#  the goal is to estimate the covariances between the medians
#   for groups 1, 2, and 5, use the command grp<-c(1,2,5)
#  before calling this function.
#
#  Missing values (values stored as NA) are not allowed.
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("The data are not stored in a matrix or list mode.")
p<-length(grp)
pm1<-p-1
for (i in 1:pm1){
ip<-i+1
if(length(x[[grp[ip]]])!=length(x[[grp[i]]]))stop("The number of observations in each group must be equal")
}
n<-length(x[[grp[1]]])
covest<-matrix(0,p,p)
for(j in 1:p){
for(k in 1:p){
if(j==k)covest[j,j]<-cov2med(x[[grp[j]]],q=q)
if(j<k){
covest[j,k]<-cov2med(x[[grp[j]]],x[[grp[k]]],q=q)
covest[k,j]<-covest[j,k]
}}}
covest
}

bootcov<-function(x,est=median,nboot=100,pr=TRUE,SEED=FALSE,...){
#
# For multivariate data, determine the squared standard errors
# and covariances when using the estimator
# est.
#
# SEED=TRUE, sets the seed of the random number generator.
#
if(SEED)set.seed(2)
if(is.list(x))x<-matl(x)
x<-elimna(x)
bvec<-matrix(NA,ncol=ncol(x),nrow=nboot)
if(pr)print("Taking Bootstrap Samples. Please wait.")
for(i in 1:nboot){
data<-sample(nrow(x),size=nrow(x),replace=TRUE)
bvec[i,]<-apply(x[data,],2,FUN=est,...)
}
covmat<-var(bvec)
covmat
}

covroc<-function(x){
#
# compute Rocke's TBS covariance matrix
#
 library(robust)
temp<-covRob(x,estim="M")
val<-temp[2]$cov
val
}

covmtrim<-function(x,tr=.2,p=length(x),grp=c(1:p)){
#
#  Estimate the covariance matrix for the sample trimmed means corresponding
#  to the data in the R variable x,
#  which is assumed to be stored in list mode or a matrix.
# (x[[1]] contains the data for group 1, x[[2]] the data for group 2, etc.)
#  The function returns a p by p matrix of covariances, the diagonal
#  elements being equal to the squared standard error of the sample
#  trimmed means, where p is the number of groups to be included.
#  By default, all the groups in x are used, but a subset of
#  the groups can be used via grp.  For example, if
#  the goal is to estimate the covariances between the sample trimmed
#  means for groups 1, 2, and 5, use the command grp<-c(1,2,5)
#  before calling this function.
#
#  The default amount of trimming is 20%
#
#  Missing values (values stored as NA) are not allowed.
#
#  This function uses winvar from chapter 2.
#
if(is.list(x))x=matl(x)
x=elimna(x)
x=listm(x)
if(!is.list(x))stop("The data are not stored in list mode or a matrix.")
p<-length(grp)
pm1<-p-1
for (i in 1:pm1){
ip<-i+1
if(length(x[[grp[ip]]])!=length(x[[grp[i]]]))stop("The number of observations in each group must be equal")
}
n<-length(x[[grp[1]]])
h<-length(x[[grp[1]]])-2*floor(tr*length(x[[grp[1]]]))
covest<-matrix(0,p,p)
covest[1,1]<-(n-1)*winvar(x[[grp[1]]],tr)/(h*(h-1))
for (j in 2:p){
jk<-j-1
covest[j,j]<-(n-1)*winvar(x[[grp[j]]],tr)/(h*(h-1))
for (k in 1:jk){
covest[j,k]<-(n-1)*wincor(x[[grp[j]]],x[[grp[k]]],tr)$cov/(h*(h-1))
covest[k,j]<-covest[j,k]
}
}
covmtrim<-covest
covmtrim
}

bwwcovm<-function(J,K,L,x,tr=.2){
#
# compute covariance matrix for a between by within by within design
#
p=J*K*L
idep=K*L
mat=matrix(0,nrow=p,ncol=p)
id=c(1:idep)
for(j in 1:J){
mat[id,id]=covmtrim(x[id],tr=tr)
id=id+idep
}
mat
}

bbwcovm<-function(J,K,L,x,tr=.2){
#
# compute covariance matrix for a between by between by within design
#
p=J*K*L
idep=L
mat=matrix(0,nrow=p,ncol=p)
id=c(1:idep)
for(j in 1:J){
for(k in 1:K){
mat[id,id]=covmtrim(x[id],tr=tr)
id=id+idep
}}
mat
}

bicovm<-function(x){
#
# compute a biweight midcovariance matrix for the vectors of
# observations in x, where x is assumed to have list mode, or
# x is an n by p matrix
#
if(is.matrix(x)){
mcov<-matrix(0,ncol(x),ncol(x))
mcor<-matrix(0,ncol(x),ncol(x))
for (i in 1:ncol(x)){
for (j in 1:ncol(x))mcov[i,j]<-bicov(x[,i],x[,j])
}
}
if(is.list(x)){
mcov<-matrix(0,length(x),length(x))
mcor<-matrix(0,length(x),length(x))
for (i in 1:length(x)){
for (j in 1:length(x))mcov[i,j]<-bicov(x[[i]],x[[j]])
}
}
for (i in 1:ncol(mcov)){
for (j in 1:ncol(mcov))mcor[i,j]<-mcov[i,j]/sqrt(mcov[i,i]*mcov[j,j])
}
list(mcov=mcov,mcor=mcor)
}

bicovM<-function(x){
M=bicovm(x)$mcov
M
}

covmve<-function(x){
library(MASS)
oldSeed <- .Random.seed
val<-cov.mve(x)
assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
list(center=val$center,cov=val$cov)
}

mvecov<-function(x){
library(MASS)
val<-cov.mve(x)
val$cov
}

covmcd<-function(x,nsamp="sample"){
#
# nsamp="best" is the default used by R,
# meaning that  the number of samples is chosen so that
# exhaustive enumeration is done up to 5000 samples
# nsamp="sample" the number of samples
#  is min(5*p, 3000)
#
library(MASS)
oldSeed <- .Random.seed
val<-cov.mcd(x,nsamp=nsamp)
assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
list(center=val$center,cov=val$cov)
}

mcdcov<-function(x,nsamp="sample"){
#
# nsamp="best" is the default used by R,
# meaning that  the number of samples is chosen so that
# exhaustive enumeration is done up to 5000 samples
# nsamp="sample" the number of samples
#  is min(5*p, 3000)
#
#library(lqs)
library(MASS)
oldSeed <- .Random.seed
val<-cov.mcd(x,nsamp=nsamp)
   assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
val$cov
}

gkcov<-function(x,y,gk.sigmamu=taulc,...){
#
# Compute robust covariance using the Gnanadesikan-Kettenring
# estimator.
# (cf. Marrona & Zomar, 2002, Technometrics
#
val<-.25*(gk.sigmamu(x+y,...)-gk.sigmamu(x-y,...))
val
}

covogk<-function(x,sigmamu=taulc,v=gkcov,n.iter=5,beta=.9,...){
#
# Compute robust (weighted) covariance matrix in Maronna and Zamar
# (2002, Technometrics, eq. 7).
#
# x is an n by p matrix
# n.iter number of iterations. 1 seems to be best
# sigmamu is any user supplied function having the form
#   sigmamu(x,mu.too=F) and which computes a robust measure of
#   of dispersion if mu.too=F. If mu.too=T, it returns
#   a robust measure of location as well.
# v is any robust covariance
#
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)  # remove any rows with missing data
temp<-ogk.pairwise(x,sigmamu=sigmamu,v=v,n.iter=n.iter,beta=beta,...)$wcovmat
temp
}

outogk<-function(x,sigmamu=taulc,v=gkcov,op=TRUE,SEED=FALSE,
beta=max(c(.95,min(c(.99,1/nrow(x)+.94)))),n.iter=1,plotit=TRUE,...){
#
# Use the ogk estimator to
# determine which points are outliers
#
#  op=T uses robust Mahalanobis distance based on
#  the OGK estimator with  beta adjusted so that
#  the outside rate per observation is approximately .05
#  under normality.
#  op=F returns the outliers based on the distances used
#  by the OGK estimator
#  (Currently, op=T seems best for detecting outliers.)
#
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

outcov<-function(x,y=NA,outfun=outogk,plotit=FALSE){
#
# Remove outliers and compute covariances
#
if(!is.na(y[1]))x<-cbind(x,y)
keep<-outfun(x,plotit=plotit)$keep
val<-var(x[keep,])
if(ncol(val)==2)val<-val[1,2]
list(cov=val)
}

covout<-function(x,y=NA,outfun=outogk,plotit=FALSE){
#
# Remove outliers and compute covariances
#
if(!is.na(y[1]))x<-cbind(x,y)
keep<-outfun(x,plotit=plotit)$keep
val<-var(x[keep,])
if(ncol(val)==2)val<-val[1,2]
val
}

skipogk<-function(x,y=NA,plotit=FALSE){
#
# Remove outliers and compute correlations
#
if(!is.na(y[1]))x<-cbind(x,y)
x<-elimna(x)
n<-nrow(x)
keep<-outogk(x,plotit=plotit)$keep
val<-cor(x[keep,])
p.value<-NA
test<-NA
crit.05<-15.49/n+2.68
vat<-val
diag(vat)<-0
test<-abs(vat*sqrt((n-2)/(1-vat^2)))
diag(test)<-NA
if(ncol(val)==2){
p.value<-c("Greater than .1")
val<-val[1,2]
test<-abs(val*sqrt((n-2)/(1-val^2)))
crit<-4.8/n+2.72
if(test>=crit)p.value<-c("Less than .1")
crit<-15.49/n+2.68
if(test>=crit)p.value<-c("Less than .05")
crit<-14.22/n+3.26
if(test>=crit)p.value<-c("Less than .025")
crit<-24.83/n+3.74
if(test>=crit)p.value<-c("Less than .01")
}
list(cor=val,test.stat=test,p.value=p.value,crit.05=crit.05)
}

covmba <- function(x, csteps = 5)
{  # gets the MBA estimator
        zx <- x
        x <- as.matrix(x)
	p <- dim(x)[2]
	##get the DGK estimator
	covs <- var(x)
	mns <- apply(x, 2, mean)	## concentrate
	for(i in 1:csteps) {
              	md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
             if(p > 1){
		mns <- apply(x[md2 <= medd2,  ], 2,
			mean)
		covs <- var(x[md2 <= medd2,  ])
             }
             if(p == 1){
		mns <- mean(x[md2 <= medd2])
		covs <- var(x[md2 <= medd2])
             }
	}
	covb <- covs
	mnb <- mns	##get the square root of det(covb)
	critb <- prod(diag(chol(covb)))
	##get the resistant estimator
	covv <- diag(p)
	med <- apply(x, 2, median)
	md2 <- mahalanobis(x, center = med, covv)
	medd2 <- median(md2)	## get the start
        if(p > 1){
	mns <- apply(x[md2 <= medd2,  ], 2, mean)
	covs <- var(x[md2 <= medd2,  ])
        }
        if(p == 1){
	mns <- mean(zx[md2 <= medd2])
	covs <- var(zx[md2 <= medd2])
        }
        ## concentrate
	for(i in 1:csteps) {
		md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
              if(p > 1){
		mns <- apply(x[md2 <= medd2,  ], 2,
			mean)
		covs <- var(x[md2 <= medd2,  ])
              }
              if(p == 1){
	        mns <- mean(zx[md2 <= medd2])
	        covs <- var(zx[md2 <= medd2])
               }
	}
	crit <- prod(diag(chol(covs)))
	if(crit < critb) {
		critb <- crit
		covb <- covs
		mnb <- mns
	}
##scale for better performance at MVN
	rd2 <- mahalanobis(x, mnb, covb)
	const <- median(rd2)/(qchisq(0.5, p))
	covb <- const * covb
	list(center = mnb, cov = covb)
}

skipcov<-function(m,cop=6,MM=FALSE,op=1,mgv.op=0,outpro.cop=3,STAND=TRUE){
#
# m is an n by p matrix
#
# Compute skipped covariance matrix
#
# op=1:
# Eliminate outliers using a projection method
# That is, first determine center of data using:
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
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
#
# op=2 use mgv (function outmgv) method to eliminate outliers
#
# Eliminate any outliers and compute means
#  using remaining data.
# mgv.op=0, mgv uses all pairwise distances to determine center of the data
# mgv.op=1 uses MVE
# mgv.op=2 uses MCD
#
temp<-NA
m<-elimna(m)
m<-as.matrix(m)
if(op==2)temp<-outmgv(m,plotit=FALSE,op=mgv.op)$keep
if(op==1)temp<-outpro(m,plotit=FALSE,MM=MM,cop=outpro.cop,STAND=STAND,pr=FALSE)$keep
val<-var(m[temp,])
val
}

mgvcov<-function(m,op=1,cov.fun=rmba,plotit=FALSE){
#
# m is an n by p matrix
#
# Compute skipped covariance matrix
# using the MGV method
#
# Eliminate any outliers and compute covariance matrix
#  using remaining data.
# op=0, mgv uses all pairwise distances to determine center of the data
# op=1 uses the function indicated by the argument
# cov.fun to determine center of the data cloud.
# default is Olive's median ball
#
#
temp<-NA
m<-elimna(m)
temp<-outmgv(m,plotit=plotit,op=op,cov.fun=cov.fun)$keep
val<-var(m[temp,])
val
}

rmba<-function(x, csteps = 5,na.rm=TRUE,plotit=FALSE)
{
# computes the reweighted MBA estimator
# Code supplied by David Olive
#
#       x is assumed to be a matrix
#
#  plotit=FALSE is used to avoid problems when this function is called
#  by other function in WRS
#
x=as.matrix(x)
if(na.rm)x=elimna(x)
	p <- dim(x)[2]
	n <- dim(x)[1]	##get the DGK estimator
	covs <- var(x)
	mns <- apply(x, 2, mean)	## concentrate
	for(i in 1:csteps) {
		md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		covs <- var(x[md2 <= medd2,  ])
	}
	covb <- covs
	mnb <- mns	##get the square root of det(covb)
	critb <- prod(diag(chol(covb)))	##get the resistant estimator
	covv <- diag(p)
	med <- apply(x, 2, median)
	md2 <- mahalanobis(x, center = med, covv)
	medd2 <- median(md2)	## get the start
	mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
	covs <- var(x[md2 <= medd2,  ])	## concentrate
	for(i in 1:csteps) {
		md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
#		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		covs <- var(x[md2 <= medd2,  ])
	}
	crit <- prod(diag(chol(covs)))
	if(crit < critb) {
		critb <- crit
		covb <- covs
		mnb <- mns
	}
##scale for better performance at MVN
	rd2 <- mahalanobis(x, mnb, covb)
	const <- median(rd2)/(qchisq(0.5, p))
	covb <- const * covb
	##reweight the above MBA estimator (mnb,covb) for efficiency
	rd2 <- mahalanobis(x, mnb, covb)
	up <- qchisq(0.975, p)
	rmnb <- apply(as.matrix(x[rd2 <= up,  ]), 2, mean)
	rcovb <- var(x[rd2 <= up,  ])
	rd2 <- mahalanobis(x, rmnb, rcovb)
	const <- median(rd2)/(qchisq(0.5, p))
	rcovb <- const * rcovb	## reweight again
	rd2 <- mahalanobis(x, rmnb, rcovb)
	up <- qchisq(0.975, p)
	rmnb <- apply(as.matrix(x[rd2 <= up,  ]), 2, mean)
	rcovb <- var(x[rd2 <= up,  ])
	rd2 <- mahalanobis(x, rmnb, rcovb)
	const <- median(rd2)/(qchisq(0.5, p))
	rcovb <- const * rcovb
cor.b=NULL
temp=outer(sqrt(diag(rcovb)),sqrt(diag(rcovb)),'*')
if(min(diag(rcovb)>0))cor.b=rcovb/temp
	list(center = rmnb, cov = rcovb, cor=cor.b)
}

tbscov <- function(x,eps=1e-3,maxiter=20,r=.45,alpha=.05){
#        Rocke's contrained s-estimator
#   returns covariance matrix only. For both locatiion and scatter, use tbs
#
#      r=.45 is the breakdown point
#      alpha=.05 is the asymptotic rejection probability.
#
if(!is.matrix(x))stop("x should be a matrix with two or more columns")
x<-elimna(x)
library(MASS)
temp<-cov.mve(x)
t1<-temp$center
s<-temp$cov
    n <- nrow(x)
    p <- ncol(x)
if(p==1)stop("x should be a matrix with two or more columns")
c1M<-cgen.bt(n,p,r,alpha,asymp=FALSE)
c1<-c1M$c1
if(c1==0)c1<-.001 #Otherwise get division by zero
M<-c1M$M
    b0 <- erho.bt(p,c1,M)
    crit <- 100
    iter <- 1
    w1d <- rep(1,n)
    w2d <- w1d
    while ((crit > eps)&(iter <= maxiter))
    {
        t.old <- t1
        s.old <- s
        wt.old <- w1d
        v.old <- w2d
        d2 <- mahalanobis(x,center=t1,cov=s)
        d <- sqrt(d2)
        k <- ksolve.bt(d,p,c1,M,b0)
        d <- d/k
        w1d <- wt.bt(d,c1,M)
        w2d <- v.bt(d,c1,M)
        t1 <- (w1d %*% x)/sum(w1d)
        s <- s*0
        for (i in 1:n)
        {
            xc <- as.vector(x[i,]-t1)
            s <- s + as.numeric(w1d[i])*(xc %o% xc)
        }
        s <- p*s/sum(w2d)
        mnorm <- sqrt(as.vector(t.old) %*% as.vector(t.old))
        snorm <- eigen(s.old)$values[1]
        crit1 <- max(abs(t1 - t.old))
#        crit <- max(crit1,crit2)
        crit <- max(abs(w1d-wt.old))/max(w1d)
        iter <- iter+1
    }
#    mnorm <- sqrt(as.vector(t1) %*% as.vector(t1))
#    snorm <- eigen(s)$values[1]
#    return(list(t1=t1,s=s))
s
}

cov.mba<-function(x,COR=FALSE){
val<-covmba2(x)$cov
if(COR){
val=val/outer(sqrt(diag(val)),sqrt(diag(val)))
}
val
}

covmba2<-function(x, csteps = 5)
{
# Perform the median ball algorithm.
#
# It returns a measure of location and scatter for the
# multivariate data in x, which is assumed to have
# p>-2 column and n rows.
#
# This code is based on a very slight modificatiion of code originally
# written by David Olive
#
x<-as.matrix(x)
if(!is.matrix(x))stop("x should be a matrix")
         p <- dim(x)[2]
#if(p==1)stop("x should be a matrix with two or more columns of variables")
         ##get the DGK estimator
         covs <- var(x)
         mns <- apply(x, 2, mean)        ## concentrate
         for(i in 1:csteps) {
                 md2 <- mahalanobis(x, mns, covs)
                 medd2 <- median(md2)
#                 mns <- apply(x[md2 <= medd2,  ], 2,
                 mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2,
                         mean)
                 covs <- var(x[md2 <= medd2,  ])
         }
         covb <- covs
         mnb <- mns      ##get the square root of det(covb)
         critb <- prod(diag(chol(covb)))
         ##get the resistant estimator
         covv <- diag(p)
         med <- apply(x, 2, median)
         md2 <- mahalanobis(x, center = med, covv)
         medd2 <- median(md2)    ## get the start
#         mns <- apply(x[md2 <= medd2,  ], 2, mean)
         mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
         covs <- var(x[md2 <= medd2,  ]) ## concentrate
         for(i in 1:csteps) {
                 md2 <- mahalanobis(x, mns, covs)
                 medd2 <- median(md2)
 #                mns <- apply(x[md2 <= medd2,  ], 2,mean)
       mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
                 covs <- var(x[md2 <= medd2,  ])
         }
         crit <- prod(diag(chol(covs)))
         if(crit < critb) {
                 critb <- crit
                 covb <- covs
                 mnb <- mns
         }
##scale for better performance at MVN
         rd2 <- mahalanobis(x, mnb, covb)
         const <- median(rd2)/(qchisq(0.5, p))
         covb <- const * covb
         list(center = mnb, cov = covb)
}

dcov<-function(m,tr=.2,dop=1,cop=2,pr=TRUE){
#
# Compute multivariate measure of scatter
# using Donoho-Gasko method.
#
# dop=1, use fdepth to compute depths
# dop=2, use fdepthv2  to compute depths
# which is slower but more accurate.
#
# cop=1, use Donoho-Gasko median in fdepth
# cop=2, use MCD in fdepth
# cop=3, use marginal medians in fdepth
# cop=4, use MVE in fdepth
#
if(tr>=.5)stop("Amount of trimming must be less than .5")
if(is.list(m))m<-matl(m)
if(!is.matrix(m))stop("Data must be stored in a matrix or in list mode.")
if(ncol(m)==1){
if(tr<.5)val<-mean(m,tr)
}
if(ncol(m)>1){
temp<-NA
if(ncol(m)!=2){
# Use approximate depth
if(dop==1)temp<-fdepth(m,plotit=FALSE,cop=cop)
if(dop==2)temp<-fdepthv2(m)
}
#  Use exact depth if ncol=2
if(ncol(m)==2){
for(i in 1:nrow(m))
temp[i]<-depth(m[i,1],m[i,2],m)
}}
mdep<-max(temp)
flag<-(temp==mdep)
flag2<-(temp>=tr)
if(sum(flag2)==0)stop("Trimmed all of the data")
if(sum(flag2)==1){
if(pr)print("Warning: Trimmed all but one point")
val<-0
}
if(sum(flag2)>1)val<-var(m[flag2,])
val
}

Scov<-function(m){
#
# Compute Davies' covariance S-estimator
#
library("rrcov")
res=CovSest(m)
center=res@center
z=res@cov
list(center=center,cov=z)
}

covloc<-function(x){
#
# Return mean and covarinace matrix
#
loc=apply(x,2,mean)
mcov=cov(x)
list(center=loc,cov=mcov)
}

cov.ogk<-function(x,y=NA,n.iter=1,sigmamu=taulc,v=gkcov,beta=.9,...){
#
# Compute robust (weighted) covariance matrix in Maronna and Zamar
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
temp
}

wincov<-function(m,tr=.2){
m=winall(m,tr=tr)$cov
m
}

longcov2mat<-function(x,Sid.col,dep.col){
#
# Have data in a matrix or data frame, x
# Sid.col indicates Subject's id
# Here, each subject has one or more rows of data
#
# In a regression setting, each subject has
# one or more covariates corresponding to columns.
# For example, two covariates might be stored in columns
# 3 and 6.
#
# Goal: For ith subject, store the covariate data in
# list mode, which is a matrix.
# So for ith subject, store covariate data in z[[i]], say, which
# contains a matrix of dimension  m by p,
# m is the number of observations for ith subject and p
# the number of covariates.
#
# dep.col, having length p, indicates columns containe the covariates
# Column Sid.col indicates the column containing subject's id
#
if(is.null(dim(x)))stop("x must be a matrix or data frame")
Sid=unique(x[,Sid.col])
res=list()
nid=length(Sid)
p=length(dep.col)# Number of covariates for each subject
n=nrow(x)
flag=(x[,Sid.col]==Sid[1])
n.each.s=sum(flag) # the number of rows for each subject
ns=n/n.each.s # the number of  subjects
if(!is.wholenumber(ns))stop("Not all S's have same number of rows of data")
for(i in 1:ns){
#res[[i]]=matrix(NA,nrow=n.each.s,ncol=p)
flag=(x[,Sid.col]==Sid[i])
res[[i]]=as.matrix(x[flag,dep.col])
}
res
}

cov.roc<-function(x){
library(robust)
temp<-covRob(x,estim='M')
val<-temp
list(center=val[3]$center,cov=val[2]$cov)
}

cov.funl<-function(m){
list(cov=m)
}

skip.cov<-function(x,cop = 6, MM = FALSE, op = 1, mgv.op = 0, outpro.cop = 3,
    STAND = FALSE){
ans=skipcov(x,cop=cop,MM=MM,op=op,mgv.op=mgv.op,outpro.cop=outpro.cop,STAND=STAND)
list(cov=ans)
}

wmean.cov<-function(x,tr=0){
#
# Compute Winsoriced mean and covariance for data in x
#
loc=apply(x,2,mean,tr=tr)
cv=wincov(x,tr=tr)
list(center=loc,cov=cv)
}

covl<-function(x){
res=cov(x)
list(cov=res)
}

cov2cor<-function(x){
#
# Convert a covariance matrix to a correlation matrix
#
p=ncol(x)
m=x
for(i in 1:p){
for(j in 1:p){
m[i,j]=m[i,j]/sqrt(x[i,i]*x[j,j])
}}
m
}

DETMCD<-function(x){
#
# HUber et al. (2012) deterministic version of the MCD estimator
#
library(DetMCD)
a=DetMCD(x)
list(center=a$center,cov=a$cov)
}

wincovN<-function(x,y=NULL,tr=0.2){
#
# Winsorized covariance rescaled to est cov under normality when there is no trimming
#
library(MASS)
e=wincor(x,y,tr=tr)$cov
if(tr==0)cterm=1
else cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
e=e/cterm
e
}

outDETMCD<-function(x,cov.fun=DETMCD,xlab='X',ylab='Y',qval=.975,
crit=NULL,KS=TRUE,plotit=FALSE,...){
#
#  Search for outliers using robust measures of location and scatter,
#  which are used to compute robust analogs of Mahalanobis distance.
#
#  x is an n by p matrix or a vector of data.
#
#  The function returns the values flagged as an outlier plus
#  the (row) number where the data point is stored.
#  If x is a vector, out.id=4 indicates that the fourth observation
#  is an outlier and outval=123 indicates that 123 is the value.
#  If x is a matrix, out.id=4 indicates that the fourth row of
#  the matrix is an outlier and outval reports the corresponding
#  values.
#
#  The function also returns the distance of the
#  points identified as outliers
#  in the variable dis.
#
#  For bivariate data, if plotit=TRUE, plot points and circle outliers.
#
#  cov.fun determines how the measure of scatter is estimated.
# The default is covDETMCD
#  Possible choices are
#  cov.mve (the MVE estimate)
#  cov.mcd (the MCD estimate)
#  covmba2 (the MBA or median ball algorithm)
#  rmba  (an adjustment of MBA suggested by D. Olive)
#  cov.roc (Rockes TBS estimator)
#
#  plotit=FALSE used to avoid problems when other functions in WRS call
#  this function
#
#  KS=TRUE: keep  the seed that was used
#
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

