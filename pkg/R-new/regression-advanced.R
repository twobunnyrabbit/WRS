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
gamplot<-function(x,y,sop=TRUE,pyhat=FALSE,eout=FALSE,xout=FALSE,outfun=out,plotit=TRUE,
xlab="X",ylab="",zlab="",theta=50,phi=25,expand=.5,scale=TRUE,ticktype="simple"){
#
# Plot regression surface using generalized additive model
#
# sop=F, use usual linear model y~x1+x2...
# sop=T, use splines
#
library(akima)
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
mulgreg<-function(x,y,cov.fun=rmba){
#
# Do Multivariate regression in Rousseeuw, Van Aelst, Van Driessen Agullo
# (2004) Technometrics, 46, 293-305
#
# (y can be multivariate)
#
library(MASS)
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
# mregdepth
# ============================================================================
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


# ============================================================================
# regpord.sub
# ============================================================================
regpord.sub<-function(isub,x,y,cov.fun){
xmat<-matrix(x[isub,],nrow(x),ncol(x))
vals<-regvarp(xmat,y[isub],cov.fun=cov.fun)
vals
}


# ============================================================================
# regvarp
# ============================================================================
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
library(akima)
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
library(akima)
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
mlrregCI<-function(x,y,nboot=300,MC=FALSE,SEED=TRUE,op.dis=TRUE){
#
#  Based on Rousseeuw et al.
#  multivariate regression estimator
#  compute p-value for each of the parameters using a percentile
#  bootstrap method.
#
if(SEED)set.seed(2)
if(MC)library(parallel)
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
mlrreg.est<-function(data,x,y){
xv=x[data,]
yv=y[data,]
vals=as.vector(mlrreg(xv,yv)$coef)
vals
}


# ============================================================================
# mlrregWtest
# ============================================================================
mlrregWtest<-function(x,y,nboot=300,MC=FALSE,SEED=TRUE){
#
#  Test hypothesis that all slopes=0  based on Rousseeuw et al.
#  multivariate regression estimator
#
#  Strategy: a variation of the wild bootstrap method, percentile version.
#
if(SEED)set.seed(2)
if(MC)library(parallel)
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
mlrreg.subest<-function(data,x){
vals=as.vector(mlrreg(x,data)$coef[-1,])
vals
}


# ============================================================================
# regmediate
# ============================================================================
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
if(MC)library(parallel)
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
if(MC)library(parallel)
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
library(MASS)
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
# qhdplotsm
# ============================================================================
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
regYciCV<-function(n,alpha=.05,nboot=1000,regfun=tsreg,SEED=TRUE,MC=FALSE,null.value=0,xout=FALSE,...){
#
# Determine a critical value for regYci
#
if(SEED)set.seed(2)
mv=NA
chk=0
if(MC)library(parallel)
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
regYciCV2G<-function(n1,n2,crit=NULL,g=0,h=0,nboot=1000,regfun=tsreg,ALL=TRUE,
alpha=.05,SEED=TRUE,MC=FALSE,null.value=0,pts=NULL,npts=100,nmiss=0,...){
n=max(n1,n2)
if(nmiss>n)stop('Number of missing values is greater than max(n1,n2)')
if(SEED)set.seed(2)
mv=NA
chk=0
if(n1!=n2)nmiss=max(c(n1,n2))-min(c(n1,n2))
if(MC)library(parallel)
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
regY2G.sub<-function(xy,regfun,null.value=0,...){
 pv=regYci2Gv2(xy[,1],xy[,2],xy[,3],xy[,4],SEED=FALSE,regfun=regfun,null.value=null.value,...)[,4]
 min(pv)
}


# ============================================================================
# regYci2Gv2
# ============================================================================
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
library(parallel)
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
regIVcom_sub<-function(slope,x,tr){
yhat=apply(t(slope*t(x)),1,sum)
str=winvar(yhat,tr=tr)
str
}


# ============================================================================
# regIVstr
# ============================================================================
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
# regcor
# ============================================================================
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


# ============================================================================
# regstr
# ============================================================================
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
regIVbinv2_sub<-function(data,x,y,pts){
x=as.matrix(x)
v=logreg.pred(x[data,],y[data],pts=pts)
v=sd(v)
}


# ============================================================================
# reg.hyp.split
# ============================================================================
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
# smean.depth
# ============================================================================
smean.depth<-function(m){
#
#  Skipped estimator based on projection for removing outliers.
#  Uses random projections
#
m=elimna(m)
id=outpro.depth(m)$keep
val=apply(m[id,],2,mean)
val
}


# ============================================================================
# smeancr.cord.oph
# ============================================================================
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


# ============================================================================
# smeancr.cord
# ============================================================================
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


# ============================================================================
# reg.resid
# ============================================================================
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
regIQR<-function(x,y,xr=x,regfun=Qreg,xout=FALSE,outfun=outpro,...){
#
#
IQR=regYhat(x,y,xr=xr,regfun=regfun,q=.75)-regYhat(x,y,xr=xr,regfun=regfun,q=.25)
IQR
}


# ============================================================================
# qinvreg
# ============================================================================
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
qinvreg.sub<-function(xy,q,pt,v){
e=reg.pred(xy[,1],xy[,2],xr=pt,regfun=Qreg,q=q,xout=FALSE)
a=abs(e-v)
a
}


# ============================================================================
# reg.con.dist
# ============================================================================
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


