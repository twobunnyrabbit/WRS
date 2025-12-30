# WRS Package - Bootstrap and Resampling Infrastructure
# Generic bootstrap, permutation, and resampling methods
#
# This module contains:
#   - Generic bootstrap functions: bootse, bootdep, bootcov
#   - BCA (Bias-Corrected Accelerated) methods
#   - Permutation test infrastructure
#   - Bootstrap helpers and subroutines
#
# Note: Specific bootstrap analyses (yuenbt, etc.) are in their
# respective modules (two-sample.R, anova.R, regression.R, etc.)


# bootse
bootse<-function(x,nboot=1000,est=median,SEED=TRUE,...){
#
#   Compute bootstrap estimate of the standard error of the
#   estimator est
#   The default number of bootstrap samples is nboot=100
#
if(SEED)set.seed(2) # set seed of random number generator so that
#   results can be duplicated.
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,est,...)
bootse<-sqrt(var(bvec))
bootse
}



# bootdse
bootdse<-function(x,y=NA,est=median,nboot=100,pr=TRUE,...){
#
# Determine standard error of difference between
# two measures of location for dependent groups.
#
if(!is.na(y[1]))x<-cbind(x,y)
x<-elimna(x)
if(pr)print("Taking Bootstrap Samples. Please wait.")
data<-matrix(sample(nrow(x),size=nrow(x)*nboot,replace=TRUE),nrow=nboot)
xmat<-matrix(x[data,1],nrow=nboot,ncol=length(x))
ymat<-matrix(x[data,2],nrow=nboot,ncol=length(x))
bvec<-apply(xmat,1,FUN=est,...)-apply(ymat,1,FUN=est,...)
se<-sqrt(var(bvec))
se
}


# bootdpci
bootdpci<-function(x,y,est=onestep,nboot=NA,alpha=.05,plotit=FALSE,dif=TRUE,BA=FALSE,SR=TRUE,...){
#
#   Use percentile bootstrap method,
#   compute a .95 confidence interval for the difference between
#   a measure of location or scale
#   when comparing two dependent groups.
#   By default, a one-step M-estimator (with Huber's psi) is used.
#   If, for example, it is desired to use a fully iterated
#   M-estimator, use fun=mest when calling this function.
#
okay=FALSE
if(identical(est,onestep))okay=TRUE
if(identical(est,mom))okay=TRUE
if(!okay)SR=FALSE
output<-rmmcppb(x,y,est=est,nboot=nboot,alpha=alpha,SR=SR,
plotit=plotit,dif=dif,BA=BA,...)$output
list(output=output)
}



# bootdep
bootdep<-function(x,tr=.2,nboot=500){
#
# x is a matrix (n by p) or has list mode
# Goal: Obtain boostrap samples and compute
# the trimmed each for each of the p variables.
# Return the bootstrap means in a matrix
#
# tr is the amount of trimming
# nboot is the number of bootstrap samples
#
if(is.matrix(x))m1<-x
if(is.list(x)){
# put the data into a matrix
m1<-matrix(NA,ncol=length(x))
for(j in 1:length(x))m1[,j]<-x[[j]]
}
data<-matrix(sample(nrow(m1),size=nrow(m1)*nboot,replace=TRUE),nrow=nboot)
bvec<-matrix(NA,ncol=ncol(m1),nrow=nboot)
for(j in 1:ncol(m1)){
temp<-m1[,j]
bvec[,j]<-apply(data, 1., bootdepsub,temp,tr)
}
# return a nboot by p matrix of bootstrap trimmed means.
bvec
}


# bootdepsub
bootdepsub<-function(isub,x,tr){
tsub<-mean(x[isub],tr)
tsub
}

# bootcov
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

  yuenbt<-function(x,y,tr=.2,alpha=.05,nboot=199,side=FALSE,nullval=0,pr=TRUE,
plotit=FALSE,op=1,SEED=TRUE){
#
#  Compute a 1-alpha confidence interval for the difference between
#  the trimmed means corresponding to two independent groups.
#  The bootstrap-t method is used.
#
#  The default amount of trimming is tr=.2
#  side=T indicates two-sided method using absolute value of the
#  test statistics within the bootstrap; otherwise the equal-tailed method
#  is used.
#
#  This function uses trimse.
#
side<-as.logical(side)
p.value<-NA
ybt<-vector(mode="numeric",length=2)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
x<-x[!is.na(x)]  # Remove missing values in x
y<-y[!is.na(y)]  # Remove missing values in y
xcen<-x-mean(x,tr)
ycen<-y-mean(y,tr)
test<-(mean(x,tr)-mean(y,tr))/sqrt(trimse(x,tr=tr)^2+trimse(y,tr=tr)^2)
datax<-matrix(sample(xcen,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(ycen,size=length(y)*nboot,replace=TRUE),nrow=nboot)
top<-apply(datax,1,mean,tr)-apply(datay,1,mean,tr)
botx<-apply(datax,1,trimse,tr)
boty<-apply(datay,1,trimse,tr)
tval<-top/sqrt(botx^2+boty^2)
if(plotit){
if(op == 1)
akerd(tval)
if(op == 2)
rdplot(tval)
}
if(side)tval<-abs(tval)
tval<-sort(tval)
icrit<-floor((1-alpha)*nboot+.5)
ibot<-floor(alpha*nboot/2+.5)
itop<-floor((1-alpha/2)*nboot+.5)
se<-sqrt((trimse(x,tr))^2+(trimse(y,tr))^2)
ybt[1]<-mean(x,tr)-mean(y,tr)-tval[itop]*se
ybt[2]<-mean(x,tr)-mean(y,tr)-tval[ibot]*se
if(side){
CI=mean(x,tr)-mean(y,tr)-tval[icrit]*se
ybt[1]<-mean(x,tr)-mean(y,tr)-tval[icrit]*se
ybt[2]<-mean(x,tr)-mean(y,tr)+tval[icrit]*se
CI[2]=ybt[2]
p.value<-(sum(abs(test)<=abs(tval)))/nboot
if(ybt[1]>nullval)p.value=min(p.value,alpha)  #very remote chance ci and p.value differ. Force them to agree.
if(ybt[2]<nullval)p.value=min(p.value,alpha)
}
if(!side){
ibot<-round(alpha*nboot/2)
itop<-nboot-ibot+1
ibot=ibot+1
crit=c(tval[ibot],tval[itop])
CI=mean(x,tr)-mean(y,tr)-tval[itop]*se
ybt=mean(x,tr)-mean(y,tr)-tval[itop]*se
ybt[2]<-mean(x,tr)-mean(y,tr)-tval[ibot]*se
CI[2]=mean(x,tr)-mean(y,tr)-tval[ibot]*se
if(test<0)G=mean(test>tval[1:nboot])
if(test>=0)G=mean(test<tval[1:nboot])
p.value=2*G
}
list(ci=CI,test.stat=test,p.value=p.value,est.1=mean(x,tr),est.2=mean(y,tr),est.dif=mean(x,tr)-mean(y,tr),
n1=length(x),n2=length(y))
}



# boot.TM
boot.TM<-function(x,nboot=599,alpha=.05,SEED=TRUE){
#
# Global test for equal M-measures of location, J independent groups
#
# This is method TM in 5th Ed of Intro to Robust Estimation and Testing
#
if(SEED)set.seed(2)
B=nboot
  if(is.matrix(x) || is.data.frame(x))xlist=listm(x)
  else xlist=x
  xlist=elimna(xlist)
  T.test<-TM(xlist)$TM
  k<-length(xlist)
  ylist<-vector(mode="list",length=k)
  TT<-numeric(B)
  b<-floor((1-alpha)*B)
  onesteps<-sapply(xlist,onestep)
  for (i in 1:B){
   j<-1
   repeat {
    ylist[[j]]<-(sample(xlist[[j]],length(xlist[[j]]),replace=T)-onesteps[j])
    if (mad(ylist[[j]])>0) j<-j+1 #MAD must be greater than zero for every bootstrap sample
    if (j>k)break
   }
   TT[i]<-TM(ylist,alpha)$TM
  }
  TT=sort(TT)
  if(T.test>=TT[b]){1} else{0}
pv=mean(T.test<=TT)
list(Est.=onesteps,p.value=pv)
}




# regboot
regboot<-function(isub,x,y,regfun,...){
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
xmat<-matrix(x[isub,],nrow(x),ncol(x))
vals<-regfun(xmat,y[isub],...)$coef
vals
}



# regbootg
regbootg<-function(isub,x,y,regfun,...){
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
#xmat<-matrix(x[isub,],nrow(x),ncol(x))
xmat<-x[isub,]
yy<-y[isub]
regboot<-regfun(xmat,y[isub])$coef
#regboot<-regboot$coef
regboot
}


# regbootMC
regbootMC<-function(data,x,y,regfun,...){
vals=regfun(x[data,],y[data],...)$coef
}


# ancboot
ancboot<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,nboot=599,pts=NA,plotit=TRUE,xout=FALSE,outfun=outpro,...){
#
# Compare two independent  groups using the ancova method
# in chapter 12 of Wilcox, 2017, Intro to Robust Estimation and Hypothesis Testing.
# No assumption is made about the form of the regression
# lines--a running interval smoother is used.
# Confidence intervals are computed using a bootstrap-t bootstrap
# method. Comparisons are made at five empirically chosen design points.
#
#  Assume data are in x1 y1 x2 and y2
#
if(is.na(pts[1])){
isub<-c(1:5)  # Initialize isub
test<-c(1:5)
m1=elimna(cbind(x1,y1))
x1=m1[,1]
y1=m1[,2]
m1=elimna(cbind(x2,y2))
x2=m1[,1]
y2=m1[,2]
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
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","ci.low","ci.hi",
"p.value"))
gv1<-vector("list")
for (i in 1:5){
j<-i+5
temp1<-y1[near(x1,x1[isub[i]],fr1)]
temp2<-y2[near(x2,x1[isub[i]],fr2)]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
mat[i,2]<-length(temp1)
mat[i,3]<-length(temp2)
gv1[[i]]<-temp1
gv1[[j]]<-temp2
}
I1<-diag(5)
I2<-0-I1
con<-rbind(I1,I2)
test<-linconb(gv1,con=con,tr=tr,nboot=nboot)
for(i in 1:5){
mat[i,1]<-x1[isub[i]]
}
mat[,4]<-test$psihat[,2]
mat[,5]<-test$test[,2]
mat[,6]<-test$psihat[,3]
mat[,7]<-test$psihat[,4]
mat[,8]<-test$test[,4]
}
if(!is.na(pts[1])){
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
if(n1[i]<=5)paste("Warning, there are",n1[i]," points corresponding to the design point X=",pts[i])
if(n2[i]<=5)paste("Warning, there are",n2[i]," points corresponding to the design point X=",pts[i])
}
mat<-matrix(NA,length(pts),9)
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","se","ci.low","ci.hi",
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
test<-linconb(gv,con=con,tr=tr,nboot=nboot)
mat[,1]<-pts
mat[,2]<-n1
mat[,3]<-n2
mat[,4]<-test$psihat[,2]
mat[,5]<-test$test[,2]
mat[,6]<-test$test[,3]
mat[,7]<-test$psihat[,3]
mat[,8]<-test$psihat[,4]
mat[,9]<-test$test[,4]
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
runmean2g(x1,y1,x2,y2,fr=fr1,est=mean,tr=tr)
}
list(output=mat,crit=test$crit)
}


# ancbootg
ancbootg<-function(x1,y1,x2,y2,pts,fr1=1,fr2=1,tr=.2,nboot=599){
#
# Compare two independent  groups using the ancova method
# in chapter 9. No assumption is made about the form of the regression
# lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#  Comparisons are made at the design points contained in the vector
#  pts
#
m1=elimna(cbind(x1,y1))
x1=m1[,1]
y1=m1[,2]
m1=elimna(cbind(x2,y2))
x2=m1[,1]
y2=m1[,2]
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),8)
dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","se","ci.low","ci.hi"))
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
test<-linconb(gv,con=con,tr=tr,nboot=nboot)
mat[,1]<-pts
mat[,2]<-n1
mat[,3]<-n2
mat[,4]<-test$psihat[,2]
mat[,5]<-test$test[,2]
mat[,6]<-test$test[,3]
mat[,7]<-test$psihat[,3]
mat[,8]<-test$psihat[,4]
list(output=mat,crit=test$crit)
}


# bca.mean
wmw.bca<-function(x,y,alpha=.05,nboot=1000,SEED=TRUE,...){
#
#  BCA confidence interval for P(X<Y), X and Y independent
#  Method: Bias corrected accelerated bootstrap
#
library(bcaboot)
if(SEED)set.seed(2)
n1=length(x)
n2=length(y)
nmax=max(n1,n2)
m=matrix(NA,nmax,2)
m[1:n1,1]=x
m[1:n2,2]=y
e=bmp(x,y)$phat
if(e==0 || e==1)ci=cid(x,y)$ci.p
else{
a=bcajack(m,B=500,wmw.est.only,n1=n1,n2=n2,alpha=alpha/2,verbose=FALSE,...)
ci=c(a$lims[1,1],a$lims[3,1])
}
list(n1=n1,n2=n2,phat=e,ci.low=ci[1],ci.upper=ci[2])
}



# wmw.bcapv
wmw.bcapv<-function(x,y,alpha=.05,nboot=1000,SEED=TRUE,...){
#
#  BCA confidence interval for P(X<Y), X and Y independent
#  Method: Bias corrected accelerated bootstrap
#
library(bcaboot)
if(SEED)set.seed(2)
ALPHA=c(seq(.001,.1,.001),seq(.011,.05,.01),seq(.11,.99,.01))
n1=length(x)
n2=length(y)
nmax=max(n1,n2)
m=matrix(NA,nmax,2)
m[1:n1,1]=x
m[1:n2,2]=y
al=length(ALPHA)
il=0
est=bmp(x,y)
a=bcajack(m,nboot,wmw.est.only,alpha=ALPHA/2,verbose=FALSE,...)
iu=nrow(a$lims)+1
for(j in 1:al){
il=il+1
iu=iu-1
pv=ALPHA[j]
ci=c(a$lims[il,1],a$lims[iu,1])
if(a$lims[il,1]>0.5 || a$lims[iu,1]<0.5)break
}
A=wmw.bca(x,y,alpha=alpha,nboot=nboot,SEED=SEED)
list(n1=n1,n2=n2,phat=est$phat,ci.low=A$ci.low,ci.upper=A$ci.upper,p.value=pv)
}



# wmw.bcapv.v2
wmw.bcapv.v2<-function(x,y,nboot=500,iter=100,SEED=TRUE,...){
#
#  BCA  p-value for H_0P (X<Y)=.5, X and Y independent
#  Method: Bias corrected accelerated bootstrap
#
library(bcaboot)
if(SEED)set.seed(2)
n1=length(x)
n2=length(y)
nmin=min(n1,n2)
m=matrix(NA,nmin,2)
est=bmp(x,y)$phat
pv.vals=NA
for(i in 1:iter){
ys=sample(y,nmin)
m[1:nmin,1]=sample(x,nmin)
m[1:nmin,2]=sample(y,nmin)
pv.vals[i]=wmw.bcapv(m[,1],m[,2],SEED=FALSE)$p.value
}
list(n1=n1,n2=n2,phat=est,p.value=mean(pv.vals))
}

# wmw.bcav2
wmw.bcav2<-function(x,y,alpha=.05,nboot=500,iter=100,SEED=TRUE,...){
#
#  BCA confidence interval for P(X<Y), X and Y independent
#  Method: Bias corrected accelerated bootstrap
#  Uses a modification when the estimate is zero or one
#
library(bcaboot)
if(SEED)set.seed(2)
n1=length(x)
n2=length(y)
nmin=min(n1,n2)
cimat=matrix(NA,iter,2)
m=matrix(NA,nmin,2)
e=bmp(x,y)$phat
if(e==1 ||e==0){
ci=cid(x,y)$ci.p
}
else{
for(i in 1:iter){
m[,1]=sample(x,nmin)
m[,2]=sample(y,nmin)
a=bcajack(m,B=nboot,wmw.est.only,alpha=alpha/2,verbose=FALSE,...)
cimat[i,]=c(a$lims[1,1],ci.upper=a$lims[3,1])
}
ci=mean(elimna(cimat[,1]))
ci[2]=mean(elimna(cimat[,2]))
}
list(n1=n1,n2=n2,phat=e,ci.low=ci[1],ci.upper=ci[2])
}


# corblp.bca.C
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


# pbcan
pbcan<-function(x,nboot=1000,grp=NA,est=onestep,...){
#
#   Test the hypothesis that J independent groups have
#   equal measures of location using the percentile bootstrap method.
#   in conjunction with a partially centering technique.
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   est is the measure of location and defaults to an M-estimator
#   ... can be used to set optional arguments associated with est
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#   Missing values are allowed.
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
if(!is.na(sum(grp))){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[1]]]
x<-xx
}
J<-length(x)
tempn<-0
vecm<-0
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
vecm[j]<-est(x[[j]],...)
}
xcen<-list()
flag<-rep(TRUE,J)
for(j in 1:J){
flag[j]<-FALSE
temp<-mean(vecm[flag])
xcen[[j]]<-x[[j]]-temp
flag[j]<-T
}
icrit<-round((1-alpha)*nboot)
bvec<-matrix(NA,nrow=J,ncol=nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
for(j in 1:J){
paste("Working on group ",j)
data<-matrix(sample(xcen[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,est,...) # Bootstrapped values for jth group
}
vvec<-NA
for(j in 1:J){
vvec[j]<-sum((bvec[j,]-vecm[j])^2)/(nboot-1)
}
dis<-NA
for(i in 1:nboot){
dis[i]<-sum((bvec[,i]-vecm)^2/vvec)
}
tvec<-sum((0-vecm)^2/vvec)
dis<-sort(dis)
print(tvec)
print(dis[icrit])
print(vecm)
sig<-1-sum((tvec>=dis))/nboot
list(p.value=sig)
}


# permg
permg<-function(x,y,alpha=.05,est=mean,nboot=1000){
#
# Do a two-sample permutation test based on means or any
# other measure of location or scale indicated by the
# argument est.
#
# The default number of permutations is nboot=1000
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
xx<-c(x,y)
dif<-est(x)-est(y)
vec<-c(1:length(xx))
v1<-length(x)+1
difb<-NA
temp2<-NA
for(i in 1:nboot){
data <- sample(xx, size = length(xx), replace = FALSE)
temp1<-est(data[c(1:length(x))])
temp2<-est(data[c(v1:length(xx))])
difb[i]<-temp1-temp2
}
difb<-sort(difb)
icl<-floor((alpha/2)*nboot+.5)
icu<-floor((1-alpha/2)*nboot+.5)
reject<-"no"
if(dif>=difb[icu] || dif <=difb[icl])reject<-"yes"
list(dif=dif,lower=difb[icl],upper=difb[icu],reject=reject)
}



# permg.t
permg.t<-function(x,y,alpha=.05,tr=0,nboot=1000,SEED=TRUE){
#
# Do a two-sample permutation test based on trimmed means using the
# Chung--Romano version of a permuation test.

# The default number of permutations is nboot=1000
#
if(SEED)set.seed(2)
x<-x[!is.na(x)]
y<-y[!is.na(y)]
xx<-c(x,y)
tval<-yuen(x,y,tr=tr)$teststat
vec<-c(1:length(xx))
v1<-length(x)+1
difb<-NA
tv<-NA
for(i in 1:nboot){
data <- sample(xx, size = length(xx), replace = FALSE)
temp1<-data[c(1:length(x))]
temp2<-data[c(v1:length(xx))]
tv[i]<-yuen(temp1,temp2,tr=tr)$teststat
}
tv<-sort(tv)
icl<-floor((alpha/2)*nboot+.5)
icu<-floor((1-alpha/2)*nboot+.5)
reject<-'no'
list(teststat=tval,lower.crit=tv[icl],upper.crit=tv[icu],reject=reject)
}


# perm.rho
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


# perm.rho.sub
perm.rho.sub<-function(x,y){
rho=cor(x,y)
n=length(x)
xbar=mean(x)
ybar=mean(y)
m22=sum((x-xbar)^2*(y-ybar)^2)/n
m20=sum((x-xbar)^2)/n
m02=sum((y-ybar)^2)/n
tau=sqrt(m22/(m20*m02))
S=sqrt(n)*rho/tau
S
}


# bwdepth.perm
linWMWMC.sub<-function(M,con){
L=apply(t(con*t(M)),1,sum)
L
}


# linWMWMC.sub2
linWMWMC.sub2<-function(L){
phat=mean(L<0)+.5*mean(L==0)
phat
}



# ridgeGnullMC.sub
ridgeGnullMC.sub<-function(x,p,regfun=regfun){
p1=p+1
v=ridgeG.sub(x[,1:p],x[,p1],regfun=regfun)$F.test
v
}



# btsqrk
btsqrk<-function(alist,alpha=0.05,tr=0.2){
#computes B2_tk test statistics for k independent samples.
#alist should be a list type object
#s's are computed by trimse which can be found in all Rallfun files written by Wilcox Rand
k<-length(alist)
# Remove any missing values in alist
for (i in 1:k){alist[[i]]<-alist[[i]][!is.na(alist[[i]])]}
zc<-qnorm(alpha/2)
e=trunc(tr*sapply(alist,length))
f<-(sapply(alist,length))-(2*e)
s=sapply(alist,trimse,tr=tr)^2
wden=sum(1/s)
w=(1/s)/wden
yplus<-sum(w*(sapply(alist,mean,trim=tr)))
tt<-((sapply(alist,mean,trim=tr))-yplus)/sqrt(s)
v<-(f-1)
z<-((4*v^2)+(5*((2*(zc^2))+3)/24))/((4*v^2)+v+(((4*(zc^2))+9)/12))*sqrt(v)*(sqrt(log(1+(tt^2/v))))
teststat<-sum(z^2)
crit<-qchisq(1-alpha,k-1)
bt2pvalue<-1-(pchisq(teststat,k-1))
list(p.value=bt2pvalue,teststat=teststat,crit=crit,e=e,f=f,s=s,w=w,tt=tt)
}


