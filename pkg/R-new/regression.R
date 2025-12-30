# WRS Package - Regression Methods
# Extracted from Rallfun-v45.R
#
# This module contains robust regression methods including:
# - Theil-Sen regression (tsreg, tshdreg, etc.)
# - Least trimmed squares (ltsreg, LMSreg, MMreg)
# - M-regression (chreg, bmreg, bireg, winreg)
# - Outlier-pruned regression (opreg, mopreg)
# - Skipped regression (snmreg)
# - Depth-based regression (depreg, mdepreg)
# - Quantile regression basics (qreg, Qreg)
# - Regression inference (regci, regtest, lintest)
# - Two-group comparisons (difreg, reg2ci)
# - One-way regression ANOVA (reg1way)
#
# Total functions: 98
# Extraction date: 2025-12-30


# ============================================================================
# lintestMC
# ============================================================================
lintestMC<-function(x,y,regfun=tsreg,nboot=500,alpha=.05,xout=FALSE,outfun=out,...){
#
# Test the hypothesis that the regression surface is a plane.
# Stute et al. (1998, JASA, 93, 141-149).
#
library(parallel)
set.seed(2)
if(identical(regfun,tshdreg))print('When using tshdreg, be sure to include RES=TRUE')
#if(identical(regfun,Qreg))print('When using Qreg, be sure to include res.vals=TRUE')
x<-as.matrix(x)
d<-ncol(x)
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
x<-as.matrix(x)
y<-temp[,d+1]
if(xout){
flag<-outfun(x)$keep
x<-x[flag,]
x<-as.matrix(x)
y<-y[flag]
}
mflag<-matrix(NA,nrow=length(y),ncol=length(y))
for (j in 1:length(y)){
for (k in 1:length(y)){
mflag[j,k]<-(sum(x[j,]<=x[k,])==ncol(x))
}
}
reg<-regfun(x,y,...)
yhat<-y-reg$residuals
print("Taking bootstrap samples, please wait.")
data<-matrix(runif(length(y)*nboot),nrow=nboot)
data<-sqrt(12)*(data-.5) # standardize the random numbers.
data=listm(t(data))
rvalb<-mclapply(data,lintests1,yhat,reg$residuals,mflag,x,regfun,mc.preschedule=TRUE,...)
# An n x nboot matrix of R values
rvalb=matl(rvalb)
rvalb<-rvalb/sqrt(length(y))
dstatb<-apply(abs(rvalb),2,max)
wstatb<-apply(rvalb^2,2,mean)
# compute test statistic
v<-c(rep(1,length(y)))
rval<-lintests1(v,yhat,reg$residuals,mflag,x,regfun,...)
rval<-rval/sqrt(length(y))
dstat<-max(abs(rval))
wstat<-mean(rval^2)
ib<-round(nboot*(1-alpha))
p.value.d<-1-sum(dstat>=dstatb)/nboot
p.value.w<-1-sum(wstat>=wstatb)/nboot
list(dstat=dstat,wstat=wstat,p.value.d=p.value.d,p.value.w=p.value.w)
}

# ============================================================================
# bireg
# ============================================================================
bireg<-function(x,y,iter=20,bend=1.28){
#
# Compute a biweight midregression equation
# The predictors are assumed to be stored in the n by p matrix x.
#
x<-as.matrix(x)
ma<-matrix(0,ncol(x),1)
m<-matrix(0,ncol(x),ncol(x))
mvals<-apply(x,2,mest,bend)
for (i in 1:ncol(x)){
ma[i,1]<-bicov(x[,i],y)
for (j in 1:ncol(x))m[i,j]<-bicov(x[,i],x[,j])
}
slope<-solve(m,ma)
b0<-mest(y,bend)-sum(slope%*%mvals)
for(it in 1:iter){
res<-y-x%*%slope-b0
for (i in 1:ncol(x))ma[i,1]<-bicov(x[,i],res)
slopeadd<-solve(m,ma)
b0add<-mest(res,bend)-sum(slopeadd%*%mvals)
if(max(abs(slopeadd),abs(b0add)) <.0001)break
slope<-slope+slopeadd
b0<-b0+b0add
}
if(max(abs(slopeadd),abs(b0add)) >=.0001)
paste("failed to converge in",iter,"iterations")
list(coef=c(b0,slope),residuals=res)
}

# ============================================================================
# chreg
# ============================================================================
chreg<-function(x,y,bend=1.345,SEED=TRUE,xout=FALSE,outfun=outpro,pr=TRUE,...){
#
# Compute Coakley Hettmansperger robust regression estimators
# JASA, 1993, 88, 872-880
#
# x is a n by p matrix containing the predictor values.
#
# No missing values are allowed
#
#  Comments in this function follow the notation used
#  by Coakley and Hettmansperger
#
library(MASS)
# with old version of R, need library(lqs) when using ltsreg
# as the initial estimate.
#
if(pr)print('If using chreg with a bootstrap method, use chregF instead')
if(SEED)set.seed(12) # Set seed so that results are always duplicated.
x<-as.matrix(x)
p<-ncol(x)
m<-elimna(cbind(x,y))
x<-m[,1:p]
p1<-p+1
y<-m[,p1]
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
x<-as.matrix(x)
cutoff<-bend
mve<-vector("list")
if(ncol(x)==1){
mve$center<-median(x)
mve$cov<-mad(x)^2
}
if(ncol(x)>=2)mve<-cov.mve(x)  # compute minimum volume ellipsoid measures of
                 # location and scale and store in mve.
reg0<-ltsReg(x,y) # compute initial regression est using least trimmed
                 # squares.
# Next, compute the rob-md2(i) values and store in rob
rob<-1  # Initialize vector rob
mx<-mve$center
rob<-mahalanobis(x,mx,mve$cov)
k21<-qchisq(.95,p)
c62<-k21/rob
vecone<-c(rep(1,length(y))) # Initialize vector vecone to 1
c30<-pmin(vecone,c62)  # mallows weights put in c30
k81<-median(abs(reg0$residuals)) # median of absolute residuals
k72<-1.4826*(1+(5/(length(y)-p-1)))*k81 # lms scale
c60<-reg0$residuals/(k72*c30) # standardized residuals
#  compute psi and store in c27
cvec<-c(rep(cutoff,length(y))) # Initialize vector cvec to cutoff
c27<-pmin(cvec,c60)
c27<-pmax(-1*cutoff,c27)  #c27 contains psi values
#
# compute B matrix and put in c66.
#  Also, transform B so that i th diag elem = 0 if c27[i] is
#  between -cutoff and cutoff, 1 otherwise.
#
c66<-ifelse(abs(c27)<=bend,1,0) # Have derivative of psi in c66
m1<-cbind(1,x)  # X matrix with col of 1's added
m2<-t(m1)   #X transpose
m5<-diag(c30) # matrix W, diagonal contains weights
m4<-diag(c66) # B matrix
m6<-m4%*%m1   # BX
m7<-m2%*%m6   # X'BX (nD=X'BX)
m8<-solve(m7)  #m8 = (X'-B-X)inverse
m9<-m8%*%m2 #m9=X prime-B-X inverse X'
m9<-m9%*%m5 # m9=X prime-B-X inverse X'W
m10<-m9%*%c27
c20<-m10*k72
c21<-reg0$coef+c20 #update initial estimate of parameters.
res<-y-m1%*%c21
list(coef=t(c21),residuals=res)
}

# ============================================================================
# regboot
# ============================================================================
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

# ============================================================================
# bmreg
# ============================================================================
bmreg<-function(x,y,iter=20,bend=2*sqrt((ncol(x)+1)/nrow(x)),xout=FALSE,outfun=outpro,...){
# compute a bounded M regression using Huber Psi and Schweppe weights.
# The predictors are assumed to be stored in the n by p matrix x.
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
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x<-as.matrix(x)
init<-lsfit(x,y)
resid<-init$residuals
x1<-cbind(1,x)
nu<-sqrt(1-hat(x1))
low<-ncol(x)+1
for(it in 1:iter){
ev<-sort(abs(resid))
scale<-median(ev[c(low:length(y))])/qnorm(.75)
rov<-(resid/scale)/nu
psi<-ifelse(abs(rov)<=bend,rov,bend*sign(rov))  # Huber Psi
wt<-nu*psi/(resid/scale)
new<-lsfit(x,y,wt)
if(max(abs(new$coef-init$coef))<.0001)break
init$coef<-new$coef
resid<-new$residuals
}
resid<-y-x1%*%new$coef
if(max(abs(new$coef-init$coef))>=.0001)
paste("failed to converge in",iter,"steps")
list(coef=new$coef,residuals=resid,w=wt)
}

# ============================================================================
# reglev
# ============================================================================
reglev<-function(x,y,plotit=TRUE,SEED=TRUE,DIS=FALSE){
#
#  Search for good and bad leverage points using the
#  Rousseuw and van Zomeren method.
#
#  x is an n by p matrix
#
#  The function returns the number of the rows in x that are identified
#  as outliers. (The row numbers are stored in outliers.)
#  It also returns the distance of the points identified as outliers
#  in the variable dis.
#
library(MASS)
xy=elimna(cbind(x,y))
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
plotit<-as.logical(plotit)
if(SEED)set.seed(12)
x<-as.matrix(x)
res<-lmsreg(x,y)$resid
sighat<-sqrt(median(res^2))
sighat<-1.4826*(1+(5/(length(y)-ncol(x)-1)))*sighat
stanres<-res/sighat
if(ncol(x)>=2)mve<-cov.mve(x)
if(ncol(x)==1){
mve<-vector("list")
mve$center<-median(x)
mve$cov<-mad(x)^2
}
dis<-mahalanobis(x,mve$center,mve$cov)
dis<-sqrt(dis)
crit<-sqrt(qchisq(.975,ncol(x)))
chk<-ifelse(dis>crit,1,0)
vec<-c(1:nrow(x))
id<-vec[chk==1]
chkreg<-ifelse(abs(stanres)>2.5,1,0)
idreg<-vec[chkreg==1]
if(plotit){
plot(dis,stanres,xlab="Robust distances",ylab="standardized residuals")
abline(-2.5,0)
abline(2.5,0)
abline(v=crit)
}
all=c(id,idreg)
ID=duplicated(all)
blp=all[ID]
vec=c(1:length(y))
nkeep=vec
if(length(blp)>0)nkeep=vec[-blp]
if(!DIS)dis=NULL
list(levpoints=id,regout=idreg,bad.lev.points=blp,keep=nkeep,dis=dis,stanres=stanres,crit=crit)
}

# ============================================================================
# winreg
# ============================================================================
winreg<-function(x,y,iter=20,tr=.2,xout=FALSE,outfun=outpro,...){
#
# Compute a Winsorized regression estimator
# The predictors are assumed to be stored in the n by p matrix x.
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
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x=as.matrix(x)
ma<-matrix(0,ncol(x),1)
m<-matrix(0,ncol(x),ncol(x))
mvals<-apply(x,2,win,tr)
for (i in 1:ncol(x)){
ma[i,1]<-wincor(x[,i],y,tr=tr)$cov
for (j in 1:ncol(x))m[i,j]<-wincor(x[,i],x[,j],tr=tr)$cov
}
slope<-solve(m,ma)
b0<-win(y,tr)-sum(slope%*%mvals)
for(it in 1:iter){
res<-y-x%*%slope-b0
for (i in 1:ncol(x))ma[i,1]<-wincor(x[,i],res,tr=tr)$cov
slopeadd<-solve(m,ma)
b0add<-win(res,tr)-sum(slopeadd%*%mvals)
if(max(abs(slopeadd),abs(b0add)) <.0001)break
slope<-slope+slopeadd
b0<-b0+b0add
}
if(max(abs(slopeadd),abs(b0add)) >=.0001)
paste("failed to converge in",iter,"iterations")
list(coef=c(b0,slope),resid=res)
}

# ============================================================================
# regpres1
# ============================================================================
regpres1<-function(isub,x,y,regfun,mval){
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
xmat<-matrix(x[isub,],mval,ncol(x))
regboot<-regfun(xmat,y[isub])
regboot<-regboot$coef
regboot
}

# ============================================================================
# hratio
# ============================================================================
hratio<-function(x,y,regfun=bmreg){
#
#   Compute a p by p matrix of half-slope ratios
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
#  OUTPUT:
#The first row reports the half-slope
#ratios when the data are divided into two groups using the first predictor.
#The first column is the half-slope ratio for the first predictor, the
#second column is the half-slope ratio for the second predictor, and so forth.
#The second row contains the half-slope ratios when the data are divided
#into two groups using the second predictor, and so on.
#
x<-as.matrix(x)
xmat<-matrix(0,nrow(x),ncol(x))
mval<-floor(length(y)/2)
mr<-length(y)-mval
xmatl<-matrix(0,mval,ncol(x))
xmatr<-matrix(0,mr,ncol(x))
hmat<-matrix(NA,ncol(x),ncol(x))
isub<-c(1:length(y))
ksub<-c(1:ncol(x))+1
for (k in 1:ncol(x)){
xord<-order(x[,k])
yord<-y[xord]
yl<-yord[isub<=mval]
yr<-yord[isub>mval]
for (j in 1:ncol(x)){
xmat[,j]<-x[xord,j]
xmatl[,j]<-xmat[isub<=mval,j]
xmatr[,j]<-xmat[isub>mval,j]
}
coefl<-regfun(xmatl,yl)$coef
coefr<-regfun(xmatr,yr)$coef
hmat[k,]<-coefr[ksub[ksub>=2]]/coefl[ksub[ksub>=2]]
}
hmat
}

# ============================================================================
# mbmreg
# ============================================================================
mbmreg<-function(x,y,iter=20,bend=2*sqrt(ncol(x)+1)/nrow(x),xout=FALSE,outfun=outpro,...){
#
# Compute a bounded M regression estimator using
# Huber Psi and Schweppe weights with
# regression outliers getting a weight of zero.
#
# This is the modified M-regression estimator in Chapter 8
#
# The predictors are assumed to be stored in the n by p matrix x.
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
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x<-as.matrix(x)
if(is.matrix(y)){
if(ncol(y)==1)y=as.vector(y)
}
x1<-cbind(1,x)
library(MASS)
reslms<-lmsreg(x,y)$resid
sighat<-sqrt(median(reslms^2))
sighat<-1.4826*(1+(5/(length(y)-ncol(x)-1)))*sighat
if(sighat==0)warning("The estimated measure of scale, based on the residuals using lms regression, is zero")
temp<-ifelse(sighat*reslms>0,abs(reslms)/sighat,0*reslms)
wt<-ifelse(temp<=2.5,1,0)
init<-lsfit(x,y,wt)
resid<-init$residuals
nu<-sqrt(1-hat(x1))
low<-ncol(x)+1
for(it in 1:iter){
ev<-sort(abs(resid))
scale<-median(ev[c(low:length(y))])/qnorm(.75)
rov<-(resid/scale)/nu
psi<-ifelse(abs(rov)<=bend,rov,bend*sign(rov))  # Huber Psi
wt<-nu*psi/(resid/scale)
wt<-ifelse(temp<=2.5,wt,0)
new<-lsfit(x,y,wt)
if(abs(max(new$coef-init$coef)<.0001))break
init$coef<-new$coef
resid<-new$residuals
}
resid<-y-x1%*%new$coef
if(abs(max(new$coef-init$coef)>=.0001))
paste("failed to converge in",iter,"steps")
list(coef=new$coef,residuals=resid,w=wt)
}

# ============================================================================
# regts1
# ============================================================================
regts1<-function(vstar,yhat,res,mflag,x,tr){
ystar<-yhat+res*vstar
bres<-ystar-mean(ystar,tr)
rval<-0
for (i in 1:nrow(x)){
rval[i]<-sum(bres[mflag[,i]])
}
rval
}

# ============================================================================
# depreg
# ============================================================================
depreg<-function(x,y,xout=FALSE,outfun=out,...){
#
# Compute the depth regression estimator.
# Only a single predictor is allowed in this version
#  Perhaps use instead
#
if(is.matrix(x)){
if(ncol(x)>=2)stop("Only a single predicor is allowed")
x<-as.vector(x)
}
xy=cbind(x,y)
xy=elimna(xy)
if(xout){
flag<-outfun(xy[,1],plotit=FALSE,...)$keep
xy<-xy[flag,]
}
x=xy[,1]
y=xy[,2]
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
slope<-v1/v2
vec3<-outer(ys,ys,"+")
vec4<-outer(xs,xs,"+")
v3<-vec3[vec2>0]
v4<-vec4[vec2>0]
deep<-NA
inter<-v3/2-slope*v4/2
temp<-matrix(c(inter,slope),ncol=2)
deep<-apply(temp,1,rdepth.orig,x,y)
best<-max(deep)
coef<-NA
coef[2]<-mean(slope[deep==best])
coef[1]<-mean(inter[deep==best])
res<-y-coef[2]*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# tsgreg
# ============================================================================
tsgreg<-function(x,y,tries=(length(y)^2-length(y))/2){
#
#
x<-as.matrix(x)
if(nrow(x)!=length(y))stop("Length of y must match the number of rows of x")
# eliminate any rows with missing values.
m1<-cbind(x,y)
m1<-elimna(m1)
x<-m1[,1:ncol(x)]
y<-m1[,ncol(x)+1]
set.seed(2)
data<-matrix(NA,ncol=ncol(x)+1,nrow=tries)
for(i in 1:tries){
data[i,]<-sample(length(y),size=ncol(x)+1,replace=FALSE)
}
bvec <- apply(data, 1,tsgregs1,x,y)
coef<-0
numzero<-0
loc<-0
for (i in 1:ncol(x)){
ip<-i+1
temp<-bvec[ip,]
loc[i]<-median(x[,i])
coef[i+1]<-median(temp[temp!=0])
numzero[i]<-length(temp[temp==0])
}
ip<-ncol(x)+1
coef[1]<-median(y)-sum(coef[2:ip]*loc)
res<-y-x %*% coef[2:ip] - coef[1]
list(coef=coef,residuals=res,numzero=numzero)
}

# ============================================================================
# tsgregs1
# ============================================================================
tsgregs1<-function(isub,x,y){
#
#  This function is used by tsgreg
#
#  Perform regression using x[isub,] to predict y[isub]
#  isub is a vector of length nsub, determined by tsgreg
#
tsgregs1<-lsfit(x[isub,],y[isub])$coef
}

# ============================================================================
# lts1reg
# ============================================================================
lts1reg<-function(x,y,tr=.2,h=NA){
#
# Compute the least trimmed squares regression estimator.
# Only a single predictor is allowed in this version
#
if(is.na(h))h<-length(x)-floor(tr * length(x))
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
slope<-v1/v2
vec3<-outer(ys,ys,"+")
vec4<-outer(xs,xs,"+")
v3<-vec3[vec2>0]
v4<-vec4[vec2>0]
val<-NA
inter<-v3/2-slope*v4/2
for(i in 1:length(slope)){
#risk<-(y[vec2>0]-slope[i]*x[vec2>0]-inter[i])^2
risk<-(y-slope[i]*x-inter[i])^2
risk<-sort(risk)
val[i]<-sum(risk[1:h])
}
best<-min(val)
coef<-NA
coef[2]<-mean(slope[val==best])
coef[1]<-mean(inter[val==best])
res<-y-coef[2]*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# twolsreg
# ============================================================================
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

# ============================================================================
# twolsregsub
# ============================================================================
twolsregsub<-function(isub, x, y)
{
        #
        #  Compute least squares estimate of the
        #  slope using x[isub] and y[isub]
        #  isub is a vector of length n,
        #  a bootstrap sample from the sequence of integers
        #  1, 2, 3, ..., n
        #
        twolsregsub<-lsfit(x[isub],y[isub])$coef[2]
        twolsregsub
}

# ============================================================================
# regi
# ============================================================================
regi<-function(x,y,z,pt=median(z),fr=.8,est=onestep,regfun=tsreg,testit=FALSE,...){
#
# split the data according to whether z is < or > pt, then
# use runmean2g to plot a smooth of the regression
# lines corresponding to these two groups.
#
m<-cbind(x,y,z)
m<-elimna(m)
x<-m[,1]
y<-m[,2]
z<-m[,3]
flag<-(z<pt)
runmean2g(x[flag],y[flag],x[!flag],y[!flag],fr=fr,est=est,...)
output<-"Done"
if(testit){
abline(regfun(x[flag],y[flag])$coef)
abline(regfun(x[!flag],y[!flag])$coef,lty=2)
output<-reg2ci(x[flag],y[flag],x[!flag],y[!flag],regfun=regfun,plotit=FALSE)
}
output
}

# ============================================================================
# linchk
# ============================================================================
linchk<-function(x,y,sp,pv=1,regfun=tsreg,plotit=TRUE,nboot=599,alpha=.05,pr=TRUE,xout=FALSE){
#
# Split the data into two groups according to whether
# predictor variable pv has a value less than sp.
# Then test the hypothesis that slope coefficients,
# based on the regression method regfun, are equal.
#
x<-as.matrix(x)
if(pr)print(paste("Splitting data using predictor", pv))
xx<-x[,pv]
flag<-(xx<=sp)
temp<-reg2ci(x[flag,],y[flag],x[!flag,],y[!flag],regfun=regfun,plotit=plotit,nboot=nboot,alpha=alpha,xout=xout)
temp
}

# ============================================================================
# lintests1
# ============================================================================
lintests1<-function(vstar,yhat,res,mflag,x,regfun,...){
ystar<-yhat+res*vstar
bres<-regfun(x,ystar,...)$residuals
rval<-0
for (i in 1:nrow(x)){
rval[i]<-sum(bres[mflag[,i]])
}
rval
}

# ============================================================================
# lsfitci
# ============================================================================
lsfitci<-function(x,y,nboot=599,alpha=.05,SEED=TRUE,xout=FALSE,outfun=out){
#
#   Compute a confidence interval for the slope parameters of
#   a linear regression equation when using the least squares estimator.
#
#   For p=1 predictor,
#   this function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
#   For p>1, a standard percentile bootstrap method is used
#   with FWE (the probability of at least one type I error)
#   controlled via the Bonferroni inequality.
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=599
#
#   SEED=T causes the seed of the random number generator to be set to 2,
#   otherwise the seed is not set.
#
#   Warning: probability coverage has been studied only when alpha=.05
#
x<-as.matrix(x)
p<-ncol(x)
pp<-p+1
temp<-elimna(cbind(x,y)) # Remove any missing values.
x<-temp[,1:p]
y<-temp[,p+1]
if(xout){
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,pp]
}
x<-as.matrix(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples; please wait")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,lsfit) # A p+1 by n matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
if(p==1){
if(alpha != .05){print("Resetting alpha to .05")
print("With p=1, unknown how to adjust confidence interval")
print("when alpha is not equal to .05.")
}
ilow<-15
ihi<-584
if(length(y) < 250){
ilow<-13
ihi<-586
}
if(length(y) < 180){
ilow<-10
ihi<-589
}
if(length(y) < 80){
ilow<-7
ihi<-592
}
if(length(y) < 40){
ilow<-6
ihi<-593
}
ilow<-round((ilow/599)*nboot)
ihi<-round((ihi/599)*nboot)
}
if(p>1){
ilow<-round(alpha*nboot/2)+1
ihi<-nboot-ilow
}
lsfitci<-matrix(0,ncol(x),2)
for(i in 1:ncol(x)){
ip<-i+1
bsort<-sort(bvec[ip,])
lsfitci[i,1]<-bsort[ilow+1]
lsfitci[i,2]<-bsort[ihi]
}
bsort<-sort(bvec[1,])
interceptci<-c(bsort[15],bsort[584])
crit.level<-NA
pmat<-NA
if(p>1){
crit.level<-alpha/p
pmat<-matrix(NA,nrow=p,ncol=2)
dimnames(pmat) <- list(NULL, c("Slope","p-value"))
for(pv in 1:p){
pmat[pv,1]<-pv
pp<-pv+1
pmat[pv,2]<-(sum(bvec[pp,]<0)+.5*sum(bvec[pp,]==0))/nboot
temp3<-1-pmat[pv,2]
pmat[pv,2]<-2*min(pmat[pv,2],temp3)
}}
list(intercept.ci=interceptci,slope.ci=lsfitci,crit.level=crit.level,
p.values=pmat)
}

# ============================================================================
# lsfitNci
# ============================================================================
lsfitNci<-function(x,y,alpha=.05){
#
# Compute confidence interval for least squares
# regression using heteroscedastic method
# recommended by Long and Ervin (2000).
#
x<-as.matrix(x)
if(nrow(x) != length(y))stop("Length of y does not match number of x values")
m<-cbind(x,y)
m<-elimna(m)
y<-m[,ncol(x)+1]
temp<-lsfit(x,y)
x<-cbind(rep(1,nrow(x)),m[,1:ncol(x)])
xtx<-solve(t(x)%*%x)
h<-diag(x%*%xtx%*%t(x))
hc3<-xtx%*%t(x)%*%diag(temp$res^2/(1-h)^2)%*%x%*%xtx
df<-nrow(x)-ncol(x)
crit<-qt(1-alpha/2,df)
al<-ncol(x)
ci<-matrix(NA,nrow=al,ncol=3)
for(j in 1:al){
ci[j,1]<-j
ci[j,2]<-temp$coef[j]-crit*sqrt(hc3[j,j])
ci[j,3]<-temp$coef[j]+crit*sqrt(hc3[j,j])
}
print("Confidence intervals for intercept followed by slopes:")
list(ci=ci,stand.errors=sqrt(diag(hc3)))
}

# ============================================================================
# taureg
# ============================================================================
taureg<-function(m,y,corfun=tau,...){
#
#    Compute Kendall's tau between y and each of the
#    p variables stored  in the n by p matrix m.
#
#    Alternative measures of correlation can be used via the
#    argument corfun. The only requirement is that the function
#    corfun returns the correlation in corfun$cor and the p-value
#    in corfun$p.value.
#
#    This function also returns the two-sided significance level
#    for all pairs of variables, plus a test of zero correlations
#    among all pairs. (See chapter 9 of Wilcox, 2005, for details.)
#
m<-as.matrix(m)
tauvec<-NA
siglevel<-NA
for (i in 1:ncol(m)){
pbc<-corfun(m[,i],y,...)
tauvec[i]<-pbc$cor
siglevel[i]<-pbc$p.value
}
list(cor=tauvec,p.value=siglevel)
}

# ============================================================================
# correg.sub
# ============================================================================
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

# ============================================================================
# correg
# ============================================================================
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

# ============================================================================
# lintest
# ============================================================================
lintest<-function(x,y,regfun=tsreg,nboot=500,alpha=.05,xout=FALSE,SEED=TRUE,
outfun=out,...){
#
# Test the hypothesis that the regression surface is a plane.
# Stute et al. (1998, JASA, 93, 141-149).
#
if(SEED)set.seed(2)
#if(identical(regfun,Qreg))print('When using Qreg, be sure to include res.vals=TRUE')
#if(identical(regfun,tshdreg))print('When using tshdreg, be sure to include RES=TRUE')
#if(identical(regfun,MMreg))print('When using MMreg, be sure to include RES=TRUE') # no longer necessary
x<-as.matrix(x)
d<-ncol(x)
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
x<-as.matrix(x)
y<-temp[,d+1]
if(xout){
flag<-outfun(x,...)$keep
x<-x[flag,]
x<-as.matrix(x)
y<-y[flag]
}
mflag<-matrix(NA,nrow=length(y),ncol=length(y))
for (j in 1:length(y)){
for (k in 1:length(y)){
mflag[j,k]<-(sum(x[j,]<=x[k,])==ncol(x))
}
}
reg<-regfun(x,y,...)
yhat<-y-reg$residuals
#print("Taking bootstrap samples, please wait.")
data<-matrix(runif(length(y)*nboot),nrow=nboot)
data<-sqrt(12)*(data-.5) # standardize the random numbers.
rvalb<-apply(data,1,lintests1,yhat,reg$residuals,mflag,x,regfun,...)
# An n x nboot matrix of R values
rvalb<-rvalb/sqrt(length(y))
dstatb<-apply(abs(rvalb),2,max)
wstatb<-apply(rvalb^2,2,mean)
# compute test statistic
v<-c(rep(1,length(y)))
rval<-lintests1(v,yhat,reg$residuals,mflag,x,regfun,...)
rval<-rval/sqrt(length(y))
dstat<-max(abs(rval))
wstat<-mean(rval^2)
ib<-round(nboot*(1-alpha))
p.value.d<-1-sum(dstat>=dstatb)/nboot
p.value.w<-1-sum(wstat>=wstatb)/nboot
list(dstat=dstat,wstat=wstat,p.value.d=p.value.d,p.value.w=p.value.w)
}

# ============================================================================
# regtest
# ============================================================================
regtest<-function(x,y,regfun=tsreg,nboot=600,alpha=.05,plotit=TRUE,
grp=c(1:ncol(x)),nullvec=c(rep(0,length(grp))),xout=FALSE,outfun=outpro,SEED=TRUE,pr=TRUE,...){
#
#  Test the hypothesis that q of the p predictors are equal to
#  some specified constants. By default, the hypothesis is that all
#  p predictors have a coefficient equal to zero.
#  The method is based on a confidence ellipsoid.
#  The critical value is determined with the percentile bootstrap method
#  in conjunction with Mahalanobis distance.
#
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
if(pr)print("Default for outfun is now outpro")
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,regfun=regfun,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x<-as.matrix(x)
if(length(grp)!=length(nullvec))stop("The arguments grp and nullvec must have the same length.")
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun) # A p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
grp<-grp+1  #Ignore the intercept.
est<-regfun(x,y)$coef
estsub<-est[grp]
bsub<-t(bvec[grp,])
if(length(grp)==1){
m1<-sum((bvec[grp,]-est)^2)/(length(y)-1)
dis<-(bsub-estsub)^2/m1
}
if(length(grp)>1){
mvec<-apply(bsub,2,FUN=mean)
m1<-var(t(t(bsub)-mvec+estsub))
dis<-mahalanobis(bsub,estsub,m1)
}
dis2<-order(dis)
dis<-sort(dis)
critn<-floor((1-alpha)*nboot)
crit<-dis[critn]
test<-mahalanobis(t(estsub),nullvec,m1)
sig.level<-1-sum(test>dis)/nboot
if(length(grp)==2 && plotit){
plot(bsub,xlab="Parameter 1",ylab="Parameter 2")
points(nullvec[1],nullvec[2],pch=0)
xx<-bsub[dis2[1:critn],]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(test=test,crit=crit,p.value=sig.level,nullvec=nullvec,est=estsub,n=length(y))
}

# ============================================================================
# reg2ci
# ============================================================================
reg2ci<-function(x,y,x1,y1,regfun=tsreg,nboot=599,alpha=.05,plotit=TRUE,SEED=TRUE,
xout=FALSE,outfun=outpro,xlab="X",ylab="Y",pr=FALSE,...){
#
#   Compute a .95 confidence interval for the difference between the
#   the intercepts and slopes corresponding to two independent groups.
#   The default regression method is Theil-Sen.
#
#   The predictor values for the first group are
#   assumed to be in the n by p matrix x.
#   The predictors for the second group are in x1
#
#   The default number of bootstrap samples is nboot=599
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
x1<-as.matrix(x1)
xx1<-cbind(x1,y1)
xx1<-elimna(xx1)
x1<-xx1[,1:ncol(x1)]
x1<-as.matrix(x1)
y1<-xx1[,ncol(x1)+1]
x=as.matrix(x)
x1=as.matrix(x1)
if(xout){
if(pr)print("outfun now defaults to outpro rather than out")
if(identical(outfun,outblp)){
flag1=outblp(x,y,plotit=FALSE)$keep
flag2=outblp(x1,y2,plotit=FALSE)$keep
}
if(!identical(outfun,outblp)){
flag1=outfun(x,plotit=FALSE)$keep
flag2=outfun(x1,plotit=FALSE)$keep
}
x=x[flag1,]
y=y[flag1]
x1=x1[flag2,]
y1=y1[flag2]
}
n=length(y)
n[2]=length(y1)
x<-as.matrix(x)
x1<-as.matrix(x1)
est1=regfun(x,y,...)$coef
est2=regfun(x1,y1,...)$coef
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun,xout=FALSE,...) # A p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data,1,regboot,x1,y1,regfun,xout=FALSE,...)
bvec<-bvec-bvec1
p1<-ncol(x)+1
regci<-matrix(0,p1,6)
dimnames(regci)<-list(NULL,
c("Parameter","ci.lower","ci.upper","p.value","Group 1","Group 2"))
ilow<-round((alpha/2)*nboot)+1
ihi<-nboot-(ilow-1)
for(i in 1:p1){
temp<-sum(bvec[i,]<0)/nboot+sum(bvec[i,]==0)/(2*nboot)
regci[i,4]<-2*min(temp,1-temp)
bsort<-sort(bvec[i,])
regci[i,2]<-bsort[ilow]
regci[i,3]<-bsort[ihi]
regci[,1]<-c(0:ncol(x))
}
regci[,5]=est1
regci[,6]=est2
if(ncol(x)==1 && plotit){
plot(c(x,x1),c(y,y1),type="n",xlab=xlab,ylab=ylab)
points(x,y)
points(x1,y1,pch="+")
abline(regfun(x,y,...)$coef)
abline(regfun(x1,y1,...)$coef,lty=2)
}
list(n=n,output=regci)
}

# ============================================================================
# regpre
# ============================================================================
regpre<-function(x,y,regfun=lsfit,error=absfun,nboot=100,adz=TRUE,
mval=round(5*log(length(y))),model=NULL,locfun=mean,pr=FALSE,
xout=FALSE,outfun=out,STAND=TRUE,
plotit=TRUE,xlab="Model Number",ylab="Prediction Error",SEED=TRUE,...){
#
#   Estimate prediction error using the regression method
#   regfun. The .632 method is used.
#   (See Efron and Tibshirani, 1993, pp. 252--254)
#
#   The predictor values are assumed to be in the n-by-p matrix x.
#   The default number of bootstrap samples is nboot=100
#
#   Prediction error is the expected value of the function error.
#   The argument error defaults to squared error.
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
#   The default value for mval, the number of observations to resample
#   for each of the B bootstrap samples is based on results by
#   Shao (JASA, 1996, 655-665). (Resampling n vectors of observations
#   model selection may not lead to the correct model as n->infinity.
#
#   The argument model should have list mode, model[[1]] indicates
#   which predictors are used in the first model. For example, storing
#   1,4 in model[[1]] means predictors 1 and 4 are being considered.
#   If model is not specified, and number of predictors is at most 5,
#   then all models are considered.
#
#   If adz=T, added to the models to be considered is where
#   all regression slopes are zero. That is, use measure of location only
#   corresponding to
#   locfun.
#
if(pr){
print("By default, least squares regression is used, ")
print("But from Wilcox, R. R. 2008, Journal of Applied Statistics, 35, 1-8")
print("Setting regfun=tsreg appears to be a better choice for general use.")
print("That is, replace least squares with the Theil-Sen estimator")
print("Note: Default for the argument error is now absfun")
print(" meaning absolute error is used")
print("To use squared error, set error=sqfun")
}
x<-as.matrix(x)
d<-ncol(x)
p1<-d+1
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
y<-temp[,d+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(!STAND)flag<-outfun(x,plotit=FALSE,...)$keep
if(STAND)flag<-outpro(x,STAND=TRUE,plotit=FALSE)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(is.null(model)){
if(d<=5)model<-modgen(d,adz=adz)
if(d>5)model[[1]]<-c(1:ncol(x))
}
mout<-matrix(NA,length(model),5,dimnames=list(NULL,c("apparent.error",
"boot.est","err.632","var.used","rank")))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=mval*nboot,replace=TRUE),nrow=nboot)
bid<-apply(data,1,idb,length(y))
#  bid is an n by nboot matrix. If the jth bootstrap sample from
#  1, ..., mval contains the value i, bid[i,j]=0; otherwise bid[i,j]=1
for (imod in 1:length(model)){
nmod=length(model[[imod]])-1
temp=c(nmod:0)
mout[imod,4]=sum(model[[imod]]*10^temp)
if(sum(model[[imod]]==0)!=1){
xx<-x[,model[[imod]]]
xx<-as.matrix(xx)
if(sum(model[[imod]]==0)!=1)bvec<-apply(data,1,regpres1,xx,y,regfun,mval,...)
# bvec is a p+1 by nboot matrix. The first row
# contains the bootstrap intercepts, the second row
# contains the bootstrap values for first predictor, etc.
if(sum(model[[imod]]==0)!=1)yhat<-cbind(1,xx)%*%bvec
if(sum(model[[imod]]==0)==1){
bvec0<-matrix(0,nrow=p1,ncol=nboot)
for(it in 1:nboot){
bvec0[1,it]<-locfun(y[data[it,]])
}
yhat<-cbind(1,x)%*%bvec0
}
# yhat is n by nboot matrix of predicted values based on
                           # bootstrap regressions.
bi<-apply(bid,1,sum) # B sub i in notation of Efron and Tibshirani, p. 253
temp<-(bid*(yhat-y))
diff<-apply(temp,1,error)
ep0<-sum(diff/bi)/length(y)
aperror<-error(regfun(xx,y,...)$resid)/length(y) # apparent error
regpre<-.368*aperror+.632*ep0
mout[imod,1]<-aperror
mout[imod,3]<-regpre
temp<-yhat-y
diff<-apply(temp,1,error)
mout[imod,2]<-sum(diff)/(nboot*length(y))
}
if(sum(model[[imod]]==0)==1){
mout[imod,3]<-locpre(y,error=error,est=locfun,SEED=SEED,mval=mval)
}}
mout[,5]=rank(mout[,3])
if(plotit)plot(c(1:nrow(mout)),mout[,3],xlab=xlab,ylab=ylab)
list(estimates=mout)
}

# ============================================================================
# regout
# ============================================================================
regout<-function(x,y,regest=stsreg,plotit=TRUE,mbox=TRUE){
#
# Check for regression outliers by fitting a
# a line to data using regest and then applying
# a boxplot rule to the residuals.
# mbox=T uses Carling's method
# mbox=F uses ideal fourths with conventional boxplot rules.
#
chk<-regest(x,y)
flag<-outbox(chk$residuals,mbox=mbox)$out.id
if(plotit){
plot(x,y)
points(x[flag],y[flag],pch="o")
abline(chk$coef)
}
list(out.id=flag)
}

# ============================================================================
# stsregp1
# ============================================================================
stsregp1<-function(x,y,sc=pbvar,xout=FALSE,outfun=out,...){
#
# Compute the S-type modification of
# the Theil-Sen regression estimator.
# Only a single predictor is allowed in this version
#
xy=elimna(cbind(x,y))
p=ncol(as.matrix(x))
if(p!=1)stop("Current version is limited to one predictor")
p1=p+1
x=xy[,1:p]
y=xy[,p1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
slope<-v1/v2
allvar<-NA
for(i in 1:length(slope))allvar[i]<-sc(y-slope[i]*x,...)
temp<-order(allvar)
coef<-0
coef[2]<-slope[temp[1]]
coef[1]<-median(y)-coef[2]*median(x)
res<-y-coef[2]*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# stsreg
# ============================================================================
stsreg<-function(x,y,xout=FALSE,outfun=outpro,iter=10,sc=pbvar,varfun=pbvar,
corfun=pbcor,plotit=FALSE,...){
#
#  Compute Theil-Sen regression estimator
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(ncol(x)==1){
temp1<-stsregp1(x,y,sc=sc)
coef<-temp1$coef
res<-temp1$res
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-tsp1reg(x[,p],y)$coef[2]
}
res<-y-x%*%temp
alpha<-median(res)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-stsregp1(x[,p],r[,p],sc=sc)$coef[2]
}
alpha<-median(y-x%*%temp)
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
yhat<-y-res
stre=NULL
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}
list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

# ============================================================================
# regbootg
# ============================================================================
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

# ============================================================================
# rregci
# ============================================================================
rregci<-function(x,y,regfun=chreg,nboot=599,alpha=.05, ...){
#
#   Compute a .95 confidence interval for each of the parameters of
#   a linear regression equation. The default regression method is
#   a bounded influence M-regression with Schweppe weights
#   (the Coakley-Hettmansperger estimator).
#
#   When using the least squares estimator, and when n<250, use
#   lsfitci instead.
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=599
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimated of
#   the first predictor, etc.
#
x<-as.matrix(x)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun,...) # A p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
p1<-ncol(x)+1
regci<-matrix(0,p1,2)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
se<-NA
pvec<-NA
for(i in 1:p1){
bsort<-sort(bvec[i,])
pvec[i]<-sum(bvec[i,]<0)/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
regci[i,1]<-bsort[ilow]
regci[i,2]<-bsort[ihi]
se[i]<-sqrt(var(bvec[i,]))
}
pvec<-2*pvec
list(regci=regci,p.value=pvec,se=se)
}

# ============================================================================
# wsp1reg
# ============================================================================
wsp1reg<-function(x,y,plotit=FALSE){
#
# Compute the Wilcoxon R estimate of the slope
# Only a single predictor is allowed in this version
#
temp<-matrix(c(x,y),ncol=2)
temp<-elimna(temp)     # Remove any pairs with missing values
x<-temp[,1]
y<-temp[,2]
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
slope<-v1/v2
tmin<-wrregfun(slope[1],x,y)
ikeep<-1
for(i in 2:length(slope)){
tryit<-wrregfun(slope[i],x,y)
if(tryit<tmin){
tmin<-tryit
ikeep<-i
}}
coef<-NA
coef[2]<-slope[ikeep]
coef[1]<-median(y-coef[2]*x)
if(plotit){
plot(x,y,xlab="X",ylab="Y")
abline(coef)
}
res<-y-coef[2]*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# wreg
# ============================================================================
wreg<-function(x,y,iter=10){
#
#  Compute Wilcoxon R estimate
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#
#  Argument iter is used when there is more than one
#  predictor and indicates maximum number
#  of iterations used by Gauss-Seidel algoritm.
#
temp<-NA
x<-as.matrix(x)
if(ncol(x)==1){
temp1<-wsp1reg(x,y)
coef<-temp1$coef
res<-temp1$res
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-wsp1reg(x[,p],y)$coef[2]
}
res<-y-x%*%temp
alpha<-median(res)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-wsp1reg(x[,p],r[,p],plotit=FALSE)$coef[2]
}
alpha<-median(y-x%*%temp)
if(max(abs(tempold-temp))<.0001)break
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
list(coef=coef,residuals=res)
}

# ============================================================================
# opreg
# ============================================================================
opreg<-function(x,y,regfun=tsreg,cop=3,MC=FALSE,varfun=pbvar,corfun=pbcor,STAND=TRUE,xout=FALSE){
#
# Do regression on points not labled outliers
# using projection-type outlier detection method
#
#  Note: argument xout is not relevant here, but is included to avoid conflicts when using regci.
#
if(MC)library(parallel)
x<-as.matrix(x)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
if(!MC)ivec<-outpro(m,plotit=FALSE,cop=cop,STAND=STAND)$keep
if(MC)ivec<-outproMC(m,plotit=FALSE,cop=cop,STAND=STAND)$keep
np1<-ncol(x)+1
coef<-regfun(m[ivec,1:ncol(x)],m[ivec,np1])$coef
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
# wrregfun
# ============================================================================
wrregfun<-function(slope,x=x,y=y){
x<-as.matrix(x)
res<-y-x%*%slope
v1<-rank(res)
v2<-sqrt(12)*(v1/(length(y)+1)-.5)
wrregfun<-sum(v2*res)
wrregfun
}

# ============================================================================
# opregpb
# ============================================================================
opregpb<-function(x,y,nboot=1000,alpha=.05,om=TRUE,ADJ=TRUE,SEED=TRUE,
nullvec=rep(0,ncol(x)+1),plotit=TRUE,opdis=2,gval=sqrt(qchisq(.95,ncol(x)+1))){
#
# generate bootstrap estimates
# use projection-type outlier detection method followed by
# TS regression.
#
# om=T and ncol(x)>1, means an omnibus test is performed,
# otherwise only individual tests of parameters are performed.
#
# opdis=2, means that Mahalanobis distance is used
# opdis=1, means projection-type distance is used
#
# gval is critical value for projection-type outlier detection
# method
#
# ADJ=T, Adjust p-values as described in Section 11.1.5 of the text.
#
if(SEED)set.seed(2)
x<-as.matrix(x)
m<-cbind(x,y)
p1<-ncol(x)+1
m<-elimna(m) # eliminate any rows with missing data
x<-m[,1:ncol(x)]
x<-as.matrix(x)
y<-m[,p1]
if(nrow(x)!=length(y))stop("Sample size of x differs from sample size of y")
if(!is.matrix(x))stop("Data should be stored in a matrix")
print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun=opreg)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
# using Hochberg method
bvec<-t(bvec)
dvec<-alpha/(c(1:ncol(x)))
test<-NA
icl0<-round(alpha*nboot/2)
icl<-round(alpha*nboot/(2*ncol(x)))
icu0<-nboot-icl0
icu<-nboot-icl
output<-matrix(0,p1,6)
vlabs="Intercept"
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(output)<-list(vlabs,c("Param.","p.value","p.crit",
"ci.lower","ci.upper","s.e."))
pval<-NA
for(i in 1:p1){
output[i,1]<-i-1
se.val<-var(bvec[,i])
temp<-sort(bvec[,i])
output[i,6]<-sqrt(se.val)
if(i==1){
output[i,4]<-temp[icl0+1]
output[i,5]<-temp[icu0]
}
if(i>1){
output[i,4]<-temp[icl+1]
output[i,5]<-temp[icu]
}
pval[i]<-sum((temp>nullvec[i]))/length(temp)
if(pval[i]>.5)pval[i]<-1-pval[i]
}
fac<-2
if(ADJ){
# Adjust p-value if n<60
nval<-length(y)
if(nval<20)nval<-20
if(nval>60)nval<-60
fac<-2-(60-nval)/40
}
pval[1]<-2*pval[1]
pval[2:p1]<-fac*pval[2:p1]
output[,2]<-pval
temp2<-order(0-pval[2:p1])
zvec<-dvec[1:ncol(x)]
sigvec<-(test[temp2]>=zvec)
output[temp2+1,3]<-zvec
output[1,3]<-NA
output[,2]<-pval
om.pval<-NA
temp<-opreg(x,y)$coef
if(om && ncol(x)>1){
temp2<-rbind(bvec[,2:p1],nullvec[2:p1])
if(opdis==1)dis<-pdis(temp2,center=temp[2:p1])
if(opdis==2){
cmat<-var(bvec[,2:p1]-apply(bvec[,2:p1],2,mean)+temp[2:p1])
dis<-mahalanobis(temp2,temp[2:p1],cmat)
}
om.pval<-sum((dis[nboot+1]<=dis[1:nboot]))/nboot
}
# do adjusted p-value
nval<-length(y)
if(nval<20)nval<-20
if(nval>60)nval<-60
adj.pval<-om.pval/2+(om.pval-om.pval/2)*(nval-20)/40
if(ncol(x)==2 && plotit){
plot(bvec[,2],bvec[,3],xlab="Slope 1",ylab="Slope 2")
temp.dis<-order(dis[1:nboot])
ic<-round((1-alpha)*nboot)
xx<-bvec[temp.dis[1:ic],2:3]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(output=output,om.pval=om.pval,adj.om.pval=adj.pval)
}

# ============================================================================
# ltsgreg
# ============================================================================
ltsgreg<-function(x, y, tr = 0.2, h = NA,xout=FALSE,outfun=outpro,...)
{
        #
        # Compute the least trimmed absolute value regression estimator.
        # The default amount of trimming is .2
        #
        #  Can simply use ltsreg with tr=amount of trimming.
        #
x<-as.matrix(x)
library(MASS)
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
p<-np-1
x<-X[,1:p]
x<-as.matrix(x)
y<-X[,np]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
X<-cbind(x,y)
}
if(is.na(h)) h <- length(y) - floor(tr * length(y))
temp<-ltsReg(y~x)$coef
START<-temp
coef<-nelderv2(X,np,FN=lts.sub,START=START,h=h,p=p)
        res <- y - x%*%coef[2:np] - coef[1]
        list(coef = coef, residuals = res)
}

# ============================================================================
# locreg
# ============================================================================
locreg<-function(x,y,pyhat=FALSE,pts=NA,np=100,plotit=TRUE,eout=FALSE,outfun=out,
xlab="X",ylab="Y",pch='.'){
#
# Compute local weighted regression with Epanechnikov kernel
#
# See Fan, Annals of Statistics, 1993, 21, 196-217.
# cf. Bjerve and Doksum, Annals of Statistics, 1993, 21, 890-902
#
# With np=100, the function plots a smooth using
# middle 80% of the x values versus y
# With np=0, it plots using all x values
# or all values in pts if values are stored in it.
# With np=0, pts=x is used.
#
# pyhat=T, the function returns the estimated y values
# corresponding to x values in pts. If pts=NA, pts=x
# is assumed.
#
m<-elimna(cbind(x,y))
if(eout){
keep<-outfun(m,plotit=FALSE)$keep
m<-m[keep,]
}
x<-m[,1]
y<-m[,2]
n<-length(x)
sig<-sqrt(var(x))
temp<-idealf(x)
iqr<-(temp$qu-temp$ql)/1.34
A<-min(c(sig,iqr))
yhat<-NA
temp<-NA
if(is.na(pts[1])){
if(np>0)pts<-seq(min(x),max(x),length=np)
if(np==0)pts<-x
}
pts<-sort(pts)
for(i in 1:length(pts)){
yhat[i]<-NA
for(j in 1:length(x)){
temp[j]<-((x[j]-pts[i])/A)^2
}
epan<-ifelse(temp<1,.75*(1-temp),0)
chkit<-sum(epan!=0)
if(chkit > 1){
vals<-lsfit(x,y,wt=epan)$coef
yhat[i]<-vals[2]*pts[i]+vals[1]
}
}
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab,pch=pch)
if(np>0){
ilow<-round(.1*np)
iup<-round(.9*np)
}
if(np==0){
ilow<-1
iup<-length(pts)
}
lines(pts[ilow:iup],yhat[ilow:iup])
}
m<-"Done"
if(pyhat)m<-yhat
m
}

# ============================================================================
# qreg.sub
# ============================================================================
qreg.sub<-function(X,theta,qval=.5){
np<-ncol(X)
p<-np-1
x<-X[,1:p]
y<-X[,np]
temp<-t(t(x)*theta[2:np])
yhat<-apply(temp,1,sum)+theta[1]
res<-y-yhat
flag<-(res<=0)
rval<-(qval-flag)*res
val<-sum(rval)
val
}

# ============================================================================
# qreg
# ============================================================================
qreg<-function(x, y,qval=.5, q=NULL,pr=FALSE,xout=FALSE, outfun=outpro,plotit=FALSE,xlab="X",ylab="Y",op=1,v2=TRUE,method='br',WARN=FALSE,...)
{
#
# Compute the quantile regression line. That is, the goal is to
# determine the qth (qval) quantile of Y given X using the
#  the Koenker-Bassett approach.
#
#  v2=T, uses the function rq in the R library quantreg
#  v2=F, uses an older and slower version
#  op=1 has to do with the old version.
#
#  method=scad, see Wu and  Liu (2009). VARIABLE SELECTION IN QUANTILE REGRESSION, Statistica Sinica 19, 801-817.
#
if(!is.null(q))qval=q
x<-as.matrix(x)
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
p<-np-1
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
if(!v2){
temp<-ltareg(x,y,0,op=op)
if(qval==.5){
coef<-temp$coef
res<-temp$res
}
if(qval!=.5){
START<-temp$coef
coef<-nelderv2(X,np,FN=qreg.sub,START=START,qval=qval)
}}
if(v2){
library(quantreg)
x<-as.matrix(x)
if(!WARN)options(warn=-1)
temp<-rq(y~x,tau=qval,method=method)
coef<-temp[1]$coefficients
if(!WARN)options(warn=0)
}
if(ncol(x)==1){
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab)
abline(coef)
}}
res <- y - x%*%coef[2:np] - coef[1]
list(coef = coef, residuals = res)
}

# ============================================================================
# kerreg
# ============================================================================
kerreg<-function(x,y,pyhat=FALSE,pts=NA,plotit=TRUE,theta=50,phi=25,expand=.5,
scale=FALSE,zscale=FALSE,eout=FALSE,xout=FALSE,outfun=out,np=100,xlab="X",ylab="Y",zlab="Z",
varfun=pbvar,e.pow=TRUE,pr=TRUE,ticktype="simple",pch='.',...){
#
# Compute local weighted regression with Epanechnikov kernel
#
# See Fan, Annals of Statistics, 1993, 21, 196-217.
# cf. Bjerve and Doksum, Annals of Statistics, 1993, 21, 890-902
#
# With a single predictor, this function calls locreg
# See locreg for information about np and plotting
#
library(akima)
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
d<-ncol(x)
np1<-d+1
m<-elimna(cbind(x,y))
if(xout && eout)stop("Can't have eout=xout=T")
if(eout){
flag<-outfun(m,plotit=FALSE,...)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
}
if(zscale){
for(j in 1:np1){
m[,j]<-(m[,j]-median(m[,j]))/mad(m[,j])
}}
x<-m[,1:d]
x<-as.matrix(x)
y<-m[,np1]
n<-nrow(x)
if(d>1){
xrem<-x
pi<-gamma(.5)^2
cd<-c(2,pi)
if(d==2)A<-1.77
if(d==3)A<-2.78
if(d>2){
for(j in 3:d)cd[j]<-2*pi*cd[j-2]/j  # p. 76
}
if(d>3)A<-(8*d*(d+2)*(d+4)*(2*sqrt(pi))^d)/((2*d+1)*cd[d])  # p. 87
hval<-A*(1/n)^(1/(d+4))  # p. 86
for(j in 1:d){
sig<-sqrt(var(x[,j]))
temp<-idealf(x[,j])
iqr<-(temp$qu-temp$ql)/1.34
A<-min(c(sig,iqr))
x[,j]<-x[,j]/A
}
xx<-cbind(rep(1,nrow(x)),x)
yhat<-NA
for(j in 1:n){
yhat[j]<-NA
temp1<-t(t(x)-x[j,])/(hval)
temp1<-temp1^2
temp1<-apply(temp1,1,FUN="sum")
temp<-.5*(d+2)*(1-temp1)/cd[d]
epan<-ifelse(temp1<1,temp,0) # Epanechnikov kernel, p. 76
chkit<-sum(epan!=0)
if(chkit >= np1){
vals<-lsfit(x,y,wt=epan)$coef
yhat[j]<-xx[j,]%*%vals
}}
if(plotit  && d==2){
if(pr){
if(!scale){
print("scale=F is specified")
print("If there is dependence, might use scale=T")
}}
m<-elimna(cbind(xrem,yhat))
xrem<-m[,1:d]
yhat<-m[,np1]
fitr<-yhat
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
}}
if(d==1){
yhat<-locreg(x[,1],y,pyhat=TRUE,np=np,plotit=plotit,pts=pts,
xlab=xlab,ylab=ylab,pch=pch)
yhat2<-locreg(x[,1],y,pyhat=TRUE,np=0,plotit=FALSE)
}
if(d>1)yhat2<-yhat
m<-NULL
#E.pow<-varfun(yhat2[!is.na(yhat2)])/varfun(y)
# Estimate of explanatory power performs poorly.
if(pyhat)m<-yhat
#list(Strength.Assoc=sqrt(E.pow),Explanatory.Power=E.pow,yhat=m)
m
}

# ============================================================================
# snmreg.sub
# ============================================================================
snmreg.sub<-function(X,theta){
np<-ncol(X)
p<-np-1
x<-X[,1:p]
y<-X[,np]
temp<-t(t(x)*theta)
yhat<-apply(temp,1,sum)
yhat<-yhat
res<-y-yhat
val<-pbvar(res)
val
}

# ============================================================================
# tstsreg
# ============================================================================
tstsreg<-function(x,y,sc=pbvar,xout=FALSE,outfun=outpro,iter=5,plotit=FALSE,...){
#
# Compute a modified Theil-Sen regression estimator.
# Use s-type initial estimate, eliminate points with
# outlying residuals, then do regular Theil-Sen
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
res=stsreg(x,y,sc=sc,iter=iter)$res
chk<-abs(res-median(res))/mad(res)
xx<-x[chk<=2,]
yy<-y[chk<=2]
temp<-tsreg(xx,yy)
list(coef=temp$coef,residuals=temp$res)
}

# ============================================================================
# tssnmreg
# ============================================================================
tssnmreg<-function(x,y,sc=pbvar,xout=FALSE,outfun=out,plotit=FALSE,...){
#
# Compute a modified Theil--Sen regression estimator.
# Use s-type initial estimate, eliminate points with
# outlying residuals, then do regular Theil--Sen
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
res=snmreg(x,y)$res
chk<-abs(res-median(res))/mad(res)
xx<-x[chk<=2,]
yy<-y[chk<=2]
temp<-tsreg(xx,yy)
list(coef=temp$coef,residuals=temp$res)
}

# ============================================================================
# gyreg
# ============================================================================
gyreg<-function(x,y,rinit=lmsreg,K=2.5){
xy=elimna(cbind(x,y))
p=ncol(as.matrix(x))
p1=p+1
x=xy[,1:p]
y=xy[,p1]
library(MASS)
res<-rinit(y~x)$res
res.scale<-abs(res)/mad(res)
flag<-(res.scale >=K)
i0<-sum(flag)
il<-length(y)-i0+1
res.sort<-sort(res.scale)
if(i0>0){
dval<-pnorm(res.sort[il:length(y)])-c(il:length(y))/length(y)
}
if(i0<=0)dval<-0
dval<-max(dval)
ndval<-floor(length(y)*dval)
if(ndval<0)ndval<-0
iup<-length(y)-ndval
rord<-order(res.scale)
flag<-rord[1:iup]
x=as.matrix(x)
temp<-lsfit(x[flag,],y[flag])
list(coef=temp$coef,res=temp$residual)
}

# ============================================================================
# regpreS
# ============================================================================
regpreS<-function(x,y,regfun=lsfit,error=absfun,nboot=100,
mval=round(5*log(length(y))),locfun=mean,pr=TRUE,
xout=FALSE,outfun=out,
plotit=TRUE,xlab="Model Number",ylab="Prediction Error",SEED=TRUE,...){
#
#   Stepwise selection of predictors based on
#   estimates of  prediction error using the regression method
#   regfun,
#   which defaults to least squares.  Prediction error
#   is estimated with .632 method.
#   (See Efron and Tibshirani, 1993, pp. 252--254)
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=100
#
#   Prediction error is the expected value of the function error.
#   The argument error defaults to absolute  error. To use
#   squared error, set error=sqfun.
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
#   The default value for mval, the number of observations to resample
#   for each of the B bootstrap samples is based on results by
#   Shao (JASA, 1996, 655-665). (Resampling n vectors of observations,
#   model selection may not lead to the correct model as n->infinity.
#
if(SEED)set.seed(2)
q=ncol(x)
qm1=q-1
x<-as.matrix(x)
d<-ncol(x)
p1<-d+1
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
y<-temp[,d+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,SEED=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
adit=NULL
pval=c(1:ncol(x))
#pval=c(1:q)
allp=pval
for(ip in 1:qm1){
model=list()
for(j  in 1:length(pval))model[[j]]=c(adit,pval[j])
temp=regpre(x,y,model=model,pr=FALSE,plotit=FALSE,adz=FALSE,regfun=regfun,
SEED=SEED)$estimates
pbest=order(temp[,5])
adit=model[[pbest[1]]]
pval=allp[-adit]
}
output=model[[pbest[1]]]
output=c(output,allp[-output])
output
}

# ============================================================================
# tsp1reg
# ============================================================================
tsp1reg<-function(x,y,plotit=FALSE,HD=FALSE,OPT=TRUE,tr=FALSE){
#
# Compute the Theil-Sen regression estimator.
# Only a single predictor is allowed in this version
#
#  OPT=TRUE, compute the intercept using median(y)-beta_1median(X)
#  OPT=FALSE compute the intercept using median of y-beta_1X
#
temp<-matrix(c(x,y),ncol=2)
temp<-elimna(temp)     # Remove any pairs with missing values
x<-temp[,1]
y<-temp[,2]
ord<-order(x)
xs<-x[ord]
ys<-y[ord]
vec1<-outer(ys,ys,"-")
vec2<-outer(xs,xs,"-")
v1<-vec1[vec2>0]
v2<-vec2[vec2>0]
if(!HD)slope<-median(v1/v2,na.rm=TRUE)
if(HD)slope<-hd(v1/v2,na.rm=TRUE,tr=tr)
if(OPT){
if(!HD)coef<-median(y,na.rm=TRUE)-slope*median(x,na.rm=TRUE)
if(HD)coef<-hd(y,na.rm=TRUE)-slope*hd(x,na.rm=TRUE,tr=tr)
}
if(!OPT){
if(!HD)coef<-median(y-slope*x,na.rm=TRUE)
if(HD)coef<-hd(y-slope*x,na.rm=TRUE,tr=tr)
}
names(coef)<-"Intercept"
coef<-c(coef,slope)
if(plotit){
plot(x,y,xlab="X",ylab="Y")
abline(coef)
}
res<-y-slope*x-coef[1]
list(coef=coef,residuals=res)
}

# ============================================================================
# mdepreg.sub
# ============================================================================
mdepreg.sub<-function(X,theta){
np<-ncol(X)
p<-np-1
x<-X[,1:p]
y<-X[,np]
temp<-t(t(x)*theta[2:np])
yhat<-apply(temp,1,sum)+theta[1]
res<-y-yhat
val<-0-mregdepth(x,res)
val
}

# ============================================================================
# poireg
# ============================================================================
poireg<-function(x,y,xout=FALSE,outfun=outpro,plotit=FALSE,xlab="X",ylab="Y",
varfun=var,YHAT=FALSE,STAND=TRUE,...){
#
# Perform Poisson regression.
# The predictors are assumed to be stored in the n by p matrix x.
# The y values are typically count data (integers).
#
# xout=T will remove outliers from among the x values and then fit
# the regression line.
#  Default:
# One predictor, a mad-median rule is used.
# With more than one, projection method is used.
#
# outfun=out will use MVE method
#
xy=elimna(cbind(x,y))
x<-as.matrix(x)
x=xy[,1:ncol(x)]
y=xy[,ncol(xy)]
x<-as.matrix(x)
if(xout){
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
temp=glm(formula=y~x,family=poisson)
init=summary(temp)
yhat=temp$coef[1]
for(j in 1:ncol(x)){
j1=j+1
yhat=yhat+temp$coef[j1]*x[,j]
}
yhat=exp(yhat)
if(plotit){
x=as.matrix(x)
if(ncol(x)>1)stop("Cannot plot with more than one predictor")
plot(x,y,xlab=xlab,ylab=ylab)
#points(x,yhat,pch=".")
xord=order(x)
lines(x[xord],yhat[xord])
init$coef
}
ex=varfun(yhat)/varfun(y)
str=sqrt(ex)
hatv=NULL
if(YHAT)hatv=yhat
list(results=init,Explanatory.Power=ex,Strength.Assoc=str,yhat=hatv)
}

# ============================================================================
# tsreg
# ============================================================================
tsreg<-function(x,y,xout=FALSE,outfun=outpro,iter=5,varfun=pbvar,tr=FALSE,do.stre=TRUE,
corfun=pbcor,plotit=FALSE,WARN=TRUE,HD=FALSE,OPT=FALSE,xlab='X',ylab='Y',...){
#
#  Compute Theil-Sen regression estimator
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#  By default, the number of iterations is
#  iter=5
#  starting with Rallfun-v35. It was 10 but this can result in high execution time and
#  does not appear to be necessary in terms of Type I errors.
# iter=1 can be unsatisfactory
#
#  OPT=TRUE, compute the intercept using median(y)-b_1median(X)
#  OPT=FALSE compute the intercept using median of y-b_1X
#
#  Starting with version Rallfun-v29, OPT=F is the default, which is consistent with
#  other R functions that have been supplied for computing the Theil--Sen estimator.
#
#do.stre=FALSE, strength not computed, can be useful when varfun can't be computed
#
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
n=nrow(x)
n.keep=n
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else flag<-outfun(x,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
n.keep=nrow(x)
}
if(ncol(x)==1){
temp1<-tsp1reg(x,y,HD=HD,OPT=OPT,tr=tr)
coef<-temp1$coef
res<-temp1$res
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-tsp1reg(x[,p],y,OPT=OPT,tr=tr)$coef[2]
}
res<-y-x%*%temp
if(!HD)alpha<-median(res)
if(HD)alpha<-hd(res,tr=tr)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-tsp1reg(x[,p],r[,p],plotit=FALSE,OPT=OPT,tr=tr)$coef[2]
}
if(!HD)alpha<-median(y-x%*%temp)
if(HD)alpha<-hd(y-x%*%temp,tr=tr)
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
yhat<-y-res
stre=e.pow=NULL
if(do.stre){
temp=varfun(y)
if(temp==0){
if(WARN)print("Warning: When computing strength of association, measure of variation=0")
}
e.pow=NULL
if(temp>0){
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}}}
if(plotit){
if(ncol(x)==1){
plot(x,y,xlab=xlab,ylab=ylab)
abline(coef)
}}
list(n=n,n.keep=n.keep,
coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

# ============================================================================
# regbootMC
# ============================================================================
regbootMC<-function(data,x,y,regfun,...){
vals=regfun(x[data,],y[data],...)$coef
}

# ============================================================================
# snmreg
# ============================================================================
snmreg<-function(x,y,SEED=TRUE,xout=FALSE,outfun=outpro,initreg=MMreg,...){
#
# Compute regression S-estimator via Nelder-Mead method
# The measure of scale is taken to be the percentage bend midvariance
#
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
x <- as.matrix(x)
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
npm1=np-1
x=X[,1:npm1]
x=as.matrix(x)
y=X[,np]
N<-np-1
#temp<-initreg(x,y,SEED=SEED)$coef
temp<-initreg(x,y)$coef
START<-temp[2:np]
temp<-nelder(X,N,FN=snmreg.sub,START=START)
alpha <- median(y - x %*% temp)
coef <- c(alpha,temp)
res <- y - x %*% temp - alpha
list(coef = coef, residuals = res)
}

# ============================================================================
# mopreg
# ============================================================================
mopreg<-function(x,y,regfun=tsreg,cop=3,KEEP=TRUE,MC=FALSE,STAND=TRUE){
#
# Do multiple (outcomes) regression on points not labled outliers
# using projection-type outlier detection method
# Arg=regfun determines regression method;
# by default, Theil-Sen is used.
#
#  KEEP=F, outliers will be eliminated
#  KEEP=T, outliers are not eliminated
# cop: see function outpro
library(MASS)
if(MC)library(parallel)
x<-as.matrix(x)
y<-as.matrix(y)
px<-ncol(x)
py<-ncol(y)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
if(KEEP)ivec<-c(1:nrow(x))
if(!KEEP){
if(!MC)ivec<-outpro(m,plotit=FALSE,cop=cop,STAND=STAND)$keep
if(MC)ivec<-outproMC(m,plotit=FALSE,cop=cop,STAND=STAND)$keep
}
np1<-ncol(x)+1
vec<-rep(1,nrow(m))
pxpy<-px+py
coef<-matrix(ncol=py,nrow=np1)
res<-matrix(ncol=py,nrow=nrow(m))
for(i in 1:py){
pv<-px+i
coef[,i]<-regfun(m[ivec,1:ncol(x)],m[ivec,pv])$coef
vec<-as.matrix(vec)
res[,i]<-m[,pv]-cbind(vec,m[,1:ncol(x)])%*%coef[,i]
}
list(coef=coef,residuals=res)
}

# ============================================================================
# mdepreg.coef
# ============================================================================
mdepreg.coef<-function(x,y,xout=FALSE,outfun=out,...){
#
# multiple depth regression
#
X<-cbind(x,y)
X<-elimna(X)
p1=ncol(X)
p=p1-1
if(xout){
flag<-outfun(X[,1:p],plotit=FALSE,...)$keep
X<-X[flag,]
}
library(mrfDepth)
a=rdepthmedian(X)$deep
list(coef=a)
}

# ============================================================================
# mdepreg
# ============================================================================
mdepreg<-function(x,y,xout=FALSE,outfun=out,RES=FALSE,...){
#
# multiple depth regression
#
X<-cbind(x,y)
X<-elimna(X)
n=nrow(X)
p1=ncol(X)
p=p1-1
if(xout){
flag<-outfun(X[,1:p],plotit=FALSE,...)$keep
X<-X[flag,]
}
n.keep=nrow(X)
library(mrfDepth)
a=rdepthmedian(X)$deepest
res=NA
if(RES)res=X[,p1]-X[,1:p]%*%a[2:p1]-a[1]
list(n=n,n.keep=n.keep,coef=a,residuals=res)
}

# ============================================================================
# mdepreg.orig
# ============================================================================
mdepreg.orig<-function(x,y,xout=FALSE,outfun=outpro){
#
# multiple depth regression
#
X<-cbind(x,y)
X<-elimna(X)
np=n.keep=ncol(X)
p=np-1
if(xout){
id=outfun(X[,1:p],plotit=FALSE)$keep
X=X[id,]
n.keep=nrow(X)
}
if(np==2){
temp=depreg(X[,1],X[,2])
coef=temp$coef
res=temp$residuals
}
if(np>2){
N<-np-1
x=X[,1:N]
y=X[,np]
START<-tsreg(x,y)$coef
coef<-nelderv2(X,np,FN=mdepreg.sub,START=START)
x <- as.matrix(x)
res <- y - x %*% coef[2:np] - coef[1]
}
list(n=n,n.keep=n.keep,coef = coef, residuals = res)
}

# ============================================================================
# MMreg
# ============================================================================
MMreg<-function(x,y,RES=TRUE,xout=FALSE,outfun=outpro,STAND=TRUE,varfun=pbvar,corfun=pbcor,WARN=FALSE,...){
#
#  Compute MM regression estimate derived by Yohai (1987)
#  simply by calling the R function lmrob
#  This function will remove leverage points when
#  xout=T
#  using the outlier detection method indicated by
#  outfun, which defaults to the projection method.
#
library('robustbase')
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else flag<-outpro(x,STAND=STAND,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(!WARN)options(warn=-1)
temp=lmrob(y~x)
if(!WARN)options(warn=0)
coef=temp$coefficients
p1=ncol(x)+1
res<-y-x%*%coef[2:p1]-coef[1]
yhat<-y-res
stre=NULL
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}
if(!RES)res=NULL
list(coef=coef,residuals=res,Strength.Assoc=stre)
}

# ============================================================================
# opregpbMC
# ============================================================================
opregpbMC<-function(x,y,nboot=1000,alpha=.05,om=TRUE,ADJ=TRUE,cop=3,SEED=TRUE,
nullvec=rep(0,ncol(x)+1),plotit=TRUE,opdis=2,gval=sqrt(qchisq(.95,ncol(x)+1))){
#
#  Same as opregpb, only this function takes advantage of a multi-core
#  processor assuming one is availabe and that the R package
#  multicore has been installed.
#
# generate bootstrap estimates
# use projection-type outlier detection method followed by
# TS regression.
#
# om=T and ncol(x)>1, means an omnibus test is performed,
# otherwise only individual tests of parameters are performed.
#
# opdis=2, means that Mahalanobis distance is used
# opdis=1, means projection-type distance is used
#
# gval is critical value for projection-type outlier detection
# method
#
# ADJ=T, Adjust p-values as described in Section 11.1.5 of the text.
#
if(SEED)set.seed(2)
library(parallel)
x<-as.matrix(x)
m<-cbind(x,y)
p1<-ncol(x)+1
m<-elimna(m) # eliminate any rows with missing data
x<-m[,1:ncol(x)]
x<-as.matrix(x)
y<-m[,p1]
if(nrow(x)!=length(y))stop("Sample size of x differs from sample size of y")
if(!is.matrix(x))stop("Data should be stored in a matrix")
print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,regboot,x,y,regfun=opregMC)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
# using Hochberg method
bvec<-t(bvec)
dvec<-alpha/(c(1:ncol(x)))
test<-NA
icl0<-round(alpha*nboot/2)
icl<-round(alpha*nboot/(2*ncol(x)))
icu0<-nboot-icl0
icu<-nboot-icl
output<-matrix(0,p1,6)
dimnames(output)<-list(NULL,c("Param.","p.value","crit.p.value",
"ci.lower","ci.upper","s.e."))
pval<-NA
for(i in 1:p1){
output[i,1]<-i-1
se.val<-var(bvec[,i])
temp<-sort(bvec[,i])
output[i,6]<-sqrt(se.val)
if(i==1){
output[i,4]<-temp[icl0+1]
output[i,5]<-temp[icu0]
}
if(i>1){
output[i,4]<-temp[icl+1]
output[i,5]<-temp[icu]
}
pval[i]<-sum((temp>nullvec[i]))/length(temp)
if(pval[i]>.5)pval[i]<-1-pval[i]
}
fac<-2
if(ADJ){
# Adjust p-value if n<60
nval<-length(y)
if(nval<20)nval<-20
if(nval>60)nval<-60
fac<-2-(60-nval)/40
}
pval[1]<-2*pval[1]
pval[2:p1]<-fac*pval[2:p1]
output[,2]<-pval
temp2<-order(0-pval[2:p1])
zvec<-dvec[1:ncol(x)]
sigvec<-(test[temp2]>=zvec)
output[temp2+1,3]<-zvec
output[1,3]<-NA
output[,2]<-pval
om.pval<-NA
temp<-opregMC(x,y)$coef
if(om && ncol(x)>1){
temp2<-rbind(bvec[,2:p1],nullvec[2:p1])
if(opdis==1)dis<-pdisMC(temp2,center=temp[2:p1])
if(opdis==2){
cmat<-var(bvec[,2:p1]-apply(bvec[,2:p1],2,mean)+temp[2:p1])
dis<-mahalanobis(temp2,temp[2:p1],cmat)
}
om.pval<-sum((dis[nboot+1]<=dis[1:nboot]))/nboot
}
# do adjusted p-value
nval<-length(y)
if(nval<20)nval<-20
if(nval>60)nval<-60
adj.pval<-om.pval/2+(om.pval-om.pval/2)*(nval-20)/40
if(ncol(x)==2 && plotit){
plot(bvec[,2],bvec[,3],xlab="Slope 1",ylab="Slope 2")
temp.dis<-order(dis[1:nboot])
ic<-round((1-alpha)*nboot)
xx<-bvec[temp.dis[1:ic],2:3]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
list(output=output,om.pval=om.pval,adj.om.pval=adj.pval)
}

# ============================================================================
# opregMC
# ============================================================================
opregMC<-function(x,y,regfun=tsreg,cop=3,fast=FALSE,pr=TRUE,prres=FALSE,STAND=TRUE,xout=FALSE){
#
# Do regression on points not labled outliers
# using projection-type outlier detection method
#
#  Note: argument xout is not relevant here, but is included to avoid conflicts when using regci.
#
library(parallel)
x<-as.matrix(x)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
ivec<-outproMC(m,plotit=FALSE,cop=cop,fast=FALSE,pr=FALSE,STAND=STAND)$keep
np1<-ncol(x)+1
coef<-regfun(m[ivec,1:ncol(x)],m[ivec,np1])$coef
vec<-rep(1,length(y))
residuals<-y-cbind(vec,x)%*%coef
if(fast && pr){
print("Intercept, followed by slopes:")
print(coef)
if(prres){
print("Residuals:")
print(residuals)
}}
list(coef=coef,residuals=residuals)
}

# ============================================================================
# bkreg
# ============================================================================
bkreg<-function(x,y,kerfun=akerd,pyhat=FALSE,plotit=TRUE,xlab="X",ylab="Y",
zlab="Z",xout=FALSE,outfun=outpro,pr=TRUE,theta=50,phi=25,duplicate="error",
expand=.5,scale=FALSE,ticktype="simple",...){
#
# Kernel estimator for binary regression.
# (See Signorini and Jones, JASA, 2004, 119-)
#
x=as.matrix(x)
p=ncol(x)
p1=p+1
xx<-elimna(cbind(x,y))
x<-xx[,1:p]
y<-xx[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
x=as.matrix(x)
flag<-(y==1)
mv=sum(flag)
nv=sum(!flag)
phat<-NA
fhat<-kerfun(x[flag,],pyhat=TRUE,plotit=FALSE,pts=x)
ghat<-kerfun(x[!flag,],pyhat=TRUE,plotit=FALSE,pts=x)
phat<-mv*fhat/(mv*fhat+nv*ghat)
if(p==1){
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab)
flag2<-order(x)
#lines(x[flag2],phat[flag2])
lines(x[flag2],phat)
}}
if(p==2){
if(plotit){
library(akima)
if(pr){
if(!scale)print("With dependence, suggest using scale=T")
}
fitr<-phat
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
}}
if(!pyhat)phat<-"Done"
phat
}

# ============================================================================
# regpreCV
# ============================================================================
regpreCV<-function(x,y,regfun=tsreg,varfun=pbvar,adz=TRUE,model=NULL,locfun=mean,
xout=FALSE,outfun=out,
plotit=TRUE,xlab="Model Number",ylab="Prediction Error",...){
#
# Estimate the prediction error using the regression method
#   regfun in conjunction with leave-one-out cross-validation
#
#   The argument model should have list mode, model[[1]] indicates
#   which predictors are used in the first model. For example, storing
#   1,4 in model[[1]] means predictors 1 and 4 are being considered.
#   If model is not specified, and number of predictors is at most 5,
#   then all models are considered.
#
#   If adz=T, added to the models to be considered is where
#   all regression slopes are zero. That is, use measure of location only
#   corresponding to
#   locfun.
#
x<-as.matrix(x)
d<-ncol(x)
p1<-d+1
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
y<-temp[,d+1]
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(is.null(model)){
if(d<=5)model<-modgen(d,adz=adz)
if(d>5)model[[1]]<-c(1:ncol(x))
}
mout<-matrix(NA,length(model),3,dimnames=list(NULL,c("est.error",
"var.used","rank")))
for (imod in 1:length(model)){
nmod=length(model[[imod]])-1
temp=c(nmod:0)
mout[imod,2]=sum(model[[imod]]*10^temp)
#
if(sum(model[[imod]]==0)!=1){
xx<-x[,model[[imod]]]
xx<-as.matrix(xx)
mout[imod,1]<-regpecv(xx,y,regfun=regfun,varfun=varfun,...)
}
#
if(sum(model[[imod]]==0)==1){
mout[imod,1]<-locCV(y,varfun=varfun,locfun=locfun)
}}
mout[,3]=rank(mout[,1])
if(plotit)plot(c(1:nrow(mout)),mout[,1],xlab=xlab,ylab=ylab)
mout
}

# ============================================================================
# regci
# ============================================================================
regci<-function(x,y,regfun=tsreg,nboot=599,alpha=.05,SEED=TRUE,pr=TRUE,null.val=NULL,method='hoch',
xout=FALSE,outfun=outpro,plotit=FALSE,xlab="Predictor 1",ylab="Predictor 2",WARNS=FALSE,LABELS=FALSE,...){
#
#   Compute a .95 confidence interval for each of the parameters of
#   a linear regression equation. The default regression method is
#   the Theil-Sen estimator.
#
#   When using the least squares estimator, and when n<250, use
#   lsfitci instead.
#
#   The predictor values are assumed to be in the n by p matrix x.
#   The default number of bootstrap samples is nboot=599
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimated of
#   the first predictor, etc.
#
#   plotit=TRUE: If there are two predictors, plot 1-alpha confidence region based
#  on the bootstrap samples.
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
if(pr)print("Default for outfun is now outpro, not out")
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
estit=regfun(x,y,...)$coef
if(is.null(null.val))null.val=rep(0,p1)
flagF=FALSE
flagF=identical(regfun,tsreg)
if(flagF){if(pr){
if(sum(duplicated(y)>0))print("Duplicate values detected; tshdreg might have more power than tsreg")
}}
nv=length(y)
x<-as.matrix(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
if(!WARNS)options(warn=-1)
bvec<-apply(data,1,regboot,x,y,regfun,xout=FALSE,...)
options(warn=0)
#Leverage points already removed.
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
regci<-matrix(0,p1,6)
vlabs="Intercept"
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
if(LABELS)vlabs[2:p1]=labels(x)[[2]]
dimnames(regci)<-list(vlabs,c("ci.low","ci.up","Estimate","S.E.","p-value",'p.adj'))
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
se<-NA
pvec<-NA
for(i in 1:p1){
bsort<-sort(bvec[i,])
#pvec[i]<-(sum(bvec[i,]<0)+.5*sum(bvec[i,]==0))/nboot
pvec[i]<-(sum(bvec[i,]<null.val[i])+.5*sum(bvec[i,]==null.val[i]))/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
regci[i,1]<-bsort[ilow]
regci[i,2]<-bsort[ihi]
se[i]<-sqrt(var(bvec[i,]))
}
if(p1==3){
if(plotit){
plot(bvec[2,],bvec[3,],xlab=xlab,ylab=ylab)
}}
regci[,3]=estit
pvec<-2*pvec
regci[,4]=se
regci[,5]=regci[,6]=pvec
regci[2:p1,6]=p.adjust(pvec[2:p1],method=method)
list(regci=regci,n=nrem,n.keep=nv)
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
# reg1way
# ============================================================================
reg1way<-function(x,y,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,outfun=outpro,STAND=TRUE,AD=FALSE,alpha=.05,pr=TRUE,...){
#
#  Test hypothesis that for two or more independent groups, all regression parameters
#  (the intercepts and slopes) are equal
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen MANOVA type test statistic.
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun, which defaults to the MVE method.
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#
if(pr){
if(!xout)print("Might want to consider xout=T to  remove leverage points")
}
if(SEED)set.seed(2)
if(!is.list(x))stop("Argument x should have list mode")
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop("Something is wrong. Number of covariates differs among the groups being compared")
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
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
}}
x=lapply(x,as.matrix)
K=p1
est=matrix(NA,nrow=J,ncol=p1)
grpnum=NULL
for(j in 1:J)grpnum[j]=paste("Group",j)
vlabs="Intercept"
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)=list(grpnum,vlabs)
ecov=list()
ecovinv=list()
W=rep(0,p1)
gmean=rep(0,p1)
for(j in 1:J){
est[j,]=regfun(x[[j]],y[[j]],xout=FALSE,...)$coef
nv.keep[j]=nrow(x[[j]])
vals=matrix(NA,nrow=nboot,ncol=p1)
data<-matrix(sample(length(y[[j]]),size=length(y[[j]])*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
bvec<-lapply(data,regbootMC,x[[j]],y[[j]],regfun,...)
# bvec is a p+1 by nboot matrix.
vals=t(matl(bvec))
ecov[[j]]=var(vals)
ecovinv[[j]]=solve(ecov[[j]])  #W_j
gmean=gmean+ecovinv[[j]]%*%est[j,]
W=W+ecovinv[[j]]
}
estall=solve(W)%*%gmean
F=0
for(k in 1:K){
for(m in 1:K){
for(j in 1:J){
F=F+ecovinv[[j]][k,m]*(est[j,k]-estall[k])*(est[j,m]-estall[m])
}}}
pvalad=NULL
# if xout=F or AD=T, compute corrected critical value, stemming from Johansen
df=K*(J-1)
if(!xout || AD){
iden=diag(p1)
Aw=0
for(j in 1:J){
temp=iden-solve(W)%*%ecovinv[[j]]
tempsq=temp%*%temp
Aw=Aw+(sum(diag(tempsq))+(sum(diag(temp)))^2)/(nv[j]-1)
}
Aw=Aw/2
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
crit=qchisq(alval[i],df)
critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
if(F<critad)break
}
pvalad=1-irem/1000
}
#
pval=1-pchisq(F,df)
crit=qchisq(1-alpha,df)
critad=NULL
if(!xout || AD)critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
est=data.frame(est)
hlabs=NULL
vlabs="Intercept"
for(j in 1:J)hlabs[j]=paste("Group",j)
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)<-list(hlabs,vlabs)
list(n=nv,n.keep=nv.keep,test.stat=F,crit.value=crit,adjusted.crit=critad,p.value=pval,adjusted.p.value=pvalad,est=est)
}

# ============================================================================
# reg1wayMC
# ============================================================================
reg1wayMC<-function(x,y,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,outfun=outpro,
STAND=TRUE,alpha=.05,pr=TRUE,AD=FALSE,...){
#
#  Test hypothesis that for two or more independent groups, all regression parameters are equal
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen MANOVA type test statistic
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#
library(parallel)
if(pr){
if(!xout)print("Might want to consider xout=T to  remove leverage points")
}
if(SEED)set.seed(2)
if(!is.list(x))stop("Argument x should have list mode")
if(pr){
if(xout)print("xout=T, so an adjusted critical is not computed and apparently not needed")
}
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop("Something is wrong. Number of covariates differs among the groups being compared")
nv=NULL
nv.keep=NULL
nv.all=NULL
p=ncol(x[[1]])
p1=p+1
for(j in 1:J){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1:p]
y[[j]]=xy[,p1]
x[[j]]=as.matrix(x[[j]])
nv.all[j]=c(nv,nrow(x[[j]]))
}
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
}}
x=lapply(x,as.matrix)
p1=ncol(x[[1]])+1
K=p1
est=matrix(NA,nrow=J,ncol=p1)
hlabs=NULL
vlabs="Intercept"
for(j in 1:J)hlabs[j]=paste("Group",j)
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)<-list(hlabs,vlabs)
nv=NA
ecov=list()
ecovinv=list()
W=rep(0,p1)
gmean=rep(0,p1)
for(j in 1:J){
est[j,]=regfun(x[[j]],y[[j]])$coef
nv.keep[j]=nrow(x[[j]])
nv[j]=nv.keep[j]
vals=matrix(NA,nrow=nboot,ncol=p1)
data<-matrix(sample(nv[j],size=nv[j]*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
bvec<-mclapply(data,regbootMC,x[[j]],y[[j]],regfun,...)
vals=t(matl(bvec))
ecov[[j]]=var(vals)
ecovinv[[j]]=solve(ecov[[j]])  #W_j
gmean=gmean+ecovinv[[j]]%*%est[j,]
W=W+ecovinv[[j]]
}
estall=solve(W)%*%gmean
F=0
for(k in 1:K){
for(m in 1:K){
for(j in 1:J){
F=F+ecovinv[[j]][k,m]*(est[j,k]-estall[k])*(est[j,m]-estall[m])
}}}
df=K*(J-1)
pvalad=NULL
# if xout=F, compute corrected critical value, stemming from Johansen
df=K*(J-1)
if(!xout || AD){
iden=diag(p1)
Aw=0
for(j in 1:J){
temp=iden-solve(W)%*%ecovinv[[j]]
tempsq=temp%*%temp
Aw=Aw+(sum(diag(tempsq))+(sum(diag(temp)))^2)/(nv[j]-1)
}
Aw=Aw/2
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
crit=qchisq(alval[i],df)
critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
if(F<critad)break
}
pavida=1-irem/1000
}
#
pval=1-pchisq(F,df)
crit=qchisq(1-alpha,df)
critad=NULL
if(AD)critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
est=data.frame(est)
list(n=nv.all,n.keep=nv.keep,test.stat=F,crit.value=crit,adjusted.crit=critad,p.value=pval,adjusted.p.value=pvalad,est=est)
}

# ============================================================================
# reg2ciMC
# ============================================================================
reg2ciMC<-function(x,y,x1,y1,regfun=tsreg,nboot=599,alpha=.05,plotit=TRUE,SEED=TRUE,
xout=FALSE,outfun=outpro,pr=FALSE,xlab='X',ylab='Y',...){
#
#   Compute a .95 confidence interval for the difference between the
#   the intercepts and slopes corresponding to two independent groups.
#   The default regression method is Theil-Sen.
#
#   Same as reg2ci, only takes advantage of a multi-core processor
#
#   The predictor values for the first group are
#   assumed to be in the n by p matrix x.
#   The predictors for the second group are in x1
#
#   The default number of bootstrap samples is nboot=599
#
#   regfun can be any R function that returns the coefficients in
#   the vector regfun$coef, the first element of which contains the
#   estimated intercept, the second element contains the estimate of
#   the first predictor, etc.
#
library(parallel)
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
x1<-as.matrix(x1)
xx1<-cbind(x1,y1)
xx1<-elimna(xx1)
x1<-xx1[,1:ncol(x1)]
x1<-as.matrix(x1)
y1<-xx1[,ncol(x1)+1]
x=as.matrix(x)
x1=as.matrix(x1)
if(xout){
if(identical(outfun,outblp)){
flag1=outblp(x,y,plotit=FALSE)$keep
flag2=outblp(x1,y2,plotit=FALSE)$keep
}
if(!identical(outfun,outblp)){
flag1=outfun(x,plotit=FALSE)$keep
flag2=outfun(x1,plotit=FALSE)$keep
}
x=x[flag1,]
y=y[flag1]
x1=x1[flag2,]
y1=y1[flag2]
}
n=length(y)
n[2]=length(y1)
x<-as.matrix(x)
x1<-as.matrix(x1)
est1=regfun(x,y)$coef
est2=regfun(x1,y1)$coef
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec<-mclapply(data,regbootMC,x,y,regfun,mc.preschedule=TRUE,xout=FALSE,...)
bvec=matl(bvec)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec1<-mclapply(data,regbootMC,x1,y1,regfun,mc.preschedule=TRUE,xout=FALSE,...)
bvec1=matl(bvec1)
bvec<-bvec-bvec1
p1<-ncol(x)+1
regci<-matrix(0,p1,6)
dimnames(regci)<-list(NULL,
c("Parameter","ci.lower","ci.upper","p.value","Group 1","Group 2"))
ilow<-round((alpha/2)*nboot)+1
ihi<-nboot-(ilow-1)
for(i in 1:p1){
temp<-sum(bvec[i,]<0)/nboot+sum(bvec[i,]==0)/(2*nboot)
regci[i,4]<-2*min(temp,1-temp)
bsort<-sort(bvec[i,])
regci[i,2]<-bsort[ilow]
regci[i,3]<-bsort[ihi]
regci[,1]<-c(0:ncol(x))
}
regci[,5]=est1
regci[,6]=est2
if(ncol(x)==1 && plotit){
plot(c(x,x1),c(y,y1),type="n",xlab=xlab,ylab=ylab)
points(x,y)
points(x1,y1,pch="+")
abline(regfun(x,y,...)$coef)
abline(regfun(x1,y1,...)$coef,lty=2)
}
list(n=n,output=regci)
}

# ============================================================================
# reg1wayISO
# ============================================================================
reg1wayISO<-function(x,y,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,outfun=outpro,STAND=TRUE,alpha=.05,pr=TRUE,...){
#
#  Test hypothesis that for two or more independent groups, all slope parameters
#  are equal.
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen MANOVA type test statistic.
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun
#
if(SEED)set.seed(2)
if(pr){
if(!xout)print("Might want to consider xout=T to  remove leverage points")
}
if(!is.list(x))stop("Argument x should have list mode")
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop("Something is wrong. Number of covariates differs among the groups being compared")
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
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
}}
x=lapply(x,as.matrix)
K=p1
est=matrix(NA,nrow=J,ncol=p1)
nv.keep=NULL
ecov=list()
ecovinv=list()
W=rep(0,p1)
gmean=rep(0,p)
for(j in 1:J){
est[j,]=regfun(x[[j]],y[[j]],xout=FALSE,...)$coef
nv.keep[j]=nrow(x[[j]])
vals=matrix(NA,nrow=nboot,ncol=p1)
data<-matrix(sample(length(y[[j]]),size=length(y[[j]])*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
bvec<-lapply(data,regbootMC,x[[j]],y[[j]],regfun,...)
# bvec is a p+1 by nboot matrix.
vals=t(matl(bvec))
ecov[[j]]=var(vals)
ecovinv[[j]]=solve(ecov[[j]])  #W_j
gmean=gmean+ecovinv[[j]][2:K,2:K]%*%est[j,2:K]
W=W+ecovinv[[j]]
}
estall=solve(W[2:K,2:K])%*%gmean
estall=c(0,estall)
F=0
for(k in 2:K){
for(m in 2:K){
for(j in 1:J){
F=F+ecovinv[[j]][k,m]*(est[j,k]-estall[k])*(est[j,m]-estall[m])
}}}
df=p*(J-1)
pvalad=NULL
AD=FALSE # Seems adjusted critical is not needed
if(AD){
iden=diag(p1)
Aw=0
for(j in 1:J){
temp=iden-solve(W)%*%ecovinv[[j]]
tempsq=temp%*%temp
Aw=Aw+(sum(diag(tempsq))+(sum(diag(temp)))^2)/(nv[j]-1)
}
Aw=Aw/2
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
crit=qchisq(alval[i],df)
critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
if(F<critad)break
}
pvalad=1-irem/1000
}
#
pval=1-pchisq(F,df)
crit=qchisq(1-alpha,df)
critad=NULL
if(AD)critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
est=data.frame(est)
hlabs=NULL
vlabs="Intercept"
for(j in 1:J)hlabs[j]=paste("Group",j)
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)<-list(hlabs,vlabs)
list(n=nv,n.keep=nv.keep,test.stat=F,crit.value=crit,p.value=pval,est=est)
}

# ============================================================================
# reg1wayISOMC
# ============================================================================
reg1wayISOMC<-function(x,y,regfun=tsreg,nboot=100,SEED=TRUE,xout=FALSE,outfun=outpro,
STAND=TRUE,alpha=.05,pr=TRUE,...){
#
#  Test hypothesis that for two or more independent groups, all regression parameters are equal
#  By default the Theil--Sen estimator is used
#
#  Strategy: Use bootstrap estimate of standard errors followed by
#  Johansen MANOVA type test statistic
#
#  x and y are assumed to have list mode having length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#
library(parallel)
if(SEED)set.seed(2)
if(pr){
if(!xout)print("Might want to consider xout=T to  remove leverage points")
}
if(!is.list(x))stop("Argument x should have list mode")
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop("Something is wrong. Number of covariates differs among the groups being compared")
nv=NULL
nv.keep=NULL
nv.all=NULL
p=ncol(x[[1]])
p1=p+1
for(j in 1:J){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1:p]
y[[j]]=xy[,p1]
x[[j]]=as.matrix(x[[j]])
nv.all[j]=c(nv,nrow(x[[j]]))
}
if(xout){
temp=lapply(x,outfun,plotit=FALSE,STAND=STAND,...)
for(j in 1:J){
x[[j]]=x[[j]][temp[[j]]$keep,]
y[[j]]=y[[j]][temp[[j]]$keep]
}}
x=lapply(x,as.matrix)
p1=ncol(x[[1]])+1
K=p1
est=matrix(NA,nrow=J,ncol=p1)
hlabs=NULL
vlabs="Intercept"
for(j in 1:J)hlabs[j]=paste("Group",j)
for(j in 2:p1)vlabs[j]=paste("Slope",j-1)
dimnames(est)<-list(hlabs,vlabs)
nv=NA
ecov=list()
ecovinv=list()
W=rep(0,p1)
gmean=rep(0,p)
for(j in 1:J){
est[j,]=regfun(x[[j]],y[[j]])$coef
nv.keep[j]=nrow(x[[j]])
nv[j]=nv.keep[j]
vals=matrix(NA,nrow=nboot,ncol=p1)
data<-matrix(sample(nv[j],size=nv[j]*nboot,replace=TRUE),ncol=nboot)
data=listm(data)
bvec<-mclapply(data,regbootMC,x[[j]],y[[j]],regfun,...)
vals=t(matl(bvec))
ecov[[j]]=var(vals)
ecovinv[[j]]=solve(ecov[[j]])  #W_j
gmean=gmean+ecovinv[[j]][2:K,2:K]%*%est[j,2:K]
W=W+ecovinv[[j]]
}
estall=solve(W[2:K,2:K])%*%gmean
estall=c(0,estall)
F=0
for(k in 2:K){
for(m in 2:K){
for(j in 1:J){
F=F+ecovinv[[j]][k,m]*(est[j,k]-estall[k])*(est[j,m]-estall[m])
}}}
df=p*(J-1)
#
df=p*(J-1)
#
pval=1-pchisq(F,df)
crit=qchisq(1-alpha,df)
critad=NULL
#if(!xout || AD)critad=crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
est=data.frame(est)
list(n=nv.all,n.keep=nv.keep,test.stat=F,crit.value=crit,p.value=pval,est=est)
}

# ============================================================================
# tsregNW
# ============================================================================
tsregNW<-function(x,y,xout=FALSE,outfun=out,iter=10,varfun=pbvar,
corfun=pbcor,plotit=FALSE,tol=.0001,...){
#
#  Compute Theil-Sen regression estimator
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#
#
x<-as.matrix(x)
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
if(ncol(x)==1){
temp1<-tsp1reg(x,y)
coef<-temp1$coef
res<-temp1$res
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-tsp1reg(x[,p],y)$coef[2]
}
res<-y-x%*%temp
alpha<-median(res)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-tsp1reg(x[,p],r[,p],plotit=FALSE)$coef[2]
}
if(max(abs(temp-tempold))<tol)break
alpha<-median(y-x%*%temp)
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
yhat<-y-res
stre=NULL
temp=varfun(y)
#if(temp==0)print('Warning: When computing strength of association, measure of variation=0')
e.pow=NULL
if(temp>0){
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}}
res=NULL
list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

# ============================================================================
# reg2cimcp
# ============================================================================
reg2cimcp<-function(x,y,regfun=tsreg,nboot=599,alpha=0.05,
SEED=TRUE,xout=FALSE,outfun=out,...){
#
# Like reg2ci only x1 etc have list mode containing
# data for J>1 groups. For all pairs of groups are compared via a
# call to reg2ci.
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
res=reg2ci(x[[j]],y[[j]],x[[k]],y[[k]],regfun=regfun,nboot=nboot,alpha=alpha,
plotit=FALSE,xout=xout,outfun=outfun,WARN=FALSE,...)
print(paste('Group', j,'Group', k))
print(res)
}}}
}

# ============================================================================
# reg1mcp
# ============================================================================
reg1mcp<-function(x,y,regfun=tsreg,SEED=TRUE,nboot=100,xout=FALSE,outfun=outpro,STAND=TRUE,alpha=.05,
pr=TRUE,MC=FALSE,...){
#
#  Perform all pairwise comparisons of intercepts among J independent groups
#  Do the same of the first slope, followed by the 2nd slope, etc.
#
#  Control FWE via Hochberg's methods for each set of
#  (J^2-J)/2 parameters. That is, control FWE for the intercepts
#  Do the same for the first slope, etc.
#
#  #  x and y are assumed to have list mode having
#  length J equal to the number of groups
#  For example, x[[1]] and y[[1]] contain the data for group 1.
#
#   xout=T will eliminate leverage points using the function outfun,
#   which defaults to the projection method.
#
#  OUTPUT:
#   n is sample size after missing values are removed
#   nv.keep is sample size after leverage points are removed.
#   output contains all pairwise comparisons
#   For each parameter, FWE is controlled using Hochberg's method
#   So by default, for the intercepts,
#   all pairwise comparisons are performed with FWE=.05
#   For the first slope, all pairwise comparisons are performed
#   with FWE=.05, etc.
#
if(SEED)set.seed(2)
if(!is.list(x))stop('Argument x should have list mode')
if(!is.list(y))stop('Argument y should have list mode')
J=length(x) # number of groups
x=lapply(x,as.matrix)
pchk=lapply(x,ncol)
temp=matl(pchk)
if(var(as.vector(temp))!=0)stop('Something is wrong. Number of covariates differs among the groups being compared')
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
outp=matrix(NA,ncol=7,nrow=nr)
x=lapply(x,as.matrix)
rlab=rep('Intercept',tot)
xx=list()
yy=list()
iall=0
ivp=c(1,tot)-tot
for(ip in 1:p){
#iv=ip-1
rlab=c(rlab,rep(paste('slope',ip),tot))
}
i=0
sk=1+tot*p
st=seq(1,sk,tot)
st=st-1
for(j in 1:J){
for(k in 1:J){
if(j < k){
i=i+1
st=st+1
xx[[1]]=x[[j]][,1:p]
xx[[2]]=x[[k]][,1:p]
yy[[1]]=y[[j]]
yy[[2]]=y[[k]]
if(!MC)temp=reg2ci(xx[[1]],yy[[1]],xx[[2]],yy[[2]],regfun=regfun)$output
if(MC)temp=reg2ci(xx[[1]],yy[[1]],xx[[2]],yy[[2]],regfun=regfun)$output
iall=iall+1
outp[iall,1]=j
outp[iall,2]=k
outp[st,3]=temp[,4]
outp[st,5]=temp[,2]
outp[st,6]=temp[,3]
}}}
for(i in 1:p1){
ivp=ivp+tot
temp2<-order(0-outp[ivp[1]:ivp[2],3])
icc=c(ivp[1]:ivp[2])
icc[temp2]=dvec
outp[ivp[1]:ivp[2],4]=icc
}
flag=(outp[,3]<=outp[,4])
outp[,7]=rep(0,nr)
outp[flag,7]=1
v=outp[1:tot,1]
vall=rep(v,p1)
outp[,1]=vall
v=outp[1:tot,2]
vall=rep(v,p1)
outp[,2]=vall
#outp[,7]=p.adjust(outp[,3],method=method)
dimnames(outp)=list(rlab,c('Group','Group','p.value','p.crit','ci.low','ci.hi','Sig'))
list(n=nv,n.keep=nv.keep,output=outp)
}

# ============================================================================
# chregF
# ============================================================================
chregF<-function(x,y,bend=1.345,SEED=FALSE,xout=FALSE,outfun=out,...){
#
# Compute Coakley Hettmansperger robust regression estimators
# JASA, 1993, 88, 872-880
#
# x is a n by p matrix containing the predictor values.
#
# No missing values are allowed
#
#  Comments in this function follow the notation used
#  by Coakley and Hettmansperger
#
library(MASS)
# with old version of R, need library(lqs) when using ltsreg
# as the initial estimate.
#
if(SEED)set.seed(12) # Set seed so that results are always duplicated.
x<-as.matrix(x)
p<-ncol(x)
m<-elimna(cbind(x,y))
x<-m[,1:p]
p1<-p+1
y<-m[,p1]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
x<-as.matrix(x)
cutoff<-bend
mve<-vector('list')
if(ncol(x)==1){
mve$center<-median(x)
mve$cov<-mad(x)^2
}
if(ncol(x)>=2)mve<-cov.mve(x)  # compute minimum volume ellipsoid measures of
                 # location and scale and store in mve.
reg0<-ltsReg(x,y) # compute initial regression est using least trimmed
                 # squares.
# Next, compute the rob-md2(i) values and store in rob
rob<-1  # Initialize vector rob
mx<-mve$center
rob<-mahalanobis(x,mx,mve$cov)
k21<-qchisq(.95,p)
c62<-k21/rob
vecone<-c(rep(1,length(y))) # Initialize vector vecone to 1
c30<-pmin(vecone,c62)  # mallows weights put in c30
k81<-median(abs(reg0$residuals)) # median of absolute residuals
k72<-1.4826*(1+(5/(length(y)-p-1)))*k81 # lms scale
c60<-reg0$residuals/(k72*c30) # standardized residuals
#  compute psi and store in c27
cvec<-c(rep(cutoff,length(y))) # Initialize vector cvec to cutoff
c27<-pmin(cvec,c60)
c27<-pmax(-1*cutoff,c27)  #c27 contains psi values
#
# compute B matrix and put in c66.
#  Also, transform B so that i th diag elem = 0 if c27[i] is
#  between -cutoff and cutoff, 1 otherwise.
#
c66<-ifelse(abs(c27)<=bend,1,0) # Have derivative of psi in c66
m1<-cbind(1,x)  # X matrix with col of 1's added
m2<-t(m1)   #X transpose
m5<-diag(c30) # matrix W, diagonal contains weights
m4<-diag(c66) # B matrix
m6<-m4%*%m1   # BX
m7<-m2%*%m6   # X'BX (nD=X'BX)
m8<-solve(m7)  #m8 = (X'-B-X)inverse
m9<-m8%*%m2 #m9=X prime-B-X inverse X'
m9<-m9%*%m5 # m9=X prime-B-X inverse X'W
m10<-m9%*%c27
c20<-m10*k72
c21<-reg0$coef+c20 #update initial estimate of parameters.
res<-y-m1%*%c21
list(coef=t(c21),residuals=res)
}

# ============================================================================
# DregGOLS
# ============================================================================
DregGOLS<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,SEED=TRUE,nboot=200,
STAND=TRUE,...){
#
#  Global test that two dependent (time 1 and time 2)
#  OLS regression lines are identical
#
#  Use a variation of Hotelling's test coupled with a bootstrap
#  estimate of the relevant covariance matrix associated with the differences
#  in the estimates of the parameters.
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
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
opf=identical(outfun,outpro)
if(!opf){
flag1=outfun(x1)$out.id
flag2=outfun(x2)$out.id
}
if(opf){
flag1=outpro(x1,STAND=STAND)$out.id
flag2=outfun(x2,STAND=STAND)$out.id
}
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
nk=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data,1,regboot,x1,y1,regfun=lsfit,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-apply(data,1,regboot,x2,y2,regfun=lsfit,...)
dif=t(bvec1-bvec2)
S=cov(dif)
est1=lsfit(x1,y1)$coef
est2=lsfit(x2,y2)$coef
est=est1-est2
k <- (nk-p1)/((nk - 1)*p1)
        stat <- k * crossprod(est, solve(S, est))[1, ]
        pvalue <- 1 - pf(stat, p1, nk - p1)
list(test.statistic = stat, degrees_of_freedom = c(p1, nk - p1), p.value =
pvalue,est.1=est1,est.2=est2,estimate.dif = est)
}

# ============================================================================
# difregOLS
# ============================================================================
difregOLS<-function(x1,y1,x2,y2,regfun=lsfit,xout=FALSE,outfun=outpro,nboot=200,
alpha=.05,SEED=TRUE,plotit=FALSE,xlab='X',ylab='Y',...){
#
# OLS regression data from two different times i.e., two dependent groups
#
#  compute confidence interval for the difference between intercepts
#  and the slopes
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
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
flag1=outfun(x1)$out.id
flag2=outfun(x2)$out.id
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
nk=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data,1,regboot,x1,y1,regfun,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-apply(data,1,regboot,x2,y2,regfun,...)
dif=t(bvec1)-t(bvec2)
est1=lsfit(x1,y1)$coef
est2=lsfit(x2,y2)$coef
estdif=est1-est2
se=apply(dif,2,sd)
pvec=NA
test=NA
test=estdif/se
df=nk-1
pvec=2*(1-pt(abs(test),df))
if(plotit){
reg2plot(x1,y1,x2,y2,xlab=xlab,ylab=ylab)
}
lvec='Intercept'
ci=matrix(NA,nrow=p1,ncol=3)
ci[,1]=c(0:p)
ci[,2]=estdif+qt(alpha/2,df)*se
ci[,3]=estdif-qt(alpha/2,df)*se
dimnames(ci)=list(NULL,c('Param','ci.low','ci.hi'))
for(j in 2:p1)lvec=c(lvec,paste('slope',j-1))
pvec=array(pvec,dimnames=lvec)
list(n=n,n.keep=nk,est.dif=estdif,est.1=est1,est.2=est2,
test.stat=test,standard.error=se,p.values=pvec,conf.intervals=ci)
}

# ============================================================================
# difregYvar
# ============================================================================
difregYvar<-function(x1,y1,x2,y2,regfun=tsreg,pts=NULL,
nboot=100,xout=FALSE,outfun=out,SEED=TRUE,...){
#
#  Estimate standard error of difference between the predicted value of Y
#  corresponding to two dependent groups using regression estimator indicated by
#  the argument
#  regfun
#  corresponding to the points in
#  pts
#  regfun defaults to tsreg, the Theil--Sen estimator
#  pts default is to use all unique points among x1 and x2
#
X=elimna(cbind(x1,y1,x2,y2))
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
x1<-as.matrix(x1)
x2=as.matrix(x2)
if(is.null(pts)){
pts=rbind(x1,x2)
pts=unique(pts)
}
pts=as.matrix(pts)
nvpts=nrow(pts)
bvec1=matrix(NA,nrow=nboot,ncol=nvpts)
bvec2=matrix(NA,nrow=nboot,ncol=nvpts)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
for(ib in 1:nboot){
bvec1[ib,]=regYsub(x1[data[ib,],],y1[data[ib,]],pts,p1=p1,regfun=regfun,...)
bvec2[ib,]=regYsub(x2[data[ib,],],y2[data[ib,]],pts,p1=p1,regfun=regfun,...)
}
bvec=bvec1-bvec2
sqsd=apply(bvec,2,var)
sqsd
}

# ============================================================================
# difreg
# ============================================================================
difreg<-function(x1,y1,x2,y2,regfun=tsreg,xout=FALSE,outfun=outpro,nboot=599,
alpha=.05,SEED=TRUE,plotit=FALSE,xlab='X',ylab='Y',pr=TRUE,...){
#
# regression data from two different times i.e., two dependent groups
#
#  compute confidence interval for the difference in the slopes
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
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
flag1=outfun(x1,...)$out.id
flag2=outfun(x2,...)$out.id
flag=unique(c(flag1,flag2))
X=X[-flag,]
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
nk=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec1<-lapply(data,regboot,x1,y1,regfun,xout=FALSE,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-lapply(data,regboot,x2,y2,regfun,xout=FALSE,...)
bvec1=matl(bvec1)
bvec2=matl(bvec2)
dif=t(bvec1)-t(bvec2)
dif.sort=apply(dif,2,sort)
pvec=NA
for(i in 1:p1){
pvec[i]<-(sum(dif[,i]<0)+.5*sum(dif[,i]==0))/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
}
pvec<-2*pvec
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=matrix(NA,nrow=p1,ncol=3)
ci[,1]=c(0:p)
for(i in 1:p1){
ci[i,2]=dif.sort[ilow,i]
ci[i,3]=dif.sort[ihi,i]
}
dimnames(ci)=list(NULL,c('Param','ci.low','ci.hi'))
if(plotit){
reg2plot(x1,y1,x2,y2,xlab=xlab,ylab=ylab,regfun=regfun,...)
}
lvec='Intercept'
for(j in 2:p1)lvec=c(lvec,paste('slope',j-1))
#pvec=array(pvec,dimnames=lvec)
est1=regfun(x1,y1,xout=FALSE,...)$coef
est2=regfun(x2,y2,xout=FALSE,...)$coef
list(n=n,n.keep=nk,param=lvec,p.values=pvec,est.grp1=est1,est.grp2=est2,conf.intervals=ci)
}

# ============================================================================
# tshdreg
# ============================================================================
tshdreg<-function(x,y,HD=TRUE,xout=FALSE,outfun=out,iter=5,varfun=pbvar,tr=FALSE,do.stre=TRUE,
corfun=pbcor,plotit=FALSE,tol=.0001,RES=TRUE,OPT=FALSE,xlab='X',ylab='Y',...){
#
#  Compute Theil-Sen regression estimator
#
#  Use back-fitting
#  when there is more than one predictor
#  and estimate intercept using Harrel-Davis estimator
#
x<-as.matrix(x)
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
if(ncol(x)==1){
temp1<-tshd(x,y,HD=HD,plotit=plotit,xlab=xlab,ylab=ylab,OPT=OPT,tr=tr)
coef<-temp1$coef
res<-y-coef[2]*x-coef[1]
}
if(ncol(x)>1){
for(p in 1:ncol(x)){
temp[p]<-tshd(x[,p],y)$coef[2]
}
res<-y-x%*%temp
alpha<-hd(res)
r<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
tempold<-temp
for(it in 1:iter){
for(p in 1:ncol(x)){
r[,p]<-y-x%*%temp-alpha+temp[p]*x[,p]
temp[p]<-tshd(x[,p],r[,p],plotit=FALSE,tr=tr)$coef[2]
}
if(max(abs(temp-tempold))<tol)break
alpha<-hd(y-x%*%temp,tr=tr)
tempold<-temp
}
coef<-c(alpha,temp)
res<-y-x%*%temp-alpha
}
yhat<-y-res
stre=NULL
e.pow=NULL
if(do.stre){
temp=varfun(y)
if(temp==0)print('Warning: When computing strength of association, measure of variation=0')
if(temp>0){
e.pow<-varfun(yhat)/varfun(y)
if(!is.na(e.pow)){
if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
e.pow=as.numeric(e.pow)
stre=sqrt(e.pow)
}}}
if(!RES)res=NULL
list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow,residuals=res)
}

# ============================================================================
# ltsreg
# ============================================================================
ltsreg<-function(x,y,tr=.5,xout=FALSE,outfun=outpro,STAND=TRUE,...){
#
# Leasts trimmed squares regression via the function ltsReg in the
# R package robustbase
#
x<-as.matrix(x)
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
library(robustbase)
temp=ltsReg(y~x,alpha=1-tr)
#coef=ltsReg(y~x)[8]$coefficients
coef=temp[8]$coefficients
res=temp[7]$raw.resid
list(coef=coef,residuals=res)
}

# ============================================================================
# ltsreg.2
# ============================================================================
ltsreg.2<-function(x,y,tr=.2,xout=FALSE,outfun=outpro,STAND=TRUE,...){
#
# Leasts trimmed squares regression via the function ltsReg in the
# R package robustbase
#
x<-as.matrix(x)
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
library(robustbase)
temp=ltsReg(y~x,alpha=1-tr)
coef=temp[8]$coefficients
res=temp[7]$raw.resid
list(coef=coef,residuals=res)
}

# ============================================================================
# DregG
# ============================================================================
DregG<-function(x1,y1,x2,y2,nullv=NULL,regfun=tshdreg,nboot=500,xout=FALSE,outfun=outpro,
SEED=TRUE,plotit=FALSE,pr=TRUE,...){
#
#  Global test that two dependent groups have identical
#  regression parameters.
#
#  Use a variation of Hotelling's test coupled with a bootstrap
#  estimate of the relevant covariance matrix associated with the differences
#  in the estimates of the parameters.#  For OLS, use DregGOLS
#
#  (plotit=F is used so that in simulations, if xout=T, the seed is not
#  set everytime outpro is called.)
#
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
if(is.null(nullv))nullv=rep(0,p1)
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
n=length(y1)
if(xout){
flag1=outfun(x1,...)$out.id
flag2=outfun(x2,...)$out.id
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
x1=as.matrix(x1)
x2=as.matrix(x2)
flagF=FALSE
flagF1=identical(regfun,tsreg)
flagF1[2]=identical(regfun,tshdreg)
#flagF1[3]=identical(regfun,tshdreg_C) obsolete,now it causes an error
if(sum(flagF1)>0)flagF=TRUE
if(!flagF){if(pr){
if(sum(duplicated(y1)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
pr=FALSE
}
if(pr){
if(sum(duplicated(y2)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
}}
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data,1,regboot,x1,y1,regfun=regfun,xout=FALSE,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-apply(data,1,regboot,x2,y2,regfun=regfun,xout=FALSE,...)
dif=t(bvec1-bvec2)
temp<-pdis(rbind(dif,nullv))
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
est1=regfun(x1,y1)$coef
est2=regfun(x2,y2)$coef
est=est1-est2
list(p.value=sig.level,est.1=est1,est.2=est2,estimate.dif = est)
}

# ============================================================================
# DregGMC
# ============================================================================
DregGMC<-function(x1,y1,x2,y2,nullv=NULL,regfun=tsreg,nboot=500,xout=FALSE,outfun=outpro,
SEED=TRUE,plotit=FALSE,pr=TRUE,...){
#
#  Global test that two dependent groups have identical
#  regression parameters.
#
#  Use a variation of Hotelling's test coupled with a bootstrap
#  estimate of the relevant covariance matrix associated with the differences
#  in the estimates of the parameters.#  For OLS, use DregGOLS
#
#  (plotit=F is used so that in simulations, if xout=T, the seed is not
#  set everytime outpro is called.)
#
flag=FALSE
library(parallel)
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
x1=as.matrix(x1)
x2=as.matrix(x2)
p=ncol(x1)
p1=p+1
p2=p+2
p3=p1+p
p4=p3+1
if(is.null(nullv))nullv=rep(0,p1)
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
n=length(y1)
if(xout){
flag1=outfun(x1,...)$out.id
flag2=outfun(x2,...)$out.id
flag=unique(c(flag1,flag2))
if(length(flag)>0)X=X[-flag,]
x1=X[,1:p]
y1=X[,p1]
x2=X[,p2:p3]
y2=X[,p4]
}
flagF=FALSE
flagF1=identical(regfun,tsreg)
flagF1[2]=identical(regfun,tshdreg)
#flagF1[3]=identical(regfun,tshdreg_C)
if(sum(flagF1)>0)flagF=TRUE
if(!flagF){
if(pr){
if(sum(duplicated(y1)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
pr=FALSE
}
if(pr){
if(sum(duplicated(y2)>0))print('Duplicate values detected; regfun=tshdreg might have more power than tsreg')
}}
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec1=mclapply(data,regbootMC,x1,y1,regfun,xout=FALSE,...)
bvec1=matl(bvec1)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2=mclapply(data,regbootMC,x2,y2,regfun,xout=FALSE,...)
bvec2=matl(bvec2)
dif=t(bvec1-bvec2)
temp<-pdisMC(rbind(dif,nullv))
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
est1=regfun(x1,y1)$coef
est2=regfun(x2,y2)$coef
est=est1-est2
list(p.value=sig.level,est.1=est1,est.2=est2,estimate.dif = est)
}

# ============================================================================
# difregMC
# ============================================================================
difregMC<-function(x1,y1,x2,y2,regfun=tsreg,xout=FALSE,outfun=outpro,nboot=599,
alpha=.05,SEED=TRUE,plotit=FALSE,xlab='X',ylab='Y',pr=TRUE,...){
#
# regression data from two different times i.e., two dependent groups
#
#  compute confidence interval for the difference in the slopes and intercepts
#
library(parallel)
if(SEED)set.seed(2)
X=elimna(cbind(x1,y1,x2,y2))
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
flag1=outfun(x1,...)$out.id
flag2=outfun(x2,...)$out.id
flag=unique(c(flag1,flag2))
X=X[-flag,]
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
nk=length(y1)
x1=as.matrix(x1)
x2=as.matrix(x2)
data<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
bvec1<-mclapply(data,regboot,x1,y1,regfun,mc.preschedule=TRUE,xout=FALSE,...)
# bvec is a p+1 by nboot matrix. The first row
#                     contains the bootstrap intercepts, the second row
#                     contains the bootstrap values for first predictor, etc.
bvec2<-mclapply(data,regboot,x2,y2,regfun,mc.preschedule=TRUE,xout=FALSE,...)
bvec1=matl(bvec1)
bvec2=matl(bvec2)
dif=t(bvec1)-t(bvec2)
dif.sort=apply(dif,2,sort)
pvec=NA
for(i in 1:p1){
pvec[i]<-(sum(dif[,i]<0)+.5*sum(dif[,i]==0))/nboot
if(pvec[i]>.5)pvec[i]<-1-pvec[i]
}
pvec<-2*pvec
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=matrix(NA,nrow=p1,ncol=3)
ci[,1]=c(0:p)
for(i in 1:p1){
ci[i,2]=dif.sort[ilow,i]
ci[i,3]=dif.sort[ihi,i]
}
dimnames(ci)=list(NULL,c('Param','ci.low','ci.hi'))
if(plotit){
reg2plot(x1,y1,x2,y2,xlab=xlab,ylab=ylab,regfun=regfun,...)
}
lvec='Intercept'
for(j in 2:p1)lvec=c(lvec,paste('slope',j-1))
#pvec=array(pvec,dimnames=lvec)
est1=regfun(x1,y1,xout=FALSE,...)$coef
est2=regfun(x2,y2,xout=FALSE,...)$coef
list(n=n,n.keep=nk,param=lvec,p.values=pvec,est.grp1=est1,est.grp2=est2,conf.intervals=ci)
}

# ============================================================================
# Qreg
# ============================================================================
Qreg<-function(x,y,q=.5,xout=FALSE,outfun=outpro,res.vals=TRUE,plotit=FALSE,xlab='X',ylab='Y',pch='*',...){
#
# Quantile regression. Like the function qreg, but avoids computational
# problems that can arise when there are tied values among the dependent
# variable
#
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
xx=as.matrix(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
if(xout){
x<-as.matrix(x)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=plotit,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
init=ols(x,y)$coef
v=optim(init,qfun,x=x,y=y,q=q,method='BFGS')$par
p1=ncol(x)+1
res=NULL
if(res.vals)res<-y-x%*%v[2:p1]-v[1]
if(ncol(x)==1){
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab,pch=pch)
abline(v)
}}
list(coef=v,residuals=res)
}

# ============================================================================
# regcits
# ============================================================================
regcits<-function(x,y,regfun=tshdreg,nboot=599,alpha=.05,SEED=TRUE,pr=TRUE,
xout=FALSE,outfun=outpro,plotit=FALSE,xlab='Predictor 1',ylab='Predictor 2',
MC=TRUE,...){
if(MC)v=regciMC(x,y,regfun=regfun,nboot=nboot,alpha=alpha,SEED=SEED,pr=pr,xout=xout,
outfun=outfun,plotit=plotit,xlab=xlab,ylab=ylab,...)
if(!MC)v=regci(x,y,regfun=regfun,nboot=nboot,alpha=alpha,SEED=SEED,pr=pr,xout=xout,
outfun=outfun,plotit=plotit,xlab=xlab,ylab=ylab,...)
v
}

# ============================================================================
# scorreg
# ============================================================================
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

# ============================================================================
# scorregciH
# ============================================================================
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
library(parallel)
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

# ============================================================================
# scorreg.sub
# ============================================================================
scorreg.sub<-function(data,xy,corfun=corfun,outfun=outfun,ALL=ALL){
p1=ncol(xy)
p=p1-1
est<-scorreg(xy[data,1:p],xy[data,p1],corfun=corfun,SEED=FALSE,ALL=ALL)$cor
est
}

# ============================================================================
# scorregci
# ============================================================================
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
library(parallel)
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

# ============================================================================
# scorreg.cr
# ============================================================================
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
library(parallel)
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

# ============================================================================
# scorreg.cr.sub
# ============================================================================
scorreg.cr.sub<-function(x,corfun,p=p,nboot=500,ALL=ALL,outfun=outfun){
p1=p+1
v=scorregci(x[,1:p],x[,p1],SEED=FALSE,corfun=corfun,nboot=nboot,crit.pv=1,pr=FALSE,hoch=TRUE,ALL=ALL,outfun=outfun)$Estimates[,2]
mp=min(as.vector(v),na.rm=T)
mp
}

# ============================================================================
# corregci
# ============================================================================
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

# ============================================================================
# corregci.sub
# ============================================================================
corregci.sub<-function(isub,x,y,corfun){
p=ncol(x)
xmat<-matrix(x[isub,],nrow(x),ncol(x))
e=NA
for(j in 1:p)e[j]=corfun(xmat[,j],y[isub])$cor
e
}

# ============================================================================
# LMSreg
# ============================================================================
LMSreg<-function(x,y,xout=FALSE,outfun=outpro,...){
#
#  Least median of squares
#
library(MASS)
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
temp<-NA
x<-as.matrix(x)
n=nrow(x)
n.keep=n
if(xout){
x<-as.matrix(x)
flag<-outfun(x,plotit=FALSE,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
n.keep=nrow(x)
}
a=lmsreg(x,y)
list(n=n,n.keep=n.keep,coef=a[3]$coefficients)
}

# ============================================================================
# Qreghat
# ============================================================================
Qreghat<-function(x,y,xr=x,q=.5,xout=FALSE,outfun=outpro,plotit.pts=FALSE){
#
#
#
xy=elimna(cbind(x,y))
xr=as.matrix(xr)
x=as.matrix(x)
p=ncol(x)
p1=p+1
vals=NA
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
xr=as.matrix(xr)
est=Qreg(x,y,q=q)$coef
if(ncol(xr)!=p)xr=t(xr)   # for a single point, need to transpose.
yhat=est[1]+xr%*%est[2:p1]
if(plotit.pts)points(xr,yhat)
yhat
}

# ============================================================================
# reg.reglev
# ============================================================================
reg.reglev<-function(x,y,plotit=TRUE,xlab='X',ylab='Y',GEN=TRUE,regfun=tsreg,outfun=outpro,pr=TRUE,...){

#
#  Remove any bad leverage points detected by
#  the fit using the estimator indicated by regun
#
# GEN=TRUE: use a generalization of the Rousseeuw van Zomeren method
# GEN=FALSE: usw the Rousseeuw van Zomeren method. Unknown when if ever this older approach
#     offers an advantage.
#
xy=elimna(cbind(x,y))
n=nrow(xy)
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
x<-as.matrix(x)
keep=c(1:n)
if(!GEN)a=reglev(x,y,plotit=FALSE,SEED=FALSE)$bad.lev.points
else a=reglev.gen(x,y,plotit=FALSE,regfun=regfun,outfun=outfun)$bad.lev
if(length(a)>0)keep=keep[-a]
nk=length(y[keep])
e=regfun(x[keep,],y[keep],...)
list(n=n,n.keep=nk,coef=e$coef)
}

