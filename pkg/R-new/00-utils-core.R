# WRS Package - Core Utilities
# Foundation functions used throughout the package
#
# This file contains the most commonly called utility functions,
# extracted from Rallfun-v45.R based on dependency analysis.
#
# Top functions by call count:
#   - elimna: Called by 928 functions (47% of package!)
#   - listm: Called by 327 functions
#   - matl: Called by 219 functions
#
# These functions MUST be loaded before other modules.


# elimna
elimna<-function(m){
#
# remove any rows of data having missing values
#
DONE=FALSE
if(is.list(m) && is.matrix(m)){
z=pool.a.list(m)
m=matrix(z,ncol=ncol(m))
DONE=TRUE
}
if(!DONE){
if(is.list(m) && is.matrix(m[[1]])){
for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
e=m
DONE=TRUE
}}
if(!DONE){
if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
e=m
DONE=TRUE
}}
if(!DONE){
m<-as.matrix(m)
ikeep<-c(1:nrow(m))
for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
e<-m[ikeep[ikeep>=1],]
}
e
}


# listm
listm<-function(x){
#
# Store the data in a matrix or data frame in a new
# R variable having list mode.
# Col 1 will be stored in y[[1]], col 2 in y[[2]], and so on.
#
if(is.null(dim(x)))stop("The argument x must be a matrix or data frame")
y<-list()
for(j in 1:ncol(x))y[[j]]<-x[,j]
y
}

m2l=listm

matrix2list=listm


# matl
matl<-function(x){
#
# take data in list mode and store it in a matrix
#
J=length(x)
nval=NA
for(j in 1:J)nval[j]=length(x[[j]])
temp<-matrix(NA,ncol=J,nrow=max(nval))
for(j in 1:J)temp[1:nval[j],j]<-x[[j]]
temp
}

list2mat=matl


# near
near<-function(x,pt,fr=1){
# determine which values in x are near pt
# based on fr * mad
if(!is.vector(x))stop('x should be a vector')
m<-mad(x)
if(m==0){
temp<-idealf(x)
m<-(temp$qu-temp$ql)/(qnorm(.75)-qnorm(.25))
}
if(m==0)m<-sqrt(winvar(x)/.4129)
if(m==0)stop("All measures of dispersion are equal to 0")
dis<-abs(x-pt)
dflag<-dis <= fr*m
dflag
}


# hd
hd<-function(x,q=.5,na.rm=TRUE,STAND=NULL,tr=FALSE){
#
#  Compute the Harrell-Davis estimate of the qth quantile
#
#  The vector x contains the data,
#  and the desired quantile is q
#  The default value for q is .5.
#
if(tr)e=thd(x,q=q)
else{
if(na.rm)x=elimna(x)
n<-length(x)
m1<-(n+1)*q
m2<-(n+1)*(1-q)
vec<-seq(along=x)
w<-pbeta(vec/n,m1,m2)-pbeta((vec-1)/n,m1,m2)  # W sub i values
y<-sort(x)
e<-sum(w*y)
}
e
}


# winvar
winvar<-function(x,tr=.2,na.rm=FALSE,STAND=NULL){
#
#  Compute the gamma Winsorized variance for the data in the vector x.
#  tr is the amount of Winsorization which defaults to .2.
#
remx=x
x<-x[!is.na(x)]
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
y<-ifelse(y<=xbot,xbot,y)
y<-ifelse(y>=xtop,xtop,y)
wv<-var(y)
if(!na.rm)if(sum(is.na(remx)>0))wv=NA
wv
}


# lplot
lplot<-function(x,y,low.span=2/3,span=.75,pyhat=FALSE,eout=FALSE,xout=FALSE,outfun=outpro,plotit=TRUE,
expand=.5,varfun=pbvar,cor.op=FALSE,cor.fun=pbcor,pr=TRUE,STR=TRUE,ZLIM=FALSE,pch=NULL,
scale=TRUE,xlab="X",ylab="Y",zlab="",theta=50,phi=25,family="gaussian",
duplicate="error",pc="*",ticktype="simple",frame=TRUE,...){
#
# Plot regression surface using LOESS
#
# NOTE: For a single independent variable, the function lplotCI will plot a confidence band
#  that allows heteroscedasticity and has simultaneous probability coverage 1-alpha
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
if(pr){
if(!xout)print('Suggest also looking at result using xout=TRUE')
}
if(!is.null(pch))pc=pch
library(stats)
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n.orig=nrow(m)
n.keep=n.orig
if(!is.matrix(x))stop("x is not a matrix")
d<-ncol(x)
if(d>=2){
library(akima)
if(ncol(x)==2 & !scale){
if(pr){
print("scale=FALSE is specified.")
print("If there is dependence, might use scale=TRUE")
print("To get a p-value, based on the measure of the")
print("strength of association based on this function,")
print("use the function lplotPV")
}}
x<-m[,1:d]
y<-m[,d+1]
if(eout & xout)stop("Can't have both eout and xout = FALSE")
if(eout){
flag<-outfun(m,plotit=FALSE,...)$keep
m<-m[flag,]
n.keep=nrow(m)
}
if(xout){
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
n.keep=nrow(m)
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
if(!ZLIM)persp(fitr,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,expand=expand,
scale=scale,ticktype=ticktype)
if(ZLIM)persp(fitr,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,expand=expand,
scale=scale,ticktype=ticktype,zlim=c(0,1))  #used by logreg.plot
}}
if(d==1){
m<-elimna(cbind(x,y))
x<-m[,1:d]
y<-m[,d+1]
if(eout && xout)stop("Cannot have both eout and xout = T")
if(eout){
flag<-outfun(m,plotit=FALSE,...)$keep
m<-m[flag,]
n.keep=nrow(m)
}
if(xout){
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
n.keep=nrow(m)
}
x<-m[,1:d]
y<-m[,d+1]
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab,pch=pc,frame=frame)
lines(lowess(x,y,f=low.span))
}
tempxy<-lowess(x,y,f=low.span)
yyy<-tempxy$y
xxx<-tempxy$x
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
if(!STR)E.power=NA
if(STR){
E.power<-1
if(!cor.op)E.power<-varfun(last[!is.na(last)])/varfun(y)
if(cor.op || E.power>=1){
if(d==1){
xord<-order(x)
E.power<-cor.fun(last,y[xord])$cor^2
}
if(d>1)E.power<-cor.fun(last,y)$cor^2
}
E.power=as.numeric(E.power)
}
if(!pyhat)last <- NULL
list(Strength.Assoc=sqrt(E.power),Explanatory.power=E.power,yhat.values=last,n=n.orig,
n.keep=n.keep)
}

# qest
qest<-function(x,q=.5,na.rm=TRUE){
#
# Compute an estimate of qth quantile
#  using a single order statistic
#
if(na.rm)x<-elimna(x)
if(q<=0 || q>=1)stop("q must be > 0 and < 1")
n<-length(x)
xsort<-sort(x)
iq <- floor(q * n + 0.5)
qest<-NA
if(iq>0 || iq<=n)qest<-xsort[iq]
qest
}

# pool.a.list
pool.a.list<-function(x){
#
# x has list mode. Pool all of the data into a single R variable.
#
if(!is.list(x))stop("x should have list mode")
pts=NULL
for(j in 1:length(x))pts=c(pts,x[[j]])
pts
}


# con2way
con2way<-function(J,K){
#
# For a  J by K ANOVA design
# create the contrast coefficients for
# doing all pairwise comparisons of
# main effects for Factor A and B and all interactions
#
JK <- J * K
Ja<-(J^2-J)/2
Ka<-(K^2-K)/2
JK<-J*K
conA<-matrix(0,nrow=JK,ncol=Ja)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j < jj){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[j,]<-1
mat[jj,]<-0-1
conA[,ic]<-t(mat)
}}}
conB<-matrix(0,nrow=JK,ncol=Ka)
ic<-0
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[,k]<-1
mat[,kk]<-0-1
conB[,ic]<-t(mat)
}}}
conAB<-matrix(0,nrow=JK,ncol=Ka*Ja)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j < jj){
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[j,k]<-1
mat[j,kk]<-0-1
mat[jj,k]<-0-1
mat[jj,kk]<-1
}
conAB[,ic]<-t(mat)
}}}}}
list(conA=conA,conB=conB,conAB=conAB)
}


# rmul
rmul<-function(n,p=2,cmat=NULL,rho=0,
mar.fun=ghdist,OP=FALSE,g=0,h=0,...){
#
# generate n observations from a p-variate dist
# By default, use normal distributions.
#
# Can generate data form a g-and-h distribution via the arguments
#  g and h
#
# To adjust rho so that Pearson = remains equal to rho after transforming, use rngh
#
# Example rmul(30,p=4,rho=.3,g=.5,h=.2) will
# generate 30 vectors from a 4-variate distribution where the marginals
# have a g-and-h distribution with g=.5 and h=.2.
#
# This function is similar to ghmul, only here, generate the marginal values
# and then transform the data to have correlation matrix cmat
#
# cmat: if specified, is the correlation matrix that is used to generate data
#
# If not specified, data are generated with a common correlation
# rho
#
#OP= TRUE:
# Method (e.g. Browne, M. W. (1968) A comparison of factor analytic
# techniques. Psychometrika, 33, 267-334.
#  Let U'U=R be the Cholesky decomposition of R. Generate independent data
#  from some dist yielding X. Then XU has population correlation matrix R
#
#  OP=FALSE, use mvrnorm to generate data then transform marginals to g-and-h distribution.
#
if(!is.null(cmat)){
if(ncol(cmat)!=p)stop('cmat: number of  columns must equal the value in the argument p')
}
if(abs(rho)>1)stop('rho must be between -1 and 1')
if(is.null(cmat)){
cmat<-matrix(rho,p,p)
diag(cmat)<-1
}
if(OP){
np<-n*p
if(identical(mar.fun,ghdist))x<-matrix(mar.fun(np,g=g,h=h),nrow=n,ncol=p)
else x<-matrix(mar.fun(np,...),nrow=n,ncol=p)
rmat<-matsqrt(cmat)
x<-x%*%rmat
}
if(!OP){
library(MASS)
x=mvrnorm(n,rep(0,p),cmat)
if(g==0)x=x*exp(h*x^2/2)
if(g>0)x=(exp(g*x)-1)*exp(h*x^2/2)/g
}
x
}



# regYhat
regYhat<-function(x,y,xr=x,regfun=tsreg,xout=FALSE,outfun=outpro,pr=FALSE,plot.pts=FALSE,pts=NULL,...){
#
#  For convenience, return estimate of Y based on data in xr (or pts) using
#  regression line based on regfun
#
xy=elimna(cbind(x,y))
x<-as.matrix(x)
xr=as.matrix(xr)
p=ncol(x)
p1=p+1
vals=NA
x<-xy[,1:p]
y<-xy[,p1]
#print(xr[1:10,])
if(xout){
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
if(!is.null(pts[1]))xr=pts
xr=as.matrix(xr)
est=regfun(x,y,...)$coef
if(ncol(xr)!=p)xr=t(xr)   # for a single point, need to transpose.
yhat=est[1]+xr%*%est[2:p1]
if(plot.pts)points(xr,yhat)
yhat
}

reg.pred<-regYhat


# pdis
pdis<-function(m,pts=m,MM=FALSE,cop=3,dop=1,center=NA,na.rm=TRUE){
#
# Compute projection distances for points in pts relative to points in m
#  That is, the projection distance from the center of m
#
#
#  MM=F  Projected distance scaled
#  using interquatile range.
#  MM=T  Scale projected distances using MAD.
#
#  There are five options for computing the center of the
#  cloud of points when computing projections:
#  cop=1 uses Donoho-Gasko median
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#  cop=5 uses skipped mean
#
m<-elimna(m) # Remove missing values
pts=elimna(pts)
m<-as.matrix(m)
nm=nrow(m)
pts<-as.matrix(pts)
if(ncol(m)>1){
if(ncol(pts)==1)pts=t(pts)
}
npts=nrow(pts)
mp=rbind(m,pts)
np1=nrow(m)+1
if(ncol(m)==1){
m=as.vector(m)
pts=as.vector(pts)
if(is.na(center[1]))center<-median(m)
dis<-abs(pts-center)
disall=abs(m-center)
temp=idealf(disall)
if(!MM){
pdis<-dis/(temp$qu-temp$ql)
}
if(MM)pdis<-dis/mad(disall)
}
else{
if(is.na(center[1])){
if(cop==1)center<-dmean(m,tr=.5,dop=dop)
if(cop==2)center<-cov.mcd(m)$center
if(cop==3)center<-apply(m,2,median)
if(cop==4)center<-cov.mve(m)$center
if(cop==5)center<-smean(m)
}
dmat<-matrix(NA,ncol=nrow(mp),nrow=nrow(mp))
for (i in 1:nrow(mp)){
B<-mp[i,]-center
dis<-NA
BB<-B^2
bot<-sum(BB)
if(bot!=0){
for (j in 1:nrow(mp)){
A<-mp[j,]-center
temp<-sum(A*B)*B/bot
dis[j]<-sqrt(sum(temp^2))
}
dis.m=dis[1:nm]
if(!MM){
#temp<-idealf(dis)
temp<-idealf(dis.m)
dmat[,i]<-dis/(temp$qu-temp$ql)
}
if(MM)dmat[,i]<-dis/mad(dis.m)
}}
pdis<-apply(dmat,1,max,na.rm=na.rm)
pdis=pdis[np1:nrow(mp)]
}
pdis
}

# yuen
yuen<-function(x,y=NULL,tr=.2,alpha=.05){
#
#  Perform Yuen's test for trimmed means on the data in x and y.
#  The default amount of trimming is 20%
#  Missing values (values stored as NA) are automatically removed.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuen$ci.
#  The p-value is returned in yuen$p.value
#
#  x, y: The data for the two groups are stored in x and y
#  tr=.2: indicates that the default amount of trimming is .2
#         tr=0 results in using the sample mean
#
#  The function returns both a confidence interval and a p-value.
#
#  For an omnibus test with more than two independent groups,
#  use t1way.
#  This function uses winvar from chapter 2.
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
if(tr==.5)stop("Using tr=.5 is not allowed; use a method designed for medians")
if(tr>.25)print("Warning: with tr>.25 type I error control might be poor")
x<-x[!is.na(x)]  # Remove any missing values in x
y<-y[!is.na(y)]  # Remove any missing values in y
h1<-length(x)-2*floor(tr*length(x))
h2<-length(y)-2*floor(tr*length(y))
q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
crit<-qt(1-alpha/2,df)
dif<-mean(x,tr)-mean(y,tr)
low<-dif-crit*sqrt(q1+q2)
up<-dif+crit*sqrt(q1+q2)
test<-abs(dif/sqrt(q1+q2))
yuen<-2*(1-pt(test,df))
list(n1=length(x),n2=length(y),est.1=mean(x,tr),est.2=mean(y,tr),ci=c(low,up),p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,crit=crit,df=df)
}


# akerd
akerd<-function(xx,hval=NA,aval=.5,op=1,fr=.8,pyhat=FALSE,pts=NA,plotit=TRUE,
xlab="",ylab="",zlab="",theta=50,phi=25,expand=.5,scale=TRUE,ticktype="simple",color='black'){
#
# Compute adaptive kernel density estimate
#
# (See Silverman, 1986)
#
# op=1 Use expected frequency as initial estimate of the density
# op=2 Univariate case only
#      Use normal kernel to get initial estimate of the density
#  ticktype="detailed" will create ticks as done for a two-dimensional plot
#
#  Note, when pyhat=T, returns estimate of density at points in pts AFTER
#  putting the points in ascending order.
#
xx=elimna(xx)
fval<-"Done"
if(is.matrix(xx)){
if(ncol(xx)>1)fval<-akerdmul(xx,pts=pts,hval=hval,aval=aval,fr=fr,pr=pyhat,
plotit=plotit,theta=theta,phi=phi,expand=expand,scale=scale,ticktype=ticktype)
plotit<-F
}
if(is.matrix(xx) && ncol(xx)==1)xx<-xx[,1]
if(!is.matrix(xx)){
x<-sort(xx)
if(op==1){
m<-mad(x)
if(m==0){
temp<-idealf(x)
m<-(temp$qu-temp$ql)/(qnorm(.75)-qnorm(.25))
}
if(m==0)m<-sqrt(winvar(x)/.4129)
if(m==0)stop("All measures of dispersion are equal to 0")
fhat <- rdplot(x,pyhat=TRUE,plotit=FALSE,fr=fr)
if(m>0)fhat<-fhat/(2*fr*m)
}
if(op==2){
init<-density(xx)
fhat <- init$y
x<-init$x
}
n<-length(x)
if(is.na(hval)){
sig<-sqrt(var(x))
temp<-idealf(x)
iqr<-(temp$qu-temp$ql)/1.34
A<-min(c(sig,iqr))
if(A==0)A<-sqrt(winvar(x))/.64
hval<-1.06*A/length(x)^(.2)
# See Silverman, 1986, pp. 47-48
}
gm<-exp(mean(log(fhat[fhat>0])))
alam<-(fhat/gm)^(0-aval)
dhat<-NA
if(is.na(pts[1]))pts<-x
pts<-sort(pts)
for(j in 1:length(pts)){
temp<-(pts[j]-x)/(hval*alam)
epan<-ifelse(abs(temp)<sqrt(5),.75*(1-.2*temp^2)/sqrt(5),0)
dhat[j]<-mean(epan/(alam*hval))
}
if(plotit){
plot(pts,dhat,type="n",ylab=ylab,xlab=xlab)
lines(pts,dhat,col=color)
}
if(pyhat)fval<-dhat
}
fval
}



# chi.int2
chi.int2 <- function(p,a,c1)
#   partial expectation d in (c1,\infty) of d^a under chi-squared p
 return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*(1-pchisq(c1^2,p+a)))

# outpro
outpro<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
xlab="VAR 1",ylab="VAR 2",STAND=TRUE,tr=.2,q=.5,pr=TRUE,...){
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
# plotit=TRUE creates a scatterplot when working with
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
#  When using cop=2, 3 or 4, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
#  STAND=T means that marginal distributions are standardized before
#  checking for outliers.
#
#  Donoho-Gasko (Tukey) median is marked with a cross, +.
#
m<-as.matrix(m)
if(pr){
if(!STAND){
if(ncol(m)>1)print("STAND=FALSE. If measures are on different scales, might want to use STAND=TRUE")
}}
library(MASS)
m=elimna(m)
m<-as.matrix(m)
nv=nrow(m)
if(ncol(m)==1){
dis<-(m-median(m,na.rm=TRUE))^2/mad(m,na.rm=TRUE)^2
dis<-sqrt(dis)
dis[is.na(dis)]=0
crit<-sqrt(qchisq(.975,1))
chk<-ifelse(dis>crit,1,0)
vec<-c(1:nrow(m))
outid<-vec[chk==1]
keep<-vec[chk==0]
}
if(ncol(m)>1){
M=m
if(STAND)m=standm(m,est=median,scat=mad)
if(is.na(gval) && cop==1)gval<-sqrt(qchisq(.95,ncol(m)))
if(is.na(gval) && cop!=1)gval<-sqrt(qchisq(.975,ncol(m)))
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
for (i in 1:nrow(m)){
B<-m[i,]-center
dis<-NA
BB<-B^2
bot<-sum(BB)
if(bot!=0){
for (j in 1:nrow(m)){
A<-m[j,]-center
temp<-sum(A*B)*B/bot
dis[j]<-sqrt(sum(temp^2))
}
temp<-idealf(dis)
if(!MM)cu<-median(dis)+gval*(temp$qu-temp$ql)
if(MM)cu<-median(dis)+gval*mad(dis)
outid<-NA
temp2<-(dis> cu)
flag[temp2]<-1
}}
if(sum(flag) == 0) outid <- NA
if(sum(flag) > 0)flag<-(flag==1)
outid <- vec[flag]
idv<-c(1:nrow(m))
keep<-idv[!flag]
if(ncol(m)==2){
if(plotit){
m=M # plot data using the original scale.
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
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}}}
list(n=nv,n.out=length(outid),out.id=outid,keep=keep)
}


# idealf
idealf<-function(x,na.rm=FALSE){
#
# Compute the ideal fourths for data in x
#
if(na.rm)x<-x[!is.na(x)]
j<-floor(length(x)/4 + 5/12)
y<-sort(x)
g<-(length(x)/4)-j+(5/12)
ql<-(1-g)*y[j]+g*y[j+1]
k<-length(x)-j+1
qu<-(1-g)*y[k]+g*y[k-1]
list(ql=ql,qu=qu)
}


# runmean2g
runmean2g<-function(x1,y1,x2,y2,fr=.8,est=tmean,xlab="X",ylab="Y",SCAT=TRUE,
sm=FALSE,nboot=40,SEED=TRUE,eout=FALSE,xout=FALSE,outfun=out,LP=FALSE,pch1='*',pch2='+',...){
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
m<-elimna(cbind(x1,y1))
if(eout && xout)stop("Not allowed to have eout=xout=T")
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m<-m[flag,]
}
x1<-m[,1]
y1<-m[,2]
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
}
m<-elimna(cbind(x2,y2))
if(eout && xout)stop("Not allowed to have eout=xout=T")
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m<-m[flag,]
}
x2<-m[,1]
y2<-m[,2]
if(xout){
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
ord1=order(x1)
x1=x1[ord1]
y1=y1[ord1]
ord2=order(x2)
x2=x2[ord2]
y2=y2[ord2]
if(!sm){
temp<-rungen(x1,y1,est=est,fr=fr,pyhat=TRUE,plotit=FALSE,xout=FALSE,eout=FALSE,LP=LP,...)
rmd1<-temp[1]$output
}
if(sm){
temp<-runmbo(x1,y1,est=est,fr=fr,pyhat=TRUE,plotit=FALSE,SEED=SEED,
nboot=nboot,eout=FALSE,xout=FALSE,...)
rmd1<-temp
}
if(!sm){
temp<-rungen(x2,y2,fr=fr,est=est,pyhat=TRUE,plotit=FALSE,xout=FALSE,eout=FALSE,LP=LP,...)
rmd2<-temp[1]$output
}
if(sm){
temp<-runmbo(x2,y2,est=est,fr=fr,pyhat=TRUE,plotit=FALSE,SEED=SEED,
nboot=nboot,eout=FALSE,xout=FALSE,...)
rmd2<-temp
}
plot(c(x1,x2),c(y1,y2),type="n",xlab=xlab,ylab=ylab)
sx1<-sort(x1)
sx2<-sort(x2)
xorder1<-order(x1)
xorder2<-order(x2)
sysm1<-rmd1[xorder1]
sysm2<-rmd2[xorder2]
if(LP){
sysm1=lplot(sx1,sysm1[xorder1],plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
sysm2=lplot(sx2,sysm2,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
if(SCAT)points(x1,y1,pch=pch1)
if(!SCAT)points(x1,y1,type="n")
if(SCAT)points(x2,y2,pch=pch2)
if(!SCAT)points(x2,y2,type="n")
lines(sx1,sysm1)
lines(sx2,sysm2,lty=2)
}

# smmcrit
smmcrit<-function(nuhat,C){
#
#  Determine the .95 quantile of the C-variate Studentized maximum
#  modulus distribution using linear interpolation on inverse
#  degrees of freedom
#  If C=1, this function returns the .975 quantile of Student's t
#  distribution.
#
if(C-round(C)!=0)stop("The number of contrasts, C, must be an  integer")
if(C>=29)stop("C must be less than or equal to 28")
if(C<=0)stop("C must be greater than or equal to 1")
if(nuhat<2)stop("The degrees of freedom must be greater than or equal to 2")
if(C==1)smmcrit<-qt(.975,nuhat)
if(C>=2){
C<-C-1
m1<-matrix(0,20,27)
m1[1,]<-c(5.57,6.34,6.89,7.31,7.65,7.93,8.17,8.83,8.57,
8.74,8.89,9.03,9.16,9.28,9.39,9.49,9.59, 9.68,
9.77,9.85,9.92,10.00,10.07,10.13,10.20,10.26,10.32)
m1[2,]<-c(3.96,4.43,4.76,5.02,5.23,5.41,5.56,5.69,5.81,
5.92,6.01,6.10,6.18,6.26,6.33,6.39,6.45,6.51,
6.57,6.62,6.67,6.71,6.76,6.80,6.84,6.88, 6.92)
m1[3,]<-c(3.38,3.74,4.01,4.20,4.37,4.50,4.62,4.72,4.82,
4.89,4.97,5.04,5.11,5.17,5.22,5.27,5.32, 5.37,
5.41,5.45,5.49,5.52,5.56,5.59,5.63,5.66,5.69)
m1[4,]<-c(3.09,3.39,3.62,3.79,3.93,4.04,4.14,4.23,4.31,
4.38,4.45,4.51,4.56,4.61,4.66,4.70,4.74,4.78,
4.82,4.85,4.89,4.92,4.95,4.98,5.00,5.03,5.06)
m1[5,]<-c(2.92,3.19,3.39,3.54,3.66,3.77,3.86,3.94,4.01,
4.07,4.13,4.18,4.23,4.28,4.32,4.36,4.39,4.43,
4.46,4.49,4.52,4.55,4.58,4.60,4.63,4.65,4.68)
m1[6,]<-c(2.80,3.06,3.24,3.38,3.49,3.59,3.67,3.74,3.80,
3.86,3.92,3.96,4.01,4.05,4.09,4.13,4.16,4.19,
4.22,4.25,4.28,4.31,4.33,4.35,4.38,4.39,4.42)
m1[7,]<-c(2.72,2.96,3.13,3.26,3.36,3.45,3.53,3.60,3.66,
3.71,3.76,3.81,3.85,3.89,3.93,3.96,3.99, 4.02,
4.05,4.08,4.10,4.13,4.15,4.18,4.19,4.22,4.24)
m1[8,]<-c(2.66,2.89,3.05,3.17,3.27,3.36,3.43,3.49,3.55,
3.60,3.65,3.69,3.73,3.77,3.80,3.84,3.87,3.89,
3.92,3.95,3.97,3.99,4.02,4.04,4.06,4.08,4.09)
m1[9,]<-c(2.61,2.83,2.98,3.10,3.19,3.28,3.35,3.41,3.47,
3.52,3.56,3.60,3.64,3.68,3.71,3.74,3.77,3.79,
3.82,3.85,3.87,3.89,3.91,3.94,3.95, 3.97,3.99)
m1[10,]<-c(2.57,2.78,2.93,3.05,3.14,3.22,3.29,3.35,3.40,
3.45,3.49,3.53,3.57,3.60,3.63,3.66,3.69,3.72,
3.74,3.77,3.79,3.81,3.83,3.85,3.87,3.89,3.91)
m1[11,]<-c(2.54,2.75,2.89,3.01,3.09,3.17,3.24,3.29,3.35,
3.39,3.43,3.47,3.51,3.54,3.57,3.60,3.63,3.65,
3.68,3.70,3.72,3.74,3.76,3.78,3.80,3.82,3.83)
m1[12,]<-c(2.49,2.69,2.83,2.94,3.02,3.09,3.16,3.21,3.26,
3.30,3.34,3.38,3.41,3.45,3.48,3.50,3.53,3.55,
3.58,3.59,3.62,3.64,3.66,3.68,3.69,3.71,3.73)
m1[13,]<-c(2.46,2.65,2.78,2.89,2.97,3.04,3.09,3.15,3.19,
3.24,3.28,3.31,3.35,3.38,3.40,3.43,3.46,3.48,
3.50,3.52,3.54,3.56,3.58,3.59,3.61,3.63,3.64)
m1[14,]<-c(2.43,2.62,2.75,2.85,2.93,2.99,3.05,3.11,3.15,
3.19,3.23,3.26,3.29,3.32,3.35,3.38,3.40,3.42,
3.44,3.46,3.48,3.50,3.52,3.54,3.55,3.57,3.58)
m1[15,]<-c(2.41,2.59,2.72,2.82,2.89,2.96,3.02,3.07,3.11,
3.15,3.19,3.22,3.25,3.28,3.31,3.33,3.36,3.38,
3.39,3.42,3.44,3.46,3.47,3.49,3.50,3.52,3.53)
m1[16,]<-c(2.38,2.56,2.68,2.77,2.85,2.91,2.97,3.02,3.06,
3.09,3.13,3.16,3.19,3.22,3.25,3.27,3.29,3.31,
3.33,3.35,3.37,3.39,3.40,3.42,3.43,3.45,3.46)
m1[17,]<-c(2.35,2.52,2.64,2.73,2.80,2.87,2.92,2.96,3.01,
3.04,3.07,3.11,3.13,3.16,3.18,3.21,3.23,3.25,
3.27,3.29,3.30,3.32,3.33,3.35,3.36,3.37,3.39)
m1[18,]<-c(2.32,2.49,2.60,2.69,2.76,2.82,2.87,2.91,2.95,
2.99,3.02,3.05,3.08,3.09,3.12,3.14,3.17, 3.18,
3.20,3.22,3.24,3.25,3.27,3.28,3.29,3.31,3.32)
m1[19,]<-c(2.29,2.45,2.56,2.65,2.72,2.77,2.82,2.86,2.90,
2.93,2.96,2.99,3.02,3.04,3.06,3.08,3.10, 3.12,
3.14,3.16,3.17,3.19,3.20,3.21,3.23,3.24,3.25)
m1[20,]<-c(2.24,2.39,2.49,2.57,2.63,2.68,2.73,2.77,2.79,
2.83,2.86,2.88,2.91,2.93,2.95,2.97,2.98, 3.01,
3.02,3.03,3.04,3.06,3.07,3.08,3.09,3.11,3.12)
if(nuhat>=200)smmcrit<-m1[20,C]
if(nuhat<200){
nu<-c(2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,24,30,40,60,200)
temp<-abs(nu-nuhat)
find<-order(temp)
if(temp[find[1]]==0)smmcrit<-m1[find[1],C]
if(temp[find[1]]!=0){
if(nuhat>nu[find[1]]){
smmcrit<-m1[find[1],C]-
(1/nu[find[1]]-1/nuhat)*(m1[find[1],C]-m1[find[1]+1,C])/
(1/nu[find[1]]-1/nu[find[1]+1])
}
if(nuhat<nu[find[1]]){
smmcrit<-m1[find[1]-1,C]-
(1/nu[find[1]-1]-1/nuhat)*(m1[find[1]-1,C]-m1[find[1],C])/
(1/nu[find[1]-1]-1/nu[find[1]])
}
}}
}
smmcrit
}


# depth
depth<-function(U,V,m){
#
#  Compute the halfspace depth of the point (u,v) for the pairs of points
#  in the n by 2 matrix m.
#
X<-m[,1]
Y<-m[,2]
FV<-NA
NUMS<-0
NUMH<-0
SDEP<-0.0
HDEP<-0.0
N<-length(X)
P<-acos(-1)
P2<-P*2.0
EPS<-0.000001
ALPHA<-NA
NT<-0
for(i in 1:nrow(m)){
               DV<-sqrt(((X[i]-U)*(X[i]-U)+(Y[i]-V)*(Y[i]-V)))
              if (DV <= EPS){
              NT<-NT+1
                              }
          else{
              XU<-(X[i]-U)/DV
              YU<-(Y[i]-V)/DV
              if (abs(XU) > abs(YU)){
                  if (X[i] >= U){
                      ALPHA[i-NT]<-asin(YU)
                      if(ALPHA[i-NT] < 0.0)
                          ALPHA[i-NT]<-P2+ALPHA[i-NT]
                                  }
                  else{
                      ALPHA[i-NT]<-P-asin(YU)
                      }
                                    }
              else{
                  if (Y[i] >= V)
                      ALPHA[i-NT]<-acos(XU)
                  else
                      ALPHA[i-NT]<-P2-acos(XU)
                  }
              if (ALPHA[i-NT] >= P2-EPS) ALPHA[i-NT]<-0.0
                }
}
NN<-N-NT
if(NN<=1){
NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
depths1(NT,3)
      if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
      NUMH<-NUMH+NT
      HDEP<-(NUMH+0.0)/(N+0.0)
      return(HDEP)
}
ALPHA<-sort(ALPHA[1:NN])
ANGLE<-ALPHA[1]-ALPHA[NN]+P2
for(i in 2:NN){
ANGLE<-max(c(ANGLE,ALPHA[i]-ALPHA[i-1]))
               }
if(ANGLE > (P+EPS)){
NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
depths1(NT,3)
      if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
      NUMH<-NUMH+NT
      HDEP<-(NUMH+0.0)/(N+0.0)
      return(HDEP)
                  }
ANGLE<-ALPHA[1]
NU<-0
for (i in 1:NN){
ALPHA[i]<-ALPHA[i]-ANGLE
if(ALPHA[i]<(P-EPS))NU<-NU+1
               }
if(NU >= NN){
NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
depths1(NT,3)
      if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
      NUMH<-NUMH+NT
      HDEP<-(NUMH+0.0)/(N+0.0)
      return(HDEP)
}
#
#  Mergesort the alpha with their antipodal angles beta,
#  and at the same time update I, F(I), and NBAD.
#
JA<-1
JB<-1
      ALPHK<-ALPHA[1]
      BETAK<-ALPHA[NU+1]-P
      NN2<-NN*2
      NBAD<-0
      I<-NU
      NF<-NN
for(J in 1:NN2){
           ADD<-ALPHK+EPS
          if (ADD < BETAK){
              NF<-NF+1
              if(JA < NN){
                  JA<-JA+1
                  ALPHK<-ALPHA[JA]
              }
              else
                  ALPHK<-P2+1.0
              }
          else{
              I<-I+1
              NN1<-NN+1
              if(I==NN1){
                  I<-1
                  NF<-NF-NN
              }
              FV[I]<-NF
              NFI<-NF-I
              NBAD<-NBAD+depths1(NFI,2)
              if(JB < NN){
                  JB<-JB+1
                  if(JB+NU <= NN)
                      BETAK<-ALPHA[JB+NU]-P
                  else
                      BETAK<-ALPHA[JB+NU-NN]+P
              }
              else
                  BETAK<-P2+1.0
          }
}
NUMS<-depths1(NN,3)-NBAD
#
#  Computation of NUMH for halfspace depth.
#
      GI<-0
      JA<-1
      ANGLE<-ALPHA[1]
      dif<-NN-FV[1]
      NUMH<-min(FV[1],dif)
for(I in 2:NN){
          AEPS<-ANGLE+EPS
          if(ALPHA[I] <= AEPS){
              JA<-JA+1
                              }
          else{
              GI<-GI+JA
              JA<-1
              ANGLE<-ALPHA[I]
              }
          KI<-FV[I]-GI
          NNKI<-NN-KI
          NUMH<-min(c(NUMH,min(c(KI,NNKI))))
   }
NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
depths1(NT,3)
      if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
      NUMH<-NUMH+NT
      HDEP<-(NUMH+0.0)/(N+0.0)
      HDEP
}


# lplot.pred
lplot.pred<-function(x,y,pts=NULL,xout=FALSE,outfun=outpro,span=2/3,family='gaussian',...){
#
#  Using loess, compute predicted values based on the data in  pts
#
x<-as.matrix(x)
d=ncol(x)
dp1=d+1
m<-elimna(cbind(x,y))
n.orig=nrow(m)
n.keep=n.orig
if(xout){
flag<-outfun(m[,1:d],plotit=FALSE,...)$keep
m<-m[flag,]
n.keep=nrow(m)
}
x<-m[,1:d]
y<-m[,dp1]
if(is.null(pts))pts=x
fit=loess(y~x,span=span,family=family)
pred=predict(fit,pts)
list(n=n.orig,n.keep=n.keep,x.used=x,yhat=pred)
}


# trimci
trimci<-function(x,tr=.2,alpha=.05,null.value=0,pr=TRUE,nullval=NULL){
#
#  Compute a 1-alpha confidence interval for the trimmed mean
#
#  The default amount of trimming is tr=.2
#
if(pr){
print("The p-value returned by this function is based on the")
print("null value specified by the argument null.value, which defaults to 0")
print('To get a measure of effect size using a Winsorized measure of scale,  use trimciv2')
}
if(!is.null(nullval))null.value=nullval
x<-elimna(x)
se<-sqrt(winvar(x,tr))/((1-2*tr)*sqrt(length(x)))
trimci<-vector(mode="numeric",length=2)
df<-length(x)-2*floor(tr*length(x))-1
trimci[1]<-mean(x,tr)-qt(1-alpha/2,df)*se
trimci[2]<-mean(x,tr)+qt(1-alpha/2,df)*se
test<-(mean(x,tr)-null.value)/se
sig<-2*(1-pt(abs(test),df))
list(estimate=mean(x,tr),ci=trimci,test.stat=test,se=se,p.value=sig,n=length(x))
}


# standm
standm<-function(x,locfun=lloc,est=mean,scat=var,...){
# standardize a matrix x
#
x=elimna(x)
x=as.matrix(x)
m1=lloc(x,est=est)
v1=apply(x,2,scat)
p=ncol(x)
for(j in 1:p)x[,j]=(x[,j]-m1[j])/sqrt(v1[j])
x
}


# con3way
con3way<-function(J,K,L){
#
# Generate all contrast coefficients for 3-way design.
# with the goal of testing all main effects and interactions
#
cj=cjMAT(J)
ck=cjMAT(K)
cl=cjMAT(L)
ij<-matrix(c(rep(1,J)),1,J)
ik<-matrix(c(rep(1,K)),1,K)
il<-matrix(c(rep(1,L)),1,L)
conA=t(kron(cj,kron(ik,il)))
conB=t(kron(ij,kron(ck,il)))
conC=t(kron(ij,kron(ik,cl)))
conAB=kron(cj,kron(ck,il))
conAC=kron(cj,kron(ik,cl))
conBC=kron(ij,kron(ck,cl))
conABC=kron(cj,kron(ck,cl))
list(conA=conA,conB=conB,conC=conC,conAB=t(conAB),conAC=t(conAC),
conBC=t(conBC),conABC=t(conABC))
}



# binmat
binmat<-function(m,col,lower, upper,INC=TRUE){
#
# pull out the rows of the matrix m based on the values in the column
# indicated by the argument
# col
# that are between lower and upper, inclusive.  Note: the built-in function findInterval could be used instead
#
#  Example: binmat(m,3,10,15) will return all rows such that the
#  values in column 3 are between 10 and 15, inclusive.
#
if(is.vector(m)){
m=as.matrix(m)
col=1
}
if(INC){
flag1=m[,col]<=upper
flag2=m[,col]>=lower
}
if(!INC){
flag1=m[,col]<upper
flag2=m[,col]>lower
}
flag=as.logical(flag1*flag2)
m[flag,]
}



# near3d
near3d<-function(x,pt,fr=.8,m){
# determine which values in x are near pt
# based on fr * cov.mve
#
# x is assumed to be an n by p matrix
# pt is a vector of length p (a point in p-space).
# m is cov.mve(x) computed by runm3d
#
library(MASS)
if(!is.matrix(x))stop("Data are not stored in a matrix.")
dis<-sqrt(mahalanobis(x,pt,m$cov))
dflag<-dis < fr
dflag
}


# kron
kron<-function(m1,m2){
#  compute the Kronecker product of the two matrices m1 and m2.
#
m1<-as.matrix(m1) # Vectors of length p are converted to a p by 1 matrix
m2<-as.matrix(m2)
kron<-vector(mode="numeric",length=0)
for(i in 1:nrow(m1)){
m3<-m1[i,1]*m2
for(j in 2:ncol(m1))m3<-cbind(m3,m1[i,j]*m2)
if(i==1)kron<-m3
if(i>=2)kron<-rbind(kron,m3)
}
kron
}


# fac2list
fac2list<-function(x,g,pr=TRUE){
#
# data are stored in x
# information about the level of the value in x is stored in g,
# which can be a matrix with up to 4 columns
#
# sort the data in x into groups based on values in g.
# store results in list mode.
#
#  Example: fac2list(m[,2],m[,4]) would sort the values
#  in column 2 of m according to the values in column 4 of m
#
g=as.data.frame(g)
ng=ncol(g)+1
xg=cbind(x,g)
xg=elimna(xg)
x=xg[,1]
x=as.matrix(x)
g=xg[,2:ng]
g=as.data.frame(g)
L=ncol(g)
g=listm(g)
for(j in 1:L)g[[j]]=as.factor(g[[j]])
g=matl(g)
Lp1=L+1
if(L>4)stop("Can have at most 4 factors")
if(L==1){
res=selby(cbind(x,g),2,1)
group.id=res$grpn
res=res$x
}
if(L>1){
res=selby2(cbind(x,g),c(2:Lp1),1)
group.id=res$grpn
res=res$x
}
if(pr)
{print("Group Levels:")
print(group.id)
}
res=lapply(res,as.numeric)
res
}


# smmcrit01
smmcrit01<-function(nuhat,C){
#
#  Determine the .99 quantile of the C-variate Studentized maximum
#  modulus distribution using linear interpolation on inverse
#  degrees of freedom
#  If C=1, this function returns the .995 quantile of Student's t
#  distribution.
#
if(C-round(C)!=0)stop("The number of contrasts, C, must be an  integer")
if(C>=29)stop("C must be less than or equal to 28")
if(C<=0)stop("C must be greater than or equal to 1")
if(nuhat<2)stop("The degrees of freedom must be greater than or equal to 2")
if(C==1)smmcrit01<-qt(.995,nuhat)
if(C>=2){
C<-C-1
m1<-matrix(0,20,27)
m1[1,]<-c(12.73,14.44,15.65,16.59,17.35,17.99,18.53,19.01,19.43,
19.81,20.15,20.46,20.75,20.99,20.99,20.99,20.99,20.99,
22.11,22.29,22.46,22.63,22.78,22.93,23.08,23.21,23.35)
m1[2,]<-c(7.13,7.91,8.48,8.92,9.28,9.58,9.84,10.06,10.27,
10.45,10.61,10.76,10.90,11.03,11.15,11.26,11.37,11.47,
11.56,11.65,11.74,11.82,11.89,11.97,12.07,12.11,12.17)
m1[3,]<-c(5.46,5.99,6.36,6.66,6.89,7.09,7.27,7.43,7.57,
7.69,7.80,7.91,8.01,8.09,8.17,8.25,8.32,8.39,
8.45,8.51,8.57,8.63,8.68,8.73,8.78,8.83,8.87)
m1[4,]<-c(4.70,5.11,5.39,5.63,5.81,5.97,6.11,6.23,6.33,
6.43,6.52,6.59,6.67,6.74,6.81,6.87,6.93,6.98,
7.03,7.08,7.13,7.17,7.21,7.25,7.29,7.33,7.36)
m1[5,]<-c(4.27,4.61,4.85,5.05,5.20,5.33,5.45,5.55,5.64,
5.72,5.79,5.86,5.93,5.99,6.04,6.09,6.14,6.18,
6.23,6.27,6.31,6.34,6.38,6.41,6.45,6.48,6.51)
m1[6,]<-c(3.99,4.29,4.51,4.68,4.81,4.93,5.03,5.12,5.19,
5.27,5.33,5.39,5.45,5.50,5.55,5.59,5.64,5.68,
5.72,5.75,5.79,5.82,5.85,5.88,5.91,5.94,5.96)
m1[7,]<-c(3.81,4.08,4.27,4.42,4.55,4.65,4.74,4.82,4.89,
4.96,5.02,5.07,5.12,5.17,5.21,5.25,5.29, 5.33,
5.36,5.39,5.43,5.45,5.48,5.51,5.54,5.56,5.59)
m1[8,]<-c(3.67,3.92,4.10,4.24,4.35,4.45,4.53,4.61,4.67,
4.73,4.79,4.84,4.88,4.92,4.96,5.01,5.04,5.07,
5.10,5.13,5.16,5.19,5.21,5.24,5.26,5.29,5.31)
m1[9,]<-c(3.57,3.80,3.97,4.09,4.20,4.29,4.37,4.44,4.50,
4.56,4.61,4.66,4.69,4.74,4.78,4.81,4.84,4.88,
4.91,4.93,4.96,4.99,5.01,5.03,5.06,5.08,5.09)
m1[10,]<-c(3.48,3.71,3.87,3.99,4.09,4.17,4.25,4.31,4.37,
4.42,4.47,4.51,4.55,4.59,4.63,4.66,4.69,4.72,
4.75,4.78,4.80,4.83,4.85,4.87,4.89,4.91,4.93)
m1[11,]<-c(3.42,3.63,3.78,3.89,.99,4.08,4.15,4.21,4.26,
4.31,4.36,4.40,4.44,4.48,4.51,4.54,4.57,4.59,
4.62,4.65,4.67,4.69,4.72,4.74,4.76,4.78,4.79)
m1[12,]<-c(3.32,3.52,3.66,3.77,3.85,3.93,3.99,.05,4.10,
4.15,4.19,4.23,4.26,4.29,4.33,4.36,4.39,4.41,
4.44,4.46,4.48,4.50,4.52,4.54,4.56,4.58,4.59)
m1[13,]<-c(3.25,3.43,3.57,3.67,3.75,3.82,3.88,3.94,3.99,
4.03,4.07,4.11,4.14,4.17,4.19,4.23,4.25,4.28,
4.29,4.32,4.34,4.36,4.38,4.39,4.42,4.43,4.45)
m1[14,]<-c(3.19,3.37,3.49,3.59,3.68,3.74,3.80,3.85,3.89,
3.94,3.98,4.01,4.04,4.07,4.10,4.13,4.15,4.18,
4.19,4.22,4.24,4.26,4.28,4.29,4.31,4.33,4.34)
m1[15,]<-c(3.15,3.32,3.45,3.54,3.62,3.68,3.74,3.79,3.83,
3.87,3.91,3.94,3.97,3.99,4.03,4.05,4.07,4.09,
4.12,4.14,4.16,4.17,4.19,4.21,4.22,4.24,4.25)
m1[16,]<-c(3.09,3.25,3.37,3.46,3.53,3.59,3.64,3.69,3.73,
3.77,3.80,3.83,3.86,3.89,3.91,3.94,3.96,3.98,
4.00,4.02,4.04,4.05,4.07,4.09,4.10,4.12,4.13)
m1[17,]<-c(3.03,3.18,3.29,3.38,3.45,3.50,3.55,3.59,3.64,
3.67,3.70,3.73,3.76,3.78,3.81,3.83,3.85,3.87,
3.89,3.91,3.92,3.94,3.95,3.97,3.98,4.00,4.01)
m1[18,]<-c(2.97,3.12,3.22,3.30,3.37,3.42,3.47,3.51,3.55,
3.58,3.61,3.64,3.66,3.68,3.71,3.73,3.75,3.76,
3.78,3.80,3.81,3.83,3.84,3.85,3.87,3.88,3.89)
m1[19,]<-c(2.91,3.06,3.15,3.23,3.29,3.34,3.38,3.42,3.46,
3.49,3.51,3.54,3.56,3.59,3.61,3.63,3.64,3.66,
3.68,3.69,3.71,3.72,3.73,3.75,3.76,3.77,3.78)
m1[20,]<-c(2.81,2.93,3.02,3.09,3.14,3.19,3.23,3.26,3.29,
3.32,3.34,3.36,3.38,3.40,.42,.44,3.45,3.47,
3.48,3.49,3.50,3.52,3.53,3.54,3.55,3.56,3.57)
if(nuhat>=200)smmcrit01<-m1[20,C]
if(nuhat<200){
nu<-c(2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,24,30,40,60,200)
temp<-abs(nu-nuhat)
find<-order(temp)
if(temp[find[1]]==0)smmcrit01<-m1[find[1],C]
if(temp[find[1]]!=0){
if(nuhat>nu[find[1]]){
smmcrit01<-m1[find[1],C]-
(1/nu[find[1]]-1/nuhat)*(m1[find[1],C]-m1[find[1]+1,C])/
(1/nu[find[1]]-1/nu[find[1]+1])
}
if(nuhat<nu[find[1]]){
smmcrit01<-m1[find[1]-1,C]-
(1/nu[find[1]-1]-1/nuhat)*(m1[find[1]-1,C]-m1[find[1],C])/
(1/nu[find[1]-1]-1/nu[find[1]])
}
}}
}
smmcrit01
}


# chi.int
chi.int <- function(p,a,c1)
#   partial expectation d in (0,c1) of d^a under chi-squared p
  return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*pchisq(c1^2,p+a) )

# psi.bt
psi.bt <- function(x,c1,M)
{
    x1 <- (x-M)/c1
    ivec1 <- (x1 < 0)
    ivec2 <- (x1 >  1)
    return(ivec1*x+(1-ivec1-ivec2)*x*(1-x1^2)^2)
}

# rdplot
rdplot<-function(x,fr=NA,plotit=TRUE,theta=50,phi=25,expand=.5,pyhat=FALSE,pts=NA,
xlab="X",ylab="",ticktype="simple"){
#
# Expected frequency curve
#
# fr controls amount of smoothing
#  theta is the azimuthal direction and phi the colatitude
#
plotit<-as.logical(plotit)
x<-elimna(x)
x<-as.matrix(x)
rmd<-NA
if(ncol(x)==1){
x=as.vector(x)
if(is.na(fr))fr<-.8
if(is.na(pts[1]))pts<-x
for(i in 1:length(pts)){
rmd[i]<-sum(near(x,pts[i],fr))
}
if(mad(x)!=0)rmd<-rmd/(2*fr*mad(x))
rmd<-rmd/length(x)
if(plotit){
plot(pts,rmd,type="n",ylab=ylab,xlab=xlab)
sx<-sort(pts)
xorder<-order(pts)
sysm<-rmd[xorder]
lines(sx,sysm)
}}
x<-as.matrix(x)
if(ncol(x)>1){
library(MASS)
if(is.na(fr))fr<-.6
m<-covmve(x)
for(i in 1:nrow(x)){
rmd[i]<-sum(near3d(x,x[i,],fr,m))
}
rmd<-rmd/nrow(x)
if(plotit && ncol(x)==2){
library(akima)
fitr<-rmd
iout<-c(1:length(fitr))
nm1<-length(fitr)-1
for(i in 1:nm1){
ip1<-i+1
for(k in ip1:length(fitr))if(sum(x[i,]==x[k,])==2)iout[k]<-0
}
fitr<-fitr[iout>=1]
mkeep<-x[iout>=1,]
fit<-interp(mkeep[,1],mkeep[,2],fitr)
persp(fit,theta=theta,phi=phi,expand=expand,xlab="Var 1",ylab="Var 2",zlab="",
ticktype=ticktype)
}
}
if(pyhat)last<-rmd
if(!pyhat)last<-"Done"
last
}

 rimul<-function(J,K,x,alpha=.05,p=J*K,grp=c(1:p),plotit=TRUE,op=4){
#
#  Rank-based multiple comparisons for all interactions
#  in J by K design. The method is based on an
#  extension of Cliff's heteroscedastic technique for
#  handling tied values and the Patel-Hoel definition of no interaction.
#
#  The familywise type I error probability is controlled by using
#  a critical value from the Studentized maximum modulus distribution.
#
#  It is assumed all groups are independent.
#
#  Missing values are automatically removed.
#
#  The default value for alpha is .05. Any other value results in using
#  alpha=.01.
#
#  Argument grp can be used to rearrange the order of the data.
#
 df=Inf
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
CCJ<-(J^2-J)/2
CCK<-(K^2-K)/2
CC<-CCJ*CCK
test<-matrix(NA,CC,8)
test.p<-matrix(NA,CC,7)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
}
mat<-matrix(grp,ncol=K,byrow=TRUE)
dimnames(test)<-list(NULL,c("Factor A","Factor A","Factor B","Factor B","delta","ci.lower","ci.upper","p.value"))
jcom<-0
crit=qsmm(1-alpha,CC,df)
#if(alpha!=.05)crit<-smmcrit01(200,CC)
alpha<-1-pnorm(crit)
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
temp1<-cid(x[[mat[j,k]]],x[[mat[j,kk]]],plotit=FALSE)
temp2<-cid(x[[mat[jj,k]]],x[[mat[jj,kk]]],plotit=FALSE)
delta<-temp2$d-temp1$d
sqse<-temp1$sqse.d+temp2$sqse.d
test[jcom,5]<-delta/2
test[jcom,6]<-delta/2-crit*sqrt(sqse/4)
test[jcom,7]<-delta/2+crit*sqrt(sqse/4)
test[jcom,8]=2*(1-pnorm(abs((delta/2)/sqrt(sqse/4))))
}}}}}}
if(J==2 & K==2){
if(plotit){
m1<-outer(x[[1]],x[[2]],FUN="-")
m2<-outer(x[[3]],x[[4]],FUN="-")
m1<-as.vector(m1)
m2<-as.vector(m2)
g2plot(m1,m2,op=op)
}}
list(test=test)
}


# binom.conf
binom.conf<-function(x = sum(y), nn = length(y),AUTO=TRUE,pr=TRUE,
method=c('AC','P','CP','KMS','WIL','SD'), y = NULL, n = NA, alpha = 0.05){
#
#
# P: Pratt's method
# AC: Agresti--Coull
# CP: Clopper--Pearson
# KMS:  Kulinskaya et al. 2008, p. 140
#  WIL:  Wilson type CI. Included for completeness; was used in simulations  relevant to binom2g
#   SD: Schilling--Doi
#
if(pr) print('Note: To perform the sign test, use the the R function signt')
if(nn<35){
if(AUTO)method='SD'
}
type=match.arg(method)
switch(type,
    P=binomci(x=x,nn=nn,y=y,n=n,alpha=alpha),
    AC=acbinomci(x=x,nn=nn,y=y,n=n,alpha=alpha),
    CP=binomCP(x=x,nn=nn,y=y,n=n,alpha=alpha),
    KMS=kmsbinomci(x=x,nn=nn,y=y,n=n,alpha=alpha),
    WIL=wilbinomci(x=x,y=y,n=nn,alpha=alpha),
    SD=binomLCO(x=x,nn=nn,y=y,alpha=alpha),
    )
}


# pdisMC
pdisMC<-function(m,MM=FALSE,cop=3,dop=1,center=NA){
#
# Compute projection distances for points in m
#
#
#
#  MM=F  Projected distance scaled
#  using interquatile range.
#  MM=T  Scale projected distances using MAD.
#
#  There are five options for computing the center of the
#  cloud of points when computing projections:
#  cop=1 uses Donoho-Gasko median
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#  cop=5 uses skipped mean
#
library(parallel)
m<-elimna(m) # Remove missing values
m<-as.matrix(m)
if(ncol(m)==1){
if(is.na(center[1]))center<-median(m)
dis<-abs(m[,1]-center)
if(!MM){
temp<-idealf(dis)
pdis<-dis/(temp$qu-temp$ql)
}
if(MM)pdis<-dis/mad(dis)
}
if(ncol(m)>1){
if(is.na(center[1])){
if(cop==1)center<-dmean(m,tr=.5,dop=dop)
if(cop==2)center<-cov.mcd(m,print=FALSE)$center
if(cop==3)center<-apply(m,2,median)
if(cop==4)center<-cov.mve(m,print=FALSE)$center
if(cop==5)center<-smean(m)
}
cenmat=matrix(rep(center,nrow(m)),ncol=ncol(m),byrow=TRUE)
Amat=m-cenmat
B=listm(t(Amat))  # so rows are now in B[[1]]...B[[n]]
dis=mclapply(B,outproMC.sub,Amat,mc.preschedule=TRUE)
if(!MM){
dmat<-mclapply(dis,IQRstand,mc.preschedule=TRUE)
}
if(MM)dmat<-mclapply(dis,MADstand,mc.preschedule=TRUE)
pdis<-apply(matl(dmat),1,max,na.rm=TRUE)
}
pdis
}

# wincor
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


# winall
winall<-function(m,tr=.2){
#
#    Compute the Winsorized correlation and covariance matrix for the
#    data in the n by p matrix m.
#
#    This function also returns the two-sided significance level
#
if(is.data.frame(m))m=as.matrix(m)
if(!is.matrix(m))stop("The data must be stored in a n by p matrix")
wcor<-matrix(1,ncol(m),ncol(m))
wcov<-matrix(0,ncol(m),ncol(m))
siglevel<-matrix(NA,ncol(m),ncol(m))
for (i in 1:ncol(m)){
ip<-i
for (j in ip:ncol(m)){
val<-wincor(m[,i],m[,j],tr)
wcor[i,j]<-val$cor
wcor[j,i]<-wcor[i,j]
if(i==j)wcor[i,j]<-1
wcov[i,j]<-val$cov
wcov[j,i]<-wcov[i,j]
if(i!=j){
siglevel[i,j]<-val$p.value
siglevel[j,i]<-siglevel[i,j]
}
}
}
cent=apply(m,2,mean,tr,na.rm=TRUE)
list(cor=wcor,cov=wcov,center=cent,p.values=siglevel)
}




# erho.bt
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

# yuend
yuend<-function(x,y,tr=.2,alpha=.05){
#
#  Compare the trimmed means of two dependent random variables
#  using the data in x and y.
#  The default amount of trimming is 20%
#
#  Any pair with a missing value is eliminated
#  The function rm2miss allows missing values.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuend$ci.
#  The significance level is returned in yuend$p.value
#
#  For inferences based on difference scores, use trimci
#
if(length(x)!=length(y))stop("The number of observations must be equal")
m<-cbind(x,y)
m<-elimna(m)
x<-m[,1]
y<-m[,2]
h1<-length(x)-2*floor(tr*length(x))
q1<-(length(x)-1)*winvar(x,tr)
q2<-(length(y)-1)*winvar(y,tr)
q3<-(length(x)-1)*wincor(x,y,tr)$cov
df<-h1-1
se<-sqrt((q1+q2-2*q3)/(h1*(h1-1)))
crit<-qt(1-alpha/2,df)
dif<-mean(x,tr)-mean(y,tr)
low<-dif-crit*se
up<-dif+crit*se
test<-dif/se
yuend<-2*(1-pt(abs(test),df))
list(ci=c(low,up),p.value=yuend,est1=mean(x,tr),est2=mean(y,tr),dif=dif,se=se,teststat=test,n=length(x),df=df)
}



# fdepth
fdepth<-function(m,pts=NA,plotit=TRUE,cop=3,center=NA,xlab="VAR 1",
ylab="VAR 2"){
#
# Determine depth of points in pts,  relative to
# points in m. If pts is not specified,
# depth of all points in m are determined.
#
# m and pts can be vectors or matrices with
# p columns (the number of variables).
#
# Determine center, for each point, draw a line
# connecting it with center, project points onto this line
# and determine depth of the projected points.
# The final depth of a point is its minimum depth
# among all projections.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data and pts=NA
#
#  There are three options for computing the center of the
#  cloud of points when computing projections, assuming center=NA:
#
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#
#  If a value for center is passed to this function,
#  this value is used to determine depths.
#
#  When plotting,
#  center is marked with a cross, +.
#
library(MASS)
if(cop!=2 && cop!=3 && cop!=4)stop("Only cop=2, 3 or 4 is allowed")
if(is.list(m))stop("Store data in a matrix; might use function listm")
m<-as.matrix(m)
pts<-as.matrix(pts)
if(!is.na(pts[1]))remm<-m
nm<-nrow(m)
nm1<-nm+1
if(!is.na(pts[1])){
if(ncol(m)!=ncol(pts))stop("Number of columns of m is not equal to number of columns for pts")
}
m<-elimna(m) # Remove missing values
m<-as.matrix(m)
if(ncol(m)==1)dep<-unidepth(as.vector(m[,1]),pts=pts)
if(ncol(m)>1){
if(is.na(center[1])){
if(cop==2){
center<-cov.mcd(m)$center
}
if(cop==4){
center<-cov.mve(m)$center
}
if(cop==3){
center<-apply(m,2,median)
}}
if(is.na(pts[1])){
mdep <- matrix(NA,nrow=nrow(m),ncol=nrow(m))
}
if(!is.na(pts[1])){
mdep <- matrix(NA,nrow=nrow(m),ncol=nrow(pts))
}
for (i in 1:nrow(m)){
B<-m[i,]-center
dis<-NA
BB<-B^2
bot<-sum(BB)
if(bot!=0){
if(is.na(pts[1])){
for (j in 1:nrow(m)){
A<-m[j,]-center
temp<-sum(A*B)*B/bot
dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
}}
if(!is.na(pts[1])){
m<-rbind(remm,pts)
for (j in 1:nrow(m)){
A<-m[j,]-center
temp<-sum(A*B)*B/bot
dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
}}
#
# For ith projection, store depths of
# points in mdep[i,]
#
if(is.na(pts[1]))mdep[i,]<-unidepth(dis)
if(!is.na(pts[1])){
mdep[i,]<-unidepth(dis[1:nm],dis[nm1:nrow(m)])
}}
if(bot==0)mdep[i,]<-rep(0,ncol(mdep))
}
dep<-apply(mdep,2,min)
if(ncol(m)==2 && is.na(pts[1])){
flag<-chull(m)
dep[flag]<-min(dep)
}
}
if(ncol(m)==2){
if(is.na(pts[1]) && plotit){
plot(m,xlab=xlab,ylab=ylab)
points(center[1],center[2],pch="+")
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
dep<-round(dep*nrow(m))/nrow(m)
dep
}


# ols
ols<-function(x,y,xout=FALSE,outfun=outpro,alpha=.05,plotit=FALSE,xlab='X',ylab='Y',zlab='Z',RES=TRUE,...){
#
# Performs OLS regression calling built-in R function.
#
# xout=T will eliminate any leverage points (outliers among x values)
# if one predictor,
# plotit=TRUE will plot the points and the regression line
#
m<-elimna(cbind(x,y))
n=nrow(m)
n.keep=n
x<-as.matrix(x)
p<-ncol(x)
pp<-p+1
x<-m[,1:p]
y<-m[,pp]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,pp]
n.keep=length(y)
}
x<-as.matrix(x)
temp<-summary(lm(y~x))
coef<-temp[4]$coefficients
CI=matrix(NA,nrow(coef),ncol=2)
CI[,1]=coef[,1]-qt(1-alpha/2,temp[10]$fstatistic[3])*coef[,2]
CI[,2]=coef[,1]+qt(1-alpha/2,temp[10]$fstatistic[3])*coef[,2]
dimnames(CI)=list(NULL,c("low.ci","up.ci"))
coef=cbind(coef,CI)
if(plotit){
if(p==1){
plot(x,y,xlab=xlab,ylab=ylab)
abline(coef[,1])
}
if(p==2){
regp2plot(x,y,regfun=ols,xlab=xlab,ylab=ylab,zlab=zlab)
}}
Ftest<-temp[10]$fstatistic
Ftest.p.value<-1-pf(Ftest[1],Ftest[2],Ftest[3])
Rval=Rsq.ols(x,y)
res=NULL
if(RES)res=y-x%*%coef[2:pp,1]-coef[1,1]
list(n=n,n.keep=n.keep,summary=coef,coef=coef[,1],F.test=temp[10]$fstatistic[1],Ftest.p.value=Ftest.p.value,
F.test.degrees.of.freedom=temp[10]$fstatistic[2:3],R.squared=Rval,residuals=as.vector(res))
}


# rplot
rplot<-function(x,y,est=tmean,scat=TRUE,fr=NA,plotit=TRUE,pyhat=FALSE,efr=.5,
theta=50,phi=25,scale=TRUE,expand=.5,SEED=TRUE,varfun=pbvar,outfun=outpro,
nmin=0,xout=FALSE,out=FALSE,eout=FALSE,xlab='X',ylab='Y',zscale=FALSE,
zlab=' ',pr=TRUE,duplicate='error',ticktype='simple',LP=TRUE,OLD=FALSE,pch='.',prm=TRUE,...){
# duplicate='error'
# In some situations where duplicate values occur, when plotting with
# two predictors, it is necessary to set duplicate='strip'
#
# LP=TRUE, the plot of the smooth is further smoothed via lplot (lowess)
# To get a plot as done with old version set
# LP=FALSE
#
#  zscale=TRUE will standardize the dependent variable when plotting with 2 independent variables.
#
# efr is the span when computing explanatory strength of association
#
# cf qplot in the R package ggplot2
#
if(pr){
if(!xout)print('Suggest also looking at result using xout=TRUE')
}
x<-as.matrix(x)
p=ncol(x)
p1=p+1
if(pr && !OLD){
print('A new estimate of the strength of the association is used by default.')
print(' To get the old estimate, set OLD=TRUE')
}
xx<-cbind(x,y)
xx<-elimna(xx)
n=nrow(xx)
if(eout){
flag=outfun(xx,plotit=FALSE,...)$keep
xx=xx[flag,]
}
if(xout){
if(identical(outfun,outblp))flag=outblp(xx[,1:p],xx[,p1],plotit=FALSE)$keep
else
flag=outfun(xx[,1:p],plotit=FALSE,...)$keep
xx=xx[flag,]
}
n.keep=nrow(xx)
x<-xx[,1:p]
x<-as.matrix(x)
p1=ncol(x)+1
y<-xx[,p1]
if(ncol(x)==1){
if(is.na(fr))fr<-.8
val<-rungen(x,y,est=est,scat=scat,fr=fr,plotit=plotit,pyhat=TRUE,
xlab=xlab,ylab=ylab,LP=LP,pch=pch,...)
val2<-rungen(x,y,est=est,fr=efr,plotit=FALSE,pyhat=TRUE,LP=FALSE,...)$output
val<-val$output
}
if(ncol(x)>1){
id=chk4binary(x)
Lid=length(id)
if(Lid>0)Stop('Binary independent variables detected, use rplotv2')
if(ncol(x)==2 && !scale){
if(pr){print('scale=FALSE is specified.')
print('If there is dependence, might want to use scale=T')
}}
if(is.na(fr))fr<-1
val<-rung3d(x,y,est=est,fr=fr,plotit=plotit,pyhat=TRUE,SEED=SEED,nmin=nmin,LP=LP,
scale=scale,phi=phi,theta=theta,expand=expand,zscale=zscale,pr=FALSE,
duplicate='error',xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,...)
}
E.power=NULL
if(OLD){
E.power=varfun(val)/varfun(y)
names(E.power)=''
if(E.power>1)E.power=.99
}
if(!OLD)E.power=smRstr(x,y,fr=fr)$str^2
stra=sqrt(E.power)
# Best correction at the moment. Not sure when or if needed.
# Maybe a correlation option is better, but need to check this.
xvals=x
if(ncol(x)==1)xvals=sort(xvals)
if(!pyhat){
val <- NULL
xvals=NULL
}
if(!prm){
stra=NULL
E.power=NULL
val=NULL
}
list(n=n,n.keep=n.keep,Strength.Assoc=stra,Explanatory.Power = E.power,  xvals=xvals,yhat = val)
}


# out
out<-function(x,cov.fun=cov.mve,xlab="X",ylab="Y",qval=.975,
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
#  Possible choices are
#  cov.mve (the MVE estimate)
#  cov.mcd (the MCD estimate)
#  covmba2 (the MBA or median ball algorithm)
#  rmba  (an adjustment of MBA suggested by D. Olive)
#  cov.roc (Rocke's TBS estimator)
#
#  plotit=FALSE used to avoid problems when other functions in WRS call
#  this function
#
#  KS=TRUE: keep  the seed that was used
#
library(MASS)
if(KS)oldSeed <- .Random.seed
set.seed(12)
if(is.data.frame(x))x=as.matrix(x)
if(is.list(x))stop("Data cannot be stored in list mode")
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
plot(x[,1],x[,2],xlab=xlab,ylab=ylab,type="n")
flag<-rep(TRUE,nrow(x))
flag[id]<-FALSE
points(x[flag,1],x[flag,2])
if(sum(!flag)>0)points(x[!flag,1],x[!flag,2],pch="*")
}}
if(!is.matrix(x))outval<-x[id]
if(is.matrix(x))outval<-x[id,]
n=nrow(as.matrix(x))
n.out=length(id)
assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
list(n=n,n.out=n.out,out.val=outval,out.id=id,keep=keep,dis=dis,crit=crit)
}


# bmp
bmp<-function(x,y,alpha=.05,crit=NA,plotit=FALSE,pop=0,fr=.8,xlab='',ylab=''){
#
# Brunner and Munzel (2000) heteroscedastic analog of WMW test.
#
# plotit=TRUE causes a plot of the difference scores to be created
x<-x[!is.na(x)]  # Remove any missing values
y<-y[!is.na(y)]
n1<-length(x)
n2<-length(y)
N<-n1+n2
n1p1<-n1+1
flag1<-c(1:n1)
flag2<-c(n1p1:N)
R<-rank(c(x,y))
R1<-mean(R[flag1])
R2<-mean(R[flag2])
Rg1<-rank(x)
Rg2<-rank(y)
S1sq<-sum((R[flag1]-Rg1-R1+(n1+1)/2)^2)/(n1-1)
S2sq<-sum((R[flag2]-Rg2-R2+(n2+1)/2)^2)/(n2-1)
sig1<-S1sq/n2^2
sig2<-S2sq/n1^2
se<-sqrt(N)*sqrt(N*(sig1/n1+sig2/n2))
bmtest<-(R2-R1)/se
phat<-(R2-(n2+1)/2)/n1
flag=TRUE
if(phat==0 || phat==1)flag=FALSE
dhat<-1-2*phat
df<-(S1sq/n2 + S2sq/n1)^2/((S1sq/n2)^2/(n1-1)+(S2sq/n1)^2/(n2-1))
sig<-2 * (1 - pt(abs(bmtest),df))
if(is.na(crit))vv<-qt(alpha/2,df)
if(!is.na(crit))vv<-crit
ci.p<-c(phat+vv*se/N,phat-vv*se/N)
ci.p[1]=max(0,ci.p[1])
ci.p[2]=min(1,ci.p[2])
dval=matrix(0,1,3)
for(i in 1:n1){
for(j in 1:n2){
id=sign(x[i]-y[j])+2
dval[1,id]=dval[1,id]+1
}}
dval=dval/(n1*n2)
dimnames(dval)<-list(NULL,c('P(X<Y)','P(X=Y)','P(X>Y)'))
if(!flag){
nm=max(c(length(x),length(y)))
if(phat==1)A=binomcipv(nm,nm,alpha=alpha)
if(phat==0)A=binomcipv(0,nm,alpha=alpha)
ci.p=A$ci
sig=A$p.value
}

if(plotit){
msave<-outer(x,y,FUN='-')
if(pop==0){
if(length(x)*length(y)>2500){
print('Product of sample sizes exceeds 2500.')
print('Execution time might be high when plotting and when using pop=1')
print('If this is case, might consider changing the argument pop or using plotit=F')
}
akerd(as.vector(msave),fr=fr)
}
if(pop==1)rdplot(as.vector(msave),fr=fr,xlab=xlab,ylab=ylab)
if(pop==2)kdplot(as.vector(msave),rval=rval,xlab=xlab,ylab=ylab)
if(pop==3)boxplot(as.vector(msave))
if(pop==4)stem(as.vector(msave))
if(pop==5)hist(as.vector(msave))
if(pop==6)skerd(as.vector(msave),xlab=xlab,ylab=ylab)
}
list(n1=n1,n2=n2,test.stat=bmtest,phat=phat,dhat=dhat,s.e.=se/N,p.value=sig,ci.p=ci.p,df=df,summary.dval=dval)
}


# trimse
trimse<-function(x,tr=.2,na.rm=FALSE){
#
#  Estimate the standard error of the gamma trimmed mean
#  The default amount of trimming is tr=.2.
#
if(na.rm)x<-x[!is.na(x)]
trimse<-sqrt(winvar(x,tr))/((1-2*tr)*sqrt(length(x)))
trimse
}


# ancova
ancova<-function(x1,y1,x2,y2,fr1=1,fr2=1,tr=.2,alpha=.05,plotit=TRUE,pts=NA,sm=FALSE,method="EP",SEED=TRUE,
pr=TRUE,xout=FALSE,outfun=out,LP=FALSE,SCAT=TRUE,xlab='X',ylab='Y',pch1='*',pch2='+',
skip.crit=FALSE,nmin=12,crit.val=1.09,...){
#
# Compare two independent  groups using the ancova method with a single covariate
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#
#  sm=TRUE will create smooths using bootstrap bagging.
#  pts can be used to specify the design points where the regression lines
#  are to be compared.
#
#  Argument method indicates which measure of effect size will be used
#  EP: explanatory measure of effect size  (default)
#  QS: quantile shift measure of effect size
#  AKP:  trimmed mean Winsorized variance analog of Cohen's d
#  WMW:  P(X<Y)
#
#
#   skip.crit: this argument is used to avoid computational issues associated with ancdet.pv
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
if(length(x1)!=length(y1))stop("x1 and y1 have different lengths")
if(length(x2)!=length(y2))stop("x2 and y2 have different lengths")
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
if(pr){
print("NOTE: Confidence intervals are adjusted to control the probability")
print("of at least one Type I error.")
#print("But p-values are not")
print('Effect size is based on the argument method, default is explanatory measure of effect size')
print('Other options: QS, quantile shift; AKP, robust analog of Cohen d; WMW, P(X<Y)')
print('KMS, robust analog of Cohen d')
}
if(xout){
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
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
isub[1]<-min(sub[vecn>=nmin])
isub[5]<-max(sub[vecn>=nmin])
isub[3]<-floor((isub[1]+isub[5])/2)
isub[2]<-floor((isub[1]+isub[3])/2)
isub[4]<-floor((isub[3]+isub[5])/2)
mat<-matrix(NA,5,14)
dimnames(mat)<-list(NULL,c("X","n1","n2","Est1","Est2","DIF","TEST","se","ci.low","ci.hi","p.value","crit.val","Effect.Size",'p.adjusted'))
for (i in 1:5){
g1<-y1[near(x1,x1[isub[i]],fr1)]
g2<-y2[near(x2,x1[isub[i]],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test<-yuen(g1,g2,tr=tr)
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
mat[i,4]<-test$est.1
mat[i,5]<-test$est.2
mat[i,6]<-test$dif
mat[i,7]<-test$teststat
mat[i,8]<-test$se
mat[i,13]=ESfun(g1,g2,method=method,pr=FALSE,SEED=SEED)
mat[i,14]=1-psmm(abs(test$teststat),5,test$df)
if(skip.crit)critv=crit.val
#if(!skip.crit){
critv<-NA
#if(alpha==.05)critv<-smmcrit(test$df,5)
#if(alpha==.01)critv<-smmcrit01(test$df,5)
#if(is.na(critv))critv<-smmval(test$df,5,alpha=alpha)
critv=qsmm(1-alpha,5,test$df)
mat[i,12]<-critv
#}
cilow<-test$dif-critv*test$se
cihi<-test$dif+critv*test$se
mat[i,9]<-cilow
mat[i,10]<-cihi
mat[i,11]<-test$p.value
}}
if(!is.na(pts[1])){
#if(!skip.crit){
#if(length(pts)>=29)stop("At most 28 points can be compared")
#}
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),14)
dimnames(mat)<-list(NULL,c("X","n1","n2","Est1","Est2","DIF","TEST","se","ci.low","ci.hi",
"p.value","crit.val","Effect.Size",'p.adjusted'))
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test<-yuen(g1,g2,tr=tr)
mat[i,1]<-pts[i]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
if(length(g1)<=5)print(paste("Warning, there are",length(g1)," points corresponding to the design point X=",pts[i]))
if(length(g2)<=5)print(paste("Warning, there are",length(g2)," points corresponding to the design point X=",pts[i]))
mat[i,4]<-test$est.1
mat[i,5]<-test$est.2
mat[i,6]<-test$dif
mat[i,7]<-test$teststat
mat[i,8]<-test$se
mat[i,13]=ESfun(g1,g2,method=method,pr=FALSE,SEED=SEED)
mat[i,14]=1-psmm(abs(test$teststat),length(pts),test$df)
if(skip.crit)critv=crit.val
if(!skip.crit){
if(length(pts)>=2)critv=qsmm(1-alpha,length(pts),test$df)    #smmcrit(test$df,length(pts))
if(length(pts)==1)critv<-qt(.975,test$df)
}
cilow<-test$dif-critv*test$se
cihi<-test$dif+critv*test$se
mat[i,9]<-cilow
mat[i,10]<-cihi
mat[i,11]<-test$p.value
mat[i,12]<-critv
}}
if(plotit){
runmean2g(x1,y1,x2,y2,fr=fr1,est=tmean,tr=tr,sm=sm,xout=FALSE,LP=LP,
SCAT=SCAT,xlab=xlab,ylab=ylab,pch1=pch1,pch2=pch2,...)
}
list(output=mat)
}

# covmtrim
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

# olshc4
olshc4<-function(x,y,alpha=.05,CN=FALSE,
xout=FALSE,outfun=outpro,HC3=FALSE,plotit=FALSE,xlab = "X", ylab = "Y", zlab = "Z",...){
#
# Compute confidence intervals via least squares
# regression using heteroscedastic method
# recommended by Cribari-Neto (2004).
# CN=F, degrees of freedom are n-p
# CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
# All indications are that CN=F is best for general use.
#
#  HC3=TRUE, will replace the HC4 estimator with the HC3 estimator.
#
x<-as.matrix(x)
pnum=ncol(x)
if(nrow(x) != length(y))stop("Length of y does not match number of x values")
m<-cbind(x,y)
m<-elimna(m)
y<-m[,ncol(x)+1]
x=m[,1:ncol(x)]
n=length(y)
nrem=n
n.keep=length(y)
x<-as.matrix(x)
if(xout){
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
x<-as.matrix(x)
x<-x[flag,]
y<-y[flag]
n.keep=length(y)
x<-as.matrix(x)
}
temp<-lsfit(x,y)
rsq=Rsq.ols(x,y)
x<-cbind(rep(1,nrow(x)),x)
xtx<-solve(t(x)%*%x)
h<-diag(x%*%xtx%*%t(x))
n<-length(h)
d<-(n*h)/sum(h)
for(i in 1:length(d)){
        d[i]<-min(4, d[i])
}
if(HC3)d=2
hc4<-xtx%*%t(x)%*%diag(temp$res^2/(1-h)^d)%*%x%*%xtx
df<-nrow(x)-ncol(x)
crit<-qt(1-alpha/2,df)
if(CN)crit=qnorm(1-alpha/2)
al<-ncol(x)
p=al-1
ci<-matrix(NA,nrow=al,ncol=6)
lab.out=rep("Slope",p)
dimnames(ci)<-list(c("(Intercept)",lab.out),c("Coef.","Estimates",
"ci.lower","ci.upper","p-value","Std.Error"))
for(j in 1:al){
ci[j,1]<-j-1
ci[j,2]<-temp$coef[j]
ci[j,3]<-temp$coef[j]-crit*sqrt(hc4[j,j])
ci[j,4]<-temp$coef[j]+crit*sqrt(hc4[j,j])
test<-temp$coef[j]/sqrt(hc4[j,j])
names(test)=NULL
ci[j,5]<-2*(1-pt(abs(test),df))
if(CN)ci[j,5]<-2*(1-pnorm(abs(test),df))
}
ci[,6]=sqrt(diag(hc4))
if(plotit){
if(pnum==1){
plot(x[,-1],y,xlab=xlab,ylab=ylab)
abline(ci[,2])
}
if(pnum==2){
regp2plot(x[,-1],y,regfun=ols,xlab=xlab,ylab=ylab,zlab=zlab)
}}
list(n=nrem,n.keep=n.keep,ci=ci, cov=hc4, test.stat=test,R.squared=rsq)
}

olsci<-olshc4


# con.all.pairs
con.all.pairs<-function(J){
#
# Compute contrast coefficients	 for doing all pairwise comparisons
#
C=(J^2-J)/2
con=matrix(0,nrow=J,ncol=C)
ic=0
for(i in 1:J){
for(j in 1:J){
if(i<j){
ic=ic+1
con[i,ic]=1
con[j,ic]=-1
}}}
con
}


# rmmcppb
rmmcppb<-function(x,y=NULL,alpha=.05,
con=0,est=onestep,plotit=FALSE,dif=TRUE,grp=NA,nboot=NA,BA=FALSE,hoch=FALSE,xlab="Group 1",ylab="Group 2",pr=TRUE,SEED=TRUE,SR=FALSE,...){
#
#   Use a percentile bootstrap method to  compare dependent groups.
#   By default,
#   compute a .95 confidence interval for all linear contrasts
#   specified by con, a J-by-C matrix, where  C is the number of
#   contrasts to be tested, and the columns of con are the
#   contrast coefficients.
#   If con is not specified, all pairwise comparisons are done.
#
#   If est=onestep or mom, method SR (see my book on robust methods)
#   is used to control the probability of at least one Type I error.
#
#   Otherwise, Hochberg is used.
#
#   dif=T indicates that difference scores are to be used
#   dif=F indicates that measure of location associated with
#   marginal distributions are used instead.
#
#   nboot is the bootstrap sample size. If not specified, a value will
#   be chosen depending on the number of contrasts there are.
#
#   x can be an n by J matrix or it can have list mode
#   for two groups, data for second group can be put in y
#   otherwise, assume x is a matrix (n by J) or has list mode.
#
#   A sequentially rejective method is used to control alpha using method SR.
#
#   Argument BA: When using dif=F, BA=T uses a correction term
#   when computing a p-value.
#
if(hoch)SR=FALSE #Assume Hochberg if hoch=TRUE even if SR=TRUE
if(SR){
okay=FALSE
if(identical(est,onestep))okay=TRUE
if(identical(est,mom))okay=TRUE
SR=okay # 'Only use method SR (argument SR=TRUE) when est=onestep or mom
}
if(dif){
if(pr){print("dif=TRUE, so analysis is done on difference scores.")
print(" Each confidence interval has probability coverage 1-alpha.")
print("Also note that a sequentially rejective method is being used")
}
temp<-rmmcppbd(x,y=y,alpha=alpha,con=con,est,plotit=plotit,grp=grp,nboot=nboot, SEED=SEED,
hoch=TRUE,...)
output<-temp$output
con<-temp$con
}
if(!dif){
if(pr){
print("dif=FALSE, so analysis is done on marginal distributions")
if(!BA){
if(identical(est,onestep))print("With M-estimator or MOM, suggest using BA=TRUE and hoch=TRUE")
if(identical(est,mom))print("With M-estimator or MOM, suggest using BA=TRUE and hoch=TRUE")
}}
if(!is.null(y[1]))x<-cbind(x,y)
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
if(is.matrix(con)){
if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
}}
if(is.list(x)){
# put the data in an n by J matrix
mat<-matl(x)
}
if(is.matrix(x) && is.matrix(con)){
if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
mat<-x
}
if(is.matrix(x))mat<-x
if(!is.na(sum(grp)))mat<-mat[,grp]
mat<-elimna(mat) # Remove rows with missing values.
x<-mat
J<-ncol(mat)
xcen<-x
for(j in 1:J)xcen[,j]<-x[,j]-est(x[,j],...)
Jm<-J-1
if(sum(con^2)==0){
d<-(J^2-J)/2
con<-matrix(0,J,d)
id<-0
for (j in 1:Jm){
jp<-j+1
for (k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
d<-ncol(con)
if(is.na(nboot)){
if(d<=4)nboot<-1000
if(d>4)nboot<-5000
}
n<-nrow(mat)
crit.vec<-alpha/c(1:d)
connum<-ncol(con)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
xbars<-apply(mat,2,est,...)
psidat<-NA
for (ic in 1:connum)psidat[ic]<-sum(con[,ic]*xbars)
psihat<-matrix(0,connum,nboot)
psihatcen<-matrix(0,connum,nboot)
bvec<-matrix(NA,ncol=J,nrow=nboot)
bveccen<-matrix(NA,ncol=J,nrow=nboot)
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(ib in 1:nboot){
bvec[ib,]<-apply(x[data[ib,],],2,est,...)
bveccen[ib,]<-apply(xcen[data[ib,],],2,est,...)
}
#
# Now have an nboot by J matrix of bootstrap values.
#
test<-1
bias<-NA
for (ic in 1:connum){
psihat[ic,]<-apply(bvec,1,bptdpsi,con[,ic])
psihatcen[ic,]<-apply(bveccen,1,bptdpsi,con[,ic])
bias[ic]<-sum((psihatcen[ic,]>0))/nboot-.5
ptemp<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
if(BA)test[ic]<-ptemp-.1*bias[ic]
if(!BA)test[ic]<-ptemp
test[ic]<-min(test[ic],1-test[ic])
test[ic]<-max(test[ic],0)  # bias corrected might be less than zero
}
test<-2*test
ncon<-ncol(con)
dvec<-alpha/c(1:ncon) # Assume Hochberg unless specified otherwise
if(SR){
if(alpha==.05){
dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
dvecba<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
dvecba<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
dvecba<-dvec
dvec[2]<-alpha
}}
if(hoch)dvec<-alpha/c(1:ncon)
dvecba<-dvec
if(plotit && ncol(bvec)==2){
z<-c(0,0)
one<-c(1,1)
plot(rbind(bvec,z,one),xlab=xlab,ylab=ylab,type="n")
points(bvec)
totv<-apply(x,2,est,...)
cmat<-var(bvec)
dis<-mahalanobis(bvec,totv,cmat)
temp.dis<-order(dis)
ic<-round((1-alpha)*nboot)
xx<-bvec[temp.dis[1:ic],]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
abline(0,1)
}
temp2<-order(0-test)
ncon<-ncol(con)
zvec<-dvec[1:ncon]
if(BA)zvec<-dvecba[1:ncon]
sigvec<-(test[temp2]>=zvec)
output<-matrix(0,connum,6)
dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.crit","ci.lower","ci.upper"))
tmeans<-apply(mat,2,est,...)
psi<-1
output[temp2,4]<-zvec
for (ic in 1:ncol(con)){
output[ic,2]<-sum(con[,ic]*tmeans)
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(psihat[ic,])
#icl<-round(output[ic,4]*nboot/2)+1 # This adjustment causes confusion; it's not based on Hochberg
icl<-round(alpha*nboot/2)+1
icu<-nboot-(icl-1)
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
}
ids=NA
num.sig=nrow(output)
ior=order(output[,3],decreasing=TRUE)
for(j in 1:nrow(output)){
if(output[ior[j],3]<=output[ior[j],4])break
else num.sig=num.sig-1
}
list(output=output,con=con,num.sig=num.sig)
}


pbos<-function(x,beta=.2){
#
#    Compute the one-step percentage bend measure of location
#
#
temp<-sort(abs(x-median(x)))
omhatx<-temp[floor((1-beta)*length(x))]
psi<-(x-median(x))/omhatx
i1<-length(psi[psi<(-1)])
i2<-length(psi[psi>1])
sx<-ifelse(psi<(-1),0,x)
sx<-ifelse(psi>1,0,sx)
pbos<-(sum(sx)+omhatx*(i2-i1))/(length(x)-i1-i2)
pbos
}

pbos<-function(x,beta=.2){
#
#    Compute the one-step percentage bend measure of location
#
#
temp<-sort(abs(x-median(x)))
omhatx<-temp[floor((1-beta)*length(x))]
psi<-(x-median(x))/omhatx
i1<-length(psi[psi<(-1)])
i2<-length(psi[psi>1])
sx<-ifelse(psi<(-1),0,x)
sx<-ifelse(psi>1,0,sx)
pbos<-(sum(sx)+omhatx*(i2-i1))/(length(x)-i1-i2)
pbos
}


# winval
winval<-function(x,tr=.2){
#
#  Winsorize the data in the vector x.
#  tr is the amount of Winsorization which defaults to .2.
#
#  This function is used by several other functions that come with this book.
#
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
winval<-ifelse(x<=xbot,xbot,x)
winval<-ifelse(winval>=xtop,xtop,winval)
winval
}

# pbvar
pbvar<-function(x,beta=.2){
#
#  Compute percentage bend midvariance
#
#   WARNING: Multiple values for beta seem to make sense
#
temp<-idealf(x)
omega<-temp$n*temp$span
x<-temp$x-median(temp$x)
n<-length(x)
m<-(1-beta)*n
m<-floor(m+.5)
w<-c(1:n)
w<-w/n
pbvar<-(n*wcov(x,w,beta)/(2*(n-m-1)))
pbvar<-pbvar[1,1]
pbvar
}
