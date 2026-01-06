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


# outproMC
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




# outproMC.sub
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

# outproMC.sub2
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

# outproad
outproad<-function(m,center=NA,plotit=TRUE,op=TRUE,MM=TRUE,cop=3,
xlab="VAR 1",ylab="VAR 2",rate=.05,iter=100,ip=6,pr=TRUE,SEED=TRUE,STAND=TRUE){
#
# Adjusts the critical value, gval used by outpro,
# so that the outside rate per observation, under normality
# is approximatley equal to the value given by the argument
# rate, which defaults to .05.
# That is, expected proportion of points declared outliers under normality
# is intended to be rate=.05
#
# When dealing with p-variate data, p>9, this adjustment can be crucial
#
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

# outproadMC
outproadMC<-function(m,center=NA,plotit=TRUE,op=TRUE,MM=TRUE,cop=3,
xlab="VAR 1",ylab="VAR 2",rate=.05,iter=100,ip=6,pr=TRUE,SEED=TRUE){
#
# Adjusts the critical value, gval used by outpro,
# so that the outside rate per observation, under normality
# is approximatley equal to the value given by the argument
# rate, which defaults to .05.
# That is, expected proportion of points declared outliers under normality
# is intended to be rate=.05
#
# When dealing with p-variate data, p>9, this adjustment can be crucial
#
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




# outpro.depth
outpro.depth<-function(x,ndir=1000,MM=FALSE,SEED=TRUE,plotit=FALSE,xlab='X',ylab='Y'){
#
#  Use projection distances to detect outliers.
# This function can handle large sample sizes
#
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


# outbox
outbox<-function(x,mbox=FALSE,gval=NA,plotit=FALSE,STAND=FALSE){
#
# This function detects outliers using the
# boxplot rule, but unlike the R function boxplot,
# the ideal fourths are used to estimate the quartiles.
#
# Setting mbox=TRUE results in using the modification
# of the boxplot rule suggested by Carling (2000).
#
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


# out3d
out3d<-function(x,outfun=outpro,xlab="Var 1",ylab="Var 2",zlab="Var 3",
reg.plane=FALSE,regfun=tsreg,COLOR=FALSE,tick.marks=TRUE,...){
#
# Create a 3D plot of points and indicate outliers with red dots.
#
#  Assumes that the package scatterplot3d has been installed.
#  If not, use the command install.packages("scatterplot3d")
#  assuming you are connected to the web.
#
# To add a regression plane, set
#  reg.plane=T.The regression method used is specified with the argument
#  regfun.
# First two columns are taken to be predictors and third column is the outcome
#
#  Package scatterplot3d is required. To install it, use the command
#  install.packages("scatterplot3d")
#  while connected to the web
#
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


# outbag
outbag<-function(x,plotit=FALSE){
#
#  Search for outliers using bagplot
#  bivariate data only.
#
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



# outmah
outmah<-function(x,qval=pnorm(3),plotit=TRUE,xlab="VAR 1",ylab="VAR 2"){
#
#  detect outliers using Mahalanobis Distance
#   For demonstration purposes only. Suggest
#   using a method that avoids masking.
#
#  In univariate case, default strategy is to use 3 standard deviation rule
#
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


# outmve
outmve<-function(x,mve.flag=TRUE,plotit=TRUE,SEED=TRUE,outsym='*'){
#
#  Search for outliers using the minimum volume ellipsoid method.
#
#  x is an n by p matrix
#
#  The function returns the number of the rows in x that are identified
#  as outliers. (The row numbers are stored in outliers.)
#  It also returns the distance of the points identified as outliers
#  in the variable dis.
#
#  If mve.flag=T, use the mve estimator, otherwise use the mcd
#
#  If plotit=TRUE, plot points and circle outliers.
#
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


# outmgv
outmgv<-function(x,y=NULL,plotit=TRUE,outfun=outbox,se=TRUE,op=1,ndir=1000,
cov.fun=rmba,xlab="X",ylab="Y",SEED=TRUE,STAND=FALSE,...){
#
# Check for outliers using mgv method
#
# NOTE: if columns of the input matrix are reordered, this can
# have an effect on the results due to rounding error when calling
# the R function eigen.
#
#  (Argument STAND is included simply to avoid programming issues when outmgv is called by other functions.)
#
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


# outmgvad
outmgvad<-function(m,center=NA,plotit=TRUE,op=1,
xlab="VAR 1",ylab="VAR 2",rate=.05,iter=100,ip=6,pr=TRUE){
#
# Adjusts the critical value, gval used by outmgv,
# so that the outside rate per observation, under normality
# is approximately equal to the value given by the argument
# rate, which defaults to .05.
# That is, expected proportion of points declared outliers under normality
# is intended to be rate=.05
#
# When dealing with p-variate data, p>9, this adjustment can be crucial
#
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



# outmgvf
outmgvf<-function(x,y=NA,plotit=TRUE,outfun=outbox,se=TRUE,ndir=1000,SEED=TRUE,...){
#
# Check for outliers using inward mgv method
# This method is faster than outmgv.
#
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


# outmgv.v2
outmgv.v2<-function(x,outfun=outbox,se=TRUE,op=1,
gval=sqrt(qchisq(.975,ncol(x))),
cov.fun=rmba,xlab="X",ylab="Y",SEED=TRUE,...){
#
# Check for outliers using mgv method
#
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


# outms
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

# outogk
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


# outcov
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


# outtbs
outtbs<-function(x,SEED=FALSE,plotit=TRUE,xlab="X",ylab="Y",...){
#
# Use the tbs estimator to
# determine which points are outliers
#
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)
temp<-out(x,cov.fun=tbs,plotit=plotit,SEED=SEED,xlab=xlab,ylab=ylab)
outid<-temp$out.id
keep<-temp$keep
list(out.id=outid,keep=keep,distances=temp$dis)
}


# outDETMCD
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


# outICS
outICS<-function(x,n.id=NULL){
#
#  Detect outliers using the ICS method in Archimbaud et al
#  CSDA
#
# Assumes the R package ICSOutlier has been installed.
#
#  To id outliers, must run this function twice.
#  First time, determine how many outliers were dectected
#  Second time, set n.id = number of outliers
#
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


# outmc
outmc<-function(x,plotit=FALSE){
#
#  Detect outliers using a modification of Carling's method
#  that takes into account skewness
#
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


# out.dummy
out.dummy<-function(x,outfun=outpro,id,plotit=FALSE,...){
#
#  When using dummy coding in regression
#
#  remove col indicated by
#  id
# then check for outliers using
# outfun
x=as.matrix(x)
if(ncol(x)==1)stop(' Should have two or more columns')
X=x[,-id]
a=outfun(X,plotit=FALSE)
a
}


# out.by.groups
out.by.groups<-function(x,grp.col,outfun=outpro,pr=TRUE,plotit=FALSE,...){
#
# divide data into groups, id outliers in each group
# return:
# keep = id rows in x not outliers
#  out.id =rows containing outliers
#
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



# out.methods
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


# outblp.HH
outblp.HH<-function(x,y,regfun=tsreg,omit.col=NULL,plotit=TRUE,xlab='X',ylab='Y'){
#
# indicates which points, if any, are bad leverage points
# using a blend of a homoscedastic and heteroscedastic methods.
#
# This approach helps to avoid issues with Type I errors when testing hpotheses
#
# If for example
# omit.col=c(1,3)
# columns 1 and 3 of x are ignored when checking for bad leverage points.
#   These columns might be, for example, dummy variables.
#
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




# depth2
depth2<-function(x,pts=NA,plotit=TRUE,xlab="VAR 1",ylab="VAR 2"){
#
#   Compute exact depths for bivariate data
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


# depthcom
depthcom<-function(x1,y1,x2,y2,est=tmean,fr=1){
temp1=depthcomsub(x1,y1,x2,y2,est=est,fr=fr)
temp2=depthcomsub(x2,y2,x1,y1,est=est,fr=fr)
dep=max(c(abs(temp1$dep1-temp1$dep2),abs(temp2$dep1-temp2$dep2)))
dep
}

# depthcomsub
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


# depthg2
depthg2<-function(x,y,alpha=.05,nboot=500,MD=FALSE,plotit=TRUE,op=FALSE,fast=FALSE,SEED=TRUE,
xlab="VAR 1",ylab="VAR 2"){
#
#   Compare two independent groups based on p measures
#   for each group.
#
#   The method is based on Tukey's depth if MD=F;
#   otherwise the Mahalanobis depth is used.
#   If p>2, then Mahalanobis depth is used automatically
#
#   The method is designed to be sensitive to differences in scale
#
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


# depths1
depths1<-function(m,j){
if(m < j)depths1<-0
else{
if(j==1)depths1<-m
if(j==2)depths1<-(m*(m-1))/2
if(j==3)depths1<-(m*(m-1)*(m-2))/6
}
depths1
}


# fdepthv2
fdepthv2<-function(m,pts=NA,plotit=TRUE){
#
# Determine depth of points in pts relative to
# points in m
#
# Draw a line between each pair of distinct points
# and determine depth of the projected points.
# The final depth of a point is its minimum depth
# among all projections.
#
# This function is slower than fdepth and requires
# space for a nc by nc matrix, nc=(n^2-n)/2.
# But it allows
# data to have a singular covariance matrix
# and it provides a more accurate approximation of
# halfspace depth.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data and pts=NA
#
#  When plotting,
#  center is marked with a cross, +.
#
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


# indepth
indepth<-function(m){
#
# Compute the inward depth of all points in m
# based on the generalized variance.
#
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


# pdepth
pdepth<-function(m,pts=m,MM=FALSE,cop=3,dop=1,center=NA, SEED=TRUE){
#
# projection depth
#
#  SEED, included for convenience when this function is used with certain classification techniques.
#
v=pdis(m,pts=pts,MM=MM,cop=cop,dop=dop,center=center)
v=1/(1+v)
v
}

# prodepth
prodepth<-function(x,pts=x,ndir=1000,SEED=TRUE){
#
#  Determine an approximation of the projection depth of
#  pts in
#  x
#  using the R package library(DepthProc)
#
#  ndir indicates how many randomly chosen projections will be used
#
#  Advantage over zoudepth, much faster execution time.
#  Should be noted, however, that using the function twice on the same
#  data generally results in different values for the depths.
#  Setting
#  SEED=TRUE
#  avoids this.
#
#
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

# rdepth.orig
rdepth.orig<-function(d, x, y, sortx = TRUE)
{
##########################################################################
# This function computes the regression depth of a line with coordinates d
# relative to the bivariate data set (x,y).
# The first component of the vector d indicates the intercept of the line,
# the second component is the slope.
#
# Input : d          : vector with two components
#         x,y        : vectors of equal length (data set)
#         sortx      : logical, to set to F if the data set (x,y) is
#                      already sorted by its x-coordinates
#
# Reference:
#           Rousseeuw, P.J. and Hubert, M. (1996),
#           Regression Depth, Technical report, University of Antwerp
#           submitted for publication.
##########################################################################
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


# resdepth
resdepth<-function(x,res)
{
##########################################################################
# This function computes the regression depth of a regression line based
# on its residuals. The fit could be, for example, a nonparametric
# regression or smooth.
#
# The algorithm is based on a simple modification of
#
#           Rousseeuw, P.J. and Hubert, M. (1996),
#           Regression Depth, Technical report, University of Antwerp
#
##########################################################################
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

# resdepth.sub
resdepth.sub<-function(x,res)
{
##########################################################################
# This function computes the regression depth of a regression line based
# on its residuals. The fit could be, for example, a nonparmatric
# regression or smooth.
#
# The algorithm is based on a simple modification of
#
#           Rousseeuw, P.J. and Hubert, M. (1996),
#           Regression Depth, Technical report, University of Antwerp
#
##########################################################################
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


# zdepth
zdepth<-function(m,pts=m,zloc=median,zscale=mad){
#
# Compute depth of points as in Zuo, Annals, 2003
#
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


# zdepth.sub
zdepth.sub<-function(x,theta,zloc=median,zscale=mad,pts=NA){
theta<-theta/sqrt(sum(theta^2))
temp<-t(t(x)*theta)
ppro<-sum(t(t(pts)*theta))
yhat<-apply(temp,1,sum)
val<-0-abs(ppro-zloc(yhat))/zscale(yhat)
val
}

zdist=zdepth


# zoudepth
zoudepth<-function(x,pts=x, zloc = median, zscale = mad, SEED=TRUE){
#
#  Determine projection depth using the R function zdepth
#  The Nelder--Mead method for finding the maximum of a function is used
#
#  SEED, included for convenience when this function is used with certain classification techniques.
#
res=1/(1+zdepth(x,pts,zloc,zscale))
res
}

# unidepth
unidepth<-function(x,pts=NA){
#
# Determine depth of points in the vector x
#
if(!is.vector(x))stop("x should be a vector")
if(is.na(pts[1]))pts<-x
pup<-apply(outer(pts,x,FUN="<="),1,sum)/length(x)
pdown<-apply(outer(pts,x,FUN="<"),1,sum)/length(x)
pdown<-1-pdown
m<-matrix(c(pup,pdown),nrow=2,byrow=TRUE)
dep<-apply(m,2,min)
dep
}


# discdepth
discdepth<-function(train=NULL,test=NULL,g,x1=NULL,x2=NULL,depthfun=prodepth,...){
#
# x1 and x2  contain the data for the two groups
# Goal, classify the values in test using depths associated with the training data
#
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


# mregdepth
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



# smean.depth
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


# pbadepth
pbadepth<-function(x,est=onestep,con=0,alpha=.05,nboot=2000,grp=NA,op=3,allp=TRUE,
MM=FALSE,MC=FALSE,cop=3,SEED=TRUE,na.rm=FALSE,...){
#
#   Test the hypothesis that C linear contrasts all have a value of zero.
#   By default, an M-estimator is used
#
#   Independent groups are assumed.
#
#   The data are assumed to be stored in x in list mode or in a matrix.
#   If stored in list mode,
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J, say.
#   If stored in a matrix, columns correspond to groups.
#
#   By default, all pairwise differences are used, but contrasts
#   can be specified with the argument con.
#   The columns of con indicate the contrast coefficients.
#   Con should have J rows, J=number of groups.
#   For example, con[,1]=c(1,1,-1,-1,0,0) and con[,2]=c(,1,-1,0,0,1,-1)
#   will test two contrasts: (1) the sum of the first
#   two measures of location is
#   equal to the sum of the second two, and (2) the difference between
#   the first two is equal to the difference between the
#   measures of location for groups 5 and 6.
#
#   The default number of bootstrap samples is nboot=2000
#
#   op controls how depth is measured
#   op=1, Mahalanobis
#   op=2, Mahalanobis based on MCD covariance matrix
#   op=3, Projection distance
#
#   MC=TRUE, use a multicore processor when op=3
#
#   for arguments MM and cop, see pdis.
#
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


# Qdepthcom
Qdepthcom<-function(x1,y1,x2,y2,qval){
temp1=Qdepthcomsub(x1,y1,x2,y2,qval)
temp2=Qdepthcomsub(x2,y2,x1,y1,qval)
dep=max(c(abs(temp1$dep1-temp1$dep2),abs(temp2$dep1-temp2$dep2)))
dep
}

# Qdepthcomsub
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



# comdepthsvm
comdepthsvm<-function(x1,x2,alpha=.05,depthfun=prodepth,
plotit=FALSE,kernel='radial',MISS=FALSE,TABLE=FALSE,...){
#
# compare two independent multivariate distributions based
# on a basic support vector machines method.
#
#  MISS=TRUE: returns the vectors misclassified.
#
# Leave-one-out cross validation is used, but see
# Shao (1993). Linear Model Selection by Cross-Validation, JASA, 88, 486--494
#
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



# aov2depth
aov2depth<-function(x1,x2,est=tmean,nboot=500,SEED=TRUE,nmin=12,CR=FALSE,
xlab=' DIF 1',ylab='DIF 2',zlab='DIF 3',alpha=.05,...){
#
# 2 by K ANOVA independent group (K levels not necessarily independent and
#                                 not completely dependent
#
#   Main effect Factor A only
#
# Strategy: Use depth of zero based on estimated
# differences for each column  of the K levels of Factor B
# That is, testing no main effects for Factor A in
# a manner that takes into account the pattern of the
# measures of location rather then simply averaging
# across columns.
#
#  x1 can be a matrix with K columns corrspoding to groups, ditto for x2
#  Or x1 and x2 can have list mode.
#   Assuming x1 and x2 contain data for indepedendent groups.
#
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


# bagdepth
bagdepth<-function(x,pts=NULL,SEED=TRUE){
#
#  Compute depth of points based in their bagdistance.
#
# requires R package mrfDepth
#
#  SEED, included for convenience when this function is used with certain classification techniques.
#
library(mrfDepth)
d=bagdistance(x,pts)$bagdistance
d=1/(d+1)
d
}


# bwdepth
bwdepth<-function(x,y,fun=prodepth,plotit=FALSE,xlab='V1',ylab='V2'){
#
# For two independent groups, let X and Y denote multivariate random variables
# This function estimates the extent the distributions overlap using the notion
# of projection distances. In effect, a nonparametric measure of effect size is
# estimated
# For identical distribution, effect size is .5. The more separated the distributions, the
# closer is the effect size to zero. Complete separation means the effect size is equal to zero
#
#
#
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


# bwdepthMC.ci
bwdepthMC.ci<-function(x,y,fun=prodepth,nboot=100,alpha=.05,MC=TRUE,
SEED=TRUE,plotit=FALSE,xlab='V1',ylab='V2'){
#
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


# bwdepth.perm
bwdepth.perm<-function(x,y,reps=500,
fun=prodepth,alpha=.05,SEED=TRUE){
#
# Permutation test of F=G, two independent multivariate distributions
#
#
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


# bwdepth.sub
bwdepth.sub<-function(id,x,y,n1,n2,fun){
n=n1+n2
np1=n1+1
e=bwdepth(x[id[1:n1],],y[id[np1:n],],fun=fun)$e
e
}


# Depth.class.bag
Depth.class.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,DIST=FALSE,nboot=100,SEED=TRUE,...){
#
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



# dis.depth.bag
dis.depth.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,nboot=100,SEED=TRUE,...){
#
#
#  g=class id
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2.
#  The function will then create appropriate labels and store them in g.
#
#  Uses data depths.
# KNNdist uses data depths, for the n1!=n2 it can be a bit biased, meaning that
# when there is no association, the probability of a correct classification will be less than .5
#
#
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



# pro.class.bag
pro.class.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,kernel='radial',nboot=20,
PR=FALSE,SEED=TRUE,...){
#
#
# pro.class: for n1!=n2 it can be a  biased, meaning that
# when there is no association, the probability of a correct classification will be less than .5
#
#  This function deals with this via bootstrap bagging
#  g=class labels
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2, which are assumed to be matrices or a data frame.
#  The function will then create appropriate labels and store them in g.
#
#
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




# pro.classPD.bag
pro.classPD.bag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,rule=1,nboot=100,SEED=TRUE,...){
#
#
#  A bagged version of pro.classPD
#
#  g=class id
#  if there are two classes and the training data are stored in  separate variables, can enter
#  the data for each class via the arguments
#  x1 and x2.
#  The function will then create appropriate labels and store them in g.
#
#
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
KNNbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,nboot=100,SEED=TRUE,...){
#
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

# LSMbag
LSMbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,
sm=TRUE,rule=.5,nboot=100,SEED=TRUE,...){
#
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




# NNbag
NNbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,nboot=100,SEED=TRUE,...){
#
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


# RFbag
RFbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,kernel='radial',nboot=100,SEED=TRUE,...){
#
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


# SVMbag
SVMbag<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,kernel='radial',nboot=100,SEED=TRUE,...){
#
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


