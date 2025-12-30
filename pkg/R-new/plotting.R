# ============================================================================
# plotting.R - Visualization and Plotting Functions
# ============================================================================
#
# Comprehensive collection of plotting and visualization functions for
# robust statistical methods.
#
# Categories:
# - Regression plotting (rplot, regplot, etc.)
# - Lowess/loess smoothing plots
# - GAM plotting
# - ANCOVA plotting
# - Group comparison plots
# - Error bar and box plots
# - Functional data plots
# - Interaction plots
# - Logistic/longitudinal regression plots
# - Depth-based plots (bagplot)
# - Distribution/density plots
# - 3D plotting
# - Specialized/utility plots
#
# Total functions: 97 (98 - 1 duplicate linplot excluded)
# ============================================================================


# ============================================================================
# REGRESSION PLOTTING
# ============================================================================

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

reg2plot<-function(x1,y1,x2,y2,regfun=tsreg,xlab="X",ylab="Y",xout=FALSE,outfun=outpro,pch1='.',pch2='+',...){
#
#  For convenience
#  plot two regression lines corresponding to two groups.
#
xy=elimna(cbind(x1,y1))
x1=xy[,1]
y1=xy[,2]
xy=elimna(cbind(x2,y2))
x2=xy[,1]
y2=xy[,2]
if(xout){
if(identical(outfun,outblp))flag=outblp(x1,y1,plotit=FALSE)$keep
else
flag<-outfun(x1,plotit=FALSE,...)$keep
x1=x1[flag]
y1=y1[flag]
if(identical(outfun,outblp))flag=outblp(x2,y2,plotit=FALSE)$keep
else
flag<-outfun(x2,plotit=FALSE,...)$keep
x2=x2[flag]
y2=y2[flag]
}
plot(c(x1,x2),c(y1,y2),type="n",xlab=xlab,ylab=ylab)
points(x1,y1,pch=pch1)
points(x2,y2,pch=pch2)
abline(regfun(x1,y1,...)$coef)
abline(regfun(x2,y2,...)$coef,lty=2)
}

reg2g.p2plot<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,xlab="Var 1",ylab="Var 2",zlab="Var 3",regfun=tsreg,COLOR=TRUE,STAND=TRUE,
tick.marks=TRUE,type="p",pr=TRUE,...){
#
# Create a 3D plot of points and plot regression surface for two groups.
#
#  Assumes that the package scatterplot3d has been installed.
#  If not, use the command install.packages("scatterplot3d")
#  assuming you are connected to the web.
#
# The regression method used is specified with the argument
#  regfun.
#
#  type="p", points will be plotted. Use type="n" to get only regression planes plotted
#
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=2)stop("Argument x1 must be stored in a matrix with 2 columns.")
if(ncol(x2)!=2)stop("Argument x2 must be stored in a matrix with 2 columns.")
xy1<-elimna(cbind(x1,y1))
xy2<-elimna(cbind(x2,y2))
if(xout){
if(!STAND)flag1=outfun(xy1[,1:2],plotit=FALSE,...)$keep
if(STAND)flag1=outpro(xy1[,1:2],plotit=FALSE,STAND=TRUE,...)$keep
if(!STAND)flag2=outfun(xy2[,1:2],plotit=FALSE,...)$keep
if(STAND)flag2=outpro(xy2[,1:2],plotit=FALSE,STAND=TRUE,...)$keep
xy1=xy1[flag1,]
xy2=xy2[flag2,]
}
x1=xy1[,1:2]
x2=xy2[,1:2]
y1=xy1[,3]
y2=xy2[,3]
library(scatterplot3d)
temp<-scatterplot3d(rbind(xy1,xy2),xlab=xlab,ylab=ylab,zlab=zlab,tick.marks=tick.marks,type=type)
vals1<-regfun(x1,y1,...)$coef
vals2<-regfun(x2,y2,...)$coef
if(COLOR){
if(pr)print("First group is blue")
temp$plane(vals1,col="blue")
temp$plane(vals2,col="red")
}
if(!COLOR){
temp$plane(vals1)
temp$plane(vals2)
}
list(coef.group.1=vals1,coef.group.2=vals2)
}

regp2plot<-function(x,y,xout=FALSE,outfun=out,xlab="Var 1",ylab="Var 2",zlab="Var 3",regfun=tsreg,COLOR=FALSE,tick.marks=TRUE,...){
#
# Create a 3D plot of points and plot regression surface.
# based on the regression estimator indicated by
#  regfun
#
#  Assumes that the package scatterplot3d has been installed.
#  If not, use the command install.packages("scatterplot3d")
#  assuming you are connected to the web.
#
# The regression method used is specified with the argument
#  regfun.
#
#  Package scatterplot3d is required. To install it, use the command
#  install.packages("scatterplot3d")
#  while connected to the web
#
x=as.matrix(x)
if(ncol(x)!=2)stop("Argument x must be stored in a matrix with 2 columns.")
xy<-elimna(cbind(x,y))
if(xout){
flag=outfun(xy[,1:2])$keep
xy=xy[flag,]
}
x=xy[,1:2]
y=xy[,3]
library(scatterplot3d)
temp<-scatterplot3d(xy,xlab=xlab,ylab=ylab,zlab=zlab,tick.marks=tick.marks)
vals<-regfun(x,y,...)$coef
if(COLOR)temp$plane(vals,col="blue")
if(!COLOR)temp$plane(vals)
}

regplot<-function(x,y,regfun=tsreg,xlab='X',ylab='Y',zlab='Z',
xout=FALSE,outfun=out,theta=50, phi=25,ticktype='simple',...){
x=as.matrix(x)
xy=elimna(cbind(x,y))
if(ncol(x)>2)stop('One or two predictors only is allowed,')
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
if(xout){
xy=cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag=outfun(x)$keep
x=xy[flag,1:p]
y=xy[flag,p1]
}
if(p==1){
plot(x,y,xlab=xlab,ylab=ylab)
abline(regfun(x,y,...)$coef)
}
if(p==2){
pyhat=regYhat(x,y,regfun=regfun,...)
temp=rplot(x,pyhat,scat=FALSE,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,pr=FALSE)
}
}

riplot<-function(x,y,adfun=adrun,plotfun=lplot,eout=FALSE,xout=TRUE,scale=FALSE){
#
# Plot used to investigate regression interaction
# (the extent a generalized additive model does not fit data).
# Compute additive fit, plot residuals
# versus x, an n by 2 matrix.
#
if(!is.matrix(x))stop(" x must be a matrix")
if(ncol(x)!=2)stop(" x must have two columns only")
yhat<-adfun(x,y,pyhat=TRUE,eout=eout,xout=xout,plotit=FALSE)
plotfun(x,y-yhat,eout=eout,xout=xout,scale=scale)
}

rplot.res<-function(x,y,pv=1,est = tmean, scat = TRUE, fr = NA, plotit = TRUE,
pyhat = FALSE, efr = 0.5, theta = 50, phi = 25, scale = TRUE,
expand = 0.5, SEED = TRUE, varfun = pbvar, outfun = outpro,STAND=TRUE,
nmin = 0, xout = FALSE, out = FALSE, eout = FALSE, xlab='X',
ylab ='Y',zscale=FALSE,zlab=' ', pr=TRUE,duplicate='error',
ticktype='simple',LP=TRUE,...){
#
# Apply rplot excluding the independent variable indicated by the argument
# pv.
# So pv=1 means will exclude the first predictor.
# Fit a smooth using the remaing variables, compute the residuals, then plot
# the smooth using the residuals as the dependent variable and
# the variables indicated by pv as the independent variables.
#
xy=na.omit(cbind(x,y))
p=ncol(x)
p1=p+1
if(xout){
flag=outfun(xy[,1:p],plotit=FALSE,STAND=STAND,...)$keep
xy=xy[flag,]
}
x=xy[,1:p]
y=xy[,p1]
res=y-rplot(x[,1:2],y,est=est,scat=scat,varfun=varfun,expand=expand,nmin=nmin,
pyhat=TRUE,plotit=FALSE,fr=fr,xout=FALSE)$yhat
outp=rplot(x[,pv],res,fr=fr,xout=FALSE,efr=efr,theta=theta,phi=phi,
scale=scale,SEED=SEED,xlab=xlab,ylab=ylab,zlab=zlab,pr=FALSE,
ticktype=ticktype,LP=LP,...)
outp
}

rplotCI<-function(x,y,tr=.2,fr=.5,p.crit=NA,plotit=TRUE,scat=TRUE,
SEED=TRUE,pyhat=FALSE,npts=25,xout=FALSE,
xlab='X',ylab='Y',low.span=2/3,nmin=12,pr=TRUE,
outfun=outpro,LPCI=FALSE,MID=TRUE,alpha=.05,pch='.',...){
#
# Confidence interval for running  interval smoother based on a trimmed mean.
# Unlike rplot, includes an approximate  confidence band having simultaneous probability
# coverage equal to 1-alpha.   More precisely, the simultaneous probability
# is for K=npts points
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing
#
# To specify the points where confidence intervals are computed,
# use rplotCIsmm
#
if(pr){
if(!LPCI)print('To get smoother plot, set LPCI=TRUE')
}
m<-cbind(x,y)
if(ncol(m)>2)stop('Only one independent variable can be used')
m<-elimna(m)
x=m[,1]
y=m[,2]
if(xout){
xy=cbind(x,y)
flag=outfun(x,plotit=FALSE)$keep
x=xy[flag,1]
y=xy[flag,2]
}
n.used=NA
n=length(y)
if(n<50)stop('Need at least n=50')
nv=c(50,  60,  70,  80, 100,
150, 200, 300, 400, 500, 600, 800, 1000)
if(npts==25) pv=c(0.004846408,
0.004553274,
0.004236101,
0.004099674,
 0.00353898,   #n=100
 0.003366984,
0.003038767,
 0.003061386,
  0.002793521,
  0.002479689,
 0.002606313,
 0.0026630370,
 0.002836043)
if(npts==10) pv=c(
0.007612451,
0.008383655,
0.006992874,
 0.0068073,
0.005733889,
0.005767139,
0.006130155,
0.005447419,
0.005555951,
0.005228471,
0.005642503,
0.005402162,
0.005569217)
FLAG=FALSE
if(npts==25 || npts==10)FLAG=TRUE
if(alpha!=.05 || !FLAG){
if(is.na(p.crit)){
print('p.crit must be estimated, execution time might be high')
print('Or use the R function rplotCIsmm')
}
p.crit=rplotCITAP.pv(n,tr=tr,fr=fr,alpha=alpha,nmin=nmin,npts=npts,nreps=nreps)
}
rem.n=n
if(n>1000)n=1000
if(is.na(p.crit))p.crit=lplot.pred(1/nv,pv,1/n)$yhat
n=rem.n
xord=order(x)
x=x[xord]
y=y[xord]
infit=rplot(x,y,tr=tr,xout=FALSE,plotit=plotit,LP=LPCI,fr=fr,pr=FALSE,pyhat=TRUE,xlab=xlab,
ylab=ylab)
rmd=infit$pyhat
m<-cbind(x,y)
if(ncol(m)>2)stop('One covariate only is allowed with this function')
m<-elimna(m)
nv=nrow(m)
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
x=m[,1]
y=m[,2]
n.keep=length(y)
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr,nmin=nmin)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
y.hat=NA
civ=matrix(NA,nrow=npts,ncol=2)
for(i in 1:length(pts)){
xx=y[near(x,pt=pts[i],fr)]
doit=trimci(xx,tr=tr,alpha=p.crit,pr=FALSE)
civ[i,]=doit$ci
y.hat[i]=doit$estimate
n.used[i]=doit$n
}
up=civ[,2]
low=civ[,1]
if(plotit){
if(LPCI){
up=lplot(pts,up,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
y.hat=lplot(pts,y.hat,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
low=lplot(pts,low,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
}
lines(pts,up,lty=2)
lines(pts,low,lty=2)
}
if(pyhat){output<-cbind(pts,y.hat,low,up,n.used)
dimnames(output)=list(NULL,c('pts','y.hat','ci.low','ci.up','n.used'))
}
if(!pyhat)output<-'Done'
list(output=output,str=infit$Strength.Assoc,n=nv,n.keep=n.keep)
}

rplotCIS<-function(x,y,tr=.2,fr=.8,plotit=TRUE,scat=TRUE,pyhat=FALSE,SEED=TRUE,dfmin=8,
eout=FALSE,xout=FALSE,xlab='x',ylab='y',outfun=out,LP=TRUE,alpha=.05,pch='.',...){
#
# A simple method for computing a confidence band based on
# running  interval smoother and a trimmed mean.
#
#  rplotCI adjusts the band so that FWE=1-alpha
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing
plotit<-as.logical(plotit)
scat<-as.logical(scat)
str=rplot(x,y,tr=tr,xout=xout,plotit=FALSE,LP=LP,fr=fr,pr=FALSE)$Strength.Assoc
m<-cbind(x,y)
if(ncol(m)>2)stop('Only one independent variable can be used')
m<-elimna(m)
nv=nrow(m)
if(eout && xout)stop('Not allowed to have eout=xout=T')
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
x=m[,1]
y=m[,2]
n.keep=length(y)
rmd<-c(1:length(x))
for(i in 1:length(x))rmd[i]<-mean(y[near(x,x[i],fr=fr)],tr=tr)
sedf=runse(x,y,fr=fr,tr=tr,pts=x,SEED=SEED)
df=sedf$df
flag=df>dfmin
se=sedf$se
low=rmd[flag]-qt(1-alpha/2,df[flag])*se[flag]
up=rmd[flag]+qt(1-alpha/2,df[flag])*se[flag]
rmd=rmd[flag]
x=x[flag]
y=y[flag]
if(plotit){
ord=order(x)
x=x[ord]
rmd=rmd[ord]
up=up[ord]
low=low[ord]
if(LP){
rmd=lplot(x,rmd,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
up=lplot(x,up,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
low=lplot(x,low,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
if(scat){
plot(c(x,x),c(y,rmd),xlab=xlab,ylab=ylab,type='n')
lines(x,up,lty=2)
lines(x,low,lty=2)
points(x,y,pch=pch)
}
if(!scat)plot(c(x,x),c(y,rmd),type='n',ylab=ylab,xlab=xlab)
points(x,rmd,type='n')
sx<-sort(x)
xorder<-order(x)
sysm<-rmd[xorder]
lines(sx,sysm)
lines(x,up,lty=2)
lines(x,low,lty=2)
}
if(pyhat){output<-cbind(x,rmd,low,up)
dimnames(output)=list(NULL,c('x','y.hat','ci.low','ci.up'))
}
if(!pyhat)output<-'Done'
list(output=output,str=str,n=nv,n.keep=n.keep)
}

rplotpbCI<-function(x,y,est=onestep,fr=1,plotit=TRUE,scat=TRUE,pyhat=FALSE,
xout=FALSE,xlab='x',ylab='y',outfun=out,LP=TRUE,alpha=.05,
nboot=500,SEED=TRUE,...){
#
# running  interval smoother based on any measure of location
# Unlike rplotCI, uses a percentile bootstrap
# method to get a confidence band
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing
plotit<-as.logical(plotit)
scat<-as.logical(scat)
m<-cbind(x,y)
if(ncol(m)>2)stop('Only one independent variable can be used')
m<-elimna(m)
x=m[,1]
y=m[,2]
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
x=m[,1]
y=m[,2]
low=rep(NA,length(y))
up=rep(NA,length(y))
rmd<-NA
for(i in 1:length(x)){
sel=y[near(x,x[i],fr)]
temp=onesampb(sel,est=est,nboot=nboot,alpha=alpha,SEED=SEED,...)
low[i]=temp$ci[1]
up[i]=temp$ci[2]
rmd[i]=temp$estimate
}
all=elimna(cbind(x,low,up,y,rmd))
x=all[,1]
low=all[,2]
up=all[,3]
y=all[,4]
rmd=all[,5]
if(plotit){
ord=order(x)
x=x[ord]
y=y[ord]
rmd=rmd[ord]
up=up[ord]
low=low[ord]
if(LP){
rmd=lplot(x,rmd,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
up=lplot(x,up,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
low=lplot(x,low,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
if(scat){
plot(c(x,x),c(y,rmd),xlab=xlab,ylab=ylab,type='n')
points(x,y)
lines(x,up,lty=2)
lines(x,low,lty=2)
}
if(!scat)plot(c(x,x),c(y,rmd),type='n',ylab=ylab,xlab=xlab)
points(x,rmd,type='n')
sx<-sort(x)
xorder<-order(x)
sysm<-rmd[xorder]
lines(sx,sysm)
lines(x,up,lty=2)
lines(x,low,lty=2)
}
if(pyhat)output<-rmd
if(!pyhat)output<-'Done'
list(output=output)
}

rplotCIM<-function(x,y,est=hd,fr=.5,p.crit=NA,plotit=TRUE,scat=TRUE,
pyhat=FALSE, pts=NA,npts=25,xout=FALSE,
xlab='X',ylab='Y',low.span=2/3,nmin=16,
outfun=out,LP=TRUE,LPCI=FALSE,MID=TRUE,alpha=.05,pch='.',...){
#
# Confidence interval for running  interval smoother based on a median
# Unlike rplot, includes a confidence band having simultaneous probability
# coverage equal to 1-alpha.
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing. If the association is relatively strong, might want to use fr=.2
#
chk=FALSE
if(identical(est,hd))chk=TRUE
if(!chk)stop('Current version, argument est must be hd')
n=length(y)
if(n<50)stop('Need at least n=50')
xord=order(x)
x=x[xord]
y=y[xord]
infit=rplot(x,y,est=est,xout=xout,plotit=plotit,LP=LP,fr=fr,pr=FALSE,pyhat=TRUE,xlab=xlab,ylab=ylab)
rmd=infit$yhat
m<-cbind(x,y)
if(ncol(m)>2)stop('One covariate only is allowed with this function')
m<-elimna(m)
nv=nrow(m)
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
x=m[,1]
y=m[,2]
n.keep=length(y)
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr,nmin=nmin)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
flag=duplicated(pts)
npts=length(pts)
civ=matrix(NA,nrow=npts,ncol=2)
for(i in 1:length(pts)){
xx=y[near(x,pt=pts[i],fr)]
civ[i,]=sint(xx,alpha=alpha/npts)
}
up=civ[!flag,2]
low=civ[!flag,1]
if(plotit){
if(LPCI){
up=lplot(pts,up,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
low=lplot(pts,low,plotit=FALSE,pyhat=TRUE,pr=FALSE,low.span=low.span)$yhat
}
pts=pts[!flag]
lines(pts,up,lty=2)
lines(pts,low,lty=2)
}
if(pyhat){output<-cbind(pts,rmd[!flag],low,up)
dimnames(output)=list(NULL,c('pts','y.hat','ci.low','ci.up'))
}
if(!pyhat)output<-'Done'
list(output=output,str=infit$Strength.Assoc,n=nv,n.keep=n.keep)
}

rplotCIsmm<-function(x,y,tr=.2,fr=.5,plotit=TRUE,scat=TRUE,pyhat=FALSE,SEED=TRUE,
dfmin=2,pts=NULL,npts=25,nmin=12,
eout=FALSE,xout=FALSE,xlab='x',ylab='y',outfun=out,LP=TRUE,MID=TRUE,alpha=.05,pch='.',...){
#
# Confidence interval for running  interval smoother based on a trimmed mean.
#
#  rplotCI will provide shorter and more accurate confidence intervals but
#  is limited to 10 or 25 points and alpha=.05.
#  This functions returns confidence intervals that are generally a bit wider
#  but it has low execution time if alpha differs from 0.5 or there is interest
#  using something other than 10 or 25 points.
#
# Unlike rplot,a confidence band based on the Studentized maximum modulus dist
# is computed,
# unless alpha is not equal to .05 or the number of confidence intervals
# is greater than npts=28, in which case the distribution of the max of npts
#  random variables is used.
#
# LP=TRUE, the plot is further smoothed via lowess
#
# fr controls amount of smoothing
#
xord=order(x)
x=x[xord]
y=y[xord]
if(!is.null(pts))pts=sort(pts)
str=rplot(x,y,tr=tr,xout=xout,plotit=FALSE,LP=LP,fr=fr,pr=FALSE)$Strength.Assoc
m<-cbind(x,y)
if(ncol(m)>2)stop('To get a smooth with more than one covariate, use rplot')
m<-elimna(m)
nv=nrow(m)
if(eout && xout)stop('Not allowed to have eout=xout=T')
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(m[,1])$keep
m<-m[flag,]
}
if(is.null(pts)){
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr,nmin=nmin)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
}
x=m[,1]
y=m[,2]
n.keep=length(y)
if(is.null(pts)){
if(!MID)pts=seq(min(x),max(x),length.out=npts)
vv=idealf(x)
if(MID)pts=seq(vv$ql,vv$qu,length.out=npts)
}
rmd=NA
for(i in 1:length(pts))rmd[i]<-mean(y[near(x,pts[i],fr)],tr=tr)
sedf=runse(x,y,fr=fr,tr=tr,pts=pts,SEED=SEED)
df=sedf$df
flag=df>dfmin
se=sedf$se
ntest=length(df[flag])
mdif=min(df[flag])
crit=NA
dfval=df[flag]
for(it in 1:ntest)crit[it]=qsmm(1-alpha,ntest,dfval[it])
low=rmd[flag]-crit*se[flag]
up=rmd[flag]+crit*se[flag]
ptsall=pts
rmdall=rmd
rmd=rmd[flag]
pts=pts[flag]
if(plotit){
ord=order(x)
x=x[ord]
y=y[ord]
if(LP){
rmd=lplot(pts,rmd,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
up=lplot(pts,up,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
low=lplot(pts,low,plotit=FALSE,pyhat=TRUE,pr=FALSE)$yhat
}
plot(c(x,pts),c(y,rmd),xlab=xlab,ylab=ylab,type='n')
if(scat)points(x,y,pch=pch)
lines(pts,up,lty=2)
lines(pts,low,lty=2)
lines(pts,rmd)
}
if(pyhat){output<-cbind(pts,rmd,low,up)
dimnames(output)=list(NULL,c('pts','y.hat','ci.low','ci.up'))
}
if(!pyhat)output<-'Done'
list(output=output,str=str,n=nv,n.keep=n.keep)
}

rplotCITAP.pv<-function(n,nreps=2000,alpha=.05,npts=25,tr=.2,fr=5,MC=FALSE,nmin=12,SEED=TRUE,LP=FALSE){
if(SEED)set.seed(2)
pvals=NA
xy=list()
for (i in 1:nreps){
xy[[i]]=rmul(n)
}
if(!MC)pvals=lapply(xy,rplotCITAP.sub,npts=npts,tr=tr,fr=fr,alpha=alpha,nmin=nmin)
if(MC){
pvals=mclapply(xy,rplotCITAP.sub,npts=npts,tr=tr,fr=fr,alpha=alpha,nmin=nmin)
}
pvals=matl(pvals)
pv=hd(pvals,alpha)
pv
}

rplotCITAP.sub<-function(xy,tr=.2,fr=NA,SEED=TRUE,nmin=12,
pts=NA,npts=25,LP=FALSE,alpha=.05,xout=FALSE,...){
#
# prediction interval running  interval smoother based on a trimmed mean.
# Unlike rplot, includes a confidence band.
#
x=xy[,1]
y=xy[,2]
xord=order(x)
x=x[xord]
y=y[xord]
infit=rplot(x,y,tr=tr,xout=xout,plotit=FALSE,LP=LP,fr=fr,pr=FALSE,pyhat=TRUE,nmin=nmin)
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
rmd=infit$pyhat
res1=ancova(x,y,x,y,pr=FALSE,plotit=FALSE,fr1=fr,fr2=fr,nmin=nmin)$output
pts=seq(res1[1,1],res1[5,1],length.out=npts)
y.hat=NA
pv=NA
civ=matrix(NA,nrow=npts,ncol=2)
for(i in 1:length(pts)){
doit=trimci(y[near(x,pts[i],fr)],tr=tr,alpha=alpha,pr=FALSE)
pv[i]=doit$p.value
}
min(pv)
}

rplotCIv2.pv<-function(n,nreps=4000,alpha=.05,tr=.2,fr=.5,
MC=TRUE,nmin=12,SEED=TRUE){
if(SEED)set.seed(2)
pvals=NA
xy=list()
for (i in 1:nreps){
xy[[i]]=rmul(n)
}
if(!MC)pvals=lapply(xy,rplotCIv2.sub,tr=tr,fr=fr,nmin=nmin)
if(MC){
pvals=mclapply(xy,rplotCIv2.sub,tr=tr,fr=fr,nmin=nmin)
}
pvals=matl(pvals)
pv=hd(pvals,alpha)
pv
}

rplotCIv2.sub<-function(xy,nmin,tr,fr){
x=xy[,1]
y=xy[,2]
n=length(y)
nv=NA
for(j in 1:n)nv[j]=sum(near(x,x[j],fr=fr))
pts=x[nv>=nmin]
n.keep=length(pts)
for(j in 1:n.keep)nv[j]=sum(near(x,x[j],fr=fr))
pts=x[nv>=nmin]
rmd=NA
for(i in 1:length(pts))rmd[i]<-trimci(y[near(x,pts[i],fr)],tr=tr,pr=FALSE)$p.value
pv=min(rmd)
pv
}

rplotCV<-function(x,y,fr=NA,varfun=pbvar,est=tmean,xout=FALSE,outfun=out,eout=FALSE,corfun=pbvar,...){
#
# Estimate prediction error based on
# a running interval smoother in conjunction with
# a leave-one-out cross validation method
#
#  varfun is the measure of variation used on the predicted Y values.
#  est is the measure of location used by the running interval smoother.
#  The estimate is returned in VAR.Y.HAT
#
m=elimna(cbind(x,y))
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m=m[flag,]
}
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=as.matrix(m[,1:p])
y=m[,p1]
vals=NA
if(is.na(fr)){
if(p==1)fr=.8
if(p>1)fr=1
}
if(xout){
keepit<-outfun(x,plotit=FALSE,...)$keep
x<-x[keepit,]
y<-y[keepit]
}
x=as.matrix(x)
for(i in 1:nrow(x)){
if(p==1)vals[i]=runhat(x[-i,],y[-i],fr=fr,est=est,pts=x[i,],...)
if(p>1)vals[i]=rung3hat(x[-i,],y[-i],fr=fr,pts=t(as.matrix(x[i,])))$rmd
}
dif=y-vals
ans=varfun(elimna(dif))
list(VAR.Y.HAT=ans)
}

rplotsm<-function(x,y,est=tmean,fr=1,plotit=TRUE,pyhat=FALSE,nboot=40,atr=0,nmin=0,pch='*',
outfun=outpro,eout=FALSE,xlab="X",ylab="Y",scat=TRUE,SEED=TRUE,expand=.5,scale=FALSE,STAND=TRUE,
varfun=pbvar,pr=TRUE,ticktype="simple",theta=50,phi=25,...){
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
if(ncol(x)==1){
val<-runmbo(x,y,est=est,scat=scat,fr=fr,plotit=plotit,pyhat=TRUE,STAND=STAND,pch=pch,
xlab=xlab,ylab=ylab,eout=eout,nboot=nboot,outfun=outfun,SEED=SEED,atr=atr,...)
}
if(ncol(x)>1){
if(ncol(x)==2 && !scale){
if(pr){
print("scale=F is specified.")
print("If there is dependence, use scale=T")
}}
if(ncol(x)>2)plotit<-F
val<-run3bo(x,y,est=est,fr=fr,nmin=nmin,plotit=plotit,pyhat=TRUE,phi=phi,
theta=theta,xlab=xlab,ylab=ylab,ticktype=ticktype,STAND=STAND,
SEED=SEED,expand=expand,scale=scale,nboot=nboot,...)
val<-val$output
}
E.power<-varfun(val[!is.na(val)])/varfun(y)
if(!pyhat)val <- NULL
E.power=as.numeric(E.power)
list(Strength.Assoc=sqrt(E.power),Explanatory.Power = E.power, yhat = val)
}

rplotN<-function(x,y,nsub=5000,est=tmean,fr=1,xout=FALSE,xlab='X',ylab='Y',zlab='',ticktype = 'simple',theta = 50, phi = 25, scale = TRUE, 
    expand = 0.5, SEED = TRUE,frame=TRUE){
#
 # Running interval smoother, good for large sample sizes or plots of the
 # regression surface without a scatter plot.
 #
 # nsub is size of the random sample of the data used to predict outcome using all of the data
 #  
 if(SEED)set.seed(2)
 x=as.matrix(x)
p=ncol(x)
p1=p+1
 xy=cbind(x,y)
 xy=elimna(xy)
 if(xout){
 flag<-outpro.depth(x,plotit=FALSE)$keep                                
x<-x[flag,]                                                                        
x<-as.matrix(x)                                                                    
y<-y[flag]                                                                         
} 
 n=nrow(xy)
 nsub=min(n,nsub)
 id=sample(n,nsub) 
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
w=rplot.pred(x,y,pts=x[id,],fr=fr)$Y.hat 
a=lplot(x[id,],w,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,frame=frame,phi=phi,theta=theta,scale=scale,pr=FALSE)   
}

rplot.bin<-function(x,y,est=mean,scat=TRUE,fr=NULL,plotit=TRUE,pyhat=FALSE,pts=x,LP=FALSE,
theta=50,phi=25,scale=TRUE,expand=.5,SEED=TRUE,
nmin=0,xout=FALSE,outfun=outpro,xlab=NULL,ylab=NULL,
zlab='P(Y=1)',pr=TRUE,duplicate='error',...){
#
#  This function applies the running interval smoother, but is designed
#  specifically for situations where y is  binary.
#
# duplicate='error'
# In some situations where duplicate values occur, when plotting with
# two predictors, it is necessary to set duplicate='strip'
#
y=chbin2num(y)
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
x<-as.matrix(x)
if(length(unique(y))!=2)stop('y is not binary')
n=length(y)
Y=rep(0,n)
flag=y==max(y)
Y[flag]=1
y=Y
x<-as.matrix(x)
if(ncol(x)==1){
if(is.null(ylab))ylab='P(Y=1)'
if(is.null(xlab))ylab='X'
if(is.null(fr))fr=.8
a=rplot(x,y,est=mean,xout=xout,outfun=outfun,fr=fr,xlab=xlab,ylab=ylab,pr=FALSE,LP=LP)
}
if(ncol(x)>1){
id=chk4binary(x)
Lid=length(id)
if(Lid>0)stop('Binary independent variables detected. Use rplot.binv2')
if(is.null(xlab))xlab='X1'
if(is.null(ylab))ylab='X2'
if(is.null(fr))fr=1.2
if(ncol(x)==2){
if(scale){
if(pr){print('scale=T is specified.')
print('If there is independence, might want to use scale=F')
a=rplot(x,y,est=mean,xout=xout,outfun=outfun,fr=fr,xlab=xlab,ylab=ylab,zlab=zlab,scale=scale,pr=FALSE)
}}}}
if(!pyhat)val <- 'DONE'
if(pyhat)val=rplot.pred(x,y,pts=pts,est=mean,xout=xout,outfun=outfun,fr=fr)
val
}

rplot.binCI<-function(x,y,pts=NULL,alpha=.05,nmin=5,xout=FALSE,outfun=outpro,fr=.5,tr.plot=FALSE,
method=NULL,plotit=TRUE,LP=TRUE,xlab='X',ylab='P(Y=1|X)',...){
#
#  An alternative to logistic regression.
#
#  For a collection of intervals among the values in
#  x, compute the probability of success and a confidence based on the corresponding y values
#
#  Default: use the deciles to define the intervals
#
#  Example: pts=c(-1,0,1,2). The intervals would be (-1,0), (0,1), (1,2).
#
y=chbin2num(y)
xx<-elimna(cbind(x,y))
x<-xx[,1]
y<-xx[,2]
if(is.null(pts)){
id=duplicated(x)
pts=x[!id]
}
else plotit=FALSE
if(xout){
m<-cbind(x,y)
if(ncol(m)!=2)stop('Only one  explanatory variable is allowed')
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1]
y<-m[,2]
}
n=length(x)
xor=order(x)
x=x[xor]
y=y[xor]
pts=sort(pts)
npts=length(pts)
#
if(length(unique(y))>2)stop('y should be binary')
#
# Determine which method will be used:
if(is.null(method)){
if(n<80)method='AC'
if(n>=80)method='CP'
}

# Next convert y to 0 and 1s  if not already 0 and 1s
yy=rep(0,n)
yc=as.character(y)
flag=yc==max(yc)
yy[flag]=1
y=yy
#
rmd<-matrix(NA,nrow=npts,ncol=7)
for(i in 1:npts){
isub=near(x,pts[i],fr)
if(sum(isub)>=nmin){
z=y[isub]
v=binom.conf(y=z,method=method,alpha=alpha,pr=FALSE)
rmd[i,1]=v$n
rmd[i,2]=min(x[isub])
rmd[i,3]=max(x[isub])
rmd[i,4]=pts[i]
rmd[i,5]=v$phat
rmd[i,6:7]=v$ci
}}
rs=elimna(rmd)
dimnames(rmd)=list(NULL,c('n','low.end','upper.end','pts','p.hat','ci.low','ci.up'))
if(plotit){
if(tr.plot){
v=quantile(rmd[,4],probs=c(.1,.9),na.rm=TRUE)
flag=(rmd[,4]>=v[1] & rmd[,4]<=v[2])
rmd=rmd[flag,]
}
ys=rs[,5]
plot(rs[,4],rs[,5],ylim=c(0,1),xlab=xlab,ylab=ylab,type='n')
if(LP){
z1=lplot.pred(rmd[,4],rmd[,5],pts=rmd[,2])$yhat
flag=z1>1
z1[flag]=1
flag=z1<0
z1[flag]=0
z2=lplot.pred(rmd[,4],rmd[,6],pts=rmd[,2])$yhat
flag=z2>1
z2[flag]=1
flag=z2<0
z2[flag]=0
z3=lplot.pred(rmd[,4],rmd[,7],pts=rmd[,2])$yhat
flag=z3>1
z3[flag]=1
flag=(z3<0)
z3[flag]=0
lines(rmd[,4],z1)
lines(rmd[,4],z2,lty=2)
lines(rmd[,4],z3,lty=2)
}
if(!LP){
lines(rmd[,4],rmd[,5])
lines(rmd[,4],rmd[,6],lty=2)
lines(rmd[,4],rmd[,7],lty=2)
}}
id=duplicated(rmd[,2:3])
rmd=elimna(rmd[!id,])
rmd
}

reg.vs.rplot<-function(x,y,xout=FALSE,fr=1,est=median,regfun=Qreg,Qreg.plot=TRUE,qv=c(.25,.75),SMQ=FALSE,
pr=TRUE,xlab='Reg.Est',ylab='Rplot.Est',pch='*'){
#
# If the linear model is correct, the plot returned here should be
# tightly clustered around a line having slope=1 and intercept=0, indicated
# by a dashed line.
#
if(pr)print('This function was updated July 2022')
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
x=xy[,1:p]
y=xy[,p1]
e1=regYhat(x,y,xout=xout,regfun=regfun)
e2=rplot.pred(x,y,xout=xout,est=est,fr=fr)$Y.hat
if(Qreg.plot){
if(!SMQ)qplotreg(e1,e2,xlab=xlab,ylab=ylab,pch=pch,qval=qv)
if(SMQ)qhdsm(e1,e2,xlab=xlab,ylab=ylab,pch=pch,qval=qv,LP=TRUE)
}
if(!(Qreg.plot)) lplot(e1,e2,xlab=xlab,ylab=ylab,pc=pch)
abline(0,1,lty=2)
}

reg.vs.lplot<-function(x,y,xout=FALSE,Qreg.plot=TRUE,qv=c(.25,.75),SMQ=FALSE,pch='*',pr=TRUE,
outfun=outpro,fr=1,est=mean,regfun=tsreg,xlab='Reg.est',ylab='Lplot.est',span=.75,...){
#
#
#
if(pr)print('This function was updated July 2022')
xy=elimna(cbind(x,y))
p1=ncol(xy)
p=p1-1
x=xy[,1:p]
y=xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
e1=regYhat(x,y,regfun=regfun)
e2=lplot.pred(x,y,,est=est,span=span)$yhat
if(Qreg.plot){
if(!SMQ)qplotreg(e1,e2,xlab=xlab,ylab=ylab,pch=pch,qval=qv)
if(SMQ)qhdsm(e1,e2,xlab=xlab,ylab=ylab,pch=pch,qval=qv,LP=TRUE)
}
if(!(Qreg.plot)) lplot(e1,e2,xlab=xlab,ylab=ylab,pc=pch)
abline(0,1,lty=2)
}

reg2difplot<-function(x1,y1,x2,y2,regfun=tsreg,pts=x1,xlab="VAR 1",ylab="VAR 2",zlab="Group 2 minus Group 1",xout=FALSE,outfun=out,ALL=TRUE,pts.out=FALSE,SCAT=FALSE,theta=50,phi=25,ticktype='simple',
pr=TRUE,...){
#
# Fit a regression model to both groups assuming have two predictors.
# Get predicted Y values based on points in pts. By default, use
# pts=x1
#
#  x1 a matrix containing two predictors
#  x2 a matrix containing two predictors
#
#  Compute differences in predicted values and plot the results as a function of the points in pts
#  pts=x1 by default.
#  ALL=T, pts is taken to be all points in x1 and x2.
#
#  pts.out=T will remove leverage points from pts.
#
if(!is.matrix(x1))stop("x1 should be a matrix")
if(!is.matrix(x2))stop("x2 should be a matrix")
if(!is.matrix(pts))stop("pts should be a matrix")
if(ncol(x1)!=2)stop("x1 should be a matrix with two columns")
if(ncol(x2)!=2)stop("x2 should be a matrix with two columns")
if(ncol(pts)!=2)stop("pts should be a matrix with two columns")
if(ALL)pts=rbind(x1,x2)
if(pts.out){
flag=outfun(pts,plotit=FALSE,...)$keep
pts=pts[flag,]
}
e1=regYhat(x1,y1,xout=xout,regfun=regfun,outfun=outfun,xr=pts,...)
e2=regYhat(x2,y2,xout=xout,regfun=regfun,outfun=outfun,xr=pts,...)
if(SCAT){
library(scatterplot3d)
scatterplot3d(cbind(pts,e2-e1),xlab=xlab,ylab=ylab,zlab=zlab)
}
if(!SCAT)rplot(pts,e2-e1,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,pr=FALSE,ticktype=ticktype,prm=FALSE)
}

qplotreg<-function(x, y,qval=c(.2,.8),q=NULL,plotit=TRUE,xlab="X",ylab="Y",xout=FALSE,
outfun=outpro,pch='*',...){
#
# Compute the quantile regression line for each of the
# quantiles indicated by qval.
# plotit=TRUE, plot the results.
#
if(!is.null(q))qval=q
xy=elimna(cbind(x,y))
if(ncol(xy)>2)stop("Only One Predictor Allowed")
x=xy[,1]
y=xy[,2]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
}
n<-length(qval)
coef<-matrix(NA,ncol=2,nrow=n)
x<-as.matrix(x)
if(ncol(x)>1)stop("This version allows one predictor only.")
if(plotit)plot(x,y,xlab=xlab,ylab=ylab,pch=pch)
for(it in 1:n){
coef[it,]<-qreg(x,y,qval=qval[it],pr=FALSE)$coef
dimnames(coef)=list(NULL,c("Inter.","Slope"))
if(plotit)abline(coef[it,1],coef[it,2])
}
coef
}

qregplots<-function(x, y,qval=.5,q=NULL,op=1,pr=FALSE,xout=FALSE,outfun=out,plotit=FALSE,xlab="X",ylab="Y",...){
#
# Compute the quantile regression line for one or more quantiles and plot the results
# That is, determine the qth (qval) quantile of Y given X using the
#  the Koenker-Bassett approach.
#
#  One predictor only is allowed
#
#  v2=T, uses the function rq in the R library quantreg
#  v2=F, uses an older and slower version
#
#  Example: qregplots(x,y,q=c(.25,.5,.75)) will plot the regression lines for
#   predicting quartiles.
#
if(!is.null(q))qval=q
x<-as.matrix(x)
if(ncol(x)!=1)stop("Current version allows only one predictor")
X<-cbind(x,y)
X<-elimna(X)
np<-ncol(X)
p<-np-1
x<-X[,1:p]
x<-as.matrix(x)
y<-X[,np]
if(xout){
x<-as.matrix(x)
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
est=matrix(NA,ncol=3,nrow=length(qval))
dimnames(est)=list(NULL,c("q","Inter","Slope"))
library(quantreg)
x<-as.matrix(x)
plot(x,y,xlab=xlab,ylab=ylab)
if(ncol(x)!=1)stop("Current version allows only one predictor")
for(j in 1:length(qval)){
coef=coefficients((rq(y~x,tau=qval[j])))
est[j,1]=qval[j]
est[j,2:3]=coef
abline(coef)
}
list(coef = est)
}

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

prplot<-function(x,y,pval=ncol(x),regfun=tsreg,fr=.8,est=tmean,op=1,
xlab="X",ylab="Residuals",xout=FALSE,outfun=out,...){
#
# Goal: check for curvature associated with predictor
# indicated by pval.
# This is done by creating a partial residual plot.
# That is subtracting out the linear prediction based
# on the other predictors and then
# smooth the result versus the predictor in the column of x indicated by pval
#
x=as.matrix(x)
p=ncol(x)
p1=p+1
temp=elimna(cbind(x,y))
x=temp[,1:p]
y=temp[,p1]
if(xout){
flag<-outfun(x,...)$keep
x<-as.matrix(x)
x<-x[flag,]
y<-y[flag]
x<-as.matrix(x)
}
if(!is.matrix(x))stop("Should have two or more variables stored in a matrix")
flag<-rep(TRUE,ncol(x))
flag[pval]<-FALSE
temp<-regfun(x[,flag],y)$residual
if(op!=1)rungen(x[,!flag],temp,est=est,fr=fr,xlab=xlab,ylab=ylab,...)
if(op==1)lplot(x[,!flag],temp,xlab=xlab,ylab=ylab)
}

MLRplot<-function(x, Y)
{
# Forward response plot and residual plot.
# R needs command "library(lqs)" if a robust estimator replaces lsfit.
# Advance the view with the right mouse button.
	x <- as.matrix(x)
	out <- lsfit(x, Y)
	cook <- ls.diag(out)$cooks
	n <- dim(x)[1]
	p <- dim(x)[2] + 1
	tem <- cook > min(0.5, (2 * p)/n)
	bhat <- out$coef
	FIT <- bhat[1] + x %*% bhat[-1]
	par(mfrow = c(2, 1))
	plot(FIT, Y)
	abline(0, 1)
	points(FIT[tem], Y[tem], pch = 15)
	identify(FIT, Y)
	title("Forward Response Plot")
	RES <- Y - FIT
	plot(FIT, RES)
	points(FIT[tem], RES[tem], pch = 15)
	identify(FIT, RES)
	title("Residual Plot")
}


# ============================================================================
# LOWESS/LOESS PLOTS
# ============================================================================

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
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n.orig=nrow(m)
n.keep=n.orig
if(!is.matrix(x))stop("x is not a matrix")
d<-ncol(x)
if(d>=2){
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

lplotv2<-function(x,y,span=.75,pyhat=FALSE,eout=FALSE,xout=FALSE,outfun=out,plotit=TRUE,
expand=.5,low.span=2/3,varfun=pbvar,cor.op=FALSE,cor.fun=pbcor,ADJ=FALSE,nboot=20,
scale=TRUE,xlab="X",ylab="Y",zlab="",theta=50,phi=25,family="gaussian",
duplicate="error",pr=TRUE,SEED=TRUE,ticktype="simple"){
#
# Plot regression surface using LOESS
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
st.adj=NULL
e.adj=NULL
if(ADJ){
if(SEED)set.seed(2)
}
si=1
x<-as.matrix(x)
if(!is.matrix(x))stop("x is not a matrix")
d<-ncol(x)
if(d>=2){
if(ncol(x)==2 && !scale){
if(pr){
print("scale=F is specified.")
print("If there is dependence, might use scale=T")
}}
m<-elimna(cbind(x,y))
x<-m[,1:d]
y<-m[,d+1]
if(eout && xout)stop("Can't have both eout and xout = F")
if(eout){
flag<-outfun(m,plotit=FALSE)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(x,plotit=FALSE)$keep
m<-m[flag,]
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
persp(fitr,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,expand=expand,
scale=scale,ticktype=ticktype)
}}
if(d==1){
m<-elimna(cbind(x,y))
x<-m[,1:d]
y<-m[,d+1]
if(eout && xout)stop("Can't have both eout and xout = F")
if(eout){
flag<-outfun(m)$keep
m<-m[flag,]
}
if(xout){
flag<-outfun(x)$keep
m<-m[flag,]
}
x<-m[,1:d]
y<-m[,d+1]
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab)
lines(lowess(x,y,f=low.span))
}
yyy<-lowess(x,y)$y
xxx<-lowess(x,y)$x
if(d==1){
ordx=order(xxx)
yord=yyy[ordx]
flag=NA
for (i in 2:length(yyy))flag[i-1]=sign(yord[i]-yord[i-1])
if(sum(flag)<0)si=-1
}
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
E.power<-1
if(!cor.op)E.power<-varfun(last[!is.na(last)])/varfun(y)
if(cor.op || E.power>=1){
if(d==1){
xord<-order(x)
E.power<-cor.fun(last,y[xord])$cor^2
}
if(d>1)E.power<-cor.fun(last,y)$cor^2
}
if(ADJ){
x=as.matrix(x)
val=NA
n=length(y)
data1<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(i in 1:nboot){
temp=lplot.sub(x[data1[i,],],y[data2[i,]],plotit=FALSE,pr=FALSE)
val[i]=temp$Explanatory.power
}
vindt=median(val)
v2indt=median(sqrt(val))
st.adj=(sqrt(E.power)-max(c(0,v2indt)))/(1-max(c(0,v2indt)))
e.adj=(E.power-max(c(0,vindt)))/(1-max(c(0,vindt)))
st.adj=max(c(0,st.adj))
e.adj=max(c(0,e.adj))
}
if(!pyhat)last <- NULL
list(Strength.Assoc=si*sqrt(E.power),Explanatory.power=E.power,
Strength.Adj=st.adj,Explanatory.Adj=e.adj,yhat.values=last)
}

lplot2g<-function(x1,y1,x2,y2,fr=.8,est=tmean,xlab="X",ylab="Y",xout=FALSE,eout=FALSE,pr=TRUE,
pch1='*',pch2='o',outfun=outpro,...){
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
if(pr){
if(!xout)print('Suggest also looking at result using xout=TRUE')
}
m1<-elimna(cbind(x1,y1))
if(eout && xout)stop("Can't have both eout and xout = F")
if(eout){
flag<-outfun(m1,plotit=FALSE,...)$keep
m1<-m1[flag,]
}
if(xout){
flag<-outfun(m1[,1],plotit=FALSE,...)$keep
m1<-m1[flag,]
}
x1<-m1[,1]
y1<-m1[,2]
m2<-elimna(cbind(x2,y2))
if(eout){
flag<-outfun(m2,plotit=FALSE,...)$keep
m2<-m2[flag,]
}
if(xout){
flag<-outfun(m2[,1],plotit=FALSE,...)$keep
m2<-m2[flag,]
}
x2<-m2[,1]
y2<-m2[,2]

flag=order(x1)
x1=x1[flag]
y1=y1[flag]
flag=order(x2)
x2=x2[flag]
y2=y2[flag]
temp1<-lplot(x1,y1,pyhat=TRUE,plotit=FALSE,pr=FALSE)$yhat.values
temp2<-lplot(x2,y2,pyhat=TRUE,plotit=FALSE,pr=FALSE)$yhat.values
plot(c(x1,x2),c(y1,y2),type="n",xlab=xlab,ylab=ylab)
points(x1,y1,pch=pch1)
points(x2,y2,pch=pch2)
lines(x1,temp1)
lines(x2,temp2,lty=2)
}

lplotCI<-function(x,y,plotit=TRUE,xlab='X',ylab='Y',p.crit=NULL,alpha=.05,span=NULL,
CIV=FALSE,xout=FALSE,outfun=outpro, pch='.',SEED=TRUE,nboot=100,pts=NULL,npts=25,nreps=2000,...){
#
#  Confidence band using LOESS
#
#  Method allows heteroscedasticity and adjust the confidence intervals
# so that the simultaneous probabillty coverage is approximately 1-alpha
#
# If CIV=FALSE and plotit=TRUE, creates a plot with the confidence intervals.
# CIV=TRUE, returns the confidence intervals for the points in pts
# pts =NULL, the function picks
# npts points, extending between M-1.5*mad(x) and M+1.5*mad(x)
#
#
# For alpha=0.05, n <=2000 execution time is low. Otherwise
# the adjusted critical value must be computed.
#
#  p.crit=NULL: If alpha=.05, determined quickly, otherwise it is computed.
#
xy=elimna(cbind(x,y))
if(ncol(xy)>2)stop('Current version limited to a single predictor variable')
if(xout){
flag<-outfun(xy[,1],plotit=FALSE,...)$keep
xy<-xy[flag,]
}
n=nrow(xy)
if(is.null(span)){
span=2/3
if(n >=300)span=.5
if(n >=800)span=.3
}
x=xy[,1]
y=xy[,2]
xord=order(x)
y=y[xord]
x=x[xord]
M=median(x)
low=M-1.5*mad(x)
up=M+1.5*mad(x)
if(is.null(pts))pts=seq(low,up,length.out=npts)
if(npts<=5)p.crit=alpha/npts
if(alpha==.05){
if(is.null(p.crit)){
if(n<30)stop('Should have n>=30')
nv=c(30,50,70,100,150, 200,300, 500, 1000, 2000)
pv=c(0.003599898, 0.002661925, 0.002399994, 0.002877103, 0.003000428, 0.003538190,
 0.003872710, 0.004396500, 0.004075000, 0.0045161)

if(npts==25){
if(n<=2000)p.crit=lplot.pred(1/nv,pv,pts=1/n)$yhat
if(n>2000)p.crit=.00452
}}}
if(is.null(p.crit)){
print('p.crit is being computed, this might take some time.')
pts.stand=NULL
if(!is.null(pts))pts.stand=(median(x)-pts)/mad(x)
p.crit=lplotbsepvv3(n,nreps=nreps,npts=npts,pts=pts.stand,alpha=alpha)
}
plx<-predict(loess(y ~ x,span=span), se=TRUE)
se=lplotse(x,y,nboot=nboot,SEED=SEED,pts=pts,span=span)
lfit=lplot.pred(x,y,pts=pts,span=2/3)$yhat
if(plotit){
plot(x,y,xlab=xlab,ylab=ylab,pch=pch)
lines(x,plx$fit)
if(is.null(p.crit))p.crit=alpha
lines(pts,lfit - qt(1-p.crit/2,plx$df)*se, lty=2)
lines(pts,lfit + qt(1-p.crit/2,plx$df)*se, lty=2)
}
ci.low=lfit - qt(1-p.crit/2,plx$df)*se
ci.up=lfit + qt(1-p.crit/2,plx$df)*se
if(!CIV)ci=NULL
if(CIV){
ci=cbind(pts,lfit,ci.low,ci.up)
dimnames(ci)=list(NULL,c('X','Y.hat','ci.low','ci.up'))
}
list(p.crit=p.crit,Conf.Intervals=ci)
}

lplotse<-function(x,y,pts=x,nboot=100,SEED=TRUE,span=2/3){
#
# compute estimae of SE
# return the values corresponding to the order x values
#
xord=order(x)
y=y[xord]
x=x[xord]
if(SEED)set.seed(2)
n=length(y)
ev=matrix(NA,nrow=nboot,ncol=length(pts))
for(i in 1:nboot){
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
ev[i,]=lplot.pred(x[data[i,]],y[data[i,]],pts=pts,span=span)$yhat
}
se=apply(ev,2,sd)
se
}

lplotPV<-function(x,y, span = 0.75, xout = FALSE,pr=TRUE,
    outfun = out,nboot=1000,SEED=TRUE,plotit=TRUE,pyhat = FALSE, expand = 0.5, low.span = 2/3,
    varfun = pbvar, cor.op = FALSE, cor.fun = pbcor, scale = FALSE,
    xlab = "X", ylab = "Y", zlab = "", theta = 50, phi = 25,
    family = "gaussian", duplicate = "error", pc = "*", ticktype = "simple",...){
#
# Compute a p-value based on the Strength of Association estimated via lplot
# If significant, conclude there is dependence.
#
if(SEED)set.seed(2)
x=as.matrix(x)
if(ncol(x)==2 && !scale){
if(pr){
print("scale=F is specified.")
print("If there is dependence, might use scale=T")
}}
vals=NA
nv=ncol(x)
m=elimna(cbind(x,y))
x<-m[,1:nv]
y<-m[,nv+1]
if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:nv]
y<-m[,nv+1]
}
x=as.matrix(x)
est=lplot(x,y,span=span,plotit=plotit,pr=FALSE, pyhat = pyhat,
    outfun = outfun, expand = expand, low.span = low.span,
    varfun = varfun, cor.op =cor.op, cor.fun = cor.fun, scale = scale,
    xlab = xlab, ylab = ylab, zlab =zlab, theta =theta, phi = phi,
    family = family, duplicate = duplicate, pc = pc, ticktype = ticktype,...)
n=nrow(x)
data1<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(i in 1:nboot){
vals[i]=lplot(x[data1[i,],],y[data2[i,]],plotit=FALSE,pr=FALSE)$Strength.Assoc
}
p=mean(est$Strength<vals)
list(p.value=p,Strength.Assoc=est$Strength.Assoc,Explanatory.power=est$Explanatory.power,yhat.values=est$yhat.values)
}

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

lplotcom2<-function(x,y,xout=FALSE,pts1=NULL,pts2=NULL,outfun=outpro,span=2/3,npts=10,tr=.2,...){
#
# For two independent variables, estimate their relative importance when using LOESS
#
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n.orig=nrow(m)
n.keep=n.orig
d<-ncol(x)
if(d!=2)stop('Current version is for two independent variables only')
if(xout){
flag<-outfun(m[,1:2],plotit=FALSE,...)$keep
m<-m[flag,]
}
n.keep=nrow(m)
M=apply(m,2,median)
SD=apply(m,2,mad)
low=M-1.5*SD
up=M+1.5*SD
if(is.null(pts1))pts1=seq(low[1],up[1],length.out=npts)
if(is.null(pts2))pts2=seq(low[2],up[2],length.out=npts)
e1=NA  #
e2=NA
for(j in 1:length(pts1)){     # Determine strength of x2 given a value stored in pts1.
v2=cbind(rep(pts1[j],n.keep),m[,2])
vals=lplot.pred(m[,1:2],m[,3],v2,span=span)$yhat
vals=elimna(vals)
nv=length(vals)
e2[j]=NA
if(nv>=10)e2[j]=winsd(vals,tr=tr,na.rm=TRUE)/winsd(m[,3],tr=tr,na.rm=TRUE)
}
for(j in 1:length(pts2)){     # Determine strength of x1 given a value stored in pts2.
v1=cbind(m[,1],rep(pts2[j],n.keep))
vals=lplot.pred(m[,1:2],m[,3],v1,span=span)$yhat
vals=elimna(vals)
nv=length(vals)
e1[j]=NA
if(nv>=10)e1[j]=winsd(vals,tr=tr,na.rm=TRUE)/winsd(m[,3],tr=tr,na.rm=TRUE)
 }
p=mean(outer(e1,e2,FUN='-')<0,na.rm=TRUE)
list(str1=e1,str2=e2,p=p,mean.str1=mean(e1),mean.str2=mean(e2))
}

lplotcom2v2<-function(x,y,xout=FALSE,pts1=NULL,pts2=NULL,outfun=outpro,span=2/3,npts=10,tr=.2,...){
#
# For two independent variables, estimate their relative importance when using LOESS
#
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n.orig=nrow(m)
n.keep=n.orig

d<-ncol(x)
if(d!=2)stop('Current version is for two independent variables only')
if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
} #$
n.keep=nrow(m)
M=apply(m,2,median)
SD=apply(m,2,mad)
low=M-1.5*SD
up=M+1.5*SD
if(is.null(pts1))pts1=seq(low[1],up[1],length.out=npts)
if(is.null(pts2))pts2=seq(low[2],up[2],length.out=npts)
e1=NA
e2=NA
for(j in 1:length(pts1)){     # Determine strength of x2 given a value stored in pts1.
v2=cbind(rep(pts1[j],n.keep),m[,2])
vals=lplot.pred(m[,1:2],m[,3],v2,span=span)$yhat
vals=elimna(vals)
nv=length(vals)
e2[j]=NA
if(nv>=10)e2[j]=winsd(vals,tr=tr,na.rm=TRUE)
}
for(j in 1:length(pts2)){     # Determine strength of x1 given a value stored in pts2.
v1=cbind(m[,1],rep(pts2[j],n.keep))
vals=lplot.pred(m[,1:2],m[,3],v1,span=span)$yhat
vals=elimna(vals)
nv=length(vals)
e1[j]=NA
if(nv>=10)e1[j]=winsd(vals,tr=tr,na.rm=TRUE)
 }
p=mean(outer(e1,e2,FUN='-')<0,na.rm=TRUE)
list(str1=e1,str2=e2,p=p,mean.str1=mean(e1),mean.str2=mean(e2))
}

lplotcomBCI<-function(x,y,xout=FALSE,pts1=NULL,pts2=NULL,p.crit=NULL,
outfun=outpro,span=2/3,npts=10,tr=.2,nboot=500,
SEED=TRUE,SEQ=FALSE,MAD.OP=FALSE,plotit=TRUE,ticktype='simple',
xlab='X1',ylab='X2',zlab='Y',reverse.x1=FALSE,reverse.x2=FALSE,pr=FALSE,
MEDIAN=FALSE,Q1=FALSE,Q2=FALSE,alpha=.05,MC=FALSE,...){
#
# For two independent variables, estimate their relative importance when using LOESS
# p.crit is the critical p-value. If not specified, the function returns the approximate 0.05 critical p-value
#
# By default, use the average of the strength of the associations, so essentially a global test based on the quartiles
#  MEDIAN=TRUE, use the median of the independent variables only.
#  Q1=TRUE, use the lower quartile of the independent variables only.
#  Q2=TRUE, use the upper quartile of the independent variables only.
#
#  ADJ.CI=TRUE: Confidence intervals are based on the critical p-value
#  otherwise use alpha
#
if(SEED)set.seed(2)
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n=nrow(m)
x=m[,1:2]
y=m[,3]
if(xout){
flag=outfun(x,plotit=FALSE)$keep
x=x[flag,]
y=y[flag]
n=nrow(x)
}
if(n<50)stop('The sample size must be greater than or equal to 50')
if(MEDIAN){
pts1=median(x[,1])
pts2=median(x[,2])
if(is.null(p.crit)){
if(n<=200)p.crit=regYhat(c(1/50,1/100,1/200),c(.114,.080,.065),1/n)
if(n>200)p.crit=.062
}
}
if(Q1){
pts1=qest(x[,1],.25)
pts2=qest(x[,2],.25)
if(is.null(p.crit)){
if(n<=200)p.crit=regYhat(c(1/50,1/100,1/200),c(.142,.095,.082),1/n)
if(n>200)p.crit=.062
}
}
if(Q2){
pts1=qest(x[,1],.75)
pts2=qest(x[,2],.75)
if(is.null(p.crit)){
if(n<=200)p.crit=regYhat(c(1/50,1/100,1/200),c(.142,.095,.082),1/n)
if(n>200)p.crit=.062
}
}
if(is.null(pts1)){
pts1=qest(x[,1],c(.25,.5,.75))
if(is.null(p.crit)){
if(n<=200)p.crit=regYhat(c(1/50,1/100,1/200),c(.082,.076,.067),1/n)
if(n>200)p.crit=.06
}
if(reverse.x1)pts1=sort(pts1,TRUE)
if(is.null(pts2))pts2=qest(x[,2],c(.25,.5,.75))
if(reverse.x2)pts2=sort(pts2,TRUE)
}
if(SEQ){
if(MAD.OP){
M=apply(m,2,median)
SD=apply(m,2,mad)

low=M-1.5*SD
up=M+1.5*SD
}
else{
low=apply(m,2,qest,.25)
hi=apply(m,2,qest,.75)
}
pts1=seq(low[1],up[1],length.out=npts)
pts2=seq(low[2],up[2],length.out=npts)
}
v1=NA
v2=NA
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
if(!MC){
for(i in 1:nboot){
ib=data[i,]
temp=lplotcom2v2(x[ib,],y[ib],pts1=pts1,pts2=pts2,npts=npts,tr=tr,span=span)
v1[i]=temp$mean.str1
v2[i]=temp$mean.str2
}}
if(MC){
data=listm(t(data))
bvec<-mclapply(data,lplotCIMCv2,x,y,pts1=pts1,pts2=pts2,npts=npts,tr=tr,span=span)
bvec=matl(bvec)  # a 2-by-nboot matrix.
dif=sort(bvec[1,]-bvec[2,])
}
if(!MC)dif=sort(v1-v2)
nbl=length(dif)
#ilow<-round((alpha/2) * nbl)
ilow<-round((p.crit/2) * nbl)
ihi<-nbl - ilow
ilow<-ilow+1
ci.low=dif[ilow]
ci.hi=dif[ihi]
pv=mean(dif<0,na.rm=TRUE)
pv=2*min(pv,1-pv)
est=lplotcom2(x,y,xout=FALSE,pts1=pts1,pts2=pts2,outfun=outfun,span=span,
npts=npts,tr=tr)
if(plotit)lplot(x,y,ticktype=ticktype,xlab=xlab,ylab=ylab,zlab=zlab,pr=pr)
list(p.crit=p.crit,p.value=pv,str.x1.given.x2=est$str1,str.x2.given.x1=est$str2,mean.str1=est$mean.str1,
mean.str2=est$mean.str2,
ci.low=ci.low,ci.hi=ci.hi,pts.x1=pts1,pts.x2=pts2)
}

lplotcomBCI9<-function(x,y,xout=FALSE,pr=TRUE,
outfun=outpro,span=2/3,npts=10,tr=.2,nboot=500,
SEED=TRUE,plotit=TRUE,ticktype='simple',ADJ.CI=TRUE,
xlab='X1',ylab='X2',zlab='Y',alpha=.05,MC=FALSE,...){
#
# For two independent variables, estimate their relative importance when using LOESS
# Focus  on the quartiles: none tests based on all possible combinations.
#
p.crit=NA
if(pr){
if(alpha!=.05){
if(pr)print('Critical p-value is taken to be the value of alpha. Unknown how to adjust  when alpha is not .05')
p.crit=alpha
}
if(ADJ.CI)print('Confidence intervals are based on the critical p-value')
}
if(SEED)set.seed(2)
x<-as.matrix(x)
m<-elimna(cbind(x,y))
n=nrow(m)
n.orig=n
x=m[,1:2]
y=m[,3]
if(xout){
flag=outfun(x,plotit=FALSE)$keep
x=x[flag,]
y=y[flag]
n=nrow(x)
}
if(is.na(p.crit)){
if(n<=100)p.crit=regYhat(c(1/50,1/100),c(.042,.025),1/n)
else p.crit=.025
}
output<-matrix(NA,nrow=9,ncol=7)
dimnames(output)=list(NULL,c('pts1','pts2','p-value','str.x1.given.x2','str.x2.given.x1','ci.low','ci.hi'))
pts1=qest(x[,1],c(.25,.5,.75))
pts2=qest(x[,2],c(.25,.5,.75))

v1=matrix(NA,nrow=nboot,ncol=3)
v2=matrix(NA,nrow=nboot,ncol=3)

data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
if(!MC){
for(i in 1:nboot){
ib=data[i,]
for(j in 1:3){
temp=lplotcom2v2(x[ib,],y[ib],pts1=pts1[j],pts2=pts2[j],npts=npts,tr=tr,span=span)
v1[i,j]=temp$mean.str1
v2[i,j]=temp$mean.str2
}}}
if(MC){
data=listm(t(data))
for(j in 1:3){
bvec<-mclapply(data,lplotCIMCv2,x,y,pts1=pts1[j],pts2=pts2[j],npts=npts,tr=tr,span=span)
bvec=matl(bvec)  # a 2-by-nboot matrix.
v1[,j]=bvec[1,]
v2[,j]=bvec[2,]
}
}
pc=matrix(NA,3,3) #rows are for pts1, columns for pts2
ic=0
for(j in 1:3){
for(k in 1:3){
est=lplotcom2(x,y,xout=FALSE,pts1=pts1[j],pts2=pts2[k],outfun=outfun,span=span,
npts=npts,tr=tr)
ic=ic+1
output[ic,1]=pts1[j]
output[ic,2]=pts2[k]
dif=sort(v1[,j]-v2[,k])
nbl=length(dif)
if(ADJ.CI)ilow<-round((p.crit/2) * nbl)
else ilow<-round((alpha/2) * nbl)
ihi<-nbl - ilow
ilow<-ilow+1
ci.low=dif[ilow]
ci.hi=dif[ihi]
pv=mean(dif<0,na.rm=TRUE)
pc[j,k]=2*min(pv,1-pv)
output[ic,3]=pc[j,k]
output[ic,6]=ci.low
output[ic,7]=ci.hi
output[ic,4]=est$mean.str1
output[ic,5]=est$mean.str2
}}
if(plotit)lplot(x,y,ticktype=ticktype,xlab=xlab,ylab=ylab,zlab=zlab)
list(n=n.orig,n.keep=n,p.crit=p.crit,output=output)
}

lplotCIMC<-function(data,x,y,pts1,pts2,npts,tr,span){
temp=lplotcom2(x[data,],y[data],pts1=pts1,pts2=pts2,npts=npts,tr=tr,span=span)
v=c(temp$mean.str1,temp$mean.str2)
}

lplotCIMCv2<-function(data,x,y,pts1,pts2,npts,tr,span){
temp=lplotcom2v2(x[data,],y[data],pts1=pts1,pts2=pts2,npts=npts,tr=tr,span=span)
v=c(temp$mean.str1,temp$mean.str2)
}

lplotbsepvv3<-function(n,nreps=2000,alpha=0.05,pts=NULL,npts=25){
#
# Determine critical p-value for lplotCI.
#
set.seed(2)
pv=NA
for(i in 1:nreps){
x=rnorm(n)
y=rnorm(n)
xord=order(x)
y=y[xord]
x=x[xord]
M=median(x)
low=M-1.5*mad(x)
up=M+1.5*mad(x)
if(is.null(pts))pts=seq(low,up,length.out=npts)
plx<-predict(loess(y ~ x), se=TRUE)
est=lplot.pred(x,y,pts=pts)$yhat
se=lplotse(x,y,SEED=FALSE,pts=pts)
test=abs(est/se)
pall=2*(1-pt(abs(test),plx$df))
pv[i]=min(elimna(pall))
}
hd(pv,alpha)
}

lplotN<-function(x,y,nsub=5000,est=tmean,fr=1,xout=FALSE,xlab='X',ylab='Y',zlab='',ticktype = 'simple',theta = 50, phi = 25, scale = TRUE,pc=' ',
    expand = 0.5, SEED = TRUE,frame=TRUE){
#
 # Running interval smoother, good for large sample sizes or plots of the
 # regression surface without a scatter plot.
 #
 # nsub is size of the random sample of the data used to predict outcome using all of the data
 #  
 if(SEED)set.seed(2)
 x=as.matrix(x)
p=ncol(x)
p1=p+1
 xy=cbind(x,y)
 xy=elimna(xy)
  if(xout){
 flag<-outpro.depth(x,plotit=FALSE)$keep                                
x<-x[flag,]                                                                        
x<-as.matrix(x)                                                                    
y<-y[flag]    
xy=cbind(x,y)                                                                     
} 
 n=nrow(xy)
 nsub=min(n,nsub)
 id=sample(n,nsub) 
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
w=lplot.pred(x,y,pts=x[id,],fr=fr)$yhat 
a=lplot(x[id,],w,xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,pc=pc,
frame=frame,phi=phi,theta=theta,scale=scale,pr=FALSE)   
}


# ============================================================================
# GAM PLOTS
# ============================================================================

gamplot<-function(x,y,sop=TRUE,pyhat=FALSE,eout=FALSE,xout=FALSE,outfun=out,plotit=TRUE,
xlab="X",ylab="",zlab="",theta=50,phi=25,expand=.5,scale=TRUE,ticktype="simple"){
#
# Plot regression surface using generalized additive model
#
# sop=F, use usual linear model y~x1+x2...
# sop=T, use splines
#
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
# ANCOVA PLOTS
# ============================================================================

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

linplot<-function(x,con=0,plotfun=akerd,nboot=800,plotit=TRUE,pyhat=FALSE,...){
#
#  plot distribtion of the linear contrast
#  c_1X_2+c_2X_2+...
#
#  con contains contrast coefficients. If not specified,
#  con<-c(1,1,...,1)
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
J<-length(x)
tempn<-0
mvec<-NA
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
}
Jm<-J-1
#
# Determine contrast matrix
# If not specified, assume distribution of the sum is to be plotted
#
if(sum(con^2)==0)con<-matrix(1,J,1)
bvec<-matrix(NA,nrow=J,ncol=nboot)
for(j in 1:J){
data<-matrix(sample(x[[j]],size=nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-data
}
bcon<-t(con)%*%bvec #ncon by nboot matrix
bcon<-as.vector(bcon)
dval<-plotfun(bcon,pyhat=pyhat,...)
dval
}

lin2plot<-function(x,con,op=4,nboot=800,plotit=TRUE,pyhat=FALSE){
#
#  plot two distribtions.
#   The first is the distribtion  of the linear contrast
#  c_1X_2+c_2X_2+... c_i>0
#  and the second is the distribution of c_1X_2+c_2X_2+... c_i<0
#
#  con contains contrast coefficients. If not specified,
#  function terminates.
#
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
J<-length(x)
if(J != length(con)){
stop("Number of contrast coefficients must equal the number of groups")
}
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
#
# Determine contrast matrix for positive contrast coefficients
#
flag<-(con<0)
con1<-con
con1[flag]<-0
# Determine contrast matrix for negative contrast coefficients
flag<-(con>0)
con2<-con
con2[flag]<-0
bvec<-matrix(NA,nrow=J,ncol=nboot)
for(j in 1:J){
data<-matrix(sample(x[[j]],size=nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-data
}
bcon1<-t(con1)%*%bvec
bcon2<-t(con2)%*%bvec
bcon1<-as.vector(bcon1)
bcon2<-as.vector(bcon2)
fval<-g2plot(bcon1,bcon2,op=op,rval=15,fr=0.8,aval=0.5,xlab="X",ylab="")
fval
}


# ============================================================================
# GROUP COMPARISON PLOTS
# ============================================================================

g2plot<-function(x1,x2,op=4,rval=15,fr=.8,aval=.5,xlab="X",ylab=""){
#
# plot estimates of the density functions for two groups.
#
# op=1: Use Rosenblatt shifted histogram
#
# op=2:
# Use kernel density estimate
# Using the built-in S+ function density,
#
# op=3: Use expected frequency curve.
#
# op=4: Use adaptive kernel estimator
#
x1<-elimna(x1)
x2<-elimna(x2)
if(op==3){
rd2plot(x1,x2,fr=fr,xlab=xlab,ylab=ylab)
print("Might consider using op=4 if graph is ragged")
}
if(op==2){
tempx<-density(x1,na.rm=TRUE,kernel="epanechnikov")
tempy<-density(x2,na.rm=TRUE,kernel="epanechnikov")
plot(c(tempx$x,tempy$x),c(tempx$y,tempy$y),type="n",xlab=xlab,ylab=ylab)
lines(tempx$x,tempx$y)
lines(tempy$x,tempy$y,lty=2)
}
if(op==1){
        y1 <- sort(x1)
        z1 <- 1
        z2 <- 1
        par(yaxt = "n")
        temp <- floor(0.01 * length(x1))
        if(temp == 0)
                temp <- 5
        ibot <- y1[temp]
        itop <- y1[floor(0.99 * length(x1))]
        xaxis1 <- seq(ibot, itop, length = rval)
        for(i in 1:rval)
                z1[i] <- kerden(x1, 0, xaxis1[i])
        y2 <- sort(x2)
         temp <- floor(0.01 * length(x2))
        if(temp == 0)
                temp <- 5
        ibot <- y2[temp]
        itop <- y2[floor(0.99 * length(x2))]
        xaxis2 <- seq(ibot, itop, length = rval)
        for(i in 1:rval)
                z2[i] <- kerden(x2, 0, xaxis2[i])
plot(c(xaxis1,xaxis2),c(z1,z2), xlab =xlab, ylab =ylab, type = "n")
lines(xaxis1,z1)
lines(xaxis2,z2,lty=2)
}
if(op==4){
x1<-sort(x1)
x2<-sort(x2)
z1<-akerd(x1,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
z2<-akerd(x2,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
plot(c(x1,x2),c(z1,z2), xlab =xlab, ylab =ylab, type = "n")
lines(x1,z1)
lines(x2,z2,lty=2)
}
}

g2plotdifxy<-function(x,y,xlab="Difference",ylab=""){
#
# Plot an estimate of the distribution of X-Y
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
m<-as.vector(outer(x,y,FUN="-"))
akerd(m,xlab=xlab,ylab=ylab)
}

g5plot<-function(x1,x2,x3=NULL,x4=NULL,x5=NULL,fr=.8,aval=.5,xlab='X',ylab='',color=rep('black',5),main=NULL,sub=NULL){
#
# plot estimates of the density functions for up to 5 groups.
# using an adaptive kernel density estimator
#
if(is.matrix(x1)||is.data.frame(x1))x1=listm(x1)
if(is.list(x1)){
x=x1
J=length(x)
ic=0
for(j in 1:J){
ic=ic+1
if(ic==1)x1=x[[1]]
if(ic==2)x2=x[[2]]
if(ic==3)x3=x[[3]]
if(ic==4)x4=x[[4]]
if(ic==5)x5=x[[5]]
}
}
x1<-elimna(x1)
x2<-elimna(x2)
x1<-sort(x1)
x2<-sort(x2)
if(!is.null(x3))x3<-sort(x3)
if(!is.null(x4))x4<-sort(x4)
if(!is.null(x5))x5<-sort(x5)
z3=NULL
z4=NULL
z5=NULL
z1<-akerd(x1,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
z2<-akerd(x2,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
if(!is.null(x3))z3=akerd(x3,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
if(!is.null(x4))z4=akerd(x4,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
if(!is.null(x5))z5=akerd(x5,aval=aval,fr=fr,pyhat=TRUE,plotit=FALSE)
plot(c(x1,x2,x3,x4,x5),c(z1,z2,z3,z4,z5), xlab =xlab, ylab =ylab, type = 'n',main=main,sub=sub)
lines(x1,z1,col=color[1])
lines(x2,z2,lty=2,col=color[2])
if(!is.null(x3))lines(x3,z3,lty=3,col=color[3])
if(!is.null(x4))lines(x4,z4,lty=4,col=color[4])
if(!is.null(x5))lines(x5,z5,lty=5,col=color[5])
}

g5.cen.plot<-function(x1, x2, x3 = NULL, x4 = NULL, x5 = NULL, fr = 0.8,
    aval = 0.5, xlab = 'X', ylab ='', color = rep('black', 5),
    main = NULL, sub = NULL,loc.fun=median){
#
#  Same a g5plot, only center the data based on the
#  measure of location indicated by the argument
#  loc.fun
#
x1=elimna(x1)
x2=elimna(x2)
x1=x1-loc.fun(x1)
x2=x2-loc.fun(x2)
if(!is.null(x3))x3=x3-loc.fun(x3)
if(!is.null(x4))x4=x4-loc.fun(x4)
if(!is.null(x5))x5=x5-loc.fun(x5)
g5plot(x1=x1, x2=x2, x3=x3, x4 = x4, x5 = x5, fr = fr,
    aval = aval, xlab = xlab, ylab =ylab, color = color,
    main = main, sub = sub)
}

gplot<-function(x,xlab="Group",ylab="",xnum=FALSE){
if(is.matrix(x))x<-listm(x)
if(!xnum)par(xaxt="n")
mval<-NA
vals<-x[[1]]
gval<-rep(1,length(x[[1]]))
for(j in 2:length(x)){
vals<-c(vals,x[[j]])
gval<-c(gval,rep(j,length(x[[j]])))
}
plot(gval,vals,xlab=xlab,ylab=ylab)
}

l2plot<-function(x1,y1,x2,y2,f=2/3,SCAT=TRUE,xlab="x",ylab="y",pch='*',
eout=FALSE,xout=FALSE,...){
#
# Plot LOESS smoother for two groups
#
# f is the span used by loess
# SCAT=F, scatterplot not created, just the regression lines
# Missing values are automatically removed.
#
m<-elimna(cbind(x1,y1))
x1<-m[,1]
y1<-m[,2]
m<-elimna(cbind(x2,y2))
x2<-m[,1]
y2<-m[,2]
plot(c(x1,x2),c(y1,y2),xlab=xlab,ylab=ylab,pch=pch)
lines(lowess(x1,y1,f=f))
lines(lowess(x2,y2,f=f))
}

loc2plot<-function(x,y,plotfun=akerd,xlab='X',ylab='',...){
#
# Plot an estimate of the distribution of X-Y
# By default,
# plotfun=akerd, meaning that a kernel adaptive estimator is used.
# Other options are:
#  skerd
# kdplot
# rdplot
#
#  See Wilcox Introduction to Robust Estimation and Hypothesis Testing
#  section 3.2 for details.
#
m=elimna(cbind(x,y))
x=m[,1]
y=m[,2]
temp=temp=as.vector(outer(x,y,FUN='-'))
plotfun(temp,xlab=xlab,ylab=ylab,...)
}

sumplot2g<-function(x,y=NULL,xlab="X",ylab="",eblabx="Groups",eblaby="",nse=2){
#
# create four plots useful when comparing two groups
# 1. error bars
# 2. boxplots
# 3. kernel density estimates
# 4 shift function
#
if(!is.null(y)){
xy=list()
xy[[1]]=x
xy[[2]]=y
}
if(is.null(y)){
if(is.matrix(x))xy=matl(x)
}
par(mfrow=c(2,2))
par(oma=c(4,0,0,0))
ebarplot(xy,xlab=eblabx,ylab=eblaby,nse=nse)
boxplot(xy)
g2plot(xy[[1]],xy[[2]])
sband(xy[[1]],xy[[2]])
par(mfrow=c(1,1))
}

difQplot<-function(x,y=NULL,xlab="Quantile",ylab="Effect Size"){
#
#  Plot that provides perspective on the degree a distribution is symmetric about zero.
#  This function plots the sum of q and 1-q quantiles. If the distributions are symmetric
#  the plot should be approximately a horizontal line. If in addition the median
# of the difference scores is zero, the horizontal line will intercept the y-axis at zero.
#
if(is.null(y))dif=x
if(!is.null(y))dif=x-y
x=elimna(x)
qd=NA
for(i in 1:99)qd[i]=hd(dif,.5-i/200)+hd(dif,.5+i/200)
plot(.5-c(1:99)/200,qd,xlab=xlab,ylab=ylab)
}

bplot<-function(x,y,q=c(.25,.5,.75),vals=NULL,FUN=lincon,...){
#
#
#  x is a vector
#
#  Split the data in y into groups based on values in x.
#  vals=NULL means quantiles of x will be used, quantiles indicated by argument
#  q
#
#  Next, compare and plot  the
# groups based using the method indicated by the argument
# FUN
#
if(!is.vector(x))stop('Argument x should be a vector')
v=split.mat(x,y,q=q,vals=vals)
a=FUN(v,...)
a
}

 outreg<-function(x,y,regfun=tsreg,outfun=outpro.depth,varfun=pbvar,corfun=pbcor,xout=FALSE,...){
#
# Do regression on points not labled outliers
# based on the  method indicated by outfun
#
#   A more general version of opreg
#
#  Possible alternative choices for outfun:
#  outpro
#  outmgv
#  outmcd
#
#
#  Note: argument xout is not relevant here, but is included to avoid conflicts when using regci.
#
x<-as.matrix(x)
m<-cbind(x,y)
m<-elimna(m) # eliminate any rows with missing data
ivec<-outfun(m,plotit=FALSE)$keep
np1<-ncol(x)+1
coef<-regfun(m[ivec,1:ncol(x)],m[ivec,np1],...)$coef
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
# ERROR BAR/BOX PLOTS
# ============================================================================

ebarplot<-function(x,y=NULL,nse=2, liw = uiw, aui=NULL, ali=aui,
err="y", tr=0,ylim=NULL, sfrac = 0.01, gap=0, add=FALSE,
col=par("col"), lwd=par("lwd"), slty=par("lty"), xlab="Group",
                    ylab=NULL, ...) {
# plots error bars using the data in
# x, which is assumed to be a matrix with J columns (J groups) or
# x has list mode.
# nse indicates how many standard errors to use when plotting.
#
# By default, means are used. To use a trimmed mean, set
# tr to some value between 0 and .5
# So tr=.2 would use a 20% trimmed mean
#
# Missing values are automatically removed.
#
if(tr==.5)stop("For medians, use ebarplot.med")
if(!is.null(y)){
if(is.matrix(x))stop("When y is given, x should not be a matrix")
if(is.list(x))stop("When y is given, x should not be in list mode")
rem=x
x=list()
x[[1]]=rem
x[[2]]=y
}
if(is.matrix(x))x<-listm(x)
mval<-NA
if(!is.list(x) && is.null(y))stop("This function assumes there
 are  two or more groups")
for(j in 1:length(x))mval[j]<-mean(x[[j]],na.rm=TRUE,tr=tr)
se<-NA
#for(j in 1:length(x))se[j]<-sqrt(var(x[[j]],na.rm=TRUE)/length(x[[j]])
for(j in 1:length(x))se[j]<-trimse(x[[j]],na.rm=TRUE,tr=tr)
uiw<-nse*se
plotCI(mval,y=NULL, uiw=uiw, liw = uiw, aui=NULL, ali=aui,
                    err="y", ylim=NULL, sfrac = 0.01, gap=0, add=FALSE,
                    col=par("col"), lwd=par("lwd"), slty=par("lty"), xlab=xlab,
                    ylab=ylab)
}

ebarplot.med<-function(x,y=NULL,alpha=.05,nse=2, liw = uiw, aui=NULL, ali=aui,
err="y", tr=0,ylim=NULL, sfrac = 0.01, gap=0, add=FALSE,
col=par("col"), lwd=par("lwd"), slty=par("lty"), xlab="Group",
                    ylab=NULL, ...) {
# plots error bars using the data in
# x, which is assumed to be a matrix with J columns (J groups) or
# x has list mode.
# nse indicates how many standard errors to use when plotting.
#
# Designed specifically for medians
# Uses distribution-free confidence intervals
#
# Missing values are automatically removed.
#
if(!is.null(y)){
if(is.matrix(x))stop("When y is given, x should not be a matrix")
if(is.list(x))stop("When y is given, x should not be in list mode")
rem=x
x=list()
x[[1]]=rem
x[[2]]=y
}
if(is.matrix(x))x<-listm(x)
mval<-NA
if(!is.list(x) && is.null(y))stop("This function assumes there
 are  two or more groups")
aui=NA
ali=NA
for(j in 1:length(x)){
mval[j]<-median(x[[j]],na.rm=TRUE)
temp=sint(x[[j]],alpha=alpha,pr=FALSE)
ali[j]=temp[1]
aui[j]=temp[2]
}

plotCI(mval,y=NULL, liw = uiw, aui=aui, ali=ali,
                    err="y", ylim=NULL, sfrac = 0.01, gap=0, add=FALSE,
                    col=par("col"), lwd=par("lwd"), slty=par("lty"), xlab=xlab,
                    ylab=ylab)
}

box_plot1<-function(x,fileout){
library("reshape")
library("tidyverse")
library("viridis")
library("ggplot2")
library("hrbrthemes")

x1<-melt(x)
#tiff(fileout)
#jpeg
plot<-ggplot(x1,aes(x=variable,y=value,fill=variable))+
 geom_boxplot() +
scale_fill_viridis(discrete = TRUE, alpha=0.6) +
geom_jitter(color="black", size=1.5, alpha=0.9) +
theme_ipsum() +
 theme(
   legend.position="none",
   plot.title = element_text(size=11)
  ) +
  ggtitle("             Prediction Errors (D)") +
 xlab("")+
 ylab("")
ggsave(fileout,plot=plot)

}

STRIPchart<-function(x,method ='overplot', jitter = 0.1, offset = 1/3,
           vertical = FALSE, group.names, add = FALSE,
           at = NULL, xlim = NULL, ylim = NULL,
           ylab = NULL, xlab = NULL, dlab ='', glab ='',
           log = '', pch = 0, col = par('fg'), cex = par('cex'),
           axes = TRUE, frame.plot = axes, ...){
#
#    Same as stripchart,	only it	accepts	a matrix, unlike stripchart, which
#    allows x to be  a data frame or list mode,	but not	a matrix.
#
if(is.matrix(x))x=listm(x)
stripchart(x,method=method,jitter=jitter,offset = offset,
           vertical = vertical, group.names=group.names, add = add,
           at =at, xlim = xlim, ylim = ylim,
           ylab = ylab, xlab = xlab, dlab = dlab, glab = glab,
           log = log, pch = pch, col = col, cex = cex,
           axes = axes, frame.plot = frame.plot, ...)
}


# ============================================================================
# FUNCTIONAL DATA PLOTS
# ============================================================================

fbplot<-function(fit,x=NULL,method='MBD',depth=NULL,plot=TRUE,prob=0.5,color=6,outliercol=2,barcol=4,fullout=FALSE, factor=1.5,xlab='Time',ylab='Y',...){

  if(is.fdSmooth(fit) | is.fdPar(fit)){ fit = fit$fd }
	if(is.fd(fit)){
    if(length(x)==0){
      x = seq(fit$basis$rangeval[1],fit$basis$rangeval[2],len=101)
    }
    fit = eval.fd(x,fit)
  }

	tp=dim(fit)[1]
	n=dim(fit)[2]
	if (length(x)==0) {x=1:tp}
  #compute band depth
  if (length(depth)==0){
	if (method=='BD2') {depth=BD2(t(fit))}
	else if (method=='BD3') {depth=BD3(t(fit))}
	else if (method=='MBD') {depth=MBD(t(fit))}
	else if (method=='Both') {depth=round(BD2(t(fit)),4)*10000+MBD(t(fit))}
  }

	dp_s=sort(depth,decreasing=TRUE)
	index=order(depth,decreasing=TRUE)
	if (plot) {
	plot(x,fit[,index[1]],lty=1,lwd=2,col=1,type='l',xlab=xlab,ylab=ylab,...)
	}
	for (pp in 1:length(prob)){
		m=ceiling(n*prob[pp])#at least 50%
		center=fit[,index[1:m]]
		out=fit[,index[(m+1):n]]
		inf=apply(center,1,min)
		sup=apply(center,1,max)

		if (prob[pp]==0.5){ #check outliers
			dist=factor*(sup-inf)
			upper=sup+dist
			lower=inf-dist
			outly=(fit<=lower)+(fit>=upper)
			outcol=colSums(outly)
			remove=(outcol>0)
			#outlier column
			colum=1:n
			outpoint=colum[remove==1]
			out=fit[,remove]
			woout=fit
			good=woout[,(remove==0),drop=FALSE]
			maxcurve=apply(good,1,max)
			mincurve=apply(good,1,min)
			if (sum(outly)>0){
				if (plot) {
				matlines(x,out,lty=2,col=outliercol,type='l',...)
				}
			}
			barval=(x[1]+x[tp])/2
			bar=which(sort(c(x,barval))==barval)[1]
			if (plot) {
			lines(c(x[bar],x[bar]),c(maxcurve[bar],sup[bar]),col=barcol,lwd=2)
		    lines(c(x[bar],x[bar]),c(mincurve[bar],inf[bar]),col=barcol,lwd=2)
			}
		}
		xx=c(x,x[order(x,decreasing=TRUE)])
		supinv=sup[order(x,decreasing=TRUE)]
		yy=c(inf,supinv)
		if (plot) {
		if (prob[pp]==0.5) {polygon(xx,yy,col=color[pp],border=barcol,lwd=2)}
		else {polygon(xx,yy,col=color[pp],border=NA)}
		}
	}
	if (plot) {
	lines(x,fit[,index[1]],lty=1,lwd=2,col=1,type='l')
	lines(x,maxcurve,col=barcol,lwd=2)
	lines(x,mincurve,col=barcol,lwd=2)
	if (fullout) {
		if (sum(outly)>0){
				if (plot) {
				matlines(x,out,lty=2,col=outliercol,type='l',...)
				}
			}
		}
	}
	return(list(depth=depth,outpoint=outpoint))
}

Flplot<-function(x,est=mean,xlab='Time',ylab='Y',plotit=TRUE){
#
#  average n curves and plot results
#
es=apply(x,2,est)
if(plotit){
plot(es,xlab=xlab,ylab=ylab,type='n')
lines(es)
}
es
}

FQplot<-function(x,xlab='Time',ylab='Y',plotit=TRUE){
#
# Compute the  median and quartiles of  n curves and plot results
#
es=apply(x,2,hd)
es1=apply(x,2,hd,q=.25)
es2=apply(x,2,hd,q=.75)
if(plotit){
plot(rep(c(1:ncol(x)),3),c(es,es1,es2),xlab=xlab,ylab=ylab,type='n')
lines(es)
lines(es1,lty=2)
lines(es2,lty=2)
}
es
}

Flplot2g<-function(x1,x2,est=mean,xlab='Time',ylab='Y',plotit=TRUE){
#
#  average n curves and plot results
#
x1=elimna(x1)
x2=elimna(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 should have the same number of columns')
x1=elimna(x1)
x2=elimna(x2)
es1=apply(x1,2,est)
es2=apply(x2,2,est)
if(plotit){
plot(rep(1:ncol(x1),2),c(es1,es2),xlab=xlab,ylab=ylab,type='n')
lines(es1)
lines(es2,lty=2)
}
list(est.1=es1,est.2=es2)
}

func.plot<-function(fit, x = NULL, method ='MBD', depth = NULL, plotit = TRUE,
    prob = 0.5, color = 6, outliercol = 2, barcol = 4, fullout = FALSE,
    factor = 1.5, xlim = c(1, nrow(fit)), ylim = c(min(fit) -
        0.5 * diff(range(fit)), max(fit) + 0.5 * diff(range(fit))),xlab='Time',ylab='Y',
    ...){
#
# functional boxplot for functional data using
# method in Sun and Genton.
#
#
# fit is assumed to be an n-by-p matrix
# n= number of subjects
# p= number points where the function has been evaluated.
#
#  rows with missing values are automatically removed.
#
library(fda)
elimna(fit)
fit=t(fit)
res=fbplot(fit, x = NULL, method = method, depth = depth, plot = plotit,
    prob =prob, color =color, outliercol =outliercol, barcol = barcol,
fullout = fullout, factor = factor, xlim =xlim, ylim = ylim, xlab=xlab,ylab=ylab,...)
res
}

spag.plot<-function(x, regfun=tsreg,type = c('l',
    'p', 'b', 'o', 'c'), legend = FALSE, trace.label = deparse(substitute(trace.factor)),
    fixed = FALSE, xlab = 'Time', ylab ='',
    xtick = FALSE, xaxt = par('xaxt'), axes = TRUE, fit.lin=FALSE,...){
#
# Create a spaghetti plot for data stored in a matrix with
# n rows and p columns. The p columns
# contain  measures taken at p times for each subject.
# This function converts x into a form that can be used by interaction.plot
#
#  fit.line=TRUE means that a linear fit is plotted.
#
#  regfun: The linear fit is  based on the regression estimator indicated by
#          regfun. The  default is Theil--Sen estimator
#
#
# type: the type of plot (see plot.default): lines or points or both.
#
x=as.matrix(x)
n=nrow(x)
p=ncol(x)
np=n*p
m=matrix(NA,nrow=np,3)
pvec=c(1:p)
ic=1-p
iu=0
for(i in 1:n){
ic=ic+p
iu=iu+p
m[ic:iu,1]=i  # create Subject id.
m[ic:iu,2]=pvec
m[ic:iu,3]=x[i,]
}
if(!fit.lin)interaction.plot(m[,2],m[,1],m[,3],xlab=xlab,ylab=ylab,legend=legend,
xtick=xtick,xaxt=xaxt,axes=axes)
if(fit.lin){
fit=by(m[,2:3],m[,1],regYval,regfun=regfun)
fit1 <- unlist(fit)
names(fit1) <- NULL
#plotting the linear fit by id
interaction.plot(m[,2],m[,1], fit1,
                  xlab=xlab, ylab=ylab, legend=legend)
}
}


# ============================================================================
# INTERACTION PLOTS
# ============================================================================

interplot<-function(J,K,x,locfun=mean,locvec=NULL,na.rm=TRUE,
g1lev=NULL,g2lev=NULL,type = c("l",
    "p", "b"), xlab = "Fac 1", ylab = "means",trace.label="Fac 2",...){
if(is.null(locvec))locvec=lloc(x,est=locfun,na.rm=na.rm)
if(is.list(locvec))locvec=as.vector(matl(locvec))
if(is.null(g1lev[1])){
g1=c(rep(1,K))
for(j in 2:J)g1=c(g1,rep(j,K))
}
if(!is.null(g1lev)){
g1=c(rep(g1lev[1],K))
for(j in 2:J)g1=c(g1,rep(g1lev[j],K))
}
g1=as.factor(g1)
if(is.null(g2lev[1]))g2=as.factor(rep(c(1:K),J))
if(!is.null(g2lev[1]))g2=as.factor(rep(g2lev,J))
g2=as.factor(g2)
interaction.plot(g1,g2,locvec, xlab = xlab, ylab = ylab,
trace.label=trace.label)
}

Qinterplot<-function(x,q=.5){
#
# Plot interactions based on quantiles estimated via the
#  Harrell--Davis estimator
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
if(length(x)!=4)stop('Should have a 2 by 2 design for a total of four groups')
qv=lapply(x,hd,q=q)
qv=as.vector(matl(qv))
interplot(2,2,locvec=qv,xlab='Fac 1',ylab=paste(q,'Quantile'),trace.label='Fac 2')
}

plot.inter<-function(x,nreps=500,SEED=TRUE,xlab='DV',ylab=''){
#
#  For a 2-by-2 design, plot the distribution of X_1-X_2 and X_3-X_4
#
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
if(J!=4)stop('Length should have four groups')
x=elimna(x)
L1=NA
L2=NA
linv=NA
for(i in 1:nreps){
for(j in 1:J)linv[j]=sample(x[[j]],1)
L1[i]=linv[1]-linv[2]
L2[i]=linv[3]-linv[4]
}
g2plot(L1,L2,xlab=xlab,ylab=ylab)
}

reg.plot.inter<-function(x,y, regfun=tsreg,
 pyhat = FALSE, eout = FALSE, xout = FALSE, outfun = out,
    plotit = TRUE, expand = 0.5, scale = TRUE, xlab = "X",
    ylab = "Y", zlab = "", theta = 50, phi = 25, family = "gaussian",
    duplicate = "error",ticktype="simple",...){
#
# Plot regression surface based on the classic interaction model:
#  usual product term
#
#   x is assumed to be a matrix with two columns (two predictors)
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
if(xout){
p=ncol(x)
p1=p+1
m<-cbind(x,y)
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}

if(!scale)print("scale=F. If there is an association, try scale=T")
if(ncol(x)!=2)stop("x should have two columns")
xx=cbind(x,x[,1]*x[,2])
temp=regfun(xx,y)
fitr=y-temp$residuals
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

ols.plot.inter<-function(x,y, pyhat = FALSE, eout = FALSE, xout = FALSE, outfun = out,
    plotit = TRUE, expand = 0.5, scale = TRUE, xlab = "X",
    ylab = "Y", zlab = "", theta = 50, phi = 25, family = "gaussian",
    duplicate = "error",ticktype="simple",...){
#
# Plot regression surface based on the classic interaction model:
#  usual product term
#
#   x is assumed to be a matrix with two columns (two predictors)
x<-as.matrix(x)
xx<-cbind(x,y)
xx<-elimna(xx)
x<-xx[,1:ncol(x)]
x<-as.matrix(x)
y<-xx[,ncol(x)+1]
if(ncol(x)!=2)stop("x should have two columns")
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:2]
y<-m[,3]
}
xx=cbind(x,x[,1]*x[,2])
temp=lsfit(xx,y)
fitr=y-temp$residuals
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


# ============================================================================
# LOGISTIC/LONGITUDINAL PLOTS
# ============================================================================

logreg.plot<-function(x,y,MLE=TRUE,ROB=FALSE,xlab=NULL,ylab=NULL,zlab='P(Z=1)',xout=FALSE,outfun=outpro,
theta=50,phi=25,duplicate="error",LP=TRUE,Lspan=.75,pyhat=FALSE,LABELS=FALSE,
WARN=FALSE,BY=TRUE,
expand=.5,scale=TRUE,fr=2,ticktype="simple",pr=TRUE,...){
#
# For one predictor, plot logistic regression line
#
#  if x is a matrix with more than one column, plot is  based on data in
#  in column 1.
#
#  MLE=T, will plot usual maximum likelihood estimate using a solid line
#  ROB=T, will plot robust estimate, which is indicated by a
#  dashed line.
#
library(robustbase)
xy=cbind(x,y)
xy=elimna(xy)
p1=ncol(xy)
if(p1>3)stop('Only one or two independent variables can be used')
if(!xout){
if(pr)print('Suggest also looking at result using xout=TRUE')
}
p=p1-1
x=xy[,1:p]
x=as.matrix(x)
y=xy[,p1]
if(xout){
m<-cbind(x,y)
flag<-outfun(x,plotit=FALSE,...)$keep
m<-m[flag,]
x<-m[,1:p]
y<-m[,p1]
}
if(p==1){
if(is.null(ylab))ylab='P(Y=1|X)'
if(is.matrix(x))x=x[,1]
xord=order(x)
xx=x[xord]
yy=y[xord]
est1=logreg(xx,yy)[1:2,1]
if(is.null(xlab))v='X'
if(is.null(ylab))ylab='P(Y=1|X)'
if(LABELS)v=labels(x)[[2]]
if(MLE){
plot(xx,yy,xlab=v[1],ylab=ylab)
phat=logreg.pred(xx,yy,xx)
lines(xx,phat)
}
if(ROB){
if(!WARN)options(warn=-1)
if(!BY)est2=wlogreg(xx,yy)$coef[1:2]
if(BY)est2=BYlogreg(xx,yy)$coef[1:2]
phat2=exp(est2[1]+est2[2]*xx)/(1+exp(est2[1]+est2[2]*xx))
lines(xx,phat2,lty=2)
phat=cbind(xx,phat2)
dimnames(phat)=list(NULL,c(v,'Y.hat'))
if(!WARN)options(warn=0)
}
}
if(p==2){
fitr=logreg.pred(x,y,x)
if(is.null(xlab))v='X'
if(is.null(ylab))v[2]='Y'
if(LABELS)v=labels(x)[[2]]
if(LP)lplot(x,fitr,xlab=v[1],ylab=v[2],zlab=xlab,z=zlab,ticktype=ticktype,theta=theta,phi=phi,pr=FALSE)
phat=cbind(x,fitr)
dimnames(phat)=list(NULL,c(v,'Y.hat'))
}
if(!pyhat)phat<-"Done"
phat
}

longreg.plot<-function(x,x.col,y.col,s.id,regfun=tsreg,scat=TRUE,xlab="X",
ylab="Y"){
#
# x is a data frame or matrix
#
# Longitudinal data: plot regression lines
#
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
if(p!=2)stop("Plot allows a single covariate only")
outmat=matrix(NA,nrow=n,ncol=p)
datx=NULL
daty=NULL
for(i in 1:n){
outmat[i,]=regfun(as.matrix(xvals[[i]]),ymat[i,])$coef
temp=as.matrix(xvals[[i]])
datx=c(datx,temp)
daty=c(daty,ymat[i,])
}
if(!scat)plot(datx,daty,type="n",xlab=xlab,ylab=ylab)
if(scat)plot(datx,daty,xlab=xlab,ylab=ylab)
for(i in 1:n)abline(outmat[i,1],outmat[i,2])
}


# ============================================================================
# DEPTH/BAGPLOT
# ============================================================================

Bagplot<-function(x,plotit=TRUE,colorbag = NULL, colorloop = NULL,
    colorchull = NULL, databag = TRUE, dataloop = TRUE, plot.fence = FALSE,type='hdepth'){
#
#  requires packages mrfDepth and ggplot2
#
#   type = measure of depth: 'hdepth' =  halfspace depth,
#   'projdepth' for projection depth and
#   'sprojdepth'  for skewness-adjusted projection depth.
#
library(mrfDepth)
library(ggplot2)
z=compBagplot(x,type=type)
bagplot(z,colorbag =colorbag, colorloop = colorloop,
    colorchull =colorchull, databag=databag, dataloop =dataloop,
plot.fence = plot.fence)
}

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
if(is.na(fr))fr<-.6
m<-covmve(x)
for(i in 1:nrow(x)){
rmd[i]<-sum(near3d(x,x[i,],fr,m))
}
rmd<-rmd/nrow(x)
if(plotit && ncol(x)==2){
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

rd2plot<-function(x,y,fr=.8,xlab="",ylab=""){
#
# Expected frequency curve
# for two groups.
#
# fr controls amount of smoothing
x<-elimna(x)
y<-elimna(y)
rmdx<-NA
rmdy<-NA
for(i in 1:length(x)){
rmdx[i]<-sum(near(x,x[i],fr))
}
for(i in 1:length(y)){
rmdy[i]<-sum(near(y,y[i],fr))
}
rmdx<-rmdx/length(x)
rmdy<-rmdy/length(y)
plot(c(x,y),c(rmdx,rmdy),type="n",ylab=ylab,xlab=xlab)
sx<-sort(x)
xorder<-order(x)
sysm<-rmdx[xorder]
lines(sx,sysm)
sy<-sort(y)
yorder<-order(y)
sysm<-rmdy[yorder]
lines(sy,sysm,lty=2)
}


# ============================================================================
# DISTRIBUTION/DENSITY PLOTS
# ============================================================================

splot<-function(x,op=TRUE,VL=FALSE,xlab="X",ylab="Rel. Freq.",frame.plot=TRUE,plotit=TRUE){
#
# Frequency plot
#
# For each unique value in x,
# the relatively frequency is determined and plotted.
#
# op=TRUE a line connecting the relative frequencies is drawn if VL=FALSE.
# VL=TRUE, a vertical line is drawn for each unique value in x;
# the height of the line indicates the relative frequency.
#
# op=FALSE. No lines are drawn
#
# The function returns the sample size as well as the frequencies
# associated with each unique value stored in x.
#
x<-x[!is.na(x)]
temp<-sort(unique(x))
freq<-NA
for(i in 1:length(temp)){
freq[i]<-sum(x==temp[i])
}
rmfreq=freq
nval=sum(freq)
freq<-freq/length(x)
tfreq<-freq
tfreq[1]<-0
tfreq[2]<-max(freq)
if(plotit){
plot(temp,tfreq,xlab=xlab,ylab=ylab,type="n",frame.plot=frame.plot)
points(temp,freq,pch="*")
if(op)
if(!VL)lines(temp,freq)
if(VL){
for(i in 1:length(temp))lines(c(temp[i],temp[i]),c(0,freq[i]))
}}
den=sum(rmfreq)
list(obs.values=temp,n=nval,frequencies=rmfreq,rel.freq=rmfreq/den)
}

splotg5<-function(x1,x2=NULL,x3=NULL,x4=NULL,x5= NULL,xlab="X",ylab="Rel. Freq."){
#
# Frequency plot for up to five variables.
#
#
freqx2=NULL
freqx3=NULL
freqx4=NULL
freqx5=NULL
x1<-x1[!is.na(x1)]
x2<-x2[!is.na(x2)]
x3<-x3[!is.na(x3)]
x4<-x4[!is.na(x4)]
x5<-x5[!is.na(x5)]

xall=c(x1,x2,x3,x4,x5)
xall=xall[!is.na(xall)]
temp=sort(unique(xall))
XL=list(x1,x2,x3,x4,x5)
NN=0
for(j in 1:5)if(!is.null(XL[[j]]))NN=NN+1
freqx1<-NA
for(i in 1:length(temp)){
freqx1[i]<-sum(x1==temp[i])
}
freqx1<-freqx1/length(x1)
if(!is.null(x2)){
freqx2<-NA
for(i in 1:length(temp)){
freqx2[i]<-sum(x2==temp[i])
}
freqx2<-freqx2/length(x2)
}
if(!is.null(x3)){
freqx3<-NA
for(i in 1:length(temp)){
freqx3[i]<-sum(x3==temp[i])
}
freqx3<-freqx3/length(x3)
}
if(!is.null(x4)){
x4<-x4[!is.na(x4)]
freqx4<-NA
for(i in 1:length(temp)){
freqx4[i]<-sum(x4==temp[i])
}
freqx4<-freqx4/length(x4)
}
if(!is.null(x5)){
x5<-x5[!is.na(x5)]
freqx5<-NA
for(i in 1:length(temp)){
freqx5[i]<-sum(x5==temp[i])
}
freqx5<-freqx5/length(x5)
}
X=rep(temp,NN)
pts=c(freqx1,freqx2,freqx3,freqx4,freqx5)
plot(X,pts,type="n",xlab=xlab,ylab=ylab)
points(X,pts)
lines(temp,freqx1)
if(NN>=2)lines(temp,freqx2,lty=2)
if(NN>=3)lines(temp,freqx3,lty=3)
if(NN>=4)lines(temp,freqx4,lty=4)
if(NN>=5)lines(temp,freqx5,lty=5)
}

kdplot<-function(x,rval=15,xlab="X",ylab="Y"){
#
#   Compute the kernel density estimator for a range of values
#   and plot results.
#
#   x contains vector of observations
#
x<-x[!is.na(x)]  #  Remove any missing values
y<-sort(x)
z<-1
temp<-floor(.01*length(x))
if(temp==0)temp<-5
ibot<-y[temp]
itop<-y[floor(.99*length(x))]
xaxis<-seq(ibot,itop,length=rval)
for(i in 1:rval)z[i]<-kerden(x,0,xaxis[i])
plot(xaxis,z,xlab=xlab,ylab=ylab)
lines(xaxis,z)
}

piplot<-function(x, y, alpha = 0.05)
{
# Makes an FY plot with prediction limits added.
	x <- as.matrix(x)
	p <- dim(x)[2] + 1
	n <- length(y)
	up <- 1:n
	low <- up
	out <- lsfit(x, y)
	tem <- ls.diag(out)
	lev <- tem$hat
	res <- out$residuals
	FIT <- y - res
	Y <- y
	corfac <- (1 + 15/n)*sqrt(n/(n - p))
	val2 <- quantile(res, c(alpha/2, 1 - alpha/2))
	#get lower and upper PI limits for each case
	for(i in 1:n) {
		val <- sqrt(1 + lev[i])
		val3 <- as.single(corfac * val2[1] * val)
		val4 <- as.single(corfac * val2[2] * val)
		up[i] <- FIT[i] + val4
		low[i] <- FIT[i] + val3
	}
	zy <- c(min(low), Y, max(up))
	zx <- c(min(FIT), FIT, max(FIT))
        #change labels so plot labels are good
        ff <- FIT
        yy <- Y
        Y <- zy
        FIT <- zx
	plot(FIT, Y, type = "n")
	points(ff, yy)
	abline(0, 1)
	points(ff, up, pch = 17)
	points(ff, low, pch = 17)
}


# ============================================================================
# 3D PLOTTING
# ============================================================================

plot3D<-function(x,y,xlab='X1',ylab='X2',zlab='Y',theta=50,phi=25,
duplicate='error',pc='*',ticktype='simple',expand=.5){
#
#  A 3D plot: supplied for convenience
#
# Example: plot a regression surface
# x and y generated from regression model with no error term.
#
x=as.matrix(x)
if(ncol(x)!=2)stop('x should have two columns only')
fitr<-interp(x[,1],x[,2],y,duplicate=duplicate)
persp(fitr,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,expand=expand,
scale=scale,ticktype=ticktype)
}


# ============================================================================
# SPECIALIZED/UTILITY PLOTS
# ============================================================================

bwiJ2plot<-function(J,K,x,fr=.8,aval=.5,xlab = 'X', ylab = '',
color = rep('black', 5),BOX=FALSE){
#
# This function is for a J by 2 between by within design
#
# Plot distribution of the difference scores for
# each of the J independpent groups
#
# x: can be a matrix, organized as expected by bwimcp
# or it can have list mode.
#
if(K!=2)stop('Should have only two dependent variables')
if(J>5)stop('Can only have five levels for the independent factor')
      if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}
dif=list()
for(j in 1:5)dif[[j]]=NULL
JK=J*K
m<-matrix(c(1:JK),J,K,byrow=TRUE)
ic<-c(-1,0)
for(j in 1:J){
ic<-ic+2
dif[[j]]=x[[ic[1]]]-x[[ic[2]]]
}
if(!BOX)g5plot(dif[[1]],dif[[2]],dif[[3]],dif[[4]],dif[[5]],fr = fr,
    aval = aval, xlab = xlab, ylab = ylab, color = color)
if(BOX)boxplot(dif)
}

dlinplot<-function(x,con,xlab='DV',ylab='',sym.test=FALSE){
#
# For dependent variables,
# determine distribution of Y_i=sum_j c_jX_j
# and then plot the distribution
#
# The function also tests the hypothesis that Y	has a median of zero.
# sym.test=TRUE: test the hypothesis that Y is symmetric.
#
#  A quantile shift measure of effect size is returned as well.
#
if(is.matrix(con)){
if(ncol(con>1))print('Warning: Argument con should be a vector. Only the first contrast coefficients are used.')
}
con=as.vector(con)
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(is.data.frame(x))x=as.matrix(x)
if(is.list(x))x=matl(x)
x=elimna(x)
n=nrow(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
if(length(con)!=ncol(x))stop('Length of con should equal number of groups')
x=elimna(x)
L=NA
linv=NA
for(i in 1:n){
L[i]=sum(con*x[i,])
}
akerd(L,xlab=xlab,ylab=ylab)
mt=sintv2(L)
sym=NULL
Q=depQS(L)
if(sym.test)sym=Dqdif(L)
list(median=mt$median,n=mt$n,ci.low=mt$ci.low,ci.up=mt$ci.up,
p.value=mt$p.value,Q.effect=Q$Q.effect,sym.test=sym)
}

 dlin.sign<-function(x,con){
 # For dependent variables,
# determine distribution of Y_i=sum_j c_jX_{ij}
# and then do a sign test
con=as.vector(con)
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(is.data.frame(x))x=as.matrix(x)
if(is.list(x))x=matl(x)
x=elimna(x)
n=nrow(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
if(length(con)!=ncol(x))stop('Length of con should equal number of groups')
x=elimna(x)
L=NA
linv=NA
for(i in 1:n){
L[i]=sum(con*x[i,])
}
a=signt(dif=L)
list(Prob_a_value_is_less_than_zerro=a$Prob_x_less_than_y,ci=a$ci,n=a$n,N=a$N,p.value=a$p.value)
 }

plot_robpca<-function(robpca.obj, classic=0, labod=3, labsd=3) {
	diagnosticplot <- !any(robpca.obj$od <= as.vector(1.E-06,mode(robpca.obj$od)))
	if(diagnosticplot == T) {
		xmax <- max(max(robpca.obj$sd), robpca.obj$cutoff$sd)
		ymax <- max(max(robpca.obj$od), robpca.obj$cutoff$od)
		plot(robpca.obj$sd, robpca.obj$od, xlab="Score distance", ylab="Orthogonal distance", xlim=c(0,xmax), ylim=c(0,ymax), type="p")
		abline(v=robpca.obj$cutoff$sd)
		abline(h=robpca.obj$cutoff$od)
		givelabel(robpca.obj, labod, labsd)
	}
	else {
		ymax <- max(max(robpca.obj$sd), robpca.obj$cutoff$sd)
		plot(robpca.obj$sd, xlab="Index", ylab="Score distance", ylim=c(0,ymax), type="p")
		abline(h=robpca.obj$cutoff$sd)
		givelabel(robpca.obj, labod=0, labsd, indexplot=1)
	}
	title("ROBPCA")
	if(classic == 1) {
		diagnosticplot <- !any(robpca.obj$classic$od <= as.vector(1.E-06,mode(robpca.obj$classic$od)))
		if(diagnosticplot == T) {
			xmax <- max(max(robpca.obj$classic$sd), robpca.obj$classic$cutoff$sd)
			ymax <- max(max(robpca.obj$classic$od), robpca.obj$classic$cutoff$od)
			plot(robpca.obj$classic$sd, robpca.obj$classic$od, xlab="Score distance", ylab="Orthogonal distance", xlim=c(0,xmax), ylim=c(0,ymax), type="p")
			abline(v=robpca.obj$classic$cutoff$sd)
			abline(h=robpca.obj$classic$cutoff$od)
			givelabel(robpca.obj$classic, labod, labsd)
		}
		else {
			ymax <- max(max(robpca.obj$classic$sd), robpca.obj$classic$cutoff$sd)
			plot(robpca.obj$classic$sd, xlab="Index", ylab="Score distance", ylim=c(0,ymax), type="p")
			abline(h=robpca.obj$classic$cutoff$sd)
			givelabel(robpca.obj$classic, labod=0, labsd, indexplot=1)
		}
		title("CPCA")
	}
	invisible(robpca.obj)
}
"givelabel"<-function(object, labod, labsd, indexplot=0) {
	if((labod == 0) && (labsd == 0)) {
		return(invisible(object))
	}
	if(indexplot != 1) {
		order.od <- order(object$od*(-1))
		order.sd <- order(object$sd*(-1))
		if(labod != 0) {
			for(i in 1:labod) {
				lab <- ifelse(is.character(names(object$od)), names(object$od[order.od[i]]), order.od[i])
				text(object$sd[order.od[i]], object$od[order.od[i]]+par("cxy")[2], labels=lab)
			}
		}
		if(labsd != 0) {
			for(i in 1:labsd) {
				lab <- ifelse(is.character(names(object$sd)), names(object$od[order.sd[i]]), order.sd[i])
				text(object$sd[order.sd[i]], object$od[order.sd[i]]+par("cxy")[2], labels=lab)
			}
		}
	}
	else {
		order.sd <- order(object$sd*(-1))
		if(labsd != 0) {
			for(i in 1:labsd) {
				lab <- ifelse(is.character(names(object$sd)), names(object$sd[order.sd[i]]), order.sd[i])
				text(order.sd[i], object$sd[order.sd[i]]+par("cxy")[2], labels=lab)
			}
		}
	}
	return(invisible(object))
}

testplot<-function(x){
plot(x[,1],x[,2])
}


plotDAP<-function(rawData,CRP,center,DataMean,name){
library("plotrix")
#rawDat<-read.table("rawData.txt",sep="\t",header=T)
#realDat<-read.table("liwData.txt",sep="\t",header=F)
colnames(rawData)<-c("X","Y")
colnames(CRP)<-c("X","Y")
colnames(center)<-c("X","Y")
colnames(DataMean)<-c("X","Y")

D1<-data.frame(X1=CRP$X-mean(CRP$X),Y1=CRP$Y-mean(CRP$Y),X=CRP$X,Y=CRP$Y)
tmp1<-D1[D1$Y1>0,]
tmp2<-D1[D1$Y1<0,]
Tmp1<-tmp1[order(tmp1$X1,decreasing=T),]
Tmp2<-tmp2[order(tmp2$X1,decreasing=F),]
Tmp<-rbind(Tmp1,Tmp2)
Tmp$a<-atan2(Tmp$Y1,Tmp$X1)
Tmp<-Tmp[order(Tmp$a,decreasing=T),]
## add an ending point as the starting point, so that the circle is complete
newDat<-rbind(Tmp,Tmp[1,])

D2<-data.frame(X1=DataMean$X-mean(DataMean$X),Y1=DataMean$Y-mean(DataMean$Y),X=DataMean$X,Y=DataMean$Y)
tmp1<-D2[D2$Y1>0,]
tmp2<-D2[D2$Y1<0,]
Tmp1<-tmp1[order(tmp1$X1,decreasing=T),]
Tmp2<-tmp2[order(tmp2$X1,decreasing=F),]
Tmp<-rbind(Tmp1,Tmp2)
Tmp$a<-atan2(Tmp$Y1,Tmp$X1)
Tmp<-Tmp[order(Tmp$a,decreasing=T),]
## add an ending point as the starting point, so that the circle is complete
newDat2<-rbind(Tmp,Tmp[1,])

R<-4
cos45<-cos(pi/4)
#tiff("DoubleAnglePlot.tiff")
plot(newDat$X,newDat$Y,type="p",xlim=c(-5,5),ylim=c(-5,5),col="blue",pch=4,cex=0,frame.plot=F,axes=FALSE,xlab="",ylab="",asp=1,main=paste(name,""))
points(rawData$X,rawData$Y,type="p",col="black",pch=19,cex=0.5)
points(newDat2$X,newDat2$Y,type="p",col="purple",pch=8,cex=0)
points(center$X,center$Y,type="p",col="red",pch=15,cex=1)
lines(newDat$X,newDat$Y,col="blue",lwd=1.3)
lines(newDat2$X,newDat2$Y,col="purple",lwd=2.0)
### lines(D[,1],D[,2],type="p",col="red",pch=0,cex=5)

draw.circle(0,0,1,lty=1,lwd=0.5)
draw.circle(0,0,2,lty=1,lwd=0.5)
draw.circle(0,0,3,lty=1,lwd=0.5)
draw.circle(0,0,4,lty=1,lwd=0.5)
segments(-1*R*cos45,-1*R*cos45,1*R*cos45,R*cos45,lty=1,lwd=0.5)
segments(-1*R*cos45,1*R*cos45,1*R*cos45,-1*R*cos45,lty=1,lwd=0.5)
segments(0,R,0,-1*R,lty=1,lwd=0.5)
segments(R,0,-1*R,0,lty=1,lwd=0.5)

R1<-R+0.6
text(R1,0,paste0("0",intToUtf8(176)))
text(R1*cos45,R1*cos45,paste0("22.5",intToUtf8(176)))
text(0, R1,paste0("45",intToUtf8(176)))
text(-1*R1*cos45,R1*cos45,paste0("67.5",intToUtf8(176)))
text(-1*R1,0,paste0("90",intToUtf8(176)))
text(-1*R1*cos45,-1*R1*cos45,paste0("112.5",intToUtf8(176)))
text(0,-1*R1,paste0("135",intToUtf8(176)))
text(R1*cos45,-1*R1*cos45,paste0("157.5",intToUtf8(176)))
#dev.off()

#legend(x="bottomleft",pch=c(15,4,8),legend=c("Centroid","Mean Convex Polygon", "Dataset Convex Polygon"),col=c("red","blue","purple"))
legend(x="bottomleft",pch=c(15,NA,NA), lty=c(NA,1,1),cex=0.88,legend=c("Centroid","Mean Convex Polygon", "Dataset Convex Polygon"),col=c("red","blue","purple"))
}

