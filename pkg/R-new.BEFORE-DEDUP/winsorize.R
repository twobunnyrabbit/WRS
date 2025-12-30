# =============================================================================
# WRS: Winsorization Methods
# =============================================================================
#
# Core winsorization functions and utilities
#
# Functions in this module:
#   - win: Core winsorization function
#   - winmean: Winsorized mean
#   - winsd, winsd05, winsdN: Winsorized standard deviations
#   - winse, winci: Standard error and confidence interval for winsorized mean
#   - winvarN: Normalized winsorized variance
#   - winsorized: General winsorization with threshold
#   - WINCOR: Winsorized correlation matrix convenience function
#
# Note: Core utilities like winvar, winval, winall are in 00-utils-core.R
#       Winsorized correlation (wincor) is in correlation.R
#       Winsorized covariance (wincov) is in covariance.R
#       Winsorized regression (winreg) is in regression.R
#
# =============================================================================

winsd<-function(x,tr=.2,na.rm=FALSE){
val=sqrt(winvar(x,tr=tr,na.rm=na.rm))
val
}

winsd05<-function(x,tr=.2,na.rm=FALSE){
val=sqrt(winvar(x,tr=tr,na.rm=na.rm))
val
}

win<-function(x,tr=.2){
#
#  Compute the gamma Winsorized mean for the data in the vector x.
#
#  tr is the amount of Winsorization
#
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
y<-ifelse(y<=xbot,xbot,y)
y<-ifelse(y>=xtop,xtop,y)
win<-mean(y)
win
}

winmean<-function(x,tr=.2,na.rm=TRUE){
if(na.rm)x=elimna(x)
winmean<-mean(winval(x,tr))
winmean
}

winvarN<-function(x,tr=.2){
#
# rescale the Winsorized variance so that it equals one for the standard
# normal distribution
#
x=elimna(x)
cterm=NULL
if(tr==0)cterm=1
if(tr==0.1)cterm=0.6786546
if(tr==0.2)cterm=0.4120867
if(is.null(cterm))cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
bot=winvar(x,tr=tr)/cterm
bot
}

winse<-function(x,tr=.2){
#
# Estimate the standard error of the Winsorized mean
#
x=elimna(x)
n=length(x)
h=n-2*floor(tr*n)
top=(n-1)*sqrt(winvar(x,tr=tr))
bot=(h-1)*sqrt(n)
se=top/bot
se
}

winci<-function(x,tr=.2,alpha=.05,null.value=0,pr=TRUE){
#
#  Compute a 1-alpha confidence interval for the Winsorized mean
#
#  The default amount of  Winsorizing is tr=.2
#
if(pr){
print("The p-value returned by the this function is based on the")
print("null value specified by the argument null.value, which defaults to 0")
}
x<-elimna(x)
se<-winse(x,tr=tr)
df<-length(x)-2*floor(tr*length(x))-1
trimci<-winmean(x,tr)-qt(1-alpha/2,df)*se
trimci[2]<-winmean(x,tr)+qt(1-alpha/2,df)*se
test<-(winmean(x,tr)-null.value)/se
sig<-2*(1-pt(abs(test),df))
list(ci=trimci,test.stat=test,p.value=sig)
}

winsorized<- function(x,a=1.5,sigma=1) {
s<-sigma
newx<-x
indp<-x>(a*s)
newx[indp]<-(a*s)
indn<- x<(a*-s)
newx[indn]<- (-a*s)
newx
}

WINCOR<-function(x,tr=.2){
#
# For convenience, compute Winsorized correlation matrix only.
#
a=winall(x,tr=tr)$cor
a
}

winsdN<-function(x,tr=.2){
#
# Rescale a Winsorized standard deviation so that it estimates
# the population standard deviation under normality.
#
x=elimna(x)
e=winsd(x,tr=tr)
if(tr==0)cterm=1
else
cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
cterm=sqrt(cterm)
e=e/cterm
e
}
