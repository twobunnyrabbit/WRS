# ==============================================================================
# Two-Sample Comparison Functions
# ==============================================================================
#
# This module contains functions for comparing two groups (independent or paired).
#
# Main functions:
#   - yuen(), yuend() - Trimmed mean comparisons (independent/dependent)
#   - pb2gen() - Percentile bootstrap for two groups
#   - wmw*() - Wilcoxon-Mann-Whitney and related tests
#   - cid*() - Cliff's analog and related methods
#   - qcomhd*() - Quantile comparisons using Harrell-Davis
#   - trimpb*() - Bootstrap trimmed mean comparisons
#   - two*() - Various two-sample comparison methods
#
# Extracted: 2025-12-30
# Functions: 103
# ==============================================================================

qcomhdMC<-function(x,y,est=hd,q=c(.1,.25,.5,.75,.9),nboot=4000,plotit=TRUE,SEED=TRUE,xlab="Group 1",ylab="Est.1-Est.2",alpha=.05,ADJ.CI=TRUE){
#
# Compare quantiles using pb2gen
# via hd estimator. Tied values are allowed.
#
#  ADJ.CI=TRUE means that the confidence intervals are adjusted based on the level used by the corresponding
#  test statistic. If a test is performed with at the .05/3 level, for example, the confidence returned has
#  1-.05/3 probability coverage.
#
# When comparing lower or upper quartiles, both power and the probability of Type I error
# compare well to other methods that have been derived.
# q: can be used to specify the quantiles to be compared
# q defaults to comparing the .1,.25,.5,.75, and .9 quantiles
#
#   Function returns p-values and critical p-values based on Hochberg's method.
#
library(parallel)
if(SEED)set.seed(2)
print('Can also use the function qcomhd with the argument MC=TRUE')
pv=NULL
output=matrix(NA,nrow=length(q),ncol=10)
dimnames(output)<-list(NULL,c("q","n1","n2","est.1","est.2","est.1_minus_est.2","ci.low","ci.up","p_crit","p-value"))
for(i in 1:length(q)){
output[i,1]=q[i]
output[i,2]=length(elimna(x))
output[i,3]=length(elimna(y))
output[i,4]=hd(x,q=q[i])
output[i,5]=hd(y,q=q[i])
output[i,6]=output[i,4]-output[i,5]
temp=qcom.sub(x,y,nboot=nboot,q=q[i],SEED=FALSE,alpha=alpha)
output[i,7]=temp$ci[1]
output[i,8]=temp$ci[2]
output[i,10]=temp$p.value
}
temp=order(output[,10],decreasing=TRUE)
zvec=alpha/c(1:length(q))
output[temp,9]=zvec
if(ADJ.CI){
for(i in 1:length(q)){
temp=pb2genMC(x,y,nboot=nboot,est=est,q=q[i],SEED=FALSE,alpha=output[i,9],pr=FALSE)
output[i,7]=temp$ci[1]
output[i,8]=temp$ci[2]
output[i,10]=temp$p.value
}
temp=order(output[,10],decreasing=TRUE)
}
output <- data.frame(output)
output$signif=rep("YES",nrow(output))
for(i in 1:nrow(output)){
if(output[temp[i],10]>output[temp[i],9])output$signif[temp[i]]="NO"
#if(output[temp[i],10]<=output[temp[i],9])break
}
if(plotit){
xax=rep(output[,4],3)
yax=c(output[,6],output[,7],output[,8])
plot(xax,yax,xlab=xlab,ylab=ylab,type="n")
points(output[,4],output[,6],pch="*")
lines(output[,4],output[,6])
points(output[,4],output[,7],pch="+")
points(output[,4],output[,8],pch="+")
}
output
}

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

twomanbt<-function(x,y,tr=.2,alpha=.05,nboot=599){
#
#   Two-sample Behrens-Fisher problem.
#
#   For each of two independent groups,
#   have p measures for each subject. The goal is to compare the
#   trimmed means of the first measure, the trimmed means for the second
#   and so on.   So there are a total of p comparisons between the two
#   groups, one for each measure.
#
#   The percentile t bootstrap method is used to
#   compute a .95 confidence interval.
#
#   By default, 20% trimming is used with B=599 bootstrap samples.
#
#   x contains the data for the first group; it
#    can be an n by J matrix or it can have list mode.
#   y contains the data for the second group.
#
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(!is.list(y) && !is.matrix(y))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
# put the data in an n by p matrix
matx<-matrix(0,length(x[[1]]),length(x))
for (j in 1:length(x))matx[,j]<-x[[j]]
}
if(is.list(y)){
# put the data in an n by p matrix
maty<-matrix(0,length(y[[1]]),length(y))
for (j in 1:length(y))maty[,j]<-y[[j]]
}
if(is.matrix(x)){
matx<-x
}
if(is.matrix(y)){
maty<-y
}
if(ncol(matx)!=ncol(maty))stop("The number of variables for group one is not equal to the number for group 2")
if(sum(is.na(mat)>=1))stop("Missing values are not allowed.")
J<-ncol(mat)
connum<-ncol(matx)
bvec<-matrix(0,connum,nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
xcen<-matrix(0,nrow(matx),ncol(matx))
ycen<-matrix(0,nrow(maty),ncol(maty))
for (j in 1:connum)xcen[,j]<-matx[,j]-mean(matx[,j],tr) #Center data
for (j in 1:connum)ycen[,j]<-maty[,j]-mean(maty[,j],tr) #Center data
print("Taking bootstrap samples. Please wait.")
bootx<-sample(nrow(matx),size=nrow(matx)*nboot,replace=TRUE)
booty<-sample(nrow(maty),size=nrow(maty)*nboot,replace=TRUE)
matval<-matrix(0,nrow=nboot,ncol=connum)
for (j in 1:connum){
datax<-matrix(xcen[bootx,j],ncol=nrow(matx))
datay<-matrix(ycen[booty,j],ncol=nrow(maty))
paste("Working on variable", j)
top<- apply(datax, 1., mean, tr) - apply(datay, 1., mean, tr)
botx <- apply(datax, 1., trimse, tr)
boty <- apply(datay, 1., trimse, tr)
matval[,j]<-abs(top)/sqrt(botx^2. + boty^2.)
}
bvec<-apply(matval,1,max)
icrit<-round((1-alpha)*nboot)
bvec<-sort(bvec)
crit<-bvec[icrit]
psihat<-matrix(0,ncol=4,nrow=connum)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol=3,nrow=connum)
dimnames(test)<-list(NULL,c("con.num","test","se"))
for(j in 1:ncol(matx)){
temp<-yuen(matx[,j],maty[,j],tr=tr)
test[j,1]<-j
test[j,2]<-abs(temp$test)
test[j,3]<-temp$se
psihat[j,1]<-j
psihat[j,2]<-mean(matx[,j],tr)-mean(maty[,j])
psihat[j,3]<-mean(matx[,j],tr)-mean(maty[,j])-crit*temp$se
psihat[j,4]<-mean(matx[,j],tr)-mean(maty[,j])+crit*temp$se
}
list(psihat=psihat,teststat=test,critical.value=crit)
}

twobici<-function(r1=sum(elimna(x)),n1=length(elimna(x)),r2=sum(elimna(y)),n2=length(elimna(y)),
x=NA,y=NA,alpha=.05){
#
# Compute confidence interval for p1-p2,
# the difference between probabilities of
# success for a two binomials using Beal's method.
#
# r is number of successes
# n is sample size
# if x contains data, r1 is taken to be the
# number of 1s in x and n1 is length(x)
#
if(length(r1)>1)stop("r1 must be a single number, not a vector")
if(length(n1)>1)stop("n1 must be a single number, not a vector")
if(length(r2)>1)stop("r2 must be a single number, not a vector")
if(!is.na(sum(r1)) || !is.na(sum(n1)) || !is.na(sum(r2)) || !is.na(sum(n2))){
if(r1<0 || n1<0)stop("Both r1 and n1 must be greater than 0")
if(r1 > n1)stop("r1 can't be greater than n1")
if(r2<0 || n2<0)stop("Both r2 and n2 must be greater than 0")
if(r2 > n2)stop("r2 can't be greater than n2")
}
if(!is.na(sum(x))){
r1<-sum(x)
n1<-length(x)
}
if(!is.na(sum(y))){
r2<-sum(y)
n2<-length(y)
}
a<-(r1/n1)+(r2/n2)
b<-(r1/n1)-(r2/n2)
u<-.25*((1/n1)+(1/n2))
v<-.25*((1/n1)-(1/n2))
V<-u*((2-a)*a-b^2)+2*v*(1-a)*b
crit<-qchisq(1-alpha/2,1)
A<-sqrt(crit*(V+crit*u^2*(2-a)*a+crit*v^2*(1-a)^2))
B<-(b+crit*v*(1-a))/(1+crit*u)
ci<-NA
ci[1]<-B-A/(1+crit*u)
ci[2]<-B+A/(1+crit*u)
p1<-r1/n1
p2<-r2/n2
list(ci=ci,p1=p1,p2=p2)
}

trimpb2<-function(x,y,tr=.2,alpha=.05,nboot=2000,WIN=FALSE,win=.1,plotit=FALSE,op=4,
SEED=TRUE){
#
#   Compute a 1-alpha confidence interval for
#   the difference between two 20% trimmed means.
#   Independent groups are assumed.
#
#   The default number of bootstrap samples is nboot=2000
#
#   tr is the amount of trimming
#
#   win is the amount of Winsorizing before bootstrapping
#   when WIN=T.
#
#   Missing values are automatically removed.
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
if(WIN){
if(win>tr)stop("Cannot Winsorize more than you trim")
if(tr < .2){print("When Winsorizing, the amount of trimming")
print("should be at least .2")
}
if(min(c(length(x),length(y))) < 15){
print ("Warning: Winsorizing with sample sizes less than 15")
print("can result in poor control over the probability of a Type I error")
}
x<-winval(x,win)
y<-winval(y,win)
}
xx<-list()
xx[[1]]<-x
xx[[2]]<-y
e1=mean(xx[[1]],tr=tr)
e2=mean(xx[[2]],tr=tr)
#est.dif<-tmean(xx[[1]],tr=tr)-tmean(xx[[2]],tr=tr)
est.dif=e1-e2
crit<-alpha/2
temp<-round(crit*nboot)
icl<-temp+1
icu<-nboot-temp
bvec<-matrix(NA,nrow=2,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:2){
data<-matrix(sample(xx[[j]],size=length(xx[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,mean,tr) # Bootstrapped trimmed means for jth group
}
top<-bvec[1,]-bvec[2,]
test<-sum(top<0)/nboot+.5*sum(top==0)/nboot
if(test > .5)test<-1-test
top<-sort(top)
ci<-NA
ci[1]<-top[icl]
ci[2]<-top[icu]
if(plotit)g2plot(bvec[1,],bvec[2,],op=op)
list(Est1=e1,Est2=e2,p.value=2*test,ci=ci,est.dif=est.dif)
}

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

wmw<-function(x,y){
#
# Do Mann-Whitney test
# Return the usual p-value followed by adjusted
# p-value using Hodges, Ramsey and Wechsler (1990) method
# (See Wilcox, 2003, p. 559.)
#
x=elimna(x)
y=elimna(y)
m<-length(x)
n<-length(y)
com<-rank(c(x,y))
xp1<-length(x)+1
x<-com[1:length(x)]
y<-com[xp1:length(com)]
u<-sum(y)-n*(n+1)/2
sigsq<-m*n*(n+m+1)/12
yv<-(u+.5-m*n/2)/sqrt(sigsq)
kv<-20*m*n*(m+n+1)/(m^2+n^2+n*m+m+n)
S<-yv^2
T1<-S-3
T2<-(155*S^2-416*S-195)/42
cv<-1+T1/kv+T2/kv^2
sighrw<-2*(1-pnorm(abs(cv*yv)))
z<-(u-(.5*m*n))/sqrt(sigsq)
sig<-2*(1-pnorm(abs(z)))
list(p.value=sig,adj.p.value=sighrw,p.hat=u/(n*m))
}

twopcor<-function(x1,y1,x2,y2,SEED=TRUE){
#
#   Compute a .95 confidence interval for
#   the difference between two Pearson
#   correlations corresponding to two independent
#   goups.
#
#   This function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
#
#   WARNING: If the number of bootstrap samples is altered, it is
#   unknown how to adjust the confidence interval when n1+n2 < 250.
#
nboot<-599  #Number of bootstrap samples
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
X<-elimna(cbind(x1,y1))
x1<-X[,1]
y1<-X[,2]
X<-elimna(cbind(x2,y2))
x2<-X[,1]
y2<-X[,2]
data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(data1,1,pcorbsub,x1,y1) # A 1 by nboot matrix.
data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
bvec2<-apply(data2,1,pcorbsub,x2,y2) # A 1 by nboot matrix.
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
r1<-cor(x1,y1)
r2<-cor(x2,y2)
ci<-c(bsort[ilow],bsort[ihi])
list(r1=r1,r2=r2,ci=ci)
}

tworhobt<-function(X1,Y1,X2,Y2,alpha=.05,nboot=499,SEED=TRUE){
#
# compare two independent correlations using a bootstrap-t method in conjunction with the HC4 estimator
#
if(SEED)set.seed(2)
r1=cor(X1,Y1)
r2=cor(X2,Y2)
n1=length(X1)
n2=length(X2)
v=NA
Nboot=nboot+1
for(i in 1:Nboot){
if(i<=nboot){
id1=sample(n1,n1,replace=TRUE)
id2=sample(n2,n2,replace=TRUE)
}
if(i==Nboot){
id1=c(1:n1)
id2=c(1:n2)
}
x1=X1[id1]
y1=Y1[id1]
x2=X2[id2]
y2=Y2[id2]
x1=(x1-mean(x1))/sd(x1)
y1=(y1-mean(y1))/sd(y1)
x2=(x2-mean(x2))/sd(x2)
y2=(y2-mean(y2))/sd(y2)
temp1=olshc4(x1,y1)
temp2=olshc4(x2,y2)
if(i<=nboot)v[i]=(temp1$ci[2,2]-r1-temp2$ci[2,2]+r2)/sqrt(temp1$ci[2,6]^2+temp2$ci[2,6]^2)
if(i==Nboot)v[i]=(temp1$ci[2,2]-temp2$ci[2,2])/sqrt(temp1$ci[2,6]^2+temp2$ci[2,6]^2)
}
ibot<-round(alpha*nboot/2)
itop<-nboot-ibot+1
ibot=ibot+1    #adjusted so that p-value and confidence interval give consistent results.
vs=sort(v[1:nboot])
crit=c(vs[ibot],vs[itop])
test=v[Nboot]
if(test<0)G=mean(test>v[1:nboot])
if(test>=0)G=mean(test<v[1:nboot])
pv=2*G
if(pv>1)pv=1
if(pv<0)pv=0
list(test=test,crit.val=crit,p.value=pv)
}

twobinom<-function(r1=sum(elimna(x)),n1=length(elimna(x)),r2=sum(elimna(y)),n2=length(elimna(y)),x=NA,y=NA,alpha=.05){
#
# Test the hypothesis that two independent binomials have equal
# probability of success using the Storer--Kim method.
#
# r1=number of successes in group 1
# n1=number of observations in group 1
#
n1p<-n1+1
n2p<-n2+1
n1m<-n1-1
n2m<-n2-1
chk<-abs(r1/n1-r2/n2)
x<-c(0:n1)/n1
y<-c(0:n2)/n2
phat<-(r1+r2)/(n1+n2)
m1<-outer(x,y,"-")
m2<-matrix(1,n1p,n2p)
flag<-(abs(m1)>=chk)
m3<-m2*flag
b1<-1
b2<-1
xv<-c(1:n1)
yv<-c(1:n2)
xv1<-n1-xv+1
yv1<-n2-yv+1
dis1<-c(1,pbeta(phat,xv,xv1))
dis2<-c(1,pbeta(phat,yv,yv1))
pd1<-NA
pd2<-NA
for(i in 1:n1)pd1[i]<-dis1[i]-dis1[i+1]
for(i in 1:n2)pd2[i]<-dis2[i]-dis2[i+1]
pd1[n1p]<-phat^n1
pd2[n2p]<-phat^n2
m4<-outer(pd1,pd2,"*")
test<-sum(m3*m4)
list(p.value=test,p1=r1/n1,p2=r2/n2,est.dif=r1/n1-r2/n2)
}

twodcor8<-function(x,y){
#
#   Compute a .95 confidence interval for
#   the difference between two dependent
#   correlations corresponding to two independent
#   goups.
#
#
# x is a matrix with two columns,
# y is a vector
#  Goal: test equality of Pearson correlation for x1, y versus x2, y.
#
# For general use, twodcor10 is probably better,
# which calls this function and estimates an adjusted p-value.
#
X<-elimna(cbind(x,y))
Z1<-(X[,1]-mean(X[,1]))/sqrt(var(X[,1]))
Z2<-(X[,2]-mean(X[,2]))/sqrt(var(X[,2]))
temp<-cor.test(Z1-Z2,X[,3])
temp<-temp[3]$p.value
list(p.value=temp)
}

twodcor10<-function(x,y,nboot=500,SEED=TRUE,alpha=.05){
#
#   Compute a .95 confidence interval for
#   the difference between two dependent
#   correlations corresponding to two independent
#   goups.
#
# x is a matrix with two columns,
# y is a vector
#  Goal: test equality of Pearson correlation for x1, y versus x2, y.
#
#   This function uses an adjusted p-value, the adjustment
#  being made assuming normality.
#
#  nboot indicates how many samples from a normal distribution
#  are used to approximate the adjustment.
#
# Simulations suggest that this fucntion
#  continues to work well under non-normality.
#
if(SEED)set.seed(2)
X<-elimna(cbind(x,y))
if(ncol(X)!=3)stop("x should be a matrix with  two columns")
n<-nrow(X)
cval<-cor(X)
nval<-(cval[1,3]+cval[2,3])/2
cmat<-bdiag(1,3,nval)
cmat[1,2]<-nval
cmat[2,1]<-nval
pval<-NA
for(i in 1:nboot){
d<-rmul(n,p=3,cmat=cmat)
pval[i]<-twodcor8(d[,1:2],d[,3])$p.value
}
pval<-sort(pval)
iv<-round(alpha*nboot)
est.p<-pval[iv]
adp<-alpha/est.p
test<-twodcor8(X[,1:2],X[,3])$p.value
p.value<-test*adp
if(p.value>1)p.value<-1
list(p.value=p.value)
}

twodcor10<-function(x,y,nboot=500,SEED=TRUE,alpha=.05){
#
#   Compute a .95 confidence interval for
#   the difference between two dependent
#   correlations corresponding to two independent
#   goups.
#
# x is a matrix with two columns,
# y is a vector
#  Goal: test equality of Pearson correlation for x1, y versus x2, y.
#
#   This function uses an adjusted p-value, the adjustment
#  being made assuming normality.
#
#  nboot indicates how many samples from a normal distribution
#  are used to approximate the adjustment.
#
# Simulations suggest that this fucntion
#  continues to work well under non-normality.
#
if(SEED)set.seed(2)
X<-elimna(cbind(x,y))
if(ncol(X)!=3)stop("x should be a matrix with  two columns")
n<-nrow(X)
cval<-cor(X)
nval<-(cval[1,3]+cval[2,3])/2
cmat<-bdiag(1,3,nval)
cmat[1,2]<-nval
cmat[2,1]<-nval
pval<-NA
for(i in 1:nboot){
d<-rmul(n,p=3,cmat=cmat)
pval[i]<-twodcor8(d[,1:2],d[,3])$p.value
}
pval<-sort(pval)
iv<-round(alpha*nboot)
est.p<-pval[iv]
adp<-alpha/est.p
test<-twodcor8(X[,1:2],X[,3])$p.value
p.value<-test*adp
if(p.value>1)p.value<-1
list(p.value=p.value)
}

twodcor8<-function(x,y){
#
#   Compute a .95 confidence interval for
#   the difference between two dependent
#   correlations corresponding to two independent
#   goups.
#
#
# x is a matrix with two columns,
# y is a vector
#  Goal: test equality of Pearson correlation for x1, y versus x2, y.
#
# For general use, twodcor10 is probably better,
# which calls this function and estimates an adjusted p-value.
#
X<-elimna(cbind(x,y))
Z1<-(X[,1]-mean(X[,1]))/sqrt(var(X[,1]))
Z2<-(X[,2]-mean(X[,2]))/sqrt(var(X[,2]))
temp<-cor.test(Z1-Z2,X[,3])
temp<-temp[3]$p.value
list(p.value=temp)
}

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

mcppb20<-function(x,crit=NA,con=0,tr=.2,alpha=.05,nboot=2000,grp=NA,WIN=FALSE,
win=.1){
#
#   Compute a 1-alpha confidence interval for a set of d linear contrasts
#   involving trimmed means using the percentile bootstrap method.
#   Independent groups are assumed.
#
#   The data are assumed to be stored in x in list mode.  Thus,
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J, say.
#
#   By default, all pairwise comparisons are performed, but contrasts
#   can be specified with the argument con.
#   The columns of con indicate the contrast coefficients.
#   Con should have J rows, J=number of groups.
#   For example, con[,1]=c(1,1,-1,-1,0,0) and con[,2]=c(,1,-1,0,0,1,-1)
#   will test two contrasts: (1) the sum of the first two trimmed means is
#   equal to the sum of the second two, and (2) the difference between
#   the first two is equal to the difference between the trimmed means of
#   groups 5 and 6.
#
#   The default number of bootstrap samples is nboot=2000
#
#
con<-as.matrix(con)
if(is.matrix(x)){
xx<-list()
for(i in 1:ncol(x)){
xx[[i]]<-x[,i]
}
x<-xx
}
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
if(!is.na(sum(grp))){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[1]]]
x<-xx
}
J<-length(x)
tempn<-0
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
}
Jm<-J-1
d<-ifelse(sum(con^2)==0,(J^2-J)/2,ncol(con))
if(is.na(crit) && tr != .2){
print("A critical value must be specified when")
stop("the amount of trimming differs from .2")
}
if(WIN){
if(tr < .2){
print("Warning: When Winsorizing, the amount")
print("of trimming should be at least .2")
}
if(win > tr)stop("Amount of Winsorizing must <= amount of trimming")
if(min(tempn) < 15){
print("Warning: Winsorizing with sample sizes")
print("less than 15 can result in poor control")
print("over the probability of a Type I error")
}
for (j in 1:J){
x[[j]]<-winval(x[[j]],win)
}
}
if(is.na(crit)){
if(d==1)crit<-alpha/2
if(d==2 && alpha==.05 && nboot==1000)crit<-.014
if(d==2 && alpha==.05 && nboot==2000)crit<-.014
if(d==3 && alpha==.05 && nboot==1000)crit<-.009
if(d==3 && alpha==.05 && nboot==2000)crit<-.0085
if(d==3 && alpha==.025 && nboot==1000)crit<-.004
if(d==3 && alpha==.025 && nboot==2000)crit<-.004
if(d==3 && alpha==.01 && nboot==1000)crit<-.001
if(d==3 && alpha==.01 && nboot==2000)crit<-.001
if(d==4 && alpha==.05 && nboot==2000)crit<-.007
if(d==5 && alpha==.05 && nboot==2000)crit<-.006
if(d==6 && alpha==.05 && nboot==1000)crit<-.004
if(d==6 && alpha==.05 && nboot==2000)crit<-.0045
if(d==6 && alpha==.025 && nboot==1000)crit<-.002
if(d==6 && alpha==.025 && nboot==2000)crit<-.0015
if(d==6 && alpha==.01 && nboot==2000)crit<-.0005
if(d==10 && alpha==.05 && nboot<=2000)crit<-.002
if(d==10 && alpha==.05 && nboot==3000)crit<-.0023
if(d==10 && alpha==.025 && nboot<=2000)crit<-.0005
if(d==10 && alpha==.025 && nboot==3000)crit<-.001
if(d==15 && alpha==.05 && nboot==2000)crit<-.0016
if(d==15 && alpha==.025 && nboot==2000)crit<-.0005
if(d==15 && alpha==.05 && nboot==5000)crit<-.0026
if(d==15 && alpha==.025 && nboot==5000)crit<-.0006
}
if(is.na(crit) && alpha==.05)crit<-0.0268660714*(1/d)-0.0003321429
if(is.na(crit))crit<-alpha/(2*d)
if(d> 10 && nboot <5000){
print("Warning: Suggest using nboot=5000")
print("when the number of contrasts exceeds 10.")
}
icl<-round(crit*nboot)+1
icu<-round((1-crit)*nboot)
if(sum(con^2)==0){
con<-matrix(0,J,d)
id<-0
for (j in 1:Jm){
jp<-j+1
for (k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
psihat<-matrix(0,ncol(con),6)
dimnames(psihat)<-list(NULL,c("con.num","psihat","se","ci.lower",
"ci.upper","p-value"))
if(nrow(con)!=length(x)){
print("The number of groups does not match")
stop("the number of contrast coefficients.")
}
bvec<-matrix(NA,nrow=J,ncol=nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
for(j in 1:J){
paste("Working on group ",j)
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,mean,tr) # Bootstrapped trimmed means for jth group
}
test<-NA
for (d in 1:ncol(con)){
top<-0
for (i in 1:J){
top<-top+con[i,d]*bvec[i,]
}
test[d]<-(sum(top>0)+.5*sum(top==0))/nboot
test[d]<-min(test[d],1-test[d])
top<-sort(top)
psihat[d,4]<-top[icl]
psihat[d,5]<-top[icu]
}
for (d in 1:ncol(con)){
psihat[d,1]<-d
testit<-lincon(x,con[,d],tr,pr=FALSE)
psihat[d,6]<-2*test[d]
psihat[d,2]<-testit$psihat[1,2]
psihat[d,3]<-testit$test[1,4]
}
list(psihat=psihat,crit.p.value=2*crit,con=con)
}

comvar2d<-function(x,y,SEED=TRUE){
#
#  Compare the variances of two dependent groups.
#
nboot<-599
m<-cbind(x,y)
m<-elimna(m) # Remove missing values
U<-m[,1]-m[,2]
V<-m[,1]+m[,2]
ci<-pcorb(U,V,SEED=SEED)$ci
list(n=nrow(m),ci=ci)
}

mulwmw<-function(m1,m2,plotit=TRUE,cop=3,alpha=.05,nboot=1000,pop=4,fr=.8,pr=FALSE,SEED=TRUE,tr=.5,NC=TRUE){
#
#
# Determine center correpsonding to two
# independent groups, project all  points onto line
# connecting the centers,
# then based on the projected distances,
# estimate p=probability that a randomly sampled
# point from group 1 is less than a point from group 2
# based on the projected distances.
#
# plotit=TRUE creates a plot of the projected data
# pop=1 plot two dotplots based on projected distances
# pop=2 boxplots
# pop=3 expected frequency curve.
# pop=4 adaptive kernel density
#
#  There are three options for computing the center of the
#  cloud of points when computing projections:
#  cop=1 uses Donoho-Gasko median
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#
#  When using cop=2 or 3, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
#  NC=F: critical values not computed
#
#  Donoho-Gasko (Tukey) median is marked with a cross, +.
#
if(is.null(dim(m1))||dim(m1)[2]<2){print("Data are assumed to be stored in")
print(" a matrix or data frame having two or more columns.")
stop(" For univariate data, use the function outbox or out")
}
m1<-elimna(m1) # Remove missing values
m2<-elimna(m2)
n1=nrow(m1)
n2=nrow(m2)
if(cop==1){
if(ncol(m1)>2){
center1<-dmean(m1,tr=.5)
center2<-dmean(m2,tr=.5)
}
if(ncol(m1)==2){
tempd<-NA
for(i in 1:nrow(m1))
tempd[i]<-depth(m1[i,1],m1[i,2],m1)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center1<-m1[flag,]
if(sum(flag)>1)center1<-apply(m1[flag,],2,mean)
for(i in 1:nrow(m2))
tempd[i]<-depth(m2[i,1],m2[i,2],m2)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center2<-m2[flag,]
if(sum(flag)>1)center2<-apply(m2[flag,],2,mean)
}}
if(cop==2){
center1<-cov.mcd(m1)$center
center2<-cov.mcd(m2)$center
}
if(cop==3){
center1<-apply(m1,2,mean,tr=tr)
center2<-apply(m2,2,mean,tr=tr)
}
if(cop==4){
center1<-smean(m1)
center2<-smean(m2)
}
center<-(center1+center2)/2
B<-center1-center2
if(sum(center1^2)<sum(center2^2))B<-(0-1)*B
BB<-B^2
bot<-sum(BB)
disx<-NA
disy<-NA
if(bot!=0){
for (j in 1:nrow(m1)){
AX<-m1[j,]-center
tempx<-sum(AX*B)*B/bot
disx[j]<-sign(sum(AX*B))*sqrt(sum(tempx^2))
}
for (j in 1:nrow(m2)){
AY<-m2[j,]-center
tempy<-sum(AY*B)*B/bot
disy[j]<-sign(sum(AY*B))*sqrt(sum(tempy^2))
}
}
if(plotit){
if(pop==1){
par(yaxt="n")
xv<-rep(2,length(disx))
yv<-rep(1,length(disy))
plot(c(disx,disy),c(xv,yv),type="n",xlab="",ylab="")
xv<-rep(1.6,length(disx))
yv<-rep(1.4,length(disy))
points(disx,xv)
points(disy,yv)
}
if(pop==2)boxplot(disx,disy)
if(pop==3)rd2plot(disx,disy,fr=fr)
if(pop==4)g2plot(disx,disy,fr=fr)
}
m<-outer(disx,disy,FUN="-")
m<-sign(m)
phat<-(1-mean(m))/2
if(bot==0)phat<-.5
if(pr)print("Computing critical values")
m1<-t(t(m1)-center1)
m2<-t(t(m2)-center2)
v1=c(NA,NA)
if(NC)v1<-mulwmwcrit(m1,m2,cop=cop,alpha=alpha,iter=nboot,pr=pr,SEED=SEED)
list(phat=phat,lower.crit=v1[1],upper.crit=v1[2],n1=n1,n2=n2)
}

mulwmwcrit<-function(mm1,mm2,plotit=TRUE,cop=3,iter=1000,alpha=.05,SEED=TRUE,pr=FALSE){
#
#
# Determine critical value for the function mulwmw
#
if(!is.matrix(mm1))stop("Data are assumed to be stored in a matrix having two or more columns. For univariate data, use the function outbox or out")
#if(is.na(SEED))set.seed(2)
#if(!is.na(SEED))set.seed(SEED)
if(SEED)set.seed(2)
val<-NA
n1<-nrow(mm1)
n2<-nrow(mm2)
for(it in 1:iter){
ivec1<-sample(c(1:n1),replace=TRUE)
ivec2<-sample(c(1:n2),replace=TRUE)
m1<-mm1[ivec1,]
m2<-mm2[ivec2,]
if(cop==1){
if(ncol(m1)>2){
center1<-dmean(m1,tr=.5)
center2<-dmean(m2,tr=.5)
}
if(ncol(m1)==2){
tempd<-NA
for(i in 1:nrow(m1))
tempd[i]<-depth(m1[i,1],m1[i,2],m1)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center1<-m1[flag,]
if(sum(flag)>1)center1<-apply(m1[flag,],2,mean)
for(i in 1:nrow(m2))
tempd[i]<-depth(m2[i,1],m2[i,2],m2)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center2<-m2[flag,]
if(sum(flag)>1)center2<-apply(m2[flag,],2,mean)
}}
if(cop==2){
center1<-cov.mcd(m1)$center
center2<-cov.mcd(m2)$center
}
if(cop==3){
center1<-apply(m1,2,median)
center2<-apply(m2,2,median)
}
center<-(center1+center2)/2
B<-center1-center2
if(sum(center1^2)>sum(center2^2))B<-(0-1)*B
BB<-B^2
bot<-sum(BB)
disx<-NA
disy<-NA
if(bot!=0){
for (j in 1:nrow(m1)){
AX<-m1[j,]-center
tempx<-sum(AX*B)*B/bot
disx[j]<-sign(sum(AX*B))*sqrt(sum(tempx^2))
}
for (j in 1:nrow(m2)){
AY<-m2[j,]-center
tempy<-sum(AY*B)*B/bot
disy[j]<-sign(sum(AY*B))*sqrt(sum(tempy^2))
}}
m<-outer(disx,disy,FUN="-")
m<-sign(m)
val[it]<-(1-mean(m))/2
if(bot==0)val[it]<-.5
if(pr)print(paste("Iteration ",it," of ",iter," complete"))
}
val<-sort(val)
low<-round(alpha*iter/2)+1
up<-iter-low
crit<-NA
crit[1]<-val[low]
crit[2]<-val[up]
crit
}

mwmw<-function(m1,m2,cop=5,pr=TRUE,plotit=TRUE,pop=1,fr=.8,op=1,dop=1){
#
# Compute measure of effect size, p,
# a multivariate analog of Wilcoxon-Mann-Whitney p
#
# When plotting:
# pop=1 Use scatterplot
# pop=2 Use expected frequency curve.
# pop=3 Use adaptive kernel density
#
# dop=1, use method A1 approximation of halfspace depth
# dop=2, use method A2 approximation of halfspace depth
#
# cop determines how center of data is determined when
# approximating halfspace depth
# cop=1, Halfspace median
# cop=2, MCD
# cop=3, marginal medians
# cop=4, MVE
# cop=5, skipped mean
#
library(akima)
if(is.null(dim(m1)))stop("m1 is not a matrix or data frame")
if(is.null(dim(m2)))stop("m2 is not a matrix or data frame")
if(ncol(m1)!=ncol(m2))stop("number of columns for m1 and m2 are not equal")
if(ncol(m1)==1)stop("Use R function cid or bmp")
nn<-min(c(nrow(m1),nrow(m2)))
mdif<-matrix(as.vector(outer(m1[,1],m2[,1],"-")),ncol=1)
for(j in 2:ncol(m1)){
mdif<-cbind(mdif,matrix(as.vector(outer(m1[,j],m2[,j],"-")),ncol=1))
}
if(op==1){
if(ncol(m1)==2)temp2<-depth2(rbind(mdif,c(rep(0,ncol(m1)))))
#if(ncol(m1)==3)temp2<-depth3(rbind(mdif,c(rep(0,ncol(m1)))))
if(ncol(m1)>2){
if(cop==1)center<-dmean(mdif,tr=.5,dop=dop)
if(cop==2)center<-cov.mcd(mdif)$center
if(cop==3)center<-apply(mdif,2,median)
if(cop==4)center<-cov.mve(mdif)$center
if(cop==5)center<-smean(mdif)
temp2<-fdepth(rbind(mdif,c(rep(0,ncol(m1)))))
}}
if(op==2){
temp2<-pdis(rbind(mdif,c(rep(0,ncol(m1)))))
temp2<-1/(temp2+1)
}
center<-dmean(mdif,tr=.5,dop=dop)
phat<-temp2[nrow(mdif)+1]/max(temp2)
# phat is relative depth of zero vector
# Determine critical value
crit<-NA
alpha<-c(.1,.05,.025,.01)
crit[1]<-1-1.6338/sqrt(nn)
crit[2]<-1-1.8556/sqrt(nn)
crit[3]<-1-2.0215/sqrt(nn)
crit[4]<-1-2.1668/sqrt(nn)
if(pr){
print("For alpha=.1,.05,.025,.01, the correspoding critical values are")
print(crit)
print("Reject if phat is less than or equal to the critical value")
}
if(plotit && ncol(m1)==2){
if(pop==2)rdplot(mdif,fr=fr)
if(pop==1){
plot(mdif[,1],mdif[,2],xlab="VAR 1",ylab="VAR 2",type="n")
points(mdif[,1],mdif[,2],pch=".")
points(center[1],center[2],pch="o")
points(0,0,pch="+")
}
if(pop==3)akerdmul(mdif,fr=fr)
}
list(phat=phat,center=center,crit.val=crit)
}

TWOpNOV<-function(x,y,HC4=FALSE,alpha=.05){
#
#   Compute a .95 confidence interval
#   for the difference between two dependent Pearson correlations,
#   non-overlapping case.
#
#    Both x and y are assumed to be matrices with two columns.
#   The function compares the correlation between x[,1] and x[,2]
#   to the correlation between y[,1] and y[,2].
#
#  For simulation results, see Wilcox (2009).
#  COMPARING PEARSON CORRELATIONS: DEALING WITH
#  HETEROSCEDASTICITY AND NON-NORMALITY, Communications in Statistics--Simulations
#   and Computations, 38, 2220-2234.
#
#
if(!HC4 && alpha!=.05)stop('For alpha not equal to .05, must use HC4=TRUE')
#if(!is.matrix(x))stop("x should be a matrix")
#if(!is.matrix(y))stop("y should be a matrix")
if(ncol(x)!=2)stop("x should be a matrix or data a frame with 2 columns")
if(ncol(y)!=2)stop("y should be a matrix or a data frame with 2 columns")
xy=elimna(cbind(x,y))
x1=xy[,1]
x2=xy[,2]
y1=xy[,3]
y2=xy[,4]
r12=cor(x1,x2)
r13=cor(x1,y1)
r14=cor(x1,y2)
r23=cor(x2,y1)
r24=cor(x2,y2)
r34=cor(y1,y2)
term1=.5*r12*r34*(r13^2+r14^2+r23^2+r24^2)
term2=r12*r13*r14+r12*r23*r24+r13*r23*r34+r14*r24*r34
corhat=(term1+r13*r24+r14*r23-term2)/((1-r12^2)*(1-r34^2))
if(!HC4)temp=pcorbv4(x1,x2,SEED=FALSE)
if(HC4)temp=pcorhc4(x1,x2,alpha=alpha)
ci12=temp$ci[1]
ci12[2]=temp$ci[2]
if(!HC4)temp=pcorbv4(y1,y2,SEED=FALSE)
if(HC4)temp=pcorhc4(y1,y2,alpha=alpha)
ci34=temp$ci[1]
ci34[2]=temp$ci[2]
terml=2*corhat*(r12-ci12[1])*(ci34[2]-r34)
termu=2*corhat*(ci12[2]-r12)*(r34-ci34[1])
L=r12-r34-sqrt((r12-ci12[1])^2+(ci34[2]-r34)^2-terml)
U=r12-r34+sqrt((r12-ci12[2])^2+(ci34[1]-r34)^2-termu)
if(ZCI){
if(is.na(L) || is.na(U))L=U=0
}
list(est.1=r12,est.2=r34,ci.lower=L,ci.upper=U)
}

TWOpov<-function(x,y,alpha=.05,CN=FALSE,BOOT=TRUE, nboot=499,SEED=TRUE,ZCI=FALSE){
#
# Comparing two dependent correlations: Overlapping case
#
# x is assumed to be a matrix with 2 columns
#
#  Compare correlation of x[,1] with y to x[,2] with y
#
#  returns a confidence stored in
#  ci
#
if(ncol(x)!=2)stop('x should be a matrix with two columns')
x1y=elimna(cbind(x[,1],y))
x2y=elimna(cbind(x[,2],y))
xx=elimna(x)
r12=cor(x1y[,1],x1y[,2])
r13=cor(x2y[,1],x2y[,2])
r23=cor(xx[,1],xx[,2])
if(!BOOT){
ci12=pcorhc4(x1y[,1],x1y[,2],alpha=alpha,CN=CN)$ci
ci13=pcorhc4(x2y[,1],x2y[,2],alpha=alpha,CN=CN)$ci
}
if(BOOT){
ci12=rhohc4bt(x1y[,1],x1y[,2],alpha=alpha,SEED=SEED,nboot=nboot)$ci
ci13=rhohc4bt(x2y[,1],x2y[,2],alpha=alpha,SEED=SEED,nboot=nboot)$ci
}
corhat=((r23-.5*r12*r13)*(1-r12^2-r13^2-r23^2)+r23^3)/((1-r12^2)*(1-r13^2))
term1=2*corhat*(r12-ci12[1])*(ci13[2]-r13)
term2=2*corhat*(r12-ci12[2])*(ci13[1]-r13)
L=r12-r13-sqrt((r12-ci12[1])^2+(ci13[2]-r13)^2-term1)
U=r12-r13+sqrt((r12-ci12[2])^2+(ci13[1]-r13)^2-term2)
if(ZCI){
if(is.na(L) || is.na(U))L=U=0
}
list(est.rho1=r12,est.rho2=r13,dif=r12-r13,ci=c(L,U))
}

trimpb<-function(x,y=NULL,tr=.2,alpha=.05,nboot=2000,WIN=FALSE,win=.1,
plotit=FALSE,pop=1,null.value=0,pr=TRUE,xlab="X",fr=NA,SEED=TRUE){
#
#   Compute a 1-alpha confidence interval for
#   a trimmed mean.
#
#   The default number of bootstrap samples is nboot=2000
#
#   win is the amount of Winsorizing before bootstrapping
#   when WIN=T.
#
#   Missing values are automatically removed.
#
#  nv is null value. That test hypothesis trimmed mean equals nv
#
#  plotit=TRUE gives a plot of the bootstrap values
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate
#  pop=3 boxplot
#  pop=4 stem-and-leaf
#  pop=5 histogram
#  pop=6 adaptive kernel density estimate.
#
#  fr controls the amount of smoothing when plotting the bootstrap values
#  via the function rdplot. fr=NA means the function will use fr=.8
#  (When plotting bivariate data, rdplot uses fr=.6 by default.)
#
#  If y is not null, the function uses x-y; so can be used for two dependent variables.
#
if(pr){
print("The p-value returned by this function is based on the")
print("null value specified by the argument null.value, which defaults to 0")
}
if(!is.null(y))x=x-y
x<-x[!is.na(x)]
if(WIN){
if(win > tr)stop("The amount of Winsorizing must be <= to the amount of trimming")
x<-winval(x,win)
}
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
bvec<-NA
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,mean,tr) # Bootstrapped trimmed means
bvec<-sort(bvec)
#p.value<-sum(bvec<null.value)/nboot
p.value<-mean(bvec<null.value)+.5*mean(bvec==null.value)
p.value<-2*min(p.value,1-p.value)
ci<-NA
ci[1]<-bvec[icl]
ci[2]<-bvec[icu]
if(plotit){
if(pop==1)rdplot(as.vector(bvec),fr=fr,xlab=xlab)
if(pop==2)kdplot(as.vector(bvec),rval=rval)
if(pop==3)boxplot(as.vector(bvec))
if(pop==4)stem(as.vector(bvec))
if(pop==5)hist(as.vector(bvec))
if(pop==6)akerd(as.vector(bvec),xlab=xlab)
}
list(estimate=mean(x,tr=tr),ci=ci,p.value=p.value)
}

yuendna<-function(x,y=NULL,tr=.2,alpha=.05){
#
#  Compare the trimmed means of two dependent random variables
#  using the data in x and y.
#  The default amount of trimming is 20%
#
# If y is not supplied, this function assumes x is a matrix with 2 columns.
#
#  pairs of observations, for which one value is missing, are NOT deleted.
#  Marginal trimmed means are compared
#  using all available data.
#
if(is.null(y)){
if(!is.matrix(x))stop("y is null and x is not a matrix")
y=x[,2]
x=x[,1]
}
if(length(x)!=length(y))stop("The number of observations must be equal")
m<-cbind(x,y)
# first eliminate any rows with both values missing.
flag=(apply(is.na(m),1,sum)==2)
m=m[!flag,]
x<-m[,1]
y<-m[,2]
flagx=is.na(y) # Indicates observed x values for which y is missing
flagy=is.na(x) # Indicates the y values for which x is missing
m<-elimna(m)   # m has data where both values are available--no missing values
n=nrow(m)
n1=sum(flagx)  # number of x values for which y is missing
n2=sum(flagy)
h=n-2*floor(tr*n)
h1=n1-2*floor(tr*n1)
h2=n2-2*floor(tr*n2)
xbarn=mean(x,tr=tr,na.rm=TRUE)
xbarn1=0
if(h1>0)xbarn1=mean(x[flagx],tr=tr)
ybarn=mean(y[!flagy],tr=tr,na.rm=TRUE)
ybarn1=0
if(h2>0)ybarn1=mean(y[flagy],tr=tr)
lam1=h/(h+h1)
lam2=h/(h+h2)
est=lam1*xbarn-lam2*ybarn+(1-lam1)*xbarn1-(1-lam2)*ybarn1
sex=trimse(elimna(x),tr=tr)
sey=trimse(elimna(y),tr=tr)
q1<-(n-1)*winvar(m[,1],tr)
q2<-(n-1)*winvar(m[,2],tr)
q3<-(n-1)*wincor(m[,1],m[,2],tr)$cov
sen=sqrt((lam1^2*q1+lam2^2*q2-2*lam1*lam2*q3)/(h*(h-1)))
SE=sqrt(sen^2+(1-lam1)^2*sex^2+(1-lam2)^2*sey^2)
test=est/SE
list(estimate=est,test=test,se=SE)
}

yuenv2<-function(x,y=NULL,tr=.2,alpha=.05,plotit=FALSE,plotfun=splot,op=TRUE, VL=TRUE,cor.op=FALSE, loc.fun=median,
xlab="Groups",ylab="",PB=FALSE,nboot=100, SEED=TRUE){
#
#  Perform Yuen's test for trimmed means on the data in x and y.
#  The default amount of trimming is 20%
#  Missing values (values stored as NA) are automatically removed.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuen$ci.
#  The significance level is returned in yuen$p.value
#
#  For an omnibus test with more than two independent groups,
#  use t1way.
#
#   Unlike the function yuen, a robust heteroscedastic measure
#   of effect size is returned.
#  PB=FALSE means that a Winsorized variation of prediction error is used to measure effect size.
#  PB=TRUE:  A percentage bend measure of variation is used instead.
#
if(tr==.5)stop("Use medpb to compare medians.")
if(tr>.5)stop("Can't have tr>.5")
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
library(MASS)
x<-x[!is.na(x)]  # Remove any missing values in x
y<-y[!is.na(y)]  # Remove any missing values in y
n1=length(x)
n2=length(y)
h1<-length(x)-2*floor(tr*length(x))
h2<-length(y)-2*floor(tr*length(y))
q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
crit<-qt(1-alpha/2,df)
m1=mean(x,tr)
m2=mean(y,tr)
mbar=(m1+m2)/2
dif=m1-m2
low<-dif-crit*sqrt(q1+q2)
up<-dif+crit*sqrt(q1+q2)
test<-abs(dif/sqrt(q1+q2))
yuen<-2*(1-pt(test,df))
xx=c(rep(1,length(x)),rep(2,length(y)))
if(h1==h2){
pts=c(x,y)
top=var(c(m1,m2))
#
if(!PB){
if(tr==0)cterm=1
if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
bot=winvar(pts,tr=tr)/cterm
e.pow=top/bot
if(!is.na(e.pow)){
if(e.pow>1){
x0=c(rep(1,length(x)),rep(2,length(y)))
y0=c(x,y)
e.pow=wincor(x0,y0,tr=tr)$cor^2
}
}
}
#
if(PB){
bot=pbvar(pts)
e.pow=top/bot
}
#
}
if(n1!=n2){
N=min(c(n1,n2))
vals=0
if(SEED)set.seed(2)
for(i in 1:nboot)vals[i]=yuen.effect(sample(x,N),sample(y,N),tr=tr)$Var.Explained
e.pow=loc.fun(vals)
}
if(plotit){
plot(xx,pts,xlab=xlab,ylab=ylab)
if(op)
points(c(1,2),c(m1,m2))
if(VL)lines(c(1,2),c(m1,m2))
}
list(ci=c(low,up),n1=n1,n2=n2,
p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,
crit=crit,df=df,Var.Explained=e.pow,Effect.Size=sqrt(e.pow))
}

yuen.effect.ci<-function(x,y,SEED=TRUE,nboot=400,tr=.2,alpha=.05){
#
# Compute a 1-alpha  confidence interval
# for a robust, heteroscedastic  measure of effect size
#  The absolute value of the measure of effect size is used.
#
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
x=elimna(x)
y=elimna(y)
bvec=0
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(x)*nboot,replace=TRUE),nrow=nboot)
for(i in 1:nboot){
bvec[i]=yuenv2(datax[i,],datay[i,],tr=tr,SEED=FALSE)$Effect.Size
}
bvec<-sort(abs(bvec))
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
ci<-NA
ci[1]<-bvec[icl]
pchk=yuen(x,y,tr=tr)$p.value
if(pchk>alpha)ci[1]=0
ci[2]<-bvec[icu]
if(ci[1]<0)ci[1]=0
es=abs(yuenv2(x,y,tr=tr)$Effect.Size)
list(CI=ci,Effect.Size=es)
}

yuen.effect<-function(x,y,tr=.2,alpha=.05,plotit=FALSE,
plotfun=splot,op=TRUE,VL=TRUE,cor.op=FALSE,
xlab="Groups",ylab="",PB=FALSE){
#
#  Same as yuen, only it computes explanatory power and the related
# measure of effect size. Only use this with n1=n2. Called by yuenv2
# which allows n1!=n2.
#
#
#  Perform Yuen's test for trimmed means on the data in x and y.
#  The default amount of trimming is 20%
#  Missing values (values stored as NA) are automatically removed.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuen$ci.
#  The p-valueis returned in yuen$p.value
#
#  For an omnibus test with more than two independent groups,
#  use t1way.
#  This function uses winvar from chapter 2.
#
if(tr==.5)stop("Use medpb to compare medians.")
if(tr>.5)stop("Can't have tr>.5")
library(MASS)
x<-x[!is.na(x)]  # Remove any missing values in x
y<-y[!is.na(y)]  # Remove any missing values in y
h1<-length(x)-2*floor(tr*length(x))
h2<-length(y)-2*floor(tr*length(y))
q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
crit<-qt(1-alpha/2,df)
m1=mean(x,tr)
m2=mean(y,tr)
mbar=(m1+m2)/2
dif=m1-m2
low<-dif-crit*sqrt(q1+q2)
up<-dif+crit*sqrt(q1+q2)
test<-abs(dif/sqrt(q1+q2))
yuen<-2*(1-pt(test,df))
xx=c(rep(1,length(x)),rep(2,length(y)))
pts=c(x,y)
top=var(c(m1,m2))
#
if(!PB){
if(tr==0)cterm=1
if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
bot=winvar(pts,tr=tr)/cterm
}
if(PB)bot=pbvar(pts)/1.06
#
e.pow=top/bot
if(e.pow>1){
x0=c(rep(1,length(x)),rep(2,length(y)))
y0=c(x,y)
e.pow=wincor(x0,y0,tr=tr)$cor^2
}
if(plotit){
plot(xx,pts,xlab=xlab,ylab=ylab)
if(op)
points(c(1,2),c(m1,m2))
if(VL)lines(c(1,2),c(m1,m2))
}
list(ci=c(low,up),p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,
crit=crit,df=df,Var.Explained=e.pow,Effect.Size=sqrt(e.pow))
}

wmwloc2<-function(x,est=median){
#
# Compute loc2dif for all pairs of groups
#
if(is.matrix(x))x=listm(x)
locvec=NULL
ic=0
J=length(x)
for(j in 1:J){
for(k in 1:J){
if (j<k){
ic=ic+1
locvec[ic]=loc2dif(x[[j]],x[[k]],est=est)
}}}
locvec
}

cidM<-function(x,nboot=1000,alpha=.05,MC=FALSE,SEED=TRUE,g=NULL,dp=NULL){
#
# Variation of Cliff method based on median of X-Y
# i.e., use p=P(X<Y) as effect size.
# test p=.5
# All pairwise comparisons performed.
# FWE controlled via Hochberg method.
# x can be a matrix (columns are groups) or have list mode
#
#   g=NULL, x is assumed to be a matrix or have list mode
#   if g is specifed, it is assumed that column g of x is
#   a factor variable and that the dependent variable of interest is in column
#   dp of x, which can be a matrix or data frame.
#
if(!is.null(g)){
if(is.null(dp))stop("Specify a value for dp, the column containing the data")
x=fac2list(x[,dp],x[,g])
}
if(SEED)set.seed(2)
if(MC)library(parallel)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x=listm(x)
chk=tlist(x)
if(chk!=0)print("Warning: tied values detected. Suggest using cidmulv2")
J=length(x)
L=(J^2-J)/2
CC=L
pvec=NA
boot=list()
MAT=matrix(NA,nrow=nboot,ncol=L)
for(i in 1:nboot){
jcom=0
for (j in 1:J){
boot[[j]]=sample(x[[j]],size=length(x[[j]]),replace=TRUE)
}
MAT[i,]=wmwloc2(boot)
}
#
pvec=NA
test<-matrix(NA,CC,8)
dimnames(test)<-list(NULL,c("Group","Group","p-value","p.crit",
"P(X<Y)","P(X=Y)","P(X>Y)","p.hat"))
dvec<-alpha/c(1:CC)
for(j in 1:J){
for(k in 1:J){
if(j<k){
jcom=jcom+1
p.value=mean(MAT[,jcom]>0)+.5*mean(MAT[,jcom]==0)
pvec[jcom]=2*min(c(p.value,1-p.value))
if(is.na(pvec[jcom]))pvec=1
test[jcom,1]<-j
test[jcom,2]<-k
test[jcom,3]<-pvec[jcom]
test[jcom,5:7]<-cid(x[[j]],x[[k]])$summary.dvals
test[jcom,8]<-test[jcom,5]+.5*test[jcom,6]
}}}
temp2<-order(0-test[,3])
test[temp2,4]=dvec
list(test=test)
}

cidmul<-function(x,alpha=.05,g=NULL,dp=NULL,pr=TRUE){
#
#  Perform Cliff's method for all pairs of J independent groups.
#  Unlike the function meemul, ties are allowed.
#  The familywise type I error probability is controlled by using
#  a critical value from the Studentized maximum modulus distribution.
#
#  The data are assumed to be stored in $x$ in list mode.
#  Length(x) is assumed to correspond to the total number of groups, J.
#  It is assumed all groups are independent.
#
#  Missing values are automatically removed.
#
#  The default value for alpha is .05. Any other value results in using
#  alpha=.01.
#
if(pr)print('cidmulv2 might provide better power')
if(!is.null(g)){
if(is.null(dp))stop("Specify a value for dp, the column containing the data")
x=fac2list(x[,dp],x[,g])
}
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
J<-length(x)
CC<-(J^2-J)/2
test<-matrix(NA,CC,7)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
}
dimnames(test)<-list(NULL,c("Group","Group","d","ci.lower","ci.upper",
"p.hat","p-value"))
jcom<-0
crit<-smmcrit(200,CC)
if(alpha!=.05)crit<-smmcrit01(200,CC)
alpha<-1-pnorm(crit)
n=matl(lapply(x,length))
for (j in 1:J){
for (k in 1:J){
if (j < k){
temp<-cid(x[[j]],x[[k]],alpha,plotit=FALSE)
temp2<-cidv2(x[[j]],x[[k]],alpha,plotit=FALSE)
jcom<-jcom+1
test[jcom,1]<-j
test[jcom,2]<-k
test[jcom,3]<-temp$d
test[jcom,4]<-temp$cl
test[jcom,5]<-temp$cu
test[jcom,6]<-temp$phat
test[jcom,7]<-temp2$p.value
}}}
list(n=n,test=test)
}

medpb2<-function(x,y=NULL,alpha=.05,nboot=2000,SEED=TRUE){
#
#   Compare 2 independent groups using medians.
#
#   A percentile bootstrap method is used, which performs well when
#   there are tied values.
#
#   The data are assumed to be stored in x and y. If y=NULL, x is assumed to have two columns.
#
#   Missing values are automatically removed.
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
x=elimna(x)
y=elimna(y)
xx<-list()
xx[[1]]<-x
xx[[2]]<-y
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
est1=median(xx[[1]])
est2=median(xx[[2]])
est.dif<-median(xx[[1]])-median(xx[[2]])
crit<-alpha/2
temp<-round(crit*nboot)
icl<-temp+1
icu<-nboot-temp
bvec<-matrix(NA,nrow=2,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:2){
data<-matrix(sample(xx[[j]],size=length(xx[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,median) # Bootstrapped medians for jth group
}
top<-bvec[1,]-bvec[2,]
test<-sum(top<0)/nboot+.5*sum(top==0)/nboot
if(test > .5)test<-1-test
top<-sort(top)
ci<-NA
ci[1]<-top[icl]
ci[2]<-top[icu]
list(n1=length(x),n2=length(y),p.value=2*test,ci=ci,est1=est1,est2=est2,
est.dif=est.dif)
}

cid<-function(x,y,alpha=.05,plotit=FALSE,pop=0,fr=.8,rval=15,xlab="",ylab=""){
#
# For two independent groups,
#  compute a confidence interval for  P(X<Y). The method stems from
#  Cliff, 1996, p. 140, eq 5.12. Tied values are allowed.
#
#  To compare the lower and upper quantiles of the distribution of D=X-Y,
#  use cbmhd.
#
#  This function also  reports a 1-alpha confidence interval for
#  P(X>Y)-P(X<Y)
#
#  Let xi_q be the qth quantile of the distribution of D, q<.5. To test xi_q +xi_1-q=0, use cbmhd
#  If nothing is going on, D is symmetric about zero, so this function tests for symmetry and provides
# some sense of how the tails of the distribution D differ.
#  qwmwhd applies the method using a range of q values
#
#
#  plotit=TRUE creates a plot of the difference scores.
#  pop=0 adaptive kernel density estimate
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate (Rosenblatt's shifted histogram)
#  pop=3 boxplot
#  pop=4 stem-and-leaf
#  pop=5 histogram
#  pop=6  kernel density estimate
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
if(length(x)*length(y)>10^6)stop('Use bmp with a large sample size. If using rimul, use ribmp instead')
m<-outer(x,y,FUN="-")
msave<-m
m<-sign(m)
d<-mean(m)
phat<-(1-d)/2
flag=TRUE
if(phat==0 || phat==1)flag=FALSE
q0<-sum(msave==0)/length(msave)
qxly<-sum(msave<0)/length(msave)
qxgy<-sum(msave>0)/length(msave)
c.sum<-matrix(c(qxly,q0,qxgy),nrow=1,ncol=3)
dimnames(c.sum)<-list(NULL,c("P(X<Y)","P(X=Y)","P(X>Y)"))
if(flag){
sigdih<-sum((m-d)^2)/(length(x)*length(y)-1)
di<-NA
for (i in 1:length(x))di[i]<-sum(x[i]>y)/length(y)-sum(x[i]<y)/length(y)
dh<-NA
for (i in 1:length(y))dh[i]<-sum(y[i]>x)/length(x)-sum(y[i]<x)/length(x)
sdi<-var(di)
sdh<-var(dh)
sh<-((length(y)-1)*sdi+(length(x)-1)*sdh+sigdih)/(length(x)*length(y))
zv<-qnorm(alpha/2)
cu<-(d-d^3-zv*sqrt(sh)*sqrt((1-d^2)^2+zv^2*sh))/(1-d^2+zv^2*sh)
cl<-(d-d^3+zv*sqrt(sh)*sqrt((1-d^2)^2+zv^2*sh))/(1-d^2+zv^2*sh)
}
if(!flag){
sh=NULL
nm=max(c(length(x),length(y)))
if(phat==1)bci=binomci(nm,nm,alpha=alpha)
if(phat==0)bci=binomci(0,nm,alpha=alpha)
}
if(plotit){
if(pop==1 || pop==0){
if(length(x)*length(y)>2500){
print("Product of sample sizes exceeds 2500.")
print("Execution time might be high when using pop=0 or 1")
print("If this is case, might consider changing the argument pop")
}}
if(pop==0)akerd(as.vector(msave),xlab=xlab,ylab=ylab)
if(pop==1)rdplot(as.vector(msave),fr=fr,xlab=xlab,ylab=ylab)
if(pop==2)kdplot(as.vector(msave),rval=rval,xlab=xlab,ylab=ylab)
if(pop==3)boxplot(as.vector(msave))
if(pop==4)stem(as.vector(msave))
if(pop==5)hist(as.vector(msave),xlab=xlab)
if(pop==6)skerd(as.vector(msave))
}
if(flag)pci=c((1-cu)/2,(1-cl)/2)
if(!flag){
pci=bci$ci
cl=1-2*pci[2]
cu=1-2*pci[1]
}
list(n1=length(x),n2=length(y),cl=cl,cu=cu,d=d,sqse.d=sh,phat=phat,summary.dvals=c.sum,ci.p=pci)
}

cidv2<-function(x,y,alpha=.05,plotit=FALSE,pop=0,fr=.8,rval=15,xlab='',ylab=''){
#
#   p-value for Cliff's analog of WMW test
#
#  To compare the lower and upper quantiles of the distribution of D=X-Y,
#  use cbmhd.
#
if(length(x)*length(y)>10^6)stop('Use bmp with a large sample size.')
nullval<-0
ci<-cid(x,y,alpha=alpha,plotit=plotit,pop=pop,fr=fr,rval=rval)
FLAG=TRUE
if(ci$phat==0 || ci$phat==1)FLAG=FALSE
if(FLAG){
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-cid(x,y,alpha=alph[i],plotit=FALSE)
if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
}
p.value<-irem/100
if(p.value<=.01){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-cid(x,y,alpha=alph[i],plotit=FALSE,xlab=xlab,ylab=ylab)
if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
}}
if(p.value<=.001){
alph<-seq(.0001,.001,.0001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-cid(x,y,alpha=alph[i],plotit=FALSE)
if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
}}
phat<-(1-ci$d)/2
pci=c((1-ci$cu)/2,(1-ci$cl)/2)
d.ci=c(ci$cl,ci$cu)
dval=cid(x,y)$summary.dvals
}
if(!FLAG){
D=bmp(x,y)
p.value=D$p.value
d.ci=NA
pci=D$ci.p
phat=D$phat
dval=ci$summary.dvals
}
list(n1=length(elimna(x)),n2=length(elimna(y)),d.hat=ci$d,d.ci=d.ci,p.value=p.value,p.hat=phat,p.ci=pci,summary.dvals=dval)
}

pb2trmcp<-function(J,K,x,grp=c(1:p),p=J*K,tr=.2,nboot=NA,alpha=.05,SEED=TRUE,pr=TRUE,
bhop=FALSE){
#
#  Perform a J by K anova using trimmed means with
#  for two independent groups using a bootstrap-t method
#
#  tr=.2 is default trimming
#
#
#  The R variable data is assumed to contain the raw
#  data stored in list mode. data[[1]] contains the data
#  for the first level of both factors: level 1,1.
#  data[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  data[[K]] is the data for level 1,K
#  data[[K+1]] is the data for level 2,1, data[2K] is level 2,K, etc.
#
#  It is assumed that data has length JK, the total number of
#  groups being tested, but a subset of the data can be analyzed
#  using grp
#
if(SEED)set.seed(2)
if(is.list(x))x<-elimna(matl(x))
if(is.matrix(x))x<-elimna(x)
data<-x
if(is.matrix(data))data<-listm(data)
if(!is.list(data))stop("Data are not stored in list mode or a matrix")
if(p!=length(data)){
print("The total number of groups, based on the specified levels, is")
print(p)
print("The number of groups stored in x is")
print(length(data))
print("Warning: These two values are not equal")
}
if(p!=length(grp))stop("Apparently a subset of the groups was specified that does not match the total number of groups indicated by the values for J and K.")
temp=con2way(J,K)
conA<-temp$conA
conB<-temp$conB
conAB<-temp$conAB
if(pr)print("Taking bootstrap samples")
Factor.A<-pbtrmcp(x,con=conA,tr=tr,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
Factor.B<-pbtrmcp(x,con=conB,tr=tr,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
Factor.AB<-pbtrmcp(x,con=conAB,tr=tr,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.AB=Factor.AB,bhop=bhop,SEED=FALSE)
}

wmwaov<-function(x,est=median,nboot=500,MC=FALSE,SEED=TRUE,MM=FALSE){
#
# Extension of WMW to J groups
# i.e., use p=P(X<Y) as effect size.
# test p_{jk}=.5 all j<k
#
if(SEED)set.seed(2)
if(MC)library(parallel)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x=listm(x)
chk=tlist(x)
if(chk!=0){
if(identical(est,median))print("Warning: tied values detected. Suggest using est=hd or the function cidmulv2")
}
J=length(x)
L=(J^2-J)/2
ic=0
pvec=NA
boot=list()
MAT=matrix(NA,nrow=nboot,ncol=L)
for(i in 1:nboot){
for (j in 1:J){
boot[[j]]=sample(x[[j]],size=length(x[[j]]),replace=TRUE)
}
MAT[i,]=wmwloc2(boot,est=est)
}
zero=rep(0,L)
bconB=rbind(MAT,zero)
if(MC)dv=pdisMC(bconB,MM=MM)
if(!MC)dv=pdis(bconB,MM=MM)
bplus<-nboot+1
p.value<-1-sum(dv[bplus]>dv[1:nboot])/nboot-.5*sum(dv[bplus]==dv[1:nboot])/nboot
p.value
}

dtrimpb<-function(x,y=NULL,alpha=.05,con=0,est=tmean,plotit=TRUE,dif=TRUE,grp=NA,
hoch=TRUE,nboot=NA,xlab="Group 1",ylab="Group 2",
pr=TRUE,SEED=TRUE,BA=FALSE,PCI=FALSE,ylab.ebar=NULL,...){
#
#   Use a percentile bootstrap method to  compare
#   trimmed means  of dependent groups.
#
#   This is essentially the function rmmcppb, but set to compare trimmed means
#   by default.
#
#   By default,
#   compute a .95 confidence interval for all linear contasts
#   specified by con, a J by C matrix, where  C is the number of
#   contrasts to be tested, and the columns of con are the
#   contrast coefficients.
#   If con is not specified, all pairwise comparisons are done.
#
#   A sequentially rejective method
#   is used to control the probability of at least one Type I error.
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
#  PCI=TRUE, if dif=TRUE and est=median, confidence intervals for difference scores are plottted
#           So this is like plotting error bars.
#          The y-axis of the plot can be specified with the argument ylab.ebar
#
if(dif){
if(pr)print("dif=T, so analysis is done on difference scores")
temp<-rmmcppbd(x,y=y,alpha=.05,con=con,est=est,plotit=plotit,grp=grp,
nboot=nboot,hoch=hoch,...)
output<-temp$output
con<-temp$con
}
if(!dif){
if(pr)print("dif=F, so analysis is done on marginal distributions")
if(!is.null(y[1]))x<-cbind(x,y)
if(is.data.frame(x))x=as.matrix(x)

if(!is.list(x) && !is.matrix(x))
stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
if(is.matrix(con)){
if(length(x)!=nrow(con))
stop("The number of rows in con is not equal to the number of groups.")
}}
if(is.list(x)){
# put the data in an n by J matrix
mat<-matl(x)
}
if(is.matrix(x) && is.matrix(con)){
if(ncol(x)!=nrow(con))
stop("The number of rows in con is not equal to the number of groups.")
mat<-x
}
if(is.matrix(x))mat<-x
if(!is.na(sum(grp)))mat<-mat[,grp]
mat<-elimna(mat) # Remove rows with missing values.
x<-mat
J<-ncol(mat)
xcen<-x
for(j in 1:J)xcen[,j]<-x[,j]-est(x[,j])
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
xbars<-apply(mat,2,est)
psidat<-NA
for (ic in 1:connum)psidat[ic]<-sum(con[,ic]*xbars)
psihat<-matrix(0,connum,nboot)
psihatcen<-matrix(0,connum,nboot)
bvec<-matrix(NA,ncol=J,nrow=nboot)
bveccen<-matrix(NA,ncol=J,nrow=nboot)
print("Taking bootstrap samples. Please wait.")
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
tval<-NA
tvalcen<-NA
for (ic in 1:connum){
psihat[ic,]<-apply(bvec,1,bptdpsi,con[,ic])
psihatcen[ic,]<-apply(bveccen,1,bptdpsi,con[,ic])
tvalcen[ic]<-sum((psihatcen[ic,]==0))/nboot
bias[ic]<-sum((psihatcen[ic,]>0))/nboot+sum((psihatcen[ic,]==0))/nboot-.5
tval[ic]<-sum((psihat[ic,]==0))/nboot
if(BA){
test[ic]<-sum((psihat[ic,]>0))/nboot+tval[ic]-.1*bias[ic]
if(test[ic]<0)test[ic]<-0
}
if(!BA)test[ic]<-sum((psihat[ic,]>0))/nboot+tval[ic]
test[ic]<-min(test[ic],1-test[ic])
}
test<-2*test
ncon<-ncol(con)
if(alpha==.05){
dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
dvecba<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
dvecba<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(hoch)dvec<-alpha/(2* c(1:ncon))
dvec<-2*dvec
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
dvecba<-dvec
dvec[1]<-alpha/2
}
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
dimnames(output)<-list(NULL,c("con.num","psihat","p-value","p.crit",
"ci.lower","ci.upper"))
tmeans<-apply(mat,2,est,...)
psi<-1
for (ic in 1:ncol(con)){
output[ic,2]<-sum(con[,ic]*tmeans)
output[ic,1]<-ic
output[ic,3]<-test[ic]
output[temp2,4]<-zvec
temp<-sort(psihat[ic,])
icl<-round(output[ic,4]*nboot/2)+1
icu<-nboot-(icl-1)
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
}
if(PCI){
if(dif){
plotCI(output[,2],ali=output[,5],aui=output[,6],xlab='Difference',ylab=ylab.ebar)
}}
num.sig<-sum(output[,3]<=output[,4])
list(output=output,con=con,num.sig=num.sig)
}

binband<-function(x,y,KMS=FALSE,alpha=.05, plotit=TRUE,xlab="X",    #ADJ.P=FALSE, old code deleted.
ylab="Rel. Freq.", method='hoch',pr=TRUE){
#
#  Comparing two independent variables in terms of their probability function.
#  For each value that occurs, say x, test P(X=x)=P(Y=x)
#  So this method is useful when dealing with highly discrete data.
#
#  If KMS=TRUE, use Kulinskaya, Morgenthaler and Staudte (2010)
#   method for comparing binomials
# Kulinskaya, E., Morgenthaler, S. and Staudte, R. (2010).
# Variance Stabilizing the Difference of two Binomial
#  Proportions. American Statistician, 64,
#  350--356 DOI:10.1198/tast.2010.09096

#  Otherwise use Storer and Kim.
#
#  method='hoch': p-values are adjusted via Hochberg's method
#
if(!KMS){
if(pr)print('To get confidence intervals, set KMS=TRUE')
}
x=elimna(x)
y=elimna(y)
vals=sort(unique(c(x,y)))
ncon=length(vals)
n1=length(x)
n2=length(y)
p.values=NA
adj=1
cv=1
if(!KMS){
output=matrix(NA,ncol=6,nrow=length(vals))
dimnames(output)=list(NULL,c("Value","p1.est","p2.est","p1-p2","p.value","p.adj"))
}
if(KMS){
output=matrix(NA,ncol=8,nrow=length(vals))
dimnames(output)=list(NULL,c("Value","p1.est","p2.est","p1-p2","ci.low","ci.up","p.value",
"p.adj"))
}
for(i in 1:length(vals)){
x1=sum(x==vals[i])
y1=sum(y==vals[i])
if(!KMS){
output[i,5]=twobinom(x1,n1,y1,n2)$p.value
output[i,2]=x1/n1
output[i,3]=y1/n2
output[i,1]=vals[i]
output[i,4]=output[i,2]-output[i,3]
}
if(KMS){
temp=bi2KMSv2(x1,n1,y1,n2)
output[i,1]=vals[i]
output[i,5]=temp$ci[1]
output[i,6]=temp$ci[2]
output[i,2]=x1/n1
output[i,3]=y1/n2
output[i,4]=output[i,2]-output[i,3]
output[i,7]=temp$p.value
}}
# Determine adjusted  critical p-value using Hochberg method
ncon=length(vals)
dvec=alpha/c(1:ncon)

if(KMS){
output[,8]=p.adjust(output[,7],method=method)
}
if(!KMS){
output[,6]=p.adjust(output[,5],method=method)
}
if(plotit)splotg2(x,y, xlab=xlab, ylab=ylab)
output
}

tworegwb<-function(x1,y1,x2,y2,nboot=599,RAD=FALSE,alpha=.05,SEED=TRUE,xout=FALSE,
outfun=out){
#
# Simple regression (one predictor)
# Test H_0: two independent groups have equal slopes.
#
xy=elimna(cbind(x1,y1))
if(ncol(xy)>2)stop("This function only allows one covariate")
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
x=c(x1,x2)
y=c(y1,y2)
g=c(rep(0,length(x1)),rep(1,length(x2)))
xgy=elimna(cbind(x,g,x*g,y))
xg=xgy[,1:3]
y=xgy[,4]
res=olswbtest(xg,y,nboot=nboot,SEED=SEED,RAD=RAD,alpha=alpha)
res[3,6]
}

disc2comSK<-function(x,y,alpha=.05,nboot=500,SEED=TRUE){
#
#  Comparing two independent variables in terms of their probability function.
#  A global test of P(X=x)=P(Y=x) for all x.
#  Appears  to have no advantage over a chi-square test done by the R function disc2com
#
#  The R function binband tests this hypothesis for each x.
#
library(mc2d)
x=elimna(x)
y=elimna(y)
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
vals=sort(unique(c(x,y)))
n1=length(x)
n2=length(y)
K=length(vals)
C1=NULL
C2=NULL
HT=NULL
for(i in 1:K){
C1[i]=sum(x==vals[i])
C2[i]=sum(y==vals[i])
HT[i]=(C1[i]+C2[i])/(n1+n2)
}
p1hat=C1/n1
p2hat=C2/n2
test=sum((p1hat-p2hat)^2)
tv=NULL
TB=NA
VP=NA
for(ib in 1:nboot){
xx=rmultinomial(n1,1,HT)
yy=rmultinomial(n2,1,HT)
B1=NA
B2=NA
BP=NA
for(i in 1:K){
B1[i]=sum(xx[,i])
B2[i]=sum(yy[,i])
}
B1hat=B1/n1
B2hat=B2/n2
TB[ib]=sum((B1hat-B2hat)^2)
}
pv=1-mean(test>TB)-.5*mean(test==TB)
list(test=test,p.value=pv)
}

wmwloc<-function(x,y,na.rm=TRUE,est=median,...){
#
# Estimate the median of the distribution of x-y
#
if(na.rm){
x<-x[!is.na(x)]
y<-y[!is.na(y)]
}
m<-outer(x,y,FUN="-")
est=est(m,na.rm=TRUE,...)
est
}

Dqcomhd<-function(x,y,est=hd,q=c(1:9)/10,nboot=2000,pr=TRUE,
plotit=FALSE,SEED=TRUE,xlab='Group 1',
ylab='Est.1-Est.2',na.rm=TRUE,alpha=rep(.05,length(q))){
#
# Compare the quantiles of the marginal distributions associated with  two dependent groups
# via hd estimator. Tied values are allowed.
#
# est=thd would use trimmed hd estimator
#
# When comparing lower or upper quartiles, both power and the probability of Type I error
# compare well to other methods have been derived.
#
#  x: data for group 1
#  y: data for group 2
#  q: the quantiles to be compared
#  nboot: Number of bootstrap samples
#
#
if(pr){
print('Note: confidence  intervals are not adjusted to control the simultaneous probability coverage')
}
if(SEED)set.seed(2)
if(na.rm){
xy=elimna(cbind(x,y))
x=xy[,1]
y=xy[,2]
}
pv=NULL
output=matrix(NA,nrow=length(q),ncol=10)
dimnames(output)<-list(NULL,c('q','n1','n2','est.1','est.2','est.1_minus_est.2','ci.low','ci.up','p-value','adj.p.value'))
for(i in 1:length(q)){
output[i,1]=q[i]
output[i,2]=length(elimna(x))
output[i,3]=length(elimna(y))
output[i,4]=hd(x,q=q[i])
output[i,5]=hd(y,q=q[i])
output[i,6]=output[i,4]-output[i,5]
if(na.rm){
temp=bootdpci(x,y,est=est,q=q[i],dif=FALSE,plotit=FALSE,pr=FALSE,nboot=nboot,alpha=alpha[i],SEED=FALSE)
output[i,7]=temp$output[1,5]
output[i,8]=temp$output[1,6]
output[i,9]=temp$output[1,3]
}
if(!na.rm){
temp=rmmismcp(x,y,est=est,q=q[i],plotit=FALSE,pr=FALSE,nboot=nboot,alpha=alpha[i],SEED=FALSE)
output[i,7]=temp$output[1,6]
output[i,8]=temp$output[1,7]
output[i,9]=temp$output[1,4]
}
}
output[,10]=p.adjust(output[,9],method='hoch')
if(plotit){
xax=rep(output[,4],3)
yax=c(output[,6],output[,7],output[,8])
plot(xax,yax,xlab=xlab,ylab=ylab,type='n')
points(output[,4],output[,6],pch='*')
lines(output[,4],output[,6])
points(output[,4],output[,7],pch='+')
lines(output[,4],output[,7],lty=2)
points(output[,4],output[,8],pch='+')
lines(output[,4],output[,8],lty=2)
}
output
}

qwmwhd<-function(x,y,q=seq(5,40,5)/100,xlab="Quantile",ylab="Sum of q and 1-q Quantiles",plotit=TRUE,alpha=.05,nboot=1000,SEED=TRUE){
#
#  Plot that provides perspective on the degree a distribution is symmetric about zero.
#  This function plots the sum of q and 1-q quantiles of the distribution of D=X-Y, X and Y independent.
#  A 1-alpha confidence interval for the sum is indicated by a +
#  If the distribution is symmetric
#  the plot should be approximately a horizontal line.
#
#  FWE is controlled via Hochberg's method, which was used to determine critical
#  p-values based on the argument
#  alpha.
#
#  Can alter the quantiles compared via the argument
#  q
#  q must be less than .5
#
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
output=matrix(NA,ncol=8,nrow=length(q))
dimnames(output)=list(NULL,c("quantile","Est.1","Est.2","SUM","ci.low","ci.up","p_crit","p-value"))
for(i in 1:length(q)){
test=cbmhd(x,y,q=q[i],plotit=FALSE,nboot=nboot,SEED=SEED)
output[i,1]=q[i]
output[i,2]=test$Est1
output[i,3]=test$Est2
output[i,4]=test$sum
output[i,8]=test$p.value
output[i,5]=test$ci[1]
output[i,6]=test$ci[2]
}
temp=order(output[,8],decreasing=TRUE)
zvec=alpha/c(1:length(q))
output[temp,7]=zvec
output <- data.frame(output)
output$signif=rep("YES",nrow(output))
for(i in 1:nrow(output)){
if(output[temp[i],8]>output[temp[i],7])output$signif[temp[i]]="NO"
if(output[temp[i],8]<=output[temp[i],7])break
}
if(plotit){
plot(rep(q,3),c(output[,4],output[,5],output[,6]),type="n",xlab=xlab,ylab=ylab)
points(q,output[,6],pch="+")
points(q,output[,5],pch="+")
points(q,output[,4],pch="*")
}
list(n=c(n1,n2),output=output)
}

yuendv2<-function(x, y, tr = 0.2, alpha = 0.05,null.value=0,pr=TRUE){
#
#  Same as yuend, only it also returns a measure of
#  effect size similar to the one used by yuenv2.
# To get a measure of effect size based on the difference scores,
#  use the function trimci or trimcipb
#  (est.dif - null.value)/sd
#  For trimmed means, sd is a Winsorized variance
#  rescaled so that it estimates the standard deviation under normality
#
if(pr)print('This version returns an effect size similar to what is used by yuenv2')
if(pr)print('To get a measure of effect size based on the difference scores, use trimciv2')
library(MASS)
if(tr<0)stop('tr must be between 0 and .5')
if(tr>.5)stop('tr must be between 0 and .5')
res=yuend(x=x,y=y,tr=tr,alpha=alpha)
#
#if(tr==0)term=1
#if(tr>0)term=sqrt(area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr)
#epow=(res$dif-null.value)*term/sqrt(winvar(x-y,tr=tr,na.rm=TRUE))
epow=yuenv2(x,y,tr=tr)$Effect.Size
list(ci=res$ci,p.value=res$p.value,est1=res$est1,est2=res$est2,dif=res$dif,se=res$se,
teststat=res$teststat,n=res$n,df=res$df,Effect.Size=epow)
}

ancovaWMW<-function(x1,y1,x2,y2,fr1=1,fr2=1,alpha=.05,sm=FALSE,est=tmean,
plotit=TRUE,pts=NA,xout=FALSE,outfun=out,LP=TRUE,...){
#
# Compare two independent  groups using the ancova method in conjunction
# with Cliff's improvement on the Wilcoxon-Mann-Whitney test.
# No parametric assumption is made about the form of
# the regression lines--a running interval smoother is used.
#
#  Assume data are in x1 y1 x2 and y2
#
#  OLD version: sm=TRUE will use bootstrap bagging when plotting the regression lines
#  The plot is based on measure of location indicated by the argument
#  est. Default is the Harrell-Davis estimate of the median.  Not working, took this out.
#
#   LP=TRUE: use running interval smoother followed by LOESS
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
flag<-outfun(x1,...)$keep
x1<-x1[flag]
y1<-y1[flag]
flag<-outfun(x2,...)$keep
x2<-x2[flag]
y2<-y2[flag]
}
dv.sum=NULL
if(is.na(pts[1])){
npt<-5
CC=5
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
dimnames(mat)<-list(NULL,c('X','n1','n2','p.hat','ci.low','ci.hi','p.value','p.crit'))
for (i in 1:5){
g1<-y1[near(x1,x1[isub[i]],fr1)]
g2<-y2[near(x2,x1[isub[i]],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test<-cidv2(g1,g2,alpha=alpha)
dv.sum=rbind(dv.sum,test$summary.dvals)
mat[i,1]<-x1[isub[i]]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
mat[i,4]<-test$p.hat
mat[i,5]<-test$p.ci[1]
mat[i,6]<-test$p.ci[2]
mat[i,7]<-test$p.value
}}
if(!is.na(pts[1])){
CC=length(pts)
n1<-1
n2<-1
vecn<-1
for(i in 1:length(pts)){
n1[i]<-length(y1[near(x1,pts[i],fr1)])
n2[i]<-length(y2[near(x2,pts[i],fr2)])
}
mat<-matrix(NA,length(pts),8)
dimnames(mat)<-list(NULL,c('X','n1','n2','p.hat','ci.low','ci.hi','p.value','p.crit'))
for (i in 1:length(pts)){
g1<-y1[near(x1,pts[i],fr1)]
g2<-y2[near(x2,pts[i],fr2)]
g1<-g1[!is.na(g1)]
g2<-g2[!is.na(g2)]
test=cidv2(g1,g2,alpha=alpha)
dv.sum=rbind(dv.sum,test$summary.dvals)
mat[i,1]<-pts[i]
mat[i,2]<-length(g1)
mat[i,3]<-length(g2)
if(length(g1)<=5)print(paste('Warning, there are',length(g1),' points corresponding to the design point X=',pts[i]))
if(length(g2)<=5)print(paste('Warning, there are',length(g2),' points corresponding to the design point X=',pts[i]))
mat[i,4]<-test$p.hat
mat[i,5]<-test$p.ci[1]
mat[i,6]<-test$p.ci[2]
mat[i,7]<-test$p.value
}}
dvec<-alpha/c(1:CC)
temp2<-order(0-mat[,6])
mat[temp2,8]=dvec
if(plotit){
runmean2g(x1,y1,x2,y2,fr=fr1,est=est,sm=sm,xout=FALSE,LP=LP,...)
}
list(output=mat,summary=dv.sum)
}

twoKlin<-function(x=NULL,x1=NULL,x2=NULL,tr=.2,alpha=.05,pr=TRUE,opt=1){
#
#  A step-down MCP based on K independent tests.
#  It is essential that the tests are independent.
#
#  Use Fisher method based on p-values coupled with Hochberg
#
# Data are assumed to be stored in two R variables, x1 and x2 or in one
#  R variable, x
#
# If stored in x1 and x2, they are assumed to be matrices with K columns
# or to have list mode, both having length K.
#
# If the data are stored in x,
# x is assumed to have 2K columns if a matrix or length 2K if it has list mode.
#
# If data are stored in x1 and x2, for each column, compute a p-value.
# That is, perform a test based on the data in column 1 of x1 and x2,
# followed by a test using the data in column 2 of x1 and x2, etc.
#
# If data are stored in x, the first test is based
# on the data in columns 1 and K+1,
# the second test is based on columns 2 and K+2, etc.
#
#  opt=1  Fisher's method
#  opt=2  Chen-Nadarajah method
#  opt=3  Max method
#
if(is.null(x[1])){
if(is.matrix(x1))x=cbind(x1,x2)
if(is.list(x1))x=c(x1,x2)
}
if(is.matrix(x))x=listm(x)
crit=NA
n1=NA
n2=NA
if(is.matrix(x) || is.data.frame(x))K2=ncol(x)
if(is.list(x))K2=length(x)
K=floor(K2/2)
if(2*K!=K2)stop('Total number of groups, K2, should be an even number')
ic=0
ic2=K
pv=NULL
for(i in 1:K){
ic=ic+1
ic2=ic2+1
testit=yuen(x[[ic]],x[[ic2]],tr=tr,alpha=alpha)
n1[ic]=testit$n1
n2[ic]=testit$n2
pv[ic]=testit$p.value
}
pick=NULL
v=order(pv)
ic=0
for(i in K:1){
K2=2*K
flag=TRUE
if(opt==1){
i2=i*2
if(i==K)res=(0-2)*sum(log(pv))  # Fisher test statistic
if(i<K)res=(0-2)*sum(log(pv[-pick]))  # Fisher test statistic
pvF=1-pchisq(res,i2)   #Fisher p-value based on all tests.
}
if(opt==2){
if(i==K)res=sum(qnorm(pv/2)^2)  # C-N test
if(i<K)res=sum(qnorm(pv[-pick]/2)^2)
pvF=1-pchisq(res,i)
}
if(opt==3){
if(i==K)res=max(pv)
if(i<K)res=max(pv[-pick])
pvF=pbeta(res,i,1)
}
if(pvF>alpha)flag=TRUE
if(pvF<=alpha/(K+1-i)){
ic=ic+1
pick=c(pick,v[ic])
flag=FALSE
if(pv[v[ic]]>alpha)flag=TRUE
}
if(flag)break
}
Decision=rep('Not Sig',length(pv))
if(!is.null(pick))Decision[pick]='Reject'
nsig=sum(length(pick))
list(n1=n1,n2=n2,p.values=pv,
Decisions=as.matrix(Decision),num.sig=nsig)
}

twobicipv<-function(r1=sum(x),n1=length(x),r2=sum(y),n2=length(y),x=NA,y=NA,alpha=.05){
#
# Compute a p-value based on Beal's method for comparing two independent
# binomials.
#
alph=seq(.001,.999,.001)
for(i in 1:length(alph)){
pv=alph[i]
chk=twobici(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=y,alpha=alph[i])$ci #$
if(chk[1]>0 && chk[2]>0)break
if(chk[1]<0 && chk[2]<0)break
}
reg=twobici(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=y,alpha=alpha)
list(p.value=pv,ci=reg$ci,p1=reg$p1,p2=reg$p2)
}

twoDcorR_sub<-function(data,x,y,corfun=wincor,...){
#
# Used by TwoDcorR
#
rv=corfun(x[data,1],y[data],...)$cor
rv[2]=corfun(x[data,2],y[data],...)$cor
rv
}

twoDcorR<-function(x,y,corfun=wincor,alpha=.05,nboot=500,SEED=TRUE,MC=FALSE,outfun=outpro,...){
#
# Comparing two robust dependent correlations: Overlapping case
# Winsorized correlation is used by default.
#
# x is assumed to be a matrix with 2 columns
#
#  Compare correlation of x[,1] with y to x[,2] with y
#
#  The confidence interval is returned in ci
#  The estimates of the correlations are returned in est.rho1 and est.rho2
#
if(nrow(x)!=length(y))stop('x and y have different sample sizes; should be equal')
if(ncol(x)!=2)stop('Argument x should have two columns')
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1:2]
y=m1[,3]
est<-cor2xy(x,y,corfun=corfun,...)$cor
r12=est[1]
r13=est[2]
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
library(parallel)
bvec<-mclapply(data,twoDcorR_sub,x,y,corfun,...)
}
if(!MC)bvec<-lapply(data,twoDcorR_sub,x,y,corfun,...)
mat=matrix(NA,nrow=nboot,ncol=2)
for(i in 1:nboot)mat[i,]=bvec[[i]]
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
bsort<-sort(mat[,1]-mat[,2])
ci12<-1
ci12[1]<-bsort[ilow]
ci12[2]<-bsort[ihi]
pv=mean(bsort<0)+.5*mean(bsort==0)
pv=2*min(c(pv,1-pv))
list(est.rho1=r12,est.rho2=r13,ci=ci12,p.value=pv)
}

twoDNOV<-function(x,y,corfun=wincor,alpha=.05,nboot=500,SEED=TRUE,MC=FALSE){
#
# Comparing two robust dependent correlations: Non-overlapping case
# Winsorized correlation is used by default.
#
# Both x and y are assumed to be a matrix with 2 columns
#
#  Compare correlation of x[,1] x[,2] to the correlation between
#   y[,1] and  y[,2]
#
if(nrow(x)!=nrow(y))stop('x and y have different sample sizes; should be equal')
m1=cbind(x,y)
if(ncol(m1)!=4)stop('Both x and y should have two columns')
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1:2]
y=m1[,3:4]
r12=corfun(x[,1],x[,2])$cor
r13=corfun(y[,1],y[,2])$cor
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#
#  If you use corfun=scor, set plotit=F
#
data<-matrix(sample(nrow(y),size=nrow(y)*nboot,replace=TRUE),nrow=nboot)
data=listm(t(data))
if(MC){
library(parallel)
bvec1<-mclapply(data,corbsub,x[,1],x[,2],corfun)
bvec2<-mclapply(data,corbsub,y[,1],y[,2],corfun)
}
if(!MC){
bvec1<-lapply(data,corbsub,x[,1],x[,2],corfun)
bvec2<-lapply(data,corbsub,y[,1],y[,2],corfun)
}
mat1=matl(bvec1)
mat2=matl(bvec2)
ihi<-floor((1-alpha/2)*nboot+.5)
ilow<-floor((alpha/2)*nboot+.5)
bsort<-sort(mat1-mat2)
ci12<-bsort[ilow]
ci12[2]<-bsort[ihi]
ci12
pv=mean(bsort<0)
pv=2*min(c(pv,1-pv))
list(est.rho1=r12,est.rho2=r13,est.dif=r12-r13,ci=ci12,p.value=pv)
}

wmwpb<-function(x,y=NULL,est=median,alpha=.05,nboot=2000,SEED=TRUE,pr=TRUE,
na.rm=TRUE,...){
#
#   Compute a bootstrap confidence interval for a
#   measure of location associated with
#   the distribution of x-y,
#   est indicates which measure of location will be used
#   x and y are possibly dependent
#
#   loc2dif.ci  computes a non-bootstrap confidence interval
#
if(is.null(y[1])){
if(!is.matrix(x) & !is.data.frame(x))stop('With y missing, x should be a matrix')
y=x[,2]
x=x[,1]
}
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data1<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-NA
for(i in 1:nboot)bvec[i]<-wmwloc(x[data1[i,]],y[data2[i,]],est=est,na.rm=na.rm,...)
bvec<-sort(bvec)
low<-round((alpha/2)*nboot)+1
up<-nboot-low
temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
sig.level<-2*(min(temp,1-temp))
estdiff=wmwloc(x,y,est=est,na.rm=na.rm,...)
list(estimate=estdiff,ci=c(bvec[low],bvec[up]),p.value=sig.level)
}

TWOpovPV<-function(x,y,alpha=.05,CN=FALSE){
#
# Comparing two dependent correlations: Overlapping case
#
# x is assumed to be a matrix with 2 columns
#
#  Compare correlation of x[,1] with y to x[,2] with y
#  returns a confidence interval stored in
#  ci
#
# This function is exactly like TWOpov, only it returns a p-value as well.
#
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-TWOpov(x,y,alpha=alph[i],CN=CN,ZCI=TRUE)$ci
if(sign(chkit[1]*chkit[2])==1)break
}
p.value<-irem/100
if(p.value<=.1){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
irem=i
p.value<-alph[i]
chkit<-TWOpov(x,y,alpha=alph[i],CN=CN,ZCI=TRUE)$ci
if(sign(chkit[1]*chkit[2])==1)break
}}
if(p.value<=.1){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-TWOpov(x,y,alpha=alph[i],CN=CN,ZCI=TRUE)$ci
if(sign(chkit[1]*chkit[2])==1)break
}}
res=TWOpov(x,y,alpha=alpha,CN=CN)
list(p.value=p.value,est.rho1=res$est.rho1,est.rho2=res$est.rho2,ci=res$ci)
}

TWOpNOVPV<-function(x,y,HC4=TRUE,alpha=.05){
#
# Comparing two dependent correlations: Non-overlapping case
#
#   Compute a .95 confidence interval
#   for the difference between two dependent Pearson correlations,
#   non-overlapping case.
#
#    Both x and y are assumed to be matrices with two columns.
#   The function compares the correlation between x[,1] and x[,2]
#   to the correlation between y[,1] and y[,2].
#
#  For simulation results, see Wilcox (2009).
#  COMPARING PEARSON CORRELATIONS: DEALING WITH
#  HETEROSCEDASTICITY AND NON-NORMALITY, Communications in Statistics--Simulations
#   and Computations, 38, 2220-2234.
#
# This function is exactly like TWOpNOV, only it returns a p-value as well.
#
#  Note: To get a p-value, HC4=TRUE must be used.
#
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}
p.value<-irem/100
if(p.value<=.1){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
p.value<-alph[i]
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}}}
p.value<-irem/100
if(p.value<=.1){
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}}
if(p.value<=.001){
alph<-seq(.0001,.001,.0001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}}
if(p.value<=.001){
alph<-seq(.0001,.001,.0001)
for(i in 1:length(alph)){
p.value<-alph[i]
chkit<-TWOpNOV(x,y,alpha=alph[i],HC4=TRUE)
chkit=c(chkit$ci.lower,chkit$ci.upper)
if(sign(chkit[1]*chkit[2])==1)break
}}
res=TWOpNOV(x,y,alpha=alpha,HC4=TRUE)
ci=c(res$ci.lower,res$ci.upper)
list(p.value=p.value,est.1=res$est.1,est.2=res$est.2,ci=ci) #ci.lower=res$ci.lower,ci.upper=res$ci.upper)
}

twohc4cor<-function(x1,y1,x2,y2,alpha=.05){
#
#   Compare two independent Pearson correlations using the HC4 method
#
#
X<-elimna(cbind(x1,y1))
x1<-X[,1]
y1<-X[,2]
X<-elimna(cbind(x2,y2))
x2<-X[,1]
y2<-X[,2]
x1=(x1-mean(x1))/sd(x1)
y1=(y1-mean(y1))/sd(y1)
x2=(x2-mean(x2))/sd(x2)
y2=(y2-mean(y2))/sd(y2)
temp1=olshc4(x1,y1)
temp2=olshc4(x2,y2)
test=(temp1$ci[2,2]-temp2$ci[2,2])/sqrt(temp1$ci[2,6]^2+temp2$ci[2,6]^2)
df=length(x1)+length(x2)-4
pv=2*(1-pt(abs(test),df))
pv
}

funyuenpb<-function(x1,x2,tr=.2,pts=NULL,npts=25,plotit=TRUE,alpha=.05,
SEED=TRUE,
nboot=2000,xlab='T',ylab='Est.dif',FBP=TRUE,method='hochberg',COLOR=TRUE){
#
#  x1 and x2 are n-by-p matrices,
#  Designed for functional data.
#  For example, p measures taken over time where  p is typically large
#
#  Goal: at speficied times, compare the two groups.
#  pts: Can specify time points where comparisons are to be made
#  if pts=NULL, pick
#  npts points evenly space between min and max time points
#
p=ncol(x1)
pm1=p-1
if(p!=ncol(x2))stop('ncol(x1) is not equal to ncol(x2)')
n1=nrow(x1)
n2=nrow(x2)
if(SEED)set.seed(2)
if(is.null(pts)){
np=round(p/npts)
if(np==0)np=1
pts=seq(2,pm1,np)
notpts=-1*length(pts)
pts=pts[c(-1,notpts)]
}
npts=length(pts)
xsub1=x1[,pts]
xsub2=x2[,pts]
res=NA
dif=NA
bvals=matrix(nrow=nboot,ncol=npts)
for(j in 1:nboot){
data1=sample(n1,size=n1,replace=TRUE)
data2<-sample(n2,size=n2,replace=TRUE)
bvals[j,]=apply(xsub1[data1,],2,tmean,tr=tr)-apply(xsub2[data2,],2,tmean,tr=tr)
}
bsort=apply(bvals,2,sort)
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
op=matrix(NA,nrow=length(pts),ncol=7)
dimnames(op)=list(NULL,c('est 1','est 2','dif','p.value',
'adjust.p.value','ci.low','ci.hi'))
op[,1]=apply(xsub1,2,tmean,tr=tr)
op[,2]=apply(xsub2,2,tmean,tr=tr)
op[,3]=op[,1]-op[,2]
bsort=apply(bvals,2,sort)
bs=bvals<0
pv=apply(bs,2,mean)
pv2=rbind(pv,1-pv)
pv2=apply(pv2,2,min)
op[,4]=2*pv2
#flag0=op[,4]==0
#op[flag0,4]=.004
op[,5]=p.adjust(op[,4],method=method)
op[,6]=bsort[icl,]
op[,7]=bsort[icu,]
if(plotit){
if(!FBP){
xlow=c(1:nrow(op))
xax=rep(c(1:nrow(op)),3)
rplot(xlow,op[,3],xlab=xlab,ylab=ylab,scat=FALSE)
plot(xax,as.vector(op[,c(3,6,7)]),type='n',xlab=xlab,ylab=ylab)
lines(xlow,op[,3])
lines(xlow,op[,6],lty=2)
lines(xlow,op[,7],lty=2)
}
if(FBP){
par(mfrow=c(1,2))
if(COLOR)FBplot(x1)
if(!COLOR)func.out(x1)
lines(medcurve(x2))
if(COLOR)FBplot(x2)
if(!COLOR)func.out(x2)
lines(medcurve(x1))
par(mfrow=c(1,1))
}}
op=cbind(pts,op)
op
}

linWMW<-function(x,con,locfun=median,nreps=100,SEED=TRUE){
#
# Determine distribution of Y_i=sum_j c_jX_j
# Then estimate P(Y<0) and measure of location
# based on
# locfun, which defaults to the median.
#
con=as.vector(con)
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
if(length(con)!=J)stop('Length of con should equal number of groups')
x=elimna(x)
nv=as.vector(matl(lapply(x,FUN='length')))
nmin=min(nv)
est=NA
p=NA
B=list()
M=matrix(NA,nrow=nmin,ncol=J)
for(i in 1:nreps){
for(j in 1:J)M[,j]=sample(x[[j]],nmin)
B[[i]]=M
}
L=lapply(B,linWMWMC.sub,con=con)
est=lapply(L,locfun)
p=lapply(L,linWMWMC.sub2)
est=as.vector(matl(est))
p=as.vector(matl(p))
list(p=mean(p),center=mean(est))
}

interWMWpb<-function(x,nreps=100,SEED=TRUE,nboot=500,alpha=.05,nmax=10^8,MC=TRUE){
#
#
#
if(MC)library(parallel)
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
J<-length(x)
if(J!=4)stop('Number of groups should be four')
nv=lapply(x,length)
y=list()
pv=NA
N=max(pool.a.list(nv))
mat=matrix(NA,nrow=N,ncol=4)
for(i in 1:nboot){
for(j in 1:4)mat[1:nv[[j]],j]=sample(x[[j]],nv[[j]],replace=TRUE)
y[[i]]=mat
}
if(!MC)pv=lapply(y,interWMWpb.lsub)
if(MC)pv=mclapply(y,interWMWpb.lsub)
pv=pool.a.list(pv)
est=interWMW(x,nreps=nreps,SEED=SEED,nmax=nmax)
pv=sort(pv)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=pv[ilow]
ci[2]=pv[ihi]
pval=mean(pv<.5)+.5*mean(pv==.5)
pval=2*min(c(pval,1-pval))
list(p.est=est$p.est,ci=ci,p.value=pval,row.results=est$results.4.rows)
}

interWMWpb.lsub<-function(x,nreps=nreps){
v=interWMW(x,nreps=nreps,SEED=FALSE)$p.est
v
}

linWMWpb<-function(x,con,nreps=100,SEED=TRUE,nboot=500,alpha=.05,MC=FALSE){
#
# Compute a confidence interval for the probability that a linear contrast
# is less than zero.
#
con=as.vector(con)
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
J<-length(x)
if(J!=length(con))stop('Number of groups should be equal to the number of rows in con')
nv=lapply(x,length)
N=max(pool.a.list(nv))
mat=matrix(NA,nrow=N,ncol=J)
y=list()
pv=NA
est=linWMW(x,con=con,nreps=nreps,SEED=SEED)$p
for(i in 1:nboot){
#for(j in 1:J)y[[j]]=sample(x[[j]],nv[[j]],replace=TRUE)
for(j in 1:J)mat[1:nv[[j]],j]=sample(x[[j]],nv[[j]],replace=TRUE)
y[[i]]=mat
}
if(!MC)pv=lapply(y,linWMWpb.lsub,con=con,nreps=nreps,SEED=SEED)
if(MC){
library(parallel)
pv=mclapply(y,linWMWpb.lsub,con=con,nreps=nreps,SEED=SEED)
}
pv=pool.a.list(pv)
pv=sort(pv)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=pv[ilow]
ci[2]=pv[ihi]
pval=mean(pv<.5)+.5*mean(pv==.5)
pval=2*min(c(pval,1-pval))
list(p.est=est,ci=ci,p.value=pval)
}

linWMWpb.lsub<-function(x,nreps=nreps,con=con,SEED=SEED){
v=linWMW(x,nreps=nreps,con=con,SEED=SEED)$p
v
}

WMW2med<-function(x,y,q){
#
# If P(X<Y)=q, determine the value delta such that
#   P(X-delta-Y<=0.0)=q
#  So an estimate of delta yields an estimate of the qth quantile
#  of the sampling distribution of  the estimator p used by Cliff's (and the WMW) method.
#
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
n=max(n1,n2)
X=matrix(NA,nrow=n+1,ncol=2)
X[1:n1,1]=x
X[1:n2,2]=y
X[n+1,]=q
v=nelder(X,1,FN=WMW2med.sub,START=0)
v
}

WMW2med.sub<-function(X,delta){
n=nrow(X)
n1m=n-1
pv=cid(X[1:n1m,1]-delta,X[1:n1m,2])$phat
dif=abs(pv-X[n,1])
dif
}

twoKgen<-function(x=NULL,x1=NULL,x2=NULL,func=cidv2,alpha=.05,pr=TRUE,opt=1){
#
#  A step-down MCP based on K independent tests.
#  It is essential that the tests are independent.
#
#  Use Fisher method based on p-values coupled with Hochberg
#
# Data are assumed to be stored in two R variables, x1 and x2 or in one
#  R variable, x
#
# If stored in x1 and x2, they are assumed to be matrices with K columns
# or to have list mode, both having length K.
#
# If the data are stored in x,
# x is assumed to have 2K columns if a matrix or length 2K if it has list mode.
#
# If data are stored in x1 and x2, for each column, compute a p-value.
# That is, perform a test based on the data in column 1 of x1 and x2,
# followed by a test using the data in column 2 of x1 and x2, etc.
#
# If data are stored in x, the first test is based
# on the data in columns 1 and K+1,
# the second test is based on columns 2 and K+2, etc.
#
#  opt=1  Fisher's method
#  opt=2  Chen-Nadarajah method
#  opt=3  Max method
#
if(is.null(x[1])){
if(is.matrix(x1))x=cbind(x1,x2)
if(is.list(x1))x=c(x1,x2)
}
if(is.matrix(x))x=listm(x)
crit=NA
n1=NA
n2=NA
if(is.matrix(x) || is.data.frame(x))K2=ncol(x)
if(is.list(x))K2=length(x)
K=floor(K2/2)
if(2*K!=K2)stop('Total number of groups, K2, should be an even number')
ic=0
ic2=K
pv=NULL
for(i in 1:K){
ic=ic+1
ic2=ic2+1
testit=func(x[[ic]],x[[ic2]],alpha=alpha)
n1[ic]=testit$n1
n2[ic]=testit$n2
pv[ic]=testit$p.value
}
pick=NULL
v=order(pv)
ic=0
for(i in K:1){
K2=2*K
flag=TRUE
if(opt==1){
i2=i*2
if(i==K)res=(0-2)*sum(log(pv))  # Fisher test statistic
if(i<K)res=(0-2)*sum(log(pv[-pick]))  # Fisher test statistic
pvF=1-pchisq(res,i2)   #Fisher p-value based on all tests.
}
if(opt==2){
if(i==K)res=sum(qnorm(pv/2)^2)  # C-N test
if(i<K)res=sum(qnorm(pv[-pick]/2)^2)
pvF=1-pchisq(res,i)
}
if(opt==3){
if(i==K)res=max(pv)
if(i<K)res=max(pv[-pick])
pvF=pbeta(res,i,1)
}
if(pvF>alpha)flag=TRUE
if(pvF<=alpha/(K+1-i)){
ic=ic+1
pick=c(pick,v[ic])
flag=FALSE
if(pv[v[ic]]>alpha)flag=TRUE
}
if(flag)break
}
Decision=rep('Not Sig',length(pv))
if(!is.null(pick))Decision[pick]='Reject'
nsig=sum(length(pick))
list(n1=n1,n2=n2,p.values=pv,
Decisions=as.matrix(Decision),num.sig=nsig)
}

interWMW<-function(x,locfun=median,nreps=200,SEED=TRUE,nmax=10^8){
#
#  Goal: estimate P(X_1-X_2 < X_3-X_4).
#
# That is, dealing with an interaction in a 2-by-2 ANOVA design based on
# a Wilcoxon--Mann--Whitney approach but allow heteroscedasticity.
#
#  Strategy: estimate the distribution of X_1-X_2, non-parametrically  do the same
#  for X_3-X_4, then estimate P(X_1-X_2< X_3-X_4)
#
#  x should be a matrix with four columns or have list mode with length=4
#
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
J<-length(x)
if(J!=4)stop('Number of groups should be four')
#nx=pool.a.list(lapply(x,FUN='length'))
LL=list()
LL[[1]]=outer(x[[1]],x[[2]],FUN='-')
LL[[2]]=outer(x[[3]],x[[4]],FUN='-')
nv=c(length(LL[[1]]),length(LL[[2]]))
ntot=nv[1]*nv[2]
if(ntot<=nmax)p=bmp(LL[[1]],LL[[2]])$phat
else{
nmin=min(nv)
est=NA
p=NA
pest=NA
B=list()
M=matrix(NA,nrow=nmin,ncol=2)
for(i in 1:nreps){
for(j in 1:2)M[,j]=sample(LL[[j]],nmin)
B[[i]]=M
pest[i]=mean(M[,1]<M[,2])
}
L=lapply(B,linWMWMC.sub,con=c(1,-1))
est=lapply(L,locfun)
p=lapply(L,linWMWMC.sub2)
est=as.vector(matl(est))
p=as.vector(matl(p))
}
#
# NOTE:
# When computing a confidence interval
#L1=outer(x[[1]],x[[2]],FUN='-')
#L2=outer(x[[3]],x[[4]],FUN='-')
#bm=bmp(L1,L2)
# does not work due to the dependence among the values in L1 as well as L2
#
row=matrix(NA,nrow=2,ncol=6)
dimnames(row)=list(c('Row 1','Row 2'),c('n1','n2','p-value','p.hat','ci.low','ci.up'))
p1=cidv2(x[[1]],x[[2]])
p2=cidv2(x[[3]],x[[4]])
temp=matl(p1)
row[1,]=c(temp[1,c(1,2,5,6,7)],temp[2,7])
temp=matl(p2)
row[2,]=c(temp[1,c(1,2,5,6,7)],temp[2,7])
list(p.est=mean(p),results.4.rows=row)
}

WMWinter.est<-function(x,iter=10,SEED=TRUE){
#
# For a 2-by-2 design,
# Estimate P(Z<Z*)
# Z = X_1-X_2
# Z* =X_3-X_4
#
#  Compute a WMW type measure of effect	size.
#
ef=NA
if(is.matrix(x))x=listm(x)
nv=as.vector(matl(lapply(x,FUN='length')))
nt=prod(nv)
if(nt>10^3){
if(SEED)set.seed(2)
Nmin1=min(c(nv[1],nv[2],100))
Nmin2=min(c(nv[3],nv[4],100))
for(i in 1:iter){
id1=sample(nv[1],Nmin1)
id2=sample(nv[2],Nmin1)
L1=outer(x[[1]][id1],x[[2]][id2],FUN='-')
id1=sample(nv[3],Nmin2)
id2=sample(nv[4],Nmin2)
L2=outer(x[[3]][id1],x[[4]][id2],FUN='-')
ef[i]=pxly(L1,L2,iter=iter,SEED=SEED)
}}
if(nt<=10^3){
L1=outer(x[[1]],x[[2]],FUN='-')
L2=outer(x[[3]],x[[4]],FUN='-')
ef=pxly(L1,L2,iter=iter,SEED=SEED)
}
ef=mean(ef)
ef
}

interWMWAP<-function(x,nreps=100,SEED=TRUE,nboot=500,alpha=.05,nmax=10^8,MC=TRUE){
#
#  Interaction in a 2-by-2 design using P(X_1-X_2<X_3-X_4)
#  Compared to WMWinterci, this method is  better at dealing with
#  small unequal sample sizes when there is heteroscedasticity.
#
#
if(MC)library(parallel)
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
J<-length(x)
if(J!=4)stop('Number of groups should be four')
TV=linWMWpb(x,con=c(1,-1,-1,1),nreps=nreps,SEED=SEED,nboot=nboot,
alpha = alpha, MC = MC)
TV
}

WMWinterci<-function(x,alpha=0.05,SW=FALSE){
#
#
# Jan De Neve & Olivier Thas (2016): A Mann--Whitney type effect measure
# of interaction for factorial designs, Communications in Statistics - Theory and Methods, DOI:
#10.1080/03610926.2016.1263739
#
# For a 2-by-2 design,
# Estimate P(Z<Z^*)
# Z = X_1-X-2
# Z^* =X_3-X_4
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
x=elimna(x)
if(SW)x=x[c(1,3,2,4)]
J<-length(x)
if(J!=4)stop('Number of groups should be four')
bhat=WMWinter.est(x)
gb=qnorm(bhat)
nx=pool.a.list(lapply(x,FUN='length'))
N=prod(nx)
I1=0
I2=0
I3=0
I4=0
for(i in 1:length(x[[1]]))I1=I1+Ifun(x[[1]][1],x[[2]],x[[3]],x[[4]],bhat)^2
for(i in 1:length(x[[2]]))I2=I2+Ifun(x[[1]],x[[2]][i],x[[3]],x[[4]],bhat)^2
for(i in 1:length(x[[3]]))I3=I3+Ifun(x[[1]],x[[2]],x[[3]][i],x[[4]],bhat)^2
for(i in 1:length(x[[4]]))I4=I4+Ifun(x[[1]],x[[2]],x[[3]],x[[4]][i],bhat)^2
deriv=1/dnorm(qnorm(bhat))
sighat=(deriv/N)^2*(I1+I2+I3+I4)
crit=qnorm(1-alpha/2)
ci=gb-crit*sqrt(sighat)
ci[2]=gb+crit*sqrt(sighat)
ci=pnorm(ci)
test=(gb-0.)/sqrt(sighat)
pv=2*(1-pnorm(abs(test)))
list(n=nx,p.hat=bhat,p.value=pv,ci=ci)
}

yuenQS<-function(x,y=NULL,tr=.2,alpha=.05, plotit=FALSE,op=TRUE,
cor.op=FALSE,loc.fun=median,pr=TRUE,xlab='X',ylab=' ' ){
#
#  Perform Yuen's test for trimmed means on the data in x and y.
#  The default amount of trimming is 20%
#  Missing values (values stored as NA) are automatically removed.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuen$ci.
#  The significance level is returned in yuen$p.value
#
#   Unlike the function yuen, a robust quantile shift measure
#   of effect size is returned.
#
if(pr){
print('Note: Under normality and homoscedasticity, Cohen d= 0, .2, .5, .8')
print('correspond approximately to Q.Effect = 0.5, 0.55, 0.65 and 0.70, respectively')
}
if(tr==.5)stop('Use medpb to compare medians.')
if(tr>.5)stop('cannot have tr>.5')
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
library(MASS)
x<-x[!is.na(x)]  # Remove any missing values in x
y<-y[!is.na(y)]  # Remove any missing values in y
n1=length(x)
n2=length(y)
h1<-length(x)-2*floor(tr*length(x))
h2<-length(y)-2*floor(tr*length(y))
q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
crit<-qt(1-alpha/2,df)
m1=mean(x,tr)
m2=mean(y,tr)
mbar=(m1+m2)/2
dif=m1-m2
low<-dif-crit*sqrt(q1+q2)
up<-dif+crit*sqrt(q1+q2)
test<-abs(dif/sqrt(q1+q2))
yuen<-2*(1-pt(test,df))
e.pow=shiftQS(x,y,tmean,tr=tr)$Q.Effect
if(plotit){
g2plot(x,y,xlab=xlab,ylab=ylab)
}
list(ci=c(low,up),n1=n1,n2=n2,
p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,
crit=crit,df=df,Q.Effect=e.pow)
}

ancdetwmw<-function(x1,y1,x2,y2,fr1=1,fr2=1,nmin=8,Ycrit=FALSE,
alpha=.05,plotit=TRUE,pts=NA,span=2/3,sm=TRUE,BOTH=TRUE,
pr=TRUE,xout=FALSE,outfun=out,MC=FALSE,
npts=25,p.crit=NULL,EST=FALSE,
SCAT=TRUE,xlab='X',ylab='P.hat',pc='.',...){
#
#  Like the function ancdet, only use analog of Wilcoxon--Mann--Whitney
#  plot=TRUE:  plot  estimates P.hat plus a
# confidence band having simultaneous probability coverage 1-alpha
#
#   Method S: choose covariate points based on nmin.
#  In contrast method Q, performed by ancdetwmwQ, picks points between the q and 1-q quantiles.
#
#
#  span = the span when using loess to plot the regression line.
#
# npts = number of  covariate values to be used
#
# sm=TRUE will smooth the plot using lowess
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
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
xor1=order(x1)
xor2=order(x2)
x1=x1[xor1]
x2=x2[xor2]
y1=y1[xor1]
y2=y2[xor2]
n1<-1
n2<-1
vecn<-1
for(i in 1:length(x1))n1[i]<-length(y1[near(x1,x1[i],fr1)])
for(i in 1:length(x1))n2[i]<-length(y2[near(x2,x1[i],fr2)])
for(i in 1:length(x1))vecn[i]<-min(n1[i],n2[i])
sub<-c(1:length(x1))
bot<-min(sub[vecn>=nmin])
itop<-max(sub[vecn>=nmin])
xbot=x1[bot]
xup=x1[itop]

if(BOTH){
vecn=1
n1=1
n2=1
for(i in 1:length(x2))n1[i]<-length(y1[near(x1,x1[i],fr1)])
for(i in 1:length(x2))n2[i]<-length(y2[near(x2,x1[i],fr2)])
for(i in 1:length(x2))vecn[i]<-min(n1[i],n2[i])
sub<-c(1:length(x2))
bot<-max(sub[vecn>=nmin])
itop<-min(sub[vecn>=nmin])
xbot[2]=x2[itop]  #CORRECT, need to switch
xup[2]=x2[bot]
}
xbot=max(xbot)
xup=min(xup)
pts=seq(xbot,xup,length.out=npts)
if(alpha!=.05)EST=TRUE
if(is.null(p.crit)){
nv=c(30,  50,  60,  70,  80, 100, 150, 200,
300, 400, 500, 600, 800)
if(Ycrit)pv=c(0.00824497,0.00581, 0.005435089, 0.004763079,
0.00416832, 0.004406774, 0.00388228,0.003812836,0.003812836,0.003453055, 0.003625061,
.003372966, 0.003350022)
if(!Ycrit)
pv=c(0.008566, # 30
0.0083847, # 50
0.006758,  # 60
0.006871,   # 70
0.006157,  # 80
0.006629, #100
0.006629, #  150
0.004681, # 200
0.004537,  # 300
0.004952, # 400
 0.004294,  # 500
 0.004288,  # 600
 0.004148)
n1= length(y1)
 n2=length(y2)
p.crit=(lplot.pred(1/nv,pv,1/n1)$yhat+lplot.pred(1/nv,pv,1/n2)$yhat)/2
p.crit=(alpha/.05)*p.crit # Crude approximation when alpha != .05, tends to be conservative.
}
temp=ancovaWMW(x1,y1,x2,y2,pts=pts,fr1=fr1,fr2=fr2,alpha=p.crit,plotit=plotit)
res=temp$output[,1:7]
if(plotit){
x=res[,1]
y=res[,4]
minx=min(x)
maxx=max(x)
plot(c(minx,maxx,x),c(0,1,y),xlab=xlab,ylab=ylab,type='n')
points(x,y,pch=pc)
if(!sm){lines(res[,1],res[,5],lty=2)
lines(res[,1],res[,6],lty=2)
lines(res[,1],res[,4])
}
if(sm){
plin=lplot.pred(res[,1],res[,4],span=span)$yhat
lines(res[,1],plin)
low.line=lplot.pred(res[,1],res[,5],span=span)$yhat
lines(res[,1],low.line,lty=2)
up.line=lplot.pred(res[,1],res[,6],span=span)$yhat
lines(res[,1],up.line,lty=2)
}

}
sig=rep(0,nrow(res))
sig[res[,7]<=p.crit]=1
sig=as.matrix(sig,ncol=1)
dimnames(sig)=list(NULL,'Sig.Dif')
res=cbind(res,sig)
list(p.crit=p.crit,output=res,summary=temp$summary,num.sig=sum(sig))
}

ancdetwmwQ<-function(x1,y1,x2,y2,fr1=1,fr2=1,nmin=8,q=.05,
alpha=.05,plotit=TRUE,pts=NA,span=2/3,sm=TRUE, xout=FALSE,outfun=out,MC=FALSE,
npts=25,p.crit=NULL,
SCAT=TRUE,xlab='X',ylab='P.hat',pc='.',...){
#
#  Like the function ancdet, only use analog of Wilcoxon--Mann--Whitney
#  plot=TRUE:  plot  estimates P.hat plus a
# confidence band having simultaneous probability coverage 1-alpha
#
#  span = the span when using loess to plot the regression line.
#
# npts = number of  covariate values to be used
#
# sm=TRUE will smooth the plot using lowess
#
#  Covariate points are chosen that lie between the q and 1-q quantiles
#
if(ncol(as.matrix(x1))>1)stop('One covariate only is allowed with this function')
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
xor1=order(x1)
xor2=order(x2)
x1=x1[xor1]
x2=x2[xor2]
y1=y1[xor1]
y2=y2[xor2]
n1<-1
n2<-1
vecn<-1

xbot=max(qest(x1,q),qest(x2,q))
xup=min(qest(x1,1-q),qest(x2,1-q))

pts=seq(xbot,xup,length.out=npts)

nchk1=0
for(i in 1:length(pts))nchk1[i]=length(y1[near(x1,pts[i],fr1)])
nchk2=0
for(i in 1:length(pts))nchk2[i]=length(y2[near(x2,pts[i],fr2)])
flag1=nchk1>=nmin
flag2=nchk2>=nmin
flag=as.logical(flag1*flag2)
pts=pts[flag]
if(is.null(p.crit)){
nv=c(30,  50,  60,  70,  80, 100, 150, 200,
300, 400, 500, 600, 800)
pv=c(0.00824497,0.00581, 0.005435089, 0.004763079,
0.00416832, 0.004406774, 0.00388228,0.003812836,0.003812836,0.003453055, 0.003625061,
.003372966, 0.003350022)
n1= length(y1)
 n2=length(y2)
p.crit=(lplot.pred(1/nv,pv,1/n1)$yhat+lplot.pred(1/nv,pv,1/n2)$yhat)/2
p.crit=(alpha/.05)*p.crit # Crude approximation when alpha != .05, tends to be conservative.
}
temp=ancovaWMW(x1,y1,x2,y2,pts=pts,fr1=fr1,fr2=fr2,alpha=p.crit,plotit=plotit)
res=temp$output
if(plotit){
x=res[,1]
y=res[,4]
minx=min(x)
maxx=max(x)
plot(c(minx,maxx,x),c(0,1,y),xlab=xlab,ylab=ylab,type='n')
points(x,y,pch=pc)
if(!sm){lines(res[,1],res[,5],lty=2)
lines(res[,1],res[,6],lty=2)
lines(res[,1],res[,4])
}
if(sm){
plin=lplot.pred(res[,1],res[,4],span=span)$yhat
lines(res[,1],plin)
low.line=lplot.pred(res[,1],res[,5],span=span)$yhat
lines(res[,1],low.line,lty=2)
up.line=lplot.pred(res[,1],res[,6],span=span)$yhat
lines(res[,1],up.line,lty=2)
}

}
sig=rep(0,nrow(res))
sig[res[,7]<=p.crit]=1
sig=as.matrix(sig,ncol=1)
dimnames(sig)=list(NULL,'Sig.Dif')
res=cbind(res,sig)
list(p.crit=p.crit,output=res,summary=temp$summary,num.sig=sum(sig),p.crit=p.crit)
}

mulwmw.dist.new<-function(m1,m2,new,cop=3){
#
#
# Determine center corresponding to two
# independent groups, project all  points onto line
# connecting the centers based on m1 and m2. Return projected distances for m1
# m2.
#  new:  new data, not known whether it came from group 1 or 2.
#  This function is used in pro.class, a classification method.
#
#
#  There are three options for computing the center of the
#  cloud of points when computing projections:
#  cop=1 uses Donoho-Gasko median
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#
#  When using cop=2 or 3, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
if(is.null(dim(m1))||dim(m1)[2]<2){stop("m1 and m2 should have two or more columns")
}
m1<-elimna(m1) # Remove missing values
m2<-elimna(m2)
new=elimna(new)
FLAG=FALSE
new=as.matrix(new)
if(ncol(new)==1){
FLAG=TRUE
new=t(new)  # If test is a vector, a single point, transpose to get correct number of columns.
new=rbind(new,new)   #avoid R from aborting.
}
n1=nrow(m1)
n2=nrow(m2)
if(cop==1){
if(ncol(m1)>2){
center1<-dmean(m1,tr=.5)
center2<-dmean(m2,tr=.5)
}
if(ncol(m1)==2){
tempd<-NA
for(i in 1:nrow(m1))
tempd[i]<-depth(m1[i,1],m1[i,2],m1)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center1<-m1[flag,]
if(sum(flag)>1)center1<-apply(m1[flag,],2,mean)
for(i in 1:nrow(m2))
tempd[i]<-depth(m2[i,1],m2[i,2],m2)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center2<-m2[flag,]
if(sum(flag)>1)center2<-apply(m2[flag,],2,mean)
}}
if(cop==2){
center1<-cov.mcd(m1)$center
center2<-cov.mcd(m2)$center
}
if(cop==3){
center1<-apply(m1,2,median)
center2<-apply(m2,2,median)
}
if(cop==4){
center1<-smean(m1)
center2<-smean(m2)
}
center<-(center1+center2)/2
B<-center1-center2
if(sum(center1^2)<sum(center2^2))B<-(0-1)*B
BB<-B^2
bot<-sum(BB)
disx<-NA
disy<-NA
dis.new=NA
if(bot!=0){
for (j in 1:nrow(m1)){
AX<-m1[j,]-center
tempx<-sum(AX*B)*B/bot
disx[j]<-sign(sum(AX*B))*sqrt(sum(tempx^2))
}
for (j in 1:nrow(m2)){
AY<-m2[j,]-center
tempy<-sum(AY*B)*B/bot
disy[j]<-sign(sum(AY*B))*sqrt(sum(tempy^2))
}
for (j in 1:nrow(new)){
AY<-new[j,]-center
tempy<-sum(AY*B)*B/bot
dis.new[j]<-sign(sum(AY*B))*sqrt(sum(tempy^2))
}}
if(FLAG)dis.new=dis.new[1]
list(dist1=disx,dist2=disy,dis.new=dis.new)
}

mulwmw.dist.new<-function(m1,m2,new,cop=3){
#
#
# Determine center corresponding to two
# independent groups, project all  points onto line
# connecting the centers based on m1 and m2. Return projected distances for m1
# m2.
#  new:  new data, not known whether it came from group 1 or 2.
#  This function is used in pro.class, a classification method.
#
#
#  There are three options for computing the center of the
#  cloud of points when computing projections:
#  cop=1 uses Donoho-Gasko median
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#
#  When using cop=2 or 3, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
if(is.null(dim(m1))||dim(m1)[2]<2){stop("m1 and m2 should have two or more columns")
}
m1<-elimna(m1) # Remove missing values
m2<-elimna(m2)
new=elimna(new)
FLAG=FALSE
new=as.matrix(new)
if(ncol(new)==1){
FLAG=TRUE
new=t(new)  # If test is a vector, a single point, transpose to get correct number of columns.
new=rbind(new,new)   #avoid R from aborting.
}
n1=nrow(m1)
n2=nrow(m2)
if(cop==1){
if(ncol(m1)>2){
center1<-dmean(m1,tr=.5)
center2<-dmean(m2,tr=.5)
}
if(ncol(m1)==2){
tempd<-NA
for(i in 1:nrow(m1))
tempd[i]<-depth(m1[i,1],m1[i,2],m1)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center1<-m1[flag,]
if(sum(flag)>1)center1<-apply(m1[flag,],2,mean)
for(i in 1:nrow(m2))
tempd[i]<-depth(m2[i,1],m2[i,2],m2)
mdep<-max(tempd)
flag<-(tempd==mdep)
if(sum(flag)==1)center2<-m2[flag,]
if(sum(flag)>1)center2<-apply(m2[flag,],2,mean)
}}
if(cop==2){
center1<-cov.mcd(m1)$center
center2<-cov.mcd(m2)$center
}
if(cop==3){
center1<-apply(m1,2,median)
center2<-apply(m2,2,median)
}
if(cop==4){
center1<-smean(m1)
center2<-smean(m2)
}
center<-(center1+center2)/2
B<-center1-center2
if(sum(center1^2)<sum(center2^2))B<-(0-1)*B
BB<-B^2
bot<-sum(BB)
disx<-NA
disy<-NA
dis.new=NA
if(bot!=0){
for (j in 1:nrow(m1)){
AX<-m1[j,]-center
tempx<-sum(AX*B)*B/bot
disx[j]<-sign(sum(AX*B))*sqrt(sum(tempx^2))
}
for (j in 1:nrow(m2)){
AY<-m2[j,]-center
tempy<-sum(AY*B)*B/bot
disy[j]<-sign(sum(AY*B))*sqrt(sum(tempy^2))
}
for (j in 1:nrow(new)){
AY<-new[j,]-center
tempy<-sum(AY*B)*B/bot
dis.new[j]<-sign(sum(AY*B))*sqrt(sum(tempy^2))
}}
if(FLAG)dis.new=dis.new[1]
list(dist1=disx,dist2=disy,dis.new=dis.new)
}

cidMULT<-function(x1,x2,alpha=.05,BMP=FALSE){
#
#
#  Deals with  two independent p-variate random variables
# for each variable, compare the two groups with Cliff's analog of the
# Wilcoxon--Mann--Whitney test.
#
p=ncol(x1)
V=list()
res=matrix(NA,nrow=p,ncol=9)
for(j in 1:p){
if(!BMP)a=cidv2(x1[,j],x2[,j])
if(BMP)a=bmp(x1[,j],x2[,j])
res[j,1]=a$n1
res[j,2]=a$n2
res[j,3]=a$p.value
if(!BMP)res[j,4]=a$p.hat
if(BMP)res[j,4]=a$phat
if(!BMP)res[j,5:6]=a$p.ci
#if(BMP)res[j,5:6]=a$ci.p
res[j,7:9]=a$summary.dvals
}
L=NULL
for(j in 1:p)L[j]=paste('Var',j)
dimnames(res)=list(L,c('n1','n2','p-value','p.hat','ci.low','ci.up','P(X<Y)','P(X=Y)','P(X>Y)'))
res
}

twoway.pool<-function(J,K,x){
#
#  For a two-way design,for each level of Factor A, pool over B
#  Do the same for Factor B
#
#  Return results in A and B having list mode,
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
JK=J*K
imat=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
A=list()
B=list()
for(j in 1:J)A[[j]]=pool.a.list(x[imat[j,]])
for(k in 1:K)B[[k]]=pool.a.list(x[imat[,k]])
list(A=A,B=B)
}

twowayESM<-function(J,K,x,fun=ES.summary,...){
#
# For each level of Factor A, pool data over levels of Factor B, compute
# measures of effect size.
#
# Do the same for Factor B
#
#  argument fun, see the function IND.PAIR.ES
#
a=twoway.pool(J,K,x)
A=IND.PAIR.ES(a$A,fun=fun,...)
B=IND.PAIR.ES(a$B,fun=fun,...)
list(Factor.A=A,Factor.B=B)
}

RMcomvar.locdis<-function(x,y,
loc.fun=median,CI=FALSE,plotit=TRUE,xlab='First Group',
ylab='Est.1 - Est.2',ylabQCOM='Est.2 - Est.1',sm=TRUE,QCOM=TRUE,q=c(.1,.25,.75,.9),MC=FALSE,nboot=2000,PR=TRUE,...){
#
#  Compare the marginal distributions of two dependent groups in terms of the
#  variation in the tails using all of the quantiles
#  after centering the data.
#
#  CI=FALSE, suppresses confidence intervals
#
if(!QCOM){
if(PR){
print('Interpretation: when using QCOM=F:  If values in  q.sig.greater are less than .5')
print('this indicates more variation in the lower tail for group 1')
print('Interpretation: If values in  q.sig.greater are greater than .5')
print('This indicates more variation in the lower tail for group 2')

print('Interpretation: If values in  q.sig.less are less than .5')
print('this indicates more variation in the upper tail for group 2')
print('Interpretation: If values in  q.sig.less are greater than .5')
print('This indicates more variation in the upper tail for group 1')
}
}
x=elimna(x)
y=elimna(y)
mx=loc.fun(x,...)
my=loc.fun(y,...)
X=x-mx
Y=y-my
if(!QCOM){
a=lband(X,Y,plotit=plotit,xlab=xlab,ylab=ylabQCOM,sm=sm,CI=CI)
if(!CI)a$m=NULL
}
else{
a=Dqcomhd(X,Y,q=q,nboot=nboot,plotit=plotit,xlab=xlab,ylab=ylab)
}
a
}

yuend.QS.ci<-function(x,y=NULL,tr=.2,alpha=.05,nboot=1000,SEED=TRUE,...){
#
# confidence interval for the quantile shift measure of effect size.
#
if(is.null(y)){
if(ncol(x)!=2)stop('If y is null, x should have	two columns')
}
if(SEED)set.seed(2)
xy=cbind(x,y)
xy=elimna(xy)
n1=nrow(xy)
v=NA
for(i in 1:nboot){
id=sample(n1,replace=TRUE)
v[i]=shiftes(xy[id,1],xy[id,2],locfun=tmean,tr=tr)$Q.Effect
}
v=sort(v)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=v[ilow]
ci[2]=v[ihi]
est=shiftes(xy[,1],xy[,2],locfun=tmean,tr=tr)$Q.Effect
list(Q.effect=est,ci=ci)
}

comvar.mcp<-function(x,method='hoch',SEED=TRUE){
#
#  Compare the variances of J indepenent variables.
#  Perform all pairwise comparisons using
#   a slight extension of  HC4 vesion of the Morgan-Pitman test
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
J=length(x)
CC=(J^2-J)/2
output<-matrix(0,CC,7)
dimnames(output)<-list(NULL,c('Var','Var','SD 1','SD 2','Dif','p.value','Adj.p.value'))
ic=0
for (j in 1:J){
for (k in 1:J){
if (j < k){
ic=ic+1
a=varcom.IND.MP(x[[j]],x[[k]],SEED=SEED)
a=pool.a.list(a)
output[ic,1]=j
output[ic,2]=k
output[ic,3:4]=a[1:2]
output[ic,3:4]=sqrt(output[ic,3:4])
output[ic,5]=output[ic,3]-output[ic,4]
output[ic,6]=a[3]
}}}
output[,7]=p.adjust(output[,6],method=method)
output
}

oph.ind.comvar<-function(x,y=NULL,method='hommel',invalid=4,SEED=TRUE,STOP=TRUE){
#
#  This function is designed specifically for dealing with
#  Prediction Error for Intraocular Lens Power Calculation
#  It is assumed that any value less than -3 diopters or greater than 3 diopters
#  is invalid.  The argument invalid can be used to change this decision rule.
#
#  Goal: compare the variances of J independent measures.
#  All pairwise comparisons are performed using
#   a slight extension of the  HC4 vesion of the Morgan-Pitman test
#
#   x can be a matrix, a data frame or it can have list mode.
#   if y is not NULL, the function assumes x is a vector
#   and the goal is to compare the variances of the data in x and y.
#
#  By default, Hochberg's method is used to control the probability of one
#  or more TypeI errors
#
if(!is.null(y))x=list(x,y)
if(is.matrix(x) || is.data.frame(x))x=listm(x)
J=length(x)
for(j in 1:J)x[[j]]=elimna(x[[j]])
for(j in 1:J){
flag=abs(elimna(x[[j]]))>invalid
if(sum(flag,na.rm=TRUE)>0){
print(paste('The value of argument invalid indicates that any value greater in absolute value than', invalid,' is invalid'))
print(paste('Variable', j, 'has one  or more invalid values'))
print('They occur in the following positions')
nr=c(1:length(x[[j]]))
print(nr[flag])
if(STOP)stop()
}
}
CC=(J^2-J)/2
output<-matrix(0,CC,7)
dimnames(output)<-list(NULL,c('Var','Var','SD 1','SD 2','Dif','p.value','Adj.p.value'))
ic=0
for (j in 1:J){
for (k in 1:J){
if (j < k){
ic=ic+1
a=varcom.IND.MP(x[[j]],x[[k]],SEED=SEED)
a=pool.a.list(a)
output[ic,1]=j
output[ic,2]=k
output[ic,3:4]=a[1:2]
output[ic,3:4]=sqrt(output[ic,3:4])
output[ic,5]=output[ic,3]-output[ic,4]
output[ic,6]=a[3]
}}}
output[,7]=p.adjust(output[,6],method=method)
output
}

oph.dep.comvar<-function(x, y=NULL, invalid=4, method='hommel',STOP=TRUE){
#
#  This function is designed specifically for dealing with
#  Prediction Error for Intraocular Lens Power Calculation
#  It is assumed that any value less than -3 diopters or greater than 3 diopters
#  is invalid.  The argument invalid can be used to change this decision rule.
#
#  Goal: compare the variances of J dependent measures.
#  All pairwise comparisons are performed using
#   a slight extension of the  HC4 version of the Morgan-Pitman test
#  Compare the variances of J dependent variables.
#  Perform all pairwise comparisons using the HC4 extension of the Morgan-Pitman test

#   x can be a matrix, a data frame or it can have list mode.
#   if y is not NULL, the function assumes x is a vector
#   and the goal is to compare the variances of the data in x and y.
#
#  By default, Hochberg's method is used to control the probability of one
#  or more Type I errors
#
if(!is.null(y))x=cbind(x,y)
if(is.list(x)){
n=pool.a.list(lapply(x,length))
if(var(n)!=0)stop('lengths have different values')
x=matl(x)
}
J=ncol(x)
flag=abs(elimna(x))>invalid
if(sum(flag,na.rm=TRUE)>0){
nr=c(1:nrow(x))
if(sum(flag)>1){
print(paste('The value of argument invalid indicates that any value greater in absolute value than', invalid,' is invalid'))
print('The following rows have invalid values')
}
if(sum(flag)==1){
print(paste('The value of argument invalid indicates that any value greater in absolute value than', invalid,' is invalid'))
print('The following row has an  invalid value')
}
irow=NA
ic=0
N=nrow(x)
for(i  in 1:N){
iflag=abs(x[i,])>invalid
if(sum(iflag,na.rm=TRUE)>0){
ic=ic+1
irow[ic]=i
}}
print(irow)
if(STOP)stop()
}
CC=(J^2-J)/2
output<-matrix(0,CC,7)
dimnames(output)<-list(NULL,c('Var','Var','SD 1','SD 2','Dif','p.value','Adj.p.value'))
ic=0
for (j in 1:J){
for (k in 1:J){
if (j < k){
ic=ic+1
a=comdvar(x[,j],x[,k])
output[ic,1]=j
output[ic,2]=k
output[ic,3]=a$est1
output[ic,4]=a$est2
output[ic,3]=sqrt(a$est1)
output[ic,4]=sqrt(a$est2)
output[ic,5]=sqrt(a$est1)- sqrt(a$est2)
output[ic,6]=a$p.value
}}}
output[,7]=p.adjust(output[,6],method=method)
output
}

two.dep.pb<-function(x,y=NULL,alpha=.05,est=tmean,plotit=FALSE,dif=TRUE,
nboot=NA,xlab='Group 1',ylab='Group 2',pr=TRUE,SEED=TRUE,...){
#
#  Two dependent groups
# Compare measures of location  via a percentile bootstrap.
#  Trimmed mean	used by	default.
#
#  Calls rmmcppb, provided for convenience
#  nboot, number of bootstrap samples defaults to 1000
#
if(pr){
if(dif)print('dif=TRUE, difference scores were used')
if(!dif)print('dif=FALSE, marginal trimmed means  were used')
}
if(is.null(y)){
if(ncol(x)!=2)stop('y is null  so x should have two columns')
}
if(!is.null(y)){
xy=cbind(x,y)
xy=elimna(xy)
x=xy[,1]
y=xy[,2]
}
e=apply(cbind(x,y),2,est,...)
a=rmmcppb(x,y,est=est,nboot=nboot,alpha=alpha,SR=FALSE,SEED=SEED,
plotit=plotit,dif=dif,BA=FALSE,pr=FALSE,...)$output
if(!dif){
output=matrix(c(e[1],e[2],a[1,2],a[1,3],a[1,5],a[1,6]),nrow=1)
dimnames(output)=list(NULL,c('Est.1','Est.2','Est.dif','p.value','ci.lower','ci.upper'))
}
if(dif){
output=matrix(c(a[1,2],a[1,3],a[1,5],a[1,6]),nrow=1)
dimnames(output)=list(NULL,c('Est.typical.dif','p.value','ci.lower','ci.upper'))
}
output
}

wmw.est.only<-function(m,M=m,n1,n2){
#
# For two independent random variables,
# estimate P(X<Y)
#
# m = matrix, ncol=2
#
if(sum(is.na(m[,1])!=n1))x1=sample(M[,1])
else
x1=sample(m[,1])
if(sum(is.na(m[,2])!=n1))x2=sample(M[,2])
else
x2=sample(m[,2])
e=bmp(x1,x2)$phat
e
}

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

wmw.RZR<-function(x,y,nboot=1000,SEED=TRUE){
#
#   Perform the Reiczigel et al. (2005) improvement of of the
# Wilcoxon--Mann--Whitney  test
#
if(SEED)set.seed(2)
val=0
n1=length(x)
n2=length(y)
xy=rank(c(x,y))
N=n1+n2
n1p1=n1+1
a=yuen(xy[1:n1],xy[n1p1:N],tr=0)$teststat
#print(a)
LOC=loc2dif(x,y)
x=x-a
y=y-a
#print(yuen(x,x2,tr=0))
bval=0
for(i in 1:nboot){
z1=sample(x,n1,replace=TRUE)
z2=sample(y,n2,replace=TRUE)
XY=rank(c(z1,z2))
bval[i]=yuen(XY[1:n1],XY[n1p1:N],tr=0)$teststat
}
#print(bval[1:10])
pv1=mean(a>bval)
pv2=mean(a<bval)

pv=2*min(pv1,pv2)
list(p.value=pv)
}

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

comvar2.astig<-function(m1,m2,method='holm'){
#
#  This function is designed specifically for dealing with astigmatism.
#
#  m= matrix or data frame, four columns
#
# compare variances of m[,1] vs m[,3] as well as m[,2] vs m[,4]
#
output=matrix(NA,2,5)
a=varcom.IND.MP(m1[,1], m2[,1])
output[1,1]=a$est1
output[1,2]=a$est2
output[1,3]=a$est1/a$est2
output[1,4]=a$p.value
a=varcom.IND.MP(m1[,2], m2[,2])
output[2,1]=a$est1
output[2,2]=a$est2
output[2,3]=a$est1/a$est2
output[2,4]=a$p.value
output[,5]=p.adjust(output[,4],method=method)
dimnames(output)=list(NULL,c('VAR 1','VAR 2','Ratio','p.value','p.adjusted') )
output
}

wmw.ancbse<-function(x1,y1,x2,y2,pts,nboot=100,SEED=TRUE,MC=FALSE,null.value=.5,xout=FALSE,
outfun=outpro,alpha=.05,...){
#
#  ANCOVA based on a heteroscedastic analog of the
#  Wilcoxon--Mann-Whitney test
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

if(is.null(pts))pts=ancova.KMS(x1,y1,x2,y2,plotit=FALSE)[,1]
if(SEED)set.seed(2)
e=wmw.anc(x1,y1,x2,y2,pts=pts)
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
if(!MC)tv=lapply(bl,wmw.ancbse.sub,pts)
else{
library(parallel)
tv=mclapply(bl,wmw.ancbse.sub,pts)
}
E=t(matl(tv))
se=apply(E,2,sd)
ci[,1]<-e-qnorm(1-alpha/2)*se
ci[,2]<-e+qnorm(1-alpha/2)*se
test<-(e-null.value)/se
sig<-2*(1-pnorm(abs(test)))
list(Est=e,SE=se,test.stat=test,conf.int=ci,p.value=sig)
}

wmw.ancbse.sub<-function(m,pts){
v=wmw.anc(m[[1]],m[[2]],m[[3]],m[[4]],pts=pts)
v
}

wmw.anc<-function(x1,y1,x2,y2,pts,xout=FALSE,outfun=outpro){
#
#  Estimate P(y1<y2)  given that the covariate X is equal to the values in pts.
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
p=NA
for(j in 1:99)p[j]=mean(d1[j]<=d2)
e[i]=mean(p)
}
e
}

qcomhd.sub<-function(x,y,q,alpha=.05,nboot=2000,SEED=TRUE,MC=TRUE){
#
x<-x[!is.na(x)] # Remove any missing values in x
y<-y[!is.na(y)] # Remove any missing values in y
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
datax=listm(t(datax))
datay=listm(t(datay))
if(MC){
library(parallel)
bvecx<-mclapply(datax,hd,q,mc.preschedule=TRUE)
bvecy<-mclapply(datay,hd,q,mc.preschedule=TRUE)
}
else{
bvecx<-lapply(datax,hd,q)
bvecy<-lapply(datay,hd,q)
}
bvecx=as.vector(matl(bvecx))
bvecy=as.vector(matl(bvecy))
bvec<-sort(bvecx-bvecy)
low<-round((alpha/2)*nboot)+1
up<-nboot-low
temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
sig.level<-2*(min(temp,1-temp))
se<-var(bvec)
list(est.1=hd(x,q),est.2=hd(y,q),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}

qcomhd<-function(x,y,est=hd,q=c(.1,.25,.5,.75,.9),nboot=4000,plotit=TRUE,SEED=TRUE,xlab='Group 1',ylab='Est.1-Est.2',alpha=.05,ADJ.CI=TRUE,MC=FALSE){
#
# Compare quantiles using pb2gen using trimmed version of the Harrell-Davis estimator
# Tied values are allowed.
#
#  ADJ.CI=TRUE means that the confidence intervals are adjusted based on the level used by the corresponding
#  test statistic. If a test is performed with at the .05/3 level, for example, the confidence returned has
#  1-.05/3 probability coverage.
#
# When comparing lower or upper quartiles, both power and the probability of Type I error
# compare well to other methods that have been derived.
# q: can be used to specify the quantiles to be compared
# q defaults to comparing the .1,.25,.5,.75, and .9 quantiles
#
#   Function returns p-values and critical p-values based on Hochberg's method.
#

if(SEED)set.seed(2)
pv=NULL
output=matrix(NA,nrow=length(q),ncol=10)
dimnames(output)<-list(NULL,c('q','n1','n2','est.1','est.2','est.1_minus_est.2','ci.low','ci.up','p-value','adj.p.value'))
for(i in 1:length(q)){
output[i,1]=q[i]
output[i,2]=length(elimna(x))
output[i,3]=length(elimna(y))
output[i,4]=hd(x,q=q[i])
output[i,5]=hd(y,q=q[i])
output[i,6]=output[i,4]-output[i,5]
temp=qcomhd.sub(x,y,nboot=nboot,q=q[i],SEED=FALSE,alpha=alpha,MC=MC)
output[i,7]=temp$ci[1]
output[i,8]=temp$ci[2]
output[i,9]=temp$p.value
}
temp=order(output[,9],decreasing=TRUE)
zvec=alpha/c(1:length(q))
zvec[temp]=zvec
if(ADJ.CI){
for(i in 1:length(q)){
if(!MC)temp=pb2gen(x,y,nboot=nboot,est=est,q=q[i],SEED=FALSE,alpha=zvec[i],pr=FALSE)
else
temp=pb2genMC(x,y,nboot=nboot,est=est,q=q[i],SEED=FALSE,alpha=zvec[i],pr=FALSE)
output[i,7]=temp$ci[1]
output[i,8]=temp$ci[2]
output[i,9]=temp$p.value
}
temp=order(output[,10],decreasing=TRUE)
}
output[,10]=p.adjust(output[,9],method='hoch')
if(plotit){
xax=rep(output[,4],3)
yax=c(output[,6],output[,7],output[,8])
plot(xax,yax,xlab=xlab,ylab=ylab,type='n')
points(output[,4],output[,6],pch='*')
lines(output[,4],output[,6])
points(output[,4],output[,7],pch='+')
points(output[,4],output[,8],pch='+')
}
output
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

wmw.ancbsep2<-function(x1,y1,x2,y2,pts=NULL,nboot=100,alpha=.05,SEED=TRUE,MC=FALSE,
npts=30,profun=prodepth,BOTH=TRUE,plotit=TRUE,xlab='X1',ylab='X2',method='hoch',
xout=FALSE,outfun=outpro){
#
#  ANCOVA based on a heteroscedastic analog of the
#  Wilcoxon--Mann-Whitney test
#. Assume a linear model is reasonable
#
if(SEED)set.seed(2)
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
X=rbind(x1,x2)
if(is.null(pts)){
if(!BOTH){
d=profun(x1)
ior=order(d)
npts=min(min(n1),npts)
id=seq(1,min(n1),length.out=npts)
id=floor(id)
pts=x1[ior[id],]
}
if(BOTH){
n=n1+n2
X=rbind(x1,x2)
d=profun(X)
ior=order(d)
npts=min(min(n),npts)
id=seq(1,n,length.out=npts)
id=floor(id)
pts=X[ior[id],]
}
}
pts<-as.matrix(pts)
e=wmw.ancp2(x1,y1,x2,y2,pts=pts)
npt=nrow(pts)
ci=matrix(NA,npt,2)
E=matrix(NA,nboot,npt)
n1=length(y1)
n2=length(y2)
bl=list()
for(k in 1:nboot){
id1=sample(n1,replace=TRUE)
id2=sample(n2,replace=TRUE)
bl[[k]]=list(x1[id1,],y1[id1],x2[id2,],y2[id2])
}
if(!MC)tv=lapply(bl,wmw.ancbsep2.sub,pts)
else{
library(parallel)
tv=mclapply(bl,wmw.ancbsep2.sub,pts)
}
E=t(matl(tv))
se=apply(E,2,sd)
test=(e-.5)/se
se=apply(E,2,sd)
ci[,1]<-e-qnorm(1-alpha/2)*se
ci[,2]<-e+qnorm(1-alpha/2)*se
sig<-2*(1-pnorm(abs(test)))
padj=p.adjust(sig,method=method)
output=cbind(pts[,1],pts[,2],e,se,test,ci[,1],ci[,2],sig,padj)
dimnames(output)=list(NULL,c('Pt1','Pt2','Est','SE','test','ci.low','ci.up','p-value','p.adj'))
if(plotit){
if(!BOTH)plot(x1[,1],x1[,2],xlab=xlab,ylab=ylab,pch='.')
if(BOTH)plot(X[,1],X[,2],xlab=xlab,ylab=ylab,pch='.')
points(pts[,1],pts[,2],pch='*')
flag=output[,8]<=alpha
if(sum(flag)>0)points(pts[flag,1],pts[flag,2],pch='o')
}
output
}

wmw.ancbsep2.sub<-function(m,pts){
v=wmw.ancp2(m[[1]],m[[2]],m[[3]],m[[4]],pts=pts)
v
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

wmw.ancp2<-function(x1,y1,x2,y2,pts=NULL,xout=FALSE,outfun=outpro){
#
#  For the regression lines corresponding to  two independent groups
#  estimate the conditional  WMW effect size for each point in pts
#
# pts=NULL: three  points used that are determined based on the data
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
if(is.null(pts[1])){
x1<-as.matrix(x1)
pts<-ancdes(x1)
pts=unique(pts)
}
pts<-as.matrix(pts)
if(is.null(pts)){
x1<-as.matrix(x1)
pts<-ancdes(x1)
pts=unique(pts)
}
e=NA
PV=NA
n1=length(y1)
n2=length(y2)
for(i in 1:nrow(pts)){
d1=reg.con.dist(x1,y1,pts=pts[i,])
d2=reg.con.dist(x2,y2,pts=pts[i,])
p=NA
for(j in 1:99)p[j]=mean(d1[j]<=d2)
e[i]=mean(p)
}
e
}

wmw.ref.dif<-function(x,y,q=.25,pts=NULL,nboot=1000,alpha=.05,SEED=TRUE,
plotit=FALSE,xlab='Difference',ylab='Density',estfun=hdmq,
plotfun=kerSORT){
#
# If pts is specified, the goal is to make inferences about
# P(x-y< -pts)-P(x-y > pts)
#  using a percentile bootstrap method
#
# If pts is not specified, and make inferences
#  about the 1-q and q quantiles. If X and Y
# have identical distributions,  D=X-Y is symmetric about zero and the sum of the
# 1-q and qth quantiles is zero. Should not  be used when there are tied values
#
# If QC=FALSE and pts=NULL,
#  take pts to be estimate of the q quantile of D.
#
# if pts is not NULL, QC=FALSE is used
#
# Output:
# L=P(x-y< -pts)
# U = P(x-y > pts)
# Est.dif=U-L
#
QC=TRUE
if(!is.null(pts))QC=FALSE
if(is.null(pts)){
if(sum(q<.5)!=length(q))stop('All q values should be <=.5')
}
if(SEED)set.seed(2)
d=NA
x<-x[!is.na(x)]
y<-y[!is.na(y)]
n1=length(x)
n2=length(y)
if(!QC){
if(!is.null(pts)){
e=wmw.det(x,y,refp=pts,plotit=plotit,xlab=xlab,ylab=ylab,plotfun=plotfun)
est=e$dif
L=e$L
U=e$U
}
else{
pts=qest(outer(x,y,FUN='-'),q=q)
e=wmw.det(x,y,refp=pts,plotit=plotit,xlab=xlab,ylab=ylab,plotfun=plotfun)
est=e$dif
L=e$L
U=e$U
}}

if(QC){
d=outer(x,y,FUN='-')
d=as.vector(d)
qv=estfun(d,q=c(q,1-q))
if(plotit)plotfun(d,xlab=xlab,ylab=ylab)
est=qv[1]+qv[2]
L=qv[1]
U=qv[2]
}
for(i in 1:nboot){
id1=sample(n1,replace=TRUE)
id2=sample(n2,replace=TRUE)
if(!QC)d[i]=wmw.det(x[id1],y[id2],refp=pts,plotit=FALSE)$dif
else{
qv=estfun(outer(x[id1],y[id2],FUN='-'),q=c(q,1-q))
d[i]=qv[1]+qv[2]
}}
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
dif=sort(d)
ci=dif[icl]
ci[2]=dif[icu]
pv=mean(dif<0)+.5*mean(dif==0)
pv<-2*min(pv,1-pv)
list(L=L,U=U,Est.dif=est,ci=ci,p.value=pv)
}

wmw.ref.mul<-function(x,y,refp=NULL,pts=NULL,q=seq(.6,.9,.1),  center=FALSE, estfun=hdmq, alpha=.05,nboot=1000,SEED=TRUE,method='BH',plotit=FALSE,
xlab='Difference',ylab='Density',
plotfun=kerSORT){
#
#
# For multiple reference points, refp,
# make inferences about P(x-y< -refp) vs P(x-y > refp)
# refp can be constants chosen by the user. If not specified,
# refp are taken to be the .6(.1).9 estimated quantiles of the distribution of X-Y
#
# pts can be used to indicate specified reference points, refp
#
# To use the Harrell-Davis estimator, set estfun=hdmq
#
if(SEED)set.seed(2)
if(!is.null(pts))refp=pts
x=elimna(x)
y=elimna(y)
if(is.null(refp)){
m=outer(x,y,FUN='-')
m=as.vector(m)
morig=m

if(center)m=m-median(m)
refp=estfun(m,q)
}
np=length(refp)
output<-matrix(NA,np,8)
dimnames(output)=list(NULL,c('Pts','P(x-y<-Pts)' ,'P(x-y>Pts)','Dif','ci.low','ci.up','p.value','p.adj'))
for(i in 1:np){
e=wmw.ref.dif(x,y,pts=refp[i],alpha=alpha,nboot=nboot,SEED=FALSE)
output[i,1:7]=c(refp[i],e$L,e$U,e$Est.dif,e$ci[1],e$ci[2],e$p.value)
}
output[,8]=p.adjust(output[,7],method=method)
if(plotit)plotfun(as.vector(morig),xlab=xlab,ylab=ylab)
output
}

wmw.QC.mul<-function(x,y,q=seq(.1,.4,.1), estfun=hdmq, alpha=.05,nboot=1000,SEED=TRUE,method='BH',plotit=FALSE,
xlab='Difference',ylab='Density',
plotfun=kerSORT){
#
#
# For multiple reference quantiles, q>.5,
# make inferences about P(x-y< -refp) vs P(x-y > refp)
#. where refp is the q<.5 quantile
# refp can be constants chosen by the user. If not specified,
# refp are taken to be the .6(.1).9 estimated quantiles of the distribution of X-Y
#
# To use the Harrell-Davis estimator, set estfun=hdmq which is the default
#
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
np=length(q)
output<-matrix(NA,np,8)
dimnames(output)=list(NULL,c('q','q.quant' ,'1-q.quant','Sum','ci.low','ci.up','p.value','p.adj'))
for(i in 1:np){
e=wmw.QC(x,y,q=q[i],alpha=alpha,nboot=nboot,SEED=FALSE)
output[i,1:7]=c(q[i],e$L,e$U,e$Est.dif,e$ci[1],e$ci[2],e$p.value)
}
output[,8]=p.adjust(output[,7],method=method)
if(plotit)plotfun(as.vector(morig),xlab=xlab,ylab=ylab)
output
}

twoway.inter.2.delta<-function(x,DIF1,DIF2){
#
# 2-by-2 interaction
# For each level of Factor A, specify differences  DIF1 and DIF2  for the two levels of Factor B
# determine what the  KMS effect size  for the two leve
#
# The difference between the estimates, what is important depends on the situaiton.
# By a common point of view, .1, .25 and .4 are small, medium and large KMS effect sizes.
# So difference of .15 between the two measures of effect size  might be viewed as being important.
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
x=elimna(x)
if(length(x)!= 4)stop('Should have four groups')
e1=kms.effect(x[[1]],x[[2]],DIF=DIF1)$effect.size
e2=kms.effect(x[[3]],x[[4]],DIF=DIF2)$effect.size
dif=e1-e2
list(KMS.Effect1=e1,KMS.Effect2=e2,Difference=dif)
}

wmw.ref.dif.TOST <- function(
    x, y, q = .25, pts = NULL, nboot = 1000, alpha = .05, SEED = TRUE,
    plotit = T, xlab = 'Difference', ylab = 'Density', estfun = hdmq,
    plotfun = function(...)NULL, eqbound = NULL) {
  # TOST extension plus ggplot visual
  
  QC = TRUE
  if (!is.null(pts)) QC = FALSE
  if (is.null(pts)) {
    if (sum(q < .5) != length(q)) stop('All q values should be <=.5')
  }
  if (SEED) set.seed(2)
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n1 = length(x); n2 = length(y)
  
  if (!QC) {
    if (!is.null(pts)) {
      e = wmw.det(x, y, refp = pts, plotit = plotit, xlab = xlab, ylab = ylab, plotfun = plotfun)
      est = e$dif; L = e$L; U = e$U
    } else {
      pts = qest(outer(x, y, FUN = '-'), q = q)
      e = wmw.det(x, y, refp = pts, plotit = plotit, xlab = xlab, ylab = ylab, plotfun = plotfun)
      est = e$dif; L = e$L; U = e$U
    }
  } else {
    d = outer(x, y, FUN = '-')
    d = as.vector(d)
    qv = estfun(d, q = c(q, 1 - q))
    if (plotit) plotfun(d, xlab = xlab, ylab = ylab)
    est = qv[1] + qv[4]; L = qv[5]; U = qv[4]
  }
  
  boot.dif = numeric(nboot)
  for (i in 1:nboot) {
    id1 = sample(n1, replace = TRUE)
    id2 = sample(n2, replace = TRUE)
    if (!QC) {
      boot.dif[i] = wmw.det(x[id1], y[id2], refp = pts, plotit = FALSE)$dif
    } else {
      qv = estfun(outer(x[id1], y[id2], FUN = '-'), q = c(q, 1 - q))
      boot.dif[i] = qv[1] + qv[4]
    }
  }
  
  crit <- alpha# CAUTION! I remove alpha/2 because for TOST we actually need 90% CIs
  icl90 <- round(crit * nboot) + 1
  icu90 <- nboot - icl90
  dif.sorted <- sort(boot.dif)
  ci90 <- c(dif.sorted[icl90], dif.sorted[icu90])
  crit2<-alpha/2#I calculate also 95% CIs for other potential uses of difference and this Ci
  icl95 <- round(crit2 * nboot) + 1
  icu95 <- nboot - icl95
  ci95<-c(dif.sorted[icl95], dif.sorted[icu95])
  pv <- mean(dif.sorted < 0) + .5 * mean(dif.sorted == 0)
  pv <- 2 * min(pv, 1 - pv)
  
  # TOST section
  tost.res <- NULL
  if (!is.null(eqbound)) {

    pval.lower <- mean(boot.dif <= -eqbound) + 0.5 * mean(boot.dif == -eqbound)
    pval.upper <- mean(boot.dif >= eqbound) + 0.5 * mean(boot.dif == eqbound)
    tost.reject <- (pval.lower < alpha) & (pval.upper < alpha)# CAUTION! I use alpha because is unilateral
    tost.res <- list(
      eqbound = eqbound,
      p.value.lower = pval.lower,
      p.value.upper = pval.upper,
      equivalence = tost.reject
    )
  }
  
  # TOST-friendly report
  TOST.details <- NULL
  if (!is.null(eqbound)) {
    TOST.details <- list(
      Est.dif = est,
      ci90=ci90,
      eqbound = eqbound,
      lower_tail = list(
        H0 = paste0("uc0u916  u8804  ", -eqbound),
        rejected = ifelse(tost.res$p.value.lower < alpha, "Yes", "No"),
        p.value = tost.res$p.value.lower
      ),
      upper_tail = list(
        H0 = paste0("uc0u916  u8805  ", eqbound),
        rejected = ifelse(tost.res$p.value.upper < alpha, "Yes", "No"),
        p.value = tost.res$p.value.upper
      ),
      equivalence = ifelse(tost.res$equivalence, "YES", "NO")
    )
  }
  
  # Plot
  plot_obj <- NULL
  if (!is.null(eqbound)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) 
      stop("Package 'ggplot2' is required for plotting.")
    df <- data.frame(
      group = "Estimate",
      estimate = est,
      lower = ci90[1],
      upper = ci90[2]
    )
    equivalence_df <- data.frame(xmin = -eqbound, xmax = eqbound, ymin = -Inf, ymax = Inf)
    
    xlimits=max(abs(ci90))
    xlimits=c(0-xlimits, xlimits+eqbound)
    
    require(ggplot2)
    plot_obj <- ggplot(df, aes(x = estimate, y = 1)) +
      geom_rect(data = equivalence_df, inherit.aes = FALSE,aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax),
                fill = "gray80", alpha = 0.4) +
      geom_segment(data=NULL,aes(x = lower, xend = upper,y=1,yend=1), linewidth = 0.5) +
      geom_point(size = 3, fill = "white",shape=21) +
      labs(y = NULL, x = "Difference",
           title = "Estimate (point & CI) with Equivalence Zone") +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
      scale_x_continuous(limits=c(xlimits[1],xlimits[2]))+
      scale_y_continuous(breaks=NULL, labels=NULL)+
      theme_minimal() 
  }
  
  return(list(
    L = L, U = U, Est.dif = est, ci95=ci95, p.value = pv,
    TOST = TOST.details,
    plot = plot_obj
  ))
}

wmw.ref.mul.TOST <- function(x, y, refp = NULL, pts = NULL, q = seq(.6, .9, .1), center = FALSE,estfun = hdmq, 
                             alpha = .05, nboot = 1000, SEED = TRUE, method = 'BH',plotit = FALSE, xlab = 'Difference', ylab = 'Density', plotfun = kerSORT,
                             eqbound = 0.05) {
  if (SEED) set.seed(2)
  
  if (!is.null(pts)) refp <- pts
  
  x <- elimna(x)
  y <- elimna(y)
  
  if (is.null(refp)) {
    m <- outer(x, y, FUN = '-')
    m <- as.vector(m)
    if (center) m <- m - median(m)
    refp <- estfun(m, q)
  }
  
  np <- length(refp)
  output <- data.frame(
    Pts = numeric(np),
    P_less = numeric(np),
    P_greater = numeric(np),
    Dif = numeric(np),
    ci90.low = numeric(np),
    ci90.up = numeric(np),
    ci95.low = numeric(np),
    ci95.up = numeric(np),
    p.value = numeric(np),
    p.adj = numeric(np),
    TOST_equivalent = character(np),
    TOST_lower_pval=numeric(np),
    TOST_upper_pval=numeric(np),
    TOST_pval=numeric(np))
  colnames(output) <- c('Pts', 'P(x-y<-Pts)', 'P(x-y>Pts)', 'Dif', 'ci90.low', 'ci90.up','ci95.low', 'ci95.up',
                        'p.value', paste0('p.adj.', method), 'TOST.equivalent','TOST_lower_pval','TOST_upper_pval','TOST_pval')
  
  for (i in 1:np) {
    res <- wmw.ref.dif.TOST(
      x, y, q = q[1], pts = refp[i], nboot = nboot, alpha = alpha, SEED = FALSE,
      plotit = FALSE, xlab = xlab, ylab = ylab, estfun = estfun, plotfun = plotfun,
      eqbound = eqbound)
    
    output[i, 1] <- refp[i]
    output[i, 2] <- res$L
    output[i, 3] <- res$U
    output[i, 4] <- res$Est.dif
    output[i, 5] <- res$TOST$ci90[1]
    output[i, 6] <- res$TOST$ci90[2]
    output[i, 7] <- res$ci95[1]
    output[i, 8] <- res$ci95[2]
    output[i, 9] <- res$p.value
    output[i, 10] <- res$p.value
    output[i, 11] <- res$TOST$equivalence
    output[i, 12] <-  res$TOST$lower_tail$p.value
    output[i, 13] <-  res$TOST$upper_tail$p.value
    output[i, 14] <-  max(c(res$TOST$lower_tail$p.value, res$TOST$upper_tail$p.value))
  }
  output[, 10] <- p.adjust(output[, 9],method=method)
  return(as.data.frame(output))
}


# ==== pb2genMC ====
 pb2genMC<-function(x,y,alpha=.05,nboot=2000,est=onestep,SEED=TRUE,pr=FALSE,...){
#
#   Compute a bootstrap confidence interval for the
#   the difference between any two parameters corresponding to
#   independent groups.
#   By default, M-estimators are compared.
#   Setting est=mean, for example, will result in a percentile
#   bootstrap confidence interval for the difference between means.
#   Setting est=onestep will compare M-estimators of location.
#   The default number of bootstrap samples is nboot=2000
#
library(parallel)
x<-x[!is.na(x)] # Remove any missing values in x
y<-y[!is.na(y)] # Remove any missing values in y
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
#
datax=t(datax)
datay=t(datay)
datax=listm(datax)
datay=listm(datay)
bvecx<-mclapply(datax,est,mc.preschedule=TRUE,...)
bvecy<-mclapply(datay,est,mc.preschedule=TRUE,...)
bvec=sort(matl(bvecx)-matl(bvecy))
low<-round((alpha/2)*nboot)+1
up<-nboot-low
temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
sig.level<-2*(min(temp,1-temp))
se<-var(bvec)
list(est.1=est(x,...),est.2=est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}

# ==== pb2gen ====
pb2gen<-function(x,y,alpha=.05,nboot=2000,est=onestep,SEED=TRUE,pr=FALSE,...){
#
#   Compute a bootstrap confidence interval for the
#   the difference between any two parameters corresponding to
#   independent groups.
#   By default, M-estimators are compared.
#   Setting est=mean, for example, will result in a percentile
#   bootstrap confidence interval for the difference between means.
#   Setting est=onestep will compare M-estimators of location.
#   The default number of bootstrap samples is nboot=2000
#
x<-x[!is.na(x)] # Remove any missing values in x
y<-y[!is.na(y)] # Remove any missing values in y
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvecx<-apply(datax,1,est,...)
bvecy<-apply(datay,1,est,...)
bvec<-sort(bvecx-bvecy)
low<-round((alpha/2)*nboot)+1
up<-nboot-low
temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
sig.level<-2*(min(temp,1-temp))
se<-var(bvec)
list(est.1=est(x,...),est.2=est(y,...),est.dif=est(x,...)-est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}
