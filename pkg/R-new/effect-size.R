# ============================================================================
# Effect Size Functions
# ============================================================================
# 
# This module contains general effect size estimation and reporting functions.
# 
# Core functions:
#   - ES.summary, ES.summary.CI - General effect size summaries
#   - qhat, qhatd - Q statistic (nonparametric effect size)
#   - akp.effect - AKP robust effect size
#   - Various factorial/interaction effect sizes
#
# Dependencies: Depends on functions from utils-core, location, bootstrap
# 
# ============================================================================


shiftdhd<-function(x,y,nboot=200,plotit=TRUE,plotop=FALSE,SEED=TRUE,pr=TRUE,xlab='x (first group)',
ylab='Delta'){
#
#   Compute confidence intervals for the difference between deciles
#   of two dependent groups. The simultaneous probability coverage is .95.
#   The Harrell-Davis estimate of the qth quantile is used.
#   The default number of bootstrap samples is nboot=200
#
#   The results are stored and returned in a 9 by 4 matrix,
#   the ith row corresponding to the i/10 quantile.
#   The first column is the lower end of the confidence interval.
#   The second column is the upper end.
#   The third column is the estimated difference between the deciles
#   (second group minus first).
#   The fourth column contains the estimated standard error.
#
#   No missing values are allowed.
#
if(pr){
print("NOTE:  if the goal is to use an alpha value different from .05")
print("use the function qcomdhd or  qdec2ci")
}
xy=elimna(cbind(x,y))
x=xy[,1]
y=xy[,2]
plotit<-as.logical(plotit)
if(SEED)set.seed(2) # set seed of random number generator so that
#   results can be duplicated.
crit<-37/length(x)^(1.4)+2.75
if(pr)print("The approximate .05 critical value is")
if(pr)print(crit)
m<-matrix(0,9,6)
if(pr)print("Taking Bootstrap Samples. Please wait.")
data<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
xmat<-matrix(x[data],nrow=nboot,ncol=length(x))
ymat<-matrix(y[data],nrow=nboot,ncol=length(x))
for (i in 1:9){
q<-i/10
bvec<-apply(xmat,1,hd,q)-apply(ymat,1,hd,q)
se<-sqrt(var(bvec))
dif<-hd(x,q)-hd(y,q)
m[i,1]=hd(x,q)
m[i,2]=hd(y,q)
m[i,3]<-dif
m[i,4]<-dif-crit*se
m[i,5]<-dif+crit*se
m[i,6]<-se
}
dimnames(m)<-list(NULL,c('est.1','est.2','est.dif','ci.lower','ci.upper','se'))
if(plotit){
if(plotop){
xaxis<-c(1:9)/10
xaxis<-c(xaxis,xaxis)
}
if(!plotop)xaxis<-c(deciles(x),deciles(x))
par(pch="+")
#yaxis<-c(m[,1],m[,2])
yaxis<-c(m[,4],m[,5])
if(!plotop)plot(xaxis,yaxis,ylab=ylab,xlab=xlab)
if(plotop)plot(xaxis,yaxis,ylab="delta",xlab="Deciles")
par(pch="*")
if(!plotop)points(deciles(x),m[,3])
if(plotop)points(c(1:9)/10,m[,3])
}
m
}


qhatds1<-function(isubx,x,y){
#
#  function used by qhat  when working on bootstrap estimates.
#
xx<-x[isubx]
yy<-y[isubx]
group<-disker(xx,yy,x,op=2)$zhat
group
}


qhatd<-function(x,y,nboot=50){
#
#   Estimate Q, a nonparametric measure of effect size, using
#   the .632 method of estimating prediction error.
#   (See Efron and Tibshirani, 1993, pp. 252--254)
#
#   The default number of bootstrap samples is nboot=50
#
#   This function is for dependent groups. For independent groups, use
#   qhati
#
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
#    data is an nboot by n matrix containing subscripts for bootstrap sample
bid<-apply(data,1,idb,length(x))
#  bid is a n by nboot matrix. If the jth bootstrap sample from
#  1, ..., n contains the value i, bid[i,j]=0; otherwise bid[i,j]=1
yhat<-apply(data,1,qhatds1,x,y)
bi<-apply(bid,1,sum) # B sub i in notation of Efron and Tibshirani, p. 253
temp<-(bid*yhat)
diff<-apply(temp,1,sum)
temp<-diff/bi
ep0<-sum(temp[!is.na(temp)])/length(y)
aperror<-disker(x,y)$phat  # apparent error
regpre<-.368*aperror+.632*ep0
list(app.error=aperror,qhat.632=regpre)
}


qhat<-function(x,y,nboot=50,op=2,SEED=TRUE,pr.track=FALSE){
#
#   Estimate Q, a nonparametric measure of effect size, using
#   the .632 method of estimating prediction error.
#   (See Efron and Tibshirani, 1993, pp. 252--254)
#
#   The default number of bootstrap samples is nboot=50
#
#   Missing values are automatically removed
#
# op=1, use Rosenblatt's shifted histogram version of kernel estimate
# op=2, use adaptive kernel estimate with initial estimate based
#       on expected frequency curve.
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
datax<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
#    datax is an nboot by n matrix containing subscripts for bootstrap sample
#    associated with first group.
datay<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
#    datay is an nboot by m matrix containing subscripts for bootstrap sample
#    associated with second group.
bidx<-apply(datax,1,idb,n=length(x))
#  bidx is a n by nboot matrix. If the jth bootstrap sample from
#  1, ..., n contains the value i, bid[i,j]=0; otherwise bid[i,j]=1
bidy<-apply(datay,1,idb,n=length(y))
temp3<-matrix(0,ncol=length(x),nrow=nboot)
temp5<-matrix(0,ncol=length(y),nrow=nboot)
for(i in 1:nboot){
temp3[i,]<-disker(x[datax[i,]],y[datay[i,]],x,op=op)$zhat
# temp3 contains vector of 0s and 1s, 1 if x[i] is
#  is classified as coming from group 1.
temp5[i,]<-disker(y[datay[i,]],x[datax[i,]],y,op=op)$zhat
if(pr.track)print(paste("Iteration ", i, "of ", nboot," is complete"))
}
temp4<-temp3*t(bidx)
temp4<-apply(temp4,2,sum)/apply(bidx,1,sum)
temp6<-temp5*t(bidy)
temp6<-apply(temp6,2,sum)/apply(bidy,1,sum)
ep0x<-mean(temp4,na.rm=TRUE)  # epsilon hat_x
aperrorx<-disker(x,y,op=op)$phat  # apparent error
regprex<-.368*aperrorx+.632*ep0x
ep0y<-mean(temp6,na.rm=TRUE)
aperrory<-disker(y,x,op=op)$phat  # apparent error
regprey<-.368*aperrory+.632*ep0y
aperror<-(length(x)*aperrorx+length(y)*aperrory)/(length(x)+length(y))
regpre<-(length(x)*regprex+length(y)*regprey)/(length(x)+length(y))
list(qhat.632=regpre)
}


akp.effect<-function(x,y,EQVAR=TRUE,tr=.2){
#
# Computes the robust effect size suggested by
#Algina, Keselman, Penfield Psych Methods, 2005, 317-328
x<-elimna(x)
y<-elimna(y)
n1<-length(x)
n2<-length(y)
s1sq=winvar(x,tr=tr)
s2sq=winvar(y,tr=tr)
spsq<-(n1-1)*s1sq+(n2-1)*s2sq
sp<-sqrt(spsq/(n1+n2-2))
cterm=1
if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
cterm=sqrt(cterm)
if(EQVAR)dval<-cterm*(tmean(x,tr)-tmean(y,tr))/sp
if(!EQVAR){
dval<-cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s1sq)
dval[2]=cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s2sq)
}
dval
}


akp.effect.ci<-function(x,y,alpha=.05,tr=.2,nboot=1000,SEED=TRUE,null.val=0){
#
# Computes the robust effect size for two-sample case using
# Algina, Keselman, Penfield Pcyh Methods, 2005, 317-328
#
#
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
be.f=NA
for(i in 1:nboot){
X=sample(x,n1,replace=TRUE)
Y=sample(y,n2,replace=TRUE)
be.f[i]=akp.effect(X,Y,tr=tr)
}
L=alpha*nboot/2
U=nboot-L
be.f=sort(be.f)
ci=be.f[L+1]
ci[2]=be.f[U]
est=akp.effect(x,y,tr=tr)
pv=mean(be.f<null.val)+mean(be.f==null.val)
pv=2*min(c(pv,1-pv))
list(akp.effect=est,ci=ci,p.value=pv)
}


qhatDEP<-function(x1,x2,depthfun=prodepth,...){
#
#  Compute apparent probability of correct classification
#
x1<-x1[!is.na(x1)]
x2<-x2[!is.na(x2)]
x1=as.matrix(x1)
x2=as.matrix(x2)
tv=c(rep(1,nrow(x1)),rep(2,nrow(x2)))
see=discdepth(x1,x2,z=rbind(x1,x2))
qhat=mean(tv==see)
qhat
}


qhatdepPB<-function(x1,x2,nboot=500,alpha=.05,depthfun=prodepth,
SEED=TRUE,...){
#
#
if(SEED)set.seed(2)
bvec=NA
x1=as.matrix(x1)
x2=as.matrix(x2)
n1=nrow(x1)
n2=nrow(x2)
for(i in 1:nboot){
dat1=sample(n1,n1,replace=TRUE)
dat2=sample(n2,n2,replace=TRUE)
bvec[i]=qhatDEP(x1[dat1,],x2[dat2,],depthfun=depthfun)
}
est=qhatDEP(x1,x2)
bvec=sort(bvec)
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
ci<-bvec[icl]
ci[2]<-bvec[icu]
list(estimate=est,ci=ci)
}


lin.ES<-function(x,con,locfun=tmean,tr=.2,nreps=200,SEED=TRUE,
MAIN=FALSE,INT=FALSE,POOL=FALSE,nmax=10^8,...){
#
#  Estimate a  quantile shift measure of effect size for a linear contrast based
#  on an estimate of the distribution of the linear contrast.
#
#  con = contrast coefficients
#  x is a matrix or has list mode.
#
#  POOL=FALSE: Estimate distribution of sum c_iX_i, c_i are contrast coefficient
#              Then compute a relative shift measure of effect size.
#  POOL=TRUE: Pool data corresponding to c_i=1, do the same for c_i=-1, then
#  use then compute a relative shift measure of effect size.
#
#
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
nv=as.vector(matl(lapply(x,FUN='length')))
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
nv=as.vector(matl(lapply(x,FUN='length')))
if(length(con)!=J)stop('Length of con should equal number of groups')
x=elimna(x)
if(!POOL){
np=prod(nv)
nmin=min(nv)
if(np>nmax){
nmin=min(c(nmin,100))
}
B=list()
M=matrix(NA,nrow=nmin,ncol=J)
for(i in 1:nreps){
for(j in 1:J)M[,j]=sample(x[[j]],nmin)
B[[i]]=M
}
L=lapply(B,linWMWMC.sub,con=con)
ef.size=NA
for(j in 1:length(L))ef.size[j]=linES.sub(L[[j]],locfun=locfun,...)
ef=mean(ef.size)
}
if(POOL){
y=list()
id1=which(con==1)
id2=which(con==-1)
v1=pool.a.list(x[id1])
v2=pool.a.list(x[id2])
if(length(v1)*length(v2)<nmax){
L=outer(pool.a.list(x[id1]),pool.a.list(x[id2]),FUN='-')
ef=linES.sub(L,locfun=locfun,...)
}
if(length(v1)*length(v2)>=nmax){
B=list()
M=matrix(NA,nrow=nmin,ncol=J)
for(i in 1:nreps){
for(j in 1:J)M[,j]=sample(x[[j]],nmin)
B[[i]]=M
}
L=lapply(B,linWMWMC.sub,con=con)
ef.size=NA
for(j in 1:length(L))ef.size[j]=linES.sub(L[[j]],locfun=locfun,...)
ef=mean(ef.size)
}}
list(Effect.Size=ef)
}


linES.sub<-function(L,locfun,...){
est=locfun(L,...)
ef.size=mean(L-est<=est)
ef.size
}


rmlinES<-function(x, con = NULL){
#
#  Dependent groups:
#  For each linear contrast, compute Algina et al. effect size based on the linear sum
#
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
J=ncol(x)
if(is.null(con)){
C=(J^2-J)/2
con=matrix(0,ncol=C,nrow=J)
ic=0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic=ic+1
con[j,ic]=1
con[k,ic]=-1
}}}}
x=elimna(x)
n=nrow(x)
ES=NA
for (d in 1:ncol(con)){
S=NA
for(i in 1:n)S[i]=sum(con[,d]*x[i,])
ES[d]=D.akp.effect(S,tr=tr)
}
list(con=con,Effect.Size=ES)
}


ES.sum.REL.MAG<-function(REL.M,n = 10000,reps=10){
#
#  Determine small medium and large equivalent measures of effect size based on the values in
#  REL.M
#
if(length(REL.M)!=3)stop('Should have three value in REL.M')
if(n>10000)n=10000
x=rnorm(n)
y=rnorm(n)
output=matrix(0,ncol=3,nrow=6)
int=matrix(NA,ncol=3,nrow=6)
dimnames(output)=list(c('AKP','EP','QS (median)','QStr','WMW','KMS'),c('S','M','L'))
for(k in 1:reps){
for(j in 1:3)int[,j]=ES.summary(x,y-REL.M[j],)[,1]
output=output+int
}
output=output/reps
output
}


ES.summary<-function(x,y,tr=.2,NULL.V=c(0,0,.5,.5,.5,0),REL.MAG=NULL, REL.M=NULL,n.est=1000000){
#
# Estimate a collection of effect sizes:
#  AKP: Homoscedastic robust analog of Cohen's d
#  EP:  Explanatory power
#  QS:  Quantile shift based on the median of the distribution of X-Y,
#  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
#  WMW:  P(X<Y)
#  KMS:  Robust heteroscedastic analog  of Cohen's d
#
#  REL.M can be used to indicate what is a small, medium and large effect when
#  using AKP under normality and homosecedasticity. An estimate of the corresponding
#  values for the other measures of effect size is estimated based on
#  n.est values sampled from a normal distribution
#
if(!is.null(REL.MAG)) REL.M=REL.MAG
SM=c(.2,.14,.55,.55,.55,.1)
MED=c(.5,.34,.64,.64,.64,.25)
LAR=c(.8,.52,.71,.71,.71,.4)
if(is.null(REL.M)){
SM=c(.2,.14,.55,.55,.55,.1)
MED=c(.5,.34,.64,.64,.64,.25)
LAR=c(.8,.52,.71,.71,.71,.4)
}
if(!is.null(REL.M)){
v=ES.sum.REL.MAG(REL.M,n=n.est)
SM=v[,1]
MED=v[,2]
LAR=v[,3]
SM[1]=REL.M[1]
SM[6]=REL.M[1]/2
MED[1]=REL.M[2]
MED[6]=REL.M[2]/2
LAR[1]=REL.M[3]
LAR[6]=REL.M[3]/2
}
a=c('AKP','EP','QS','QStr','WMW','KMS')
output=matrix(NA,ncol=5,nrow=6)
for(j in 1:6){
output[j,1]=ESfun(x,y,method=a[j],tr=tr,pr=FALSE)
}
output[,2:5]=cbind(NULL.V,SM,MED,LAR)
if(output[1,1]<0)output[1,3:5]=-1*output[1,3:5]
for(j in 3:5){
if(!is.na(output[j,1])){
if(output[j,1]<.5) {
dif=output[j,3:5]-.5
output[j,3:5]=.5-dif
}}}
if(output[6,1]<0)output[6,3:5]=-1*output[6,3:5]
dimnames(output)=list(c('AKP','EP','QS (median)','QStr','WMW','KMS'),c('Est','NULL','S','M','L'))
output
}


ES.summary.CI<-function(x,y,tr=.2,QSfun=median,SEED=TRUE,alpha=.05,nboot=2000,method='hoch',
NULL.V=c(0,0,.5,.5,.5,0),REL.MAG=NULL,REL.M=NULL,n.est=1000000){
#
# Estimate a collection of effect sizes:
#  AKP: Homoscedastic robust analog of Cohen's d
#  EP:  Explanatory power
#  QS:  Quantile shift based on the median of the distribution of X-Y,
#  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
#  WMW:  P(X<Y)
#  KMS:  Robust heteroscedastic analog  of Cohen's d
#
a=c('AKP','EP','QS','QStr','WMW','KMS')
output=matrix(NA,ncol=9,nrow=6)
dimnames(output)=list(c('AKP','EP','QS (median)','QStr','WMW','KMS'),c('Est','NULL','S','M','L','ci.low','ci.up','p.value','adj.p.value'))
for(j in 1:6){
output[j,1]=ESfun(x,y,QSfun=QSfun,method=a[j],tr=tr,pr=FALSE,SEED=SEED)
}
if(!is.null(REL.MAG))REL.M=REL.MAG
if(is.null(REL.M)){
SM=c(.2,.14,.55,.55,.55,.1)
MED=c(.5,.34,.64,.64,.64,.25)
LAR=c(.8,.52,.71,.71,.71,.4)
}
if(!is.null(REL.M)){
v=ES.sum.REL.MAG(REL.M,n=n.est)
SM=v[,1]
MED=v[,2]
LAR=v[,3]
SM[1]=REL.M[1]
SM[6]=REL.M[1]/2
MED[1]=REL.M[2]
MED[6]=REL.M[2]/2
LAR[1]=REL.M[3]
LAR[6]=REL.M[3]/2
}
output[,2:5]=cbind(NULL.V,SM,MED,LAR)
if(output[1,1]<0)output[1,3:5]=-1*output[1,3:5]
for(j in 3:4){
if(output[j,1]<.5){
dif=output[j,3:5]-.5
output[j,3:5]=.5-dif
}}
if(output[5,1]>.5 & output[5,5]< .5){
output[5,3:5]=1-output[5,3:5]
}
if(output[5,1]<.5 & output[5,5]> .5){
output[5,3:5]=1-output[5,3:5]
}

a=akp.effect.ci(x,y,tr=tr,alpha=alpha,nboot=nboot,SEED=SEED)
output[1,6:7]=a$ci
output[1,8]=a$p.value
a=EPci(x,y,tr=tr,alpha=alpha,SEED=SEED,nboot=nboot)
output[2,6:7]=a$ci
output[2,8]=yuen(x,y,tr=tr)$p.value
a=shiftPBci(x,y,locfun=QSfun,alpha=alpha,nboot=nboot,SEED=SEED)
output[3,6:7]=a$ci
output[3,8]=a$p.value
a=shiftPBci(x,y,locfun=tmean,alpha=alpha,nboot=nboot,SEED=SEED)
output[4,6:7]=a$ci
output[4,8]=a$p.value
a=cidv2(x,y,alpha=alpha)
output[5,6:7]=a$p.ci
output[5,8]=a$p.value
a=KMS.ci(x,y,alpha=alpha,nboot=nboot,SEED=SEED)
output[6,6:7]=a$ci
output[6,8]=a$p.value
if(output[6,1]<0)output[6,3:5]=-1*output[6,3:5]
output[,9]=p.adjust(output[,8],method=method)
output
}


MUL.ES.sum<-function(x1,x2){
#
# For multivariate data, compute measures of effect for	each variable
#
V=list()
p=ncol(x1)
for(j in 1:p)V[[j]]=ES.summary(x1[,j],x2[,j])
V
}


IND.PAIR.ES<-function(x,con=NULL,fun=ES.summary,tr=.2,...){
#
#  J independent groups
#  For each column of a specified matrix of linear contrast
#  coefficients, pool the data having a contrast coefficient of 1 into one
#  group, do the same for contrast coefficient of -1,
#  then estimate measures of effect size.
#
#  Default, all pairwise comparisons using all of the measures of effect size
#  via ES.summary
#
#  To get individual measures of effect size, use
#  fun=ESfun and include the argument
#  method.
#  Example: fun=ESfun, method='EP' does explanatory power.
# Choices for method are:
#  EP: Explanatory power
#  QS: Quantile shift based on the medians
#  QStr: Based on trimmed means
#  AKP: Robust analog of Cohen's d
#  WMW:   P(X<Y)
#  KMS:    Heteroscedastic analog of Cohen's d
#
#  Returns contrast coefficients followed by effect size measures for each column
#  of con.
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
CON=list()
J=length(x)
if(is.null(con))conmat=con.all.pairs(J)
else conmat=con
P=ncol(conmat)
ic=0
for(j in 1:P){
id1=which(conmat[,j]==1)
id2=which(conmat[,j]==-1)
ic=ic+1
if(length(id1)==1)CON[[ic]]=fun(x[[id1]],x[[id2]],tr=tr,...)
else{
z1=pool.a.list(x[id1])
z2=pool.a.list(x[id2])
CON[[ic]]=fun(z1,z2,...)
}}
list(con=conmat,effect.size=CON)
}


RCES<-function(J,K,x,fun=ES.summary,tr=.2,...){
#
# For each level of Factor A, estimate effect sizes for all pairwise comparisons
#  among the levels of B
#
#  Do the same of Factor B
#
#  argument fun, see the function IND.PAIR.ES
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
JK=J*K
imat=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
A=list()
B=list()
for(j in 1:J){
z=x[imat[j,]]
A[[j]]=IND.PAIR.ES(z,tr=tr,fun=fun,...)
}
for(k in 1:K){
z=x[imat[,k]]
B[[k]]=IND.PAIR.ES(z,tr=tr,fun=fun,...)
}
list(A=A,B=B)
}


inter.ES<-function(x,method='EP',iter=5,SEED=TRUE,tr=.2,pr=FALSE){
#
#
# Measures of effect size for an  interaction in a 2-by-2 design
# For  level 1 of Factor A, estimate the distribution of the
#  the typical difference for levels 1  and 2 Factor B
#  Do the same for level 2 of Factor A, and compute  a measure of
#  effect size based on these two distributions.
#
# Choices for the argument method:
# 'DNT', 'PH,`EP',`QS',`QStr',`AKP',`WMW',`KMS'
# DNT= De Neve and Thas  method
# PH=Patel--Hoel
# EP=explanatory power,
#  QS= quantile shift (median,
#  QStr= quantile shift (trimmed mean) ,
#  AKP =trimmed mean version of Cohen's d,
# WMW= Wilcoxon type measure
# KMS=heteroscedastic analog of Cohen's d
#
#
if(is.matrix(x))x=listm(x)
ef=NA
if(length(x)!=4)stop('Limited to a two-by-two design')
x=elimna(x)
FLAG=TRUE
if(method=='DNT'){
ef=WMWinter.est(x,iter=iter,SEED=SEED)
FLAG=FALSE
}
if(method=='PH'){
ef=rimul(2,2,x)$test[,5]
FLAG=FALSE
}
if(FLAG){
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
ef[i]=ESfun(L1,L2,method=method,tr=tr,pr=pr,SEED=FALSE)
}}
else{
L1=outer(x[[1]],x[[2]],FUN='-')
L2=outer(x[[3]],x[[4]],FUN='-')
ef=ESfun(L1,L2,method=method,tr=tr,pr=pr)
}
}
ef=mean(ef)
ef
}


inter.TDES.sub<-function(x,method='QS',iter=5,SEED=TRUE,tr=.2,pr=FALSE,switch=FALSE){
#
#
# Measures of effect size for an  interaction in a 2-by-2 design
# For  level 1 of Factor A, estimate the distribution of the
#  the typical difference for levels 1  and 2 Factor B
#  Do the same for level 2 of Factor A, and compute  a measure of
#  effect size based on these two distributions.
#
#  swithch=TRUE, interchange the rows and columns
#
# Choices for the argument method:
# 'DNT',`EP',`QS',`QStr',`AKP',`KMS'
# DNT= De Neve and Thas P(X_1-X_2 < X_3-X_4) so a WMW-type measure
# EP=explanatory power,
#  QS= quantile shift (median,
#  QStr= quantile shift (trimmed mean) ,
#  AKP =trimmed mean version of Cohen's d,
# KMS=heteroscedastic analog of Cohen's d
#
#
if(is.matrix(x))x=listm(x)
if(switch)x=x[1,3,2,4]
ef=NA
if(length(x)!=4)stop('Limited to a two-by-two design')
x=elimna(x)
FLAG=TRUE
if(method=='DNT'){
ef=WMWinter.est(x,iter=iter,SEED=SEED)
FLAG=FALSE
}
if(FLAG){
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
ef[i]=ESfun(L1,L2,method=method,tr=tr,pr=pr,SEED=FALSE)
}}
else{
L1=outer(x[[1]],x[[2]],FUN='-')
L2=outer(x[[3]],x[[4]],FUN='-')
ef=ESfun(L1,L2,method=method,tr=tr,pr=pr)
}
}
ef=mean(ef)
ef
}


interES.2by2<-function(x,tr=.2,SW=FALSE){
#
# Estimate a collection of effect sizes
# for the first row of a 2-by-2 design
# do the same for the second row
# return estimates of the differences
#
#  AKP: Homoscedastic robust analog of Cohen's d
#  EP:  Explanatory power
#  QS:  Quantile shift based on the median of the distribution of X-Y,
#  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
#  KMS:  Robust heteroscedastic analog  of Cohen's d
#  PH: Patel--Hoel, uses Cliff'a analog of Wilcoxon--Mann--Whitney
#
#  switch=TRUE: reverses rows and columns

if(is.matrix(x)  || is.data.frame(x))x=listm(x)
if(SW)x=x[c(1,3,2,4)]
J=length(x)
if(J!=4)stop('Should have four groups; designed for a 2-by-2 ANOVA only')
a=c('AKP','EP','QS','QStr','KMS','WMW')
output=matrix(NA,ncol=4,nrow=6)
output[,1]=c(0.0,0.0,0.5,0.5,0.0,0.5)
for(j in 1:6){
output[j,2]=ESfun(x[[1]],x[[2]],method=a[j],tr=tr,pr=FALSE)
output[j,3]=ESfun(x[[3]],x[[4]],method=a[j],tr=tr,pr=FALSE)
output[j,4]=output[j,2]-output[j,3]
}
dimnames(output)=list(c('AKP','EP','QS (median)','QStr',
'KMS','PH'),c('NULL','Est 1','Est 2','Diff'))
output
}


interJK.ESmul<-function(J,K,x,method='QS',tr=.2,SEED=TRUE){
#
#  Compute measures of effect size for interactions associated with
#  in J-by-K design.
#  This is done for all relevant tetrad cells using interES.2by2
#  Missing values are automatically removed.
#
#  Methods, see the R function ESfun
#  Defaults to quantile shfit
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
CCJ<-(J^2-J)/2
CCK<-(K^2-K)/2
CC<-CCJ*CCK
JK=J*K
test<-matrix(NA,CC,7)
x=elimna(x)
mat=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
dimnames(test)<-list(NULL,c("Factor A","Factor A","Factor B","Factor B","Effect Size  1","Effect Size  2","Diff"))
jcom<-0
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
id1=mat[j,k]
id2=mat[j,kk]
a=ESfun(x[[id1]],x[[id2]],method=method,tr=tr,pr=FALSE,SEED=SEED)
id1=mat[jj,k]
id2=mat[jj,kk]
b=ESfun(x[[id1]],x[[id2]],method=method,tr=tr,pr=FALSE,SEED=SEED)
test[jcom,5:7]<-c(a,b,a-b)
}}}}}}
list(EFFECT.est=test)
}


LCES<-function(x,con,nreps=200,tr=.2,SEED=TRUE){
#
# For each column of con,  compute four measures of effect size:
# quantile shift based on median
# quantile shift based on a trimmed mean
# AKP generalization of Cohen's d
# SIGN: analog of the sign test.
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
x=elimna(x)
con=as.matrix(con)
d=ncol(con)
mat=matrix(NA,nrow=4,ncol=d)
LAB=NULL
for(i in 1:d){
LAB[i]=paste('Con',i)
mat[1,i]=lin.ES(x,as.vector(con[,i]),locfun=median,nreps=nreps)$Effect.Size
mat[2,i]=lin.ES(x,as.vector(con[,i]),locfun=mean,nreps=nreps,tr=tr)$Effect.Size
mat[3,i]=lin.akp(x,con[,i],locfun=mean,nreps=nreps,tr=tr)$Effect.Size
mat[4,i]=linsign(x,con[,i],nreps=nreps)
}
mat=cbind(c(0.5,0.5,0.0,0.5),mat)
LAB=c('NULL',LAB)
dimnames(mat)=list(c('QS','Qstr','AKP','SIGN'),LAB)
list(EST=mat,con=con)
}


qno.est<-function(x,q=.5){
#
#  Estimate of the qth quantile
#  In some situations, offers a distinct advantage over the Harrell-Davis estimator when
# comparing extreme quantiles and distributions have heavy tails.
#
n<-length(x)
x<-sort(x)
s<-numeric()
ifelse(n>2, {for(g in 1:(n-2)){
s[g]<-x[g+1]*(dbinom(g, size=n, prob=q)*(1-q)+dbinom(g+1, size=n, prob=q)*q)
}
sum(s,na.rm = TRUE)
t1<-(2*dbinom(0, size=n, prob=q)*q+dbinom(1, size=n, prob=q)*q)*x[1]
t2<-(2*(1-q)*dbinom(n, size=n, prob=q)+dbinom(n-1, size=n, prob=q)*(1-q))*x[n]
t3<-dbinom(0, size=n, prob=q)*(2-3*q)*x[2]-dbinom(0, size=n, prob=q)*(1-q)*x[3]-
dbinom(n, size=n, prob=q)*q*x[n-2]+dbinom(n, size=n, prob=q)*(3*q-1)*x[n-1]
quan<-sum(s,na.rm = T)+t1+t2+t3},
ifelse(n==2,{quan <- (1-q)*x[1]+q*x[2]},quan<-x))
quan
}


bw.es.A<-function(J,K,x,tr=.2,pr=TRUE,fun=ES.summary,...){
#
#  Between-by-within design.
#
# Assumed that fun function has an argument tr
#
#Using REL.M, can change default values for small, medium and large
# Example REL.M=c(.1,.3,.5)
#
# For each level of Factor B,  compute  effect sizes
# for all pairs of  levels of  Factor A .
#
#  The R variable x is assumed to contain the raw
#  x stored in list mode. x[[1]] contains the x
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the x for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the x for level 1,K
#  x[[K+1]] is the x for level 2,1, x[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested, but a subset of the x can be analyzed
#  using grp
#
#
if(pr){
if(!identical(fun,ES.summary.CI))print('To get confidence intervals, set the argument fun=ES.summary.CI')
}
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK=J*K
mat=matrix(c(1:JK),J,K,byrow=TRUE)
B=list()
for(k in 1:K){
B[[k]]=IND.PAIR.ES(x[mat[,k]],fun=fun,tr=tr,...)
if(k==1){
if(pr){
print('B[[1]] contains pairwise measures of effect size for all levels of Factor A')
print(' and level 1 of Factor B')
print(' B[[2]]   contains pairwise measures of effect size for all levels of Factor A')
print('and level 2 of Factor B')
}}}
list(B=B)
}


dep.ES.summary<-function(x,y=NULL,tr=.2, alpha=.05, REL.MAG=NULL,SEED=TRUE,nboot=2000){
#
#
# For two dependent groups,
# compute confidence intervals for four measures of effect size based on difference scores:
#
#  AKP: robust standardized difference similar to  Cohen's d
#  QS:  Quantile shift based on the median of the distribution of difference scores,
#  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
#  SIGN:  P(X<Y), probability that for a random pair, the first is less than the second.
#
#  y=NULL: Assume difference scores stored in x or only a single distribution is being used.
#
#  REL.MAG: suppose it is decided that AKP values .1, .3 and .5 are viewed as small, medium and large under normality Then
#
# REL.MAG=c(.1,.3,.5)  the default,
# means that the function will compute the corresponding values for the measures of effect size used here.

ecom=c(0.10,  0.54,  0.54,  0.46, 0.30,  0.62,
   0.62,  0.38, 0.50,  0.69,  0.69,  0.31)
REL.EF=matrix(NA,4,3)

if(!is.null(y))x=x-y
x=elimna(x)
n=length(x)
output=matrix(NA,ncol=5,nrow=4)
dimnames(output)=list(c('AKP','QS (median)','QStr','SIGN'),c('NULL','Est','S','M','L'))
output[1,1:2]=c(0,D.akp.effect(x,tr=tr))
output[2,1:2]=c(0.5,depQS(x)$Q.effect)
output[3,1:2]=c(0.5,depQS(x,locfun=mean,tr=tr)$Q.effect)
output[4,1:2]=c(0.5,mean(x[x!=0]<0))
if(is.null(REL.MAG)){
REL.MAG=c(.1,.3,.5)
REL.EF=matrix(ecom,4,3)
}
else{
if(length(REL.MAG)!=3)stop('REL.MAG should have three values')
if(SEED)set.seed(2)
z=rnorm(1000000)
v=dep.ES.summary.sub(z+REL.MAG[1],tr=tr)
REL.EF[,1]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[2],tr=tr)
REL.EF[,2]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[3],tr=tr)
REL.EF[,3]=v[,2]
}
if(output[1,2]<0)REL.EF[1,]=0-REL.EF[1,]
if(output[2,2]<0.5)REL.EF[2,]=.5-(REL.EF[2,]-.5)
if(output[3,2]<0.5)REL.EF[3,]=.5-(REL.EF[3,]-.5)
if(output[4,2]>0.5)REL.EF[4,]=.5-(REL.EF[4,]-.5)
output[,3:5]=REL.EF
output
}


dep.ES.summary.sub<-function(x,y=NULL,tr=.2){
#
#
#  Used  to determine equivalent effect size based on specified standard deviations
#

if(!is.null(y))x=x-y
output=matrix(NA,ncol=2,nrow=4)
dimnames(output)=list(c('AKP','QS (median)','QStr','SIGN'),c('NULL','Est'))
output[1,1:2]=c(0,D.akp.effect(x))
output[2,1:2]=c(0.5,depQS(x)$Q.effect)
output[3,1:2]=c(0.5,depQS(x,locfun=mean,tr=tr)$Q.effect)
output[4,1:2]=c(0.5,mean(x[x!=0]<0))
output
}


dep.ES.summary.CI<-function(x,y=NULL,tr=.2, alpha=.05, REL.MAG=NULL,SEED=TRUE,nboot=1000,AUTO=FALSE){
#
#
# For two dependent groups,
# compute confidence intervals for four measures of effect size based on difference scores:
#
#  AKP: robust standardized difference similar to  Cohen's d
#  QS:  Quantile shift based on the median of the distribution of difference scores,
#  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
#  SIGN:  P(X<Y), probability that for a random pair, the first is less than the second.
#
#  OPT=TRUE: No effect, difference scores are symmetric about zero.
#  Under normality, suppose a shift of .2, .5 and .8 standard deviation of the difference score
#  is considered small, medium and large. The corresponding values for QS and SIGN are printed.
#
#  y=NULL: Assume difference scores stored in x or only a single distribution is being used.
#
#  REL.MAG: suppose it is decided that AKP values .1, .3 and .5 are viewed as small, medium and large under normality Then
#
#  if OPT=T and
# REL.MAG=c(.1,.3,.5)  the default,
# means that the function will compute the corresponding values for the measures of effect size used here.

ecom=c(0.10,  0.54,  0.54,  0.46, 0.30,  0.62,
   0.62,  0.38, 0.50,  0.69,  0.69,  0.31)
REL.EF=matrix(NA,4,3)

if(!is.null(y))x=x-y
x=elimna(x)
n=length(x)
output=matrix(NA,ncol=8,nrow=4)
dimnames(output)=list(c('AKP','QS (median)','QStr','SIGN'),c('NULL','Est','S','M','L','ci.low','ci.up','p.value'))
output[1,1:2]=c(0,D.akp.effect(x))
output[2,1:2]=c(0.5,depQS(x)$Q.effect)
output[3,1:2]=c(0.5,depQS(x,locfun=mean,tr=tr)$Q.effect)
output[4,1:2]=c(0.5,mean(x[x!=0]<0))
if(is.null(REL.MAG)){
REL.MAG=c(.1,.3,.5)
REL.EF=matrix(ecom,4,3)
}
else{
if(length(REL.MAG)!=3)stop('REL.MAG should have three values')
if(SEED)set.seed(2)
z=rnorm(1000000)
v=dep.ES.summary.sub(z+REL.MAG[1],tr=tr)
REL.EF[,1]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[2],tr=tr)
REL.EF[,2]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[3],tr=tr)
REL.EF[,3]=v[,2]
}
if(output[1,2]<0)REL.EF[1,]=0-REL.EF[1,]
if(output[2,2]<0.5)REL.EF[2,]=.5-(REL.EF[2,]-.5)
if(output[3,2]<0.5)REL.EF[3,]=.5-(REL.EF[3,]-.5)
if(output[4,2]>0.5)REL.EF[4,]=.5-(REL.EF[4,]-.5)
output[,3:5]=REL.EF
a=D.akp.effect.ci(x,alpha=alpha,SEED=SEED,tr=tr,nboot=nboot)
output[1,6:7]=a$ci
output[1,8]=a$p.value
#output[1,6:7]=D.akp.effect.ci(x,alpha=alpha,SEED=SEED,tr=tr,nboot=nboot)$ci
a=depQSci(x,alpha=alpha,SEED=SEED,nboot=nboot)
output[2,6:7]=a$ci
output[2,8]=a$p.value
a=depQSci(x,locfun=tmean,alpha=alpha, SEED=SEED,tr=tr,nboot=nboot)
output[3,6:7]=a$ci
output[3,8]=a$p.value
Z=sum(x<0)
nm=length(x[x!=0])
a=binom.conf.pv(Z,nm,alpha=alpha,AUTO=AUTO,pr=FALSE)
output[4,6:7]=a$ci
output[4,8]=a$p.value
output
}


bw.es.B<-function(J,K,x,tr=.2,POOL=FALSE,OPT=FALSE,CI=FALSE,SEED=TRUE,REL.MAG=NULL,pr=TRUE){
#
#  Between-by-within design.
#
# For each level of Factor A,  compute  effect sizes
# for all j<m levels of  Factor B .
#
# If POOL=TRUE, pool the data over levels of Factor B
#
# If the data are stored in a matrix or data frame, it is converted to list mode.
#   Once in list,
# x[[1]] contains the x
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the x for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the x for level 1,K
#  x[[K+1]] is the x for level 2,1, x[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested, but a subset of the x can be analyzed
#  using grp
#
#
if(pr){
if(!CI)print('To get confidence intervals, set the argument CI=TRUE')
}
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK=J*K
mat=matrix(c(1:JK),J,K,byrow=TRUE)
if(POOL){
Y=list()
for(k in 1:K){
Y[[k]]=pool.a.list(x[mat[,k]])
}
A=DEP.PAIR.ES(Y,tr=tr,SEED=SEED,REL.MAG=REL.MAG,CI=CI)
}
if(!POOL){
if(OPT){
REL.EF=matrix(NA,nrow=4,ncol=3)
if(!is.null(REL.MAG))
if(length(REL.MAG)!=3)stop('REL.MAG should have three values')
if(SEED)set.seed(2)
z=rnorm(1000000)
v=dep.ES.summary.sub(z+REL.MAG[1],tr=tr)
REL.EF[,1]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[2],tr=tr)
REL.EF[,2]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[3],tr=tr)
REL.EF[,3]=v[,2]
}
A=list()
for(j in 1:J){
A[[j]]=DEP.PAIR.ES(x[mat[j,]],tr=tr,SEED=SEED,REL.MAG=REL.MAG,CI=CI)
}}
list(A=A)
}


bw.es.I<-function(J,K,x,tr=.2,OPT=FALSE,SEED=TRUE,CI=FALSE, alpha=.05, REL.MAG=NULL,pr=TRUE){
#
#  Between-by-within design.
#
# Effect size based on an interaction: compare difference scores
#
# If POOL=TRUE, pool the data over levels of Factor B
#
# If the data are stored in a matrix or data frame, it is converted to list mode.
#   Once in list,
# x[[1]] contains the x
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the x for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the x for level 1,K
#  x[[K+1]] is the x for level 2,1, x[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested, but a subset of the x can be analyzed
#  using grp
#
#
if(pr){
if(!CI)print('To get confidence intervals, set the argument CI=TRUE')
}
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK=J*K
mat=matrix(c(1:JK),J,K,byrow=TRUE)
Int=list()
mat=con2way(J,K)$conAB
NC=ncol(mat)
for(j in 1:NC){
id=which(mat[,j]!=0)
dif1=x[[id[1]]]-x[[id[2]]]
dif2=x[[id[3]]]-x[[id[4]]]
if(!CI)Int[[j]]=ES.summary(dif1,dif2,tr=tr)
if(CI)Int[[j]]=ES.summary.CI(dif1,dif2,tr=tr,alpha=alpha)
}
list(Interaction.ES=Int,con=mat)
}


ww.es<-function(J,K,x,tr=.2,CI=FALSE,SEED=TRUE,REL.MAG=NULL){
#
#  within-by-within design.
#
# For each level of Factor A,  compute  effect sizes
# for all j<m levels of  Factor B .
#
#
# If the data are stored in a matrix or data frame, it is converted to list mode.
#   Once in list,
# x[[1]] contains the x
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the x for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the x for level 1,K
#  x[[K+1]] is the x for level 2,1, x[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested
#
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK=J*K
mat=matrix(c(1:JK),J,K,byrow=TRUE)
if(!is.null(REL.MAG)){
REL.EF=matrix(NA,nrow=4,ncol=3)
if(length(REL.MAG)!=3)stop('REL.MAG should have three values')
if(SEED)set.seed(2)
z=rnorm(1000000)
v=dep.ES.summary.sub(z+REL.MAG[1],tr=tr)
REL.EF[,1]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[2],tr=tr)
REL.EF[,2]=v[,2]
v=dep.ES.summary.sub(z+REL.MAG[3],tr=tr)
REL.EF[,3]=v[,2]
}
A=list()
B=list()
AB=list()
for(j in 1:J)A[[j]]=DEP.PAIR.ES(x[mat[j,]],tr=.2,SEED=SEED,REL.MAG=REL.MAG,CI=CI)
for(k in 1:K)B[[k]]=DEP.PAIR.ES(x[mat[,k]],tr=.2,SEED=SEED,REL.MAG=REL.MAG,CI=CI)
MAT=con2way(J,K)$conAB
for(k in 1:ncol(MAT)){
id=which(MAT[,k]!=0)
dif1=x[[id[1]]]-x[[id[2]]]
dif2=x[[id[3]]]-x[[id[4]]]
if(!CI)AB[[k]]=ES.summary(dif1,dif2)
if(CI)AB[[k]]=ES.summary.CI(dif1,dif2)
}
list(Factor.A=A,Factor.B=B,INT=AB,conA=con.all.pairs(J),conB=con.all.pairs(K),AB=MAT)
}


deplin.ES.summary.CI<-function(x,con=NULL,tr=.2,REL.MAG=NULL,SEED=TRUE,nboot=1000){
#
#  For J dependent variables,
#  compute four measures of effect size based on a linear contrast of the J variables specified by the argument
#  con
#
#  Generalizes dep.ES.summary.CI
#  Example:
#  If x is a matrix with two columns and con=c(1,-1), get the same results as dep.ES.summary.CI
#
#  By default, do all pairwise comparisons
#
#  Measures of effect size:
#
#  AKP: robust standardized difference similar to  Cohen's d
#  QS:  Quantile shift based on the median of the distribution of difference scores,
#  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
#  SIGN:  P(X<Y), probability that for a random pair, the first is less than the second.
#
#  If there is  no effect, distribution of the  linear contrast is  symmetric about zero.
#  Under normality, suppose a shift of .2, .5 and .8 standard deviation of the difference score
#  is considered small, medium and large. The corresponding values for QS and SIGN are printed.
#
#
#  REL.MAG: suppose it is decided that AKP values .1, .3 and .5 are viewed as small, medium and large under normality Then
# REL.MAG=c(.1,.3,.5)  the default,
# means that the function will compute the corresponding values for the measures of effect size used here.
#
CON=list()
if(is.list(x))x=matl(x)
x=elimna(x)
n=length(x)
J=ncol(x)
if(is.null(con))con=con.all.pairs(J)
for(k in 1:ncol(con)){
X=con[1,k]*x[,1]
for(j in 2:J)X=X+con[j,k]*x[,j]
output=dep.ES.summary.CI(X,tr=tr,REL.MAG=REL.MAG,SEED=SEED,nboot=nboot)
CON[[k]]=output
}
list(con=con, output=CON)
}


DEP.PAIR.ES<-function(x,tr=.2,CI=FALSE,
REL.MAG=NULL,
SEED=TRUE){
#
#  J dependent groups
#  For each pair of groups, compute a collection of  effect sizes
#  based on difference scores
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
J=length(x)
conmat=con.all.pairs(J)
ES=list()
Ncon=ncol(conmat)
ic=0

grp.com=matrix(NA,Ncon,2)
for(j in 1:Ncon){
id1=which(conmat[,j]==1)
id2=which(conmat[,j]==-1)
grp.com[j,]=c(id1,id2)
if(!CI)ES[[j]]=dep.ES.summary(x[[id1]],y=x[[id2]],tr=tr,REL.MAG=REL.MAG,SEED=SEED)
if(CI)ES[[j]]=dep.ES.summary.CI(x[[id1]],y=x[[id2]],tr=tr,REL.MAG=REL.MAG,SEED=SEED)
}
list(groups.compared=grp.com,Effect.size.estimates=ES)
}


wwlin.es<-function(J,K,x,tr = 0.2, REL.MAG = NULL, SEED = TRUE, nboot = 1000){
#
#  # For within-by-within
#
# Effect sizes based on linear sum of the random variables.
# Simplest case, compute effect sizes on x-y, difference scores
#
con=con2way(J,K)
A=deplin.ES.summary.CI(x,con=con$conA,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
B=deplin.ES.summary.CI(x,con=con$conB,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
AB=deplin.ES.summary.CI(x,con=con$conAB,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
list(Factor.A=A,Factor.B=B,Interactions=AB)
}


wwwlin.es<-function(J,K,L,x,tr = 0.2, REL.MAG = NULL, SEED = TRUE, nboot = 1000){
#
# For within-by-within-by-within
#
# Effect sizes based on linear sum of the random variables.
# Simplest case, compute effect sizes based on x-y, difference scores
#
con=con3way(J,K,L)
A=deplin.ES.summary.CI(x,con=con$conA,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
B=deplin.ES.summary.CI(x,con=con$conB,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
C=deplin.ES.summary.CI(x,con=con$conC,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
AB=deplin.ES.summary.CI(x,con=con$conAB,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
AC=deplin.ES.summary.CI(x,con=con$conAC,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
BC=deplin.ES.summary.CI(x,con=con$conBC,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
ABC=deplin.ES.summary.CI(x,con=con$conABC,tr=tr,SEED=SEED,REL.MAG=REL.MAG,nboot=nboot)
list(Factor.A=A,Factor.B=B,Factor.C=C,Inter.AB=AB,Inte.AC=AC,Inter.BC=BC,Inter.ABC=ABC)
}


bwwA.es<-function(J,K,L,x,fun=KMS.ci,nboot=1000,...){
#
#  For every two levels  of Factor A, compute  KMS shift effect size  assuming 20% trim
#  and do this for each 
#  level of Factors B and C.
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
JKL=J*K*L
KL=K*L
id=matrix(c(1:JKL),ncol=KL,byrow=TRUE)
con=con.all.pairs(J)
A=list()
N=(J^2-J)/2
for(i in 1:N){
w=which(con[,i]!=0)
isel=id[w,]
A[[i]]=bwwA.es.sub(2,K,L,x[isel],fun=fun,nboot=nboot,...)
}
list(A=A,con=con)
}


bwwA.es.sub<-function(J,K,L,x,fun=QSci,nboot=1000,...){
#
#  Effect sizes for the between factor, computed
#  for each level of the within factors
#
if(J!=2)stop('Must have J=2')
JKL=J*K*L
KL=K*L
ic=0
LOW=1
UP=1+KL
id=matrix(c(1:JKL),ncol=KL,byrow=TRUE)
ES=matrix(NA,nrow=KL,ncol=7)
for(k in 1:K){
for(l in 1:L){
ic=ic+1
ES[ic,1]=k
ES[ic,2]=l
isel=id[,ic]
d=IND.PAIR.ES(x[id[,ic]],fun=fun,...)$effect.size
temp=c(d[[1]]$n1,d[[1]]$n2,d[[1]]$effect.size,d[[1]][4]$ci[1],d[[1]][4]$ci[2])
ES[ic,3:7]=temp
}}
dimnames(ES)=list(NULL,c('B.Level','C.Level','n1','n2','Effect.Size','ci.low','ci.up'))
ES
}


BEST.cell<-function(x,alpha=.05,LARGEST=TRUE,method='AC',p.crit=NULL,AUTO=FALSE,iter=2000,SEED=TRUE,pr=TRUE){
#
#  For a multinomial distribution, can a decision be made about
#  about which cell has the highest probability?
#
#  PV if specified, is a N by iter matrix of p-values that can be computed via best.cell.crit
#  N=number of cells
# x  Assumed to contain the  cell frequencies
#
if(pr)print('Confidence intervals are based on the critical p-values')
if(SEED)set.seed(2)
x=elimna(x)
n=sum(x)
NCELL=length(x)
NCm1=NCELL-1
xor=order(x,decreasing = LARGEST)
IND.pv=NA
ic=0
CI=matrix(NA,nrow=NCm1,ncol=2)
for(j in 2:NCELL){
ic=ic+1
IND.pv[ic]=cell.com.pv(x,xor[1],xor[j])
}
if(is.null(p.crit))p.crit=best.cell.crit(n,NCELL,LARGEST=LARGEST,iter=iter,AUTO=FALSE,SEED=SEED)
output=matrix(NA,nrow=NCm1,8)
output[,1]=rep(x[xor[1]]/n,NCm1)
output[,2]=xor[2:NCELL]
output[,3]=x[xor[2:NCELL]]/n
output[,4]=output[,1]-output[,3]
ic=0
for(j in 2:NCELL){
ic=ic+1
CI=cell.com(x,xor[1],xor[j],AUTO=AUTO,method=method,alpha=p.crit[ic])
output[ic,5:6]=CI$ci
}
output[,7]=IND.pv
output[,8]=p.crit
dimnames(output)=list(NULL,c('Largest.Est','CELL','Est','Dif','ci.low','ci.up','p.value','p.crit'))
flag=IND.pv<=p.crit
id=output[flag,2]
setClass('BIN',slots=c('Cell.with.largest.estimate','Larger.than','n','output'))
put=new('BIN',Cell.with.largest.estimate=xor[1],Larger.than=id,n=n,output=output)

if(!LARGEST){
dimnames(output)=list(NULL,c('Smallest.Est','CELL','Est','Dif','ci.low','ci.up','p.value','p.crit'))
setClass('BIN',slots=c('Cell.with.smallest.estimate','smaller.than','n','output'))
put=new('BIN',Cell.with.smallest.estimate=xor[1],smaller.than=id,n=n,output=output)
}
put
}


KMS.ES.M<-function(x,y){
#
# Computes the robust effect size using a simple generalization of the method in
# Kulinskaya, E., Morgenthaler, S. & Staudte, R. (2008).
# Meta Analysis: A guide to calibrating and combining statistical evidence  p. 177
# based	on an M-estimator and percentage bend variance
#Cohen d=.2, .5 .8 correspond to .1, .25 and .4') (KMS p. 180)

x<-elimna(x)
y<-elimna(y)
n1<-length(x)
n2<-length(y)
N=n1+n2
q=n1/N
s1sq=pbvar(x)
s2sq=pbvar(y)
t1=onestep(x)
t2=onestep(y)
top=q*s1sq+(1-q)*s2sq
bot=q*(1-q)
sigsq=top/bot  #  Quantity in brackets KMS p. 176 eq 21.1
varrho=s2sq/s1sq
d1=(t1-t2)/sqrt(sigsq)
list(effect.size=d1,Cohen.d.equiv=2*d1)
}


ES.summary.sub<-function(x,n1,n2){
id1=c(1:n1)
n1p=n1+1
N=n1+n2
id2=c(n1p:N)
a=ES.summary.CI(x[id1],x[id2],SEED=F)[,8]
}

