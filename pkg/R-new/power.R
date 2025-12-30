# power.R - Power Analysis Functions
# Part of WRS package refactoring
#
# Functions for statistical power analysis and sample size estimation.
# These functions estimate statistical power for various robust methods
# including trimmed means, ANOVA, ANCOVA, and chi-square tests.
#
# Dependencies: trimse, winvar, yuen, elimna, tsreg, pbcor, outmgvf, pbvar, ancmg1
# External packages: pwr (for power.chisq.test)

# ============================================================================
# Power Analysis Functions (10 functions)
# ============================================================================

epow<-function(x,y,pcor=FALSE,regfun=tsreg,corfun=pbcor,outkeep=FALSE,outfun=outmgvf,varfun=pbvar,op=TRUE){
#
# Estimate the explanatory power between x and y
#
xx<-elimna(cbind(x,y))
pval<-1
if(is.matrix(x))pval<-ncol(x)
pp<-pval+1
x<-xx[,1:pval]
y<-xx[,pp]
x<-as.matrix(x)
flag<-rep(TRUE,nrow(x))
temp<-regfun(x,y)
ip<-ncol(x)+1
yhat<-y-temp$res
if(!outkeep){
temp<-outfun(cbind(x,y),plotit=FALSE)$out.id
flag[temp]<-FALSE
}
epow1<-varfun(yhat[flag])/varfun(y[flag])
if(pcor)epow2<-cor(yhat[flag],y[flag])^2
if(!pcor)epow2<-corfun(yhat[flag],y[flag])$cor^2
if(op)est<-epow2
if(!op)est<-epow1
est
}

pow1<-function(n,Del,alpha){
#
#  Determine power of Student's T in the
#  one-sided, one-sample case where
#
#  n=sample size
#  Del=(mu0-mu1)/sigma
#  alpha=Type I error probability
#  mu0 is hypothesized value
#  mu1 is some non-null value for the mean.
#
Del<-abs(Del)
if(alpha<=0 || alpha>=1)stop("alpha must be between 0 and 1")
K11<-1-alpha
K5<-sqrt(n)*Del
#  Next, use the Kraemer-Paik (1979, Technometrics, 21, 357-360)
#  approximation of the noncentral T.
K6<-n-1
K14<-qt(K11,K6)
K7<-K14*sqrt(1+K5*K5/K6)
K8<-K5*sqrt(1+K14*K14/K6)
K9<-K7-K8
pow1<-1-pt(K9,K6)
pow1
}

pow2an<-function(x,y,ci=FALSE,plotit=TRUE,nboot=800){
#
# Do a power analysis when comparing the 20% trimmed
# means of two independent groups with the percentile
# bootstrap method.
#
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
lp<-NA
se<-yuen(x,y)$se
gval<-NA
dv<-seq(0,3.5*se,length=15)
for(i in 1:length(dv)){
gval[i]<-powest(x,y,dv[i],se)
}
if(!ci){
if(plotit){
plot(dv,gval,type="n",xlab="delta",ylab="power")
lines(dv,gval)
}}
if(ci){
print("Taking bootstrap samples. Please wait.")
datax <- matrix(sample(x, size = length(x) * nboot, replace = TRUE),
                nrow = nboot)
datay <- matrix(sample(y, size = length(y) * nboot, replace = TRUE),
                nrow = nboot)
pboot<-matrix(NA,ncol=15,nrow=nboot)
for(i in 1:nboot){
se<-yuen(datax[i,],datay[i,])$se
for(j in 1:length(dv)){
pboot[i,j]<-powest(x,y,dv[j],se)
}}
ll<-floor(.05*nboot+.5)
for(i in 1:15){
temp<-sort(pboot[,i])
lp[i]<-temp[ll]
}
plot(c(dv,dv),c(gval,lp),type="n",xlab="delta",ylab="power")
lines(dv,gval)
lines(dv,lp,lty=2)
}
list(delta=dv,power=gval,lowp=lp)
}

power.chisq.test<-function(w = NULL, N = NULL, df = NULL,
 sig.level = 0.05, power = NULL){
library(pwr)
res=pwr.chisq.test(w=w,N=N,df=df,sig.level=sig.level, power = power)
res
}

powest<-function(x=NA,y=NA,delta=0,se=NA,wv1=NA,wv2=NA,n1=NA,n2=NA){
#
# wv1 = Winsorized variance for group 1
# wv2 = Winsorized variance for group 2
#
# Only 20% trimming is allowed.
#
tr<-.2
if(is.na(se)){
if(is.na(wv1)){
h1 <- length(x) - 2 * floor(tr * length(x))
h2 <- length(y) - 2 * floor(tr * length(y))
q1 <- ((length(x) - 1) * winvar(x, tr))/(h1 * (h1 - 1))
q2 <- ((length(y) - 1) * winvar(y, tr))/(h2 * (h2 - 1))
}
if(!is.na(wv1)){
if(is.na(n1))stop("Need to specify sample size for group 1")
if(is.na(n2))stop("Need to specify sample size for group 2")
h1<-n1-2*floor(tr*n1)
h2<-n2-2*floor(tr*n2)
q1<-(n1-1)*wv1/(h1*(h1-1))
q2<-(n2-1)*wv2/(h2*(h2-1))
}
se<-sqrt(q1+q2)
}
ygam<-sqrt(2*.01155)*c(0:35)/8
pow<-c(500.0,540.0,607.0, 706.0, 804.0,981.0,1176.0,1402.0,1681.0, 2008.0,
   2353.0, 2769.0, 3191.0, 3646.0, 4124.0, 4617.0, 5101.0, 5630.0,
   6117.0, 6602.0, 7058.0, 7459.0, 7812.0, 8150.0, 8479.0, 8743.0,
   8984.0, 9168.0, 9332.0, 9490.0, 9607.0, 9700.0, 9782.0, 9839.0,
   9868.0)/10000
flag<-(delta==0 & se==0)
if(flag)powest<-.05
else{
chk<-floor(8*delta/se)+1
chk1<-chk+1
gval<-delta/se
d1<-(gval-(chk-1)/8)*8
if(chk > length(pow))powest<-1
if(chk == length(pow))pow[chk1]<-1
if(chk <= length(pow))
powest<-pow[chk]+d1*(pow[chk1]-pow[chk])
}
powest
}

powt1an<-function(x,ci=FALSE,plotit=TRUE,nboot=800){
#
# Do a power analysis for the one-sample case with 20% trimmed
# mean and when the percentile bootstrap is to be used to test
# hypoltheses.
#
x<-x[!is.na(x)]
lp<-NA
se<-trimse(x)
gval<-NA
dv<-seq(0,3.5*se,length=15)
for(i in 1:length(dv)){
gval[i]<-powest(x,rep(0,5),dv[i],se)
}
if(!ci){
if(plotit){
plot(dv,gval,type="n",xlab="delta",ylab="power")
lines(dv,gval)
}}
if(ci){
set.seed(2)
print("Taking bootstrap samples. Please wait.")
datax <- matrix(sample(x, size = length(x) * nboot, replace = TRUE),
                nrow = nboot)
pboot<-matrix(NA,nrow=nboot,ncol=length(dv))
for(i in 1:nboot){
se<-trimse(datax[i,])
for(j in 1:length(dv)){
pboot[i,j]<-powest(x,rep(0,5),dv[j],se)
}}
ll<-floor(.05*nboot+.5)
for(i in 1:15){
temp<-sort(pboot[,i])
lp[i]<-temp[ll]
}
plot(c(dv,dv),c(gval,lp),type="n",xlab="delta",ylab="power")
lines(dv,gval)
lines(dv,lp,lty=2)
}
list(delta=dv,power=gval,lowp=lp)
}

powt1est<-function(x,delta=0,ci=FALSE,nboot=800){
#
# Estimate power for a given value of delta
#
# Only 20% trimming is allowed.
#
temp1<-powest(x,rep(0,5),delta,se=trimse(x))
if(ci){
set.seed(2)
pboot<-NA
datay<-rep(0,5)
print("Taking bootstrap samples. Please wait.")
datax <- matrix(sample(x, size = length(x) * nboot, replace = TRUE
                        ), nrow = nboot)
for(i in 1:nboot) {
se <- trimse(datax[i,  ])
pboot[i] <- powest(x, rep(0,5), delta, se)
}
temp <- sort(pboot)
}
ll<-floor(0.05 * nboot + 0.5)
list(est.power=temp1,ci=temp[ll])
}

z.power<-function(n,alpha=.05,del=NULL,var=NULL){
 q=qnorm(1-alpha/2)
 sig=sqrt(var)
 p1=pnorm(0-q-(sqrt(n)*del)/sig)
 p2=1-pnorm(q-(sqrt(n)*del)/sig)
 p=p1+p2
 list(power=p)
 }
