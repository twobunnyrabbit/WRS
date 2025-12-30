# zzz-internal.R
# 
# Internal utility functions and edge cases for the WRS package
#
# This module contains the final remaining functions that were not 
# extracted to other modules during the refactoring process. These are
# typically:
#   - Internal helper functions with narrow use cases
#   - Specialized functions with unique dependencies
#   - Functions that don't fit cleanly into other modules
#
# All functions in this module are marked as internal (@keywords internal)
# to indicate they are not part of the primary user-facing API.
#
# Functions in this module (3 total):
#   1. wlogregv2       - Bianco and Yohai logistic regression estimator
#   2. best.cell.crit  - Critical p-values for multinomial cell comparisons  
#   3. bestPB.DO       - Identify group with largest location measure
#


# Function 1: wlogregv2
# Bianco and Yohai (1996) logistic regression estimator
# @keywords internal

wlogregv2<-function(x0,y,initwml=FALSE,const=0.5,kmax=1e3,maxhalf=10)
{
#  Computation of the estimator of Bianco and Yohai (1996) in logistic regression
#  -------------
# This is a slightly modified version of code due to
#  Christophe Croux, Gentiane Haesbroeck, and Kristel Joossens
# (Here initwml defaults to F
#
#  This program computes the estimator of Bianco and Yohai in
#  logistic regression. By default, an intercept term is included
#  and p parameters are estimated.
#
#  For more details we refer to
#     Croux, C., and Haesbroeck, G. (2003), ``Implementing the Bianco and Yohai
#     estimator for Logistic Regression'',
#     Computational Statistics and Data Analysis, 44, 273-295
#
#Input:
#-------
# x0= n x (p-1) matrix containing the explanatory variables;
# y= n-vector containing binomial response (0 or 1);
#
# initwml= logical value for selecting one of the two possible methods for computing
#          the initial value of the optimization process. If initwml=T (default), a
#          weighted ML estimator is computed with weights derived from the MCD estimator
#          computed on the explanatory variables. If initwml=F, a classical ML fit is perfomed.
#          When the explanatory variables contain binary observations, it is recommended
#          to set initwml to F or to modify the code of the algorithm to compute the weights
#          only on the continuous variables.
# const= tuning constant used in the computation of the estimator (default=0.5);
# kmax= maximum number of iterations before convergence (default=1000);
# maxhalf= max number of step-halving (default=10).
#
# Example:
# x0=matrix(rnorm(100,1))
# y0=numeric(runif(100)>0.5)
# BYlogreg(x0,y)
#
#Output:
#--------
# list with
# 1st component: T or F if convergence achieved or not
# 2nd component: value of the objective function at the minimum
# p next components: estimates for the parameters.
# p last components: standard errors of the parameters (if first component is T)

x0=as.matrix(x0)
#  n=nrow(x0)
  p=ncol(x0)+1
p0=p-1
  #Smallest value of the scale parameter before implosion
  sigmamin=1e-4

# eliminate any rows with missing values
zz=elimna(cbind(x,y))
x=as.matrix(zz[,1:p0])
y=zz[,p]
n=nrow(x)
#  x=as.matrix(cbind(rep(1,n),x0))
  x=as.matrix(cbind(rep(1,n),x))
  y=as.numeric(y)

  # Computation of the initial value of the optimization process
  if (initwml==TRUE)
  {
    hp=floor(n*(1-0.25))+1
    mcdx=cov.mcd(x0, quantile.used =hp,method="mcd")
    rdx=sqrt(mahalanobis(x0,center=mcdx$center,cov=mcdx$cov))
    vc=sqrt(qchisq(0.975,p-1))
    wrd=(rdx<=vc)
    gstart=glm(y~x0,family=binomial,subset=wrd)$coef
  }
else {gstart=glm(y~x0,family=binomial)$coef}
  sigmastart=1/sqrt(sum(gstart^2))
  xistart=gstart*sigmastart
  stscores=x %*% xistart
sigma1=sigmastart
  #Initial value for the objective function
  oldobj=mean(phiBY3(stscores/sigmastart,y,const))
  kstep=jhalf=1
  while ((kstep < kmax) & (jhalf<maxhalf)){
unisig <- function(sigma)
{ mean(phiBY3(stscores/sigma,y,const))}
optimsig=nlminb(sigma1,unisig,lower=0)
sigma1=optimsig$par
    if (sigma1<sigmamin) {print("Explosion");kstep=kmax
      } else {
      gamma1=xistart/sigma1
      scores=stscores/sigma1
      newobj=mean(phiBY3(scores,y,const))
      oldobj=newobj
      gradBY3=colMeans((derphiBY3(scores,y,const)%*%matrix(1,ncol=p))*x)
      h=-gradBY3+((gradBY3 %*% xistart) *xistart)
      finalstep=h/sqrt(sum(h^2))
      xi1=xistart+finalstep
      xi1=xi1/(sum(xi1^2))
      scores1=(x%*%xi1)/sigma1
      newobj=mean(phiBY3(scores1,y,const))

      ####stephalving
      hstep=jhalf=1
      while  ((jhalf <=maxhalf) & (newobj>oldobj)){
        hstep=hstep/2
        xi1=xistart+finalstep*hstep
        xi1=xi1/sqrt(sum(xi1^2))
        scores1=x%*%xi1/sigma1
        newobj=mean(phiBY3(scores1,y,const))
        jhalf=jhalf+1
      }
    CONV=FALSE
      if ((jhalf==maxhalf+1) & (newobj>oldobj)) {CONV=TRUE
        } else {
        jhalf=1
        xistart=xi1
        oldobj=newobj
        stscores=x%*% xi1
        kstep=kstep+1
      }
    }
  }

  if (kstep == kmax) {
CONV=FALSE #    print("No convergence")
    result=list(convergence=FALSE,objective=0,coef=t(rep(NA,p)))
    } else {
    gammaest=xistart/sigma1
    stander=sterby3(x0,y,const,gammaest)
    result=list(convergence=CONV,coef=t(gammaest),sterror=stander)
  }
  return(result)
}


# Function 2: best.cell.crit
# Compute critical p-values for multinomial cell comparisons
# @keywords internal

best.cell.crit<-function(N,ncell,LARGEST=TRUE,iter=1000,alpha=.05,SEED=TRUE,AUTO=FALSE){
#
#
# N sample size
# ncell number of cells
#
if(SEED)set.seed(2)
NCm1=ncell-1
pv=NA
a=rmultinom(iter,N,rep(1/ncell,ncell))
pv.mat=apply(a,2,best.cell.sub,AUTO=AUTO,LARGEST=LARGEST)
init=apply(pv.mat,1,qest,alpha)
pv.mat=t(pv.mat)  # For simplicity when using extant code related to this function
z=optim(0,anc.best.fun,init=init,iter=iter,rem=pv.mat,Jm1=NCm1,
alpha=alpha,method='Brent',lower=0,upper=1)
p.crit=z$par*init
p.crit
}


# Function 3: bestPB.DO
# Determine which group has largest measure of location
# @keywords internal

bestPB.DO<-function(x,est=tmean,nboot=NA,SEED=TRUE,pr=TRUE,...){
#
#  Determine whether it is reasonable to
#  decide which group has largest measure of location
#
#
if(is.matrix(x)||is.data.frame(x))x<-listm(x)
J=length(x)
e=lapply(x,tmean,tr)
e=pool.a.list(e)
id=which(e==max(e))
id=id[1]
J=J-1
x=Z
e=lapply(x,tmean,tr)
e=pool.a.list(e)
pv=NA
for(j in 1:J)pv[j]=trimpb(x[[j]],tr=tr,pr=FALSE)$p.value
pv=max(pv)
if(!dif){
e=lapply(x,tmean,tr)
e=pool.a.list(e)
CON=conCON(J,id)$conCON
a=linconpb(x,con=CON,est=est,nboot=nboot,SEED=SEED,pr=FALSE,...)
pv=max(a$output[,3])
}
list(Best.Group=id,Est.=e,p.value=pv)
}
