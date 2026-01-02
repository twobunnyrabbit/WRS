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


#' Bianco and Yohai Robust Logistic Regression Estimator
#'
#' @description
#' Computes the robust logistic regression estimator proposed by Bianco and Yohai (1996).
#' This is a modified version of code by Christophe Croux, Gentiane Haesbroeck, and
#' Kristel Joossens. The estimator is resistant to outliers in the predictor space and
#' provides robust alternatives to maximum likelihood estimation.
#'
#' @param x0 An \eqn{n \times (p-1)} matrix of explanatory variables (predictors).
#' @param y An \eqn{n}-vector containing the binomial response (0 or 1).
#' @param initwml Logical. If \code{TRUE}, uses weighted ML estimator for initial values
#'   with weights derived from MCD estimator on predictors. If \code{FALSE} (default),
#'   uses classical ML fit. Set to \code{FALSE} when predictors contain binary variables.
#' @param const Tuning constant for the estimator. Default is 0.5.
#' @param kmax Maximum number of iterations before convergence. Default is 1000.
#' @param maxhalf Maximum number of step-halving iterations. Default is 10.
#'
#' @details
#' This function implements the Bianco and Yohai (1996) estimator for logistic regression,
#' which provides robustness against outliers in the predictor space. By default, an
#' intercept term is included and \eqn{p} parameters are estimated.
#'
#' The algorithm uses an iterative optimization process with step-halving to ensure
#' convergence. Initial values can be obtained either from weighted ML (using MCD-based
#' weights) or classical ML estimation.
#'
#' **Note**: When explanatory variables contain binary observations, it is recommended
#' to set \code{initwml = FALSE}.
#'
#' @return A list with components:
#' \item{convergence}{Logical indicating whether convergence was achieved.}
#' \item{coef}{Vector of parameter estimates (length \eqn{p}).}
#' \item{sterror}{Standard errors of the parameters (only if convergence achieved).}
#' \item{objective}{Value of the objective function at the minimum (if no convergence).}
#'
#' @references
#' Bianco, A.M. and Yohai, V.J. (1996). Robust estimation in the logistic regression
#' model. In Robust Statistics, Data Analysis, and Computer Intensive Methods
#' (ed. H. Rieder), 17-34. Springer.
#'
#' Croux, C. and Haesbroeck, G. (2003). Implementing the Bianco and Yohai estimator
#' for Logistic Regression. Computational Statistics and Data Analysis, 44, 273-295.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' x0 <- matrix(rnorm(100), ncol=1)
#' y <- as.numeric(runif(100) > 0.5)
#'
#' # Fit robust logistic regression
#' result <- wlogregv2(x0, y)
#' result$coef  # Parameter estimates
#' }
#'
#' @keywords internal
#' @export
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


#' Critical P-Values for Multinomial Cell Comparisons
#'
#' @description
#' Computes critical p-values for identifying the best (largest or smallest) cell
#' in multinomial data. Uses simulation to determine appropriate critical values
#' for controlling familywise error rate when identifying extreme cells.
#'
#' @param N Sample size (total number of observations).
#' @param ncell Number of cells (categories) in the multinomial distribution.
#' @param LARGEST Logical. If \code{TRUE} (default), identifies largest cell.
#'   If \code{FALSE}, identifies smallest cell.
#' @param iter Number of simulation iterations. Default is 1000.
#' @param alpha Familywise error rate. Default is 0.05.
#' @param SEED Logical. If \code{TRUE} (default), sets random seed to 2 for
#'   reproducibility.
#' @param AUTO Logical. Passed to internal helper function. Default is \code{FALSE}.
#'
#' @details
#' This function simulates from a multinomial distribution with equal cell probabilities
#' and determines critical p-values for testing whether the observed "best" cell
#' (largest or smallest) is significantly different from what would be expected by chance.
#'
#' The function uses optimization to find critical values that control the familywise
#' error rate at the specified alpha level when making comparisons involving the
#' extreme cell.
#'
#' @return A numeric vector of critical p-values for each of the \eqn{ncell - 1}
#'   comparisons involving the best cell.
#'
#' @seealso \code{\link{best.cell.sub}} (internal helper function)
#'
#' @keywords internal
#' @export
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


#' Identify Group with Largest Location Measure
#'
#' @description
#' Determines which group has the largest measure of location and tests whether
#' this identification is statistically reasonable using bootstrap methods.
#' Performs pairwise comparisons to validate the identified best group.
#'
#' @param x A matrix, data frame, or list of vectors, one for each group.
#'   If matrix or data frame, groups are defined by columns.
#' @param est The location estimator to use. Default is \code{tmean} (trimmed mean).
#' @param nboot Number of bootstrap samples. If \code{NA}, uses default value.
#' @param SEED Logical. If \code{TRUE} (default), sets random seed for reproducibility.
#' @param pr Logical. If \code{TRUE} (default), prints progress/results.
#' @param ... Additional arguments passed to the estimator function and comparison tests.
#'
#' @details
#' This function identifies which group has the largest location measure and then
#' performs statistical tests to determine if this identification is reasonable.
#' It uses bootstrap methods via \code{\link{linconpb}} to test contrasts comparing
#' the best group against all others.
#'
#' **Note**: This function appears to have dependencies on internal variables
#' that may not be properly defined in all contexts. Use with caution.
#'
#' @return A list with components:
#' \item{Best.Group}{Index of the group with the largest location estimate.}
#' \item{Est.}{Vector of location estimates for each group.}
#' \item{p.value}{Maximum p-value from the pairwise comparisons testing the best group.}
#'
#' @seealso \code{\link{linconpb}}, \code{\link{tmean}}, \code{\link{trimpb}}
#'
#' @keywords internal
#' @export
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
