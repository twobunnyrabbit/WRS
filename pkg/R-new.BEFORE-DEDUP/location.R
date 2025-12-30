# WRS Package - Robust Location Estimators
# Methods for estimating central tendency and location
#
# This module contains:
#   - M-estimators: mest, mom, onestep
#   - Harrell-Davis quantile estimators: hd variants
#   - Trimmed means and variants
#   - Bivariate and multivariate location estimators
#   - Location comparison and difference functions
#
# See: Wilcox (2022), Introduction to Robust Estimation


# mest
mest<-function(x,bend=1.28,na.rm=FALSE){
#
#  Compute M-estimator of location using Huber's Psi.
#  The default bending constant is 1.28
#
if(na.rm)x<-x[!is.na(x)]
if(mad(x)==0)stop("MAD=0. The M-estimator cannot be computed.")
y<-(x-median(x))/mad(x)  #mad in splus is madn in the book.
A<-sum(hpsi(y,bend))
B<-length(x[abs(y)<=bend])
mest<-median(x)+mad(x)*A/B
repeat{
y<-(x-mest)/mad(x)
A<-sum(hpsi(y,bend))
B<-length(x[abs(y)<=bend])
newmest<-mest+mad(x)*A/B
if(abs(newmest-mest) <.0001)break
mest<-newmest
}
mest
}



# mestci
mestci<-function(x,alpha=.05,nboot=4000,bend=1.28,os=FALSE,pr=TRUE){
#
#   Compute a bootstrap, .95 confidence interval for the
#   M-estimator of location based on Huber's Psi.
#   The default percentage bend is bend=1.28
#   The default number of bootstrap samples is nboot=4000
#
#   By default, the fully iterated M-estimator is used. To use the
#   one-step M-estimator instead, set os=TRUE
#
os<-as.logical(os)
if(pr){
if(length(x) <=19)
print("The number of observations is less than 20.")
print("This function might fail due to division by zero,")
print("which in turn causes an error in function hpsi")
print("having to do with a missing value.")
}
set.seed(1) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
if(!os)bvec<-apply(data,1,mest,bend)
if(os)bvec<-apply(data,1,onestep,bend)
bvec<-sort(bvec)
low<-round((alpha/2)*nboot)
up<-nboot-low
low<-low+1
list(ci=c(bvec[low],bvec[up]))
}



# mestse
mestse<-function(x,bend=1.28,op=2){
#
#   Estimate the standard error of M-estimator using Huber's Psi
#   using estimate of influence function
#
n<-length(x)
mestse<-sqrt(sum((ifmest(x,bend,op=2)^2))/(n*(n-1)))
mestse
}


# mestseb
mestseb<-function(x,nboot=1000,bend=1.28,SEED=TRUE){
#
#   Compute bootstrap estimate of the standard error of the
#   M-estimators with Huber's Psi.
#   The default percentage bend is bend=1.28
#   The default number of bootstrap samples is nboot=100
#
if(SEED)set.seed(1) # set seed of random number generator so that
#   results can be duplicated.
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,mest,bend=bend)
mestseb<-sqrt(var(bvec))
mestseb
}


# mom
mom<-function(x,bend=2.24,na.rm=TRUE){
#
#  Compute MOM-estimator of location.
#  The default bending constant is 2.24
#
if(na.rm)x<-x[!is.na(x)] #Remove missing values
flag1<-(x>median(x)+bend*mad(x))
flag2<-(x<median(x)-bend*mad(x))
flag<-rep(TRUE,length(x))
flag[flag1]<-FALSE
flag[flag2]<-FALSE
mom<-mean(x[flag])
mom
}



# momci
momci<-function(x,alpha=.05,nboot=2000,bend=2.24,SEED=TRUE,null.value=0){
#
#   Compute a bootstrap, .95 confidence interval for the
#   MOM-estimator of location based on Huber's Psi.
#   The default number of bootstrap samples is nboot=500
#
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
est=mom(x,bend=bend)
bvec<-apply(data,1,mom,bend)
bvec<-sort(bvec)
low<-round((alpha/2)*nboot)
up<-nboot-low
low<-low+1
pv=mean(bvec>null.value)+.5*mean(bvec==null.value)
pv=2*min(c(pv,1-pv))
list(ci=c(bvec[low],bvec[up]),p.value=pv,est.mom=est)
}



# onestep
onestep<-function(x,bend=1.28,na.rm=FALSE,MED=TRUE){
#
#  Compute one-step M-estimator of location using Huber's Psi.
#  The default bending constant is 1.28
#
#  MED=TRUE: initial estimate is the median
#  Otherwise use modified one-step M-estimator
#
if(na.rm)x<-x[!is.na(x)]
if(MED)init.loc=median(x)
if(!MED)init.loc=mom(x,bend=bend)
y<-(x-init.loc)/mad(x)  #mad in splus is madn in the book.
A<-sum(hpsi(y,bend))
B<-length(x[abs(y)<=bend])
onestep<-median(x)+mad(x)*A/B
onestep
}



# tmean
tmean<-function(x,tr=.2,na.rm=FALSE,STAND=NULL){
if(na.rm)x<-x[!is.na(x)]
val<-mean(x,tr)
val
}


# hdci
hdci<-function(x,q=.5,alpha=.05,nboot=100,SEED=TRUE,pr=TRUE){
#
#   Compute a 1-alpha confidence for qth quantile using the
#   Harrell-Davis estimator in conjunction with the
#   bootstrap estimate of the standard error.
#
#   The default quantile is .5.
#   The default value for alpha is .05.
#
if(alpha!=.05)stop("Use the function qcipb. Generally works well even when alpha is not equal to .05")
x=elimna(x)
if(pr){
if(sum(duplicated(x)>0))print("Duplicate values detected; use hdpb")
}
se<-hdseb(x,q,nboot,SEED=SEED)
crit<-.5064/(length(x)^(.25))+1.96
if(q<=.2 || q>=.8){
if(length(x) <=20)crit<-(-6.23)/length(x)+5.01
}
if(q<=.1 || q>=.9){
if(length(x) <=40)crit<-36.2/length(x)+1.31
}
if(length(x)<=10){
print("The number of observations is less than 11.")
print("Accurate critical values have not been determined for this case.")
}
low<-hd(x,q)-crit*se
hi<-hd(x,q)+crit*se
list(ci=c(low,hi),crit=crit,se=se)
}


# hdpb
hdpb<-function(x,est=hd,alpha=.05,nboot=2000,SEED=TRUE,nv=0,...){
#
#   Compute a bootstrap, .95 confidence interval for the
#   measure of location corresponding to the argument est.
#   By default, the Harrell-Davis estimator is used
#
#   The default number of bootstrap samples is nboot=2000
#
#   The parameter q determines the quantile estimated via the function hd
#    This function is the same as onesampb, only for convenience it defaults
#   to using an estimate of the median.
#
#    nv=null value when  computing a p-value
#
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
x=elimna(x)
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,est,...)
bvec<-sort(bvec)
low<-round((alpha/2)*nboot)
up<-nboot-low
low<-low+1
pv=mean(bvec>nv)+.5*mean(bvec==nv)
pv=2*min(c(pv,1-pv))
estimate=est(x,...)
list(ci=c(bvec[low],bvec[up]),n=length(x),estimate=estimate,p.value=pv)
}


# hdseb
hdseb<-function(x,q=.5,nboot=100,SEED=TRUE){
#
#   Compute bootstrap estimate of the standard error of the
#   Harrell-Davis estimator of the qth quantile.
#   The default quantile is the median, q=.5
#   The default number of bootstrap samples is nboot=100
#
if(SEED)set.seed(2) # set seed of random number generator so that
#   results can be duplicated.
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,hd,q)
hdseb<-sqrt(var(bvec))
hdseb
}


# hdno
hdno<-function(x,q=.5){
#
# Use hd when .1<1<.9,
# otherwise use no estimator
#
x=elimna(x)
if(q<=.1 || q>=.9)e=qno.est(x,q)
else
e=hd(x,q)
e
}


# hdmq
hdmq<-function(x,q=.5,tr=FALSE){
#
#
# Estimate one or more quantiles.
e=NA
nq=length(q)
if(!tr)for(i in 1:nq)e[i]=hd(x,q[i])
if(tr)for(i in 1:nq)e[i]=thd(x,q[i])
e
}


# hdep
hdep <- function(PNT, X, NDIR=100, EPS=10E-7, SEED=NA, PRINT=F )
{
#========================
 #
 # X  - A numeric matrix with N rows and NP columns
 # PNT  - A numeric vector representing a point in the same space as
 # defined by X, so
 #     length of T has to equal to NP.
 # NDIR - A number of samples to draw
 # EPS  - Precision.
 # SEED - If specified, sets the seen of the random
 # number generator
 # PRINT - Default=F. If T, prints warning messages, such as Eigenvectors
 # are 0.
 #
 #=============================
  # SUBROUTINES
 #  DEP
 #  Reduce
 #

   DEP  <- function( X, PNT, NDIR, EPS=10E-8, PRINT=F )
   {
    #=================================================

    #
    # X  - A numeric matrix with N rows and NP columns
    # PNT  - A numeric vector representing a point in the same space as
    # defined by X, so
    #     length of T has to equal to NP.
    # NDIR - A number of samples to draw
    # EPS  - Precision.
    #
    #==================================
    # ---------------------------------------
    # Initialize Number of singular samples
    # ---------------------------------------
     NSIN <- 0
     N  <- nrow( X )
     NP  <- ncol( X )

    # ---------------------------------------
    # Intitialize Halfspace Depth at random
    # seed
    # ---------------------------------------
     NDEP <- N
     for( NRAN in 1:NDIR )
     {
       foundSingular <- F
      # --------------------------------------- ---
      # Draw a random sample of size NP without
      #    replacement
      # ------------------------------------------
       JSAMP <- sample( 1:N, size=NP, replace=FALSE )

      # ------------------------------------------
      #   Compute covariance matrix of the sample
      # ------------------------------------------
       sX  <- matrix( X[JSAMP,], nrow=NP, ncol=NP )
       COV  <- var( sX )

      # ------------------------------------------
      # Computing Eigen Values And Eigen Vector
      # for COV matrix
      # ------------------------------------------
       resEigen <- eigen( COV )
       Eval  <- resEigen[[1]]
       Evec  <- resEigen[[2]]

       if (Eval[NP] > EPS)
       {
        NSIN <- NSIN + 1
        foundSingular <- T
        if (PRINT)
         paste( "ERROR: No Eigenvalue = 0 for sample", NRAN)
        next
       }

      # ------------------------------------------
      # Need to test for singularity
      # ------------------------------------------
       if (Eval[NP-1] <= EPS)
       {
        NSIN <- NSIN + 1
       }

      # ------------------------------------------
      # Projecting all pints on line through
      # theta with direction given by the eigen
      # vector of the smallest eigenvalue, i.e.,
      # the direction orthogonal on the hyperplane
      # given by the NP-subset.
      # Compute the one-dimensional halfspace depth
      # of theta on this line.
      # ------------------------------------------
      # in Splus the smallest eigenvalue is the
      # last one and corresponding vector is the
      # last one, hence Eval[NP] is the smallest
      # and Evec[,NP] is the corresponding vector
      # ------------------------------------------
       eigenVec <- Evec[,NP]
       NT   <- sum( ifelse( eigenVec <= EPS, 1, 0 ) )
       KT   <- sum( ifelse( eigenVec > EPS, PNT * eigenVec, 0 ) )
       if (NT == NP)
       {
        NSIN <- NSIN + 1
        foundSingular <- T
        if (PRINT)
         paste( "      ERROR: Eigenvector = 0 for sample", NRAN )
        if (foundSingular) next             # Do next Sample
       }
       K  <- X %*% eigenVec
       K  <- K - KT
       NUMH <- sum( ifelse( K > EPS, 1, 0 ) )
       NT  <- sum( ifelse( abs(K) <= EPS, 1, 0 ) )
      # -------------------------------------------
      # If all projections collapse with theta,
      # return to reduce the dimension
      # -------------------------------------------
       if (NT == N)
       {
        NSIN <- -1
        return( list( NDEP=NDEP, NSIN=NSIN, EVEC=Evec ) ) # Will need
#Eigen Vector matrix to reduce dimension
       }

      # -------------------------------------------
      # Update halfspace depth
      # -------------------------------------------
       NDEP <- min( NDEP, min( NUMH+NT,N-NUMH ) )
     }

     return( list( NDEP=NDEP, NSIN=NSIN, EVEC=Evec ) )
   }

   #================================================
   Reduce <- function( X, PNT, Evec )
   {
    Det <- det(Evec)
    if (Det==0)
    {
     return( list( X=X, PNT=PNT, DET=Det ) )
    }
    NP <- ncol(X)

    # ---------------------------------------
    # Compute (NP-1)-dimentional coordinates
    # for all points and theta
    # ---------------------------------------
    RedEvec <- matrix(Evec[,1:(NP-1)],nrow=NP,ncol=(NP-1)) # Reducing
    #  dimension by removing the last dimension with 0 variance.
    PNT   <- PNT %*% RedEvec
    X   <- X %*% RedEvec
    if (!is.matrix(X)) X <- matrix(X,ncol=(NP-1))
    return( list( X=X, PNT=PNT, DET=Det ) )
   }

#
# PROGRAM BEGINS
#
  if (!is.na(SEED)) set.seed( SEED )
 # ---------------------------------------
 # Initialize Number of singular samples
 # ---------------------------------------
  Nsin <- 0

  X  <- as.matrix( X )
  N  <- nrow( X )
  NP  <- ncol( X )

if (length(PNT) != NP){print("Length of 'PNT' has to equal to")
stop("number of columns in X !!!   " )
}

 # ---------------------------------------
 # Handle special case where N=1
 # ---------------------------------------
  if (N==1)
  {
   NDEP <- ifelse( abs(X[1,]-PNT) > EPS, 0, 1 )  # if any dimension
#  different from point PNT, NDEP=0, else = 1
   NDEP <- min( NDEP )
   DEPTH <- NDEP/ N
   return( DEPTH )
  }

 # ---------------------------------------
 # Handle special case where NP=1
 # ---------------------------------------
 repeat #+++++++++++++++++++++++++++++++++
 {
 # In this case depth is equal to number of points <= to T
  if (NP==1)
  {
   MORE <- sum( ifelse( X[,1] >= (PNT-EPS), 1, 0 ) )
   LESS <- sum( ifelse( X[,1] <= (PNT+EPS), 1, 0 ) )
   NDEP <- min( LESS, MORE )
   DEPTH <- NDEP / N
   return( DEPTH )
  }

 # ---------------------------------------
 # General Case, call function DEP
 # ---------------------------------------
  if (N > NP)
  {
   RES  <- DEP( X=X, PNT=PNT, NDIR=NDIR, EPS=EPS, PRINT=PRINT )
   NDEP <- RES$NDEP
   NSIN <- RES$NSIN
   EVEC <- RES$EVEC
  }
  else
  {
   NSIN <- -1  # Needs to reduce dimensions
   EVEC <- eigen( var( X ) )[[2]]  # Getting eigenvector
  }

 # ---------------------------------------
 # If all points and theta are identified
 # as lying on the same hyperplane, reduce
 # the dimension of the data set by projection
 # on that hyperplane, and compute the depth
 # on the reduced data set
 # ---------------------------------------
  if (NSIN == -1)
  {
   NSIN <- 0
   if (PRINT) print( "      Direction with zero variance detected" )
   RED  <- Reduce( X=X, PNT=PNT, Evec=EVEC )
   X  <- RED$X
   PNT  <- RED$PNT
   Det  <- RED$DET
   if (Det==0)
   {
print("\n\n\t DIMENSION REDUCTION TERMINATED\n\t EIGENVECTORS ARE NOT")
stop("INDEPENDENT\n\n" )
   }
   NP  <- ncol(X)
   if (PRINT) paste("     Dimension reduced to", NP )
  }
  else
  {
   break   # No need to reduce dimension of X and hence no need to
#return, breaks 'repeat' loop
  }
 } # End repeat+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


 DEPTH <- NDEP / N
 return( DEPTH )
}




# IQRhd
IQRhd<-function(x){
e=hd(x,.75)-hd(x,.25)
e
}


# medhd2g
medhd2g<-function(x, y, alpha = 0.05, nboot = 2000,SEED=TRUE,pr=TRUE, ...){
#
# Compare medians via the Harrell-Davis estimator
#
res=pb2gen(x,y,alpha=alpha,nboot=2000,est=hd,SEED=SEED,pr=pr, ...)
res
}


# Dqcihd
Dqcihd<-function(x,y,alpha=.05,q=c(1:9/10),
plotit=FALSE,pop=0,fr=.8,rval=15,xlab="",ylab="",nboot=600,SEED=TRUE){
#
#  Compute a confidence interval for the quantiles for D=X-Y, X and Y independent.
#
#  The Harrell-Davis estimator is used
#  If the distribution of X and Y are identical, then in particular the
#  distribution of D=X-Y is symmetric about zero.
#
#  plotit=TRUE causes a plot of the difference scores to be created
#  pop=0 adaptive kernel density estimate
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate (Rosenblatt's shifted histogram)
#  pop=3 boxplot
#  pop=4 stem-and-leaf
#  pop=5 histogram
#
if(SEED)set.seed(2)
x<-x[!is.na(x)]
y<-y[!is.na(y)]
n1=length(x)
n2=length(y)
m<-outer(x,y,FUN="-")
est=NA
for(i in 1:length(q))est[i]=hd(m,q=q[i])
data1<-matrix(sample(n1,size=n1*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n2,size=n2*nboot,replace=TRUE),nrow=nboot)
bvec=matrix(NA,nrow=nboot,ncol=length(q))
for(i in 1:nboot){
mb=outer(x[data1[i,]],y[data2[i,]],"-")
for(j in 1:length(q))
bvec[i,j]=hd(mb,q=q[j])
}
p=NA
ci=matrix(NA,nrow=length(q),ncol=2)
for(j in 1:length(q)){
p[j]=mean(bvec[,j]>0)+.5*mean(bvec[,j]==0)
p[j]=2*min(c(p[j],1-p[j]))
sbv=sort(bvec[,j])
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci[j,1]=sbv[ilow]
ci[j,2]=sbv[ihi]
}
if(plotit){
if(pop==1 || pop==0){
if(length(x)*length(y)>2500){
print("Product of sample sizes exceeds 2500.")
print("Execution time might be high when using pop=0 or 1")
print("If this is case, might consider changing the argument pop")
print("pop=2 might be better")
}}
MM=as.vector(m)
if(pop==0)akerd(MM,xlab=xlab,ylab=ylab)
if(pop==1)rdplot(MM,fr=fr,xlab=xlab,ylab=ylab)
if(pop==2)kdplot(MM,rval=rval,xlab=xlab,ylab=ylab)
if(pop==3)boxplot(MM)
if(pop==4)stem(MM)
if(pop==5)hist(MM,xlab=xlab)
if(pop==6)skerd(MM)
}
output=cbind(as.matrix(q),as.matrix(est),ci,as.matrix(p))
dimnames(output)=list(NULL,c("Quantile","Estimates","ci.low","ci.up","p-value"))
output
}

# Dqcomhd
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




# cbmhd
cbmhd<-function(x,y,qest=hd,alpha=.05,q=.25,plotit=FALSE,pop=0,fr=.8,rval=15,xlab="",ylab="",nboot=600,SEED=TRUE){
#
#  Compute a confidence interval for the sum of the qth and (1-q)th quantiles
#  of the distribution of D=X-Y, where X and Y are two independent random variables.
#  The Harrell-Davis estimator is used
#  If the distribution of X and Y are identical, then in particular the
#  distribution of D=X-Y is symmetric about zero.
#
#  plotit=TRUE causes a plot of the difference scores to be created
#  pop=0 adaptive kernel density estimate
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate (Rosenblatt's shifted histogram)
#  pop=3 boxplot
#  pop=4 stem-and-leaf
#  pop=5 histogram
#
if(SEED)set.seed(2)
if(q>=.5)stop("q should be less than .5")
if(q<=0)stop("q should be greater than 0")
x<-x[!is.na(x)]
y<-y[!is.na(y)]
n1=length(x)
n2=length(y)
m<-outer(x,y,FUN="-")
q2=1-q
est1=qest(m,q)
est2=qest(m,q2)
data1<-matrix(sample(n1,size=n1*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n2,size=n2*nboot,replace=TRUE),nrow=nboot)
bvec=NA
for(i in 1:nboot){
mb=outer(x[data1[i,]],y[data2[i,]],"-")
bvec[i]=qest(mb,q)+qest(mb,q2)
}
p=mean(bvec>0)+.5*mean(bvec==0)
p=2*min(c(p,1-p))
sbv=sort(bvec)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=sbv[ilow]
ci[2]=sbv[ihi]
if(plotit){
if(pop==1 || pop==0){
if(length(x)*length(y)>2500){
print("Product of sample sizes exceeds 2500.")
print("Execution time might be high when using pop=0 or 1")
print("If this is case, might consider changing the argument pop")
print("pop=2 might be better")
}}
MM=as.vector(m)
if(pop==0)akerd(MM,xlab=xlab,ylab=ylab)
if(pop==1)rdplot(MM,fr=fr,xlab=xlab,ylab=ylab)
if(pop==2)kdplot(MM,rval=rval,xlab=xlab,ylab=ylab)
if(pop==3)boxplot(MM)
if(pop==4)stem(MM)
if(pop==5)hist(MM,xlab=xlab)
if(pop==6)skerd(MM)
}
list(q=q,Est1=est1,Est2=est2,sum=est1+est2,ci=ci,p.value=p)
}


# cbmhdMC
cbmhdMC<-function(x,y,alpha=.05,q=.25,plotit=FALSE,pop=0,fr=.8,rval=15,xlab="",ylab="",nboot=600,SEED=TRUE){
#
#  Compute a confidence interval for the sum of the qth and (1-q)th quantiles
#  of the distribution of D=X-Y, where X and Y are two independent random variables.
#  The Harrell-Davis estimator is used
#  If the distribution of X and Y are identical, then in particular the
#  distribution of D=X-Y is symmetric about zero.
#
#  plotit=TRUE causes a plot of the difference scores to be created
#  pop=0 adaptive kernel density estimate
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate (Rosenblatt's shifted histogram)
#  pop=3 boxplot
#  pop=4 stem-and-leaf
#  pop=5 histogram
#
if(SEED)set.seed(2)
if(q>=.5)stop("q should be less than .5")
if(q<=0)stop("q should be greater than 0")
x<-x[!is.na(x)]
y<-y[!is.na(y)]
n1=length(x)
n2=length(y)
m<-outer(x,y,FUN="-")
q2=1-q
est1=hd(m,q)
est2=hd(m,q2)
data1<-matrix(sample(n1,size=n1*nboot,replace=TRUE),nrow=nboot)
data2<-matrix(sample(n2,size=n2*nboot,replace=TRUE),nrow=nboot)
data=cbind(data1,data2)
data=listm(t(data))
bvec=NA
bvec<-mclapply(data,cbmhd_subMC,x=x,y=y,q=q,q2=q2,n1=n1,n2=n2,mc.preschedule=TRUE)
bvec=list2vec(bvec)
p=mean(bvec>0)+.5*mean(bvec==0)
p=2*min(c(p,1-p))
sbv=sort(bvec)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=sbv[ilow]
ci[2]=sbv[ihi]
if(plotit){
if(pop==1 || pop==0){
if(length(x)*length(y)>2500){
print("Product of sample sizes exceeds 2500.")
print("Execution time might be high when using pop=0 or 1")
print("If this is case, might consider changing the argument pop")
print("pop=2 might be better")
}}
MM=as.vector(m)
if(pop==0)akerd(MM,xlab=xlab,ylab=ylab)
if(pop==1)rdplot(MM,fr=fr,xlab=xlab,ylab=ylab)
if(pop==2)kdplot(MM,rval=rval,xlab=xlab,ylab=ylab)
if(pop==3)boxplot(MM)
if(pop==4)stem(MM)
if(pop==5)hist(MM,xlab=xlab)
if(pop==6)skerd(MM)
}
list(q=q,Est1=est1,Est2=est2,sum=est1+est2,ci=ci,p.value=p)
}


# Dcbmhd
Dcbmhd<-function(x=NULL,y=NULL,d=NULL,qest=hd,alpha=.05,q=.25,plotit=FALSE,pop=0,
fr=.8,rval=15,xlab='',ylab='',nboot=600,SEED=TRUE){
#
#
#  Compute a confidence interval for the sum of the qth and (1-q)th quantiles
#  of the distribution of D=X-Y, where X and Y are two
#  dependent random variables.
#  The Harrell-Davis estimator is used
#  If the distribution of X and Y are identical, then in particular the
#  distribution of D=X-Y is symmetric about zero.
#
#  plotit=TRUE causes a plot of the difference scores to be created
#  pop=0 adaptive kernel density estimate
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate (Rosenblatt's shifted histogram)
#  pop=3 boxplot
#  pop=4 stem-and-leaf
#  pop=5 histogram
#
if(SEED)set.seed(2)
if(q>=.5)stop('q should be less than .5')
if(q<=0)stop('q should be greater than 0')
if(is.null(d))d=x-y
if(is.null(d))stop('Apparently x or y contain no data')
d=elimna(d)
n=length(d)
q2=1-q
est1=qest(d,q)
est2=qest(d,q2)
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
bvec=NA
for(i in 1:nboot){
bvec[i]=qest(d[data[i,]],q)+qest(d[data[i,]],q2)
}
p=mean(bvec>0)+.5*mean(bvec==0)
p=2*min(c(p,1-p))
sbv=sort(bvec)
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
ci=sbv[ilow]
ci[2]=sbv[ihi]
if(plotit){
if(pop==1 || pop==0){
if(length(x)*length(y)>2500){
print('Product of sample sizes exceeds 2500.')
print('Execution time might be high when using pop=0 or 1')
print('If this is case, might consider changing the argument pop')
print('pop=2 might be better')
}}
MM=d
if(pop==0)akerd(MM,xlab=xlab,ylab=ylab)
if(pop==1)rdplot(MM,fr=fr,xlab=xlab,ylab=ylab)
if(pop==2)kdplot(MM,rval=rval,xlab=xlab,ylab=ylab)
if(pop==3)boxplot(MM)
if(pop==4)stem(MM)
if(pop==5)hist(MM,xlab=xlab)
if(pop==6)skerd(MM)
}
list(q=q,n=n,Est1=est1,Est2=est2,sum=est1+est2,ci=ci,p.value=p)
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


# ghtrim
ghtrim<-function(tr=.2,g=0,h=0){
#
#  Compute trimmed mean of a g-and-h distribution.
#
#
if(g==0)val=0
if(g>0){
low=qnorm(tr)
up=-1*low
val=integrate(ftrim,low,up,tr=tr,g=g,h=h)$value
val=val/(1-2*tr)
}
val
}


# ftrim
ftrim<-function(z,tr,g,h){
gz=(exp(g*z)-1)*exp(h*z^2/2)/g
res=dnorm(z)*gz
res
}


# btrim
btrim<-function(x,tr=.2,grp=NA,g=NULL,dp=NULL,nboot=599,SEED=TRUE){
#
#   Test the hypothesis of equal trimmed means, corresponding to J independent
#   groups, using a bootstrap-t method.
#
#   The data are assumed to be stored in x in list mode
#   or in a matrix. In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, columns correspond to groups.
#
#   grp is used to specify some subset of the groups, if desired.
#   By default, all J groups are used.
#   g=NULL, x is assumed to be a matrix or have list mode
#
#   if g is specifed, it is assumed that column g of x is
#   a factor variable and that the dependent variable of interest is in column
#   dp of x, which can be a matrix or data frame.
#
#   The default number of bootstrap samples is nboot=599
#
if(!is.null(g)){
if(is.null(dp))stop("Specify a value for dp, the column containing the data")
x=fac2list(x[,dp],x[,g])
}
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
if(is.na(grp[1]))grp<-c(1:length(x))
J<-length(grp)
nval=NA
x=lapply(x,elimna)
nval=lapply(x,length)
xbar=lapply(x,mean,tr=tr)
bvec<-array(0,c(J,2,nboot))
hval<-vector("numeric",J)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
for(j in 1:J){
hval[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
   # hval is the number of observations in the jth group after trimming.
print(paste("Working on group ",grp[j]))
xcen<-x[[grp[j]]]-mean(x[[grp[j]]],tr)
data<-matrix(sample(xcen,size=length(x[[grp[j]]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,,]<-apply(data,1,trimparts,tr) # A 2 by nboot matrix. The first row
#                     contains the bootstrap trimmed means, the second row
#                     contains the bootstrap squared standard errors.
}
m1<-bvec[,1,]  # J by nboot matrix containing the bootstrap trimmed means
m2<-bvec[,2,]  # J by nboot matrix containing the bootstrap sq standard errors
wvec<-1/m2  # J by nboot matrix of w values
uval<-apply(wvec,2,sum)  # Vector having length nboot
blob<-wvec*m1
xtil<-apply(blob,2,sum)/uval # nboot vector of xtil values
blob1<-matrix(0,J,nboot)
for (j in 1:J)blob1[j,]<-wvec[j,]*(m1[j,]-xtil)^2
avec<-apply(blob1,2,sum)/(length(x)-1)
blob2<-(1-wvec/uval)^2/(hval-1)
cvec<-apply(blob2,2,sum)
cvec<-2*(length(x)-2)*cvec/(length(x)^2-1)
testb<-avec/(cvec+1)
#            A vector of length nboot containing bootstrap test values
ct<-sum(is.na(testb))
if(ct>0)print("Some bootstrap estimates of the test statistic could not be computed")
test<-t1way(x,tr=tr,grp=grp)
pval<-sum(test$TEST<=testb)/nboot
#
# Determine explanatory effect size
#
e.pow=t1wayv2(x)$Effect.Size
list(test=test$TEST,p.value=pval,Explanatory.Power=e.pow,
Effect.Size=e.pow)
}



# bbtrim
bbtrim<-function(J,K,x,grp=c(1:p),p=J*K,tr=.2,nboot=600,alpha=.05,pr=FALSE){
#
#  Perform a J by K anova using trimmed means with
#  for independent groups using a bootstrap-t method
#
#  tr=.2 is default trimming
#
#
#  The R variable x is assumed to contain the raw
#  data stored in list mode or a matrix with JK columns.
#  If in list mode, x[[1]] contains the data
#  for the first level of both factors: level 1,1.
#  data[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the data for level 1,K
#  x[[K+1]] is the data for level 2,1, x[2K] is level 2,K, etc.
#
#  It is assumed that data has length JK, the total number of
#  groups being tested, but a subset of the data can be analyzed
#  using grp
#
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
Factor.A<-linconb(x,con=conA,tr=tr,alpha=alpha,nboot=nboot,pr=pr)
Factor.B<-linconb(x,con=conB,tr=tr,alpha=alpha,nboot=nboot,pr=pr)
Factor.AB<-linconb(x,con=conAB,tr=tr,alpha=alpha,nboot=nboot,pr=pr)
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.AB=Factor.AB,pr=pr)
}


# bwtrim
bwtrim<-function(J,K,data,tr=.2,grp=c(1:p),p=J*K,MAT=FALSE,grpc=1,coln=c(2:3)){
#  Perform a J-by-K anova on trimmed means with
#  repeated measures on the second factor. That is, a split-plot design
#  is assumed, with the first factor consisting of independent groups.
#
# If the data are stored in a matrix or data frame, it is converted to list mode.
#   Once in list mode,
#   data[[1]] contains the data
#  for the first level of both factors: level 1,1.
#  data[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  data[[K]] is the data for level 1,K
#  data[[K+1]] is the data for level 2,1, data[2K] is level 2,K, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that data has length JK, the total number of
#  groups being tested, but a subset of the data can be analyzed
#  using grp
#
#  If the between groups are denoted by groups numbers stored in a column
# of dat, you can set MAT=T, which will store the data in the format
# expected by this function
#
#  Example, grpc=1 means group id numbers are in col 1.
#  coln=c(3:6) means the within variables are stored in col 3-6.
#
#  Or you can use the function selbybw to sort the data.
#
if(is.data.frame(data))data=as.matrix(data)
if(MAT)
data=selbybw(data,grpc=grpc,coln=coln)$x
x<-data
       if(is.matrix(x) || is.data.frame(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                data <- y
        }
if(!is.list(data))stop("Data are not stored in list mode or a matrix")
if(p!=length(data)){
print("The total number of groups, based on the specified levels, is")
print(p)
print("The number of groups in data is")
print(length(data))
print("Warning: These two values are not equal")
}
if(p!=length(grp))stop("Apparently a subset of the groups was specified that does not match the total number of groups indicated by the values for J and K.")
tmeans<-0
h<-0
v<-matrix(0,p,p)
klow<-1-K
kup<-0
for (i in 1:p)tmeans[i]<-mean(data[[grp[i]]],tr,na.rm=TRUE)
for (j in 1:J){
h[j]<-length(data[[grp[j]]])-2*floor(tr*length(data[[grp[j]]]))
#    h is the effective sample size for the jth level of factor A
#   Use covmtrim to determine blocks of squared standard errors and
#   covariances.
klow<-klow+K
kup<-kup+K
sel<-c(klow:kup)
v[sel,sel]<-covmtrim(data[grp[klow:kup]],tr)
}
ij<-matrix(c(rep(1,J)),1,J)
ik<-matrix(c(rep(1,K)),1,K)
jm1<-J-1
cj<-diag(1,jm1,J)
for (i in 1:jm1)cj[i,i+1]<-0-1
km1<-K-1
ck<-diag(1,km1,K)
for (i in 1:km1)ck[i,i+1]<-0-1
#  Do test for factor A
cmat<-kron(cj,ik)  # Contrast matrix for factor A
Qa<-johansp(cmat,tmeans,v,h,J,K)
# Do test for factor B
cmat<-kron(ij,ck)  # Contrast matrix for factor B
Qb<-johansp(cmat,tmeans,v,h,J,K)
# Do test for factor A by B interaction
cmat<-kron(cj,ck)  # Contrast matrix for factor A by B
Qab<-johansp(cmat,tmeans,v,h,J,K)
list(Qa=Qa$teststat,Qa.p.value=Qa$p.value,
Qb=Qb$teststat,Qb.p.value=Qb$p.value,
Qab=Qab$teststat,Qab.p.value=Qab$p.value)
}



# bbwtrim
bbwtrim<-function(J,K,L,data,tr=.2,alpha=.05,p=J*K*L){
#  Perform a between-between-within (three-way) anova on trimmed means where
#
#  JK independent groups, L dependent groups
#
#  The variable data is assumed to contain the raw
#  data stored in list mode. data[[1]] contains the data
#  for the first level of all three factors: level 1,1,1.
#  data][2]] is assumed to contain the data for level 1 of the
#  first two factors and level 2 of the third factor: level 1,1,2
#  data[[L]] is the data for level 1,1,L
#  data[[L+1]] is the data for level 1,2,1. data[[2L]] is level 1,2,L.
#  data[[KL+1]] is level 2,1,1, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that data has length JKL, the total number of
#  groups being tested.
#
if(is.data.frame(data)) data <- as.matrix(data)
if(is.list(data))data=bbwna(J,K,L,data)
if(is.matrix(data))data=bbwmatna(J,K,L,data)
grp=c(1:p)
data=bbwna(J,K,L,data)
if(!is.list(data))stop("Data are not stored in list mode")
if(p!=length(data)){
print("The total number of groups, based on the specified levels, is")
print(p)
print("The number of groups is")
print(length(data))
print("Warning: These two values are not equal")
}
tmeans<-0
h<-0
v<-0
for (i in 1:p){
tmeans[i]<-mean(data[[grp[i]]],tr)
h[i]<-length(data[[grp[i]]])-2*floor(tr*length(data[[grp[i]]]))
#    h is the effective sample size
}
v=bbwcovm(J,K,L,data,tr=tr)
ij<-matrix(c(rep(1,J)),1,J)
ik<-matrix(c(rep(1,K)),1,K)
il<-matrix(c(rep(1,L)),1,L)
jm1<-J-1
cj<-diag(1,jm1,J)
cj<-diag(1,jm1,J)
for (i in 1:jm1)cj[i,i+1]<-0-1
km1<-K-1
ck<-diag(1,km1,K)
for (i in 1:km1)ck[i,i+1]<-0-1
lm1<-L-1
cl<-diag(1,lm1,L)
for (i in 1:lm1)cl[i,i+1]<-0-1
#  Do test for factor A
cmat<-kron(cj,kron(ik,il))  # Contrast matrix for factor A
Qa=bwwtrim.sub(cmat, tmeans, v, h,p)
Qa.siglevel <- 1 - pf(Qa, J - 1, 999)
# Do test for factor B
cmat<-kron(ij,kron(ck,il))  # Contrast matrix for factor B
Qb=bwwtrim.sub(cmat, tmeans, v, h,p)
 Qb.siglevel <- 1 - pf(Qb, K - 1, 999)
# Do test for factor C
cmat<-kron(ij,kron(ik,cl))  # Contrast matrix for factor C
Qc<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qc.siglevel <- 1 - pf(Qc, L - 1, 999)
# Do test for factor A by B interaction
cmat<-kron(cj,kron(ck,il))  # Contrast matrix for factor A by B
Qab<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qab.siglevel <- 1 - pf(Qab, (J - 1) * (K - 1), 999)
# Do test for factor A by C interaction
cmat<-kron(cj,kron(ik,cl))  # Contrast matrix for factor A by C
Qac<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qac.siglevel <- 1 - pf(Qac, (J - 1) * (L - 1), 999)
# Do test for factor B by C interaction
cmat<-kron(ij,kron(ck,cl))  # Contrast matrix for factor B by C
Qbc<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qbc.siglevel <- 1 - pf(Qbc, (K - 1) * (L - 1), 999)
# Do test for factor A by B by C interaction
cmat<-kron(cj,kron(ck,cl))  # Contrast matrix for factor A by B by C
Qabc<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qabc.siglevel <-1-pf(Qabc,(J-1)*(K-1)*(L-1), 999)
list(Qa=Qa,Qa.p.value=Qa.siglevel,Qb=Qb,Qb.p.value=Qb.siglevel,
Qc=Qc,Qc.p.value=Qc.siglevel,Qab=Qab,Qab.p.value=Qab.siglevel,
Qac=Qac,Qac.p.value=Qac.siglevel,Qbc=Qbc,Qbc.p.value=Qbc.siglevel,
Qabc=Qabc,Qabc.p.value=Qabc.siglevel)
}



# bwtrimbt
bwtrimbt<-function(J,K,x,tr=.2,alpha=.05,JK=J*K,grp=c(1:JK),nboot=599,
SEED=TRUE,monitor=FALSE){
#
# A bootstrap-t for performing a split-plot design
# with trimmed means.
# By default, 20% trimming is used with B=599 bootstrap samples.
#
#
#  The R variable x is assumed to contain the raw
#  data stored in list mode or in a matrix.
#  If in list mode, x[[1]] contains the data
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second: level 1,2
#  x[[K]] is the data for level 1,K
#  x[[K+1]] is the data for level 2,1, x[[2K]] is level 2,K, etc.
#
#  If the data are in a matrix, column 1 is assumed to
#  correspond to x[[1]], column 2 to x[[2]], etc.
#
#  When in list mode x is assumed to have length JK, the total number of
#  groups being tested, but a subset of the data can be analyzed
#  using grp
#
if(SEED)set.seed(2)
if(is.matrix(x)) {
y <- list()
ik=0
il=c(1:K)-K
for(j in 1:J){
il=il+K
zz=x[,il]
zz=elimna(zz)
for(k in 1:K){
ik=ik+1
y[[ik]]=zz[,k]
}}
                x <- y
}
JK<-J*K
data<-list()
xcen<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
xcen[[j]]<-data[[j]]-mean(data[[j]],tr) # Centered data
}
x<-data
#
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
# Next determine the n_j values
nvec<-NA
jp<-1-K
for(j in 1:J){
jp<-jp+K
nvec[j]<-length(x[[j]])
}
blist<-list()
print("Taking bootstrap samples. Please wait.")
testmat<-matrix(NA,ncol=3,nrow=nboot)
for(iboot in 1:nboot){
iv<-0
for(j in 1:J){
temp<-sample(nvec[j],replace = T)
for(k in 1:K){
iv<-iv+1
tempx<-xcen[[iv]]
blist[[iv]]<-tempx[temp]
}}
if(monitor)print(paste("Bootstrap iteration" ,iboot, "is complete"))
btest<-tsplit(J,K,blist,tr)
testmat[iboot,1]<-btest$Qa
testmat[iboot,2]<-btest$Qb
testmat[iboot,3]<-btest$Qab
}
lcrit<-round((1-alpha)*nboot)
temp3<-sort(testmat[,1])
crit.Qa<-temp3[lcrit]
temp3<-sort(testmat[,2])
crit.Qb<-temp3[lcrit]
temp3<-sort(testmat[,3])
crit.Qab<-temp3[lcrit]
temp4<-tsplit(J,K,x,tr=tr)
list(Qa=temp4$Qa,Qb=temp4$Qb,Qab=temp4$Qab,crit.Qa=crit.Qa,crit.Qb=crit.Qb,crit.Qab=crit.Qab)
}

# bbwtrimbt
bbwtrimbt<-function(J, K, L, x, tr = 0.2, JKL = J * K*L, con = 0,
 alpha = 0.05, grp =c(1:JKL), nboot = 599, SEED = TRUE, ...)
{
        #
        # Bootstrap-t for omniubs tests associated with
        # all main effects and interactions
        # for a between-by-between-within design.
        #
        #  The R variable x is assumed to contain the raw
        #  data stored in list mode or in a matrix.
        #  If in list mode, x[[1]] contains the data
        #  for the first level for all three factors: levels 1,1,1.
        #  x[[2]] is assumed to contain the data for level 1 of the
        #  first factor and second factor and level 2 of the third: level 1,1,2
        #  x[[K]] is the data for level 1,1,L
        #  x[[L+1]] is the data for level 1,2,1, x[[2K]] is level 1,2,2, etc.
        #
        #  If the data are in a matrix, column 1 is assumed to
        #  correspond to x[[1]], column 2 to x[[2]], etc.
        #
        #  When in list mode x is assumed to have length JK, the total number
        #  groups being tested, but a subset of the data can be analyzed
        #  using grp
        #
if(is.data.frame(x))data=as.matrix(x)
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x)) y[[j]] <- x[, j]
                x <- y
}

#        conM = con3way(J,K,L)
 p <- J*K*L
if(p>length(x))stop("JKL is less than the Number of groups")
JK=J*K
        v <- matrix(0, p, p)
        data <- list()
xx=list()
        for(j in 1:length(x)) {
                data[[j]] <- x[[grp[j]]]
xx[[j]]=x[[grp[j]]] # save input data
                # Now have the groups in proper order.
                data[[j]] = data[[j]] - mean(data[[j]], tr = tr)
        }
#ilow=0-L
#iup=0
#for(j in 1:JK){
#ilow <- ilow + L
# iup = iup + L
#sel <- c(ilow:iup)
#xx[sel]=listm(elimna(matl(xx[sel])))
# v[sel, sel] <- covmtrim(xx[sel], tr)
#                }
test.stat=bbwtrim(J,K,L,xx,tr=tr)
        x <- data  # Centered data
#        jp <- 1 - K
#        kv <- 0
        if(SEED)
                set.seed(2)
        # set seed of random number generator so that
        #             results can be duplicated.
        testA = NA
        testB = NA
        testC=NA
        testAB = NA
        testAC = NA
        testBC = NA
        testABC = NA
        bsam = list()
        bdat = list()
aboot=NA
bboot=NA
cboot=NA
abboot=NA
acboot=NA
bcboot=NA
abcboot=NA
nvec=NA
        for(j in 1:JK){
                nvec[j] = length(x[[j]])
        for(ib in 1:nboot) {
                ilow <- 1 - L
                iup = 0
 for(j in 1:JK) {
 ilow <- ilow + L
 iup = iup + L
nv=length(x[[ilow]])
 bdat[[j]] = sample(nv, size = nv, replace =TRUE)
for(k in ilow:iup){
 bsam[[k]] = x[[k]][bdat[[j]]]
}
}
temp=bbwtrim(J,K,L,bsam,tr=tr)
aboot[ib]=temp$Qa
bboot[ib]=temp$Qb
cboot[ib]=temp$Qc
acboot[ib]=temp$Qac
bcboot[ib]=temp$Qbc
abboot[ib]=temp$Qab
abcboot[ib]=temp$Qabc
}}
pbA=NA
pbB=NA
pbC=NA
pbAB=NA
pbAC=NA
pbBC=NA
pbABC=NA
pbA=mean(test.stat$Qa[1,1]<aboot)
pbB=mean(test.stat$Qb[1,1]<bboot)
pbC=mean(test.stat$Qc[1,1]<cboot)
pbAB=mean(test.stat$Qab[1,1]<abboot)
pbAC=mean(test.stat$Qac[1,1]<acboot)
pbBC=mean(test.stat$Qbc[1,1]<bcboot)
pbABC=mean(test.stat$Qabc[1,1]<abcboot)
list(p.value.A=pbA,p.value.B=pbB,p.value.C=pbC,p.value.AB=pbAB,
p.value.AC=pbAC,p.value.BC=pbBC,p.value.ABC=pbABC)
}

# bwwtrim
bwwtrim<-function(J,K,L,data,tr=.2,grp=c(1:p),alpha=.05,p=J*K*L){
#  Perform a between by within by within (three-way) anova
#  on trimmed means where
#
#  J independent groups, KL dependent groups
#
#  The variable data is assumed to contain the raw
#  data stored in list mode. data[[1]] contains the data
#  for the first level of all three factors: level 1,1,1.
#  data][2]] is assumed to contain the data for level 1 of the
#  first two factors and level 2 of the third factor: level 1,1,2
#  data[[L]] is the data for level 1,1,L
#  data[[L+1]] is the data for level 1,2,1. data[[2L]] is level 1,2,L.
#  data[[KL+1]] is level 2,1,1, etc.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that data has length JKL, the total number of
#  groups being tested.
#
if(is.data.frame(data))data=as.matrix(data)
if(is.list(data))data=bwwna(J,K,L,data) # remove missing values
if(is.matrix(data))data=bwwmatna(J,K,L,data) #remove missing values
#                                     and convert to list mode
if(!is.list(data))stop("The data are not stored in list mode or a matrix")
if(p!=length(data)){
print("The total number of groups, based on the specified levels, is")
print(p)
print("The number of groups in data is")
print(length(data))
print("Warning: These two values are not equal")
}
tmeans<-0
h<-0
v<-0
for (i in 1:p){
tmeans[i]<-mean(data[[grp[i]]],tr)
h[i]<-length(data[[grp[i]]])-2*floor(tr*length(data[[grp[i]]]))
#    h is the effective sample size
}
v=bwwcovm(J,K,L,data,tr=tr)
ij<-matrix(c(rep(1,J)),1,J)
ik<-matrix(c(rep(1,K)),1,K)
il<-matrix(c(rep(1,L)),1,L)
jm1<-J-1
cj<-diag(1,jm1,J)
cj<-diag(1,jm1,J)
for (i in 1:jm1)cj[i,i+1]<-0-1
km1<-K-1
ck<-diag(1,km1,K)
for (i in 1:km1)ck[i,i+1]<-0-1
lm1<-L-1
cl<-diag(1,lm1,L)
for (i in 1:lm1)cl[i,i+1]<-0-1
#  Do test for factor A
cmat<-kron(cj,kron(ik,il))  # Contrast matrix for factor A
Qa=bwwtrim.sub(cmat, tmeans, v, h,p)
Qa.siglevel <- 1 - pf(Qa, J - 1, 999)
# Do test for factor B
cmat<-kron(ij,kron(ck,il))  # Contrast matrix for factor B
Qb=bwwtrim.sub(cmat, tmeans, v, h,p)
 Qb.siglevel <- 1 - pf(Qb, K - 1, 999)
# Do test for factor C
cmat<-kron(ij,kron(ik,cl))  # Contrast matrix for factor C
Qc<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qc.siglevel <- 1 - pf(Qc, L - 1, 999)
# Do test for factor A by B interaction
cmat<-kron(cj,kron(ck,il))  # Contrast matrix for factor A by B
Qab<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qab.siglevel <- 1 - pf(Qab, (J - 1) * (K - 1), 999)
# Do test for factor A by C interaction
cmat<-kron(cj,kron(ik,cl))  # Contrast matrix for factor A by C
Qac<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qac.siglevel <- 1 - pf(Qac, (J - 1) * (L - 1), 999)
# Do test for factor B by C interaction
cmat<-kron(ij,kron(ck,cl))  # Contrast matrix for factor B by C
Qbc<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qbc.siglevel <- 1 - pf(Qbc, (K - 1) * (L - 1), 999)
# Do test for factor A by B by C interaction
cmat<-kron(cj,kron(ck,cl))  # Contrast matrix for factor A by B by C
Qabc<-bwwtrim.sub(cmat, tmeans, v, h,p)
Qabc.siglevel <-1-pf(Qabc,(J-1)*(K-1)*(L-1), 999)
list(Qa=Qa,Qa.p.value=Qa.siglevel,Qb=Qb,Qb.p.value=Qb.siglevel,
Qc=Qc,Qc.p.value=Qc.siglevel,Qab=Qab,Qab.p.value=Qab.siglevel,
Qac=Qac,Qac.p.value=Qac.siglevel,Qbc=Qbc,Qbc.p.value=Qbc.siglevel,
Qabc=Qabc,Qabc.p.value=Qabc.siglevel)
}



# bwwtrimbt
bwwtrimbt<-function(J, K, L, x, tr = 0.2, JKL = J * K*L, con = 0,
 alpha = 0.05, grp =c(1:JKL), nboot = 599, SEED = TRUE, ...)
{
        #
        # A bootstrap-t for omnibus tests associated with
        # all main effects and interactions.
        # in a between-by-within-within design.
        #
        #  The R variable x is assumed to contain the raw
        #  data stored in list mode or in a matrix.
        #  If in list mode, x[[1]] contains the data
        #  for the first level of both factors: level 1,1.
        #  x[[2]] is assumed to contain the data for level 1 of the
        #  first factor and level 2 of the second: level 1,2
        #  x[[K]] is the data for level 1,K
        #  x[[K+1]] is the data for level 2,1, x[[2K]] is level 2,K, etc.
        #
        #  If the data are in a matrix, column 1 is assumed to
        #  correspond to x[[1]], column 2 to x[[2]], etc.
        #
        #  When in list mode x is assumed to have length JK, the total number
        #  groups being tested, but a subset of the data can be analyzed
        #  using grp
        #
if(is.data.frame(x))data=as.matrix(x)
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x)) y[[j]] <- x[, j]
                x <- y
}

        conM = con3way(J,K,L)
 p <- J*K*L
if(p>length(x))stop("JKL is less than the Number of groups")
JK=J*K
KL=K*L
        v <- matrix(0, p, p)
        data <- list()
xx=list()
        for(j in 1:length(x)) {
                data[[j]] <- x[[grp[j]]]
xx[[j]]=x[[grp[j]]] # save input data
                # Now have the groups in proper order.
                data[[j]] = data[[j]] - mean(data[[j]], tr = tr)
        }
        x <- data  # centered data xx has original data
test=bwwtrim(J,K,L,xx,tr=tr)
        if(SEED)
                set.seed(2)
        # set seed of random number generator so that
        #             results can be duplicated.
        bsam = list()
        bdat = list()
aboot=NA
bboot=NA
cboot=NA
abboot=NA
acboot=NA
bcboot=NA
abcboot=NA
        for(ib in 1:nboot) {
                ilow <- 1 - KL
                iup = 0
 for(j in 1:J) {
 ilow <- ilow + KL
 iup = iup + KL
nv=length(x[[ilow]])
 bdat[[j]] = sample(nv, size = nv, replace =TRUE)
for(k in ilow:iup){
 bsam[[k]] = x[[k]][bdat[[j]]]
}
}
temp=bwwtrim(J,K,L,bsam,tr=tr)
aboot[ib]=temp$Qa
bboot[ib]=temp$Qb
cboot[ib]=temp$Qc
acboot[ib]=temp$Qac
bcboot[ib]=temp$Qbc
abboot[ib]=temp$Qab
abcboot[ib]=temp$Qabc
        }
pbA=NA
pbB=NA
pbC=NA
pbAB=NA
pbAC=NA
pbBC=NA
pbABC=NA
pbA=mean(test$Qa[1,1]<aboot)
pbB=mean(test$Qb[1,1]<bboot)
pbC=mean(test$Qc[1,1]<cboot)
pbAB=mean(test$Qab[1,1]<abboot)
pbAC=mean(test$Qac[1,1]<acboot)
pbBC=mean(test$Qbc[1,1]<bcboot)
pbABC=mean(test$Qabc[1,1]<abcboot)
list(p.value.A=pbA,p.value.B=pbB,p.value.C=pbC,p.value.AB=pbAB,
p.value.AC=pbAC,p.value.BC=pbBC,p.value.ABC=pbABC)

}

# dtrimpb
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

# dtrimQS
dtrimQS<-function(x,y=NULL,tr=.2,pr=TRUE){
#
#  Trimmed mean based on difference scores
#  Also returns quantile shift measure of location
#
#
if(pr){
print('Note: Under normality and homoscedasticity, Cohen d= 0, .2, .5, .8')
print('correspond approximately to Q.effect = 0.5, 0.55, 0.65 and 0.70, respectively')
}
if(!is.null(y))L=x-y
else L=x
L=elimna(L)
output=trimci(L,tr=tr,pr=FALSE)
ef=depQS(L,locfun=mean,tr=tr)
list(ci=output$ci,estimate=output$estimate,test=output$test.stat,
se=output$se,p.value=output$p.value,n=output$n,Q.effect=ef$Q.effect)
}

# dlintrim
dlintrim<-function(x,con,SEED=TRUE,xlab='DV',ylab='',sym.test=FALSE,
plotit=TRUE,tr=.2,PB=FALSE){
#
# For dependent variables,
# determine distribution of Y_i=sum_j c_jX_j
# and then plot the distribution and test the hypothesis that
# Y has a trimmed mean equal to zero.
#
#  If x is a matrix, columns correspond to groups.
#
con=as.vector(con)
if(sum(con)!=0)stop('Contrast coefficients must sum to zero')
if(SEED)set.seed(2)
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
if(!PB)mt=trimci(L,tr=tr,pr=FALSE)
if(PB)mt=trimpb(L,tr=tr,pr=FALSE)
sym=NULL
Q=depQS(L,locfun=mean,tr=tr)
if(sym.test)sym=Dqdif(L)
list(trim.mean=mt$estimaten,n=mt$n,ci=mt$ci,
p.value=mt$p.value,Q.effect=Q$Q.effect,sym.test=sym)
}


# ghtrim
ghtrim<-function(tr=.2,g=0,h=0){
#
#  Compute trimmed mean of a g-and-h distribution.
#
#
if(g==0)val=0
if(g>0){
low=qnorm(tr)
up=-1*low
val=integrate(ftrim,low,up,tr=tr,g=g,h=h)$value
val=val/(1-2*tr)
}
val
}


# bbtrimQS
bbtrimQS<-function(J,K,x, con = 0, tr = 0.2,
alpha = 0.05, pr = TRUE, crit = NA, SEED = TRUE, INT = FALSE, locfun = tmean){
if(pr){
print('Note: confidence intervals are adjusted to control FWE')
print('But p-values are not adjusted to control FWE')
print('Adjusted p-values can be computed with the R function p.adjusted')
print('Under normality and homoscedasticity, Cohen d= 0, .2, .5, .8')
print('correspond approximately to Q.effect = 0.5, 0.55, 0.65 and 0.70, respectively')
}
CON=con2way(J,K)
A=linconQS(x,con=CON$conA,tr=tr,alpha=alpha,pr=FALSE,crit=crit,SEED=SEED,INT=FALSE,locfun=locfun)
B=linconQS(x,con=CON$conB,tr=tr,alpha=alpha,pr=FALSE,crit=crit,SEED=SEED,INT=FALSE,locfun=locfun)
AB=linconQS(x,con=CON$conAB,tr=tr,alpha=alpha,pr=FALSE,crit=crit,SEED=SEED,INT=TRUE,locfun=locfun)
list(A=A,B=B,AB=AB)
}


# biloc
biloc<-function(x){
#
# compute biweight measure of location
# This function is used by relplot
#
m<-median(x)
u<-abs((x-m)/(9*mad(x)*qnorm(.75)))
top<-sum((x[u<=1]-m)*(1-u[u<=1]^2)^2)
bot<-sum((1-u[u<=1]^2)**2)
bi<-m+top/bot^2
bi
}


# dep.loc.summary
dep.loc.summary<-function(x,y){
#
# Estimate the measures of location based on difference scores:
#  MEAN:
#  MEAN20: 20% mean
#  MED:  median
#  OS: One-step M-estimator
#  Mdif: median of typical difference
#
d=NULL
chk=ncol(x)
if(is.vector(x))d=x-y
if(!is.null(chk)){
if(dim(x)[2]==2)d=x[,1]-x[,2]
}
if(is.null(d))stop('x and y should be vectors, or x should have two columns')
d=elimna(d)
output=matrix(NA,ncol=1,nrow=4)
output[1,1]=tmean(d,tr=0)
output[2,1]=tmean(d,tr=0.2)
output[3,1]=median(d)
output[4,1]=onestep(d)
dimnames(output)=list(c('MEAN','MEAN.20%','MEDIAN','M-EST'),c('Est'))
output
}




# loc2dif
loc2dif<-function(x,y=NULL,est=median,na.rm=TRUE,plotit=FALSE,xlab="",ylab="",...){
#
# Compute a measure of location associated with the
# distribution of x-y, the typical difference between two randomly sampled values.
#  The measure of location is indicated by the argument
# est.
# x and y are paired data or independent variables having the same length.
#  If x and y have different lengths, use the function wmwloc
#
# Advantage of this estimator: relatively high efficiency even under normality versus
# using sample means.
#
if(is.null(y)){
if(ncol(x)!=2)stop("x should be an n-by-2 matrix")
y=x[,2]
x=x[,1]
if(na.rm)m=elimna(cbind(x,y))
x=m[,1]
y=m[,2]
}
x=elimna(x)
y=elimna(y)
temp=as.vector(outer(x,y,FUN="-"))
val<-est(temp,na.rm=TRUE,...)
if(plotit)akerd(temp,xlab=xlab,ylab=ylab)
val
}


# loc2difpb
loc2difpb<-function(x,y,est=median,alpha=.05,nboot=2000,SEED=TRUE){
#
# A percentile bootstrap
# confidence interval for the median of D=X-Y
#
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
v=NA
for(i in 1:nboot){
X=sample(x,n1,replace=TRUE)
Y=sample(y,n2,replace=TRUE)
v[i]=loc2dif(X,Y)
}
pv=mean(v<0)+.5*mean(v==0)
pv=2*min(c(pv,1-pv))
vs=sort(v)
crit<-alpha/2
icl<-round(crit*nboot)+1
icu<-nboot-icl
ci=vs[icl]
ci[2]=vs[icu]
list(ci=ci,p.value=pv)
}


# loc2gmulpb
loc2gmulpb<-function(x,y,est=tmean,nboot=2000,alpha = 0.05,SEED=TRUE,...){
#
# For two independent  p-variate distributions,  apply yuen to each column of data
# FWE controlled with Hochberg's method
#
# x and y are  matrices having p columns. (Can have list mode as well.)
#
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
if(!is.matrix(y))y<-matl(y)
if(!is.matrix(y))stop('Data must be stored in a matrix or in list mode.')

J<-ncol(x)
if(J!=ncol(y))stop('x and y should have the same number of columns')

xbar<-vector('numeric',J)
ncon<-J
dvec<-alpha/c(1:ncon)
psihat<-matrix(0,J,4)
dimnames(psihat)<-list(NULL,c('Variable','difference','ci.lower','ci.upper'))
test<-matrix(0,J,4)
dimnames(test)<-list(NULL,c('Variable','p.value','p.crit','se'))
temp1<-NA
nval=NULL
for (d in 1:J){
psihat[d,1]<-d
#temp=yuen(x[,d],y[,d],tr=tr)
temp=pb2gen(x[,d],y[,d],est=est,SEED=SEED,...)
test[d,1]<-d
#test[d,2]<-temp$teststat
test[d,2]=temp$p.value
test[d,4]<-sqrt(temp$sq.se)
psihat[d,2]<-temp$est.dif
psihat[d,3]<-temp$ci[1]
psihat[d,4]<-temp$ci[2]
}
temp1=test[,2]
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
test[temp2,3]<-zvec
num.sig=sum(test[,2]<=test[,3])
list(n=c(nrow(x),nrow(y)),test=test,psihat=psihat,num.sig=num.sig)
}



# loc2plot
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

# mdifloc
mdifloc<-function(x,y,est=tukmed,...){
#
# Compute multivariate measure of location associated
# with the distribution of x-y
#
# By default, use Tukey's median.
#
x<-as.matrix(x)
y<-as.matrix(y)
FLAG<-F
if(ncol(x)!=ncol(y))stop("x and y should have the same number of columns")
if(ncol(x)==1 && ncol(y)==1)FLAG<-T
if(FLAG)val<-loc2dif(x,y,est=est,...)
if(!FLAG){
J<-(ncol(x)^2-ncol(x))/2
mat<-matrix(NA,ncol=ncol(x),nrow=nrow(x)*nrow(y))
for(j in 1:ncol(x))mat[,j]<-as.vector(outer(x[,j], y[,j], FUN = "-"))
val<-est(mat,...)
}
val
}

# M2m.loc
M2m.loc<-function(m,grpc,col.dat,locfun=tmean,...){
#
# m is a matrix or data frame.
# Compute a measure of location for each of several categories, with
# categories indicated by the values in the column of m given by the
# argument grpc.
# The argument grpc can have up to 4 values, which correspond to factors.
#
#  col.dat indicates the column of m containing the outcome measure
# of interest.
# locfun indicates the measure of location, which defaults to the 20%
# trimmed mean.
#
#  Example,
# M2m.loc(x,c(1,4),5,locfun=mean)
# indicates that there are 2 factors, with levels of the factors indicated
# by the values in columns 1 and 4 of the matrix x. For each combination
# of levels,
# locfun=mean
# indicates that the sample mean will be computed.
#
flagit=F
if(is.null(dim(m)))stop("Data must be stored in a matrix or data frame")
if(is.na(grpc[1]))stop("The argument grpc is not specified")
if(is.na(col.dat[1]))stop("The argument col.dat is not specified")
if(length(grpc)>4)stop("grpc must have length <= 4")
m=as.data.frame(m)
if(length(grpc)==1){
p1=ncol(m)+1
dum=rep(1,nrow(m))
flagit=T
m=cbind(m,dum)
grpc=c(NULL,gprc,p1)
cat1<-sort(unique(m[,grpc[1]]))
M=NULL
for (ig1 in 1:length(cat1)){
flag1=(m[,grpc[1]]==cat1[ig1])
flag=(flag1==1)
msub=as.data.frame(m[flag,])
loc=locfun(m[flag,col.dat],...)
M=rbind(M,as.data.frame(cbind(msub[1,grpc],loc)))
}
M=M[,c(1,3)]
}
if(length(grpc)==2){
cat1<-sort(unique(m[,grpc[1]]))
cat2<-sort(unique(m[,grpc[2]]))
M=NULL
for (ig1 in 1:length(cat1)){
for (ig2 in 1:length(cat2)){
flag1=(m[,grpc[1]]==cat1[ig1])
flag2=(m[,grpc[2]]==cat2[ig2])
flag=(flag1*flag2==1)
msub=m[flag,]
loc=locfun(m[flag,col.dat],...)
M=rbind(M,as.data.frame(cbind(msub[1,grpc],loc)))
}}}
if(length(grpc)==3){
cat1<-sort(unique(m[,grpc[1]]))
cat2<-sort(unique(m[,grpc[2]]))
cat3<-sort(unique(m[,grpc[3]]))
M=NULL
for (ig1 in 1:length(cat1)){
for (ig2 in 1:length(cat2)){
for (ig3 in 1:length(cat3)){
flag1=(m[,grpc[1]]==cat1[ig1])
flag2=(m[,grpc[2]]==cat2[ig2])
flag3=(m[,grpc[3]]==cat3[ig3])
flag=(flag1*flag2*flag3==1)
msub=m[flag,]
loc=locfun(m[flag,col.dat],...)
M=rbind(M,as.data.frame(cbind(msub[1,grpc],loc)))
}}}}
if(length(grpc)==4){
cat1<-sort(unique(m[,grpc[1]]))
cat2<-sort(unique(m[,grpc[2]]))
cat3<-sort(unique(m[,grpc[3]]))
cat4<-sort(unique(m[,grpc[4]]))
M=NULL
for (ig1 in 1:length(cat1)){
for (ig2 in 1:length(cat2)){
for (ig3 in 1:length(cat3)){
for (ig4 in 1:length(cat4)){
flag1=(m[,grpc[1]]==cat1[ig1])
flag2=(m[,grpc[2]]==cat2[ig2])
flag3=(m[,grpc[3]]==cat3[ig3])
flag4=(m[,grpc[4]]==cat4[ig4])
flag=(flag1*flag2*flag3*flag4==1)
msub=m[flag,]
loc=locfun(m[flag,col.dat],...)
M=rbind(M,as.data.frame(cbind(msub[1,grpc],loc)))
}}}}}
if(flagit)M=M[,c(1,3)]
M
}

# mul.loc2g
mul.loc2g<-function(m1,m2,nullv=rep(0,ncol(m1)),locfun=smean,alpha=.05,SEED=TRUE,
nboot=500,...){
#
# m is an n by p matrix
#
# For two independent groups,
# test hypothesis that multivariate measures of locations  are all equal.
# This is done by computing a confidence interval for each of the
# p variables under study using a percentile bootstrap method.
#
#
if(ncol(m1) != ncol(m2)){
stop('Number of variables in group 1 is not equal to the number in group 2.')
}
p=ncol(m1)
nb1=nboot+1
if(SEED)set.seed(2)
m1<-elimna(m1)
m2<-elimna(m2)
n1<-nrow(m1)
n2<-nrow(m2)
ci=matrix(NA,nrow=2,ncol=p)
bvec=matrix(NA,nrow=nboot,ncol=p)
val<-matrix(0,ncol=ncol(m1),nrow=nboot)
est1=locfun(m1,...)
est2=locfun(m2,...)
if(is.list(est1)){
est1=est1$center
est2=est2$center
}
for(j in 1: nboot){
id1<-sample(n1,size=n1,replace=TRUE)
id2<-sample(n2,size=n2,replace=TRUE)
v1<-locfun(m1[id1,],...)
v2<-locfun(m2[id2,],...)
if(is.list(v1)){
val[j,]=(v1$center<v2$center)
bvec[j,]=v1$center-v2$center
}
else{
 val[j,]=(v1<v2)
bvec[j,]=v1-v2
}}
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
pv=apply(val,2,mean)
mat=rbind(pv,1-pv)
pv=2*apply(mat,2,min)

for(k in 1:p){
temp=sort(bvec[,k])
ci[1,k]=temp[ilow]
ci[2,k]=temp[ihi]
}
L=NA
for(i in 1:p)L[i]=paste('v',i)
dimnames(ci)=list(c('Low','Up'),L)
list(est1=est1,est2=est2,ci=ci,p.value=pv)
}


# Dmul.loc2g
Dmul.loc2g<-function(m1,m2=NULL,nullv=rep(0,ncol(m1)),locfun=smean,alpha=.05,SEED=TRUE,
nboot=500,...){
#
#  Inferences about a multivariate measure of location. A percentile bootstrap method
#  is used for each of the p variables.
#
#  By default, a skipped estimator is used where outliers are identified and removed via projection method.
#  For a large number of variables, say > 8, might include the argument
#  op=3 and
#  outfun=outproadj
#  But this will increase execution time considerably.
#
# m1  is an n by p matrix
#
# For single distribution, m2=NULL,
# test the hypothesis  that  the measure of location estimated by
# locfun is equal to the value specified by
# nullv.
#
#  For two dependent groups, meaning that
#  m2 is not mull
# test hypothesis that the difference scores have measures of locations that  are equal to nullv
#
# This is done by computing a confidence interval for each of the
# p variables under study using a percentile bootstrap.
#
#
if(is.null(m2))D=m1
else{
if(ncol(m1) != ncol(m2))stop('Number of variables in group 1 does not equal the number in group 2.')
D=m1-m2
}
names(D)=NULL
p=ncol(D)
nb1=nboot+1
if(SEED)set.seed(2)
D<-elimna(D)
n<-nrow(D)
val<-matrix(0,ncol=p,nrow=nboot)
bvec=matrix(NA,nboot,p)
ci=matrix(NA,2,p)
est=locfun(D,...)
if(is.list(est))est=est$center
for(j in 1: nboot){
id<-sample(n,size=n,replace=TRUE)
v1<-locfun(D[id,],...)
if(is.list(v1)){
val[j,]=(v1$center<nullv)
bvec[j,]=v1$center
}
else{
 val[j,]=(v1<nullv)
 bvec[j,]=v1
}}
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
pv=apply(val,2,mean)
mat=rbind(pv,1-pv)
pv=2*apply(mat,2,min)
for(k in 1:p){
temp=sort(bvec[,k])
ci[1,k]=temp[ilow]
ci[2,k]=temp[ihi]
}
L=NA
for(i in 1:p)L[i]=paste('v',i)
dimnames(ci)=list(c('Low','Up'),L)
pad=p.adjust(pv,method='holm')
list(est=est,ci=ci,p.values=pv,adjusted.p.value=pad)
}



# loc.dif.summary
loc.dif.summary<-function(x,y){
#
# Estimate the difference between a collection of measures of location::
#  MEAN:
#  MEAN20: 20% mean
#  MED:  median
#  OS: One-step M-estimator
#  Mdif: median of typical difference
#
x=elimna(x)
y=elimna(y)
output=matrix(NA,ncol=1,nrow=5)
output[1,1]=tmean(x,tr=0)-tmean(y,tr=0)
output[2,1]=tmean(x,tr=0.2)-tmean(y,tr=0.2)
output[3,1]=median(x)-median(y)
output[4,1]=onestep(x)-onestep(y)
output[5,1]=wmwloc(x,y)
dimnames(output)=list(c('MEAN','MEAN.20%','MEDIAN','M-EST','Mdif'),c('Est'))
output
}



# L1median
L1median <- function(X, tol = 1e-08, maxit = 200, m.init = apply(X, 2, median),
                     trace = FALSE)
{
  ## L1MEDIAN calculates the multivariate L1 median
  ## I/O: mX=L1median(X,tol);
  ##
  ## X  : the data matrix
  ## tol: the convergence criterium:
  ##      the iterative process stops when ||m_k - m_{k+1}|| < tol.
  ## maxit: maximum number of iterations
  ## init.m: starting value for m; typically coordinatewise median
  ##
  ## Ref: Hossjer and Croux (1995)
  ##  "Generalizing Univariate Signed Rank Statistics for Testing
  ##   and Estimating a Multivariate Location Parameter";
  ##   Non-parametric Statistics, 4, 293-308.
  ##
  ## Implemented by Kristel Joossens
  ## Many thanks to Martin Maechler for improving the program!

  ## slightly faster version of 'sweep(x, 2, m)':
  centr <- function(X,m) X - rep(m, each = n)
  ## computes objective function in m based on X and a:
  mrobj <- function(X,m) sum(sqrt(rowSums(centr(X,m)^2)))

  d <- dim(X); n <- d[1]; p <- d[2]
  m <- m.init
  if(!is.numeric(m) || length(m) != p)
      stop("'m.init' must be numeric of length p =", p)
  k <- 1
  if(trace) nstps <- 0
  while (k <= maxit) {
    mold <- m
    obj.old <- if(k == 1) mrobj(X,mold) else obj
    X. <- centr(X, m)
    Xnorms <- sqrt(rowSums(X. ^ 2))
    inorms <- order(Xnorms)
    dx <- Xnorms[inorms] # smallest first, i.e., 0's if there are
    X  <- X [inorms,]
    X. <- X.[inorms,]
    ## using 1/x weighting {MM: should this be generalized?}
    w <- ## (0 norm -> 0 weight) :
        if (all(dn0 <- dx != 0))  1/dx
        else c(rep.int(0, length(dx)- sum(dn0)), 1/dx[dn0])
    delta <- colSums(X. * rep(w,p)) / sum(w)
    nd <- sqrt(sum(delta^2))

    maxhalf <- if (nd < tol) 0 else ceiling(log2(nd/tol))
    m <- mold + delta    # computation of a new estimate
    ## If step 'delta' is too far, we try halving the stepsize
    nstep <- 0
    while ((obj <- mrobj(X, m)) >= obj.old && nstep <= maxhalf) {
      nstep <- nstep+1
      m <- mold + delta/(2^nstep)
    }
    if(trace) {
        if(trace >= 2)
            cat(sprintf("k=%3d obj=%19.12g m=(",k,obj),
                paste(formatC(m),collapse=","),
                ")", if(nstep) sprintf(" nstep=%2d halvings",nstep) else "",
                "\n", sep="")
        nstps[k] <- nstep
    }
    if (nstep > maxhalf) { ## step halving failed; keep old
        m <- mold
        ## warning("step halving failed in ", maxhalf, " steps")
        break
      }
    k <- k+1
  }
  if (k > maxit) warning("iterations did not converge in ", maxit, " steps")
  if(trace == 1)
      cat("needed", k, "iterations with a total of",
          sum(nstps), "stepsize halvings\n")
  return(m)
}

# bmean
bmean<-function(x,na.rm=TRUE){
#
#  Compute a skipped estimator of location.
#  where outliers are flagged based on a boxplot rule
#
if(na.rm)x<-x[!is.na(x)] #Remove missing values
flag<-outbox(x)$keep
es<-mean(x[flag])
es
}


# dmean
dmean<-function(m,tr=.2,dop=1,cop=2){
#
# Compute multivariate measure of location
# using Donoho-Gasko method.
#
# dop=1, use fdepth to compute depths
# dop=2, use fdepthv2  to compute depths
#
# cop=1, Tukey median; can't be used here.
# cop=2, use MCD in fdepth
# cop=3, use marginal medians in fdepth
# cop=4, use MVE in fdepth
#
if(is.list(m))m<-matl(m)
if(!is.matrix(m))stop("Data must be stored in a matrix or in list mode.")
if(ncol(m)==1){
if(tr==.5)val<-median(m)
if(tr>.5)stop("Amount of trimming must be at most .5")
if(tr<.5)val<-mean(m,tr)
}
if(ncol(m)>1){
temp<-NA
if(ncol(m)!=2){
# Use approximate depth
if(dop==1)temp<-fdepth(m,plotit=FALSE,cop=cop)
if(dop==2)temp<-fdepthv2(m)
}
#  Use exact depth if ncol=2
if(ncol(m)==2){
for(i in 1:nrow(m))
temp[i]<-depth(m[i,1],m[i,2],m)
}
mdep<-max(temp)
flag<-(temp==mdep)
if(tr==.5){
if(sum(flag)==1)val<-m[flag,]
if(sum(flag)>1)val<-apply(m[flag,],2,mean)
}
if(tr<.5){
flag2<-(temp>=tr)
if(sum(flag2)==0 && sum(flag)>1)val<-apply(as.matrix(m[flag,]),2,mean)
if(sum(flag2)==0 && sum(flag)==1)val=m[flag,]
if(sum(flag2)==1)val<-m[flag2,]
if(sum(flag2)>1)val<-apply(m[flag2,],2,mean)
}}
val
}


# mmean
mmean<-function(x,est=tmean,...){
center<-NA
if(is.list(x))center=lapply(x,est,...)
if(is.matrix(x))center<-apply(x,2,est,...)
center
}



# ghmean
ghmean<-function(g,h){
#
#Compute the mean and variance of a g-and-h distribution
#
val=0
if(h==0){
if(g>0){
val=(exp(g^2/2)-1)/g
val2=(1-2*exp(g^2/2)+exp(2*g^2))/g^2
val2=val2-val^2
}}
if(g>0 & h!=0){
if(h<1)
val=(exp(g^2/(2*(1-h)))-1)/(g*sqrt(1-h))
val2=NA
if(h>0){
if(h<.5)
val2=(exp(2*g^2/(1-2*h))-2*exp(g^2/(2*(1-2*h)))+1)/(g^2*sqrt(1-2*h))-
(exp(g^2/(2*(1-h)))-1)^2/(g^2*(1-h))
}}
if(g==0){
val=0
val2=1/(1-2*h)^1.5   #Headrick et al. (2008)
}
list(mean=val,variance=val2)
}


# mgvmean
mgvmean<-function(m,op=0,outfun=outbox,se=TRUE){
#
# m is an n by p matrix
#
# Compute a multivariate skipped measure of location
# using the MGV method
#
# Eliminate outliers using MGV method
#
# op=0 pairwise distances of points
# op=1 MVE distances
# op=2 MCD distances
#
# outfun indicates outlier rule to be applied to
# the MGV distances.
# By default, use boxplot rule
#
# Eliminate any outliers and compute means
#  using remaining data.
#
m<-elimna(m)
temp<-outmgv(m,op=op,plotit=FALSE)$keep
val<-apply(m[temp,],2,mean)
val
}


# harmonic.mean
harmonic.mean<-function(x)1/mean(1/x)

rplotv2<-
function(x,y,est=tmean,scat=TRUE,fr=NA,plotit=TRUE,pyhat=FALSE,efr=.5,pch1='*',pch2='.',
theta=50,phi=25,scale=TRUE,expand=.5,SEED=TRUE,varfun=pbvar,outfun=outpro,
nmin=0,xout=FALSE,out=FALSE,eout=FALSE,xlab='X',ylab='Y',zscale=FALSE,
zlab=' ',pr=TRUE,duplicate='error',ticktype='simple',LP=TRUE,OLD=FALSE,pch='.',prm=TRUE,...){
#
#  Like rplot but can handle one or two binary independent variables,
#  at least one non-binary independent variable is required.
#
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
xx<-cbind(x,y)
xx<-elimna(xx)
n=nrow(xx)
if(eout){
flag=outfun(xx,plotit=FALSE,...)$keep
xx=xx[flag,]
}
if(xout){
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
xvals=list()
id=chk4binary(x)
Lid=length(id)
if(Lid==ncol(x))stop('All independent variables are binary, a smoother is inappropriate')
if(Lid>2)stop('Can have a most two binary independent variables')
val=list()
if(Lid==1){
xval=list()
yhat=list()
if(is.na(fr))fr=.8
irow0=which(x[,id]==0)
val[[1]]=rplot(x[irow0,-id],y[irow0],pyhat=TRUE,plotit=FALSE,est=est,xlab=xlab,ylab=ylab,pr=FALSE)$yhat
irow1=which(x[,id]==1)
#print(x[irow1,-id])
val[[2]]=rplot(x[irow1,-id],y[irow1],pyhat=TRUE,plotit=FALSE,est=est,xlab=xlab,ylab=ylab,pr=FALSE)$yhat
rplot2g(x[irow0,-id],y[irow0],x[irow1,-id],y[irow1],est=est,xlab=xlab,ylab=ylab,fr=fr,pch1=pch1,pch2=pch2) #,xout=xout,SEED=SEED)
xvals[[1]]=x[irow0,-id]
xvals[[2]]=x[irow1,-id]
}
if(Lid==2){
if(ncol(x)>3)stop(' With two binary IVs, current version limited to a third continuous IV')
xval=NULL
yhat=NULL
xval=list()
yhat=list()
xy=cbind(x[,id],x[,-id],y)
v=bin2binary.IV(xy)
val=rplot4g(v,est=est,xlab=xlab,ylab=ylab,pyhat=pyhat)
if(pyhat){
xvals=list()
val=list()
for(j in 1:4)xvals[[j]]=v[[j]][,3]
for(j in 1:4)val[[j]]=v[[j]][,4]
}
}

if(Lid==0){
if(pr && !OLD){
print('A new estimate of the strength of the association is used by default.')
print(' To get the old estimate, set OLD=TRUE')
}
if(ncol(x)==2 && !scale){
if(pr){print('scale=FALSE is specified.')
print('If there is dependence, might want to use scale=T')
}}
if(is.na(fr))fr<-1
val<-rung3d(x,y,est=est,fr=fr,plotit=plotit,pyhat=TRUE,SEED=SEED,nmin=nmin,LP=LP,
scale=scale,phi=phi,theta=theta,expand=expand,zscale=zscale,pr=FALSE,
duplicate='error',xlab=xlab,ylab=ylab,zlab=zlab,ticktype=ticktype,...)
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
}}}
list(n=n,n.keep=n.keep,xvals=xvals,yhat = val)
}


# lognormal.mean
lognormal.mean<-function(mu=0,sig.sq=1)exp(mu+sig.sq/2)


# lognormal.mom
lognormal.mom<-function(mu=0,sig.sq=1){
me=lognormal.mean(mu,sig.sq)
V=lognormal.var(mu,sig.sq)
sk=lognormal.skew(mu,sig.sq)
ku=lognormal.kurt(mu,sig.sq)
list(mean=me,var=V,skew=sk,kurtosis=ku)
}


# bca.mean
bca.mean<-function(x,alpha=.05,SEED=TRUE){
#
# BCA confidence interval for the mean
#
# Not recommended, commonly used, but provided only to illustrate
# situations where it fails.
#
library(bcaboot)
if(SEED)set.seed(2)
x=elimna(x)
x=as.matrix(x)
a=bcajack2(x,1000,mean,alpha=alpha/2)
ci=c(a$lims[1,1],a$lims[3,1])
be=a$stats[1,1]
sd=a$stats[1,2]
list(Est=be,SD=sd, CI=ci)
}


# bptdmean
bptdmean<-function(isub,x,tr){
#
#  Compute  trimmed means
#  when comparing dependent groups.
#  By default, 20% trimmed means are used.
#  isub is a vector of length n,
#  a bootstrap sample from the sequence of integers
#  1, 2, 3, ..., n
#
#  This function is used by bptd.
#
bptdmean<-mean(x[isub],tr)
bptdmean
}



# center.m
center.m<-function(x,est=tmean,...){
x=elimna(x)
x=as.matrix(x)
m1=lloc(x,est=est)
p=ncol(x)
for(j in 1:p)x[,j]=(x[,j]-m1[j])
x
}




# funloc
funloc<-function(x,tr=.2,pts=NULL,npts=25,plotit=TRUE,alpha=.05,nv=rep(0,ncol(x)),
xlab='T',ylab='Est.',FBP=TRUE,method='hochberg',COLOR=TRUE){
#
#  x1 and x2 are n-by-p matrices,
#  Designed for functional data.
#  For example, p measures taken over time where  p is typically large
#
# nv is the null value when testing some hypothesis
#
#  Goal: at speficied times, compute measures of location and confidence intervals.
#  pts: Can specify time points where comparisons are to be made
#  if pts=NULL, pick
#  npts points evenly space between min and max time points
#
#  FBP=T: creates a  functional boxplot
#  FBP=F: plot an estimate of the typical value plus 1-alpha confidence intervals.
#
p=ncol(x)
pm1=p-1
if(is.null(pts)){
np=round(p/npts)
if(np==0)np=1
pts=seq(2,pm1,np)
}
res=NA
dif=NA
op=matrix(NA,nrow=length(pts),ncol=6)
dimnames(op)=list(NULL,c('est.','s.e.','p.value',
'adjust.p.value','ci.low','ci.hi'))
for(i in 1:length(pts)){
z=trimci(x[,i],tr=tr,null.value=nv[i],pr=FALSE)
op[i,1]=z$estimate
op[i,2]=z$se
op[i,3]=z$p.value
op[i,5]=z$ci[1]
op[i,6]=z$ci[2]
}
op[,4]=p.adjust(op[,3],method=method)
if(plotit){
if(!FBP){
xlow=c(1:nrow(op))
xax=rep(c(1:nrow(op)),3)
plot(xax,as.vector(op[,c(3,5,6)]),type='n',xlab=xlab,ylab=ylab)
lines(xlow,op[,1])
lines(xlow,op[,5],lty=2)
lines(xlow,op[,6],lty=2)
}
if(FBP){
if(COLOR)FBplot(x)
if(!COLOR)func.out(x)
}}
op=cbind(pts,op)
op

}



# funlocpb
funlocpb<-function(x,est=tmean,nboot=2000,SEED=TRUE,
pts=NULL,npts=25,plotit=TRUE,alpha=.05,nv=rep(0,ncol(x)),
xlab='T',ylab='Est.',FBP=TRUE,method='hochberg',COLOR=TRUE,...){
#
#  x1 and x2 are n-by-p matrices,
#  Designed for functional data.
#  For example, p measures taken over time where  p is typically large
#
# nv is the null value when testing some hypothesis
#
#  Goal: at speficied times, compute measures of location and confidence intervals.
#  pts: Can specify time points where comparisons are to be made
#  if pts=NULL, pick
#  npts points evenly space between min and max time points
#
#  FBP=T: creates a  functional boxplot
#  FBP=F: plot an estimate of the typical value plus 1-alpha confidence intervals.
#
p=ncol(x)
pm1=p-1
if(is.null(pts)){
np=round(p/npts)
if(np==0)np=1
pts=seq(2,pm1,np)
}
res=NA
dif=NA
op=matrix(NA,nrow=length(pts),ncol=5)
dimnames(op)=list(NULL,c('est.','p.value',
'adjust.p.value','ci.low','ci.hi'))
x=elimna(x)
n=nrow(x)
for(i in 1:length(pts)){
z=onesampb(x[,i],est=est,nboot=nboot,alpha=alpha,SEED=SEED,nv=nv,...)
op[i,1]=z$estimate
op[i,2]=z$p.value
op[i,4]=z$ci[1]
op[i,5]=z$ci[2]
}
op[,3]=p.adjust(op[,2],method=method)
if(plotit){
if(!FBP){
xlow=c(1:nrow(op))
xax=rep(c(1:nrow(op)),3)
plot(xax,as.vector(op[,c(3,5,6)]),type='n',xlab=xlab,ylab=ylab)
lines(xlow,op[,1])
lines(xlow,op[,5],lty=2)
lines(xlow,op[,6],lty=2)
}
if(FBP){
if(COLOR)FBplot(x)
if(!COLOR)func.out(x)
}}
op=cbind(pts,op)
op

}



# locpre
locpre<-function(y,est=mean,error=sqfun,nboot=100,SEED=TRUE,pr=TRUE,mval=round(5*log(length(y)))){
#
#   Estimate the prediction error using a measure of location
#   given by the argument
#   est
#
#   The .632 method is used.
#   (See Efron and Tibshirani, 1993, pp. 252--254)
#
#   Prediction error is the expected value of the function error.
#   The argument error defaults to squared error.
#
#   est can be any R function that returns a measure of location
#
#   The default value for mval, the number of observations to resample
#   for each of the B bootstrap samples is based on results by
#   Shao (JASA, 1996, 655-665). (Resampling n vectors of observations
#   model selection may not lead to the correct model as n->infinity.
#
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data<-matrix(sample(length(y),size=mval*nboot,replace=TRUE),nrow=nboot)
bid<-apply(data,1,idb,length(y))
#  bid is an n by nboot matrix. If the jth bootstrap sample from
#  1, ..., mval contains the value i, bid[i,j]=0; otherwise bid[i,j]=1
#
yhat<-apply(data,1,locpres1,y,est=est)
# yhat is nboot vector
# containing the bootstrap estimates
#
yhat<-matrix(yhat,nrow=length(y),ncol=nboot) # convert to n x nboot matrix
bi<-apply(bid,1,sum) # B sub i in notation of Efron and Tibshirani, p. 253
temp<-(bid*(yhat-y))
diff<-apply(temp,1,error)
ep0<-sum(diff/bi)/length(y)
aperror<-error(y-est(y))/length(y) # apparent error
val<-.368*aperror+.632*ep0
val
}



# locpres1
locpres1<-function(isub,x,est){
#
#  Compute a measure of location x[isub]
#  isub is a vector of length mval,
#  a bootstrap sample from the sequence of integers
#  1, 2, 3, ..., n
#
#  mval is the sample size
#  of the bootstrap sample, where mval<n is used to get
#  consistency when choosing the correct model.

#  This function is used by other functions when computing
#  bootstrap estimates.
#
regboot<-est(x[isub])
regboot
}


# locreg
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


# locvar
locvar<-function(x,y,pyhat=FALSE,pts=x,plotit=TRUE){
#
# For each x, estimate VAR(y|x)
# with the method used by Bjerve and Doksum
# i.e., use Fan's kernel regression method.
#
yhat<-locreg(x,y,pyhat=TRUE,plotit=FALSE,pts=x)
val<-locreg(x,(y-yhat)^2,pyhat=pyhat,pts=pts,plotit=plotit)
val
}


# locvarsm
locvarsm<-function(x,y,pyhat=FALSE,pts=x,plotit=TRUE,nboot=40,RNA=TRUE,xlab="X",
ylab="VAR(Y|X)",op=2,xout=TRUE,eout=FALSE,pr=TRUE,fr=.6,scat=TRUE,outfun=out,SEED=TRUE){
#
# For each x, estimate VAR(y|x) using bootstrap bagging.
# with
# op=1 uses Fan's kernel method plus bootstrap bagging.
# op=2 uses running interval smoother plus bootstrap bagging
#
# xout=T eliminates points where there are outliers among x values
#        this option applies only when using op=2 and when using
#        running interval smoother.
# eout=T eliminates outliers among cloud of all data.
#
if(SEED)set.seed(2)
temp<-cbind(x,y)
temp<-elimna(temp)
x<-temp[,1]
y<-temp[,2]
if(op==2){
if(pr){
print("Running interval method plus bagging has been chosen")
print("op=1 will use Fan's method plus bagging")
}}
if(op==1){
if(pr){
print("Fan's method plus bagging has been chosen (cf. Bjerve and Doksum)")
print("op=2 will use running interval plus bagging")
}
mat <- matrix(NA, nrow = nboot, ncol = nrow(temp))
for(it in 1:nboot) {
idat <- sample(c(1:length(y)), replace = T)
xx <- temp[idat, 1]
yy <- temp[idat, 2]
mat[it,  ] <- locvar(xx,yy,pts=x,pyhat=TRUE,plotit=FALSE)
}
rmd<-apply(mat,2,mean)
 if(plotit) {
plot(c(x, x), c(y, rmd), type = "n", xlab = xlab, ylab= ylab)
sx <- sort(x)
xorder <- order(x)
sysm <- rmd[xorder]
lines(sx, sysm)
}

output<-"Done"
if(pyhat)output <- rmd
}
if(op==2){
output<-runmbo(x,y,fr=fr,est=var,xlab=xlab,ylab=ylab,pyhat=pyhat,eout=eout,
xout=xout,RNA=RNA,plotit=plotit,scat=scat,nboot=nboot,outfun=outfun,SEED=SEED)
}
output
}


# locCV
locCV<-function(y,varfun=pbvar,locfun=median,...){
vals=NA
n=length(y)
est=locfun(y)
for(i in 1:n)vals[i]=y[i]-locfun(y[-i],...)
res=varfun(vals)
res
}



# llocv2
llocv2<-function(x,est=median,...){
if(!is.list(x))val<-est(x,...)
if(is.list(x)){
val<-NA
for(i in 1:length(x))val[i]<-est(x[[i]],...)
}
if(is.matrix(x))val<-apply(x,2,est,...)
list(center=val)
}

# covloc
covloc<-function(x){
#
# Return mean and covarinace matrix
#
loc=apply(x,2,mean)
mcov=cov(x)
list(center=loc,cov=mcov)
}

# ogk.center
ogk.center<-function(x,beta=.9,...){
#
# Compute OGK multivariate measure of location
#
center=ogk(x,beta=beta,...)$center
center
}

# ghdist
ghdist<-function(n,g=0,h=0){
#
# generate n observations from a g-and-h dist.
#
x<-rnorm(n)
if (g>0){
ghdist<-(exp(g*x)-1)*exp(h*x^2/2)/g
}
if(g==0)ghdist<-x*exp(h*x^2/2)
ghdist
}



# meancr.cord.oph
meancr.cord.oph<-function(m,nullv=rep(0,ncol(m)),cop=3,MM=FALSE,SEED=TRUE,tr=0,
nboot=500,plotit=TRUE,MC=FALSE,xlab="VAR 1",ylab="VAR 2",STAND=TRUE){
#
# m is an n by p matrix
#
# Test hypothesis that the means
# are all equal to the null value, which defaults to zero.
# The level of the test is .05.
#
# Eliminate outliers using a projection method
# That is, determine center of data using:
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
# cop=4 MVE
#
# For each point
# consider the line between it and the center
# project all points onto this line, and
# check for outliers using
#
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# Eliminate any outliers and compute means
#  using remaining data.
#
if(SEED)set.seed(2)
m<-elimna(m)
n<-nrow(m)
#est=smean(m,MC=MC,cop=cop,STAND=STAND)
est=apply(m,2,mean,tr=tr)
crit.level<-.05
if(n<=120)crit.level<-.045
if(n<=80)crit.level<-.04
if(n<=60)crit.level<-.035
if(n<=40)crit.level<-.03
if(n<=30)crit.level<-.025
if(n<=20)crit.level<-.02
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
val<-matrix(NA,ncol=ncol(m),nrow=nboot)
for(j in 1: nboot){
mm<-m[data[j,],]
#val[j,]<-smean(mm,MC=MC,cop=cop,STAND=STAND)
val[j,]=apply(mm,2,mean)
}
if(!MC)temp<-pdis(rbind(val,nullv),center=est)
if(MC)temp<-pdisMC(rbind(val,nullv),center=est)
sig.level<-sum(temp[nboot+1]<temp[1:nboot])/nboot
list(p.value=sig.level,boot.vals=val,center=est)
}



# hpsi
hpsi<-function(x,bend=1.28){
#
#   Evaluate Huber`s Psi function for each value in the vector x
#   The bending constant defaults to 1.28.
#
hpsi<-ifelse(abs(x)<=bend,x,bend*sign(x))
hpsi
}



# lloc
 lloc<-function(x,est=tmean,...){
if(is.data.frame(x)){
x=as.matrix(x)
x=apply(x,2,as.numeric) # earlier versions of R require this command
}
if(!is.list(x))val<-est(x,...)
if(is.list(x))val=lapply(x,est,...)
if(is.matrix(x))val<-apply(x,2,est,...)
val
}
