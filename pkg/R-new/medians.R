# Median-based Methods
#
# Functions for median-based inference, tests, and comparisons
# Extracted from WRS package Rallfun-v45.R
#
# Created: 2025-12-30 


# ==== msmed ====
#' Test Linear Contrasts Using Medians
#'
#' @description
#' Performs hypothesis tests for a set of linear contrasts using sample medians
#' with heteroscedastic standard errors based on the McKean-Schrader method.
#' Tests whether median-based linear contrasts equal zero for independent groups.
#'
#' @param x Data in matrix or list mode. If a matrix, columns represent groups.
#'   If a list, each element is a group. Can also be a vector if y is specified.
#' @param y Optional second group when comparing exactly two groups. Default is NA.
#' @param con A J × d matrix of contrast coefficients where J is the number of groups
#'   and d is the number of contrasts. If con=0 (default), all pairwise comparisons
#'   are performed.
#' @param alpha Nominal alpha level for confidence intervals (default: 0.05).
#'
#' @return A list with components:
#' \item{test}{Matrix with columns: Group(s), test statistic, critical value,
#'   standard error, p-value, and adjusted p-value.}
#' \item{psihat}{Matrix with columns: Group(s)/contrast number, estimate,
#'   CI lower bound, CI upper bound.}
#'
#' @details
#' The function uses sample medians and McKean-Schrader standard errors to test
#' linear contrasts among independent groups. When con=0, all pairwise comparisons
#' are performed automatically. The test statistic follows an asymptotic normal
#' distribution under the null hypothesis.
#'
#' When tied values are detected, a warning is issued suggesting the use of
#' \code{medpb} (bootstrap method) instead, as ties can affect the validity
#' of the McKean-Schrader standard errors.
#'
#' For multiple comparisons, the critical value is adjusted using the maximum
#' modulus distribution, and both unadjusted and adjusted p-values are provided.
#'
#' @references
#' McKean, J. W., & Schrader, R. M. (1984). A comparison of methods for
#' studentizing the sample median. Communications in Statistics-Simulation
#' and Computation, 13(6), 751-773.
#'
#' @seealso \code{\link{medpb}}, \code{\link{med2g}}, \code{\link{med1way}}
#'
#' @examples
#' # Compare three independent groups
#' set.seed(123)
#' x1 <- rnorm(20, mean = 0)
#' x2 <- rnorm(20, mean = 0.5)
#' x3 <- rnorm(20, mean = 1)
#' data <- list(x1, x2, x3)
#'
#' # All pairwise comparisons
#' msmed(data)
#'
#' # Custom contrast: compare group 1 vs average of groups 2 and 3
#' con <- matrix(c(1, -0.5, -0.5), ncol = 1)
#' msmed(data, con = con)
#'
#' @export
msmed<-function(x,y=NA,con=0,alpha=.05){
#
# Test a set of linear contrasts using Medians
#
#  The data are assumed to be stored in $x$ in a matrix or in list mode.
#  Length(x) is assumed to correspond to the total number of groups, J
#  It is assumed all groups are independent.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  If con is not specified, all pairwise comparisons are made.
#
#  Missing values are automatically removed.
#
if(!is.na(y[1])){
xx<-list()
xx[[1]]<-x
xx[[2]]<-y
if(is.matrix(x) || is.list(x))stop("When y is speficied, x should not have list mode or be a matrix")
x<-xx
}
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-length(x)
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
if(sum(duplicated(val)>0)){
print(paste("Warning: Group",j, "has tied values. Might want to used medpb"))
}
x[[j]]<-val[xx]  # Remove missing values
xbar[j]<-median(x[[j]])
w[j]<-msmedse(x[[j]])^2 # Squared standard error.
}
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0){
CC<-(J^2-J)/2
psihat<-matrix(0,CC,5)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
test<-matrix(NA,CC,7)
dimnames(test)<-list(NULL,c("Group","Group","test","crit","se","p.value",'p.adjusted'))
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
test[jcom,3]<-abs(xbar[j]-xbar[k])/sqrt(w[j]+w[k])
test[jcom,6]<-2*(1-pt(test[jcom,3],999))
test[jcom,7]=1-psmm(abs(test[jcom,3]),CC,500)
sejk<-sqrt(w[j]+w[k])
test[jcom,5]<-sejk
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-(xbar[j]-xbar[k])
crit<-NA
if(CC==1)crit<-qnorm(1-alpha/2)
if(CC>1){
crit=qsmm(1-alpha,CC,500)
}
test[jcom,4]<-crit
psihat[jcom,4]<-psihat[jcom,3]-crit*test[jcom,5]
psihat[jcom,5]<-psihat[jcom,3]+crit*test[jcom,5]
}}}}
if(sum(con^2)>0){
if(nrow(con)!=length(x))warning("The number of groups does not match the number of contrast coefficients.")
psihat<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol(con),6)
dimnames(test)<-list(NULL,c("con.num","test","crit","se","p.value",'p.adjusted'))
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-sqrt(sum(con[,d]^2*w))
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
test[d,5]<-2*(1-pt(abs(test[d,2]),999))
test[d,6]=1-psmm(abs(test[d,2]),ncol(con),500)
crit<-NA
if(CC==1)crit<-qnorm(1-alpha/2)
if(CC>1)crit=qsmm(1-alpha,CC,500)
test[d,3]<-crit
test[d,4]<-sejk
psihat[d,3]<-psihat[d,2]-crit*sejk
psihat[d,4]<-psihat[d,2]+crit*sejk
}}
list(test=test,psihat=psihat)
}

# ==== med2way ====
#' Two-Way ANOVA for Medians with Independent Groups
#'
#' @description
#' Performs a J × K (two-way) ANOVA for medians with all groups independent,
#' allowing for heteroscedasticity. Tests main effects for factors A and B,
#' as well as the A × B interaction.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in list mode. Length should equal J×K, with x[[1]] containing
#'   data for level (1,1), x[[2]] for level (1,2), etc.
#' @param grp Vector specifying the order of groups if data are not in the
#'   default order. For example, grp=c(2,4,3,1) for a 2×2 design indicates
#'   x[[2]] is level (1,1), x[[4]] is level (1,2), etc.
#' @param p Total number of groups (should equal J×K). Used for validation.
#' @param ADJ.P.VALUE Logical; if TRUE (default), compute adjusted p-values
#'   via simulation. If FALSE, use asymptotic chi-square/F distributions.
#' @param iter Number of simulation iterations for adjusted p-values when
#'   ADJ.P.VALUE=TRUE (default: 5000).
#' @inheritParams common_params
#'
#' @return A list with components:
#' \item{test.A}{Test statistic for Factor A main effect.}
#' \item{p.val.A}{P-value for Factor A.}
#' \item{test.B}{Test statistic for Factor B main effect.}
#' \item{p.val.B}{P-value for Factor B.}
#' \item{test.AB}{Test statistic for A × B interaction.}
#' \item{p.val.AB}{P-value for interaction.}
#'
#' @details
#' This function extends \code{\link{med1way}} to two-way factorial designs
#' with independent groups. It uses McKean-Schrader standard errors for medians
#' and constructs test statistics that allow for heteroscedasticity.
#'
#' For Factor A and B main effects, weighted test statistics are constructed
#' similar to one-way ANOVA for medians. The interaction test uses a chi-square
#' statistic based on deviations from the additive model.
#'
#' When \code{ADJ.P.VALUE=TRUE}, p-values are computed via simulation to account
#' for the sampling distribution of the test statistics under the null hypothesis.
#' When FALSE, asymptotic distributions are used (F for main effects, chi-square
#' for interaction).
#'
#' A warning is issued if tied values are detected, suggesting the use of
#' \code{m2way} (bootstrap method) instead.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{med1way}}, \code{\link{m2way}}, \code{\link{t2way}}
#'
#' @examples
#' # Two-way ANOVA for medians
#' set.seed(123)
#' # 2x3 design: 2 levels of A, 3 levels of B
#' data <- list(
#'   rnorm(15, mean = 0),    # A1,B1
#'   rnorm(15, mean = 0.5),  # A1,B2
#'   rnorm(15, mean = 0.3),  # A1,B3
#'   rnorm(15, mean = 1),    # A2,B1
#'   rnorm(15, mean = 1.5),  # A2,B2
#'   rnorm(15, mean = 1.2)   # A2,B3
#' )
#'
#' med2way(J = 2, K = 3, x = data, iter = 2000)
#'
#' @export
med2way<-function(J,K,x,grp=c(1:p),p=J*K, ADJ.P.VALUE=TRUE, iter=5000,SEED=TRUE){
#
#  Perform a J by K (two-way) anova on  medians where
#  all jk groups are independent.
#
#  The argument x is assumed to contain the raw
#  data stored in list mode.
#  If grp is unspecified, it is assumed x[[1]] contains the data
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second factor: level 1,2
#  x[[j+1]] is the data for level 2,1, etc.
#  If the data are in wrong order, grp can be used to rearrange the
#  groups. For example, for a two by two design, grp<-c(2,4,3,1)
#  indicates that the second group corresponds to level 1,1;
#  group 4 corresponds to level 1,2; group 3 is level 2,1;
#  and group 1 is level 2,2.
#
#  It is assumed that the input variable x has length JK, the total number of
#  groups being tested. If not, a warning message is printed.
#
if(L.ties(x))print("There are tied values, suggest using the function m2way instead")
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data are not stored in a matrix or in list mode")
if(p!=length(x)){
print("Warning: The number of groups in your data is not equal to JK")
}
xbar<-0
h<-0
d<-0
R<-0
W<-0
d<-0
r<-0
w<-0
nuhat<-0
omegahat<-0
DROW<-0
DCOL<-0
xtil<-matrix(0,J,K)
aval<-matrix(0,J,K)
for (j in 1:p){
xbar[j]<-median(x[[grp[j]]])
h[j]<-length(x[[grp[j]]])
d[j]<-msmedse(x[[grp[j]]])^2
}
d<-matrix(d,J,K,byrow=TRUE)
xbar<-matrix(xbar,J,K,byrow=TRUE)
h<-matrix(h,J,K,byrow=TRUE)
for(j in 1:J){
R[j]<-sum(xbar[j,])
nuhat[j]<-(sum(d[j,]))^2/sum(d[j,]^2/(h[j,]-1))
r[j]<-1/sum(d[j,])
DROW[j]<-sum(1/d[j,])
}
for(k in 1:K){
W[k]<-sum(xbar[,k])
omegahat[k]<-(sum(d[,k]))^2/sum(d[,k]^2/(h[,k]-1))
w[k]<-1/sum(d[,k])
DCOL[k]<-sum(1/d[,k])
}
D<-1/d
for(j in 1:J){
for(k in 1:K){
xtil[j,k]<-sum(D[,k]*xbar[,k]/DCOL[k])+sum(D[j,]*xbar[j,]/DROW[j])-
sum(D*xbar/sum(D))
aval[j,k]<-(1-D[j,k]*(1/sum(D[j,])+1/sum(D[,k])-1/sum(D)))^2/(h[j,k]-3)
}
}
Rhat<-sum(r*R)/sum(r)
What<-sum(w*W)/sum(w)
Ba<-sum((1-r/sum(r))^2/nuhat)
Bb<-sum((1-w/sum(w))^2/omegahat)
Va<-sum(r*(R-Rhat)^2)/((J-1)*(1+2*(J-2)*Ba/(J^2-1)))
Vb<-sum(w*(W-What)^2)/((K-1)*(1+2*(K-2)*Bb/(K^2-1)))
sig.A<-1-pf(Va,J-1,9999999)
sig.B<-1-pf(Vb,K-1,9999999)
# Next, do test for interactions
Vab<-sum(D*(xbar-xtil)^2)
dfinter<-(J-1)*(K-1)
sig.AB<-1-pchisq(Vab,dfinter)
if(ADJ.P.VALUE){
a=med2way.crit(J,K,h,iter=iter,SEED=SEED)
sig.A=mean(Va<=a$A.dist)
sig.B=mean(Vb<=a$B.dist)
sig.AB=mean(Vab<=a$AB.dist)
}
list(test.A=Va,p.val.A=sig.A,test.B=Vb,p.val.B=sig.B,test.AB=Vab,p.val.AB=sig.AB)
}

# ==== med2way.sub ====
#' Helper Function for Two-Way Median ANOVA (Internal)
#'
#' @description
#' Internal helper function used by \code{\link{med2way.crit}} to compute
#' test statistics for two-way ANOVA on medians. Not intended for direct use.
#'
#' @inheritParams med2way
#'
#' @return A list with test statistics and p-values for Factor A, Factor B,
#'   and A × B interaction (same structure as \code{\link{med2way}}).
#'
#' @details
#' This function performs the same computations as \code{\link{med2way}} but
#' without p-value adjustments or simulation options. It is called repeatedly
#' during simulation to determine critical values for the adjusted p-values.
#'
#' @seealso \code{\link{med2way}}, \code{\link{med2way.crit}}
#'
#' @keywords internal
med2way.sub<-function(J,K,x,grp=c(1:p),p=J*K){
#
#  Perform a J by K (two-way) anova on  medians where
#  all jk groups are independent.
#
#  The argument x is assumed to contain the raw
#  data stored in list mode.
#  If grp is unspecified, it is assumed x[[1]] contains the data
#  for the first level of both factors: level 1,1.
#  x[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second factor: level 1,2
#  x[[j+1]] is the data for level 2,1, etc.
#  If the data are in wrong order, grp can be used to rearrange the
#  groups. For example, for a two by two design, grp<-c(2,4,3,1)
#  indicates that the second group corresponds to level 1,1;
#  group 4 corresponds to level 1,2; group 3 is level 2,1;
#  and group 1 is level 2,2.
#
#  It is assumed that the input variable x has length JK, the total number of
#  groups being tested. If not, a warning message is printed.
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data are not stored in a matrix or in list mode')
if(p!=length(x)){
print('Warning: The number of groups in your data is not equal to JK')
}
xbar<-0
h<-0
d<-0
R<-0
W<-0
d<-0
r<-0
w<-0
nuhat<-0
omegahat<-0
DROW<-0
DCOL<-0
xtil<-matrix(0,J,K)
aval<-matrix(0,J,K)
for (j in 1:p){
xbar[j]<-median(x[[grp[j]]])
h[j]<-length(x[[grp[j]]])
d[j]<-msmedse(x[[grp[j]]])^2
}
d<-matrix(d,J,K,byrow=TRUE)
xbar<-matrix(xbar,J,K,byrow=TRUE)
h<-matrix(h,J,K,byrow=TRUE)
for(j in 1:J){
R[j]<-sum(xbar[j,])
nuhat[j]<-(sum(d[j,]))^2/sum(d[j,]^2/(h[j,]-1))
r[j]<-1/sum(d[j,])
DROW[j]<-sum(1/d[j,])
}
for(k in 1:K){
W[k]<-sum(xbar[,k])
omegahat[k]<-(sum(d[,k]))^2/sum(d[,k]^2/(h[,k]-1))
w[k]<-1/sum(d[,k])
DCOL[k]<-sum(1/d[,k])
}
D<-1/d
for(j in 1:J){
for(k in 1:K){
xtil[j,k]<-sum(D[,k]*xbar[,k]/DCOL[k])+sum(D[j,]*xbar[j,]/DROW[j])-
sum(D*xbar/sum(D))
aval[j,k]<-(1-D[j,k]*(1/sum(D[j,])+1/sum(D[,k])-1/sum(D)))^2/(h[j,k]-3)
}
}
Rhat<-sum(r*R)/sum(r)
What<-sum(w*W)/sum(w)
Ba<-sum((1-r/sum(r))^2/nuhat)
Bb<-sum((1-w/sum(w))^2/omegahat)
Va<-sum(r*(R-Rhat)^2)/((J-1)*(1+2*(J-2)*Ba/(J^2-1)))
Vb<-sum(w*(W-What)^2)/((K-1)*(1+2*(K-2)*Bb/(K^2-1)))
sig.A<-1-pf(Va,J-1,9999999)
sig.B<-1-pf(Vb,K-1,9999999)
# Next, do test for interactions
Vab<-sum(D*(xbar-xtil)^2)
dfinter<-(J-1)*(K-1)
sig.AB<-1-pchisq(Vab,dfinter)
list(test.A=Va,p.val.A=sig.A,test.B=Vb,p.val.B=sig.B,test.AB=Vab,p.val.AB=sig.AB)
}

# ==== med2way.crit ====
#' Estimate Null Distribution for Two-Way Median ANOVA (Internal)
#'
#' @description
#' Internal helper function that simulates the null distribution for
#' \code{\link{med2way}} test statistics via Monte Carlo. Used to compute
#' adjusted p-values. Not intended for direct use.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param n Matrix (J×K) of sample sizes for each group.
#' @param iter Number of simulation iterations.
#' @inheritParams common_params
#'
#' @return A list with components:
#' \item{A.dist}{Vector of simulated test statistics for Factor A under the null.}
#' \item{B.dist}{Vector of simulated test statistics for Factor B under the null.}
#' \item{AB.dist}{Vector of simulated test statistics for A × B interaction under the null.}
#'
#' @details
#' For each iteration, data are generated from a standard normal distribution
#' with sample sizes matching the original data. Test statistics are computed
#' using \code{\link{med2way.sub}}. The resulting empirical distributions are
#' used by \code{\link{med2way}} to compute adjusted p-values that account for
#' the sampling distribution under the null hypothesis.
#'
#' @seealso \code{\link{med2way}}, \code{\link{med2way.sub}}
#'
#' @keywords internal
med2way.crit<-function(J,K,n,iter,SEED=TRUE){
#
# Estimate the null distribution for med2way
#
x=list()
p=J*K
A.dist=NA
B.dist=NA
AB.dist=NA
for(i in 1:iter){
for(j in 1:p)x[[j]]=rmul(n[j])
a=med2way.sub(J,K,x)
A.dist[i]=a$test.A
B.dist[i]=a$test.B
AB.dist[i]=a$test.AB
}
list(A.dist=A.dist,B.dist=B.dist,AB.dist=AB.dist)
}

# ==== med1way ====
#' Heteroscedastic One-Way ANOVA for Medians
#'
#' @description
#' Performs a one-way ANOVA for comparing medians across J independent groups,
#' allowing for heteroscedasticity (unequal variances). Uses a test statistic
#' based on weighted medians with McKean-Schrader standard errors.
#'
#' @inheritParams common_params
#' @param x Data in matrix or list mode. If a matrix, columns represent groups.
#'   If a list, each element is a group.
#' @param crit Optional critical value. If NA (default), the critical value
#'   is determined via simulation for the specified alpha level. Specifying
#'   a value reduces execution time but no p-value is computed.
#' @param iter Number of simulation iterations for determining the critical
#'   value when crit=NA (default: 5000).
#' @param pr Logical; if TRUE, print progress messages (default: TRUE).
#'
#' @return A list with components:
#' \item{TEST}{Test statistic value.}
#' \item{crit.val}{Critical value for the test.}
#' \item{p.value}{P-value (only when crit=NA).}
#'
#' @details
#' The test statistic is a weighted sum of squared deviations of group medians
#' from the weighted grand median, divided by J-1. Weights are the inverse
#' squared McKean-Schrader standard errors.
#'
#' When \code{crit=NA}, the critical value is determined via simulation under
#' the null hypothesis of equal medians, assuming normality with group-specific
#' variances estimated from the data. The simulation approach accounts for the
#' dependence of the critical value on sample sizes.
#'
#' The test allows for heteroscedasticity, making it more robust than classical
#' ANOVA when variances are unequal across groups.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{med2way}}, \code{\link{msmed}}, \code{\link{t1way}}
#'
#' @examples
#' # One-way ANOVA for medians
#' set.seed(123)
#' g1 <- rnorm(20, mean = 0, sd = 1)
#' g2 <- rnorm(25, mean = 0.5, sd = 1.5)
#' g3 <- rnorm(30, mean = 1, sd = 2)
#' data <- list(g1, g2, g3)
#'
#' med1way(data, iter = 2000)
#'
#' # Compare subset of groups
#' med1way(data, grp = c(1, 3))
#'
#' @export
med1way<-function(x,grp=NA,alpha=.05,crit=NA,iter=5000,SEED=TRUE,pr=TRUE){
#
#  A heteroscedastic one-way ANOVA for medians.
#
#  If
#  crit=NA, an appropriate critical value is determined
#  for the alpha value specified and
#  a p-value is returned. This is done via a simulation assuming normality.
#  iter: the number of replications in the simulation.
#
#  If a value for crit is specified, it is used as the critical
#  value, but no p-value is reported. Specifying a value for
#  crit reduces execution time.
#  With crit=NA, the critical value is a function of the sample
#  sizes and is determined by calling the function med1way.crit.
#
#  The data are assumed to be stored in $x$ in list mode.
#  length(x) is assumed to correspond to the total number of groups.
#  By default, the null hypothesis is that all groups have a common mean.
#  To compare a subset of the groups, use grp to indicate which
#  groups are to be compared. For example, if you type the
#  command grp<-c(1,3,4), and then execute this function, groups
#  1, 3, and 4 will be compared with the remaining groups ignored.
#
#  Missing values are automatically removed.
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x)){
y<-list()
for(j in 1:ncol(x))y[[j]]<-x[,j]
x<-y
}
if(is.na(sum(grp[1])))grp<-c(1:length(x))
if(!is.list(x))stop("Data are not stored in a matrix or in list mode.")
J<-length(grp)  # The number of groups to be compared
n<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
w[j]<-1/msmedse(x[[grp[j]]])^2
xbar[j]<-median(x[[grp[j]]])
n[j]<-length(x[[grp[j]]])
}
pval<-NA
u<-sum(w)
xtil<-sum(w*xbar)/u
TEST<-sum(w*(xbar-xtil)^2)/(J-1)
if(is.na(crit)){
temp<-med1way.crit(n,alpha,SEED=SEED,iter=iter,TEST=TEST)
crit.val<-temp$crit.val
}
if(!is.na(crit))crit.val<-crit
list(TEST=TEST,crit.val=crit.val,p.value=temp$p.value)
}

# ==== med1way.crit ====
#' Determine Critical Value for One-Way Median ANOVA (Internal)
#'
#' @description
#' Internal helper function that determines critical values for \code{\link{med1way}}
#' via simulation under normality. Used to account for the dependence of critical
#' values on sample sizes. Not intended for direct use.
#'
#' @param n Vector of sample sizes for each group.
#' @inheritParams common_params
#' @param iter Number of simulation iterations (default: 1000).
#' @param TEST Optional observed test statistic value. If provided, a p-value
#'   is computed (default: NA).
#'
#' @return A list with components:
#' \item{crit.val}{Critical value for the specified alpha level.}
#' \item{p.value}{P-value if TEST is specified, otherwise NA.}
#'
#' @details
#' For each iteration, data are generated from standard normal distributions
#' with sample sizes matching the input. The test statistic from \code{\link{med1way}}
#' is computed for each simulated dataset. The critical value is the (1-alpha)
#' quantile of the simulated distribution.
#'
#' If an observed test statistic is provided via TEST, a p-value is computed
#' as the proportion of simulated statistics that are less than or equal to
#' the observed value.
#'
#' @seealso \code{\link{med1way}}
#'
#' @keywords internal
med1way.crit<-function(n,alpha=.05,iter=1000,TEST=NA,SEED=TRUE){
#
#  Determine the critical value for the function
#  med1way, assuming normality, based on the sample
#  sizes in n.
#
J<-length(n)
x<-list()
w<-vector("numeric",J)
xbar<-vector("numeric",J)
if(SEED)set.seed(2)
chk<-NA
grp<-c(1:J)
for (it in 1:iter){
for(j in 1:J){
x[[j]]<-rnorm(n[j])
w[j]<-1/msmedse(x[[grp[j]]])^2
xbar[j]<-median(x[[grp[j]]])
n[j]<-length(x[[grp[j]]])
}
u<-sum(w)
xtil<-sum(w*xbar)/u
chk[it]<-sum(w*(xbar-xtil)^2)/(J-1)
}
chk<-sort(chk)
iv<-round((1-alpha)*iter)
crit.val<-chk[iv]
pval<-NA
if(!is.na(TEST))pval<-sum((TEST<=chk))/iter
list(crit.val=crit.val,p.value=pval)
}

# ==== bpmed ====
#' Bootstrap Percentile Method for Median Contrasts (Bonett-Price SE)
#'
#' @description
#' Tests linear contrasts among medians of independent groups using the
#' Bonett-Price standard error estimator. Similar to \code{\link{msmed}} but
#' uses a different method for estimating standard errors.
#'
#' @inheritParams common_params
#' @param x Data in matrix or list mode. If a matrix, columns represent groups.
#'   If a list, each element is a group.
#' @param con A J × d matrix of contrast coefficients. If con=0 (default),
#'   all pairwise comparisons are performed.
#'
#' @return A list with components:
#' \item{test}{Matrix with columns: Group(s), test statistic, critical value,
#'   and standard error.}
#' \item{psihat}{Matrix with columns: Group(s)/contrast number, estimate,
#'   CI lower bound, CI upper bound.}
#'
#' @details
#' This function is similar to \code{\link{msmed}} but uses the Bonett-Price
#' method for estimating the standard error of the median, implemented in
#' \code{\link{bpmedse}}. The Bonett-Price method may provide better coverage
#' in small samples compared to the McKean-Schrader method.
#'
#' Test statistics are computed as the difference in medians divided by the
#' pooled standard error. For multiple comparisons, critical values are
#' adjusted using the studentized maximum modulus distribution.
#'
#' Missing values are automatically removed before analysis.
#'
#' @references
#' Price, R. M., & Bonett, D. G. (2001). Estimating the variance of the
#' sample median. Journal of Statistical Computation and Simulation,
#' 68(3), 295-305.
#'
#' @seealso \code{\link{msmed}}, \code{\link{bpmedse}}, \code{\link{medpb}}
#'
#' @examples
#' # Compare three groups using Bonett-Price SE
#' set.seed(123)
#' x1 <- rnorm(20, mean = 0)
#' x2 <- rnorm(20, mean = 0.5)
#' x3 <- rnorm(20, mean = 1)
#' data <- list(x1, x2, x3)
#'
#' # All pairwise comparisons
#' bpmed(data)
#'
#' # Custom contrast
#' con <- matrix(c(1, -0.5, -0.5), ncol = 1)
#' bpmed(data, con = con)
#'
#' @export
bpmed<-function(x,con=0,alpha=.05){
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-length(x)
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
xbar[j]<-median(x[[j]])
w[j]<-bpmedse(x[[j]])^2 # Squared standard error.
}
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0){
CC<-(J^2-J)/2
psihat<-matrix(0,CC,5)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
test<-matrix(NA,CC,5)
dimnames(test)<-list(NULL,c("Group","Group","test","crit","se"))
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
test[jcom,3]<-abs(xbar[j]-xbar[k])/sqrt(w[j]+w[k])
sejk<-sqrt(w[j]+w[k])
test[jcom,5]<-sejk
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-(xbar[j]-xbar[k])
crit<-NA
if(CC==1)crit<-qnorm(1-alpha/2)
if(CC>1){
#if(alpha==.05)crit<-smmcrit(500,CC)
#if(alpha==.01)crit<-smmcrit01(500,CC)
#if(is.na(crit))warning("Can only be used with alpha=.05 or .01")
crit=qsmm(1-alpha,CC,500)
}
test[jcom,4]<-crit
psihat[jcom,4]<-psihat[jcom,3]-crit*test[jcom,5]
psihat[jcom,5]<-psihat[jcom,3]+crit*test[jcom,5]
}}}}
if(sum(con^2)>0){
if(nrow(con)!=length(x))warning("The number of groups does not match the number of contrast coefficients.")
psihat<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c("con.num","test","crit","se","df"))
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-sqrt(sum(con[,d]^2*w))
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
crit<-NA
if(CC==1)crit<-qnorm(1-alpha/2)
#if(alpha==.05)crit<-smmcrit(500,ncol(con))
#if(alpha==.01)crit<-smmcrit01(500,ncol(con))
crit=qsmm(1-alpha,ncol(con),500)
test[d,3]<-crit
test[d,4]<-sejk
psihat[d,3]<-psihat[d,2]-crit*sejk
psihat[d,4]<-psihat[d,2]+crit*sejk
}}
list(test=test,psihat=psihat)
}

# ==== bpmedse ====
#' Bonett-Price Standard Error for the Median
#'
#' @description
#' Computes the standard error of the sample median using the method
#' recommended by Price and Bonett (2001).
#'
#' @param x Numeric vector of data.
#'
#' @return Numeric standard error estimate for the median.
#'
#' @details
#' The Bonett-Price method estimates the standard error by constructing
#' a confidence interval for the median and converting it to a standard
#' error estimate. It uses the binomial distribution to determine the
#' appropriate order statistics, then computes:
#'
#' SE = (U - L) / (2 * z),
#'
#' where U and L are the upper and lower confidence bounds from order
#' statistics, and z is the normal critical value corresponding to the
#' confidence level.
#'
#' This method tends to provide better coverage in small to moderate
#' samples compared to asymptotic methods.
#'
#' Missing values are automatically removed before computation.
#'
#' @references
#' Price, R. M., & Bonett, D. G. (2001). Estimating the variance of the
#' sample median. Journal of Statistical Computation and Simulation,
#' 68(3), 295-305.
#'
#' @seealso \code{\link{bpmed}}, \code{\link{msmedse}}
#'
#' @examples
#' # Standard error for median
#' set.seed(123)
#' x <- rnorm(50, mean = 5, sd = 2)
#' bpmedse(x)
#'
#' @export
bpmedse<-function(x){
#
# compute standard error of the median using method
# recommended by Price and Bonett (2001)
#
y<-sort(x)
n<-length(x)
av<-round((n+1)/2-sqrt(n))
if(av==0)av<-1
avm<-av-1
astar<-pbinom(avm,n,.5)  #alpha*/2
zval<-qnorm(1-astar)
top<-n-av+1
sqse<-((y[top]-y[av])/(2*zval))^2 # The sq. standard error
se<-sqrt(sqse)
se
}

# ==== exmed ====
#' Exact Median Contrasts with Simulation-Based Critical Values
#'
#' @description
#' Tests linear contrasts among medians of independent groups using simulation
#' to determine critical values that provide exact familywise error control
#' under normality (when iteration count is sufficiently large).
#'
#' @inheritParams common_params
#' @param x Data in matrix or list mode. If a matrix, columns represent groups.
#'   If a list, each element is a group. Can also be a vector if y is specified.
#' @param y Optional second group when comparing exactly two groups. Default is NA.
#' @param con A J × d matrix of contrast coefficients. If con=0 (default),
#'   all pairwise comparisons are performed.
#' @param iter Number of simulation iterations for determining critical values
#'   (default: 1000). Larger values provide more accurate control of Type I error.
#' @param se.fun Function for estimating standard errors. Default is
#'   \code{\link{bpmedse}} (Bonett-Price method). Can also use
#'   \code{\link{msmedse}} (McKean-Schrader method).
#'
#' @return A list with components:
#' \item{test}{Matrix with columns: Group(s), test statistic, critical value,
#'   standard error, and p-value (unadjusted).}
#' \item{psihat}{Matrix with columns: Group(s)/contrast number, estimate,
#'   CI lower bound (adjusted for FWE), CI upper bound (adjusted for FWE).}
#'
#' @details
#' This function provides exact control over the familywise error rate (FWE)
#' under normality by simulating the null distribution of the test statistics.
#' For each iteration, data are generated from standard normal distributions
#' with sample sizes matching the original data, and the maximum absolute
#' test statistic is recorded. The (1-alpha) quantile of this distribution
#' is used as the critical value.
#'
#' \strong{Advantages}:
#' \itemize{
#'   \item Exact FWE control under normality (assuming sufficient iterations)
#'   \item Accounts for dependence among test statistics
#'   \item Provides both adjusted CIs and unadjusted p-values
#' }
#'
#' The confidence intervals are simultaneous intervals with coverage 1-alpha.
#' Individual p-values are provided but are not adjusted for multiplicity.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{msmed}}, \code{\link{bpmed}}, \code{\link{medpb}}
#'
#' @examples
#' # Compare three groups with exact FWE control
#' set.seed(123)
#' x1 <- rnorm(20, mean = 0)
#' x2 <- rnorm(20, mean = 0.5)
#' x3 <- rnorm(20, mean = 1)
#' data <- list(x1, x2, x3)
#'
#' # All pairwise comparisons (500 iterations for speed)
#' exmed(data, iter = 500)
#'
#' # Using McKean-Schrader SE instead
#' exmed(data, se.fun = msmedse, iter = 500)
#'
#' @export
exmed<-function(x,y=NA,con=0,alpha=.05,iter=1000,se.fun=bpmedse,SEED=TRUE){
#
# Test a set of linear contrasts using medians
#
# Get exact control over type I errors under normality, provided
# iter is sufficietly large.
# iter determines number of replications used in a simulation
# to determine critical value.
#
# se.fun indicates method used to estimate standard errors.
# default is the method used by Bonett and Price (2002)
# To use the McKean-Shrader method,
# set se.fun=msmedse
#
#  The data are assumed to be stored in $x$ in a matrix or in list mode.
#  Length(x) is assumed to correspond to the total number of groups, J
#  It is assumed all groups are independent.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  If con is not specified, all pairwise comparisons are made.
#
#  Missing values are automatically removed.
#
#  Function returns the critical value used so that FWE=alpha
#  (under the column crit)
#  p-values are determined for each test but are not adjusted so
#  that FWE=alpha.
#  The confidence intervals are adjusted so that the simultaneous
#  probability coverage is 1-alpha.
#
if(!is.na(y[1])){
xx<-list()
xx[[1]]<-x
xx[[2]]<-y
if(is.matrix(x) || is.list(x))stop("When y is speficied, x should not have list mode or be a matrix")
x<-xx
}
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-length(x)
h<-vector("numeric",J)
w<-vector("numeric",J)
nval<-vector("numeric",J)
xbar<-vector("numeric",J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
xbar[j]<-median(x[[j]])
nval[j]<-length(x[[j]])
# w[j]<-msmedse(x[[j]])^2
 w[j]<-se.fun(x[[j]])^2
}
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0){
CC<-(J^2-J)/2
psihat<-matrix(0,CC,5)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c("Group","Group","test","crit","se","p.value"))
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
test[jcom,3]<-abs(xbar[j]-xbar[k])/sqrt(w[j]+w[k])
# Next determine p-value for each individual test
temp<-msmedsub(c(nval[j],nval[k]),se.fun=se.fun,SEED=SEED,iter=iter)
test[jcom,6]<-sum((test[jcom,3]<=temp))/iter
sejk<-sqrt(w[j]+w[k])
test[jcom,5]<-sejk
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-(xbar[j]-xbar[k])
# Determine critical value for controlling FWE
temp<-msmedsub(nval,se.fun=se.fun,SEED=SEED,iter=iter)
ic<-round((1-alpha)*iter)
crit<-temp[ic]
test[jcom,4]<-crit
psihat[jcom,4]<-psihat[jcom,3]-crit*test[jcom,5]
psihat[jcom,5]<-psihat[jcom,3]+crit*test[jcom,5]
}}}}
if(sum(con^2)>0){
if(nrow(con)!=length(x))warning("The number of groups does not match the number of contrast coefficients.")
psihat<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c("con.num","test","crit","se","p.value"))
# Determine critical value that controls FWE
temp<-msmedsub(nval,con=con,se.fun=se.fun,SEED=SEED,iter=iter)
ic<-round((1-alpha)*iter)
crit<-temp[ic]
for (d in 1:ncol(con)){
flag<-(con[,d]==0)
nvec<-nval[!flag]
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-sqrt(sum(con[,d]^2*w))
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
#  Determine p-value for individual (dth) test
temp<-msmedsub(nvec,iter=iter,se.fun=se.fun,SEED=SEED)
test[d,3]<-crit
test[d,4]<-sejk
test[d,5]<-sum(abs((test[d,2])<=temp))/iter
psihat[d,3]<-psihat[d,2]-crit*sejk
psihat[d,4]<-psihat[d,2]+crit*sejk
}}
list(test=test,psihat=psihat)
}

# ==== msmedsub ====
#' Simulate Studentized Critical Values for Median Contrasts (Internal)
#'
#' @description
#' Internal helper function for \code{\link{exmed}} that simulates studentized
#' critical values for median contrasts under normality and homoscedasticity.
#' Not intended for direct use.
#'
#' @param n Vector of sample sizes for each group.
#' @param con A J × d matrix of contrast coefficients. If con=0 (default),
#'   all pairwise comparisons are performed.
#' @inheritParams common_params
#' @param se.fun Function for estimating standard errors (default: \code{\link{bpmedse}}).
#' @param iter Number of simulation iterations (default: 1000).
#'
#' @return A sorted vector of maximum absolute studentized test statistics
#'   from the simulated null distribution.
#'
#' @details
#' For each iteration, data are generated from standard normal distributions
#' with the specified sample sizes. Medians and standard errors are computed,
#' then test statistics are calculated for each contrast. The maximum absolute
#' test statistic is recorded. The returned vector can be used to determine
#' critical values at any alpha level.
#'
#' @seealso \code{\link{exmed}}
#'
#' @keywords internal
msmedsub<-function(n,con=0,alpha=.05,se.fun=bpmedse,iter=1000,SEED=TRUE){
#
# Determine a Studentized critical value, assuming normality
# and homoscedasticity, for the function msmedv2
#
# Goal: Test a set of linear contrasts using medians
#
#  The data are assumed to be stored in $x$ in a matrix or in list mode.
#  Length(x) is assumed to correspond to the total number of groups, J
#  It is assumed all groups are independent.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  If con is not specified, all pairwise comparisons are made.
#
if(SEED)set.seed(2)
con<-as.matrix(con)
J<-length(n)
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
x<-list()
test<-NA
testmax<-NA
for (it in 1:iter){
for(j in 1:J){
x[[j]]<-rnorm(n[j])
xbar[j]<-median(x[[j]])
 w[j]<-se.fun(x[[j]])^2
}
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0){
CC<-(J^2-J)/2
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
test[jcom]<-abs(xbar[j]-xbar[k])/sqrt(w[j]+w[k])
}}}}
if(sum(con^2)>0){
for (d in 1:ncol(con)){
sejk<-sqrt(sum(con[,d]^2*w))
test[d]<-sum(con[,d]*xbar)/sejk
}}
testmax[it]<-max(abs(test))
}
testmax<-sort(testmax)
testmax
}

# ==== medpb.old ====
#' Old Version of Bootstrap Median Comparisons (Deprecated)
#'
#' @description
#' Deprecated version of \code{\link{medpb}} using Rom's method for FWER control.
#' Retained for backward compatibility. Use \code{\link{medpb}} instead, which
#' offers more flexible FWER control methods.
#'
#' @inheritParams medpb
#'
#' @return Same as \code{\link{medpb}}.
#'
#' @details
#' This is the original implementation of \code{medpb} that used Rom's method
#' for familywise error rate control. The current version (\code{\link{medpb}})
#' defaults to Hochberg's method but offers multiple options via the \code{method}
#' parameter.
#'
#' The main difference is that this version does not include adjusted p-values
#' in the output, while the current version does.
#'
#' @seealso \code{\link{medpb}}
#'
#' @keywords internal
medpb.old<-function(x,alpha=.05,nboot=NA,grp=NA,est=median,con=0,bhop=FALSE,
SEED=TRUE,...){
#
#   Multiple comparisons for  J independent groups using medians.
#
#   A percentile bootstrap method is used. FWE is controlled with Rom's method.
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   est is the measure of location and defaults to the median
#   ... can be used to set optional arguments associated with est
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#
#   con can be used to specify linear contrasts; see the function lincon
#
#   Missing values are allowed.
#
con<-as.matrix(con)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
if(!is.na(sum(grp))){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
x<-xx
}
J<-length(x)
tempn<-0
mvec<-NA
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
mvec[j]<-est(temp,...)
}
Jm<-J-1
#
# Determine contrast matrix
#
if(sum(con^2)==0){
ncon<-(J^2-J)/2
con<-matrix(0,J,ncon)
id<-0
for (j in 1:Jm){
jp<-j+1
for (k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
ncon<-ncol(con)
dvec<-alpha/c(1:ncon)
if(nrow(con)!=J)stop("Something is wrong with con; the number of rows does not match the number of groups.")
#  Determine nboot if a value was not specified
if(is.na(nboot)){
nboot<-5000
if(J <= 8)nboot<-4000
if(J <= 3)nboot<-2000
}
# Determine critical values
if(!bhop){
if(alpha==.05){
dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
}
}
if(bhop)dvec<-(ncon-c(1:ncon)+1)*alpha/ncon
bvec<-matrix(NA,nrow=J,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
for(j in 1:J){
#print(paste("Working on group ",j))
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,est,...) # Bootstrapped values for jth group
}
test<-NA
bcon<-t(con)%*%bvec #ncon by nboot matrix
tvec<-t(con)%*%mvec
for (d in 1:ncon){
tv<-sum(bcon[d,]==0)/nboot
test[d]<-sum(bcon[d,]>0)/nboot+.5*tv
if(test[d]> .5)test[d]<-1-test[d]
}
test<-2*test
output<-matrix(0,ncon,6)
dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.crit","ci.lower","ci.upper"))
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,4]<-zvec
icl<-round(dvec[ncon]*nboot/2)+1
icu<-nboot-icl-1
for (ic in 1:ncol(con)){
output[ic,2]<-tvec[ic,]
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(bcon[ic,])
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
num.sig<-sum(output[,3]<=output[,4])
list(output=output,con=con,num.sig=num.sig)
}

# ==== medpb ====
#' Multiple Comparisons for Independent Groups Using Bootstrap Medians
#'
#' @description
#' Performs multiple comparisons among J independent groups using medians
#' (or other location estimators) with percentile bootstrap inference.
#' Controls familywise error rate (FWER) using adjustable methods.
#'
#' @inheritParams common_params
#' @param x Data in matrix or list mode. If a matrix, columns represent groups.
#'   If a list, each element is a group.
#' @param est Measure of location to use (default: median). Can be any function
#'   that computes a location estimate.
#' @param ... Optional arguments passed to the estimator function \code{est}.
#' @param con A J × d matrix of contrast coefficients. If con=0 (default),
#'   all pairwise comparisons are performed.
#' @param bhop Logical; if TRUE, use Benjamini-Hochberg procedure for FDR control
#'   instead of FWER control (default: FALSE).
#' @param method Method for FWER control. Options: 'hoch' (Hochberg, default),
#'   'rom' (Rom), 'holm' (Holm), 'BY' (Benjamini-Yekutieli).
#'
#' @return A list with components:
#' \item{output}{Matrix with columns: Group(s), p-value, comparison estimate,
#'   CI lower bound, CI upper bound, and adjusted p-value.}
#' \item{con}{The contrast matrix used.}
#' \item{num.sig}{Number of significant comparisons at the alpha level.}
#'
#' @details
#' This function uses percentile bootstrap to perform multiple comparisons among
#' independent groups. For each comparison, bootstrap samples are drawn independently
#' from each group, and the difference in location estimates is computed.
#'
#' The p-value is computed by determining the proportion of bootstrap replicates
#' where the estimated difference is on the opposite side of zero compared to
#' the observed difference.
#'
#' FWER control methods:
#' \itemize{
#'   \item Hochberg (default): Step-up procedure, more powerful than Holm
#'   \item Rom: Improved version of Holm's procedure
#'   \item Holm: Step-down procedure
#'   \item BY: Benjamini-Yekutieli for dependent tests
#' }
#'
#' When \code{bhop=TRUE}, the Benjamini-Hochberg FDR procedure is used instead.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{msmed}}, \code{\link{med2g}}, \code{\link{linconb}}
#'
#' @examples
#' # Compare three groups
#' set.seed(123)
#' x1 <- rnorm(20, mean = 0)
#' x2 <- rnorm(20, mean = 0.5)
#' x3 <- rnorm(20, mean = 1)
#' data <- list(x1, x2, x3)
#'
#' # All pairwise comparisons with median
#' medpb(data, nboot = 1000)
#'
#' # Custom contrast with trimmed mean
#' con <- matrix(c(1, -0.5, -0.5), ncol = 1)
#' medpb(data, con = con, est = mean, trim = 0.2, nboot = 1000)
#'
#' @export
medpb<-function(x,alpha=.05,nboot=NA,grp=NA,est=median,con=0,bhop=FALSE,method='hoch',
SEED=TRUE,...){
#
#   Multiple comparisons for  J independent groups using medians.
#
#   A percentile bootstrap method.
#   FWE controlled via argument method
#   method =hoch  Hochberg;s method is used by default
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   est is the measure of location and defaults to the median
#   ... can be used to set optional arguments associated with est
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#
#   con can be used to specify linear contrasts; see the function lincon
#
#   Missing values are allowed.
#
con<-as.matrix(con)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in list mode or in matrix mode.')
if(!is.na(sum(grp))){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
x<-xx
}
J<-length(x)
tempn<-0
mvec<-NA
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
mvec[j]<-est(temp,...)
}
Jm<-J-1
#
# Determine contrast matrix
#
if(sum(con^2)==0){
ncon<-(J^2-J)/2
con<-matrix(0,J,ncon)
id<-0
for (j in 1:Jm){
jp<-j+1
for (k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
ncon<-ncol(con)
dvec<-alpha/c(1:ncon)
if(nrow(con)!=J)stop('Something is wrong with con; the number of rows does not match the number of groups.')
#  Determine nboot if a value was not specified
if(is.na(nboot)){
nboot<-5000
if(J <= 8)nboot<-4000
if(J <= 3)nboot<-2000
}
# Determine critical values
if(!bhop){
if(alpha==.05){
dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
}
}
if(bhop)dvec<-(ncon-c(1:ncon)+1)*alpha/ncon
bvec<-matrix(NA,nrow=J,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:J){
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,est,...) # Bootstrapped values for jth group
}
test<-NA
bcon<-t(con)%*%bvec #ncon by nboot matrix
tvec<-t(con)%*%mvec
for (d in 1:ncon){
tv<-sum(bcon[d,]==0)/nboot
test[d]<-sum(bcon[d,]>0)/nboot+.5*tv
if(test[d]> .5)test[d]<-1-test[d]
}
test<-2*test
output<-matrix(0,ncon,7)
dimnames(output)<-list(NULL,c('con.num','psihat','p.value','p.crit','ci.lower','ci.upper','adj.p.value'))
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,4]<-zvec
icl<-round(dvec[ncon]*nboot/2)+1
icu<-nboot-icl-1
for (ic in 1:ncol(con)){
output[ic,2]<-tvec[ic,]
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(bcon[ic,])
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
num.sig<-sum(output[,3]<=output[,4])
output[,7]=p.adjust(output[,3],method=method)

list(output=output,con=con,num.sig=num.sig)
}

# ==== med2g ====
#' Compare Medians of Two Independent Groups Using Bootstrap
#'
#' @description
#' Tests the hypothesis that two independent groups have equal medians using
#' a percentile bootstrap method. Provides a bootstrap p-value and confidence
#' interval for the difference in medians.
#'
#' @inheritParams common_params
#' @param x Numeric vector for first group.
#' @param y Numeric vector for second group.
#'
#' @return A list with components:
#' \item{p.value}{Bootstrap p-value for test of equal medians.}
#' \item{est.1}{Sample median of first group.}
#' \item{est.2}{Sample median of second group.}
#' \item{est.dif}{Difference in sample medians (group 1 - group 2).}
#' \item{ci.low}{Lower confidence limit for the difference.}
#' \item{ci.up}{Upper confidence limit for the difference.}
#'
#' @details
#' The function generates bootstrap samples from each group independently and
#' computes the median for each sample. The p-value is computed by determining
#' the proportion of bootstrap replicates where group 1 median exceeds group 2
#' median, adjusting for ties and converting to a two-sided test.
#'
#' The confidence interval is obtained from the empirical distribution of
#' bootstrap differences using the percentile method.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{msmed}}, \code{\link{medpb}}, \code{\link{msmedci}}
#'
#' @examples
#' # Compare medians of two groups
#' set.seed(123)
#' x <- rnorm(25, mean = 0)
#' y <- rnorm(30, mean = 0.5)
#' med2g(x, y, nboot = 1000)
#'
#' @export
med2g<-function(x,y,alpha=.05,nboot=2000,SEED=TRUE,...){
#
#   Compare medians of two independent groups using percentile bootstrap
#
#   Missing values are allowed.
#
x<-elimna(x)
y<-elimna(y)
mvec<-NA
mvec[1]<-median(x)
mvec[2]<-median(y)
bvec<-NA
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
datay<-matrix(sample(y,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec1<-apply(datax,1,median) # Bootstrapped values for jth group
bvec2<-apply(datay,1,median) # Bootstrapped values for jth group
test<-sum((bvec1>bvec2))/nboot
tv<-sum(bvec1==bvec2)/nboot
test<-test+.5*tv
if(test> .5)test<-1-test
test<-2*test
dvec<-sort(bvec1-bvec2)
icl<-round(alpha*nboot/2)+1
icu<-nboot-icl-1
cilow<-dvec[icl]
ciup<-dvec[icu]
list(p.value=test,est.1=mvec[1],est.2=mvec[2],est.dif=mvec[1]-mvec[2],ci.low=cilow,ci.up=ciup)
}

# ==== msmedci ====
#' Confidence Interval for the Median Using McKean-Schrader Method
#'
#' @description
#' Computes a confidence interval for the population median using the
#' McKean-Schrader standard error. Also provides a hypothesis test for
#' whether the median equals a specified null value.
#'
#' @inheritParams common_params
#' @param x Numeric vector of data.
#' @param nullval Null hypothesis value for the median (default: 0).
#'
#' @return A list with components:
#' \item{test}{Test statistic for H0: median = nullval.}
#' \item{ci.low}{Lower confidence limit for the median.}
#' \item{ci.hi}{Upper confidence limit for the median.}
#' \item{p.value}{Two-sided p-value for the test.}
#' \item{median}{Sample median.}
#'
#' @details
#' The function uses the McKean-Schrader method to estimate the standard error
#' of the sample median. The test statistic is computed as (median - nullval) / SE,
#' and p-values are based on the standard normal distribution.
#'
#' Confidence intervals are constructed as median ± critical value × SE, where
#' the critical value is from the standard normal distribution.
#'
#' This approach is appropriate when the sample median is approximately normally
#' distributed, which holds for reasonably large samples from continuous distributions.
#'
#' @references
#' McKean, J. W., & Schrader, R. M. (1984). A comparison of methods for
#' studentizing the sample median. Communications in Statistics-Simulation
#' and Computation, 13(6), 751-773.
#'
#' @seealso \code{\link{medcipb}}, \code{\link{msmed}}
#'
#' @examples
#' # CI for median
#' set.seed(123)
#' x <- rnorm(50, mean = 5, sd = 2)
#' msmedci(x)
#'
#' # Test if median equals 5
#' msmedci(x, nullval = 5)
#'
#' @export
msmedci<-function(x,alpha=.05,nullval=0){
#
# Confidence interval for the median
#
se<-msmedse(x)
est<-median(x)
ci.low<-est-qnorm(1-alpha/2)*se
ci.hi<-est+qnorm(1-alpha/2)*se
test<-(est-nullval)/se
p.value<-2*(1-pnorm(abs(test)))
list(test=test,ci.low=ci.low,ci.hi=ci.hi,p.value=p.value,median=est)
}

# ==== medcipb ====
#' Bootstrap Confidence Interval for the Median
#'
#' @description
#' Computes a percentile bootstrap confidence interval for the median of a
#' single variable. Optionally performs a hypothesis test if a null value
#' is specified.
#'
#' @inheritParams common_params
#' @param x Numeric vector of data.
#' @param null.val Optional null hypothesis value for the median. If specified,
#'   a bootstrap p-value is computed (default: NA).
#'
#' @return A list with components:
#' \item{Est.}{Sample median.}
#' \item{ci.low}{Lower confidence limit.}
#' \item{ci.up}{Upper confidence limit.}
#' \item{p.value}{Bootstrap p-value (only if null.val is specified).}
#'
#' @details
#' The function generates bootstrap samples by resampling with replacement
#' from the original data. For each bootstrap sample, the median is computed.
#'
#' The confidence interval is obtained using the percentile method: the
#' alpha/2 and 1-alpha/2 quantiles of the bootstrap distribution.
#'
#' If a null value is specified, the p-value is computed by determining
#' the proportion of bootstrap medians greater than the null value, adjusting
#' for ties, and converting to a two-sided test.
#'
#' Missing values are automatically removed before analysis.
#'
#' For the Harrell-Davis estimator version, see \code{\link{hdpb}}.
#'
#' @seealso \code{\link{msmedci}}, \code{\link{hdpb}}, \code{\link{onesampb}}
#'
#' @examples
#' # Bootstrap CI for median
#' set.seed(123)
#' x <- rnorm(30, mean = 5, sd = 2)
#' medcipb(x, nboot = 1000)
#'
#' # Test if median equals 5
#' medcipb(x, null.val = 5, nboot = 1000)
#'
#' @export
medcipb<-function(x,alpha=.05,null.val=NA,nboot=500,SEED=TRUE,...){
#
#   Bootstrap confidence interval for the median of single variable.
#   The usual sample median is used. hdpb uses the Harrell--Davis estimator
#   Missing values are allowed.
#
x<-elimna(x)
est=median(x)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,median) # Bootstrapped values
test<-NULL
if(!is.na(null.val)){
tv<-sum(bvec==null.val)/nboot
test<-sum(bvec>null.val)/nboot+.5*tv
if(test> .5)test<-1-test
test<-2*test
}
bvec<-sort(bvec)
icl<-round(alpha*nboot/2)+1
icu<-nboot-icl-1
cilow<-bvec[icl]
ciup<-bvec[icu]
list(Est.=est,ci.low=cilow,ci.up=ciup,p.value=test)
}

# ==== medr ====
#' Robust Median Comparisons Using Depth-Based Methods
#'
#' @description
#' Tests the hypothesis that pairwise differences among groups have location
#' equal to zero using a depth-based multivariate approach. Provides a global
#' test that accounts for the dependence among all pairwise comparisons.
#'
#' @inheritParams common_params
#' @param x Data in matrix or list mode. If a matrix, columns represent groups.
#'   If a list, each element is a group.
#' @param est Measure of location to use (default: median).
#' @param ... Optional arguments passed to the estimator function \code{est}.
#' @param op Method for computing depth (default: 1). Options:
#'   \itemize{
#'     \item 1: Mahalanobis depth
#'     \item 2: Mahalanobis depth based on MCD covariance matrix
#'     \item 3: Projection depth
#'     \item 4: Projection depth using FORTRAN implementation
#'   }
#' @param MM Logical; for projection depth, use modified algorithm (see \code{\link{pdis}}).
#' @param cop Integer; for projection depth, method for computing center (see \code{\link{pdis}}).
#' @param pr Logical; if TRUE (default), print progress messages.
#'
#' @return A list with components:
#' \item{p.value}{Overall p-value for the global test.}
#' \item{output}{Matrix showing group pairs and estimated pairwise differences.}
#'
#' @details
#' This function treats the vector of all pairwise differences as a multivariate
#' observation and uses data depth to construct a global test. The null hypothesis
#' is that all pairwise differences have median (or other location measure) equal
#' to zero.
#'
#' \strong{Procedure}:
#' \enumerate{
#'   \item Compute all pairwise differences in location estimates
#'   \item Generate bootstrap samples and compute pairwise differences
#'   \item Use depth-based methods to compare observed differences to bootstrap null distribution
#' }
#'
#' The test accounts for the dependence among pairwise comparisons by treating
#' them jointly in a multivariate space. Depth is measured relative to the
#' bootstrap distribution centered at the observed differences.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{medpb}}, \code{\link{pdis}}, \code{\link{mahalanobis}}
#'
#' @examples
#' # Global test for median differences
#' set.seed(123)
#' x1 <- rnorm(20, mean = 0)
#' x2 <- rnorm(20, mean = 0)
#' x3 <- rnorm(20, mean = 0)
#' data <- list(x1, x2, x3)
#'
#' # Using Mahalanobis depth
#' medr(data, nboot = 500)
#'
#' @export
medr<-function(x,est=median,alpha=.05,nboot=500,grp=NA,op=1,MM=FALSE,cop=3,pr=TRUE,
SEED=TRUE,...){
#
#   Test the hypothesis that the distribution for each pairwise
#   difference has a measure of location = 0
#   By default, the  median is used
#
#   The default number of bootstrap samples is nboot=500
#
#   op controls how depth is measured
#   op=1, Mahalanobis
#   op=2, Mahalanobis based on MCD covariance matrix
#   op=3, Projection distance
#   op=4, Projection distance using FORTRAN version
#
#   for arguments MM and cop, see pdis.
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x)){
xx<-list()
for(i in 1:ncol(x)){
xx[[i]]<-x[,i]
}
x<-xx
}
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
if(!is.na(grp)){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[1]]]
x<-xx
}
J<-length(x)
mvec<-NA
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
mvec[j]<-est(temp,...)
}
Jm<-J-1
d<-(J^2-J)/2
data<-list()
bvec<-matrix(NA,ncol=d,nrow=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print("Taking bootstrap samples. Please wait.")
for(it in 1:nboot){
for(j in 1:J)data[[j]]<-sample(x[[j]],size=length(x[[j]]),replace=TRUE)
dval<-0
for(j in 1:J){
for(k in 1:J){
if(j<k){
dval<-dval+1
bvec[it,dval]<-loc2dif(data[[j]],data[[k]],est=est,...)
}}}}
output<-matrix(NA,nrow=d,ncol=3)
dimnames(output)<-list(NULL,c("Group","Group","psihat"))
tvec<-NA
dval<-0
for(j in 1:J){
for(k in 1:J){
if(j<k){
dval<-dval+1
output[dval,1]<-j
output[dval,2]<-k
tvec[dval]<-loc2dif(x[[j]],x[[k]],est=est,...)
output[dval,3]<-tvec[dval]
}}}
tempcen<-apply(bvec,1,mean)
vecz<-rep(0,d)
smat<-var(bvec-tempcen+tvec)
temp<-bvec-tempcen+tvec
bcon<-rbind(bvec,vecz)
if(op==1)dv<-mahalanobis(bcon,tvec,smat)
if(op==2){
smat<-cov.mcd(temp)$cov
dv<-mahalanobis(bcon,tvec,smat)
}
if(op==3){
print("Computing p-value. Might take a while with op=3")
dv<-pdis(bcon,MM=MM,cop=cop,center=tvec)
}
if(op==4)dv<-pdis.for(bcon,MM=MM,cop=cop,pr=FALSE,center=tvec)
bplus<-nboot+1
sig.level<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
if(op==4)print(sig.level)

list(p.value=sig.level,output=output)
}

# ==== medind ====
#' Independence Test Based on Quantile Regression for Median
#'
#' @description
#' Tests whether a regression relationship exists between predictors and a
#' response by testing if the regression surface is a flat horizontal plane.
#' Uses a quantile regression approach with the median (or other quantiles).
#'
#' @param x Matrix or data frame of predictor variables (can have 1-8 columns).
#' @param y Numeric vector of response variable.
#' @param qval Quantile to use (default: 0.5 for median). Can also use 0.25 or 0.75.
#' @inheritParams common_params
#' @param nboot Number of simulations for computing p-value when com.pval=TRUE
#'   (default: 1000).
#' @param com.pval Logical; if TRUE, compute p-value via simulation. If FALSE
#'   (default), use approximate critical values from stored tables.
#' @param pr Logical; if TRUE (default), print informational messages.
#' @param xout Logical; if TRUE, remove outliers from predictors before testing
#'   (default: FALSE).
#' @param outfun Function for detecting outliers when xout=TRUE (default: \code{\link{out}}).
#' @param ... Additional arguments passed to \code{outfun}.
#'
#' @return A list with components:
#' \item{test.stat}{The test statistic value.}
#' \item{crit.value}{Critical value for the test.}
#' \item{p.value}{P-value (only if com.pval=TRUE).}
#' \item{Decision}{Character string: "Reject" or "Fail To Reject".}
#'
#' @details
#' This function tests the null hypothesis that the predictors have no
#' relationship with the response at the specified quantile level. It is
#' based on a modification of the method by He and Zhu (2003).
#'
#' \strong{Null hypothesis}: The qval-th quantile of Y is constant regardless
#' of the values of X (i.e., the quantile regression surface is flat).
#'
#' Approximate critical values are available for:
#' \itemize{
#'   \item Sample sizes: 10 ≤ n ≤ 400
#'   \item Number of predictors: p = 1, ..., 8
#'   \item Quantiles: 0.25, 0.5, 0.75
#'   \item Alpha levels: 0.01, 0.025, 0.05, 0.10
#' }
#'
#' If data characteristics fall outside these ranges, set com.pval=TRUE to
#' compute a p-value via simulation (which increases computation time).
#'
#' Missing values are automatically removed before analysis.
#'
#' @references
#' He, X., & Zhu, L. X. (2003). A lack-of-fit test for quantile regression.
#' Journal of the American Statistical Association, 98(464), 1013-1022.
#'
#' @seealso \code{\link{qreg}}, \code{\link{medcipb}}
#'
#' @examples
#' # Test for median regression relationship
#' set.seed(123)
#' n <- 50
#' x <- rnorm(n)
#' y <- 2 * x + rnorm(n)  # Linear relationship
#'
#' # Should reject (relationship exists)
#' medind(x, y, qval = 0.5)
#'
#' # No relationship
#' y2 <- rnorm(n)
#' medind(x, y2, qval = 0.5)
#'
#' @export
medind<-function(x,y,qval=.5,nboot=1000,com.pval=FALSE,SEED=TRUE,alpha=.05,pr=TRUE,
xout=FALSE,outfun=out,...){
#
# Test the hypothesis that the regression surface is a flat
# horizontal plane.
# The method is based on a modification of a method derived by
#  He and Zhu 2003, JASA, 98, 1013-1022.
# Here, resampling is avoided using approximate critical values if
# com.pval=F
#
#  critical values are available for 10<=n<=400, p=1,...,8 and
#  quantiles
#  qval=.25,.5, .75.
#
#  To get a p-value, via simulations, set  com.pval=T
#  nboot is number of simulations used to determine the p-value.
#
if(pr){
if(!com.pval)print("To get a p-value, set com.pval=T")
print("Reject if the test statistic exceeds the critical value")
if(length(y)>400)print("If n>400, current version requires com.pval=TRUE, resulting in high execution time")
}
#store.it=F
x<-as.matrix(x)
p<-ncol(x)
pp1<-p+1
p.val<-NULL
crit.val<-NULL
yx<-elimna(cbind(y,x)) #Eliminate missing values.
y<-yx[,1]
x<-yx[,2:pp1]
x<-as.matrix(x)
if(xout){
flag<-outfun(x,...)$keep
x<-x[flag,]
y<-y[flag]
}
n<-length(y)
if(n>400)com.pval=T
if(qval==.5){
resmat1=matrix(c( 0.0339384580, 0.044080032, 0.050923441, 0.064172557,
 0.0153224731, 0.021007108, 0.027687963, 0.032785044,
 0.0106482053, 0.014777728, 0.018249546, 0.023638611,
 0.0066190573, 0.009078091, 0.011690825, 0.014543009,
 0.0031558563, 0.004374515, 0.005519069, 0.007212951,
 0.0015448987, 0.002231473, 0.002748314, 0.003725916,
 0.0007724197, 0.001021767, 0.001370776, 0.001818037),ncol=4,nrow=7,byrow=TRUE)
resmat2=matrix(c(
 0.052847794, 0.061918744, 0.071346969, 0.079163419,
 0.021103277, 0.027198076, 0.031926052, 0.035083610,
 0.013720585, 0.018454145, 0.022177381, 0.026051716,
 0.008389969, 0.010590374, 0.012169233, 0.015346065,
 0.004261627, 0.005514060, 0.007132021, 0.008416836,
 0.001894753, 0.002416311, 0.003085230, 0.003924706,
 0.001045346, 0.001347837, 0.001579373, 0.001864344),ncol=4,nrow=7,byrow=TRUE)
resmat3=matrix(c(
0.071555715, 0.082937665, 0.089554679, 0.097538044,
0.031060795, 0.035798539, 0.043862556, 0.053712151,
0.019503635, 0.023776479, 0.027180121, 0.030991367,
0.011030001, 0.013419347, 0.015557409, 0.017979524,
0.005634478, 0.006804788, 0.007878358, 0.008807657,
0.002552182, 0.003603778, 0.004275965, 0.005021989,
0.001251044, 0.001531919, 0.001800608, 0.002037870),ncol=4,nrow=7,byrow=TRUE)
resmat4=matrix(c(
0.093267532, 0.101584002, 0.108733965, 0.118340448,
0.038677863, 0.045519806, 0.051402903, 0.060097046,
0.024205231, 0.029360145, 0.034267265, 0.039381482,
0.013739157, 0.015856343, 0.018065898, 0.019956084,
0.006467562, 0.007781030, 0.009037972, 0.010127143,
0.003197162, 0.003933525, 0.004656625, 0.005929469,
0.001652690, 0.001926060, 0.002363874, 0.002657071),ncol=4,nrow=7,byrow=TRUE)
resmat5=matrix(c(
0.117216934, 0.124714114, 0.129458602, 0.136456163,
0.048838630, 0.055608712, 0.060580045, 0.067943676,
0.030594644, 0.035003872, 0.040433885, 0.047648696,
0.016940240, 0.019527491, 0.022047442, 0.025313443,
0.008053039, 0.009778574, 0.011490394, 0.013383628,
0.003760567, 0.004376294, 0.005097890, 0.005866240,
0.001894616, 0.002253522, 0.002612405, 0.002938808),ncol=4,nrow=7,byrow=TRUE)
resmat6=matrix(c(
0.136961531, 0.144120225, 0.149003907, 0.152667432,
0.055909481, 0.062627211, 0.069978086, 0.081189957,
0.034634825, 0.040740587, 0.044161376, 0.047722045,
0.020165417, 0.023074738, 0.025881208, 0.028479913,
0.009436297, 0.011246968, 0.013220963, 0.015100546,
0.004644596, 0.005334418, 0.006040595, 0.007237195,
0.002277590, 0.002635712, 0.002997398, 0.003669488),ncol=4,nrow=7,byrow=TRUE)
resmat7=matrix(c(
 0.156184672, 0.163226643, 0.171754686, 0.177142753,
 0.070117003, 0.077052773, 0.082728047, 0.090410797,
 0.041774517, 0.047379662, 0.053101833, 0.057674454,
 0.023384451, 0.026014421, 0.029609042, 0.032619018,
 0.010856382, 0.012567043, 0.013747870, 0.016257014,
 0.005164004, 0.006131755, 0.006868101, 0.008351046,
 0.002537642, 0.003044154, 0.003623654, 0.003974469),ncol=4,nrow=7,byrow=TRUE)
resmat8=matrix(c(
0.178399742, 0.180006714, 0.193799396, 0.199585892,
0.078032767, 0.085624186, 0.091511226, 0.102491785,
0.045997886, 0.052181615, 0.057362163, 0.062630424,
0.025895739, 0.029733034, 0.033764463, 0.037873655,
0.012195876, 0.013663248, 0.015487587, 0.017717864,
0.005892418, 0.006876488, 0.007893475, 0.008520783,
0.002839731, 0.003243909, 0.003738571, 0.004124057),ncol=4,nrow=7,byrow=TRUE)
crit5=array(cbind(resmat1,resmat2,resmat3,resmat4,resmat5,resmat6,resmat7,
resmat8),c(7,4,8))
flag=TRUE
crit.val=NULL
if(p > 8)flag=FALSE
if(n<10 || n>=400)flag=FALSE
aval<-c(.1,.05,.025,.01)
aokay<-duplicated(c(alpha,aval))
if(sum(aokay)==0)flag=FALSE
if(flag){
nalpha=c(0:4)
asel=c(0,aval)
ialpha=nalpha[aokay]
critit=crit5[,ialpha,p]
nvec<-c(10,20,30,50,100,200,400)
nval<-duplicated(c(n,nvec))
nval<-nval[2:8]
if(sum(nval)>0)crit.val<-critit[nval]
loc<-rank(c(n,nvec))
xx<-c(1/nvec[loc[1]-1]^1.5,1/nvec[loc[1]]^1.5)
yy<-c(critit[loc[1]-1],critit[loc[1]])
icoef<-tsp1reg(xx,yy)$coef
crit.val<-icoef[1]+icoef[2]/n^1.5
}}
mqval<-min(c(qval,1-qval))
if(mqval==.25){
resmat1=matrix(c(
 0.029933486, 0.0395983678, 0.054087714, 0.062961453,
 0.011122294, 0.0149893431, 0.018154062, 0.022685244,
 0.009207200, 0.0113020766, 0.014872309, 0.019930730,
 0.004824185, 0.0070402246, 0.010356886, 0.013176896,
 0.002370379, 0.0033146605, 0.004428004, 0.005122988,
 0.001106460, 0.0016110185, 0.001984450, 0.002650256,
 0.000516646, 0.0006796144, 0.000868751, 0.001202042),ncol=4,nrow=7,byrow=TRUE)
resmat2=matrix(c(
0.0448417783, 0.0602598211, 0.066001091, 0.087040667,
0.0173410522, 0.0224713157, 0.027370822, 0.033435727,
0.0121205549, 0.0150409465, 0.018938516, 0.022643559,
0.0064894201, 0.0084611518, 0.010700320, 0.013232000,
0.0029734778, 0.0040641310, 0.004911086, 0.005769038,
0.0015149104, 0.0020584993, 0.002582982, 0.003114029,
0.0007984207, 0.0009929547, 0.001182739, 0.001398774),ncol=4,nrow=7,byrow=TRUE)
resmat3=matrix(c(
0.0636530860, 0.072974943, 0.083840562, 0.097222407,
0.0216586978, 0.027436566, 0.031875356, 0.036830302,
0.0152898678, 0.018964066, 0.021728817, 0.028959751,
0.0083568493, 0.010071525, 0.012712862, 0.015254576,
0.0039033578, 0.004764140, 0.005577071, 0.006660322,
0.0019139215, 0.002343152, 0.002833612, 0.003465269,
0.0009598105, 0.001146689, 0.001355930, 0.001547572),ncol=4,nrow=7,byrow=TRUE)
resmat4=matrix(c(
 0.085071252, 0.095947936, 0.104197413, 0.118449765,
 0.029503024, 0.034198704, 0.039543410, 0.045043759,
 0.019203266, 0.022768842, 0.026886843, 0.033481535,
 0.011440493, 0.013555017, 0.016138970, 0.018297815,
 0.004863139, 0.005756305, 0.007385239, 0.009114958,
 0.002635144, 0.003111160, 0.003769051, 0.004215897,
 0.001188837, 0.001435179, 0.001727871, 0.001956372),ncol=4,nrow=7,byrow=TRUE)
resmat5=matrix(c(
0.102893512, 0.114258558, 0.122545016, 0.130222265,
0.036733497, 0.042504996, 0.048663576, 0.055456582,
0.024192946, 0.028805967, 0.032924489, 0.038209545,
0.012663224, 0.014635216, 0.017275594, 0.019736410,
0.006105572, 0.007310803, 0.008960242, 0.009745320,
0.003067163, 0.003614637, 0.003997615, 0.004812373,
0.001441008, 0.001732819, 0.002078651, 0.002307551),ncol=4,nrow=7,byrow=TRUE)
resmat6=matrix(c(
0.117642769, 0.126566104, 0.133106804, 0.142280074,
0.044309420, 0.049731991, 0.053912739, 0.060512997,
0.028607224, 0.033826020, 0.038616476, 0.043546500,
0.015445120, 0.017557181, 0.020040720, 0.022747707,
0.007334749, 0.008406468, 0.009392098, 0.010919651,
0.003352200, 0.003814582, 0.004380562, 0.005252154,
0.001703698, 0.002001713, 0.002338651, 0.002772864),ncol=4,nrow=7,byrow=TRUE)
resmat7=matrix(c(
0.106573121, 0.113058950, 0.117388191, 0.121286795,
0.052170054, 0.058363322, 0.064733684, 0.069749344,
0.030696897, 0.035506926, 0.039265698, 0.044437674,
0.016737307, 0.019605734, 0.021253610, 0.022922988,
0.007767232, 0.009231789, 0.010340874, 0.011471110,
0.003998261, 0.004590177, 0.005506926, 0.006217415,
0.001903372, 0.002174748, 0.002519055, 0.002858655),ncol=4,nrow=7,byrow=TRUE)
resmat8=matrix(c(
 0.119571179, 0.126977461, 0.130120853, 0.133258294,
 0.059499563, 0.067185338, 0.071283297, 0.079430577,
 0.034310968, 0.039827130, 0.044451690, 0.048512464,
 0.018599530, 0.021093909, 0.023273085, 0.027471116,
 0.009135712, 0.010901687, 0.012288682, 0.013729545,
 0.004382249, 0.005191810, 0.005598429, 0.006484433,
 0.002196973, 0.002525918, 0.002818550, 0.003242426),ncol=4,nrow=7,byrow=TRUE)
crit5=array(cbind(resmat1,resmat2,resmat3,resmat4,resmat5,resmat6,resmat7,
resmat8),c(7,4,8))
flag=TRUE
crit.val=NULL
if(p > 8)flag=FALSE
if(n<10 || n>=400)flag=FALSE
aval<-c(.1,.05,.025,.01)
aokay<-duplicated(c(alpha,aval))
if(sum(aokay)==0)flag=FALSE
if(flag){
nalpha=c(0:4)
asel=c(0,aval)
ialpha=nalpha[aokay]
critit=crit5[,ialpha,p]
nvec<-c(10,20,30,50,100,200,400)
nval<-duplicated(c(n,nvec))
nval<-nval[2:8]
if(sum(nval)>0)crit.val<-critit[nval,p]
loc<-rank(c(n,nvec))
xx<-c(1/nvec[loc[1]-1]^1.5,1/nvec[loc[1]]^1.5)
yy<-c(critit[loc[1]-1],critit[loc[1]])
icoef<-tsp1reg(xx,yy)$coef
crit.val<-icoef[1]+icoef[2]/n^1.5
}}
if(is.null(crit.val))com.pval=TRUE
# no critical value found, so a p-value will be computed
# the code for checking the file medind.crit, which appears
# next, is not working yet.
gdot<-cbind(rep(1,n),x)
gdot<-ortho(gdot)
x<-gdot[,2:pp1]
x<-as.matrix(x)
coef<-NULL
if(qval==.5)coef<-median(y)
if(qval==.25)coef<-idealf(y)$ql
if(qval==.75)coef<-idealf(y)$qu
if(is.null(coef))coef<-qest(y,q=qval)
res<-y-coef
psi<-NA
psi<-ifelse(res>0,qval,qval-1)
rnmat<-matrix(0,nrow=n,ncol=pp1)
ran.mat<-apply(x,2,rank)
flagvec<-apply(ran.mat,1,max)
for(j in 1:n){
flag<-ifelse(flagvec<=flagvec[j],TRUE,FALSE)
flag<-as.numeric(flag)
rnmat[j,]<-apply(flag*psi*gdot,2,sum)
}
rnmat<-rnmat/sqrt(n)
temp<-matrix(0,pp1,pp1)
for(i in 1:n)temp<-temp+rnmat[i,]%*%t(rnmat[i,])
temp<-temp/n
test<-max(eigen(temp)$values)
if(com.pval){
if(SEED)set.seed(2)
p.val<-0
rem<-0
for(i in 1:nboot){
yboot<-rnorm(n)
if(p==1)xboot<-rnorm(n)
if(p>1)xboot<-rmul(n,p=p)
temp3<-medindsub(x,yboot,qval=qval)
if(test>=temp3)p.val<-p.val+1
rem[i]<-temp3
}
ic10<-round(.9*nboot)
ic05<-round(.95*nboot)
ic025<-round(.975*nboot)
ic001<-round(.99*nboot)
rem<-sort(rem)
p.val<-1-p.val/nboot
# now remember the critical values by storing them in "medind.crit"
#if(store.it)
#write(c(n,p,qval,rem[ic10],rem[ic05],rem[ic025],rem[ic001]),"medind.crit",
#append=T,ncolumns=7)
print("The .1, .05, .025 and .001 critical values are:")
print(c(rem[ic10],rem[ic05],rem[ic025],rem[ic001]))
crit.val<-rem[ic05]
}
names(crit.val)=""
Decision="Fail To Reject"
if(test>=crit.val)Decision="Reject"
list(test.stat=test,crit.value=crit.val,p.value=p.val,Decision=Decision)
}

# ==== medindsub ====
#' Helper Function for Median Independence Test (Internal)
#'
#' @description
#' Internal helper function for \code{\link{medind}} that computes the test
#' statistic for a single dataset. Used during simulation. Not intended for
#' direct use.
#'
#' @param x Matrix of predictor variables.
#' @param y Numeric vector of response variable.
#' @param qval Quantile level (default: 0.5 for median).
#'
#' @return Numeric test statistic value.
#'
#' @details
#' Computes the test statistic for the quantile regression independence test
#' without storing critical values or computing p-values. Called repeatedly
#' during simulation in \code{\link{medind}}.
#'
#' @seealso \code{\link{medind}}
#'
#' @keywords internal
medindsub<-function(x,y,qval=.5){
#
x<-as.matrix(x)
n<-length(y)
p<-ncol(x)
pp1<-p+1
tvec<-c(qval,0-qval,1-qval,qval-1)
pval<-c((1-qval)/2,(1-qval)/2,qval/2,qval/2)
gdot<-cbind(rep(1,n),x)
gdot<-ortho(gdot)
x<-gdot[,2:pp1]
x<-as.matrix(x)
if(qval==.5)coef<-median(y)
if(qval!=.5)coef<-qest(y)
res<-y-coef
psi<-NA
psi<-ifelse(res>0,qval,qval-1)
rnmat<-matrix(0,nrow=n,ncol=pp1)
ran.mat<-apply(x,2,rank)
flagvec<-apply(ran.mat,1,max)
for(j in 1:n){
flag<-ifelse(flagvec>=flagvec[j],T,F)
rnmat[j,]<-apply(flag*psi*gdot,2,sum)
}
rnmat<-rnmat/sqrt(n)
temp<-matrix(0,pp1,pp1)
for(i in 1:n)temp<-temp+rnmat[i,]%*%t(rnmat[i,])
temp<-temp/n
test<-max(eigen(temp)$values)
test
}

# ==== msmedse ====
#' McKean-Schrader Standard Error for the Median
#'
#' @description
#' Computes the standard error of the sample median using the method
#' recommended by McKean and Schrader (1984).
#'
#' @param x Numeric vector of data.
#'
#' @return Numeric standard error estimate for the median.
#'
#' @details
#' The McKean-Schrader method estimates the standard error by computing
#' an empirical confidence interval for the median and converting it to
#' a standard error. Specifically, it uses the 99% confidence interval
#' bounds and divides by the appropriate normal critical value.
#'
#' A warning is issued if tied values are detected, as ties can make the
#' standard error estimate highly inaccurate, even with large samples.
#' In such cases, consider using bootstrap methods instead.
#'
#' Missing values are automatically removed before computation.
#'
#' @references
#' McKean, J. W., & Schrader, R. M. (1984). A comparison of methods for
#' studentizing the sample median. Communications in Statistics-Simulation
#' and Computation, 13(6), 751-773.
#'
#' @seealso \code{\link{msmedci}}, \code{\link{msmed}}
#'
#' @examples
#' # Standard error for median
#' set.seed(123)
#' x <- rnorm(50, mean = 5, sd = 2)
#' msmedse(x)
#'
#' @keywords internal
#' @export
msmedse<-function(x){
#
# Compute  standard error of the median using method
# recommended by McKean and Shrader (1984).
#
x=elimna(x)
chk=sum(duplicated(x))
if(chk>0){
print("WARNING: tied values detected.")
print("Estimate of standard error might be highly inaccurate, even with n large")
}
y<-sort(x)
n<-length(x)
av<-round((n+1)/2-qnorm(.995)*sqrt(n/4))
if(av==0)av<-1
top<-n-av+1
sqse<-((y[top]-y[av])/(2*qnorm(.995)))^2
sqse<-sqrt(sqse)
sqse
}

# ==== med2mcp ====
#' Multiple Comparisons for Two-Way ANOVA Using Medians
#'
#' @description
#' Performs multiple comparisons for a J × K factorial design using medians
#' with percentile bootstrap methods. Tests main effects for factors A and B,
#' as well as interaction contrasts.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in list mode. Length should equal J×K, with x[[1]] containing
#'   data for level (1,1), x[[2]] for level (1,2), etc.
#' @param grp Vector specifying the order of groups if data are not in the
#'   default order.
#' @param p Total number of groups (should equal J×K). Used for validation.
#' @inheritParams common_params
#' @param bhop Logical; if TRUE, use Benjamini-Hochberg FDR control instead of
#'   FWER control (default: FALSE).
#' @param pr Logical; if TRUE (default), print progress messages.
#'
#' @return A list with components:
#' \item{Factor.A}{Results for Factor A main effect comparisons (from \code{medpb}).}
#' \item{Factor.B}{Results for Factor B main effect comparisons (from \code{medpb}).}
#' \item{Factor.AB}{Results for A × B interaction contrasts (from \code{medpb}).}
#' \item{bhop}{Value of bhop parameter used.}
#' \item{SEED}{Value of SEED parameter used.}
#'
#' @details
#' This function performs bootstrap multiple comparisons for two-way factorial
#' designs with independent groups. It generates appropriate contrast matrices
#' for Factor A main effects, Factor B main effects, and A × B interactions
#' using \code{\link{con2way}}, then applies \code{\link{medpb}} to perform
#' the comparisons.
#'
#' For each set of contrasts, percentile bootstrap confidence intervals and
#' p-values are computed. FWER control is applied separately to each family
#' of tests (Factor A, Factor B, and interaction).
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{med2way}}, \code{\link{medpb}}, \code{\link{con2way}}
#'
#' @examples
#' # Two-way multiple comparisons
#' set.seed(123)
#' # 2x3 design
#' data <- list(
#'   rnorm(15, mean = 0),    # A1,B1
#'   rnorm(15, mean = 0.5),  # A1,B2
#'   rnorm(15, mean = 0.3),  # A1,B3
#'   rnorm(15, mean = 1),    # A2,B1
#'   rnorm(15, mean = 1.8),  # A2,B2  (interaction effect)
#'   rnorm(15, mean = 1.2)   # A2,B3
#' )
#'
#' result <- med2mcp(J = 2, K = 3, x = data, nboot = 1000)
#'
#' # View Factor A comparisons
#' result$Factor.A
#'
#' # View interaction comparisons
#' result$Factor.AB
#'
#' @export
med2mcp<-function(J,K,x,grp=c(1:p),p=J*K,nboot=NA,alpha=.05,SEED=TRUE,pr=TRUE,
bhop=FALSE){
#
#  Perform multiple comparisons for  J by K anova using medians with
#   using a percentile bootstrap method
#
#
#  The R variable data is assumed to contain the raw
#  data stored in a matrix or in list mode.
#  If stored in list mode, data[[1]] contains the data
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
Factor.A<-medpb(x,con=conA,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
Factor.B<-medpb(x,con=conB,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
Factor.AB<-medpb(x,con=conAB,alpha=alpha,nboot=nboot,bhop=bhop,SEED=FALSE)
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.AB=Factor.AB,bhop=bhop,SEED=FALSE)
}

# ==== dmedpb ====
#' Bootstrap Comparisons for Medians of Dependent Groups
#'
#' @description
#' Performs multiple comparisons for medians (or other location measures) of
#' dependent (repeated measures) groups using percentile bootstrap. Can analyze
#' either difference scores or marginal distributions.
#'
#' @inheritParams common_params
#' @param x Data in matrix or list mode. If a matrix, rows are subjects and
#'   columns are repeated measures. For two groups only, can be a vector with
#'   y specified.
#' @param y Optional second group for two-group comparison (default: NULL).
#' @param est Measure of location to use (default: median).
#' @param ... Optional arguments passed to the estimator function \code{est}.
#' @param con A J × d matrix of contrast coefficients. If con=0 (default),
#'   all pairwise comparisons are performed.
#' @param dif Logical; if TRUE, analyze difference scores. If FALSE (default),
#'   analyze marginal distributions.
#' @param hoch Logical; if TRUE (default), use Hochberg's method for FWER control.
#' @param plotit Logical; if TRUE (default), create plots of results.
#' @param xlab Label for x-axis in plots (default: "Group 1").
#' @param ylab Label for y-axis in plots (default: "Group 2").
#' @param ylab.ebar Label for y-axis in error bar plot.
#' @param pr Logical; if TRUE (default), print progress messages.
#' @param BA Logical; related to Bland-Altman plot option.
#' @param PCI Logical; if TRUE and dif=TRUE with est=median, plot confidence
#'   intervals for difference scores (default: FALSE).
#' @param EBAR Alias for PCI (default: PCI).
#'
#' @return A list with components:
#' \item{output}{Matrix with columns: comparison info, p-value, estimate,
#'   CI lower bound, CI upper bound, critical p-value.}
#' \item{con}{The contrast matrix used.}
#' \item{num.sig}{Number of significant comparisons.}
#'
#' @details
#' This function extends \code{\link{medpb}} to repeated measures designs.
#' It wraps \code{\link{rmmcppb}} with defaults set for median comparisons
#' and handles tied values appropriately.
#'
#' \strong{Analysis modes}:
#' \itemize{
#'   \item \code{dif=TRUE}: Analyzes difference scores using \code{rmmcppbd}.
#'     Appropriate when interest is in within-subject changes.
#'   \item \code{dif=FALSE}: Analyzes marginal distributions. Tests whether
#'     location measures differ while accounting for dependence.
#' }
#'
#' The bootstrap procedure resamples subjects (rows) with replacement, preserving
#' the dependence structure within subjects.
#'
#' FWER is controlled using Hochberg's sequentially rejective procedure by default,
#' which is more powerful than Holm's method while maintaining strong FWER control.
#'
#' Missing values are automatically removed (listwise deletion).
#'
#' @seealso \code{\link{rmmcppb}}, \code{\link{rmmcppbd}}, \code{\link{medpb}}
#'
#' @examples
#' # Two dependent groups
#' set.seed(123)
#' n <- 30
#' pre <- rnorm(n, mean = 100, sd = 15)
#' post <- pre + rnorm(n, mean = 5, sd = 10)  # Correlated
#'
#' # Using difference scores
#' dmedpb(pre, post, dif = TRUE, nboot = 1000, plotit = FALSE)
#'
#' # Using marginal distributions
#' dmedpb(pre, post, dif = FALSE, nboot = 1000, plotit = FALSE)
#'
#' # Multiple dependent groups
#' data <- matrix(rnorm(30 * 4), ncol = 4)
#' data[, 2] <- data[, 1] + rnorm(30, mean = 2, sd = 5)
#' data[, 3] <- data[, 1] + rnorm(30, mean = 4, sd = 5)
#' data[, 4] <- data[, 1] + rnorm(30, mean = 6, sd = 5)
#'
#' dmedpb(data, dif = TRUE, nboot = 1000, plotit = FALSE)
#'
#' @export
dmedpb<-function(x,y=NULL,alpha=.05,con=0,est=median,plotit=TRUE,dif=FALSE,grp=NA,
hoch=TRUE,nboot=NA,xlab="Group 1",ylab="Group 2",ylab.ebar=NULL,
pr=TRUE,SEED=TRUE,BA=FALSE,PCI=FALSE,EBAR=PCI,...){
#
#   Use a percentile bootstrap method to  compare
#   medians of dependent groups.
#
#   This is essentially the function rmmcppb, but set to compare medians
#   by default.
#   And it is adjusted to handle tied values.
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
#   otherwise, assume x is a matrix (n-by-J) or has list mode.
#
#  PCI=TRUE, if dif=TRUE and est=median, confidence intervals for difference scores are plottted
#           So this is like plotting error bars.
#
#
if(dif){
if(pr){
print("dif=T, so analysis is done on difference scores.")
print(" Each confidence interval has probability coverage 1-alpha.")
print(" Also note a sequentially rejective method is being used.")
}
temp<-rmmcppbd(x,y=y,alpha=alpha,con=con,est=est,plotit=plotit,grp=grp,
nboot=nboot,hoch=hoch,...)
output<-temp$output
con<-temp$con
}
if(!dif){
if(pr){
print("dif=F, so analysis is done on marginal distributions.")
print(" Each confidence interval has probability coverage 1-alpha.")
print(" Also note that a sequentially rejective method is being used")
}
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
if(pr)print("Taking bootstrap samples. Please wait.")
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
icl=round(alpha*nboot/2)+1
icu<-nboot-(icl-1)
cimat=matrix(NA,nrow=connum,ncol=2)
for (ic in 1:connum){
psihat[ic,]<-apply(bvec,1,bptdpsi,con[,ic])
tp=sort(psihat[ic,])
cimat[ic,1]=tp[icl]
cimat[ic,2]=tp[icu]
psihatcen[ic,]<-apply(bveccen,1,bptdpsi,con[,ic])
tvalcen[ic]<-sum((psihatcen[ic,]==0))/nboot
bias[ic]<-sum((psihatcen[ic,]>0))/nboot+sum((psihatcen[ic,]==0))/nboot-.5
tval[ic]<-sum((psihat[ic,]==0))/nboot
if(BA){
test[ic]<-sum((psihat[ic,]>0))/nboot+tval[ic]-.1*bias[ic]
if(test[ic]<0)test[ic]<-0
}
if(!BA)test[ic]<-sum((psihat[ic,]>0))/nboot+.5*tval[ic]
test[ic]<-min(test[ic],1-test[ic])
}
test<-2*test
ncon<-ncol(con)
dvec<-alpha/c(1:ncon)
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
#dvec[1]<-alpha/2
}
if(!EBAR){
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
}}
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
#print(psihat[ic,])
#icl=round(alpha*nboot/2)+1
#icl<-round(output[ic,4]*nboot/2)+1 # This adjustment causes confusion, it is not based on Hochberg
#icu<-nboot-(icl-1)
#output[ic,5]<-temp[icl]
#output[ic,6]<-temp[icu]
output[ic,5:6]<-cimat[ic,]
}
}
num.sig=nrow(output)
ior=order(output[,3],decreasing=TRUE)
for(j in 1:nrow(output)){
if(output[ior[j],3]<=output[ior[j],4])break
else num.sig=num.sig-1
}
num.sig<-sum(output[,3]<=output[,4])
#if(nrow(output)>1)ids=which(output[,3]<=output[,4])
if(EBAR ){
if(identical(est,median)){
if(dif){
plotCI(output[,2],ali=output[,5],aui=output[,6],xlab='Difference',ylab=ylab.ebar)
}}}
list(output=output,con=con,num.sig=num.sig)
}

# ==== med.effect ====
#' Robust Effect Size Based on Median and Robust Variance
#'
#' @description
#' Computes a robust analog of Cohen's d using the median (or Harrell-Davis
#' estimator) and a robust measure of variance (percentage bend midvariance
#' by default).
#'
#' @param x Numeric vector for first group.
#' @param y Numeric vector for second group.
#' @param HD Logical; if TRUE (default), use Harrell-Davis estimator for median.
#'   If FALSE, use sample median.
#' @param eq.var Logical; if FALSE (default), use explanatory measure of effect
#'   size (variance of estimates / variance of pooled data). If TRUE, use analog
#'   of Cohen's d (difference / pooled standard deviation).
#' @param nboot Number of bootstrap samples when sample sizes are unequal and
#'   eq.var=FALSE (default: 100).
#' @param loc.fun Location function to summarize bootstrap effect sizes when
#'   sample sizes are unequal (default: median).
#' @param varfun Variance function to use (default: pbvar, percentage bend variance).
#'
#' @return Numeric effect size value.
#'
#' @details
#' This function provides two types of robust effect sizes:
#'
#' \strong{Standardized difference (eq.var=TRUE)}:
#' The difference in medians divided by the pooled robust standard deviation,
#' analogous to Cohen's d but using robust estimates.
#'
#' \strong{Explanatory measure (eq.var=FALSE)}:
#' The square root of the ratio of the variance of the two location estimates
#' to the robust variance of the pooled data. This measures the proportion of
#' variability in the pooled data explained by group differences.
#'
#' When sample sizes are equal and eq.var=FALSE, the effect size is computed
#' directly. When sample sizes are unequal, bootstrap resampling is used to
#' create equal-sized samples, and the median effect size across bootstrap
#' samples is reported.
#'
#' Missing values are automatically removed before analysis.
#'
#' @references
#' Wilcox, R. R. (2022). Introduction to Robust Estimation and Hypothesis
#' Testing (5th ed.). Academic Press.
#'
#' @seealso \code{\link{MED.ES}}, \code{\link{akp.effect}}, \code{\link{pbvar}}
#'
#' @examples
#' # Effect size using medians
#' set.seed(123)
#' x <- rnorm(30, mean = 0, sd = 1)
#' y <- rnorm(30, mean = 0.8, sd = 1)
#'
#' # Explanatory measure
#' med.effect(x, y, HD = FALSE, eq.var = FALSE)
#'
#' # Cohen's d analog
#' med.effect(x, y, HD = FALSE, eq.var = TRUE)
#'
#' # Using Harrell-Davis estimator
#' med.effect(x, y, HD = TRUE)
#'
#' @export
med.effect<-function(x,y,HD=TRUE,eq.var=FALSE,nboot=100,loc.fun=median,varfun=pbvar){
#
#  Compute robust analog of Cohen's d using
#  the median and the percentage bend midvariance
#
#  HD=TRUE, use Harrell-Davis estimator
#  HD=FALSE, use usual sample median
#
#  eq.var=FALSE, use explanatory measure of effect size
#  eq.var=TRUE, use analog of Cohen's d.
#
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
if(HD){
e1=hd(x)
e2=hd(y)
}
if(!HD){
e1=median(x)
e2=median(y)
}
if(eq.var){
s1sq=varfun(x)
s2sq=varfun(y)
spsq<-(n1-1)*s1sq+(n2-1)*s2sq
sp<-sqrt(spsq/(n1+n2-2))
dval=(e1-e2)/sp
}
if(!eq.var){
n1=length(x)
n2=length(x)
if(n1==n2)dval=var(c(e1,e2))/varfun(c(x,y))
if(n1!=n2){
N=min(c(n1,n2))
vals=0
for(i in 1:nboot)vals[i]=med.effect.sub(sample(x,N),sample(y,N),HD=HD,varfun=varfun)
dval=loc.fun(vals)
}
dval=sqrt(dval)
}
dval
}

# ==== med.effect.sub ====
#' Helper Function for Median Effect Size (Internal)
#'
#' @description
#' Internal helper function for \code{\link{med.effect}} when sample sizes
#' are unequal. Computes effect size for resampled equal-sized groups.
#' Not intended for direct use.
#'
#' @param x Numeric vector for first group.
#' @param y Numeric vector for second group.
#' @param HD Logical; if TRUE, use Harrell-Davis estimator. If FALSE, use median.
#' @param varfun Variance function to use (default: \code{\link{pbvar}}).
#'
#' @return Numeric variance ratio (not square-rooted).
#'
#' @details
#' Computes the ratio of the variance of the two location estimates to the
#' variance of the pooled data. Called repeatedly by \code{\link{med.effect}}
#' when sample sizes are unequal, with equal-sized bootstrap samples.
#'
#' @seealso \code{\link{med.effect}}
#'
#' @keywords internal
med.effect.sub<-function(x,y,HD,varfun=pbvar){
if(HD){
e1=hd(x)
e2=hd(y)
}
if(!HD){
e1=median(x)
e2=median(y)
}
val=var(c(e1,e2))/varfun(c(x,y))
val
}

# ==== medcurve ====
#' Median Curve for Functional Data
#'
#' @description
#' Computes the median curve for functional data by identifying the curve
#' with maximum depth according to the functional boxplot method.
#'
#' @param x Matrix or data frame where each row is a functional observation
#'   (curve) and columns represent time points or locations.
#'
#' @return Numeric vector representing the median curve (same length as number
#'   of columns in x).
#'
#' @details
#' The median curve is defined as the observation with the highest depth
#' according to the band depth criterion used in functional boxplots
#' (see \code{\link{FBplot}}).
#'
#' If multiple curves tie for maximum depth, the result is the pointwise
#' average of those curves.
#'
#' This provides a robust center for functional data that is less sensitive
#' to outlying curves than the pointwise mean.
#'
#' @references
#' Sun, Y., & Genton, M. G. (2011). Functional boxplots. Journal of
#' Computational and Graphical Statistics, 20(2), 316-334.
#'
#' @seealso \code{\link{FBplot}}, \code{\link{dmedian}}
#'
#' @examples
#' # Median curve for functional data
#' set.seed(123)
#' # Generate 20 curves with 50 time points each
#' t_points <- seq(0, 1, length.out = 50)
#' curves <- matrix(nrow = 20, ncol = 50)
#' for (i in 1:20) {
#'   curves[i, ] <- sin(2 * pi * t_points) + rnorm(50, sd = 0.2)
#' }
#'
#' med_curve <- medcurve(curves)
#' plot(t_points, med_curve, type = "l", lwd = 2)
#'
#' @export
medcurve<-function(x){
#
# returns the median curve for functional data
#
chk=FBplot(x,plot=FALSE)$depth
id=which(chk==max(chk))
if(length(id)==1)est=x[id,]
if(length(id)>1)est=apply(x[id,],2,mean)
est
}

# ==== dmedian ====
#' Depth-Based Multivariate Median (Deepest Point)
#'
#' @description
#' Computes the multivariate median as the deepest point in the data cloud
#' according to a specified depth function. Provides a unique median for
#' multivariate data.
#'
#' @param x Matrix or data frame with two or more columns (multivariate data).
#' @param depfun Depth function to use (default: \code{\link{pdepth}} for
#'   projection depth). Another option is \code{\link{zdepth}} (zonoid depth).
#' @param ... Additional arguments passed to the depth function.
#'
#' @return A list with component:
#' \item{center}{The observation(s) with maximum depth (multivariate median).}
#'
#' @details
#' The depth-based median is defined as the data point with the highest depth
#' value. Unlike coordinate-wise medians (which may not correspond to any
#' actual data point), this approach identifies an actual observation from
#' the dataset.
#'
#' \strong{Depth functions}:
#' \itemize{
#'   \item \code{pdepth}: Projection depth (default) - computationally intensive
#'     but provides good affine-equivariant properties
#'   \item \code{zdepth}: Zonoid depth - faster but less robust
#' }
#'
#' For continuous variables, this typically returns a unique median. For
#' discrete or tied data, multiple points may achieve maximum depth.
#'
#' This median has good robustness properties and is affine-equivariant
#' (for appropriate depth functions), making it suitable for multivariate
#' location estimation.
#'
#' @references
#' Liu, R. Y., Parelius, J. M., & Singh, K. (1999). Multivariate analysis
#' by data depth: Descriptive statistics, graphics and inference. The Annals
#' of Statistics, 27(3), 783-858.
#'
#' @seealso \code{\link{pdepth}}, \code{\link{zdepth}}, \code{\link{medcurve}}
#'
#' @examples
#' # Depth-based median for bivariate data
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 2)
#'
#' # Using projection depth (default)
#' result <- dmedian(x)
#' result$center
#'
#' # Visualize
#' plot(x, main = "Bivariate Data with Depth Median")
#' points(result$center[1], result$center[2], pch = 19, col = "red", cex = 2)
#'
#' @export
dmedian<-function(x,depfun=pdepth,...){
#
# Compute the median based on the deepest point for the multivariate
# data in x
#
#  For continuous variables, this function returns a unique median
#
#  Projection distances are used by default.
# Another option is depfun=zdepth
#
if(is.null(dim(x)) || dim(x)==1)stop('x should be a matrix with two or more columns')
val=depfun(x,...)
id=which(val==max(val))
list(center=x[id,])
}

# ==== medpb.es ====
#' Multiple Comparisons for Medians with Effect Sizes
#'
#' @description
#' Performs multiple comparisons for medians (or other location estimators)
#' across J independent groups using percentile bootstrap, with shift-type
#' effect sizes (Q statistics) reported for each comparison.
#'
#' @inheritParams common_params
#' @param x Data in matrix or list mode. If a matrix, columns represent groups.
#'   If a list, each element is a group.
#' @param est Measure of location to use (default: median).
#' @param ... Optional arguments passed to the estimator function \code{est}.
#' @param con A J × d matrix of contrast coefficients. If con=0 (default),
#'   all pairwise comparisons are performed.
#' @param bhop Logical; if TRUE, use Benjamini-Hochberg FDR control instead
#'   of FWER control (default: FALSE).
#' @param INT Logical; advanced option for interaction contrasts (default: FALSE).
#'
#' @return A list with components:
#' \item{output}{Matrix with columns: Group(s), p-value, estimate difference,
#'   CI lower bound, CI upper bound, effect size (Q), and critical p-value.}
#' \item{con}{The contrast matrix used.}
#' \item{num.sig}{Number of significant comparisons.}
#'
#' @details
#' This function extends \code{\link{medpb}} by adding shift-type effect sizes
#' to the output. The effect size Q is a probability-based measure:
#'
#' For two groups with difference D = X - Y, let M be the population median of D,
#' and let F be the distribution of D - M. Then Q = F(M).
#'
#' \strong{Interpretation of Q}:
#' \itemize{
#'   \item Q = 0.5: No effect (median difference is zero)
#'   \item Q > 0.5: Positive shift (group 1 tends to be larger)
#'   \item Q < 0.5: Negative shift (group 2 tends to be larger)
#'   \item Distance from 0.5 indicates magnitude of effect
#' }
#'
#' This effect size represents how far the distribution has shifted from the
#' null hypothesis of no difference, expressed as a quantile of the centered
#' distribution.
#'
#' FWER control uses Rom's method by default. When \code{bhop=TRUE}, the
#' Benjamini-Hochberg FDR procedure is used instead.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{medpb}}, \code{\link{qhat}}, \code{\link{med.effect}}
#'
#' @examples
#' # Compare three groups with effect sizes
#' set.seed(123)
#' x1 <- rnorm(20, mean = 0)
#' x2 <- rnorm(20, mean = 0.8)
#' x3 <- rnorm(20, mean = 1.5)
#' data <- list(x1, x2, x3)
#'
#' result <- medpb.es(data, nboot = 1000)
#'
#' # View results with effect sizes
#' result$output
#'
#' @export
medpb.es<-function(x,alpha=.05,nboot=NA,grp=NA,est=median,con=0,bhop=FALSE,
SEED=TRUE,INT=FALSE,...){
#
#   Multiple comparisons for  J independent groups using medians.
#
#   A percentile bootstrap method with Rom's method is used.
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   est is the measure of location and defaults to the median
#   ... can be used to set optional arguments associated with est
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#
#   con can be used to specify linear contrasts; see the function lincon
#
#   Missing values are allowed.
#
#  A shift-type measure of effect size,Q,  is reported. No effect, Q=.5
#  For two groups, let D=X-Y, let M be the  population median of D.
#  Let F be the  distribution D-M. Then
#  Q=F(M). If the median of D is M, there is no effect.
#  Q represents a shift in location  to some relatively high or low quantile associated with $F_0$
#
con<-as.matrix(con)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in list mode or in matrix mode.')
if(!is.na(sum(grp))){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
x<-xx
}
J<-length(x)
tempn<-0
mvec<-NA
for(j in 1:J){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
tempn[j]<-length(temp)
x[[j]]<-temp
mvec[j]<-est(temp,...)
}
Jm<-J-1
#
# Determine contrast matrix
#
if(sum(con^2)==0){
ncon<-(J^2-J)/2
con<-matrix(0,J,ncon)
id<-0
for (j in 1:Jm){
jp<-j+1
for (k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
ncon<-ncol(con)
dvec<-alpha/c(1:ncon)
if(nrow(con)!=J)stop('Something is wrong with con; the number of rows does not match the number of groups.')
#  Determine nboot if a value was not specified
if(is.na(nboot)){
nboot<-5000
if(J <= 8)nboot<-4000
if(J <= 3)nboot<-2000
}
# Determine critical values
if(!bhop){
if(alpha==.05){
dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
}
}
if(bhop)dvec<-(ncon-c(1:ncon)+1)*alpha/ncon
bvec<-matrix(NA,nrow=J,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:J){
#print(paste('Working on group ',j))
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,est,...) # Bootstrapped values for jth group
}
test<-NA
bcon<-t(con)%*%bvec #ncon by nboot matrix
tvec<-t(con)%*%mvec
for (d in 1:ncon){
tv<-sum(bcon[d,]==0)/nboot
test[d]<-sum(bcon[d,]>0)/nboot+.5*tv
if(test[d]> .5)test[d]<-1-test[d]
}
test<-2*test
output<-matrix(0,ncon,7)
dimnames(output)<-list(NULL,c('con.num','psihat','p.value','p.crit','ci.lower','ci.upper','Q.effect'))
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,4]<-zvec
icl<-round(dvec[ncon]*nboot/2)+1
icu<-nboot-icl-1
output[,2]=tvec
for (ic in 1:ncol(con)){
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(bcon[ic,])
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
see=lin.ES(x,con=con[,ic],locfun=median)
if(!INT)output[ic,7]=lin.ES(x,con=con[,ic],locfun=median)$Effect.Size
if(!INT)output[ic,7]=lin.ES(x,con=con[,ic],locfun=median)$Effect.Size
if(INT){
id=which(con[,ic]!=0)
output[ic,7]=interQS(x[id],locfun=median,SEED=SEED)$Q.Effect
}
}
num.sig<-sum(output[,3]<=output[,4])
list(output=output,con=con,num.sig=num.sig)
}

# ==== dlinmed ====
#' Distribution and Test for Linear Combination of Dependent Medians
#'
#' @description
#' For dependent (repeated measures) variables, computes a linear combination
#' of the variables, plots its distribution, and tests whether the median
#' of the linear combination equals zero.
#'
#' @param x Data in matrix or list mode. If a matrix, columns are dependent
#'   variables (repeated measures).
#' @param con Vector of contrast coefficients. Must sum to zero.
#' @param xlab Label for x-axis in plot (default: 'DV').
#' @param ylab Label for y-axis in plot (default: '').
#' @param sym.test Logical; if TRUE, also perform symmetry test for the
#'   distribution (default: FALSE).
#' @param plotit Logical; if TRUE (default), plot the distribution.
#'
#' @return A list with components:
#' \item{median}{Sample median of the linear combination.}
#' \item{n}{Sample size.}
#' \item{ci.low}{Lower confidence limit for the median.}
#' \item{ci.up}{Upper confidence limit for the median.}
#' \item{p.value}{P-value for test that median equals zero.}
#' \item{Q.effect}{Effect size based on quantile shift.}
#' \item{sym.test}{Results of symmetry test (only if sym.test=TRUE).}
#'
#' @details
#' This function creates a new variable Y = sum(con * X) where X represents
#' the dependent variables and con are the contrast coefficients. It then:
#'
#' \enumerate{
#'   \item Plots the distribution of Y using adaptive kernel density estimation
#'   \item Tests H0: median(Y) = 0 using \code{\link{sintv2}}
#'   \item Computes a Q effect size measuring the quantile shift
#'   \item Optionally tests for symmetry of the distribution
#' }
#'
#' The contrast coefficients must sum to zero for the test to be meaningful
#' as a comparison among the dependent variables.
#'
#' Missing values are removed (listwise deletion).
#'
#' @seealso \code{\link{sintv2}}, \code{\link{depQS}}, \code{\link{Dqdif}},
#'   \code{\link{dmedpb}}
#'
#' @examples
#' # Test difference between two time points
#' set.seed(123)
#' n <- 30
#' time1 <- rnorm(n, mean = 100)
#' time2 <- time1 + rnorm(n, mean = 5, sd = 10)  # Correlated
#' data <- cbind(time1, time2)
#'
#' # Test time2 - time1
#' con <- c(-1, 1)
#' result <- dlinmed(data, con)
#'
#' # Three time points: test linear trend
#' time3 <- time1 + rnorm(n, mean = 10, sd = 10)
#' data3 <- cbind(time1, time2, time3)
#' con_trend <- c(-1, 0, 1)  # Linear trend contrast
#' dlinmed(data3, con_trend)
#'
#' @export
dlinmed<-function(x,con,xlab='DV',ylab='',sym.test=FALSE,plotit=TRUE){
#
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
if(plotit)akerd(L,xlab=xlab,ylab=ylab)
mt=sintv2(L)
sym=NULL
Q=depQS(L)
if(sym.test)sym=Dqdif(L)
list(median=mt$median,n=mt$n,ci.low=mt$ci.low,ci.up=mt$ci.up,
p.value=mt$p.value,Q.effect=Q$Q.effect,sym.test=sym)
}

# ==== wwmed ====
#' Two-Way Within-Within ANOVA for Medians with Multiple Comparisons
#'
#' @description
#' Performs all multiple comparisons for a J × K (two-way) within-subjects
#' (repeated measures) design using medians. Tests main effects for both
#' factors and interaction effects.
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data in matrix mode. Rows are subjects, columns are the J×K
#'   repeated measures in standard order.
#' @inheritParams common_params
#'
#' @return A list with components:
#' \item{Factor_A}{Results for Factor A main effect comparisons (from \code{\link{sintv2mcp}}).}
#' \item{Factor_B}{Results for Factor B main effect comparisons (from \code{\link{sintv2mcp}}).}
#' \item{Factor_AB}{Results for A × B interaction contrasts (from \code{\link{sintv2mcp}}).}
#'
#' @details
#' This function is a wrapper for \code{\link{sintv2mcp}} that generates
#' appropriate contrast matrices for a two-way within-subjects design using
#' \code{\link{con2way}}, then performs multiple comparisons on the medians
#' of linear combinations of the dependent variables.
#'
#' Each subject is measured under all J×K conditions. The standard order for
#' columns is: (1,1), (1,2), ..., (1,K), (2,1), (2,2), ..., (J,K).
#'
#' For each comparison, a sign test confidence interval is computed for the
#' median of the contrast. FWER is controlled separately for each family of
#' tests (Factor A, Factor B, and interaction).
#'
#' Missing values are handled via listwise deletion.
#'
#' @seealso \code{\link{sintv2mcp}}, \code{\link{con2way}}, \code{\link{wwwmed}},
#'   \code{\link{dmedpb}}
#'
#' @examples
#' # Two-way within-subjects design
#' set.seed(123)
#' n <- 20
#' # 2x2 design: 4 measurements per subject
#' data <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
#' # Add some main effects and interaction
#' data[, 2] <- data[, 2] + 0.5  # B effect
#' data[, 3] <- data[, 3] + 0.8  # A effect
#' data[, 4] <- data[, 4] + 1.5  # A + B + AB
#'
#' result <- wwmed(J = 2, K = 2, x = data)
#'
#' # View Factor A comparisons
#' result$Factor_A
#'
#' @export
wwmed<-function(J,K,x,alpha=.05){
#
# Do all multiple comparisons for a within-by-within design
# using medians
#
conM=con2way(J,K)
A=sintv2mcp(x,con=conM$conA,alpha=alpha)
B=sintv2mcp(x,con=conM$conB,alpha=alpha)
AB=sintv2mcp(x,con=conM$conAB,alpha=alpha)
list(Factor_A=A,Factor_B=B,Factor_AB=AB)
}

# ==== wwwmed ====
#' Three-Way Within-Within-Within ANOVA for Medians with Multiple Comparisons
#'
#' @description
#' Performs all multiple comparisons for a J × K × L (three-way) within-subjects
#' (repeated measures) design using medians. Tests main effects, two-way
#' interactions, and three-way interaction.
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param L Number of levels for Factor C (within-subjects).
#' @param x Data in matrix mode. Rows are subjects, columns are the J×K×L
#'   repeated measures in standard order.
#' @inheritParams common_params
#'
#' @return A list with components:
#' \item{Factor_A}{Results for Factor A main effect comparisons.}
#' \item{Factor_B}{Results for Factor B main effect comparisons.}
#' \item{Factor_C}{Results for Factor C main effect comparisons.}
#' \item{Factor_AB}{Results for A × B interaction contrasts.}
#' \item{Factor_AC}{Results for A × C interaction contrasts.}
#' \item{Factor_BC}{Results for B × C interaction contrasts.}
#' \item{Factor_ABC}{Results for A × B × C three-way interaction contrasts.}
#'
#' @details
#' This function extends \code{\link{wwmed}} to three-way within-subjects
#' designs. It uses \code{\link{con3way}} to generate appropriate contrast
#' matrices for all main effects, two-way interactions, and the three-way
#' interaction, then applies \code{\link{sintv2mcp}} to perform multiple
#' comparisons on medians.
#'
#' Each subject is measured under all J×K×L conditions. The standard order for
#' columns follows nested loops: varying C fastest, then B, then A.
#'
#' For each comparison, a sign test confidence interval is computed for the
#' median of the contrast. FWER is controlled separately for each family of
#' tests.
#'
#' Missing values are handled via listwise deletion.
#'
#' @seealso \code{\link{wwmed}}, \code{\link{sintv2mcp}}, \code{\link{con3way}},
#'   \code{\link{dmedpb}}
#'
#' @examples
#' # Three-way within-subjects design
#' set.seed(123)
#' n <- 15
#' # 2x2x2 design: 8 measurements per subject
#' data <- matrix(rnorm(n * 8), nrow = n, ncol = 8)
#' # Add some main effects
#' data[, 5:8] <- data[, 5:8] + 1  # Factor A effect
#'
#' result <- wwwmed(J = 2, K = 2, L = 2, x = data)
#'
#' # View Factor A comparisons
#' result$Factor_A
#'
#' @export
wwwmed<-function(J,K,L,x,alpha=.05){
#
# Do all multiple comparisons for a within-by-within design
# using medians
#
conM=con3way(J,K,L)
A=sintv2mcp(x,con=conM$conA,alpha=alpha)
B=sintv2mcp(x,con=conM$conB,alpha=alpha)
C=sintv2mcp(x,con=conM$conC,alpha=alpha)
AB=sintv2mcp(x,con=conM$conAB,alpha=alpha)
AC=sintv2mcp(x,con=conM$conAC,alpha=alpha)
BC=sintv2mcp(x,con=conM$conBC,alpha=alpha)
ABC=sintv2mcp(x,con=conM$conABC,alpha=alpha)
list(Factor_A=A,Factor_B=B,Factor_C=C,Factor_AB=AB,Factor_AC=AC,Factor_BC=BC,Factor_ABC=ABC)
}

# ==== runstest.med ====
#' Runs Test Based on Median
#'
#' @description
#' Performs a runs test for randomness based on whether values are below or
#' above the sample median. Tests the null hypothesis that the sequence is
#' random.
#'
#' @param x Numeric vector representing a time series or sequence.
#'
#' @return P-value from the runs test.
#'
#' @details
#' This function converts the data into a binary sequence based on whether
#' each value is below (0) or above (1) the sample median. It then applies
#' the runs test using \code{runs.test} from the tseries package.
#'
#' The runs test examines whether the number of runs (sequences of consecutive
#' 0s or 1s) is consistent with randomness. Too few runs suggests positive
#' autocorrelation (clustering of similar values), while too many runs suggests
#' negative autocorrelation (alternating pattern).
#'
#' \strong{Null hypothesis}: The sequence is random (values above and below
#' the median occur in random order).
#'
#' \strong{Alternative hypothesis}: The sequence is not random (there is a
#' pattern or autocorrelation).
#'
#' Missing values are automatically removed before analysis.
#'
#' @note Requires the tseries package to be installed.
#'
#' @seealso \code{runs.test} in package tseries
#'
#' @examples
#' # Random sequence (should not reject)
#' set.seed(123)
#' x_random <- rnorm(50)
#' runstest.med(x_random)
#'
#' # Sequence with positive autocorrelation (should reject)
#' x_clustered <- rep(c(rep(-1, 5), rep(1, 5)), 5)
#' runstest.med(x_clustered)
#'
#' @export
runstest.med<-function(x){
#
# runs test based on whether values are	< or > than the	median
#
library(tseries)
x=elimna(x)
n=length(x)
g=rep(1,n)
flag=x<median(x)
g[flag]=0
p=runs.test(as.factor(g))
p[3]
}

# ==== oph.astig.medianconvexpoly ====
#' Bivariate Median Confidence Regions for Ophthalmology Astigmatism Data
#'
#' @description
#' Specialized function for computing bootstrap confidence regions for bivariate
#' medians of prediction errors in astigmatism studies. Handles multiple
#' measurement methods/formulas simultaneously, designed specifically for
#' ophthalmology applications.
#'
#' @param m Matrix or data frame with an even number of columns. Each pair of
#'   consecutive columns represents one measurement formula/method (first two
#'   columns = formula 1, next two = formula 2, etc.).
#' @inheritParams common_params
#' @param plotit Logical; if TRUE (default), create plots of confidence regions.
#' @param xlab Label for x-axis in plots (default: 'V1').
#' @param ylab Label for y-axis in plots (default: 'V2').
#' @param nboot Number of bootstrap samples (default: 500).
#' @param MC Logical; if TRUE, use parallel processing via \code{\link{pdisMC}}
#'   for projection depth (default: FALSE).
#'
#' @return A list with components:
#' \item{centers}{List of bivariate median coordinates for each formula, with
#'   sample sizes.}
#' \item{conf.region.points}{List of matrices containing (x,y) coordinates of
#'   confidence region boundaries for each formula.}
#' \item{p.values}{List of p-values testing if median equals (0,0) for each formula.}
#'
#' @details
#' This function is tailored for ophthalmology research on astigmatism, where
#' prediction errors from different formulas are compared. Each formula produces
#' bivariate errors (e.g., prediction error in two meridians).
#'
#' \strong{Procedure for each formula}:
#' \enumerate{
#'   \item Compute bivariate median using \code{\link{smeancr.cord.oph}}
#'   \item Generate bootstrap samples and compute bootstrap medians
#'   \item Construct confidence region using projection depth and convex hull
#'   \item Test if (0,0) lies within the confidence region
#' }
#'
#' The confidence level is adjusted based on sample size to improve coverage:
#' smaller samples use lower nominal alpha for better actual coverage.
#'
#' Special handling for zeros: If more than one value in a column equals zero,
#' small jitter is added (divided by 100) to avoid computational issues while
#' keeping values very close to zero.
#'
#' When multiple formulas are analyzed, results are displayed in a 2×2 grid
#' of plots.
#'
#' @note This is a specialized function for ophthalmology applications. For
#' general bivariate median inference, see \code{\link{smeancr}} or
#' \code{\link{dmedian}}.
#'
#' @seealso \code{\link{smeancr.cord.oph}}, \code{\link{pdis}}, \code{\link{dmedian}}
#'
#' @examples
#' # Simulated astigmatism prediction errors for 2 formulas
#' set.seed(123)
#' n <- 50
#' # Formula 1: two error components
#' formula1_v1 <- rnorm(n, mean = 0.1, sd = 0.5)
#' formula1_v2 <- rnorm(n, mean = -0.05, sd = 0.4)
#' # Formula 2: two error components
#' formula2_v1 <- rnorm(n, mean = -0.2, sd = 0.6)
#' formula2_v2 <- rnorm(n, mean = 0.15, sd = 0.5)
#'
#' # Combine into matrix
#' data <- cbind(formula1_v1, formula1_v2, formula2_v1, formula2_v2)
#'
#' # Compute confidence regions
#' result <- oph.astig.medianconvexpoly(data, nboot = 200, plotit = FALSE)
#'
#' # View median for first formula
#' result$centers[[1]]
#'
#' @export
oph.astig.medianconvexpoly<-function(m,alpha=.05,plotit=TRUE,xlab='V1',ylab='V2',nboot=500,MC=FALSE,
SEED=TRUE){
#
# This function is designed to compute confidence region for  the median of bivariate data.
#  The function is designed  for dealing with
#  prediction errors when dealing with astigmatism.
#
# alpha=.05 means that the function
#  determine the 1-.05=.95 confidence region for the median of the data cloud.
#  This is done for each formula
#  Assume m is a matrix or data frame having
#  J columns. First two columns first formula, next two columns next formula..
#
#  So J should be an even integer
#
J=ncol(m)
N=J/2
if(N != floor(N))stop('Should have an even number of columns')
# Check of zeros and jitter if more than two. Jitter values are divided by 100 to be
# sure they are sufficiently close to zero.
for(j in 1:J){
id=which(m[,j]==0)
if(length(id>1)){
d=jitter(m[id,j])/100
m[id,j]=d
}}
if(plotit){
plot.new()
if(N>1)par(mfrow=c(2,2))
}
region=list()
centers=list()
val=list()
pv=list()
CENTERS=list()
id=c(-1,0)
for(j in 1:N){
id=id+2
a=smeancr.cord.oph(m[,id],SEED=SEED,plotit=FALSE,xlab=xlab,ylab=ylab,nboot=nboot)
centers[[j]]=a$center
region[[j]]=a$conf.region.points
val[[j]]=a$boot.vals
centers[[j]]=a$center
n=nrow(elimna(m[,id]))
n=as.integer(n)
CENTERS[[j]]=c(a$center,n)
names(CENTERS[[j]])=c('V1','V2','N')
pv[[j]]=a$p.value
}
VAL=val
#if(N>1)par(mfrow=c(2,2))
id=c(-1,0)
for(j in 1:N){
id=id+2
n=nrow(m[,id])
crit.level<-.05
if(n<=120)crit.level<-.045
if(n<=80)crit.level<-.04
if(n<=60)crit.level<-.035
if(n<=40)crit.level<-.03
if(n<=30)crit.level<-.025
if(n<=20)crit.level<-.02
ic<-round((1-crit.level)*nboot)
val=VAL[[j]]
est=centers[[j]]
temp3<-est
ic<-round((1-crit.level)*nboot)
if(!MC)temp<-pdis(val,center=est)
if(MC)temp<-pdisMC(val,center=est)
temp.dis<-order(temp)
xx<-val[temp.dis[1:ic],]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
region[[j]]=xx[temp,]
if(plotit){
#if(N>1)par(mfrow=c(2,2))
plot(val[,1],val[,2],xlab=xlab,ylab=ylab)
points(temp3[1],temp3[2],pch="+")
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}
}
par(mfrow=c(1,1))
list(centers=CENTERS,conf.region.points=region,p.values=pv)
}

# ==== MED.ES ====
#' One-Sample Effect Size Based on Median
#'
#' @description
#' Computes a one-sample effect size analogous to Cohen's d, using the median
#' and either median absolute deviation (MAD) or Winsorized standard deviation
#' (rescaled to estimate the population standard deviation under normality).
#'
#' @param x Numeric vector of data.
#' @param tr Trimming proportion for Winsorized standard deviation if MAD equals
#'   zero (default: 0.25).
#' @param null.val Null hypothesis value to subtract from the median (default: 0).
#' @param est Estimator for location (default: median).
#'
#' @return Numeric effect size value: (est(x) - null.val) / scale, where scale
#'   is MAD if non-zero, otherwise Winsorized SD.
#'
#' @details
#' This function computes a standardized effect size for a single sample, measuring
#' how many standard deviations the location estimate is from the null value.
#'
#' The scale estimate is MAD (median absolute deviation) rescaled by a factor
#' that makes it consistent for estimating the standard deviation under normality.
#' If MAD equals zero (can occur with highly discrete data or many ties), the
#' Winsorized standard deviation is used instead, also rescaled for consistency
#' with the normal distribution.
#'
#' This provides a robust alternative to (mean - null.val) / SD for effect
#' size calculation, resistant to outliers and heavy-tailed distributions.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{med.effect}}, \code{\link{akp.effect}}, \code{\link{mad}}
#'
#' @examples
#' # One-sample effect size
#' set.seed(123)
#' x <- rnorm(50, mean = 5, sd = 2)
#'
#' # Effect size relative to 0
#' MED.ES(x, null.val = 0)
#'
#' # Effect size relative to 5
#' MED.ES(x, null.val = 5)
#'
#' # Using trimmed mean instead of median
#' MED.ES(x, null.val = 0, est = function(x) mean(x, trim = 0.2))
#'
#' @export
MED.ES<-function(x,tr=.25,null.val=0,est=median){
#
# One-sample effect size analog of Cohen's d based on the median
# and either MAD or Winsorized standard deviation rescaled to estimate the standard deviation when
# sampling from a normal distribution
#
x=elimna(x)
e=est(x)
bot=mad(x)
if(bot==0)bot=winsdN(x,tr=tr)
if(bot==0)stop('Both measures of scale are equal to zero')
es=(e-null.val)/bot
es
}

# ==== oph.dep.comMedAE ====
