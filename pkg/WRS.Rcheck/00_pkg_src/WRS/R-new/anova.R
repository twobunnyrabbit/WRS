# ==============================================================================
# ANOVA and Related Functions
# ==============================================================================
#
# This module contains functions for ANOVA analyses including:
#   - One-way ANOVA: t1way(), t1waybt(), MEDanova()
#   - Two-way ANOVA: t2way(), bwtrim(), bdm2way()
#   - Three-way ANOVA: t3way()
#   - Repeated measures ANOVA: rmanova(), rmanovab()
#   - Bootstrap ANOVA: pbanova(), bbwtrim()
#   - Specialized ANOVA: anova1(), anova.nestA(), etc.
#
# Extracted: 2025-12-30
# Functions: 57
# ==============================================================================

#' Rust-Fligner Rank-Based ANOVA
#'
#' Performs a Rust-Fligner ANOVA using ranks for comparing J independent groups.
#' This is a robust nonparametric alternative to traditional ANOVA.
#'
#' @param x Data in list mode (each element is a group) or a matrix/data frame
#'   where columns correspond to groups.
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{test}{The chi-square test statistic.}
#'   \item{p.value}{The p-value for the test.}
#'   \item{df}{Degrees of freedom (J-1).}
#'
#' @details
#' The Rust-Fligner test is a rank-based procedure for testing equality of
#' distributions across J independent groups. Missing values are automatically
#' removed. The test assumes tied values occur with probability zero; a warning
#' is issued if ties are detected.
#'
#' The test statistic has an asymptotic chi-square distribution with J-1 degrees
#' of freedom.
#'
#' @seealso \code{\link{t1way}}, \code{\link{MEDanova}}
#'
#' @export
#' @examples
#' # Compare three groups
#' set.seed(123)
#' x1 <- rnorm(20, mean=0)
#' x2 <- rnorm(20, mean=0.5)
#' x3 <- rnorm(20, mean=1)
#' rfanova(list(x1, x2, x3))
rfanova<-function(x,grp=0){
#
#  Perform Rust-Fligner anova using ranks.
#  x is assumed to have list mode. x[[j]] contains data for jth group.
#  or x is a matrix with columns corresponding to groups.
#
#  missing values are automatically removed.
#  Tied values are a ssumed to occur with probability zero.
#
if(!is.list(x))x<-listm(x)
chk=tlist(x)
if(chk!=0)print("Warning: tied values detected")
xall<-NA
if(sum(grp)==0){
J<-length(x)
grp<-c(1:J)
}
if(sum(grp)>0)J<-length(grp)
nval<-1
nrat<-1
nmax<-0
rbar<-1
mrbar<-0
for (j in grp){
temp<-x[[j]]
temp<-temp[!is.na(temp)] #Missing values are removed.
nrat[j]<-(length(temp)-1)/length(temp)
nval[j]<-length(temp)
if(j==grp[1])xall<-temp
if(j!=grp[1])xall<-c(xall,temp)
if(length(temp)>nmax)nmax<-length(temp)
}
pv<-array(NA,c(J,nmax,J))
tv<-matrix(NA,J,nmax)
rv<-matrix(0,J,nmax)
for (i in 1:J){
data<-x[[i]]
data<-data[!is.na(data)]
for (j in 1:length(data)){
tempr<-data[j]-xall
rv[i,j]<-length(tempr[tempr>=0])
for (l in 1:J){
templ<-x[[l]]
templ<-templ[!is.na(templ)]
temp<-data[j]-templ
pv[i,j,l]<-length(temp[temp>=0])
}
tv[i,j]<-sum(pv[i,j,])-pv[i,j,i]
}
rbar[i]<-sum(rv[i,])/nval[i]
mrbar<-mrbar+sum(rv[i,])
}
amat<-matrix(0,J,J)
for(i in 1:J){
temptv<-tv[i,]
temptv<-temptv[!is.na(temptv)]
amat[i,i]<-(length(temptv)-1)*var(temptv)
for (l in 1:J){
tempp<-pv[l,,i]
tempp<-tempp[!is.na(tempp)]
if(l!=i){
amat[i,i]<-amat[i,i]+(length(tempp)-1)*var(tempp)
}}
for (j in 1:J){
if(j>i){
for (l in 1:J){
temp1<-pv[l,,i]
temp2<-pv[l,,j]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
#if(i!=l && l!=j)amat[i,j]<-(length(temp1)-1)*var(temp1,temp2)
if(i!=l && l!=j)amat[i,j]<-amat[i,j]+(length(temp1)-1)*var(temp1,temp2)
}
temp1<-pv[i,,j]
temp2<-tv[i,]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
amat[i,j]<-amat[i,j]-(length(temp1)-1)*var(temp1,temp2)
temp1<-pv[j,,i]
temp2<-tv[j,]
temp1<-temp1[!is.na(temp1)]
temp2<-temp2[!is.na(temp2)]
amat[i,j]<-amat[i,j]-(length(temp1)-1)*var(temp1,temp2)
}
amat[j,i]<-amat[i,j]
}}
N<-sum(nval)
amat<-amat/N^3
amati<-ginv(amat)
uvec<-1
mrbar<-mrbar/N
for (i in 1:J)uvec[i]<-nval[i]*(rbar[i]-mrbar)/(N*(N+1))
testv<-N*prod(nrat)*uvec%*%amati%*%uvec
test<-testv[1,1]
df<-J-1
siglevel<-1-pchisq(test,df)
list(test=test,p.value=siglevel,df=df)
}

#' Agresti-Pendergast Rank Test for Dependent Groups
#'
#' Performs an Agresti-Pendergast rank test for comparing J dependent (repeated
#' measures) groups. This is a nonparametric alternative to repeated measures ANOVA.
#'
#' @param data Data in matrix form (n by J) or in list mode where each element
#'   contains data for one group. For dependent groups, rows represent matched
#'   observations.
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{FTEST}{The F test statistic.}
#'   \item{df1}{Numerator degrees of freedom (J-1).}
#'   \item{df2}{Denominator degrees of freedom ((J-1)(n-1)).}
#'   \item{p.value}{The p-value for the test.}
#'
#' @details
#' The Agresti-Pendergast test uses ranks to test for differences among J
#' dependent groups. It is based on an F-statistic with (J-1) and (J-1)(n-1)
#' degrees of freedom.
#'
#' For n <= 20, the function suggests using \code{bprm} instead for better
#' small sample performance.
#'
#' @seealso \code{\link{rmanova}}, \code{\link{bprm}}
#'
#' @export
#' @examples
#' # Three dependent measurements
#' set.seed(123)
#' time1 <- rnorm(15, mean=10)
#' time2 <- time1 + rnorm(15, mean=0.5, sd=0.5)
#' time3 <- time1 + rnorm(15, mean=1, sd=0.5)
#' apanova(cbind(time1, time2, time3))
apanova<-function(data,grp=0){
#
#  Perform Agresti-Pendergast rank test for J dependent groups
#  The data are assumed to be stored in an n by J matrix or
#  in list mode. In the latter case, length(data)=J.
#
if(is.list(data)){
x<-matrix(0,length(data[[1]]),length(data))
for (j in 1:length(data))x[,j]<-data[[j]]
}
if(is.matrix(data))x<-data
if(sum(grp==0))grp<-c(1:ncol(x))
x<-x[,grp]
J<-ncol(x)
n<-nrow(x)
if(n<=20)print("With n<=20, suggest using bprm")
rm<-matrix(rank(x),n,J)
rv<-apply(rm,2,mean)
sm<-(n-1)*winall(rm,tr=0)$cov/(n-J+1)
jm1<-J-1
cv<-diag(1,jm1,J)
for (i in 2:J){
k<-i-1
cv[k,i]<--1
}
cr<-cv%*%rv
ftest<-n*t(cr)%*%solve(cv%*%sm%*%t(cv))%*%cr/(J-1)
df1<-J-1
df2<-(J-1)*(n-1)
siglevel<-1-pf(ftest,df1,df2)
list(FTEST=ftest,df1=df1,df2=df2,p.value=siglevel)
}

#' Heteroscedastic One-Way Repeated Measures ANOVA for Trimmed Means
#'
#' Performs a heteroscedastic one-way repeated measures ANOVA comparing trimmed
#' means across J dependent groups. Uses adjusted degrees of freedom to handle
#' heteroscedasticity and violations of sphericity.
#'
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{num.groups}{Number of groups being compared.}
#'   \item{test}{The test statistic.}
#'   \item{df}{Vector of degrees of freedom c(df1, df2).}
#'   \item{p.value}{The p-value for the test.}
#'   \item{tmeans}{Vector of trimmed means for each group.}
#'   \item{ehat}{Epsilon-hat estimate for sphericity adjustment.}
#'   \item{etil}{Epsilon-tilde estimate (adjusted version of epsilon-hat).}
#'
#' @details
#' This function implements a heteroscedastic repeated measures ANOVA for trimmed
#' means. It does not assume equal variances or sphericity. The degrees of freedom
#' are adjusted using epsilon-hat and epsilon-tilde corrections based on the
#' Winsorized covariance matrix.
#'
#' Data can be provided as an n by J matrix (rows are subjects, columns are
#' conditions) or in list mode. Use the \code{grp} argument to select a subset
#' of groups for comparison.
#'
#' @seealso \code{\link{rmanovab}}, \code{\link{apanova}}, \code{\link{bprm}}
#'
#' @export
#' @examples
#' # Three time points
#' set.seed(123)
#' n <- 20
#' time1 <- rnorm(n, mean=10)
#' time2 <- time1 + rnorm(n, mean=0.5, sd=1.2)
#' time3 <- time1 + rnorm(n, mean=1.5, sd=1.5)
#' rmanova(cbind(time1, time2, time3), tr=0.2)
rmanova<-function(x,tr=.2,grp=c(1:length(x))){
#
#  A heteroscedastic one-way repeated measures ANOVA for trimmed means.
#
#  The data are assumed to be stored in $x$ which can
#  be either an n by J matrix, or an R variable having list mode.
#  If the data are stored in list mode,
#  length(x) is assumed to correspond to the total number of groups.
#  By default, the null hypothesis is that all group have a common mean.
#  To compare a subset of the groups, use grp to indicate which
#  groups are to be compared. For example, if you type the
#  command grp<-c(1,3,4), and then execute this function, groups
#  1, 3, and 4 will be compared with the remaining groups ignored.
#
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
J<-length(grp)  # The number of groups to be compared
#print("The number of groups to be compared is")
#print(J)
m1<-matrix(x[[grp[1]]],length(x[[grp[1]]]),1)
for(i in 2:J){     # Put the data into an n by J matrix
m2<-matrix(x[[grp[i]]],length(x[[i]]),1)
m1<-cbind(m1,m2)
}
}
if(is.matrix(x)){
if(length(grp)<ncol(x))m1<-as.matrix(x[,grp])
if(length(grp)>=ncol(x))m1<-as.matrix(x)
J<-ncol(x)
#print("The number of groups to be compared is")
#print(J)
}
#
#  Raw data are now in the matrix m1
#
m2<-matrix(0,nrow(m1),ncol(m1))
xvec<-1
g<-floor(tr*nrow(m1))  #2g is the number of observations trimmed.
for(j in 1:ncol(m1)){  # Putting Winsorized values in m2
m2[,j]<-winval(m1[,j],tr)
xvec[j]<-mean(m1[,j],tr)
}
xbar<-mean(xvec)
qc<-(nrow(m1)-2*g)*sum((xvec-xbar)^2)
m3<-matrix(0,nrow(m1),ncol(m1))
m3<-sweep(m2,1,apply(m2,1,mean))  # Sweep out rows
m3<-sweep(m3,2,apply(m2,2,mean))  # Sweep out columns
m3<-m3+mean(m2)  # Grand Winsorized mean swept in
qe<-sum(m3^2)
test<-(qc/(qe/(nrow(m1)-2*g-1)))
#
#  Next, estimate the adjusted degrees of freedom
#
v<-winall(m1,tr=tr)$cov
vbar<-mean(v)
vbard<-mean(diag(v))
vbarj<-1
for(j in 1:J){
vbarj[j]<-mean(v[j,])
}
A<-J*J*(vbard-vbar)^2/(J-1)
B<-sum(v*v)-2*J*sum(vbarj^2)+J*J*vbar^2
ehat<-A/B
etil<-(nrow(m2)*(J-1)*ehat-2)/((J-1)*(nrow(m2)-1-(J-1)*ehat))
etil<-min(1.,etil)
df1<-(J-1)*etil
df2<-(J-1)*etil*(nrow(m2)-2*g-1)
siglevel<-1-pf(test,df1,df2)
list(num.groups=J,test=test,df=c(df1,df2),p.value=siglevel,tmeans=xvec,ehat=ehat,etil=etil)
}

#' Bootstrap-t for Repeated Measures ANOVA with Trimmed Means
#'
#' Performs a bootstrap-t test for comparing trimmed means of J dependent
#' (repeated measures) groups. Provides an alternative to \code{rmanova} with
#' better control of Type I error rates.
#'
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{teststat}{The test statistic.}
#'   \item{crit}{The critical value from the bootstrap distribution.}
#'   \item{p.value}{The bootstrap p-value.}
#'
#' @details
#' This function uses the bootstrap-t method to test for differences among J
#' dependent groups using trimmed means. Data are centered before bootstrap
#' resampling to impose the null hypothesis.
#'
#' The bootstrap samples are generated by resampling rows (subjects) with
#' replacement. A fixed seed (set.seed(2)) is used for reproducibility.
#'
#' @seealso \code{\link{rmanova}}, \code{\link{t1waybt}}
#'
#' @export
#' @examples
#' # Three dependent measurements
#' set.seed(123)
#' n <- 20
#' time1 <- rnorm(n, mean=10)
#' time2 <- time1 + rnorm(n, mean=0.5, sd=1)
#' time3 <- time1 + rnorm(n, mean=1, sd=1)
#' rmanovab(cbind(time1, time2, time3), tr=0.2, nboot=500)
rmanovab<-function(x,tr=.2,alpha=.05,grp=0,nboot=599){
#
#   A bootstrap-t for comparing the trimmed means of dependent groups.
#   By default, 20% trimming is used with B=599 bootstrap samples.
#
#   The optional argument grp is used to select a subset of the groups
#   and exclude the rest.
#
#   x can be an n by J matrix or it can have list mode
#
if(is.data.frame(x))x=as.matrix(x)
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x))x=matl(x)
#{
#if(sum(grp)==0)grp<-c(1:length(x))
# put the data in an n by J matrix
#mat<-matrix(0,length(x[[1]]),length(grp))
#for (j in 1:length(grp))mat[,j]<-x[[grp[j]]]
#}
if(is.matrix(x)){
if(sum(grp)==0)grp<-c(1:ncol(x))
mat<-x[,grp]
}
mat=elimna(mat)
J<-ncol(mat)
connum<-(J^2-J)/2
bvec<-matrix(0,connum,nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(nrow(mat),size=nrow(mat)*nboot,replace=TRUE),nrow=nboot)
xcen<-matrix(0,nrow(mat),ncol(mat))
for (j in 1:J)xcen[,j]<-mat[,j]-mean(mat[,j],tr) #Center data
bvec<-apply(data,1,tsubrmanovab,xcen,tr)
# bvec is vector of nboot  bootstrap test statistics.
icrit<-round((1-alpha)*nboot)
bvec<-sort(bvec)
crit<-bvec[icrit]
test<-rmanova(mat,tr,grp)$test
pv=mean(test<=bvec)
list(teststat=test,crit=crit,p.value=pv)
}

#' Compute Test Statistic for Bootstrap Sample (Internal Helper)
#'
#' Computes the test statistic for a bootstrap sample when comparing dependent
#' groups using trimmed means. Used internally by \code{rmanovab}.
#'
#' @param isub Vector of bootstrap indices (integers from 1 to n).
#' @param x Matrix of data (centered).
#' @inheritParams common-params
#'
#' @return The test statistic value for the bootstrap sample.
#'
#' @keywords internal
tsubrmanovab<-function(isub,x,tr){
#
#  Compute test statistic for trimmed means
#  when comparing dependent groups.
#  By default, 20% trimmed means are used.
#  isub is a vector of length n,
#  a bootstrap sample from the sequence of integers
#  1, 2, 3, ..., n
#
#  This function is used by rmanovab
#
tsub<-rmanovab1(x[isub,],tr=tr)$test
tsub
}

#' Repeated Measures ANOVA Helper Function (Internal)
#'
#' Internal version of \code{rmanova} used by bootstrap procedures. Computes
#' the test statistic without additional output formatting.
#'
#' @inheritParams rmanova
#'
#' @return A list with components:
#'   \item{test}{The test statistic.}
#'   \item{df}{Vector of degrees of freedom c(df1, df2).}
#'   \item{p.value}{The p-value.}
#'   \item{tmeans}{Vector of trimmed means.}
#'   \item{ehat}{Epsilon-hat estimate.}
#'   \item{etil}{Epsilon-tilde estimate.}
#'
#' @keywords internal
rmanovab1<-function(x,tr=.2,grp=c(1:length(x))){
#
#  A heteroscedastic one-way repeated measures ANOVA for trimmed means.
#
#  The data are assumed to be stored in $x$ which can
#  be either an n by J matrix, or an R variable having list mode.
#  If the data are stored in list mode,
#  length(x) is assumed to correspond to the total number of groups.
#  By default, the null hypothesis is that all group have a common mean.
#  To compare a subset of the groups, use grp to indicate which
#  groups are to be compared. For example, if you type the
#  command grp<-c(1,3,4), and then execute this function, groups
#  1, 3, and 4 will be compared with the remaining groups ignored.
#
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
J<-length(grp)  # The number of groups to be compared
m1<-matrix(x[[grp[1]]],length(x[[grp[1]]]),1)
for(i in 2:J){     # Put the data into an n by J matrix
m2<-matrix(x[[grp[i]]],length(x[[i]]),1)
m1<-cbind(m1,m2)
}
}
if(is.matrix(x)){
if(length(grp)<ncol(x))m1<-as.matrix(x[,grp])
if(length(grp)>=ncol(x))m1<-as.matrix(x)
J<-ncol(x)
}
#
#  Raw data are now in the matrix m1
#
m2<-matrix(0,nrow(m1),ncol(m1))
xvec<-1
g<-floor(tr*nrow(m1))  #2g is the number of observations trimmed.
for(j in 1:ncol(m1)){  # Putting Winsorized values in m2
m2[,j]<-winval(m1[,j],tr)
xvec[j]<-mean(m1[,j],tr)
}
xbar<-mean(xvec)
qc<-(nrow(m1)-2*g)*sum((xvec-xbar)^2)
m3<-matrix(0,nrow(m1),ncol(m1))
m3<-sweep(m2,1,apply(m2,1,mean))  # Sweep out rows
m3<-sweep(m3,2,apply(m2,2,mean))  # Sweep out columns
m3<-m3+mean(m2)  # Grand Winsorized mean swept in
qe<-sum(m3^2)
test<-(qc/(qe/(nrow(m1)-2*g-1)))
#
#  Next, estimate the adjusted degrees of freedom
#
v<-winall(m1)$cov
vbar<-mean(v)
vbard<-mean(diag(v))
vbarj<-1
for(j in 1:J){
vbarj[j]<-mean(v[j,])
}
A<-J*J*(vbard-vbar)^2/(J-1)
B<-sum(v*v)-2*J*sum(vbarj^2)+J*J*vbar^2
ehat<-A/B
etil<-(nrow(m2)*(J-1)*ehat-2)/((J-1)*(nrow(m2)-1-(J-1)*ehat))
etil<-min(1.,etil)
df1<-(J-1)*etil
df2<-(J-1)*etil*(nrow(m2)-2*g-1)
siglevel<-1-pf(test,df1,df2)
list(test=test,df=c(df1,df2),p.value=siglevel,tmeans=xvec,ehat=ehat,etil=etil)
}

#' Bishop-Dudewicz ANOVA - Stage 1 Sample Size Determination
#'
#' Determines the required sample size for each group in the first stage of
#' the Bishop-Dudewicz two-stage ANOVA procedure to achieve specified power.
#'
#' @inheritParams common-params
#' @param power Desired power (1 - beta) for the test (default: 0.9).
#' @param delta Effect size parameter - the minimum detectable difference.
#'
#' @return A list with components:
#'   \item{N}{Vector of required sample sizes for each group.}
#'   \item{d}{Computed effect size parameter.}
#'   \item{crit}{Critical value for the test.}
#'
#' @details
#' The Bishop-Dudewicz procedure is a two-stage sampling method that guarantees
#' a specified power for detecting differences among J groups. Stage 1 uses
#' pilot data to determine how many additional observations are needed in each
#' group.
#'
#' The required sample size for each group is based on the variance estimate
#' from the pilot data and the specified effect size \code{delta}.
#'
#' @seealso \code{\link{bdanova2}}, \code{\link{anova1}}
#'
#' @export
#' @examples
#' # Pilot data from three groups
#' set.seed(123)
#' x1 <- rnorm(10, mean=0, sd=1)
#' x2 <- rnorm(10, mean=0.5, sd=1.2)
#' x3 <- rnorm(10, mean=1, sd=0.8)
#' bdanova1(list(x1, x2, x3), delta=0.5, power=0.9)
bdanova1<-function(x,alpha=.05,power=.9,delta=NA){
#
#  Do the first stage of a Bishop-Dudewicz ANOVA method.
#  That is, based on the data in x
#  determine N_j, the number of observations needed
#  in the jth group to achieve power 1-beta.
#
#  The argument x is assumed to have list mode or the
#  data is assumed to be stored in an n by J matrix
#
if(is.na(delta))stop("A value for delta was not specified")
if(!is.list(x)){
if(!is.matrix(x))stop("Data must be stored in matrix or in list mode")
}
y<-x
if(is.list(y))y=matl(y)
x<-list()
for(j in 1:ncol(y))x[[j]]<-elimna(y[,j])
nvec<-NA
svec<-NA
J<-length(x)
for(j in 1:length(x)){
nvec[j]<-length(x[[j]])
svec[j]<-var(x[[j]])
}
nu<-nvec-1
nu1<-sum(1/(nu-2))
nu1<-J/nu1+2
A<-(J-1)*nu1/(nu1-2)
B<-(nu1^2/J)*(J-1)/(nu1-2)
C<-3*(J-1)/(nu1-4)
D<-(J^2-2*J+3)/(nu1-2)
E<-B*(C+D)
M<-(4*E-2*A^2)/(E-A^2-2*A)
L<-A*(M-2)/M
f<-qf(1-alpha,L,M)
crit<-L*f
b<-(nu1-2)*crit/nu1
zz<-qnorm(power)
A<-.5*(sqrt(2)*zz+sqrt(2*zz^2+4*(2*b-J+2)))
B<-A^2-b
d<-((nu1-2)/nu1)*delta/B
N<-NA
for(j in 1:length(x)){
N[j]<-max(c(nvec[j]+1,floor(svec[j]/d)+1))
}
list(N=N,d=d,crit=crit)
}

#' Conventional One-Way ANOVA
#'
#' Performs a standard one-way ANOVA F-test comparing means across J independent
#' groups. This is the classic parametric ANOVA assuming normality and equal
#' variances.
#'
#' @param x Data in list mode (each element is a group), matrix, or data frame
#'   where columns correspond to groups.
#'
#' @return A list with components:
#'   \item{F.test}{The F test statistic.}
#'   \item{p.value}{The p-value for the test.}
#'   \item{df1}{Numerator degrees of freedom (J-1).}
#'   \item{df2}{Denominator degrees of freedom (N-J).}
#'   \item{MSBG}{Mean square between groups.}
#'   \item{MSWG}{Mean square within groups.}
#'
#' @details
#' This function implements the traditional one-way ANOVA F-test. It assumes:
#' \itemize{
#'   \item Normality within each group
#'   \item Equal variances across groups (homoscedasticity)
#'   \item Independent observations
#' }
#'
#' For robust alternatives that don't require these assumptions, see
#' \code{\link{t1way}}, \code{\link{MEDanova}}, or \code{\link{pbanova}}.
#'
#' Missing values are automatically removed.
#'
#' @seealso \code{\link{t1way}}, \code{\link{MEDanova}}, \code{\link{pbanova}}
#'
#' @export
#' @examples
#' # Three groups
#' set.seed(123)
#' x1 <- rnorm(20, mean=10, sd=2)
#' x2 <- rnorm(20, mean=12, sd=2)
#' x3 <- rnorm(20, mean=11, sd=2)
#' anova1(list(x1, x2, x3))
anova1<-function(x){
#
# conventional one-way anova
#
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
x=elimna(x)
A<-0
B<-0
C<-0
N<-0
for(j in 1:length(x)){
N<-N+length(x[[j]])
A<-A+sum(x[[j]]^2)
B<-B+sum(x[[j]])
C<-C+(sum(x[[j]]))^2/length(x[[j]])
}
SST<-A-B^2/N
SSBG<-C-B^2/N
SSWG<-A-C
nu1<-length(x)-1
nu2<-N-length(x)
MSBG<-SSBG/nu1
MSWG<-SSWG/nu2
FVAL<-MSBG/MSWG
pvalue<-1-pf(FVAL,nu1,nu2)
list(F.test=FVAL,p.value=pvalue,df1=nu1,df2=nu2,MSBG=MSBG,MSWG=MSWG)
}

#' Bishop-Dudewicz ANOVA - Stage 2 Test
#'
#' Performs the second stage of the Bishop-Dudewicz two-stage ANOVA procedure,
#' combining pilot data with additional observations to test for group differences
#' with guaranteed power.
#'
#' @param x1 Pilot data from stage 1 (list mode or matrix).
#' @param x2 Additional data collected in stage 2 (list mode or matrix). If NULL,
#'   only stage 1 sample size calculations are returned.
#' @inheritParams bdanova1
#'
#' @return A list with components:
#'   \item{test.stat}{The test statistic (if x2 is provided).}
#'   \item{crit}{The critical value.}
#'   \item{N}{Required sample sizes (if x2 is NULL).}
#'   \item{d}{Effect size parameter (if x2 is NULL).}
#'
#' @details
#' If \code{x2} is NULL, the function returns the required sample sizes from
#' stage 1 (same as \code{bdanova1}).
#'
#' If \code{x2} is provided, the function combines the pilot data (x1) with the
#' new observations (x2) and performs the test. The function checks that the
#' sample sizes in x2 meet the requirements determined in stage 1.
#'
#' @seealso \code{\link{bdanova1}}, \code{\link{anova1}}
#'
#' @export
#' @examples
#' # Stage 1: Pilot data
#' set.seed(123)
#' x1_pilot <- list(rnorm(10), rnorm(10, 0.5), rnorm(10, 1))
#' stage1 <- bdanova2(x1_pilot, delta=0.5)
#'
#' # Stage 2: Collect additional data and test
#' # x2_additional <- list(rnorm(stage1$N[1]-10), rnorm(stage1$N[2]-10), ...)
#' # bdanova2(x1_pilot, x2_additional, delta=0.5)
bdanova2<-function(x1,x2=NULL,alpha=.05,power=.9,delta){
#
#  Do the second stage of the Bishop-Duewicz ANOVA
#
if(is.null(x2[1])){
stage1=bdanova1(x1,alpha=alpha,power=power,delta=delta)
return(list(N=stage1$N,d=stage1$d,crit=stage1$crit))
}
if(!is.null(x2[1])){
if(is.na(delta))stop("A value for delta was not specified")
if(!is.list(x1)){
if(!is.matrix(x1))stop("Data must be stored in a matrix or in list mode")
y<-x1
x1<-list()
for(j in 1:ncol(y))x1[[j]]<-y[,j]
}
if(is.na(delta))stop("A value for delta was not specified")
if(!is.list(x2)){
if(!is.matrix(x2))stop("Data must be stored in matrix or in list mode")
y<-x2
x2<-list()
for(j in 1:ncol(y))x2[[j]]<-y[,j]
}
if(length(x1)!=length(x2))stop("Length of x1 does not match the length of x2")
TT<-NA
U<-NA
J<-length(x1)
nvec<-NA
nvec2<-NA
svec<-NA
for(j in 1:length(x1)){
nvec[j]<-length(x1[[j]])
nvec2[j]<-length(x2[[j]])
svec[j]<-var(x1[[j]])
TT[j]<-sum(x1[[j]])
U[j]<-sum(x2[[j]])
}
temp<-bdanova1(x1,alpha=alpha,power=power,delta=delta)
need<-temp$N-nvec
#for(j in 1:length(x1))print(c(nvec2[j],need[j]))
for(j in 1:length(x1))if(nvec2[j]<need[j]){
print(paste("Warning: For Group", j))
print("The first stage analysis based on bdanova1 reports that a larger")
print(" sample is required than what was found in the argument x2")
}
b<-sqrt(nvec*((nvec+nvec2)*temp$d-svec)/(nvec2*svec))
b<-(b+1)/(nvec+nvec2)
xtil<-TT*(1-nvec2*b)/nvec+b*U
ftil<-sum((xtil-mean(xtil))^2)/temp$d
return(list(test.stat=ftil,crit=temp$crit))
}
}

#' Heteroscedastic One-Way ANOVA for Trimmed Means (Welch-Type)
#'
#' Performs a heteroscedastic one-way ANOVA for trimmed means using a
#' generalization of Welch's method. Does not assume equal variances across
#' groups. This is the primary robust ANOVA function in the WRS package.
#'
#' @param x Data in matrix form (columns = groups), list mode, or a vector
#'   (when used with \code{IV}).
#' @param MAT Logical. If TRUE, \code{x} is a matrix with group indicators in
#'   one column and data in another (default: FALSE).
#' @param lev.col Column number indicating group membership when MAT=TRUE
#'   (default: 1).
#' @param var.col Column number containing the data values when MAT=TRUE
#'   (default: 2).
#' @param IV Vector of group identifiers. When specified, \code{x} should be
#'   a vector of all data values.
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{TEST}{The test statistic.}
#'   \item{nu1}{Numerator degrees of freedom.}
#'   \item{nu2}{Denominator degrees of freedom.}
#'   \item{n}{Vector of sample sizes for each group.}
#'   \item{p.value}{The p-value for the test.}
#'
#' @details
#' This function implements a heteroscedastic one-way ANOVA for trimmed means,
#' generalizing Welch's test to trimmed means. It does not assume equal variances
#' and uses adjusted degrees of freedom.
#'
#' The default 20% trimming provides good power while maintaining robustness.
#' Setting tr=0 performs Welch's test on means (though this is not recommended
#' for typical use).
#'
#' WARNING: Do not use this function to compare medians (tr=0.5); use
#' \code{\link{MEDanova}} or \code{\link{Qanova}} instead.
#'
#' @seealso \code{\link{t1wayv2}}, \code{\link{t1waybt}}, \code{\link{t1wayF}},
#'   \code{\link{MEDanova}}, \code{\link{pbanova}}
#'
#' @export
#' @examples
#' # Three groups in list mode
#' set.seed(123)
#' x1 <- rnorm(20, mean=0, sd=1)
#' x2 <- rnorm(20, mean=0.5, sd=2)
#' x3 <- rnorm(20, mean=1, sd=1.5)
#' t1way(list(x1, x2, x3), tr=0.2)
#'
#' # Using a data frame with grouping variable
#' df <- data.frame(score=c(x1,x2,x3), group=rep(1:3, each=20))
#' t1way(df$score, IV=df$group)
t1way<-function(x,tr=.2,grp=NA,MAT=FALSE,lev.col=1,var.col=2,IV=NULL,pr=TRUE){
#
#  A heteroscedastic one-way ANOVA for trimmed means
#  using a generalization of Welch's method.
#
#  The data are assumed to be stored in $x$ in a matrix or in list mode.
#
# MAT=F, if x is a matrix, columns correspond to groups.
# if MAT=T, assumes argument
# lev.col
# indicates which column of x denotes the groups. And
#  var.col indicates the column where the data are stored.
#
# if x has list mode:
#  length(x) is assumed to correspond to the total number of groups.
#  By default, the null hypothesis is that all groups have a common mean.
#  To compare a subset of the groups, use grp to indicate which
#  groups are to be compared. For example, if you type the
#  command grp<-c(1,3,4), and then execute this function, groups
#  1, 3, and 4 will be compared with the remaining groups ignored.
#
#  IV, if specified, taken to be the independent variable
#      That is, the group id values
#      and x is assumed to be a vector containing all of the data
#
#  Missing values are automatically removed.
#
if(is.data.frame(x))x=as.matrix(x)
if(tr==.5)print("Warning: Comparing medians should not be done with this function")
if(!is.null(IV[1])){
if(pr)print("Assuming x is a vector containing all of the data, the dependent variable")
xi=elimna(cbind(x,IV))
x=fac2list(xi[,1],xi[,2])
}
if(MAT){
if(!is.matrix(x))stop("With MAT=T, data must be stored in a matrix")
if(length(lev.col)!=1)stop("Argument lev.col should have 1 value")
temp=selby(x,lev.col,var.col)
x=temp$x
grp2=rank(temp$grpn)
x=x[grp2]
}
if(is.matrix(x))x<-listm(x)
#nv=lapply(x,length)
if(is.na(sum(grp[1])))grp<-c(1:length(x))
if(!is.list(x))stop("Data are not stored in a matrix or in list mode.")
J<-length(grp)
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
nv=NA
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
nv[j]=length(x[[j]])
h[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
   # h is the number of observations in the jth group after trimming.
if(winvar(x[[grp[j]]],tr)==0)print(paste('The Winsorized variance is zero for group',j))
w[j]<-h[j]*(h[j]-1)/((length(x[[grp[j]]])-1)*winvar(x[[grp[j]]],tr))
xbar[j]<-mean(x[[grp[j]]],tr)
}
u<-sum(w)
xtil<-sum(w*xbar)/u
A<-sum(w*(xbar-xtil)^2)/(J-1)
B<-2*(J-2)*sum((1-w/u)^2/(h-1))/(J^2-1)
TEST<-A/(B+1)
nu1<-J-1
nu2<-1./(3*sum((1-w/u)^2/(h-1))/(J^2-1))
sig<-1-pf(TEST,nu1,nu2)
list(TEST=TEST,nu1=nu1,nu2=nu2,n=nv,p.value=sig)
}

#' Heteroscedastic Three-Way ANOVA for Trimmed Means (Version 2 with P-Values)
#'
#' Performs a J by K by L three-way ANOVA on trimmed means for independent
#' groups. Does not assume equal variances. Tests main effects and all
#' interactions. This version computes p-values (unlike \code{t3way} which
#' only provides critical values).
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param L Number of levels for Factor C.
#' @param x Data in list mode (length J*K*L), matrix (columns = groups), or
#'   matrix in long format (when MAT=TRUE).
#' @param p Total number of groups (default: J*K*L).
#' @param MAT Logical. If TRUE, \code{x} is a matrix with factor levels in
#'   columns specified by \code{lev.col} and data in \code{var.col}.
#' @param lev.col Column numbers for the three factors when MAT=TRUE (default: c(1,2,3)).
#' @param var.col Column number for the outcome variable when MAT=TRUE (default: 4).
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{Qa}{Test statistic for main effect of Factor A.}
#'   \item{Qa.p.value}{P-value for Factor A.}
#'   \item{Qb}{Test statistic for main effect of Factor B.}
#'   \item{Qb.p.value}{P-value for Factor B.}
#'   \item{Qc}{Test statistic for main effect of Factor C.}
#'   \item{Qc.p.value}{P-value for Factor C.}
#'   \item{Qab}{Test statistic for A by B interaction.}
#'   \item{Qab.p.value}{P-value for A by B interaction.}
#'   \item{Qac}{Test statistic for A by C interaction.}
#'   \item{Qac.p.value}{P-value for A by C interaction.}
#'   \item{Qbc}{Test statistic for B by C interaction.}
#'   \item{Qbc.p.value}{P-value for B by C interaction.}
#'   \item{Qabc}{Test statistic for A by B by C interaction.}
#'   \item{Qabc.p.value}{P-value for three-way interaction.}
#'
#' @details
#' This function performs a heteroscedastic three-way ANOVA using trimmed means.
#' It does not assume equal variances across groups.
#'
#' Data ordering in list mode: Groups are ordered with the third factor (C)
#' varying fastest, then B, then A. For example:
#' \itemize{
#'   \item x[[1]] = (1,1,1), x[[2]] = (1,1,2), ..., x[[L]] = (1,1,L)
#'   \item x[[L+1]] = (1,2,1), ..., x[[2*L]] = (1,2,L)
#'   \item x[[K*L+1]] = (2,1,1), etc.
#' }
#'
#' @seealso \code{\link{t3way}}, \code{\link{t2way}}, \code{\link{t2wayv2}}
#'
#' @export
#' @examples
#' # 2x2x2 design with 8 groups
#' set.seed(123)
#' x <- vector("list", 8)
#' for(i in 1:8) x[[i]] <- rnorm(15, mean=i*0.2)
#' t3wayv2(2, 2, 2, x, tr=0.2)
t3wayv2<-function(J,K,L,x,tr=.2,grp=c(1:p),alpha=.05,p=J*K*L,MAT=FALSE,
lev.col=c(1:3),var.col=4,pr=TRUE){
#  Perform a J by K by L (three-way) anova on trimmed means where
#  all JKL groups are independent.
#
#  Same as t3way, only computes p-values
#
# if MAT=F (default)
#  The R variable data is assumed to contain the raw
#  data stored in list mode. data[[1]] contains the data
#  for the first level of all three factors: level 1,1,1.
#  data][2]] is assumed to contain the data for level 1 of the
#  first two factors and level 2 of the third factor: level 1,1,2
#  data[[L]] is the data for level 1,1,L
#  data[[L+1]] is the data for level 1,2,1. data[[2L]] is level 1,2,L.
#  data[[KL+1]] is level 2,1,1, etc.
#
#  MAT=T, assumes data are stored in matrix with 3 columns indicating
#  levels of the three factors.
#  That is, this function calls selby2 for you.
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that data has length JKL, the total number of
#  groups being tested.
#
if(is.data.frame(x))x=as.matrix(x)
data=x  #Yes, odd code
if(MAT){
if(!is.matrix(data))stop("With MAT=T, data must be a matrix")
if(length(lev.col)!=3)stop("Argument lev.col should have 3 values")
temp=selby2(data,lev.col,var.col)
lev1=length(unique(temp$grpn[,1]))
lev2=length(unique(temp$grpn[,2]))
lev3=length(unique(temp$grpn[,3]))
gv=apply(temp$grpn,2,rank)
gvad=100*gv[,1]+10*gv[,2]+gv[,3]
grp=rank(gvad)
if(pr){
print(paste("Factor 1 has", lev1, "levels"))
print(paste("Factor 2 has", lev2, "levels"))
print(paste("Factor 3 has", lev3, "levels"))
}
if(J!=lev1)warning("J is being reset to the number of levels found")
if(K!=lev2)warning("K is being reset to the number of levels found")
if(L!=lev3)warning("K is being reset to the number of levels found")
J=lev1
K=lev2
L=lev2
data=temp$x
}
if(is.matrix(data))data=listm(data)
if(!is.list(data))stop("Data is not stored in list mode")
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
v[i]<-(length(data[[grp[i]]])-1)*winvar(data[[grp[i]]],tr)/(h[i]*(h[i]-1))
#    v contains the squared standard errors
}
v<-diag(v,p,p)   # Put squared standard errors in a diag matrix.
ij<-matrix(c(rep(1,J)),1,J)
ik<-matrix(c(rep(1,K)),1,K)
il<-matrix(c(rep(1,L)),1,L)
jm1<-J-1
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
Qa <- johan(cmat, tmeans, v, h, alpha)
Qa.pv=t3pval(cmat, tmeans, v, h)
# Do test for factor B
cmat<-kron(ij,kron(ck,il))  # Contrast matrix for factor B
Qb<-johan(cmat,tmeans,v,h,alpha)
Qb.pv=t3pval(cmat, tmeans, v, h)
# Do test for factor C
cmat<-kron(ij,kron(ik,cl))  # Contrast matrix for factor C
Qc<-johan(cmat,tmeans,v,h,alpha)
Qc.pv=t3pval(cmat, tmeans, v, h)
# Do test for factor A by B interaction
cmat<-kron(cj,kron(ck,il))  # Contrast matrix for factor A by B
Qab<-johan(cmat,tmeans,v,h,alpha)
Qab.pv=t3pval(cmat, tmeans, v, h)
# Do test for factor A by C interaction
cmat<-kron(cj,kron(ik,cl))  # Contrast matrix for factor A by C
Qac<-johan(cmat,tmeans,v,h,alpha)
Qac.pv=t3pval(cmat, tmeans, v, h)
# Do test for factor B by C interaction
cmat<-kron(ij,kron(ck,cl))  # Contrast matrix for factor B by C
Qbc<-johan(cmat,tmeans,v,h,alpha)
Qbc.pv=t3pval(cmat, tmeans, v, h)
# Do test for factor A by B by C interaction
cmat<-kron(cj,kron(ck,cl))  # Contrast matrix for factor A by B by C
Qabc<-johan(cmat,tmeans,v,h,alpha)
Qabc.pv=t3pval(cmat, tmeans, v, h)
list(Qa=Qa$teststat,Qa.crit=Qa$crit,Qa.p.value=Qa.pv,
Qb=Qb$teststat,Qb.crit=Qb$crit,Qb.p.value=Qb.pv,
Qc=Qc$teststat,Qc.crit=Qc$crit,Qc.p.value=Qc.pv,
Qab=Qab$teststat,Qab.crit=Qab$crit,Qab.p.value=Qab.pv,
Qac=Qac$teststat,Qac.crit=Qac$crit,Qac.p.value=Qac.pv,
Qbc=Qbc$teststat,Qbc.crit=Qbc$crit,Qbc.p.value=Qbc.pv,
Qabc=Qabc$teststat,Qabc.crit=Qabc$crit,Qabc.p.value=Qabc.pv)
}

#' Heteroscedastic Two-Way ANOVA for Trimmed Means
#'
#' Performs a J by K two-way ANOVA on trimmed means for independent groups.
#' Does not assume equal variances. Tests main effects and interactions.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in list mode (length J*K), matrix (columns = groups), or
#'   vector (when used with IV1 and IV2).
#' @param p Total number of groups (default: J*K).
#' @param MAT Logical. If TRUE, \code{x} is a matrix with factor levels in
#'   columns specified by \code{lev.col} and data in \code{var.col}.
#' @param lev.col Column numbers for the two factors when MAT=TRUE (default: c(1,2)).
#' @param var.col Column number for the outcome variable when MAT=TRUE (default: 3).
#' @param IV1 Vector of first factor levels (when x is a vector of data).
#' @param IV2 Vector of second factor levels (when x is a vector of data).
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{Qa}{Test statistic for main effect of Factor A.}
#'   \item{Qa.crit}{Critical value for Factor A.}
#'   \item{Qa.p.value}{P-value for Factor A.}
#'   \item{Qb}{Test statistic for main effect of Factor B.}
#'   \item{Qb.crit}{Critical value for Factor B.}
#'   \item{Qb.p.value}{P-value for Factor B.}
#'   \item{Qab}{Test statistic for A by B interaction.}
#'   \item{Qab.crit}{Critical value for interaction.}
#'   \item{Qab.p.value}{P-value for interaction.}
#'
#' @details
#' This function performs a heteroscedastic two-way ANOVA using trimmed means.
#' It uses the Johansen procedure and does not assume equal variances across
#' groups.
#'
#' Data can be provided in several formats:
#' \itemize{
#'   \item List mode: x[[1]] = group (1,1), x[[2]] = group (1,2), ..., x[[K]] = group (1,K), x[[K+1]] = group (2,1), etc.
#'   \item Matrix: Columns correspond to groups in the same order
#'   \item Long format: Use MAT=TRUE with factor levels in lev.col and data in var.col
#'   \item Vectors: Provide data in x with factor levels in IV1 and IV2
#' }
#'
#' @seealso \code{\link{t2wayv2}}, \code{\link{t3way}}, \code{\link{t2waybt}},
#'   \code{\link{bwtrim}}, \code{\link{sppbb}}
#'
#' @export
#' @examples
#' # 2x3 design with 6 groups
#' set.seed(123)
#' x <- vector("list", 6)
#' for(i in 1:6) x[[i]] <- rnorm(15, mean=i*0.3)
#' t2way(2, 3, x, tr=0.2)
#'
#' # Using IV vectors
#' n <- 30
#' data <- rnorm(n)
#' factor1 <- rep(1:2, each=15)
#' factor2 <- rep(1:3, times=10)
#' t2way(2, 3, data, IV1=factor1, IV2=factor2)
t2way<-function(J,K,x,tr=.2,grp=c(1:p),p=J*K,MAT=FALSE,
lev.col=c(1:2),var.col=3,pr=TRUE,IV1=NULL,IV2=NULL){
#  Perform a J by K  (two-way) ANOVA on trimmed means where
#  all groups are independent.
#
#  The R variable x is assumed to contain the raw
#  data stored in list mode, or a matrix with columns
#  corresponding to groups. If stored in list mode, x[[1]] contains the data
#  for the first level of both factors: level 1,1,.
#  x[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second factor: level 1,2
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that x has length JK, the total number of
#  groups being tested.
#
#  MAT=T, assumes x are stored in matrix with 3 columns
#  with two of the columns indicated by the argument
#  lev.col
#  specifying the columns of x containing the values of the
#  levels of the two factors.
#  The outcome variable is in column
#  var.col
#  which defaults to column 3
#  That is, this function calls selby2 for you.
#
#  IV1 and IV2: if specified, taken to be the independent variable
#      That is, the group id values
#      and x is assumed to be a vector containing all of the data
#  EXAMPLE: t2way(x=data,IV1=iv1,IV2=iv2)
#  would do a two-way ANOVA based on group id's in iv1 and iv2 and
#  dependent variable data
#
if(is.data.frame(x))data=as.matrix(x)
if(tr==.5){
print("For medians, use med2way if there are no ties")
print("With ties, use linear contrasts in conjunction with medpb")
stop("")
}
if(MAT){
if(!is.matrix(x))stop("With MAT=T, data must be a matrix")
if(length(lev.col)!=2)stop("Argument lev.col should have 3 values")
temp=selby2(x,lev.col,var.col)
lev1=length(unique(temp$grpn[,1]))
lev2=length(unique(temp$grpn[,2]))
gv=apply(temp$grpn,2,rank)
gvad=10*gv[,1]+gv[,2]
grp=rank(gvad)
if(pr){
print(paste("Factor 1 has", lev1, "levels"))
print(paste("Factor 2 has", lev2, "levels"))
}
if(J!=lev1)warning("J is being reset to the number of levels found")
if(K!=lev2)warning("K is being reset to the number of levels found")
J=lev1
K=lev2
x=temp$x
}
if(!is.null(IV1[1])){
if(is.null(IV2[1]))stop("IV2 is NULL")
if(pr)print("Assuming data is a vector containing all of the data; the dependent variable")
xi=elimna(cbind(x,IV1,IV2))
J=length(unique(xi[,2]))
K=length(unique(xi[,3]))
x=fac2list(xi[,1],xi[,2:3])
}
if(is.matrix(x))x=listm(x)
if(!is.list(x))stop("Data are not stored in list mode")
if(p!=length(x)){
print("The total number of groups, based on the specified levels, is")
print(p)
print("The number of groups is")
print(length(x))
print("Warning: These two values are not equal")
}
tmeans<-0
h<-0
v<-0
for (i in 1:p){
x[[grp[i]]]=elimna(x[[grp[i]]])
tmeans[i]<-mean(x[[grp[i]]],tr)
h[i]<-length(x[[grp[i]]])-2*floor(tr*length(x[[grp[i]]]))
#    h is the effective sample size
if(winvar(x[[grp[i]]],tr)==0)print(paste('The Winsorized variance is zero for group',i))
v[i]<-(length(x[[grp[i]]])-1)*winvar(x[[grp[i]]],tr)/(h[i]*(h[i]-1))
#    v contains the squared standard errors
}
v<-diag(v,p,p)   # Put squared standard errors in a diag matrix.
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
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
Qa<-johan(cmat,tmeans,v,h,alval[i])
if(i==1)dfA=Qa$df
if(Qa$teststat>Qa$crit)break
}
A.p.value=irem/1000
# Do test for factor B
cmat<-kron(ij,ck)  # Contrast matrix for factor B
for(i in 1:999){
irem<-i
Qb<-johan(cmat,tmeans,v,h,alval[i])
if(i==1)dfB=Qb$df
if(Qb$teststat>Qb$crit)break
}
B.p.value=irem/1000
# Do test for factor A by B interaction
cmat<-kron(cj,ck)  # Contrast matrix for factor A by B
for(i in 1:999){
irem<-i
Qab<-johan(cmat,tmeans,v,h,alval[i])
if(i==1)dfAB=Qab$df
if(Qab$teststat>Qab$crit)break
}
AB.p.value=irem/1000
tmeans=matrix(tmeans,J,K,byrow=TRUE)
list(Qa=Qa$teststat,A.p.value=A.p.value, df.A=dfA,
Qb=Qb$teststat,B.p.value=B.p.value,df.B=dfB,
Qab=Qab$teststat,AB.p.value=AB.p.value,df.AB=dfAB,means=tmeans)
}

#' Percentile Bootstrap One-Way ANOVA for Trimmed Means
#'
#' Tests the hypothesis that J independent groups have equal trimmed means
#' using a percentile bootstrap method. Provides robust alternative to
#' parametric ANOVA.
#'
#' @param win Amount of Winsorizing when WIN=TRUE (default: 0.1). Must be
#'   less than or equal to \code{tr}.
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{test}{The test statistic.}
#'   \item{crit.val}{Critical value from the bootstrap distribution.}
#'   \item{p.value}{Bootstrap p-value.}
#'   \item{con}{Contrast matrix used for testing.}
#'   \item{num.sig}{Number of significant contrasts.}
#'
#' @details
#' This function uses the percentile bootstrap to test for equality of trimmed
#' means across J independent groups. It generates bootstrap samples from each
#' group and computes the distribution of the test statistic under the null
#' hypothesis.
#'
#' The number of bootstrap samples is determined automatically if not specified:
#' 5000 for J >= 10, otherwise 2000.
#'
#' If WIN=TRUE, data are Winsorized before bootstrap resampling. This can
#' improve performance but should only be used with tr >= 0.2 and n >= 15.
#'
#' @seealso \code{\link{t1way}}, \code{\link{t1waybt}}, \code{\link{MEDanova}},
#'   \code{\link{pbanovag}}
#'
#' @export
#' @examples
#' # Three groups
#' set.seed(123)
#' x1 <- rnorm(25, mean=0, sd=1)
#' x2 <- rnorm(25, mean=0.5, sd=1.2)
#' x3 <- rnorm(25, mean=1, sd=0.8)
#' pbanova(list(x1, x2, x3), tr=0.2, nboot=1000)
pbanova<-function(x,tr=.2,alpha=.05,nboot=NA,grp=NA,WIN=FALSE,win=.1){
#
#   Test the hypothesis that J independent groups have
#   equal trimmed means using the percentile bootstrap method.
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   tr is the amount of trimming
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#   WIN=T means data are Winsorized before taking bootstraps by the
#   amount win.
#
#   Missing values are allowed.
#
if(is.matrix(x))x<-listm(x)
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
if(WIN){
if(tr < .2){print("Warning: When Winsorizing,")
print("the amount of trimming should be at least.2")
}
if(win > tr)stop("Amount of Winsorizing must be <= amount of trimming")
if(min(tempn) < 15){
print("Warning: Winsorizing with sample sizes less than 15")
print("can result in poor control over the probability of a Type I error")
}
for (j in 1:J){
x[[j]]<-winval(x[[j]],win)
}
}
con<-matrix(0,J,J-1)
for (j in 1:Jm){
jp<-j+1
con[j,j]<-1
con[jp,j]<-0-1
}
#  Determine nboot if a value was not specified
if(is.na(nboot)){
nboot<-5000
if(J <= 8)nboot<-4000
if(J <= 3)nboot<-2000
}
# Determine critical values
if(alpha==.05){
dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(Jm > 10){
avec<-.05/c(11:Jm)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(Jm > 10){
avec<-.01/c(11:Jm)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:Jm)
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
for (d in 1:Jm){
dp<-d+1
test[d]<-sum(bvec[d,]>bvec[dp,])/nboot
if(test[d]> .5)test[d]<-1-test[d]
}
test<-(0-1)*sort(-2*test)
sig<-sum((test<dvec[1:Jm]))
if(sig>0)print("Significant result obtained: Reject")
if(sig==0)print("No significant result obtained: Fail to reject")
list(test.vec=test,crit.vec=dvec[1:Jm])
}

#' Generalized Percentile Bootstrap ANOVA for Any Estimator
#'
#' Tests the hypothesis that J independent groups have equal measures of
#' location (or scale) using the percentile bootstrap method. Generalizes
#' \code{pbanova} to work with any estimator function.
#'
#' @param est Estimator function to use. Default is \code{onestep} (one-step
#'   M-estimator). Other options include \code{mean}, \code{median}, \code{tmean},
#'   \code{mom}, \code{mest}, \code{mad}, \code{winvar}, etc.
#' @inheritParams common-params
#'
#' @return A list with components:
#'   \item{test}{The test statistic.}
#'   \item{crit.val}{Critical value from the bootstrap distribution.}
#'   \item{p.value}{Bootstrap p-value.}
#'   \item{con}{Contrast matrix used for testing.}
#'   \item{num.sig}{Number of significant contrasts.}
#'
#' @details
#' This function extends \code{pbanova} to work with any estimator function,
#' not just trimmed means. The estimator can be any function that takes a
#' vector and returns a scalar.
#'
#' Additional arguments to the estimator can be passed via \code{...}. For
#' example, to use a trimmed mean with 10\% trimming, use
#' \code{est=mean, trim=0.1}.
#'
#' The number of bootstrap samples is determined automatically if not specified:
#' 5000 for J >= 10, otherwise 2000.
#'
#' @seealso \code{\link{pbanova}}, \code{\link{t1way}}, \code{\link{MEDanova}}
#'
#' @export
#' @examples
#' # Compare groups using one-step M-estimator
#' set.seed(123)
#' x1 <- rnorm(25, mean=0, sd=1)
#' x2 <- rnorm(25, mean=0.5, sd=1.2)
#' x3 <- rnorm(25, mean=1, sd=0.8)
#' pbanovag(list(x1, x2, x3), est=onestep, nboot=1000)
#'
#' # Compare medians
#' pbanovag(list(x1, x2, x3), est=median, nboot=1000)
#'
#' # Compare MADs (median absolute deviation)
#' pbanovag(list(x1, x2, x3), est=mad, nboot=1000)
pbanovag<-function(x,alpha=.05,nboot=NA,grp=NA,est=onestep,...){
#
#   Test the hypothesis that J independent groups have
#   equal measures of location using the percentile bootstrap method.
#   (Robust measures of scale can be compared as well.)
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   est is the measure of location and defaults to a M-estimator
#   ... can be used to set optional arguments associated with est
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#   Missing values are allowed.
#
con<-as.matrix(con)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
if(!is.na(sum(grp))){
# Only analyze specified groups.
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
icl<-ceiling(crit*nboot)
icu<-ceiling((1-crit)*nboot)
con<-matrix(0,J,J-1)
for (j in 1:Jm){
jp<-j+1
con[j,j]<-1
con[jp,j]<-0-1
}
#  Determine nboot if a value was not specified
if(is.na(nboot)){
nboot<-5000
if(J <= 8)nboot<-4000
if(J <= 3)nboot<-2000
}
# Determine critical values
if(alpha==.05){
dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(Jm > 10){
avec<-.05/c(11:Jm)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(Jm > 10){
avec<-.01/c(11:Jm)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:Jm)
bvec<-matrix(NA,nrow=J,ncol=nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
for(j in 1:J){
paste("Working on group ",j)
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,est,...) # Bootstrapped trimmed means for jth group
}
test<-NA
for (d in 1:Jm){
dp<-d+1
test[d]<-sum(bvec[d,]>bvec[dp,])/nboot
if(test[d]> .5)test[d]<-1-test[d]
}
test<-(0-1)*sort(-2*test)
sig<-sum((test<dvec[1:Jm]))
if(sig>0)print("Significant result obtained: Reject")
if(sig==0)print("No significant result obtained: Fail to reject")
list(test.vec=test,crit.vec=dvec[1:Jm])
}

#' Heteroscedastic One-Way Random Effects ANOVA for Trimmed Means
#'
#' Performs a heteroscedastic one-way random effects ANOVA for trimmed means.
#' This function tests whether all groups have a common trimmed mean, allowing
#' for unequal variances across groups.
#'
#' @param x Data in matrix or list format. If matrix, groups correspond to columns.
#'   If list, each element contains data for one group.
#' @param tr Proportion of trimming (default is 0.2 for 20% trimming).
#' @param grp Vector specifying which groups to compare. If NA (default), all groups
#'   are compared. Use grp=c(1,3,4) to compare only groups 1, 3, and 4.
#'
#' @return A list with components:
#'   \item{teststat}{The test statistic D}
#'   \item{df}{Degrees of freedom (numerator and denominator)}
#'   \item{p.value}{p-value for the test}
#'   \item{rho}{Estimated Winsorized intraclass correlation}
#'   \item{num.groups}{Number of groups being compared}
#'
#' @details
#' The function uses a heteroscedastic approach that does not assume equal variances
#' across groups. It computes trimmed means and Winsorized variances for each group,
#' then performs an F-test with adjusted degrees of freedom.
#'
#' The null hypothesis is that all groups have a common trimmed mean. The test
#' statistic follows an F-distribution under the null hypothesis.
#'
#' @export
#' @examples
#' # Example with matrix data
#' x <- matrix(rnorm(60), ncol=3)
#' rananova(x, tr=0.2)
#'
#' # Example with list data
#' x <- list(rnorm(20), rnorm(25), rnorm(30))
#' rananova(x)
#'
#' # Compare only groups 1 and 3
#' rananova(x, grp=c(1,3))
rananova<-function(x,tr=.2,grp=NA){
#
#  A heteroscedastic one-way random effects ANOVA for trimmed means.
#
#  The data are assumed to be stored in a matrix on in list mode.
#  If in list mode,
#  Length(x) is assumed to correspond to the total number of groups.
#  If the data are stored in a matrix, groups correspond to columns.
#  By default, the null hypothesis is that all group have a common mean.
#  To compare a subset of the groups, use grp to indicate which
#  groups are to be compared. For example, if you type the
#  command grp<-c(1,3,4), and then execute this function, groups
#  1, 3, and 4 will be compared with the remaining groups ignored.
#
if(is.matrix(x))x<-listm(x)
if(is.na(grp[1]))grp<-c(1:length(x))
if(!is.list(x))stop("Data are not stored in a matrix or in list mode")
J<-length(grp)  # The number of groups to be compared
#if(pr)print("The number of groups to be compared is")
#print(J)
h<-1
xbar<-1
ybar<-1
wvar<-1
ell<-0
for(j in 1:J){
ell[j]<-length(x[[grp[j]]])/(length(x[[grp[j]]])+1)
h[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
   # h is the number of observations in the jth group after trimming.
ybar[j]<-winmean(x[[grp[j]]],tr)
xbar[j]<-mean(x[[grp[j]]],tr)
wvar[j]<-winvar(x[[grp[j]]],tr)
}
q<-NA
bsst<-var(xbar)
for (j in 1:J)q[j]<-(length(x[[grp[j]]]-1)-1)*wvar[j]/(h[j]*(h[j]-1))
wssw<-mean(q)
D<-bsst/wssw
g<-q/J
nu1<-((J-1)*sum(q))^2/((sum(q))^2+(J-2)*J*sum(q^2))
nu2<-(sum(J*q))^2/sum((J*q)^2/(h-1))
sig<-1-pf(D,nu1,nu2)
# Next, estimate the Winsorized intraclass correlation
sighat<-mean(ell*(ybar-(sum(ell*ybar)/sum(ell)))^2)
rho<-sighat/(sighat+winmean(wvar,tr))
list(teststat=D,df=c(nu1,nu2),p.value=sig,rho=rho,num.groups=J)
}

#' Two-Way ANOVA for Independent Groups Using Percentile Bootstrap
#'
#' Performs a two-way ANOVA for independent groups based on robust measures
#' of location and a percentile bootstrap method. Tests main effects and
#' interaction using bootstrap hypothesis testing.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in matrix or list format. If list, x[[1]] contains data for
#'   level (1,1), x[[2]] for level (1,2), etc. Groups are ordered by columns
#'   of the cross-classification table.
#' @param alpha Significance level (default is 0.05).
#' @param nboot Number of bootstrap samples. If NA, defaults to 5000 (or 2000
#'   if the number of contrasts is <= 4).
#' @param grp Optional vector to rearrange groups. For example, grp=c(2,4,3,1)
#'   for a 2x2 design indicates group 2 is (1,1), group 4 is (1,2), group 3 is
#'   (2,1), and group 1 is (2,2).
#' @param est Estimator to use (default is onestep M-estimator).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list with components:
#'   \item{testA}{Test statistics for Factor A main effect}
#'   \item{crit.vecA}{Critical values for Factor A}
#'   \item{testB}{Test statistics for Factor B main effect}
#'   \item{crit.vecB}{Critical values for Factor B}
#'   \item{testAB}{Test statistics for interaction}
#'   \item{crit.vecAB}{Critical values for interaction}
#'
#' @details
#' This function performs a heteroscedastic two-way ANOVA using a percentile
#' bootstrap approach. It tests the main effects of both factors and their
#' interaction. The default estimator is the one-step M-estimator, but other
#' robust estimators can be specified.
#'
#' Missing values are automatically removed. The random seed is set to 2 for
#' reproducibility. The function prints messages indicating significant results
#' for each effect.
#'
#' @export
#' @examples
#' # Create data for 2x3 design
#' set.seed(123)
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(20, mean=i)
#'
#' # Run two-way ANOVA
#' result <- t2waypbg(J=2, K=3, x=x, nboot=500)
t2waypbg<-function(J,K,x,alpha=.05,nboot=NA,grp=NA,est=onestep,...){
#
#   Two-way ANOVA for independent groups based on
#   robust measures of location
#   and a percentile bootstrap method.

#   The data are assumed to be stored in x in list mode or in a matrix.
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
#   Missing values are automatically removed.
#
if(is.data.frame(x))x=as.matrix(x)
JK<-J*K
if(is.matrix(x))x<-listm(x)
if(!is.na(grp)){
yy<-x
for(j in 1:length(grp))
x[[j]]<-yy[[grp[j]]]
}
if(!is.list(x))stop("Data must be stored in list mode or a matrix.")
for(j in 1:JK){
xx<-x[[j]]
x[[j]]<-xx[!is.na(xx)] # Remove any missing values.
}
#
# Create the three contrast matrices
#
       ij <- matrix(c(rep(1, J)), 1, J)
        ik <- matrix(c(rep(1, K)), 1, K)
       jm1 <- J - 1
        cj <- diag(1, jm1, J)
        for(i in 1:jm1)
                cj[i, i + 1] <- 0 - 1
        km1 <- K - 1
        ck <- diag(1, km1, K)
        for(i in 1:km1)
                ck[i, i + 1] <- 0 - 1
conA<-t(kron(cj,ik))
conB<-t(kron(ij,ck))
conAB<-t(kron(cj,ck))
ncon<-max(nrow(conA),nrow(conB),nrow(conAB))
if(JK!=length(x)){
print("Warning: The number of groups does not match")
print("the number of contrast coefficients.")
}
if(is.na(nboot)){
nboot<-5000
if(ncon<=4)nboot<-2000
}
m1<-matrix(0,nrow=JK,ncol=nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
print("Taking bootstrap samples. Please wait.")
for(j in 1:JK){
paste("Working on group ",j)
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
m1[j,]<-apply(data,1,est,...)
}
bootA<-matrix(0,ncol(conA),nboot)
bootB<-matrix(0,ncol(conB),nboot)
bootAB<-matrix(0,ncol(conAB),nboot)
testA<-NA
testB<-NA
testAB<-NA
testvecA<-NA
testvecB<-NA
testvecAB<-NA
for (d in 1:ncol(conA)){
bootA[d,]<-apply(m1,2,trimpartt,conA[,d])
# A vector of length nboot containing psi hat values
# corresponding to the dth linear contrast
testA[d]<-sum((bootA[d,]>0))/nboot
testA[d]<-min(testA[d],1-testA[d])
}
for (d in 1:ncol(conB)){
bootB[d,]<-apply(m1,2,trimpartt,conB[,d])
# A vector of length nboot containing psi hat values
# corresponding to the dth linear contrast
testB[d]<-sum((bootB[d,]>0))/nboot
testB[d]<-min(testB[d],1-testB[d])
}
for (d in 1:ncol(conAB)){
bootAB[d,]<-apply(m1,2,trimpartt,conAB[,d])
# A vector of length nboot containing psi hat values
# corresponding to the dth linear contrast
testAB[d]<-sum((bootAB[d,]>0))/nboot
testAB[d]<-min(testAB[d],1-testAB[d])
}
#
#  Determine critical value
#
Jm<-J-1
Km<-K-1
JKm<-(J-1)*(K-1)
dvecA <- alpha/c(1:Jm)
dvecB <- alpha/c(1:Km)
dvecAB <- alpha/c(1:JKm)
testA<-(0 - 1) * sort(-2 * testA)
testB<-(0 - 1) * sort(-2 * testB)
testAB<-(0 - 1) * sort(-2 * testAB)
sig <- sum((testA < dvecA[1:Jm]))
if(sig > 0)
print("Significant result obtained for Factor A: Reject")
if(sig == 0)
print("No significant result Factor A: Fail to reject")
sig <- sum((testB < dvecB[1:Km]))
if(sig > 0)
print("Significant result obtained for Factor B: Reject")
if(sig == 0)
print("No significant result Factor B: Fail to reject")
sig <- sum((testAB < dvecAB[1:JKm]))
if(sig > 0)
print("Significant Interaction: Reject")
if(sig == 0)
print("No significant Interaction: Fail to reject")
list(testA=testA,crit.vecA=dvecA,testB=testB,crit.vecB=dvecB,testAB=testAB,crit.vecAB=dvecAB)
}

#' Brunner-Dette-Munk Rank-Based One-Way ANOVA
#'
#' Performs the Brunner, Dette, and Munk rank-based ANOVA for comparing
#' independent groups. This is a heteroscedastic rank-based method that does
#' not assume equal variances.
#'
#' @param x Data in matrix or list format. If matrix, groups correspond to columns.
#'   If list, each element contains data for one group.
#' @param grp Vector specifying which groups to compare. If NA (default), all groups
#'   are compared.
#'
#' @return Output from bdms1 function containing test statistics and p-values.
#'
#' @details
#' Implements the Brunner, Dette, and Munk (1997) rank-based ANOVA method.
#' This approach is robust to heteroscedasticity and does not require normality
#' assumptions. The method uses a multivariate rank transformation approach.
#'
#' @references
#' Brunner, E., Dette, H., & Munk, A. (1997). Box-type approximations in
#' nonparametric factorial designs. Journal of the American Statistical
#' Association, 92, 1494-1502.
#'
#' @export
#' @examples
#' # Example with matrix data
#' x <- matrix(rnorm(60), ncol=3)
#' bdm(x)
#'
#' # Example with list data
#' x <- list(rnorm(20), rnorm(25), rnorm(30))
#' bdm(x)
bdm<-function(x,grp=NA){
#
# Perform the Brunner, Dette, Munk rank-based ANOVA
# (JASA, 1997, 92, 1494--1502)
#
# x can be a matrix with columns corresponding to groups
# or it can have list mode.
#
if(is.matrix(x))x<-listm(x)
J<-length(x)
xx<-list()
if(is.na(grp[1]))grp<-c(1:J)
for(j in 1:J)xx[[j]]<-x[[grp[j]]]
Ja<-matrix(1,J,J)
Ia<-diag(1,J)
Pa<-Ia-Ja/J
cona<-Pa
outA<-bdms1(xx,cona)
outA
}

#' Choi-Marden Multivariate One-Way Rank-Based ANOVA
#'
#' Performs the Choi and Marden multivariate one-way rank-based ANOVA for
#' comparing groups on multiple dependent measures. This method handles
#' J independent groups with K dependent measures per group.
#'
#' @param J Number of independent groups.
#' @param K Number of dependent measures.
#' @param x Data in matrix or list format. Data should contain J*K groups
#'   representing the J groups measured on K dependent variables.
#' @param grp Vector specifying group ordering (default is 1:(J*K)).
#' @param JK Total number of groups (default is J*K).
#'
#' @return A list with components:
#'   \item{test.stat}{The Choi-Marden test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p.value}{p-value for the test}
#'
#' @details
#' Implements the Choi and Marden (1997) multivariate rank-based ANOVA.
#' This method is designed for multivariate data where each group is measured
#' on multiple dependent variables. The test uses spatial ranks and does not
#' require multivariate normality.
#'
#' Missing values are automatically removed. The test statistic follows an
#' approximate chi-squared distribution with K*(J-1) degrees of freedom.
#'
#' @references
#' Choi, K., & Marden, J. (1997). An approach to multivariate rank tests in
#' multivariate analysis of variance. Journal of the American Statistical
#' Association, 92, 1581-1590.
#'
#' @export
#' @examples
#' # Example: 3 groups, 2 dependent measures
#' set.seed(123)
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(20)
#' cmanova(J=3, K=2, x=x)
cmanova<-function(J,K,x,grp=c(1:JK),JK=J*K){
#
# Perform the Choi and Marden
# multivariate one-way rank-based ANOVA
# (Choi and Marden, JASA, 1997, 92, 1581-1590.
#
# x can be a matrix with columns corresponding to groups
# or it can have list mode.
#
# Have a J by K design with J independent levels and K dependent
# measures
#
#
x=elimna(x)
if(is.matrix(x))x<-listm(x)
xx<-list()
nvec<-NA
jk<-0
for(j in 1:J){
for(k in 1:K){
jk<-jk+1
xx[[jk]]<-x[[grp[jk]]]
if(k==1)nvec[j]<-length(xx[[jk]])
}}
N<-sum(nvec)
RVALL<-matrix(0,nrow=N,K)
x<-xx
jk<-0
rmean<-matrix(NA,nrow=J,ncol=K)
for(j in 1:J){
RV<-matrix(0,nrow=nvec[j],ncol=K)
jk<-jk+1
temp1<-matrix(x[[jk]],ncol=1)
for(k in 2:K){
jk<-jk+1
temp1<-cbind(temp1,x[[jk]])
}
X<-temp1
if(j==1)XALL<-X
if(j>1)XALL<-rbind(XALL,X)
n<-nvec[j]
for(i in 1:n){
for (ii in 1:n){
temp3<-sqrt(sum((X[i,]-X[ii,])^2))
if(temp3 != 0)RV[i,]<-RV[i,]+(X[i,]-X[ii,])/temp3
}
RV[i,]<-RV[i,]/nvec[j]
if(j==1 && i==1)sighat<-RV[i,]%*%t(RV[i,])
if(j>1 || i>1)sighat<-sighat+RV[i,]%*%t(RV[i,])
}
}
# Assign ranks to pooled data and compute R bar for each group
for(i in 1:N){
for (ii in 1:N){
temp3<-sqrt(sum((XALL[i,]-XALL[ii,])^2))
if(temp3 != 0)RVALL[i,]<-RVALL[i,]+(XALL[i,]-XALL[ii,])/temp3
}
RVALL[i,]<-RVALL[i,]/N
}
bot<-1-nvec[1]
top<-0
for(j in 1:J){
bot<-bot+nvec[j]
top<-top+nvec[j]
flag<-c(bot:top)
rmean[j,]<-apply(RVALL[flag,],2,mean)
}
sighat<-sighat/(N-J)
shatinv<-solve(sighat)
KW<-0
for(j in 1:J){
KW<-KW+nvec[j]*t(rmean[j,])%*%shatinv%*%rmean[j,]
}
df<-K*(J-1)
sig.level<-1-pchisq(KW,df)
list(test.stat=KW[1,1],df=df,p.value=sig.level)
}

#' Split-Plot Bootstrap Test for Main Effects in Dependent Groups
#'
#' Performs a percentile bootstrap test for main effects among dependent groups
#' in a split-plot design. Tests whether all pairwise differences have a
#' typical value of zero.
#'
#' @param J Number of levels for between-subjects factor.
#' @param K Number of levels for within-subjects factor.
#' @param x Data in matrix or list format. If list, x[[1]] contains data for
#'   level (1,1), x[[2]] for level (1,2), ..., x[[K]] for level (1,K),
#'   x[[K+1]] for level (2,1), etc. If matrix, columns correspond to groups.
#' @param est Estimator to use (default is tmean for trimmed mean).
#' @param JK Total number of groups (default is J*K).
#' @param grp Vector specifying group ordering (default is 1:JK).
#' @param nboot Number of bootstrap samples (default is 500).
#' @param SEED Logical indicating whether to set random seed (default is TRUE).
#' @param pr Logical indicating whether to print messages (default is TRUE).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list with components:
#'   \item{p.value}{p-values for the tests}
#'   \item{center}{Estimated centers for the groups}
#'
#' @details
#' This function analyzes split-plot designs where subjects are measured under
#' multiple conditions (within-subjects factor) and there may be different groups
#' of subjects (between-subjects factor). The analysis is based on all pairs of
#' difference scores.
#'
#' Missing values are automatically removed. The random seed is set to 2 if
#' SEED=TRUE for reproducibility.
#'
#' @export
#' @examples
#' # Create data for 2x3 split-plot design
#' set.seed(123)
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(20)
#'
#' # Run split-plot bootstrap test
#' result <- sppbb(J=2, K=3, x=x, nboot=100)
sppbb<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),nboot=500,SEED=TRUE,pr=TRUE,...){
#
# A percentile bootstrap for main effects
# among dependent groups in a split-plot design
# The analysis is done based on all pairs
# of difference scores. The null hypothesis is that
# all such differences have a typical value of zero.
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
if(pr)print('As of Oct, 2014, the argument est defaults to tmean')
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}

JK<-J*K
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
jp<-1-K
kv<-0
kv2<-0
for(j in 1:J){
jp<-jp+K
xmat<-matrix(NA,ncol=K,nrow=length(x[[jp]]))
for(k in 1:K){
kv<-kv+1
xmat[,k]<-x[[kv]]
}
xmat<-elimna(xmat)
for(k in 1:K){
kv2<-kv2+1
x[[kv2]]<-xmat[,k]
}}
xx<-x
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
# Next determine the n_j values
nvec<-NA
jp<-1-K
for(j in 1:J){
jp<-jp+K
nvec[j]<-length(x[[jp]])
}
#
# Now stack the data in an N by K matrix
#
x<-matrix(NA,nrow=nvec[1],ncol=K)
#
for(k in 1:K)x[,k]<-xx[[k]]
kc<-K
for(j in 2:J){
temp<-matrix(NA,nrow=nvec[j],ncol=K)
for(k in 1:K){
kc<-kc+1
temp[,k]<-xx[[kc]]
}
x<-rbind(x,temp)
}
# Now call function rmdzero to do the analysis
temp<-rmdzero(x,est=est,nboot=nboot,...)
list(p.value=temp$p.value,center=temp$center)
}

#' Two-Way ANOVA on Trimmed Means Without p-values
#'
#' Performs a J by K two-way ANOVA on trimmed means where all JK groups are
#' independent. Returns test statistics and critical values but does not
#' compute p-values for all tests (interaction uses critical value instead).
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in matrix or list format. Groups are ordered by columns of
#'   the cross-classification table. x[[1]] is level (1,1), x[[2]] is (1,2), etc.
#' @param tr Proportion of trimming (default is 0.2 for 20% trimming).
#' @param grp Vector to rearrange group order if needed (default is 1:p).
#' @param alpha Significance level for interaction critical value (default is 0.05).
#' @param p Total number of groups (default is J*K).
#'
#' @return A list with components:
#'   \item{Qa}{Test statistic for Factor A main effect}
#'   \item{sig.A}{p-value for Factor A}
#'   \item{Qb}{Test statistic for Factor B main effect}
#'   \item{sig.B}{p-value for Factor B}
#'   \item{Qab}{Test statistic for interaction}
#'   \item{critinter}{Adjusted critical value for interaction}
#'
#' @details
#' This function performs a heteroscedastic two-way ANOVA using trimmed means.
#' Main effects are tested using F-tests with adjusted degrees of freedom.
#' The interaction is tested using a chi-squared statistic with an adjusted
#' critical value rather than a p-value.
#'
#' Missing values are automatically removed. The function uses Winsorized
#' variances and does not assume equal variances across groups.
#'
#' @export
#' @examples
#' # Create data for 2x3 design
#' set.seed(123)
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(20, mean=i)
#'
#' # Run two-way ANOVA
#' result <- t2way.no.p(J=2, K=3, x=x)
t2way.no.p<-function(J,K,x,tr=.2,grp=c(1:p),alpha=.05,p=J*K){
#  Perform a J by K (two-way) anova on trimmed means where
#  all jk groups are independent.
#
#  The R variable x is assumed to contain the raw
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
#  The default amount of trimming is tr=.2
#
#  It is assumed that the input variable x has length JK, the total number of
#  groups being tested. If not, a warning message is printed.
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data are not stored in a matrix or in list mode")
if(p!=length(x)){
print("Warning: The number of groups in your data is not equal to JK")
}
for(j in 1:p)x[[j]]<-elimna(x[[j]])
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
xbar[j]<-mean(x[[grp[j]]],tr)
h[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
d[j]<-(length(x[[grp[j]]])-1)*winvar(x[[grp[j]]],tr)/(h[j]*(h[j]-1))
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
nu2<-(J^2-1)/(3*Ba)
sig.A<-1-pf(Va,J-1,nu2)
nu2<-(K^2-1)/(3*Bb)
sig.B<-1-pf(Vb,K-1,nu2)
# Next, do test for interactions
Vab<-sum(D*(xbar-xtil)^2)
dfinter<-(J-1)*(K-1)
crit<-qchisq(1-alpha,dfinter)
hc<-(crit/(2*dfinter))*(1+(3*crit)/(dfinter+2))*sum(aval)
adcrit<-crit+hc
list(Qa=Va,sig.A=sig.A,Qb=Vb,sig.B=sig.B,Qab=Vab,critinter=adcrit)
}

#' Two-Way ANOVA on Trimmed Means Using Bootstrap-t Method
#'
#' Performs a two-way ANOVA based on trimmed means using a bootstrap-t method
#' for hypothesis testing. Tests main effects and interaction for independent
#' groups.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in matrix or list format. Groups are ordered as described in t2way.
#' @param tr Proportion of trimming (default is 0.2 for 20% trimming).
#' @param grp Vector to rearrange group order if needed (default is 1:p).
#' @param p Total number of groups (default is J*K).
#' @param nboot Number of bootstrap samples (default is 599).
#' @param SEED Logical indicating whether to set random seed (default is TRUE).
#'
#' @return A list with components:
#'   \item{A.p.value}{Bootstrap p-value for Factor A main effect}
#'   \item{B.p.value}{Bootstrap p-value for Factor B main effect}
#'   \item{AB.p.value}{Bootstrap p-value for interaction}
#'
#' @details
#' This function uses a bootstrap-t approach to test hypotheses in a two-way
#' ANOVA with independent groups. Data are centered by subtracting trimmed means,
#' then bootstrap samples are drawn from the centered data.
#'
#' The random seed is set to 2 if SEED=TRUE for reproducibility. The function
#' prints progress messages during bootstrap sampling.
#'
#' @export
#' @examples
#' # Create data for 2x3 design
#' set.seed(123)
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(20, mean=i)
#'
#' # Run bootstrap two-way ANOVA
#' result <- t2waybt(J=2, K=3, x=x, nboot=100)
t2waybt<-function(J,K,x,tr=.2,grp=c(1:p),p=J*K,nboot=599,SEED=TRUE){
#
#   Two-way ANOVA based on trimmed means and a bootstrap-t method
#
#   The data are assumed to be stored as described in the function t2way
#
#   The default number of bootstrap samples is nboot=599
#
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
# compute test statistics:
tests=t2way.no.p(J=J,K=K,x,tr=tr,grp=grp)
TA=NULL
TB=NULL
TAB=NULL
data=list()
xcen=list()
for(j in 1:length(x))xcen[[j]]<-x[[j]]-mean(x[[j]],tr)
print("Taking bootstrap samples. Please wait.")
for(b in 1:nboot){
for(j in 1:length(x))data[[j]]<-sample(xcen[[j]],size=length(x[[j]]),replace=TRUE)
bt=t2way.no.p(J,K,data,tr=tr,grp=grp)
TA[b]=bt$Qa
TB[b]=bt$Qb
TAB[b]=bt$Qab
}
pA<-sum(tests$Qa<=TA)/nboot
pB<-sum(tests$Qb<=TB)/nboot
pAB<-sum(tests$Qab<=TAB)/nboot
list(A.p.value=pA,B.p.value=pB,AB.p.value=pAB)
}

#' Three-Way ANOVA on Trimmed Means for Independent Groups
#'
#' Performs a J by K by L three-way ANOVA on trimmed means where all JKL groups
#' are independent. Tests all main effects, two-way interactions, and the
#' three-way interaction.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param L Number of levels for Factor C.
#' @param x Data in list or matrix format. If list, x[[1]] contains data for
#'   level (1,1,1), x[[2]] for (1,1,2), ..., x[[L]] for (1,1,L), x[[L+1]] for
#'   (1,2,1), etc.
#' @param tr Proportion of trimming (default is 0.2 for 20% trimming).
#' @param grp Vector to rearrange group order if needed (default is 1:p).
#' @param alpha Significance level (default is 0.05).
#' @param p Total number of groups (default is J*K*L).
#' @param MAT Logical; if TRUE, x is a matrix with columns for factor levels
#'   (default is FALSE).
#' @param lev.col Column numbers indicating factor levels when MAT=TRUE
#'   (default is c(1:3)).
#' @param var.col Column number for the dependent variable when MAT=TRUE
#'   (default is 4).
#' @param pr Logical indicating whether to print messages (default is TRUE).
#' @param IV1 Optional vector specifying levels of Factor A (for alternative input).
#' @param IV2 Optional vector specifying levels of Factor B (for alternative input).
#' @param IV3 Optional vector specifying levels of Factor C (for alternative input).
#'
#' @return A list with components for each effect:
#'   \item{Qa, Qb, Qc}{Test statistics for main effects A, B, C}
#'   \item{Qa.crit, Qb.crit, Qc.crit}{Critical values for main effects}
#'   \item{A.p.value, B.p.value, C.p.value}{p-values for main effects}
#'   \item{Qab, Qac, Qbc}{Test statistics for two-way interactions}
#'   \item{Qab.crit, Qac.crit, Qbc.crit}{Critical values for two-way interactions}
#'   \item{AB.p.value, AC.p.value, BC.p.value}{p-values for two-way interactions}
#'   \item{Qabc}{Test statistic for three-way interaction}
#'   \item{Qabc.crit}{Critical value for three-way interaction}
#'   \item{ABC.p.value}{p-value for three-way interaction}
#'
#' @details
#' This function performs a heteroscedastic three-way ANOVA using trimmed means.
#' It computes test statistics and p-values for all main effects, all two-way
#' interactions, and the three-way interaction.
#'
#' The function can accept data in multiple formats: list mode, matrix with
#' columns for groups, or matrix with factor level columns. When using IV1, IV2,
#' and IV3, x should be a vector of dependent variable values.
#'
#' @export
#' @examples
#' # Create data for 2x2x2 design
#' set.seed(123)
#' x <- list()
#' for(i in 1:8) x[[i]] <- rnorm(20)
#'
#' # Run three-way ANOVA
#' result <- t3way(J=2, K=2, L=2, x=x)
t3way<-function(J,K,L,x,tr=.2,grp=c(1:p),alpha=.05,p=J*K*L,MAT=FALSE,
lev.col=c(1:3),var.col=4,pr=TRUE,IV1=NULL,IV2=NULL,IV3=NULL){
#  Perform a J by K by L (three-way) anova on trimmed means where
#  all JKL groups are independent.
#
#  The R variable data is assumed to contain the raw
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
#  MAT=T, assumes data are stored in matrix with 3 columns indicating
#  levels of the three factors.
#  That is, this function calls selby2 for you.
#
if(is.data.frame(x))x=as.matrix(x)
if(!is.null(IV1[1])){
if(is.null(IV2[1]))stop("IV2 is NULL")
if(is.null(IV3[1]))stop("IV3 is NULL")
if(pr)print("Assuming x is a vector containing all of the data; the dependent variable")
xi=elimna(cbind(x,IV1,IV2,IV3))
x=fac2list(xi[,1],xi[,2:4])
J=length(unique(IV1))
K=length(unique(IV2))
L=length(unique(IV3))
p=J*K*L
}
data=x
if(MAT){
if(!is.matrix(data))stop("With MAT=T, data must be a matrix")
if(length(lev.col)!=3)stop("Argument lev.col should have 3 values")
temp=selby2(data,lev.col,var.col)
lev1=length(unique(temp$grpn[,1]))
lev2=length(unique(temp$grpn[,2]))
lev3=length(unique(temp$grpn[,3]))
gv=apply(temp$grpn,2,rank)
gvad=100*gv[,1]+10*gv[,2]+gv[,3]
grp=rank(gvad)
if(pr){
print(paste("Factor 1 has", lev1, "levels"))
print(paste("Factor 2 has", lev2, "levels"))
print(paste("Factor 3 has", lev3, "levels"))
}
if(J!=lev1)warning("J is being reset to the number of levels found")
if(K!=lev2)warning("K is being reset to the number of levels found")
if(L!=lev3)warning("K is being reset to the number of levels found")
J=lev1
K=lev2
L=lev3
data=temp$x
}
if(is.matrix(data))data=listm(data)
if(!is.list(data))stop("Data are not stored in list mode")
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
v[i]<-(length(data[[grp[i]]])-1)*winvar(data[[grp[i]]],tr)/(h[i]*(h[i]-1))
#    v contains the squared standard errors
}
v<-diag(v,p,p)   # Put squared standard errors in a diag matrix.
ij<-matrix(c(rep(1,J)),1,J)
ik<-matrix(c(rep(1,K)),1,K)
il<-matrix(c(rep(1,L)),1,L)
jm1<-J-1
cj<-diag(1,jm1,J)
for (i in 1:jm1)cj[i,i+1]<-0-1
km1<-K-1
ck<-diag(1,km1,K)
for (i in 1:km1)ck[i,i+1]<-0-1
lm1<-L-1
cl<-diag(1,lm1,L)
for (i in 1:lm1)cl[i,i+1]<-0-1
alval<-c(1:999)/1000
#  Do test for factor A
cmat<-kron(cj,kron(ik,il))  # Contrast matrix for factor A
Qa<-johan(cmat,tmeans,v,h,alpha)
A.p.value=t3pval(cmat,tmeans,v,h)
# Do test for factor B
cmat<-kron(ij,kron(ck,il))  # Contrast matrix for factor B
Qb<-johan(cmat,tmeans,v,h,alpha)
B.p.value=t3pval(cmat,tmeans,v,h)
# Do test for factor C
cmat<-kron(ij,kron(ik,cl))  # Contrast matrix for factor C
#Qc<-johan(cmat,tmeans,v,h,alpha)
for(i in 1:999){
irem<-i
Qc<-johan(cmat,tmeans,v,h,alval[i])
if(Qc$teststat>Qc$crit)break
}
C.p.value=irem/1000
# Do test for factor A by B interaction
cmat<-kron(cj,kron(ck,il))  # Contrast matrix for factor A by B
for(i in 1:999){
irem<-i
Qab<-johan(cmat,tmeans,v,h,alval[i])
if(Qab$teststat>Qab$crit)break
}
AB.p.value=irem/1000
# Do test for factor A by C interaction
cmat<-kron(cj,kron(ik,cl))  # Contrast matrix for factor A by C
for(i in 1:999){
irem<-i
Qac<-johan(cmat,tmeans,v,h,alval[i])
if(Qac$teststat>Qac$crit)break
}
AC.p.value=irem/1000
#Qac<-johan(cmat,tmeans,v,h,alpha)
# Do test for factor B by C interaction
cmat<-kron(ij,kron(ck,cl))  # Contrast matrix for factor B by C
#Qbc<-johan(cmat,tmeans,v,h,alpha)
for(i in 1:999){
irem<-i
Qbc<-johan(cmat,tmeans,v,h,alval[i])
if(Qbc$teststat>Qbc$crit)break
}
BC.p.value=irem/1000
# Do test for factor A by B by C interaction
cmat<-kron(cj,kron(ck,cl))  # Contrast matrix for factor A by B by C
#Qabc<-johan(cmat,tmeans,v,h,alpha)
for(i in 1:999){
irem<-i
Qabc<-johan(cmat,tmeans,v,h,alval[i])
if(Qabc$teststat>Qabc$crit)break
}
ABC.p.value=irem/1000
list(Qa=Qa$teststat,Qa.crit=Qa$crit,A.p.value=A.p.value,
Qb=Qb$teststat,Qb.crit=Qb$crit,
B.p.value=B.p.value,
Qc=Qc$teststat,Qc.crit=Qc$crit,C.p.value=C.p.value,
Qab=Qab$teststat,Qab.crit=Qab$crit,
AB.p.value=AB.p.value,
Qac=Qac$teststat,Qac.crit=Qac$crit,AC.p.value=AC.p.value,
Qbc=Qbc$teststat,Qbc.crit=Qbc$crit,
BC.p.value=BC.p.value,
Qabc=Qabc$teststat,Qabc.crit=Qabc$crit,ABC.p.value=ABC.p.value)
}

#' Brunner-Dette-Munk Rank-Based Two-Way ANOVA
#'
#' Performs the Brunner, Dette, and Munk rank-based ANOVA for a J by K
#' two-way independent groups design. Tests main effects and interaction
#' using rank-based methods.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in matrix or list format. If matrix, groups correspond to columns.
#'   If list, each element contains data for one group.
#' @param grp Vector specifying group ordering (default is 1:p).
#' @param p Total number of groups (default is J*K).
#'
#' @return A list with components:
#'   \item{p.valueA}{p-value for Factor A main effect}
#'   \item{p.valueB}{p-value for Factor B main effect}
#'   \item{p.valueAB}{p-value for interaction}
#'   \item{Relative.Effects}{Matrix of relative treatment effects}
#'   \item{A.F}{F statistic for Factor A}
#'   \item{B.F}{F statistic for Factor B}
#'   \item{AB.F}{F statistic for interaction}
#'
#' @details
#' Implements the Brunner, Dette, and Munk (1997) rank-based two-way ANOVA.
#' This method is robust to heteroscedasticity and non-normality.
#'
#' @references
#' Brunner, E., Dette, H., & Munk, A. (1997). Box-type approximations in
#' nonparametric factorial designs. Journal of the American Statistical
#' Association, 92, 1494-1502.
#'
#' @export
#' @examples
#' # Create data for 2x3 design
#' set.seed(123)
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(20)
#' bdm2way(J=2, K=3, x=x)
bdm2way<-function(J,K,x,grp=c(1:p),p=J*K){
#
# Perform the Brunner, Dette, Munk rank-based ANOVA
# (JASA, 1997, 92, 1494--1502)
# for a J by K independent groups design.
#
# x can be a matrix with columns corresponding to groups
# or it can have list mode.
#
if(is.matrix(x))x<-listm(x)
xx<-list()
for(j in 1:p)xx[[j]]<-x[[grp[j]]]
Ja<-matrix(1,J,J)
Ia<-diag(1,J)
Pa<-Ia-Ja/J
Jb<-matrix(1,K,K)
Ib<-diag(1,K)
Pb<-Ib-Jb/K
cona<-kron(Pa,Jb/K)
conb<-kron(Ja/J,Pb)
conab<-kron(Pa,Pb)
outA<-bdms1(xx,cona)
releff=matrix(outA$q.hat,nrow=J,ncol=K,byrow=TRUE)
outB<-bdms1(xx,conb)
outAB<-bdms1(xx,conab)
#  Could report degrees of freedom, but they are meaningless in terms of understanding the data.
list(p.valueA=outA$p.value,p.valueB=outB$p.value, p.valueAB=outAB$p.value,
Relative.Effects=releff,A.F=outA$F,B.F=outB$F,AB.F=outAB$F)
}

#' Heteroscedastic One-Way ANOVA for Trimmed Means with Effect Size
#'
#' @description
#' Performs a heteroscedastic one-way ANOVA for trimmed means using a
#' generalization of Welch's method. Unlike \code{t1way}, this function also
#' computes explanatory power and related effect sizes. Only use this function
#' with equal sample sizes; use \code{t1wayv2} in general, which calls this
#' function when sample sizes are equal.
#'
#' @param x Data matrix (with columns as groups if MAT=FALSE) or list where each
#'   element is a vector for a group. Can also be a matrix with grouping column
#'   if MAT=TRUE.
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param grp Numeric vector indicating which groups to compare. NA (default)
#'   means all groups.
#' @param MAT Logical. If FALSE (default), x is matrix with columns as groups or
#'   a list. If TRUE, x is a matrix with separate columns for group levels and values.
#' @param lev.col Column number indicating group levels when MAT=TRUE (default: 1)
#' @param var.col Column number indicating data values when MAT=TRUE (default: 2)
#'
#' @return A list with components:
#'   \item{TEST}{Test statistic (F-statistic)}
#'   \item{nu1}{Numerator degrees of freedom}
#'   \item{nu2}{Denominator degrees of freedom}
#'   \item{p.value}{p-value for the test}
#'   \item{Var.Explained}{Proportion of variance explained (R-squared analog)}
#'   \item{Effect.Size}{Effect size (square root of variance explained)}
#'
#' @details
#' The function uses Welch's heteroscedastic method for trimmed means. The
#' explanatory effect size is computed as the ratio of between-group variance
#' to total winsorized variance. If this ratio exceeds 1, the function computes
#' it as the squared correlation between group means and individual values.
#'
#' Missing values are automatically removed. This function is designed for equal
#' sample sizes; for unequal sample sizes, use \code{t1wayv2}.
#'
#' @seealso \code{\link{t1way}}, \code{\link{t1wayv2}}
#'
#' @export
t1way.effect<-function(x,tr=.2,grp=NA,MAT=FALSE,lev.col=1,var.col=2){
if(MAT){
if(!is.matrix(x))stop("With MAT=T, data must be stored in a matrix")
if(length(lev.col)!=1)stop("Argument lev.col should have 1 value")
temp=selby(x,lev.col,var.col)
x=temp$x
grp2=rank(temp$grpn)
x=x[grp2]
}
if(is.matrix(x))x<-listm(x)
grp=c(1:length(x))
if(is.na(sum(grp[1])))grp<-c(1:length(x))
if(!is.list(x))stop("Data are not stored in a matrix or in list mode.")
J<-length(grp)
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
pts=NULL
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
pts=c(pts,val)
x[[j]]<-val[xx]  # missing values have been removed
h[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
   # h is the number of observations in the jth group after trimming.
w[j]<-h[j]*(h[j]-1)/((length(x[[grp[j]]])-1)*winvar(x[[grp[j]]],tr))
xbar[j]<-mean(x[[grp[j]]],tr)
}
u<-sum(w)
xtil<-sum(w*xbar)/u
A<-sum(w*(xbar-xtil)^2)/(J-1)
B<-2*(J-2)*sum((1-w/u)^2/(h-1))/(J^2-1)
TEST<-A/(B+1)
nu1<-J-1
nu2<-1./(3*sum((1-w/u)^2/(h-1))/(J^2-1))
sig<-1-pf(TEST,nu1,nu2)
#
# Determine explanatory effect size
#
top=var(xbar)
bot=winvarN(pts,tr=tr)
if(bot>0)e.pow=top/bot
if(bot==0)e.pow=1
if(e.pow>=1){
v1=NULL
v2=NULL
for(j in 1:J){
v1=c(v1,rep(xbar[j],length(x[[j]])))
v2=c(v2,x[[j]])
}
e.pow=wincor(v1,v2,tr=tr)$cor^2
}
list(TEST=TEST,nu1=nu1,nu2=nu2,p.value=sig,Var.Explained=e.pow,
Effect.Size=sqrt(e.pow))
}

#' Power Analysis for ANOVA F Test
#'
#' Determines sample sizes or power when using the ANOVA F test. Given any three
#' of the four parameters (groups, n, delta, sig.level, power), computes the
#' fourth parameter.
#'
#' @param groups Number of groups (must be specified).
#' @param n Sample size per group (NULL if to be calculated).
#' @param delta Cohen's effect size: the sum of squared deviations among the
#'   means divided by the within-group variance (NULL if to be calculated).
#' @param sig.level Significance level (default is 0.05).
#' @param power Statistical power (NULL if to be calculated).
#'
#' @return A list with components:
#'   \item{groups}{Number of groups}
#'   \item{n}{Sample size per group}
#'   \item{delta}{Cohen's effect size}
#'   \item{sig.level}{Significance level}
#'   \item{power}{Statistical power}
#'
#' @details
#' This function is a wrapper for power.anova.test that uses Cohen's delta
#' as the effect size measure. Excluding the groups parameter, exactly one
#' of the other parameters must be NULL, and the function will solve for that
#' parameter.
#'
#' @export
#' @examples
#' # Find required sample size
#' anova_power(groups=3, delta=0.5, sig.level=0.05, power=0.80)
#'
#' # Find power given sample size
#' anova_power(groups=3, n=20, delta=0.5, sig.level=0.05)
anova_power<-function(groups=NULL,n=NULL,delta=NULL,sig.level=0.05,power=NULL){
#
# Determine sample sizes or power when using the ANOVA F test.
#
#  groups is the number of groups and must be specified.
#
#  delta is Cohen's effect size: the sum of the
#  squared devaitions among the means divided by the within group variance.
#
#  Excluding groups, all but one of the NULL arguments must be specified.
#  The function determines the value for the one argument that is NULL
#
#library(stats)
if(is.null(groups))stop("Need to specify the number of groups")
within.var=1
between.var=delta/(groups-1)
res=power.anova.test(groups=groups,n=n, between.var=between.var,
within.var=within.var,sig.level=sig.level,power=power)
list(groups=res[1]$groups,n=res[2]$n,delta=delta,
sig.level=res[5]$sig.level,power=res[6]$power)
}

#' Heteroscedastic One-Way ANOVA with Factor Variable
#'
#' @description
#' Performs a heteroscedastic one-way ANOVA for trimmed means using Welch's method.
#' Similar to \code{t1way}, but designed to work with data stored in a data frame
#' format with a separate factor variable for group identification.
#'
#' @param x Numeric vector containing the data to be analyzed
#' @param fac Factor vector indicating group membership for each observation in x
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param nboot Number of bootstrap samples (default: 100, used when EP=TRUE)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#' @param EP Logical. If TRUE, compute explanatory effect size measure (default: FALSE)
#' @param pr Logical. If TRUE (default), print message about EP option
#'
#' @return A list with components:
#'   \item{TEST}{Test statistic (F-statistic)}
#'   \item{nu1}{Numerator degrees of freedom}
#'   \item{nu2}{Denominator degrees of freedom}
#'   \item{p.value}{p-value for the test}
#'   \item{Var.Explained}{Proportion of variance explained (NA if EP=FALSE)}
#'   \item{Effect.Size}{Effect size, square root of variance explained (NA if EP=FALSE)}
#'
#' @details
#' This function is a convenient wrapper for \code{t1way} that accepts data in
#' data frame format. For example, if \code{dat} is a data frame with column 1
#' containing outcome measures and column 2 containing a factor variable for groups,
#' use \code{t1wayF(dat[,1], dat[,2])}.
#'
#' Missing values are automatically removed. When EP=TRUE, the function computes
#' the explanatory effect size by calling \code{t1wayv2}.
#'
#' @seealso \code{\link{t1way}}, \code{\link{t1wayv2}}, \code{\link{t1way.effect}}
#'
#' @export
t1wayF<-function(x,fac,tr=.2,nboot=100,SEED=TRUE,EP=FALSE,pr=TRUE){
if(!EP){
if(pr)print('To get an estimate of the explanatory measure of effect size, set EP=TRUE')
}
if(SEED)set.seed(2)
x=fac2list(x,fac)
J<-length(x)
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
pts=NULL
nval=0
for(j in 1:J)x[[j]]=elimna(x[[j]])
for(j in 1:J){
val<-x[[j]]
val<-elimna(val)
nval[j]=length(val)
pts=c(pts,val)
x[[j]]<-val # missing values have been removed
h[j]<-length(x[[j]])-2*floor(tr*length(x[[j]]))
   # h is the number of observations in the jth group after trimming.
w[j]<-h[j]*(h[j]-1)/((length(x[[j]])-1)*winvar(x[[j]],tr))
xbar[j]<-mean(x[[j]],tr)
}
u<-sum(w)
xtil<-sum(w*xbar)/u
A<-sum(w*(xbar-xtil)^2)/(J-1)
B<-2*(J-2)*sum((1-w/u)^2/(h-1))/(J^2-1)
TEST<-A/(B+1)
nu1<-J-1
nu2<-1./(3*sum((1-w/u)^2/(h-1))/(J^2-1))
sig<-1-pf(TEST,nu1,nu2)
#
# Determine explanatory effect size
#
e.pow=NA
if(EP)e.pow=t1wayv2(x,tr=tr)$Var.Explained
list(TEST=TEST,nu1=nu1,nu2=nu2,p.value=sig,Var.Explained=e.pow,
Effect.Size=sqrt(e.pow))
}

#' Bootstrap-t One-Way ANOVA for Trimmed Means
#'
#' @description
#' Tests the hypothesis of equal trimmed means across J independent groups using
#' a bootstrap-t method. Provides better small sample performance compared to
#' asymptotic methods.
#'
#' @param x Data in list mode (each element is a group) or matrix (columns are groups)
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param grp Numeric vector specifying subset of groups to compare. NA (default)
#'   means all groups.
#' @param nboot Number of bootstrap samples (default: 599)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#'
#' @return A list with components:
#'   \item{test}{Test statistic}
#'   \item{p.value}{Bootstrap p-value}
#'   \item{crit}{Critical value for the test}
#'
#' @details
#' The bootstrap-t method centers each group by its trimmed mean and resamples
#' from the centered distributions. This approach provides better control of
#' Type I error rates in small samples compared to asymptotic methods.
#'
#' Missing values are automatically removed. If data is in matrix or data frame
#' format, columns correspond to groups.
#'
#' @seealso \code{\link{t1way}}, \code{\link{t1wayv2}}, \code{\link{t1waybtv2}}
#'
#' @export
t1waybt<-function(x,tr=.2,grp=NA,nboot=599,SEED=TRUE){

if(is.matrix(x)||is.data.frame(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
if(!is.na(grp[1]))J=length(grp)
if(is.na(grp[1])){
J<-length(x)
grp<-c(1:J)
}
for(j in 1:J){
temp<-x[[j]]
x[[j]]<-temp[!is.na(temp)] # Remove any missing values.
}
bvec<-array(0,c(J,2,nboot))
hval<-vector("numeric",J)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:J){
hval[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
   # hval is the number of observations in the jth group after trimming.
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
if(ct>0){
print("Some bootstrap estimates of the test statistic could not be computed")
print("Effective number of bootstrap samples was")
print(sum(!is.na(testb)))
}
test<-t1wayv2(x,tr=tr,grp=grp)
pval<-mean(test$TEST<=testb,na.rm=TRUE)
list(test=test$TEST,p.value=pval,Var.Explained=test$Var.Explained,Effect.Size=test$Effect.Size)
}

#' Bootstrap-t One-Way ANOVA for Trimmed Means (Version 2)
#'
#' @description
#' Tests the hypothesis of equal trimmed means across J independent groups using
#' a bootstrap-t method. This version supports both traditional list/matrix input
#' and data frame input with factor variable specification.
#'
#' @param x Data in list mode (each element is a group), matrix (columns are groups),
#'   or data frame (when g and dp are specified)
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param grp Numeric vector specifying subset of groups to compare. NA (default)
#'   means all groups.
#' @param g Column number in x containing the grouping factor (NULL by default).
#'   When specified, x must be a data frame.
#' @param dp Column number in x containing the dependent variable (NULL by default).
#'   Required when g is specified.
#' @param nboot Number of bootstrap samples (default: 599)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#'
#' @return A list with components:
#'   \item{test}{Test statistic}
#'   \item{p.value}{Bootstrap p-value}
#'   \item{crit}{Critical value for the test}
#'
#' @details
#' This is an enhanced version of \code{t1waybt} that provides additional input
#' flexibility. When g is NULL, x is treated as a matrix or list. When g is
#' specified, column g of x is treated as a factor variable and column dp
#' contains the dependent variable.
#'
#' The function prints progress messages during bootstrap sampling. Missing
#' values are automatically removed.
#'
#' @seealso \code{\link{t1waybt}}, \code{\link{t1way}}, \code{\link{t1wayv2}}
#'
#' @export
t1waybtv2<-function(x,tr=.2,grp=NA,g=NULL,dp=NULL,nboot=599,SEED=TRUE){
if(!is.null(g)){
if(is.null(dp))stop("Specify a value for dp, the column containing the data")
x=fac2list(x[,dp],x[,g])
}
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
e.pow=t1wayv2(x)$Var.Explained
list(test=test$TEST,p.value=pval,Explanatory.Power=e.pow,
Effect.Size=sqrt(e.pow))
}

#' Two-Way Independent Groups ANOVA for Trimmed Means (Version 2)
#'
#' @description
#' Performs a J by K two-way ANOVA on trimmed means for completely independent
#' groups designs. This version includes support for data frame input with
#' factor variables.
#'
#' @param J Number of levels for first factor
#' @param K Number of levels for second factor
#' @param data Data in list mode (length J*K) or matrix (columns are groups).
#'   Can also be a data frame when g and dp are specified.
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param grp Numeric vector specifying which groups to include (default: all groups)
#' @param p Total number of groups, J*K (default: J*K)
#' @param g Vector of length 2 specifying columns in data containing the two
#'   factor variables (NULL by default). When specified, data must be a data frame.
#' @param dp Column number in data containing the dependent variable (NULL by default).
#'   Required when g is specified.
#' @param pr Logical. If TRUE (default), print warning messages
#'
#' @return A list with components:
#'   \item{Qa}{Test statistic for Factor A main effect}
#'   \item{Qb}{Test statistic for Factor B main effect}
#'   \item{Qab}{Test statistic for A*B interaction}
#'   \item{A.p.value}{p-value for Factor A}
#'   \item{B.p.value}{p-value for Factor B}
#'   \item{AB.p.value}{p-value for interaction}
#'   \item{A.df}{Degrees of freedom for Factor A}
#'   \item{B.df}{Degrees of freedom for Factor B}
#'   \item{AB.df}{Degrees of freedom for interaction}
#'
#' @details
#' Data organization: In list mode, data[[1]] contains observations for level (1,1),
#' data[[2]] for level (1,2), etc. The function tests main effects for both factors
#' and their interaction using Welch-type statistics for trimmed means.
#'
#' Missing values are automatically removed. When g is specified, it should be a
#' vector of length 2 indicating which columns contain the factor variables.
#'
#' @seealso \code{\link{t2way}}, \code{\link{t3way}}
#'
#' @export
t2wayv2<-function(J,K,data,tr=.2,grp=c(1:p),p=J*K,g=NULL,dp=NULL,pr=TRUE){
if(!is.null(g[1])){
if(length(g)!=2)stop("Argument g should have two values")
if(is.null(dp[1]))
stop("Specify a value for dp, the column containing the data")
data=fac2list(data[,dp],data[,g])
}
if(is.matrix(data))data=listm(data)
if(!is.list(data))stop("Data are not stored in list mode")
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
data[[grp[i]]]=elimna(data[[grp[i]]])
tmeans[i]<-mean(data[[grp[i]]],tr)
h[i]<-length(data[[grp[i]]])-2*floor(tr*length(data[[grp[i]]]))
#    h is the effective sample size
   if(winvar(data[[grp[i]]],tr)==0)print(paste('The Winsorized variance is zero for group',i))
v[i]<-(length(data[[grp[i]]])-1)*winvar(data[[grp[i]]],tr)/(h[i]*(h[i]-1))
#    v contains the squared standard errors
}
v<-diag(v,p,p)   # Put squared standard errors in a diag matrix.
ij<-matrix(c(rep(1,J)),1,J)
ik<-matrix(c(rep(1,K)),1,K)
jm1<-J-1
cj<-diag(1,jm1,J)
for (i in 1:jm1)cj[i,i+1]<-0-1
km1<-K-1
ck<-diag(1,km1,K)
for (i in 1:km1)ck[i,i+1]<-0-1
#  Do test for factor A
#cmat<-kron(cj,kron(ik,il))  # Contrast matrix for factor A
cmat<-kron(cj,ik)  # Contrast matrix for factor A
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
Qa<-johan(cmat,tmeans,v,h,alval[i])
if(Qa$teststat>Qa$crit)break
}
A.p.value=irem/1000
# Do test for factor B
cmat<-kron(ij,ck)  # Contrast matrix for factor B
for(i in 1:999){
irem<-i
Qb<-johan(cmat,tmeans,v,h,alval[i])
if(Qb$teststat>Qb$crit)break
}
B.p.value=irem/1000
# Do test for factor A by B interaction
cmat<-kron(cj,ck)  # Contrast matrix for factor A by B
for(i in 1:999){
irem<-i
Qab<-johan(cmat,tmeans,v,h,alval[i])
if(Qab$teststat>Qab$crit)break
}
AB.p.value=irem/1000
tmeans=matrix(tmeans,J,K,byrow=TRUE)
list(Qa=Qa$teststat,A.p.value=A.p.value,
Qb=Qb$teststat,B.p.value=B.p.value,
Qab=Qab$teststat,AB.p.value=AB.p.value,means=tmeans)
}

#' Percentile Bootstrap for Interactions in Split-Plot Design
#'
#' @description
#' Tests for interactions in a J by K split-plot (mixed) design using percentile
#' bootstrap. The analysis examines whether differences among dependent groups
#' (Factor B) vary across levels of the independent factor (Factor A).
#'
#' @param J Number of levels for Factor A (between-subjects factor)
#' @param K Number of levels for Factor B (within-subjects factor)
#' @param x Data in list mode (length J*K) or matrix (columns are groups).
#'   For list mode: x[[1]] through x[[K]] are data for first level of Factor A
#'   at all K levels of Factor B, x[[K+1]] through x[[2K]] are for second level
#'   of Factor A, etc.
#' @param est Estimator function to use (default: tmean for trimmed mean)
#' @param JK Total number of groups, J*K (default: J*K)
#' @param grp Numeric vector specifying which groups to include (default: all groups)
#' @param nboot Number of bootstrap samples (default: 500)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#' @param pr Logical. If TRUE (default), print informational messages
#' @param ... Additional arguments passed to the estimator function
#'
#' @return A list with components:
#'   \item{test}{Test statistic for interaction}
#'   \item{p.value}{Bootstrap p-value}
#'   \item{crit.value}{Critical value}
#'
#' @details
#' This function analyzes split-plot designs by computing difference scores among
#' all pairs of dependent groups (Factor B levels) and testing whether these
#' differences vary across the independent factor (Factor A levels).
#'
#' Data organization: For each subject, measurements at all K levels of Factor B
#' are in the same row. Missing values are handled by listwise deletion within
#' each level of Factor A.
#'
#' As of October 2014, the default estimator is tmean (trimmed mean).
#'
#' @seealso \code{\link{sppba}}, \code{\link{spp}}
#'
#' @export
sppbi<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),nboot=500,SEED=TRUE,pr=TRUE,...){
if(pr)print('As of Oct. 2014, argument est defaults to tmean')
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}

JK<-J*K
MJ<-(J^2-J)/2
MK<-(K^2-K)/2
JMK<-J*MK
Jm<-J-1
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
jp<-1-K
kv<-0
kv2<-0
for(j in 1:J){
jp<-jp+K
xmat<-matrix(NA,ncol=K,nrow=length(x[[jp]]))
for(k in 1:K){
kv<-kv+1
xmat[,k]<-x[[kv]]
}
xmat<-elimna(xmat)
for(k in 1:K){
kv2<-kv2+1
x[[kv2]]<-xmat[,k]
}}
xx<-x
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
# Next determine the n_j values
nvec<-NA
jp<-1-K
for(j in 1:J){
jp<-jp+K
nvec[j]<-length(x[[jp]])
}
#
# Now take bootstrap samples from jth level
# of Factor A and average K  corresponding estimates
# of location.
#
bloc<-matrix(NA,ncol=J,nrow=nboot)
#print("Taking bootstrap samples. Please wait.")
mvec<-NA
it<-0
for(j in 1:J){
paste("Working on level ",j," of Factor A")
x<-matrix(NA,nrow=nvec[j],ncol=MK)
#
im<-0
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
im<-im+1
kp<-j*K+k-K
kpp<-j*K+kk-K
x[,im]<-xx[[kp]]-xx[[kpp]]
it<-it+1
mvec[it]<-est(x[,im],...)
}}}
data<-matrix(sample(nvec[j],size=nvec[j]*nboot,replace=TRUE),nrow=nboot)
bvec<-matrix(NA,ncol=MK,nrow=nboot)
mat<-listm(x)
for(k in 1:MK){
temp<-x[,k]
bvec[,k]<-apply(data,1,rmanogsub,temp,est,...) # An nboot by MK matrix
}
if(j==1)bloc<-bvec
if(j>1)bloc<-cbind(bloc,bvec)
}
#
MJMK<-MJ*MK
con<-matrix(0,nrow=JMK,ncol=MJMK)
cont<-matrix(0,nrow=J,ncol=MJ)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j<jj){
ic<-ic+1
cont[j,ic]<-1
cont[jj,ic]<-0-1
}}}
tempv<-matrix(0,nrow=MK-1,ncol=MJ)
con1<-rbind(cont[1,],tempv)
for(j in 2:J){
con2<-rbind(cont[j,],tempv)
con1<-rbind(con1,con2)
}
con<-con1
if(MK>1){
for(k in 2:MK){
con1<-push(con1)
con<-cbind(con,con1)
}}
bcon<-t(con)%*%t(bloc) #C by nboot matrix
tvec<-t(con)%*%mvec
tvec<-tvec[,1]
tempcen<-apply(bcon,1,mean)
vecz<-rep(0,ncol(con))
bcon<-t(bcon)
temp=bcon
for(ib in 1:nrow(temp))temp[ib,]=temp[ib,]-tempcen+tvec
smat<-var(temp)
if(sum(is.na(smat))==0){
chkrank<-qr(smat)$rank
bcon<-rbind(bcon,vecz)
if(chkrank==ncol(smat))dv<-mahalanobis(bcon,tvec,smat)
if(chkrank<ncol(smat)){
smat<-ginv(smat)
dv<-mahalanobis(bcon,tvec,smat,inverted=TRUE)
}}
if(sum(is.na(smat))>0)print('Computational Problem. Try est=tmean or use function spmcpi or tsplitbt')
bplus<-nboot+1
sig.level<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
list(p.value=sig.level,psihat=tvec,con=con)
}

#' Percentile Bootstrap for Main Effects in Split-Plot Design
#'
#' @description
#' Tests for main effects of the between-subjects factor (Factor A) in a J by K
#' split-plot design using percentile bootstrap. Two analysis approaches are available.
#'
#' @param J Number of levels for Factor A (between-subjects factor)
#' @param K Number of levels for Factor B (within-subjects factor)
#' @param x Data in list mode (length J*K) or matrix (columns are groups).
#'   Organization same as \code{sppbi}: x[[1]] through x[[K]] are for first level
#'   of Factor A, x[[K+1]] through x[[2K]] for second level, etc.
#' @param est Estimator function to use (default: tmean for trimmed mean)
#' @param JK Total number of groups, J*K (default: J*K)
#' @param grp Numeric vector specifying which groups to include (default: all groups)
#' @param avg Logical. If TRUE (default), average K location measures for each level
#'   of Factor A and test all pairwise differences. If FALSE, test K simultaneous
#'   equalities (one for each level of Factor B).
#' @param nboot Number of bootstrap samples (default: 500)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#' @param MC Logical. If TRUE, use parallel processing via mclapply (default: FALSE)
#' @param MDIS Logical. Related to multivariate analysis (default: FALSE)
#' @param pr Logical. If TRUE (default), print informational messages
#' @param ... Additional arguments passed to the estimator function
#'
#' @return A list with components:
#'   \item{p.value}{Bootstrap p-value for Factor A main effect}
#'   \item{psihat}{Estimated contrasts}
#'   \item{con}{Contrast matrix used}
#'
#' @details
#' When avg=TRUE, the function averages the K location measures for each level of
#' Factor A and tests whether all pairwise differences equal zero. When avg=FALSE,
#' it tests K simultaneous hypotheses: for each level k of Factor B, whether
#' theta_1k = ... = theta_Jk.
#'
#' Missing values are handled by listwise deletion within each level of Factor A.
#' As of October 2014, the default estimator is tmean (trimmed mean).
#'
#' @seealso \code{\link{sppbi}}, \code{\link{spp}}
#'
#' @export
sppba<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),avg=TRUE,nboot=500,SEED=TRUE,
MC=FALSE,MDIS=FALSE,pr=TRUE,...){
if(pr)print('As of Oct. 2014 the argument est defaults to tmean')
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}

JK<-J*K
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
jp<-1-K
kv<-0
kv2<-0
for(j in 1:J){
jp<-jp+K
xmat<-matrix(NA,ncol=K,nrow=length(x[[jp]]))
for(k in 1:K){
kv<-kv+1
xmat[,k]<-x[[kv]]
}
xmat<-elimna(xmat)
for(k in 1:K){
kv2<-kv2+1
x[[kv2]]<-xmat[,k]
}
}
xx<-x
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
# Next determine the n_j values
nvec<-NA
jp<-1-K
for(j in 1:J){
jp<-jp+K
nvec[j]<-length(x[[jp]])
}
#
# Now take bootstrap samples from jth level
# of Factor A.
#
bloc<-matrix(NA,nrow=J,ncol=nboot)
#print("Taking bootstrap samples. Please wait.")
mvec<-NA
ik<-0
for(j in 1:J){
paste("Working on level ",j," of Factor A")
x<-matrix(NA,nrow=nvec[j],ncol=K)
#
for(k in 1:K){
ik<-ik+1
x[,k]<-xx[[ik]]
if(!avg)mvec[ik]<-est(xx[[ik]],...)
}
tempv<-apply(x,2,est,...)
data<-matrix(sample(nvec[j],size=nvec[j]*nboot,replace=TRUE),nrow=nboot)
bvec<-matrix(NA,ncol=K,nrow=nboot)
for(k in 1:K){
temp<-x[,k]
bvec[,k]<-apply(data,1,rmanogsub,temp,est,...) # An nboot by K matrix
}
if(avg){
mvec[j]<-mean(tempv)
bloc[j,]<-apply(bvec,1,mean)
}
if(!avg){
if(j==1)bloc<-bvec
if(j>1)bloc<-cbind(bloc,bvec)
}
}
if(avg){
d<-(J^2-J)/2
con<-matrix(0,J,d)
id<-0
Jm<-J-1
for (j in 1:Jm){
jp<-j+1
for(k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
if(!avg){
MJK<-K*(J^2-J)/2 # NUMBER OF COMPARISONS
JK<-J*K
MJ<-(J^2-J)/2
cont<-matrix(0,nrow=J,ncol=MJ)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j<jj){
ic<-ic+1
cont[j,ic]<-1
cont[jj,ic]<-0-1
}}}
tempv<-matrix(0,nrow=K-1,ncol=MJ)
con1<-rbind(cont[1,],tempv)
for(j in 2:J){
con2<-rbind(cont[j,],tempv)
con1<-rbind(con1,con2)
}
con<-con1
if(K>1){
for(k in 2:K){
con1<-push(con1)
con<-cbind(con,con1)
}}}
if(!avg)bcon<-t(con)%*%t(bloc) #C by nboot matrix
if(avg)bcon<-t(con)%*%(bloc)
tvec<-t(con)%*%mvec
tvec<-tvec[,1]
tempcen<-apply(bcon,1,mean)
vecz<-rep(0,ncol(con))
bcon<-t(bcon)
temp=bcon
for(ib in 1:nrow(temp))temp[ib,]=temp[ib,]-tempcen+tvec
bcon<-rbind(bcon,vecz)
if(!MDIS){
if(!MC)dv=pdis(bcon,center=tvec,na.rm=FALSE)
if(MC)dv=pdisMC(bcon,center=tvec)
lbcon=length(elimna(bcon))
bplus<-nboot+1
if(lbcon<bplus){
print(paste('Effective value for nboot is', lbcon-1))
nboot=lbcon-1
}
}
if(MDIS){
smat<-var(temp)
bcon<-rbind(bcon,vecz)
chkrank<-qr(smat)$rank
if(chkrank==ncol(smat))dv<-mahalanobis(bcon,tvec,smat)
if(chkrank<ncol(smat)){
smat<-ginv(smat)
dv<-mahalanobis(bcon,tvec,smat,inverted=TRUE)
}}
bplus<-nboot+1
sig.level<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
list(p.value=sig.level,psihat=tvec,con=con)
}

#' Multivariate ANOVA for Trimmed Means (Two Groups)
#'
#' @description
#' Performs multivariate ANOVA (MANOVA) for two independent groups using trimmed
#' means. Based on a generalization of the Yanagihara and Yuan (2005) approach
#' to the multivariate Behrens-Fisher problem.
#'
#' @param x1 Matrix where rows are observations and columns are variables for group 1
#' @param x2 Matrix where rows are observations and columns are variables for group 2
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#'
#' @return A list with components:
#'   \item{test.stat}{Test statistic (F-statistic analog)}
#'   \item{p.value}{p-value for the test}
#'
#' @details
#' This function implements equation (2.7) from Yanagihara and Yuan (2005),
#' "Three approximate solutions to the multivariate Behrens-Fisher problem,"
#' Communications in Statistics--Simulation and Computation, 34, 975-988.
#'
#' The method extends the multivariate Behrens-Fisher test to trimmed means,
#' allowing for heterogeneous covariance matrices between groups. Missing values
#' are removed via listwise deletion.
#'
#' @references
#' Yanagihara, H., & Yuan, K. H. (2005). Three approximate solutions to the
#' multivariate Behrens-Fisher problem. Communications in Statistics--Simulation
#' and Computation, 34, 975-988.
#'
#' @seealso \code{\link{MULtr.anova}}
#'
#' @export
YYmanova<-function(x1,x2,tr=.2){
x1=elimna(x1)
x2=elimna(x2)
s1=winall(x1,tr=tr)$cov
s2=winall(x2,tr=tr)$cov
n1=nrow(x1)
n2=nrow(x2)
n=n1+n2
g1=floor(n1*tr)
g2=floor(n2*tr)
h1=n1-2*g1
h2=n2-2*g2
h=h1+h2
sbar=n2*s1/n+n1*s2/n
sbarinv=solve(sbar)
psi1=n2^2*(n-2)*(sum(diag(s1%*%sbarinv)))^2/(n^2*(n1-1))+
n1^2*(n-2)*(sum(diag(s2%*%sbarinv)))^2/(n^2*(n2-1))
psi2=n2^2*(n-2)*(sum(diag(s1%*%sbarinv%*%s1%*%sbarinv)))/(n^2*(n1-1))+
n1^2*(n-2)*(sum(diag(s2%*%sbarinv%*%s2%*%sbarinv)))/(n^2*(n2-1))
p=ncol(x1)
theta1=(p*psi1+(p-2)*psi2)/(p*(p+2))
theta2=(psi1+2*psi2)/(p*(p+2))
nuhat=(h-2-theta1)^2/((h-2)*theta2-theta1)
xb1=apply(x1,2,mean,tr=tr)
xb2=apply(x2,2,mean,tr=tr)
dif=xb1-xb2
dif=as.matrix(dif)
Ttest=t(dif)%*%solve((n1-1)*s1/(h1*(h1-1))+(n2-1)*s2/(h2*(h2-1)))%*%dif
TF=(n-2-theta1)*Ttest/((n-2)*p)
pv=1-pf(TF,p,nuhat)
list(test.stat=TF,p.value=pv)
}

#' Multivariate ANOVA for Trimmed Means Using Johansen's Method
#'
#' @description
#' Performs multivariate ANOVA (MANOVA) for J independent groups using trimmed
#' means based on Johansen's method. Tests whether J groups have identical
#' multivariate trimmed mean vectors.
#'
#' @param x Data in list mode where x[[j]] is an n_j by p matrix for group j,
#'   or a matrix when J and p are specified (stored as expected by bwtrim)
#' @param J Number of groups (NULL by default, required if x is a matrix)
#' @param p Number of variables (NULL by default, required if x is a matrix)
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param alpha Significance level for critical value computation (default: 0.05)
#'
#' @return A list with components:
#'   \item{test.stat}{Test statistic}
#'   \item{crit.value}{Critical value at specified alpha level}
#'
#' @details
#' When x is a list, length(x) = J (number of groups) and x[[j]] is an n_j by p
#' matrix where p is the number of variables. When x is a matrix, both J and p
#' must be specified.
#'
#' The method uses winsorized covariance matrices and weighted grand means. The
#' test statistic follows an approximate chi-squared distribution with adjustment
#' for trimmed data.
#'
#' To obtain a p-value, use the function \code{MULAOVp}. Missing values are
#' removed via listwise deletion within each group.
#'
#' @seealso \code{\link{YYmanova}}, \code{\link{MULAOVp}}
#'
#' @export
MULtr.anova<-function(x,J=NULL,p=NULL,tr=.2,alpha=.05){
if(is.matrix(x) || is.data.frame(x)){
if(is.null(J) && is.null(p))stop("Specify J or P")
x=MAT2list(x,p=p,J=J)
}
x=lapply(x,as.matrix)
x=lapply(x,elimna)
p=ncol(x[[1]])
iden=diag(p)
J=length(x)
tvec=list()
nval=lapply(x,nrow)
Rtil=lapply(x,wincov,tr=tr)
tvec=lapply(x,mmean,tr=tr)
g=list()
gmean=rep(0,p) # grand mean eventually
groupm=list()
Wsum=matrix(0,ncol=p,nrow=p)
W=list()
f=0
Aw=0
for(j in 1:J){
dimnames(x[[j]])=list(NULL,NULL)
tvec[[j]]=as.matrix(tvec[[j]])
g[[j]]=floor(nval[[j]]*tr)
Rtil[[j]]=Rtil[[j]]*(nval[[j]]-1)/((nval[[j]]-2*g[[j]])*(nval[[j]]-2*g[[j]]-1))
f[j]=nval[[j]]-2*g[[j]]-1
W[[j]]=solve(Rtil[[j]])
groupm[[j]]=apply(x[[j]],2,tmean,tr=tr)
Wsum=Wsum+W[[j]]
gmean=gmean+W[[j]]%*%tvec[[j]]
}
Wsuminv=solve(Wsum)
for(j in 1:J){
temp=iden-Wsuminv%*%W[[j]]
tempsq=temp%*%temp
Aw=Aw+(sum(diag(tempsq))+(sum(diag(temp)))^2)/f[j]
}
Aw=Aw/2
gmean=as.matrix(gmean)
gmean=solve(Wsum)%*%gmean # Final weighted grand mean
df=p*(J-1)
crit<-qchisq(1-alpha,df)
crit<-crit+(crit/(2*df))*(Aw+3*Aw*crit/(df+2))
test=0
for(k in 1:p){
for(m in 1:p){
for(j in 1:J){
test=test+W[[j]][k,m]*(groupm[[j]][m]-gmean[m])*(groupm[[j]][k]-gmean[k])
}}}
list(test.stat=test,crit.value=crit)
}

#' Heteroscedastic One-Way ANOVA with Effect Size (Version 2)
#'
#' @description
#' Performs a heteroscedastic one-way ANOVA for trimmed means using Welch's method
#' and computes explanatory power and effect size. This is the recommended general
#' version that handles both equal and unequal sample sizes. For unequal n, it
#' calls \code{t1way.effect}.
#'
#' @param x Data in multiple formats: matrix (columns as groups if MAT=FALSE),
#'   list (each element is a group), data frame, or vector (when IV is specified)
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param grp Numeric vector specifying subset of groups to compare. NA (default)
#'   means all groups.
#' @param MAT Logical. If FALSE (default), matrix columns are groups. If TRUE,
#'   use lev.col and var.col to specify group and data columns.
#' @param lev.col Column number indicating group levels when MAT=TRUE (default: 1)
#' @param var.col Column number indicating data values when MAT=TRUE (default: 2)
#' @param nboot Number of bootstrap samples for unequal n case (default: 100)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#' @param pr Logical. If TRUE (default), print informational messages
#' @param IV Vector of group identifiers when x is a vector of all data (NULL by default)
#' @param loc.fun Location function for computing group centers (default: median)
#'
#' @return A list with components:
#'   \item{TEST}{Test statistic (F-statistic)}
#'   \item{nu1}{Numerator degrees of freedom}
#'   \item{nu2}{Denominator degrees of freedom}
#'   \item{p.value}{p-value for the test}
#'   \item{Var.Explained}{Proportion of variance explained (R-squared analog)}
#'   \item{Effect.Size}{Effect size (square root of variance explained)}
#'
#' @details
#' This is the recommended general-purpose function for one-way ANOVA with effect
#' sizes. For equal sample sizes, the effect size is computed as the ratio of
#' between-group to total winsorized variance. For unequal sample sizes, bootstrap
#' resampling equalizes sample sizes before computing effect size.
#'
#' When IV is specified, x is treated as a vector containing all data and IV
#' contains the corresponding group identifiers. Missing values are automatically
#' removed.
#'
#' @seealso \code{\link{t1way}}, \code{\link{t1way.effect}}, \code{\link{t1wayF}}
#'
#' @export
t1wayv2<-function(x,tr=.2,grp=NA,MAT=FALSE,lev.col=1,var.col=2,nboot=100,SEED=TRUE,pr=TRUE,IV=NULL,loc.fun=median){
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(MAT){
if(!is.matrix(x))stop("With MAT=T, data must be stored in a matrix")
if(length(lev.col)!=1)stop("Argument lev.col should have 1 value")
temp=selby(x,lev.col,var.col)
x=temp$x
grp2=rank(temp$grpn)
x=x[grp2]
}
if(!is.null(IV[1])){
if(pr)print("Assuming x is a vector containing all of the data, the dependent variable")
xi=elimna(cbind(x,IV))
x=fac2list(xi[,1],xi[,2])
}
if(is.matrix(x))x<-listm(x)
if(is.na(sum(grp[1])))grp<-c(1:length(x))
if(!is.list(x))stop("Data are not stored in a matrix or in list mode.")
J<-length(grp)
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
pts=NULL
nval=0
for(j in 1:J)x[[j]]=elimna(x[[j]])
for(j in 1:J){
val<-x[[j]]
val<-elimna(val)
nval[j]=length(val)
pts=c(pts,val)
x[[j]]<-val # missing values have been removed
h[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
   # h is the number of observations in the jth group after trimming.
   if(winvar(x[[grp[j]]],tr)==0)print(paste('The Winsorized variance is zero for group',j))
w[j]<-h[j]*(h[j]-1)/((length(x[[grp[j]]])-1)*winvar(x[[grp[j]]],tr))
xbar[j]<-mean(x[[grp[j]]],tr)
}
u<-sum(w)
xtil<-sum(w*xbar)/u
A<-sum(w*(xbar-xtil)^2)/(J-1)
B<-2*(J-2)*sum((1-w/u)^2/(h-1))/(J^2-1)
TEST<-A/(B+1)
nu1<-J-1
nu2<-1./(3*sum((1-w/u)^2/(h-1))/(J^2-1))
sig<-1-pf(TEST,nu1,nu2)
nv=lapply(x,length)
#
# Determine explanatory effect size
#
chkn=var(nval)
if(chkn==0){
top=var(xbar)
bot=winvarN(pts,tr=tr)
e.pow=top/bot
}
if(chkn!=0){
vals=0
N=min(nval)
xdat=list()
for(i in 1:nboot){
for(j in 1:J){
xdat[[j]]=sample(x[[j]],N)
}
vals[i]=t1way.effect(xdat,tr=tr)$Var.Explained
}
e.pow=loc.fun(vals,na.rm=TRUE)
}
list(TEST=TEST,nu1=nu1,nu2=nu2,n=nv,p.value=sig,Explanatory.Power=e.pow,
Effect.Size=sqrt(e.pow))
}

#' One-Way ANOVA for Medians Using Projection-Based Depth
#'
#' Tests the global hypothesis that J independent groups have equal medians.
#' Performs well even when there are tied values, which is common with medians.
#'
#' @inheritParams common-params
#'
#' @return A list with components from \code{\link{pbadepth}}:
#'   \item{test}{The test statistic.}
#'   \item{p.value}{Bootstrap p-value.}
#'   \item{est.1}{Estimates for first group vs all others.}
#'   \item{est.2}{Estimates for second group vs all others.}
#'   \item{...}{Additional estimates for each group.}
#'
#' @details
#' This function is a wrapper for \code{pbadepth} that uses the Harrell-Davis
#' estimator for medians. It is specifically designed to handle tied values,
#' which frequently occur when comparing medians.
#'
#' The test uses projection-based depth and bootstrap resampling to generate
#' the null distribution. The \code{op} parameter controls the number of
#' outliers removed from each tail before analysis.
#'
#' @seealso \code{\link{Qanova}}, \code{\link{pbadepth}}, \code{\link{t1way}},
#'   \code{\link{pbanova}}
#'
#' @export
#' @examples
#' # Three groups with tied values
#' set.seed(123)
#' x1 <- round(rnorm(30, mean=0))
#' x2 <- round(rnorm(30, mean=0.5))
#' x3 <- round(rnorm(30, mean=1))
#' MEDanova(list(x1, x2, x3), nboot=500)
MEDanova<-function(x,op=3,nboot=600,MC=FALSE,SEED=TRUE){
#
# Test global hypothesis that J independent groups
# have equal medians.
# Performs well when there are tied values.
#
# Basically, use pbadepth in conjunction with the Harrell--Davis
# estimator.
#
output=pbadepth(x,est=hd,allp=TRUE,SEED=SEED,op=op,nboot=nboot,MC=MC)
output
}

#' One-Way ANOVA for Quantiles Using Projection-Based Depth
#'
#' Tests the global hypothesis that J independent groups have equal quantiles
#' (default: medians). A generalization of \code{MEDanova} that works with any
#' quantile. Performs well even when there are tied values.
#'
#' @inheritParams common-params
#'
#' @return A list with components from \code{\link{pbadepth}}:
#'   \item{test}{The test statistic.}
#'   \item{p.value}{Bootstrap p-value.}
#'   \item{est.1}{Estimates for first group vs all others.}
#'   \item{est.2}{Estimates for second group vs all others.}
#'   \item{...}{Additional estimates for each group.}
#'
#' @details
#' This function is a wrapper for \code{pbadepth} that uses the Harrell-Davis
#' estimator for the qth quantile. It is specifically designed to handle tied
#' values, which can occur when comparing quantiles.
#'
#' The test uses projection-based depth and bootstrap resampling. A warning is
#' issued if any group has fewer than 20 observations, as Type I error control
#' may be compromised with very small samples.
#'
#' @seealso \code{\link{MEDanova}}, \code{\link{pbadepth}}, \code{\link{t1way}}
#'
#' @export
#' @examples
#' # Compare medians (q=0.5) for three groups
#' set.seed(123)
#' x1 <- round(rnorm(30, mean=0))
#' x2 <- round(rnorm(30, mean=0.5))
#' x3 <- round(rnorm(30, mean=1))
#' Qanova(list(x1, x2, x3), q=0.5, nboot=500)
#'
#' # Compare 75th percentiles
#' Qanova(list(x1, x2, x3), q=0.75, nboot=500)
Qanova<-function(x,q=.5,op=3,nboot=2000,MC=FALSE,SEED=TRUE){
#
# Test global hypothesis that J independent groups
# have equal medians.
# Performs well when there are tied values.
#
# Basically, use pbadepth in conjunction with the Harrell--Davis
# estimator.
#
#
if(is.matrix(x)	|| is.data.frame(x))x=listm(x)
chkcar=NA
for(j in 1:length(x))chkcar[j]=length(unique(x[[j]]))
if(min(chkcar<20)){
print('Warning: Sample size is less than')
print('20 for one or more groups. Type I error might not be controlled')
}
output=pbadepth(x,est=hd,q=q,allp=TRUE,SEED=SEED,op=op,nboot=nboot,MC=MC,na.rm=TRUE)
output
}

#' Bootstrap One-Way ANOVA for Trimmed Means Using Squared Ranks
#'
#' @description
#' Performs one-way ANOVA for trimmed means with independent groups using a
#' bootstrap method based on squared ranks. Method studied by Ozdemir et al.
#'
#' @param x Data in list mode (each element is a group) or matrix (columns are groups)
#' @param alpha Significance level (default: 0.05)
#' @param nboot Number of bootstrap samples (default: 599)
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#'
#' @return A list with components:
#'   \item{test.stat}{Test statistic}
#'   \item{crit.value}{Bootstrap critical value at specified alpha}
#'   \item{p.value}{Bootstrap p-value}
#'
#' @details
#' This function uses the squared rank method with bootstrap resampling. Data
#' are centered by their group trimmed means and bootstrap samples are drawn
#' from the centered distributions. The test statistic is computed using
#' \code{btsqrk}.
#'
#' Missing values are automatically removed. The method provides robust inference
#' for comparing trimmed means across groups.
#'
#' @references
#' Ozdemir et al. (method reference for squared rank bootstrap approach)
#'
#' @seealso \code{\link{t1waybt}}, \code{\link{t1way}}, \code{\link{btsqrk}}
#'
#' @export
t1waybtsqrk<-function(x,alpha=.05,nboot=599,tr=0.2,SEED=TRUE){
if(SEED)set.seed(2)
B=nboot
if(is.matrix(x))x=listm(x)
x=lapply(x,elimna)
T.test<-btsqrk(x,alpha=alpha,tr=tr)$teststat
means<-c()
ylist<-list(0)
TT<-c()
b<-floor((1-alpha)*B)
means<-sapply(x,mean,tr)
k<-length(x)
for (i in 1:B)
{
	for (j in 1:k)
	{ylist[[j]]<-(sample(x[[j]],length(x[[j]]),replace=TRUE)-means[j])}
	TT<-c(TT,btsqrk(ylist,alpha,tr)$teststat)
}
TT=sort(TT)
pval<-mean(T<=TT,na.rm=TRUE)
list(test.stat=T.test,crit.value=TT[b],p.value=pval)
}

#' ANOVA for Discrete/Categorical Data (Multinomial Distributions)
#'
#' @description
#' Tests the global hypothesis that two or more independent groups have identical
#' discrete (multinomial) distributions. Uses a generalization of the Storer-Kim
#' method with bootstrap resampling.
#'
#' @param x Data in matrix form (rows are observations, columns are groups) or
#'   list mode (each element is a group)
#' @param nboot Number of bootstrap samples (default: 500)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#'
#' @return A list with components:
#'   \item{test}{Test statistic}
#'   \item{p.value}{Bootstrap p-value}
#'
#' @details
#' This method is designed for situations where sample sizes are relatively small
#' and data are discrete/categorical. It can detect differences that might be
#' missed when using only measures of location. Type I error control is excellent,
#' even with n=10.
#'
#' The test uses multinomial resampling under the null hypothesis to generate
#' the bootstrap distribution. For each category value, the function computes
#' the variance of proportions across groups and sums these variances.
#'
#' As an alternative, consider using chi-squared test via \code{disc2.chi.sq}.
#' This function requires the mc2d package for multinomial resampling.
#'
#' @seealso \code{\link{discANOVA.sub}}, \code{\link{disc2.chi.sq}}
#'
#' @export
discANOVA<-function(x,nboot=500,SEED=TRUE){
if(is.matrix(x) || is.data.frame(x))x=listm(x)
library(mc2d)
if(SEED)set.seed(2)
vals=lapply(x,unique)
vals=sort(elimna(list2vec(vals)))
K=length(unique(vals))
n=lapply(x,length)
n=list2vec(n)
J=length(x)
step1=discANOVA.sub(x)
test=step1$test
C1=step1$C1
HT=NULL
for(i in 1:K)HT[i]=mean(C1[i,])
tv=NULL
TB=NA
VP=NA
B1hat=NA
xx=list()
for(ib in 1:nboot){
xx=list()
for(j in 1:J){
temp=rmultinomial(n[j],1,HT)
xx[[j]]=which(temp[1,]==1)
for(i in 2:n[j])xx[[j]][i]=which(temp[i,]==1)
}
TB[ib]=discANOVA.sub(xx)$test
}
pv=1-mean(test>TB)-.5*mean(test==TB)
list(test=test,p.value=pv)
}

#' Helper Function for discANOVA
#'
#' @description
#' Internal helper function for \code{discANOVA} that computes the test statistic
#' and proportion matrix for discrete/categorical data.
#'
#' @param x Data in list mode where each element is a group
#'
#' @return A list with components:
#'   \item{test}{Test statistic (sum of variances across categories)}
#'   \item{C1}{K by J matrix of proportions, where K is number of unique values
#'     and J is number of groups}
#'
#' @details
#' This function computes proportions for each unique value in each group,
#' then calculates the variance of proportions across groups for each unique
#' value. The test statistic is the sum of these variances.
#'
#' Missing values are removed. This is an internal function primarily called
#' by \code{discANOVA}.
#'
#' @seealso \code{\link{discANOVA}}
#'
#' @export
discANOVA.sub<-function(x){
x=lapply(x,elimna)
vals=lapply(x,unique)
vals=sort(elimna(unique(list2vec(vals))))
n=lapply(x,length)
n=list2vec(n)
K=length(vals)
J=length(x)
C1=matrix(0,nrow=K,ncol=J)
for(j in 1:J){
for(i in 1:K){
C1[i,j]=C1[i,j]+sum(x[[j]]==vals[i])
}
C1[,j]=C1[,j]/n[j]
}
test=0
for(i in 1:K)test=test+var(C1[i,])
list(test=test,C1=C1)
}

#' Nested ANOVA for Trimmed Means (Design A)
#'
#' @description
#' Performs a J-by-K nested ANOVA for trimmed means. Factor B is nested within
#' Factor A. For each level of Factor A, computes trimmed means for each level
#' of Factor B, then performs ANOVA on these trimmed means.
#'
#' @param x Data in list mode with length J. Each x[[j]] is a matrix with n_j
#'   rows and K columns, where rows are observations and columns are levels of
#'   Factor B nested within level j of Factor A.
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#'
#' @return Result from \code{t1way}: a list with test statistic, degrees of freedom,
#'   and p-value
#'
#' @details
#' This function implements a nested design where Factor B is nested within
#' Factor A. The strategy is: for each fixed level of Factor A (j=1,...,J),
#' compute the trimmed mean for each of the K levels of Factor B. These J*K
#' trimmed means become the units of analysis in a one-way ANOVA.
#'
#' Data structure: x must be a list of length J. Each element x[[j]] is a matrix
#' where rows are observations and columns correspond to the K levels of Factor B.
#'
#' @seealso \code{\link{anova.nestAP}}, \code{\link{t1way}}
#'
#' @export
anova.nestA<-function(x,tr=.2){
if(!is.list(x))stop('x should have list mode')
y=list()
J=length(x)
for(j in 1:J)y[[j]]=apply(x[[j]],2,tmean,tr=tr)
res=t1way(y,tr=tr)
res
}

#' Nested ANOVA for Trimmed Means with Pooled Data (Design AP)
#'
#' @description
#' Performs a J-by-K nested ANOVA for trimmed means where Factor B is nested
#' within Factor A. Pools all data within each level of Factor A (across all
#' levels of Factor B) before performing ANOVA.
#'
#' @param x Data in list mode with length J. Each x[[j]] is a matrix with n_j
#'   rows and K columns, where rows are observations and columns are levels of
#'   Factor B nested within level j of Factor A.
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#'
#' @return Result from \code{t1way}: a list with test statistic, degrees of freedom,
#'   and p-value
#'
#' @details
#' Unlike \code{anova.nestA}, which computes trimmed means separately for each
#' level of Factor B, this function pools all observations within each level of
#' Factor A (vectorizing the matrix across all K levels of Factor B) and then
#' performs a one-way ANOVA on these J pooled groups.
#'
#' Data structure: x must be a list of length J. Each element x[[j]] is a matrix
#' that gets converted to a vector for the analysis.
#'
#' @seealso \code{\link{anova.nestA}}, \code{\link{t1way}}
#'
#' @export
anova.nestAP<-function(x,tr=.2){
if(!is.list(x))stop('x should have list mode')
y=list()
J=length(x)
for(j in 1:J)y[[j]]=as.vector(x[[j]])
res=t1way(y,tr=tr)
res
}

#' Two-Way ANOVA for Medians with Tied Values
#'
#' @description
#' Performs a J by K two-way ANOVA for medians using projection-based depth.
#' Handles tied values well, which commonly occur with median-based analyses.
#'
#' @param J Number of levels for Factor A
#' @param K Number of levels for Factor B
#' @param x Data in list mode (length J*K) or matrix (columns are groups)
#' @param alpha Significance level (default: 0.05)
#' @param nboot Number of bootstrap samples (default: 2000)
#' @param MC Logical. If TRUE, use parallel processing via mclapply (default: FALSE)
#'
#' @return A list with components:
#'   \item{Fac.A}{Results for Factor A main effect (from pbadepth)}
#'   \item{Fac.B}{Results for Factor B main effect (from pbadepth)}
#'   \item{Fac.AB}{Results for A*B interaction (from pbadepth)}
#'
#' @details
#' This function uses the Harrell-Davis estimator for medians with projection-based
#' depth testing. It generates contrasts for main effects and interaction using
#' \code{con2way}, then tests each effect using \code{pbadepth}.
#'
#' A warning is issued if any group has fewer than 20 observations, as Type I
#' error control may be compromised. Total number of groups must equal J*K.
#'
#' @seealso \code{\link{Q3anova}}, \code{\link{pbadepth}}, \code{\link{con2way}}
#'
#' @export
Q2anova<-function(J,K,x,alpha=.05,nboot=2000,MC=FALSE){
if(is.matrix(x)|| is.data.frame(x))x=listm(x)
if(J*K != length(x))stop('Total number of groups is not equal to JK')
chkcar=NA
for(j in 1:length(x))chkcar[j]=length(unique(x[[j]]))
if(min(chkcar)<20){
print('Warning: Sample size is less than')
print('20 for one more groups. Type I error might not be controlled')
}
con=con2way(J,K)
A=pbadepth(x,est=hd,con=con$conA,alpha=alpha,nboot=nboot,MC=MC)
B=pbadepth(x,est=hd,con=con$conB,alpha=alpha,nboot=nboot,MC=MC)
AB=pbadepth(x,est=hd,con=con$conAB,alpha=alpha,nboot=nboot,MC=MC)
list(Fac.A=A,Fac.B=B,Fac.AB=AB)
}

#' Three-Way ANOVA for Medians with Tied Values
#'
#' @description
#' Performs a J by K by L three-way ANOVA for medians using projection-based depth.
#' Handles tied values well, which commonly occur with median-based analyses.
#'
#' @param J Number of levels for Factor A
#' @param K Number of levels for Factor B
#' @param L Number of levels for Factor C
#' @param x Data in list mode (length J*K*L) or matrix (columns are groups)
#' @param alpha Significance level (default: 0.05)
#' @param nboot Number of bootstrap samples (default: 600)
#' @param MC Logical. If TRUE, use parallel processing via mclapply (default: FALSE)
#'
#' @return A list with components:
#'   \item{Fac.A}{Results for Factor A main effect}
#'   \item{Fac.B}{Results for Factor B main effect}
#'   \item{Fac.C}{Results for Factor C main effect}
#'   \item{Fac.AB}{Results for A*B interaction}
#'   \item{Fac.AC}{Results for A*C interaction}
#'   \item{Fac.BC}{Results for B*C interaction}
#'   \item{Fac.ABC}{Results for A*B*C three-way interaction}
#'
#' @details
#' This function extends \code{Q2anova} to three-way designs. It uses the
#' Harrell-Davis estimator for medians with projection-based depth testing.
#' Contrasts for all main effects and interactions are generated using
#' \code{con3way}, then tested using \code{pbadepth}.
#'
#' A warning is issued if any group has fewer than 20 observations, as Type I
#' error control may be compromised. Total number of groups must equal J*K*L.
#'
#' @seealso \code{\link{Q2anova}}, \code{\link{pbadepth}}, \code{\link{con3way}}
#'
#' @export
Q3anova<-function(J,K,L,x,alpha=.05,nboot=600,MC=FALSE){
if(is.matrix(x)|| is.data.frame(x))x=listm(x)
if(J*K*L != length(x))stop('Total number of groups is not equal to JKL')
chkcar=NA
for(j in 1:length(x))chkcar[j]=length(unique(x[[j]]))
if(min(chkcar)<20){
print('Warning: Sample size is less than')
print('20 for one more groups. Type I error might not be controlled')
}
con=con3way(J,K,L)
A=pbadepth(x,est=hd,con=con$conA,alpha=alpha,nboot=nboot,MC=MC)
B=pbadepth(x,est=hd,con=con$conB,alpha=alpha,nboot=nboot,MC=MC)
C=pbadepth(x,est=hd,con=con$conC,alpha=alpha,nboot=nboot,MC=MC)
AB=pbadepth(x,est=hd,con=con$conAB,alpha=alpha,nboot=nboot,MC=MC)
AC=pbadepth(x,est=hd,con=con$conAC,alpha=alpha,nboot=nboot,MC=MC)
BC=pbadepth(x,est=hd,con=con$conBC,alpha=alpha,nboot=nboot,MC=MC)
ABC=pbadepth(x,est=hd,con=con$conABC,alpha=alpha,nboot=nboot,MC=MC)
list(Fac.A=A,Fac.B=B,Fac.C=C,Fac.AB=AB,Fac.AC=AC,Fac.BC=BC,Fac.ABC=ABC)
}

#' Brunner-Dette-Munk Test for One-Way Design
#'
#' @description
#' Tests the null hypothesis from Brunner et al. (2016) for a one-way design
#' using the rankFD package. Tests for differences in distributions across
#' independent groups.
#'
#' @param x Data in list mode (each element is a group) or matrix (columns are groups)
#'
#' @return A list with components:
#'   \item{test.stat}{ANOVA-type test statistic}
#'   \item{df1}{Numerator degrees of freedom}
#'   \item{df2}{Denominator degrees of freedom}
#'   \item{p.value}{p-value for the test}
#'   \item{q.hat}{Estimated relative effects from \code{bdm}}
#'
#' @details
#' This function implements the method from Brunner et al. (2016) for testing
#' the null hypothesis H0p (purely nonparametric hypothesis). It uses rank-based
#' methods via the rankFD package.
#'
#' The function converts list data to a data frame with group indicators and
#' calls \code{rankFD} with hypothesis='H0p'. It also computes relative effects
#' using the \code{bdm} function.
#'
#' @references
#' Brunner et al. (2016). [Complete reference for the BDM method]
#'
#' @seealso \code{\link{bdm}}
#'
#' @export
bdmP<-function(x){
if(is.matrix(x))x=listm(x)
J=length(x)
library(rankFD)
nv=lapply(x,length)
nv=as.vector(matl(nv))
z=x[[1]]
g=rep(1,nv[1])
for(j in 2:J){
z=c(z,x[[j]])
g=c(g,rep(j,nv[j]))
}
xg=cbind(z,g)
xg=as.data.frame(xg)
res=rankFD(z~g,data=xg,hypothesis = 'H0p')
w=as.vector(res$ANOVA.Type.Statistic)
list(test.stat=w[1],df1=w[2],df2=w[3],p.value=w[4],q.hat=bdm(x)$q.hat)
}

#' Bootstrap Brown-Forsythe ANOVA
#'
#' @description
#' Performs one-way ANOVA using a bootstrap version of the Brown-Forsythe test.
#' Provides robust inference for heteroscedastic data.
#'
#' @param x Data in list mode (each element is a group), matrix, or data frame
#' @param nboot Number of bootstrap samples (default: 1000)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#'
#' @return A list with components:
#'   \item{test.stat}{Brown-Forsythe test statistic}
#'   \item{p.value}{Bootstrap p-value}
#'
#' @details
#' This function implements a bootstrap version of the Brown-Forsythe ANOVA,
#' which is robust to heterogeneous variances. Data are centered by their group
#' means, and bootstrap samples are drawn from the centered distributions.
#'
#' The test statistic is computed using \code{BFANOVA}. Bootstrap resampling
#' generates the null distribution, and the p-value is the proportion of bootstrap
#' statistics at least as large as the observed statistic.
#'
#' Missing values are automatically removed.
#'
#' @seealso \code{\link{BFANOVA}}, \code{\link{t1way}}
#'
#' @export
BFBANOVA<-function(x,nboot=1000,SEED=TRUE){
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
J=length(x)
for(j in 1:J)x[[j]]=elimna(x[[j]])
TV=BFANOVA(x)$test.statistic
ylist<-list()
dat=list()
TT<-NA
#means<-sapply(x,mean)
for (j in 1:J)ylist[[j]]<-x[[j]]-mean(x[[j]])
for (i in 1:nboot){
for(j in 1:J)dat[[j]]=sample(ylist[[j]],length(ylist[[j]]),replace=TRUE)
TT[i]<-BFANOVA(dat)$test.statistic
}
pval<-mean(TV<=TT,na.rm=TRUE)
list(test.stat=TV,p.value=pval)
}

#' Brown-Forsythe ANOVA
#'
#' @description
#' Performs Brown-Forsythe one-way ANOVA, which is robust to heterogeneous
#' variances. Can be generalized to trimmed means (though performance of trimmed
#' version is not well established).
#'
#' @param x Data in list mode (each element is a group), matrix, or data frame
#'
#' @return A list with components:
#'   \item{test.statistic}{Brown-Forsythe F-statistic}
#'   \item{df1}{Numerator degrees of freedom}
#'   \item{df2}{Denominator degrees of freedom}
#'   \item{p.value}{p-value from F-distribution}
#'
#' @details
#' The Brown-Forsythe test is a modification of the classical ANOVA F-test that
#' provides better control of Type I error rates when variances are heterogeneous.
#' It uses weighted group means with weights inversely proportional to group
#' variances.
#'
#' The test statistic is: FS = sum(n*(ybar-YB)^2) / sum((1-n/N)*v)
#' where n are group sizes, ybar are group means, YB is the grand mean, v are
#' group variances, and N is total sample size.
#'
#' Note: While this can be generalized to trimmed means, the performance of such
#' generalization is not well established.
#'
#' Missing values are automatically removed.
#'
#' @seealso \code{\link{BFBANOVA}}, \code{\link{t1way}}
#'
#' @export
BFANOVA<-function(x){
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
J=length(x)
for(j in 1:J)x[[j]]=elimna(x[[j]])
xall=elimna(as.vector(matl(x)))
YB=mean(xall)
ybar=lapply(x,mean)
v=lapply(x,var)
n=lapply(x,length)
ybar=matl(ybar)
n=as.vector(matl(n))
v=as.vector(matl(v))
ybar=as.vector(matl(ybar))
w=n/v
N=sum(n)
top=sum(n*(ybar-YB)^2)
bot=sum((1-n/N)*v)
FS=top/bot
df1=J-1
fv=NA
for(j in 1:J)fv[j]=(1-n[j]/N)*v[j]
df2=1/sum(fv^2/(n-1))
pv=1-pf(FS,df1,df2)
list(test.statistic=FS,df1=df1,df2=df2,p.value=pv)
}

#' Confidence Interval for Explanatory Effect Size in One-Way ANOVA
#'
#' @description
#' Computes a bootstrap confidence interval for the explanatory measure of effect
#' size in one-way ANOVA. Provides interval estimate for the proportion of
#' variance explained.
#'
#' @param x Data in list mode (each element is a group), matrix, or data frame
#' @param alpha Significance level for confidence interval (default: 0.05 for 95% CI)
#' @param tr Proportion to trim from each end (default: 0, no trimming)
#' @param nboot Number of bootstrap samples (default: 500)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#' @param ITER Number of iterations for unequal sample sizes (default: 5).
#'   Passed to yuenv2 for iteration to get estimate.
#' @param adj Logical. If TRUE (default), apply bias correction adjustment for
#'   number of groups (available for J <= 8)
#' @param ... Additional arguments (currently unused)
#'
#' @return A list with components:
#'   \item{Effect.Size}{Point estimate of effect size (possibly adjusted)}
#'   \item{ci}{Vector of length 2 containing lower and upper confidence limits}
#'
#' @details
#' This function uses bootstrap resampling to construct a confidence interval
#' for the explanatory effect size from \code{t1wayv2}. If the omnibus test is
#' not significant at level alpha, the lower confidence limit is set to 0.
#'
#' When adj=TRUE and J <= 8, a bias correction factor is applied based on the
#' number of groups. The correction factors are: c(1, 1.268757, 1.467181,
#' 1.628221, 1.763191, 1.856621, 1.993326) for J-1 = 1 through 7.
#'
#' Missing values are automatically removed.
#'
#' @seealso \code{\link{t1wayv2}}, \code{\link{KS.ANOVA.ES}}
#'
#' @export
t1way.EXES.ci<-function(x,alpha=.05,tr=0,nboot=500,SEED=TRUE,ITER=5,adj=TRUE,...){
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x=listm(x)
J=length(x)
n=lapply(x,length)
if(SEED)set.seed(2)
chk=t1wayv2(x,tr=tr,SEED=FALSE)
v=list()
val=NA
x=elimna(x)
for(i in 1:nboot){
for(j in 1:J)v[[j]]=sample(x[[j]],replace=TRUE)
val[i]=t1wayv2(v,tr=tr,nboot=ITER,SEED=FALSE)$Effect.Size
}
ilow<-round((alpha/2) * nboot)
ihi<-nboot - ilow
ilow<-ilow+1
val=sort(val)
ci=val[ilow]
ci[2]=val[ihi]
if(chk$p.value>alpha)ci[1]=0
if(adj){
fix=c(1, 1.268757, 1.467181, 1.628221, 1.763191, 1.856621, 1.993326)
if(J>8)print('No adjustment available when J>8')
J1=J-1
if(j<=8){
chk$Effect.Size=fix[J1]*chk$Effect.Size
ci=fix[J1]*ci
}
}
list(Effect.Size=chk$Effect.Size,ci=ci)
}

#' Kulinskaya-Staudte Effect Size for One-Way ANOVA
#'
#' @description
#' Computes a measure of effect size for J independent groups based on a variation
#' of the Kulinskaya and Staudte (2005) method. Provides a standardized measure
#' of group differences.
#'
#' @param x Data in list mode (each element is a group), matrix, or data frame
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param adj Logical. If TRUE (default), rescale so that under a location shift
#'   in one group only (assuming normality), the magnitude approximately matches
#'   the magnitude when J=2 groups
#'
#' @return Numeric effect size value
#'
#' @details
#' The effect size is computed as: es = (J/2) * sqrt(sum(q*(m-gmean)^2/v))
#' where q are sample size proportions, m are trimmed means, gmean is the grand
#' mean, and v are winsorized variances.
#'
#' When adj=TRUE and J <= 8, a scaling adjustment is applied to make effect sizes
#' comparable across different numbers of groups. The adjustment factors are:
#' c(1, 0.3832, 0.3086, 0.2689, 0.2427, 0.2230, 0.2067) for J-1 = 1 through 7.
#'
#' Missing values are automatically removed.
#'
#' @references
#' Kulinskaya, E., & Staudte, R. G. (2005). [Complete reference]
#'
#' @seealso \code{\link{t1way.EXES.ci}}, \code{\link{anova.KMS.ND}}
#'
#' @export
KS.ANOVA.ES<-function(x,tr=0.2,adj=TRUE){
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x=listm(x)
J=length(x)
for(j in 1:J)x[[j]]=elimna(x[[j]])
n=lapply(x,length)
n=pool.a.list(n)
q=n/sum(n)
v=pool.a.list(lapply(x,winvarN,tr=tr))  #standard error
m=pool.a.list(lapply(x,tmean,tr))
gmean=sum(q*m)/sum(q)
term=q*(m-gmean)^2/v
es=sqrt(sum(term))
es=J*es/2
fix=c(1, 0.3832,  0.3086, 0.2689,  0.2427,  0.2230,  0.2067)
J1=J-1
if(adj){
if(J>8)print('No adjustment available when J>8')
else es=fix[J1]*es
}
es
}

#' Generate Null Distribution for KMS ANOVA
#'
#' @description
#' Generates the null distribution for the Kulinskaya-Staudte (KMS) ANOVA effect
#' size measure. Used for hypothesis testing with \code{anova.KMS}.
#'
#' @param n Vector of sample sizes for each group
#' @param tr Proportion to trim from each end (default: 0.2, meaning 20% from each tail)
#' @param iter Number of iterations for simulation (default: 5000)
#' @param nulldist Placeholder parameter for null distribution (currently unused in function body)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#'
#' @return Numeric vector of length iter containing simulated effect sizes under
#'   the null hypothesis
#'
#' @details
#' This function simulates data under the null hypothesis of no group differences
#' by generating samples from a g-and-h distribution (with g=0.75) for each group.
#' For each iteration, it computes the KMS effect size using \code{KS.ANOVA.ES}.
#'
#' The resulting vector can be used to determine critical values or p-values for
#' observed KMS effect sizes. The g-and-h distribution with g=0.75 provides a
#' flexible distribution for the null.
#'
#' @seealso \code{\link{KS.ANOVA.ES}}, \code{\link{anova.KMS}}
#'
#' @export
anova.KMS.ND<-function(n,tr=.2,iter=5000,nulldist=nulldist,SEED=TRUE){
if(SEED)set.seed(2)
v=NA
dat=list()
J=length(n)  # number of groups
for(i in 1:iter){
for(L in 1:J)dat[[L]]=ghdist(n[L],g=0.75)
v[i]=KS.ANOVA.ES(dat,tr=tr)
}
v
}

#' Plot KMS Effect Size Curve for 2x2 Design with Covariate
#'
#' @description
#' For a 2-by-2 design with a covariate, plots the KMS measure of effect size
#' associated with the two levels of the first factor as a function of the
#' covariate. Visualizes interaction effects.
#'
#' @param x List of length 4 containing covariate values for each of the four groups
#' @param y List of length 4 containing response values for each of the four groups
#' @param pts Vector of covariate values at which to evaluate effect size
#'   (NULL = use deciles from 10th to 90th percentile)
#' @param SW Logical. If TRUE, switch rows and columns (default: FALSE)
#' @param npts Number of points to use if pts is NULL (default: 15)
#' @param xlab Label for x-axis (default: 'X')
#' @param ylab Label for y-axis (default: 'Effect.Size')
#'
#' @return Creates a plot (no return value)
#'
#' @details
#' This function compares KMS effect sizes for the two levels of the first factor
#' in a 2x2 design, conditional on covariate values. For each point on the covariate,
#' it computes effect sizes using \code{ancova.KMS} and plots the difference.
#'
#' The function expects exactly 4 groups in both x and y. Missing values are
#' removed pairwise. When SW=TRUE, the groups are reordered as c(1,3,2,4) to
#' switch the factor levels.
#'
#' The plot shows v1-v2 where v1 and v2 are effect sizes for the two levels of
#' Factor A.
#'
#' @seealso \code{\link{t2way.KMS.interbt}}, \code{\link{ancova.KMS}}
#'
#' @export
t2way.KMS.curve<-function(x,y,pts=NULL,SW=FALSE,npts=15,xlab='X',ylab='Effect.Size'){

if(is.matrix(x) || is.data.frame(x))x=listm(x)
if(is.matrix(y) || is.data.frame(y))y=listm(y)
if(length(x)!=4)stop('Should have four groups exactly. Fix argument x')
if(length(y)!=4)stop('Should have four groups exactly. Fix argument y')
for(j in 1:4){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1]
y[[j]]=xy[,2]
}
n=lapply(y,length)
n=as.vector(matl(n))
adj=4.4/n+1.00086
flag=n>150
adj[flag]=1
adj=mean(adj)
nmax=max(n)
if(is.null(pts)){
xlow=max(matl((lapply(x,qest,.1))))
xhi=min(matl((lapply(x,qest,.9))))
pts=seq(xlow,xhi,length=npts)
}
nv=lapply(x,length)
if(SW){
x=x[c(1,3,2,4)]
y=y[c(1,3,2,4)]
}
v1=ancova.KMS(x[[1]],y[[1]],x[[2]],y[[2]],pts=pts,plotit=FALSE)[,2]
v2=ancova.KMS(x[[3]],y[[3]],x[[4]],y[[4]],pts=pts,plotit=FALSE)[,2]
v=v1-v2
plot(pts,v,xlab=xlab,ylab=ylab,type='n')
lines(pts,v)
}

#' Bootstrap Test for KMS Interaction Effect in 2x2 Design with Covariate
#'
#' @description
#' For a 2-by-2 design with a covariate, tests for interaction effects by comparing
#' KMS effect sizes associated with the two levels of the first factor using
#' bootstrap resampling.
#'
#' @param x List of length 4 containing covariate values for each of the four groups
#' @param y List of length 4 containing response values for each of the four groups
#' @param pts Vector of covariate values at which to evaluate effect size
#'   (NULL = use deciles from 10th to 90th percentile)
#' @param alpha Significance level (default: 0.05)
#' @param nboot Number of bootstrap samples (default: 100)
#' @param MC Logical. If TRUE, use parallel processing via mclapply (default: FALSE)
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility
#' @param SW Logical. If TRUE, switch rows and columns (default: FALSE)
#'
#' @return A list with components:
#'   \item{test}{Test statistic}
#'   \item{p.value}{Bootstrap p-value}
#'   \item{crit.value}{Bootstrap critical value at alpha level}
#'   \item{pts}{Covariate points used in analysis}
#'
#' @details
#' This function tests for interactions in a 2x2 design with a covariate by
#' comparing KMS effect sizes. For each covariate value, it computes the difference
#' in effect sizes between the two levels of Factor A using \code{ancova.KMS}.
#'
#' Bootstrap resampling generates the null distribution. The function expects
#' exactly 4 groups in both x and y. Missing values are removed pairwise.
#' When SW=TRUE, groups are reordered as c(1,3,2,4).
#'
#' The test uses the helper function \code{t2way.KMS.inter.sub} for bootstrap
#' computations.
#'
#' @seealso \code{\link{t2way.KMS.curve}}, \code{\link{t2way.KMS.inter.sub}}, \code{\link{ancova.KMS}}
#'
#' @export
t2way.KMS.interbt<-function(x,y,pts=NULL,alpha=.05,nboot=100,MC=FALSE,SEED=TRUE,SW=FALSE){

if(SEED)set.seed(2)
if(is.matrix(x) || is.data.frame(x))x=listm(x)
if(is.matrix(y) || is.data.frame(y))y=listm(y)
if(length(x)!=4)stop('Should have four groups exactly. Fix argument x')
if(length(y)!=4)stop('Should have four groups exactly. Fix argument y')
for(j in 1:4){
xy=elimna(cbind(x[[j]],y[[j]]))
x[[j]]=xy[,1]
y[[j]]=xy[,2]
}
n=lapply(y,length)
nn=as.vector(matl(n))
adj=4.4/nn+1.00086
flag=nn>150
adj[flag]=1
adj=mean(adj)
nmax=max(nn)
if(is.null(pts)){
xlow=max(matl((lapply(x,qest,.1))))
xhi=min(matl((lapply(x,qest,.9))))
pts=c(xlow,(xlow+xhi)/2,xhi)
}
nv=lapply(x,length)
if(SW){
x=x[c(1,3,2,4)]
y=y[c(1,3,2,4)]
}
npts=length(pts)
MAT1=matrix(NA,nmax,4)
MAT2=matrix(NA,nmax,4)
dat1=list()
dat2=list()
for(i in 1:nboot){
id1=sample(nv[[1]],replace=TRUE)
id2=sample(nv[[2]],replace=TRUE)
MAT1[1:nv[[1]],1:2]=cbind(x[[1]][id1],y[[1]][id1])
MAT1[1:nv[[2]],3:4]=cbind(x[[2]][id2],y[[2]][id2])
dat1[[i]]=MAT1
id1=sample(nv[[3]],replace=TRUE)
id2=sample(nv[[4]],replace=TRUE)
MAT2[1:nv[[3]],1:2]=cbind(x[[3]][id1],y[[3]][id1])
MAT2[1:nv[[4]],3:4]=cbind(x[[4]][id2],y[[4]][id2])
dat2[[i]]=MAT2
}
if(MC){
a1=mclapply(dat1,t2way.KMS.inter.sub,pts=pts)
a2=mclapply(dat2,t2way.KMS.inter.sub,pts=pts)
}
if(!MC){
a1=lapply(dat1,t2way.KMS.inter.sub,pts=pts)
a2=lapply(dat2,t2way.KMS.inter.sub,pts=pts)
}
a1=t(matl(a1))
a2=t(matl(a2))

zq=qnorm(1-alpha/2)
sqse1=NA
sqse2=NA
for(j in 1:npts)sqse1[j]=var(a1[,j])
for(j in 1:npts)sqse2[j]=var(a2[,j])
Results=matrix(NA,npts,10)
Results[,1:2]=ancova.KMS(x[[1]],y[[1]],x[[2]],y[[2]],pts=pts,plotit=FALSE)
Results[,3]=ancova.KMS(x[[3]],y[[3]],x[[4]],y[[4]],pts=pts,plotit=FALSE)[,2]
Results[,4]=Results[,2]-Results[,3]
Results[,5]=adj*Results[,4]/sqrt(sqse1+sqse2)       # n=40 1.15 work well get .054
pv=2*(1-pnorm(abs(Results[,5])))
Results[,6]=pv
Results[,7]=Results[,4]-zq*sqrt(sqse1+sqse2)/adj
Results[,8]=Results[,4]+zq*sqrt(sqse1+sqse2)/adj
Results[,9]=sqrt(sqse1+sqse2)/adj
Results[,10]=p.adjust(pv,method='hoch')
dimnames(Results)=list(NULL,c('pts','Est1','Est2','Dif','Test.Stat','p.value','ci.low','ci.up','SE','p.adjusted'))
n=matl(nv)
n=as.vector(n)
list(n=nn,Results=Results)
}

#' Helper Function for Two-Way KMS Interaction Bootstrap
#'
#' @description
#' Internal helper function for \code{t2way.KMS.interbt} that extracts KMS effect
#' sizes from ANCOVA analysis. Used in bootstrap computations.
#'
#' @param z Matrix with 4 columns: covariate and response for two groups
#' @param pts Vector of covariate values at which to evaluate effect size
#'
#' @return Numeric vector of effect size estimates (second column from ancova.KMS results)
#'
#' @details
#' This function is called internally by \code{t2way.KMS.interbt} during bootstrap
#' resampling. It expects z to be a matrix where columns 1 and 2 contain the
#' covariate and response for the first group, and columns 3 and 4 contain the
#' covariate and response for the second group.
#'
#' The function calls \code{ancova.KMS} with plotit=FALSE and extracts the second
#' column of results, which contains the effect size estimates.
#'
#' @seealso \code{\link{t2way.KMS.interbt}}, \code{\link{ancova.KMS}}
#'
#' @export
t2way.KMS.inter.sub<-function(z,pts){
ancova.KMS(z[,1],z[,2],z[,3],z[,4],pts=pts,plotit=FALSE)[,2]
}

