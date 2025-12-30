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

t1way.effect<-function(x,tr=.2,grp=NA,MAT=FALSE,lev.col=1,var.col=2){
#
# Same as t1way, but computes explanatory power and related effect size
# Only use this function with = n's
# Use t1wayv2 in general, which calls this function when sample sizes differ.
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
#  Missing values are automatically removed.
#
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

t1wayF<-function(x,fac,tr=.2,nboot=100,SEED=TRUE,EP=FALSE,pr=TRUE){
#
# Same a t1way, but now the data are assumed to be
# stored in a matrix or data frame where one of the columns contain
# the data to be analyzed and another column contains the group
# identification.
#
# For example, if dat is a data frame, with column 1 containing
# the outcome measures of interest, and column 2 is a factor variable
# indicating to  which group a value in column 1 belongs, then
#  t1wayF(dat[,1],dat[,2])
#  will test the hypothesis that all J groups have identical
#  trimmed means.
#
#  Missing values are automatically removed.
#
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

t1waybt<-function(x,tr=.2,grp=NA,nboot=599,SEED=TRUE){
#
#   Test the hypothesis of equal trimmed means, corresponding to J independent
#   groups, using a bootstrap-t bootstrap method.
#
#   The data are assumed to be stored in x in list mode
#   or in a matrix. In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, columns correspond to groups.
#
#   grp is used to specify some subset of the groups, if desired.
#   By default, all J groups are used.
#
#   The default number of bootstrap samples is nboot=599
#

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

t1waybtv2<-function(x,tr=.2,grp=NA,g=NULL,dp=NULL,nboot=599,SEED=TRUE){
#
#   Test the hypothesis of equal trimmed mdeans, corresponding to J independent
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

t2wayv2<-function(J,K,data,tr=.2,grp=c(1:p),p=J*K,g=NULL,dp=NULL,pr=TRUE){
#  Perform a J by K  (two-way) anova on trimmed means where
#  all groups are independent.
#
#  The R variable data is assumed to contain the raw
#  data stored in list mode, or a matrix with columns
#  corresponding to groups. If stored in list mode, data[[1]] contains the data
#  for the first level of all three factors: level 1,1,.
#  data[[2]] is assumed to contain the data for level 1 of the
#  first factor and level 2 of the second factor: level 1,2
#
#  The default amount of trimming is tr=.2
#
#  It is assumed that data has length JK, the total number of
#  groups being tested.
#
#   g=NULL, x is assumed to be a matrix or have list mode
#
#   if g is specifed, it is assumed that column g of x is
#   a factor variable and that the dependent variable of interest is in column
#   dp of x, which can be a matrix or data frame.
#
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

sppbi<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),nboot=500,SEED=TRUE,pr=TRUE,...){
#
# A percentile bootstrap for interactions
# in a split-plot design.
# The analysis is done by taking difference scores
# among all pairs of dependent groups and seeing whether
# these differences differ across levels of Factor A.
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
#
#  When in list mode x is assumed to have length JK, the total number of
#  groups being tested, but a subset of the data can be analyzed
#  using grp
#
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

sppba<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),avg=TRUE,nboot=500,SEED=TRUE,
MC=FALSE,MDIS=FALSE,pr=TRUE,...){
#
# A percentile bootstrap for main effects
# among independent groups in a split-plot design
#
# avg=T: The analysis is done by averaging K measures of
# location for  each level of Factor A,
# and then comparing averages by testing the hypothesis
# that all pairwise differences are equal to zero.
#
# avg=F:   The analysis is done by testing whether $K$ equalities are
# simultaneously true. For kth level of Factor B, the kth equality is
# theta_{1k}= ... theta_{Jk}, k=1,...,K.
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

YYmanova<-function(x1,x2,tr=.2){
#
# Do MANOVA using generalization of
# Yanagihara, H. \& Yuan, K. H. (2005).
# Three approximate solutions to the
#  multivariate Behrens-Fisher problem.  Communications in Statistics--
# Simulation and Computation, 34, 975--988; see their eq. (2.7).
#
#  x1 and x2 are assumed to be matrices
#
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

bwtrimbt<-function(J,K,x,tr=.2,JK=J*K,grp=c(1:JK),nboot=599,
SEED=TRUE){
#
# A bootstrap-t for performing a split-plot design
# with trimmed means.
# By default, 20% trimming is used with B=599 bootstrap samples.
#
#
#  The R variable x is assumed to contain the raw
#  data stored in list mode or in a matrix or a data frame
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
if(is.data.frame(x) || is.matrix(x)) {
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
btest<-tsplit(J,K,blist,tr)
testmat[iboot,1]<-btest$Qa
testmat[iboot,2]<-btest$Qb
testmat[iboot,3]<-btest$Qab
}
test=tsplit(J,K,x,tr=tr)
pbA=mean(test$Qa[1]<testmat[,1])
pbB=mean(test$Qb[1]<testmat[,2])
pbAB=mean(test$Qab[1]<testmat[,3])
list(p.value.A=pbA,p.value.B=pbB,p.value.AB=pbAB)
}

MULtr.anova<-function(x,J=NULL,p=NULL,tr=.2,alpha=.05){
#
#  Do Multivariate ANOVA with trimmed means using
#  Johansen's method
#
#  x is assumed to have list mode with length(x)=J=number of groups and
#  x[[j]] is an n_j-by-p  matrix, p is the number of variables.
#
#  x can also be a matrix when J and p are specified. It is assumed the data are stored in
#   a matrix in the same manner expected by bwtrim.
#
#  To get a p-value, use the function MULAOVp
#
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

t1wayv2<-function(x,tr=.2,grp=NA,MAT=FALSE,lev.col=1,var.col=2,nboot=100,SEED=TRUE,pr=TRUE,IV=NULL,loc.fun=median){
#
# Same a t1way, but computes explanatory power and related effect size
#
# For n1!=n2, this function calls t1way.effect.
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
#  IV, if specified, taken to be the independent variable
#      That is, the group id values
#      and x is assumed to be a vector containing all of the data
#
# if x has list mode:
#  length(x) is assumed to correspond to the total number of groups.
#  By default, the null hypothesis is that all groups have a common mean.
#  To compare a subset of the groups, use grp to indicate which
#  groups are to be compared. For example, if you type the
#  command grp<-c(1,3,4), and then execute this function, groups
#  1, 3, and 4 will be compared with the remaining groups ignored.
#
#  Missing values are automatically removed.
#
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

t1waybtsqrk<-function(x,alpha=.05,nboot=599,tr=0.2,SEED=TRUE){
#
#  One-way ANOVA for trimmed means, independent groups.
#  Uses a method studied by Ozdemir et al.
#
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

discANOVA<-function(x,nboot=500,SEED=TRUE){
#
#  Test the global hypothesis that for two or more independent groups,
#  the corresponding discrete distributions are identical.
#  That is, test the hypothesis that independent groups have identical
#  multinomial distributions. A generalization of the Storer--Kim method is used.
#
#  Could also use a chi-squared test via the function: disc2.chi.sq
#
#  The method is designed for situations where the
#  sample size is relatively small. The method can be sensitive to
#  differences that are missed using a measure of location.
#
#  Control over the Type I error probability is excellent, even when n=10
#
#  x is a matrix with n rows and J columns
#  or it can have list mode.
#
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

discANOVA.sub<-function(x){
#
#
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

anova.nestA<-function(x,tr=.2){
#
# J-by-K  nested ANOVA
#  x is assumed to have list mode with length J.
#  x[[j]] is assumed to be a matrix with n_j rows and K columns
#   j=1, ..., J
#
# Strategy: For fixed level of factor A compute trimmed mean for each
# level of factor B and use these trimmed means as the unit of analysis
#  That is, perform an ANOVA using these trimmed means
#
if(!is.list(x))stop('x should have list mode')
y=list()
J=length(x)
for(j in 1:J)y[[j]]=apply(x[[j]],2,tmean,tr=tr)
res=t1way(y,tr=tr)
res
}

anova.nestAP<-function(x,tr=.2){
#
# J-by-K  nested ANOVA
#  x is assumed to have list mode with length J.
#  x[[j]] is assumed to be a matrix with n_j rows and K columns
#   j=1, ..., J
#
#  pool data for each level of A and do anova
#
if(!is.list(x))stop('x should have list mode')
y=list()
J=length(x)
for(j in 1:J)y[[j]]=as.vector(x[[j]])
res=t1way(y,tr=tr)
res
}

Q2anova<-function(J,K,x,alpha=.05,nboot=2000,MC=FALSE){
#
#  Two-way ANOVA for medians, tied values allowed.
#
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

Q3anova<-function(J,K,L,x,alpha=.05,nboot=600,MC=FALSE){
#
#  Three-way ANOVA for medians, tied values allowed.
#
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

bdmP<-function(x){
#
# Test the null hypothesis in Bruner et al. (2016)
# for a one-way design
#
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

BFBANOVA<-function(x,nboot=1000,SEED=TRUE){
#
#  One-way ANOVA bootstrap version of Brown-Forsyhte
#
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

BFANOVA<-function(x){
#
#  Brown-Forsyhte ANOVA, generalized to trimmed means
#
#  (Not known whether a generalization to trimmed means performs relatively well.)
#
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

t1way.EXES.ci<-function(x,alpha=.05,tr=0,nboot=500,SEED=TRUE,ITER=5,adj=TRUE,...){
#
# Confidence interval for explanatory measure of effect size
#
#  ITER:  yuenv2, for unequal sample sizes.  iterates to get estimate
#
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

KS.ANOVA.ES<-function(x,tr=0.2,adj=TRUE){
#
#  J independent groups. Compute a measure of effect size
#  based on a slight variation of a method in
#  Kulinskaya and Staudte (2005)
#
# adj=TRUE:rescale so that under a shift of location in one group only, assuming normality, the
# magnitude matches, approximately, the magnitude when there are J=2 groups
#
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

anova.KMS.ND<-function(n,tr=.2,iter=5000,nulldist=nulldist,SEED=TRUE){
#
# Null distribution for anova.KMS
#
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

t2way.KMS.curve<-function(x,y,pts=NULL,SW=FALSE,npts=15,xlab='X',ylab='Effect.Size'){
#
# For a 2-by-2 design, compare
# KMS measure of effect size  associated with the two levels of the first factor
#    plots  an interaction effect when there is a covariate.
#
#  SW=TRUE, switches rows and column

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

t2way.KMS.interbt<-function(x,y,pts=NULL,alpha=.05,nboot=100,MC=FALSE,SEED=TRUE,SW=FALSE){
#
# For a 2-by-2 design, compare
# KMS measure of effect size  associated with the two levels of the first factor
#    to get an interaction effect when there is a covariate.
#
#  SW=TRUE, switches rows and column

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

t2way.KMS.inter.sub<-function(z,pts){
ancova.KMS(z[,1],z[,2],z[,3],z[,4],pts=pts,plotit=FALSE)[,2]
}

