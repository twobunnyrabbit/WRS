# mcp.R
#
# Multiple Comparisons Procedures (MCP) Functions
# Part of WRS package refactoring
#
# This module contains functions for:
# - Contrast matrix generators (con1way, con2way, con3way)
# - Linear contrasts for independent groups (lincon*, linconb, linconpb, linconbt)
# - Linear contrasts for dependent groups (lindep*, pairdepb)
# - MCP for factorial designs (mcp2a, mcp2atm, mcp3atm, mcp3med, rm3mcp)
# - Between-within designs (bwmcp, bwwmcp, bbwmcp, bwrmcp, bwimcp, bwbmcp)
# - Bootstrap MCP (mcppb, tmcppb, bmcppb, pbmcp, bbmcppb, bwmcppb, wwmcppb, etc.)
# - Split-plot MCP (spmcpa, spmcpi, spmcpb, sppba, sppbb, sppbi)
# - Quantile-based MCP (qdmcp, qdmcpdif)
# - MCP with effect sizes (ES, EP, QS variants)
# - Specialized MCP (stepmcp, signmcp, discmcp, sintmcp, anctsmcp, skmcp)
# - P-value adjustment utilities (mcpPV, mcpKadjp)
# - Nested and other designs (mcp.nestAP, binmcp, etc.)
#
# Extracted: 2025-12-30
# Number of functions: 102 main functions + 5 helper functions (.sub)

################################################################################
# HELPER FUNCTIONS
################################################################################

tsub<-function(isub,x,y,tr){
#
#  Compute test statistic for trimmed means
#  when comparing dependent groups.
#  By default, 20% trimmed means are used.
#  isub is a vector of length n,
#  a bootstrap sample from the sequence of integers
#  1, 2, 3, ..., n
#
#  This function is used by ydbt
#
tsub<-yuend(x[isub],y[isub],tr=tr)$teststat
tsub
}


################################################################################
# CONTRAST MATRIX GENERATORS
################################################################################

con1way<-function(J){
#
#   Create contrast coefficients for all pairwise comparisons
Ja=(J^2-J)/2
con<-matrix(0,J,Ja)
id<-0
for (j in 1:J){
for(k in 1:J){
if(j<k){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
con
}


linconm<-function(x,con=0,est=onestep,alpha=.05,nboot=500,pr=TRUE,...){
#
#   Compute a 1-alpha confidence interval for a set of d linear contrasts
#   involving M-estimators using a bootstrap method. (See Chapter 6.)
#   Independent groups are assumed.
#
#   The data are assumed to be stored in x in list mode.  Thus,
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J, say.
#
#   con is a J by d matrix containing the contrast coefficents of interest.
#   If unspecified, all pairwise comparisons are performed.
#   For example, con[,1]=c(1,1,-1,-1,0,0) and con[,2]=c(,1,-1,0,0,1,-1)
#   will test two contrasts: (1) the sum of the first two measures of location is
#   equal to the sum of the second two, and (2) the difference between
#   the first two is equal to the difference between the measure of location for
#   groups 5 and 6.
#
#   The default number of bootstrap samples is nboot=399
#
#   This function uses the function trimpartt written for this
#   book.
#
#
#
#
if(pr){
print("Note: confidence intervals are adjusted to control FWE")
print("But p-values are not adjusted to control FWE")
}
if(is.matrix(x))x<-listm(x)
con<-as.matrix(con)
if(!is.list(x))stop("Data must be stored in list mode.")
J<-length(x)
Jm<-J-1
d<-(J^2-J)/2
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
if(nrow(con)!=length(x))stop("The number of groups does not match the number of contrast coefficients.")
m1<-matrix(0,J,nboot)
m2<-1 # Initialize m2
mval<-1
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:J){
mval[j]<-est(x[[j]],...)
xcen<-x[[j]]-est(x[[j]],...)
data<-matrix(sample(xcen,size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
m1[j,]<-apply(data,1,est,...) # A J by nboot matrix.
m2[j]<-var(m1[j,])
}
boot<-matrix(0,ncol(con),nboot)
bot<-1
for (d in 1:ncol(con)){
top<-apply(m1,2,trimpartt,con[,d])
#            A vector of length nboot containing psi hat values
consq<-con[,d]^2
bot[d]<-trimpartt(m2,consq)
boot[d,]<-abs(top)/sqrt(bot[d])
}
testb<-apply(boot,2,max)
ic<-floor((1-alpha)*nboot)
testb<-sort(testb)
psihat<-matrix(0,ncol(con),6)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper","se","p.value"))
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,2]<-trimpartt(mval,con[,d])
psihat[d,3]<-psihat[d,2]-testb[ic]*sqrt(bot[d])
psihat[d,4]<-psihat[d,2]+testb[ic]*sqrt(bot[d])
psihat[d,5]<-sqrt(bot[d])
pval<-mean((boot[d,]<abs(psihat[d,2])/psihat[d,5]))
psihat[d,6]<-1-pval
}
list(psihat=psihat,crit=testb[ic],con=con)
}


lincon.old<-function(x,con=0,tr=.2,alpha=.05,pr=TRUE,crit=NA,SEED=TRUE,KB=FALSE){
#
#  A heteroscedastic test of d linear contrasts using trimmed means.
#
#  The data are assumed to be stored in $x$ in list mode, a matrix
#  or a data frame. If in list mode,
#  length(x) is assumed to correspond to the total number of groups.
#  It is assumed all groups are independent.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  If con is not specified, all pairwise comparisons are made.
#
#  Missing values are automatically removed.
#
#  To apply the Kaiser-Bowden method, use the function kbcon
#
if(tr==.5)stop("Use the R function medpb to compare medians")
if(is.data.frame(x))x=as.matrix(x)
if(KB)stop("Use the function kbcon")
flag<-TRUE
if(alpha!= .05 && alpha!=.01)flag<-FALSE
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-length(x)
sam=NA
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
sam[j]=length(x[[j]])
h[j]<-length(x[[j]])-2*floor(tr*length(x[[j]]))
   # h is the number of observations in the jth group after trimming.
w[j]<-((length(x[[j]])-1)*winvar(x[[j]],tr))/(h[j]*(h[j]-1))
xbar[j]<-mean(x[[j]],tr)
}
if(sum(con^2)==0){
CC<-(J^2-J)/2
if(CC>28)print("For faster execution time but less power, use kbcon")
psihat<-matrix(0,CC,8)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper",
"p.value","Est.1","Est.2"))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c("Group","Group","test","crit","se","df"))
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
df<-(w[j]+w[k])^2/(w[j]^2/(h[j]-1)+w[k]^2/(h[k]-1))
test[jcom,6]<-df
psihat[jcom,6]<-2*(1-pt(test[jcom,3],df))
psihat[jcom,7]=xbar[j]
psihat[jcom,8]=xbar[k]
if(!KB){
if(CC>28)flag=FALSE
if(flag){
if(alpha==.05)crit<-smmcrit(df,CC)
if(alpha==.01)crit<-smmcrit01(df,CC)
}
if(!flag || CC>28)crit<-smmvalv2(dfvec=rep(df,CC),alpha=alpha,SEED=SEED)
}
if(KB)crit<-sqrt((J-1)*(1+(J-2)/df)*qf(1-alpha,J-1,df))
test[jcom,4]<-crit
psihat[jcom,4]<-(xbar[j]-xbar[k])-crit*sejk
psihat[jcom,5]<-(xbar[j]-xbar[k])+crit*sejk
}}}}
if(sum(con^2)>0){
if(nrow(con)!=length(x)){
stop("The number of groups does not match the number of contrast coefficients.")
}
psihat<-matrix(0,ncol(con),5)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper",
"p.value"))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c("con.num","test","crit","se","df"))
df<-0
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-sqrt(sum(con[,d]^2*w))
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
df<-(sum(con[,d]^2*w))^2/sum(con[,d]^4*w^2/(h-1))
if(flag){
if(alpha==.05)crit<-smmcrit(df,ncol(con))
if(alpha==.01)crit<-smmcrit01(df,ncol(con))
}
if(!flag)crit<-smmvalv2(dfvec=rep(df,ncol(con)),alpha=alpha,SEED=SEED)
test[d,3]<-crit
test[d,4]<-sejk
test[d,5]<-df
psihat[d,3]<-psihat[d,2]-crit*sejk
psihat[d,4]<-psihat[d,2]+crit*sejk
psihat[d,5]<-2*(1-pt(abs(test[d,2]),df))
}
}
if(pr){
print("Note: confidence intervals are adjusted to control FWE")
print("But p-values are not adjusted to control FWE")
print('Adjusted p-values can be computed with the R function p.adjust')
}
list(n=sam,test=test,psihat=psihat)
}


lincon.pool<-function(x,con=0,tr=.2,alpha=.05,POOL=FALSE){
#
#  Same as lincon but with a pooling option that is used when
#  dealing with main effects in a two-way and three-way designs
#
#  See, for example, the function twowayA.poolB
#

if(tr==.5)stop('Use the R function medpb to compare medians')
if(is.data.frame(x))x=as.matrix(x)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
if(sum(con^2)>0){
if(POOL){
ic=0
y=list()
nc=ncol(con)
nc2=nc*2
Ncon=matrix(0,nrow=nc2,ncol=nc)
for(k in 1:nc){
id1=which(con[,k]==1)
id2=which(con[,k]==-1)
ic=ic+1
print(ic)
Ncon[ic,k]=1
y[[ic]]=pool.a.list(x[id1])
ic=ic+1
Ncon[ic,k]=-1
y[[ic]]=pool.a.list(x[id2])
}
res=lincon(y,con=Ncon,tr=tr)
print(Ncon)
}}
if(!POOL)res=lincon(x,con=con,tr=tr,alpha=alpha)
res
}


linconb<-function(x,con=0,tr=.2,alpha=.05,nboot=599,pr=FALSE,SEED=TRUE,method='holm'){
#
#   Compute a 1-alpha confidence interval for a set of d linear contrasts
#   involving trimmed means using the bootstrap-t bootstrap method.
#   Independent groups are assumed.
#
#   The data are assumed to be stored in x in list mode.  Thus,
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J, say.
#
#   Missing values are automatically removed.
#
#   con is a J by d matrix containing the contrast coefficents of interest.
#   If unspecified, all pairwise comparisons are performed.
#   For example, con[,1]=c(1,1,-1,-1,0,0) and con[,2]=c(,1,-1,0,0,1,-1)
#   will test two contrasts: (1) the sum of the first two trimmed means is
#   equal to the sum of the second two, and (2) the difference between
#   the first two is equal to the difference between the trimmed means of
#   groups 5 and 6.
#
#   The default number of bootstrap samples is nboot=599
#
#   This function uses functions trimparts and trimpartt written for this
#   book.
#
#
#
#
if(is.data.frame(x))x=as.matrix(x)
con<-as.matrix(con)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
for(j in 1:J){
xx<-x[[j]]
x[[j]]<-xx[!is.na(xx)] # Remove any missing values.
}
Jm<-J-1
d<-(J^2-J)/2
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
if(nrow(con)!=length(x))stop('The number of groups does not match the number of contrast coefficients.')
bvec<-array(0,c(J,2,nboot))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
if(pr)print('Taking bootstrap samples. Please wait.')
nsam=matl(lapply(x,length))
for(j in 1:J){
paste('Working on group ',j)
xcen<-x[[j]]-mean(x[[j]],tr)
data<-matrix(sample(xcen,size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,,]<-apply(data,1,trimparts,tr) # A 2 by nboot matrix. The first row
#                     contains the bootstrap trimmed means, the second row
#                     contains the bootstrap squared standard errors.
}
m1<-bvec[,1,]  # J by nboot matrix containing the bootstrap trimmed means
m2<-bvec[,2,]  # J by nboot matrix containing the bootstrap sq. se.
boot<-matrix(0,ncol(con),nboot)
for (d in 1:ncol(con)){
top<-apply(m1,2,trimpartt,con[,d])
#            A vector of length nboot containing psi hat values
consq<-con[,d]^2
bot<-apply(m2,2,trimpartt,consq)
boot[d,]<-abs(top)/sqrt(bot)
}
testb<-apply(boot,2,max)
ic<-floor((1-alpha)*nboot)
testb<-sort(testb)
psihat<-matrix(0,ncol(con),4)
test<-matrix(0,ncol(con),5)
dimnames(psihat)<-list(NULL,c('con.num','psihat','ci.lower','ci.upper'))
dimnames(test)<-list(NULL,c('con.num','test','se','p.value','p.adjusted'))
for (d in 1:ncol(con)){
test[d,1]<-d
psihat[d,1]<-d
testit<-lincon(x,con[,d],tr,pr=FALSE)
test[d,2]<-testit$test[1,2]
pval<-mean((abs(testit$test[1,2])<boot[d,]))
test[d,4]<-pval
psihat[d,3]<-testit$psihat[1,2]-testb[ic]*testit$test[1,4]
psihat[d,4]<-testit$psihat[1,2]+testb[ic]*testit$test[1,4]
psihat[d,2]<-testit$psihat[1,2]
test[d,3]<-testit$test[1,4]
}
test[,5]=p.adjust(test[,4],method=method)
list(n=nsam,psihat=psihat,test=test,crit=testb[ic],con=con)
}


linconbt<-function(x,con=0,tr=.2,alpha=.05,nboot=599,pr=FALSE,SEED=TRUE,method='holm'){
#
#   Compute a 1-alpha confidence interval for a set of d linear contrasts
#   involving trimmed means using the bootstrap-t bootstrap method.
#   Independent groups are assumed.
#
#   The data are assumed to be stored in x in list mode.  Thus,
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J, say.
#
#   Missing values are automatically removed.
#
#   con is a J by d matrix containing the contrast coefficents of interest.
#   If unspecified, all pairwise comparisons are performed.
#   For example, con[,1]=c(1,1,-1,-1,0,0) and con[,2]=c(,1,-1,0,0,1,-1)
#   will test two contrasts: (1) the sum of the first two trimmed means is
#   equal to the sum of the second two, and (2) the difference between
#   the first two is equal to the difference between the trimmed means of
#   groups 5 and 6.
#
#   The default number of bootstrap samples is nboot=599
#
#   This function uses functions trimparts and trimpartt written for this
#   book.
#
#
#
#
if(is.data.frame(x))x=as.matrix(x)
con<-as.matrix(con)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
for(j in 1:J){
xx<-x[[j]]
x[[j]]<-xx[!is.na(xx)] # Remove any missing values.
}
Jm<-J-1
d<-(J^2-J)/2
FLAG=FALSE
if(sum(con^2)==0){
FLAG=TRUE
con=con.all.pairs(J)
}
if(nrow(con)!=length(x))stop('The number of groups does not match the number of contrast coefficients.')
bvec<-array(0,c(J,2,nboot))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
nsam=matl(lapply(x,length))
for(j in 1:J){
xcen<-x[[j]]-mean(x[[j]],tr)
data<-matrix(sample(xcen,size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,,]<-apply(data,1,trimparts,tr) # A 2 by nboot matrix. The first row
#                     contains the bootstrap trimmed means, the second row
#                     contains the bootstrap squared standard errors.
}
m1<-bvec[,1,]  # J by nboot matrix containing the bootstrap trimmed means
m2<-bvec[,2,]  # J by nboot matrix containing the bootstrap sq. se.
boot<-matrix(0,ncol(con),nboot)
for (d in 1:ncol(con)){
top<-apply(m1,2,trimpartt,con[,d])
#            A vector of length nboot containing psi hat values
consq<-con[,d]^2
bot<-apply(m2,2,trimpartt,consq)
boot[d,]<-abs(top)/sqrt(bot)
}
testb<-apply(boot,2,max)
ic<-floor((1-alpha)*nboot)
ic.crit=ic
testb<-sort(testb)
psihat<-matrix(0,ncol(con),4)
test<-matrix(0,ncol(con),5)
dimnames(psihat)<-list(NULL,c('con.num','psihat','ci.lower','ci.upper'))
dimnames(test)<-list(NULL,c('con.num','test','se','p.value','p.adjusted'))
for (d in 1:ncol(con)){
test[d,1]<-d
psihat[d,1]<-d
testit<-lincon(x,con[,d],tr,pr=FALSE)
test[d,2]<-testit$test[1,2]
pval<-mean((abs(testit$test[1,2])<boot[d,]))
test[d,4]<-pval
psihat[d,3]<-testit$psihat[1,2]-testb[ic]*testit$test[1,4]
psihat[d,4]<-testit$psihat[1,2]+testb[ic]*testit$test[1,4]
psihat[d,2]<-testit$psihat[1,2]
test[d,3]<-testit$test[1,4]
}
test[,5]=p.adjust(test[,4],method=method)
if(FLAG){
rem=psihat
Trem=test
# For all pairwise comparisons, adjust format of output to make it easier to read
CC<-(J^2-J)/2
psihat<-matrix(0,CC,9)
dimnames(psihat)<-list(NULL,c('Group','Group','Est.1-Est.2','ci.lower','ci.upper',
'p.value','Est.1','Est.2','adj.p.value'))
test<-matrix(NA,CC,4)
dimnames(test)<-list(NULL,c('Group','Group','test','se'))
ic=0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic=ic+1
psihat[ic,1]=j
psihat[ic,2]=k
psihat[ic,3:5]=rem[ic,2:4]
psihat[ic,6]=Trem[ic,4]
psihat[ic,7]=tmean(x[[j]],tr=tr)
psihat[ic,8]=tmean(x[[k]],tr=tr)
psihat[ic,9]=Trem[ic,5]
test[ic,1]=j
test[ic,2]=k
test[ic,3]=Trem[ic,2]
test[ic,4]=Trem[ic,3]
}}}
}
list(n=nsam,psihat=psihat,test=test,crit=testb[ic.crit],con=con)
}


linconpb<-function(x,alpha=.05,nboot=NA,grp=NA,est=tmean,con=0,method='holm',bhop=FALSE,SEED=TRUE,...){
#
#   Multiple comparisons for  J independent groups using trimmed means
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
#   Missing values are allowed.
#
con<-as.matrix(con)
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
if(bhop)method='BH'
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
if(nrow(con)!=J)stop('Something is wrong with con; the number of rows does not match the number of groups.')
#  Determine nboot if a value was not specified
if(is.na(nboot)){
nboot<-5000
if(J <= 8)nboot<-4000
if(J <= 3)nboot<-2000
}
# Determine critical values
if(method!='BH'){
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
if(method=='BH')dvec<-(ncon-c(1:ncon)+1)*alpha/ncon
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
dimnames(output)<-list(NULL,c('con.num','psihat','p.value','p.crit','ci.lower','ci.upper','p.adjusted'))
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
output[,7]=p.adjust(output[,3],method=method)
num.sig<-sum(output[,3]<=output[,4])
list(output=output,con=con,num.sig=num.sig)
}


linconMpb<-function(x,alpha=.05,nboot=1000,grp=NA,est=tmean,con=0,bhop=FALSE,
SEED=TRUE,PDIS=FALSE,J=NULL,p=NULL,...){
#
#   Multiple comparisons for  J independent groups using trimmed means
#   with multivariate data for each group.
#
#   A percentile bootstrap method with Rom's method is used.
#
#   The data are assumed to be stored in x
#   which  has list mode,
#   x[[1]] contains the data for the first group in the form of a
#   matrix, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#
#   est is the measure of location and defaults to the median
#   ... can be used to set optional arguments associated with est
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#   Missing values are automatically removed.
#
con<-as.matrix(con)
if(is.matrix(x) || is.data.frame(x)){
if(is.null(J) && is.null(p))stop("Specify J or P")
x=MAT2list(x,p=p,J=J)
}
if(!is.list(x))stop("Data must be stored in list mode.")
if(!is.na(sum(grp))){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
x<-xx
}
J<-length(x)
nullvec=rep(0,ncol(x[[1]]))
bplus=nboot+1
tempn<-0
mvec<-list
for(j in 1:J){
x[[j]]<-elimna(x[[j]])
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
if(nrow(con)!=J)stop("Something is wrong with con; the number of rows does not match the number of groups.")
# Determine critical levels
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
bvec<-array(NA,c(J,nboot,ncol(x[[1]])))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
nvec=lapply(x,nrow)
for(j in 1:J){
data<-matrix(sample(nvec[[j]],size=nvec[[j]]*nboot,replace=TRUE),nrow=nboot)
bvec[j,,]<-apply(data,1,linconMpb.sub,x[[j]],est,...) # Bootstrapped values for jth group
}
test<-NA
for (d in 1:ncon){
tv=matrix(0,nboot,ncol(x[[1]])) #nboot by p matrix reflecting Psi hat
estit=rep(0,ncol(x[[1]]))
for(j in 1:J){
tv=tv+con[j,d]*bvec[j,,]
estit=estit+con[j,d]*apply(x[[j]],2,est,...)
}
if(!PDIS)m1=cov(tv)
tv=rbind(tv,nullvec)
if(!PDIS)dv=mahalanobis(tv,center=estit,m1)
if(PDIS)dv=pdis(tv,center=estit) # projection distances
test[d]=1-sum(dv[bplus]>=dv[1:nboot])/nboot
}
output<-matrix(0,ncon,3)
dimnames(output)<-list(NULL,c("con.num","p.value","p.crit"))
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,3]<-zvec
for (ic in 1:ncol(con)){
output[ic,1]<-ic
output[ic,2]<-test[ic]
}
num.sig<-sum(output[,2]<=output[,3])
list(output=output,con=con,num.sig=num.sig)
}


linconSpb<-function(x,alpha=.05,nboot=1000,grp=NA,est=smean,con=0,bhop=FALSE,
SEED=TRUE,PDIS=FALSE,J=NULL,p=NULL,...){
#
#   Multiple comparisons for  J independent groups
#   with multivariate data for each group.
#   That is, linear contrasts relevant to MANOVA can be tested.
#   The method can handle
#   multivariate measures of location that take into account
#   the overall structure of the data, as opposed to using, for example
#   the marginal trimmed means, which is done by default when using
#   linconMpb.
#   The argument
#
#   est=smean,
#
#   means that by default the skipped measure of location, based on
#   on projection method for detecting outliers, is used.
#
#   Mahalanobis distance is used to compute a p-value, but projection
#   distances could be used by setting PDIS=T.
#
#   A percentile bootstrap method with Rom's method is used.
#
#   alpha=.05 means the probability of one or more type I errors is .05.
#
#   The data are assumed to be stored in x
#   which  has list mode,
#   x[[1]] contains the data for the first group in the form of a
#   matrix, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#
#   est is the measure of location and defaults to the median
#   ... can be used to set optional arguments associated with est
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#   Missing values are automatically removed.
#
if(is.matrix(x) || is.data.frame(x)){
if(is.null(J) && is.null(p))stop("Specify J or P")
x=MAT2list(x,p=p,J=J)
}
con<-as.matrix(con)
if(!is.list(x))stop("Data must be stored in list mode.")
if(!is.na(sum(grp))){  # Only analyze specified groups.
xx<-list()
for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
x<-xx
}
J<-length(x)
nullvec=rep(0,ncol(x[[1]]))
bplus=nboot+1
tempn<-0
mvec<-list
for(j in 1:J){
x[[j]]<-elimna(x[[j]])
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
if(nrow(con)!=J)stop("Something is wrong with con; the number of rows does not match the number of groups.")
# Determine critical levels
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
bvec<-array(NA,c(J,nboot,ncol(x[[1]])))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
nvec=lapply(x,nrow)
for(j in 1:J){
data<-matrix(sample(nvec[[j]],size=nvec[[j]]*nboot,replace=TRUE),nrow=nboot)
bvec[j,,]<-apply(data,1,linconSpb.sub,x[[j]],est,...) # Bootstrapped values for jth group
}
test<-NA
for (d in 1:ncon){
tv=matrix(0,nboot,ncol(x[[1]])) #nboot by p matrix reflecting Psi hat
estit=rep(0,ncol(x[[1]]))
for(j in 1:J){
tv=tv+con[j,d]*bvec[j,,]
estit=estit+con[j,d]*est(x[[j]],...)
}
if(!PDIS)m1=cov(tv)
tv=rbind(tv,nullvec)
if(!PDIS)dv=mahalanobis(tv,center=estit,m1)
if(PDIS)dv=pdis(tv,center=estit) # projection distances
test[d]=1-sum(dv[bplus]>=dv[1:nboot])/nboot
}
output<-matrix(0,ncon,3)
dimnames(output)<-list(NULL,c("con.num","p.value","p.crit"))
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,3]<-zvec
for (ic in 1:ncol(con)){
output[ic,1]<-ic
output[ic,2]<-test[ic]
}
num.sig<-sum(output[,2]<=output[,3])
list(output=output,con=con,num.sig=num.sig)
}


linconEP<-function(x,con=0,tr=.2,alpha=.05,pr=TRUE,crit=NA,SEED=TRUE,INT=FALSE,nreps=200,POOL=FALSE){
#
#
#  This function is used when  estimating effect size via
#  a variation of explanatory power.
#
#  It is restricted to the usual main effects and interactions in a two-way design.
# This function is used by bbmcpEP.
#
#  con: used to indicate main effects and is passed to this function via bbmcpEP
#
#  POOL=TRUE: For the usual main effects  in a two-way where
#  for a fixed level of Factor A, say, one can simply pool the data over the
#  levels of Factor A. POOL=TRUE means that data with contrast coefficients
#  = 1 are pooled, the same is for data with contrast coefficients
#  = -1 and the resulting two groups are compared.
#
#  A heteroscedastic test of d linear contrasts using trimmed means.
#
#  The data are assumed to be stored in $x$ in list mode, a matrix
#  or a data frame. If in list mode,
#  length(x) is assumed to correspond to the total number of groups.
#  It is assumed all groups are independent.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  If con is not specified, all pairwise comparisons are made.
#
#  Missing values are automatically removed.
#
#
if(tr==.5)stop('Use the R function medpb to compare medians')
if(is.data.frame(x))x=as.matrix(x)
flag<-TRUE
if(alpha!= .05 && alpha!=.01)flag<-FALSE
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
J<-length(x)
sam=NA
h<-vector('numeric',J)
w<-vector('numeric',J)
xbar<-vector('numeric',J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
sam[j]=length(x[[j]])
h[j]<-length(x[[j]])-2*floor(tr*length(x[[j]]))
   # h is the number of observations in the jth group after trimming.
w[j]<-((length(x[[j]])-1)*winvar(x[[j]],tr))/(h[j]*(h[j]-1))
xbar[j]<-mean(x[[j]],tr)
}
if(sum(con^2)==0){
CC<-(J^2-J)/2
psihat<-matrix(0,CC,7)
dimnames(psihat)<-list(NULL,c('Group','Group','psihat','ci.lower','ci.upper',
'p.value','Effect.Size'))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c('Group','Group','test','crit','se','df'))
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
df<-(w[j]+w[k])^2/(w[j]^2/(h[j]-1)+w[k]^2/(h[k]-1))
test[jcom,6]<-df
psihat[jcom,6]<-2*(1-pt(test[jcom,3],df))
psihat[jcom,7]=yuenv2(x[[j]],x[[k]])$Effect.Size
if(CC>28)flag=FALSE
if(flag){
if(alpha==.05)crit<-smmcrit(df,CC)
if(alpha==.01)crit<-smmcrit01(df,CC)
if(!flag || CC>28)crit<-smmvalv2(dfvec=rep(df,CC),alpha=alpha,SEED=SEED)
}
test[jcom,4]<-crit
psihat[jcom,4]<-(xbar[j]-xbar[k])-crit*sejk
psihat[jcom,5]<-(xbar[j]-xbar[k])+crit*sejk
}}}}
if(sum(con^2)>0){
if(nrow(con)!=length(x)){
stop('The number of groups does not match the number of contrast coefficients.')
}
psihat<-matrix(0,ncol(con),6)
dimnames(psihat)<-list(NULL,c('con.num','psihat','ci.lower','ci.upper',
'p.value','Effect.Size'))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c('con.num','test','crit','se','df'))
df<-0
for (d in 1:ncol(con)){
if(POOL){
id1=which(con[,d]==1)
id2=which(con[,d]==-1)
y1=pool.a.list(x[id1])
y2=pool.a.list(x[id2])
xx=list(y1,y2)
conP=matrix(c(1,-1))
}
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-sqrt(sum(con[,d]^2*w))
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
df<-(sum(con[,d]^2*w))^2/sum(con[,d]^4*w^2/(h-1))
if(flag){
if(alpha==.05)crit<-smmcrit(df,ncol(con))
if(alpha==.01)crit<-smmcrit01(df,ncol(con))
}
if(!flag)crit<-smmvalv2(dfvec=rep(df,ncol(con)),alpha=alpha,SEED=SEED)
test[d,3]<-crit
test[d,4]<-sejk
test[d,5]<-df
if(!POOL)temp=linEP(x,con[,d],tr=tr,nreps=nreps,SEED=SEED)
if(POOL)temp=linEP(xx,conP,tr=tr,nreps=nreps,SEED=SEED)
if(!INT){
psihat[d,6]=linEP(x,con[,d],tr=tr,nreps=nreps,SEED=SEED)$Effect.Size
}
if(INT){
id=con[,d]!=0
psihat[d,6]=Inter.EP(x[id],tr=tr,nreps=nreps,SEED=SEED)$Effect.Size
}
psihat[d,3]<-psihat[d,2]-crit*sejk
psihat[d,4]<-psihat[d,2]+crit*sejk
psihat[d,5]<-2*(1-pt(abs(test[d,2]),df))
}
}
list(n=sam,test=test,psihat=psihat)
}


linconES<-function(x,con=0,tr=.2,alpha=.05,pr=TRUE,crit=NA,SEED=TRUE,INT=FALSE,
locfun=tmean){
#
#  Like the function lincon, only
#  this function  estimates effect size via
#  quantile shift perspective.
#
#  A heteroscedastic test of d linear contrasts using trimmed means.
#
#  The data are assumed to be stored in $x$ in list mode, a matrix
#  or a data frame. If in list mode,
#  length(x) is assumed to correspond to the total number of groups.
#  It is assumed all groups are independent.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  If con is not specified, all pairwise comparisons are made.
#
#  Missing values are automatically removed.
#
#
if(tr==.5)stop('Use the R function medpb to compare medians')
if(is.data.frame(x))x=as.matrix(x)
flag<-TRUE
if(alpha!= .05 && alpha!=.01)flag<-FALSE
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
con<-as.matrix(con)
J<-length(x)
sam=NA
h<-vector('numeric',J)
w<-vector('numeric',J)
xbar<-vector('numeric',J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
sam[j]=length(x[[j]])
h[j]<-length(x[[j]])-2*floor(tr*length(x[[j]]))
   # h is the number of observations in the jth group after trimming.
w[j]<-((length(x[[j]])-1)*winvar(x[[j]],tr))/(h[j]*(h[j]-1))
xbar[j]<-mean(x[[j]],tr)
}
if(sum(con^2)==0){
CC<-(J^2-J)/2
psihat<-matrix(0,CC,8)
dimnames(psihat)<-list(NULL,c('Group','Group','psihat','ci.lower','ci.upper',
'p.value','Q.effect','Rel.Q'))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c('Group','Group','test','crit','se','df'))
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
df<-(w[j]+w[k])^2/(w[j]^2/(h[j]-1)+w[k]^2/(h[k]-1))
test[jcom,6]<-df
psihat[jcom,6]<-2*(1-pt(test[jcom,3],df))
psihat[jcom,7]=lin.ES(x[c(j,k)],con=c(1,-1))$Effect.Size
psihat[jcom,8]=(psihat[jcom,7]-.5)/.5
if(CC>28)flag=FALSE
if(flag){
if(alpha==.05)crit<-smmcrit(df,CC)
if(alpha==.01)crit<-smmcrit01(df,CC)
}
if(!flag || CC>28)crit<-smmvalv2(dfvec=rep(df,CC),alpha=alpha,SEED=SEED)
test[jcom,4]<-crit
psihat[jcom,4]<-(xbar[j]-xbar[k])-crit*sejk
psihat[jcom,5]<-(xbar[j]-xbar[k])+crit*sejk
}}}}
if(sum(con^2)>0){
if(nrow(con)!=length(x)){
stop('The number of groups does not match the number of contrast coefficients.')
}
psihat<-matrix(0,ncol(con),6)
dimnames(psihat)<-list(NULL,c('con.num','psihat','ci.lower','ci.upper',
'p.value','Q.effect'))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c('con.num','test','crit','se','df'))
df<-0
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-sqrt(sum(con[,d]^2*w))
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
df<-(sum(con[,d]^2*w))^2/sum(con[,d]^4*w^2/(h-1))
if(flag){
if(alpha==.05)crit<-smmcrit(df,ncol(con))
if(alpha==.01)crit<-smmcrit01(df,ncol(con))
}
if(!flag)crit<-smmvalv2(dfvec=rep(df,ncol(con)),alpha=alpha,SEED=SEED)
test[d,3]<-crit
test[d,4]<-sejk
test[d,5]<-df
if(!INT)psihat[d,6]=lin.ES(x,con[,d],tr=tr,nreps=nreps,SEED=SEED)$Effect.Size
if(INT){
id=con[,d]!=0
psihat[d,6]=interQS(x[id],nreps=nreps,locfun=locfun,SEED=SEED)$Q.Effect
}
psihat[d,3]<-psihat[d,2]-crit*sejk
psihat[d,4]<-psihat[d,2]+crit*sejk
psihat[d,5]<-2*(1-pt(abs(test[d,2]),df))
}
}
if(pr){
print('Note: confidence intervals are adjusted to control FWE')
print('But p-values are not adjusted to control FWE')
print('Adjusted p-values can be computed with the R function p.adjusted')
print('Under normality and homoscedasticity, Cohen d= .2, .5, .8')
print('corresponds approximately  to Rel.Q = 0.55, 0.65 and 0.70, respectively')
}
list(n=sam,test=test,psihat=psihat)
}


linconQS<-function(x,con=0,tr=.2,alpha=.05,pr=TRUE,crit=NA,SEED=TRUE,INT=FALSE,
locfun=tmean){
#
#
#  This function is used when  estimating effect size via
#  quantile shift perspective.
#
#  A heteroscedastic test of d linear contrasts using trimmed means.
#
#  The data are assumed to be stored in x in list mode, a matrix
#  or a data frame. If in list mode,
#  length(x) is assumed to correspond to the total number of groups.
#  It is assumed all groups are independent.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  If con is not specified, all pairwise comparisons are made.
#
#  Missing values are automatically removed.
#
#
if(tr==.5)stop('Use the R function medpb to compare medians')
if(is.data.frame(x))x=as.matrix(x)
flag<-TRUE
if(alpha!= .05 && alpha!=.01)flag<-FALSE
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
con<-as.matrix(con)
J<-length(x)
sam=NA
h<-vector('numeric',J)
w<-vector('numeric',J)
xbar<-vector('numeric',J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
sam[j]=length(x[[j]])
h[j]<-length(x[[j]])-2*floor(tr*length(x[[j]]))
   # h is the number of observations in the jth group after trimming.
w[j]<-((length(x[[j]])-1)*winvar(x[[j]],tr))/(h[j]*(h[j]-1))
xbar[j]<-mean(x[[j]],tr)
}
if(sum(con^2)==0){
CC<-(J^2-J)/2
psihat<-matrix(0,CC,8)
dimnames(psihat)<-list(NULL,c('Group','Group','psihat','ci.lower','ci.upper',
'p.value','Q.effect','Rel.Q'))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c('Group','Group','test','crit','se','df'))
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
df<-(w[j]+w[k])^2/(w[j]^2/(h[j]-1)+w[k]^2/(h[k]-1))
test[jcom,6]<-df
psihat[jcom,6]<-2*(1-pt(test[jcom,3],df))
psihat[jcom,7]=lin.ES(x[c(j,k)],con=c(1,-1))$Effect.Size
psihat[jcom,8]=(psihat[jcom,7]-.5)/.5
if(CC>28)flag=FALSE
if(flag){
if(alpha==.05)crit<-smmcrit(df,CC)
if(alpha==.01)crit<-smmcrit01(df,CC)
}
if(!flag || CC>28)crit<-smmvalv2(dfvec=rep(df,CC),alpha=alpha,SEED=SEED)
test[jcom,4]<-crit
psihat[jcom,4]<-(xbar[j]-xbar[k])-crit*sejk
psihat[jcom,5]<-(xbar[j]-xbar[k])+crit*sejk
}}}}
if(sum(con^2)>0){
if(nrow(con)!=length(x)){
stop('The number of groups does not match the number of contrast coefficients.')
}
psihat<-matrix(0,ncol(con),6)
dimnames(psihat)<-list(NULL,c('con.num','psihat','ci.lower','ci.upper',
'p.value','Q.effect'))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c('con.num','test','crit','se','df'))
df<-0
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-sqrt(sum(con[,d]^2*w))
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
df<-(sum(con[,d]^2*w))^2/sum(con[,d]^4*w^2/(h-1))
if(flag){
if(alpha==.05)crit<-smmcrit(df,ncol(con))
if(alpha==.01)crit<-smmcrit01(df,ncol(con))
}
if(!flag)crit<-smmvalv2(dfvec=rep(df,ncol(con)),alpha=alpha,SEED=SEED)
test[d,3]<-crit
test[d,4]<-sejk
test[d,5]<-df
if(!INT)psihat[d,6]=lin.ES(x,con[,d],tr=tr,nreps=nreps,SEED=SEED)$Effect.Size
if(INT){
id=con[,d]!=0
psihat[d,6]=interQS(x[id],nreps=nreps,locfun=locfun,SEED=SEED)$Q.Effect
}
psihat[d,3]<-psihat[d,2]-crit*sejk
psihat[d,4]<-psihat[d,2]+crit*sejk
psihat[d,5]<-2*(1-pt(abs(test[d,2]),df))
}
}
if(pr){
print('Note: confidence intervals are adjusted to control FWE')
print('But p-values are not adjusted to control FWE')
print('Adjusted p-values can be computed with the R function p.adjusted')
print('Under normality and homoscedasticity, Cohen d= 0, .2, .5, .8')
print('correspond approximately to Q.effect = 0.5, 0.55, 0.65 and 0.70, respectively')
}
list(n=sam,test=test,psihat=psihat)
}


lincon.bin<-function(r,n,con=NULL,alpha=.05,null.value=0,x=NULL,method='KMS',binCI=acbinomci){
#
#  r: number of successes for J independent groups
#  n: corresponding sample sizes
#
#  Compute confidence interval for a linear combination of independent binomials
#
# For 3 or more groups use:
# A note on confidence interval estimation for a linear function
# of binomial proportion.
#  Zou, G. Y., Huang, W. & Zheng, X (2009) CSDA, 53, 1080-1085
#  Otherwise, use binom2g with
#  method='KMS' or 'SK'
#
#  con: contrast coeffiients
#  if NULL, all pairwise comparisons are performed.
#
#  x: if not NULL, taken to be a matrix containing 0s and 1s, columns correspond to groups
#  r and n are computed using the data in  x
#
#  binCI, required by Zou et al. methods, defaults to Agresti--Coull
#  Other choices for binCI:
#  binomci:  Pratt's method
#  binomCP:  Clopper--Pearson
# kmsbinomci:  Kulinskaya et al
#  wilbinomci:  Wilson
#  binomLCO:  Schilling--Doi
#
if(!is.null(x)){
r=apply(x,2,sum)
n=rep(nrow(x),ncol(x))
}
J=length(r)
est=matrix(NA,nrow=J,ncol=3)
for(j in 1:J){
v=binCI(r[j],n[j],alpha=alpha)
est[j,]=c(v$phat,v$ci)
}
if(!is.null(con))con=as.matrix(con)
if(is.null(con))con=con.all.pairs(J)
NT=ncol(con)
L=NA
U=NA
EST=NA
PV=NA
for(k in 1:NT){
if(sum(con[,k]!=0)==2){
id1=which(con[,k]==1)
id2=which(con[,k]==-1)
a=binom2g(r[id1],n[id1],r[id2],n[id2],method=method)
L[k]=ifelse(method!='SK',a$ci[1],NA)
U[k]=ifelse(method!='SK',a$ci[2],NA )
PV[k]=a$p.value
EST[k]=a$est.dif
}
if(sum(con[,k]!=0)>2){
mat=cbind(con[,k]*est[,2],con[,k]*est[,3])
LM=apply(mat,1,min)
UM=apply(mat,1,max)
term1=sum(con[,k]*est[,1])
EST[k]=term1
term2=sqrt(sum((con[,k]*est[,1]-LM)^2))
term3=sqrt(sum((con[,k]*est[,1]-UM)^2))
L[k]=term1-term2
U[k]=term1+term3
PV[k]=lincon.binPV(r,n,con=con[,k],nullval=null.value,binCI=acbinomci)$p.value
}
}
adj=p.adjust(PV,method='hoch')
CI=cbind(EST,L,U,PV,adj)
dimnames(CI)=list(NULL,c('Est','ci.low','ci.hi','p-value','Adjusted p.value'))
list(p.hat=est[,1],CI=CI,con=con)
}


lincon.binPV<-function(r,n,con=NULL,alpha=.05,nullval=0,binCI=acbinomci){
#
#  Compare two binomials using the method in Zou et al.2009 CSDA.
#
#  x and y are vectors of 1s and 0s.
#  Or can use the argument
#  r1 = the number of successes observed among group 1
#  r2 = the number of successes observed among group 2
#  n1 = sample size for group 1
#  n2 = sample size for group 2
#
#  nullval is the hypothesized value of the linear contrast
#
#  binCI defaults to Agresti--Coull
#  Other choices for binCI:
#  binomci:  Pratt's method
#  binomCP:  Clopper--Pearson
# kmsbinomci:  Kulinskaya et al
#  wilbinomci:  Wilson
#  binomLCO: Schilling--Doi
#
ci=lincon.bin.sub(r=r,n=n,alpha=alpha,con=con,binCI=binCI)
p.value=1
p1.hat=r[1]/n1
p2.hat=r[2]/n2
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-lincon.bin.sub(r=r,n=n,alpha=alph[i],con=con,binCI=binCI)$CI[2:3]
if(chkit[1]>nullval || chkit[2]<nullval)break
}
p.value<-irem/100
iup<-(irem+1)/100
alph<-seq(.001,iup,.001)
for(i in 1:length(alph)){
irem<-i
p.value<-alph[i]
chkit<-lincon.bin.sub(r=r,n=n,alpha=alph[i],con=con,binCI=binCI)$CI[2:3]
if(chkit[1]>nullval || chkit[2]<nullval)break
}
list(n=n,p.est=ci$CI[1],ci=ci$CI[2:3],p.value=p.value)
}



################################################################################
# LINEAR CONTRASTS - DEPENDENT GROUPS
################################################################################

lindep<-function(x,con,cmat,alpha=.05,tr=.2){
#
#  Compute a test statistic based on the
#  linear contrast coefficients in con and the covariance matrix
#  cmat.
#
#  The data are assumed to be stored in x in list mode
#  or a matrix with columns correpsonding to groups.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  d=number of linear contrasts
#
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-length(x)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
for(j in 1:J){
xbar[j]<-mean(x[[j]],tr=tr,na.rm=TRUE)
}
ncon<-ncol(con)
psihat<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","se","test"))
w<-cmat
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
cvec<-as.matrix(con[,d])
sejk<-sqrt(t(cvec)%*%w%*%cvec)
psihat[d,3]<-sejk
psihat[d,4]<-psihat[d,2]/sejk
}
list(test.stat=psihat)
}


lindepbt<-function(x, con = NULL, tr = 0.2, alpha = 0.05,nboot=599,dif=TRUE,method='holm',
SEED=TRUE){
#
# MCP on trimmed means with FWE controlled with Rom's method
# Using a bootstrap-t method.
#
#  dif=T, difference scores are used. And for linear contrasts a simple
#  extension is used.
#
#  dif=F, hypotheses are tested based on the marginal trimmed means.
#
if(SEED)set.seed(2)
if(is.data.frame(x))x=as.matrix(x)
if(is.list(x))x=matl(x)
if(is.null(con))con=con.all.pairs(ncol(x))   # all pairwise
x=elimna(x)
n=nrow(x)
flagcon=FALSE
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-ncol(x)
xbar<-vector("numeric",J)
nval<-nrow(x)
h1<-nrow(x)-2*floor(tr*nrow(x))
df<-h1-1
xbar=apply(x,2,mean,tr=tr)
if(sum(con^2!=0))CC<-ncol(con)
ncon<-CC
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
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
if(nrow(con)!=ncol(x))warning("The number of groups does not match the number
 of contrast coefficients.")
ncon<-ncol(con)
psihat<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol(con),6)
dimnames(test)<-list(NULL,c("con.num","test","p.value","p.crit","se",'p.adjusted'))
temp1<-NA
for (d in 1:ncol(con)){
psihat[d,1]<-d
#
#  !dif  Use marginal trimmed means
#
if(!dif){
psihat[d,2]<-sum(con[,d]*xbar)
#
#
sejk<-0
for(j in 1:J){
for(k in 1:J){
djk<-(nval-1)*wincor(x[,j],x[,k], tr)$cov/(h1*(h1-1))
sejk<-sejk+con[j,d]*con[k,d]*djk
}}
sejk<-sqrt(sejk)
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
test[d,5]<-sejk
#
# now use bootstrap-t to determine p-value
#
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
xcen=x
for(j in 1:ncol(x))xcen[,j]=xcen[,j]-tmean(x[,j],tr=tr)
bvec=apply(data,1,lindep.sub,xcen,con[,d],tr)
bsort<-sort(abs(bvec))
ic<-round((1-alpha)*nboot)
ci<-0
psihat[d,3]<-psihat[d,2]-bsort[ic]*test[d,5]
psihat[d,4]<-psihat[d,2]+bsort[ic]*test[d,5]
p.value<-mean(abs(test[d,2])<=abs(bvec))
temp1[d]=p.value
}
if(dif){
for(j in 1:J){
if(j==1)dval<-con[j,d]*x[,j]
if(j>1)dval<-dval+con[j,d]*x[,j]
}
temp=trimcibt(dval,tr=tr,alpha=alpha,nboot=nboot,pr=FALSE)
temp1[d]<-temp$p.value #trimci(dval,tr=tr,pr=FALSE)$p.value
test[d,1]<-d
test[d,2]=temp$test.stat
test[d,5]<-trimse(dval,tr=tr)
psihat[d,2]<-mean(dval,tr=tr)
psihat[d,3]<-temp$ci[1] #psihat[,2]-qt(1-test[,4]/2,df)*test[,5]
psihat[d,4]<-temp$ci[2] #psihat[,2]+qt(1-test[,4]/2,df)*test[,5]
}}
#
#   d ends here
#
test[,3]<-temp1
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2,3]>=zvec)
test[temp2,4]<-zvec
if(flagcon)num.sig<-sum(test[,4]<=test[,5])
if(!flagcon){num.sig<-sum(test[,3]<=test[,4])
test[,6]=p.adjust(test[,3],method=method)
}
list(test=test,psihat=psihat,con=con,num.sig=num.sig)
}


pairdepb<-function(x,tr=.2,alpha=.05,grp=0,nboot=599){
#
#   Using the percentile t bootstrap method,
#   compute a .95 confidence interval for all pairwise differences between
#   the trimmed means of dependent groups.
#   By default, 20% trimming is used with B=599 bootstrap samples.
#
#   x can be an n by J matrix or it can have list mode
#
if(is.data.frame(x)) x <- as.matrix(x)
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
if(sum(grp)==0)grp<-c(1:length(x))
# put the data in an n by J matrix
mat<-matrix(0,length(x[[1]]),length(grp))
for (j in 1:length(grp))mat[,j]<-x[[grp[j]]]
}
if(is.matrix(x)){
if(sum(grp)==0)grp<-c(1:ncol(x))
mat<-x[,grp]
}
if(sum(is.na(mat)>=1))stop("Missing values are not allowed.")
J<-ncol(mat)
connum<-(J^2-J)/2
bvec<-matrix(0,connum,nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(nrow(mat),size=nrow(mat)*nboot,replace=TRUE),nrow=nboot)
xcen<-matrix(0,nrow(mat),ncol(mat))
for (j in 1:J)xcen[,j]<-mat[,j]-mean(mat[,j],tr) #Center data
it<-0
for (j in 1:J){
for (k in 1:J){
if(j<k){
it<-it+1
bvec[it,]<-apply(data,1,tsub,xcen[,j],xcen[,k],tr)
# bvec is a connum by nboot matrix containing the bootstrap test statistics.
}}}
bvec<-abs(bvec)  #Doing two-sided confidence intervals
icrit<-round((1-alpha)*nboot)
critvec<-apply(bvec,2,max)
critvec<-sort(critvec)
crit<-critvec[icrit]
psihat<-matrix(0,connum,5)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
test<-matrix(NA,connum,4)
dimnames(test)<-list(NULL,c("Group","Group","test","se"))
it<-0
for (j in 1:J){
for (k in 1:J){
if(j<k){
it<-it+1
estse<-yuend(mat[,j],mat[,k])$se
dif<-mean(mat[,j],tr)-mean(mat[,k],tr)
psihat[it,1]<-grp[j]
psihat[it,2]<-grp[k]
psihat[it,3]<-dif
psihat[it,4]<-dif-crit*estse
psihat[it,5]<-dif+crit*estse
test[it,1]<-grp[j]
test[it,2]<-grp[k]
test[it,3]<-yuend(mat[,j],mat[,k])$teststat
test[it,4]<-estse
}}}
list(test=test,psihat=psihat,crit=crit)
}



################################################################################
# CORE MCP FUNCTIONS - ONE-WAY DESIGNS
################################################################################

mcppb<-function(x,crit=NA,con=0,tr=.2,alpha=.05,nboot=2000,grp=NA,WIN=FALSE,
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
#   Or the data can be stored in a matrix with J columns
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
if(is.na(crit) && tr != .2)stop("A critical value must be specified when
the amount of trimming differs from .2")
if(WIN){
if(tr < .2)warning("When Winsorizing, the amount of trimming should be at least
.2")
if(win > tr)stop("Amount of Winsorizing must <= amount of trimming")
if(min(tempn) < 15){warning("Winsorizing with sample sizes less than 15 can")
warning(" result in poor control over the probability of a Type I error")
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
if(d> 10 && nboot <5000)warning("Suggest using nboot=5000 when the number
of contrasts exceeds 10.")
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
"ci.upper","p.value"))
if(nrow(con)!=length(x))stop("The number of groups does not match the number
 of  contrast coefficients.")
bvec<-matrix(NA,nrow=J,ncol=nboot)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:J){
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,mean,tr) # Bootstrapped trimmed means for jth group
}
test<-NA
for (d in 1:ncol(con)){
top<-0
for (i in 1:J){
top<-top+con[i,d]*bvec[i,]
}
test[d]<-sum((top>0))/nboot
test[d]<-min(test[d],1-test[d])
top<-sort(top)
psihat[d,4]<-top[icl]
psihat[d,5]<-top[icu]
}
for (d in 1:ncol(con)){
psihat[d,1]<-d
testit<-lincon(x,con[,d],tr,pr=FALSE)
psihat[d,6]<-test[d]
psihat[d,2]<-testit$psihat[1,2]
psihat[d,3]<-testit$test[1,4]
}
print("Reminder: To control FWE, reject if the p-value is less than")
print("the crit.p.value listed in the output.")
list(psihat=psihat,crit.p.value=crit,con=con)
}


tmcppb<-function(x,alpha=.05,nboot=NA,grp=NA,est=tmean,con=0,bhop=FALSE,SEED=TRUE,...){
#
#   Multiple comparisons for  J independent groups using trimmed means
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
#   Missing values are allowed.
#
con<-as.matrix(con)
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


bmcppb<-function(x,alpha=.05,nboot=NA,grp=NA,est=tmean,con=0,bhop=FALSE,SEED=TRUE,
...){
#
#   Multiple comparisons for  J independent groups using trimmed means
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
#   Missing values are allowed.
#
con<-as.matrix(con)
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


pbmcp<-function(x,alpha=.05,nboot=NA,grp=NA,est=onestep,con=0,bhop=FALSE,
SEED=TRUE,...){
#
#   Multiple comparisons for  J independent groups.
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   est is the measure of location and defaults to an M-estimator
#   ... can be used to set optional arguments associated with est
#
#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#   Missing values are allowed.
#
okay=FALSE
if(identical(est,onestep))okay=TRUE
if(identical(est,mom))okay=TRUE
if(!okay)stop('For estimators other than onestep and mom, use linconpb')
con<-as.matrix(con)
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
nmax=max(tempn)
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
if(nrow(con)!=J){
stop("Something is wrong with con; the number of rows does not match the number of groups.")
}
#  Determine nboot if a value was not specified
if(is.na(nboot)){
nboot<-5000
if(J <= 8)nboot<-4000
if(J <= 3)nboot<-2000
}
# Determine critical values
if(!bhop){
if(!identical(est,onestep))print('When est is not equal to  onestep, suggest using bhop=TRUE')
if(alpha==.05){
dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(nmax>=100)dvec[1]=.01
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
dvec[1]<-alpha/2
}
dvec<-2*dvec
}
if(nmax>80){
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
#paste("Working on group ",j)
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,est,...) # Bootstrapped values for jth group
}
chkna=sum(is.na(bvec))
if(chkna>0){
print("Bootstrap estimates of location could not be computed")
print("This can occur when using an M-estimator")
print("Might try est=tmean")
}
test<-NA
bcon<-t(con)%*%bvec #ncon by nboot matrix
tvec<-t(con)%*%mvec
for (d in 1:ncon){
test[d]<-(sum(bcon[d,]>0)+.5*sum(bcon[d,]==0))/nboot
if(test[d]> .5)test[d]<-1-test[d]
}
test<-2*test
output<-matrix(0,ncon,6)
dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.crit",
"ci.lower","ci.upper"))
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


pbtrmcp<-function(x,alpha=.05,nboot=NA,grp=NA,con=0,bhop=FALSE,tr=.2,SEED=TRUE){
#
#   Multiple comparisons for  J independent groups based on trimmed means.
#   using a percentile bootstrap method
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#

#   The argument grp can be used to analyze a subset of the groups
#   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
#
#   Missing values are allowed.
#
stop('Old function for trimmed means. Use bmcppb. (The function tmcppb gives the same results as bmcppb)')
con<-as.matrix(con)
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
mvec[j]<-tmean(temp,tr=tr)
}
nmax=max(tempn)
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
dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
dvec[1]<-alpha/2
}
dvec<-2*dvec
}
if(bhop)dvec<-(ncon-c(1:ncon)+1)*alpha/ncon
bvec<-matrix(NA,nrow=J,ncol=nboot)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
for(j in 1:J){
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,tmean,tr=tr) # Bootstrapped values for jth group
}
test<-NA
bcon<-t(con)%*%bvec #ncon by nboot matrix
tvec<-t(con)%*%mvec
for (d in 1:ncon){
test[d]<-sum(bcon[d,]>0)/nboot
if(test[d]> .5)test[d]<-1-test[d]
}
test<-2*test
output<-matrix(0,ncon,6)
dimnames(output)<-list(NULL,c("con.num","psihat","sig.test","sig.crit","ci.lower","ci.upper"))
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


mcpOV<-function(x,alpha=.05,nboot=500,grp=NA,est=smean,con=0,bhop=FALSE,SEED=TRUE,
...){
#
#   Multiple comparisons for  J independent groups using
#   using a multivariate O-type measure of location
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
#   Missing values are allowed.
#
if(is.data.frame(x))x=as.matrix(x)
con<-as.matrix(con)
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
Jm<-J-1
x=elimna(matl(x))
n=nrow(x)
mvec=est(x)
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
if(nrow(con)!=J)stop("Something is wrong with con; the number of rows does not match the number of groups.")
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
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(i in 1:nboot)bvec[,i]=est(x[data[i,],])
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



################################################################################
# TWO-WAY FACTORIAL DESIGNS
################################################################################

mcp2a<-function(J,K,x,est=mom,con=NULL,alpha=.05,nboot=NA,grp=NA,...){
#
# Do all pairwise comparisons of
# main effects for Factor A and B and all interactions
#
        #  The data are assumed to be stored in x
        #  in list mode or in a matrix.
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
        if(identical(est,median))print('Warning: med2mcp is a better when using the usual sample median')
        JK <- J * K
        if(is.matrix(x))
                x <- listm(x)
        if(!is.na(grp)) {
                yy <- x
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop("Data must be stored in list mode or a matrix.")
mvec<-NA
  tempn=0
        for(j in 1:JK) {
                xx <- x[[j]]
                x[[j]] <- xx[!is.na(xx)]
                mvec[j]<-est(x[[j]],...)
tempn[j]=length(x[[j]])
        }
nmax=max(tempn)
        #
        # Create the three contrast matrices
        #
        if(JK != length(x))
                warning("The number of groups does not match the number of contrast coefficients.")
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#  Determine nboot if a value was not specified
if(is.na(nboot)){
nboot<-5000
if(J <= 8)nboot<-4000
if(J <= 3)nboot<-2000
}
bvec<-matrix(NA,nrow=JK,ncol=nboot)
for(j in 1:JK){
data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,]<-apply(data,1,est,...) # J by nboot matrix, jth row contains
#                          bootstrapped  estimates for jth group
}
outvec<-list()
if(!is.null(con))stop('Use linconm when specifying the linear contrast coefficients')
temp3<-con2way(J,K)
for(jj in 1:3){
con<-temp3[[jj]]
con<-as.matrix(con)
ncon<-ncol(con)
# Determine critical values
if(alpha==.05){
dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(nmax>80){
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
}}
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
test<-NA
bcon<-t(con)%*%bvec #ncon by nboot matrix
tvec<-t(con)%*%mvec
for (d in 1:ncon){
test[d]<-sum(bcon[d,]>0)/nboot
if(test[d]> .5)test[d]<-1-test[d]
}
output<-matrix(0,ncon,6)
dimnames(output)<-list(NULL,c("con.num","psihat","sig.test","sig.crit","ci.lower","ci.upper"))
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
zvec[ddd:ncon]<-dvec[ddd]
}
output[temp2,4]<-zvec
icl<-round(dvec[ncon]*nboot)+1
icu<-nboot-icl-1
for (ic in 1:ncol(con)){
output[ic,2]<-tvec[ic,]
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(bcon[ic,])
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
outvec[[jj]]<-output
}
list(FactorA=outvec[[1]],FactorB=outvec[[2]],Interactions=outvec[[3]],
conA=temp3[[1]],conB=temp3[[2]],conAB=temp3[[3]])
}


mcp2atm<-function(J,K,x,tr=.2,alpha=.05,grp=NA,op=FALSE,pr=TRUE){
#
#  Test all linear contrasts associated with
# main effects for Factor A and B and all interactions based on trimmed means
# By default,
# tr=.2, meaning 20% trimming is used.
#
#   bbmcpEP has an option for pooling over the levels of the factors.
#
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
        JK <- J * K
        if(is.matrix(x))
                x <- listm(x)
        if(!is.na(grp[1])) {
                yy <- x
                x<-list()
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop("Data must be stored in list mode or a matrix.")
        for(j in 1:JK) {
                xx <- x[[j]]
                x[[j]] <- xx[!is.na(xx)] # Remove missing values
        }
        #

        if(JK != length(x))
                warning("The number of groups does not match the number of contrast coefficients.")
for(j in 1:JK){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
        # Create the three contrast matrices
temp<-con2way(J,K)
conA<-temp$conA
conB<-temp$conB
conAB<-temp$conAB
if(!op){
Factor.A<-lincon(x,con=conA,tr=tr,alpha=alpha,pr=pr)
Factor.B<-lincon(x,con=conB,tr=tr,alpha=alpha,pr=FALSE)
Factor.AB<-lincon(x,con=conAB,tr=tr,alpha=alpha,pr=FALSE)
}
All.Tests<-NA
if(op){
Factor.A<-NA
Factor.B<-NA
Factor.AB<-NA
con<-cbind(conA,conB,conAB)
All.Tests<-lincon(x,con=con,tr=tr,alpha=alpha)
}
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.AB=Factor.AB,All.Tests=All.Tests,conA=conA,conB=conB,conAB=conAB)
}


bbmcp=mcp2atm



################################################################################
# THREE-WAY FACTORIAL DESIGNS
################################################################################

mcp3atm<-function(J,K,L, x,tr=.2,con=0,alpha=.05,grp=NA,op=FALSE,pr=TRUE){
#
# Do all pairwise comparisons of
# main effects for Factor A and B and C and all interactions
# based on trimmed means
#
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
        JKL <- J*K*L
        if(is.matrix(x))
                x <- listm(x)
        if(!is.na(grp[1])) {
                yy <- x
                x<-list()
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop("Data must be stored in list mode or a matrix.")
        for(j in 1:JKL) {
                xx <- x[[j]]
                x[[j]] <- xx[!is.na(xx)] # Remove missing values
        }
        #

        if(JKL != length(x))
                warning("The number of groups does not match the number of contrast coefficients.")
for(j in 1:JKL){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
        # Create the three contrast matrices
temp<-con3way(J,K,L)
conA<-temp$conA
conB<-temp$conB
conC<-temp$conC
conAB<-temp$conAB
conAC<-temp$conAC
conBC<-temp$conBC
conABC<-temp$conABC
if(!op){
Factor.A<-lincon(x,con=conA,tr=tr,alpha=alpha,pr=pr)
Factor.B<-lincon(x,con=conB,tr=tr,alpha=alpha,pr=pr)
Factor.C<-lincon(x,con=conC,tr=tr,alpha=alpha,pr=pr)
Factor.AB<-lincon(x,con=conAB,tr=tr,alpha=alpha,pr=pr)
Factor.AC<-lincon(x,con=conAC,tr=tr,alpha=alpha,pr=pr)
Factor.BC<-lincon(x,con=conBC,tr=tr,alpha=alpha,pr=pr)
Factor.ABC<-lincon(x,con=conABC,tr=tr,alpha=alpha,pr=pr)
}
All.Tests<-NA
if(op){
Factor.A<-NA
Factor.B<-NA
Factor.C<-NA
Factor.AB<-NA
Factor.AC<-NA
Factor.BC<-NA
Factor.ABC<-NA
con<-cbind(conA,conB,conB,conAB,conAC,conBC,conABC)
All.Tests<-lincon(x,con=con,tr=tr,alpha=alpha,,pr=pr)
}
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.C=Factor.C,
Factor.AB=Factor.AB,Factor.AC=Factor.AC,Factor.BC=Factor.BC,
Factor.ABC=Factor.ABC,All.Tests=All.Tests,conA=conA,conB=conB,conC=conC,
conAB=conAB,conAC=conAC,conBC=conBC,conABC=conABC)
}


mcp3med<-function(J,K,L, x,tr=.2,con=0,alpha=.05,grp=NA,op=FALSE){
#
# Do all pairwise comparisons of
# main effects for Factor A and B and C and all interactions
# based on trimmed means
#
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
        JKL <- J*K*L
        if(is.matrix(x))
                x <- listm(x)
        if(!is.na(grp[1])) {
                yy <- x
                x<-list()
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop("Data must be stored in list mode or a matrix.")
        for(j in 1:JKL) {
                xx <- x[[j]]
                x[[j]] <- xx[!is.na(xx)] # Remove missing values
        }
        #

        if(JKL != length(x))
                warning("The number of groups does not match the number of contrast coefficients.")
for(j in 1:JKL){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
        # Create the three contrast matrices
temp<-con3way(J,K,L)
conA<-temp$conA
conB<-temp$conB
conC<-temp$conC
conAB<-temp$conAB
conAC<-temp$conAC
conBC<-temp$conBC
conABC<-temp$conABC
if(!op){
Factor.A<-msmed(x,con=conA,alpha=alpha)
Factor.B<-msmed(x,con=conB,alpha=alpha)
Factor.C<-msmed(x,con=conC,alpha=alpha)
Factor.AB<-msmed(x,con=conAB,alpha=alpha)
Factor.AC<-msmed(x,con=conAC,alpha=alpha)
Factor.BC<-msmed(x,con=conBC,alpha=alpha)
Factor.ABC<-msmed(x,con=conABC,alpha=alpha)
}
All.Tests<-NA
if(op){
Factor.A<-NA
Factor.B<-NA
Factor.C<-NA
Factor.AB<-NA
Factor.AC<-NA
Factor.BC<-NA
Factor.ABC<-NA
con<-cbind(conA,conB,conB,conAB,conAC,conBC,conABC)
All.Tests<-msmed(x,con=con,alpha=alpha)
}
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.C=Factor.C,
Factor.AB=Factor.AB,Factor.AC=Factor.AC,Factor.BC=Factor.BC,
Factor.ABC=Factor.ABC,All.Tests=All.Tests,conA=conA,conB=conB,conC=conC,
conAB=conAB,conAC=conAC,conBC=conBC,conABC=conABC)
}


rm3mcp<-function(J,K,L, x,tr=.2,alpha=.05,dif=TRUE,op=FALSE,grp=NA){
#
# MULTIPLE COMPARISONS FOR A 3-WAY within-by-within-by within ANOVA
# Do all multiple comparisons associated with
# main effects for Factor A and B and C and all interactions
# based on trimmed means
#
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
        JKL <- J*K*L
        if(is.matrix(x))
                x <- listm(x)
        if(!is.na(grp[1])) {
                yy <- x
                x<-list()
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop("Data must be stored in list mode or a matrix.")
        for(j in 1:JKL) {
                xx <- x[[j]]
                x[[j]] <- xx[!is.na(xx)] # Remove missing values
        }
        #

        if(JKL != length(x))
                warning("The number of groups does not match the number of contrast coefficients.")
for(j in 1:JKL){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
        # Create the three contrast matrices
temp<-con3way(J,K,L)
conA<-temp$conA
conB<-temp$conB
conC<-temp$conC
conAB<-temp$conAB
conAC<-temp$conAC
conBC<-temp$conBC
conABC<-temp$conABC
Factor.A<-rmmcp(x,con=conA,tr=tr,alpha=alpha,dif=dif)
Factor.B<-rmmcp(x,con=conB,tr=tr,alpha=alpha,dif=dif)
Factor.C<-rmmcp(x,con=conC,tr=tr,alpha=alpha,dif=dif)
Factor.AB<-rmmcp(x,con=conAB,tr=tr,alpha=alpha,dif=dif)
Factor.AC<-rmmcp(x,con=conAC,tr=tr,alpha=alpha,dif=dif)
Factor.BC<-rmmcp(x,con=conBC,tr=tr,alpha=alpha,dif=dif)
Factor.ABC<-rmmcp(x,con=conABC,tr=tr,alpha=alpha,dif=dif)
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.C=Factor.C,
Factor.AB=Factor.AB,Factor.AC=Factor.AC,Factor.BC=Factor.BC,
Factor.ABC=Factor.ABC,conA=conA,conB=conB,conC=conC,
conAB=conAB,conAC=conAC,conBC=conBC,conABC=conABC)
}



################################################################################
# REPEATED MEASURES MCP
################################################################################

rmmcp<-function(x, y=NULL,con = 0, tr = 0.2, alpha = 0.05,dif=TRUE,hoch=TRUE,na.rm=TRUE){
#
# MCP on trimmed means with FWE controlled with Hochberg's method
#  hoch=FALSE, will use Rom's method if alpha=.05 or .01 and number of tests is <=10
#
# Note: confidence intervals are adjusted based on the corresponding critical p-value.
#
if(!is.null(y))x=cbind(x,y)
flagcon=FALSE
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-ncol(x)
xbar<-vector("numeric",J)
x<-elimna(x)  # Remove missing values
nval<-nrow(x)
h1<-nrow(x)-2*floor(tr*nrow(x))
df<-h1-1
for(j in 1: J)xbar[j]<-mean(x[,j],tr)
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0)CC<-(J^2-J)/2
ncon<-CC
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
if(hoch)dvec<-alpha/c(1:ncon)
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
if(sum(con^2)==0){
flagcon<-TRUE
psihat<-matrix(0,CC,5)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c("Group","Group","test","p.value","p.crit","se"))
temp1<-0
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
q1<-(nrow(x)-1)*winvar(x[,j],tr)
q2<-(nrow(x)-1)*winvar(x[,k],tr)
q3<-(nrow(x)-1)*wincor(x[,j],x[,k],tr)$cov
sejk<-sqrt((q1+q2-2*q3)/(h1*(h1-1)))
if(!dif){
test[jcom,6]<-sejk
test[jcom,3]<-(xbar[j]-xbar[k])/sejk
temp1[jcom]<-2 * (1 - pt(abs(test[jcom,3]), df))
test[jcom,4]<-temp1[jcom]
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-(xbar[j]-xbar[k])
}
if(dif){
dv<-x[,j]-x[,k]
test[jcom,6]<-trimse(dv,tr)
temp<-trimci(dv,alpha=alpha/CC,pr=FALSE,tr=tr)
test[jcom,3]<-temp$test.stat
temp1[jcom]<-temp$p.value
test[jcom,4]<-temp1[jcom]
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-mean(dv,tr=tr)
psihat[jcom,4]<-temp$ci[1]
psihat[jcom,5]<-temp$ci[2]
}
}}}
if(hoch)dvec<-alpha/c(1:ncon)
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2,4]>=zvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
#zvec[ddd:ncon]<-dvec[ddd]
# Redo confidence intervals so that they are consistent with p-values:
}
test[temp2,5]<-zvec
if(!dif){
psihat[,4]<-psihat[,3]-qt(1-test[,5]/2,df)*test[,6]
psihat[,5]<-psihat[,3]+qt(1-test[,5]/2,df)*test[,6]
}}
if(sum(con^2)>0){
if(nrow(con)!=ncol(x))warning("The number of groups does not match the number
 of contrast coefficients.")
ncon<-ncol(con)
psihat<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c("con.num","test","p.value","p.crit","se"))
temp1<-NA
for (d in 1:ncol(con)){
psihat[d,1]<-d
if(!dif){
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-0
for(j in 1:J){
for(k in 1:J){
djk<-(nval-1)*wincor(x[,j],x[,k], tr)$cov/(h1*(h1-1))
sejk<-sejk+con[j,d]*con[k,d]*djk
}}
sejk<-sqrt(sejk)
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
test[d,5]<-sejk
temp1[d]<-2 * (1 - pt(abs(test[d,2]), df))
}
if(dif){
for(j in 1:J){
if(j==1)dval<-con[j,d]*x[,j]
if(j>1)dval<-dval+con[j,d]*x[,j]
}
temp1[d]<-trimci(dval,tr=tr,pr=FALSE)$p.value
test[d,1]<-d
test[d,2]<-trimci(dval,tr=tr,pr=FALSE)$test.stat
test[d,5]<-trimse(dval,tr=tr)
psihat[d,2]<-mean(dval,tr=tr)
}}
test[,3]<-temp1
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2,3]>=zvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
}
test[temp2,4]<-zvec
psihat[,3]<-psihat[,2]-qt(1-test[,4]/2,df)*test[,5]
psihat[,4]<-psihat[,2]+qt(1-test[,4]/2,df)*test[,5]
}
num.sig=nrow(test)
if(flagcon){
ior=order(test[,5],decreasing=TRUE)
for(j in 1:nrow(test)){
if(test[ior[j],4]<=test[ior[j],5])break
else num.sig=num.sig-1
}
}
if(!flagcon){
ior=order(test[,4],decreasing=TRUE)
for(j in 1:nrow(test)){
if(test[ior[j],3]<=test[ior[j],4])break
else num.sig=num.sig-1
}
}
list(n=nval,test=test,psihat=psihat,con=con,num.sig=num.sig)
}


wmcp<-rmmcp


rmmcppbd<-function(x,y=NULL,alpha=.05,con=0,est=onestep,plotit=TRUE,grp=NA,nboot=NA,
hoch=TRUE,SEED=TRUE,...){
#
#   Use a percentile bootstrap method to  compare dependent groups
#   based on difference scores.
#   By default,
#   compute a .95 confidence interval for all linear contrasts
#   specified by con, a J by C matrix, where  C is the number of
#   contrasts to be tested, and the columns of con are the
#   contrast coefficients.
#   If con is not specified, all pairwise comparisons are done.
#
#   By default, one-step M-estimator is used
#    and a sequentially rejective method
#   is used to control the probability of at least one Type I error.
#
#   nboot is the bootstrap sample size. If not specified, a value will
#   be chosen depending on the number of contrasts there are.
#
#   x can be an n by J matrix or it can have list mode
#   for two groups, data for second group can be put in y
#   otherwise, assume x is a matrix (n by J) or has list mode.
#
#   A sequentially rejective method is used to control alpha.
#   If n>=80, hochberg's method is used.
#
if(!is.null(y[1]))x<-cbind(x,y)
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
if(is.matrix(con)){
if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
}}
if(is.list(x)){
# put the data in an n by J matrix
mat<-matl(x)
}
if(is.matrix(x) && is.matrix(con)){
if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
mat<-x
}
if(is.matrix(x))mat<-x
if(!is.na(sum(grp)))mat<-mat[,grp]
x<-mat
mat<-elimna(mat) # Remove rows with missing values.
x<-mat
J<-ncol(mat)
n=nrow(mat)
if(n>=80)hoch=TRUE
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
nboot<-5000
if(d<=10)nboot<-3000
if(d<=6)nboot<-2000
if(d<=4)nboot<-1000
}
n<-nrow(mat)
crit.vec<-alpha/c(1:d)
connum<-ncol(con)
# Create set of differences based on contrast coefficients
xx<-x%*%con
xx<-as.matrix(xx)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
psihat<-matrix(0,connum,nboot)
bvec<-matrix(NA,ncol=connum,nrow=nboot)
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
# data is an nboot by n matrix
if(ncol(xx)==1){
for(ib in 1:nboot)psihat[1,ib]<-est(xx[data[ib,]],...)
}
if(ncol(xx)>1){
for(ib in 1:nboot)psihat[,ib]<-apply(elimna(xx[data[ib,],]),2,est,...)
}
#
# Now have an nboot by connum matrix of bootstrap values.
#
test<-1
icl<-round(alpha*nboot/2)+1
icu<-nboot-icl-1
cimat=matrix(NA,nrow=connum,ncol=2)
for (ic in 1:connum){
test[ic]<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
test[ic]<-min(test[ic],1-test[ic])
temp=sort(psihat[ic,])
cimat[ic,1]=temp[icl]
cimat[ic,2]=temp[icu]
}
test<-2*test
ncon<-ncol(con)
if(alpha==.05){
dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
dvec[2]<-alpha/2
}
if(hoch)dvec<-alpha/(2*c(1:ncon))
dvec<-2*dvec
if(plotit && connum==1){
plot(c(psihat[1,],0),xlab="",ylab="Est. Difference")
points(psihat[1,])
abline(0,0)
}
temp2<-order(0-test)
ncon<-ncol(con)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output<-matrix(0,connum,6)
dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.crit","ci.lower","ci.upper"))
tmeans<-apply(xx,2,est,...)
psi<-1
output[temp2,4]<-zvec
for (ic in 1:ncol(con)){
output[ic,2]<-tmeans[ic]
output[ic,1]<-ic
output[ic,3]<-test[ic]
output[ic,5:6]<-cimat[ic,]
}
num.sig<-sum(output[,3]<=output[,4])
list(output=output,con=con,num.sig=num.sig)
}


rmmcppbtm<-function(x,alpha=.05,con=0,tr=.2,grp=NA,nboot=NA){
#
#   Using the percentile bootstrap method,
#   compute a .95 confidence interval for all linear contasts
#   specified by con, a J by C matrix, where  C is the number of
#   contrasts to be tested, and the columns of con are the
#   contrast coefficients.
#
#   The trimmed means of dependent groups are being compared.
#   By default, 20% trimming is used.
#
#   nboot is the bootstrap sample size. If not specified, a value will
#   be chosen depending on the number of contrasts there are.
#
#   x can be an n by J matrix or it can have list mode
#
#   For alpha=.05, some critical values have been
#   determined via simulations and are used by this function;
#   otherwise an approximation is used.
#
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
if(is.matrix(con)){
if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
}}
if(is.list(x)){
# put the data in an n by J matrix
mat<-matrix(0,length(x[[1]]),length(x))
for (j in 1:length(x))mat[,j]<-x[[j]]
}
if(is.matrix(x) && is.matrix(con)){
if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
mat<-x
}
if(is.matrix(x))mat<-x
if(!is.na(sum(grp)))mat<-mat[,grp]
mat<-elimna(mat) # Remove rows with missing values.
J<-ncol(mat)
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
if(is.na(crit) && tr != .2){
print("A critical value must be specified when")
stop("the amount of trimming differs from .2")
}
if(is.na(nboot)){
if(d<=3)nboot<-1000
if(d==6)nboot<-2000
if(d==10)nboot<-4000
if(d==15)nboot<-8000
if(d==21)nboot<-8000
if(d==28)nboot<-10000
}
n<-nrow(mat)
crit<-NA
if(alpha==.05){
if(d==1)crit<-alpha/2
if(d==3){
crit<-.004
if(n>=15)crit<-.006
if(n>=30)crit<-.007
if(n>=40)crit<-.008
if(n>=100)crit<-.009
}
if(d==6){
crit<-.001
if(n>=15)crit<-.002
if(n>=20)crit<-.0025
if(n>=30)crit<-.0035
if(n>=40)crit<-.004
if(n>=60)crit<-.0045
}
if(d==10){
crit<-.00025
if(n>=15)crit<-.00125
if(n>=20)crit<-.0025
}
if(d==15){
crit<-.0005
if(n>=20)crit<-.0010
if(n>=30)crit<-.0011
if(n>=40)crit<-.0016
if(n>=100)crit<-.0019
}
if(d==21){
crit<-.00025
if(n>=20)crit<-.00037
if(n>=30)crit<-.00075
if(n>=40)crit<-.00087
if(n>=60)crit<-.00115
if(n>=100)crit<-.00125
}
if(d==28){
crit<-.0004
if(n>=30)crit<-.0006
if(n>=60)crit<-.0008
if(n>=100)crit<-.001
}
}
if(is.na(crit)){
crit<-alpha/(2*d)
if(n<20)crit<-crit/2
if(n<=10)crit<-crit/2
}
icl<-ceiling(crit*nboot)+1
icu<-ceiling((1-crit)*nboot)
connum<-ncol(con)
set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
# data is an nboot by n matrix
xbars<-matrix(0,nboot,ncol(mat))
psihat<-matrix(0,connum,nboot)
print("Taking bootstrap samples. Please wait.")
bvec<-bootdep(mat,tr,nboot)
#
# Now have an nboot by J matrix of bootstrap values.
#
test<-1
for (ic in 1:connum){
psihat[ic,]<-apply(bvec,1,bptdpsi,con[,ic])
test[ic]<-sum((psihat[ic,]>0))/nboot
test[ic]<-min(test[ic],1-test[ic])
}
print("Reminder: Test statistic must be less than critical value in order to reject.")
output<-matrix(0,connum,5)
dimnames(output)<-list(NULL,c("con.num","psihat","test","ci.lower","ci.upper"))
tmeans<-apply(mat,2,mean,trim=tr)
psi<-1
for (ic in 1:ncol(con)){
output[ic,2]<-sum(con[,ic]*tmeans)
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(psihat[ic,])
output[ic,4]<-temp[icl]
output[ic,5]<-temp[icu]
}
list(output=output,crit=crit,con=con)
}


rmmcppbv2<-function(x,y=NULL,alpha=.05,NA.RM=TRUE,con=NULL,method='hoch',
est=tmean,plotit=FALSE,dif=TRUE,grp=NA,nboot=NA,BA=FALSE,xlab="Group 1",ylab="Group 2",pr=TRUE,SEED=TRUE,SR=FALSE,...){
#
#   Use a percentile bootstrap method to  compare dependent groups.
#
# Unlike rmmcppb, this function defaults to comparing the marginal trimmed means
#  that uses all of the data that is not missing.
# For example, if col 1 has missing values but the other columns have no missing values,
# all of the data in the other columns are used.
# Using mmcppb, any row with missing values is removed.
#
#   By default,
#   compute a .95 confidence interval for all linear contrasts
#   specified by con, a J-by-C matrix, where  C is the number of
#   contrasts to be tested, and the columns of con are the
#   contrast coefficients.
#   If con is not specified, all pairwise comparisons are done.
#
#   FWE is controlled using  Hochberg's method.
#. Setting method ='fdr' would control the false discovery rate
#
#   dif=TRUE indicates that difference scores are to be used
#   dif=FALSE indicates that measure of location associated with
#   marginal distributions are used instead.
#
#   nboot is the bootstrap sample size. If not specified, a value will
#   be chosen depending on the number of contrasts there are.
#
#   x can be an n by J matrix or it can have list mode
#   for two groups, data for second group can be put in y
#   otherwise, assume x is a matrix (n by J) or has list mode.
#
#   Hochberg's method is used to control FWE
#
if(is.list(x))mat<-matl(x)  # put the data in an n by J matrix
J=ncol(x)
if(is.null(con))con=con.all.pairs(J)
if(dif){
if(pr){print("dif=TRUE, so analysis is done on difference scores.")
print(" Each confidence interval has probability coverage 1-alpha.")
print("Hochberg's method used to control the FWE rate")
}
temp<-rmmcppbd(x,y=y,alpha=alpha,con=con,est,plotit=plotit,grp=grp,nboot=nboot, SEED=SEED,
hoch=TRUE,...)
output<-temp$output
con<-temp$con
}
if(!dif){
if(pr){
print("dif=FALSE, so analysis is done on marginal distributions")
}}
if(!is.null(y[1]))x<-cbind(x,y)
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
if(is.matrix(con)){
if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
}}

if(is.matrix(x) && is.matrix(con)){
if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
mat<-x
}
if(is.matrix(x))mat<-x
if(!is.na(sum(grp)))mat<-mat[,grp]
x<-mat
J<-ncol(mat)
xcen<-x
for(j in 1:J)xcen[,j]<-x[,j]-est(x[,j],na.rm=NA.RM,...)
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
connum<-ncol(con)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
xbars<-apply(mat,2,est,na.rm=NA.RM,...)
psidat<-NA
for (ic in 1:connum)psidat[ic]<-sum(con[,ic]*xbars)
psihat<-matrix(0,connum,nboot)
psihatcen<-matrix(0,connum,nboot)
bvec<-matrix(NA,ncol=J,nrow=nboot)
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(ib in 1:nboot){
bvec[ib,]<-apply(x[data[ib,],],2,est,na.rm=NA.RM,...)
}
#
# Now have an nboot by J matrix of bootstrap values.
#
test<-1
for (ic in 1:connum){
psihat[ic,]<-apply(bvec,1,bptdpsi,con[,ic])
ptemp<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
test[ic]<-ptemp
test[ic]<-2*min(test[ic],1-test[ic])
}
ncon<-ncol(con)
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
ncon<-ncol(con)
output<-matrix(0,connum,6)
dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.adj","ci.lower","ci.upper"))
tmeans<-apply(mat,2,est,na.rm=NA.RM,...)
psi<-1
for (ic in 1:ncol(con)){
output[ic,2]<-sum(con[,ic]*tmeans)
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(psihat[ic,])
icl<-round(alpha*nboot/2)
icu<-nboot-icl
icl=icl+1
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
ids=NA
output[,4]=p.adjust(output[,3],method=method)
num.sig=sum(output[,4]<=alpha)
if(is.na(output[1,3])){
if(pr)print('Evidently, one or more groups have too many missing values')
}
list(output=output,con=con,num.sig=num.sig)
}


rmmcpv2<-function(x, y=NULL,con = 0, tr = 0.2, alpha = 0.05,dif=TRUE,
hoch=TRUE,na.rm=TRUE,nmin=5){
#
# MCP on trimmed means with FWE controlled with Hochberg's method
#  hoch=FALSE, will use Rom's method if alpha=.05 or .01 and number of tests is <=10
#
# Note: confidence intervals are adjusted based on the corresponding critical p-value.
#
if(!is.null(y))x=cbind(x,y)
flagcon=FALSE
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-ncol(x)
nval<-nrow(x)
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0)CC<-(J^2-J)/2
ncon<-CC
#if(alpha==.05){
#dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
#if(ncon > 10){
#avec<-.05/c(11:ncon)
#dvec<-c(dvec,avec)
#}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(hoch)dvec<-alpha/c(1:ncon)
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
if(sum(con^2)==0){
flagcon<-TRUE
psihat<-matrix(0,CC,5)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c("Group","Group","test","p.value","p.adj","se"))
temp1<-0
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
y=elimna(x[,c(j,k)])
if(is.null(dim(y)))y=matrix(c(1,1),nrow=1)
if(nrow(y)<=nmin)print(paste('Skipping group', j, ' and group', k, 'due to small sample size'))
if(nrow(y)>nmin){
h1<-nrow(y)-2*floor(tr*nrow(y))
df<-h1-1
xbar=mean(y[,1],tr=tr)
xbar[2]=mean(y[,2],tr=tr)
q1<-(nrow(y)-1)*winvar(y[,1],tr)
q2<-(nrow(y)-1)*winvar(y[,2],tr)
q3<-(nrow(y)-1)*wincor(y[,1],y[,2],tr)$cov
sejk<-sqrt((q1+q2-2*q3)/(h1*(h1-1)))
if(!dif){
test[jcom,6]<-sejk
test[jcom,3]<-(xbar[1]-xbar[2])/sejk
temp1[jcom]<-2 * (1 - pt(abs(test[jcom,3]), df))
test[jcom,4]<-temp1[jcom]
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-(xbar[1]-xbar[2])
}
if(dif){
dv<-y[,1]-y[,2]
test[jcom,6]<-trimse(dv,tr)
temp<-trimci(dv,alpha=alpha/CC,pr=FALSE,tr=tr)
test[jcom,3]<-temp$test.stat
temp1[jcom]<-temp$p.value
test[jcom,4]<-temp1[jcom]
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-mean(dv,tr=tr)
psihat[jcom,4]<-temp$ci[1]
psihat[jcom,5]<-temp$ci[2]
}
}}}}
if(hoch)dvec<-alpha/c(1:ncon)
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2,4]>=zvec)
sigvec=elimna(sigvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
}
test[temp2,5]<-zvec
if(!dif){
psihat[,4]<-psihat[,3]-qt(1-test[,5]/2,df)*test[,6]
psihat[,5]<-psihat[,3]+qt(1-test[,5]/2,df)*test[,6]
}}
if(sum(con^2)>0){
xbar=apply(x,2,mean,tr=tr)
if(nrow(con)!=ncol(x))warning("The number of groups does not match the number
 of contrast coefficients.")
ncon<-ncol(con)
psihat<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c("con.num","test","p.value","p.adj","se"))
temp1<-NA
for (d in 1:ncol(con)){
psihat[d,1]<-d
if(!dif){
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-0
for(j in 1:J){
for(k in 1:J){
djk<-(nval-1)*wincor(x[,j],x[,k], tr)$cov/(h1*(h1-1))
sejk<-sejk+con[j,d]*con[k,d]*djk
}}
sejk<-sqrt(sejk)
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
test[d,5]<-sejk
temp1[d]<-2 * (1 - pt(abs(test[d,2]), df))
}
if(dif){
for(j in 1:J){
if(j==1)dval<-con[j,d]*x[,j]
if(j>1)dval<-dval+con[j,d]*x[,j]
}
temp1[d]<-trimci(dval,tr=tr,pr=FALSE)$p.value
test[d,1]<-d
test[d,2]<-trimci(dval,tr=tr,pr=FALSE)$test.stat
test[d,5]<-trimse(dval,tr=tr)
psihat[d,2]<-mean(dval,tr=tr)
}}
test[,3]<-temp1
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2,3]>=zvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
}
test[temp2,4]<-zvec
psihat[,3]<-psihat[,2]-qt(1-test[,4]/2,df)*test[,5]
psihat[,4]<-psihat[,2]+qt(1-test[,4]/2,df)*test[,5]
}
num.sig=nrow(test)
#if(flagcon){
#ior=order(test[,5],decreasing=TRUE)
#for(j in 1:nrow(test)){
#if(test[ior[j],4]<=test[ior[j],5])break
#else num.sig=num.sig-1
#}
#}
#
test[,5]=p.adjust(test[,4],method='hoch')
num.sig=sum(elimna(test[,5])<=alpha)
list(n=nval,test=test,psihat=psihat,con=con,num.sig=num.sig)
}


rmmcpES<-function(x, con = 0, tr = 0.2, alpha = 0.05,dif=TRUE,hoch=TRUE,pr=TRUE){
#
#  Like rmmcp,only a robust version of Cohen's d is included.
#  Designed only for all pairwise comparisons.
#
if(con!=0)stop('This function is for all pairwise comparisons only')
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
a=rmmcp(x,con=con,tr=tr,alpha=alpha,dif=dif,hoch=hoch)
test=a$test
J=ncol(x)
CC=(J^2-J)/2
psihat<-matrix(0,CC,6)
dimnames(psihat)<-list(NULL,c('Group','Group','psihat','ci.lower','ci.upper','Effect.Size'))
psihat[,1:5]=a$psihat
ic=0
for (j in 1:J){
for (k in 1:J){
if (j < k){
ic=ic+1
if(dif)psihat[ic,6]=D.akp.effect(x[,j],x[,k],tr=tr)
if(!dif){
psihat[ic,6]=yuendv2(x[,j],x[,k],tr=tr)$Effect.Size
if(pr)print('Note: With dif=FALSE, explanatory measure of effect size is used')
}
}}}
list(test=test,psihat=psihat)
}


rmmcpQS<-function(x, y=NULL,con = 0, tr = 0.2, alpha = 0.05,dif=TRUE,hoch=TRUE,locfun=tmean,...){
#
#  Same as rmmcp, only includes quantile shift measure of effect size based on difference scores.
# MCP on trimmed means with FWE controlled with Hochberg's method
#  hoch=FALSE, will use Rom's method if alpha=.05 or .01 and number of tests is <=10
#
# Note: confidence intervals are adjusted based on the corresponding critical p-value.
#
if(!is.null(y))x=cbind(x,y)
flagcon=FALSE
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
con<-as.matrix(con)
J<-ncol(x)
xbar<-vector('numeric',J)
x<-elimna(x)  # Remove missing values
nval<-nrow(x)
h1<-nrow(x)-2*floor(tr*nrow(x))
df<-h1-1
for(j in 1: J)xbar[j]<-mean(x[,j],tr)
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0)CC<-(J^2-J)/2
ncon<-CC
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
if(hoch)dvec<-alpha/c(1:ncon)
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
if(sum(con^2)==0){
flagcon<-TRUE
psihat<-matrix(0,CC,6)
dimnames(psihat)<-list(NULL,c('Group','Group','psihat','ci.lower','ci.upper','Q.effect'))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c('Group','Group','test','p.value','p.crit','se'))
temp1<-0
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
q1<-(nrow(x)-1)*winvar(x[,j],tr)
q2<-(nrow(x)-1)*winvar(x[,k],tr)
q3<-(nrow(x)-1)*wincor(x[,j],x[,k],tr)$cov
sejk<-sqrt((q1+q2-2*q3)/(h1*(h1-1)))
if(!dif){
test[jcom,6]<-sejk
test[jcom,3]<-(xbar[j]-xbar[k])/sejk
temp1[jcom]<-2 * (1 - pt(abs(test[jcom,3]), df))
test[jcom,4]<-temp1[jcom]
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-(xbar[j]-xbar[k])
}
if(dif){
dv<-x[,j]-x[,k]
test[jcom,6]<-trimse(dv,tr)
temp<-trimci(dv,alpha=alpha/CC,pr=FALSE,tr=tr)
test[jcom,3]<-temp$test.stat
temp1[jcom]<-temp$p.value
test[jcom,4]<-temp1[jcom]
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-mean(dv,tr=tr)
psihat[jcom,4]<-temp$ci[1]
psihat[jcom,5]<-temp$ci[2]
psihat[jcom,6]=depQS(x[,j],x[,k],locfun=locfun,...)$Q.effect
}
}}}
if(hoch)dvec<-alpha/c(1:ncon)
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2,4]>=zvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
}
test[temp2,5]<-zvec
if(!dif){
psihat[,4]<-psihat[,3]-qt(1-test[,5]/2,df)*test[,6]
psihat[,5]<-psihat[,3]+qt(1-test[,5]/2,df)*test[,6]
psihat[,6]=depQS(x[,j],x[,k],locfun=locfun,...)$Q.effect
}}

if(sum(con^2)>0){
if(nrow(con)!=ncol(x))warning('The number of groups does not match the number
 of contrast coefficients.')
ncon<-ncol(con)
psihat<-matrix(0,ncol(con),5)
dimnames(psihat)<-list(NULL,c('con.num','psihat','ci.lower','ci.upper','Q.effect'))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c('con.num','test','p.value','p.crit','se'))
temp1<-NA
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,5]=lindQS(x,con[,d],locfun=locfun,...)$Q.effect
if(!dif){
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-0
for(j in 1:J){
for(k in 1:J){
djk<-(nval-1)*wincor(x[,j],x[,k], tr)$cov/(h1*(h1-1))
sejk<-sejk+con[j,d]*con[k,d]*djk
}}
sejk<-sqrt(sejk)
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
test[d,5]<-sejk
temp1[d]<-2 * (1 - pt(abs(test[d,2]), df))
}
if(dif){
for(j in 1:J){
if(j==1)dval<-con[j,d]*x[,j]
if(j>1)dval<-dval+con[j,d]*x[,j]
}
temp1[d]<-trimci(dval,tr=tr,pr=FALSE)$p.value
test[d,1]<-d
test[d,2]<-trimci(dval,tr=tr,pr=FALSE)$test.stat
test[d,5]<-trimse(dval,tr=tr)
psihat[d,2]<-mean(dval,tr=tr)
}}
test[,3]<-temp1
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2,3]>=zvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
}
test[temp2,4]<-zvec
psihat[,3]<-psihat[,2]-qt(1-test[,4]/2,df)*test[,5]
psihat[,4]<-psihat[,2]+qt(1-test[,4]/2,df)*test[,5]
}
num.sig=nrow(test)
if(flagcon){
ior=order(test[,5],decreasing=TRUE)
for(j in 1:nrow(test)){
if(test[ior[j],4]<=test[ior[j],5])break
else num.sig=num.sig-1
}
}
if(!flagcon){
ior=order(test[,4],decreasing=TRUE)
for(j in 1:nrow(test)){
if(test[ior[j],3]<=test[ior[j],4])break
else num.sig=num.sig-1
}
}
list(n=nval,test=test,psihat=psihat,con=con,num.sig=num.sig)
}


rmmismcp<-function(x,y=NA,alpha=.05,con=0,est=tmean,plotit=TRUE,grp=NA,nboot=500,
SEED=TRUE,xlab="Group 1",ylab="Group 2",pr=FALSE,...){
#
#   Use a percentile bootstrap method to  compare  marginal measures of location for dependent groups.
#   Missing values are allowed; vectors of observations that contain
#   missing values are not simply removed as done by rmmcppb.
#   Only marginal measures of location are compared,
#   The function computes a .95 confidence interval for all linear contrasts
#   specified by con, a J by C matrix, where  C is the number of
#   contrasts to be tested, and the columns of con are the
#   contrast coefficients.
#   If con is not specified, all pairwise comparisons are done.
#
#   By default, a 20% trimmed is used and a sequentially rejective method
#   is used to control the probability of at least one Type I error.
#
#   nboot is the bootstrap sample size.
#
#   x can be an n by J matrix or it can have list mode
#   for two groups, data for second group can be put in y
#   otherwise, assume x is a matrix (n by J) or has list mode.
#
#
if(!is.na(y[1]))x<-cbind(x,y)
if(is.list(x))x=matl(x)
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
if(is.matrix(con)){
if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
}}
if(is.list(x)){
# put the data in an n by J matrix
mat<-matl(x)
}
if(is.matrix(x) && is.matrix(con)){
if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
mat<-x
}
J<-ncol(x)
Jm<-J-1
flag.con=F
if(sum(con^2)==0){
flag.con=T
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
n<-nrow(x)
crit.vec<-alpha/c(1:d)
connum<-ncol(con)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
xbars<-apply(x,2,est,na.rm=TRUE)
psidat<-NA
bveccen<-matrix(NA,ncol=J,nrow=nboot)
for (ic in 1:connum)psidat[ic]<-sum(con[,ic]*xbars)
psihat<-matrix(0,connum,nboot)
psihatcen<-matrix(0,connum,nboot)
bvec<-matrix(NA,ncol=J,nrow=nboot)
if(pr)print("Taking bootstrap samples. Please wait.")
data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
for(ib in 1:nboot){
bvec[ib,]<-apply(x[data[ib,],],2,est,na.rm=TRUE,...)
}
#
# Now have an nboot by J matrix of bootstrap measures of location.
#
test<-1
for (ic in 1:connum){
for(ib in 1:nboot){
psihat[ic,ib]=sum(con[,ic]*bvec[ib,])
}
matcon=c(0,psihat[ic,])
dis=mean((psihat[ic,]<0))+.5*mean((psihat[ic,]==0))
test[ic]<-2*min(c(dis,1-dis)) # the p-value
}
ncon<-ncol(con)
dvec<-alpha/c(1:ncon)
if(plotit && ncol(bvec)==2){
z<-c(0,0)
one<-c(1,1)
plot(rbind(bvec,z,one),xlab=xlab,ylab=ylab,type="n")
points(bvec)
totv<-apply(x,2,est,na.rm=TRUE,...)
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
sigvec<-(test[temp2]>=zvec)
output<-matrix(0,connum,6)
dimnames(output)<-list(NULL,c("con.num","psihat","p.value",
"crit.sig","ci.lower","ci.upper"))
tmeans<-apply(x,2,est,na.rm=TRUE,...)
psi<-1
output[temp2,4]<-zvec
for (ic in 1:ncol(con)){
output[ic,2]<-sum(con[,ic]*tmeans)
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(psihat[ic,])
icl<-round(output[ic,4]*nboot/2)+1
icu<-nboot-(icl-1)
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
if(!flag.con){
}
if(flag.con){
CC=(J^2-J)/2
test<-matrix(NA,CC,7)
dimnames(test)<-list(NULL,c("Group","Group","psi.hat","p.value","p.crit",
"ci.low","ci.upper"))
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
test[jcom,1]=j
test[jcom,2]=k
test[jcom,3:5]=output[jcom,2:4]
test[jcom,6:7]=output[jcom,5:6]
con=NULL
}}}}
if(!flag.con)test=output
#num.sig<-sum(output[,4]<=output[,5])
if(flag.con)num.sig<-sum(test[,4]<=test[,5])
if(!flag.con)num.sig<-sum(test[,3]<=test[,4])
list(output=test,con=con,num.sig=num.sig)
}



################################################################################
# BETWEEN-WITHIN DESIGNS
################################################################################

bwmcp<-function(J, K, x, tr = 0.2, JK = J * K, con = 0,
 alpha = 0.05, grp =c(1:JK), nboot = 599, method='hoch',SEED = TRUE, ...)
{
        #
        # A bootstrap-t for multiple comparisons among
        # for all main effects and interactions.
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
x=y
}

        conM = con2way(J, K)
 p <- J * K
        v <- matrix(0, p, p)
        data <- list()
xx=list()
        for(j in 1:length(x)) {
                data[[j]] <- x[[grp[j]]]
xx[[j]]=x[[grp[j]]] # save input data
                # Now have the groups in proper order.
                data[[j]] = data[[j]] - mean(data[[j]], tr = tr,na.rm=TRUE)  #centered data for bootstrapping
        }
ilow=1-K
iup=0
for(j in 1:J){
ilow <- ilow + K
 iup = iup + K
sel <- c(ilow:iup)
xx[sel]=listm(elimna(matl(xx[sel])))
 v[sel, sel] <- covmtrim(xx[sel], tr)
                }
A=lindep(xx,conM$conA,cmat=v,tr=tr)$test.stat
B=lindep(xx,conM$conB,cmat=v,tr=tr)$test.stat
AB=lindep(xx,conM$conAB,cmat=v,tr=tr)$test.stat
        x <- data
        jp <- 1 - K
        kv <- 0
        if(SEED)
                set.seed(2)
        # set seed of random number generator so that
        #             results can be duplicated.
        # Next determine the n_j values
        nvec <- NA
        testA = NA
        testB = NA
        testAB = NA
        bsam = list()
        bdat = list()
aboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conA))
bboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conB))
abboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conAB))
#        for(j in 1:J)
#                nvec[j] = length(x[[j]])
        for(ib in 1:nboot) {
                ilow <- 1 - K
                iup = 0
 for(j in 1:J) {
 ilow <- ilow + K
 iup = iup + K
nv=length(xx[[ilow]])
bdat[[j]] = sample(nv, size = nv, replace =TRUE)
for(k in ilow:iup){
# bsam[[k]] = xx[[k]][bdat[[j]]]
bsam[[k]] = data[[k]][bdat[[j]]]  # Use centered data to determine critical value.
}
 }
ilow=0-K
iup=0
for(j in 1:J){
ilow <- ilow + K
 iup = iup + K
sel <- c(ilow:iup)
 v[sel, sel] <- covmtrim(bsam[sel], tr)
                }

temp=abs(lindep(bsam,conM$conA, cmat=v,tr=tr)$test.stat[,4])
aboot[ib,]=temp
testA[ib] = max(temp)
temp=abs(lindep(bsam,conM$conB,cmat=v,tr=tr)$test.stat[,4])
bboot[ib,]=temp
testB[ib]= max(temp)
temp=abs(lindep(bsam,conM$conAB,cmat=v,tr=tr)$test.stat[,4])
testAB[ib] = max(temp)
abboot[ib,]=temp
        }
pbA=NA
pbB=NA
pbAB=NA
for(j in 1:ncol(aboot))pbA[j]=mean((abs(A[j,4])<aboot[,j]))
for(j in 1:ncol(bboot))pbB[j]=mean((abs(B[j,4])<bboot[,j]))
for(j in 1:ncol(abboot))pbAB[j]=mean((abs(AB[j,4])<abboot[,j]))
        critA = sort(testA)
        critB = sort(testB)
        critAB = sort(testAB)
        ic <- floor((1 - alpha) * nboot)
        critA = critA[ic]
        critB = critB[ic]
        critAB = critAB[ic]
cr=matrix(critA,ncol=1,nrow=nrow(A))
dimnames(cr)<-list(NULL,c('adj.p.value'))
A=cbind(A,cr)
pv=matrix(pbA,ncol=1,nrow=nrow(A))
dimnames(pv)<-list(NULL,c('p.value'))
A=cbind(A,pv)
cr=matrix(critB,ncol=1,nrow=nrow(B))
dimnames(cr)<-list(NULL,c('adj.p.value'))
B=cbind(B,cr)
pv=matrix(pbB,ncol=1,nrow=nrow(B))
dimnames(pv)<-list(NULL,c('p.value'))
B=cbind(B,pv)
cr=matrix(critAB,ncol=1,nrow=nrow(AB))
dimnames(cr)<-list(NULL,c('adj.p.value'))
AB=cbind(AB,cr)
pv=matrix(pbAB,ncol=1,nrow=nrow(AB))
dimnames(pv)<-list(NULL,c('p.value'))
AB=cbind(AB,pv)
# Update the output using adjusted p-values
A[,5]=p.adjust(A[,6],method=method)
B[,5]=p.adjust(B[,6],method=method)
AB[,5]=p.adjust(AB[,6],method=method)
A=A[,c(1,2,3,4,6,5)]
B=B[,c(1,2,3,4,6,5)]
AB=AB[,c(1,2,3,4,6,5)]
list(Fac.A=A,Fac.B=B,Fac.AB=AB,contrast.coef=conM)
}


bwwmcp<-function(J, K, L, x, tr = 0.2, JKL = J * K*L, con = 0,
 alpha = 0.05, grp =c(1:JKL), nboot = 599, SEED = TRUE, ...)
{
        #
        # A bootstrap-t for multiple comparisons
#       all main effects and interactions.
#       a between-by-within-within design.
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x)) y[[j]] <- x[, j]
                x <- y
}

        conM = con3way(J,K,L)
 p <- J*K*L
if(p>length(x))stop('JKL is less than the Number of groups')
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
ilow=1-KL
iup=0
for(j in 1:J){
ilow <- ilow + KL
 iup = iup + KL
sel <- c(ilow:iup)
xx[sel]=listm(elimna(matl(xx[sel])))
 v[sel, sel] <- covmtrim(xx[sel], tr)
                }
A=lindep(xx,conM$conA,cmat=v,tr=tr)$test.stat
B=lindep(xx,conM$conB,cmat=v,tr=tr)$test.stat
C=lindep(xx,conM$conC,cmat=v,tr=tr)$test.stat
AB=lindep(xx,conM$conAB,cmat=v,tr=tr)$test.stat
AC=lindep(xx,conM$conAC,cmat=v,tr=tr)$test.stat
BC=lindep(xx,conM$conBC,cmat=v,tr=tr)$test.stat
ABC=lindep(xx,conM$conABC,cmat=v,tr=tr)$test.stat
        x <- data
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
aboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conA))
bboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conB))
cboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conC))
abboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conAB))
acboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conAC))
bcboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conBC))
abcboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conABC))
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
ilow=1-KL
iup=0
for(j in 1:J){
ilow <- ilow + KL
 iup = iup + KL
sel <- c(ilow:iup)
 v[sel, sel] <- covmtrim(bsam[sel], tr)
                }
temp=abs(lindep(bsam,conM$conA, cmat=v,tr=tr)$test.stat[,4])
aboot[ib,]=temp
testA[ib] = max(temp)
temp=abs(lindep(bsam,conM$conB,cmat=v,tr=tr)$test.stat[,4])
bboot[ib,]=temp
testB[ib]= max(temp)

temp=abs(lindep(bsam,conM$conC,cmat=v,tr=tr)$test.stat[,4])
cboot[ib,]=temp
testC[ib]= max(temp)

temp=abs(lindep(bsam,conM$conAC,cmat=v,tr=tr)$test.stat[,4])
acboot[ib,]=temp
testAC[ib]= max(temp)

temp=abs(lindep(bsam,conM$conBC,cmat=v,tr=tr)$test.stat[,4])
bcboot[ib,]=temp
testBC[ib]= max(temp)

temp=abs(lindep(bsam,conM$conAB,cmat=v,tr=tr)$test.stat[,4])
testAB[ib] = max(temp)
abboot[ib,]=temp

temp=abs(lindep(bsam,conM$conABC,cmat=v,tr=tr)$test.stat[,4])
abcboot[ib,]=temp
testABC[ib]= max(temp)

        }
pbA=NA
pbB=NA
pbC=NA
pbAB=NA
pbAC=NA
pbBC=NA
pbABC=NA
for(j in 1:ncol(aboot))pbA[j]=mean((abs(A[j,4])<aboot[,j]))
for(j in 1:ncol(bboot))pbB[j]=mean((abs(B[j,4])<bboot[,j]))
for(j in 1:ncol(cboot))pbC[j]=mean((abs(C[j,4])<cboot[,j]))
for(j in 1:ncol(abboot))pbAB[j]=mean((abs(AB[j,4])<abboot[,j]))
for(j in 1:ncol(acboot))pbAC[j]=mean((abs(AC[j,4])<acboot[,j]))
for(j in 1:ncol(bcboot))pbBC[j]=mean((abs(BC[j,4])<bcboot[,j]))
for(j in 1:ncol(abcboot))pbABC[j]=mean((abs(ABC[j,4])<abcboot[,j]))
        critA = sort(testA)
        critB = sort(testB)
        critC = sort(testC)
        critAB = sort(testAB)
        critAC = sort(testAC)
        critBC = sort(testBC)
        critABC = sort(testABC)
        ic <- floor((1 - alpha) * nboot)
        critA = critA[ic]
        critB = critB[ic]
        critC = critC[ic]
        critAB = critAB[ic]
        critAC = critAC[ic]
        critBC = critBC[ic]
        critABC = critABC[ic]
critA=matrix(critA,ncol=1,nrow=nrow(A))
dimnames(critA)=list(NULL,c('crit.val'))
p.value=pbA
A=cbind(A,critA,p.value)


critB=matrix(critB,ncol=1,nrow=nrow(B))
dimnames(critB)=list(NULL,c('crit.val'))
p.value=pbB
B=cbind(B,critB,p.value)

critC=matrix(critC,ncol=1,nrow=nrow(C))
dimnames(critC)=list(NULL,c('crit.val'))
p.value=pbC
C=cbind(C,critC,p.value)

critAB=matrix(critAB,ncol=1,nrow=nrow(AB))
dimnames(critAB)=list(NULL,c('crit.val'))
p.value=pbAB
AB=cbind(AB,critAB,p.value)

critAC=matrix(critAC,ncol=1,nrow=nrow(AC))
dimnames(critAC)=list(NULL,c('crit.val'))
p.value=pbAC
AC=cbind(AC,critAC,p.value)


critBC=matrix(critBC,ncol=1,nrow=nrow(BC))
dimnames(critBC)=list(NULL,c('crit.val'))
p.value=pbBC
BC=cbind(BC,critBC,p.value)

critABC=matrix(critABC,ncol=1,nrow=nrow(ABC))
dimnames(critABC)=list(NULL,c('crit.val'))
p.value=pbABC
ABC=cbind(ABC,critABC,p.value)

list(Fac.A=A,Fac.B=B,Fac.C=C,Fac.AB=AB,Fac.AC=AC,Fac.BC=BC,Fac.ABC=ABC)
}


bbwmcp<-function(J, K, L, x, tr = 0.2, JKL = J * K*L, con = 0,
 alpha = 0.05, grp =c(1:JKL), nboot = 599, SEED = TRUE, ...)
{
        #
        # A bootstrap-t for multiple comparisons among
        # all main effects and interactions
#         for a between-by-between-within design.
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x)) y[[j]] <- x[, j]
                x <- y
}

        conM = con3way(J,K,L)
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
ilow=1-L
iup=0
for(j in 1:JK){
ilow <- ilow + L
 iup = iup + L
sel <- c(ilow:iup)
xx[sel]=listm(elimna(matl(xx[sel])))
 v[sel, sel] <- covmtrim(xx[sel], tr)
                }
A=lindep(xx,conM$conA,cmat=v,tr=tr)$test.stat
B=lindep(xx,conM$conB,cmat=v,tr=tr)$test.stat
C=lindep(xx,conM$conC,cmat=v,tr=tr)$test.stat
AB=lindep(xx,conM$conAB,cmat=v,tr=tr)$test.stat
AC=lindep(xx,conM$conAC,cmat=v,tr=tr)$test.stat
BC=lindep(xx,conM$conBC,cmat=v,tr=tr)$test.stat
ABC=lindep(xx,conM$conABC,cmat=v,tr=tr)$test.stat
        x <- data
        jp <- 1 - K
        kv <- 0
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
aboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conA))
bboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conB))
cboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conC))
abboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conAB))
acboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conAC))
bcboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conBC))
abcboot=matrix(NA,nrow=nboot,ncol=ncol(conM$conABC))
#        for(j in 1:JK)
#                nvec[j] = length(x[[j]])
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
ilow=1-L
iup=0
for(j in 1:JK){
ilow <- ilow + L
 iup = iup + L
sel <- c(ilow:iup)
 v[sel, sel] <- covmtrim(bsam[sel], tr)
                }
temp=abs(lindep(bsam,conM$conA, cmat=v,tr=tr)$test.stat[,4])
aboot[ib,]=temp
testA[ib] = max(temp)
temp=abs(lindep(bsam,conM$conB,cmat=v,tr=tr)$test.stat[,4])
bboot[ib,]=temp
testB[ib]= max(temp)

temp=abs(lindep(bsam,conM$conC,cmat=v,tr=tr)$test.stat[,4])
cboot[ib,]=temp
testC[ib]= max(temp)

temp=abs(lindep(bsam,conM$conAC,cmat=v,tr=tr)$test.stat[,4])
acboot[ib,]=temp
testAC[ib]= max(temp)

temp=abs(lindep(bsam,conM$conBC,cmat=v,tr=tr)$test.stat[,4])
bcboot[ib,]=temp
testBC[ib]= max(temp)

temp=abs(lindep(bsam,conM$conAB,cmat=v,tr=tr)$test.stat[,4])
testAB[ib] = max(temp)
abboot[ib,]=temp

temp=abs(lindep(bsam,conM$conABC,cmat=v,tr=tr)$test.stat[,4])
abcboot[ib,]=temp
testABC[ib]= max(temp)

        }
pbA=NA
pbB=NA
pbC=NA
pbAB=NA
pbAC=NA
pbBC=NA
pbABC=NA
for(j in 1:ncol(aboot))pbA[j]=mean((abs(A[j,4])<aboot[,j]))
for(j in 1:ncol(bboot))pbB[j]=mean((abs(B[j,4])<bboot[,j]))
for(j in 1:ncol(cboot))pbC[j]=mean((abs(C[j,4])<cboot[,j]))
for(j in 1:ncol(abboot))pbAB[j]=mean((abs(AB[j,4])<abboot[,j]))
for(j in 1:ncol(acboot))pbAC[j]=mean((abs(AC[j,4])<acboot[,j]))
for(j in 1:ncol(bcboot))pbBC[j]=mean((abs(BC[j,4])<bcboot[,j]))
for(j in 1:ncol(abcboot))pbABC[j]=mean((abs(ABC[j,4])<abcboot[,j]))
        critA = sort(testA)
        critB = sort(testB)
        critC = sort(testC)
        critAB = sort(testAB)
        critAC = sort(testAC)
        critBC = sort(testBC)
        critABC = sort(testABC)
        ic <- floor((1 - alpha) * nboot)
        critA = critA[ic]
        critB = critB[ic]
        critC = critC[ic]
        critAB = critAB[ic]
        critAC = critAC[ic]
        critBC = critBC[ic]
        critABC = critABC[ic]
critA=matrix(critA,ncol=1,nrow=nrow(A))
dimnames(critA)=list(NULL,c("crit.val"))
p.value=pbA
A=cbind(A,critA,p.value)

critB=matrix(critB,ncol=1,nrow=nrow(B))
dimnames(critB)=list(NULL,c("crit.val"))
p.value=pbB
B=cbind(B,critB,p.value)

critC=matrix(critC,ncol=1,nrow=nrow(C))
dimnames(critC)=list(NULL,c("crit.val"))
p.value=pbC
C=cbind(C,critC,p.value)

critAB=matrix(critAB,ncol=1,nrow=nrow(AB))
dimnames(critAB)=list(NULL,c("crit.val"))
p.value=pbAB
AB=cbind(AB,critAB,p.value)

critAC=matrix(critAC,ncol=1,nrow=nrow(AC))
dimnames(critAC)=list(NULL,c("crit.val"))
p.value=pbAC
AC=cbind(AC,critAC,p.value)


critBC=matrix(critBC,ncol=1,nrow=nrow(BC))
dimnames(critBC)=list(NULL,c("crit.val"))
p.value=pbBC
BC=cbind(BC,critBC,p.value)

critABC=matrix(critABC,ncol=1,nrow=nrow(ABC))
dimnames(critABC)=list(NULL,c("crit.val"))
p.value=pbABC
ABC=cbind(ABC,critABC,p.value)

list(Fac.A=A,Fac.B=B,Fac.C=C,Fac.AB=AB,Fac.AC=AC,Fac.BC=BC,Fac.ABC=ABC)
}


bwrmcp<-function(J,K,x,grp=NA,alpha=.05,bhop=TRUE){
#
# Do all pairwise comparisons of
# main effects for Factor A and B and all interactions
# using a rank-based method that tests for equal distributions.
#
#  A between by within subjects design is assumed.
#  Levels of Factor A are assumed to be independent and
#  levels of Factor B are dependent.
#
#  The data are assumed to be stored in x in list mode or in a matrix.
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
 if(is.list(x))xrem=matl(x)
        JK <- J * K
        if(is.matrix(x)){
                xrem=x
                x <- listm(x)
}

        if(!is.na(grp[1])) {
                yy <- x
                x<-list()
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop("Data must be stored in list mode or a matrix.")
#        for(j in 1:JK) {
#                xx <- x[[j]]
#                x[[j]] <- xx[!is.na(xx)] # Remove missing values
#        }
        #
if(JK != length(x))warning("The number of groups does not match the number of contrast coefficients.")
for(j in 1:JK){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
#
CC<-(J^2-J)/2
# Determine critical values
ncon<-CC*(K^2-K)/2
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
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
}
if(bhop)dvec<-(ncon-c(1:ncon)+1)*alpha/ncon
Fac.A<-matrix(0,CC,5)
dimnames(Fac.A)<-list(NULL,c("Level","Level","test.stat","p-value","sig.crit"))
mat<-matrix(c(1:JK),ncol=K,byrow=TRUE)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j < jj){
ic<-ic+1
Fac.A[ic,1]<-j
Fac.A[ic,2]<-jj
datsub=xrem[,c(mat[j,],mat[jj,])]
datsub=elimna(datsub)
#temp<-bwrank(2,K,elimna(x[,c(mat[j,],mat[jj,])]))
temp<-bwrank(2,K,datsub)
Fac.A[ic,3]<-temp$test.A
Fac.A[ic,4]<-temp$p.value.A
}}}
temp2<-order(0-Fac.A[,4])
Fac.A[temp2,5]<-dvec[1:length(temp2)]
CCB<-(K^2-K)/2
ic<-0
Fac.B<-matrix(0,CCB,5)
dimnames(Fac.B)<-list(NULL,c("Level","Level","test.stat","p-value","sig.crit"))
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
Fac.B[ic,1]<-k
Fac.B[ic,2]<-kk
mat1<-cbind(mat[,k],mat[,kk])
rv=c(mat1[,1],mat1[,2])
datsub=elimna(xrem[,sort(c(mat1[,1],mat1[,2]))])
temp<-bwrank(J,2,datsub)
Fac.B[ic,3]<-temp$test.B
Fac.B[ic,4]<-temp$p.value.B
}}}
temp2<-order(0-Fac.B[,4])
Fac.B[temp2,5]<-dvec[1:length(temp2)]
CCI<-CC*CCB
Fac.AB<-matrix(0,CCI,7)
dimnames(Fac.AB)<-list(NULL,c("Lev.A","Lev.A","Lev.B","Lev.B","test.stat","p-value","sig.crit"))
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j < jj){
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
Fac.AB[ic,1]<-j
Fac.AB[ic,2]<-jj
Fac.AB[ic,3]<-k
Fac.AB[ic,4]<-kk
val<-c(mat[j,k],mat[j,kk],mat[jj,k],mat[jj,kk])
val<-sort(val)
datsub=elimna(xrem[,val])
temp<-bwrank(2,2,datsub)
Fac.AB[ic,5]<-temp$test.AB
#Fac.AB[ic,6]<-temp$sig.AB
Fac.AB[ic,6]<-temp$p.value.AB
}}}}}}
temp2<-order(0-Fac.AB[,6])
Fac.AB[temp2,7]<-dvec[1:length(temp2)]
list(Factor.A=Fac.A,Factor.B=Fac.B,Factor.AB=Fac.AB)
}


bwamcp<-function(J,K,x,tr=.2,JK=J*K,grp=c(1:JK),alpha=.05,op=TRUE){
#
# All pairwise comparisons among levels of Factor A
# in a split-plot design using trimmed means.
#
# Data among dependent groups are pooled for each level
# of Factor A.
# Then this function calls lincon.
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
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}

JK<-J*K
if(!op){
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
data<-list()
jp<-1-K
kv<-0
for(j in 1:J){
jp<-jp+K
for(k in 1:K){
kv<-kv+1
if(k==1)temp<-x[[jp]]
if(k>1)temp<-c(temp,x[[kv]])
}
data[[j]]<-temp
}
print("Group numbers refer to levels of Factor A")
temp<-lincon(data,tr=tr,alpha=alpha)
}
if(op){
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
}}
print("Contrast Matrix Used:")
print(con)
temp<-lincon(x,con=con,tr=tr,alpha=alpha)
}
temp
}


bwimcp<-function(J,K,x,tr=.2,JK=J*K,grp=c(1:JK),alpha=.05,method='hoch'){
#
# Multiple comparisons for interactions
# in a split-plot design.
# The analysis is done by taking difference scores
# among all pairs of dependent groups and
# determining which of
# these differences differ across levels of Factor A
# using trimmed means.
#
#  FWE is controlled via Hochberg's method
# To adjusted p-values, use the function p.adjust
#
# For MOM or M-estimators, use spmcpi which uses a bootstrap method
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
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}

JK<-J*K
if(JK!=length(x))stop('Something is wrong. Expected ',JK,' groups but x contains ', length(x), 'groups instead.')
MJ<-(J^2-J)/2
MK<-(K^2-K)/2
JMK<-J*MK
MJMK<-MJ*MK
Jm<-J-1
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
output<-matrix(0,MJMK,7)
dimnames(output)<-list(NULL,c('A','A','B','B','psihat','p.value','p.crit'))
jp<-1-K
kv<-0
kv2<-0
test<-NA
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
m<-matrix(c(1:JK),J,K,byrow=T)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j<jj){
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
output[ic,1]<-j
output[ic,2]<-jj
output[ic,3]<-k
output[ic,4]<-kk
x1<-x[[m[j,k]]]-x[[m[j,kk]]]
x2<-x[[m[jj,k]]]-x[[m[jj,kk]]]
temp<-yuen(x1,x2)
output[ic,5]<-mean(x1,tr)-mean(x2,tr)
test[ic]<-temp$p.value
output[ic,6]<-test[ic]
}}}}}}
output[,7]<-p.adjust(output[,6],method=method)
sigvec=mean(output[,7] <= alpha)
output
}


bwimcpES<-function(J,K,x,tr=.2,JK=J*K,grp=c(1:JK),CI=FALSE, REL.M=NULL,alpha=.05,SEED=TRUE){
#
#  Same as bwimcp only several measures of effect size reported as well via the R function
#  bw.es.I
#
#
# Multiple comparisons for interactions
# in a split-plot design.
# The analysis is done by taking difference scores
# among all pairs of dependent groups and
# determining which of
# these differences differ across levels of Factor A
# using trimmed means.
#
#  FWE is controlled via Hochberg's method
# To adjusted p-values, use the function p.adjust
#
# For MOM or M-estimators, use spmcpi which uses a bootstrap method
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
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}

JK<-J*K
if(JK!=length(x))stop('Something is wrong. Expected ',JK,' groups but x contains ', length(x), 'groups instead.')
MJ<-(J^2-J)/2
MK<-(K^2-K)/2
JMK<-J*MK
MJMK<-MJ*MK
Jm<-J-1
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
output<-matrix(0,MJMK,7)
dimnames(output)<-list(NULL,c('A','A','B','B','psihat','p.value','p.crit'))
jp<-1-K
kv<-0
kv2<-0
test<-NA
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
m<-matrix(c(1:JK),J,K,byrow=T)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j<jj){
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
output[ic,1]<-j
output[ic,2]<-jj
output[ic,3]<-k
output[ic,4]<-kk
x1<-x[[m[j,k]]]-x[[m[j,kk]]]
x2<-x[[m[jj,k]]]-x[[m[jj,kk]]]
temp<-yuenv2(x1,x2,SEED=SEED)
output[ic,5]<-mean(x1,tr)-mean(x2,tr)
test[ic]<-temp$p.value
output[ic,6]<-test[ic]
}}}}}}
ncon<-length(test)
dvec<-alpha/c(1:ncon)
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,7]<-zvec
output[,7]<-output[,7]
INT=bw.es.I(J,K,x,CI=CI,tr=tr,alpha=alpha,REL.MAG=REL.M)$Interaction.ES
list(output=output,Effect.Sizes=INT)
}


bwimcpAKP<-function(J,K,x,tr=.2,JK=J*K,grp=c(1:JK),alpha=.05,SEED=TRUE){
#
#  Same as bwimcp only include Algina et al.  measures of effect size.
#
#
# Multiple comparisons for interactions
# in a split-plot design.
# The analysis is done by taking difference scores
# among all pairs of dependent groups and
# determining which of
# these differences differ across levels of Factor A
# using trimmed means.
#
#  FWE is controlled via Hochberg's method
# To adjusted p-values, use the function p.adjust
#
# For MOM or M-estimators, use spmcpi which uses a bootstrap method
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
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}

JK<-J*K
if(JK!=length(x))stop('Something is wrong. Expected ',JK,' groups but x contains ', length(x), 'groups instead.')
MJ<-(J^2-J)/2
MK<-(K^2-K)/2
JMK<-J*MK
MJMK<-MJ*MK
Jm<-J-1
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
output<-matrix(0,MJMK,8)
dimnames(output)<-list(NULL,c('A','A','B','B','psihat','p.value','p.crit','Effect.Size'))
jp<-1-K
kv<-0
kv2<-0
test<-NA
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
m<-matrix(c(1:JK),J,K,byrow=TRUE)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j<jj){
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
output[ic,1]<-j
output[ic,2]<-jj
output[ic,3]<-k
output[ic,4]<-kk
x1<-x[[m[j,k]]]-x[[m[j,kk]]]
x2<-x[[m[jj,k]]]-x[[m[jj,kk]]]
temp<-yuenv2(x1,x2,SEED=SEED)
output[ic,5]<-mean(x1,tr)-mean(x2,tr)
test[ic]<-temp$p.value
output[ic,6]<-test[ic]
output[ic,8]<-akp.effect(x1,x2)
}}}}}}
ncon<-length(test)
dvec<-alpha/c(1:ncon)
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,7]<-zvec
output[,7]<-output[,7]
output
}


bwbmcp<-function(J,K,x,tr=.2,JK=J*K,grp=c(1:JK),con=0,alpha=.05,dif=TRUE,pool=FALSE,hoch=TRUE,pr=TRUE){
#
# All pairwise comparisons among levels of Factor B
# in a split-plot design using trimmed means.
#
#  If pool=F, levels of Factor A are ignored.
#
# This function calls rmmcp.
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
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}
if(pr)print("Group numbers refer to levels of Factor B")
JK<-J*K
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
if(pool){
data<-list()
m1<-matrix(c(1:JK),J,K,byrow=TRUE)
for(k in 1:K){
for(j in 1:J){
flag<-m1[j,k]
if(j==1)temp<-x[[flag]]
if(j>1){
temp<-c(temp,x[[flag]])
}}
data[[k]]<-temp
}
POOLED<-rmmcp(data,con=con,tr=tr,alpha=alpha,dif=dif,hoch=hoch)
PSI=NULL
A=NULL
}
if(!pool){
A=list()
PSI=list()
mat<-matrix(c(1:JK),ncol=K,byrow=TRUE)
for(j in 1:J){
data<-list()
ic<-0
for(k in 1:K){
ic<-ic+1
data[[ic]]<-x[[mat[j,k]]]
}
temp=rmmcp(data,con=con,tr=tr,alpha=alpha,dif=dif,hoch=hoch)
A[[j]]=temp$test
PSI[[j]]=temp$psihat
POOLED=NULL
}}
list(TESTS.4.EACH.LEVEL.OF.A=A,PSIHAT.4.EACH.LEVEL.OF.A=PSI,POOLED.RESULTS=POOLED)
}


bwmedimcp<-function(J,K,x,JK=J*K,grp=c(1:JK),alpha=.05){
#
# Multiple comparisons for interactions
# in a split-plot design.
# The analysis is done by taking difference scores
# among all pairs of dependent groups and
# determining which of
# these differences differ across levels of Factor A
# using trimmed means.
#
# For MOM or M-estimators, use spmcpi which uses a bootstrap method
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
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}

JK<-J*K
if(JK!=length(x))stop("Something is wrong. Expected ",JK," groups but x contains ", length(x), "groups instead.")
MJ<-(J^2-J)/2
MK<-(K^2-K)/2
JMK<-J*MK
MJMK<-MJ*MK
Jm<-J-1
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x<-data
output<-matrix(0,MJMK,7)
dimnames(output)<-list(NULL,c("A","A","B","B","psihat","sig","crit.sig"))
jp<-1-K
kv<-0
kv2<-0
test<-NA
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
m<-matrix(c(1:JK),J,K,byrow=TRUE)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j<jj){
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
output[ic,1]<-j
output[ic,2]<-jj
output[ic,3]<-k
output[ic,4]<-kk
x1<-x[[m[j,k]]]-x[[m[j,kk]]]
x2<-x[[m[jj,k]]]-x[[m[jj,kk]]]
temp<-qdtest(x1,x2)
output[ic,5]<-median(x1)-median(x2)
test[ic]<-temp$p.value
output[ic,6]<-test[ic]
}}}}}}
ncon<-length(test)
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
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
for (ic in 1:ncol(con)){
output[temp2,7]<-zvec
}
output
}


bwmedbmcp<-function(J,K,x,JK=J*K,grp=c(1:JK),con=0,alpha=.05,dif=FALSE,pool=FALSE,bop=FALSE,nboot=100,SEED=TRUE){
#
# All pairwise comparisons among levels of Factor B
# in a split-plot design using trimmed means.
#
# Data are pooled for each level
# of Factor B.
# bop=T, use bootstrap estimates of standard errors.
# FWE controlled with Rom's method
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
if(pool){
data<-list()
m1<-matrix(c(1:JK),J,K,byrow=TRUE)
for(k in 1:K){
for(j in 1:J){
flag<-m1[j,k]
if(j==1)temp<-x[[flag]]
if(j>1){
temp<-c(temp,x[[flag]])
}}
data[[k]]<-temp
}
print("Group numbers refer to levels of Factor B")
if(!dif)temp<-lincdm(data,con=con,alpha=alpha,nboot=nboot,mop=bop)
if(dif)temp<-qdmcpdif(data,con=con,alpha=alpha)
return(temp)
}
if(!pool){
mat<-matrix(c(1:JK),ncol=K,byrow=TRUE)
for(j in 1:J){
data<-list()
ic<-0
for(k in 1:K){
ic<-ic+1
data[[ic]]<-x[[mat[j,k]]]
}
print(paste("For level ", j, " of Factor A:"))
if(!dif)temp<-lincdm(data,con=con,alpha=alpha,nboot=nboot,mop=bop)
if(dif)temp<-qdmcpdif(data,con=con,alpha=alpha)
print(temp$test)
print(temp$psihat)
}}
}


bwmcpAKP<-function(J,K,x,tr=.2,pr=TRUE){
#
#  Compute Algina et al measure of effect size for all pairwise comparisons
#  in a between-by-within design
#
if(pr){
print('A[[1]] contains the estimated effect size for level 1 of Factor A;')
print(' all pairwise comparisons over Factor B')
print('A[[2]] contains results for level 2, etc.')
}
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
JK=J*K
ID=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
A=list()
for (j in 1:J)A[[j]]=wmcpAKP(x[ID[j,]])
B=list()
for(k in 1:K)B[[k]]=bmcpAKP(x[ID[,k]],tr=tr)
AB=bwimcpAKP(J,K,x)[,c(1:4,8)]
list(Factor.A=A,Factor.B=B,interactions=AB)
}


bwmcpORD<-function(J,K,x,pr=TRUE){
#
#  Compute sign test for within and Cliff's generalization for between factors
#  in a between-by-within design
#
if(pr){
print('A[[1]] contains the estimated effect size for level 1 of Factor A;')
print(' all pairwise comparisons over Factor B')
print('A[[2]] contains results for level 2, etc.')
}
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
JK=j*K
ID=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
A=list()
for (j in 1:J)A[[j]]=signmcp(x[ID[j,]])
B=list()
for(k in 1:K)B[[k]]=cidmulv2(x[ID[,k]])
AB=bwiPATEL(J,K,x)
list(Factor.A=A,Factor.B=B,interactions=AB)
}


bwmcpQS<-function(J,K,x,locfun=median,pr=TRUE){
#
#  Compute quantile shift measure of effect size for all pairwise comparisons
#  in a between-by-within design
#
if(pr){
print('A[[1]] contains results for level 1 of Factor A;')
print(' all pairwise comparisons over Factor B')
print('A[[2]] contains results for level 2, etc.')
print('	')
print('Note: Under normality and homoscedasticity, Cohen d= 0, .2, .5, .8')
print('correspond approximately to Q.effect = 0.5, 0.55, 0.65 and 0.70, respectively')
}
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
JK=J*K
ID=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
A=list()
for (j in 1:J)A[[j]]=wmcpQS(x[ID[j,]],locfun=locfun)
B=list()
for(k in 1:K)B[[k]]=bmcpQS(x[ID[,k]],locfun=locfun)
AB=bwiQS(J,K,x)[,5]
list(Factor.A=A,Factor.B=B,interactions=AB)
}



################################################################################
# BOOTSTRAP MCP FOR MIXED DESIGNS
################################################################################

bwmcppb<-function(J, K, x, est=tmean,JK = J * K,
 alpha = 0.05, grp =c(1:JK), nboot = 500, bhop=TRUE,SEED = TRUE,...)
{
        #
        # A percentile bootstrap for multiple comparisons
        # for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
        #  bhop=TRUE, use Benjaminin--Hochberg. When using a one-step M-estimator
        #  and the sample sizes are small, say less than 80, bhop=TRUE is a bit better.
        #
con=con2way(J,K)
A=bwmcppb.sub(J=J, K=K, x, est=est,JK = J * K,con=con$conA,
 alpha = alpha, grp =c(1:JK), nboot = nboot, bhop=bhop,SEED = SEED,...)
B=bwmcppb.sub(J=J, K=K, x, est=est,JK = J * K,con=con$conB,
 alpha = alpha, grp =c(1:JK), nboot = nboot, bhop=bhop,SEED = SEED,...)
AB=bwmcppb.sub(J=J, K=K, x, est=est,JK = J * K,con=con$conAB,
 alpha = alpha, grp =c(1:JK), nboot = nboot, bhop=bhop,SEED = SEED,...)
list(Fac.A=A,Fac.B=B,Fac.AB=AB)
}


bwmcppb.sub<-function(J, K, x, est=tmean, JK = J * K, con = 0,method='hoch',
 alpha = 0.05, grp =c(1:JK), nboot = 500, bhop=TRUE,SEED = TRUE, ...){
        #
        # A percentile bootstrap for multiple comparisons
        #  for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
x=y
}
        nvec <- NA
for(j in 1:length(x))nvec[j]=length(x[[j]])
nmax=max(nvec)
ncon=ncol(con)
 p <- J * K
        v <- matrix(0, p, p)
        data <- list()
xx=list()
        for(j in 1:length(x)) {
#                data[[j]] <- x[[grp[j]]]
xx[[j]]=x[[grp[j]]] # save input data
#                # Now have the groups in proper order.
#                data[[j]] = data[[j]] - mean(data[[j]], tr = tr)
        }
ilow=1-K
iup=0
for(j in 1:J){
ilow <- ilow + K
 iup = iup + K
sel <- c(ilow:iup)
xx[sel]=listm(elimna(matl(xx[sel])))
                }
        jp <- 1 - K
        kv <- 0
        if(SEED)
                set.seed(2)
        # set seed of random number generator so that
        #             results can be duplicated.
        # Next determine the n_j values
        testA = NA
        bsam = list()
        bdat = list()
aboot=matrix(NA,nrow=nboot,ncol=ncol(con))
tvec=NA
x=xx
tvec=linhat(x,con,est=est,...)
        for(ib in 1:nboot) {
                ilow <- 1 - K
                iup = 0
 for(j in 1:J) {
 ilow <- ilow + K
 iup = iup + K
nv=length(x[[ilow]])
bdat[[j]] = sample(nv, size = nv, replace =TRUE)
for(k in ilow:iup){
 bsam[[k]] = x[[k]][bdat[[j]]]
}
}
ilow=1-K
iup=0
aboot[ib,]=linhat(bsam,con=con,est=est,...)
}
pbA=NA
for(j in 1:ncol(aboot)){
pbA[j]=mean(aboot[,j]>0)
pbA[j]=2*min(c(pbA[j],1-pbA[j]))
}
# Determine critical values
if(!bhop)dvec=alpha/c(1:ncol(con))
if(bhop)dvec<-(ncol(con)-c(1:ncol(con))+1)*alpha/ncol(con)
outputA<-matrix(0,ncol(con),6)
dimnames(outputA)<-list(NULL,c("con.num","psihat","p.value","adj.p.value",
"ci.lower","ci.upper"))
test=pbA
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
outputA[temp2,4]<-zvec
icl<-round(dvec[ncon]*nboot/2)+1
icu<-nboot-icl-1
outputA[,2]<-tvec
for (ic in 1:ncol(con)){
outputA[ic,1]<-ic
outputA[ic,3]<-test[ic]
temp<-sort(aboot[,ic])
outputA[ic,5]<-temp[icl]
outputA[ic,6]<-temp[icu]
}
outputA[,4]=p.adjust(outputA[,3],method=method)
outputA
}


bwmcppb.adj<-function(J, K, x, est=tmean,JK = J * K,method='hoch',
 alpha = 0.05, grp =c(1:JK), nboot = 500, bhop=TRUE,SEED = TRUE,...)
{
        #
        # A percentile bootstrap for multiple comparisons
        # for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
        #  bhop=TRUE, use Benjaminin--Hochberg. When using a one-step M-estimator
        #  and the sample sizes are small, say less than 80, bhop=TRUE is a bit better.
        #
con=con2way(J,K)
A=bwmcppb.sub(J=J, K=K, x, est=est,JK = J * K,con=con$conA,
 alpha = alpha, grp =c(1:JK), nboot = nboot, bhop=bhop,SEED = SEED,...)
B=bwmcppb.sub(J=J, K=K, x, est=est,JK = J * K,con=con$conB,
 alpha = alpha, grp =c(1:JK), nboot = nboot, bhop=bhop,SEED = SEED,...)
AB=bwmcppb.sub(J=J, K=K, x, est=est,JK = J * K,con=con$conAB,
 alpha = alpha, grp =c(1:JK), nboot = nboot, bhop=bhop,SEED = SEED,...)
 A[,4]=p.adjust(A[,3],method=method)
 B[,4]=p.adjust(B[,3],method=method)
 AB[,4]=p.adjust(AB[,3],method=method)
 dimnames(A)=list(NULL,c('con.num', 'psihat','p.value' , 'adjusted.p',' ci.lower' ,'ci.upper'))
 dimnames(B)=list(NULL,c('con.num', 'psihat', 'p.value' , 'adjusted.p',' ci.lower' ,'ci.upper'))
 dimnames(AB)=list(NULL,c('con.num', 'psihat', 'p.value' , 'adjusted.p',' ci.lower' ,'ci.upper'))
list(Fac.A=A,Fac.B=B,Fac.AB=AB)
}


bwwmcppb<-function(J, K,L, x, est=tmean,JKL = J * K*L,
 alpha = 0.05, grp =c(1:JKL), nboot = 2000, bhop=FALSE,SEED = TRUE,...)
{
        #
        # A percentile bootstrap for multiple comparisons
        # for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
con=con3way(J,K,L)
A=bwwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conA,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
B=bwwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conB,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
C=bwwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
AB=bwwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conAB,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
AC=bwwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conAC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
BC=bwwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conBC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
ABC=bwwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conABC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
list(Fac.A=A,Fac.B=B,Fac.C=C,Fac.AB=AB,Fac.AC=AC,Fac.BC=BC,Fac.ABC=ABC)
}


bwwmcppb.sub<-function(J, K,L, x, est=tmean, JKL = J * K*L, con = 0,
 alpha = 0.05, grp =c(1:JKL), nboot = 500, bhop=FALSE,SEED = TRUE, ...){
        #
        # A percentile bootstrap for multiple comparisons
        # for all main effects and interactions.
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
#
#  J independent groups, KL dependent groups
#

       #  If the data are in a matrix, column 1 is assumed to
        #  correspond to x[[1]], column 2 to x[[2]], etc.
        #
        #  When in list mode x is assumed to have length JK, the total number
        #  groups being tested, but a subset of the data can be analyzed
        #  using grp
        #
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
x=y
}
#        nvec <- NA
#for(j in 1:length(x))nvec[j]=length(x[[j]])
ncon=ncol(con)
 p <- J*K*L
if(p>length(x))stop("JKL is less than the Number of groups")
JK=J*K
KL=K*L
#        v <- matrix(0, p, p)
        data <- list()
xx=list()
        for(j in 1:length(x)) {
#                data[[j]] <- x[[grp[j]]]
xx[[j]]=x[[grp[j]]] # save input data
#                # Now have the groups in proper order.
        }
ilow=1-KL
iup=0
for(j in 1:J){
ilow <- ilow + KL
 iup = iup + KL
sel <- c(ilow:iup)
xx[sel]=listm(elimna(matl(xx[sel])))
}

        jp <- 1 - KL
        kv <- 0
        if(SEED)
                set.seed(2)
        # set seed of random number generator so that
        #             results can be duplicated.
        # Next determine the n_j values
        testA = NA
        bsam = list()
        bdat = list()
aboot=matrix(NA,nrow=nboot,ncol=ncol(con))
tvec=NA
x=xx
tvec=linhat(x,con,est=est,...)
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
ilow=1-KL
iup=0
aboot[ib,]=linhat(bsam,con=con,est=est,...)
}
pbA=NA
for(j in 1:ncol(aboot)){
pbA[j]=mean(aboot[,j]>0)
pbA[j]=2*min(c(pbA[j],1-pbA[j]))
}
# Determine critical values
if(!bhop){
if(alpha==.05){
dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncol(con) > 10){
avec<-.05/c(11:(ncol(con)))
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(con > 10){
avec<-.01/c(11:ncol(con))
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncol(con))
}
}
if(bhop)dvec<-(ncol(con)-c(1:ncol(con))+1)*alpha/ncol(con)
outputA<-matrix(0,ncol(con),6)
dimnames(outputA)<-list(NULL,c("con.num","psihat","p.value","p.crit",
"ci.lower","ci.upper"))
test=pbA
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
outputA[temp2,4]<-zvec
icl<-round(dvec[ncon]*nboot/2)+1
icu<-nboot-icl-1
outputA[,2]<-tvec
for (ic in 1:ncol(con)){
outputA[ic,1]<-ic
outputA[ic,3]<-test[ic]
temp<-sort(aboot[,ic])
outputA[ic,5]<-temp[icl]
outputA[ic,6]<-temp[icu]
}
outputA
}


bbmcppb<-function(J, K, x, est=tmean,JK = J*K,
 alpha = 0.05, grp =c(1:JK), nboot = 2000, bhop=FALSE,SEED = TRUE,...)
{
#
#  BETWEEN-BY-BETWEEN DESIGN
#
        # A percentile bootstrap for multiple comparisons
        #  for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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

if(!is.null(dim(x)))x= listm(x)
x=elimna(x)
n=lapply(x,length)
con=con2way(J,K)
A=bbmcppb.sub(J=J, K=K, x, est=est,con=con$conA,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
B=bbmcppb.sub(J=J, K=K, x, est=est,con=con$conB,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
AB=bbmcppb.sub(J=J, K=K, x, est=est,con=con$conAB,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
list(n=n,Fac.A=A,Fac.B=B,Fac.AB=AB)
}


bbwmcppb<-function(J, K,L, x, est=tmean,JKL = J * K*L,
 alpha = 0.05, grp =c(1:JKL), nboot = 2000, bhop=FALSE,SEED = TRUE,...)
{
#
#  BETWEEN-BETWEEN-WITHIN DESIGN
#
        # A percentile bootstrap for multiple comparisons
        #  for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
        #
        #  The R variable x is assumed to contain the raw
        #  data stored in list mode or in a matrix.
        #  If in list mode, x[[1]] contains the data
        #  for the first level of all three factors: level 1,1,1.
        #  x[[2]] is assumed to contain the data for level 1 of the
        #  first two factors and level 3 of the third: level 1,1,2
        #  x[[K]] is the data for level 1,1,K
        #  x[[K+1]] is the data for level 1,2,1, x[[2K]] is level 1,2,K, etc.
        #
        #  If the data are in a matrix, column 1 is assumed to
        #  correspond to x[[1]], column 2 to x[[2]], etc.
        #
        #  When in list mode x is assumed to have length JKL, the total number
        #  groups being tested, but a subset of the data can be analyzed
        #  using grp
        #
con=con3way(J,K,L)
A=bbwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conA,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
B=bbwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conB,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
C=bbwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
AB=bbwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conAB,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
AC=bbwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conAC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
BC=bbwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conBC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
ABC=bbwmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conABC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
list(Fac.A=A,Fac.B=B,Fac.C=C,Fac.AB=AB,Fac.AC=AC,Fac.BC=BC,Fac.ABC=ABC)
}


bbwmcppb.sub<-function(J, K,L, x, est=tmean, JKL = J * K*L, con = 0,
 alpha = 0.05, grp =c(1:JKL), nboot = 500, bhop=FALSE,SEED = TRUE, ...){
        #
        # A percentile bootstrap for multiple comparisons among
        # multiple comparisons for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
        #
        #  The R variable x is assumed to contain the raw
        #  data stored in list mode or in a matrix.
        #  If in list mode, x[[1]] contains the data
        #  for the first level of both factors: level 1,1,1.
        #  x[[2]] is assumed to contain the data for levels 1,1,2, etc.
        #
#
#  JK independent groups, L dependent groups
#

       #  If the data are in a matrix, column 1 is assumed to
        #  correspond to x[[1]], column 2 to x[[2]], etc.
        #
        #  When in list mode x is assumed to have length JKL, the total number
        #  groups being tested, but a subset of the data can be analyzed
        #  using grp
        #
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
x=y
}
ncon=ncol(con)
 p <- J*K*L
if(p>length(x))stop("JKL is less than the Number of groups")
JK=J*K
KL=K*L
        data <- list()
xx=list()
        for(j in 1:length(x)) {
xx[[j]]=x[[grp[j]]] # save input data
#                # Now have the groups in proper order.
        }
ilow=1-L
iup=0
for(j in 1:JK){
ilow <- ilow + L
 iup = iup + L
sel <- c(ilow:iup)
xx[sel]=listm(elimna(matl(xx[sel])))
}

        jp <- 1 - L
        kv <- 0
        if(SEED)
                set.seed(2)
        # set seed of random number generator so that
        #             results can be duplicated.
        # Next determine the n_j values
        testA = NA
        bsam = list()
        bdat = list()
aboot=matrix(NA,nrow=nboot,ncol=ncol(con))
tvec=NA
x=xx
tvec=linhat(x,con,est=est,...)
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
ilow=0-L
iup=0
aboot[ib,]=linhat(bsam,con=con,est=est,...)
}
pbA=NA
for(j in 1:ncol(aboot)){
pbA[j]=mean(aboot[,j]>0)
pbA[j]=2*min(c(pbA[j],1-pbA[j]))
}
# Determine critical values
if(!bhop)dvec=alpha/c(1:ncol(con))
if(bhop)dvec<-(ncol(con)-c(1:ncol(con))+1)*alpha/ncol(con)
outputA<-matrix(0,ncol(con),6)
dimnames(outputA)<-list(NULL,c("con.num","psihat","p.value","p.crit",
"ci.lower","ci.upper"))
test=pbA
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
outputA[temp2,4]<-zvec
icl<-round(dvec[ncon]*nboot/2)+1
icu<-nboot-icl-1
outputA[,2]<-tvec
for (ic in 1:ncol(con)){
outputA[ic,1]<-ic
outputA[ic,3]<-test[ic]
temp<-sort(aboot[,ic])
outputA[ic,5]<-temp[icl]
outputA[ic,6]<-temp[icu]
}
outputA
}


bbbmcppb<-function(J, K,L, x, est=tmean,JKL = J * K*L,
 alpha = 0.05, grp =c(1:JKL), nboot = 2000, bhop=FALSE,SEED = TRUE,...)
{
#
#  BETWEEN-BETWEEN-BETWEEN DESIGN
#
        # A percentile bootstrap for multiple comparisons among
        # multiple comparisons for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
con=con3way(J,K,L)
A=bbbmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conA,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
B=bbbmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conB,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
C=bbbmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
AB=bbbmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conAB,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
AC=bbbmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conAC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
BC=bbbmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conBC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
ABC=bbbmcppb.sub(J=J, K=K,L=L, x, est=est,con=con$conABC,
 alpha = alpha, nboot = nboot, bhop=bhop,SEED = SEED,grp=grp,...)
list(Fac.A=A,Fac.B=B,Fac.C=C,Fac.AB=AB,Fac.AC=AC,Fac.BC=BC,Fac.ABC=ABC)
}


 bbbmcppb.sub<-function(J, K,L, x, est=tmean, JKL = J * K*L, con = 0,
 alpha = 0.05, grp =c(1:JKL), nboot = 500, bhop=FALSE,SEED = TRUE, ...){
#
#   between-by-between-by-between design
#
        #
        # A percentile bootstrap for
        # multiple comparisons for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using and appropriate linear contrast.
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
#
        #  If the data are in a matrix, column 1 is assumed to
        #  correspond to x[[1]], column 2 to x[[2]], etc.
        #
        #  When in list mode x is assumed to have length JKL, the total number
        #  groups being tested, but a subset of the data can be analyzed
        #  using grp
        #
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
x=y
}
ncon=ncol(con)
 p <- J*K*L
JKL=p
if(p>length(x))stop('JKL is less than the Number of groups')
JK=J*K
KL=K*L
        data <- list()
xx=list()
        for(j in 1:length(x)) {
xx[[j]]=x[[grp[j]]] # save input data
#                # Now have the groups in proper order.
        }
for(j in 1:p){
xx[j]=elimna(xx[j])
}
  crit=alpha/2
 icl<-round(crit*nboot)+1
icu<-nboot-icl
        if(SEED)
                set.seed(2)
        # set seed of random number generator so that
        #             results can be duplicated.
        # Next determine the n_j values
        testA = NA
        bsam = list()
        bdat = list()
aboot=matrix(NA,nrow=nboot,ncol=ncol(con))
tvec=NA
tvec=linhat(x,con,est=est,...)
        for(ib in 1:nboot) {
 for(j in 1:JKL) {
nv=length(x[[j]])
bdat[[j]] = sample(nv, size = nv, replace =TRUE)
bsam[[j]] = x[[j]][bdat[[j]]]
}
aboot[ib,]=linhat(bsam,con=con,est=est,...)
}
pbA=NA
for(j in 1:ncol(aboot)){
pbA[j]=mean(aboot[,j]>0)
pbA[j]=2*min(c(pbA[j],1-pbA[j]))
}
outputA<-matrix(0,ncol(con),6)
dimnames(outputA)<-list(NULL,c('con.num','psihat','p.value','p.adjust',
'ci.lower','ci.upper'))
test=pbA
outputA[,2]<-tvec
for (ic in 1:ncol(con)){
outputA[ic,1]<-ic
outputA[ic,3]<-test[ic]
temp<-sort(aboot[,ic])
outputA[ic,5]<-temp[icl]
outputA[ic,6]<-temp[icu]
}
outputA[,4]=p.adjust(outputA[,3],method='hoch')
outputA
}


wwwmcppb<-function(J,K,L,x, alpha = 0.05, con = 0,est=tmean, plotit = FALSE,
    dif = TRUE, grp = NA, nboot = NA, BA = TRUE, hoch = TRUE, xlab = "Group 1",
    ylab = "Group 2", pr = TRUE, SEED = TRUE,...){
#
# Do all multiple comparisons for a within-by-within-by-within design.
# using a percentile bootstrap method and trimmed means
#
if(pr){
print('This new version includes the option to use difference scores and defaults to dif=TRUE')
print('Number of bootstrap samples differs from the old version')
print('To use the old version, use wwwmcppb.OLD')
}
conM=con3way(J,K,L)
A=rmmcppb(x,con=conM$conA,alpha=alpha,dif=dif,plotit=plotit,est=est,
nboot=nboot,BA=BA,hoch=hoch,SEED=SEED,xlab=xlab,ylab=ylab,pr=pr,...)
B=rmmcppb(x,con=conM$conB,alpha=alpha,dif=dif,
plotit=plotit,est=est,nboot=nboot,BA=BA,hoch=hoch,
SEED=SEED,xlab=xlab,ylab=ylab,pr=FALSE,...)
C=rmmcppb(x,con=conM$conC,alpha=alpha,dif=dif,
plotit=plotit,est=est,nboot=nboot,BA=BA,hoch=hoch,
SEED=SEED,xlab=xlab,ylab=ylab,pr=FALSE,...)
AB=rmmcppb(x,con=conM$conAB,alpha=alpha,dif=dif,plotit=plotit,est=est,
nboot=nboot,BA=BA,hoch=hoch,SEED=SEED,xlab=xlab,ylab=ylab,pr=FALSE,...)
AC=rmmcppb(x,con=conM$conAC,alpha=alpha,dif=dif,plotit=plotit,est=est,
nboot=nboot,BA=BA,hoch=hoch,SEED=SEED,xlab=xlab,ylab=ylab,pr=FALSE,...)
BC=rmmcppb(x,con=conM$conBC,alpha=alpha,dif=dif,plotit=plotit,est=est,
nboot=nboot,BA=BA,hoch=hoch,SEED=SEED,xlab=xlab,ylab=ylab,pr=FALSE,...)
ABC=rmmcppb(x,con=conM$conABC,alpha=alpha,dif=dif,plotit=plotit,est=est,
nboot=nboot,BA=BA,hoch=hoch,SEED=SEED,xlab=xlab,ylab=ylab,pr=FALSE,...)
list(Factor_A=A,Factor_B=B,Factor_AB=AB,Factor_AC=AC,Factor_BC=BC,Factor_ABC=ABC)
}


wwwmcppb.sub<-function(J, K,L, x, est=tmean, JKL = J * K*L, con = 0,
 alpha = 0.05, grp =c(1:JKL), nboot = 500, bhop=FALSE,SEED = TRUE, ...){
        #
        # A percentile bootstrap for multiple comparisons among
        # multiple comparisons for all main effects and interactions
        # The analysis is done by generating bootstrap samples and
        # using an appropriate linear contrast.
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
#
#  within-by-within-by-within design
#
#  JKL dependent groups
#

       #  If the data are in a matrix, column 1 is assumed to
        #  correspond to x[[1]], column 2 to x[[2]], etc.
        #
        #  When in list mode x is assumed to have length JK, the total number
        #  groups being tested, but a subset of the data can be analyzed
        #  using grp
        #
        if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
x=y
}
#        nvec <- NA
#for(j in 1:length(x))nvec[j]=length(x[[j]])
ncon=ncol(con)
 p <- J*K*L
JKL=p
if(p>length(x))stop("JKL is less than the Number of groups")
JK=J*K
KL=K*L
#        v <- matrix(0, p, p)
        data <- list()
xx=list()
        for(j in 1:length(x)) {
#                data[[j]] <- x[[grp[j]]]
xx[[j]]=x[[grp[j]]] # save input data
#                # Now have the groups in proper order.
        }
        if(SEED)
                set.seed(2)
        # set seed of random number generator so that
        #             results can be duplicated.
        # Next determine the n_j values
        testA = NA
        bsam = list()
        bdat = list()
aboot=matrix(NA,nrow=nboot,ncol=ncol(con))
tvec=NA
x=xx
tvec=linhat(x,con,est=est,...)
nv=length(x[[1]])
        for(ib in 1:nboot) {
bdat[[j]] = sample(nv, size = nv, replace =TRUE)
for(k in 1:JKL) bsam[[k]] = x[[k]][bdat[[j]]]
aboot[ib,]=linhat(bsam,con=con,est=est,...)
}
pbA=NA
for(j in 1:ncol(aboot)){
pbA[j]=mean(aboot[,j]>0)
pbA[j]=2*min(c(pbA[j],1-pbA[j]))
}
# Determine critical values
if(!bhop){
if(alpha==.05){
dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncol(con) > 10){
avec<-.05/c(11:(ncol(con)))
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(con > 10){
avec<-.01/c(11:ncol(con))
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncol(con))
}
}
if(bhop)dvec<-(ncol(con)-c(1:ncol(con))+1)*alpha/ncol(con)
outputA<-matrix(0,ncol(con),6)
dimnames(outputA)<-list(NULL,c("con.num","psihat","p.value","p.crit",
"ci.lower","ci.upper"))
test=pbA
temp2<-order(0-test)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
outputA[temp2,4]<-zvec
icl<-round(dvec[ncon]*nboot/2)+1
icu<-nboot-icl-1
outputA[,2]<-tvec
for (ic in 1:ncol(con)){
outputA[ic,1]<-ic
outputA[ic,3]<-test[ic]
temp<-sort(aboot[,ic])
outputA[ic,5]<-temp[icl]
outputA[ic,6]<-temp[icu]
}
outputA
}

wwwmcppbtr<-function(J,K,L, x,tr=.2,alpha=.05,dif=TRUE,op=FALSE,grp=NA,nboot=2000,SEED=TRUE,pr=TRUE){
#
#  Based on a percentile bootstrap method.
#
#  dif=TRUE: use a  linear  combination of the variables, test the hypothesis that the trimmed mean is zero
#  dif=FALSE: Use the marginal trimmed means instead.
#
# MULTIPLE COMPARISONS FOR A 3-WAY within-by-within-by within ANOVA
# Do all multiple comparisons associated with
# main effects for Factor A and B and C and all interactions
# based on trimmed means
#
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
        JKL <- J*K*L
        if(is.matrix(x))
                x <- listm(x)
        if(!is.na(grp[1])) {
                yy <- x
                x<-list()
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop("Data must be stored in list mode or a matrix.")
        for(j in 1:JKL) {
                xx <- x[[j]]
                x[[j]] <- xx[!is.na(xx)] # Remove missing values
        }
        #

        if(JKL != length(x))
                warning("The number of groups does not match the number of contrast coefficients.")
for(j in 1:JKL){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
        # Create the three contrast matrices
temp<-con3way(J,K,L)
conA<-temp$conA
conB<-temp$conB
conC<-temp$conC
conAB<-temp$conAB
conAC<-temp$conAC
conBC<-temp$conBC
conABC<-temp$conABC
Factor.A<-rmmcppb(x,con=conA,est=tmean,tr=tr,alpha=alpha,dif=dif,nboot=nboot,SEED=SEED,pr=pr)
Factor.B<-rmmcppb(x,con=conB,est=tmean,tr=tr,alpha=alpha,dif=dif,nboot=nboot,SEED=SEED,pr=FALSE)
Factor.C<-rmmcppb(x,con=conC,est=tmean,tr=tr,alpha=alpha,dif=dif,nboot=nboot,SEED=SEED,pr=FALSE)
Factor.AB<-rmmcppb(x,con=conAB,est=tmean,tr=tr,alpha=alpha,dif=dif,nboot=nboot,SEED=SEED,pr=FALSE)
Factor.AC<-rmmcppb(x,con=conAC,est=tmean,tr=tr,alpha=alpha,dif=dif,nboot=nboot,SEED=SEED,pr=FALSE)
Factor.BC<-rmmcppb(x,con=conBC,est=tmean,tr=tr,alpha=alpha,dif=dif,nboot=nboot,SEED=SEED,pr=FALSE)
Factor.ABC<-rmmcppb(x,con=conABC,est=tmean,tr=tr,alpha=alpha,dif=dif,nboot=nboot,SEED=SEED,pr=FALSE)
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.C=Factor.C,
Factor.AB=Factor.AB,Factor.AC=Factor.AC,Factor.BC=Factor.BC,
Factor.ABC=Factor.ABC,conA=conA,conB=conB,conC=conC,
conAB=conAB,conAC=conAC,conBC=conBC,conABC=conABC)
}

wmcppb<-function(x, y=NULL,alpha = 0.05, con = 0,est=tmean, plotit = FALSE,
    dif = TRUE, grp = NA, nboot = NA, BA = TRUE, hoch = TRUE, xlab = "Group 1",
    ylab = "Group 2", pr = TRUE, SEED = TRUE, ...){
#
# Do all multiple comparisons for a repeated measures design.
# using a percentile bootstrap method and trimmed means
#
if(!is.null(y))x=cbind(x,y)
A=rmmcppb(x,con=con,alpha=alpha,dif=dif,plotit=plotit,est=est,
nboot=nboot,BA=BA,hoch=hoch,SEED=SEED,xlab=xlab,ylab=ylab,pr=pr,...)
A
}


wwmcppb<-function(J,K,x, alpha = 0.05, con = 0,est=tmean, plotit = FALSE,
    dif = TRUE, grp = NA, nboot = NA, BA = TRUE, hoch = TRUE, xlab = "Group 1",
    ylab = "Group 2", pr = TRUE, SEED = TRUE,...){
#
# Do all multiple comparisons for a within-by-within design.
# using a percentile bootstrap method and trimmed means
#
conM=con2way(J,K)
A=rmmcppb(x,con=conM$conA,alpha=alpha,dif=dif,plotit=plotit,est=est,
nboot=nboot,BA=BA,hoch=hoch,SEED=SEED,xlab=xlab,ylab=ylab,pr=pr,...)
B=rmmcppb(x,con=conM$conB,alpha=alpha,dif=dif,
plotit=plotit,est=est,nboot=nboot,BA=BA,hoch=hoch,
SEED=SEED,xlab=xlab,ylab=ylab,pr=pr,...)
AB=rmmcppb(x,con=conM$conAB,alpha=alpha,dif=dif,plotit=plotit,est=est,
nboot=nboot,BA=BA,hoch=hoch,SEED=SEED,xlab=xlab,ylab=ylab,pr=pr,...)
list(Factor_A=A,Factor_B=B,Factor_AB=AB)
}


wwmcpbt<-function(J,K,x, tr=.2, dif=TRUE, alpha = 0.05, nboot = 599){
#
# Do multiple comparisons for a within-by-within design.
# using a bootstrap-t method and trimmed means.
# All linear contrasts relevant to main effects and interactions
# are tested.
#
#
conM=con2way(J,K)
A=lindepbt(x,con=conM$conA,alpha=alpha,tr=tr,dif=dif,nboot=nboot)
B=lindepbt(x,con=conM$conB,alpha=alpha,tr=tr,dif=dif,nboot=nboot)
AB=lindepbt(x,con=conM$conAB,alpha=alpha,tr=tr,dif=dif,nboot=nboot)
list(Factor_A=A,Factor_B=B,Factor_AB=AB)
}


wwmcpES<-function(J,K,x,tr=.2,alpha=.05,dif=TRUE){
#
# Do all multiple comparisons for a within-by-within design
# using trimmed means
#
stop('Use ww.es instead')
conM=con2way(J,K)
A=rmmcpES(x,con=conM$conA,tr=tr,alpha=alpha,dif=dif)
B=rmmcpES(x,con=conM$conB,tr=tr,alpha=alpha,dif=dif)
AB=rmmcpES(x,con=conM$conAB,tr=tr,alpha=alpha,dif=dif)
list(Factor_A=A,Factor_B=B,Factor_AB=AB)
}


wwmcpQS<-function(J,K,x,tr=.2,alpha=.05,dif=TRUE,locfun=tmean,...){
#
# Do all multiple comparisons for a within-by-within design
# using trimmed means
#
conM=con2way(J,K)
A=rmmcpQS(x,con=conM$conA,tr=tr,alpha=alpha,dif=dif,locfun=locfun)
B=rmmcpQS(x,con=conM$conB,tr=tr,alpha=alpha,dif=dif,locfun=locfun)
AB=rmmcpQS(x,con=conM$conAB,tr=tr,alpha=alpha,dif=dif,locfun=locfun)
list(Factor_A=A,Factor_B=B,Factor_AB=AB)
}


wwwmcpQS<-function(J,K,L, x,tr=.2,alpha=.05,dif=TRUE,op=FALSE,grp=NA,locfun=tmean,...){
#
# MULTIPLE COMPARISONS FOR A 3-WAY within-by-within-by within ANOVA
#
#  Include quantile shift measure of location.
#
# Do all multiple comparisons associated with
# main effects for Factor A and B and C and all interactions
# based on trimmed means
#
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
        JKL <- J*K*L
        if(is.matrix(x))
                x <- listm(x)
        if(!is.na(grp[1])) {
                yy <- x
                x<-list()
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop("Data must be stored in list mode or a matrix.")
        for(j in 1:JKL) {
                xx <- x[[j]]
                x[[j]] <- xx[!is.na(xx)] # Remove missing values
        }
        #

        if(JKL != length(x))
                warning("The number of groups does not match the number of contrast coefficients.")
for(j in 1:JKL){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
        # Create the three contrast matrices
temp<-con3way(J,K,L)
conA<-temp$conA
conB<-temp$conB
conC<-temp$conC
conAB<-temp$conAB
conAC<-temp$conAC
conBC<-temp$conBC
conABC<-temp$conABC
Factor.A<-rmmcpQS(x,con=conA,tr=tr,alpha=alpha,dif=dif,locfun=locfun,...)
Factor.B<-rmmcpQS(x,con=conB,tr=tr,alpha=alpha,dif=dif,locfun=locfun,...)
Factor.C<-rmmcpQS(x,con=conC,tr=tr,alpha=alpha,dif=dif,locfun=locfun,...)
Factor.AB<-rmmcpQS(x,con=conAB,tr=tr,alpha=alpha,dif=dif,locfun=locfun,...)
Factor.AC<-rmmcpQS(x,con=conAC,tr=tr,alpha=alpha,dif=dif,locfun=locfun,...)
Factor.BC<-rmmcpQS(x,con=conBC,tr=tr,alpha=alpha,dif=dif,locfun=locfun,...)
Factor.ABC<-rmmcpQS(x,con=conABC,tr=tr,alpha=alpha,dif=dif,locfun=locfun,...)
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.C=Factor.C,
Factor.AB=Factor.AB,Factor.AC=Factor.AC,Factor.BC=Factor.BC,
Factor.ABC=Factor.ABC,conA=conA,conB=conB,conC=conC,
conAB=conAB,conAC=conAC,conBC=conBC,conABC=conABC)
}



################################################################################
# SPLIT-PLOT DESIGNS
################################################################################

spmcpa<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),avg=FALSE,alpha=.05,SEED=TRUE,
nboot=NA,pr=TRUE,...){
#
# A percentile bootstrap for multiple comparisons among
# all main effects for independent groups in a split-plot design
# The analysis is done by generating bootstrap samples and
# using an appropriate linear contrast.
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
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}
if(pr)print("As of Sept. 2005, est defaults to tmean")
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
d<-ncol(con)
if(is.na(nboot)){
if(d<=4)nboot<-1000
if(d>4)nboot<-5000
}
#
# Now take bootstrap samples from jth level
# of Factor A and average K  corresponding estimates
# of location.
#
bloc<-matrix(NA,nrow=J,ncol=nboot)
print("Taking bootstrap samples. Please wait.")
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
mat<-listm(x)
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
if(avg)bloc<-t(bloc)
connum<-d
psihat<-matrix(0,connum,nboot)
test<-1
for (ic in 1:connum){
psihat[ic,]<-apply(bloc,1,bptdpsi,con[,ic])
#test[ic]<-sum((psihat[ic,]>0))/nboot
test[ic]<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
test[ic]<-min(test[ic],1-test[ic])
}
ncon<-ncol(con)
if(alpha==.05){
dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
dvec[1]<-alpha/2
}
temp2<-order(0-test)
ncon<-ncol(con)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output<-matrix(0,connum,6)
dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.sig","ci.lower","ci.upper"))
tmeans<-mvec
psi<-1
output[temp2,4]<-zvec
for (ic in 1:ncol(con)){
output[ic,2]<-sum(con[,ic]*tmeans)
output[ic,1]<-ic
output[ic,3]<-test[ic]
temp<-sort(psihat[ic,])
temp3<-round(output[ic,4]*nboot)+1
icl<-round(dvec[ncon]*nboot)+1
icu<-nboot-(icl-1)
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
output[,3]<-2*output[,3]
output[,4]<-2*output[,4]
num.sig<-sum(output[,3]<=output[,4])
list(output=output,con=con,num.sig=num.sig)
}


spmcpi<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),alpha=.05,nboot=NA,
SEED=TRUE,pr=TRUE,SR=FALSE,...){
#
# Multiple comparisons for interactions
# in a split-plot design.
# The analysis is done by taking difference scores
# among all pairs of dependent groups and
# determining which of
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
#  When in list mode x is assumed to have length JK, the total number of
#  groups being tested, but a subset of the data can be analyzed
#  using grp
#
       if(is.matrix(x)) {
                y <- list()
                for(j in 1:ncol(x))
                        y[[j]] <- x[, j]
                x <- y
}
if(pr)print("As of Sept. 2005, est defaults to tmean")
JK<-J*K
if(JK!=length(x)){
print("Something is wrong.")
paste(" Expected ",JK," groups but x contains ", length(x), "groups instead.")
stop()
}
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
d<-ncol(con)
if(is.na(nboot)){
if(d<=4)nboot<-1000
if(d>4)nboot<-5000
}
connum<-d
psihat<-matrix(0,connum,nboot)
#
# Now take bootstrap samples from jth level
# of Factor A and average K  corresponding estimates
# of location.
#
bloc<-matrix(NA,ncol=J,nrow=nboot)
print("Taking bootstrap samples. Please wait.")
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
test<-1
for (ic in 1:connum){
psihat[ic,]<-apply(bloc,1,bptdpsi,con[,ic])
#test[ic]<-sum((psihat[ic,]>0))/nboot
test[ic]<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
test[ic]<-min(test[ic],1-test[ic])
}
ncon<-ncol(con)
dvec<-alpha/c(1:ncon)
if(SR){
okay=FALSE
if(identical(est,onestep))okay=TRUE
if(identical(est,mom))okay=TRUE
if(!okay)stop('For estimators other than onestep and mom, use SR=FALSE')
if(alpha==.05){
dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
if(ncon > 10){
avec<-.05/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha==.01){
dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
if(ncon > 10){
avec<-.01/c(11:ncon)
dvec<-c(dvec,avec)
}}
if(alpha != .05 && alpha != .01){
dvec<-alpha/c(1:ncon)
dvec[1]<-alpha/2
}}
temp2<-order(0-test)
ncon<-ncol(con)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output<-matrix(0,connum,6)
dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.crit","ci.lower","ci.upper"))
tmeans<-mvec
psi<-1
for (ic in 1:ncol(con)){
output[ic,2]<-sum(con[,ic]*tmeans)
output[ic,1]<-ic
output[ic,3]<-test[ic]
output[temp2,4]<-zvec
temp<-sort(psihat[ic,])
icl<-round(dvec[ncon]*nboot)+1
icu<-nboot-(icl-1)
output[ic,5]<-temp[icl]
output[ic,6]<-temp[icu]
}
output[,3]<-2*output[,3]
if(SR)output[,4]<-2*output[,4]
num.sig<-sum(output[,3]<=output[,4])
list(output=output,con=con,num.sig=num.sig)
}


spmcpb<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),dif=TRUE,alpha=.05,SEED=TRUE,
nboot=NA,...){
#
# A percentile bootstrap for all pairwise
# multiple comparisons
# among dependent groups in a split-plot design
#
#  Levels of A are ignored.
#
# If dif=T, the analysis is done based on all pairs
# of difference scores.
# Otherwise, marginal measures of location are used.
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
 if(is.matrix(x) || is.data.frame(x))x=listm(x)
JK<-J*K
data<-list()
for(j in 1:length(x)){
data[[j]]<-x[[grp[j]]] # Now have the groups in proper order.
}
x=data
x=pool.fun(J,K,x)
temp<-rmmcppb(x,est=est,nboot=nboot,dif=dif,alpha=alpha,plotit=FALSE,SEED=SEED,...)
list(output=temp$output,con=temp$con,num.sig=temp$num.sig)
}


spmcpbA<-function(J,K,x,est=tmean,JK=J*K,grp=c(1:JK),dif=TRUE,alpha=.05,
nboot=NA,SEED=TRUE,...){
#
# For each level of Factor A
# use a percentile bootstrap for all pairwise
# multiple comparisons
# among dependent groups in a split-plot design
#
#
# If dif=T, the analysis is done based on all pairs
# of difference scores.
# Otherwise, marginal measures of location are used.
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

       if(is.matrix(x) || is.data.frame(x)) {
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
A=list()
imat=matrix(c(1:JK),nrow=J,byrow=TRUE)
for(j in 1:J){
# Now call function rmmcppb to do the analysis
id=imat[j,]
A[[j]]<-rmmcppb(x[id],est=est,pr=FALSE,nboot=nboot,dif=dif,alpha=alpha,plotit=FALSE,SEED=SEED,...)
}
list(A.Level=A)
}


qdmcp<-function(x,alpha = 0.05,bop=FALSE,nboot=100,pr=TRUE,q=.5,SEED=TRUE){
#
# For dependent groups,
# Perform all pairwise comparisons
# using quantiles estimated with a single order statistic.
#  FWE controlled with Rom's method
#
if(is.data.frame(x))x=as.matrix(x)
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
J<-ncol(x)
xbar<-vector("numeric",J)
x<-elimna(x)  # Remove missing values
df<-nrow(x)-1
nval<-nrow(x)
for(j in 1: J){
if(!bop)xbar[j]<-qest(x[,j],q=q)
if(bop)xbar[j]<-median(x[,j])
}
CC<-(J^2-J)/2
ncon<-CC
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
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
psihat<-matrix(0,CC,5)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c("Group","Group","test","p-value","p.crit","se"))
if(bop)se.val<-bootdse(x,nboot=nboot,pr=pr)
temp1<-0
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
if(!bop)temp<-qdtest(x[,j],x[,k],q=q,bop=bop)
if(bop)temp<-qdtest(x[,j],x[,k],se.val=se.val[jcom])
sejk<-temp$se
test[jcom,6]<-sejk
test[jcom,3]<-temp$test.stat
test[jcom,4]<-temp$p.value
if(length(x[,j])<20)test[jcom,4]<-mrm1way(x[,c(j,k)],q=q,SEED=SEED)$p.value
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-(xbar[j]-xbar[k])
}}}
temp1<-test[,4]
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
test[temp2,5]<-zvec
psihat[,4]<-psihat[,3]-qt(1-test[,5]/2,df)*test[,6]
psihat[,5]<-psihat[,3]+qt(1-test[,5]/2,df)*test[,6]
num.sig<-sum(test[,4]<=test[,5])
list(test=test,psihat=psihat,num.sig=num.sig)
}


qdmcpdif<-function(x, con = 0,alpha = 0.05){
#
# MCP with medians on difference scores
# FWE controlled with Rom's method
#
if(is.data.frame(x))x=as.matrix(x)
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-ncol(x)
xbar<-vector("numeric",J)
x<-elimna(x)  # Remove missing values
nval<-nrow(x)
h1<-nrow(x)
df<-h1-1
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0)CC<-(J^2-J)/2
ncon<-CC
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
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
if(sum(con^2)==0){
psihat<-matrix(0,CC,5)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
test<-matrix(NA,CC,5)
dimnames(test)<-list(NULL,c("Group","Group","p-value","p.crit","se"))
temp1<-0
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
dv<-x[,j]-x[,k]
test[jcom,5]<-msmedse(dv)
temp<-sintv2(dv,alpha=alpha/CC)
temp1[jcom]<-temp$p.value
test[jcom,3]<-temp$p.value
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-median(dv)
psihat[jcom,4]<-temp$ci.low
psihat[jcom,5]<-temp$ci.up
}}}
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
zvec[ddd:ncon]<-dvec[ddd]
}
test[temp2,4]<-zvec
}
if(sum(con^2)>0){
if(nrow(con)!=ncol(x))print("WARNING: The number of groups does not match the number of contrast coefficients.")
ncon<-ncol(con)
psihat<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
test<-matrix(0,ncol(con),4)
dimnames(test)<-list(NULL,c("con.num","sig","crit.sig","se"))
temp1<-NA
for (d in 1:ncol(con)){
psihat[d,1]<-d
for(j in 1:J){
if(j==1)dval<-con[j,d]*x[,j]
if(j>1)dval<-dval+con[j,d]*x[,j]
}
temp3<-sintv2(dval)
temp1[d]<-temp3$p.value
test[d,1]<-d
test[d,4]<-msmedse(dval)
psihat[d,2]<-median(dval)
psihat[d,3]<-temp3$ci.low
psihat[d,4]<-temp3$ci.up
}
test[,2]<-temp1
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
print(c(ncon,zvec))
sigvec<-(test[temp2,2]>=zvec)
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
zvec[ddd:ncon]<-dvec[ddd]
}
test[temp2,3]<-zvec
}
if(sum(con^2)==0)num.sig<-sum(psihat[,4]>0)+ sum(psihat[,5]<0)
if(sum(con^2)>0)num.sig<-sum(psihat[,3]>0)+ sum(psihat[,4]<0)
list(test=test,psihat=psihat,con=con,num.sig=num.sig)
}


mwwmcp<-function(J,K,x,grp=c(1:p),p=J*K,q=.5,bop=FALSE,alpha=.05,nboot=100,
SEED=TRUE){
#
#  For a J by K anova using quantiles with
#  repeated measures on both factors,
#  Perform all multiple comparisons for main effects
#  and interactions.
#
#  q=.5 by default meaning medians are compared
#  bop=F means bootstrap option not used;
#  with bop=T, function uses usual medians rather
#  rather than a single order statistic to estimate median
#  in conjunction with a bootstrap estimate of covariances
#  among sample medians.
#
#  The R variable data is assumed to contain the raw
#  data stored in a matrix or in list mode.
#  When in list mode data[[1]] contains the data
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
Qa<-NA
Qab<-NA
if(is.data.frame(x))x=as.matrix(x)
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
tmeans<-0
        # Create the three contrast matrices
        #
Ja<-(J^2-J)/2
Ka<-(K^2-K)/2
JK<-J*K
conA<-matrix(0,nrow=JK,ncol=Ja)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j < jj){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[j,]<-1
mat[jj,]<-0-1
conA[,ic]<-t(mat)
}}}
conB<-matrix(0,nrow=JK,ncol=Ka)
ic<-0
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[,k]<-1
mat[,kk]<-0-1
conB[,ic]<-t(mat)
}}}
conAB<-matrix(0,nrow=JK,ncol=Ka*Ja)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j < jj){
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[j,k]<-1
mat[j,kk]<-0-1
mat[jj,k]<-0-1
mat[jj,kk]<-1
}
conAB[,ic]<-t(mat)
}}}}}
Qa<-lincdm(x,con=conA,alpha=alpha,mop=bop,nboot=nboot,SEED=SEED)
# Do test for factor B
Qb<-lincdm(x,con=conB,alpha=alpha,mop=bop,nboot=nboot,SEED=SEED)
# Do test for factor A by B interaction
Qab<-lincdm(x,con=conAB,alpha=alpha,mop=bop,nboot=nboot,SEED=SEED)
list(Qa=Qa,Qb=Qb,Qab=Qab)
}


twwmcp<-function(J,K,x,grp=c(1:p),p=J*K,tr=.2,alpha=.05,dif=FALSE){
#
#  For a J by K anova using quantiles with
#  repeated measures on both factors,
#  Perform all multiple comparisons for main effects
#  and interactions.
#
#  tr=.2. default trimming
#  bop=F means bootstrap option not used;
#  with bop=T, function uses usual medians rather
#  rather than a single order statistic to estimate median
#  in conjunction with bootstrap estimate of covariances
#  among the sample medians.
#
#  The R variable data is assumed to contain the raw
#  data stored in a matrix or in list mode.
#  When in list mode data[[1]] contains the data
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
Qa<-NA
Qab<-NA
if(is.list(x))x<-elimna(matl(x))
if(is.data.frame(x))x=as.matrix(x)
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
tmeans<-0
temp<-con2way(J,K) # contrasts matrices stored in temp
Qa<-rmmcp(x,con=temp$conA,alpha=alpha,dif=dif,tr=tr)
# Do test for factor B
Qb<-rmmcp(x,con=temp$conB,alpha=alpha,dif=dif,tr=tr)
# Do test for factor A by B interaction
Qab<-rmmcp(x,con=temp$conAB,alpha=alpha,dif=dif,tr=tr)
list(Qa=Qa,Qb=Qb,Qab=Qab)
}


tkmcp<-function(x,alpha=.05,ind.pval=TRUE){
#
# conventional Tukey-Kramer multiple comparison procedure
# for all pairiwise comparisons.
#
#  ind.pval=T, computes p-value for each individual test
#  ind.pval=F computes p-value based on controlling the
#  familywise error rate. (The alpha level at which the
#  Tukey-Kramer test would reject.)
#
if(is.matrix(x))x<-listm(x)
J<-length(x)
A<-0
B<-0
C<-0
N<-0
for(j in 1:J){
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
numcom<-length(x)*(length(x)-1)/2
output<-matrix(nrow=numcom,ncol=7)
dimnames(output)<-list(NULL,c("Group","Group","t.test","est.difference",
"ci.lower","ci.upper","p.value"))
ic<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
ic<-ic+1
output[ic,1]<-j
output[ic,2]<-k
dif<-mean(x[[j]])-mean(x[[k]])
output[ic,3]<-abs(dif)/sqrt(MSWG*(1/length(x[[j]])+1/length(x[[k]]))/2)
output[ic,4]<-dif
crit<-qtukey(1-alpha,length(x),nu2)
output[ic,5]<-dif-crit*sqrt(MSWG*(1/length(x[[j]])+1/length(x[[k]]))/2)
output[ic,6]<-dif+crit*sqrt(MSWG*(1/length(x[[j]])+1/length(x[[k]]))/2)
if(!ind.pval)output[ic,7]<-1-ptukey(output[ic,3],length(x),nu2)
if(ind.pval)output[ic,7]<-2*(1-pt(output[ic,3],nu2))
}}}
output
}



################################################################################
# EFFECT SIZE MCP
################################################################################

esmcp<-function(x,tr=.2,method='EP',pr=FALSE,SEED=TRUE){
#
#  A one-way design is assumed, independent groups
#
#  For all pairs of groups,
#  compute a heteroscedastic robust measure of effect size
#  using one of six methods, indicated by the argument method:
#
#  EP: explanatory measure of effect size
#  QS: a median-type quantile shift measure of effect size
#  QStr:  a trimmed mean  quantile shift measure of effect size
#  AKP:  trimmed mean Winsorized variance analog of Cohen's d
#  WMW:  P(X<Y)
#  KMS: Kulinskaya et al. method
#
#
#  OUTPUT: effect size for all pairs of groups
#
if(is.data.frame(x))x=as.matrix(x)
if(SEED)set.seed(2)
if(is.matrix(x))x<-listm(x)
J<-length(x)
JALL=(J^2-J)/2
est=matrix(NA,JALL,3)
dimnames(est)=list(NULL,c("Group","Group","Effect Size"))
ic=0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic=ic+1
est[ic,1]=j
est[ic,2]=k
est[ic,3]=ESfun(x[[j]],x[[k]],tr=tr,method=method,pr=pr)
}}}
list(Estimates=est)
}


bbmcpEP<-function(J,K,x,tr=.2,alpha=.05,grp=NA,op=FALSE,nreps=200,SEED=TRUE,pr=TRUE,POOL=TRUE){
#
#  Test all linear contrasts associated with
# main effects for Factor A and B and all interactions based on trimmed means
# By default,
# tr=.2, meaning 20% trimming is used.
#
#  This function is the same as bbmpc, only it also reports a measures of effect size
#  based on explanatory power.
#
#
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
        JK <- J * K
        if(is.matrix(x))
                x <- listm(x)
        if(!is.na(grp[1])) {
                yy <- x
                x<-list()
                for(j in 1:length(grp))
                        x[[j]] <- yy[[grp[j]]]
        }
        if(!is.list(x))
                stop('Data must be stored in list mode or a matrix.')
        for(j in 1:JK) {
                xx <- x[[j]]
                x[[j]] <- xx[!is.na(xx)] # Remove missing values
        }
        #

        if(JK != length(x))
                warning('The number of groups does not match the number of contrast coefficients.')
for(j in 1:JK){
temp<-x[[j]]
temp<-temp[!is.na(temp)] # Remove missing values.
x[[j]]<-temp
}
        # Create the three contrast matrices
temp<-con2way(J,K)
conA<-temp$conA
conB<-temp$conB
conAB<-temp$conAB
if(!op){
Factor.A<-linconEP(x,con=conA,tr=tr,nreps=nreps,INT=FALSE,pr=FALSE,POOL=POOL)
Factor.B<-linconEP(x,con=conB,tr=tr,nreps=nreps,INT=FALSE,pr=FALSE,POOL=POOL)
Factor.AB<-linconEP(x,con=conAB,tr=tr,nreps=nreps,INT=TRUE,pr=FALSE,POOL=FALSE)
}
All.Tests<-NA
if(op){
Factor.A<-NA
Factor.B<-NA
Factor.AB<-NA
con<-cbind(conA,conB,conAB)
All.Tests<-lincon(x,con=con,tr=tr,alpha=alpha)
}
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.AB=Factor.AB,All.Tests=All.Tests,conA=conA,conB=conB,conAB=conAB)
}


bbmcpQS<-function(J,K,x,locfun,nreps=100,SEED=TRUE,POOL=TRUE,pr=TRUE){
#
# For independent groups,
#  compute quantile shift measure of effect size for all main effects and interactions.
#
#   To get an explanatory measure of effect size, use bbmcpEP
#
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
        JK <- J * K
        if(is.matrix(x) || is.data.frame(x))
                x <- listm(x)
        if(!is.list(x))
                stop('Data must be stored in list mode or a matrix.')
        if(JK != length(x))
                print('Warning: JK does not match the number of groups.')
x=elimna(x)  # Remove missing values.
DONE=FALSE
if(J==2 & K==2){
Factor.A=list()
Factor.A[[1]]=linconQS(x[1:2],pr=FALSE)
Factor.A[[2]]=linconQS(x[3:4],pr=FALSE)
Factor.B=list()
Factor.B[[1]]=linconQS(x[c(1,3)],pr=FALSE)
Factor.B[[2]]=linconQS(x[c(2,4)],pr=FALSE)
Factor.AB=linconQS(x,con=c(1,-1,-1,1),INT=TRUE,pr=FALSE)
DONE=TRUE
}
temp<-con2way(J,K)
conA<-temp$conA
conB<-temp$conB
conAB<-temp$conAB
if(!DONE){
        # Create the three contrast matrices
if(!POOL){  # For each level of Factor A, compute effect size
#          for all  pairwise comparisons among the  levels of B
ID=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
Factor.A=list()
for(j in 1:J){
id=as.vector(ID[j,])
Factor.A[[j]]=linconQS(x[id],pr=FALSE)
}
Factor.B=list()
ID=t(ID)
for(k in 1:K){
id=as.vector(ID[k,])
Factor.B[[k]]=linconQS(x[id],pr=FALSE)
}
}}
# Do interactions
Factor.AB=list()
for(l in 1:ncol(conAB)){
#id=which(conAB[,l]!=0)
Factor.AB[[l]]=linconQS(x,con=conAB[,l],INT=TRUE,pr=FALSE)
}
#
if(POOL){
ID=matrix(c(1:JK),nrow=J,ncol=K,byrow=TRUE)
LEV.A=list()
for(j in 1:J){
id=as.vector(ID[j,])
LEV.A[[j]]=pool.a.list(x[id])
}
Factor.A=linconQS(LEV.A,pr=FALSE)
ID=t(ID)
LEV.B=list()
for(k in 1:K){
id=as.vector(ID[k,])
LEV.B[[k]]=pool.a.list(x[id])
}
Factor.B=linconQS(LEV.B,pr=FALSE)
}
if(pr){
print('The columns of conAB contain the contrast coefficients for the interactions.')
print('For example, the output in FactorAB[[1]] are the results based')
print('on the contrast coefficients in column 1')
print('which is the interaction for the first two rows and the first two columns')
print('  ')
print('Note: Under normality and homoscedasticity, Cohen d= 0, .2, .5, .8')
print('correspond approximately to Q.effect = 0.5, 0.55, 0.65 and 0.70, respectively')
}
if(!POOL){
print('Factor.A: for each row of 1st factor, perform all pairwise')
print(' among the levels of Factor B and store the results in Factor.A')
print('Do the same for the second factor and store the results in Factor.B')
}
if(POOL){
print('Factor.A: for each row of 1st factor, pool the data over the levels')
print('of Factor B. Then do all pairwise comparisons and store the results in Factor.A')
print('Do the same for the second factor and store the results in Factor.B')
}
list(Factor.A=Factor.A,Factor.B=Factor.B,Factor.AB=Factor.AB,conAB=conAB)
}


bbdetmcp<-function(J,K,x,tr=0.2){
#
# For each level of Factor A, do all pairiwise comparisons
# among levels of B and store results in A in list mode.
#
# For each level of Factor B, do all pairiwise comparisons
# among levels of A and store results in B in list mode.
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
A=list()
B=list()
JK=J*K
mat<-matrix(c(1:JK),ncol=K,byrow=TRUE)
for(j in 1:J)A[[j]]=lincon(x[mat[j,]],tr=tr,pr=FALSE)
for(k in 1:K)B[[k]]=lincon(x[mat[,k]],tr=tr,pr=FALSE)
list(Levels.of.A=A,Level.of.B=B)
}


bbdetmcpQS<-function(J,K,x,tr=0.2){
#
# For each level of Factor A, do all pairiwise comparisons
# among levels of B and store results in A in list mode.
#
# For each level of Factor B, do all pairiwise comparisons
# amonglevels of A andstore results in B in list mode.
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
A=list()
B=list()
JK=J*K
mat<-matrix(c(1:JK),ncol=K,byrow=TRUE)
for(j in 1:J)A[[j]]=linconQS(x[mat[j,]],tr=tr,pr=FALSE)
for(k in 1:K)B[[k]]=linconQS(x[mat[,k]],tr=tr,pr=FALSE)
list(Levels.of.A=A,Level.of.B=B)
}


bmcpAKP<-function(x,tr=.2){
#
#  Compute quantile shift measure of effect size for pairs of J dependent groups
#
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
J=length(x)
C=(J^2-J)/2
A=matrix(NA,nrow=C,ncol=3)
dimnames(A)=list(NULL,c('Group','Group','Effect.Size'))
ic=0
for (j in 1:J){
for (k in 1:J){
if(j<k){
ic=ic+1
A[ic,1]=j
A[ic,2]=k
A[ic,3]=akp.effect(x[[j]],x[[k]],tr=tr)
}}}
A
}


bmcpQS<-function(x,locfun=median,...){
#
#  Compute quantile shift measure of effect size for all pairs of J independent groups
#
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
J=length(x)
C=(J^2-J)/2
A=matrix(NA,nrow=C,ncol=3)
dimnames(A)=list(NULL,c('Group','Group','Effect.Size'))
ic=0
for (j in 1:J){
for (k in 1:J){
if(j<k){
ic=ic+1
A[ic,1]=j
A[ic,2]=k
A[ic,3]=shiftes(x[[j]],x[[k]],locfun=locfun,...)$Q.Effect
}}}
A
}


wmcpAKP<-function(x,tr=.2){
#
#  Compute Algina et al measure of effect size for pairs of J dependent groups
#
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
J=length(x)
C=(J^2-J)/2
A=matrix(NA,nrow=C,ncol=3)
dimnames(A)=list(NULL,c('Group','Group','Effect.Size'))
ic=0
for (j in 1:J){
for (k in 1:J){
if(j<k){
ic=ic+1
A[ic,1]=j
A[ic,2]=k
A[ic,3]=D.akp.effect(x[[j]],x[[k]],tr=tr)
}}}
A
}


wmcpQS<-function(x,locfun=median){
#
#  Compute quantile shift measure of effect size for pairs of J dependent groups
#
if(is.matrix(x) || is.data.frame(x))x<-listm(x)
J=length(x)
C=(J^2-J)/2
A=matrix(NA,nrow=C,ncol=3)
dimnames(A)=list(NULL,c('Group','Group','Q.Effect'))
ic=0
for (j in 1:J){
for (k in 1:J){
if(j<k){
ic=ic+1
A[ic,1]=j
A[ic,2]=k
A[ic,3]=depQS(x[[j]],x[[k]],locfun=locfun)$Q.effect
}}}
A
}



################################################################################
# SPECIALIZED MCP FUNCTIONS
################################################################################

stepmcp<-function(x,tr=.2,alpha=.05){
#
# Step-down MCP method based on trimmed means
#
# x is assumed to have list mode, or a matrix or data with J columns
#   J=number of groups.
#
if(is.matrix(x) || is.data.frame(x))x=listm(x)
J=length(x)
if(J<3)stop('For two groups only, use yuen')
if(J>5)stop('Currently limited to at most five groups')
com=modgen(J)
jp1=J+1
mout=matrix(NA,nrow=length(com),ncol=3,
dimnames=list(NULL,c('Groups','p-value','p.crit')))
mout[,3]=alpha
jm2=J-2
com=com[jp1:length(com)]
mout=mout[jp1:nrow(mout),]
for(i in 1:length(com)){
nmod=length(com[[i]])-1
temp=c(nmod:0)
mout[i,1]=sum(com[[i]]*10^temp)
temp=t1way(x[com[[i]]],tr=tr)$p.value
pnum=length(com[[i]])
pe=1-(1-alpha)^(pnum/J)
if(length(com[[i]])<=jm2)mout[i,3]=pe
mout[i,2]=t1way(x[com[[i]]],tr=tr)$p.value
}
mout
}


signmcp<-function(x,y = NULL, alpha = 0.05, method='AC' , AUTO=TRUE,Method="hochberg"){
#
#  Dependent groups
#  Perform sign test for all pairwise differences
#
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
J<-ncol(x)
nval<-nrow(x)
ncon<-(J^2-J)/2
dvec<-alpha/c(1:ncon)
psihat<-matrix(NA,ncon,9)
dimnames(psihat)<-list(NULL,c("Group","Group","n","N","Prob_x_less_than_y","ci.lower","ci.upper",
"p.value","p.adjusted"))
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
temp=signt(x[,j],x[,k],alpha=alpha,method=method,AUTO=AUTO)
psihat[jcom,1]<-j
psihat[jcom,2]<-k
psihat[jcom,3]<-temp$n
psihat[jcom,4]<-temp$N
psihat[jcom,5]<-temp$Prob_x_less_than_y
psihat[jcom,6:7]=temp$ci
if(method!='SD')psihat[jcom,8]=temp$p.value
}}}
if(method!='SD')psihat[,9]=p.adjust(psihat[,8],method=Method)
list(output=psihat)
}


discmcp<-function(x,alpha=.05,nboot=500,SEED=TRUE,...){
#
#   Multiple comparisons for  J independent groups
#   having discrete distributions.
#   The method is based on a chi-squared test for each pair of groups to be compared
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   Missing values are allowed.
#
# Probability of one or more Type I errors controlled using Hochberg's method.
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in list mode or in matrix mode.')
J<-length(x)
ncon=(J^2-J)/2
Jm<-J-1
#
# Determine critical values
dvec=alpha/c(1:ncon)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
output<-matrix(NA,nrow=ncon,ncol=4)
dimnames(output)<-list(NULL,c('Group','Group','p.value','p.crit'))
ic=0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic=ic+1
output[ic,1]=j
output[ic,2]=k
output[ic,3]=disc2com(x[[j]],x[[k]],B=nboot)$p.value
}}}
temp2<-order(0-output[,3])
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,4]<-zvec
num.sig<-sum(output[,3]<=output[,4])
list(output=output,num.sig=num.sig)
}


sintmcp<-function(x, con=0, alpha=0.05){
#
#  Dependent groups
#  Multiple comparisons using medians on difference scores
#
flagcon=F
if(!is.matrix(x))x<-matl(x)
if(!is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
con<-as.matrix(con)
J<-ncol(x)
#xbar<-NULL
x<-elimna(x)  # Remove missing values
nval<-nrow(x)
if(sum(con^2!=0))CC<-ncol(con)
if(sum(con^2)==0)CC<-(J^2-J)/2
ncon<-CC
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
if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
if(sum(con^2)==0){
flagcon<-TRUE
psihat<-matrix(0,CC,7)
dimnames(psihat)<-list(NULL,c('Group','Group','psihat','ci.lower','ci.upper','p.value','p.crit'))
temp1<-0
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
dv<-x[,j]-x[,k]
temp=sintv2(dv,pr=FALSE)
temp1[jcom]<-temp$p.value
psihat[jcom,1]<-j
psihat[jcom,2]<-k
psihat[jcom,3]<-median(dv)
psihat[jcom,4]<-temp$ci.low
psihat[jcom,5]<-temp$ci.up
psihat[jcom,6]<-temp$p.value
}}}
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(psihat[temp2,6]>=zvec)
dd=0
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
}
psihat[temp2,7]<-zvec
}
if(sum(con^2)>0){
if(nrow(con)!=ncol(x))warning('The number of groups does not match the number
 of contrast coefficients.')
ncon<-ncol(con)
psihat<-matrix(0,ncol(con),6)
dimnames(psihat)<-list(NULL,c('con.num','psihat','ci.lower','ci.upper','p.value','p.crit'))
temp1<-NA
for (d in 1:ncol(con)){
psihat[d,1]<-d
for(j in 1:J){
if(j==1)dval<-con[j,d]*x[,j]
if(j>1)dval<-dval+con[j,d]*x[,j]
}
temp=sintv2(dval,pr=FALSE)
temp1[d]=temp$p.value
psihat[d,5]=temp$p.value
psihat[d,2]<-median(dval)
psihat[d,3]<-temp$ci.low
psihat[d,4]<-temp$ci.up
}
temp2<-order(0-temp1)
zvec<-dvec[1:ncon]
sigvec<-(psihat[temp2,5]>=zvec)
psihat[temp2,6]<-zvec
dd=0
if(sum(sigvec)<ncon){
dd<-ncon-sum(sigvec) #number that are sig.
ddd<-sum(sigvec)+1
}
}
list(output=psihat,con=con,num.sig=dd)
}


skmcp<-function(x,alpha=.05){
#
#   Multiple comparisons for  J independent groups
#   and binary data.
#   The method is based on the Storer--Kim
#   method for comparing independent binomials.
#
#   The data are assumed to be stored in x
#   which either has list mode or is a matrix.  In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, the columns of the matrix correspond
#   to groups.
#
#   Missing values are allowed.
#
# Probability of one or more Type I errors controlled using Hochberg's method.
#
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop('Data must be stored in list mode or in matrix mode.')
J<-length(x)
ncon=(J^2-J)/2
Jm<-J-1
#
# Determine critical values
dvec=alpha/c(1:ncon)
output<-matrix(NA,nrow=ncon,ncol=4)
dimnames(output)<-list(NULL,c('Group','Group','p.value','p.crit'))
ic=0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic=ic+1
output[ic,1]=j
output[ic,2]=k
output[ic,3]=twobinom(x=x[[j]],y=x[[k]])$p.value
}}}
temp2<-order(0-output[,3])
zvec<-dvec[1:ncon]
sigvec<-(test[temp2]>=zvec)
output[temp2,4]<-zvec
num.sig<-sum(output[,3]<=output[,4])
list(output=output,num.sig=num.sig)
}


mcp.nestAP<-function(x,tr=.2){
#
# Nested ANOVA
#
# Strategy: for each level of factor A, pool the data
# and then perform the analysis
#
# x is assumed to have list mode with length J,
# the number of independent groups.
#
# x[[1]] contains an n by K matrix, the nested data
# for the first level of the first factor.
# x[[2]] contains an n by K matrix, the nested data
# for the second level of the first factor, etc.
#
 xx=list()
for(j in 1: length(x))xx[[j]]=as.vector(x[[j]])
results=lincon(xx,tr=tr)
results
}


binmcp<-function(x,n,p.crit=NULL,alpha=.05,iter=2000,SEED=TRUE){
#
#
#  x is a vector containing the number of successes.
#  n is a vector indicating the sample sizes.
#
#  p.crit:  If NULL, critical p-values are determined so that that FWE is alpha
#  This is done using a simulation to determine the null distribution based on
#  iter=5000 replications.
#
#  All pairwise comparisons are done.
#
#   Confidence intervals having simultaneous probability coverage 1-alpha
#  using the adjusted level.
#
J=length(x)
A=(J^2-J)/2
if(J<2)stop('Should have 2 or more groups')
Jm1=J-1
est=x/n
pvec=NA
init=rep(alpha,Jm1)
if(is.null(p.crit)){
phat=sum(x)/sum(n)
pv.mat=binmcp.crit(phat,n=n,iter=iter,SEED=SEED)
}
p.crit=qest(pv.mat,alpha)
output<-matrix(NA,A,9)
dimnames(output)=list(NULL,c('Grp','Grp','Est 1','Est 2','Dif','ci.low','ci.up','p.value','p.crit'))
p.crit=p.crit/A
ic=0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic=ic+1
a=bi2KMSv2(x[j],n[j],x[k],n[k],alpha=p.crit)
output[ic,]=c(j,k,a$p1, a$p2,a$est.dif,a$ci[1],a$ci[2],a$p.value,p.crit)
}}}
list(n=n,output=output)
}


binmcp.crit<-function(p,n,iter=5000,SEED=TRUE){
#
#
#
if(SEED)set.seed(2)
J=length(n)  #Number of groups
A=(J^2-J)/2
pv.mat=matrix(NA,iter,A)
for(i in 1:iter){
x=rbinom(J,n,p)
pv.mat[i,]=binmcp.sub(x,n)
}
pv.mat
}



################################################################################
# P-VALUE ADJUSTMENT UTILITIES
################################################################################

mcpPV<-function(pv,alpha=.05,opt=1){
#
#  pv: A collection of p-values based on independent tests
#
#  Perform the step-down multiple comparison method in
#  Wilcox, R. R. \& Clark, F. (in press).
#   Robust multiple comparisons based on combined
#    probabilities from independent tests.  Journal of Data Science
#  based on K independent p-values
#
#  opt=1  Fisher's method
#  opt=2  Chen-Nadarajah method
#  opt=3  Max method
#
K=length(pv)
pick=NULL
v=order(pv)
ic=0
for(i in K:1){
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
list(p.values=pv,
Decisions=as.matrix(Decision),num.sig=nsig)
}


mcpKadjp <- function (p, k=1, proc = c('Holm'), rawp=p) {
#
#  MCP method based on results in
#
# Keselman, H. J., Miller, C. E., & Holland, B. (2011).
# Many tests of significance: New methods for controlling Type I errors.
# Psychological Methods, 16, 420-431.
#
# Also see
# Keselman, H. J., & Miller, C. E. (2012).
# Correction to many tests of significance:
# New methods for controlling Type I errors. Psychological Methods, 17(4), 679.
#
# p: The p-values to be adjusted.
# k: The  value for k-FWER
#  proc: indicates the method to be used. Choices are:
#' Holm'
# 'Hochberg',
#' 'RS', 	Romano-Shaikh procedure
 # 'Sarkar',
 # 'BH' , Benjamini--Hochberg
#
##  Generalized Hochberg is valid under MTP2 condition of the joint null
##     distribution of the p-values
##  Sarkar procedure is only valid for independent test statistics
#
D1 <- function(k=1, s=1000) {
#To calculate D1 values for Romano-Shaikh procedure
  alpha <- NULL
  for (i in 1:s) {
    if (i <= k) alpha[i]=k/s
    else alpha[i]=k/(s+k-i)
  }
  S <- NULL
  S[1:k] <- 0
  for (I in (k+1):s) {
    tmp <- NULL
    tmp[1:(k-1)] <- 0
    tmp[k]=I*alpha[s-I+k]/k
    for (j in (k+1):I) tmp[j]=I*(alpha[s-I+j]-alpha[s-I+j-1])/j
    S[I] <- sum(tmp)
    if (S[I] < S[I-1]) break

    maxI <- I
    maxS <- round(S[I],4)
  }
  return(list(S=S, maxI=maxI, maxS=maxS))
}
#modified the function mt.rawp2adjp from MTP package to k-FWER procedures
    m <- length(rawp)
    n <- length(proc)
    index <- order(rawp)
    spval <- rawp[index]
    adjp <- matrix(0, m, n + 1)
    dimnames(adjp) <- list(NULL, c('rawp', proc))
    adjp[, 1] <- spval

#################### Calculate adjusted p-values  ######################

#generalized Holm procedure based on Lehmann and Romano (2005)
    if (is.element('Holm', proc)) {
        crit <- sapply(c(rep(k,k-1),k:m), function(i)
                       k/(m-i+k))
        tmp <- 1/crit*spval
        tmp[tmp>1] <- 1
        for (i in 2:m) tmp[i] <- max(tmp[i-1], tmp[i])
        adjp[, 'Holm'] <- tmp
    }
#generalized Hochberg procedure (Step-up version of Lehmann and Romano)
    if (is.element('Hochberg', proc)) {
        crit <- sapply(c(rep(k,k-1),k:m), function(i)
                       k/(m-i+k))
        tmp <- 1/crit*spval
        tmp[tmp>1] <- 1
        for (i in (m-1):1) tmp[i] <- min(tmp[i+1], tmp[i])
        adjp[, 'Hochberg'] <- tmp
    }
#generalized Hochberg procedure based on Romano and Shaikh(2006)
    if (is.element('RS', proc)) {
        d <- D1(k,m)$maxS
        crit <- sapply(c(rep(k,k-1),k:m), function(i)
                       k/(m-i+k)/d)
        tmp <- 1/crit*spval
        tmp[tmp>1] <- 1
        for (i in (m-1):1) tmp[i] <- min(tmp[i+1], tmp[i])
        adjp[, 'RS'] <- tmp
    }
#generalized Hochberg procedure based on Sarkar(2008)
    ### Only for independent case  ###
    if (is.element('Sarkar', proc)) {
        crit <- sapply(c(rep(k,k-1),k:m), function(i)
                       (prod((1:k)/(m-i+(1:k)))))
        tmp <- 1/crit*(spval^k)
        tmp[tmp>1] <- 1
        for (i in (m-1):1) tmp[i] <- min(tmp[i+1], tmp[i])
       # Next line used to protect against possibility of adjp<rawp
        tmp <- pmax(tmp, spval)
        adjp[, 'Sarkar'] <- tmp
      }
#Benjamini and Hochberg (1995) procedure
        if (is.element('BH', proc)) {
        crit <- sapply(c(1:m), function(i)
                       i/m)
        tmp <- 1/crit*spval
        tmp[tmp>1] <- 1
        for (i in (m-1):1) tmp[i] <- min(tmp[i+1], tmp[i])
        adjp[, 'BH'] <- tmp
    }
### The following line returns original order of p-values
    adjp <- adjp[order(index),]
    return(adjp)
}


