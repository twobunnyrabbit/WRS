# WRS Package - Classification and Machine Learning Methods
# Extracted from Rallfun-v45.R
#
# This module contains classification and machine learning methods including:
# - K-means clustering (Kmeans, Kmeans.grp)
# - K-nearest neighbors (KNN, KNNv2, KNNdist)
# - Depth-based classification (Depth.class)
# - Logistic regression (class.logR)
# - Random forests (class.forest)
# - Gradient boosting (class.gbm)
# - Neural networks (NN.class)
# - Bagging (CLASS.BAG)
# - AdaBoost (class.ada)
# - General classifier (CLASS.fun)
# - Classification error estimation (class.error.com, class.error.loo, class.error.CM, class.error.CP, class.uni.error)
# - Ridge regression (ridge.test, ridge.est.k, ridge.Liu, ridge.Gtest, ridge.G.sub, ridgeGnull, ridgeGnullMC)
# - LASSO regression (lasso.est, lasso.rep)
#
# Total functions: 27
# Extraction date: 2025-12-30


# ============================================================================
# class.ada
# ============================================================================
class.ada<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,sm=FALSE,fr=2,SEED=NULL,baselearner='btree'){
#
# Do classification using adaboost
#
#   baselearner='btree': Stumps
#                         bbs: Splines
#                         bols:  linear models
#
# train is the training set
# test is the test data
# g contains labels for the data in the training set,
#
# This function removes the need to call library mboost.
#
#  SEED=NULL, used for convenience when called by other functions that expect SEED
#
library(mboost)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}
g=as.numeric(as.vector(g))
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
}
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(0,n1),rep(1,n2))
train=rbind(x1,x2)
if(is.null(test))stop('Argument test is null, contains  no data')
if(is.vector(test))stop('Argument test is a vector, should contain two or more variables')
g=as.factor(g)
ddata=data.frame(g,train)
d=mboost(g~., data=ddata,family=AdaExp(),baselearner=baselearner)
test=data.frame(test)
e=predict(d,newdata=test)
res=rep(1,nrow(test))
flag=e>0
res[flag]=2
res
}

# ============================================================================
# CLASS.BAG
# ============================================================================
CLASS.BAG<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,SEED=TRUE,kernel='radial',nboot=100,
method=c('KNN','SVM','DIS','DEP','PRO','NN','RF','ADA','LSM'),depthfun=prodepth,baselearner ='bbs',sm=TRUE,rule=.5,...){
#
# A collection of classification methods  for which the error rate is not
#  impacted by unequal sample sizes.
#  Bagged version of various classification methods is used.
#
#   For methods that do not require bagging, see UB.class.
#
#  KNN: calls KNNbag: a robust analog of the K nearest neighbor method
#  SVM:  a type of bagging method used in conjunction with support vector machine
#  DIS:  Points classified based on their depths
#  DEP: Uses depths as suggested by Makinde and Fasoranbaku (2018). JAS
#  PRO: project the points onto a line connecting the centers of the data clouds.
#       Then use estimate of the pdf for each group to make a decision about future points.
#  NN: Neural network
#  RF: Random forest
#  ADA: adaboost method
#  LSM: Uses a smoother designed for a binary dependent variable. sm=FALSE, uses logistic  regression
#
type=match.arg(method)
switch(type,
      KNN=KNNbag(train=train,test=test,g=g,x1=x1,x2=x2,SEED=SEED,nboot=nboot),
    SVM=SVMbag(train=train,test=test,g=g,x1=x1,x2=x2,depthfun=depthfun,SEED=SEED,nboot=nboot,...),
    DIS=dis.depth.bag(train=train,test=test,g=g,x1=x1,x2=x2,SEED=SEED,nboot=nboot),
    DEP=Depth.class.bag(train=train,test=test,g=g,x1=x1,x2=x2,SEED=SEED,nboot=nboot),
    PRO=pro.class.bag(train=train,test=test,g=g,x1=x1,x2=x2,SEED=SEED,nboot=nboot,...),
    NN=NNbag(train=train,test=test,g=g,x1=x1,x2=x2,SEED=SEED,...),
    RF=RFbag(train=train,test=test,g=g,x1=x1,x2=x2,SEED=SEED,kernel=kernel,nboot=nboot,...),
    ADA=class.ada.bag(train=train,test=test,g=g,x1=x1,x2=x2,baselearner=baselearner,SEED=SEED,nboot=nboot,...),
    LSM=LSMbag(train=train,test=test,g=g,x1=x1,x2=x2,sm=sm,rule=rule,SEED=SEED,nboot=nboot,...),
     )
}

# ============================================================================
# class.error.CM
# ============================================================================
class.error.CM<-function(x1=NULL,x2=NULL,train=NULL,g=NULL,method='KNN',nboot=100,EN=TRUE,FAST=TRUE,
AUC=FALSE,SEED=TRUE,...){
#
#  For a classification methods indicated by the argument
#  method
# use cross validation leaving one out.
#
#. Return a confusion matrix: unconditional, To get a conditional result use class.error.CP
#
#  The data for the two groups can be entered via  the arguments
#  x1 and x2
#  or
#  store all of the data in the argument train in which case g specifies the group
# AUC=TRUE, returns auc. Default is FALSE because conditions can be created where
#    Error: $ operator is invalid for atomic vectors
#
# Current choices available:
#  KNN: Nearest neighbor using robust depths
#  DIS:  Points classified based on their depths
#  DEP: Uses depths as suggested by Makinde and Fasoranbaku (2018). JAS
#  SVM: support vector machine
#  RF: Random forest
#  NN: neural network
#  ADA: ada boost
#  PRO: project the points onto a line connecting the centers of the data clouds.
#       Then use estimate of the pdf for each group to make a decision about future points.
#  LSM: smooth version of logistic regression when sm=TRUE; otherwise use logistic regression.
#
#  Returns confusion matrix
#
#
# method='KNN'     is default
#
#  nboot=number of samples
#
if(length(method)!=1)stop('Only one method at a time is allowed')
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g))if(dim(g)>1)stop('Argument g should be a vector')
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
}
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
if(is.null(x1))stop('Something is wrong, no data in x1')
if(is.null(x2))stop('Something is wrong, no data in x2')
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
dimnames(x1)=list(NULL,NULL) # can be necessary to eliminate labels  to avoid an error in randomForest.
dimnames(x2)=list(NULL,NULL)
n1=nrow(x1)
n2=nrow(x2)
ns1=min(n1,nboot)
ns2=min(n2,nboot)
mn=min(ns1,ns2)
CM=matrix(0,2,2)
isub1=sample(c(1:ns1))
isub2=sample(c(1:ns2))
A1=NULL
A2=NULL
ic1=0
ic2=0
for(k in 1:mn){
N1=isub1[k]
N2=isub2[k]
train1=x1[-N1,]
train2=x2[-N2,]
test=rbind(x1[N1,],x2[N2,])
a=CLASS.fun(x1=train1,x2=train2,test=test,method=method,...)
a1=a[1]
a2=a[2]
A1[k]=a1
A2[k]=a2
if(a1==1)CM[1,1]=CM[1,1]+1 #true = 1 pred 1
else
CM[1,2]=CM[1,2]+1
if(a2==2)CM[2,2]=CM[2,2]+1    #true =2 and pred 2
else
CM[2,1]=CM[2,1]+1
}
FREQ=CM
CM=CM/(2*nboot)
F=matrix(NA,3,3)
dimnames(F)=list(c('True 1','True 2','Sum'),c('Pred 1','Pred 2','Sum'))
F[1,1]=FREQ[1,1]
F[1,2]=FREQ[1,2]
F[2,1]=FREQ[2,1]
F[2,2]=FREQ[2,2]
F[1,3]=F[1,1]+F[1,2]
F[2,3]=F[2,1]+F[2,2]
F[3,1]=F[1,1]+F[2,1]
F[3,2]=F[1,2]+F[2,2]
F[3,3]=F[1,3]+F[2,3]
RES=F/F[3,3]
auroc=NULL
if(AUC){
library(ROCR)
PRED=c(A1,A2)
LABS=c(rep(1,length(A1)),rep(2,length(A2)))
pred=prediction(PRED,LABS)
perf=performance(pred, "auc")
 auroc<- perf@y.values[[1]]
}
dimnames(RES)=list(c('True 1','True 2','Sum'),c('Pred 1','Pred 2','Sum'))
list(C.MAT=RES,COUNTS=F,AUC=auroc[[1]])
}

# ============================================================================
# class.error.com
# ============================================================================
class.error.com<-function(x1=NULL,x2=NULL,train=NULL,g=NULL,method=NULL,
pro.p=.8,nboot=100,EN=FALSE,FAST=TRUE,
SEED=TRUE,...){
#
#  For two classification methods indicated by the arguments
#  class.fun1  and
#  class.fun2
# use cross validation coupled with resampling  to estimate the probability that of a correct classification.
#
#  The data for the two groups can be entered via  the arguments
#  x1 and x2
#  or
#  store all of the data in the argument train in which case g specifies the group
#
#  Choices for these two arguments:
#
# Current choices available:
#  KNN: Nearest neighbor using robust depths
#  DIS:  Points classified based on their depths
#  DEP: Uses depths as suggested by Makinde and Fasoranbaku (2018). JAS
#  SVM: support vector machine
#  RF: Random forest
#  NN: neural network
#  ADA: ada boost
#  PRO: project the points onto a line connecting the centers of the data clouds.
#       Then use estimate of the pdf for each group to make a decision about future points.
#  LSM: smooth version of logistic regression when sm=TRUE; otherwise use logistic regression.
#
#   EN=TRUE; use equal samples for the test data to deal with classification bias
#      Otherwise, the  ratio of the sample sizes is n1/n2
#
# method=NULL  All of the methods listed above will be compared if
#  FAST=FALSE
#   For method 'PRO', execution time might take several minutes if the sample sizes are large
#  For this reason.PRO is is not used if FAST=TRUE
#
#  pro.p=.8  means 80% of the data will be used as training data
#  nboot=number of bootstrap samples

#  Returns estimate of the error rate plus
#  FP (false positive): average proportion of values in x1 erroneously classified as coming from x2
#   Example, x1 contains no fracture, x2 contains fractures.
#  FN (false negative): average proportion of values in x2 erroneously classified as coming from x1
#
#
#
if(is.null(method))method=c('KNN','DIS','DEP','SVM','RF','NN','ADA','PRO','LSM')
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g))if(dim(g)>1)stop('Argument g should be a vector')
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
}
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
if(is.null(x1))stop('Something is wrong, no data in x1')
if(is.null(x2))stop('Something is wrong, no data in x2')
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
dimnames(x1)=list(NULL,NULL) # can be necessary to eliminate labels  to avoid an error in randomForest.
dimnames(x2)=list(NULL,NULL)
n1=nrow(x1)
n2=nrow(x2)
ns1=round(pro.p*n1)
ns2=round(pro.p*n2)
if(EN)ns1=ns2=min(c(ns1,ns2))
P1hat=NA
P2hat=NA
Av=NA
Bv=NA
Cv=NA
Dv=NA

J=length(method)
TE=matrix(NA,nrow=nboot,ncol=J)
FP=matrix(NA,nrow=nboot,ncol=J)
FN=matrix(NA,nrow=nboot,ncol=J)

for(k in 1:nboot){
N1=sample(n1,ns1)
N2=sample(n2,ns2)
test1=x1[-N1,]
test2=x2[-N2,]
for(j in 1:J){
a1=CLASS.fun(x1=x1[N1,],x2=x2[N2,],test=test1,method=method[j],...)
a2=CLASS.fun(x1=x1[N1,],x2=x2[N2,],test=test2,method=method[j],...)
flag1=a1!=1 # ID  False negatives e..g., method 1 predict no fracture but fracture occurred. So !flag1 is correct decision
flag2=a2!=2 # ID  False positives e..g., predict fracture but no fracture occurred.
flag=c(flag1,flag2)   #Overall mistakes
TE[k,j]=mean(flag)
FN[k,j]=mean(flag1)
FP[k,j]=mean(flag2)
}}

ERR=matrix(NA,nrow=3,ncol=J)
dimnames(ERR)=list(c('TE','FP','FN'),method)
v=apply(TE,2,mean)
ERR[1,]=v
v=apply(FP,2,mean)
ERR[2,]=v
v=apply(FN,2,mean)
ERR[3,]=v

list(Error.rates=ERR)
}

# ============================================================================
# class.error.CP
# ============================================================================
class.error.CP<-function(x1=NULL,x2=NULL,train=NULL,g=NULL,method='KNN',nboot=100,EN=TRUE,FAST=TRUE,
AUC=FALSE,SEED=TRUE,...){
#
#. Requires ROCR pacakage
#
#  For a classification methods indicated by the argument
#  method
# use cross validation leaving one out.
#
#. Return a confusion matrix
#
#
#  The data for the two groups can be entered via  the arguments
#  x1 and x2
#  or
#  store all of the data in the argument train in which case g specifies the group
#
# AUC=TRUE, returns auc. Default is FALSE because conditions can be created where
#    Error: $ operator is invalid for atomic vectors
#
# Current choices available:
#  KNN: Nearest neighbor using robust depths
#  DIS:  Points classified based on their depths
#  DEP: Uses depths as suggested by Makinde and Fasoranbaku (2018). JAS
#  SVM: support vector machine
#  RF: Random forest
#  NN: neural network
#  ADA: ada boost
#  PRO: project the points onto a line connecting the centers of the data clouds.
#       Then use estimate of the pdf for each group to make a decision about future points.
#  LSM: smooth version of logistic regression when sm=TRUE; otherwise use logistic regression.
#
#  Returns confusion matrix
#
#
# method='KNN'     is default
#
#  nboot=number of samples
#
if(length(method)!=1)stop('Only one method at a time is allowed')
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g))if(dim(g)>1)stop('Argument g should be a vector')
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
}
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
if(is.null(x1))stop('Something is wrong, no data in x1')
if(is.null(x2))stop('Something is wrong, no data in x2')
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
dimnames(x1)=list(NULL,NULL) # can be necessary to eliminate labels  to avoid an error in randomForest.
dimnames(x2)=list(NULL,NULL)
n1=nrow(x1)
n2=nrow(x2)
ns1=min(n1,nboot)
ns2=min(n2,nboot)
mn=min(ns1,ns2)
CM=matrix(0,2,2)
isub1=sample(c(1:ns1))
isub2=sample(c(1:ns2))
A1=NULL
A2=NULL
ic1=0
ic2=0
for(k in 1:mn){
N1=isub1[k]
N2=isub2[k]
train1=x1[-N1,]
train2=x2[-N2,]
test=rbind(x1[N1,],x2[N2,])
a=CLASS.fun(x1=train1,x2=train2,test=test,method=method,...)
a1=a[1]
a2=a[2]
A1[k]=a1
A2[k]=a2
if(a1==1)CM[1,1]=CM[1,1]+1 #true = 1 pred 1
else
CM[1,2]=CM[1,2]+1
if(a2==2)CM[2,2]=CM[2,2]+1    #true =2 and pred 2
else
CM[2,1]=CM[2,1]+1
}
FREQ=CM
CM=CM/(2*nboot)
F=matrix(NA,2,3)
dimnames(F)=list(c('True 1','True 2'),c('Pred 1','Pred 2','Sum'))
F[1,1]=FREQ[1,1]
F[1,2]=FREQ[1,2]
F[2,1]=FREQ[2,1]
F[2,2]=FREQ[2,2]
F[1,3]=F[1,1]+F[1,2]
F[2,3]=F[2,1]+F[2,2]
RES=F
RES[1,1]=F[1,1]/(F[1,1]+F[1,2])
RES[1,2]=F[1,2]/(F[1,1]+F[1,2])
RES[2,1]=F[2,1]/(F[2,1]+F[2,2])
RES[2,2]=F[2,2]/(F[2,1]+F[2,2])
RES[,3]=1
auroc=NULL
if(AUC){
library(ROCR)
PRED=c(A1,A2)
LABS=c(rep(1,length(A1)),rep(2,length(A2)))
pred=prediction(PRED,LABS)
perf=performance(pred, "auc")
 auroc<- perf@y.values[[1]]
 }
list(C.MAT=RES,COUNTS=F,AUC=auroc[[1]])
}

# ============================================================================
# class.error.loo
# ============================================================================
class.error.loo<-function(x1=NULL,x2=NULL,train=NULL,g=NULL,method=NULL,nboot=100,EN=TRUE,FAST=TRUE,
SEED=TRUE,...){
#
#  For one or more classification methods indicated by the argument
#  method
# use cross validation leaving one out.
#
#  The data for the two groups can be entered via  the arguments
#  x1 and x2
#  or
#  store all of the data in the argument train in which case g specifies the group
#
# Current choices available:
#  KNN: Nearest neighbor using robust depths
#  DIS:  Points classified based on their depths
#  DEP: Uses depths as suggested by Makinde and Fasoranbaku (2018). JAS
#  SVM: support vector machine
#  RF: Random forest
#  NN: neural network
#  ADA: ada boost
#  PRO: project the points onto a line connecting the centers of the data clouds.
#       Then use estimate of the pdf for each group to make a decision about future points.
#  LSM: smooth version of logistic regression when sm=TRUE; otherwise use logistic regression.
#
#  Returns estimate of the error rates plus
#  FP: average proportion of values in x1 erroneously classified as coming from x2
#   Example, x1 contains no fracture, x2 contains fractures.
# FN: average proportion of values in x2 erroneously classified as coming from x1
#
#   EN=TRUE; use equal samples for the test data to deal with classification bias
#      Otherwise, the  ratio of the sample sizes is n1/n2
#
# method=NULL  All of the methods listed above will be compared if
#  FAST=FALSE
#   For method 'PRO', execution time might take several minutes if the sample sizes are large
#  For this reason.PRO is is not used if FAST=TRUE
#
#  pro.p=.8  means 80% of the data will be used as training data
#  nboot=number of bootstrap samples
#
if(is.null(method)){
method=c('KNN','DIS','DEP','SVM','RF','NN','ADA','PRO','LSM','GBT')
if(FAST)method=c('KNN','DIS','DEP','SVM','RF','NN','ADA','LSM')
}
if(SEED)set.seed(2)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g))if(dim(g)>1)stop('Argument g should be a vector')
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
}
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
if(is.null(x1))stop('Something is wrong, no data in x1')
if(is.null(x2))stop('Something is wrong, no data in x2')
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
dimnames(x1)=list(NULL,NULL) # can be necessary to eliminate labels  to avoid an error in randomForest.
dimnames(x2)=list(NULL,NULL)
n1=nrow(x1)
n2=nrow(x2)
ns1=min(n1,nboot)
ns2=min(n2,nboot)
nm=min(ns1,ns2)
if(EN)ns1=ns2=min(c(ns1,ns2))
P1hat=NA
P2hat=NA
Av=NA
Bv=NA
Cv=NA
Dv=NA

J=length(method)
TE=matrix(NA,nrow=nm,ncol=J)
FP=matrix(NA,nrow=nm,ncol=J)
FN=matrix(NA,nrow=nm,ncol=J)
TP=matrix(NA,nrow=nm,ncol=J)
TN=matrix(NA,nrow=nm,ncol=J)

isub1=sample(c(1:nm))
isub2=sample(c(1:nm))
for(k in 1:nm){
N1=isub1[k]
N2=isub2[k]
train1=x1[-N1,]
train2=x2[-N2,]
test=rbind(x1[N1,],x2[N2,])
for(j in 1:J){
a=CLASS.fun(x1=train1,x2=train2,test=test,method=method[j],...)
a1=a[1]
a2=a[2]
flag1=a1!=1 # ID  False negatives e..g., method 1 predict no fracture but fracture occurred. So !flag1 is correct decision
flag2=a2!=2 # ID  False positives e..g., predict fracture but no fracture occurred.
flag=c(flag1,flag2)   #Overall mistakes
TE[k,j]=mean(flag)
FN[k,j]=mean(flag1)
FP[k,j]=mean(flag2)
flag3=a1==1
flag4=a2==2
TP[k,j]=mean(flag3)  #method 1 predict  fracture and fracture occurred
TN[k,j]=mean(flag4) #method 1 predict no fracture and  no fracture occurred
}}
ERR=matrix(NA,nrow=5,ncol=J)
dimnames(ERR)=list(c('TE','FP','FN','TP','TN'),method)
#dimnames(CAT)=list(c('TRUE 1','TRUE 2'),c('PRED 1','PRED2'))
v=apply(TE,2,mean)
ERR[1,]=v
v=apply(FP,2,mean)
ERR[2,]=v
v=apply(FN,2,mean)
ERR[3,]=v
v=apply(TP,2,mean)
ERR[4,]=v
v=apply(TN,2,mean)
ERR[5,]=v
list(Error.rates=ERR)
}

# ============================================================================
# class.forest
# ============================================================================
class.forest<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,sm=FALSE,fr=2,SEED=NULL){
#
# Do classification using random forest
#
#
# train is the training set
# test is the test data
# g contains labels for the data in the training set,
#
#  Alternatively, store the data for the two groups in
#  x1 and x2, in which case the function creates labels, i.e., no need to specify train and g.
#
#  SEED=NULL, used for convenience when called by other functions that expect SEED
#
if(is.null(test))stop('Argument test is null, contains  no data')
library(randomForest)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}
g=as.numeric(as.vector(g))
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
}
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(0,n1),rep(1,n2))
train=rbind(x1,x2)
if(is.vector(test))stop('Argument test is a vector, should contain two or more variables')
g=as.factor(g)
train=as.data.frame(train)
d=randomForest(g~., train)
test=as.data.frame(test)
e=predict(d,newdata=test)
res=as.numeric(as.vector(e))+1
res
}

# ============================================================================
# CLASS.fun
# ============================================================================
CLASS.fun<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,
method=c('KNN','DIS','DEP','SVM','RF','NN','ADA','PRO','LSM','GBT'),
depthfun=prodepth,kernel='radial',baselearner='btree',sm=TRUE,rule=.5,...){
#
# A collection of classification methods:
#
#  KNN: Nearest neighbor using robust depths
#  DIS:  Points classified based on their depths
#  DEP: Uses depths as suggested by Makinde and Fasoranbaku (2018). JAS
#  SVM: support vector machine
#  RF: Random forest
#  NN: neural network
#  ADA: ada boost
#  PRO: project the points onto a line connecting the centers of the data clouds.
#       Then use estimate of the pdf for each group to make a decision about future points.
# GBT: Gradient boosted trees: requires R package gbm
#
type=match.arg(method)
switch(type,
    KNN=KNNdist(train=train,test=test,g=g,x1=x1,x2=x2,depthfun=depthfun),
    DIS=discdepth(train=train,test=test,g=g,x1=x1,x2=x2),
    DEP=Depth.class(train=train,test=test,g=g,x1=x1,x2=x2),
    SVM=SVM(train=train,test=test,g=g,x1=x1,x2=x2,kernel=kernel),
    RF=class.forest(train=train,test=test,g=g,x1=x1,x2=x2),
    NN=NN.class(train=train,test=test,g=g,x1=x1,x2=x2),
    ADA=class.ada(train=train,test=test,g=g,x1=x1,x2=x2,baselearner=baselearner),
    PRO=pro.classPD(train=train,test=test,g=g,x1=x1,x2=x2),
    LSM=class.logR(train=train,test=test,g=g,x1=x1,x2=x2,sm=sm,rule=rule),
    GBT=class.gbm(train=train,test=test,g=g,x1=x1,x2=x2),
    )
}

# ============================================================================
# class.gbm
# ============================================================================
class.gbm<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,SEED=NULL,n.trees=100){
#
# Do classification using boosting via the R package gdm
#
#
# train is the training set
# test is the test data
# g contains labels for the data in the training set,
#
#  Alternatively, store the data for the two groups in
#  x1 and x2, in which case the function creates labels, i.e., no need to specify train and g.
#
# This function removes the need to call library gbm.
#
#  SEED=NULL, used for convenience when called by other functions that expect SEED
#
if(is.null(test))stop('Argument test is null, contains  no data')
library(gbm)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
}

x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(0,n1),rep(1,n2))
train=rbind(x1,x2)
if(is.vector(test))stop('Argument test is a vector, should contain two or more variables')
if(ncol(train)!=ncol(test))stop('training data and test data have different number of columns')
data=data.frame(g,train)
a=gbm(g~., data=data,distribution ='bernoulli',n.trees=n.trees)
test=data.frame(test)
e=predict(a,newdata=test,n.trees=n.trees)
res=rep(1,nrow(test))
flag=e>0
res[flag]=2
res
}

# ============================================================================
# class.logR
# ============================================================================
class.logR<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,sm=TRUE,fr=2,rule=.5,SEED=NULL){
#
# Do classification using logistic or a smoother
#  sm=TRUE: a smoother will be used with the span taken to be
#  fr is this span
#  sm=FALSE: use logistic regression.
#
#
# train is the training set
# test is the test data
# g contains labels for the data in the training set,
#
# This function removes the need to call library class.
# For more information, use the command ?knn
#
#  SEED=NULL, used for convenience when called by other functions that expect SEED
#
if(!is.null(train)){
train=as.matrix(train)
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g))if(dim(g)>1)stop('Argument g should be a vector')
flag=g==min(g)
x1=train[flag,]
x2=train[!flag,]
}
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
if(is.null(x1))stop('Something is wrong, no data in x1')
if(is.null(x2))stop('Something is wrong, no data in x2')
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
n1=nrow(x1)
n2=nrow(x2)
train=rbind(x1,x2)
g=c(rep(1,n1),rep(2,n2))
if(!sm)e=logreg.pred(train,g,test)
if(sm)e=logSMpred(train,g,test,fr=fr)
if(is.null(test))stop('Argument test is null, contains  no data')
test=as.matrix(test)
res=rep(1,nrow(test))
flag=e>rule
res[flag]=2
res
}

# ============================================================================
# class.uni.error
# ============================================================================
class.uni.error<-function(x,y,xy=NULL,g=NULL){
#
#  For  univariate data, estimate prediction error using
# a kernel density estimator.
#  Returns conditional   estimate of the error rates
#  Example: Given that a values is from group 1  erroneously classify it as coming from
#  group 2
#
# Also returns unconditional probabilities
# Example, the probability that a randomly sample subject is in group 1 and
# is classified as being in group 1
#
#
if(!is.null(xy)){
xy=fac2list(xy,g)
x=xy[[1]]
y=xy[[2]]
}
x=elimna(x)
y=elimna(y)
n1=length(x)
n2=length(y)
xsort=sort(x)
ysort=sort(y)
n1p=n1+1
N=n1+n2
if(is.null(x1))stop('Something is wrong, no data in x1')
if(is.null(x2))stop('Something is wrong, no data in x2')
UC=matrix(NA,2,2)
CP=matrix(NA,2,2)
d1<-akerd(x,pts=xsort,pyhat=TRUE,plotit=FALSE)
d2<-akerd(y,pts=xsort,pyhat=TRUE,plotit=FALSE)
e1=d1>d2 #means predict group 1 for data in group 1
D1<-akerd(x,pts=ysort,pyhat=TRUE,plotit=FALSE)
D2<-akerd(y,pts=ysort,pyhat=TRUE,plotit=FALSE)
e2=D2>D1 # means predict group 2 for data in group 2
CP[1,1]=mean(e1==1)
CP[1,2]=1-CP[1,1]
CP[2,2]=mean(e2==0)
CP[2,1]=1-CP[2,2]
#
e3=c(e1,e2)
UC[1,1]=sum(e3[1:n1]==1)
UC[1,2]=sum(e3[1:n1]==0)
UC[2,1]=sum(e3[n1p:N]==0)
UC[2,2]=sum(e3[n1p:N]==1)
UC=UC/sum(UC)
dimnames(CP)=list(c('True 1','True 2'),c('Pred 1','Pred 2'))
dimnames(UC)=list(c('True 1','True 2'),c('Pred 1','Pred 2'))
list(Conditional.prob=CP, Unconditiion.prob=UC,prob.correct.decision=UC[1,1]+UC[2,2])
}

# ============================================================================
# Depth.class
# ============================================================================
Depth.class<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,depthfun=prodepth,DIST=TRUE,SEED=TRUE,...){
#
# Do classification using max depths or max depth distribution as suggested by
# Makinde and Fasoranbaku (2018). JAS
#
#  depthfun indicates how the depth of a point is computed.
#  By default, projection depth is used. depthfun=zonoid would use zonoid depth
#
# train is the training set
# test is the test data
# g:  labels for the data in the training set.
#
#  depthfun must be a function having the form depthfun(x,pts).
#  That is, compute depth for the points in pts relative to points in x.
#
#
library(class)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}
Train=cbind(train,g)
Train=elimna(Train)
p=ncol(train)
p1=p+1
train=Train[,1:p]
g=Train[,p1]
flag=g==min(g)
x1=Train[flag,1:p]
x2=Train[!flag,1:p]
}
x1=elimna(x1)
x2=elimna(x2)
x1=as.matrix(x1)
x2=as.matrix(x2)
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(0,n1),rep(1,n2))
train=rbind(x1,x2)
train=elimna(train)
test=elimna(test)
train=as.matrix(train)
test=as.matrix(test)
if(ncol(test)==1)test=t(test)  # If test is a vector, a single point, transpose to get correct number of columns.
if(ncol(test)!=ncol(train))stop('test and train do not have the same number of columns')
ntest=nrow(test)
P=ncol(train)
P1=P+1
xall=as.data.frame(matrix(NA,nrow=nrow(train),ncol=P1))
xall[,1:P]=train
xall[,P1]=g
xall=elimna(xall)
x1=xall[,1:P]
xall=as.matrix(xall)
g=as.vector(xall[,P1])
ids=unique(g)  # Number of categories
x2=elimna(test)
x1=as.matrix(x1)
x2=as.matrix(x2)
n=nrow(x1)
n2=nrow(x2)
p=length(ids)
D=matrix(NA,nrow=p,ncol=n2)
for(i in 1:p){
flag=g==ids[i]
D[i,]=depthfun(as.matrix(x1[flag,]),pts=x2,SEED=SEED,...) # depth of test points relative to train data in cat i
}
if(!DIST)res=apply(D,2,which.max)
if(DIST){
res=NA
all.dep=list()
for(i in 1:p){
flag=g==ids[i]
all.dep[[i]]=depthfun(x1[flag,],pts=x1[flag,],...) #Have depth of all point in the training set for class i
}
for(j in 1:ntest){
dt=NA
cum.depth=NA
for(i in 1:p){
flag=g==ids[i]
chkpt=matrix(x2[j,],nrow=1)
#dt[i]=depthfun(x1[flag,],pts=x2[j,],...)
dt[i]=depthfun(x1[flag,],chkpt,...)
cum.depth[i]=mean(all.dep[[i]]<=dt[i])
}
chkit=which(cum.depth==max(cum.depth))
res[j]=chkit[1]
}}
res
}

# ============================================================================
# Kmeans
# ============================================================================
Kmeans<-function(x,k,xout=FALSE,outfun=outpro){
#
# Do K means cluster analysis, outliers removed
#  if xout=T
x=elimna(x)
x=as.matrix(x)
if(xout){
flag<-outfun(x,plotit=FALSE)$keep
x<-x[flag,]
}
x=as.matrix(x)
res=kmeans(x,k)
res
}

# ============================================================================
# Kmeans.grp
# ============================================================================
Kmeans.grp<-function(x,k,y,xout=FALSE,outfun=outpro){
#
#
# x and y are assumed to be matrices, each having n rows.
# This function create k groups based on data in x
#  Then it sorts the data in y into k groups based as indicated
#  by k-means cluter method applied to x
#
x=elimna(x)
x=as.matrix(x)
y=as.matrix(y)
if(xout){
flag<-outfun(x,plotit=FALSE)$keep
x<-x[flag,]
y=y[flag,]
}
y=as.matrix(y)
res=kmeans(x,k)
grpid=res[[1]]
grpdat=list()
for(i in 1:k){
flag=(grpid==i)
grpdat[[i]]=y[flag,]
}
grpdat
}

# ============================================================================
# KNN
# ============================================================================
KNN<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,k=3,prob=TRUE,SEED=NULL){
#
# Do classification using the kNN method
#
# k: number of nearest neighbors
#
# train is the training set
# test is the test data
# g contains labels for the data in the training set,
# If data for the two groups are stored in
#  x1
# and
# x2,
#  the function creates labels for you.
#
# This function removes the need to call library class.
# For more information, use the command ?knn
#
#  SEED=NULL, used for convenience when called by other functions that expect SEED
#
library(class)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}}
if(!is.null(x1)){
if(!is.null(x2)){
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(0,n1),rep(1,n2))
train=rbind(x1,x2)
}}
res=knn(train,test,cl=as.factor(g),k=k,prob=prob)
res
}

# ============================================================================
# KNNdist
# ============================================================================
KNNdist<-function(train=NULL,test=NULL,g,k=3,x1=NULL,x2=NULL,prob=FALSE,plotit=FALSE,SEED=NULL,
xlab='Group 1',ylab='Group 2',depthfun=prodepth,...){
#
# Do classification using depths and the kNN method.
# Points are transformed to their depth in each group and knn is applied
# using the resulting depth values.
# See Li et al., 2012, DD-classifier: nonparametric classification
#  procedure based on DD-plot. Journal of the American Statistical Association,
#  107, 737--753
#
#  depthfun indicates how the depth of a point is computed.
#  By default, projection distances are used.
#
# train is the training set
# test is the test data
# g:  labels for the data in the training set.
#
#  depthfun must be a function having the form depthfun(x,pts).
#  That is, compute depth for the points in pts relative to points in x.
#
#  SEED: not used here,  included for convenience when this function is called by other functions
#
library(class)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}}
if(!is.null(x1)){
if(!is.null(x2)){
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(0,n1),rep(1,n2))
train=rbind(x1,x2)
}}
g=as.numeric(as.vector(g))
train=elimna(train)
if(is.null(test))stop('Argument test is null,  no data')
test=elimna(test)
train=as.matrix(train)
test=as.matrix(test)
if(ncol(train)!=ncol(test))stop('The first two arguments, train and test, should have the same number of coluimns')
P=ncol(train)
P1=P+1
xall=as.data.frame(matrix(NA,nrow=nrow(train),ncol=P1))
xall[,1:P]=train
xall[,P1]=g
xall=elimna(xall)
x1=xall[,1:P]
xall=as.matrix(xall)
x1=as.matrix(x1)
g=as.vector(xall[,P1])
ids=unique(g)
x2=elimna(test)
x1=as.matrix(x1)
x2=as.matrix(x2)
n=nrow(x1)
n2=nrow(x2)
p=length(ids)
d=matrix(NA,nrow=n,ncol=p)
D=matrix(NA,nrow=n2,ncol=p)
for(i in 1:length(ids)){
flag=g==ids[i]
d[,i]=depthfun(as.matrix(x1[flag,]),pts=x1,...)
D[,i]=depthfun(as.matrix(x1[flag,]),pts=x2,...)
}
res=NULL
res=knn(d,D,cl=as.factor(g),k=k,prob=prob)
if(plotit){
if(p==2){
plot(d[,1],d[,2],xlab=xlab,ylab=ylab,type='n')
flag=g==ids[1]
points(d[flag,1],d[flag,2])
points(d[!flag,1],d[!flag,2],pch='*')
}}
res=as.numeric(res)
res
}

# ============================================================================
# KNNv2
# ============================================================================
KNNv2<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,k=3,prob=TRUE,SEED=NULL){
#
# Do classification using the kNN method
#
#
# train is the training set
# test is the test data
# g contains labels for the data in the training set,
#
# This function removes the need to call library class.
# For more information, use the command ?knn
#
#  SEED=NULL, used for convenience when called by other functions that expect SEED
#
library(class)
if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}}
if(!is.null(x1)){
if(!is.null(x2)){
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(0,n1),rep(1,n2))
train=rbind(x1,x2)
}}
res=knn(train,test,cl=as.factor(g),k=k,prob=prob)
res=as.numeric(as.vector(res))+1
res
}

# ============================================================================
# lasso.est
# ============================================================================
lasso.est<-function(x,y,xout=FALSE,STAND=TRUE,outfun=outpro,regout=FALSE,lam=NULL,...){
#
#  Lasso regression via the R package glmnet.
#  This function includes the option of eliminating leverage points
#  This function is for convenience, it returns the estimates of the
#  coefficients only. The R function cv.glmnet provides more complete details
#  and includes other options.
#
#   xout=TRUE eliminate leverage points with the function
#   outfun
#   regout=TRUE eliminate regression outliers with the function elo
#
library(glmnet)
x<-as.matrix(x)
p1<-ncol(x)+1
p<-ncol(x)
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
nrem=length(y)
if(regout){
flag=elo(x,y,outfun=outfun,lev=TRUE,reg=xout)$keep
xy<-xy[flag,]
x<-xy[,1:p]
y<-xy[,p1]
xout=FALSE
}
if(xout){
flag<-outfun(x,plotit=FALSE,...)$keep
xy<-xy[flag,]
x<-xy[,1:p]
y<-xy[,p1]
}
if(STAND)x=standm(x)
z=cv.glmnet(x, y, family='gaussian',lam=lam)
e=coef(z,s=z$lambda.min)
e=as.vector(e)
list(coef=e,lambda.min.used=z$lambda.min)
}

# ============================================================================
# lasso.rep
# ============================================================================
lasso.rep<-function(x,y,lasso.fun=lasso.est,xout=FALSE,outfun=outpro,iter=20,EST=FALSE,SEED=TRUE){
#
#  R function lasso.est can result in different results when
#  applied repeatedly using the same data.
#  lasso.est uses the R package glmnet
#
#  choices for lasso.fun:
#  H.lasso
#  OS.lasso
#  RA.lasso  same as LADlasso.Z
#
#  This function
#  returns the smallest absolute value among the estimates after
#  call lasso.est
#  iter times.
#  So the goal is to determine whether zero tends to occur.
#
x=as.matrix(x)
p=ncol(x)
p1=p+1
xy<-cbind(x,y)
xy<-elimna(xy)
x<-xy[,1:p]
y<-xy[,p1]
if(xout){
flag<-outfun(x,plotit=FALSE)$keep
xy<-xy[flag,]
x<-xy[,1:p]
y<-xy[,p1]
}
e1=lasso.fun(x,y)$coef
pl=length(e1)
v=matrix(NA,iter,pl)
v[1,]=e1
for(i in 2:iter)v[i,]=lasso.fun(x,y)$coef
e=apply(abs(v),2,min)
all.e=NULL
if(EST)all.e=apply(v,2,sort)
list(coef=e,all.est=all.e)
}

# ============================================================================
# NN.class
# ============================================================================
NN.class<-function(train=NULL,test=NULL,g=NULL,x1=NULL,x2=NULL,SEED=NULL){
#
# Do classification using the neural network method via the R package neuralnet.
#  This function provides another way of applying this approach using R commands
#  consistent with other classification methods in Rallfun
#
#
# train is the training set
# test is the test data
# g contains labels for the data in the training set,
#
#  Alternatively, store the data for the two groups in
#  x1 and x2, in which case the function creates labels, i.e., no need to specify train and g.
#
#  SEED=NULL, used for convenience when called by other functions that expect SEED
#
if(is.null(test))stop('Argument test is NULL, contains  no data')
library(neuralnet)

if(!is.null(train)){
if(is.null(g))stop('Argument g, group ids, must be specified')
if(is.matrix(g)){
if(dim(g)>1)stop('Argument g should be a vector')
}}
if(!is.null(x1)){
if(!is.null(x2)){
if(ncol(x1)!=ncol(x2))stop('x1 and x2 have different number of columns')
x1=elimna(x1)
x2=elimna(x2)
n1=nrow(x1)
n2=nrow(x2)
g=c(rep(0,n1),rep(1,n2))
train=rbind(x1,x2)
}}
g=as.numeric(as.vector(g))
train=elimna(train)
if(is.null(test))stop('Argument test is null,  no data')
test=elimna(test)
train=as.matrix(train)
test=as.matrix(test)


if(is.vector(test))stop('Argument test is a vector, should contain two or more variables')
#  Next, store data as expected by neuralnet
ddata=as.matrix(cbind(g,train))
dimnames(ddata)=list(NULL,NULL)
p1=ncol(ddata)
p=p1-1
ddata=as.data.frame(ddata)
if(p>8)stop('Current version limited to 8 independent variables')
if(p==2)d=neuralnet(V1~V2+V3,data=ddata)
if(p==3)d=neuralnet(V1~V2+V3+V4,data=ddata)
if(p==4)d=neuralnet(V1~V2+V3+V4+V5,data=ddata)
if(p==5)d=neuralnet(V1~V2+V3+V4+V5+V6,data=ddata)
if(p==6)d=neuralnet(V1~V2+V3+V4+V5+V6+V7,data=ddata)
if(p==7)d=neuralnet(V1~V2+V3+V4+V5+V6+V7+V8,data=ddata)
if(p==8)d=neuralnet(V1~V2+V3+V4+V5+V6+V7+V8+V9,data=ddata)
e=predict(d,newdata=test)
res=rep(1,nrow(test))
flag=e>.5
res[flag]=2
res
}

# ============================================================================
# ridge.est.k
# ============================================================================
ridge.est.k<-function(x,y,STAND=TRUE,locfun=median,scat=madsq,regfun=ols,...){
#
# Estimate the bias parameter k as suggested by Kibria
#
x<-as.matrix(x)
if(nrow(x) != length(y))stop("Length of y does not match number of x values")
m<-cbind(x,y)
m<-elimna(m)
p=ncol(x)
if(p==1)stop('There should be at least two independent variables')
p1=p+1
x=m[,1:p]
y=m[,p1]
n=nrow(x)
df=n-p1
if(STAND){
x=standm(x[,1:2],locfun=median,scat=madsq)
y=y-locfun(y)
}
est=regfun(x,y,...)
C=t(x)%*%x
sighat=sum(est$res^2)/(n-p1)
cr=eigen(C)
D=as.matrix(cr$vectors)          
LAM=diag(cr$values)                                    
XD=x%*%D  #  X^* in Kibria notation                                  
ALPHA=solve(LAM)%*%t(XD)%*%as.matrix(y)
k=sighat/((prod(ALPHA^2))^(1/p)) # Recommended by Kibria
k
}

# ============================================================================
# ridgeGnull
# ============================================================================
ridgeGnull<-function(n,p,regfun=MMreg,iter=5000,SEED=TRUE,MSF=TRUE){
#
# Determine null distribution of ridge.Gtest
#
if(SEED)set.seed(45)
fv=NA
p1=p+1
for(i in 1:iter){
x=rmul(n,p1)
fv[i]=ridgeG.sub(x[,1:p],x[,p1],regfun=regfun,MSF=MSF)$F.test  # stored in ridgeG_sub_chk.tex
}
fv
}

# ============================================================================
# ridgeGnullMC
# ============================================================================
ridgeGnullMC<-function(n,p,regfun=MMreg,iter=5000,SEED=TRUE){
#
# Determine null distribution of ridge.Gtest
#
if(SEED)set.seed(45)
fv=NA
p1=p+1
a=list()
library(parallel)
for(i in 1:iter)a[[i]]=rmul(n,p1)
fv=mclapply(a,ridgeGnullMC.sub,p=p,regfun=regfun)
fv=matl(fv)
fv=as.vector(fv)
fv
}

# ============================================================================
# ridgeG.sub
# ============================================================================
ridgeG.sub<-function(x,y,k=NULL,regfun=tsreg,outfun=outpro,STAND=FALSE,MSF=TRUE,
locfun=mean,scat=var,...){
#
# Test the hypothesis that all slope parameters are zero
#  using a robust analog of a ridge estimator.
#
#
#
# STAND=TRUE: x is standardized and y is centered, based on measures of location
#             and scatter indicatd by
#             locfun and scat
#             locfun=median would use the median and scat=madsq uses MAD.
#
x<-as.matrix(x)
if(nrow(x) != length(y))stop("Length of y does not match number of x values")
m<-cbind(x,y)
m<-elimna(m)
p=ncol(x)
if(p==1)stop('There should be at least two independent variables')
p1=p+1
x=m[,1:p]
y=m[,p1]
n=nrow(x)
if(STAND){
x=standm(x,locfun=locfun,scat=scat)
y=y-locfun(y)
}
if(is.null(k)){
if(!MSF)k=ridge.est.k(x,y)
else{
ires=regfun(x,y)$residuals
sigh=sqrt(sum(ires^2)/(n-p-1))
k=p^(1+1/p)*sigh
}}
est=rob.ridge(x,y,k=k,Regfun=regfun)$coef
est=as.matrix(est)
x<-cbind(rep(1,nrow(x)),x[,1:ncol(x)])
res<-y-x%*%est
p=ncol(x)
kmat=matrix(0,p,p)
diag(kmat)=k
xtx<-solve(t(x)%*%x+kmat)
h<-diag(x%*%xtx%*%t(x))
hc3<-xtx%*%t(x)%*%diag(as.vector(res^2/(1-h)^2))%*%x%*%xtx
slopes=as.matrix(est[2:p1])
Ssq=hc3[2:p1,2:p1]
f.test=t(slopes)%*%solve(Ssq)%*%slopes
f.test=(n-p)*f.test/((n-1)*p)
list(Ridge.est=est,F.test=f.test)
}

# ============================================================================
# ridge.Gtest
# ============================================================================
ridge.Gtest<-function(x,y,k=NULL,regfun=tsreg,xout=FALSE,outfun=outpro,STAND=FALSE,PV=FALSE,iter=5000,
locfun=mean,scat=var,MC=FALSE,MSF=TRUE,...){
#
# Test the hypothesis that all slope parameters are zero.
#  using a robust analog of a ridge estimator.
#
#
# PV=TRUE: computes a p-value at the expense of higher execution time.
# Otherwise, simply test at the 0.05 level using an approximate critical value.
#
#
# STAND=TRUE: x is standardized and y is centered, based on measures of location
#             and scatter indicated by
#             locfun and scat
#             locfun=median would use the median and scat=madsq uses MAD.
#
x<-as.matrix(x)
if(nrow(x) != length(y))stop("Length of y does not match number of x values")
m<-cbind(x,y)
m<-elimna(m)
p=ncol(x)
if(p==1)stop('There should be at least two independent variables')
p1=p+1
x=m[,1:p]
y=m[,p1]
if(xout){
flag<-outfun(x,plotit=FALSE)$keep
x<-x[flag,]
x<-as.matrix(x)
y<-y[flag]
}
n=nrow(x)
if(STAND){
if(!PV)print('Suggest using PV=TRUE when STAND=TRUE')
x=standm(x,locfun=locfun,scat=scat)
y=y-locfun(y)
}
if(is.null(k)){
if(!MSF)k=ridge.est.k(x,y)
else{
ires=ols(x,y)$residuals
sigh=sqrt(sum(ires^2)/(n-p-1))
k=p^(1+1/p)*sigh
}}
est=rob.ridge(x,y,k=k,Regfun=regfun,MSF=MSF)$coef
est=as.matrix(est)
x<-cbind(rep(1,nrow(x)),x[,1:ncol(x)])
res<-y-x%*%est
p=ncol(x)
kmat=matrix(0,p,p)
diag(kmat)=k
xtx<-solve(t(x)%*%x+kmat)
h<-diag(x%*%xtx%*%t(x))
hc3<-xtx%*%t(x)%*%diag(as.vector(res^2/(1-h)^2))%*%x%*%xtx
slopes=as.matrix(est[2:p1])
Ssq=hc3[2:p1,2:p1]
f.test=t(slopes)%*%solve(Ssq)%*%slopes
f.test=(n-p)*f.test/((n-1)*p)
crit.val=NULL
if(!PV){
if(n<20)crit.val=2.56
if(n>500)crit.val=2.24
if(is.null(crit.val)){
nx=c(20,30,40,
50,75,100,
200,500)
ny=c(2.56,2.51,2.40,
2.41,2.34,2.30,
 2.29,2.24)
options(warn=-1)
if(n<=500)crit.f=lplot.pred(1/nx,ny,1/n)$yhat
options(warn=0)
}}
f.test=as.vector(f.test)
pv=NULL
if(PV){
crit.f=NULL
if(!MC)v=ridgeGnull(n,p,regfun=regfun,iter=iter)
if(MC)v=ridgeGnullMC(n,p,regfun=regfun,iter=iter)
pv=mean(v>=f.test)
}
list(n=n,Ridge.est=est,F.test=f.test,critical.05.value=crit.f,p.value=pv)
}

# ============================================================================
# ridge.Liu
# ============================================================================
ridge.Liu<-function(x,y,Regfun=ltsreg,xout=FALSE,outfun=outpro,
plotit=FALSE,d=NULL,tr=.2,
STAND=TRUE,INT=FALSE,locfun=median,...){
#
#     Betul Kan, Ozlem Alpu & Berna Yazici (2013).
#     Robust ridge and robust Liu estimator for regression based on the LTS estimator
#      Journal of Applied Statistics Volume 40,  Issue 3
#
x=as.matrix(x)
xy=elimna(cbind(x,y))
x=as.matrix(x)
p=ncol(x)
p1=p+1
x=xy[,1:p]
y=xy[,p1]
x=as.matrix(x)
if(xout){
flag<-outfun(x,plotit=plotit)$keep
x<-x[flag,]
x<-as.matrix(x)
y<-y[flag]
}
n=nrow(x)
if(STAND){
x=standm(x)
y=(y-mean(y))/sd(y)
}
#  First compute a robust estimate
if(identical(Regfun,ltsreg))EST=ltsreg(x,y,tr=tr)
else EST=Regfun(x,y,...)
sig.hat=sum(EST$residuals^2)/(n-p)
C=t(x)%*%x
cr=eigen(C)
LAM=diag(cr$values)
dlam=cr$values
# Estimate d as in Alheety and Kibria (2009)
if(is.null(d)){
dtop=sum(dlam*(EST$coef[2:p1]^2)/(dlam*(dlam+1)^2))
dbot=sum((sig.hat+dlam*EST$coef[2:p1]^2)/(dlam*(dlam+1)^2))
d=dtop/dbot
}
Ident=diag(1,nrow=p,ncol=p)
#
# Compute the Kan et al.  estimate
Liu.est=solve(LAM+Ident)%*%(LAM+d*Ident)%*%as.matrix(EST$coef[2:p1])
res=y-x%*%Liu.est
Inter=locfun(res)
f.est=c(Inter,as.vector(Liu.est))
list(coef=f.est,d=d)
}

# ============================================================================
# ridge.test
# ============================================================================
ridge.test<-function(x,y,k=NULL,alpha=.05,pr=TRUE,xout=FALSE,outfun=outpro,STAND=TRUE,method='hoch',
locfun=mean,scat=var,MSF=TRUE,...){
#
#
#  Using a ridge estimator
# test the hypothesis of a zero slope for each of p independent variables.
#  If the smallest p-value is less than or equal to alpha, reject the corresponding hypothesis
# as well as the hypothesis that all slopes are zero.
#  But no other slopes can be declared significant due to the bias associated with the
#  ridge estimator.
#
#  The method uses an analog of the heteroscedastic method
# recommended by Long and Ervin (2000).
#  p-values are adjusted to control the probability of one or more Type I errors;'
#  Hochberg's method is used by default.
#
#  Advantage:
#  Power tends to be at least as high as OLS and potentially much higher
#  But when the null hypothsis is false, confidence intervals can be highly
#  inaccurate.
#
#
# STAND=TRUE: x is standardized and y is centered, based on measues of location
#             and scatter indicatd by
#             locfun and scat
#             locfun=median would use the median and scat=madsq uses MAD.
#    For n<=40, this helps control the Type I error probability.
#
x<-as.matrix(x)
if(nrow(x) != length(y))stop('Length of y does not match number of x values')
m<-cbind(x,y)
m<-elimna(m)
p=ncol(x)
if(p==1)stop('There should be at least two independent variables')
p1=p+1
x=m[,1:p]
y=m[,p1]
n=nrow(x)
df=n-p1
if(xout){
if(identical(outfun,outblp))flag=outblp(x,y,plotit=FALSE)$keep
else
 flag<-outfun(x,plotit=FALSE)$keep
x<-x[flag,]
x<-as.matrix(x)
y<-y[flag]
n=nrow(x)
df=n-p1
}
if(STAND){
x=standm(x,locfun=locfun,scat=scat)
y=y-locfun(y)
}
if(is.null(k)){
if(!MSF)k=ridge.est.k(x,y)
else{
ires=ols(x,y)$residuals
sigh=sqrt(sum(ires^2)/(n-p-1))
k=p^(1+1/p)*sigh
}
}
x1=cbind(rep(1,n),x)
ols.est=ols(x,y)$coef
est=ols.ridge(x,y,k=k)$coef
x<-cbind(rep(1,nrow(x)),x[,1:ncol(x)])
res<-y-x1%*%est
p=ncol(x)
kmat=matrix(0,p,p)
diag(kmat)=k
xtx<-solve(t(x)%*%x+kmat)
h<-diag(x%*%xtx%*%t(x))
hc3<-xtx%*%t(x)%*%diag(as.vector(res^2/(1-h)^2))%*%x%*%xtx
df<-nrow(x)-ncol(x)
crit<-qt(1-alpha/2,df)
al<-ncol(x)
ci<-matrix(NA,nrow=p,ncol=6)
se=sqrt(diag(hc3))
p=p-1
for(j in 2:p1){
ci[j,1]=se[j]
ci[j,2]=est[j]/sqrt(hc3[j,j])
ci[j,3]=2*(1-pt(abs(ci[j,2]),df))
ci[j,5]=est[j]
ci[j,6]=ols.est[j]
}
ci[,4]=p.adjust(ci[,3],method=method)
sig='No slope is signficant'
if(sum(ci[2:al,4]<=alpha)>0){
id=which(ci[2:al,4]==min(ci[2:al,4]))
sig=paste('Slope',id, 'is signficant')
}
ci=ci[2:p1,]  # Eliminate results related to the intercept, not relevant
vlabs=NA
for(j in 1:p)vlabs[j]=paste('Slope',j)
dimnames(ci)=list(vlabs,c('s.e.','test.stat','p-value','Adjusted p','Ridge.Est.','OLS.est'))
list(output=ci,Sig.slope=sig)
}

