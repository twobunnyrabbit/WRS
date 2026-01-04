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

#' AdaBoost Classification
#'
#' Performs binary classification using the AdaBoost (Adaptive Boosting) algorithm
#' via the mboost package.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param sm Logical. Currently not used (included for consistency with other functions).
#' @param fr Numeric. Currently not used (included for consistency with other functions).
#' @param SEED Seed value for reproducibility. If NULL, no seed is set. Included for
#'   compatibility when called by other functions.
#' @param baselearner Character string specifying the base learner:
#'   \itemize{
#'     \item `"btree"`: Decision stumps (default)
#'     \item `"bbs"`: Splines
#'     \item `"bols"`: Linear models
#'   }
#'
#' @return A numeric vector of predicted class labels (1 or 2) for the test data.
#'
#' @details
#' AdaBoost is an ensemble learning method that combines multiple weak classifiers
#' to create a strong classifier. This function uses the exponential loss function
#' (AdaExp family in mboost).
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from training and test data.
#'
#' @note This function requires the `mboost` package to be installed.
#'
#' @seealso \code{\link{CLASS.BAG}} for bagged AdaBoost, \code{\link{CLASS.fun}}
#' for general classification dispatcher, \code{\link{class.forest}} for random forests
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # Classify using decision stumps
#' class.ada(x1 = x1, x2 = x2, test = test)
#'
#' # Using train/g specification
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' class.ada(train = train, g = g, test = test)
#'
#' # Using spline base learners
#' class.ada(x1 = x1, x2 = x2, test = test, baselearner = "bbs")
#' }
#'
#' @export
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

#' Bagged Classification Methods
#'
#' Performs binary classification using bagged (bootstrap aggregated) versions of
#' various classification methods. Bagging helps reduce the impact of unequal sample
#' sizes on error rates.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param SEED Logical. If TRUE, sets random seed to 2 for reproducibility.
#' @param kernel Character string specifying SVM kernel type (default: "radial").
#'   Used only when method = "SVM".
#' @param nboot Number of bootstrap samples (default: 100).
#' @param method Character string specifying the classification method:
#'   \itemize{
#'     \item `"KNN"`: Robust K-nearest neighbors using depths (default)
#'     \item `"SVM"`: Support vector machine with bagging
#'     \item `"DIS"`: Depth-based classification (distance method)
#'     \item `"DEP"`: Depth-based (Makinde & Fasoranbaku, 2018)
#'     \item `"PRO"`: Projection method using PDF estimates
#'     \item `"NN"`: Neural network
#'     \item `"RF"`: Random forest
#'     \item `"ADA"`: AdaBoost
#'     \item `"LSM"`: Logistic smoother (or logistic regression if sm=FALSE)
#'   }
#' @param depthfun Function for computing data depth (default: \code{prodepth} for
#'   projection depth). Used for depth-based methods.
#' @param baselearner Character string specifying AdaBoost base learner (default: "bbs").
#'   Used only when method = "ADA".
#' @param sm Logical. If TRUE (default), uses smoother for LSM method. If FALSE,
#'   uses standard logistic regression.
#' @param rule Numeric threshold for classification rule (default: 0.5). Used for LSM method.
#' @param ... Additional arguments passed to specific classification methods.
#'
#' @return A numeric vector of predicted class labels (1 or 2) for the test data.
#'
#' @details
#' Bagging (Bootstrap Aggregating) creates multiple bootstrap samples from the training
#' data, fits a classifier to each sample, and aggregates predictions. This approach
#' reduces variance and makes classification more robust to unequal sample sizes.
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Different methods have different computational requirements and may require specific
#' R packages. For non-bagged versions of some methods, see \code{\link{CLASS.fun}}.
#'
#' @references
#' Makinde, O. S., & Fasoranbaku, O. A. (2018). On classification based on depth functions.
#' Journal of Applied Statistics, 45(8), 1541-1555.
#'
#' @seealso \code{\link{CLASS.fun}} for non-bagged classification methods,
#' \code{\link{class.ada}} for AdaBoost details, \code{\link{class.forest}} for
#' random forests, \code{\link{Depth.class}} for depth-based classification
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # Bagged KNN classification
#' CLASS.BAG(x1 = x1, x2 = x2, test = test, method = "KNN", nboot = 100)
#'
#' # Bagged random forest
#' CLASS.BAG(x1 = x1, x2 = x2, test = test, method = "RF", nboot = 50)
#'
#' # Using train/g specification with depth-based method
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' CLASS.BAG(train = train, g = g, test = test, method = "DEP")
#' }
#'
#' @export
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

#' Classification Error Estimation with Unconditional Confusion Matrix
#'
#' Estimates classification error rates using cross-validation (leave-one-out) and
#' returns an unconditional confusion matrix.
#'
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param train Training data matrix (alternative to x1/x2 specification).
#' @param g Group labels for training data (required if train is provided).
#' @param method Character string specifying classification method (default: "KNN").
#'   See \code{\link{CLASS.fun}} for available methods.
#' @param nboot Number of bootstrap samples for cross-validation (default: 100).
#' @param EN Logical. If TRUE (default), uses equal sample sizes for test data to
#'   address classification bias.
#' @param FAST Logical. If TRUE (default), excludes computationally intensive methods.
#' @param AUC Logical. If TRUE, computes area under ROC curve (requires ROCR package).
#'   Default is FALSE to avoid potential errors.
#' @param SEED Logical. If TRUE (default), sets random seed to 2 for reproducibility.
#' @param ... Additional arguments passed to the classification method.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `C.MAT`: Unconditional confusion matrix (proportions)
#'     \item `COUNTS`: Confusion matrix with raw counts
#'     \item `AUC`: Area under ROC curve (if AUC=TRUE), otherwise NULL
#'   }
#'
#' @details
#' This function estimates classification error rates using cross-validation with
#' leave-one-out sampling. It returns an **unconditional** confusion matrix showing
#' the joint probabilities of true and predicted classes.
#'
#' The confusion matrix shows:
#' - Rows: True class membership
#' - Columns: Predicted class membership
#'
#' For a **conditional** confusion matrix (showing probabilities conditional on true
#' class), use \code{\link{class.error.CP}} instead.
#'
#' @note Requires the `ROCR` package if AUC=TRUE.
#'
#' @seealso \code{\link{class.error.CP}} for conditional confusion matrix,
#' \code{\link{class.error.loo}} for leave-one-out error rates,
#' \code{\link{class.error.com}} for comparing methods, \code{\link{CLASS.fun}}
#' for classification methods
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Estimate error with KNN
#' class.error.CM(x1 = x1, x2 = x2, method = "KNN", nboot = 50)
#'
#' # With AUC calculation
#' class.error.CM(x1 = x1, x2 = x2, method = "RF", AUC = TRUE)
#' }
#'
#' @export
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

#' Compare Classification Methods Using Cross-Validation
#'
#' Compares multiple classification methods using cross-validation coupled with
#' resampling to estimate error rates (total error, false positives, false negatives).
#'
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param method Character vector specifying classification methods to compare. Available
#'   methods: 'KNN' (K-nearest neighbors), 'DIS' (depth-based), 'DEP' (Makinde & Fasoranbaku
#'   depth), 'SVM' (support vector machine), 'RF' (random forest), 'NN' (neural network),
#'   'ADA' (AdaBoost), 'PRO' (projection-based), 'LSM' (logistic/smooth logistic regression).
#'   If NULL, all methods are compared (except 'PRO' if FAST=TRUE).
#' @param pro.p Proportion of data to use for training (default: 0.8 = 80%).
#' @param nboot Number of bootstrap samples for cross-validation (default: 100).
#' @param EN Logical. If TRUE, use equal sample sizes for both groups to deal with
#'   classification bias. If FALSE, use the natural ratio n1/n2 (default: FALSE).
#' @param FAST Logical. If TRUE (default), excludes 'PRO' method which can be slow for
#'   large samples.
#' @param SEED Logical or numeric. If TRUE, sets seed to 2 for reproducibility. If numeric,
#'   sets that seed value (default: TRUE).
#' @param ... Additional arguments passed to individual classification methods.
#'
#' @return A list with component:
#'   \item{Error.rates}{A 3 x J matrix (J = number of methods) with rows:
#'     \itemize{
#'       \item TE: Total Error rate (proportion of misclassifications)
#'       \item FP: False Positive rate (proportion of x1 classified as x2)
#'       \item FN: False Negative rate (proportion of x2 classified as x1)
#'     }}
#'
#' @details
#' This function performs cross-validation by:
#' 1. Randomly sampling `pro.p` proportion of data for training
#' 2. Using the remaining data for testing
#' 3. Repeating `nboot` times
#' 4. Computing average error rates across all bootstrap samples
#'
#' **Error rate definitions:**
#' - **Total Error (TE)**: Overall proportion of misclassifications
#' - **False Positive (FP)**: Proportion of group 1 observations misclassified as group 2
#'   (e.g., predict fracture when no fracture occurred)
#' - **False Negative (FN)**: Proportion of group 2 observations misclassified as group 1
#'   (e.g., predict no fracture when fracture occurred)
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from data.
#'
#' @note
#' - The 'PRO' method can take several minutes for large samples
#' - Different methods may require different packages to be installed
#'
#' @seealso \code{\link{class.error.loo}} for leave-one-out cross-validation,
#' \code{\link{class.error.CP}} for confusion matrix with single method,
#' \code{\link{CLASS.fun}} for general classification dispatcher
#'
#' @examples
#' \dontrun{
#' # Generate two groups
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Compare all fast methods
#' class.error.com(x1 = x1, x2 = x2)
#'
#' # Compare specific methods
#' class.error.com(x1 = x1, x2 = x2, method = c('KNN', 'RF', 'SVM'))
#'
#' # Use 90% of data for training
#' class.error.com(x1 = x1, x2 = x2, pro.p = 0.9)
#'
#' # Use equal sample sizes for testing
#' class.error.com(x1 = x1, x2 = x2, EN = TRUE)
#' }
#'
#' @export
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

#' Classification Error with Confusion Matrix (Leave-One-Out)
#'
#' Estimates classification error for a single classification method using leave-one-out
#' cross-validation and returns a confusion matrix showing classification performance.
#'
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param method Character string specifying single classification method. Available
#'   methods: 'KNN' (K-nearest neighbors, default), 'DIS' (depth-based), 'DEP' (Makinde &
#'   Fasoranbaku depth), 'SVM' (support vector machine), 'RF' (random forest), 'NN' (neural
#'   network), 'ADA' (AdaBoost), 'PRO' (projection-based), 'LSM' (logistic regression).
#' @param nboot Number of samples to use for leave-one-out cross-validation (default: 100).
#'   Limited by min(n1, n2).
#' @param EN Logical. If TRUE (default), use equal sample sizes for both groups to deal
#'   with classification bias.
#' @param FAST Logical. Currently not used (included for consistency).
#' @param AUC Logical. If TRUE, computes area under ROC curve. Default is FALSE because
#'   certain conditions can create errors. Requires ROCR package.
#' @param SEED Logical or numeric. If TRUE, sets seed to 2 for reproducibility. If numeric,
#'   sets that seed value (default: TRUE).
#' @param ... Additional arguments passed to the classification method.
#'
#' @return A list with components:
#'   \item{C.MAT}{Confusion matrix showing proportions: rows are true classes, columns
#'     are predicted classes}
#'   \item{COUNTS}{Confusion matrix showing counts instead of proportions}
#'   \item{AUC}{Area under ROC curve (if AUC=TRUE), otherwise NULL}
#'
#' @details
#' This function performs leave-one-out cross-validation:
#' 1. Randomly selects one observation from each group
#' 2. Trains classifier on remaining data
#' 3. Predicts class for the two held-out observations
#' 4. Repeats `nboot` times (or min(n1, n2) times, whichever is smaller)
#' 5. Summarizes results in a confusion matrix
#'
#' **Confusion matrix interpretation:**
#' - Rows: True class labels
#' - Columns: Predicted class labels
#' - Diagonal elements: Correct classifications
#' - Off-diagonal elements: Misclassifications
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from data.
#'
#' @note
#' - Only one method can be analyzed at a time
#' - ROCR package required if AUC=TRUE
#' - Different methods may require different packages to be installed
#'
#' @seealso \code{\link{class.error.loo}} for multiple methods with leave-one-out,
#' \code{\link{class.error.com}} for cross-validation comparison,
#' \code{\link{CLASS.fun}} for general classification dispatcher
#'
#' @examples
#' \dontrun{
#' # Generate two groups
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Get confusion matrix for KNN
#' result <- class.error.CP(x1 = x1, x2 = x2)
#' result$C.MAT  # Proportion matrix
#' result$COUNTS # Count matrix
#'
#' # Try different method
#' class.error.CP(x1 = x1, x2 = x2, method = 'RF')
#'
#' # Include AUC
#' class.error.CP(x1 = x1, x2 = x2, AUC = TRUE)
#' }
#'
#' @export
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

#' Classification Error Using Leave-One-Out Cross-Validation
#'
#' Estimates classification error rates (total error, false positives, false negatives,
#' true positives, true negatives) for one or more classification methods using
#' leave-one-out cross-validation.
#'
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param method Character vector specifying classification methods to compare. Available
#'   methods: 'KNN' (K-nearest neighbors), 'DIS' (depth-based), 'DEP' (Makinde & Fasoranbaku
#'   depth), 'SVM' (support vector machine), 'RF' (random forest), 'NN' (neural network),
#'   'ADA' (AdaBoost), 'PRO' (projection-based), 'LSM' (logistic/smooth logistic regression),
#'   'GBT' (gradient boosting). If NULL, all methods are compared (except 'PRO' and 'GBT'
#'   if FAST=TRUE).
#' @param nboot Number of samples to use for leave-one-out cross-validation (default: 100).
#'   Limited by min(n1, n2).
#' @param EN Logical. If TRUE (default), use equal sample sizes for both groups to deal
#'   with classification bias.
#' @param FAST Logical. If TRUE (default), excludes 'PRO' and 'GBT' methods which can be
#'   slow for large samples.
#' @param SEED Logical or numeric. If TRUE, sets seed to 2 for reproducibility. If numeric,
#'   sets that seed value (default: TRUE).
#' @param ... Additional arguments passed to individual classification methods.
#'
#' @return A list with component:
#'   \item{Error.rates}{A 5 x J matrix (J = number of methods) with rows:
#'     \itemize{
#'       \item TE: Total Error rate (proportion of misclassifications)
#'       \item FP: False Positive rate (proportion of x1 classified as x2)
#'       \item FN: False Negative rate (proportion of x2 classified as x1)
#'       \item TP: True Positive rate (proportion of x1 correctly classified as x1)
#'       \item TN: True Negative rate (proportion of x2 correctly classified as x2)
#'     }}
#'
#' @details
#' This function performs leave-one-out cross-validation:
#' 1. Randomly selects one observation from each group
#' 2. Trains classifier on remaining data
#' 3. Predicts class for the two held-out observations
#' 4. Repeats `nboot` times (or min(n1, n2) times, whichever is smaller)
#' 5. Computes average error rates across all samples
#'
#' **Error rate definitions:**
#' - **Total Error (TE)**: Overall proportion of misclassifications
#' - **False Positive (FP)**: Proportion of group 1 observations misclassified as group 2
#' - **False Negative (FN)**: Proportion of group 2 observations misclassified as group 1
#' - **True Positive (TP)**: Proportion of group 1 observations correctly classified as group 1
#' - **True Negative (TN)**: Proportion of group 2 observations correctly classified as group 2
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from data.
#'
#' @note
#' - The 'PRO' method can take several minutes for large samples
#' - Different methods may require different packages to be installed
#'
#' @seealso \code{\link{class.error.com}} for cross-validation with training/test split,
#' \code{\link{class.error.CP}} for confusion matrix with single method,
#' \code{\link{CLASS.fun}} for general classification dispatcher
#'
#' @examples
#' \dontrun{
#' # Generate two groups
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Compare all fast methods
#' class.error.loo(x1 = x1, x2 = x2)
#'
#' # Compare specific methods
#' class.error.loo(x1 = x1, x2 = x2, method = c('KNN', 'RF', 'SVM'))
#'
#' # Use equal sample sizes for testing
#' class.error.loo(x1 = x1, x2 = x2, EN = TRUE)
#' }
#'
#' @export
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

#' Random Forest Classification
#'
#' Performs binary classification using random forest via the randomForest package.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param sm Logical. Currently not used (included for consistency with other functions).
#' @param fr Numeric. Currently not used (included for consistency with other functions).
#' @param SEED Seed value for reproducibility. If NULL, no seed is set. Included for
#'   compatibility when called by other functions.
#'
#' @return A numeric vector of predicted class labels (1 or 2) for the test data.
#'
#' @details
#' Random forests are ensemble learning methods that construct multiple decision trees
#' during training and output the class that is the mode of the classes from individual trees.
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from training and test data.
#'
#' @note This function requires the `randomForest` package to be installed.
#'
#' @seealso \code{\link{CLASS.fun}} for general classification dispatcher,
#' \code{\link{CLASS.BAG}} for bagged random forest, \code{\link{class.ada}}
#' for AdaBoost
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # Random forest classification
#' class.forest(x1 = x1, x2 = x2, test = test)
#'
#' # Using train/g specification
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' class.forest(train = train, g = g, test = test)
#' }
#'
#' @export
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

#' General Classification Method Dispatcher
#'
#' Dispatcher function providing access to multiple binary classification methods.
#' This is a unified interface for various classification algorithms.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param method Character string specifying the classification method:
#'   \itemize{
#'     \item `"KNN"`: K-nearest neighbors using robust depths (default)
#'     \item `"DIS"`: Depth-based classification (distance method)
#'     \item `"DEP"`: Depth-based (Makinde & Fasoranbaku, 2018)
#'     \item `"SVM"`: Support vector machine
#'     \item `"RF"`: Random forest
#'     \item `"NN"`: Neural network
#'     \item `"ADA"`: AdaBoost
#'     \item `"PRO"`: Projection method using PDF estimates
#'     \item `"LSM"`: Logistic smoother (or logistic regression if sm=FALSE)
#'     \item `"GBT"`: Gradient boosted trees
#'   }
#' @param depthfun Function for computing data depth (default: \code{prodepth} for
#'   projection depth). Used for depth-based methods.
#' @param kernel Character string specifying SVM kernel type (default: "radial").
#' @param baselearner Character string specifying AdaBoost base learner (default: "btree").
#' @param sm Logical. If TRUE (default), uses smoother for LSM method. If FALSE,
#'   uses standard logistic regression.
#' @param rule Numeric threshold for classification rule (default: 0.5). Used for LSM method.
#' @param ... Additional arguments passed to specific classification methods.
#'
#' @return A numeric vector of predicted class labels (1 or 2) for the test data.
#'
#' @details
#' This function provides a unified interface to multiple classification methods.
#' It dispatches to the appropriate classification function based on the `method` argument.
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Different methods may require different R packages to be installed. The function
#' handles loading required packages automatically.
#'
#' @references
#' Makinde, O. S., & Fasoranbaku, O. A. (2018). On classification based on depth functions.
#' Journal of Applied Statistics, 45(8), 1541-1555.
#'
#' @seealso \code{\link{CLASS.BAG}} for bagged versions, \code{\link{class.ada}} for
#' AdaBoost, \code{\link{class.forest}} for random forests, \code{\link{Depth.class}}
#' for depth-based classification, \code{\link{KNN}} for K-nearest neighbors
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # K-nearest neighbors classification
#' CLASS.fun(x1 = x1, x2 = x2, test = test, method = "KNN")
#'
#' # Random forest classification
#' CLASS.fun(x1 = x1, x2 = x2, test = test, method = "RF")
#'
#' # Depth-based classification
#' CLASS.fun(x1 = x1, x2 = x2, test = test, method = "DEP")
#'
#' # Using train/g specification with SVM
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' CLASS.fun(train = train, g = g, test = test, method = "SVM")
#' }
#'
#' @export
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

#' Gradient Boosted Classification
#'
#' Performs binary classification using gradient boosted trees via the gbm package.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param SEED Seed value for reproducibility. If NULL, no seed is set. Included for
#'   compatibility when called by other functions.
#' @param n.trees Number of trees to fit (default: 100).
#'
#' @return A numeric vector of predicted class labels (1 or 2) for the test data.
#'
#' @details
#' Gradient boosting is an ensemble method that builds trees sequentially, with each
#' tree attempting to correct errors made by previous trees. This implementation uses
#' the Bernoulli distribution for binary classification.
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from training and test data.
#'
#' @note This function requires the `gbm` package to be installed.
#'
#' @seealso \code{\link{CLASS.fun}} for general classification dispatcher,
#' \code{\link{class.forest}} for random forests, \code{\link{class.ada}} for AdaBoost
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # Gradient boosted classification with 100 trees
#' class.gbm(x1 = x1, x2 = x2, test = test, n.trees = 100)
#'
#' # Using train/g specification with 200 trees
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' class.gbm(train = train, g = g, test = test, n.trees = 200)
#' }
#'
#' @export
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

#' Logistic Regression or Smoother Classification
#'
#' Performs binary classification using either logistic regression or a smoother for
#' binary outcomes.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param sm Logical. If TRUE (default), uses a smoother designed for binary outcomes.
#'   If FALSE, uses standard logistic regression.
#' @param fr Numeric. Span parameter for the smoother (default: 2). Used only when sm=TRUE.
#' @param rule Numeric threshold for classification (default: 0.5). Predicted probabilities
#'   greater than `rule` are classified as group 2, otherwise group 1.
#' @param SEED Seed value for reproducibility. If NULL, no seed is set. Included for
#'   compatibility when called by other functions.
#'
#' @return A numeric vector of predicted class labels (1 or 2) for the test data.
#'
#' @details
#' This function provides two classification approaches:
#' - When `sm=TRUE`: Uses a smoothing method designed for binary dependent variables
#' - When `sm=FALSE`: Uses standard logistic regression
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from training and test data.
#'
#' The classification rule assigns observations to group 2 if the predicted probability
#' exceeds `rule`, and to group 1 otherwise.
#'
#' @seealso \code{\link{CLASS.fun}} for general classification dispatcher,
#' \code{\link{CLASS.BAG}} for bagged version
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # Classification using smoother
#' class.logR(x1 = x1, x2 = x2, test = test, sm = TRUE)
#'
#' # Classification using logistic regression
#' class.logR(x1 = x1, x2 = x2, test = test, sm = FALSE)
#'
#' # Using different classification threshold
#' class.logR(x1 = x1, x2 = x2, test = test, rule = 0.6)
#' }
#'
#' @export
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

#' Univariate Classification Error Using Kernel Density Estimation
#'
#' For univariate (one-dimensional) data, estimates classification error using
#' kernel density estimation. Returns both conditional and unconditional error
#' probabilities.
#'
#' @param x Numeric vector of observations from group 1.
#' @param y Numeric vector of observations from group 2.
#' @param xy Alternative data input: matrix or data frame with all observations
#'   (requires `g` to specify group labels). If NULL, `x` and `y` must be provided.
#' @param g Group labels when using `xy` input. Should be a vector with two unique
#'   values indicating the two groups.
#'
#' @return A list with components:
#'   \item{Conditional.prob}{A 2x2 matrix of conditional probabilities: rows are true
#'     classes, columns are predicted classes. Entry [i,j] is P(predict j | true i).}
#'   \item{Unconditiion.prob}{A 2x2 matrix of unconditional (joint) probabilities.
#'     Entry [i,j] is P(true i AND predict j).}
#'   \item{prob.correct.decision}{Probability of correct classification (sum of
#'     diagonal elements of unconditional probability matrix).}
#'
#' @details
#' This function uses kernel density estimation (via \code{\link{akerd}}) to estimate
#' the probability density function for each group. For each observation, classification
#' is based on which group's density is higher at that point.
#'
#' **Conditional probabilities** answer questions like:
#' - P(classify as group 2 | actually group 1): Conditional false positive rate
#' - P(classify as group 1 | actually group 1): Conditional true positive rate
#'
#' **Unconditional probabilities** answer questions like:
#' - P(actually group 1 AND correctly classified as group 1)
#' - P(actually group 1 BUT misclassified as group 2)
#'
#' The data can be specified in two ways:
#' 1. Using `x` and `y`: Provide the two groups separately as numeric vectors
#' 2. Using `xy` and `g`: Provide all data in `xy` with group labels in `g`
#'
#' Missing values are automatically removed.
#'
#' @note
#' This function is designed specifically for univariate (one-dimensional) data.
#' For multivariate data, use other classification methods like \code{\link{CLASS.fun}}.
#'
#' @seealso \code{\link{akerd}} for adaptive kernel density estimation,
#' \code{\link{class.error.com}} for multivariate classification error,
#' \code{\link{CLASS.fun}} for general classification
#'
#' @examples
#' \dontrun{
#' # Generate univariate data from two groups
#' set.seed(123)
#' x <- rnorm(100, mean = 0)
#' y <- rnorm(100, mean = 1.5)
#'
#' # Estimate classification error
#' result <- class.uni.error(x, y)
#'
#' # View conditional probabilities
#' result$Conditional.prob
#'
#' # View unconditional probabilities
#' result$Unconditiion.prob
#'
#' # Overall correct classification probability
#' result$prob.correct.decision
#'
#' # Alternative input format
#' xy <- c(x, y)
#' g <- c(rep(1, length(x)), rep(2, length(y)))
#' class.uni.error(xy = xy, g = g)
#' }
#'
#' @export
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

#' Depth-Based Classification
#'
#' Performs binary classification using data depth, as suggested by Makinde and
#' Fasoranbaku (2018).
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param depthfun Function for computing data depth (default: \code{prodepth} for
#'   projection depth). Must have the form `depthfun(x, pts)` where it computes
#'   depth for points in `pts` relative to points in `x`.
#' @param DIST Logical. If TRUE (default), uses depth distribution method. If FALSE,
#'   uses maximum depth method.
#' @param SEED Logical. If TRUE (default), sets random seed for reproducibility.
#' @param ... Additional arguments passed to the depth function.
#'
#' @return A numeric vector of predicted class labels for the test data. The labels
#'   correspond to the unique values in the training group labels.
#'
#' @details
#' This function implements depth-based classification using two approaches:
#'
#' - **Maximum Depth (DIST=FALSE)**: Classifies each test point to the group for which
#'   it has the maximum depth.
#' - **Depth Distribution (DIST=TRUE)**: Classifies based on the cumulative distribution
#'   of depths, comparing where the test point's depth falls in each group's depth
#'   distribution. Recommended method.
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from training and test data.
#'
#' Common depth functions include:
#' - `prodepth`: Projection depth (default, robust)
#' - `zonoid`: Zonoid depth
#' - Other depth functions following the same interface
#'
#' @references
#' Makinde, O. S., & Fasoranbaku, O. A. (2018). On classification based on depth functions.
#' Journal of Applied Statistics, 45(8), 1541-1555.
#'
#' @note This function requires the `class` package to be installed.
#'
#' @seealso \code{\link{CLASS.fun}} for general classification dispatcher,
#' \code{\link{KNNdist}} for KNN with depths, \code{\link{CLASS.BAG}} for bagged version
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # Depth distribution method (recommended)
#' Depth.class(x1 = x1, x2 = x2, test = test, DIST = TRUE)
#'
#' # Maximum depth method
#' Depth.class(x1 = x1, x2 = x2, test = test, DIST = FALSE)
#'
#' # Using train/g specification
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' Depth.class(train = train, g = g, test = test)
#' }
#'
#' @export
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

#' K-Means Clustering with Optional Outlier Removal
#'
#' Performs K-means clustering with an option to remove outliers before clustering.
#'
#' @param x Data matrix where rows are observations and columns are variables.
#' @param k Number of clusters.
#' @param xout Logical. If TRUE, outliers are removed before clustering using the
#'   function specified in `outfun`. Default is FALSE.
#' @param outfun Function for detecting outliers (default: \code{outpro} for projection-based
#'   outlier detection). Only used when xout=TRUE.
#'
#' @return A kmeans object containing:
#'   \itemize{
#'     \item `cluster`: A vector of integers (from 1:k) indicating the cluster to which each point is allocated
#'     \item `centers`: A matrix of cluster centers
#'     \item `totss`: The total sum of squares
#'     \item `withinss`: Vector of within-cluster sum of squares, one component per cluster
#'     \item `tot.withinss`: Total within-cluster sum of squares
#'     \item `betweenss`: The between-cluster sum of squares
#'     \item `size`: The number of points in each cluster
#'   }
#'
#' @details
#' This function wraps the standard `kmeans` function with an option to remove outliers
#' before clustering. Removing outliers can improve clustering when the data contains
#' extreme values that would otherwise distort the cluster centers.
#'
#' Missing values are automatically removed before analysis.
#'
#' @seealso \code{\link{Kmeans.grp}} for K-means with group sorting,
#' \code{\link{kmeans}} for the underlying clustering function
#'
#' @examples
#' \dontrun{
#' # Generate data with 3 clusters
#' set.seed(123)
#' x <- rbind(
#'   matrix(rnorm(100, mean = 0), ncol = 2),
#'   matrix(rnorm(100, mean = 3), ncol = 2),
#'   matrix(rnorm(100, mean = 6), ncol = 2)
#' )
#'
#' # K-means clustering without outlier removal
#' res1 <- Kmeans(x, k = 3, xout = FALSE)
#' res1$centers
#'
#' # K-means clustering with outlier removal
#' res2 <- Kmeans(x, k = 3, xout = TRUE)
#' res2$centers
#' }
#'
#' @export
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

#' K-Means Clustering with Group Sorting
#'
#' Performs K-means clustering on one set of variables (x) and sorts a second set
#' of variables (y) into groups based on the cluster assignments.
#'
#' @param x Data matrix used for clustering. Rows are observations, columns are
#'   clustering variables.
#' @param k Number of clusters.
#' @param y Data matrix to be sorted into groups. Must have the same number of rows as x.
#' @param xout Logical. If TRUE, outliers are removed from both x and y before clustering
#'   using the function specified in `outfun`. Default is FALSE.
#' @param outfun Function for detecting outliers (default: \code{outpro} for projection-based
#'   outlier detection). Only used when xout=TRUE.
#'
#' @return A list of length k, where each element contains the rows of y corresponding
#'   to observations assigned to that cluster.
#'
#' @details
#' This function is useful when you want to cluster observations based on one set of
#' variables (x) but then examine patterns in a different set of variables (y) within
#' those clusters.
#'
#' The function:
#' 1. Applies K-means clustering to x
#' 2. Identifies which cluster each observation belongs to
#' 3. Sorts the corresponding rows of y into k groups based on cluster membership
#'
#' Missing values are automatically removed before analysis. When outliers are removed
#' (xout=TRUE), the same observations are removed from both x and y to maintain alignment.
#'
#' @seealso \code{\link{Kmeans}} for K-means clustering without group sorting,
#' \code{\link{kmeans}} for the underlying clustering function
#'
#' @examples
#' \dontrun{
#' # Generate clustering variables
#' set.seed(123)
#' x <- rbind(
#'   matrix(rnorm(100, mean = 0), ncol = 2),
#'   matrix(rnorm(100, mean = 3), ncol = 2)
#' )
#'
#' # Generate response variables to sort
#' y <- matrix(rnorm(100), ncol = 2)
#'
#' # Cluster based on x, sort y into groups
#' groups <- Kmeans.grp(x, k = 2, y)
#'
#' # Examine y values in each cluster
#' lapply(groups, colMeans)
#' }
#'
#' @export
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

#' K-Nearest Neighbors Classification
#'
#' Performs binary classification using the K-nearest neighbors (KNN) method via
#' the class package.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param k Number of nearest neighbors to use for classification (default: 3).
#' @param prob Logical. If TRUE (default), returns probability estimates. If FALSE,
#'   returns only class labels.
#' @param SEED Seed value for reproducibility. If NULL, no seed is set. Included for
#'   compatibility when called by other functions.
#'
#' @return A factor vector of predicted class labels for the test data, with attributes
#'   containing probability estimates if prob=TRUE.
#'
#' @details
#' K-nearest neighbors is a non-parametric classification method that assigns a test
#' observation to the most common class among its k nearest training observations.
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately (function creates labels automatically)
#'
#' Missing values are automatically removed from training and test data.
#'
#' The choice of k affects the bias-variance tradeoff:
#' - Small k: More flexible, lower bias, higher variance
#' - Large k: Less flexible, higher bias, lower variance
#'
#' @note This function requires the `class` package to be installed.
#'
#' @seealso \code{\link{KNNdist}} for KNN using data depths, \code{\link{KNNv2}} for
#' alternative KNN implementation, \code{\link{CLASS.fun}} for general classification
#' dispatcher, \code{\link{knn}} for the underlying function
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # KNN classification with k=3
#' KNN(x1 = x1, x2 = x2, test = test, k = 3)
#'
#' # KNN classification with k=5
#' KNN(x1 = x1, x2 = x2, test = test, k = 5)
#'
#' # Using train/g specification
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' KNN(train = train, g = g, test = test, k = 3)
#' }
#'
#' @export
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

#' K-Nearest Neighbors Classification Using Data Depths
#'
#' Performs binary classification using the K-nearest neighbors (KNN) method applied
#' to data depths. Points are transformed to their depth in each group (DD-plot) and
#' KNN is applied using the resulting depth coordinates.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param k Number of nearest neighbors to use (default: 3).
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param prob Logical. If TRUE, returns class probabilities instead of class labels
#'   (default: FALSE).
#' @param plotit Logical. If TRUE, creates a DD-plot (depth-depth plot) showing the
#'   transformation (default: FALSE).
#' @param SEED Seed value for reproducibility. If NULL, no seed is set. Included for
#'   compatibility when called by other functions.
#' @param xlab Label for x-axis in DD-plot (default: 'Group 1').
#' @param ylab Label for y-axis in DD-plot (default: 'Group 2').
#' @param depthfun Function for computing data depth (default: \code{prodepth} for
#'   projection depth). Must have the form `depthfun(x, pts)` where it computes
#'   depth for points in `pts` relative to points in `x`.
#' @param ... Additional arguments passed to the depth function.
#'
#' @return A numeric vector of predicted class labels for the test data (or class
#'   probabilities if `prob=TRUE`). The labels correspond to the unique values in
#'   the training group labels.
#'
#' @details
#' This function implements the DD-classifier (Depth-Depth classifier) approach:
#' 1. For each point, compute its depth in group 1
#' 2. For each point, compute its depth in group 2
#' 3. Transform each point to 2D coordinates (depth1, depth2)
#' 4. Apply standard KNN classification in this depth space
#'
#' This transformation can improve classification performance by reducing the
#' dimensionality to 2D (regardless of original dimension) and making the decision
#' boundary more linear in the transformed space.
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from data.
#'
#' Common depth functions include:
#' - `prodepth`: Projection depth (default, robust)
#' - `zonoid`: Zonoid depth
#' - Other depth functions following the interface `depthfun(x, pts)`
#'
#' @references
#' Li, J., Cuesta-Albertos, J. A., & Liu, R. Y. (2012). DD-classifier: Nonparametric
#' classification procedure based on DD-plot. Journal of the American Statistical
#' Association, 107(498), 737-753.
#'
#' @note This function requires the `class` package to be installed.
#'
#' @seealso \code{\link{KNN}} for standard KNN, \code{\link{Depth.class}} for
#' depth-based classification, \code{\link{CLASS.fun}} for general classification
#' dispatcher
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # KNN classification using depths
#' KNNdist(x1 = x1, x2 = x2, test = test, k = 5)
#'
#' # With DD-plot visualization
#' KNNdist(x1 = x1, x2 = x2, test = test, plotit = TRUE)
#'
#' # Using train/g specification
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' KNNdist(train = train, g = g, test = test, k = 5)
#' }
#'
#' @export
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

#' K-Nearest Neighbors Classification (Version 2)
#'
#' Performs binary classification using the K-nearest neighbors (KNN) method. This is
#' an alternative version of \code{\link{KNN}} that returns numeric class labels (1, 2)
#' instead of factor levels.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param k Number of nearest neighbors to use (default: 3).
#' @param prob Logical. If TRUE (default), includes class probabilities in the result.
#' @param SEED Seed value for reproducibility. If NULL, no seed is set. Included for
#'   compatibility when called by other functions.
#'
#' @return A numeric vector of predicted class labels (1 or 2) for the test data.
#'
#' @details
#' This function is a wrapper for \code{knn} from the class package. The key difference
#' from \code{\link{KNN}} is that it returns numeric class labels (1, 2) instead of
#' factor levels, which can be more convenient for downstream processing.
#'
#' **KNN algorithm:**
#' 1. For each test point, find the k nearest training points (by Euclidean distance)
#' 2. Assign the class based on majority vote among the k neighbors
#' 3. If `prob=TRUE`, also compute the proportion of votes for the winning class
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from training and test data.
#'
#' @note
#' This function requires the `class` package to be installed. For more information
#' about the underlying KNN implementation, see `?knn` from the class package.
#'
#' @seealso \code{\link{KNN}} for standard KNN, \code{\link{KNNdist}} for KNN with
#' data depths, \code{\link{CLASS.fun}} for general classification dispatcher
#'
#' @examples
#' \dontrun{
#' # Generate training data
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # KNN classification with k=3
#' KNNv2(x1 = x1, x2 = x2, test = test, k = 3)
#'
#' # Try different k values
#' KNNv2(x1 = x1, x2 = x2, test = test, k = 5)
#' KNNv2(x1 = x1, x2 = x2, test = test, k = 7)
#'
#' # Using train/g specification
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' KNNv2(train = train, g = g, test = test, k = 3)
#' }
#'
#' @export
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

#' LASSO Regression with Outlier Detection
#'
#' Performs LASSO (Least Absolute Shrinkage and Selection Operator) regression
#' via the glmnet package, with options for eliminating leverage points and
#' regression outliers before fitting.
#'
#' @param x Predictor matrix (n x p). Can be a matrix or data frame.
#' @param y Response vector of length n.
#' @param xout Logical. If TRUE, eliminate leverage points using `outfun` before
#'   fitting (default: FALSE).
#' @param STAND Logical. If TRUE (default), standardize predictors before fitting.
#' @param outfun Function for detecting outliers/leverage points (default: \code{outpro}).
#'   Only used if `xout=TRUE` or `regout=TRUE`.
#' @param regout Logical. If TRUE, eliminate regression outliers using the \code{elo}
#'   function before fitting (default: FALSE). This option removes both leverage points
#'   and regression outliers.
#' @param lam Optional vector of lambda values for cross-validation. If NULL, glmnet
#'   chooses the lambda sequence automatically.
#' @param ... Additional arguments passed to the outlier detection function `outfun`.
#'
#' @return A list with components:
#'   \item{coef}{Vector of estimated coefficients (including intercept). Length is
#'     p+1 where the first element is the intercept.}
#'   \item{lambda.min.used}{The lambda value that minimizes cross-validation error.}
#'
#' @details
#' This function is a convenience wrapper for \code{cv.glmnet} from the glmnet package,
#' providing simplified output (coefficients only) and built-in outlier removal options.
#'
#' **LASSO regression** performs variable selection and regularization by adding an L1
#' penalty to the least squares objective. This shrinks some coefficients exactly to zero,
#' effectively performing automatic variable selection.
#'
#' **Outlier removal options:**
#' - `xout=TRUE`: Remove leverage points (outliers in predictor space) using `outfun`
#' - `regout=TRUE`: Remove both leverage points AND regression outliers (outliers in
#'   residual space) using the \code{elo} function
#' - If `regout=TRUE`, the `xout` parameter is ignored (elo handles both types)
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' - This function requires the `glmnet` package to be installed
#' - LASSO estimates can vary slightly across runs due to random aspects of cross-validation.
#'   Use \code{\link{lasso.rep}} for more stable estimates via repeated fitting.
#' - For more complete output and options, use `cv.glmnet` directly from glmnet package
#'
#' @seealso \code{\link{lasso.rep}} for repeated LASSO to identify stable coefficients,
#' \code{\link{ridge.test}} for ridge regression, \code{\link{elo}} for regression
#' outlier detection
#'
#' @examples
#' \dontrun{
#' # Generate data
#' set.seed(123)
#' x <- matrix(rnorm(100 * 5), ncol = 5)
#' y <- x[,1] + 2*x[,2] + rnorm(100)
#'
#' # Basic LASSO
#' result <- lasso.est(x, y)
#' result$coef  # Coefficients (intercept + 5 slopes)
#'
#' # With leverage point removal
#' lasso.est(x, y, xout = TRUE)
#'
#' # With regression outlier removal
#' lasso.est(x, y, regout = TRUE)
#'
#' # Without standardization
#' lasso.est(x, y, STAND = FALSE)
#' }
#'
#' @export
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

#' Repeated LASSO Regression for Stable Coefficient Identification
#'
#' Performs LASSO regression repeatedly and identifies stable coefficients by
#' finding the minimum absolute value across repeated fits. Helps determine which
#' coefficients consistently tend toward zero (indicating they should be excluded).
#'
#' @param x Predictor matrix (n x p). Can be a matrix or data frame.
#' @param y Response vector of length n.
#' @param lasso.fun LASSO function to use (default: \code{lasso.est}). Alternative
#'   choices include `H.lasso`, `OS.lasso`, `RA.lasso` (same as `LADlasso.Z`).
#' @param xout Logical. If TRUE, eliminate leverage points using `outfun` before
#'   fitting (default: FALSE).
#' @param outfun Function for detecting outliers/leverage points (default: \code{outpro}).
#'   Only used if `xout=TRUE`.
#' @param iter Number of times to repeat the LASSO estimation (default: 20).
#' @param EST Logical. If TRUE, return all `iter` estimates sorted by value for each
#'   coefficient. If FALSE (default), return only the minimum absolute values.
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#'
#' @return A list with components:
#'   \item{coef}{Vector of minimum absolute coefficient values across all iterations.
#'     Length is p+1 where the first element is the intercept. Coefficients that
#'     consistently shrink to zero will have values close to zero.}
#'   \item{all.est}{Matrix of all estimates (if EST=TRUE), otherwise NULL. Each column
#'     corresponds to a coefficient, rows are sorted estimates across iterations.}
#'
#' @details
#' LASSO estimates can vary across runs due to:
#' - Random aspects of cross-validation in glmnet
#' - Numerical precision in optimization
#'
#' This function addresses this variability by:
#' 1. Fitting the LASSO model `iter` times (default: 20)
#' 2. For each coefficient, finding the minimum absolute value across all fits
#' 3. Coefficients that are consistently shrunk to zero will have minimum values near zero
#'
#' **Interpretation:**
#' - If `min(|coef|)` is close to zero across iterations, the coefficient should likely
#'   be excluded from the model
#' - If `min(|coef|)` is substantially different from zero, the coefficient is stable
#'   and important for prediction
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' This function requires the `glmnet` package (or the package for the chosen `lasso.fun`)
#' to be installed.
#'
#' @seealso \code{\link{lasso.est}} for single LASSO fit, \code{\link{ridge.test}}
#' for ridge regression
#'
#' @examples
#' \dontrun{
#' # Generate data with some irrelevant predictors
#' set.seed(123)
#' x <- matrix(rnorm(100 * 5), ncol = 5)
#' # Only x[,1] and x[,2] are important
#' y <- x[,1] + 2*x[,2] + rnorm(100)
#'
#' # Repeated LASSO to identify stable coefficients
#' result <- lasso.rep(x, y, iter = 20)
#' result$coef  # Minimum absolute values
#'
#' # Get all estimates to see variability
#' result2 <- lasso.rep(x, y, iter = 20, EST = TRUE)
#' result2$all.est  # Each column shows sorted estimates for one coefficient
#'
#' # With leverage point removal
#' lasso.rep(x, y, xout = TRUE)
#' }
#'
#' @export
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

#' Neural Network Classification
#'
#' Performs binary classification using neural networks via the neuralnet package.
#'
#' @param train Training data matrix or data frame (optional if x1 and x2 are provided).
#' @param test Test data matrix or data frame containing observations to classify.
#' @param g Group labels for training data (required if train is provided). Should be
#'   a vector with two unique values indicating the two groups.
#' @param x1 Training data for group 1 (alternative to train/g specification).
#' @param x2 Training data for group 2 (alternative to train/g specification).
#' @param SEED Seed value for reproducibility. If NULL, no seed is set. Included for
#'   compatibility when called by other functions.
#'
#' @return A numeric vector of predicted class labels (1 or 2) for the test data.
#'
#' @details
#' This function implements neural network classification using the neuralnet package.
#' Predictions are based on a threshold of 0.5: predicted probabilities greater than
#' 0.5 are classified as group 2, otherwise group 1.
#'
#' The training data can be specified in two ways:
#' 1. Using `train` and `g`: Provide all training data in `train` with group labels in `g`
#' 2. Using `x1` and `x2`: Provide the two groups separately
#'
#' Missing values are automatically removed from training and test data.
#'
#' **Current limitation**: The function supports up to 8 independent variables (predictors).
#'
#' @note This function requires the `neuralnet` package to be installed.
#'
#' @seealso \code{\link{CLASS.fun}} for general classification dispatcher,
#' \code{\link{CLASS.BAG}} for bagged version
#'
#' @examples
#' \dontrun{
#' # Generate training data (2 predictors)
#' set.seed(123)
#' x1 <- matrix(rnorm(50), ncol = 2)
#' x2 <- matrix(rnorm(50, mean = 2), ncol = 2)
#'
#' # Generate test data
#' test <- matrix(rnorm(20, mean = 1), ncol = 2)
#'
#' # Neural network classification
#' NN.class(x1 = x1, x2 = x2, test = test)
#'
#' # Using train/g specification
#' train <- rbind(x1, x2)
#' g <- c(rep(1, nrow(x1)), rep(2, nrow(x2)))
#' NN.class(train = train, g = g, test = test)
#' }
#'
#' @export
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

#' Estimate Ridge Parameter Using Kibria's Method
#'
#' Estimates the ridge regression bias parameter k using the method recommended
#' by Kibria (2003).
#'
#' @param x Predictor matrix (n x p). Must have at least 2 predictors.
#' @param y Response vector of length n.
#' @param STAND Logical. If TRUE (default), standardize predictors and center response
#'   using robust location and scale estimators.
#' @param locfun Location estimator for standardization (default: median).
#' @param scat Scale estimator for standardization (default: madsq - median absolute
#'   deviation squared).
#' @param regfun Regression function to use for initial fit (default: ols). Can be
#'   any robust regression function.
#' @param ... Additional arguments passed to the regression function.
#'
#' @return A numeric value: the estimated ridge parameter k.
#'
#' @details
#' The ridge parameter k controls the amount of shrinkage in ridge regression. Larger
#' k values lead to more shrinkage (more bias, less variance). This function implements
#' Kibria's recommended estimator:
#'
#' k = sigma^2 / (product of alpha^2)^(1/p)
#'
#' where alpha are coefficients in the transformed space and sigma^2 is the residual
#' variance estimate.
#'
#' Missing values are automatically removed.
#'
#' @note Requires at least 2 predictors (p ≥ 2).
#'
#' @references
#' Kibria, B. M. G. (2003). Performance of some new ridge regression estimators.
#' Communications in Statistics-Simulation and Computation, 32(2), 419-435.
#'
#' @seealso \code{\link{ridge.test}} for ridge regression tests,
#' \code{\link{ridge.Gtest}} for global ridge tests
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(100 * 3), ncol = 3)
#' y <- x[,1] + 2*x[,2] + rnorm(100)
#'
#' # Estimate ridge parameter
#' k <- ridge.est.k(x, y)
#' k
#' }
#'
#' @keywords internal
#' @export
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

#' Critical Value for Global Ridge Regression Test
#'
#' Determines the critical value for testing the global null hypothesis (all slopes = 0)
#' in ridge regression via simulation.
#'
#' @param n Sample size.
#' @param p Number of predictors.
#' @param regfun Regression function to use (default: MMreg - M-estimator).
#' @param iter Number of simulation iterations (default: 5000).
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#' @param MSF Logical. If TRUE (default), use MSE-based standard errors.
#'
#' @return The 95th percentile of the simulated null distribution (critical value).
#'
#' @details
#' This function simulates the null distribution of the global test statistic under
#' the assumption that all slopes are zero. It generates data from the null model,
#' fits ridge regression, and computes the test statistic `iter` times, then returns
#' the 95th percentile as the critical value for a 0.05 level test.
#'
#' @seealso \code{\link{ridge.Gtest}} for global ridge test, \code{\link{ridgeGnullMC}}
#' for parallel version
#'
#' @examples
#' \dontrun{
#' # Get critical value for n=50, p=3
#' crit <- ridgeGnull(n = 50, p = 3, iter = 1000)
#' }
#'
#' @keywords internal
#' @export
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

#' Critical Value for Global Ridge Test (Parallel Version)
#'
#' Parallel version of \code{\link{ridgeGnull}} using multicore processing.
#'
#' @param n Sample size.
#' @param p Number of predictors.
#' @param regfun Regression function (default: MMreg).
#' @param iter Number of iterations (default: 5000).
#' @param SEED Logical. Set random seed (default: TRUE).
#'
#' @return The 95th percentile critical value.
#'
#' @seealso \code{\link{ridgeGnull}} for sequential version
#'
#' @keywords internal
#' @export
ridgeGnullMC<-function(n,p,regfun=MMreg,iter=5000,SEED=TRUE){
#
# Determine null distribution of ridge.Gtest
#
if(SEED)set.seed(45)
fv=NA
p1=p+1
a=list()
for(i in 1:iter)a[[i]]=rmul(n,p1)
fv=mclapply(a,ridgeGnullMC.sub,p=p,regfun=regfun)
fv=matl(fv)
fv=as.vector(fv)
fv
}

# ============================================================================
# ridgeG.sub
# ============================================================================

#' Bootstrap Helper for Global Ridge Test
#'
#' Helper function for \code{\link{ridge.Gtest}} that computes the test statistic
#' for one bootstrap sample.
#'
#' @param x Predictor matrix.
#' @param y Response vector.
#' @param k Ridge parameter.
#' @param regfun Regression function (default: tsreg).
#' @param outfun Outlier detection function (default: outpro).
#' @param STAND Logical. Standardize predictors (default: FALSE).
#' @param MSF Logical. Use MSE for standard errors (default: TRUE).
#' @param locfun Location estimator (default: mean).
#' @param scat Scale estimator (default: var).
#' @param ... Additional arguments.
#'
#' @return Test statistic value.
#'
#' @seealso \code{\link{ridge.Gtest}}
#'
#' @keywords internal
#' @export
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

#' Global Test for Ridge Regression
#'
#' Tests the global null hypothesis that all slope parameters are zero using
#' ridge regression with a robust analog.
#'
#' @param x Predictor matrix (n x p).
#' @param y Response vector.
#' @param k Ridge parameter. If NULL, estimated via \code{\link{ridge.est.k}}.
#' @param regfun Regression function (default: tsreg - Theil-Sen).
#' @param xout Logical. Remove leverage points (default: FALSE).
#' @param outfun Outlier detection function (default: outpro).
#' @param STAND Logical. Standardize predictors (default: FALSE).
#' @param PV Logical. If TRUE, compute p-value via simulation (slower). If FALSE,
#'   test at 0.05 level using approximate critical value (default: FALSE).
#' @param iter Number of iterations for p-value computation (default: 5000).
#' @param locfun Location estimator (default: mean).
#' @param scat Scale estimator (default: var).
#' @param MC Logical. Use parallel processing (default: FALSE).
#' @param MSF Logical. Use MSE for standard errors (default: TRUE).
#' @param ... Additional arguments.
#'
#' @return A list with test statistic, critical value, and decision.
#'
#' @details
#' Tests H0: all slopes = 0 using ridge regression. Can use simulation to compute
#' exact p-values or use an approximate critical value for faster computation.
#'
#' @seealso \code{\link{ridge.test}} for individual slope tests,
#' \code{\link{ridgeGnull}} for critical value computation
#'
#' @export
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

#' Liu Ridge Estimator
#'
#' Computes the Liu estimator, a modification of ridge regression that can have
#' smaller MSE than standard ridge regression.
#'
#' @param x Predictor matrix (n x p).
#' @param y Response vector.
#' @param Regfun Robust regression function (default: ltsreg - Least Trimmed Squares).
#' @param xout Logical. Remove leverage points (default: FALSE).
#' @param outfun Outlier detection function (default: outpro).
#' @param plotit Logical. Create diagnostic plot (default: FALSE).
#' @param d Liu parameter. If NULL (default), estimated using Alheety & Kibria method.
#' @param tr Trimming proportion (default: 0.2).
#' @param STAND Logical. Standardize predictors (default: TRUE).
#' @param INT Logical. Include intercept in estimation (default: FALSE).
#' @param locfun Location estimator (default: median).
#' @param ... Additional arguments passed to regression function.
#'
#' @return A list with components:
#'   \item{coef}{Vector of Liu estimator coefficients.}
#'   \item{d}{Liu parameter d used in estimation.}
#'
#' @details
#' The Liu estimator is a biased estimator that modifies ridge regression:
#'
#' Liu.est = (X'X + I)^(-1) (X'X + dI) beta.hat
#'
#' where d is the Liu parameter. When properly chosen, this can have smaller MSE
#' than both OLS and ridge regression.
#'
#' This implementation uses robust regression (LTS by default) for the initial
#' estimate and Alheety & Kibria's method for estimating d.
#'
#' Missing values are automatically removed.
#'
#' @references
#' Kan, B., Alpu, O., & Yazici, B. (2013). Robust ridge and robust Liu estimator
#' for regression based on the LTS estimator. Journal of Applied Statistics, 40(3), 644-655.
#'
#' @seealso \code{\link{ridge.test}} for ridge regression tests,
#' \code{\link{ltsreg}} for LTS regression
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(100 * 3), ncol = 3)
#' y <- x[,1] + 2*x[,2] + rnorm(100)
#'
#' # Liu estimator
#' ridge.Liu(x, y)
#'
#' # Specify d parameter
#' ridge.Liu(x, y, d = 0.5)
#' }
#'
#' @export
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

#' Ridge Regression Hypothesis Tests with Heteroscedastic Standard Errors
#'
#' Tests hypotheses about individual regression slopes using ridge regression with
#' heteroscedastic-consistent standard errors. Controls family-wise error rate when
#' testing multiple slopes.
#'
#' @param x Predictor matrix (n x p). Can be a matrix or data frame.
#' @param y Response vector of length n.
#' @param k Ridge parameter. If NULL (default), estimated using Kibria's method via
#'   \code{\link{ridge.est.k}}.
#' @param alpha Significance level for hypothesis tests (default: 0.05).
#' @param pr Logical. If TRUE (default), print results to console.
#' @param xout Logical. If TRUE, eliminate leverage points using `outfun` before
#'   fitting (default: FALSE).
#' @param outfun Function for detecting outliers/leverage points (default: \code{outpro}).
#' @param STAND Logical. If TRUE (default), standardize predictors and center response.
#' @param method Method for p-value adjustment: 'hoch' (Hochberg, default), 'bh'
#'   (Benjamini-Hochberg), or 'rom' (Rom's method).
#' @param locfun Location estimator for standardization (default: mean).
#' @param scat Scale estimator for standardization (default: var).
#' @param MSF Logical. If TRUE (default), use MSE to estimate standard errors. If FALSE,
#'   use residual-based method.
#' @param ... Additional arguments passed to the outlier detection function.
#'
#' @return A list with components:
#'   \item{coef}{Vector of ridge regression coefficient estimates.}
#'   \item{TEST}{Test statistics for each slope.}
#'   \item{p.value}{Unadjusted p-values for each slope.}
#'   \item{p.crit}{Critical p-values after adjustment for multiple testing.}
#'   \item{num.sig}{Number of slopes significantly different from zero.}
#'   \item{k.used}{Ridge parameter value used in the analysis.}
#'
#' @details
#' This function performs ridge regression and tests H0: beta_j = 0 for each predictor j.
#'
#' **Ridge regression** adds an L2 penalty to handle multicollinearity, but introduces
#' bias. This function uses heteroscedastic-consistent standard errors (Long & Ervin, 2000)
#' to compute test statistics and p-values.
#'
#' **Multiple testing procedure:**
#' 1. Compute p-value for each slope
#' 2. Adjust p-values using the specified method to control family-wise error rate
#' 3. If smallest p-value ≤ alpha, reject that hypothesis AND the global null (all slopes = 0)
#' 4. Due to ridge bias, no other slopes can be declared significant
#'
#' **Advantages:**
#' - Power often at least as high as OLS, potentially much higher
#' - Handles multicollinearity better than OLS
#'
#' **Limitations:**
#' - Confidence intervals can be highly inaccurate when null is false
#' - Only the most significant slope can be formally tested
#'
#' Missing values are automatically removed.
#'
#' @note This function requires at least 2 predictors (p ≥ 2).
#'
#' @references
#' Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent standard
#' errors in the linear regression model. The American Statistician, 54(3), 217-224.
#'
#' @seealso \code{\link{ridge.est.k}} for estimating ridge parameter,
#' \code{\link{ridge.Gtest}} for global test variant, \code{\link{lasso.est}}
#' for LASSO regression
#'
#' @examples
#' \dontrun{
#' # Generate data with multicollinearity
#' set.seed(123)
#' x1 <- rnorm(100)
#' x2 <- x1 + rnorm(100, sd = 0.3)
#' x3 <- rnorm(100)
#' x <- cbind(x1, x2, x3)
#' y <- x1 + 2*x2 + rnorm(100)
#'
#' # Ridge regression tests
#' ridge.test(x, y)
#'
#' # Specify ridge parameter
#' ridge.test(x, y, k = 0.5)
#'
#' # With outlier removal
#' ridge.test(x, y, xout = TRUE)
#' }
#'
#' @export
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

