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

#' Compute Bootstrap Test Statistic for Dependent Groups (Internal)
#'
#' Internal helper function that computes the test statistic for comparing trimmed
#' means of two dependent groups using a bootstrap sample. Used by bootstrap
#' procedures like \code{ydbt}.
#'
#' @param isub Vector of bootstrap sample indices (integers from 1 to n)
#' @param x Numeric vector for first group
#' @param y Numeric vector for second group (same length as x)
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming)
#'
#' @return Test statistic from \code{\link{yuend}} for the bootstrap sample
#'
#' @details
#' This function is used internally by bootstrap procedures to compute the test
#' statistic for each bootstrap sample. It:
#' \itemize{
#'   \item Takes a bootstrap sample using indices in \code{isub}
#'   \item Applies \code{\link{yuend}} to compare trimmed means of the sampled
#'     paired data
#'   \item Returns the test statistic
#' }
#'
#' The function is designed to be called repeatedly in bootstrap loops to generate
#' the null distribution of the test statistic.
#'
#' @keywords internal
#' @seealso \code{\link{yuend}} for the underlying test,
#'   \code{ydbt} for bootstrap procedure that uses this
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

#' Generate Contrast Matrix for All Pairwise Comparisons (One-Way Design)
#'
#' Creates a contrast coefficient matrix for performing all pairwise comparisons
#' in a one-way ANOVA design. Each column represents one pairwise comparison.
#'
#' @param J Number of groups to compare
#'
#' @return A J × C matrix where C = J(J-1)/2 is the number of pairwise comparisons.
#'   Each column contains the contrast coefficients for one pairwise comparison,
#'   with 1 for one group, -1 for the compared group, and 0 for all others.
#'
#' @details
#' This function generates contrast coefficients for all possible pairwise
#' comparisons among J groups. The resulting matrix has J rows (one per group)
#' and J(J-1)/2 columns (one per pairwise comparison).
#'
#' For example, with J=3 groups, there are 3 pairwise comparisons:
#' \itemize{
#'   \item Group 1 vs Group 2: (1, -1, 0)
#'   \item Group 1 vs Group 3: (1, 0, -1)
#'   \item Group 2 vs Group 3: (0, 1, -1)
#' }
#'
#' The contrast matrix is used with functions like \code{\link{linconb}},
#' \code{\link{linconpb}}, and other MCP procedures.
#'
#' @examples
#' # Generate contrasts for 4 groups
#' con <- con1way(4)
#' # Results in 4×6 matrix with 6 pairwise comparisons
#'
#' # Use with linear contrast tests
#' # data_list <- list(group1_data, group2_data, group3_data, group4_data)
#' # linconb(data_list, con=con)
#'
#' @family multiple comparison procedures
#' @family contrast generators
#' @export
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


#' Linear Contrasts for M-estimators with Bootstrap
#'
#' Computes confidence intervals for linear contrasts involving M-estimators
#' (robust measures of location) using bootstrap methodology. Tests multiple
#' contrasts among independent groups while controlling family-wise error (FWE).
#'
#' @inheritParams linconb
#' @param est M-estimator function to use (default: \code{onestep}, one-step M-estimator)
#' @param nboot Number of bootstrap samples (default: 500)
#' @param ... Additional arguments passed to the estimator function
#'
#' @return A list with components:
#'   \item{psihat}{Matrix with columns: con.num, psihat (contrast estimate),
#'     ci.lower, ci.upper, se (standard error), p.value}
#'   \item{crit}{Critical value used for simultaneous confidence intervals}
#'   \item{con}{The contrast matrix used}
#'
#' @details
#' This function performs simultaneous inference for multiple linear contrasts
#' of M-estimators using a percentile bootstrap method. The confidence intervals
#' are adjusted to control the family-wise error rate across all contrasts.
#'
#' M-estimators are robust alternatives to the mean that downweight outliers.
#' The default is the one-step M-estimator based on Huber's Psi function.
#'
#' \strong{Contrast Specification:}
#' If \code{con} is not specified, all pairwise comparisons are performed.
#' Otherwise, \code{con} should be a J × d matrix where J is the number of
#' groups and d is the number of contrasts. Each column specifies one contrast.
#'
#' The bootstrap procedure resamples centered data from each group and computes
#' the maximum absolute test statistic across all contrasts to determine the
#' critical value that controls FWE.
#'
#' @note
#' Confidence intervals are adjusted to control FWE, but p-values are not
#' adjusted. Use \code{p.adjust()} for adjusted p-values if needed.
#'
#' @examples
#' # Compare 4 groups using one-step M-estimator
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1), rnorm(20, 1.5))
#'
#' # All pairwise comparisons
#' result <- linconm(x)
#'
#' # Custom contrasts
#' con <- matrix(c(1, -1, 0, 0,   # Group 1 vs 2
#'                 1, 0, -1, 0),  # Group 1 vs 3
#'               nrow=4, ncol=2)
#' result2 <- linconm(x, con=con, nboot=1000)
#'
#' @family multiple comparison procedures
#' @family robust M-estimator methods
#' @seealso \code{\link{linconb}} for trimmed means, \code{\link{linconpb}} for percentile bootstrap
#' @export
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


#' Linear Contrasts for Independent Groups (Old Version - Deprecated)
#'
#' @description
#' \strong{DEPRECATED:} This is an older version of \code{\link{lincon}}.
#' Use \code{\link{lincon}} instead for current functionality.
#'
#' Performs heteroscedastic tests of linear contrasts using trimmed means.
#' This older version uses different critical value tables than the current version.
#'
#' @inheritParams lincon
#' @param crit Critical value (default: NA, computed automatically)
#' @param KB Logical; use Kaiser-Bowden method. Now deprecated - use \code{kbcon} instead
#'
#' @return List with components similar to \code{\link{lincon}}
#'
#' @details
#' This function is retained for backward compatibility but is deprecated.
#' It uses older critical value tables and methods.
#'
#' Key differences from \code{\link{lincon}}:
#' \itemize{
#'   \item Uses different critical value lookup tables
#'   \item Slower for C > 28 comparisons
#'   \item Kaiser-Bowden method option (now handled by \code{kbcon})
#' }
#'
#' @note
#' \strong{Deprecated:} Use \code{\link{lincon}} for current analyses.
#'
#' @seealso \code{\link{lincon}} for the current version,
#'   \code{\link{kbcon}} for Kaiser-Bowden method
#' @keywords internal
#' @export
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


#' Linear Contrasts with Pooling Option (Deprecated)
#'
#' @description
#' \strong{DEPRECATED:} This function extends \code{\link{lincon.old}} with
#' a pooling option for factorial designs. Use modern factorial MCP functions instead.
#'
#' Tests linear contrasts with an option to pool data across factor levels,
#' primarily used for main effects in two-way and three-way designs.
#'
#' @inheritParams lincon
#' @param POOL Logical; if TRUE, prints pooling information for factorial designs
#'
#' @return List with results similar to \code{\link{lincon}}
#'
#' @details
#' This is an older utility function that adds pooling functionality to
#' \code{lincon.old}. When POOL=TRUE, it prints information about how data
#' is pooled for analyzing main effects in factorial designs.
#'
#' @note
#' \strong{Deprecated:} Use modern factorial MCP functions like \code{\link{mcp2a}},
#' \code{\link{mcp3atm}}, or related functions for factorial designs.
#'
#' @seealso \code{\link{lincon}} for current version,
#'   \code{\link{mcp2a}} for two-way factorial MCP,
#'   \code{\link{mcp3atm}} for three-way factorial MCP
#' @keywords internal
#' @export
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


#' Linear Contrasts for Trimmed Means with Bootstrap-t
#'
#' Computes confidence intervals and tests for multiple linear contrasts involving
#' trimmed means using the bootstrap-t method. Provides simultaneous inference
#' with family-wise error rate (FWE) control for independent groups.
#'
#' @param x Data in list mode where \code{x[[1]]} contains data for group 1, etc.
#'   Can also be a matrix or data frame (converted to list mode)
#' @inheritParams ancova
#' @param con Contrast matrix (J × d) where J = number of groups, d = number of contrasts.
#'   Each column specifies one contrast. If 0 or unspecified, all pairwise
#'   comparisons are performed. Default: 0
#' @inheritParams yuen
#' @param nboot Number of bootstrap samples (default: 599)
#' @param pr Logical; if TRUE, prints progress messages (default: FALSE)
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#' @param method P-value adjustment method passed to \code{p.adjust()} (default: 'holm').
#'   Options include 'holm', 'hochberg', 'bonferroni', 'BH', 'BY', 'fdr'
#'
#' @return A list with components:
#'   \item{n}{Sample sizes for each group}
#'   \item{psihat}{Matrix with columns: con.num, psihat (contrast estimate),
#'     ci.lower, ci.upper}
#'   \item{test}{Matrix with columns: con.num, test (test statistic), se, p.value, p.adjusted}
#'   \item{crit}{Critical value for simultaneous CIs controlling FWE}
#'   \item{con}{The contrast matrix used}
#'
#' @details
#' This function performs simultaneous inference for multiple linear contrasts of
#' trimmed means using the bootstrap-t method. The bootstrap-t approach tends to
#' provide better control of Type I error rates than percentile bootstrap methods,
#' especially for small to moderate sample sizes.
#'
#' \strong{Bootstrap-t Method:}
#' The bootstrap-t uses a studentized bootstrap where the test statistic is
#' standardized by its bootstrap standard error. This provides better accuracy
#' than raw bootstrap percentiles, particularly when the sampling distribution
#' is skewed or has heavy tails.
#'
#' \strong{Contrast Specification:}
#' Contrasts are linear combinations of group means specified by coefficients.
#' For example, with 4 groups:
#' \itemize{
#'   \item Simple contrast: (1, -1, 0, 0) compares groups 1 and 2
#'   \item Complex contrast: (1, 1, -1, -1) compares average of groups 1-2 with 3-4
#' }
#'
#' If \code{con=0}, all J(J-1)/2 pairwise comparisons are automatically generated.
#'
#' \strong{Multiple Comparison Adjustment:}
#' \itemize{
#'   \item Confidence intervals are adjusted via critical value to control FWE
#'   \item P-values are adjusted using the method specified in \code{method}
#' }
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' For comparing medians (tr=0.5), use \code{\link{medpb}} instead.
#' For M-estimators, see \code{\link{linconm}}.
#' For percentile bootstrap, see \code{\link{linconpb}}.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'   Academic Press. Chapter 7.
#'
#' @examples
#' # Compare 4 groups with 20% trimmed means
#' set.seed(123)
#' g1 <- rnorm(20)
#' g2 <- rnorm(20, mean=0.5)
#' g3 <- rnorm(20, mean=1)
#' g4 <- rnorm(20, mean=1.5)
#' x <- list(g1, g2, g3, g4)
#'
#' # All pairwise comparisons
#' result <- linconb(x, tr=0.2, nboot=599)
#'
#' # Custom contrasts
#' con <- matrix(c(1, -1, 0, 0,    # Group 1 vs 2
#'                 1, 0, -1, 0,    # Group 1 vs 3
#'                 0, 1, 0, -1),   # Group 2 vs 4
#'               nrow=4, ncol=3)
#' result2 <- linconb(x, con=con, tr=0.1, nboot=1000)
#'
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @family trimmed mean methods
#' @seealso \code{\link{linconpb}}, \code{\link{linconbt}}, \code{\link{linconm}}
#' @export
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


#' Linear Contrasts for Trimmed Means with Bootstrap-t (Alternative Implementation)
#'
#' An alternative implementation of bootstrap-t linear contrasts for trimmed means.
#' Similar to \code{\link{linconb}} but uses a slightly different computational
#' approach. Provides simultaneous inference with FWE control for independent groups.
#'
#' @inheritParams linconb
#'
#' @return A list with components:
#'   \item{n}{Sample sizes for each group}
#'   \item{psihat}{Matrix with columns: con.num, psihat (contrast estimate),
#'     ci.lower, ci.upper}
#'   \item{test}{Matrix with columns: con.num, test (test statistic), se, p.value, p.adjusted}
#'   \item{crit}{Critical value for simultaneous CIs controlling FWE}
#'   \item{con}{The contrast matrix used}
#'
#' @details
#' This function is an alternative implementation of the bootstrap-t method for
#' linear contrasts. It uses \code{con.all.pairs()} to generate all pairwise
#' contrasts when \code{con=0}, which may produce slightly different ordering
#' compared to \code{\link{linconb}}.
#'
#' The bootstrap-t method provides studentized bootstrap confidence intervals,
#' which tend to have better coverage properties than percentile bootstrap methods,
#' especially for small to moderate sample sizes.
#'
#' \strong{Differences from linconb:}
#' \itemize{
#'   \item Uses \code{con.all.pairs()} for generating pairwise contrasts
#'   \item May have slightly different computational details
#'   \item Generally produces similar results to \code{linconb}
#' }
#'
#' In most cases, \code{\link{linconb}} and \code{linconbt} will give very similar
#' results. The choice between them is largely a matter of preference or specific
#' computational considerations.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' For most applications, \code{\link{linconb}} is recommended as the primary
#' function. This function is retained for compatibility and as an alternative
#' implementation.
#'
#' @examples
#' # Compare 3 groups with 20% trimmed means
#' set.seed(123)
#' x <- list(rnorm(15), rnorm(15, 0.5), rnorm(15, 1))
#'
#' # All pairwise comparisons
#' result <- linconbt(x, tr=0.2, nboot=599)
#'
#' # Custom contrast: Group 1 vs average of groups 2 and 3
#' con <- matrix(c(2, -1, -1), nrow=3, ncol=1)
#' result2 <- linconbt(x, con=con, tr=0.1)
#'
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @family trimmed mean methods
#' @seealso \code{\link{linconb}}, \code{\link{linconpb}}
#' @export
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


#' Linear Contrasts with Percentile Bootstrap and Rom's Method
#'
#' Performs multiple comparisons for independent groups using percentile bootstrap
#' with Rom's sequentially rejective method for controlling family-wise error rate.
#' Provides a flexible framework for robust comparisons using any location estimator.
#'
#' @param x Data in list mode (each element is one group) or a matrix/data frame
#'   (columns are groups). Can contain missing values.
#' @inheritParams ancova
#' @param nboot Number of bootstrap samples. If NA, automatically determined:
#'   5000 for J>8 groups, 4000 for J=4-8, 2000 for J≤3. Default: NA
#' @param grp Optional vector specifying subset of groups to analyze.
#'   Example: \code{grp=c(1,3,5)} compares only groups 1, 3, and 5. Default: NA
#' @param est Measure of location function (default: \code{tmean}, trimmed mean).
#'   Can be any location estimator function
#' @param con Contrast matrix (J × d). If 0 or unspecified, all pairwise
#'   comparisons are performed. Default: 0
#' @param method P-value adjustment method for Rom's procedure. Options: 'holm',
#'   'hochberg', 'BH', etc. Default: 'holm'
#' @param bhop Logical; if TRUE, uses Benjamini-Hochberg procedure (same as method='BH').
#'   Default: FALSE
#' @param SEED Logical; if TRUE, sets random seed for reproducibility. Default: TRUE
#' @param ... Additional arguments passed to the estimator function \code{est}
#'
#' @return A matrix with columns:
#'   \item{con.num}{Contrast number}
#'   \item{psihat}{Estimated contrast value}
#'   \item{p.value}{Bootstrap p-value}
#'   \item{p.crit}{Rom's critical p-value for this contrast}
#'   \item{ci.lower}{Lower confidence limit}
#'   \item{ci.upper}{Upper confidence limit}
#'   \item{p.adjusted}{Adjusted p-value using specified method}
#'
#' @details
#' This function uses percentile bootstrap to perform multiple comparisons while
#' controlling the family-wise error rate via Rom's sequentially rejective method.
#' Rom's method is more powerful than Bonferroni while still strongly controlling FWE.
#'
#' \strong{Percentile Bootstrap:}
#' The percentile bootstrap resamples from each group independently and computes
#' the bootstrap distribution of each contrast. P-values are determined from the
#' proportion of bootstrap samples where the contrast has the opposite sign.
#'
#' \strong{Rom's Sequentially Rejective Method:}
#' Rom's method provides adjusted critical values that are less conservative than
#' Bonferroni. Contrasts are ordered by p-value and tested sequentially using
#' adjusted alpha levels that depend on the number of remaining tests.
#'
#' \strong{Contrast Specification:}
#' If \code{con=0}, all J(J-1)/2 pairwise comparisons are automatically generated.
#' Custom contrasts can be specified as a J × d matrix where each column defines
#' one contrast.
#'
#' \strong{Location Estimator:}
#' The default estimator is the trimmed mean (\code{tmean}), but any robust
#' location estimator can be used (e.g., \code{median}, \code{onestep}, \code{mom}).
#' Additional arguments for the estimator (e.g., \code{tr} for trim proportion)
#' can be passed via \code{...}.
#'
#' Missing values are automatically removed from each group before analysis.
#'
#' @references
#' Rom, D.M. (1990). A sequentially rejective test procedure based on a modified
#'   Bonferroni inequality. \emph{Biometrika}, \strong{77}, 663-665.
#'
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing (5th ed.).
#'   Academic Press. Chapter 7.
#'
#' @examples
#' # Compare 4 groups using trimmed means
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1), rnorm(20, 1.5))
#'
#' # All pairwise comparisons with 20% trimming
#' result <- linconpb(x, est=tmean, tr=0.2, nboot=2000)
#'
#' # Compare only groups 1, 2, and 4
#' result2 <- linconpb(x, grp=c(1,2,4), est=tmean, tr=0.2)
#'
#' # Use median instead of trimmed mean
#' result3 <- linconpb(x, est=median, nboot=2000)
#'
#' # Custom contrasts
#' con <- matrix(c(1, -1, 0, 0,    # Group 1 vs 2
#'                 1, 0, 0, -1),   # Group 1 vs 4
#'               nrow=4, ncol=2)
#' result4 <- linconpb(x, con=con, est=tmean, tr=0.1)
#'
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @family trimmed mean methods
#' @seealso \code{\link{linconb}} for bootstrap-t method, \code{\link{tmcppb}} for alternative implementation
#' @export
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


#' Multiple Comparisons for Independent Groups with Multivariate Data
#'
#' Performs multiple comparisons for J independent groups with multivariate data using
#' percentile bootstrap with Rom's method. Tests linear contrasts using marginal
#' estimators (default: trimmed means) and Mahalanobis or projection distances.
#'
#' @param x Data in list mode or matrix. If list, \code{x[[j]]} contains the multivariate
#'   data for group j (matrix form). If matrix, converted via \code{MAT2list}.
#' @param alpha Family-wise Type I error rate (default: 0.05). Rom's method controls FWER.
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param grp Vector specifying subset of groups to compare (default: NA uses all groups).
#'   Example: \code{grp=c(1,3,5)} compares groups 1, 3, and 5.
#' @param est Measure of location for marginal distributions (default: \code{tmean}).
#' @param con Contrast matrix (default: 0 generates all pairwise comparisons).
#'   Rows correspond to groups, columns to contrasts.
#' @param bhop Logical. If TRUE, use Benjamini-Hochberg method instead of Rom's (default: FALSE).
#' @param SEED Logical. If TRUE, set random seed to 2 for reproducibility (default: TRUE).
#' @param PDIS Logical. If TRUE, use projection distances instead of Mahalanobis (default: FALSE).
#' @param J Optional: number of groups (required if x is matrix).
#' @param p Optional: number of variables (required if x is matrix).
#' @param ... Additional arguments passed to \code{est}.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: contrast number, p-value, critical p-value.}
#'   \item{con}{Contrast matrix used.}
#'   \item{num.sig}{Number of significant contrasts.}
#'
#' @details
#' This function tests linear contrasts for multivariate data by applying the estimator
#' \code{est} to each marginal distribution separately. The bootstrap method generates
#' samples, computes the contrast estimates, and uses Mahalanobis distance (or projection
#' distance if \code{PDIS=TRUE}) to compute p-values. Rom's method adjusts critical values
#' to control family-wise error rate.
#'
#' For overall multivariate location (accounting for data structure), see \code{linconSpb}
#' which uses methods like \code{smean}.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{linconSpb}}, \code{\link{linconpb}}, \code{\link{con1way}}
#' @examples
#' \dontrun{
#' # Three groups with bivariate data
#' set.seed(123)
#' x <- list(
#'   matrix(rnorm(30), ncol=2),
#'   matrix(rnorm(30, mean=0.5), ncol=2),
#'   matrix(rnorm(30), ncol=2)
#' )
#' result <- linconMpb(x, nboot=500)
#' print(result$output)
#' }
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


#' Multivariate Linear Contrasts with Skipped Estimator
#'
#' Performs multiple comparisons for J independent groups with multivariate data using
#' multivariate measures of location that account for the overall data structure.
#' Uses percentile bootstrap with Rom's method.
#'
#' @param x Data in list mode or matrix. If list, \code{x[[j]]} contains the multivariate
#'   data for group j (matrix form). If matrix, converted via \code{MAT2list}.
#' @param alpha Family-wise Type I error rate (default: 0.05). Rom's method controls FWER.
#' @param nboot Number of bootstrap samples (default: 1000).
#' @param grp Vector specifying subset of groups to compare (default: NA uses all groups).
#'   Example: \code{grp=c(1,3,5)} compares groups 1, 3, and 5.
#' @param est Multivariate measure of location (default: \code{smean}, the skipped mean
#'   based on projection outlier detection).
#' @param con Contrast matrix (default: 0 generates all pairwise comparisons).
#'   Rows correspond to groups, columns to contrasts.
#' @param bhop Logical. If TRUE, use Benjamini-Hochberg method instead of Rom's (default: FALSE).
#' @param SEED Logical. If TRUE, set random seed to 2 for reproducibility (default: TRUE).
#' @param PDIS Logical. If TRUE, use projection distances instead of Mahalanobis (default: FALSE).
#' @param J Optional: number of groups (required if x is matrix).
#' @param p Optional: number of variables (required if x is matrix).
#' @param ... Additional arguments passed to \code{est}.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: contrast number, p-value, critical p-value.}
#'   \item{con}{Contrast matrix used.}
#'   \item{num.sig}{Number of significant contrasts.}
#'
#' @details
#' Unlike \code{linconMpb} which applies estimators to marginal distributions separately,
#' this function uses multivariate measures of location that account for the overall
#' structure of the data. The default \code{est=smean} uses the skipped mean based on
#' projection method for outlier detection, which can provide more powerful tests when
#' data have multivariate outliers.
#'
#' Mahalanobis distance is used by default to compute p-values, but projection distances
#' can be used by setting \code{PDIS=TRUE}.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{linconMpb}}, \code{\link{smean}}, \code{\link{linconpb}}
#' @examples
#' \dontrun{
#' # Three groups with bivariate data
#' set.seed(123)
#' x <- list(
#'   matrix(rnorm(30), ncol=2),
#'   matrix(rnorm(30, mean=0.5), ncol=2),
#'   matrix(rnorm(30), ncol=2)
#' )
#' result <- linconSpb(x, nboot=500)
#' print(result$output)
#' }
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


#' Linear Contrasts with Explanatory Power Effect Sizes
#'
#' Computes confidence intervals and tests for linear contrasts of trimmed means
#' with explanatory power-based effect size estimates. Designed for factorial
#' designs (especially two-way) where pooling over factor levels may be appropriate.
#'
#' @param x Data in list mode, matrix, or data frame. For list mode: \code{x[[1]]}
#'   is group 1, \code{x[[2]]} is group 2, etc. Matrix/data frame: columns are groups
#' @param con Contrast matrix (J × d) where J = number of groups, d = number of contrasts.
#'   Each column specifies one contrast. If 0 (default), all pairwise comparisons performed
#' @inheritParams yuen
#' @inheritParams ancova
#' @param pr Logical; if TRUE, prints explanatory messages (default: TRUE)
#' @param crit Critical value for confidence intervals. If NA (default), computed automatically
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#' @param INT Logical; if TRUE, computes interaction effect sizes (default: FALSE)
#' @param nreps Number of replications for effect size computation (default: 200)
#' @param POOL Logical; if TRUE, pools data with matching contrast coefficients
#'   (useful for main effects in factorial designs, default: FALSE)
#'
#' @return A list with components:
#'   \item{n}{Sample sizes for each group}
#'   \item{test}{Matrix with columns: Group/con.num, Group/con.num, test, crit, se, df}
#'   \item{psihat}{Matrix with columns: Group/con.num, Group/con.num, psihat,
#'     ci.lower, ci.upper, p.value, Effect.Size (explanatory power estimate)}
#'
#' @details
#' This function performs linear contrast tests using trimmed means while computing
#' effect sizes based on explanatory power, a measure of the proportion of variance
#' in the outcome explained by group differences.
#'
#' \strong{Explanatory Power Effect Size:}
#' Explanatory power provides an alternative to standardized mean differences,
#' estimating how much of the outcome variability is attributable to group
#' membership. It is computed via the \code{\link{linEP}} function (for main effects)
#' or \code{\link{Inter.EP}} (for interactions when \code{INT=TRUE}).
#'
#' \strong{Pooling in Factorial Designs:}
#' When \code{POOL=TRUE}, groups with contrast coefficient +1 are pooled together,
#' and groups with coefficient -1 are pooled together, then the two pooled groups
#' are compared. This is appropriate for testing main effects in balanced
#' factorial designs where you want to combine data across levels of other factors.
#'
#' \strong{Contrast Specification:}
#' \itemize{
#'   \item If \code{con=0}, all J(J-1)/2 pairwise comparisons are performed
#'   \item Each column of \code{con} defines one linear combination of group means
#'   \item Confidence intervals control family-wise error rate (FWE)
#'   \item P-values are not adjusted for multiplicity
#' }
#'
#' This function is primarily used internally by \code{\link{bbmcpEP}} for
#' two-way factorial designs but can be called directly for custom analyses.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' For comparing medians (tr=0.5), use \code{\link{medpb}} instead.
#' For contrasts without effect sizes, use \code{\link{linconb}}.
#' For quantile shift effect sizes, see \code{\link{linconES}} or \code{\link{linconQS}}.
#'
#' @examples
#' # Two-way design: 2 × 3 factorial
#' set.seed(123)
#' # Generate data for 6 groups (2 levels × 3 levels)
#' x <- list(
#'   rnorm(20, 0), rnorm(20, 0.3), rnorm(20, 0.6),  # Factor A level 1
#'   rnorm(20, 0.5), rnorm(20, 0.8), rnorm(20, 1.1)  # Factor A level 2
#' )
#'
#' # All pairwise comparisons with effect sizes
#' result <- linconEP(x, tr=0.2, nreps=100)
#'
#' # Custom contrast: Factor A main effect (pooling over Factor B)
#' # Groups 1-3 vs Groups 4-6
#' con <- matrix(c(1, 1, 1, -1, -1, -1), nrow=6, ncol=1)
#' result2 <- linconEP(x, con=con, POOL=TRUE, nreps=100)
#'
#' @family multiple comparison procedures
#' @family trimmed mean methods
#' @family effect size methods
#' @seealso \code{\link{linconES}}, \code{\link{linconQS}}, \code{\link{bbmcpEP}}
#' @export
linconEP<-function(x,con=0,tr=.2,alpha=.05,pr=TRUE,crit=NA,SEED=TRUE,INT=FALSE,nreps=200,POOL=FALSE){
#
#  (Internal comments preserved for reference)
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


#' Linear Contrasts with Quantile Shift Effect Sizes
#'
#' Computes confidence intervals and tests for linear contrasts of trimmed means
#' with quantile-based effect size estimates. Provides Q-effect and relative
#' Q-effect measures that quantify distributional shifts between groups.
#'
#' @param x Data in list mode, matrix, or data frame. For list mode: \code{x[[1]]}
#'   is group 1, \code{x[[2]]} is group 2, etc. Matrix/data frame: columns are groups
#' @param con Contrast matrix (J × d) where J = number of groups, d = number of contrasts.
#'   Each column specifies one contrast. If 0 (default), all pairwise comparisons performed
#' @inheritParams yuen
#' @inheritParams ancova
#' @param pr Logical; if TRUE, prints explanatory messages (default: TRUE)
#' @param crit Critical value for confidence intervals. If NA (default), computed automatically
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#' @param INT Logical; if TRUE, computes interaction effect sizes using \code{\link{interQS}}
#'   (default: FALSE)
#' @param locfun Location function for effect size computation (default: \code{tmean})
#'
#' @return A list with components:
#'   \item{n}{Sample sizes for each group}
#'   \item{test}{Matrix with columns: Group/con.num, Group/con.num, test, crit, se, df}
#'   \item{psihat}{For pairwise: Group, Group, psihat, ci.lower, ci.upper, p.value,
#'     Q.effect, Rel.Q. For custom contrasts: con.num, psihat, ci.lower, ci.upper,
#'     p.value, Q.effect}
#'
#' @details
#' This function performs linear contrast tests using trimmed means while computing
#' quantile-based effect sizes that describe the magnitude and direction of
#' distributional shifts.
#'
#' \strong{Quantile Shift Effect Sizes:}
#' \itemize{
#'   \item \strong{Q.effect}: Proportion of group 1 scores below the median of group 2
#'     (or linear combination for custom contrasts). Q = 0.5 indicates no shift,
#'     Q > 0.5 indicates group 1 tends to have smaller values
#'   \item \strong{Rel.Q}: Relative quantile effect = (Q - 0.5) / 0.5, rescaled to
#'     range [-1, 1]. Only reported for pairwise comparisons
#' }
#'
#' Under normality and homoscedasticity, Cohen's d values of 0, 0.2, 0.5, 0.8
#' correspond approximately to Q.effect values of 0.5, 0.55, 0.65, 0.70, respectively.
#'
#' \strong{Contrast Specification:}
#' \itemize{
#'   \item If \code{con=0}, all J(J-1)/2 pairwise comparisons are performed
#'   \item Each column of \code{con} defines one linear combination of group means
#'   \item Confidence intervals control family-wise error rate (FWE)
#'   \item P-values are not adjusted; use \code{p.adjust()} if needed
#' }
#'
#' Effect sizes are computed via \code{\link{lin.ES}} for main effects or
#' \code{\link{interQS}} for interactions when \code{INT=TRUE}.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' For comparing medians (tr=0.5), use \code{\link{medpb}} instead.
#' For explanatory power effect sizes, see \code{\link{linconEP}}.
#' \code{\link{linconQS}} provides identical quantile-based effect sizes.
#'
#' @examples
#' # Compare 4 groups with quantile effect sizes
#' set.seed(123)
#' x <- list(
#'   rnorm(25, 0),
#'   rnorm(25, 0.5),
#'   rnorm(25, 1.0),
#'   rnorm(25, 1.5)
#' )
#'
#' # All pairwise comparisons
#' result <- linconES(x, tr=0.2)
#'
#' # Custom contrast: groups 1+2 vs 3+4
#' con <- matrix(c(1, 1, -1, -1), nrow=4, ncol=1)
#' result2 <- linconES(x, con=con, tr=0.1)
#'
#' @family multiple comparison procedures
#' @family trimmed mean methods
#' @family effect size methods
#' @seealso \code{\link{linconQS}}, \code{\link{linconEP}}, \code{\link{lin.ES}}
#' @export
linconES<-function(x,con=0,tr=.2,alpha=.05,pr=TRUE,crit=NA,SEED=TRUE,INT=FALSE,
locfun=tmean){
#
#  Like the function lincon, only
#  this function  estimates effect size via
#  quantile shift perspective.
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


#' Linear Contrasts with Quantile Shift Effect Sizes (Alternative Implementation)
#'
#' Computes confidence intervals and tests for linear contrasts of trimmed means
#' with quantile-based effect size estimates. Functionally identical to
#' \code{\link{linconES}}, providing Q-effect and relative Q-effect measures.
#'
#' @param x Data in list mode, matrix, or data frame. For list mode: \code{x[[1]]}
#'   is group 1, \code{x[[2]]} is group 2, etc. Matrix/data frame: columns are groups
#' @param con Contrast matrix (J × d) where J = number of groups, d = number of contrasts.
#'   Each column specifies one contrast. If 0 (default), all pairwise comparisons performed
#' @inheritParams yuen
#' @inheritParams ancova
#' @param pr Logical; if TRUE, prints explanatory messages (default: TRUE)
#' @param crit Critical value for confidence intervals. If NA (default), computed automatically
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#' @param INT Logical; if TRUE, computes interaction effect sizes using \code{\link{interQS}}
#'   (default: FALSE)
#' @param locfun Location function for effect size computation (default: \code{tmean})
#'
#' @return A list with components:
#'   \item{n}{Sample sizes for each group}
#'   \item{test}{Matrix with columns: Group/con.num, Group/con.num, test, crit, se, df}
#'   \item{psihat}{For pairwise: Group, Group, psihat, ci.lower, ci.upper, p.value,
#'     Q.effect, Rel.Q. For custom contrasts: con.num, psihat, ci.lower, ci.upper,
#'     p.value, Q.effect}
#'
#' @details
#' This function is an alternative implementation of \code{\link{linconES}}, providing
#' the same quantile shift effect size computations. The functions differ only in
#' minor implementation details but produce equivalent results.
#'
#' \strong{Quantile Shift Effect Sizes:}
#' \itemize{
#'   \item \strong{Q.effect}: Proportion of group 1 scores below the median of group 2.
#'     Q = 0.5 indicates no distributional shift
#'   \item \strong{Rel.Q}: Relative quantile effect = (Q - 0.5) / 0.5, rescaled to [-1, 1]
#' }
#'
#' \strong{Interpretation Guidelines:}
#' Under normality and homoscedasticity:
#' \itemize{
#'   \item Cohen's d = 0 ≈ Q.effect = 0.50 (no effect)
#'   \item Cohen's d = 0.2 ≈ Q.effect = 0.55 (small effect)
#'   \item Cohen's d = 0.5 ≈ Q.effect = 0.65 (medium effect)
#'   \item Cohen's d = 0.8 ≈ Q.effect = 0.70 (large effect)
#' }
#'
#' \strong{Multiple Comparison Control:}
#' \itemize{
#'   \item Confidence intervals are adjusted to control family-wise error rate (FWE)
#'   \item P-values are NOT adjusted - use \code{p.adjust()} if needed
#' }
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' This function provides identical effect sizes to \code{\link{linconES}}.
#' Use either function based on preference; results will be the same.
#' For comparing medians (tr=0.5), use \code{\link{medpb}} instead.
#'
#' @examples
#' # Compare 3 groups with quantile effect sizes
#' set.seed(123)
#' x <- list(
#'   rnorm(30, 0),
#'   rnorm(30, 0.4),
#'   rnorm(30, 0.8)
#' )
#'
#' # All pairwise comparisons
#' result <- linconQS(x, tr=0.2)
#'
#' # Custom contrast with median as location function
#' con <- matrix(c(1, -0.5, -0.5), nrow=3, ncol=1)
#' result2 <- linconQS(x, con=con, tr=0.1, locfun=median)
#'
#' @family multiple comparison procedures
#' @family trimmed mean methods
#' @family effect size methods
#' @seealso \code{\link{linconES}}, \code{\link{linconEP}}, \code{\link{lin.ES}}
#' @export
linconQS<-function(x,con=0,tr=.2,alpha=.05,pr=TRUE,crit=NA,SEED=TRUE,INT=FALSE,
locfun=tmean){
#
#  This function is used when  estimating effect size via
#  quantile shift perspective.
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


#' Linear Contrasts for Independent Binomial Proportions
#'
#' Computes confidence intervals and tests for linear contrasts of independent
#' binomial proportions. Supports multiple methods for inference on linear
#' combinations of success probabilities.
#'
#' @param r Vector of number of successes for J independent groups
#' @param n Vector of corresponding sample sizes
#' @param con Contrast coefficient vector or matrix. If NULL, performs all
#'   pairwise comparisons (default: NULL)
#' @param alpha Significance level (default: 0.05)
#' @param null.value Null hypothesis value for the contrast (default: 0)
#' @param x Optional raw binary data in list mode. If provided, r and n are
#'   computed from x (default: NULL)
#' @param method Method for computing CI: 'KMS' (Kulinskaya-Morgenthaler-Staudte,
#'   default) or other methods
#' @param binCI Function for computing binomial confidence intervals
#'   (default: \code{acbinomci})
#'
#' @return List with components:
#'   \item{psihat}{Contrast estimate(s)}
#'   \item{ci}{Confidence interval(s)}
#'   \item{p.value}{P-value(s) for testing null.value}
#'   \item{test.stat}{Test statistic(s)}
#'
#' @details
#' This function performs inference on linear contrasts of binomial proportions,
#' such as differences or weighted combinations of success probabilities across
#' J independent groups.
#'
#' \strong{Methods:}
#' \itemize{
#'   \item KMS: Kulinskaya-Morgenthaler-Staudte method (recommended)
#'   \item Other methods available via the binCI parameter
#' }
#'
#' When con=NULL, all pairwise differences are tested. Custom contrasts can
#' be specified via the con parameter.
#'
#' @examples
#' \dontrun{
#' # Compare 3 groups
#' r <- c(15, 22, 18)  # Number of successes
#' n <- c(30, 30, 30)  # Sample sizes
#'
#' # All pairwise comparisons
#' result <- lincon.bin(r, n)
#'
#' # Custom contrast: group 1 vs average of groups 2 and 3
#' con <- c(1, -0.5, -0.5)
#' result2 <- lincon.bin(r, n, con=con)
#' }
#'
#' @family multiple comparison procedures
#' @family binomial methods
#' @seealso \code{\link{binmcp}} for multiple comparisons with FWE control,
#'   \code{\link{bi2KMSv2}} for two-group comparison,
#'   \code{\link{lincon.binPV}} for p-value version
#' @export
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


#' Compute P-value for Linear Contrast of Binomial Proportions
#'
#' Computes p-value for testing a linear contrast of independent binomial
#' proportions using the Zou et al. (2009) method. Typically called by
#' \code{\link{lincon.bin}} for computing p-values.
#'
#' @param r Vector of number of successes for J independent groups
#' @param n Vector of corresponding sample sizes
#' @param con Contrast coefficient vector. If NULL, computes for all pairwise
#'   differences (default: NULL)
#' @param alpha Significance level (default: 0.05)
#' @param nullval Hypothesized value of the linear contrast (default: 0)
#' @param binCI Function for computing binomial confidence intervals
#'   (default: \code{acbinomci})
#'
#' @return List with components:
#'   \item{test.stat}{Test statistic}
#'   \item{p.value}{P-value for testing the contrast equals nullval}
#'   \item{est}{Estimated contrast value}
#'
#' @details
#' This function uses the method from Zou, Huang, & Zheng (2009) for testing
#' linear contrasts of binomial proportions. It computes confidence intervals
#' for each proportion using the specified binCI function, then combines them
#' to test the linear contrast.
#'
#' The test statistic is based on whether the null value falls within the
#' confidence interval for the linear contrast. The p-value is computed
#' accordingly.
#'
#' @note
#' This is typically called internally by \code{\link{lincon.bin}}. Users
#' should generally use \code{\link{lincon.bin}} directly.
#'
#' @references
#' Zou, G. Y., Huang, W., & Zheng, X. (2009). A note on confidence interval
#' estimation for a linear function of binomial proportions. Computational
#' Statistics & Data Analysis, 53, 1080-1085.
#'
#' @examples
#' \dontrun{
#' # Compute p-value for difference between two groups
#' r <- c(15, 22)
#' n <- c(30, 30)
#' con <- c(1, -1)  # Group 1 - Group 2
#'
#' result <- lincon.binPV(r, n, con=con, nullval=0)
#' print(result$p.value)
#' }
#'
#' @family binomial methods
#' @seealso \code{\link{lincon.bin}} for the main function,
#'   \code{\link{binmcp}} for multiple comparisons
#' @keywords internal
#' @export
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

#' Linear Contrasts for Dependent Groups
#'
#' Computes test statistics for linear contrasts involving trimmed means from
#' dependent (repeated measures) groups using a supplied covariance matrix.
#' Provides basic infrastructure for testing custom contrasts in dependent designs.
#'
#' @param x Data in list mode (each element is one repeated measure) or a matrix
#'   (columns are repeated measures). Missing values handled via \code{na.rm=TRUE}
#' @param con Contrast matrix (J × d) where J = number of groups/measures,
#'   d = number of contrasts. Each column specifies one contrast. Required parameter.
#' @param cmat Covariance matrix (J × J) for the contrasts. Typically obtained
#'   from winsorized covariance or other robust covariance estimator
#' @inheritParams ancova
#' @inheritParams yuen
#'
#' @return A list with component:
#'   \item{test.stat}{Matrix with columns: con.num (contrast number),
#'     psihat (estimated contrast), se (standard error), test (test statistic)}
#'
#' @details
#' This function provides the computational core for testing linear contrasts
#' in dependent (repeated measures) designs. Unlike independent group contrasts,
#' dependent contrasts must account for correlations between measurements.
#'
#' \strong{Covariance Matrix:}
#' The \code{cmat} parameter should contain the covariance matrix for the
#' measurements. For robust methods, this is typically a winsorized covariance
#' matrix obtained from \code{wincov()} or similar functions.
#'
#' \strong{Test Statistic:}
#' For each contrast with coefficients c, the test statistic is:
#' \deqn{T = \hat{\psi} / \sqrt{c' \Sigma c}}
#' where \eqn{\hat{\psi}} is the estimated contrast and \eqn{\Sigma} is the
#' covariance matrix.
#'
#' \strong{Contrast Specification:}
#' Contrasts must sum to zero for proper interpretation in repeated measures
#' designs. For example:
#' \itemize{
#'   \item Time 1 vs Time 2: (1, -1, 0, 0)
#'   \item Average of Times 1-2 vs Times 3-4: (0.5, 0.5, -0.5, -0.5)
#' }
#'
#' This function computes test statistics but does not provide p-values or
#' critical values. For complete inference with FWE control, use higher-level
#' functions like \code{\link{lindepbt}} or \code{\link{rmmcp}}.
#'
#' @note
#' This is a lower-level function. Most users should use \code{\link{lindepbt}}
#' for bootstrap-t inference or \code{\link{rmmcp}} for general repeated measures
#' multiple comparisons.
#'
#' @examples
#' # Four repeated measures
#' set.seed(123)
#' time1 <- rnorm(20)
#' time2 <- time1 + rnorm(20, 0.2, 0.5)
#' time3 <- time1 + rnorm(20, 0.5, 0.5)
#' time4 <- time1 + rnorm(20, 0.8, 0.5)
#' x <- cbind(time1, time2, time3, time4)
#'
#' # Define contrasts
#' con <- matrix(c(1, -1, 0, 0,      # Time 1 vs 2
#'                 0, 1, -1, 0,      # Time 2 vs 3
#'                 1, 0, 0, -1),     # Time 1 vs 4
#'               nrow=4, ncol=3)
#'
#' # Compute winsorized covariance
#' cmat <- wincov(x, tr=0.2)$cov
#'
#' # Test contrasts
#' result <- lindep(x, con=con, cmat=cmat, tr=0.2)
#'
#' @family multiple comparison procedures
#' @family dependent groups
#' @family trimmed mean methods
#' @seealso \code{\link{lindepbt}} for bootstrap-t inference, \code{\link{rmmcp}} for general MCP
#' @export
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


#' Linear Contrasts for Dependent Groups with Bootstrap-t
#'
#' Computes confidence intervals and tests for linear contrasts involving
#' trimmed means of dependent (repeated measures) groups using the percentile
#' bootstrap-t method. Controls family-wise error rate using Rom's method.
#'
#' @param x Data in matrix form (n × J, columns are groups/measures) or list mode.
#'   Missing values are automatically removed
#' @param con Contrast matrix (J × d) where J = number of groups, d = number of contrasts.
#'   If NULL (default), all pairwise comparisons are performed
#' @inheritParams yuen
#' @inheritParams ancova
#' @param nboot Number of bootstrap samples (default: 599)
#' @param dif Logical; if TRUE (default), uses difference scores for contrasts.
#'   If FALSE, tests based on marginal trimmed means with covariance structure
#' @param method P-value adjustment method for \code{p.adjust()} (default: 'holm').
#'   Options: 'holm', 'hochberg', 'bonferroni', 'BH', 'BY', 'fdr'
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#'
#' @return A list with components:
#'   \item{test}{Matrix with columns: con.num, test (test statistic), p.value,
#'     p.crit (Rom's critical value), se, p.adjusted}
#'   \item{psihat}{Matrix with columns: con.num, psihat (contrast estimate),
#'     ci.lower, ci.upper}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts}
#'
#' @details
#' This function performs simultaneous inference for multiple linear contrasts
#' involving trimmed means from dependent groups. It uses a bootstrap-t approach
#' combined with Rom's sequentially rejective method to control the family-wise
#' error rate.
#'
#' \strong{Bootstrap-t Procedure:}
#' The bootstrap-t method provides better accuracy than simple percentile bootstrap,
#' especially for small to moderate sample sizes, by studentizing the bootstrap
#' distribution.
#'
#' \strong{Rom's Method:}
#' Rom's sequentially rejective procedure is a more powerful alternative to
#' Bonferroni correction for controlling FWE. It uses sequential critical values
#' that depend on the number and ordering of p-values.
#'
#' \strong{Difference Score vs Marginal Approach:}
#' \itemize{
#'   \item \code{dif=TRUE}: Forms difference scores directly from the contrast
#'     coefficients (e.g., for contrast (1,-1,0), computes X1-X2). This is simpler
#'     and recommended for most applications
#'   \item \code{dif=FALSE}: Tests based on marginal trimmed means while accounting
#'     for the covariance structure. More complex but preserves the marginal
#'     interpretation
#' }
#'
#' \strong{Contrast Specification:}
#' If \code{con=NULL}, generates all J(J-1)/2 pairwise comparisons automatically.
#' For custom contrasts, each column of \code{con} specifies coefficients for
#' one linear combination of group means.
#'
#' The bootstrap resamples entire rows (subjects) to preserve within-subject
#' correlations across repeated measures.
#'
#' @note
#' For pairwise comparisons only, \code{\link{pairdepb}} provides a simpler interface.
#' For independent groups, use \code{\link{linconb}} or \code{\link{linconbt}}.
#'
#' @examples
#' # Four repeated measures with custom contrasts
#' set.seed(123)
#' n <- 30
#' t1 <- rnorm(n, 0)
#' t2 <- t1 + rnorm(n, 0.3, 0.4)
#' t3 <- t1 + rnorm(n, 0.5, 0.4)
#' t4 <- t1 + rnorm(n, 0.8, 0.4)
#' data <- cbind(t1, t2, t3, t4)
#'
#' # All pairwise comparisons
#' result <- lindepbt(data, tr=0.2, nboot=999)
#'
#' # Custom contrasts: linear trend and quadratic trend
#' con <- matrix(c(-3, -1, 1, 3,      # Linear trend
#'                  1, -1, -1, 1),    # Quadratic trend
#'               nrow=4, ncol=2)
#' result2 <- lindepbt(data, con=con, dif=TRUE, nboot=999)
#'
#' # Using marginal means approach
#' result3 <- lindepbt(data, con=con, dif=FALSE, nboot=599)
#'
#' @family multiple comparison procedures
#' @family dependent groups
#' @family bootstrap methods
#' @family trimmed mean methods
#' @seealso \code{\link{pairdepb}}, \code{\link{lindep}}, \code{\link{rmmcp}}
#' @export
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


#' Pairwise Comparisons for Dependent Groups with Bootstrap-t
#'
#' Computes simultaneous confidence intervals for all pairwise differences
#' between trimmed means of dependent (repeated measures) groups using the
#' percentile bootstrap-t method. Controls family-wise error rate across
#' all pairwise comparisons.
#'
#' @param x Data in matrix form (n × J, columns are groups/measures) or list mode.
#'   Missing values are NOT allowed
#' @inheritParams yuen
#' @inheritParams ancova
#' @param grp Optional vector specifying subset of groups/measures to compare.
#'   If 0 (default), all groups are compared. Example: \code{grp=c(1,3,4)}
#'   compares only columns 1, 3, and 4
#' @param nboot Number of bootstrap samples (default: 599)
#'
#' @return A list with components:
#'   \item{psihat}{Matrix with columns: Group (first group number),
#'     Group (second group number), psihat (trimmed mean difference),
#'     ci.lower, ci.upper}
#'   \item{test}{Matrix with columns: Group, Group, test (test statistic), se}
#'   \item{crit}{Critical value used for simultaneous confidence intervals}
#'
#' @details
#' This function performs all pairwise comparisons among J dependent groups
#' (repeated measures) using a bootstrap-t method. The procedure controls the
#' family-wise error rate, meaning the probability of making at least one
#' Type I error across all C = J(J-1)/2 pairwise comparisons is ≤ α.
#'
#' \strong{Bootstrap-t Procedure:}
#' \enumerate{
#'   \item Center each group's data by subtracting its trimmed mean
#'   \item For each bootstrap sample, resample rows (subjects) with replacement
#'   \item Compute test statistics for all pairwise differences
#'   \item Critical value is the (1-α) quantile of the maximum absolute
#'     test statistic across bootstrap samples
#' }
#'
#' \strong{Dependent Groups:}
#' Unlike independent groups comparisons, this function resamples entire rows
#' (subjects) to preserve the correlation structure across repeated measures.
#' This is essential for proper inference in within-subjects designs.
#'
#' \strong{Missing Values:}
#' This function does NOT allow missing values. All subjects must have complete
#' data across all measures. Use list-wise deletion before calling this function
#' if needed.
#'
#' \strong{Trimmed Means:}
#' The default 20% trimming (\code{tr=0.2}) removes the highest and lowest 20%
#' of values from each group before computing means, providing robustness against
#' outliers while maintaining good efficiency.
#'
#' @note
#' For dependent groups with missing values, consider using \code{\link{rmmcp}}
#' which handles missing data. For more flexible contrast structures beyond
#' pairwise comparisons, use \code{\link{lindepbt}}.
#'
#' @examples
#' # Four repeated measures (e.g., time points)
#' set.seed(123)
#' n <- 25
#' time1 <- rnorm(n)
#' time2 <- time1 + rnorm(n, 0.3, 0.5)
#' time3 <- time1 + rnorm(n, 0.6, 0.5)
#' time4 <- time1 + rnorm(n, 0.9, 0.5)
#' data <- cbind(time1, time2, time3, time4)
#'
#' # All pairwise comparisons
#' result <- pairdepb(data, tr=0.2, nboot=999)
#'
#' # Compare only times 1, 2, and 4
#' result2 <- pairdepb(data, grp=c(1,2,4), tr=0.1, nboot=599)
#'
#' # List mode data
#' data_list <- list(time1, time2, time3, time4)
#' result3 <- pairdepb(data_list, tr=0.2)
#'
#' @family multiple comparison procedures
#' @family dependent groups
#' @family bootstrap methods
#' @family trimmed mean methods
#' @seealso \code{\link{lindepbt}} for custom contrasts, \code{\link{rmmcp}} for handling missing data
#' @export
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

#' Multiple Comparisons with Percentile Bootstrap (General MCP)
#'
#' Computes simultaneous confidence intervals for linear contrasts involving
#' trimmed means using the percentile bootstrap method. Provides a general
#' framework for multiple comparisons among independent groups with optional
#' Winsorization.
#'
#' @param x Data in list mode (each element is one group) or matrix (columns are groups).
#'   Missing values are allowed and will be removed
#' @param crit Critical value for confidence intervals. If NA (default), values are
#'   automatically selected based on the number of contrasts, alpha level, and nboot.
#'   Must be specified when \code{tr} ≠ 0.2
#' @param con Contrast matrix (J × d). If 0 (default), all pairwise comparisons
#'   are performed. Each column specifies one contrast
#' @inheritParams yuen
#' @inheritParams ancova
#' @param nboot Number of bootstrap samples (default: 2000)
#' @param grp Optional vector specifying subset of groups to analyze.
#'   Example: \code{grp=c(1,3,5)} compares only groups 1, 3, and 5. Default: NA
#' @param WIN Logical; if TRUE, applies Winsorization before trimming for
#'   additional robustness. Default: FALSE
#' @param win Amount of Winsorization when \code{WIN=TRUE}. Must be ≤ \code{tr}.
#'   Default: 0.1
#'
#' @return A list with components:
#'   \item{psihat}{Matrix with columns: psihat (contrast estimate), ci.lower, ci.upper}
#'   \item{test}{Matrix with columns: test (test statistic), p.value, crit.value}
#'   \item{crit}{Critical value used for confidence intervals}
#'   \item{con}{The contrast matrix used}
#'   \item{n}{Sample sizes for each group}
#'
#' @details
#' This function provides a flexible percentile bootstrap framework for multiple
#' comparisons among independent groups. It can handle any set of linear contrasts
#' and optionally combines Winsorization with trimming for enhanced robustness.
#'
#' \strong{Percentile Bootstrap Method:}
#' \enumerate{
#'   \item Bootstrap samples are drawn independently from each group
#'   \item Trimmed means (and optionally Winsorized values) are computed
#'   \item The bootstrap distribution of each contrast is obtained
#'   \item Critical values are determined to control family-wise error
#' }
#'
#' \strong{Critical Value Selection:}
#' When \code{crit=NA}, the function uses pre-determined critical values for
#' standard scenarios (tr=0.2, specific combinations of d, alpha, and nboot).
#' For other scenarios, you must specify \code{crit} manually. This can be
#' obtained via simulation for the specific design.
#'
#' \strong{Winsorization Option:}
#' When \code{WIN=TRUE}, extreme values are first Winsorized (replaced with
#' less extreme values) before trimming. This provides additional robustness
#' but requires larger sample sizes (n ≥ 15 recommended). The amount of
#' Winsorization (\code{win}) should be less than or equal to the trim level
#' (\code{tr}), with \code{tr} ≥ 0.2 recommended.
#'
#' \strong{Contrast Specification:}
#' If \code{con=0}, all J(J-1)/2 pairwise comparisons are generated automatically.
#' Custom contrasts should be specified as a J × d matrix where each column
#' represents one contrast.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' This function requires careful selection of the critical value when using
#' non-standard trimming levels. For tr=0.2 with standard parameters, critical
#' values are provided automatically.
#'
#' @examples
#' # Compare 4 groups with default settings
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1), rnorm(20, 1.5))
#'
#' # All pairwise comparisons
#' result <- mcppb(x, tr=0.2, nboot=2000)
#'
#' # Custom contrasts
#' con <- matrix(c(1, -1, 0, 0,      # Group 1 vs 2
#'                 1, 0, 0, -1),     # Group 1 vs 4
#'               nrow=4, ncol=2)
#' result2 <- mcppb(x, con=con, tr=0.2, nboot=2000)
#'
#' # With Winsorization for extra robustness
#' result3 <- mcppb(x, tr=0.2, WIN=TRUE, win=0.1, nboot=2000)
#'
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @family trimmed mean methods
#' @seealso \code{\link{linconpb}} for Rom's method, \code{\link{tmcppb}} for alternative implementation
#' @export
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


#' Multiple Comparisons with Percentile Bootstrap (Trimmed Means)
#'
#' Performs all pairwise comparisons for J independent groups using a percentile
#' bootstrap method with Rom's sequentially rejective procedure to control
#' family-wise error rate. Default estimator is the trimmed mean.
#'
#' @param x Data in list mode or matrix. For list: \code{x[[1]]} is group 1,
#'   \code{x[[2]]} is group 2, etc. For matrix: columns are groups
#' @inheritParams ancova
#' @param nboot Number of bootstrap samples. If NA (default), automatically
#'   determined: 5000 for J>8, 4000 for 3<J≤8, 2000 for J≤3
#' @param grp Optional vector specifying subset of groups to compare. Example:
#'   \code{grp=c(1,3,5)} compares only groups 1, 3, and 5
#' @param est Measure of location (default: \code{tmean}, trimmed mean).
#'   Other options: \code{median}, \code{mean}, \code{mom}, etc.
#' @param con Contrast matrix (J × d). If 0 (default), all pairwise comparisons performed
#' @param bhop Logical; if TRUE, uses Benjamini-Hochberg procedure instead of Rom's
#'   method (default: FALSE)
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#' @param ... Additional arguments passed to the estimator function
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat (contrast estimate),
#'     p.value, p.crit (Rom's critical value), ci.lower, ci.upper}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts at level α}
#'
#' @details
#' This function implements a percentile bootstrap approach for multiple comparisons
#' using trimmed means (or other location estimators) with Rom's sequentially
#' rejective method for family-wise error control.
#'
#' \strong{Percentile Bootstrap:}
#' Bootstrap samples are drawn independently from each group. For each contrast,
#' the p-value is twice the proportion of bootstrap samples where the contrast
#' has the opposite sign from the observed contrast (two-tailed test).
#'
#' \strong{Rom's Method:}
#' A sequentially rejective procedure that is uniformly more powerful than
#' Bonferroni correction while still controlling FWE at level α. It uses
#' critical values that depend on the number of contrasts and their ordering
#' by p-value.
#'
#' \strong{Benjamini-Hochberg Option:}
#' Setting \code{bhop=TRUE} switches to the Benjamini-Hochberg step-up procedure,
#' which controls the false discovery rate (FDR) instead of FWE. This is more
#' powerful but allows some false positives.
#'
#' \strong{Automatic Bootstrap Samples:}
#' The default number of bootstrap samples balances accuracy and computation time
#' based on the number of groups. More groups require fewer bootstrap samples
#' per comparison due to multiple comparison adjustments.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' For M-estimators (onestep, mom), use \code{\link{pbmcp}} instead.
#' For general estimators without bootstrap, see \code{\link{linconpb}}.
#' \code{\link{bmcppb}} is functionally identical to this function.
#'
#' @examples
#' # Compare 4 groups using 20% trimmed means
#' set.seed(123)
#' x <- list(
#'   rnorm(25, 0),
#'   rnorm(25, 0.5),
#'   rnorm(25, 1.0),
#'   rnorm(25, 1.5)
#' )
#'
#' # All pairwise comparisons
#' result <- tmcppb(x, alpha=0.05, nboot=2000)
#'
#' # Compare only groups 1, 2, and 4
#' result2 <- tmcppb(x, grp=c(1,2,4), nboot=2000)
#'
#' # Using median with Benjamini-Hochberg
#' result3 <- tmcppb(x, est=median, bhop=TRUE, nboot=2000)
#'
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @family trimmed mean methods
#' @seealso \code{\link{bmcppb}}, \code{\link{pbmcp}}, \code{\link{linconpb}}, \code{\link{mcppb}}
#' @export
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


#' Multiple Comparisons with Percentile Bootstrap (Alternative Name)
#'
#' Performs all pairwise comparisons for J independent groups using a percentile
#' bootstrap method. This function is functionally identical to \code{\link{tmcppb}}.
#'
#' @inheritParams tmcppb
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat, p.value, p.crit, ci.lower, ci.upper}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts}
#'
#' @details
#' This function is an alias for \code{\link{tmcppb}}. It provides the same
#' percentile bootstrap multiple comparison procedure with Rom's method for
#' controlling family-wise error rate.
#'
#' See \code{\link{tmcppb}} for complete documentation and usage details.
#'
#' @note
#' Use either \code{\link{tmcppb}} or \code{bmcppb} based on preference - they
#' produce identical results. For M-estimators, see \code{\link{pbmcp}}.
#'
#' @examples
#' # See tmcppb for examples
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1))
#' result <- bmcppb(x, nboot=2000)
#'
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @seealso \code{\link{tmcppb}}, \code{\link{pbmcp}}
#' @export
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


#' Multiple Comparisons with Percentile Bootstrap (M-Estimators)
#'
#' Performs all pairwise comparisons for J independent groups using percentile
#' bootstrap with robust M-estimators (onestep or MOM). Uses Rom's method or
#' Benjamini-Hochberg for multiplicity control.
#'
#' @param x Data in list mode or matrix. For list: \code{x[[1]]} is group 1,
#'   \code{x[[2]]} is group 2, etc. For matrix: columns are groups
#' @inheritParams ancova
#' @param nboot Number of bootstrap samples. If NA (default), automatically
#'   determined: 5000 for J>8, 4000 for 3<J≤8, 2000 for J≤3
#' @param grp Optional vector specifying subset of groups to compare
#' @param est M-estimator to use (default: \code{onestep}, one-step M-estimator).
#'   Also accepts \code{mom} (modified one-step). For other estimators, use
#'   \code{\link{linconpb}} instead
#' @param con Contrast matrix (J × d). If 0 (default), all pairwise comparisons performed
#' @param bhop Logical; if TRUE, uses Benjamini-Hochberg instead of Rom's method
#'   (default: FALSE)
#' @param SEED Logical; if TRUE, sets random seed (default: TRUE)
#' @param ... Additional arguments passed to the estimator
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat, p.value, p.crit, ci.lower, ci.upper}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts}
#'
#' @details
#' This function specializes in robust M-estimators for multiple comparisons,
#' providing resistance to outliers while maintaining good statistical power.
#'
#' \strong{M-Estimators:}
#' \itemize{
#'   \item \strong{onestep}: One-step M-estimator based on Huber's Psi function.
#'     Robust to outliers, good statistical efficiency
#'   \item \strong{mom}: Modified one-step M-estimator. Alternative robust estimator
#'     with slightly different properties
#' }
#'
#' For other location estimators (trimmed means, median, etc.), use
#' \code{\link{linconpb}} or \code{\link{tmcppb}}.
#'
#' \strong{Bootstrap Procedure:}
#' Independent bootstrap resampling from each group. P-values computed as twice
#' the minimum of P(contrast > 0) and P(contrast < 0), accounting for ties.
#'
#' \strong{Multiplicity Control:}
#' \itemize{
#'   \item \code{bhop=FALSE}: Rom's sequentially rejective method (controls FWE)
#'   \item \code{bhop=TRUE}: Benjamini-Hochberg procedure (controls FDR)
#' }
#'
#' \strong{Sample Size Adjustments:}
#' Critical values are adjusted based on maximum sample size. For nmax > 80,
#' uses less conservative critical values appropriate for larger samples.
#'
#' If M-estimator bootstrap fails (produces NAs), the function suggests trying
#' \code{est=tmean} as an alternative.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' This function is restricted to \code{onestep} and \code{mom} estimators.
#' For trimmed means, use \code{\link{tmcppb}} or \code{\link{bmcppb}}.
#' For general estimators, use \code{\link{linconpb}}.
#'
#' @examples
#' # Compare 4 groups using one-step M-estimator
#' set.seed(123)
#' x <- list(
#'   rnorm(25, 0),
#'   rnorm(25, 0.5),
#'   rnorm(25, 1.0),
#'   rnorm(25, 1.5)
#' )
#'
#' # All pairwise comparisons
#' result <- pbmcp(x, nboot=2000)
#'
#' # Using MOM estimator
#' result2 <- pbmcp(x, est=mom, nboot=2000)
#'
#' # With Benjamini-Hochberg FDR control
#' result3 <- pbmcp(x, bhop=TRUE, nboot=2000)
#'
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @family robust M-estimator methods
#' @seealso \code{\link{tmcppb}}, \code{\link{linconpb}}, \code{\link{linconm}}
#' @export
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


#' Deprecated: Percentile Bootstrap MCP for Trimmed Means
#'
#' This function is deprecated. Use \code{\link{bmcppb}} or \code{\link{tmcppb}} instead,
#' which provide the same functionality with improved implementation.
#'
#' @param x Data in list mode or matrix. Columns correspond to groups.
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param nboot Number of bootstrap samples (default: NA, auto-determined).
#' @param grp Vector specifying subset of groups to compare (default: NA uses all groups).
#' @param con Contrast matrix (default: 0 generates all pairwise comparisons).
#' @param bhop Logical. If TRUE, use Benjamini-Hochberg method instead of Rom's (default: FALSE).
#' @param tr Trim proportion (default: 0.2 for 20% trimming).
#' @param SEED Logical. If TRUE, set random seed for reproducibility (default: TRUE).
#'
#' @return This function stops with an error message directing users to \code{bmcppb}.
#'
#' @keywords internal
#' @family MCP functions
#' @seealso \code{\link{bmcppb}}, \code{\link{tmcppb}}
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


#' Multiple Comparisons with O-Type Multivariate Location Estimator
#'
#' Performs multiple comparisons for J independent groups using multivariate O-type
#' measures of location. Uses percentile bootstrap with Rom's method.
#'
#' @inheritParams linconpb
#' @param est Multivariate measure of location (default: \code{smean}, the skipped mean).
#' @param ... Additional arguments passed to \code{est}.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: contrast number, estimate, p-value, critical p-value,
#'     confidence interval bounds.}
#'   \item{con}{Contrast matrix used.}
#'   \item{num.sig}{Number of significant contrasts.}
#'
#' @details
#' This function tests linear contrasts using multivariate O-type location estimators.
#' The default \code{est=smean} uses the skipped mean based on projection outlier detection,
#' which provides a robust multivariate location estimate that accounts for the overall
#' data structure. Rom's method controls the family-wise error rate.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{smean}}, \code{\link{linconpb}}, \code{\link{linconSpb}}
#' @examples
#' \dontrun{
#' # Three groups, univariate data
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, mean=0.5), rnorm(20))
#' result <- mcpOV(x, nboot=500)
#' print(result$output)
#' }
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

#' Multiple Comparisons for Two-Way Factorial Design (General Estimators)
#'
#' Performs all pairwise comparisons for main effects of Factor A, Factor B,
#' and all interaction contrasts in a two-way factorial design using percentile
#' bootstrap with a general measure of location (default: MOM M-estimator).
#'
#' @param J Number of levels for Factor A
#' @param K Number of levels for Factor B
#' @param x Data in list mode or matrix with J×K groups. For list mode: \code{x[[1]]}
#'   is cell (1,1), \code{x[[2]]} is cell (1,2), ..., \code{x[[K+1]]} is cell (2,1), etc.
#'   For matrix: columns correspond to cells
#' @param est Measure of location (default: \code{mom}, modified one-step M-estimator).
#'   Other options: \code{onestep}, \code{median}, \code{mean}, etc.
#' @param con Contrast matrix. If NULL (default), all pairwise comparisons are performed.
#'   If specified, use \code{\link{linconm}} instead
#' @inheritParams ancova
#' @param nboot Number of bootstrap samples. If NA (default), automatically determined:
#'   5000 for J>8, 4000 for 3<J≤8, 2000 for J≤3
#' @param grp Optional vector to reorder groups. Example: \code{grp=c(2,4,3,1)} for 2×2
#'   design reorders cells to (1,1)=x[[2]], (1,2)=x[[4]], (2,1)=x[[3]], (2,2)=x[[1]]
#' @param ... Additional arguments passed to the estimator function
#'
#' @return A list with components:
#'   \item{FactorA}{Matrix with columns: con.num, psihat, sig.test, sig.crit,
#'     ci.lower, ci.upper for Factor A main effect comparisons}
#'   \item{FactorB}{Same structure for Factor B main effects}
#'   \item{Interactions}{Same structure for all interaction contrasts}
#'   \item{conA}{Contrast matrix for Factor A}
#'   \item{conB}{Contrast matrix for Factor B}
#'   \item{conAB}{Contrast matrix for interactions}
#'
#' @details
#' This function implements a comprehensive multiple comparison procedure for
#' two-way factorial designs using robust measures of location with bootstrap
#' inference and Rom's method for FWE control.
#'
#' \strong{Data Organization:}
#' Data must be ordered by Factor B within Factor A. For a 2×3 design:
#' \itemize{
#'   \item \code{x[[1]]}: A1,B1 (first level of A, first level of B)
#'   \item \code{x[[2]]}: A1,B2
#'   \item \code{x[[3]]}: A1,B3
#'   \item \code{x[[4]]}: A2,B1
#'   \item \code{x[[5]]}: A2,B2
#'   \item \code{x[[6]]}: A2,B3
#' }
#'
#' \strong{Contrast Generation:}
#' Automatically creates contrast matrices via \code{\link{con2way}} for:
#' \itemize{
#'   \item All pairwise comparisons of Factor A marginal means
#'   \item All pairwise comparisons of Factor B marginal means
#'   \item All pairwise comparisons of cell means (interaction contrasts)
#' }
#'
#' \strong{Bootstrap Procedure:}
#' Uses percentile bootstrap with Rom's sequentially rejective method for
#' controlling family-wise error rate. P-values are computed as the proportion
#' of bootstrap samples where the contrast differs from zero.
#'
#' \strong{M-estimators:}
#' The default MOM (modified one-step M-estimator) provides robustness to
#' outliers while maintaining good efficiency. For median comparisons with
#' large sample sizes, consider \code{\link{med2mcp}} for better performance.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' For trimmed means, use \code{\link{mcp2atm}} instead.
#' To specify custom contrasts, use \code{\link{linconm}}.
#'
#' @examples
#' # 2 × 3 factorial design with MOM estimator
#' set.seed(123)
#' # Generate data: Factor A (2 levels) × Factor B (3 levels)
#' x <- list(
#'   rnorm(20, 0.0), rnorm(20, 0.3), rnorm(20, 0.6),  # A1: B1, B2, B3
#'   rnorm(20, 0.5), rnorm(20, 0.8), rnorm(20, 1.1)   # A2: B1, B2, B3
#' )
#'
#' # All main effects and interaction comparisons
#' result <- mcp2a(J=2, K=3, x=x, nboot=2000)
#'
#' # Using median instead of MOM
#' result2 <- mcp2a(J=2, K=3, x=x, est=median, nboot=2000)
#'
#' @family multiple comparison procedures
#' @family factorial design methods
#' @family robust M-estimator methods
#' @seealso \code{\link{mcp2atm}}, \code{\link{mcp3atm}}, \code{\link{con2way}}, \code{\link{med2mcp}}
#' @export
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


#' Multiple Comparisons for Two-Way Factorial Design (Trimmed Means)
#'
#' Tests all linear contrasts for main effects of Factor A, Factor B, and
#' interactions in a two-way factorial design using trimmed means with
#' heteroscedastic methods and FWE control.
#'
#' @param J Number of levels for Factor A
#' @param K Number of levels for Factor B
#' @param x Data in list mode or matrix with J×K groups. For list mode: \code{x[[1]]}
#'   is cell (1,1), \code{x[[2]]} is cell (1,2), ..., \code{x[[K+1]]} is cell (2,1), etc.
#'   For matrix: columns correspond to cells
#' @inheritParams yuen
#' @inheritParams ancova
#' @param grp Optional vector to reorder groups (default: NA, sequential order assumed)
#' @param op Logical; if TRUE, performs omnibus test across all contrasts simultaneously.
#'   If FALSE (default), tests Factor A, Factor B, and interactions separately
#' @param pr Logical; if TRUE (default), prints explanatory messages
#'
#' @return A list with components:
#'   \item{Factor.A}{Results for Factor A main effects (from \code{\link{lincon}})}
#'   \item{Factor.B}{Results for Factor B main effects}
#'   \item{Factor.AB}{Results for interaction contrasts}
#'   \item{All.Tests}{If \code{op=TRUE}, omnibus results across all contrasts}
#'   \item{conA}{Contrast matrix for Factor A}
#'   \item{conB}{Contrast matrix for Factor B}
#'   \item{conAB}{Contrast matrix for interactions}
#'
#' @details
#' This function provides comprehensive multiple comparison analysis for two-way
#' factorial designs using trimmed means. It automatically generates and tests
#' all relevant contrasts while controlling the family-wise error rate.
#'
#' \strong{Data Organization:}
#' Data should be ordered by Factor B within Factor A levels. For a 2×2 design:
#' x[[1]]=A1B1, x[[2]]=A1B2, x[[3]]=A2B1, x[[4]]=A2B2.
#'
#' \strong{Contrast Testing:}
#' Uses \code{\link{lincon}} for heteroscedastic inference on trimmed means:
#' \itemize{
#'   \item \code{op=FALSE}: Separate analyses for Factor A, Factor B, and interactions
#'   \item \code{op=TRUE}: Single omnibus test combining all contrasts
#' }
#'
#' \strong{Trimmed Means:}
#' Default 20% trimming (\code{tr=0.2}) provides robustness to outliers while
#' maintaining good statistical power. Adjust \code{tr} based on anticipated
#' outlier prevalence.
#'
#' \strong{Pooling Option:}
#' For pooling data across factor levels (useful for main effects),
#' consider \code{\link{bbmcpEP}} which includes a pooling option.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' This function is also available as \code{\link{bbmcp}} (alias).
#' For general measures of location (M-estimators), use \code{\link{mcp2a}}.
#' For three-way designs, see \code{\link{mcp3atm}}.
#'
#' @examples
#' # 2 × 3 factorial design
#' set.seed(123)
#' x <- list(
#'   rnorm(20, 0), rnorm(20, 0.3), rnorm(20, 0.6),    # A1: B1, B2, B3
#'   rnorm(20, 0.5), rnorm(20, 0.8), rnorm(20, 1.1)   # A2: B1, B2, B3
#' )
#'
#' # Separate tests for each effect
#' result <- mcp2atm(J=2, K=3, x=x, tr=0.2)
#'
#' # Omnibus test across all contrasts
#' result2 <- mcp2atm(J=2, K=3, x=x, tr=0.2, op=TRUE)
#'
#' # With 10% trimming
#' result3 <- mcp2atm(J=2, K=3, x=x, tr=0.1, pr=FALSE)
#'
#' @family multiple comparison procedures
#' @family factorial design methods
#' @family trimmed mean methods
#' @seealso \code{\link{mcp2a}}, \code{\link{mcp3atm}}, \code{\link{lincon}}, \code{\link{bbmcpEP}}
#' @export
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

#' Multiple Comparisons for Three-Way Factorial Design (Trimmed Means)
#'
#' Tests all linear contrasts for main effects of Factors A, B, and C, plus
#' all two-way and three-way interactions in a three-way factorial design
#' using trimmed means with heteroscedastic methods.
#'
#' @param J Number of levels for Factor A
#' @param K Number of levels for Factor B
#' @param L Number of levels for Factor C
#' @param x Data in list mode or matrix with J×K×L groups. Groups ordered by
#'   C within B within A: x[[1]]=A1B1C1, x[[2]]=A1B1C2, ..., x[[L+1]]=A1B2C1, etc.
#' @inheritParams yuen
#' @param con Reserved for future use (currently unused, default: 0)
#' @inheritParams ancova
#' @param grp Optional vector to reorder groups (default: NA, sequential order assumed)
#' @param op Logical; if TRUE, performs omnibus test across all contrasts.
#'   If FALSE (default), tests each effect separately
#' @param pr Logical; if TRUE (default), prints explanatory messages
#'
#' @return A list with components:
#'   \item{Factor.A}{Results for Factor A main effects}
#'   \item{Factor.B}{Results for Factor B main effects}
#'   \item{Factor.C}{Results for Factor C main effects}
#'   \item{Factor.AB}{Results for A×B interaction contrasts}
#'   \item{Factor.AC}{Results for A×C interaction contrasts}
#'   \item{Factor.BC}{Results for B×C interaction contrasts}
#'   \item{Factor.ABC}{Results for three-way A×B×C interaction contrasts}
#'   \item{All.Tests}{If \code{op=TRUE}, omnibus results across all contrasts}
#'   \item{conA, conB, conC}{Contrast matrices for main effects}
#'   \item{conAB, conAC, conBC, conABC}{Contrast matrices for interactions}
#'
#' @details
#' This function extends \code{\link{mcp2atm}} to three-way factorial designs,
#' providing comprehensive analysis of all main effects and interaction contrasts.
#'
#' \strong{Data Organization:}
#' Data must be ordered by Factor C within Factor B within Factor A. For a 2×2×2 design:
#' \itemize{
#'   \item x[[1]]: A1,B1,C1
#'   \item x[[2]]: A1,B1,C2
#'   \item x[[3]]: A1,B2,C1
#'   \item x[[4]]: A1,B2,C2
#'   \item x[[5]]: A2,B1,C1
#'   \item x[[6]]: A2,B1,C2
#'   \item x[[7]]: A2,B2,C1
#'   \item x[[8]]: A2,B2,C2
#' }
#'
#' \strong{Contrast Generation:}
#' Automatically creates contrast matrices via \code{\link{con3way}} for:
#' \itemize{
#'   \item 3 sets of main effect contrasts (A, B, C)
#'   \item 3 sets of two-way interaction contrasts (A×B, A×C, B×C)
#'   \item 1 set of three-way interaction contrasts (A×B×C)
#' }
#'
#' \strong{Analysis Options:}
#' \itemize{
#'   \item \code{op=FALSE}: Separate hypothesis tests for each effect type
#'   \item \code{op=TRUE}: Single omnibus test combining all contrasts
#' }
#'
#' Each set of contrasts controls family-wise error rate via \code{\link{lincon}}.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' Three-way designs can produce a large number of contrasts. Consider whether
#' all interactions are theoretically meaningful before running this analysis.
#' For two-way designs, use \code{\link{mcp2atm}}.
#'
#' @examples
#' # 2 × 2 × 2 factorial design
#' set.seed(123)
#' x <- list(
#'   rnorm(15, 0.0), rnorm(15, 0.2),  # A1,B1: C1, C2
#'   rnorm(15, 0.3), rnorm(15, 0.5),  # A1,B2: C1, C2
#'   rnorm(15, 0.4), rnorm(15, 0.6),  # A2,B1: C1, C2
#'   rnorm(15, 0.7), rnorm(15, 0.9)   # A2,B2: C1, C2
#' )
#'
#' # Separate tests for each effect
#' result <- mcp3atm(J=2, K=2, L=2, x=x, tr=0.2)
#'
#' # Omnibus test
#' result2 <- mcp3atm(J=2, K=2, L=2, x=x, tr=0.2, op=TRUE)
#'
#' @family multiple comparison procedures
#' @family factorial design methods
#' @family trimmed mean methods
#' @seealso \code{\link{mcp2atm}}, \code{\link{mcp2a}}, \code{\link{con3way}}
#' @export
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


#' Three-Way Median-Based Multiple Comparisons
#'
#' Performs all pairwise comparisons for main effects and interactions in a three-way
#' design using median-based methods. Tests Factor A, B, C main effects and all
#' two-way and three-way interactions.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param L Number of levels for Factor C.
#' @param x Data in list mode or matrix. If matrix, columns correspond to groups.
#'   Groups ordered as: (1,1,1), (1,1,2), ..., (1,1,L), (1,2,1), ..., (J,K,L).
#' @param tr Trim proportion (default: 0.2, currently unused - uses medians).
#' @param con Contrast matrix (default: 0 generates appropriate contrasts).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param grp Vector for reordering groups if data not in standard order (default: NA).
#' @param op Logical. If FALSE (default), return separate results for each effect.
#'   If TRUE, combine all tests and return in \code{All.Tests}.
#'
#' @return A list with components:
#'   \item{Factor.A}{Results for Factor A main effect (unless \code{op=TRUE}).}
#'   \item{Factor.B}{Results for Factor B main effect (unless \code{op=TRUE}).}
#'   \item{Factor.C}{Results for Factor C main effect (unless \code{op=TRUE}).}
#'   \item{Factor.AB}{Results for A×B interaction (unless \code{op=TRUE}).}
#'   \item{Factor.AC}{Results for A×C interaction (unless \code{op=TRUE}).}
#'   \item{Factor.BC}{Results for B×C interaction (unless \code{op=TRUE}).}
#'   \item{Factor.ABC}{Results for A×B×C interaction (unless \code{op=TRUE}).}
#'   \item{All.Tests}{Combined results for all tests (only if \code{op=TRUE}).}
#'   \item{conA, conB, conC, conAB, conAC, conBC, conABC}{Contrast matrices used.}
#'
#' @details
#' This function uses \code{msmed} to perform median-based comparisons for all
#' main effects and interactions in a three-way design. Results are organized by
#' effect unless \code{op=TRUE}, which combines all tests into a single output.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{msmed}}, \code{\link{con3way}}, \code{\link{mcp3atm}}
#' @examples
#' \dontrun{
#' # 2x2x2 design
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20), rnorm(20, mean=0.5), rnorm(20),
#'           rnorm(20), rnorm(20), rnorm(20), rnorm(20, mean=0.5))
#' result <- mcp3med(J=2, K=2, L=2, x=x)
#' print(result$Factor.A)
#' }
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


#' Three-Way Repeated Measures Multiple Comparisons
#'
#' Performs all pairwise comparisons for main effects and interactions in a three-way
#' repeated measures (within-by-within-by-within) design using trimmed means. Tests
#' Factor A, B, C main effects and all two-way and three-way interactions.
#'
#' @inheritParams rmmcp
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param L Number of levels for Factor C.
#' @param x Data in list mode or matrix with J×K×L columns. Each column (or list element)
#'   contains repeated measures for one condition. Groups ordered as: (1,1,1), (1,1,2),
#'   ..., (1,1,L), (1,2,1), ..., (J,K,L).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param dif Logical. If TRUE (default), test hypotheses based on difference scores
#'   (more powerful when normality holds). If FALSE, use marginal measures.
#' @param op Logical. Currently not used (reserved for future functionality).
#' @param grp Vector for reordering groups if data not in standard order (default: NA).
#'
#' @return A list with components:
#'   \item{Factor.A}{Results from \code{\link{rmmcp}} for Factor A main effect.}
#'   \item{Factor.B}{Results from \code{\link{rmmcp}} for Factor B main effect.}
#'   \item{Factor.C}{Results from \code{\link{rmmcp}} for Factor C main effect.}
#'   \item{Factor.AB}{Results from \code{\link{rmmcp}} for A×B interaction.}
#'   \item{Factor.AC}{Results from \code{\link{rmmcp}} for A×C interaction.}
#'   \item{Factor.BC}{Results from \code{\link{rmmcp}} for B×C interaction.}
#'   \item{Factor.ABC}{Results from \code{\link{rmmcp}} for A×B×C three-way interaction.}
#'   \item{conA, conB, conC}{Contrast matrices for main effects.}
#'   \item{conAB, conAC, conBC, conABC}{Contrast matrices for interactions.}
#'
#' @details
#' This function generates appropriate contrast matrices using \code{\link{con3way}} and
#' applies \code{\link{rmmcp}} to test each main effect and interaction. The FWE
#' (family-wise error) rate is controlled using either Hochberg's or Rom's method,
#' depending on the \code{dif} parameter setting in \code{\link{rmmcp}}.
#'
#' Missing values are automatically removed. For three-way designs, this provides
#' a comprehensive analysis of all effects using robust trimmed means.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{rmmcp}}, \code{\link{con3way}}, \code{\link{mcp3atm}}, \code{\link{mcp3med}}
#' @examples
#' \dontrun{
#' # 2×2×2 repeated measures design
#' set.seed(123)
#' n <- 20
#' x <- matrix(rnorm(n * 8), ncol=8)
#' x[,4] <- x[,4] + 0.5  # Add effect to condition (1,2,2)
#' result <- rm3mcp(J=2, K=2, L=2, x=x)
#' print(result$Factor.A)
#' print(result$Factor.ABC)
#' }
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

#' Multiple Comparisons for Repeated Measures (Trimmed Means)
#'
#' Performs multiple comparisons among repeated measures (dependent groups) using
#' trimmed means. Family-wise error rate is controlled using Hochberg's or Rom's method.
#'
#' @param x Data matrix (n × J) or list of J dependent groups. Each row represents
#'   one subject measured across J time points or conditions.
#' @param y Optional second group (for J=2 case). If provided, data combined with x.
#' @param con Contrast matrix (J × C) where C is the number of contrasts. If 0
#'   (default), all pairwise comparisons are performed. Each column specifies one contrast.
#' @inheritParams common_params
#' @param dif Logical. If TRUE (default), use difference scores for inference.
#'   If FALSE, use marginal trimmed means with dependency adjustment.
#' @param hoch Logical. If TRUE (default), use Hochberg's method for FWE control.
#'   If FALSE and alpha ∈ {.05, .01} with ≤10 contrasts, Rom's method is used.
#' @param na.rm Logical. If TRUE (default), remove rows with missing values.
#'
#' @return A list with components:
#'   \item{n}{Sample size (number of subjects after removing missing values)}
#'   \item{test}{Matrix with columns: Group (2 cols for pairwise), test (statistic),
#'     p.value, p.crit (adjusted critical p-value), se (standard error)}
#'   \item{psihat}{Matrix with columns: Group (2 cols for pairwise), psihat (contrast estimate),
#'     ci.lower, ci.upper (simultaneous confidence intervals adjusted for FWE)}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts at level alpha}
#'
#' @details
#' This function performs simultaneous inference for multiple contrasts involving
#' repeated measures (dependent/paired data) using trimmed means to achieve robustness.
#'
#' \strong{Difference Scores vs. Marginal Means:}
#' \itemize{
#'   \item If \code{dif=TRUE} (default and recommended): Each contrast is computed
#'     as a weighted sum of the original measurements (creating a difference score),
#'     then inference is based on the trimmed mean of these difference scores.
#'   \item If \code{dif=FALSE}: Inference is based on marginal trimmed means, adjusting
#'     standard errors for the dependency structure using winsorized covariances.
#' }
#'
#' \strong{FWE Control:}
#' The family-wise error rate is controlled using either:
#' \itemize{
#'   \item Hochberg's (1988) sequentially rejective method (default): Conservative,
#'     works for any alpha and number of contrasts.
#'   \item Rom's (1990) method: Slightly more powerful than Hochberg for alpha=.05
#'     or .01 with ≤10 contrasts. Set \code{hoch=FALSE} to use.
#' }
#'
#' Confidence intervals are adjusted to be consistent with the critical p-values
#' from the chosen FWE control method.
#'
#' \strong{Contrast Specification:}
#' All pairwise comparisons (default) or custom contrasts can be specified. For
#' custom contrasts, each column of \code{con} should sum to zero (for comparing
#' groups) and specify the weights for each group in the linear combination.
#'
#' @note
#' Missing values are handled by complete-case deletion (removing entire rows).
#' For marginal comparisons allowing missing values, see \code{\link{rmmismcp}}.
#'
#' @examples
#' # Example: 4 repeated measures on same subjects
#' set.seed(123)
#' # Simulate data: 20 subjects, 4 time points
#' x <- matrix(rnorm(80), nrow=20, ncol=4)
#' x[,2] <- x[,2] + 0.5  # Add effect at time 2
#' x[,4] <- x[,4] + 0.8  # Larger effect at time 4
#'
#' # All pairwise comparisons with 20% trimming
#' result <- rmmcp(x)
#'
#' # Custom contrasts: Compare time 1 with average of times 2-4
#' con <- matrix(c(1, -1/3, -1/3, -1/3), nrow=4)
#' result2 <- rmmcp(x, con=con, tr=0.1)
#'
#' # Use Rom's method for potentially more power
#' result3 <- rmmcp(x, alpha=0.05, hoch=FALSE)
#'
#' @family multiple comparison procedures
#' @family repeated measures methods
#' @references
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance.
#'   \emph{Biometrika}, \emph{75}, 800-802.
#'
#' Rom, D. M. (1990). A sequentially rejective test procedure based on a modified
#'   Bonferroni inequality. \emph{Biometrika}, \emph{77}, 663-665.
#'
#' Wilcox, R. R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#'   (5th Ed.). Academic Press. Chapter 8.
#'
#' @seealso \code{\link{rmmcppbd}} for bootstrap version with M-estimators,
#'   \code{\link{rmmcpv2}} for improved missing value handling,
#'   \code{\link{rmmcpES}} to include effect sizes,
#'   \code{\link{rmmcpQS}} to include quantile shift effect sizes
#' @export
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


#' @rdname rmmcp
#' @export
wmcp<-rmmcp


#' Repeated Measures MCP with Percentile Bootstrap (Difference Scores)
#'
#' Performs multiple comparisons for repeated measures using percentile bootstrap
#' methodology applied to difference scores. Supports robust estimators (default:
#' one-step M-estimator) with Hochberg's method for FWE control.
#'
#' @param x Data matrix (n × J) or list of J dependent groups. Each row represents
#'   one subject measured across J time points or conditions.
#' @param y Optional second group (for J=2 case). If provided, data combined with x.
#' @param con Contrast matrix (J × C) where C is the number of contrasts. If 0
#'   (default), all pairwise comparisons are performed.
#' @inheritParams common_params
#' @param est Robust estimator function (default: \code{onestep} for one-step M-estimator).
#'   Can also use \code{tmean}, \code{median}, \code{mom}, etc.
#' @param plotit Logical. Currently not implemented (for future plotting functionality).
#' @param grp Optional vector specifying subset of groups to analyze.
#' @param nboot Number of bootstrap samples. If NA (default), chosen based on number
#'   of contrasts: 5000 for >10, 3000 for 7-10, 2000 for 5-6, 1000 for ≤4 contrasts.
#' @param hoch Logical. If TRUE (default when n<80), use Hochberg's method.
#'   If FALSE, use sequentially rejective method. Automatically TRUE if n≥80.
#' @inheritParams common_params
#' @param ... Additional arguments passed to the estimator function
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num (contrast number), psihat (estimate),
#'     ci.lower, ci.upper (bootstrap CI), p.value, p.adj (adjusted p-value)}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts at level alpha}
#'
#' @details
#' This function uses the percentile bootstrap method to compare dependent groups
#' based on \strong{difference scores}. For each contrast specified by the columns
#' of \code{con}, the function:
#' \enumerate{
#'   \item Forms difference scores: d = X × con (weighted combination of measurements)
#'   \item Computes the estimator on each difference score vector
#'   \item Resamples subjects (rows) with replacement to get bootstrap distribution
#'   \item Constructs percentile bootstrap confidence intervals
#'   \item Controls FWE using Hochberg's or sequentially rejective method
#' }
#'
#' \strong{Bootstrap Sample Size:}
#' The default \code{nboot} is adaptive based on the number of contrasts:
#' \itemize{
#'   \item ≤4 contrasts: 1000 bootstrap samples
#'   \item 5-6 contrasts: 2000 samples
#'   \item 7-10 contrasts: 3000 samples
#'   \item >10 contrasts: 5000 samples
#' }
#'
#' \strong{FWE Control:}
#' Family-wise error rate is controlled using Hochberg's (1988) sequentially rejective
#' method (automatically used if n≥80). This adjusts p-values to maintain the Type I
#' error rate across all contrasts.
#'
#' @note
#' Unlike \code{\link{rmmcp}} which uses trimmed means, this function allows any
#' robust estimator through the \code{est} argument. The percentile bootstrap provides
#' valid inference without normality assumptions.
#'
#' Rows with missing values are removed before analysis (complete-case deletion).
#'
#' @examples
#' # Example: 3 time points, 30 subjects
#' set.seed(42)
#' x <- matrix(rnorm(90), nrow=30, ncol=3)
#' x[,2] <- x[,2] + 0.5
#' x[,3] <- x[,3] + 1.0
#'
#' # All pairwise comparisons with one-step M-estimator
#' result <- rmmcppbd(x, nboot=1000)
#'
#' # Using trimmed mean instead
#' result2 <- rmmcppbd(x, est=tmean, tr=0.2, nboot=1000)
#'
#' # Custom contrast: Time 1 vs. average of Times 2-3
#' con <- matrix(c(1, -0.5, -0.5), nrow=3)
#' result3 <- rmmcppbd(x, con=con, nboot=2000)
#'
#' @family multiple comparison procedures
#' @family repeated measures methods
#' @family bootstrap methods
#' @references
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance.
#'   \emph{Biometrika}, \emph{75}, 800-802.
#'
#' Wilcox, R. R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#'   (5th Ed.). Academic Press.
#'
#' @seealso \code{\link{rmmcp}} for parametric version with trimmed means,
#'   \code{\link{rmmcppbtm}} for bootstrap trimmed means specifically,
#'   \code{\link{rmmismcp}} for allowing missing values
#' @export
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


#' Repeated Measures MCP with Percentile Bootstrap for Trimmed Means
#'
#' Performs multiple comparisons for repeated measures using percentile bootstrap
#' methodology specifically for trimmed means. This is a specialized version of
#' \code{\link{rmmcppbd}} optimized for trimmed means.
#'
#' @param x Data matrix (n × J) where rows are subjects and columns are repeated
#'   measurements, or data in list mode.
#' @param con Contrast matrix (J × C) where C is the number of contrasts. If 0
#'   (default), all pairwise comparisons are performed.
#' @inheritParams common_params
#' @param grp Optional vector specifying subset of groups to analyze.
#' @param nboot Number of bootstrap samples. If NA (default), chosen adaptively
#'   based on number of contrasts (same strategy as \code{\link{rmmcppbd}}).
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat (trimmed mean difference),
#'     p.value, p.crit (critical p-value), ci.lower, ci.upper}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts}
#'
#' @details
#' This function is equivalent to \code{rmmcppbd(x, est=tmean, tr=tr, ...)} but
#' optimized for the specific case of trimmed means. It:
#' \enumerate{
#'   \item Forms difference scores based on contrast coefficients
#'   \item Computes trimmed means of difference scores
#'   \item Uses percentile bootstrap to construct confidence intervals
#'   \item Controls FWE using Hochberg's method
#' }
#'
#' The bootstrap resamples subjects (complete rows) to preserve the dependency
#' structure among repeated measurements.
#'
#' @note
#' For general robust estimators, use \code{\link{rmmcppbd}}. This function is
#' specifically for trimmed means and may be slightly more efficient.
#'
#' @examples
#' # 4 time points, 25 subjects
#' set.seed(99)
#' x <- matrix(rnorm(100), 25, 4)
#' x[,3] <- x[,3] + 0.8
#'
#' # All pairwise comparisons with 20% trimming
#' result <- rmmcppbtm(x, nboot=1000)
#'
#' # Compare time 1 vs. average of times 2-4
#' con <- matrix(c(1, -1/3, -1/3, -1/3), nrow=4)
#' result2 <- rmmcppbtm(x, con=con, tr=0.1, nboot=2000)
#'
#' @family multiple comparison procedures
#' @family repeated measures methods
#' @family bootstrap methods
#' @seealso \code{\link{rmmcppbd}} for general estimators,
#'   \code{\link{rmmcp}} for parametric version
#' @export
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


#' Repeated Measures MCP with Bootstrap (Version 2, Improved Missing Value Handling)
#'
#' Performs multiple comparisons for repeated measures using percentile bootstrap.
#' Unlike \code{\link{rmmcppbd}}, this version allows analysis of marginal distributions
#' using all available (non-missing) data for each variable.
#'
#' @param x Data matrix (n × J) or list of J dependent groups. Missing values (NA) allowed.
#' @param y Optional second group (for J=2 case). If provided, data combined with x.
#' @param con Contrast matrix (J × C) where C is the number of contrasts. If NULL
#'   (default), all pairwise comparisons are performed.
#' @inheritParams common_params
#' @param NA.RM Logical. If TRUE (default), handle missing values by analyzing marginal
#'   distributions (each variable uses all its non-missing values).
#' @param method Method for p-value adjustment: 'hoch' (Hochberg, default) or 'fdr'
#'   (false discovery rate control).
#' @param est Robust estimator function (default: \code{tmean} for trimmed mean).
#' @param plotit Logical. If TRUE, create plots (currently not fully implemented).
#' @param dif Logical. If TRUE (default), use difference scores. If FALSE, analyze
#'   marginal distributions (each using its own available data).
#' @param grp Optional vector specifying subset of groups to analyze.
#' @param nboot Number of bootstrap samples. If NA, chosen adaptively.
#' @param BA Logical. Beh rens-Fisher approach (currently not implemented).
#' @param xlab,ylab Labels for plotting (if plotit=TRUE).
#' @param pr Logical. If TRUE (default), print informational messages.
#' @inheritParams common_params
#' @param SR Logical. Studentized range method (currently not implemented).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat (estimate), ci.lower,
#'     ci.upper, p.value, p.adj (adjusted p-value based on method)}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts at level alpha}
#'
#' @details
#' This is an improved version of \code{\link{rmmcppbd}} with better handling of
#' missing values:
#'
#' \strong{Key Differences from rmmcppbd:}
#' \itemize{
#'   \item \strong{Missing value handling}: When \code{dif=FALSE}, each marginal
#'     distribution uses all its available (non-missing) data. For example, if
#'     column 1 has missing values but columns 2-3 don't, columns 2-3 use all
#'     their data while column 1 uses only non-missing observations.
#'   \item \code{rmmcppbd} removes entire rows with any missing values (listwise deletion).
#'   \item Supports FDR control in addition to FWE control (via \code{method='fdr'}).
#' }
#'
#' \strong{Analysis Modes:}
#' \itemize{
#'   \item \code{dif=TRUE} (default): Forms difference scores for each contrast,
#'     then analyzes these (calls \code{rmmcppbd} internally).
#'   \item \code{dif=FALSE}: Analyzes marginal distributions, allowing each
#'     variable to use all its non-missing data.
#' }
#'
#' @note
#' For complete data or when listwise deletion is acceptable, \code{\link{rmmcppbd}}
#' is simpler. Use this function when you want to maximize use of available data
#' with missing values.
#'
#' @examples
#' # Example with missing values
#' set.seed(55)
#' x <- matrix(rnorm(80), nrow=20, ncol=4)
#' x[,2] <- x[,2] + 0.5
#' # Introduce some missing values
#' x[c(1,5,10), 1] <- NA
#' x[c(3,8), 3] <- NA
#'
#' # Analyze using all available data for each variable
#' result <- rmmcppbv2(x, dif=FALSE, nboot=1000)
#'
#' # Using difference scores (like rmmcppbd)
#' result2 <- rmmcppbv2(x, dif=TRUE, nboot=1000)
#'
#' # Control FDR instead of FWE
#' result3 <- rmmcppbv2(x, dif=FALSE, method='fdr', nboot=1000)
#'
#' @family multiple comparison procedures
#' @family repeated measures methods
#' @family bootstrap methods
#' @seealso \code{\link{rmmcppbd}} for standard bootstrap MCP,
#'   \code{\link{rmmismcp}} for another approach to missing values
#' @export
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


#' Repeated Measures MCP (Version 2, Improved Missing Value Handling)
#'
#' Improved version of \code{\link{rmmcp}} with better handling of missing values.
#' Analyzes each pairwise comparison using only the complete pairs for that comparison,
#' rather than removing all rows with any missing values.
#'
#' @param x Data matrix (n × J) or list of J dependent groups. Missing values allowed.
#' @param y Optional second group (for J=2 case). If provided, data combined with x.
#' @param con Contrast matrix (J × C) where C is the number of contrasts. If 0
#'   (default), all pairwise comparisons are performed.
#' @inheritParams common_params
#' @param dif Logical. If TRUE (default), use difference scores for inference.
#'   If FALSE, use marginal trimmed means with dependency adjustment.
#' @param hoch Logical. If TRUE (default), use Hochberg's method for FWE control.
#' @param na.rm Logical. If TRUE (default), handle missing values pairwise.
#' @param nmin Minimum sample size for a pairwise comparison (default: 5).
#'   Comparisons with fewer than \code{nmin} complete pairs are skipped.
#'
#' @return A list with components:
#'   \item{n}{Total sample size (may vary across comparisons with missing data)}
#'   \item{test}{Matrix with columns: Group (2 cols), test (statistic),
#'     p.value, p.adj (adjusted p-value), se}
#'   \item{psihat}{Matrix with columns: Group (2 cols), psihat (estimate),
#'     ci.lower, ci.upper}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts at level alpha}
#'
#' @details
#' This function improves upon \code{\link{rmmcp}} by using \strong{pairwise deletion}
#' of missing values instead of listwise deletion:
#'
#' \strong{Key Improvement:}
#' \itemize{
#'   \item \code{rmmcp}: Removes entire rows if any measurement is missing (listwise
#'     deletion). If one column has missing values, all data from that row is lost.
#'   \item \code{rmmcpv2}: Each pairwise comparison uses only the complete pairs
#'     for those two variables (pairwise deletion). Maximizes use of available data.
#' }
#'
#' For example, with 3 time points where time 1 has missing values:
#' \itemize{
#'   \item Time 1 vs Time 2: Uses all rows with both time 1 and time 2 observed
#'   \item Time 2 vs Time 3: Uses all rows with both time 2 and time 3 observed
#'     (includes rows where time 1 is missing!)
#' }
#'
#' \strong{Minimum Sample Size:}
#' Comparisons with fewer than \code{nmin} complete pairs are automatically skipped
#' (with a printed message). Default is 5 pairs minimum.
#'
#' \strong{Adjusted P-values:}
#' Unlike \code{rmmcp} which outputs critical p-values (\code{p.crit}), this function
#' outputs adjusted p-values (\code{p.adj}) using Hochberg's method. Reject if p.adj ≤ alpha.
#'
#' @note
#' When data are complete (no missing values), results should be very similar to
#' \code{\link{rmmcp}}. This version is primarily useful when missing data are present.
#'
#' @examples
#' # Data with missing values
#' set.seed(66)
#' x <- matrix(rnorm(100), nrow=25, ncol=4)
#' x[,2] <- x[,2] + 0.7
#' # Introduce missing values in different columns
#' x[c(1,5,10), 1] <- NA   # Missing in column 1
#' x[c(3,8,15), 3] <- NA   # Missing in column 3
#'
#' # Pairwise deletion: each comparison uses its available data
#' result_v2 <- rmmcpv2(x)
#'
#' # Compare with listwise deletion (removes rows 1,3,5,8,10,15)
#' result_orig <- rmmcp(x)
#'
#' # Set higher minimum sample size
#' result_strict <- rmmcpv2(x, nmin=20)
#'
#' @family multiple comparison procedures
#' @family repeated measures methods
#' @seealso \code{\link{rmmcp}} for standard version (listwise deletion),
#'   \code{\link{rmmcppbv2}} for bootstrap version with improved missing value handling
#' @export
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


#' Repeated Measures MCP with Effect Sizes
#'
#' Performs multiple comparisons for repeated measures using trimmed means and
#' includes robust effect sizes (AKP robust Cohen's d or explanatory measure).
#' This extends \code{\link{rmmcp}} by adding effect size estimates.
#'
#' @param x Data matrix (n × J) where rows are subjects and columns are repeated
#'   measurements, or data in list mode.
#' @param con Must be 0 (default). Currently only supports all pairwise comparisons.
#' @inheritParams common_params
#' @param dif Logical. If TRUE (default), compute AKP robust Cohen's d based on
#'   difference scores. If FALSE, use explanatory measure of effect size.
#' @param hoch Logical. If TRUE (default), use Hochberg's method for FWE control.
#' @param pr Logical. If TRUE (default), print informational messages.
#'
#' @return A list with components:
#'   \item{test}{Matrix with columns: Group (2 cols), test (statistic),
#'     p.value, p.crit (critical p-value), se}
#'   \item{psihat}{Matrix with columns: Group (2 cols), psihat (trimmed mean difference),
#'     ci.lower, ci.upper, Effect.Size (robust effect size)}
#'
#' @details
#' This function extends \code{\link{rmmcp}} by computing and reporting robust
#' effect sizes for each pairwise comparison. The effect size measure depends
#' on the \code{dif} parameter:
#'
#' \strong{Effect Size Measures:}
#' \itemize{
#'   \item \code{dif=TRUE} (default): AKP robust Cohen's d based on difference scores.
#'     Computed via \code{D.akp.effect()}. This is a robust version of Cohen's d
#'     that uses trimmed means and winsorized variances.
#'   \item \code{dif=FALSE}: Explanatory measure of effect size from \code{yuendv2()}.
#'     This quantifies the proportion of variance explained, robust to outliers.
#'     A message is printed noting this choice.
#' }
#'
#' All other aspects (hypothesis tests, confidence intervals, FWE control) are
#' identical to \code{\link{rmmcp}}.
#'
#' @note
#' Currently restricted to all pairwise comparisons (\code{con=0}). Custom contrasts
#' are not supported. An error is raised if \code{con} is not 0.
#'
#' Effect sizes provide a standardized measure of the magnitude of differences,
#' complementing the hypothesis tests and confidence intervals.
#'
#' @examples
#' # 4 time points, 30 subjects
#' set.seed(777)
#' x <- matrix(rnorm(120), 30, 4)
#' x[,2] <- x[,2] + 0.5
#' x[,4] <- x[,4] + 1.2
#'
#' # All pairwise comparisons with robust effect sizes
#' result <- rmmcpES(x)
#' print(result$psihat)  # Includes Effect.Size column
#'
#' # Using explanatory measure of effect size
#' result2 <- rmmcpES(x, dif=FALSE)
#'
#' # With different trimming
#' result3 <- rmmcpES(x, tr=0.1)
#'
#' @family multiple comparison procedures
#' @family repeated measures methods
#' @family effect size methods
#' @seealso \code{\link{rmmcp}} for version without effect sizes,
#'   \code{\link{rmmcpQS}} for quantile shift effect sizes,
#'   \code{\link{D.akp.effect}} for AKP effect size details
#' @export
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


#' Repeated Measures MCP with Quantile Shift Effect Sizes
#'
#' Performs multiple comparisons for repeated measures using trimmed means and
#' includes quantile shift measures of effect size. This extends \code{\link{rmmcp}}
#' by adding distribution-based effect size estimates.
#'
#' @param x Data matrix (n × J) where rows are subjects and columns are repeated
#'   measurements, or data in list mode.
#' @param y Optional second group (for J=2 case). If provided, data combined with x.
#' @param con Contrast matrix (J × C) where C is the number of contrasts. If 0
#'   (default), all pairwise comparisons are performed.
#' @inheritParams common_params
#' @param dif Logical. If TRUE (default), compute quantile shift on difference scores.
#'   If FALSE, use marginal distributions.
#' @param hoch Logical. If TRUE (default), use Hochberg's method for FWE control.
#'   If FALSE and alpha ∈ {.05, .01} with ≤10 contrasts, Rom's method is used.
#' @param locfun Location function for quantile shift computation (default: \code{tmean}).
#' @param ... Additional arguments passed to \code{locfun}.
#'
#' @return A list with components:
#'   \item{n}{Sample size (after removing missing values)}
#'   \item{test}{Matrix with columns: Group (2 cols), test (statistic),
#'     p.value, p.crit (critical p-value), se}
#'   \item{psihat}{Matrix with columns: Group (2 cols), psihat (trimmed mean difference),
#'     ci.lower, ci.upper, Q.effect (quantile shift effect size)}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts}
#'
#' @details
#' This function extends \code{\link{rmmcp}} by computing quantile shift (Q) effect
#' sizes, which provide a distribution-based measure of how much one distribution
#' is shifted relative to another.
#'
#' \strong{Quantile Shift Effect Size:}
#' The Q effect size estimates what proportion of observations from one distribution
#' would need to be shifted to make the distributions identical. Computed via:
#' \itemize{
#'   \item \code{dif=TRUE}: \code{depQS()} for dependent groups using difference scores
#'   \item \code{dif=FALSE}: Still uses \code{depQS()} but applied to marginal distributions
#' }
#'
#' Values range from 0 (identical distributions) to 1 (completely separated distributions).
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item Q ≈ 0.1-0.2: Small effect
#'   \item Q ≈ 0.3-0.4: Medium effect
#'   \item Q ≥ 0.5: Large effect
#' }
#'
#' All hypothesis testing aspects are identical to \code{\link{rmmcp}}.
#'
#' @note
#' For custom contrasts (not just pairwise), quantile shift is computed via
#' \code{lindQS()} for each contrast.
#'
#' @examples
#' # 3 time points, 25 subjects
#' set.seed(888)
#' x <- matrix(rnorm(75), 25, 3)
#' x[,2] <- x[,2] + 0.6
#' x[,3] <- x[,3] + 1.1
#'
#' # All pairwise comparisons with quantile shift effect sizes
#' result <- rmmcpQS(x)
#' print(result$psihat)  # Includes Q.effect column
#'
#' # Using median instead of trimmed mean
#' result2 <- rmmcpQS(x, locfun=median)
#'
#' # Compare marginal distributions
#' result3 <- rmmcpQS(x, dif=FALSE)
#'
#' @family multiple comparison procedures
#' @family repeated measures methods
#' @family effect size methods
#' @references
#' Wilcox, R. R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing}
#'   (5th Ed.). Academic Press.
#' @seealso \code{\link{rmmcp}} for version without effect sizes,
#'   \code{\link{rmmcpES}} for AKP robust Cohen's d effect sizes,
#'   \code{\link{depQS}} for quantile shift details
#' @export
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


#' Repeated Measures MCP with Improved Missing Value Support
#'
#' Performs multiple comparisons for repeated measures using percentile bootstrap,
#' with special handling for missing values. Unlike \code{\link{rmmcppbd}}, this
#' function compares marginal measures of location without removing rows with missing values.
#'
#' @param x Data matrix (n × J) or list of J dependent groups. Missing values allowed.
#' @param y Optional second group (for J=2 case). If provided, data combined with x.
#' @param con Contrast matrix (J × C) where C is the number of contrasts. If 0
#'   (default), all pairwise comparisons are performed.
#' @inheritParams common_params
#' @param est Robust estimator function (default: \code{tmean} for 20% trimmed mean).
#' @param plotit Logical. If TRUE, create plots (currently not fully implemented).
#' @param grp Optional vector specifying subset of groups to analyze.
#' @param nboot Number of bootstrap samples (default: 500).
#' @inheritParams common_params
#' @param xlab,ylab Labels for plotting (if plotit=TRUE).
#' @param pr Logical. If TRUE, print progress/diagnostic messages (default: FALSE).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat (estimate), ci.lower,
#'     ci.upper, p.value, p.adj (adjusted p-value)}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant contrasts at level alpha}
#'
#' @details
#' This function performs bootstrap-based multiple comparisons for repeated measures
#' with a focus on handling missing values effectively.
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Missing value handling}: Compares \strong{marginal} measures of
#'     location. Each variable uses all its non-missing observations. Unlike
#'     \code{rmmcppbd}, rows with missing values are NOT removed.
#'   \item \strong{Bootstrap resampling}: Resamples subjects (rows) to preserve
#'     dependency structure, then computes the estimator on each variable's
#'     available data.
#'   \item \strong{Flexible estimators}: Supports any robust estimator (default
#'     is trimmed mean with 20% trimming).
#' }
#'
#' \strong{Comparison with Other Functions:}
#' \itemize{
#'   \item \code{rmmcppbd}: Removes all rows with any missing values (listwise deletion),
#'     analyzes difference scores
#'   \item \code{rmmcppbv2}: When \code{dif=FALSE}, similar approach but different implementation
#'   \item \code{rmmismcp}: Focuses on marginal comparisons, maximizes use of available data
#' }
#'
#' @note
#' This function is specifically designed for scenarios where:
#' \enumerate{
#'   \item Missing data are present
#'   \item You want to compare marginal distributions (not difference scores)
#'   \item You want to maximize use of all available observations
#' }
#'
#' For complete data or when difference scores are preferred, use \code{\link{rmmcppbd}}.
#'
#' @examples
#' # Data with missing values
#' set.seed(111)
#' x <- matrix(rnorm(80), 20, 4)
#' x[,2] <- x[,2] + 0.5
#' x[,4] <- x[,4] + 1.0
#' # Introduce missing values
#' x[c(1,5,8), 1] <- NA
#' x[c(2,9), 3] <- NA
#' x[c(4,7,10), 4] <- NA
#'
#' # Compare marginal distributions using all available data
#' result <- rmmismcp(x, nboot=1000)
#'
#' # Use median instead of trimmed mean
#' result2 <- rmmismcp(x, est=median, nboot=1000)
#'
#' # Custom contrasts
#' con <- matrix(c(1, -1, 0, 0,   # Time 1 vs 2
#'                 0, 0, 1, -1),  # Time 3 vs 4
#'               nrow=4, ncol=2)
#' result3 <- rmmismcp(x, con=con, nboot=1000)
#'
#' @family multiple comparison procedures
#' @family repeated measures methods
#' @family bootstrap methods
#' @seealso \code{\link{rmmcppbd}} for standard bootstrap MCP,
#'   \code{\link{rmmcppbv2}} for another approach to missing values,
#'   \code{\link{rmmcpv2}} for parametric version with pairwise deletion
#' @export
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

#' Multiple Comparisons for Between-Within Two-Way Design (Bootstrap-t)
#'
#' Bootstrap-t method for all pairwise comparisons of main effects and interactions
#' in a two-factor design with one between-subjects factor (A) and one within-subjects
#' factor (B). Uses trimmed means and dependent contrasts for within-subjects comparisons.
#'
#' @param J Number of levels for Factor A (between-subjects)
#' @param K Number of levels for Factor B (within-subjects)
#' @param x Data in list mode or matrix. In list mode, \code{x[[1]]} contains level 1,1;
#'   \code{x[[2]]} contains level 1,2; ...; \code{x[[K]]} contains level 1,K;
#'   \code{x[[K+1]]} contains level 2,1; etc. In matrix mode, columns correspond to groups.
#' @inheritParams tmean
#' @param JK Total number of groups (J × K), used for validation
#' @param con Unused (kept for compatibility)
#' @inheritParams linconb
#' @inheritParams mcppb
#' @param method P-value adjustment method passed to \code{\link[stats]{p.adjust}}
#'   (default: \code{'hoch'} for Hochberg)
#' @inheritParams yuen
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{Fac.A}}{Matrix with test statistics for Factor A main effects (pairwise comparisons)}
#'     \item{\code{Fac.B}}{Matrix with test statistics for Factor B main effects (pairwise comparisons)}
#'     \item{\code{Fac.AB}}{Matrix with test statistics for A×B interaction contrasts}
#'     \item{\code{contrast.coef}}{List containing contrast matrices (conA, conB, conAB) from \code{\link{con2way}}}
#'   }
#'   Each matrix contains: estimate, standard error, test statistic, critical value,
#'   p-value, and adjusted p-value.
#'
#' @details
#' This function performs multiple comparisons for a between-within (mixed) two-way design:
#' \itemize{
#'   \item \strong{Factor A (between-subjects)}: Independent groups, bootstrap resampling within each level
#'   \item \strong{Factor B (within-subjects)}: Dependent measurements, uses dependent contrasts with covariance matrix
#'   \item \strong{Contrasts}: Automatically generates all pairwise comparisons using \code{\link{con2way}}
#'   \item \strong{Method}: Bootstrap-t with trimmed means (default 20\% trimming)
#'   \item \strong{FWE control}: Rom's method for dependent contrasts + p-value adjustment via \code{method}
#' }
#'
#' The bootstrap procedure:
#' 1. Centers data within each group using trimmed means
#' 2. Resamples within each level of the between-subjects factor
#' 3. Computes bootstrap critical values for each effect separately
#' 4. Applies additional p-value adjustment (Hochberg by default)
#'
#' For within-subjects comparisons, the function uses \code{\link{lindep}} which
#' accounts for correlations via the covariance matrix estimated with \code{\link{covmtrim}}.
#'
#' @note
#' \itemize{
#'   \item Missing values are removed pairwise within each level of Factor A
#'   \item Data should have equal sample sizes within each level of Factor A
#'   \item Set \code{SEED=TRUE} (default) for reproducible results
#' }
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing} (5th ed.).
#' Academic Press. Section on mixed designs and multiple comparisons.
#'
#' @seealso
#' \code{\link{bwwmcp}} for three-way between-within-within designs,
#' \code{\link{bbwmcp}} for between-between-within designs,
#' \code{\link{bwbmcp}} for alternative Rom's method implementation,
#' \code{\link{con2way}} for contrast generation,
#' \code{\link{lindep}} for dependent contrasts
#'
#' @family between-within MCP
#' @family multiple comparison procedures
#' @export
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


#' Multiple Comparisons for Between-Within-Within Three-Way Design (Bootstrap-t)
#'
#' Bootstrap-t method for all pairwise comparisons of main effects and interactions
#' in a three-factor design with one between-subjects factor (A) and two within-subjects
#' factors (B and C). Extends \code{\link{bwmcp}} to three-way designs.
#'
#' @param J Number of levels for Factor A (between-subjects)
#' @param K Number of levels for Factor B (within-subjects)
#' @param L Number of levels for Factor C (within-subjects)
#' @param x Data in list mode or matrix with J×K×L groups. Groups are ordered as:
#'   level 1,1,1; 1,1,2; ...; 1,1,L; 1,2,1; ...; 1,K,L; 2,1,1; ...; J,K,L
#' @inheritParams tmean
#' @param JKL Total number of groups (J × K × L), used for validation
#' @param con Unused (kept for compatibility)
#' @inheritParams linconb
#' @inheritParams mcppb
#' @inheritParams yuen
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{Fac.A}}{Matrix for Factor A main effects}
#'     \item{\code{Fac.B}}{Matrix for Factor B main effects}
#'     \item{\code{Fac.C}}{Matrix for Factor C main effects}
#'     \item{\code{Fac.AB}}{Matrix for A×B interaction}
#'     \item{\code{Fac.AC}}{Matrix for A×C interaction}
#'     \item{\code{Fac.BC}}{Matrix for B×C interaction}
#'     \item{\code{Fac.ABC}}{Matrix for A×B×C three-way interaction}
#'   }
#'   Each matrix contains: estimate, standard error, test statistic, critical value,
#'   and p-value for all pairwise comparisons.
#'
#' @details
#' This function performs multiple comparisons for a three-way design where:
#' \itemize{
#'   \item \strong{Factor A}: Between-subjects (independent groups)
#'   \item \strong{Factor B}: Within-subjects (repeated measures)
#'   \item \strong{Factor C}: Within-subjects (repeated measures)
#'   \item \strong{Contrasts}: Automatically generated using \code{\link{con3way}}
#'   \item \strong{Method}: Bootstrap-t with trimmed means
#'   \item \strong{FWE control}: Rom's method for dependent comparisons
#' }
#'
#' The bootstrap procedure resamples within each level of the between-subjects factor
#' while preserving the within-subjects dependencies through the covariance matrix.
#'
#' @note
#' \itemize{
#'   \item This function tests 7 families of hypotheses (3 main effects, 3 two-way interactions, 1 three-way interaction)
#'   \item Requires equal sample sizes within each level of Factor A
#'   \item Can be computationally intensive with large \code{nboot} values
#'   \item Set \code{SEED=TRUE} (default) for reproducibility
#' }
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing} (5th ed.).
#' Academic Press.
#'
#' @seealso
#' \code{\link{bwmcp}} for two-way between-within designs,
#' \code{\link{bbwmcp}} for between-between-within designs,
#' \code{\link{con3way}} for three-way contrast generation,
#' \code{\link{wwwmcppb}} for within-within-within designs
#'
#' @family between-within MCP
#' @family multiple comparison procedures
#' @export
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


#' Multiple Comparisons for Between-Between-Within Three-Way Design (Bootstrap-t)
#'
#' Bootstrap-t method for all pairwise comparisons of main effects and interactions
#' in a three-factor design with two between-subjects factors (A and B) and one
#' within-subjects factor (C).
#'
#' @param J Number of levels for Factor A (between-subjects)
#' @param K Number of levels for Factor B (between-subjects)
#' @param L Number of levels for Factor C (within-subjects)
#' @param x Data in list mode or matrix with J×K×L groups. Groups are ordered as:
#'   level 1,1,1; 1,1,2; ...; 1,1,L; 1,2,1; ...; 1,K,L; 2,1,1; ...; J,K,L
#' @inheritParams tmean
#' @param JKL Total number of groups (J × K × L)
#' @param con Unused (kept for compatibility)
#' @inheritParams linconb
#' @inheritParams mcppb
#' @inheritParams yuen
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{Fac.A}}{Matrix for Factor A main effects}
#'     \item{\code{Fac.B}}{Matrix for Factor B main effects}
#'     \item{\code{Fac.C}}{Matrix for Factor C main effects}
#'     \item{\code{Fac.AB}}{Matrix for A×B interaction}
#'     \item{\code{Fac.AC}}{Matrix for A×C interaction}
#'     \item{\code{Fac.BC}}{Matrix for B×C interaction}
#'     \item{\code{Fac.ABC}}{Matrix for A×B×C three-way interaction}
#'   }
#'   Each matrix contains: estimate, standard error, test statistic, critical value,
#'   and p-value for all pairwise comparisons.
#'
#' @details
#' This function performs multiple comparisons for a three-way design where:
#' \itemize{
#'   \item \strong{Factor A}: Between-subjects (independent)
#'   \item \strong{Factor B}: Between-subjects (independent)
#'   \item \strong{Factor C}: Within-subjects (repeated measures)
#'   \item \strong{Contrasts}: Automatically generated using \code{\link{con3way}}
#'   \item \strong{Method}: Bootstrap-t with trimmed means
#'   \item \strong{FWE control}: Rom's method for within-subjects comparisons
#' }
#'
#' The bootstrap procedure resamples within each of the J×K between-subjects cells
#' (combinations of Factors A and B) while preserving dependencies for the within-subjects
#' factor through the covariance matrix.
#'
#' @note
#' \itemize{
#'   \item Requires equal sample sizes within each of the J×K between-subjects cells
#'   \item Tests 7 families of hypotheses (3 main effects, 3 two-way interactions, 1 three-way)
#'   \item Set \code{SEED=TRUE} (default) for reproducible results
#' }
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing} (5th ed.).
#' Academic Press.
#'
#' @seealso
#' \code{\link{bwmcp}} for two-way between-within designs,
#' \code{\link{bwwmcp}} for between-within-within designs,
#' \code{\link{con3way}} for three-way contrast generation
#'
#' @family between-within MCP
#' @family multiple comparison procedures
#' @export
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


#' Multiple Comparisons for Between-Within Design Using Rank-Based Method
#'
#' Performs all pairwise comparisons of main effects and interactions for a
#' between-within two-way design using a rank-based method that tests for equal
#' distributions (not just equal medians).
#'
#' @param J Number of levels for Factor A (between-subjects)
#' @param K Number of levels for Factor B (within-subjects)
#' @param x Data in list mode or matrix. See Details for data structure.
#' @inheritParams mcppb
#' @inheritParams linconb
#' @param bhop Logical; if \code{TRUE}, uses Benjamini-Hochberg procedure for FWE control.
#'   If \code{FALSE}, uses Hochberg's method.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{Fac.A}}{Matrix of test statistics for Factor A main effect comparisons}
#'     \item{\code{Fac.B}}{Matrix of test statistics for Factor B main effect comparisons}
#'     \item{\code{Fac.AB}}{Matrix of test statistics for A×B interaction comparisons}
#'   }
#'   Each matrix contains Wilcoxon-type test statistics and adjusted p-values.
#'
#' @details
#' This function uses a rank-based (distribution-free) method to compare groups:
#' \itemize{
#'   \item \strong{Factor A comparisons}: Uses \code{\link{cidv2}} (for independent groups)
#'   \item \strong{Factor B comparisons}: Uses \code{\link{wmwpb}} (for dependent groups)
#'   \item \strong{Interaction comparisons}: Tests differences in Factor B effects across Factor A levels
#'   \item \strong{Null hypothesis}: Equal distributions (not just equal medians/means)
#'   \item \strong{FWE control}: Benjamini-Hochberg or Hochberg adjustment
#' }
#'
#' **Data structure**: Data in \code{x} should be ordered as:
#' \itemize{
#'   \item \code{x[[1]]}: Level 1,1 (Factor A level 1, Factor B level 1)
#'   \item \code{x[[2]]}: Level 1,2 (Factor A level 1, Factor B level 2)
#'   \item ...
#'   \item \code{x[[K]]}: Level 1,K
#'   \item \code{x[[K+1]]}: Level 2,1
#'   \item etc.
#' }
#'
#' Use \code{grp} to reorder groups if needed. For example, \code{grp=c(2,4,3,1)} for
#' a 2×2 design maps: group 2 → 1,1; group 4 → 1,2; group 3 → 2,1; group 1 → 2,2.
#'
#' Missing values are automatically removed pairwise.
#'
#' @note
#' This is a distribution-free alternative to \code{\link{bwmcp}} that does not assume
#' specific distributional forms. It is robust to outliers and heavy-tailed distributions.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing} (5th ed.).
#' Academic Press.
#'
#' @seealso
#' \code{\link{bwmcp}} for trimmed mean version,
#' \code{\link{cidv2}} for independent group rank tests,
#' \code{\link{wmwpb}} for dependent group rank tests
#'
#' @family between-within MCP
#' @family multiple comparison procedures
#' @family rank-based methods
#' @export
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

#' Multiple Comparisons for Split-Plot Designs - Factor A Main Effects
#'
#' Performs all pairwise comparisons among levels of Factor A (between-subjects)
#' in a split-plot design using trimmed means. Data among dependent groups
#' (Factor B levels) are pooled for each level of Factor A.
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (within-subjects factor).
#' @param x Data in matrix or list format. If a matrix, columns represent groups.
#'   If a list, \code{x[[1]]} contains data for level (1,1), \code{x[[2]]} for
#'   level (1,2), ..., \code{x[[K]]} for level (1,K), \code{x[[K+1]]} for
#'   level (2,1), etc. Length should be J*K.
#' @param tr Proportion to trim (default: 0.2 for 20% trimming).
#' @param JK Total number of groups (default: J*K).
#' @param grp Vector indicating which groups to include (default: all J*K groups).
#' @param alpha Family-wise error rate (default: 0.05).
#' @param op Logical. If TRUE, uses omnibus pooling approach with contrasts.
#'   If FALSE, pools data for each level of A and calls \code{\link{lincon}}
#'   (default: TRUE).
#'
#' @return When \code{op=FALSE}, returns output from \code{\link{lincon}}.
#'   When \code{op=TRUE}, returns output from contrast-based approach with
#'   components similar to \code{\link{bwmcp}}.
#'
#' @details
#' This function compares the J levels of Factor A (between-subjects factor) in
#' a split-plot design by pooling data across the K levels of Factor B
#' (within-subjects factor) for each level of A.
#'
#' Two approaches are available:
#' \itemize{
#'   \item When \code{op=FALSE}: Pools data for each level of Factor A, then
#'     performs pairwise comparisons using \code{\link{lincon}}.
#'   \item When \code{op=TRUE}: Uses omnibus pooling with a contrast matrix
#'     approach, accounting for the repeated measures structure.
#' }
#'
#' Rom's method is used to control the family-wise error rate.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{spmcpa}}, \code{\link{bwmcp}}, \code{\link{lincon}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 3x4 split-plot design (3 groups, 4 time points)
#' set.seed(123)
#' J <- 3
#' K <- 4
#' n <- 20
#' x <- list()
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     x[[(j-1)*K + k]] <- rnorm(n, mean = j)  # Differ on Factor A
#'   }
#' }
#'
#' # Compare Factor A levels (pooled across Factor B)
#' result <- bwamcp(J, K, x)
#' }
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


#' Multiple Comparisons for Interactions in Between-Within Design
#'
#' Tests all pairwise interaction contrasts in a between-within (split-plot) design
#' by computing difference scores for all pairs of dependent groups (within-subjects)
#' and testing whether these differences vary across levels of the between-subjects factor.
#'
#' @param J Number of levels for Factor A (between-subjects)
#' @param K Number of levels for Factor B (within-subjects)
#' @inheritParams bwmcp
#' @inheritParams tmean
#' @inheritParams linconb
#' @param method P-value adjustment method passed to \code{\link[stats]{p.adjust}}
#'   (default: \code{'hoch'} for Hochberg)
#'
#' @return A list with test results for interaction comparisons (format from \code{\link{lincon.old}}):
#'   matrix with columns for contrasts, estimates, test statistics, critical values,
#'   and p-values.
#'
#' @details
#' This function focuses specifically on interaction contrasts:
#' \itemize{
#'   \item \strong{Method}: For each pair of within-subjects levels, computes difference scores
#'   \item \strong{Test}: Compares these difference scores across between-subjects groups using trimmed means
#'   \item \strong{Interpretation}: Significant results indicate the pattern of within-subjects
#'     differences varies across between-subjects groups
#'   \item \strong{FWE control}: Hochberg's method (or other methods via \code{method})
#' }
#'
#' For M-estimators or MOM (modified one-step M-estimator), use \code{\link{spmcpi}}
#' which implements a bootstrap method.
#'
#' @note
#' This function tests only interaction contrasts. For main effects, use \code{\link{bwmcp}}.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing} (5th ed.).
#' Academic Press.
#'
#' @seealso
#' \code{\link{bwmcp}} for complete between-within analysis (main effects and interactions),
#' \code{\link{spmcpi}} for bootstrap version with M-estimators,
#' \code{\link{lincon.old}} for underlying linear contrast method
#'
#' @family between-within MCP
#' @family multiple comparison procedures
#' @export
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


#' Between-Within Interaction Contrasts with Effect Sizes
#'
#' Performs multiple comparisons for interactions in a split-plot (between-within) design
#' using trimmed means, with additional effect size measures. Tests whether difference
#' scores between dependent groups differ across levels of the between-subjects factor.
#'
#' @inheritParams bwimcp
#' @param CI Logical. If TRUE, compute confidence intervals for effect sizes (default: FALSE).
#' @param REL.M Relative magnitude parameter for effect size computation (default: NULL).
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: A, A, B, B, psihat, p.value, p.crit.
#'     Each row represents a test comparing difference scores (B - B') between two
#'     levels of Factor A.}
#'   \item{Effect.Sizes}{Effect size measures for interactions from \code{\link{bw.es.I}}.}
#'
#' @details
#' This function extends \code{\link{bwimcp}} by also computing effect sizes via
#' \code{\link{bw.es.I}}. It takes difference scores among all pairs of within-subjects
#' groups and tests which differences vary across between-subjects groups. FWE is controlled
#' using Hochberg's method.
#'
#' For robust M-estimators or MOM, use \code{\link{spmcpi}} which employs bootstrap methods.
#' The function removes missing values listwise within each level of Factor A.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bwimcp}}, \code{\link{bw.es.I}}, \code{\link{spmcpi}}, \code{\link{bwimcpAKP}}
#' @examples
#' \dontrun{
#' # 2×3 split-plot design
#' set.seed(123)
#' J <- 2  # 2 between-subjects groups
#' K <- 3  # 3 within-subjects conditions
#' n <- 20
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(n)
#' result <- bwimcpES(J=2, K=3, x=x)
#' print(result$output)
#' print(result$Effect.Sizes)
#' }
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


#' Between-Within Interaction Contrasts with AKP Effect Sizes
#'
#' Performs multiple comparisons for interactions in a split-plot (between-within) design
#' using trimmed means, with Algina-Keselman-Penfield (AKP) robust effect sizes analogous
#' to Cohen's d. Tests whether difference scores between dependent groups differ across
#' levels of the between-subjects factor.
#'
#' @inheritParams bwimcp
#' @param SEED Logical. If TRUE (default), set random seed for reproducibility.
#'
#' @return A matrix with columns:
#'   \describe{
#'     \item{A, A}{Indices of the two levels of Factor A being compared.}
#'     \item{B, B}{Indices of the two levels of Factor B (within-subjects) being compared.}
#'     \item{psihat}{Estimated difference in trimmed means between groups.}
#'     \item{p.value}{P-value for the test.}
#'     \item{p.crit}{Critical p-value after FWE adjustment.}
#'     \item{Effect.Size}{AKP robust effect size (trimmed mean difference / Winsorized SD).}
#'   }
#'
#' @details
#' This function extends \code{\link{bwimcp}} by including AKP effect sizes for each
#' interaction contrast. The AKP effect size is a robust analog of Cohen's d using
#' trimmed means and Winsorized standard deviations, computed via \code{\link{D.akp.effect}}.
#'
#' The analysis takes difference scores among all pairs of within-subjects groups and
#' tests which differences vary across between-subjects groups. FWE is controlled using
#' Hochberg's method.
#'
#' For robust M-estimators or MOM, use \code{\link{spmcpi}} which employs bootstrap methods.
#' Missing values are removed listwise within each level of Factor A.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bwimcp}}, \code{\link{bwimcpES}}, \code{\link{D.akp.effect}}, \code{\link{akp.effect}}
#' @examples
#' \dontrun{
#' # 2×3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(n)
#' result <- bwimcpAKP(J=2, K=3, x=x)
#' print(result)
#' }
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


#' Multiple Comparisons for Factor B in Between-Within Design (Rom's Method)
#'
#' Performs all pairwise comparisons among levels of the within-subjects factor (Factor B)
#' in a between-within design using trimmed means and Rom's method for FWE control.
#'
#' @param J Number of levels for Factor A (between-subjects)
#' @param K Number of levels for Factor B (within-subjects)
#' @inheritParams bwmcp
#' @inheritParams tmean
#' @inheritParams rmmcp
#' @param pool Logical; if \code{TRUE}, pools data across levels of Factor A before comparison.
#'   If \code{FALSE} (default), analyzes Factor B comparisons separately for each level of Factor A.
#' @param hoch Logical; if \code{TRUE} (default), uses Hochberg's method for FWE control.
#'   If \code{FALSE}, uses Rom's method.
#' @param pr Logical; if \code{TRUE} (default), prints results to console.
#'
#' @return Results from \code{\link{rmmcp}} (Rom's method for repeated measures comparisons):
#'   matrix with contrasts, estimates, test statistics, critical values, and p-values
#'   for all pairwise comparisons of Factor B.
#'
#' @details
#' This function focuses on comparisons of the within-subjects factor (Factor B):
#' \itemize{
#'   \item If \code{pool=FALSE}: Performs Factor B comparisons within each level of Factor A
#'   \item If \code{pool=TRUE}: Pools all data across Factor A levels before comparing Factor B levels
#'   \item \strong{Method}: Uses \code{\link{rmmcp}} with trimmed means
#'   \item \strong{FWE control}: Rom's method (accounts for dependencies in repeated measures)
#' }
#'
#' This is useful when the primary interest is in the within-subjects factor and you want
#' to either:
#' 1. Test Factor B effects separately for each between-subjects group, or
#' 2. Test overall Factor B effects ignoring the between-subjects factor
#'
#' @note
#' For complete factorial analysis including interactions, use \code{\link{bwmcp}} instead.
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing} (5th ed.).
#' Academic Press.
#'
#' @seealso
#' \code{\link{bwmcp}} for complete between-within analysis,
#' \code{\link{rmmcp}} for underlying repeated measures comparisons,
#' \code{\link{bwimcp}} for interaction comparisons
#'
#' @family between-within MCP
#' @family multiple comparison procedures
#' @export
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


#' Median-Based Between-Within Interaction Comparisons
#'
#' Performs multiple comparisons for interactions in a split-plot (between-within) design
#' using medians instead of trimmed means. Tests whether difference scores (medians)
#' between pairs of within-subjects conditions differ across between-subjects groups.
#'
#' @param J Number of levels for Factor A (between-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data in list mode or matrix. Each list element (or column) contains data
#'   for one cell of the J×K design. Groups ordered as: (1,1), (1,2), ..., (1,K),
#'   (2,1), ..., (J,K).
#' @param JK Total number of groups (should equal J×K, default: J*K).
#' @param grp Vector for reordering groups if data not in standard order
#'   (default: c(1:JK)).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#'
#' @return A matrix with columns:
#'   \describe{
#'     \item{A, A}{Indices of the two levels of Factor A (between) being compared.}
#'     \item{B, B}{Indices of the two levels of Factor B (within) being compared.}
#'     \item{test}{Test statistic (absolute value of the test).}
#'     \item{p.value}{P-value for the test.}
#'     \item{p.crit}{Critical p-value after FWE adjustment.}
#'   }
#'
#' @details
#' This is a median-based analog of \code{\link{bwimcp}}. For each pair of
#' between-subjects groups (levels of A) and each pair of within-subjects conditions
#' (levels of B), it computes difference scores and tests whether the median difference
#' varies across the two A groups using a robust test.
#'
#' The function controls the family-wise error rate using Hochberg's method. Missing
#' values are removed listwise within each level of Factor A.
#'
#' For trimmed mean-based analysis, use \code{\link{bwimcp}} instead.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bwimcp}}, \code{\link{bwmedbmcp}}, \code{\link{med2g}}
#' @examples
#' \dontrun{
#' # 2×3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(n)
#' result <- bwmedimcp(J=2, K=3, x=x)
#' print(result)
#' }
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


#' Median-Based Between-Within Factor B Comparisons
#'
#' Performs pairwise comparisons among levels of Factor B (within-subjects) in a
#' split-plot (between-within) design using medians. Data can be pooled across levels
#' of Factor A or analyzed separately.
#'
#' @param J Number of levels for Factor A (between-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data in list mode or matrix. Each list element (or column) contains data
#'   for one cell of the J×K design.
#' @param JK Total number of groups (should equal J×K, default: J*K).
#' @param grp Vector for reordering groups if data not in standard order
#'   (default: c(1:JK)).
#' @param con Contrast matrix (default: 0 generates all pairwise comparisons).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param dif Logical. If FALSE (default), use marginal measures. If TRUE, use
#'   difference scores (currently FALSE is recommended).
#' @param pool Logical. If TRUE, pool data across Factor A levels before comparing
#'   Factor B levels. If FALSE (default), analyze Factor A levels separately.
#' @param bop Logical. If TRUE, use bootstrap estimates of standard errors.
#'   If FALSE (default), use analytic approach (default: FALSE).
#' @param nboot Number of bootstrap samples if \code{bop=TRUE} (default: 100).
#' @param SEED Logical. If TRUE (default), set random seed for bootstrap reproducibility.
#'
#' @return A list with components:
#'   \item{TESTS.4.EACH.LEVEL.OF.A}{List of test results for each level of Factor A
#'     (if \code{pool=FALSE}), or NULL if pooled.}
#'   \item{PSIHAT.4.EACH.LEVEL.OF.A}{List of estimated contrasts for each level of
#'     Factor A (if \code{pool=FALSE}), or NULL if pooled.}
#'   \item{POOLED.RESULTS}{Results when data pooled across Factor A (if \code{pool=TRUE}),
#'     or NULL otherwise.}
#'
#' @details
#' This is a median-based analog of \code{\link{bwbmcp}}. It compares levels of the
#' within-subjects factor (B) using medians. The FWE is controlled using Rom's method.
#'
#' When \code{pool=FALSE}, performs separate analyses for each level of Factor A.
#' When \code{pool=TRUE}, pools data across A levels and performs a single analysis.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bwbmcp}}, \code{\link{bwmedimcp}}, \code{\link{msmed}}
#' @examples
#' \dontrun{
#' # 2×3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(n)
#' # Pooled analysis
#' result <- bwmedbmcp(J=2, K=3, x=x, pool=TRUE)
#' print(result$POOLED.RESULTS)
#' }
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


#' Between-Within MCP with AKP Effect Sizes
#'
#' Computes Algina-Keselman-Penfield (AKP) robust effect sizes for all pairwise
#' comparisons in a between-by-within (split-plot) design. Provides effect sizes
#' for Factor A (between), Factor B (within), and interactions.
#'
#' @param J Number of levels for Factor A (between-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data in list mode or matrix with J×K groups.
#' @param tr Trim proportion for trimmed means (default: 0.2).
#' @param pr Logical. If TRUE (default), print descriptive messages about output structure.
#'
#' @return A list with components:
#'   \item{Factor.A}{List with J elements. \code{Factor.A[[j]]} contains AKP effect
#'     sizes for all pairwise comparisons among levels of Factor B at level j of
#'     Factor A (via \code{\link{wmcpAKP}}).}
#'   \item{Factor.B}{List with K elements. \code{Factor.B[[k]]} contains AKP effect
#'     sizes for all pairwise comparisons among levels of Factor A at level k of
#'     Factor B (via \code{\link{bmcpAKP}}).}
#'   \item{interactions}{Matrix of interaction effect sizes from \code{\link{bwimcpAKP}}.}
#'
#' @details
#' The AKP effect size is a robust analog of Cohen's d using trimmed means and
#' Winsorized standard deviations. This function provides a comprehensive set of
#' effect sizes for a split-plot design, computing them for:
#' \itemize{
#'   \item Within-subjects comparisons at each level of the between factor
#'   \item Between-subjects comparisons at each level of the within factor
#'   \item All interaction contrasts
#' }
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bwmcpQS}}, \code{\link{wmcpAKP}}, \code{\link{bmcpAKP}}, \code{\link{bwimcpAKP}}
#' @examples
#' \dontrun{
#' # 2×3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(n)
#' result <- bwmcpAKP(J=2, K=3, x=x)
#' print(result$Factor.A[[1]])  # Within-subjects ES at level 1 of A
#' }
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


#' Ordinal Between-Within MCP
#'
#' Performs ordinal (distribution-free) multiple comparisons in a between-by-within
#' (split-plot) design. Uses sign tests for within-subjects comparisons and Cliff's
#' analog for between-subjects comparisons. No parametric assumptions required.
#'
#' @param J Number of levels for Factor A (between-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data in list mode or matrix with J×K groups.
#' @param pr Logical. If TRUE (default), print descriptive messages about output structure.
#'
#' @return A list with components:
#'   \item{Factor.A}{List with J elements. \code{Factor.A[[j]]} contains sign test
#'     results for all pairwise comparisons among levels of Factor B at level j of
#'     Factor A.}
#'   \item{Factor.B}{List with K elements. \code{Factor.B[[k]]} contains Cliff's analog
#'     results for all pairwise comparisons among levels of Factor A at level k of
#'     Factor B.}
#'
#' @details
#' This function provides a completely distribution-free approach to MCP in split-plot
#' designs. For within-subjects comparisons (Factor B), it uses sign tests on paired
#' differences. For between-subjects comparisons (Factor A), it uses Cliff's analog,
#' which estimates P(X < Y).
#'
#' This method makes no assumptions about the shape of the distributions and is
#' particularly useful when data are ordinal or when parametric assumptions are violated.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bwmcpQS}}, \code{\link{bwmcpAKP}}, \code{\link{cid}}
#' @examples
#' \dontrun{
#' # 2×3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(n)
#' result <- bwmcpORD(J=2, K=3, x=x)
#' print(result$Factor.A[[1]])
#' }
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


#' Between-Within MCP with Quantile Shift Effect Sizes
#'
#' Computes quantile shift effect sizes for all pairwise comparisons in a
#' between-by-within (split-plot) design. The quantile shift measures the extent
#' to which distributions differ in location, providing robust effect sizes.
#'
#' @param J Number of levels for Factor A (between-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data in list mode or matrix with J×K groups.
#' @param locfun Location function to use for quantile shift (default: median).
#'   Can be any location estimator.
#' @param pr Logical. If TRUE (default), print descriptive messages and interpretive
#'   notes about effect sizes.
#'
#' @return A list with components:
#'   \item{Factor.A}{List with J elements. \code{Factor.A[[j]]} contains quantile shift
#'     effect sizes for all pairwise comparisons among levels of Factor B at level j
#'     of Factor A (via \code{\link{wmcpQS}}).}
#'   \item{Factor.B}{List with K elements. \code{Factor.B[[k]]} contains quantile shift
#'     effect sizes for all pairwise comparisons among levels of Factor A at level k
#'     of Factor B (via \code{\link{bmcpQS}}).}
#'   \item{interactions}{Vector of quantile shift effect sizes for interactions.}
#'
#' @details
#' The quantile shift effect size measures how much the distribution of one group
#' is shifted relative to another. It ranges from 0 to 1, with 0.5 indicating no
#' effect. Under normality and homoscedasticity, Cohen's d values of 0, 0.2, 0.5,
#' and 0.8 correspond approximately to quantile shifts of 0.50, 0.55, 0.65, and 0.70.
#'
#' This function provides a comprehensive set of robust effect sizes for a split-plot
#' design, covering within-subjects, between-subjects, and interaction contrasts.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bwmcpAKP}}, \code{\link{bwmcpORD}}, \code{\link{wmcpQS}}, \code{\link{bmcpQS}}, \code{\link{shiftes}}
#' @examples
#' \dontrun{
#' # 2×3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(i in 1:6) x[[i]] <- rnorm(n)
#' result <- bwmcpQS(J=2, K=3, x=x)
#' print(result$Factor.A[[1]])  # Within-subjects QS at level 1 of A
#' }
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

#' Multiple Comparisons for Between-Within Design (Percentile Bootstrap)
#'
#' Percentile bootstrap method for all pairwise comparisons of main effects and
#' interactions in a between-within two-way design. More flexible than \code{\link{bwmcp}}
#' as it allows any location estimator (not just trimmed means).
#'
#' @param J Number of levels for Factor A (between-subjects)
#' @param K Number of levels for Factor B (within-subjects)
#' @inheritParams bwmcp
#' @inheritParams linconpb
#' @inheritParams mcppb
#' @param bhop Logical; if \code{TRUE} (default), uses Benjamini-Hochberg procedure.
#'   If \code{FALSE}, uses Hochberg's method. Recommended \code{TRUE} for M-estimators with n < 80.
#' @inheritParams yuen
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{Fac.A}}{Matrix of results for Factor A main effect comparisons}
#'     \item{\code{Fac.B}}{Matrix of results for Factor B main effect comparisons}
#'     \item{\code{Fac.AB}}{Matrix of results for A×B interaction comparisons}
#'   }
#'   Each matrix contains: estimate, confidence interval, p-value, and adjusted p-value.
#'
#' @details
#' This function extends \code{\link{bwmcp}} by using percentile bootstrap instead of bootstrap-t:
#' \itemize{
#'   \item \strong{Method}: Percentile bootstrap with any robust estimator
#'   \item \strong{Default estimator}: Trimmed mean (\code{\link{tmean}})
#'   \item \strong{Alternatives}: Can use MOM (\code{\link{mom}}), M-estimator (\code{\link{onestep}}), median, etc.
#'   \item \strong{FWE control}: Benjamini-Hochberg or Hochberg adjustment
#'   \item \strong{Bootstrap}: Resamples within each level of Factor A (between-subjects)
#' }
#'
#' For Factor B (within-subjects) and interaction comparisons, the method accounts for
#' dependencies by resampling entire rows (subjects) within each level of Factor A.
#'
#' @note
#' \itemize{
#'   \item Set \code{bhop=TRUE} when using M-estimators with small samples (n < 80)
#'   \item Slower than \code{\link{bwmcp}} but more flexible (any estimator)
#'   \item Set \code{SEED=TRUE} (default) for reproducible results
#' }
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing} (5th ed.).
#' Academic Press.
#'
#' @seealso
#' \code{\link{bwmcp}} for bootstrap-t version,
#' \code{\link{bwwmcppb}} for three-way between-within-within designs,
#' \code{\link{linconpb}} for underlying percentile bootstrap method
#'
#' @family between-within MCP
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @export
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


#' Percentile Bootstrap MCP for Between-Within Designs (Internal Helper)
#'
#' Internal helper function for performing percentile bootstrap multiple comparisons
#' in between-within (mixed) factorial designs. Called by \code{bwmcppb} to analyze
#' specific contrasts.
#'
#' @param J Number of levels of Factor A (between-subjects factor)
#' @param K Number of levels of Factor B (within-subjects factor)
#' @param x Data in list mode or matrix. For list mode: x[[1]] = level (1,1),
#'   x[[2]] = level (1,2), ..., x[[K]] = level (1,K), x[[K+1]] = level (2,1), etc.
#' @param est Estimator function (default: \code{tmean} for trimmed mean)
#' @param JK Total number of groups (J*K)
#' @param con Contrast matrix for the comparisons of interest
#' @param method P-value adjustment method (default: 'hoch' for Hochberg)
#' @param alpha Family-wise error rate (default: 0.05)
#' @param grp Group indices to analyze (default: all groups 1:JK)
#' @param nboot Number of bootstrap samples (default: 500)
#' @param bhop Logical; if TRUE use Benjamini-Hochberg; if FALSE use specified method
#'   (default: TRUE)
#' @param SEED Logical; if TRUE sets random seed for reproducibility (default: TRUE)
#' @param ... Additional arguments passed to estimator function
#'
#' @return Matrix with columns: con.num, psihat (contrast estimate), p.value,
#'   adj.p.value, ci.lower, ci.upper
#'
#' @details
#' This is an internal subroutine used by \code{bwmcppb} to handle the bootstrap
#' analysis for a specific set of contrasts (main effects or interactions).
#'
#' The bootstrap procedure:
#' \itemize{
#'   \item Resamples subjects within each level of the between-subjects factor
#'   \item Computes the contrast estimates for each bootstrap sample
#'   \item Determines p-values and confidence intervals based on the bootstrap
#'     distribution
#' }
#'
#' @keywords internal
#' @seealso \code{\link{bwmcppb}} for the main user function
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


#' Percentile Bootstrap MCP for Between-Within Designs with Custom Adjustment (Internal)
#'
#' Internal helper that performs bootstrap MCP for between-within designs with
#' custom p-value adjustment applied to all contrasts together. Similar to
#' \code{bwmcppb.sub} but adjusts p-values across main effects and interactions.
#'
#' @inheritParams bwmcppb.sub
#' @param method P-value adjustment method applied globally (default: 'hoch')
#'
#' @return List with components Fac.A, Fac.B, Fac.AB, each containing a matrix
#'   of results with adjusted p-values
#'
#' @details
#' Differs from \code{bwmcppb.sub} by applying the p-value adjustment method
#' (e.g., Hochberg) to all contrasts together rather than separately for each
#' effect type.
#'
#' @keywords internal
#' @seealso \code{\link{bwmcppb.sub}}, \code{\link{bwmcppb}}
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


#' Multiple Comparisons for Between-Within-Within Design (Percentile Bootstrap)
#'
#' Percentile bootstrap method for all pairwise comparisons in a three-way
#' between-within-within design. Extends \code{\link{bwmcppb}} to three factors.
#'
#' @param J Number of levels for Factor A (between-subjects)
#' @param K Number of levels for Factor B (within-subjects)
#' @param L Number of levels for Factor C (within-subjects)
#' @param x Data in list mode or matrix with J×K×L groups
#' @inheritParams linconpb
#' @param JKL Total number of groups (J × K × L)
#' @inheritParams linconb
#' @inheritParams mcppb
#' @param bhop Logical; if \code{TRUE}, uses Benjamini-Hochberg procedure.
#'   If \code{FALSE} (default), uses Hochberg's method.
#' @inheritParams yuen
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{Fac.A}}{Results for Factor A main effect comparisons}
#'     \item{\code{Fac.B}}{Results for Factor B main effect comparisons}
#'     \item{\code{Fac.C}}{Results for Factor C main effect comparisons}
#'     \item{\code{Fac.AB}}{Results for A×B interaction}
#'     \item{\code{Fac.AC}}{Results for A×C interaction}
#'     \item{\code{Fac.BC}}{Results for B×C interaction}
#'     \item{\code{Fac.ABC}}{Results for A×B×C three-way interaction}
#'   }
#'   Each component contains estimates, confidence intervals, p-values, and adjusted p-values.
#'
#' @details
#' This function performs percentile bootstrap comparisons for a three-way design:
#' \itemize{
#'   \item \strong{Design}: One between-subjects factor (A), two within-subjects factors (B, C)
#'   \item \strong{Method}: Percentile bootstrap with any robust estimator
#'   \item \strong{Default estimator}: Trimmed mean, but can use MOM, M-estimator, median, etc.
#'   \item \strong{Bootstrap}: Resamples within each level of Factor A (between-subjects)
#'   \item \strong{FWE control}: Benjamini-Hochberg or Hochberg adjustment
#'   \item \strong{Default nboot}: 2000 (higher than two-way due to complexity)
#' }
#'
#' The method preserves dependencies for within-subjects factors by resampling entire
#' subjects within each level of the between-subjects factor.
#'
#' @note
#' \itemize{
#'   \item Tests 7 families of hypotheses (3 main effects, 3 two-way interactions, 1 three-way)
#'   \item Computationally intensive; consider using fewer bootstrap samples for initial exploration
#'   \item Set \code{SEED=TRUE} (default) for reproducibility
#' }
#'
#' @references
#' Wilcox, R.R. (2022). \emph{Introduction to Robust Estimation and Hypothesis Testing} (5th ed.).
#' Academic Press.
#'
#' @seealso
#' \code{\link{bwmcppb}} for two-way between-within design,
#' \code{\link{bwwmcp}} for bootstrap-t version,
#' \code{\link{bbwmcppb}} for between-between-within design
#'
#' @family between-within MCP
#' @family multiple comparison procedures
#' @family bootstrap methods
#' @export
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


#' Bootstrap MCP Helper for Between-Within-Within Designs (Internal)
#'
#' Internal helper function for \code{\link{bwwmcppb}} that performs bootstrap
#' analysis for specific contrasts in three-way between-within-within designs.
#'
#' @inheritParams bwmcppb.sub
#' @param L Number of levels for Factor C (third factor, within-subjects)
#' @param JKL Total number of groups (J × K × L)
#'
#' @return Matrix with contrast results including estimates, p-values, and
#'   confidence intervals
#'
#' @keywords internal
#' @seealso \code{\link{bwwmcppb}}
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


#' Between-By-Between Design Percentile Bootstrap MCP
#'
#' Performs percentile bootstrap multiple comparisons for a two-way between-subjects
#' design. Tests Factor A and B main effects and A×B interaction using any estimator.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data in list mode or matrix. Groups ordered as: (1,1), (1,2), ..., (1,K),
#'   (2,1), ..., (J,K).
#' @param est Measure of location (default: \code{tmean} for 20% trimmed mean).
#' @param JK Total number of groups, J×K (default: computed as J*K).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param grp Vector specifying group order (default: 1:JK).
#' @param nboot Number of bootstrap samples (default: 2000).
#' @param bhop Logical. If TRUE, use Benjamini-Hochberg method instead of Rom's (default: FALSE).
#' @param SEED Logical. If TRUE, set random seed for reproducibility (default: TRUE).
#' @param ... Additional arguments passed to \code{est}.
#'
#' @return A list with components:
#'   \item{n}{Sample sizes for each group.}
#'   \item{Fac.A}{Results for Factor A main effect.}
#'   \item{Fac.B}{Results for Factor B main effect.}
#'   \item{Fac.AB}{Results for A×B interaction.}
#'
#' @details
#' This function wraps \code{bbmcppb.sub} to test all main effects and interactions
#' in a two-way between-subjects design using percentile bootstrap. Rom's method
#' controls the family-wise error rate.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bwmcppb}}, \code{\link{bbbmcppb}}, \code{\link{con2way}}
#' @examples
#' \dontrun{
#' # 2x3 between-subjects design
#' set.seed(123)
#' x <- list(rnorm(15), rnorm(15), rnorm(15),
#'           rnorm(15, mean=0.5), rnorm(15), rnorm(15))
#' result <- bbmcppb(J=2, K=3, x=x, nboot=500)
#' print(result$Fac.A$output)
#' }
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


#' Bootstrap MCP for Between-Between-Within Three-Way Designs (Internal)
#'
#' Internal function for percentile bootstrap multiple comparisons in between-between-within
#' factorial designs (two between-subjects factors, one within-subjects factor).
#'
#' @inheritParams bwwmcppb
#'
#' @return List with 7 components: Fac.A, Fac.B, Fac.C (main effects), Fac.AB,
#'   Fac.AC, Fac.BC (two-way interactions), Fac.ABC (three-way interaction)
#'
#' @details
#' Design structure: Factors A and B are between-subjects, Factor C is within-subjects.
#' Bootstrap resamples independent observations for the between-subjects structure.
#'
#' @keywords internal
#' @seealso \code{\link{bwwmcppb}} for between-within-within design
#' @export
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


#' Bootstrap MCP Helper for Between-Between-Within Designs (Internal)
#'
#' Internal helper for \code{\link{bbwmcppb}} that performs bootstrap analysis
#' for specific contrasts in between-between-within three-way designs.
#'
#' @inheritParams bwwmcppb.sub
#'
#' @return Matrix with contrast results
#'
#' @keywords internal
#' @seealso \code{\link{bbwmcppb}}
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


#' Between-Between-Between Design Percentile Bootstrap MCP
#'
#' Performs percentile bootstrap multiple comparisons for a three-way between-subjects
#' design. Tests Factors A, B, and C main effects and all interactions using any estimator.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param L Number of levels for Factor C.
#' @param x Data in list mode or matrix. Groups ordered as: (1,1,1), (1,1,2), ...,
#'   (1,1,L), (1,2,1), ..., (J,K,L).
#' @param est Measure of location (default: \code{tmean} for 20% trimmed mean).
#' @param JKL Total number of groups, J×K×L (default: computed as J*K*L).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param grp Vector specifying group order (default: 1:JKL).
#' @param nboot Number of bootstrap samples (default: 2000).
#' @param bhop Logical. If TRUE, use Benjamini-Hochberg method instead of Rom's (default: FALSE).
#' @param SEED Logical. If TRUE, set random seed for reproducibility (default: TRUE).
#' @param ... Additional arguments passed to \code{est}.
#'
#' @return A list with components:
#'   \item{Fac.A}{Results for Factor A main effect.}
#'   \item{Fac.B}{Results for Factor B main effect.}
#'   \item{Fac.C}{Results for Factor C main effect.}
#'   \item{Fac.AB}{Results for A×B interaction.}
#'   \item{Fac.AC}{Results for A×C interaction.}
#'   \item{Fac.BC}{Results for B×C interaction.}
#'   \item{Fac.ABC}{Results for A×B×C three-way interaction.}
#'
#' @details
#' This function wraps \code{bbbmcppb.sub} to test all main effects and interactions
#' in a three-way between-subjects design using percentile bootstrap. Rom's method
#' controls the family-wise error rate.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{bbmcppb}}, \code{\link{wwwmcppb}}, \code{\link{con3way}}
#' @examples
#' \dontrun{
#' # 2x2x2 between-subjects design
#' set.seed(123)
#' x <- list(rnorm(15), rnorm(15), rnorm(15), rnorm(15),
#'           rnorm(15), rnorm(15), rnorm(15, mean=0.5), rnorm(15))
#' result <- bbbmcppb(J=2, K=2, L=2, x=x, nboot=500)
#' print(result$Fac.A$output)
#' }
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


#' Within-Within-Within Design Percentile Bootstrap MCP
#'
#' Performs percentile bootstrap multiple comparisons for a three-way repeated measures
#' (within-subjects) design. Tests Factors A, B, and C main effects and all interactions.
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param L Number of levels for Factor C (within-subjects).
#' @param x Data matrix or list. Columns/elements correspond to J×K×L repeated measures.
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param con Contrast matrix (default: 0 generates appropriate contrasts via \code{con3way}).
#' @param est Measure of location (default: \code{tmean} for 20% trimmed mean).
#' @param plotit Logical. If TRUE, create plots (default: FALSE).
#' @param dif Logical. If TRUE, use difference scores for dependent groups (default: TRUE,
#'   recommended for repeated measures).
#' @param grp Vector for reordering groups if needed (default: NA).
#' @param nboot Number of bootstrap samples (default: NA, auto-determined).
#' @param BA Logical. Bias adjustment parameter (default: TRUE).
#' @param hoch Logical. If TRUE, use Hochberg's method; if FALSE, use Rom's (default: TRUE).
#' @param xlab Label for x-axis in plots (default: "Group 1").
#' @param ylab Label for y-axis in plots (default: "Group 2").
#' @param pr Logical. If TRUE, print informational messages (default: TRUE).
#' @param SEED Logical. If TRUE, set random seed for reproducibility (default: TRUE).
#' @param ... Additional arguments passed to \code{est}.
#'
#' @return A list with components:
#'   \item{Factor_A}{Results for Factor A main effect.}
#'   \item{Factor_B}{Results for Factor B main effect.}
#'   \item{Factor_C}{Results for Factor C main effect (note: missing from current code).}
#'   \item{Factor_AB}{Results for A×B interaction.}
#'   \item{Factor_AC}{Results for A×C interaction.}
#'   \item{Factor_BC}{Results for B×C interaction.}
#'   \item{Factor_ABC}{Results for A×B×C three-way interaction.}
#'
#' @details
#' This function wraps \code{rmmcppb} to test all main effects and interactions
#' in a three-way repeated measures design. The use of difference scores (\code{dif=TRUE})
#' is recommended for within-subjects designs as it accounts for correlations.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{rmmcppb}}, \code{\link{wwmcppb}}, \code{\link{con3way}}
#' @examples
#' \dontrun{
#' # 2x2x2 repeated measures design
#' set.seed(123)
#' x <- matrix(rnorm(80), ncol=8)  # 10 subjects, 8 conditions
#' result <- wwwmcppb(J=2, K=2, L=2, x=x, nboot=500, pr=FALSE)
#' print(result$Factor_A$output)
#' }
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


#' Bootstrap MCP Helper for Within-Within-Within Designs (Internal)
#'
#' Internal helper for within-within-within repeated measures designs. Performs
#' bootstrap analysis for specific contrasts when all three factors are within-subjects.
#'
#' @inheritParams bwwmcppb.sub
#'
#' @return Matrix with contrast results
#'
#' @details
#' All three factors are within-subjects (repeated measures). Bootstrap resamples
#' entire subjects to preserve the dependency structure across all factors.
#'
#' @keywords internal
#' @seealso \code{\link{wwwmcppb}}, \code{\link{wwwmcppbtr}}
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

#' Bootstrap MCP for Within-Within-Within Designs with Trimmed Means (Internal)
#'
#' Internal function for percentile bootstrap multiple comparisons in within-within-within
#' (fully repeated measures) designs using trimmed means.
#'
#' @param J Number of levels for Factor A (within-subjects)
#' @param K Number of levels for Factor B (within-subjects)
#' @param L Number of levels for Factor C (within-subjects)
#' @param x Data in matrix or list mode
#' @param tr Trimming proportion (default: 0.2)
#' @param alpha Family-wise error rate (default: 0.05)
#' @param dif Logical; use difference scores (default: TRUE)
#' @param op Logical; produce plots (default: FALSE)
#' @param grp Group indices (default: NA for all groups)
#' @param nboot Number of bootstrap samples (default: 2000)
#' @param SEED Logical; set random seed (default: TRUE)
#' @param pr Logical; print messages (default: TRUE)
#'
#' @return List with results for all main effects, two-way interactions, and
#'   three-way interaction
#'
#' @details
#' Wrapper function that calls \code{\link{rmmcppb}} for each effect type
#' (main effects, interactions) in a within-within-within design using trimmed means.
#'
#' @keywords internal
#' @seealso \code{\link{wwwmcppb.sub}}, \code{\link{rmmcppb}}
#' @export
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

#' Repeated Measures Percentile Bootstrap Multiple Comparisons
#'
#' Performs percentile bootstrap multiple comparisons for a one-way repeated measures
#' (within-subjects) design. Wrapper for \code{\link{rmmcppb}}.
#'
#' @param x Data matrix or data frame. Each column represents a repeated measure.
#' @param y Optional second variable to bind with x (default: NULL).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param con Contrast matrix (default: 0 generates all pairwise comparisons).
#' @param est Measure of location (default: \code{tmean} for 20% trimmed mean).
#' @param plotit Logical. If TRUE, create plots (default: FALSE).
#' @param dif Logical. If TRUE, use difference scores for dependent groups (default: TRUE,
#'   recommended for repeated measures).
#' @param grp Vector for reordering groups if needed (default: NA).
#' @param nboot Number of bootstrap samples (default: NA, auto-determined).
#' @param BA Logical. Bias adjustment parameter (default: TRUE).
#' @param hoch Logical. If TRUE, use Hochberg's method; if FALSE, use Rom's (default: TRUE).
#' @param xlab Label for x-axis in plots (default: "Group 1").
#' @param ylab Label for y-axis in plots (default: "Group 2").
#' @param pr Logical. If TRUE, print informational messages (default: TRUE).
#' @param SEED Logical. If TRUE, set random seed for reproducibility (default: TRUE).
#' @param ... Additional arguments passed to \code{est}.
#'
#' @return Results from \code{rmmcppb}, a list with components including:
#'   \item{output}{Matrix with test results for each contrast.}
#'   \item{con}{Contrast matrix used.}
#'   \item{num.sig}{Number of significant contrasts.}
#'
#' @details
#' This is a simple wrapper for \code{rmmcppb}. The use of difference scores
#' (\code{dif=TRUE}) is recommended for within-subjects designs as it properly
#' accounts for correlations between repeated measures.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{rmmcppb}}, \code{\link{wwmcppb}}, \code{\link{rmmcp}}
#' @examples
#' \dontrun{
#' # One-way repeated measures with 4 conditions
#' set.seed(123)
#' x <- matrix(rnorm(40), ncol=4)  # 10 subjects, 4 conditions
#' result <- wmcppb(x, nboot=500, pr=FALSE)
#' print(result$output)
#' }
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


#' Within-By-Within Design Percentile Bootstrap MCP
#'
#' Performs percentile bootstrap multiple comparisons for a two-way repeated measures
#' (within-by-within) design. Tests Factor A and B main effects and A×B interaction.
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data matrix or list. Columns/elements correspond to J×K repeated measures.
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param con Contrast matrix (default: 0 generates appropriate contrasts via \code{con2way}).
#' @param est Measure of location (default: \code{tmean} for 20% trimmed mean).
#' @param plotit Logical. If TRUE, create plots (default: FALSE).
#' @param dif Logical. If TRUE, use difference scores for dependent groups (default: TRUE,
#'   recommended for repeated measures).
#' @param grp Vector for reordering groups if needed (default: NA).
#' @param nboot Number of bootstrap samples (default: NA, auto-determined).
#' @param BA Logical. Bias adjustment parameter (default: TRUE).
#' @param hoch Logical. If TRUE, use Hochberg's method; if FALSE, use Rom's (default: TRUE).
#' @param xlab Label for x-axis in plots (default: "Group 1").
#' @param ylab Label for y-axis in plots (default: "Group 2").
#' @param pr Logical. If TRUE, print informational messages (default: TRUE).
#' @param SEED Logical. If TRUE, set random seed for reproducibility (default: TRUE).
#' @param ... Additional arguments passed to \code{est}.
#'
#' @return A list with components:
#'   \item{Factor_A}{Results for Factor A main effect.}
#'   \item{Factor_B}{Results for Factor B main effect.}
#'   \item{Factor_AB}{Results for A×B interaction.}
#'
#' @details
#' This function wraps \code{rmmcppb} to test all main effects and interactions
#' in a two-way repeated measures design. The use of difference scores (\code{dif=TRUE})
#' is recommended for within-subjects designs as it accounts for correlations.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{rmmcppb}}, \code{\link{wmcppb}}, \code{\link{wwwmcppb}}, \code{\link{con2way}}
#' @examples
#' \dontrun{
#' # 3x2 repeated measures design
#' set.seed(123)
#' x <- matrix(rnorm(60), ncol=6)  # 10 subjects, 6 conditions (3x2)
#' result <- wwmcppb(J=3, K=2, x=x, nboot=500, pr=FALSE)
#' print(result$Factor_A$output)
#' }
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


#' Within-By-Within Design Bootstrap-t MCP
#'
#' Performs bootstrap-t multiple comparisons for a two-way repeated measures
#' (within-by-within) design using trimmed means. Tests Factor A and B main
#' effects and A×B interaction.
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data matrix or list. Columns/elements correspond to J×K repeated measures.
#' @param tr Trim proportion (default: 0.2 for 20% trimming).
#' @param dif Logical. If TRUE (default), use difference scores for dependent groups
#'   (recommended for repeated measures).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param nboot Number of bootstrap samples (default: 599).
#'
#' @return A list with components:
#'   \item{Factor_A}{Results from \code{\link{lindepbt}} for Factor A main effect.}
#'   \item{Factor_B}{Results from \code{\link{lindepbt}} for Factor B main effect.}
#'   \item{Factor_AB}{Results from \code{\link{lindepbt}} for A×B interaction.}
#'
#' @details
#' This function uses the bootstrap-t method via \code{\link{lindepbt}} to test
#' all main effects and interactions in a two-way repeated measures design.
#' Contrast matrices are generated using \code{\link{con2way}}.
#'
#' The bootstrap-t method can provide better control of Type I error rates than
#' percentile bootstrap methods in some situations. The use of difference scores
#' (\code{dif=TRUE}) accounts for correlations between repeated measures.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{lindepbt}}, \code{\link{wwmcppb}}, \code{\link{rmmcp}}, \code{\link{con2way}}
#' @examples
#' \dontrun{
#' # 2×3 repeated measures design
#' set.seed(123)
#' x <- matrix(rnorm(60), ncol=6)  # 10 subjects, 6 conditions (2×3)
#' result <- wwmcpbt(J=2, K=3, x=x, nboot=500)
#' print(result$Factor_A$output)
#' }
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


#' Within-By-Within Design Effect Size MCP (Deprecated)
#'
#' **Deprecated**: This function has been replaced by \code{ww.es}. Please use
#' that function instead.
#'
#' @param J Number of levels for Factor A.
#' @param K Number of levels for Factor B.
#' @param x Data matrix or list.
#' @param tr Trim proportion (default: 0.2).
#' @param alpha Significance level (default: 0.05).
#' @param dif Logical. Use difference scores (default: TRUE).
#'
#' @return Throws an error directing users to \code{ww.es}.
#'
#' @details
#' This function is deprecated. Use \code{ww.es} for within-by-within
#' effect size comparisons.
#'
#' @keywords internal
#' @export
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


#' Within-By-Within Design Quantile Shift MCP
#'
#' Performs multiple comparisons for a two-way repeated measures (within-by-within)
#' design with quantile shift effect sizes. Tests Factor A and B main effects and
#' A×B interaction.
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data matrix or list. Columns/elements correspond to J×K repeated measures.
#' @param tr Trim proportion (default: 0.2, used if \code{locfun=tmean}).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param dif Logical. If TRUE (default), use difference scores for dependent groups
#'   (recommended for repeated measures).
#' @param locfun Location function for quantile shift (default: \code{tmean}).
#'   Can be any location estimator like \code{median}, \code{mean}, etc.
#' @param ... Additional arguments passed to \code{locfun}.
#'
#' @return A list with components:
#'   \item{Factor_A}{Results from \code{\link{rmmcpQS}} for Factor A main effect,
#'     including p-values and quantile shift effect sizes.}
#'   \item{Factor_B}{Results from \code{\link{rmmcpQS}} for Factor B main effect.}
#'   \item{Factor_AB}{Results from \code{\link{rmmcpQS}} for A×B interaction.}
#'
#' @details
#' This function uses \code{\link{rmmcpQS}} to test all main effects and interactions
#' in a two-way repeated measures design, providing quantile shift effect sizes for
#' each contrast. Contrast matrices are generated using \code{\link{con2way}}.
#'
#' Quantile shift effect sizes measure how much one distribution is shifted relative
#' to another, ranging from 0 to 1 with 0.5 indicating no effect. Under normality,
#' Cohen's d = 0, 0.2, 0.5, 0.8 correspond approximately to Q = 0.50, 0.55, 0.65, 0.70.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{rmmcpQS}}, \code{\link{wwmcppb}}, \code{\link{wwwmcpQS}}, \code{\link{shiftes}}
#' @examples
#' \dontrun{
#' # 2×3 repeated measures design
#' set.seed(123)
#' x <- matrix(rnorm(60), ncol=6)  # 10 subjects, 6 conditions (2×3)
#' result <- wwmcpQS(J=2, K=3, x=x)
#' print(result$Factor_A)
#' }
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


#' Three-Way Within-Subjects Design Quantile Shift MCP
#'
#' Performs multiple comparisons for a three-way repeated measures (within-by-within-by-within)
#' design with quantile shift effect sizes. Tests all main effects and interactions for
#' Factors A, B, and C.
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param L Number of levels for Factor C (within-subjects).
#' @param x Data matrix or list. Columns/elements correspond to J×K×L repeated measures.
#' @param tr Trim proportion (default: 0.2, used if \code{locfun=tmean}).
#' @param alpha Family-wise Type I error rate (default: 0.05).
#' @param dif Logical. If TRUE (default), use difference scores for dependent groups.
#' @param op Logical. Print operation messages (default: FALSE).
#' @param grp Optional vector to reorder groups if data not in expected order (default: NA).
#' @param locfun Location function for quantile shift (default: \code{tmean}).
#' @param ... Additional arguments passed to \code{locfun}.
#'
#' @return A list with components:
#'   \item{Factor.A, Factor.B, Factor.C}{Results from \code{\link{rmmcpQS}} for main effects.}
#'   \item{Factor.AB, Factor.AC, Factor.BC}{Results for two-way interactions.}
#'   \item{Factor.ABC}{Results for three-way interaction.}
#'   \item{conA, conB, conC, conAB, conAC, conBC, conABC}{Contrast matrices used.}
#'
#' @details
#' This function performs all pairwise comparisons for main effects and interactions
#' in a three-way repeated measures design using quantile shift effect sizes.
#' Contrast matrices are generated using \code{\link{con3way}}.
#'
#' Quantile shift effect sizes measure distributional shift, ranging from 0 to 1
#' with 0.5 indicating no effect. Under normality, Cohen's d values of 0, 0.2, 0.5, 0.8
#' correspond approximately to Q values of 0.50, 0.55, 0.65, 0.70.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{rmmcpQS}}, \code{\link{wwmcpQS}}, \code{\link{rm3mcp}}, \code{\link{wwwmcppb}}
#' @examples
#' \dontrun{
#' # 2×2×2 repeated measures design
#' set.seed(123)
#' x <- matrix(rnorm(80), ncol=8)  # 10 subjects, 8 conditions (2×2×2)
#' result <- wwwmcpQS(J=2, K=2, L=2, x=x)
#' print(result$Factor.A)
#' print(result$Factor.ABC)
#' }
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

#' Multiple Comparisons for Split-Plot Designs - Main Effects (Factor A)
#'
#' Performs percentile bootstrap multiple comparisons among all pairwise main
#' effects for Factor A (between-subjects factor) in a split-plot design.
#' Uses appropriate linear contrasts with Rom's method for FWE control.
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (within-subjects factor).
#' @param x Data in matrix or list format. If a matrix, columns represent groups.
#'   If a list, \code{x[[1]]} contains data for level (1,1), \code{x[[2]]} for
#'   level (1,2), ..., \code{x[[K]]} for level (1,K), \code{x[[K+1]]} for
#'   level (2,1), etc. Length should be J*K.
#' @param est Measure of location (default: \code{\link{tmean}}).
#' @param JK Total number of groups (default: J*K).
#' @param grp Vector indicating which groups to include (default: all J*K groups).
#' @param avg Logical. If TRUE, averages across Factor B levels before comparing
#'   Factor A levels. If FALSE, uses all J*K groups separately (default: FALSE).
#' @param alpha Family-wise error rate (default: 0.05).
#' @param SEED Logical. If TRUE, sets random seed for reproducibility (default: TRUE).
#' @param nboot Number of bootstrap samples. If NA, defaults based on number of
#'   contrasts: 1000 if ≤4 contrasts, 5000 if >4 (default: NA).
#' @param pr Logical. Print informational messages (default: TRUE).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat, p.value, p.sig,
#'     ci.lower, ci.upper}
#'   \item{con}{Contrast matrix used}
#'   \item{num.sig}{Number of significant comparisons}
#'
#' @details
#' For split-plot (mixed) designs with Factor A (between-subjects) and Factor B
#' (within-subjects), this function performs all pairwise comparisons among the
#' J levels of Factor A. The bootstrap accounts for the correlation structure
#' within subjects. Rom's method controls the family-wise error rate.
#'
#' When \code{avg=TRUE}, the function averages the K repeated measures for each
#' subject before comparing Factor A levels. When \code{avg=FALSE}, all J*K cell
#' means are used with appropriate contrasts.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{spmcpi}}, \code{\link{spmcpb}}, \code{\link{bwmcp}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x3 split-plot design (2 groups, 3 time points)
#' set.seed(123)
#' J <- 2  # Between-subjects factor
#' K <- 3  # Within-subjects factor
#' n <- 20
#' x <- list()
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     x[[(j-1)*K + k]] <- rnorm(n, mean = j + k)
#'   }
#' }
#'
#' # Compare main effects of Factor A
#' result <- spmcpa(J, K, x, nboot = 500)
#' result$num.sig
#' }
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

#' Multiple Comparisons for Split-Plot Designs - Interactions
#'
#' Performs percentile bootstrap multiple comparisons for interaction contrasts
#' in a split-plot design. Tests which pairwise differences among dependent groups
#' (Factor B) differ across levels of Factor A (between-subjects factor).
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (within-subjects factor).
#' @param x Data in matrix or list format. If a matrix, columns represent groups.
#'   If a list, \code{x[[1]]} contains data for level (1,1), \code{x[[2]]} for
#'   level (1,2), ..., \code{x[[K]]} for level (1,K), \code{x[[K+1]]} for
#'   level (2,1), etc. Length should be J*K.
#' @param est Measure of location (default: \code{\link{tmean}}).
#' @param JK Total number of groups (default: J*K).
#' @param grp Vector indicating which groups to include (default: all J*K groups).
#' @param alpha Family-wise error rate (default: 0.05).
#' @param nboot Number of bootstrap samples. If NA, defaults based on number of
#'   contrasts: 1000 if ≤4 contrasts, 5000 if >4 (default: NA).
#' @param SEED Logical. If TRUE, sets random seed for reproducibility (default: TRUE).
#' @param pr Logical. Print informational messages (default: TRUE).
#' @param SR Logical. Internal parameter (default: FALSE).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat, p.value, p.sig,
#'     ci.lower, ci.upper}
#'   \item{con}{Contrast matrix used for interaction tests}
#'   \item{num.sig}{Number of significant interaction contrasts}
#'
#' @details
#' For split-plot (mixed) designs, this function tests interactions by computing
#' difference scores among all pairs of dependent groups (Factor B levels) and
#' determining which of these differences differ across levels of Factor A.
#'
#' The analysis uses a percentile bootstrap that properly accounts for the
#' correlation structure of repeated measures. Rom's method controls the
#' family-wise error rate across all interaction contrasts.
#'
#' Number of contrasts tested: J(J-1)/2 × K(K-1)/2 (all interaction contrasts).
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{spmcpa}}, \code{\link{spmcpb}}, \code{\link{bwimcp}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     # Create interaction effect
#'     x[[(j-1)*K + k]] <- rnorm(n, mean = j*k)
#'   }
#' }
#'
#' # Test for interaction contrasts
#' result <- spmcpi(J, K, x, nboot = 500)
#' result$num.sig
#' }
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

#' Multiple Comparisons for Split-Plot Designs - Within-Subjects Pairwise
#'
#' Performs percentile bootstrap pairwise multiple comparisons among all
#' dependent groups (Factor B levels) in a split-plot design, pooling across
#' levels of Factor A.
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (within-subjects factor).
#' @param x Data in matrix or list format. If a matrix, columns represent groups.
#'   If a list, \code{x[[1]]} contains data for level (1,1), \code{x[[2]]} for
#'   level (1,2), ..., \code{x[[K]]} for level (1,K), \code{x[[K+1]]} for
#'   level (2,1), etc. Length should be J*K.
#' @param est Measure of location (default: \code{\link{tmean}}).
#' @param JK Total number of groups (default: J*K).
#' @param grp Vector indicating which groups to include (default: all J*K groups).
#' @param dif Logical. If TRUE, uses difference scores. If FALSE, uses marginal
#'   measures of location (default: TRUE).
#' @param alpha Family-wise error rate (default: 0.05).
#' @param SEED Logical. If TRUE, sets random seed for reproducibility (default: TRUE).
#' @param nboot Number of bootstrap samples (default: NA, determined automatically).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: con.num, psihat, p.value, p.sig,
#'     ci.lower, ci.upper}
#'   \item{con}{Contrast matrix used}
#'   \item{num.sig}{Number of significant comparisons}
#'
#' @details
#' This function performs all pairwise comparisons among the K levels of Factor B
#' (within-subjects factor), pooling data across all J levels of Factor A.
#' Levels of Factor A are ignored in the analysis.
#'
#' When \code{dif=TRUE}, the analysis is based on all pairs of difference scores,
#' which is generally more powerful. When \code{dif=FALSE}, marginal measures of
#' location are compared directly.
#'
#' The function internally calls \code{\link{rmmcppb}} after pooling the data
#' appropriately. Rom's method controls the family-wise error rate.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{spmcpa}}, \code{\link{spmcpi}}, \code{\link{spmcpbA}},
#'   \code{\link{rmmcppb}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     x[[(j-1)*K + k]] <- rnorm(n, mean = k)
#'   }
#' }
#'
#' # Compare all pairs of Factor B levels (pooled across Factor A)
#' result <- spmcpb(J, K, x, nboot = 500)
#' result$num.sig
#' }
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

#' Multiple Comparisons for Split-Plot Designs - Within Each Level of Factor A
#'
#' Performs percentile bootstrap pairwise multiple comparisons among dependent
#' groups (Factor B levels) separately for each level of Factor A in a
#' split-plot design.
#'
#' @param J Number of levels for Factor A (between-subjects factor).
#' @param K Number of levels for Factor B (within-subjects factor).
#' @param x Data in matrix or list format. If a matrix, columns represent groups.
#'   If a list, \code{x[[1]]} contains data for level (1,1), \code{x[[2]]} for
#'   level (1,2), ..., \code{x[[K]]} for level (1,K), \code{x[[K+1]]} for
#'   level (2,1), etc. Length should be J*K.
#' @param est Measure of location (default: \code{\link{tmean}}).
#' @param JK Total number of groups (default: J*K).
#' @param grp Vector indicating which groups to include (default: all J*K groups).
#' @param dif Logical. If TRUE, uses difference scores. If FALSE, uses marginal
#'   measures of location (default: TRUE).
#' @param alpha Family-wise error rate (default: 0.05).
#' @param nboot Number of bootstrap samples (default: NA, determined automatically).
#' @param SEED Logical. If TRUE, sets random seed for reproducibility (default: TRUE).
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return A list containing J sublists (one per level of Factor A), each with:
#'   \item{output}{Matrix with columns: con.num, psihat, p.value, p.sig,
#'     ci.lower, ci.upper}
#'   \item{con}{Contrast matrix used}
#'   \item{num.sig}{Number of significant comparisons}
#'
#' @details
#' Unlike \code{\link{spmcpb}}, which pools across all levels of Factor A, this
#' function performs separate analyses for each level of Factor A. For each of
#' the J levels of Factor A, all pairwise comparisons among the K levels of
#' Factor B are tested.
#'
#' This approach is useful when you want to examine the pattern of differences
#' among repeated measures separately for each between-subjects group.
#'
#' The function calls \code{\link{rmmcppb}} separately for each level of Factor A.
#' Rom's method controls the family-wise error rate within each set of comparisons.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{spmcpb}}, \code{\link{spmcpa}}, \code{\link{spmcpi}},
#'   \code{\link{rmmcppb}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x3 split-plot design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 20
#' x <- list()
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     # Different patterns for each group
#'     x[[(j-1)*K + k]] <- rnorm(n, mean = j + k)
#'   }
#' }
#'
#' # Compare Factor B levels separately within each level of Factor A
#' result <- spmcpbA(J, K, x, nboot = 500)
#' # result is a list with J elements
#' result[[1]]$num.sig  # Significant comparisons in first group
#' result[[2]]$num.sig  # Significant comparisons in second group
#' }
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

#' Multiple Comparisons Using Quantiles - Dependent Groups
#'
#' Performs all pairwise comparisons among dependent groups using quantile
#' estimates. By default compares medians. Family-wise error rate controlled
#' with Rom's method.
#'
#' @param x Data matrix (rows are subjects, columns are groups) or list of vectors.
#' @param alpha Family-wise error rate (default: 0.05).
#' @param bop Logical. If TRUE, uses usual median with bootstrap estimate of
#'   standard error. If FALSE, uses single order statistic estimator (default: FALSE).
#' @param nboot Number of bootstrap samples when \code{bop=TRUE} (default: 100).
#' @param pr Logical. Print progress messages (default: TRUE).
#' @param q Quantile to estimate (default: 0.5 for median).
#' @param SEED Logical. If TRUE, sets random seed for small samples (default: TRUE).
#'
#' @return A list with components:
#'   \item{test}{Matrix with columns: Group, Group, test, p-value, p.crit, se}
#'   \item{psihat}{Matrix with columns: Group, Group, psihat, ci.lower, ci.upper}
#'   \item{num.sig}{Number of significant pairwise comparisons}
#'
#' @details
#' For J dependent groups, this function performs all J(J-1)/2 pairwise
#' comparisons of quantiles. Rom's method is used to control the family-wise
#' error rate across all comparisons.
#'
#' When \code{bop=FALSE} (default), quantiles are estimated using a single order
#' statistic, which provides distribution-free inference. When \code{bop=TRUE},
#' the usual sample median is used with a bootstrap estimate of the covariance
#' matrix among medians.
#'
#' For small samples (n < 20), the function uses \code{\link{mrm1way}} to
#' compute p-values via permutation tests for better control of Type I error.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{qdmcpdif}}, \code{\link{mrm1way}}, \code{\link{qest}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare medians across 4 time points
#' set.seed(123)
#' J <- 4
#' n <- 25
#' x <- matrix(rnorm(n*J), ncol = J)
#' # Add differences to some time points
#' x[,3] <- x[,3] + 0.5
#' x[,4] <- x[,4] + 1
#'
#' # All pairwise median comparisons
#' result <- qdmcp(x)
#' result$test
#' result$num.sig
#' }
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

#' Multiple Comparisons Using Medians on Difference Scores
#'
#' Performs multiple comparisons using medians computed on difference scores
#' for dependent groups. Family-wise error rate controlled with Rom's method.
#'
#' @param x Data matrix (rows are subjects, columns are groups) or list of vectors.
#' @param con Contrast matrix (default: 0 for all pairwise comparisons). Each column
#'   specifies one contrast.
#' @param alpha Family-wise error rate (default: 0.05).
#'
#' @return A list with components:
#'   \item{test}{Matrix with columns: Group, Group, p-value, p.crit, se (for pairwise)}
#'   \item{psihat}{Matrix with columns: Group, Group, psihat, ci.lower, ci.upper (for pairwise)}
#'   \item{con.results}{Results for custom contrasts (if con specified)}
#'   \item{num.sig}{Number of significant comparisons}
#'
#' @details
#' This function performs multiple comparisons on medians of difference scores.
#' For pairwise comparisons (default when \code{con=0}), all pairs of difference
#' scores are formed and their medians compared using the sign test approach via
#' \code{\link{sintv2}}.
#'
#' For custom contrasts, the function forms weighted difference scores according
#' to the contrast coefficients and tests whether the median differs from zero.
#'
#' Rom's method controls the family-wise error rate. Standard errors are computed
#' using \code{\link{msmedse}} (McKean-Schrader estimate).
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{qdmcp}}, \code{\link{sintv2}}, \code{\link{msmedse}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare medians of difference scores across 3 time points
#' set.seed(123)
#' n <- 30
#' x <- matrix(rnorm(n*3), ncol = 3)
#' x[,2] <- x[,2] + 0.5  # Add effect at time 2
#'
#' # All pairwise comparisons
#' result <- qdmcpdif(x)
#' result$num.sig
#' }
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

#' Multiple Comparisons for Within-Within Designs Using Quantiles
#'
#' Performs all multiple comparisons for main effects and interactions in a
#' J×K repeated measures design using quantiles (default: medians). Both
#' factors are within-subjects (dependent measures).
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data in matrix (rows=subjects, columns=J*K groups) or list format.
#'   If a list, \code{x[[1]]} is level (1,1), \code{x[[2]]} is level (1,2), ...,
#'   \code{x[[K]]} is level (1,K), \code{x[[K+1]]} is level (2,1), etc.
#' @param grp Vector specifying which groups to include (default: all J*K groups).
#' @param p Total number of groups (default: J*K).
#' @param q Quantile to compare (default: 0.5 for median).
#' @param bop Logical. If TRUE, uses usual median with bootstrap covariance
#'   estimate. If FALSE, uses single order statistic (default: FALSE).
#' @param alpha Family-wise error rate (default: 0.05).
#' @param nboot Number of bootstrap samples when \code{bop=TRUE} (default: 100).
#' @param SEED Logical. If TRUE, sets random seed (default: TRUE).
#'
#' @return A list with components:
#'   \item{Qa}{Results for Factor A main effect comparisons}
#'   \item{Qb}{Results for Factor B main effect comparisons}
#'   \item{Qab}{Results for A×B interaction comparisons}
#'
#' @details
#' This function tests three sets of hypotheses in a J×K within-within design:
#' \itemize{
#'   \item All pairwise comparisons for Factor A main effects: J(J-1)/2 tests
#'   \item All pairwise comparisons for Factor B main effects: K(K-1)/2 tests
#'   \item All interaction contrasts: J(J-1)/2 × K(K-1)/2 tests
#' }
#'
#' Uses \code{\link{lincdm}} internally to perform the tests. Rom's method
#' controls the family-wise error rate within each family of tests.
#'
#' When \code{bop=FALSE}, uses distribution-free inference via single order
#' statistics. When \code{bop=TRUE}, uses bootstrap estimate of the covariance
#' matrix among sample medians.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{twwmcp}}, \code{\link{lincdm}}, \code{\link{qdmcp}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x3 within-within design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 25
#' x <- matrix(rnorm(n * J * K), ncol = J*K)
#' # Add main effects and interaction
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     col <- (j-1)*K + k
#'     x[,col] <- x[,col] + j + k + 0.5*j*k
#'   }
#' }
#'
#' result <- mwwmcp(J, K, x, nboot = 100)
#' result$Qa$num.sig  # Factor A comparisons
#' result$Qb$num.sig  # Factor B comparisons
#' result$Qab$num.sig  # Interaction comparisons
#' }
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

#' Multiple Comparisons for Within-Within Designs Using Trimmed Means
#'
#' Performs all multiple comparisons for main effects and interactions in a
#' J×K repeated measures design using trimmed means. Both factors are
#' within-subjects (dependent measures).
#'
#' @param J Number of levels for Factor A (within-subjects).
#' @param K Number of levels for Factor B (within-subjects).
#' @param x Data in matrix (rows=subjects, columns=J*K groups) or list format.
#'   If a list, \code{x[[1]]} is level (1,1), \code{x[[2]]} is level (1,2), ...,
#'   \code{x[[K]]} is level (1,K), \code{x[[K+1]]} is level (2,1), etc.
#' @param grp Vector specifying which groups to include (default: all J*K groups).
#' @param p Total number of groups (default: J*K).
#' @param tr Proportion to trim (default: 0.2 for 20% trimming).
#' @param alpha Family-wise error rate (default: 0.05).
#' @param dif Logical. If TRUE, uses difference scores for main effects. If FALSE,
#'   uses marginal trimmed means (default: FALSE).
#'
#' @return A list with components:
#'   \item{Qa}{Results for Factor A main effect comparisons}
#'   \item{Qb}{Results for Factor B main effect comparisons}
#'   \item{Qab}{Results for A×B interaction comparisons}
#'
#' @details
#' This function tests three sets of hypotheses in a J×K within-within design
#' using trimmed means:
#' \itemize{
#'   \item All pairwise comparisons for Factor A main effects: J(J-1)/2 tests
#'   \item All pairwise comparisons for Factor B main effects: K(K-1)/2 tests
#'   \item All interaction contrasts: J(J-1)/2 × K(K-1)/2 tests
#' }
#'
#' Uses \code{\link{lindepbt}} internally to perform the trimmed mean comparisons
#' with bootstrap-t method. Rom's method or Hochberg's method controls the
#' family-wise error rate within each family of tests.
#'
#' When \code{dif=TRUE}, main effect comparisons are based on difference scores,
#' which can be more powerful. When \code{dif=FALSE}, marginal trimmed means
#' are compared.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @seealso \code{\link{mwwmcp}}, \code{\link{lindepbt}}, \code{\link{rmmcp}}
#'
#' @export
#' @examples
#' \dontrun{
#' # 2x3 within-within design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' n <- 25
#' x <- matrix(rnorm(n * J * K), ncol = J*K)
#' # Add effects
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     col <- (j-1)*K + k
#'     x[,col] <- x[,col] + j*0.5 + k*0.3
#'   }
#' }
#'
#' result <- twwmcp(J, K, x)
#' result$Qa$num.sig  # Factor A comparisons
#' result$Qb$num.sig  # Factor B comparisons
#' result$Qab$num.sig  # Interaction comparisons
#' }
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

#' Tukey-Kramer Multiple Comparison Procedure
#'
#' Performs the conventional Tukey-Kramer multiple comparison procedure for
#' all pairwise comparisons among independent groups. Uses the studentized
#' range distribution to control family-wise error rate.
#'
#' @param x Data in matrix format (groups in columns) or list format.
#' @param alpha Family-wise error rate (default: 0.05).
#' @param ind.pval Logical. If TRUE, computes individual p-values for each test.
#'   If FALSE, computes p-values based on controlling the family-wise error rate
#'   using the studentized range distribution (default: TRUE).
#'
#' @return Matrix with columns:
#'   \item{Group}{First group number}
#'   \item{Group}{Second group number}
#'   \item{t.test}{Test statistic (studentized range statistic)}
#'   \item{est.difference}{Estimated mean difference}
#'   \item{ci.lower}{Lower confidence limit}
#'   \item{ci.upper}{Upper confidence limit}
#'   \item{p.value}{P-value (individual or FWE-adjusted based on \code{ind.pval})}
#'
#' @details
#' This is the standard parametric Tukey-Kramer procedure assuming normality
#' and homogeneity of variance. It uses pooled variance from one-way ANOVA and
#' the studentized range distribution to construct simultaneous confidence
#' intervals for all pairwise comparisons.
#'
#' The procedure controls the family-wise error rate at level alpha across all
#' J(J-1)/2 pairwise comparisons.
#'
#' For robust alternatives, see \code{\link{linconb}}, \code{\link{linconpb}},
#' or \code{\link{mcppb}}.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' Tukey, J. W. (1953). The problem of multiple comparisons. Unpublished manuscript.
#'
#' @seealso \code{\link{linconb}}, \code{\link{mcppb}}, \code{\link{pairwise.t.test}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare 4 independent groups
#' set.seed(123)
#' x <- list(
#'   rnorm(20, mean = 0),
#'   rnorm(20, mean = 0.5),
#'   rnorm(20, mean = 1),
#'   rnorm(20, mean = 1.5)
#' )
#'
#' # All pairwise comparisons
#' result <- tkmcp(x, alpha = 0.05)
#' result
#' }
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

#' Effect Size Multiple Comparisons for One-Way Designs
#'
#' Computes heteroscedastic robust effect sizes for all pairwise comparisons
#' in a one-way independent groups design. Supports six different effect size
#' measures.
#'
#' @param x Data in list or matrix format. If list, \code{x[[j]]} contains data
#'   for group j. If matrix, columns correspond to groups.
#' @param tr Trim proportion for trimmed mean methods (default: 0.2).
#' @param method Effect size method to use (default: 'EP'). Options:
#'   \itemize{
#'     \item \code{'EP'}: Explanatory measure of effect size
#'     \item \code{'QS'}: Median-type quantile shift measure
#'     \item \code{'QStr'}: Trimmed mean quantile shift measure
#'     \item \code{'AKP'}: Trimmed mean Winsorized variance analog of Cohen's d
#'     \item \code{'WMW'}: P(X < Y) probabilistic superiority measure
#'     \item \code{'KMS'}: Kulinskaya et al. method
#'   }
#' @param pr Logical. Print messages during computation (default: FALSE).
#' @param SEED Logical. Set random seed for reproducibility (default: TRUE).
#'
#' @return A list with component:
#'   \item{Estimates}{Matrix with columns: Group, Group, Effect Size.
#'     Each row shows the effect size comparing one pair of groups.}
#'
#' @details
#' This function computes robust heteroscedastic effect sizes for all pairwise
#' comparisons among J independent groups. Unlike classical Cohen's d, these
#' methods account for unequal variances and are robust to outliers.
#'
#' The 'EP' (explanatory power) method measures the proportion of variance
#' explained. The 'QS' methods measure distributional shift (0-1 scale).
#' The 'AKP' method is a robust analog of Cohen's d using trimmed means
#' and Winsorized variances.
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @export
#' @family MCP functions
#' @family effect size functions
#' @seealso \code{\link{linconEP}}, \code{\link{linconES}}, \code{\link{akp.effect}}
#' @examples
#' \dontrun{
#' # Compare three groups with explanatory power effect sizes
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, mean=0.5), rnorm(20, mean=1))
#' result <- esmcp(x, method='EP')
#' print(result$Estimates)
#'
#' # Use AKP robust Cohen's d analog
#' result2 <- esmcp(x, method='AKP')
#' print(result2$Estimates)
#' }
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


#' Two-Way Between-Subjects MCP with Explanatory Power Effect Sizes
#'
#' Performs multiple comparisons for a two-way independent groups (between-between)
#' design with explanatory power effect sizes. Tests Factor A and B main effects
#' and A×B interaction using trimmed means.
#'
#' @param J Number of levels for Factor A (between-subjects).
#' @param K Number of levels for Factor B (between-subjects).
#' @param x Data in list or matrix format. If list, \code{x[[1]]} contains data
#'   for cell (1,1), \code{x[[2]]} for cell (1,2), etc. Length should be J×K.
#' @param tr Trim proportion (default: 0.2, 20% trimming).
#' @param alpha Significance level (default: 0.05).
#' @param grp Optional vector to reorder groups if data not in expected order (default: NA).
#' @param op Logical. If TRUE, return all tests combined; if FALSE (default),
#'   return separate results for each factor and interaction.
#' @param nreps Number of replications for effect size estimation (default: 200).
#' @param SEED Logical. Set random seed for reproducibility (default: TRUE).
#' @param pr Logical. Print informational messages (default: TRUE).
#' @param POOL Logical. Use pooled variance for main effects (default: TRUE).
#'
#' @return A list with components:
#'   \item{Factor.A}{Results from \code{\link{linconEP}} for Factor A main effect.}
#'   \item{Factor.B}{Results from \code{\link{linconEP}} for Factor B main effect.}
#'   \item{Factor.AB}{Results from \code{\link{linconEP}} for A×B interaction.}
#'   \item{All.Tests}{If \code{op=TRUE}, combined results for all tests (default: NA).}
#'   \item{conA, conB, conAB}{Contrast matrices used.}
#'
#' @details
#' This function extends \code{bbmcp} by adding explanatory power (EP) effect sizes,
#' which measure the proportion of variance explained by group differences.
#' Contrast matrices are generated using \code{\link{con2way}}.
#'
#' Explanatory power ranges from 0 to 1, with larger values indicating stronger
#' effects. Unlike R², EP uses robust trimmed means and accounts for heteroscedasticity.
#'
#' @export
#' @family MCP functions
#' @seealso \code{\link{linconEP}}, \code{\link{bbmcpQS}}, \code{\link{mcp2a}}, \code{\link{esmcp}}
#' @examples
#' \dontrun{
#' # 2×3 between-subjects design
#' set.seed(123)
#' J <- 2
#' K <- 3
#' x <- list()
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     x[[(j-1)*K + k]] <- rnorm(20, mean = j + k)
#'   }
#' }
#' result <- bbmcpEP(J, K, x)
#' print(result$Factor.A)
#' }
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


#' Two-Way Between-Subjects MCP with Quantile Shift Effect Sizes
#'
#' Performs multiple comparisons for a two-way independent groups (between-between)
#' design with quantile shift effect sizes. Computes effect sizes for Factor A
#' and B main effects and A×B interaction.
#'
#' @param J Number of levels for Factor A (between-subjects).
#' @param K Number of levels for Factor B (between-subjects).
#' @param x Data in list or matrix format. If list, \code{x[[1]]} contains data
#'   for cell (1,1), \code{x[[2]]} for cell (1,2), etc. Length should be J×K.
#' @param locfun Location function for quantile shift computation (e.g., \code{median}, \code{tmean}).
#' @param nreps Number of replications for effect size estimation (default: 100).
#' @param SEED Logical. Set random seed for reproducibility (default: TRUE).
#' @param POOL Logical. Pool data across levels when computing main effects (default: TRUE).
#' @param pr Logical. Print informational messages (default: TRUE).
#'
#' @return A list with components:
#'   \item{Factor.A}{Quantile shift results for Factor A main effect.
#'     If \code{POOL=FALSE}, a list of results for each level of A.}
#'   \item{Factor.B}{Quantile shift results for Factor B main effect.
#'     If \code{POOL=FALSE}, a list of results for each level of B.}
#'   \item{Factor.AB}{List of quantile shift results for each A×B interaction contrast.}
#'
#' @details
#' This function computes quantile shift effect sizes for all pairwise comparisons
#' in a two-way factorial design. Quantile shift measures distributional differences,
#' ranging from 0 to 1 with 0.5 indicating no effect.
#'
#' Under normality and equal variances, Cohen's d values of 0, 0.2, 0.5, 0.8
#' correspond approximately to quantile shift values of 0.50, 0.55, 0.65, 0.70.
#'
#' When \code{POOL=TRUE}, data are pooled across factor levels before computing
#' main effects. When \code{POOL=FALSE}, separate comparisons are made at each
#' level of the other factor.
#'
#' @export
#' @family MCP functions
#' @family effect size functions
#' @seealso \code{\link{bbmcpEP}}, \code{\link{linconQS}}, \code{\link{esmcp}}
#' @examples
#' \dontrun{
#' # 2×3 between-subjects design with median-based quantile shift
#' set.seed(123)
#' J <- 2
#' K <- 3
#' x <- list()
#' for(j in 1:J) {
#'   for(k in 1:K) {
#'     x[[(j-1)*K + k]] <- rnorm(20, mean = j + k)
#'   }
#' }
#' result <- bbmcpQS(J, K, x, locfun=median)
#' print(result$Factor.A)
#' }
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


#' Between-by-Between Factorial MCP for Two-Way Designs
#'
#' Performs multiple comparisons for a two-way between-subjects (independent groups)
#' factorial design. For each level of Factor A, compares all pairs of Factor B levels,
#' and vice versa.
#'
#' @param J Number of levels for Factor A
#' @param K Number of levels for Factor B
#' @param x Data in list mode (length J*K), matrix, or data frame with J*K columns.
#'   Groups ordered as: A1B1, A1B2, ..., A1BK, A2B1, A2B2, ..., AJBK
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming)
#'
#' @return A list with components:
#'   \item{Levels.of.A}{List of length J. Each element contains MCP results
#'     (from \code{\link{lincon}}) comparing all pairs of Factor B levels
#'     at that level of Factor A}
#'   \item{Level.of.B}{List of length K. Each element contains MCP results
#'     comparing all pairs of Factor A levels at that level of Factor B}
#'
#' @details
#' This function performs a "simple effects" analysis for a two-way factorial design
#' with independent groups. It uses trimmed means and the linear contrast function
#' \code{\link{lincon}} to control family-wise error rate.
#'
#' \strong{Simple Effects Strategy:}
#' \itemize{
#'   \item For Factor A: At each level of A (j=1,...,J), compare all pairs
#'     of B levels using trimmed means
#'   \item For Factor B: At each level of B (k=1,...,K), compare all pairs
#'     of A levels using trimmed means
#' }
#'
#' The function expects data organized in a specific order. If using a matrix or
#' data frame, columns should represent groups in row-major order of the J×K
#' factorial structure.
#'
#' @note
#' This function uses the "det" (deterministic/trimmed mean) approach. For a
#' quantile-based version, see \code{\link{bbdetmcpQS}}.
#'
#' @examples
#' \dontrun{
#' # 2x3 factorial design (2 levels of A, 3 levels of B)
#' set.seed(123)
#' # Generate data for 6 groups (A1B1, A1B2, A1B3, A2B1, A2B2, A2B3)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1),
#'           rnorm(20, 0.3), rnorm(20, 0.8), rnorm(20, 1.3))
#'
#' result <- bbdetmcp(J=2, K=3, x=x, tr=0.2)
#'
#' # Results for Factor A (compare B levels within each A level)
#' result$Levels.of.A[[1]]  # Compare B1, B2, B3 at level A1
#' result$Levels.of.A[[2]]  # Compare B1, B2, B3 at level A2
#'
#' # Results for Factor B (compare A levels within each B level)
#' result$Level.of.B[[1]]   # Compare A1, A2 at level B1
#' }
#'
#' @family multiple comparison procedures
#' @family factorial design methods
#' @seealso \code{\link{bbdetmcpQS}} for quantile-based version,
#'   \code{\link{lincon}} for the underlying MCP method,
#'   \code{\link{mcp2a}} for main effects and interactions
#' @export
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


#' Between-by-Between Factorial MCP Using Quantile Shift
#'
#' Performs multiple comparisons for a two-way between-subjects factorial design
#' using quantile shift measures instead of trimmed means. For each level of Factor A,
#' compares all pairs of Factor B levels, and vice versa.
#'
#' @param J Number of levels for Factor A
#' @param K Number of levels for Factor B
#' @param x Data in list mode (length J*K), matrix, or data frame with J*K columns.
#'   Groups ordered as: A1B1, A1B2, ..., A1BK, A2B1, A2B2, ..., AJBK
#' @param tr Proportion of trimming used in the quantile shift calculation (default: 0.2)
#'
#' @return A list with components:
#'   \item{Levels.of.A}{List of length J. Each element contains MCP results
#'     (from \code{\link{linconQS}}) comparing all pairs of Factor B levels
#'     at that level of Factor A using quantile shift}
#'   \item{Level.of.B}{List of length K. Each element contains MCP results
#'     comparing all pairs of Factor A levels at that level of Factor B}
#'
#' @details
#' This function is the quantile shift (QS) version of \code{\link{bbdetmcp}}.
#' It performs simple effects analysis using \code{\link{linconQS}}, which is
#' based on quantile shifts rather than differences in trimmed means.
#'
#' Quantile shift measures provide information about distributional differences
#' beyond just location shifts, making them useful when distributions may differ
#' in shape or spread.
#'
#' @examples
#' \dontrun{
#' # 2x3 factorial design
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1),
#'           rnorm(20, 0.3), rnorm(20, 0.8), rnorm(20, 1.3))
#'
#' result <- bbdetmcpQS(J=2, K=3, x=x, tr=0.2)
#'
#' # Compare with trimmed mean version
#' result_tm <- bbdetmcp(J=2, K=3, x=x, tr=0.2)
#' }
#'
#' @family multiple comparison procedures
#' @family factorial design methods
#' @family quantile-based methods
#' @seealso \code{\link{bbdetmcp}} for trimmed mean version,
#'   \code{\link{linconQS}} for the underlying quantile shift MCP method
#' @export
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


#' Effect Sizes for All Pairwise Comparisons (Independent Groups, AKP Method)
#'
#' Computes the Algina-Keselman-Penfield (AKP) robust effect size measure for all
#' pairwise comparisons among independent groups. The AKP effect size is based on
#' trimmed means and Winsorized variances.
#'
#' @param x Data in list mode (one element per group), or a matrix/data frame
#'   with groups as columns
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming)
#'
#' @return A matrix with J(J-1)/2 rows (one per pairwise comparison) and 3 columns:
#'   \item{Group}{Index of first group in the comparison}
#'   \item{Group}{Index of second group in the comparison}
#'   \item{Effect.Size}{AKP effect size for the comparison}
#'
#' @details
#' The AKP effect size is a robust analog of Cohen's d that uses trimmed means
#' for location and a pooled Winsorized variance for scale. It is computed using
#' the \code{\link{akp.effect}} function for each pair of groups.
#'
#' This measure is appropriate for independent groups and provides a standardized
#' effect size that is resistant to outliers and non-normality.
#'
#' The effect size formula is:
#' \deqn{d_{AKP} = \frac{\bar{X}_{t1} - \bar{X}_{t2}}{s_w}}
#' where \eqn{\bar{X}_{t}} are trimmed means and \eqn{s_w} is the pooled
#' Winsorized standard deviation.
#'
#' @note
#' For dependent (paired) groups, use \code{\link{wmcpAKP}} instead.
#' For quantile shift effect sizes, see \code{\link{bmcpQS}}.
#'
#' @examples
#' \dontrun{
#' # Compare 4 independent groups
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1), rnorm(20, 1.5))
#'
#' # Compute AKP effect sizes for all pairs
#' effects <- bmcpAKP(x, tr=0.2)
#' print(effects)
#'
#' # Results show Group 1 vs 2, 1 vs 3, 1 vs 4, 2 vs 3, 2 vs 4, 3 vs 4
#' }
#'
#' @family multiple comparison procedures
#' @family effect size measures
#' @seealso \code{\link{akp.effect}} for single comparison,
#'   \code{\link{wmcpAKP}} for dependent groups,
#'   \code{\link{bmcpQS}} for quantile shift effect sizes,
#'   \code{\link{lincon}} for hypothesis tests with trimmed means
#' @export
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


#' Effect Sizes for All Pairwise Comparisons (Independent Groups, Quantile Shift)
#'
#' Computes quantile shift effect size measures for all pairwise comparisons among
#' independent groups. Quantile shifts provide information about distributional
#' differences beyond location.
#'
#' @param x Data in list mode (one element per group), or a matrix/data frame
#'   with groups as columns
#' @param locfun Location estimator function (default: \code{median})
#' @param ... Additional arguments passed to \code{\link{shiftes}}
#'
#' @return A matrix with J(J-1)/2 rows (one per pairwise comparison) and 3 columns:
#'   \item{Group}{Index of first group in the comparison}
#'   \item{Group}{Index of second group in the comparison}
#'   \item{Effect.Size}{Quantile shift effect size (Q.Effect) for the comparison}
#'
#' @details
#' This function computes the quantile shift effect size using \code{\link{shiftes}}
#' for each pair of independent groups. The quantile shift measures how much one
#' distribution is shifted relative to another across their entire range.
#'
#' Unlike traditional effect sizes that focus only on central tendency, quantile
#' shift provides a more complete picture of distributional differences, detecting
#' shifts in any part of the distribution.
#'
#' The default uses the median as the location estimator, but any robust estimator
#' can be specified via \code{locfun}.
#'
#' @note
#' For dependent groups, use \code{\link{wmcpQS}} instead.
#' For AKP effect sizes based on trimmed means, see \code{\link{bmcpAKP}}.
#'
#' @examples
#' \dontrun{
#' # Compare 4 independent groups
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1), rnorm(20, 1.5))
#'
#' # Compute quantile shift effect sizes for all pairs
#' effects <- bmcpQS(x)
#' print(effects)
#'
#' # Use trimmed mean instead of median
#' effects_tm <- bmcpQS(x, locfun=mean, tr=0.2)
#' }
#'
#' @family multiple comparison procedures
#' @family effect size measures
#' @family quantile-based methods
#' @seealso \code{\link{shiftes}} for single comparison,
#'   \code{\link{wmcpQS}} for dependent groups,
#'   \code{\link{bmcpAKP}} for AKP effect sizes,
#'   \code{\link{linconQS}} for hypothesis tests using quantile shift
#' @export
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


#' Effect Sizes for All Pairwise Comparisons (Dependent Groups, AKP Method)
#'
#' Computes the Algina-Keselman-Penfield (AKP) robust effect size measure for all
#' pairwise comparisons among dependent (repeated measures) groups. Uses the dependent
#' groups version of the AKP effect size.
#'
#' @param x Data in matrix or list mode. If matrix, rows are subjects and columns
#'   are repeated measures. If list, each element contains one group's data
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming)
#'
#' @return A matrix with J(J-1)/2 rows (one per pairwise comparison) and 3 columns:
#'   \item{Group}{Index of first group in the comparison}
#'   \item{Group}{Index of second group in the comparison}
#'   \item{Effect.Size}{Dependent groups AKP effect size for the comparison}
#'
#' @details
#' This function computes a robust effect size for dependent (paired) data using
#' the \code{\link{D.akp.effect}} function for each pair of groups. The dependent
#' groups AKP measure accounts for the correlation between measurements.
#'
#' The effect size is based on the trimmed mean of the difference scores divided
#' by an appropriate Winsorized standard deviation that accounts for pairing.
#'
#' This is the within-subjects analog of \code{\link{bmcpAKP}} and is appropriate
#' for repeated measures designs where the same subjects are measured under
#' multiple conditions.
#'
#' @note
#' For independent groups, use \code{\link{bmcpAKP}} instead.
#' For quantile shift effect sizes with dependent groups, see \code{\link{wmcpQS}}.
#'
#' @examples
#' \dontrun{
#' # Compare 4 dependent (repeated measures) groups
#' set.seed(123)
#' # Generate correlated data (subjects x conditions)
#' n <- 20
#' x <- matrix(NA, n, 4)
#' for(i in 1:n) {
#'   base <- rnorm(1)
#'   x[i,] <- base + rnorm(4, mean=c(0, 0.3, 0.6, 0.9), sd=0.5)
#' }
#'
#' # Compute dependent groups AKP effect sizes for all pairs
#' effects <- wmcpAKP(x, tr=0.2)
#' print(effects)
#' }
#'
#' @family multiple comparison procedures
#' @family effect size measures
#' @family dependent groups methods
#' @seealso \code{\link{D.akp.effect}} for single comparison,
#'   \code{\link{bmcpAKP}} for independent groups,
#'   \code{\link{wmcpQS}} for quantile shift effect sizes,
#'   \code{\link{lindep}} for hypothesis tests with dependent groups
#' @export
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


#' Effect Sizes for All Pairwise Comparisons (Dependent Groups, Quantile Shift)
#'
#' Computes quantile shift effect size measures for all pairwise comparisons among
#' dependent (repeated measures) groups. Quantile shifts provide information about
#' distributional differences beyond location.
#'
#' @param x Data in matrix or list mode. If matrix, rows are subjects and columns
#'   are repeated measures. If list, each element contains one group's data
#' @param locfun Location estimator function (default: \code{median})
#'
#' @return A matrix with J(J-1)/2 rows (one per pairwise comparison) and 3 columns:
#'   \item{Group}{Index of first group in the comparison}
#'   \item{Group}{Index of second group in the comparison}
#'   \item{Q.Effect}{Quantile shift effect size for the comparison}
#'
#' @details
#' This function computes the quantile shift effect size using \code{\link{depQS}}
#' for each pair of dependent groups. The quantile shift for dependent data is
#' computed on the difference scores and measures distributional shifts beyond
#' simple location differences.
#'
#' Unlike traditional effect sizes that focus only on central tendency, quantile
#' shift provides information about how the entire distribution of differences
#' is shifted, which can be valuable for detecting effects in different parts of
#' the distribution.
#'
#' This is the within-subjects analog of \code{\link{bmcpQS}} and is appropriate
#' for repeated measures designs.
#'
#' @note
#' For independent groups, use \code{\link{bmcpQS}} instead.
#' For AKP effect sizes with dependent groups, see \code{\link{wmcpAKP}}.
#'
#' @examples
#' \dontrun{
#' # Compare 4 dependent (repeated measures) groups
#' set.seed(123)
#' # Generate correlated data (subjects x conditions)
#' n <- 20
#' x <- matrix(NA, n, 4)
#' for(i in 1:n) {
#'   base <- rnorm(1)
#'   x[i,] <- base + rnorm(4, mean=c(0, 0.3, 0.6, 0.9), sd=0.5)
#' }
#'
#' # Compute quantile shift effect sizes for all pairs
#' effects <- wmcpQS(x)
#' print(effects)
#' }
#'
#' @family multiple comparison procedures
#' @family effect size measures
#' @family dependent groups methods
#' @family quantile-based methods
#' @seealso \code{\link{depQS}} for single comparison,
#'   \code{\link{bmcpQS}} for independent groups,
#'   \code{\link{wmcpAKP}} for AKP effect sizes,
#'   \code{\link{lindep}} for hypothesis tests with dependent groups
#' @export
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

#' Step-Down Multiple Comparison Procedure for Trimmed Means
#'
#' Performs a step-down multiple comparison procedure for testing hypotheses about
#' trimmed means across J independent groups. Uses a modified Bonferroni approach
#' with step-down critical values.
#'
#' @param x Data in list mode (one element per group), or a matrix/data frame
#'   with groups as columns
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming)
#' @param alpha Family-wise error rate (default: 0.05)
#'
#' @return A matrix with columns:
#'   \item{Groups}{Numeric code identifying which groups are being compared}
#'   \item{p-value}{P-value from \code{\link{t1way}} for this subset of groups}
#'   \item{p.crit}{Critical p-value for this comparison using step-down method}
#'
#' @details
#' The step-down method tests nested hypotheses, starting with all J groups and
#' progressively testing smaller subsets. For each subset, it uses the heteroscedastic
#' one-way ANOVA for trimmed means (\code{\link{t1way}}) and compares the p-value
#' to an adjusted critical value.
#'
#' The critical values are determined using a step-down modification of Bonferroni's
#' method that accounts for the number of groups in each comparison. This can provide
#' more power than standard multiple comparison procedures.
#'
#' \strong{Limitations:}
#' \itemize{
#'   \item Currently limited to at most 5 groups
#'   \item Requires at least 3 groups (use \code{\link{yuen}} for 2 groups)
#' }
#'
#' The step-down procedure controls the family-wise error rate at level alpha
#' across all tested hypotheses.
#'
#' @note
#' For more than 5 groups, use \code{\link{lincon}} or other MCP methods.
#'
#' @examples
#' \dontrun{
#' # Compare 4 groups using step-down method
#' set.seed(123)
#' x <- list(rnorm(20), rnorm(20, 0.5), rnorm(20, 1), rnorm(20, 1.5))
#'
#' result <- stepmcp(x, tr=0.2, alpha=0.05)
#' print(result)
#'
#' # Identify significant differences
#' sig <- result[result[,"p-value"] <= result[,"p.crit"], ]
#' }
#'
#' @family multiple comparison procedures
#' @family trimmed mean methods
#' @seealso \code{\link{t1way}} for the omnibus test,
#'   \code{\link{lincon}} for general linear contrasts,
#'   \code{\link{yuen}} for two-group comparison
#' @export
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


#' Sign Test Multiple Comparisons for Dependent Groups
#'
#' Performs sign tests for all pairwise comparisons among J dependent (repeated
#' measures) groups. Uses the sign test to compare medians of difference scores,
#' with p-value adjustment to control family-wise error.
#'
#' @param x Data matrix where rows are subjects and columns are repeated measures,
#'   or data in list mode (one element per group)
#' @param y Not used; included for compatibility (default: NULL)
#' @param alpha Nominal Type I error rate (default: 0.05)
#' @param method Method for computing confidence interval in sign test:
#'   'AC' (Agresti-Coull), 'SD' (score), or other methods supported by
#'   \code{\link{signt}} (default: 'AC')
#' @param AUTO Logical; if TRUE, automatically selects method based on sample size
#'   (default: TRUE)
#' @param Method P-value adjustment method for multiple comparisons. Any method
#'   accepted by \code{\link{p.adjust}}: "hochberg" (default), "holm", "bonferroni",
#'   "BH", "BY", "fdr", etc.
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: Group, Group (indices of compared groups),
#'     n (sample size after removing ties), N (original sample size),
#'     Prob_x_less_than_y (estimated probability), ci.lower, ci.upper,
#'     p.value (unadjusted), p.adjusted (adjusted using specified Method)}
#'
#' @details
#' This function performs all J(J-1)/2 pairwise sign tests for J dependent groups.
#' The sign test is a nonparametric test that makes minimal distributional assumptions
#' and is based on the median of the difference scores.
#'
#' For each pair of groups, it computes:
#' \itemize{
#'   \item The sign test comparing the two groups using \code{\link{signt}}
#'   \item A confidence interval for the probability P(X < Y)
#'   \item Unadjusted and adjusted p-values
#' }
#'
#' P-values are adjusted using the method specified in \code{Method} to control
#' the family-wise error rate. The default Hochberg method is slightly more
#' powerful than the Bonferroni correction.
#'
#' @note
#' The sign test is less powerful than tests based on trimmed means when the
#' distribution is approximately symmetric, but more robust to extreme outliers.
#'
#' @examples
#' \dontrun{
#' # Compare 4 dependent groups using sign test
#' set.seed(123)
#' # Generate correlated data (subjects x conditions)
#' n <- 25
#' x <- matrix(NA, n, 4)
#' for(i in 1:n) {
#'   base <- rnorm(1)
#'   x[i,] <- base + rnorm(4, mean=c(0, 0.3, 0.6, 0.9), sd=0.5)
#' }
#'
#' result <- signmcp(x, alpha=0.05)
#' print(result$output)
#'
#' # Identify significant comparisons (using adjusted p-values)
#' sig <- result$output[result$output[,"p.adjusted"] < 0.05, ]
#' }
#'
#' @family multiple comparison procedures
#' @family dependent groups methods
#' @family nonparametric methods
#' @seealso \code{\link{signt}} for single sign test comparison,
#'   \code{\link{lindep}} for trimmed means with dependent groups,
#'   \code{\link{p.adjust}} for adjustment methods
#' @export
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


#' Multiple Comparisons for Independent Groups with Discrete Distributions
#'
#' Performs all pairwise comparisons for J independent groups when data have
#' discrete distributions. Uses a chi-squared based test for each comparison
#' with Hochberg's method to control family-wise error rate.
#'
#' @param x Data in list mode (one element per group) or matrix (groups as columns).
#'   Missing values are allowed
#' @param alpha Family-wise error rate (default: 0.05)
#' @param nboot Number of bootstrap samples for determining critical values
#'   (default: 500)
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#' @param ... Additional arguments (currently unused)
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: Group, Group (indices of compared groups),
#'     p.value (from chi-squared test), p.crit (critical p-value using Hochberg)}
#'   \item{num.sig}{Number of significant pairwise comparisons}
#'
#' @details
#' This function is designed for comparing distributions when the data are discrete
#' (e.g., count data, ordinal data, Likert scales). It uses \code{\link{disc2com}}
#' to perform a chi-squared test for each pair of groups.
#'
#' The chi-squared test detects any distributional differences between groups,
#' not just location shifts. This makes it appropriate for discrete data where
#' traditional mean or median comparisons may not be suitable.
#'
#' \strong{Multiple Comparison Adjustment:}
#' Hochberg's step-up method is used to control the family-wise error rate. This
#' method is uniformly more powerful than the Bonferroni correction while still
#' controlling Type I error.
#'
#' Critical p-values are determined by ordering the observed p-values and comparing
#' them to adjusted alpha levels: alpha/1, alpha/2, ..., alpha/C where C is the
#' number of comparisons.
#'
#' @note
#' For continuous data, use \code{\link{lincon}} or related functions instead.
#' The bootstrap is used within \code{disc2com} to estimate the p-value.
#'
#' @examples
#' \dontrun{
#' # Compare 3 groups with discrete data (e.g., Likert scale 1-5)
#' set.seed(123)
#' x <- list(
#'   sample(1:5, 30, replace=TRUE, prob=c(0.1, 0.2, 0.4, 0.2, 0.1)),
#'   sample(1:5, 30, replace=TRUE, prob=c(0.2, 0.3, 0.3, 0.1, 0.1)),
#'   sample(1:5, 30, replace=TRUE, prob=c(0.1, 0.1, 0.2, 0.3, 0.3))
#' )
#'
#' result <- discmcp(x, alpha=0.05, nboot=1000)
#' print(result$output)
#' print(paste("Number of significant differences:", result$num.sig))
#' }
#'
#' @family multiple comparison procedures
#' @family nonparametric methods
#' @family discrete data methods
#' @seealso \code{\link{disc2com}} for two-group comparison,
#'   \code{\link{lincon}} for continuous data,
#'   \code{\link{skmcp}} for binary data
#' @export
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


#' Inferences on Medians for Dependent Groups Using Multiple Comparisons
#'
#' Performs multiple comparisons for dependent (repeated measures) groups using
#' medians. Tests linear contrasts or all pairwise comparisons of medians for
#' difference scores with Hochberg's method to control family-wise error.
#'
#' @param x Data matrix where rows are subjects and columns are repeated measures,
#'   or data in list mode (one element per group)
#' @param con Contrast matrix (J × d) where J is number of groups and d is number
#'   of contrasts. If not specified (con=0), all pairwise comparisons are performed
#' @param alpha Family-wise error rate (default: 0.05)
#'
#' @return A list with components:
#'   \item{output}{Matrix with results. For pairwise comparisons: Group, Group,
#'     psihat (median of differences), ci.lower, ci.upper, p.value, p.crit.
#'     For custom contrasts: con.num, psihat, ci.lower, ci.upper, p.value, p.crit}
#'   \item{con}{The contrast matrix used}
#'   \item{num.sig}{Number of significant comparisons}
#'
#' @details
#' This function performs inference on medians for dependent groups. For each
#' pairwise comparison or linear contrast, it computes the median of the
#' (weighted) difference scores and a confidence interval using \code{\link{sintv2}}.
#'
#' \strong{For Pairwise Comparisons} (con=0):
#' Computes all J(J-1)/2 pairwise differences and tests whether each median
#' difference is zero.
#'
#' \strong{For Custom Contrasts:}
#' For each contrast column, computes weighted difference scores using the
#' contrast coefficients and tests whether the median is zero.
#'
#' \strong{Multiple Comparison Adjustment:}
#' Uses Hochberg's step-up method with special critical values optimized for
#' alpha = 0.05 or 0.01. For other alpha levels, uses standard Bonferroni adjustment.
#'
#' Missing values are automatically removed before analysis.
#'
#' @note
#' For trimmed means instead of medians, use \code{\link{lindep}}.
#' For independent groups, use \code{\link{discmcp}} or other methods.
#'
#' @examples
#' \dontrun{
#' # Compare 4 dependent groups using medians
#' set.seed(123)
#' # Generate correlated data (subjects x conditions)
#' n <- 25
#' x <- matrix(NA, n, 4)
#' for(i in 1:n) {
#'   base <- rnorm(1)
#'   x[i,] <- base + rnorm(4, mean=c(0, 0.3, 0.6, 0.9), sd=0.5)
#' }
#'
#' # All pairwise comparisons
#' result <- sintmcp(x, alpha=0.05)
#' print(result$output)
#'
#' # Custom contrasts: (1+2)/2 vs (3+4)/2 and 1 vs 2
#' con <- matrix(c(1, 1, -1, -1,
#'                 1, -1, 0, 0) / 2, nrow=4, ncol=2)
#' result2 <- sintmcp(x, con=con, alpha=0.05)
#' }
#'
#' @family multiple comparison procedures
#' @family dependent groups methods
#' @family median-based methods
#' @seealso \code{\link{sintv2}} for single median test,
#'   \code{\link{lindep}} for trimmed means with dependent groups,
#'   \code{\link{signmcp}} for sign test version
#' @export
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


#' Multiple Comparisons for Independent Binomial Proportions (Storer-Kim Method)
#'
#' Performs all pairwise comparisons for J independent groups with binary data
#' using the Storer-Kim method for comparing binomial proportions. Uses Hochberg's
#' method to control family-wise error rate.
#'
#' @param x Data in list mode (one element per group) or matrix (groups as columns).
#'   Data should be binary (0/1). Missing values are allowed
#' @param alpha Family-wise error rate (default: 0.05)
#'
#' @return A list with components:
#'   \item{output}{Matrix with columns: Group, Group (indices of compared groups),
#'     p.value (from Storer-Kim test), p.crit (critical p-value using Hochberg)}
#'   \item{num.sig}{Number of significant pairwise comparisons}
#'
#' @details
#' This function is designed for comparing proportions when the data are binary
#' (e.g., success/failure, yes/no, 0/1). It uses the Storer-Kim method via
#' \code{\link{twobinom}} to test each pair of groups.
#'
#' The Storer-Kim method provides accurate inference for binomial proportions
#' without relying on large-sample normal approximations. It is particularly
#' useful when sample sizes are small or when proportions are near 0 or 1.
#'
#' \strong{Multiple Comparison Adjustment:}
#' Hochberg's step-up method is used to control the family-wise error rate. This
#' method is uniformly more powerful than the Bonferroni correction while still
#' controlling Type I error at the nominal alpha level.
#'
#' Critical p-values are determined by ordering the observed p-values and comparing
#' them to adjusted alpha levels: alpha/C, alpha/(C-1), ..., alpha/1 where C is
#' the number of comparisons.
#'
#' @note
#' For continuous data, use \code{\link{lincon}} or related functions.
#' For general discrete data (not just binary), use \code{\link{discmcp}}.
#'
#' @examples
#' \dontrun{
#' # Compare 3 groups with binary data
#' set.seed(123)
#' x <- list(
#'   rbinom(30, 1, 0.3),  # Group 1: 30% success rate
#'   rbinom(30, 1, 0.5),  # Group 2: 50% success rate
#'   rbinom(30, 1, 0.7)   # Group 3: 70% success rate
#' )
#'
#' result <- skmcp(x, alpha=0.05)
#' print(result$output)
#' print(paste("Number of significant differences:", result$num.sig))
#' }
#'
#' @family multiple comparison procedures
#' @family nonparametric methods
#' @family binomial methods
#' @seealso \code{\link{twobinom}} for two-group Storer-Kim test,
#'   \code{\link{binmcp}} for alternative binomial MCP method,
#'   \code{\link{discmcp}} for general discrete data
#' @export
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


#' Multiple Comparisons for Nested Designs with Trimmed Means
#'
#' Performs multiple comparisons for nested (hierarchical) ANOVA designs using
#' trimmed means. Tests all pairwise comparisons among the levels of the main
#' factor (Factor A) by pooling nested observations.
#'
#' @param x Data in list mode with length J (number of levels of Factor A).
#'   Each element x[[j]] is a matrix (n × K) containing the nested data for
#'   level j of Factor A, where K is the number of nested observations per subject
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming)
#'
#' @return Results from \code{\link{lincon}}: a list containing pairwise comparisons
#'   of the J levels of Factor A using trimmed means, with adjusted confidence
#'   intervals and p-values
#'
#' @details
#' This function handles nested (hierarchical) designs where observations are
#' nested within levels of a main factor. The strategy is to pool all nested
#' observations for each level of Factor A into a single vector, then perform
#' all pairwise comparisons using \code{\link{lincon}}.
#'
#' \strong{Data Structure:}
#' \itemize{
#'   \item J levels of Factor A
#'   \item For each level j, there is an n × K matrix
#'   \item n = number of subjects at this level
#'   \item K = number of nested observations per subject
#' }
#'
#' For example, in a study with 3 schools (Factor A), each school has multiple
#' classrooms, and each classroom has multiple students. The function pools all
#' students within each school and compares schools.
#'
#' The function uses trimmed means to provide robust inference that is resistant
#' to outliers and violations of normality.
#'
#' @examples
#' \dontrun{
#' # Nested design: 3 schools, each with 4 classrooms of 10 students
#' set.seed(123)
#' x <- list()
#' for(j in 1:3) {
#'   # Each school has different mean performance
#'   school_effect <- (j-2) * 0.5
#'   x[[j]] <- matrix(rnorm(40, mean=school_effect), nrow=10, ncol=4)
#' }
#'
#' # Compare schools using trimmed means
#' result <- mcp.nestAP(x, tr=0.2)
#' print(result)
#' }
#'
#' @family multiple comparison procedures
#' @family nested and hierarchical designs
#' @family trimmed mean methods
#' @seealso \code{\link{lincon}} for the underlying MCP method,
#'   \code{\link{t1way}} for omnibus test
#' @export
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


#' Multiple Comparisons for Binomial Proportions
#'
#' Performs all pairwise comparisons of binomial proportions for J independent
#' groups. Uses simulation to determine critical p-values that control family-wise
#' error rate, and provides simultaneous confidence intervals for differences.
#'
#' @param x Vector of length J containing the number of successes in each group
#' @param n Vector of length J containing the sample sizes for each group
#' @param p.crit Optional vector of critical p-values. If NULL (default), critical
#'   values are determined via simulation to control FWE at level alpha
#' @param alpha Family-wise error rate (default: 0.05)
#' @param iter Number of simulation iterations for determining critical p-values
#'   (default: 2000). Used only if p.crit is NULL
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#'
#' @return A list with components:
#'   \item{n}{Sample sizes for each group}
#'   \item{output}{Matrix with columns: Grp, Grp (indices of compared groups),
#'     Est 1, Est 2 (estimated proportions), Dif (difference in proportions),
#'     ci.low, ci.up (simultaneous confidence interval), p.value, p.crit
#'     (critical p-value for this comparison)}
#'
#' @details
#' This function compares binomial proportions across J independent groups using
#' the Kulinskaya-Morgenthaler-Staudte (KMS) method via \code{\link{bi2KMSv2}}
#' for each pairwise comparison.
#'
#' \strong{Critical P-value Determination:}
#' When \code{p.crit=NULL}, the function uses simulation to determine critical
#' p-values that control the family-wise error rate. Under the null hypothesis
#' of equal proportions, it:
#' \itemize{
#'   \item Estimates the common proportion as the pooled success rate
#'   \item Generates \code{iter} datasets under this null
#'   \item Computes p-values for all pairwise comparisons
#'   \item Determines the quantile that controls FWE at level alpha
#' }
#'
#' \strong{Simultaneous Confidence Intervals:}
#' Confidence intervals are adjusted to have simultaneous coverage probability
#' of 1-alpha across all J(J-1)/2 pairwise comparisons.
#'
#' @note
#' The KMS method provides better performance than traditional normal approximations,
#' especially with small sample sizes or extreme proportions.
#'
#' @examples
#' \dontrun{
#' # Compare success rates across 4 groups
#' x <- c(15, 22, 18, 25)  # Number of successes
#' n <- c(30, 30, 30, 30)  # Sample sizes
#'
#' result <- binmcp(x, n, alpha=0.05, iter=5000)
#' print(result$output)
#'
#' # Identify significant differences
#' sig <- result$output[result$output[,"p.value"] < result$output[,"p.crit"], ]
#' }
#'
#' @family multiple comparison procedures
#' @family binomial methods
#' @seealso \code{\link{bi2KMSv2}} for two-group comparison,
#'   \code{\link{binmcp.crit}} for critical value simulation,
#'   \code{\link{skmcp}} for Storer-Kim method
#' @export
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


#' Simulate Critical P-values for Binomial Multiple Comparisons
#'
#' Simulates the null distribution of p-values for multiple comparisons of binomial
#' proportions to determine critical values that control family-wise error rate.
#' This is a helper function used by \code{\link{binmcp}}.
#'
#' @param p Common proportion under the null hypothesis (pooled proportion)
#' @param n Vector of sample sizes for each group
#' @param iter Number of simulation iterations (default: 5000)
#' @param SEED Logical; if TRUE, sets random seed for reproducibility (default: TRUE)
#'
#' @return A matrix with \code{iter} rows and J(J-1)/2 columns, where J is the
#'   number of groups. Each row contains the p-values from all pairwise comparisons
#'   for one simulated dataset under the null hypothesis.
#'
#' @details
#' This function generates the null distribution of p-values for all pairwise
#' comparisons when all groups have the same binomial proportion \code{p}.
#'
#' For each iteration:
#' \itemize{
#'   \item Generates J binomial samples with sizes \code{n} and probability \code{p}
#'   \item Computes p-values for all J(J-1)/2 pairwise comparisons using
#'     \code{binmcp.sub}
#'   \item Stores the p-values in one row of the output matrix
#' }
#'
#' The resulting matrix can be used to determine critical p-values that control
#' the family-wise error rate at a desired level alpha. Typically, one would use
#' the alpha-quantile of the maximum or minimum p-values across comparisons.
#'
#' @note
#' This is an internal helper function. Users should typically use \code{\link{binmcp}}
#' which calls this function automatically when \code{p.crit=NULL}.
#'
#' @examples
#' \dontrun{
#' # Simulate null distribution for 4 groups
#' n <- c(30, 30, 30, 30)
#' p <- 0.5  # Common proportion under null
#'
#' # Generate null distribution
#' pv_matrix <- binmcp.crit(p, n, iter=5000)
#'
#' # Determine critical value for alpha=0.05
#' # (this is done automatically within binmcp)
#' }
#'
#' @family multiple comparison procedures
#' @family binomial methods
#' @keywords internal
#' @seealso \code{\link{binmcp}} for the main function that uses this
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

#' Multiple Comparison Procedure Using Combined P-values
#'
#' Performs a step-down multiple comparison procedure that combines p-values from
#' independent tests using Fisher's method, Chen-Nadarajah method, or the maximum
#' p-value method. Controls family-wise error rate.
#'
#' @param pv Vector of p-values from independent tests
#' @param alpha Family-wise error rate (default: 0.05)
#' @param opt Method for combining p-values:
#'   \itemize{
#'     \item 1 = Fisher's method (default)
#'     \item 2 = Chen-Nadarajah method
#'     \item 3 = Maximum p-value method
#'   }
#'
#' @return A list with components:
#'   \item{p.values}{Original p-values}
#'   \item{Decisions}{Vector of decisions: "Reject" or "Not Sig" for each test}
#'   \item{num.sig}{Number of significant tests}
#'
#' @details
#' This function implements the step-down multiple comparison method described in
#' Wilcox, R. R. & Clark, F. Robust multiple comparisons based on combined
#' probabilities from independent tests. Journal of Data Science.
#'
#' \strong{Combining Methods:}
#' \describe{
#'   \item{Fisher's method (opt=1)}{Combines p-values using -2*sum(log(p)).
#'     Under the null, this follows a chi-squared distribution with 2K degrees
#'     of freedom (K = number of tests)}
#'   \item{Chen-Nadarajah method (opt=2)}{Uses sum(qnorm(p/2)^2), which follows
#'     a chi-squared distribution with K degrees of freedom}
#'   \item{Maximum method (opt=3)}{Uses the maximum p-value, compared to a
#'     beta(K, 1) distribution}
#' }
#'
#' \strong{Step-down Procedure:}
#' The method tests hypotheses in a step-down fashion, starting with all K tests
#' and progressively removing the most significant ones. This can provide more
#' power than standard multiple comparison adjustments.
#'
#' The combined p-value is compared to alpha/(K+1-i) at each step i, providing
#' control of the family-wise error rate.
#'
#' @note
#' This method assumes the tests are independent. Results may be conservative
#' if tests are positively correlated.
#'
#' @examples
#' \dontrun{
#' # Combine p-values from 5 independent tests
#' pv <- c(0.01, 0.03, 0.15, 0.08, 0.25)
#'
#' # Using Fisher's method
#' result1 <- mcpPV(pv, alpha=0.05, opt=1)
#' print(result1$Decisions)
#'
#' # Using Chen-Nadarajah method
#' result2 <- mcpPV(pv, alpha=0.05, opt=2)
#' print(result2$Decisions)
#' }
#'
#' @family multiple comparison procedures
#' @family p-value adjustment methods
#' @seealso \code{\link{p.adjust}} for standard adjustment methods,
#'   \code{\link{mcpKadjp}} for k-FWER procedures
#' @export
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


#' K-FWER Multiple Comparison Procedures with P-value Adjustment
#'
#' Adjusts p-values for multiple comparisons using k-FWER (k familywise error rate)
#' controlling procedures. Provides more liberal control than standard FWER by
#' allowing up to k false rejections.
#'
#' @param p Vector of raw p-values to be adjusted
#' @param k The k value for k-FWER control (default: 1, equivalent to standard FWER).
#'   Setting k>1 allows more liberal testing by permitting up to k false rejections
#' @param proc Vector of procedure names to use. Choices are:
#'   \itemize{
#'     \item 'Holm' - Generalized Holm procedure (default)
#'     \item 'Hochberg' - Generalized Hochberg procedure (requires MTP2 condition)
#'     \item 'RS' - Romano-Shaikh procedure
#'     \item 'Sarkar' - Sarkar procedure (requires independent test statistics)
#'     \item 'BH' - Benjamini-Hochberg (FDR control)
#'   }
#' @param rawp Original raw p-values (default: same as p). Used for maintaining order
#'
#' @return A matrix with adjusted p-values. First column contains raw p-values,
#'   subsequent columns contain adjusted p-values for each requested procedure
#'
#' @details
#' This function implements k-FWER controlling procedures based on:
#' \itemize{
#'   \item Keselman, H. J., Miller, C. E., & Holland, B. (2011). Many tests of
#'     significance: New methods for controlling Type I errors. Psychological
#'     Methods, 16, 420-431.
#'   \item Keselman, H. J., & Miller, C. E. (2012). Correction to many tests of
#'     significance. Psychological Methods, 17(4), 679.
#' }
#'
#' \strong{k-FWER Control:}
#' Instead of controlling the probability of making any false rejections (FWER),
#' k-FWER controls the probability of making more than k false rejections. This
#' provides:
#' \itemize{
#'   \item More power than standard FWER when some false positives are tolerable
#'   \item Less conservative p-value adjustments
#'   \item Flexibility in balancing Type I and Type II errors
#' }
#'
#' \strong{Procedure Validity Conditions:}
#' \itemize{
#'   \item Generalized Hochberg: Valid under MTP2 condition of the joint null
#'     distribution of p-values
#'   \item Sarkar: Valid only for independent test statistics
#'   \item Holm: Valid under general dependence
#' }
#'
#' @note
#' When k=1, these procedures reduce to their standard FWER-controlling versions.
#' For k>1, they provide more liberal testing while still maintaining control
#' over the number of false rejections.
#'
#' @examples
#' \dontrun{
#' # Standard FWER control (k=1)
#' pvals <- c(0.001, 0.01, 0.03, 0.05, 0.10, 0.20)
#' result1 <- mcpKadjp(pvals, k=1, proc=c('Holm', 'Hochberg'))
#' print(result1)
#'
#' # More liberal k-FWER control (k=2, allow up to 2 false rejections)
#' result2 <- mcpKadjp(pvals, k=2, proc=c('Holm', 'RS'))
#' print(result2)
#' }
#'
#' @family multiple comparison procedures
#' @family p-value adjustment methods
#' @seealso \code{\link{p.adjust}} for standard methods,
#'   \code{\link{mcpPV}} for combined p-value methods
#' @export
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
#' Calculate D1 Values for Romano-Shaikh Procedure (Internal)
#'
#' Internal helper function to calculate D1 values used in the Romano-Shaikh
#' k-FWER procedure. This function is defined within \code{\link{mcpKadjp}}.
#'
#' @param k The k value for k-FWER
#' @param s Number of hypotheses (default: 1000)
#'
#' @return A list with components:
#'   \item{S}{Vector of cumulative values}
#'   \item{maxI}{Index of maximum value}
#'   \item{maxS}{Maximum S value (rounded to 4 decimals)}
#'
#' @keywords internal
#' @noRd
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


