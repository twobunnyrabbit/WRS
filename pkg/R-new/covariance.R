# covariance.R
#
# Covariance and Scatter Matrix Estimation Functions
# Part of WRS package refactoring
#
# This module contains functions for:
# - Robust covariance estimation (OGK, MVE, MCD, S-estimators, MBA)
# - Winsorized and trimmed covariance matrices
# - Skipped covariance (outlier removal)
# - Bootstrap covariance estimation
# - Distance covariance and ROC-based covariance
# - Mixed design covariance (between-within)
# - Median-based covariance methods
# - Utility functions for covariance matrices
#
# Extracted: 2025-12-30
# Number of functions: 43


#' Compute Biweight Midcovariance Between Two Variables
#'
#' Calculates the biweight midcovariance, a robust measure of covariation between
#' two variables that is resistant to outliers.
#'
#' @param x Numeric vector for first variable
#' @param y Numeric vector for second variable (must have same length as x)
#'
#' @return Numeric value representing the biweight midcovariance between x and y
#'
#' @details
#' The biweight midcovariance is a robust alternative to the classical covariance
#' that downweights the influence of outliers. The function:
#' \itemize{
#'   \item Centers both variables using the median
#'   \item Scales using MAD (median absolute deviation)
#'   \item Applies biweight weighting to points within threshold (|u| <= 1)
#'   \item Computes weighted covariance with down-weighting based on distance from center
#' }
#'
#' Points more than 9 * qnorm(.75) * MAD from the median receive zero weight.
#' This provides resistance to outliers while maintaining efficiency at normal
#' distributions.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{bicovm}} for biweight midcovariance matrix,
#'   \code{\link{wincov}} for winsorized covariance
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @examples
#' \dontrun{
#' # Compare robust vs classical covariance
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100)
#'
#' # Classical covariance
#' cov(x, y)
#'
#' # Robust biweight covariance
#' bicov(x, y)
#'
#' # Add outliers - biweight is more resistant
#' x2 <- c(x, 10, -10)
#' y2 <- c(y, 10, -10)
#' cov(x2, y2)    # Affected by outliers
#' bicov(x2, y2)  # More resistant
#' }
bicov<-function(x,y){
#
# compute biweight midcovariance of x and y
#
mx<-median(x)
my<-median(y)
ux<-abs((x-mx)/(9*qnorm(.75)*mad(x)))
uy<-abs((y-my)/(9*qnorm(.75)*mad(y)))
aval<-ifelse(ux<=1,1,0)
bval<-ifelse(uy<=1,1,0)
top<-sum(aval*(x-mx)*(1-ux^2)^2*bval*(y-my)*(1-uy^2)^2)
top<-length(x)*top
botx<-sum(aval*(1-ux^2)*(1-5*ux^2))
boty<-sum(bval*(1-uy^2)*(1-5*uy^2))
bi<-top/(botx*boty)
bi
}

#' Compute Skipped Covariance Matrix Using Projection Method
#'
#' Calculates a covariance matrix after removing multivariate outliers detected
#' using projection-based outlier detection (\code{\link{outpro}}).
#'
#' @param m Numeric matrix (n × p) of observations
#' @param STAND Logical; if TRUE (default), project data onto a sphere before
#'   computing projection distances
#'
#' @return A p × p covariance matrix computed from non-outlying observations
#'
#' @details
#' This function implements a skipped covariance estimator that:
#' \enumerate{
#'   \item Identifies multivariate outliers using the projection method
#'     (\code{\link{outpro}}):
#'     \itemize{
#'       \item Computes Donoho-Gasko median as center
#'       \item For each point, projects all data onto line connecting that point
#'         to the median
#'       \item Checks for outliers on each projection using boxplot rule
#'       \item Declares a point an outlier if it's outlying on any projection
#'     }
#'   \item Removes identified outliers
#'   \item Computes classical covariance matrix on remaining data
#' }
#'
#' The \code{STAND} parameter controls whether data are projected onto a sphere
#' before computing distances, which can improve detection with high-dimensional
#' data.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{outpro}} for the outlier detection method,
#'   \code{\link{skipcov}} for alternative skipped covariance methods,
#'   \code{\link{covout}} for covariance with general outlier removal
#'
#' @examples
#' \dontrun{
#' # Bivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(200), ncol = 2)
#' x <- rbind(x, c(5, 5), c(-5, -5))  # Add outliers
#'
#' # Classical covariance (affected by outliers)
#' cov(x)
#'
#' # Skipped covariance (outliers removed)
#' mscov(x)
#'
#' # Compare number of outliers detected
#' sum(!outpro(x, plotit = FALSE)$keep)
#' }
mscov<-function(m,STAND=TRUE){
#
# m is an n by p matrix
#
# Compute a skipped covariance matrix
#
# Eliminate outliers using a projection method
# That is, compute Donoho-Gasko median, for each point
# consider the line between it and the median,
# project all points onto this line, and
# check for outliers using a boxplot rule.
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# Eliminate any outliers and compute covariances
#  using remaining data.
#
m<-elimna(m)
temp<-outpro(m,plotit=FALSE,STAND=STAND)$keep
mcor<-var(m[temp,])
mcor
}

#' Estimate Covariance Between Dependent Order Statistics
#'
#' Computes the estimated covariance between order statistics (typically medians)
#' of two dependent variables, useful for constructing standard errors in paired
#' quantile comparisons.
#'
#' @param x Numeric vector for first variable
#' @param y Numeric vector for second variable (must have same length as x);
#'   if NA, function returns squared standard error for x
#' @param q Quantile to estimate (default: 0.5 for median)
#'
#' @return Estimated covariance between the q-th order statistics of x and y.
#'   If y = NA, returns squared standard error of the q-th quantile of x.
#'
#' @details
#' This function estimates the covariance between dependent order statistics,
#' particularly useful when comparing quantiles of paired data. The method:
#' \itemize{
#'   \item Computes the probability that both x and y fall in the same quadrant
#'     relative to their q-th quantiles
#'   \item Estimates marginal densities at the quantiles using \code{\link{akerd}}
#'   \item Combines these to form asymptotic covariance estimate
#' }
#'
#' When y = NA or x = y, the function returns the squared standard error of a
#' single order statistic, equivalent to calling \code{qse(x, q = q, op = 3)^2}.
#'
#' This is used internally by \code{\link{covmmed}} to construct covariance
#' matrices for median-based inference with dependent groups.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covmmed}} for covariance matrix of sample medians,
#'   \code{\link{qse}} for quantile standard errors,
#'   \code{\link{akerd}} for adaptive kernel density estimation
#'
#' @examples
#' \dontrun{
#' # Paired data
#' x <- rnorm(50)
#' y <- x + rnorm(50, sd = 0.5)
#'
#' # Covariance of sample medians
#' cov2med(x, y)
#'
#' # Squared standard error of median
#' cov2med(x)
#'
#' # For upper quartile
#' cov2med(x, y, q = 0.75)
#' }
cov2med<-function(x,y=NA,q=.5){
#
# Estimate the covariance between two dependent
# order statistics
# By default, q=.5 meaning that an estimate of
# of covariance is made when a single order statistic
# is used to estimate the median.
# y=NA, function returns squared standard error.
#
if(is.na(y[1]))val<-qse(x,q=q,op=3)^2
if(!is.na(y[1])){
if(sum((x-y)^2)==0)val<-qse(x,q=q,op=3)^2
if(sum((x-y)^2)>0){
n<-length(x)
m<-floor(q*n+.5)
yord<-sort(y)
flag<-(y<=yord[m])
xord<-sort(x)
xq<-xord[m]
yord<-sort(y)
yq<-yord[m]
flag1<-(x<=xq)
flag2<-(y<=yq)
A<-mean(flag1*flag2)
flag1<-(x<=xq)
flag2<-(y>yq)
B<-mean(flag1*flag2)
flag1<-(x>xq)
flag2<-(y<=yq)
C1<-mean(flag1*flag2)
flag1<-(x>xq)
flag2<-(y>yq)
D1<-mean(flag1*flag2)
fx<-akerd(x,pts=xq,plotit=FALSE,pyhat=TRUE)
fy<-akerd(y,pts=yq,plotit=FALSE,pyhat=TRUE)
v1<-(q-1)^2*A
v2<-(q-1)*q*B
v3<-(q-1)*q*C1
v4<-q*q*D1
val<-((v1+v2+v3+v4)/(fx*fy))/n
}}
val
}

#' Estimate Covariance Matrix for Sample Medians in Dependent Groups
#'
#' Computes the covariance matrix for sample medians (or other quantiles) across
#' multiple dependent groups, accounting for within-subject correlation.
#'
#' @param x Data in list or matrix format. If list, \code{x[[i]]} contains data
#'   for group i. If matrix, columns represent groups.
#' @param p Number of groups (default: length of x)
#' @param grp Vector specifying which groups to include (default: all groups).
#'   For example, \code{grp = c(1, 2, 5)} uses only groups 1, 2, and 5.
#' @param q Quantile to estimate (default: 0.5 for median)
#'
#' @return A p × p covariance matrix where p is the number of groups. Diagonal
#'   elements are squared standard errors of the quantiles; off-diagonal elements
#'   are covariances between quantiles of different groups.
#'
#' @details
#' This function estimates the covariance matrix for quantiles across dependent
#' (within-subject) groups. It:
#' \itemize{
#'   \item Uses \code{\link{cov2med}} to estimate covariances between order
#'     statistics of paired groups
#'   \item Requires equal sample sizes across all groups (paired design)
#'   \item Diagonal elements = \code{cov2med(x[[i]])} (squared SE)
#'   \item Off-diagonal elements = \code{cov2med(x[[i]], x[[j]])} (covariance)
#' }
#'
#' The resulting covariance matrix is used for inference on contrasts of medians
#' in repeated measures designs. Missing values are not allowed; all groups must
#' have complete paired observations.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{cov2med}} for pairwise covariance estimation,
#'   \code{\link{covmtrim}} for trimmed mean covariance matrix,
#'   \code{\link{Dqcomhd}} for quantile comparisons using this covariance
#'
#' @examples
#' \dontrun{
#' # Three dependent groups (e.g., three time points)
#' x <- list(
#'   rnorm(30),
#'   rnorm(30) + 0.5,
#'   rnorm(30) + 1.0
#' )
#'
#' # Covariance matrix of sample medians
#' covmmed(x)
#'
#' # Use only groups 1 and 3
#' covmmed(x, grp = c(1, 3))
#'
#' # For upper quartile instead of median
#' covmmed(x, q = 0.75)
#' }
covmmed<-function(x,p=length(x),grp=c(1:p),q=.5){
#
#  Estimate the covariance matrix for the sample medians
#  based on a SINGLE order statistic, using
#  the data in the R variable x.
# (x[[1]] contains the data for group 1, x[[2]] the data for group 2, etc.)
#  The function returns a p by p matrix of covariances, the diagonal
#  elements being equal to the squared standard error of the sample
#  trimmed means, where p is the number of groups to be included.
#  By default, all the groups in x are used, but a subset of
#  the groups can be used via grp.  For example, if
#  the goal is to estimate the covariances between the medians
#   for groups 1, 2, and 5, use the command grp<-c(1,2,5)
#  before calling this function.
#
#  Missing values (values stored as NA) are not allowed.
#
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("The data are not stored in a matrix or list mode.")
p<-length(grp)
pm1<-p-1
for (i in 1:pm1){
ip<-i+1
if(length(x[[grp[ip]]])!=length(x[[grp[i]]]))stop("The number of observations in each group must be equal")
}
n<-length(x[[grp[1]]])
covest<-matrix(0,p,p)
for(j in 1:p){
for(k in 1:p){
if(j==k)covest[j,j]<-cov2med(x[[grp[j]]],q=q)
if(j<k){
covest[j,k]<-cov2med(x[[grp[j]]],x[[grp[k]]],q=q)
covest[k,j]<-covest[j,k]
}}}
covest
}

#' Compute Rocke's TBS (Translated Biweight S-estimator) Covariance Matrix
#'
#' Calculates a robust covariance matrix using Rocke's TBS M-estimator, which
#' combines high breakdown point with good efficiency.
#'
#' @param x Numeric matrix (n × p) of observations
#'
#' @return A p × p robust covariance matrix
#'
#' @details
#' This function wraps the \code{covRob} function from the \pkg{robust} package
#' with \code{estim = "M"} to compute Rocke's TBS covariance estimator. The TBS
#' method:
#' \itemize{
#'   \item Achieves high breakdown point (resists contamination)
#'   \item Maintains good efficiency at normal distributions
#'   \item Uses constrained M-estimation with biweight loss function
#' }
#'
#' Requires the \pkg{robust} package to be installed.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{tbscov}} for TBS scatter estimator only,
#'   \code{\link{cov.roc}} for wrapper returning both center and covariance,
#'   \code{\link{Scov}} for Davies' S-estimator
#'
#' @references
#' Rocke, D.M. (1996). Robustness properties of S-estimators of multivariate
#' location and shape in high dimension. The Annals of Statistics, 24, 1327-1345.
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x <- rbind(x, matrix(5, nrow = 5, ncol = 3))  # Add outliers
#'
#' # Classical covariance
#' cov(x)
#'
#' # Rocke's robust covariance
#' covroc(x)
#' }
covroc<-function(x){
#
# compute Rocke's TBS covariance matrix
#
 library(robust)
temp<-covRob(x,estim="M")
val<-temp[2]$cov
val
}

#' Compute Covariance Matrix for Between-Within-Within Design
#'
#' Calculates trimmed covariance matrix for a three-way mixed design with one
#' between-subjects factor and two within-subjects factors.
#'
#' @param J Number of levels for between-subjects factor
#' @param K Number of levels for first within-subjects factor
#' @param L Number of levels for second within-subjects factor
#' @param x List of dependent variables in proper order (length J*K*L)
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming)
#'
#' @return A (J*K*L) × (J*K*L) block-diagonal covariance matrix. Non-zero blocks
#'   correspond to within-subjects covariances within each level of the between
#'   factor.
#'
#' @details
#' This function constructs a covariance matrix for analyzing three-way mixed
#' designs where:
#' \itemize{
#'   \item Factor 1 is between-subjects (J levels)
#'   \item Factor 2 is within-subjects (K levels)
#'   \item Factor 3 is within-subjects (L levels)
#' }
#'
#' The resulting matrix has block-diagonal structure because subjects in different
#' between-groups are independent. Within each block (size K*L × K*L), the
#' trimmed covariance matrix is computed using \code{\link{covmtrim}}.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{bbwcovm}} for between-between-within design,
#'   \code{\link{covmtrim}} for trimmed covariance estimation,
#'   \code{\link{bwwmcp}} for multiple comparisons in this design
#'
#' @examples
#' \dontrun{
#' # 2x3x2 design (2 between, 3x2 within)
#' # Generate data for 12 conditions
#' set.seed(123)
#' x <- lapply(1:12, function(i) rnorm(20))
#'
#' # Compute covariance matrix
#' C <- bwwcovm(J = 2, K = 3, L = 2, x = x)
#' dim(C)  # 12 x 12 matrix
#' }
bwwcovm<-function(J,K,L,x,tr=.2){
#
# compute covariance matrix for a between by within by within design
#
p=J*K*L
idep=K*L
mat=matrix(0,nrow=p,ncol=p)
id=c(1:idep)
for(j in 1:J){
mat[id,id]=covmtrim(x[id],tr=tr)
id=id+idep
}
mat
}

#' Compute Covariance Matrix for Between-Between-Within Design
#'
#' Calculates trimmed covariance matrix for a three-way mixed design with two
#' between-subjects factors and one within-subjects factor.
#'
#' @param J Number of levels for first between-subjects factor
#' @param K Number of levels for second between-subjects factor
#' @param L Number of levels for within-subjects factor
#' @param x List of dependent variables in proper order (length J*K*L)
#' @param tr Proportion of trimming (default: 0.2 for 20% trimming)
#'
#' @return A (J*K*L) × (J*K*L) block-diagonal covariance matrix. Non-zero blocks
#'   correspond to within-subjects covariances within each combination of the
#'   between factors.
#'
#' @details
#' This function constructs a covariance matrix for analyzing three-way mixed
#' designs where:
#' \itemize{
#'   \item Factor 1 is between-subjects (J levels)
#'   \item Factor 2 is between-subjects (K levels)
#'   \item Factor 3 is within-subjects (L levels)
#' }
#'
#' The resulting matrix has block-diagonal structure with J*K blocks (one per
#' combination of between-group levels), each of size L × L. Within each block,
#' the trimmed covariance matrix is computed using \code{\link{covmtrim}}.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{bwwcovm}} for between-within-within design,
#'   \code{\link{covmtrim}} for trimmed covariance estimation,
#'   \code{\link{bbwmcp}} for multiple comparisons in this design
#'
#' @examples
#' \dontrun{
#' # 2x2x3 design (2x2 between, 3 within)
#' # Generate data for 12 conditions
#' set.seed(123)
#' x <- lapply(1:12, function(i) rnorm(20))
#'
#' # Compute covariance matrix
#' C <- bbwcovm(J = 2, K = 2, L = 3, x = x)
#' dim(C)  # 12 x 12 matrix
#' }
bbwcovm<-function(J,K,L,x,tr=.2){
#
# compute covariance matrix for a between by between by within design
#
p=J*K*L
idep=L
mat=matrix(0,nrow=p,ncol=p)
id=c(1:idep)
for(j in 1:J){
for(k in 1:K){
mat[id,id]=covmtrim(x[id],tr=tr)
id=id+idep
}}
mat
}

#' Compute Biweight Midcovariance and Correlation Matrices
#'
#' Calculates robust biweight midcovariance and midcorrelation matrices for
#' multivariate data, providing resistance to outliers.
#'
#' @param x Data in matrix (n × p) or list format. If list, each element
#'   contains observations for one variable.
#'
#' @return A list with two components:
#'   \item{mcov}{p × p biweight midcovariance matrix}
#'   \item{mcor}{p × p biweight midcorrelation matrix}
#'
#' @details
#' This function computes robust covariance and correlation matrices using the
#' biweight midcovariance (\code{\link{bicov}}) for all pairs of variables. The
#' biweight estimator:
#' \itemize{
#'   \item Uses median-based centering for robustness
#'   \item Downweights outliers using biweight loss function
#'   \item Maintains high efficiency at normal distributions
#' }
#'
#' The correlation matrix is derived from the covariance matrix by standardizing:
#' \code{cor[i,j] = cov[i,j] / sqrt(cov[i,i] * cov[j,j])}.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{bicov}} for pairwise biweight covariance,
#'   \code{\link{bicovM}} for covariance matrix only,
#'   \code{\link{wincov}} for winsorized covariance
#'
#' @examples
#' \dontrun{
#' # Multivariate data
#' set.seed(123)
#' x <- matrix(rnorm(150), ncol = 3)
#'
#' # Classical covariance and correlation
#' cov(x)
#' cor(x)
#'
#' # Robust biweight versions
#' result <- bicovm(x)
#' result$mcov   # Biweight covariance
#' result$mcor   # Biweight correlation
#'
#' # With outliers - biweight is more resistant
#' x2 <- rbind(x, c(5, 5, 5))
#' bicovm(x2)$mcor
#' }
bicovm<-function(x){
#
# compute a biweight midcovariance matrix for the vectors of
# observations in x, where x is assumed to have list mode, or
# x is an n by p matrix
#
if(is.matrix(x)){
mcov<-matrix(0,ncol(x),ncol(x))
mcor<-matrix(0,ncol(x),ncol(x))
for (i in 1:ncol(x)){
for (j in 1:ncol(x))mcov[i,j]<-bicov(x[,i],x[,j])
}
}
if(is.list(x)){
mcov<-matrix(0,length(x),length(x))
mcor<-matrix(0,length(x),length(x))
for (i in 1:length(x)){
for (j in 1:length(x))mcov[i,j]<-bicov(x[[i]],x[[j]])
}
}
for (i in 1:ncol(mcov)){
for (j in 1:ncol(mcov))mcor[i,j]<-mcov[i,j]/sqrt(mcov[i,i]*mcov[j,j])
}
list(mcov=mcov,mcor=mcor)
}

#' Extract Biweight Midcovariance Matrix Only
#'
#' Wrapper function that returns only the biweight midcovariance matrix from
#' \code{\link{bicovm}}, discarding the correlation matrix.
#'
#' @param x Data in matrix (n × p) or list format
#'
#' @return A p × p biweight midcovariance matrix
#'
#' @details
#' This is a convenience function equivalent to \code{bicovm(x)$mcov}. Use this
#' when only the covariance matrix is needed and the correlation matrix can be
#' discarded.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{bicovm}} for both covariance and correlation matrices,
#'   \code{\link{bicov}} for pairwise biweight covariance
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(150), ncol = 3)
#'
#' # These are equivalent:
#' bicovM(x)
#' bicovm(x)$mcov
#' }
bicovM<-function(x){
M=bicovm(x)$mcov
M
}

#' Compute MVE (Minimum Volume Ellipsoid) Robust Covariance Matrix
#'
#' Calculates robust location and scatter estimates using the Minimum Volume
#' Ellipsoid (MVE) estimator, which finds the ellipsoid of smallest volume
#' covering half the data.
#'
#' @param x Numeric matrix (n × p) of observations
#'
#' @return A list with two components:
#'   \item{center}{p-dimensional robust location estimate (MVE center)}
#'   \item{cov}{p × p robust scatter matrix (MVE covariance)}
#'
#' @details
#' The MVE estimator provides a high breakdown point estimator (50%) that:
#' \itemize{
#'   \item Finds the ellipsoid of minimum volume containing at least h = ⌈(n+p+1)/2⌉
#'     observations
#'   \item Uses the sample mean and covariance of those h observations
#'   \item Rescales to achieve consistency at multivariate normal
#' }
#'
#' This function wraps \code{MASS::cov.mve} while preserving the random number
#' generator state. The MVE has lower efficiency than MCD but is sometimes more
#' robust in small samples.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{mvecov}} for covariance matrix only,
#'   \code{\link{covmcd}} for MCD estimator (generally preferred),
#'   \code{\link{covmba}} for MBA estimator
#'
#' @references
#' Rousseeuw, P.J. and van Zomeren, B.C. (1990). Unmasking multivariate outliers
#' and leverage points. Journal of the American Statistical Association, 85, 633-639.
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)  # Add outliers
#'
#' # Classical estimates (affected by outliers)
#' colMeans(x)
#' cov(x)
#'
#' # Robust MVE estimates
#' result <- covmve(x)
#' result$center  # Robust center
#' result$cov     # Robust covariance
#' }
covmve<-function(x){
oldSeed <- .Random.seed
val<-cov.mve(x)
assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
list(center=val$center,cov=val$cov)
}

#' Extract MVE Covariance Matrix Only
#'
#' Convenience function that returns only the MVE covariance matrix, discarding
#' the location estimate.
#'
#' @param x Numeric matrix (n × p) of observations
#'
#' @return A p × p MVE robust covariance matrix
#'
#' @details
#' This function is a wrapper for \code{MASS::cov.mve} that extracts only the
#' covariance component. Unlike \code{\link{covmve}}, it does not preserve the
#' random seed state.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covmve}} for both center and covariance,
#'   \code{\link{mcdcov}} for MCD covariance only
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#'
#' # MVE covariance only
#' mvecov(x)
#'
#' # Equivalent to:
#' covmve(x)$cov
#' }
mvecov<-function(x){
val<-cov.mve(x)
val$cov
}

#' Compute MCD (Minimum Covariance Determinant) Robust Covariance Matrix
#'
#' Calculates robust location and scatter estimates using the Minimum Covariance
#' Determinant (MCD) estimator, which finds the subset of observations with
#' smallest covariance determinant.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param nsamp Sampling method for candidate subsets:
#'   \itemize{
#'     \item \code{"sample"} (default): Uses min(5*p, 3000) random subsets
#'     \item \code{"best"}: Exhaustive enumeration up to 5000 subsets
#'   }
#'
#' @return A list with two components:
#'   \item{center}{p-dimensional robust location estimate (MCD center)}
#'   \item{cov}{p × p robust scatter matrix (MCD covariance)}
#'
#' @details
#' The MCD estimator provides a high breakdown point estimator (50%) that:
#' \itemize{
#'   \item Finds h observations (h ≈ n/2) whose covariance matrix has smallest
#'     determinant
#'   \item Uses the sample mean and covariance of those h observations
#'   \item Applies consistency and finite-sample correction factors
#' }
#'
#' This function wraps \code{MASS::cov.mcd} while preserving the random number
#' generator state. MCD generally has better efficiency than MVE and is the
#' preferred high breakdown estimator for most applications.
#'
#' The \code{nsamp} parameter controls the computational effort. Use \code{"best"}
#' for small datasets where exhaustive search is feasible, or \code{"sample"} for
#' larger datasets.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{mcdcov}} for covariance matrix only,
#'   \code{\link{DETMCD}} for deterministic MCD variant,
#'   \code{\link{covmve}} for MVE estimator,
#'   \code{\link{covmba}} for MBA estimator
#'
#' @references
#' Rousseeuw, P.J. and Van Driessen, K. (1999). A fast algorithm for the minimum
#' covariance determinant estimator. Technometrics, 41, 212-223.
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)  # Add outliers
#'
#' # Robust MCD estimates
#' result <- covmcd(x)
#' result$center  # Robust center
#' result$cov     # Robust covariance
#'
#' # Use exhaustive search for small dataset
#' result2 <- covmcd(x, nsamp = "best")
#' }
covmcd<-function(x,nsamp="sample"){
#
# nsamp="best" is the default used by R,
# meaning that  the number of samples is chosen so that
# exhaustive enumeration is done up to 5000 samples
# nsamp="sample" the number of samples
#  is min(5*p, 3000)
#
oldSeed <- .Random.seed
val<-cov.mcd(x,nsamp=nsamp)
assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
list(center=val$center,cov=val$cov)
}

#' Extract MCD Covariance Matrix Only
#'
#' Convenience function that returns only the MCD covariance matrix, discarding
#' the location estimate.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param nsamp Sampling method: \code{"sample"} (default) or \code{"best"}
#'
#' @return A p × p MCD robust covariance matrix
#'
#' @details
#' This function wraps \code{MASS::cov.mcd} and extracts only the covariance
#' component. It preserves the random seed state. See \code{\link{covmcd}} for
#' details on the MCD estimator and the \code{nsamp} parameter.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covmcd}} for both center and covariance,
#'   \code{\link{mvecov}} for MVE covariance only
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#'
#' # MCD covariance only
#' mcdcov(x)
#'
#' # Equivalent to:
#' covmcd(x)$cov
#' }
mcdcov<-function(x,nsamp="sample"){
#
# nsamp="best" is the default used by R,
# meaning that  the number of samples is chosen so that
# exhaustive enumeration is done up to 5000 samples
# nsamp="sample" the number of samples
#  is min(5*p, 3000)
#
#library(lqs)
oldSeed <- .Random.seed
val<-cov.mcd(x,nsamp=nsamp)
   assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
val$cov
}

#' Compute Gnanadesikan-Kettenring Robust Covariance
#'
#' Calculates robust covariance between two variables using the Gnanadesikan-
#' Kettenring estimator based on robust scale estimates of sums and differences.
#'
#' @param x Numeric vector for first variable
#' @param y Numeric vector for second variable (same length as x)
#' @param gk.sigmamu Function to compute robust scale (default: \code{\link{taulc}})
#' @param ... Additional arguments passed to \code{gk.sigmamu}
#'
#' @return Robust covariance estimate between x and y
#'
#' @details
#' The Gnanadesikan-Kettenring (GK) estimator computes covariance as:
#' \deqn{cov(x,y) = 0.25 * [scale(x+y)^2 - scale(x-y)^2]}
#'
#' This extends the familiar identity for variances to use robust scale estimates
#' instead of standard deviations. The method:
#' \itemize{
#'   \item Works with any robust scale estimator (default: tau scale)
#'   \item Inherits robustness properties of the scale estimator
#'   \item Forms the basis for OGK covariance matrix estimation
#' }
#'
#' This is used as a building block in \code{\link{covogk}} for multivariate
#' robust covariance estimation.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covogk}} for OGK covariance matrix using this estimator,
#'   \code{\link{taulc}} for default robust scale function,
#'   \code{\link{bicov}} for alternative robust covariance
#'
#' @references
#' Gnanadesikan, R. and Kettenring, J.R. (1972). Robust estimates, residuals, and
#' outlier detection with multiresponse data. Biometrics, 28, 81-124.
#'
#' Maronna, R.A. and Zamar, R.H. (2002). Robust estimates of location and
#' dispersion for high-dimensional datasets. Technometrics, 44, 307-317.
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100)
#'
#' # GK robust covariance
#' gkcov(x, y)
#'
#' # Compare to classical covariance
#' cov(x, y)
#'
#' # With outliers
#' x2 <- c(x, 10)
#' y2 <- c(y, 10)
#' gkcov(x2, y2)  # More resistant
#' cov(x2, y2)    # Affected by outliers
#' }
gkcov<-function(x,y,gk.sigmamu=taulc,...){
#
# Compute robust covariance using the Gnanadesikan-Kettenring
# estimator.
# (cf. Marrona & Zomar, 2002, Technometrics
#
val<-.25*(gk.sigmamu(x+y,...)-gk.sigmamu(x-y,...))
val
}

#' Compute OGK (Orthogonalized Gnanadesikan-Kettenring) Robust Covariance Matrix
#'
#' Calculates a robust covariance matrix using the Orthogonalized Gnanadesikan-
#' Kettenring (OGK) estimator, which combines pairwise robust covariances with
#' orthogonalization for computational efficiency.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param sigmamu Function to compute robust univariate scale and location
#'   (default: \code{\link{taulc}})
#' @param v Function to compute pairwise robust covariance (default: \code{\link{gkcov}})
#' @param n.iter Number of iterations (default: 5; often 1 is sufficient)
#' @param beta Quantile for reweighting (default: 0.9)
#' @param ... Additional arguments passed to \code{sigmamu} and \code{v}
#'
#' @return A p × p robust covariance matrix
#'
#' @details
#' The OGK estimator is a computationally efficient robust covariance estimator
#' that:
#' \enumerate{
#'   \item Computes robust pairwise covariances using the GK estimator
#'   \item Performs eigenvalue decomposition to orthogonalize variables
#'   \item Applies robust scale estimation to orthogonalized components
#'   \item Transforms back to original coordinates
#'   \item Optionally reweights observations for efficiency
#' }
#'
#' The method is particularly effective for high-dimensional data where
#' algorithms like MCD become computationally prohibitive. It has lower breakdown
#' point than MCD but better computational scaling.
#'
#' This function wraps \code{ogk.pairwise} and extracts the weighted covariance
#' matrix component.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{cov.ogk}} for wrapper with different interface,
#'   \code{\link{gkcov}} for pairwise GK covariance,
#'   \code{\link{skipogk}} for OGK-based skipped correlation,
#'   \code{\link{covmcd}} for higher breakdown alternative
#'
#' @references
#' Maronna, R.A. and Zamar, R.H. (2002). Robust estimates of location and
#' dispersion for high-dimensional datasets. Technometrics, 44, 307-317.
#'
#' @examples
#' \dontrun{
#' # High-dimensional data
#' set.seed(123)
#' x <- matrix(rnorm(500), ncol = 10)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 10)  # Add outliers
#'
#' # Classical covariance
#' cov(x)
#'
#' # OGK robust covariance
#' covogk(x)
#'
#' # Single iteration often sufficient
#' covogk(x, n.iter = 1)
#' }
covogk<-function(x,sigmamu=taulc,v=gkcov,n.iter=5,beta=.9,...){
#
# Compute robust (weighted) covariance matrix in Maronna and Zamar
# (2002, Technometrics, eq. 7).
#
# x is an n by p matrix
# n.iter number of iterations. 1 seems to be best
# sigmamu is any user supplied function having the form
#   sigmamu(x,mu.too=F) and which computes a robust measure of
#   of dispersion if mu.too=F. If mu.too=T, it returns
#   a robust measure of location as well.
# v is any robust covariance
#
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)  # remove any rows with missing data
temp<-ogk.pairwise(x,sigmamu=sigmamu,v=v,n.iter=n.iter,beta=beta,...)$wcovmat
temp
}

#' Compute Covariance After Removing Outliers
#'
#' Calculates classical covariance matrix after detecting and removing multivariate
#' outliers using a specified outlier detection method.
#'
#' @param x Numeric matrix (n × p) or vector. If vector and y is provided,
#'   combined into matrix.
#' @param y Optional numeric vector (same length as x)
#' @param outfun Function for outlier detection (default: \code{\link{outogk}})
#' @param plotit Logical; if TRUE, create diagnostic plot via \code{outfun}
#'
#' @return If bivariate (p=2), returns scalar covariance. Otherwise, returns
#'   p × p covariance matrix computed from non-outlying observations.
#'
#' @details
#' This function implements a skipped covariance estimator:
#' \enumerate{
#'   \item Detects multivariate outliers using \code{outfun}
#'   \item Removes flagged outliers
#'   \item Computes classical covariance on remaining data
#' }
#'
#' The default outlier detector (\code{\link{outogk}}) uses OGK-based Mahalanobis
#' distances. Other options include \code{\link{outpro}} for projection-based
#' detection or \code{\link{outmgv}} for MGV method.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{skipcov}} for skipped covariance with more options,
#'   \code{\link{outogk}} for default outlier detection,
#'   \code{\link{mscov}} for projection-based skipped covariance
#'
#' @examples
#' \dontrun{
#' # Bivariate data with outliers
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100)
#' x <- c(x, 5, 5)
#' y <- c(y, 5, -5)
#'
#' # Classical covariance (affected)
#' cov(x, y)
#'
#' # Skipped covariance (outliers removed)
#' covout(x, y)
#'
#' # Multivariate version
#' X <- cbind(x, y)
#' covout(X, plotit = TRUE)
#' }
covout<-function(x,y=NA,outfun=outogk,plotit=FALSE){
#
# Remove outliers and compute covariances
#
if(!is.na(y[1]))x<-cbind(x,y)
keep<-outfun(x,plotit=plotit)$keep
val<-var(x[keep,])
if(ncol(val)==2)val<-val[1,2]
val
}

#' Compute Skipped Correlation Using OGK Outlier Detection
#'
#' Calculates correlation matrix after removing outliers detected using OGK-based
#' Mahalanobis distances, with optional p-values for bivariate case.
#'
#' @param x Numeric matrix (n × p) or vector
#' @param y Optional numeric vector (same length as x)
#' @param plotit Logical; if TRUE, create diagnostic plot
#'
#' @return A list with components:
#'   \item{cor}{Correlation matrix computed from non-outlying observations}
#'   \item{test.stat}{Test statistic(s) for correlations}
#'   \item{p.value}{P-value (bivariate case only)}
#'   \item{crit.05}{Critical value for 0.05 significance level}
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Detects outliers using \code{\link{outogk}} (OGK-based Mahalanobis
#'     distances)
#'   \item Computes correlation on remaining observations
#'   \item For bivariate data, provides test of H₀: ρ = 0 with critical values
#' }
#'
#' The test statistic is \code{r*sqrt((n-2)/(1-r^2))} where r is the correlation
#' and n is the sample size after outlier removal. Critical values are based on
#' empirical relationships for skipped correlations.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{outogk}} for outlier detection method,
#'   \code{\link{covout}} for skipped covariance,
#'   \code{\link{skipcov}} for more general skipped covariance
#'
#' @references
#' Wilcox, R.R. (2022). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @examples
#' \dontrun{
#' # Bivariate data with outliers
#' set.seed(123)
#' x <- rnorm(100)
#' y <- 0.6 * x + rnorm(100)
#' x <- c(x, 5, -5)
#' y <- c(y, -5, 5)
#'
#' # Classical correlation
#' cor.test(x, y)
#'
#' # Skipped correlation
#' result <- skipogk(x, y)
#' result$cor
#' result$p.value
#'
#' # Multivariate version
#' X <- matrix(rnorm(300), ncol = 3)
#' skipogk(X)$cor
#' }
skipogk<-function(x,y=NA,plotit=FALSE){
#
# Remove outliers and compute correlations
#
if(!is.na(y[1]))x<-cbind(x,y)
x<-elimna(x)
n<-nrow(x)
keep<-outogk(x,plotit=plotit)$keep
val<-cor(x[keep,])
p.value<-NA
test<-NA
crit.05<-15.49/n+2.68
vat<-val
diag(vat)<-0
test<-abs(vat*sqrt((n-2)/(1-vat^2)))
diag(test)<-NA
if(ncol(val)==2){
p.value<-c("Greater than .1")
val<-val[1,2]
test<-abs(val*sqrt((n-2)/(1-val^2)))
crit<-4.8/n+2.72
if(test>=crit)p.value<-c("Less than .1")
crit<-15.49/n+2.68
if(test>=crit)p.value<-c("Less than .05")
crit<-14.22/n+3.26
if(test>=crit)p.value<-c("Less than .025")
crit<-24.83/n+3.74
if(test>=crit)p.value<-c("Less than .01")
}
list(cor=val,test.stat=test,p.value=p.value,crit.05=crit.05)
}

#' Compute MBA (Median Ball Algorithm) Robust Covariance Estimator
#'
#' Calculates robust location and scatter estimates using Olive's Median Ball
#' Algorithm (MBA), which finds the subset with smallest scatter around its median.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param csteps Number of concentration steps (default: 5)
#'
#' @return A list with two components:
#'   \item{center}{p-dimensional robust location estimate}
#'   \item{cov}{p × p robust scatter matrix}
#'
#' @details
#' The MBA estimator provides a computationally efficient robust estimator that:
#' \enumerate{
#'   \item Uses the marginal medians as starting point
#'   \item Concentrates by iteratively selecting observations with smallest
#'     Mahalanobis distances and recomputing center/scatter
#'   \item Compares two candidate solutions (DGK-based and median-based starts)
#'   \item Selects the solution with smallest determinant
#'   \item Rescales for consistency at multivariate normal
#' }
#'
#' MBA is faster than MCD/MVE and works well in moderate dimensions. It provides
#' a good balance between robustness and computational efficiency.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{rmba}} for reweighted MBA (generally preferred),
#'   \code{\link{cov.mba}} for covariance only,
#'   \code{\link{covmba2}} for original implementation,
#'   \code{\link{covmcd}} for MCD estimator
#'
#' @references
#' Olive, D.J. (2004). A resistant estimator of multivariate location and
#' dispersion. Computational Statistics & Data Analysis, 46, 93-102.
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)
#'
#' # MBA estimates
#' result <- covmba(x)
#' result$center
#' result$cov
#' }
covmba <- function(x, csteps = 5)
{  # gets the MBA estimator
        zx <- x
        x <- as.matrix(x)
	p <- dim(x)[2]
	##get the DGK estimator
	covs <- var(x)
	mns <- apply(x, 2, mean)	## concentrate
	for(i in 1:csteps) {
              	md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
             if(p > 1){
		mns <- apply(x[md2 <= medd2,  ], 2,
			mean)
		covs <- var(x[md2 <= medd2,  ])
             }
             if(p == 1){
		mns <- mean(x[md2 <= medd2])
		covs <- var(x[md2 <= medd2])
             }
	}
	covb <- covs
	mnb <- mns	##get the square root of det(covb)
	critb <- prod(diag(chol(covb)))
	##get the resistant estimator
	covv <- diag(p)
	med <- apply(x, 2, median)
	md2 <- mahalanobis(x, center = med, covv)
	medd2 <- median(md2)	## get the start
        if(p > 1){
	mns <- apply(x[md2 <= medd2,  ], 2, mean)
	covs <- var(x[md2 <= medd2,  ])
        }
        if(p == 1){
	mns <- mean(zx[md2 <= medd2])
	covs <- var(zx[md2 <= medd2])
        }
        ## concentrate
	for(i in 1:csteps) {
		md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
              if(p > 1){
		mns <- apply(x[md2 <= medd2,  ], 2,
			mean)
		covs <- var(x[md2 <= medd2,  ])
              }
              if(p == 1){
	        mns <- mean(zx[md2 <= medd2])
	        covs <- var(zx[md2 <= medd2])
               }
	}
	crit <- prod(diag(chol(covs)))
	if(crit < critb) {
		critb <- crit
		covb <- covs
		mnb <- mns
	}
##scale for better performance at MVN
	rd2 <- mahalanobis(x, mnb, covb)
	const <- median(rd2)/(qchisq(0.5, p))
	covb <- const * covb
	list(center = mnb, cov = covb)
}

#' Compute Skipped Covariance Matrix with Flexible Outlier Detection
#'
#' Calculates covariance matrix after removing outliers using projection-based or
#' MGV outlier detection with multiple options for robust centering.
#'
#' @param m Numeric matrix (n × p) of observations
#' @param cop Center option for projection method (default: 6):
#'   \itemize{
#'     \item 1: Donoho-Gasko median
#'     \item 2: MCD center
#'     \item 3: Marginal medians
#'     \item 4: MVE center
#'     \item 5: TBS center
#'     \item 6: RMBA (Olive's median ball algorithm)
#'   }
#' @param MM Logical; if FALSE (default), use boxplot rule for projection outliers;
#'   if TRUE, use MAD-median rule
#' @param op Outlier detection method (default: 1):
#'   \itemize{
#'     \item 1: Projection method (\code{\link{outpro}})
#'     \item 2: MGV method (\code{\link{outmgv}})
#'   }
#' @param mgv.op MGV center option (when op=2):
#'   \itemize{
#'     \item 0: Use all pairwise distances
#'     \item 1: MVE center
#'     \item 2: MCD center
#'   }
#' @param outpro.cop Center option passed to \code{outpro} (default: 3)
#' @param STAND Logical; if TRUE (default), standardize data before computing
#'   projection distances
#'
#' @return A p × p covariance matrix computed from non-outlying observations
#'
#' @details
#' This function provides flexible skipped covariance estimation with multiple
#' outlier detection options:
#'
#' \strong{Projection method (op=1):}
#' \itemize{
#'   \item Projects data onto lines through each point and robust center
#'   \item Flags outliers on each projection using boxplot or MAD rule
#'   \item Declares a point outlying if flagged on any projection
#' }
#'
#' \strong{MGV method (op=2):}
#' \itemize{
#'   \item Uses generalized variances from subsets of observations
#'   \item More computationally intensive but can be more powerful
#' }
#'
#' After outlier removal, classical covariance is computed on remaining data.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{outpro}} for projection-based outlier detection,
#'   \code{\link{outmgv}} for MGV outlier detection,
#'   \code{\link{mscov}} for simpler skipped covariance,
#'   \code{\link{mgvcov}} for MGV-specific interface,
#'   \code{\link{skip.cov}} for wrapper returning list format
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)
#'
#' # Skipped covariance with RMBA center
#' skipcov(x)
#'
#' # Use MCD center instead
#' skipcov(x, cop = 2)
#'
#' # Use MGV outlier detection
#' skipcov(x, op = 2)
#' }
skipcov<-function(m,cop=6,MM=FALSE,op=1,mgv.op=0,outpro.cop=3,STAND=TRUE){
#
# m is an n by p matrix
#
# Compute skipped covariance matrix
#
# op=1:
# Eliminate outliers using a projection method
# That is, first determine center of data using:
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)
#
# For each point
# consider the line between it and the center,
# project all points onto this line, and
# check for outliers using
#
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
#
# op=2 use mgv (function outmgv) method to eliminate outliers
#
# Eliminate any outliers and compute means
#  using remaining data.
# mgv.op=0, mgv uses all pairwise distances to determine center of the data
# mgv.op=1 uses MVE
# mgv.op=2 uses MCD
#
temp<-NA
m<-elimna(m)
m<-as.matrix(m)
if(op==2)temp<-outmgv(m,plotit=FALSE,op=mgv.op)$keep
if(op==1)temp<-outpro(m,plotit=FALSE,MM=MM,cop=outpro.cop,STAND=STAND,pr=FALSE)$keep
val<-var(m[temp,])
val
}

#' Compute Skipped Covariance Using MGV Outlier Detection
#'
#' Calculates covariance matrix after removing outliers detected using the
#' Minimum Generalized Variance (MGV) method.
#'
#' @param m Numeric matrix (n × p) of observations
#' @param op Center determination method (default: 1):
#'   \itemize{
#'     \item 0: Use all pairwise distances
#'     \item 1: Use robust center from \code{cov.fun}
#'   }
#' @param cov.fun Function to compute robust center/covariance (default: \code{\link{rmba}})
#' @param plotit Logical; if TRUE, create diagnostic plot
#'
#' @return A p × p covariance matrix computed from non-outlying observations
#'
#' @details
#' This function detects multivariate outliers using the MGV (Minimum Generalized
#' Variance) method via \code{\link{outmgv}}, then computes classical covariance
#' on the remaining observations.
#'
#' The MGV method identifies outliers by examining how much the generalized
#' variance (determinant of covariance matrix) changes when observations are removed.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{outmgv}} for MGV outlier detection,
#'   \code{\link{skipcov}} for more general skipped covariance,
#'   \code{\link{rmba}} for default robust center estimation
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)
#'
#' # MGV-based skipped covariance
#' mgvcov(x)
#'
#' # Use all pairwise distances
#' mgvcov(x, op = 0)
#' }
mgvcov<-function(m,op=1,cov.fun=rmba,plotit=FALSE){
#
# m is an n by p matrix
#
# Compute skipped covariance matrix
# using the MGV method
#
# Eliminate any outliers and compute covariance matrix
#  using remaining data.
# op=0, mgv uses all pairwise distances to determine center of the data
# op=1 uses the function indicated by the argument
# cov.fun to determine center of the data cloud.
# default is Olive's median ball
#
#
temp<-NA
m<-elimna(m)
temp<-outmgv(m,plotit=plotit,op=op,cov.fun=cov.fun)$keep
val<-var(m[temp,])
val
}

#' Compute Reweighted MBA (Median Ball Algorithm) Robust Covariance
#'
#' Calculates robust location, scatter, and correlation using the reweighted
#' Median Ball Algorithm (RMBA), which improves efficiency over MBA through
#' iterative reweighting.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param csteps Number of concentration steps (default: 5)
#' @param na.rm Logical; if TRUE (default), remove rows with missing values
#' @param plotit Logical; if FALSE (default), suppress plotting (used internally)
#'
#' @return A list with three components:
#'   \item{center}{p-dimensional robust location estimate}
#'   \item{cov}{p × p robust scatter matrix}
#'   \item{cor}{p × p robust correlation matrix (NULL if any variances <= 0)}
#'
#' @details
#' The RMBA improves upon the basic MBA (\code{\link{covmba}}) through reweighting:
#' \enumerate{
#'   \item Computes initial MBA estimate
#'   \item Identifies observations within chi-square 97.5% quantile
#'   \item Recomputes center/scatter using only those observations
#'   \item Rescales for consistency at multivariate normal
#'   \item Repeats reweighting step for additional efficiency
#' }
#'
#' The reweighting increases statistical efficiency while maintaining good
#' robustness properties. This is the recommended MBA variant for most applications.
#'
#' Code contributed by David J. Olive.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covmba}} for basic MBA,
#'   \code{\link{covmba2}} for alternative MBA implementation,
#'   \code{\link{cov.mba}} for covariance only
#'
#' @references
#' Olive, D.J. (2004). A resistant estimator of multivariate location and
#' dispersion. Computational Statistics & Data Analysis, 46, 93-102.
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)
#'
#' # RMBA estimates (preferred)
#' result <- rmba(x)
#' result$center
#' result$cov
#' result$cor
#' }
rmba<-function(x, csteps = 5,na.rm=TRUE,plotit=FALSE)
{
# computes the reweighted MBA estimator
# Code supplied by David Olive
#
#       x is assumed to be a matrix
#
#  plotit=FALSE is used to avoid problems when this function is called
#  by other function in WRS
#
x=as.matrix(x)
if(na.rm)x=elimna(x)
	p <- dim(x)[2]
	n <- dim(x)[1]	##get the DGK estimator
	covs <- var(x)
	mns <- apply(x, 2, mean)	## concentrate
	for(i in 1:csteps) {
		md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		covs <- var(x[md2 <= medd2,  ])
	}
	covb <- covs
	mnb <- mns	##get the square root of det(covb)
	critb <- prod(diag(chol(covb)))	##get the resistant estimator
	covv <- diag(p)
	med <- apply(x, 2, median)
	md2 <- mahalanobis(x, center = med, covv)
	medd2 <- median(md2)	## get the start
	mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
	covs <- var(x[md2 <= medd2,  ])	## concentrate
	for(i in 1:csteps) {
		md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
#		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		covs <- var(x[md2 <= medd2,  ])
	}
	crit <- prod(diag(chol(covs)))
	if(crit < critb) {
		critb <- crit
		covb <- covs
		mnb <- mns
	}
##scale for better performance at MVN
	rd2 <- mahalanobis(x, mnb, covb)
	const <- median(rd2)/(qchisq(0.5, p))
	covb <- const * covb
	##reweight the above MBA estimator (mnb,covb) for efficiency
	rd2 <- mahalanobis(x, mnb, covb)
	up <- qchisq(0.975, p)
	rmnb <- apply(as.matrix(x[rd2 <= up,  ]), 2, mean)
	rcovb <- var(x[rd2 <= up,  ])
	rd2 <- mahalanobis(x, rmnb, rcovb)
	const <- median(rd2)/(qchisq(0.5, p))
	rcovb <- const * rcovb	## reweight again
	rd2 <- mahalanobis(x, rmnb, rcovb)
	up <- qchisq(0.975, p)
	rmnb <- apply(as.matrix(x[rd2 <= up,  ]), 2, mean)
	rcovb <- var(x[rd2 <= up,  ])
	rd2 <- mahalanobis(x, rmnb, rcovb)
	const <- median(rd2)/(qchisq(0.5, p))
	rcovb <- const * rcovb
cor.b=NULL
temp=outer(sqrt(diag(rcovb)),sqrt(diag(rcovb)),'*')
if(min(diag(rcovb)>0))cor.b=rcovb/temp
	list(center = rmnb, cov = rcovb, cor=cor.b)
}

#' Compute Rocke's TBS S-Estimator Covariance Matrix
#'
#' Calculates robust covariance using Rocke's constrained S-estimator (TBS),
#' which achieves high breakdown point with good efficiency.
#'
#' @param x Numeric matrix (n × p) with p >= 2
#' @param eps Convergence tolerance (default: 1e-3)
#' @param maxiter Maximum number of iterations (default: 20)
#' @param r Breakdown point (default: 0.45 for 45% breakdown)
#' @param alpha Asymptotic rejection probability (default: 0.05)
#'
#' @return A p × p robust covariance matrix
#'
#' @details
#' The TBS (Translated Biweight S-estimator) is a constrained S-estimator that:
#' \itemize{
#'   \item Starts with MVE estimate as initial values
#'   \item Iteratively reweights using biweight loss function
#'   \item Uses iteratively reweighted algorithm to minimize robust scale
#'   \item Achieves specified breakdown point while maintaining efficiency
#' }
#'
#' This function returns only the covariance matrix. For both location and
#' scatter, use the \code{tbs} function instead.
#'
#' The algorithm stops when the change in estimates is less than \code{eps} or
#' after \code{maxiter} iterations.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covroc}} for wrapper using robust package,
#'   \code{\link{Scov}} for Davies' S-estimator,
#'   \code{\link{covmcd}} for MCD estimator
#'
#' @references
#' Rocke, D.M. (1996). Robustness properties of S-estimators of multivariate
#' location and shape in high dimension. The Annals of Statistics, 24, 1327-1345.
#'
#' @examples
#' \dontrun{
#' # Multivariate data
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)
#'
#' # TBS covariance
#' tbscov(x)
#'
#' # With higher breakdown point
#' tbscov(x, r = 0.49)
#' }
tbscov <- function(x,eps=1e-3,maxiter=20,r=.45,alpha=.05){
#        Rocke's contrained s-estimator
#   returns covariance matrix only. For both locatiion and scatter, use tbs
#
#      r=.45 is the breakdown point
#      alpha=.05 is the asymptotic rejection probability.
#
if(!is.matrix(x))stop("x should be a matrix with two or more columns")
x<-elimna(x)
temp<-cov.mve(x)
t1<-temp$center
s<-temp$cov
    n <- nrow(x)
    p <- ncol(x)
if(p==1)stop("x should be a matrix with two or more columns")
c1M<-cgen.bt(n,p,r,alpha,asymp=FALSE)
c1<-c1M$c1
if(c1==0)c1<-.001 #Otherwise get division by zero
M<-c1M$M
    b0 <- erho.bt(p,c1,M)
    crit <- 100
    iter <- 1
    w1d <- rep(1,n)
    w2d <- w1d
    while ((crit > eps)&(iter <= maxiter))
    {
        t.old <- t1
        s.old <- s
        wt.old <- w1d
        v.old <- w2d
        d2 <- mahalanobis(x,center=t1,cov=s)
        d <- sqrt(d2)
        k <- ksolve.bt(d,p,c1,M,b0)
        d <- d/k
        w1d <- wt.bt(d,c1,M)
        w2d <- v.bt(d,c1,M)
        t1 <- (w1d %*% x)/sum(w1d)
        s <- s*0
        for (i in 1:n)
        {
            xc <- as.vector(x[i,]-t1)
            s <- s + as.numeric(w1d[i])*(xc %o% xc)
        }
        s <- p*s/sum(w2d)
        mnorm <- sqrt(as.vector(t.old) %*% as.vector(t.old))
        snorm <- eigen(s.old)$values[1]
        crit1 <- max(abs(t1 - t.old))
#        crit <- max(crit1,crit2)
        crit <- max(abs(w1d-wt.old))/max(w1d)
        iter <- iter+1
    }
#    mnorm <- sqrt(as.vector(t1) %*% as.vector(t1))
#    snorm <- eigen(s)$values[1]
#    return(list(t1=t1,s=s))
s
}

#' Extract MBA Covariance or Correlation Matrix
#'
#' Convenience wrapper for \code{\link{covmba2}} that extracts covariance or
#' correlation matrix only.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param COR Logical; if TRUE, return correlation matrix; if FALSE (default),
#'   return covariance matrix
#'
#' @return A p × p covariance matrix (if COR=FALSE) or correlation matrix (if COR=TRUE)
#'
#' @details
#' This is a wrapper around \code{\link{covmba2}} that allows easy extraction of
#' either covariance or correlation matrix from the MBA estimator.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covmba2}} for full MBA implementation,
#'   \code{\link{covmba}} for MBA with both starts,
#'   \code{\link{rmba}} for reweighted MBA
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#'
#' # MBA covariance
#' cov.mba(x)
#'
#' # MBA correlation
#' cov.mba(x, COR = TRUE)
#' }
cov.mba<-function(x,COR=FALSE){
val<-covmba2(x)$cov
if(COR){
val=val/outer(sqrt(diag(val)),sqrt(diag(val)))
}
val
}

#' Compute MBA Covariance Matrix (Alternative Implementation)
#'
#' Alternative implementation of the Median Ball Algorithm covariance estimator,
#' used internally by \code{\link{cov.mba}}.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param csteps Number of concentration steps (default: 5)
#'
#' @return A list with two components:
#'   \item{center}{p-dimensional robust location estimate}
#'   \item{cov}{p × p robust scatter matrix}
#'
#' @details
#' This function implements the basic MBA algorithm similar to \code{\link{covmba}}.
#' See \code{\link{covmba}} for algorithm details. For most applications,
#' \code{\link{rmba}} (reweighted MBA) is preferred for better efficiency.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covmba}} for MBA with both starts,
#'   \code{\link{rmba}} for reweighted MBA (preferred),
#'   \code{\link{cov.mba}} for wrapper extracting covariance only
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#' result <- covmba2(x)
#' result$center
#' result$cov
#' }
covmba2<-function(x, csteps = 5)
{
# Perform the median ball algorithm.
#
# It returns a measure of location and scatter for the
# multivariate data in x, which is assumed to have
# p>-2 column and n rows.
#
# This code is based on a very slight modificatiion of code originally
# written by David Olive
#
x<-as.matrix(x)
if(!is.matrix(x))stop("x should be a matrix")
         p <- dim(x)[2]
#if(p==1)stop("x should be a matrix with two or more columns of variables")
         ##get the DGK estimator
         covs <- var(x)
         mns <- apply(x, 2, mean)        ## concentrate
         for(i in 1:csteps) {
                 md2 <- mahalanobis(x, mns, covs)
                 medd2 <- median(md2)
#                 mns <- apply(x[md2 <= medd2,  ], 2,
                 mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2,
                         mean)
                 covs <- var(x[md2 <= medd2,  ])
         }
         covb <- covs
         mnb <- mns      ##get the square root of det(covb)
         critb <- prod(diag(chol(covb)))
         ##get the resistant estimator
         covv <- diag(p)
         med <- apply(x, 2, median)
         md2 <- mahalanobis(x, center = med, covv)
         medd2 <- median(md2)    ## get the start
#         mns <- apply(x[md2 <= medd2,  ], 2, mean)
         mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
         covs <- var(x[md2 <= medd2,  ]) ## concentrate
         for(i in 1:csteps) {
                 md2 <- mahalanobis(x, mns, covs)
                 medd2 <- median(md2)
 #                mns <- apply(x[md2 <= medd2,  ], 2,mean)
       mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
                 covs <- var(x[md2 <= medd2,  ])
         }
         crit <- prod(diag(chol(covs)))
         if(crit < critb) {
                 critb <- crit
                 covb <- covs
                 mnb <- mns
         }
##scale for better performance at MVN
         rd2 <- mahalanobis(x, mnb, covb)
         const <- median(rd2)/(qchisq(0.5, p))
         covb <- const * covb
         list(center = mnb, cov = covb)
}

#' Compute Donoho-Gasko Depth-Based Scatter Matrix
#'
#' Calculates a robust scatter matrix using depth-based trimming, where observations
#' are weighted according to their halfspace depth.
#'
#' @param m Numeric matrix (n × p) or list of vectors
#' @param tr Proportion of trimming based on depth (default: 0.2). Must be < 0.5.
#' @param dop Depth computation method (default: 1):
#'   \itemize{
#'     \item 1: Use \code{\link{fdepth}} (faster approximate depth for p > 2)
#'     \item 2: Use \code{\link{fdepthv2}} (slower but more accurate)
#'   }
#' @param cop Center option for depth computation (default: 2):
#'   \itemize{
#'     \item 1: Donoho-Gasko median
#'     \item 2: MCD center
#'     \item 3: Marginal medians
#'     \item 4: MVE center
#'   }
#' @param pr Logical; if TRUE (default), print warning if only one point remains
#'
#' @return For p=1, returns trimmed mean. For p>1, returns p × p scatter matrix
#'   computed from observations with depth >= tr.
#'
#' @details
#' This function computes a depth-based scatter estimator:
#' \enumerate{
#'   \item Computes halfspace depth for each observation
#'   \item Trims observations with lowest depths (proportion \code{tr})
#'   \item Computes classical covariance on remaining high-depth observations
#' }
#'
#' For bivariate data (p=2), exact depth is computed using \code{\link{depth}}.
#' For higher dimensions, approximate depth is used via \code{\link{fdepth}} or
#' \code{\link{fdepthv2}}.
#'
#' Depth-based trimming provides affine-equivariant robust estimation.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{fdepth}} for depth computation,
#'   \code{\link{depth}} for exact bivariate depth,
#'   \code{\link{covmcd}} for MCD covariance
#'
#' @references
#' Donoho, D.L. and Gasko, M. (1992). Breakdown properties of location estimates
#' based on halfspace depth and projected outlyingness. The Annals of Statistics,
#' 20, 1803-1827.
#'
#' @examples
#' \dontrun{
#' # Multivariate data
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)
#'
#' # Depth-based scatter (20% trimming)
#' dcov(x)
#'
#' # More aggressive trimming
#' dcov(x, tr = 0.3)
#'
#' # Use exact depth (slower)
#' dcov(x, dop = 2)
#' }
dcov<-function(m,tr=.2,dop=1,cop=2,pr=TRUE){
#
# Compute multivariate measure of scatter
# using Donoho-Gasko method.
#
# dop=1, use fdepth to compute depths
# dop=2, use fdepthv2  to compute depths
# which is slower but more accurate.
#
# cop=1, use Donoho-Gasko median in fdepth
# cop=2, use MCD in fdepth
# cop=3, use marginal medians in fdepth
# cop=4, use MVE in fdepth
#
if(tr>=.5)stop("Amount of trimming must be less than .5")
if(is.list(m))m<-matl(m)
if(!is.matrix(m))stop("Data must be stored in a matrix or in list mode.")
if(ncol(m)==1){
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
}}
mdep<-max(temp)
flag<-(temp==mdep)
flag2<-(temp>=tr)
if(sum(flag2)==0)stop("Trimmed all of the data")
if(sum(flag2)==1){
if(pr)print("Warning: Trimmed all but one point")
val<-0
}
if(sum(flag2)>1)val<-var(m[flag2,])
val
}

#' Compute Davies' S-Estimator for Covariance
#'
#' Calculates robust location and scatter using Davies' S-estimator, which
#' minimizes a robust scale measure of multivariate residuals.
#'
#' @param m Numeric matrix (n × p) of observations
#'
#' @return A list with two components:
#'   \item{center}{p-dimensional robust location estimate}
#'   \item{cov}{p × p robust scatter matrix}
#'
#' @details
#' This function wraps \code{CovSest} from the \pkg{rrcov} package to compute
#' Davies' S-estimator. The S-estimator:
#' \itemize{
#'   \item Minimizes a robust M-scale of Mahalanobis distances
#'   \item Achieves high breakdown point (up to 50%)
#'   \item Has high computational complexity (slower than MCD)
#' }
#'
#' Requires the \pkg{rrcov} package to be installed.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{tbscov}} for Rocke's TBS S-estimator,
#'   \code{\link{covmcd}} for MCD estimator (generally faster),
#'   \code{\link{covmve}} for MVE estimator
#'
#' @references
#' Davies, P.L. (1987). Asymptotic behaviour of S-estimates of multivariate
#' location parameters and dispersion matrices. The Annals of Statistics, 15,
#' 1269-1292.
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)
#'
#' # Davies' S-estimator
#' result <- Scov(x)
#' result$center
#' result$cov
#' }
Scov<-function(m){
#
# Compute Davies' covariance S-estimator
#
library("rrcov")
res=CovSest(m)
center=res@center
z=res@cov
list(center=center,cov=z)
}

#' Compute OGK Covariance Matrix (Alternative Interface)
#'
#' Alternative interface to \code{\link{covogk}} that works with vectors or
#' matrices and wraps the \code{ogk} function.
#'
#' @param x Numeric matrix (n × p) or vector. If vector and y provided, combined
#'   into matrix.
#' @param y Optional numeric vector (same length as x)
#' @param n.iter Number of iterations (default: 1, which is often sufficient)
#' @param sigmamu Function to compute robust scale and location (default: \code{\link{taulc}})
#' @param v Function to compute pairwise robust covariance (default: \code{\link{gkcov}})
#' @param beta Quantile for reweighting (default: 0.9)
#' @param ... Additional arguments passed to \code{sigmamu} and \code{v}
#'
#' @return A p × p robust covariance matrix
#'
#' @details
#' This is an alternative interface to the OGK (Orthogonalized Gnanadesikan-
#' Kettenring) estimator. See \code{\link{covogk}} for detailed algorithm description.
#'
#' This function wraps the internal \code{ogk} function and extracts the covariance
#' component. It accepts both vector pairs and matrices.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covogk}} for primary OGK interface,
#'   \code{\link{gkcov}} for pairwise GK covariance,
#'   \code{\link{skipogk}} for OGK-based skipped correlation
#'
#' @examples
#' \dontrun{
#' # High-dimensional data
#' set.seed(123)
#' x <- matrix(rnorm(500), ncol = 10)
#'
#' # OGK covariance
#' cov.ogk(x)
#'
#' # Bivariate version
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100)
#' cov.ogk(x, y)
#' }
cov.ogk<-function(x,y=NA,n.iter=1,sigmamu=taulc,v=gkcov,beta=.9,...){
#
# Compute robust (weighted) covariance matrix in Maronna and Zamar
# (2002, Technometrics, eq. 7).
#
# n.iter number of iterations. 1 seems to be best
# sigmamu computes a robust measure of location and scale for
#  data stored in a single vector.
#  v robust correlation coefficient
#  estloc, a robust measure of location
#
if(!is.na(y[1]))x<-cbind(x,y)
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)
n<-nrow(x)
p<-ncol(x)
val<-matrix(NA,p,p)
temp<-ogk(x,sigmamu=sigmamu,v=v,n.iter=n.iter,beta=beta,...)$cov
temp
}

#' Compute Winsorized Covariance Matrix
#'
#' Calculates covariance matrix after winsorizing each variable, providing
#' robustness to outliers while maintaining all observations.
#'
#' @param m Numeric matrix (n × p) of observations
#' @param tr Proportion of winsorization (default: 0.2 for 20% winsorization on
#'   each tail)
#'
#' @return A p × p winsorized covariance matrix
#'
#' @details
#' This function winsorizes each variable separately using \code{\link{winall}},
#' then computes the classical covariance matrix of the winsorized data.
#' Winsorization replaces extreme values with less extreme values rather than
#' removing them entirely.
#'
#' For each variable, the lowest \code{tr} proportion of values are set equal to
#' the (\code{tr}*n)th order statistic, and the highest \code{tr} proportion are
#' set equal to the ((1-\code{tr})*n)th order statistic.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{wincovN}} for normalized winsorized covariance,
#'   \code{\link{wmean.cov}} for winsorized mean and covariance together,
#'   \code{\link{winall}} for underlying winsorization,
#'   \code{\link{bicov}} for biweight covariance
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1, ] <- c(10, 10, 10)  # Add outlier
#'
#' # Classical covariance
#' cov(x)
#'
#' # Winsorized covariance
#' wincov(x)
#' }
wincov<-function(m,tr=.2){
m=winall(m,tr=tr)$cov
m
}

#' Convert Long Format Data to Covariance Matrix Format
#'
#' Transforms longitudinal data in long format (one row per observation) to list
#' format where each subject's covariate data is stored as a matrix.
#'
#' @param x Data frame or matrix in long format
#' @param Sid.col Column number containing subject IDs
#' @param dep.col Vector of column numbers containing the covariates/dependent
#'   variables
#'
#' @return A list where each element contains a matrix of covariate data for one
#'   subject. Each matrix has dimensions (m × p) where m is the number of
#'   observations for that subject and p is the number of covariates.
#'
#' @details
#' This function is useful for converting repeated measures or longitudinal data
#' from "long" format to the list-of-matrices format required by many WRS
#' covariance and ANOVA functions.
#'
#' Requirements:
#' \itemize{
#'   \item All subjects must have the same number of observations
#'   \item Subject IDs in \code{Sid.col} identify which rows belong to each subject
#'   \item \code{dep.col} specifies which columns contain the variables of interest
#' }
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covmtrim}} for trimmed covariance on dependent groups,
#'   \code{\link{bwwcovm}} for mixed design covariance matrices
#'
#' @examples
#' \dontrun{
#' # Create long format data
#' # 3 subjects, 2 time points, 2 variables
#' df <- data.frame(
#'   subj = rep(1:3, each = 2),
#'   time = rep(1:2, 3),
#'   var1 = rnorm(6),
#'   var2 = rnorm(6)
#' )
#'
#' # Convert to list format
#' # Extract variables from columns 3-4 for each subject in column 1
#' result <- longcov2mat(df, Sid.col = 1, dep.col = c(3, 4))
#' # result[[1]] is 2x2 matrix for subject 1
#' # result[[2]] is 2x2 matrix for subject 2
#' # result[[3]] is 2x2 matrix for subject 3
#' }
longcov2mat<-function(x,Sid.col,dep.col){
#
# Have data in a matrix or data frame, x
# Sid.col indicates Subject's id
# Here, each subject has one or more rows of data
#
# In a regression setting, each subject has
# one or more covariates corresponding to columns.
# For example, two covariates might be stored in columns
# 3 and 6.
#
# Goal: For ith subject, store the covariate data in
# list mode, which is a matrix.
# So for ith subject, store covariate data in z[[i]], say, which
# contains a matrix of dimension  m by p,
# m is the number of observations for ith subject and p
# the number of covariates.
#
# dep.col, having length p, indicates columns containe the covariates
# Column Sid.col indicates the column containing subject's id
#
if(is.null(dim(x)))stop("x must be a matrix or data frame")
Sid=unique(x[,Sid.col])
res=list()
nid=length(Sid)
p=length(dep.col)# Number of covariates for each subject
n=nrow(x)
flag=(x[,Sid.col]==Sid[1])
n.each.s=sum(flag) # the number of rows for each subject
ns=n/n.each.s # the number of  subjects
if(!is.wholenumber(ns))stop("Not all S's have same number of rows of data")
for(i in 1:ns){
#res[[i]]=matrix(NA,nrow=n.each.s,ncol=p)
flag=(x[,Sid.col]==Sid[i])
res[[i]]=as.matrix(x[flag,dep.col])
}
res
}

#' Compute Rocke's ROC Covariance with Center and Scatter
#'
#' Wrapper for Rocke's TBS M-estimator that returns both robust center and
#' covariance matrix.
#'
#' @param x Numeric matrix (n × p) of observations
#'
#' @return A list with two components:
#'   \item{center}{p-dimensional robust location estimate}
#'   \item{cov}{p × p robust scatter matrix}
#'
#' @details
#' This function wraps \code{covRob} from the \pkg{robust} package with
#' \code{estim = 'M'}. See \code{\link{covroc}} for details on the TBS estimator.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covroc}} for covariance only,
#'   \code{\link{tbscov}} for direct TBS implementation
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#' result <- cov.roc(x)
#' result$center
#' result$cov
#' }
cov.roc<-function(x){
library(robust)
temp<-covRob(x,estim='M')
val<-temp
list(center=val[3]$center,cov=val[2]$cov)
}

#' Wrap Covariance Matrix in List Format
#'
#' Simple wrapper that places a covariance matrix in list format for compatibility
#' with functions expecting list input.
#'
#' @param m Covariance matrix (p × p)
#'
#' @return A list with one component \code{cov} containing the input matrix
#'
#' @keywords internal
cov.funl<-function(m){
list(cov=m)
}

#' Wrapper for Skipped Covariance Returning List Format
#'
#' Wrapper for \code{\link{skipcov}} that returns the covariance matrix in list
#' format for compatibility with functions expecting that structure.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param cop Center option (default: 6 for RMBA)
#' @param MM Logical; projection outlier rule
#' @param op Outlier detection method (1: projection, 2: MGV)
#' @param mgv.op MGV center option
#' @param outpro.cop Center option for outpro
#' @param STAND Logical; standardize before projection
#'
#' @return A list with one component \code{cov} containing the skipped covariance matrix
#'
#' @details
#' See \code{\link{skipcov}} for detailed parameter descriptions and algorithm details.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{skipcov}} for direct skipped covariance
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#' result <- skip.cov(x)
#' result$cov
#' }
skip.cov<-function(x,cop = 6, MM = FALSE, op = 1, mgv.op = 0, outpro.cop = 3,
    STAND = FALSE){
ans=skipcov(x,cop=cop,MM=MM,op=op,mgv.op=mgv.op,outpro.cop=outpro.cop,STAND=STAND)
list(cov=ans)
}

#' Compute Winsorized Mean and Covariance Matrix
#'
#' Calculates both winsorized mean vector and winsorized covariance matrix,
#' returning them together in list format.
#'
#' @param x Numeric matrix (n × p) of observations
#' @param tr Proportion of winsorization (default: 0 for no winsorization)
#'
#' @return A list with two components:
#'   \item{center}{p-dimensional vector of winsorized means}
#'   \item{cov}{p × p winsorized covariance matrix}
#'
#' @details
#' This function computes winsorized location and scatter estimates:
#' \itemize{
#'   \item The center is computed by applying winsorized mean to each variable
#'   \item The covariance is computed via \code{\link{wincov}}
#' }
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{wincov}} for winsorized covariance only,
#'   \code{\link{winmean}} for univariate winsorized mean
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#' result <- wmean.cov(x, tr = 0.2)
#' result$center  # Winsorized means
#' result$cov     # Winsorized covariance
#' }
wmean.cov<-function(x,tr=0){
#
# Compute Winsoriced mean and covariance for data in x
#
loc=apply(x,2,mean,tr=tr)
cv=wincov(x,tr=tr)
list(center=loc,cov=cv)
}

#' Wrap Classical Covariance in List Format
#'
#' Simple wrapper that computes classical covariance and returns it in list format
#' for compatibility with functions expecting list input.
#'
#' @param x Numeric matrix (n × p) of observations
#'
#' @return A list with one component \code{cov} containing the classical covariance matrix
#'
#' @export
#' @family covariance estimation
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#' result <- covl(x)
#' result$cov  # Same as cov(x)
#' }
covl<-function(x){
res=cov(x)
list(cov=res)
}

#' Convert Covariance Matrix to Correlation Matrix
#'
#' Transforms a covariance matrix into the corresponding correlation matrix by
#' standardizing.
#'
#' @param x A p × p covariance matrix
#'
#' @return A p × p correlation matrix
#'
#' @details
#' Computes correlations as \code{cor[i,j] = cov[i,j] / sqrt(cov[i,i] * cov[j,j])}.
#'
#' Note: Base R also has \code{cov2cor} which accomplishes the same task. This
#' function provides an alternative implementation within the WRS package.
#'
#' @export
#' @family covariance estimation
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(300), ncol = 3)
#' C <- cov(x)
#' R <- cov2cor(C)  # Convert to correlation
#' all.equal(diag(R), rep(1, 3))  # Diagonal should be 1
#' }
cov2cor<-function(x){
#
# Convert a covariance matrix to a correlation matrix
#
p=ncol(x)
m=x
for(i in 1:p){
for(j in 1:p){
m[i,j]=m[i,j]/sqrt(x[i,i]*x[j,j])
}}
m
}

#' Compute Deterministic MCD Covariance Estimator
#'
#' Calculates robust location and scatter using the deterministic MCD algorithm,
#' which avoids random subsampling for reproducible results.
#'
#' @param x Numeric matrix (n × p) of observations
#'
#' @return A list with two components:
#'   \item{center}{p-dimensional robust location estimate}
#'   \item{cov}{p × p robust scatter matrix}
#'
#' @details
#' This function wraps \code{DetMCD} from the \pkg{DetMCD} package, which implements
#' a deterministic version of the MCD (Minimum Covariance Determinant) estimator.
#' Unlike the random-start MCD in \code{\link{covmcd}}, this version:
#' \itemize{
#'   \item Produces deterministic, reproducible results
#'   \item Doesn't require setting random seeds
#'   \item May be slower for large datasets
#' }
#'
#' Requires the \pkg{DetMCD} package to be installed.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{covmcd}} for standard (random-start) MCD,
#'   \code{\link{mcdcov}} for MCD covariance only
#'
#' @references
#' Hubert, M., Rousseeuw, P.J., and Verdonck, T. (2012). A deterministic algorithm
#' for robust location and scatter. Journal of Computational and Graphical
#' Statistics, 21, 618-637.
#'
#' @examples
#' \dontrun{
#' # Multivariate data with outliers
#' set.seed(123)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:5, ] <- matrix(5, nrow = 5, ncol = 3)
#'
#' # Deterministic MCD (reproducible without seed)
#' result <- DETMCD(x)
#' result$center
#' result$cov
#' }
DETMCD<-function(x){
#
# HUber et al. (2012) deterministic version of the MCD estimator
#
library(DetMCD)
a=DetMCD(x)
list(center=a$center,cov=a$cov)
}

#' Compute Normalized Winsorized Covariance Matrix
#'
#' Calculates winsorized covariance matrix rescaled to estimate the classical
#' covariance under normality when no trimming is used.
#'
#' @param x Numeric matrix or vector
#' @param y Optional numeric vector (if x is a vector)
#' @param tr Proportion of winsorization (default: 0.2)
#'
#' @return Normalized winsorized covariance matrix (or scalar if bivariate)
#'
#' @details
#' This function computes the winsorized covariance via \code{\link{wincor}} and
#' then rescales it to be consistent for the classical covariance under normality.
#' The rescaling factor accounts for the efficiency loss due to winsorization.
#'
#' When \code{tr = 0} (no winsorization), the function returns the classical
#' covariance.
#'
#' @export
#' @family covariance estimation
#' @seealso \code{\link{wincov}} for unscaled winsorized covariance,
#'   \code{\link{wincor}} for winsorized correlation
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100)
#'
#' # Standard winsorized covariance
#' wincov(cbind(x, y), tr = 0.2)[1, 2]
#'
#' # Normalized version
#' wincovN(x, y, tr = 0.2)
#' }
wincovN<-function(x,y=NULL,tr=0.2){
#
# Winsorized covariance rescaled to est cov under normality when there is no trimming
#
e=wincor(x,y,tr=tr)$cov
if(tr==0)cterm=1
else cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
e=e/cterm
e
}

