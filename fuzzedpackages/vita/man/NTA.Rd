\name{NTA}
\alias{NTA}
\alias{NTA.default}
\alias{print.NTA}
\title{
Novel testing approach
}
\description{
Calculates the p-values for each permutation variable importance measure, based on the empirical null distribution
from non-positive importance values as described in Janitza et al. (2015).
}
\usage{
\method{NTA}{default}(PerVarImp)
\method{print}{NTA}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PerVarImp}{ permutation variable importance measures in a vector.}
  \item{x}{for the print method, an \code{NTA} object }
  \item{\dots}{optional parameters for \code{\link[base]{print}} }
}
\details{
The observed non-positive permutation variable importance values are used to approximate the distribution of
variable importance for non-relevant variables. The null distribution Fn0 is computed by mirroring the
non-positive variable importance values on the y-axis. Given the approximated null importance distribution,
the p-value is the probability of observing the \code{original PerVarImp} or a larger value. This testing
approach is suitable for data with large number of variables without any effect.

\code{PerVarImp} should be computed based on the hold-out permutation variable importance measures. If using
standard variable importance measures the results may be biased.

This function has not been tested for regression tasks so far, so this routine is meant for the expert user
only and its current state is rather experimental.
}
\value{
\item{PerVarImp }{ the orginal permutation variable importance measures. }
\item{M}{The non-positive variable importance values with the mirrored values on the y-axis.}
\item{pvalue }{the p-value is the probability of observing the \code{orginal PerVarImp} or a
               larger value, given the approximated null importance distribution.}
}
\references{
Janitza S, Celik E, Boulesteix A-L, (2015),\emph{ A computationally fast variable importance test for random forest for high dimensional data},Technical Report 185, University of Munich, <\url{http://nbn-resolving.de/urn/resolver.pl?urn=nbn:de:bvb:19-epub-25587-4}>
}

\seealso{
\code{\link[vita]{CVPVI}},\code{\link[randomForest]{importance}}, \code{\link[randomForest]{randomForest}}
}
\examples{
##############################
#      Classification        #
##############################
## Simulating data
X = replicate(100,rnorm(200))
X= data.frame( X) #"X" can also be a matrix
z  = with(X,2*X1 + 3*X2 + 2*X3 + 1*X4 -
            2*X5 - 2*X6 - 2*X7 + 1*X8 )
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = as.factor(rbinom(200,1,pr))
##################################################################
# cross-validated permutation variable importance

cv_vi = CVPVI(X,y,k = 2,mtry = 3,ntree = 500,ncores = 2)
##################################################################
#compare them with the original permutation variable importance
library("randomForest")
cl.rf = randomForest(X,y,mtry = 3,ntree = 500, importance = TRUE)
##################################################################
# Novel Test approach
cv_p = NTA(cv_vi$cv_varim)
summary(cv_p,pless = 0.1)
pvi_p = NTA(importance(cl.rf, type=1, scale=FALSE))
summary(pvi_p)


###############################
#      Regression             #
###############################
##################################################################
## Simulating data:
X = replicate(100,rnorm(200))
X = data.frame( X) #"X" can also be a matrix
y = with(X,2*X1 + 2*X2 + 2*X3 + 1*X4 - 2*X5 - 2*X6 - 1*X7 + 2*X8 )

##################################################################
# cross-validated permutation variable importance
cv_vi = CVPVI(X,y,k = 2,mtry = 3,ntree = 500,ncores = 2)
##################################################################
#compare them with the original permutation variable importance
reg.rf = randomForest(X,y,mtry = 3,ntree = 500, importance = TRUE)
##################################################################
# Novel Test approach (not tested for regression so far!)
cv_p = NTA(cv_vi$cv_varim)
summary(cv_p,pless = 0.1)
pvi_p = NTA(importance(reg.rf, type=1, scale=FALSE))
summary(pvi_p)
}

