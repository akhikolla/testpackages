\name{partial.mahalanobis}
\alias{partial.mahalanobis}
\title{Partial Square Mahalanobis Distance}
\description{
   Computes the partial square Mahalanobis distance for all observations in x. 
   Let \eqn{\mathbf{x} = (x_{i1},...,x_{ip})'} be a p-dimensional random vector and 
   \eqn{\mathbf{u} = (u_{i1},...,u_{ip})'} be a p-dimensional vectors of zeros and ones indicating which entry is missing:
   0 as missing and 1 as observed. Then partial mahalanobis distance is given by:
   
   \eqn{ d(\mathbf{x, u, m, \Sigma}) = (\mathbf{x^{(u)}} - \mathbf{m^{(u)}})' (\mathbf{\Sigma^{(u)}})^{-1}(\mathbf{x^{(u)}} - \mathbf{m^{(u)}})}.
   
   With no missing data, this function is equivalent to mahalanobis distance.
}
\usage{partial.mahalanobis(x, mu, Sigma)}
\arguments{
  \item{x}{a matrix or data frame. May contain missing values, but cannot contain columns with completely missing entries.}
  \item{mu}{location estimate}
  \item{Sigma}{scatter estimate. Must be positive definite}
}
\value{
An S4 object of class \code{\link{CovRobMiss-class}}.
The output S4 object contains the following slots:

  \tabular{ll}{
	\code{mu} \tab Estimated location. Can be accessed via \code{\link{getLocation}}. \cr
	\code{S} \tab Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. \cr 
	\code{pmd} \tab Squared partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. \cr
	\code{pmd.adj} \tab Adjusted squared partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. \cr
	\code{pu} \tab Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. \cr
	\code{call} \tab Object of class \code{"language"}. Not meant to be accessed. \cr
	\code{x} \tab Input data matrix. Not meant to be accessed. \cr
	\code{p} \tab Column dimension of input data matrix. Not meant to be accessed. \cr
	\code{estimator}\tab Character string of the name of the estimator used. Not meant to be accessed. \cr
	}
}
\examples{
\dontrun{
## suppose we would like to compute pmd for an MLE
x <- matrix(rnorm(1000),100,10)
U <- matrix(rbinom(1000,1,0.1),100,10)
x <- x * ifelse(U==1,NA,1)
## compute MLE (i.e. EM in this case)
res <- CovEM(x)
## compute pmd
res.pmd <- partial.mahalanobis(x, mu=getLocation(res), S=getScatter(res))
summary(res.pmd)
plot(res.pmd, which="index")
}
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}}



