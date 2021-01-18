\name{CovEM}
\alias{CovEM}
\title{Gaussian MLE of mean and covariance}
\description{Computes the Gaussian MLE via EM-algorithm for missing data.}
\usage{CovEM(x, tol=0.001, maxiter=1000)}
\arguments{
  \item{x}{a matrix or data frame. May contain missing values, but cannot contain columns with completely missing entries.}
  \item{tol}{tolerance level for the maximum relative change of the estimates. Default is 0.001.}
  \item{maxiter}{maximum iteration for the EM algorithm. Default is 1000.}
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
\author{Mike Danilov, Andy Leung \email{andy.leung@stat.ubc.ca}}



