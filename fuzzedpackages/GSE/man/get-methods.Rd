\name{get-methods}
\docType{methods}
\alias{getDist}
\alias{getDist-methods} 
\alias{getDistAdj}
\alias{getDistAdj-methods} 
\alias{getDim}
\alias{getDim-methods} 
\alias{getLocation}
\alias{getLocation-methods} 
\alias{getScatter}
\alias{getScatter-methods} 
\alias{getMissing}
\alias{getMissing-methods} 
\alias{getOutliers}
\alias{getOutliers-methods} 
\alias{getScale}
\alias{getScale-methods} 
\alias{getFiltDat}
\alias{getFiltDat-methods}
\title{Accessor methods to the essential slots of classes CovRobMiss, TSGS, GSE, emve, and HuberPairwise}
\description{
    Accessor methods to the slots of objects of classes \code{CovRobMiss}, \code{TSGS}, \code{GSE}, \code{emve}, and \code{HuberPairwise}
}
\usage{
getLocation(object)
getScatter(object)
getDist(object)
getDistAdj(object)
getDim(object)
getMissing(object)
getOutliers(object, cutoff)
getScale(obj)
getFiltDat(object)
} 
\arguments{
\item{obj,object}{an object of any of the following classes
  \code{\link{CovRobMiss-class}}, \code{\link{GSE-class}},
  \code{\link{emve-class}}, and \code{\link{HuberPairwise-class}}. For
  function \code{getScale} the package defines a method for objects of class \code{\link{GSE-class}}
 objects are allowed. 
}
    \item{cutoff}{optional argument for \code{getOutliers} - 
	quantiles of chi-square to be used as a threshold for outliers detection, defaults to 0.99}
}
\details{
  \describe{
  	\item{getLocation}{\code{signature(object = "CovRobMiss")}: return the estimated location vector}
	\item{getScatter}{\code{signature(object = "CovRobMiss", cutoff = "numeric")}: return the estimated scatter matrix}
	\item{getDist}{\code{signature(object = "CovRobMiss")}: return the squared partial Mahalanobis distances }
	\item{getDistAdj}{\code{signature(object = "CovRobMiss")}: return the adjusted squared partial Mahalanobis distances }
	\item{getDim}{\code{signature(object = "CovRobMiss")}: return the dimension of observed entries for each case }
	\item{getMissing}{\code{signature(object = "CovRobMiss")}: return the case number with completely missing data, if any}
	\item{getOutliers}{\code{signature(object = "CovRobMiss", cutoff = "numeric")}: return the case number(s) adjusted squared 
		distances above \code{(1 - cutoff)}th quantile of chi-square p-degrees of freedom.}	
 	\item{getScale}{\code{signature(object = "CovRobMissSc")}: return either the estimated generalized S-scale or MVE-scale. See \code{\link{GSE}} and \code{\link{emve}} for details.}
	\item{getFiltDat}{\code{signature(object = "TSGS")}: return filtered data matrix from the first step of 2SGS.}
	}
}
\examples{
\dontrun{
data(boston)
res <- GSE(boston)

## extract estimated location
getLocation(res)

## extract estimated scatter
getScatter(res)

## extract estimated adjusted distances
getDistAdj(res)

## extract outliers
getOutliers(res)
}}
\keyword{get}
\keyword{methods}

