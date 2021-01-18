\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-method}
\alias{plot,CovRobMiss-method}
\alias{plot,CovRobMiss,missing-method}
\title{Plot methods for objects of class 'CovRobMiss'}
\description{
    Plot methods for objects of class 'CovRobMiss'. The following plots are available:
	
	- chi-square qqplot for adjusted square partial Mahalanobis distances
	
	- index plot for adjusted square partial Mahalanobis distances

	- distance-distance plot comparing the adjusted distances based on classical MLE and robust estimators
	
	Cases with completely missing data will be dropped out. Outliers are identified using some pre-specific 
	cutoff value, for instance 99\% quantile of chi-square with p degrees of freedom, where p is the column
	dimension of the data. Identified outliers can also be retrieved using \code{\link{getOutliers}} with
	an optional argument of \code{cutoff}, ranged from 0 to 1.
}
\usage{
	\S4method{plot}{CovRobMiss}(x, which = c("all","distance","qqchi2", "dd"), 
		which = c("all", "distance", "qqchisq", "dd"),
		ask = (which=="all" && dev.interactive(TRUE)),
		cutoff = 0.99, xlog10 = FALSE, ylog10 = FALSE)
} 
\arguments{
 \item{x}{an object of class \code{"CovRobMiss"}}
  \item{which}{ Which plot to show? Default is \code{which}="all". }
  \item{ask}{ logical; if 'TRUE', the user is \emph{ask}ed before each plot, see 'par(ask=.)'. 
	Default is \code{ask = which=="all" && dev.interactive()}.  }
  \item{cutoff}{ The quantile cutoff for the distances. Default is 0.99. }
  \item{xlog10}{ Base-10 logged x-axis? Default is \code{FALSE}. }
  \item{ylog10}{ Base-10 logged y-axis? Default is \code{FALSE}. }
}
\examples{
\dontrun{
data(boston)
res <- GSE(boston)

## plot all graphs
plot(res)

## plot individuals plots
plot(res, which="qqchisq")
plot(res, which="index")
plot(res, which="dd")

## control the coordinates, e.g. log10 transform the y-axis
plot(res, which="qqchisq", xlog10=TRUE, ylog10=TRUE)
plot(res, which="index", ylog10=TRUE)
plot(res, which="dd", xlog10=TRUE, ylog10=TRUE)
}}
\keyword{methods}

