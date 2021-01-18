\name{HuberPairwise-class}
\docType{class}
\alias{HuberPairwise-class}
\title{ Quadrant Covariance and Huberized Pairwise Scatter }
\description{ Class of Quadrant Covariance and Huberized Pairwise Scatter. It has the superclass of \code{CovRobMiss}. }
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("HuberPairwise", ...)},
    but the best way of creating \code{HuberPairwise} objects is a call to the function
    \code{HuberPairwise} which serves as a constructor.
}
\section{Slots}{
	\describe{
	\item{\code{mu}}{Estimated location. Can be accessed via \code{\link{getLocation}}. }
	\item{\code{S}}{Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. } 
	\item{\code{pmd}}{Squared partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. }
	\item{\code{pmd.adj}}{Adjusted squared partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. }
	\item{\code{pu}}{Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. }
	\item{\code{R}}{Estimated correlation matrix. Not meant to be accessed.}
	\item{\code{call}}{Object of class \code{"language"}. Not meant to be accessed. }
	\item{\code{x}}{Input data matrix. Not meant to be accessed. }
	\item{\code{p}}{Column dimension of input data matrix. Not meant to be accessed. }
	\item{\code{estimator}}{Character string of the name of the estimator used. Not meant to be accessed. }
	}	
}
\section{Extends}{ Class \code{"\linkS4class{CovRobMiss}"}, directly. }
\section{Methods}{ No methods defined with class "HuberPairwise" in the signature. }
\author{ Andy Leung \email{andy.leung@stat.ubc.ca} }
\seealso{
   \code{\link{HuberPairwise}}, \code{\link{CovRobMiss-class}}
}
\keyword{classes}


