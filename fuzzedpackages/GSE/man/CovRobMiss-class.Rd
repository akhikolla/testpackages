\name{CovRobMiss-class}
\docType{class}
\alias{CovRobMiss-class}
\alias{show,CovRobMiss-method}
\alias{summary,CovRobMiss-method}
\alias{getDist,CovRobMiss-method}
\alias{getDistAdj,CovRobMiss-method}
\alias{getDim,CovRobMiss-method}
\alias{getLocation,CovRobMiss-method}
\alias{getScatter,CovRobMiss-method}
\alias{getMissing,CovRobMiss-method}
\alias{getOutliers,CovRobMiss-method}
\title{ Class "CovRobMiss" -- a superclass for the robust estimates of location and scatter for missing data }
\description{ The Superclass of all the objects output from the various robust estimators of location 
	and scatter for missing data, which includes Generalized S-estimator \code{\link{GSE}}, Extended Minimum Volumn 
	Ellipsoid \code{\link{emve}}, and Huberized Pairwise \code{\link{HuberPairwise}}. It can also be constructed
	using the code \code{\link{partial.mahalanobis}}. }
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovRobMiss", ...)},
    but the best way of creating \code{CovRobMiss} objects is a call to either of the folowing 
	functions:\code{GSE}, \code{emve}, \code{HuberPairwise}, and \code{partial.mahalanobis}, 
	which all serve as a constructor.
}
\section{Slots}{
	\describe{
	\item{\code{mu} }{ Estimated location. Can be accessed via \code{\link{getLocation}}. }
	\item{\code{S} }{ Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. } 
	\item{\code{pmd} }{ Square partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. }
	\item{\code{pmd.adj} }{ Adjusted square partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. }
	\item{\code{pu} }{ Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. }
	\item{\code{call} }{ Object of class \code{"language"}. Not meant to be accessed. }
	\item{\code{x} }{ Input data matrix. Not meant to be accessed. }
	\item{\code{p} }{ Column dimension of input data matrix. Not meant to be accessed. }
	\item{\code{estimator}}{ Character string of the name of the estimator used. Not meant to be accessed. }
	}
}
\section{Methods}{ 
  \describe{
    \item{show}{\code{signature(object = "CovRobMiss")}: display the object }
    \item{summary}{\code{signature(object = "CovRobMiss")}: calculate summary information }
	\item{plot}{\code{signature(object = "CovRobMiss", cutoff = "numeric")}: plot the object. See \code{\link{plot}}}
	\item{getDist}{\code{signature(object = "CovRobMiss")}: return the squared partial Mahalanobis distances }
	\item{getDistAdj}{\code{signature(object = "CovRobMiss")}: return the adjusted squared partial Mahalanobis distances }
	\item{getDim}{\code{signature(object = "CovRobMiss")}: return the dimension of observed entries for each case }
	\item{getLocation}{\code{signature(object = "CovRobMiss")}: return the estimated location vector}
	\item{getScatter}{\code{signature(object = "CovRobMiss", cutoff = "numeric")}: return the estimated scatter matrix}
	\item{getMissing}{\code{signature(object = "CovRobMiss")}: return the case number with completely missing data, if any}
	\item{getOutliers}{\code{signature(object = "CovRobMiss", cutoff = "numeric")}: return the case number(s) adjusted squared 
		distances above \code{(1 - cutoff)}th quantile of chi-square p-degrees of freedom.}	
	}
}
\author{ Andy Leung \email{andy.leung@stat.ubc.ca} }
\seealso{
   \code{\link{GSE}}, \code{\link{emve}}, \code{\link{HuberPairwise}}, \code{\link{partial.mahalanobis}}
}
\keyword{classes}

