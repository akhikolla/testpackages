\name{CovRobMissSc-class}
\docType{class}
\alias{CovRobMissSc-class}
\alias{getScale,CovRobMissSc-method}
\title{ Class "CovRobMissSc" -- a subclass of "CovRobMiss" with scale estimate }
\description{ The Superclass of the \code{\link{GSE-class}} and \code{\link{emve-class}} objects. }
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovRobMissSc", ...)},
    but the best way of creating \code{CovRobMissSc} objects is a call to either of the folowing 
	functions:\code{GSE} or \code{emve}.
}
\section{Slots}{
    \describe{
	\item{\code{mu} }{ Estimated location. Can be accessed via \code{\link{getLocation}}. }
	\item{\code{S} }{ Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. } 
	\item{\code{sc} }{ Estimated M-scale (either GS-scale or MVE-scale). Can be accessed via \code{\link{getScale}}. } 	
	\item{\code{pmd} }{ Square partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. }
	\item{\code{pmd.adj} }{ Adjusted square partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. }
	\item{\code{pu} }{ Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. }
	\item{\code{call} }{ Object of class \code{"language"}. Not meant to be accessed. }
	\item{\code{x} }{ Input data matrix. Not meant to be accessed. }
	\item{\code{p} }{ Column dimension of input data matrix. Not meant to be accessed. }
	\item{\code{estimator}}{ Character string of the name of the estimator used. Not meant to be accessed. }
	}
}
\section{Extends}{ Class \code{"\linkS4class{CovRobMiss}"}, directly. }
\section{Methods}{ 
  \describe{
	In addition to methods inheritedfrom the class "CovRobMiss":
	\item{getScale}{\code{signature(object = "CovRobMissSc")}}: return the GS-scale or MVE-scale of the best candidate.  
	}
}
\author{ Andy Leung \email{andy.leung@stat.ubc.ca} }
\seealso{
   \code{\link{GSE}}, \code{\link{CovRobMiss-class}}
}
\keyword{classes}

