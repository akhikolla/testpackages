\name{TSGS-class}
\docType{class}
\alias{TSGS-class}
\alias{getFiltDat,TSGS-method}
\title{Two-Step Generalized S-Estimator for cell- and case-wise outliers}
\description{Class of Two-Step Generalized S-Estimator. It has the superclass of \code{GSE}. }
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("TSGS", ...)},
    but the best way of creating \code{TSGS} objects is a call to the function
    \code{TSGS} which serves as a constructor.
}
\section{Slots}{
    \describe{
	\item{\code{mu} }{ Estimated location. Can be accessed via \code{\link{getLocation}}. }
	\item{\code{S} }{ Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. } 
	\item{\code{xf} }{ Filtered data matrix from the first step of 2SGS. Can be accessed via \code{\link{getFiltDat}}. }
	}
}
\section{Extends}{ Class \code{"\linkS4class{GSE}"}, directly. }
\section{Methods}{ 
	In addition to the methods defined in the superclass "GSE", the following methods are also defined:
	  \describe{
 	\item{getFiltDat}{\code{signature(object = "TSGS")}: return the filtered data matrix.}
	}
}
\author{ Andy Leung \email{andy.leung@stat.ubc.ca} }
\seealso{
   \code{\link{TSGS}}, \code{\link{GSE}}, \code{\link{GSE-class}}
}
\keyword{classes}

