\name{SummaryCovGSE-class}
\docType{class}
\alias{SummaryCovGSE-class}
\alias{show,SummaryCovGSE-method}

\title{Class "SummaryCovGSE" - displaying summary of "CovRobMiss" objects }
\description{ Displays summary information for \code{\link{CovRobMiss-class}} objects }
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("SummaryCovGSE", ...)}.
}
\section{Slots}{
    \describe{
        \item{\code{obj}:}{\code{\link{CovRobMiss-class}} object}
        \item{\code{evals}:}{Eigenvalues and eigenvectors of the covariance or correlation matrix}
    }
}
\section{Methods}{ 
  \describe{
    \item{show}{\code{signature(object = "SummaryCovGSE")}: display the object }
  }
}
\author{ Andy Leung \email{andy.leung@stat.ubc.ca} }
\keyword{classes}

