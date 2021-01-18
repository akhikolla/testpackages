\name{Idtlda-class}
\Rdversion{1.1}
\docType{class}
\alias{Idtlda-class}
\alias{predict,Idtlda-method}
\alias{show,Idtlda-method}
\alias{CovCase,Idtlda-method}

\title{Class "Idtlda"}
\description{Idtlda contains the results of Linear Discriminant Analysis for the interval data
}

\section{Slots}{
  \describe{
    \item{\code{prior}:}{Prior probabilities of class membership; if unspecified, the class proportions for the training set are used; if present, the probabilities should be specified in the order of the factor levels.}
    \item{\code{means}:}{Matrix with the mean vectors for each group}
    \item{\code{scaling}:}{Matrix which transforms observations to discriminant functions, normalized so that within groups covariance matrix is spherical.}
    \item{\code{N}:}{Number of observations}
    \item{\code{CovCase}:}{Configuration case of the variance-covariance matrix: Case 1 through Case 4}
  }
}

\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "Idtlda")}: Classifies interval-valued observations in conjunction with lda. }
    \item{show}{\code{signature(object = "Idtlda")}: show S4 method for the IDdtlda-class}
    \item{CovCase}{\code{signature(object = "Idtlda")}: Returns the configuration case of the variance-covariance matrix}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Duarte Silva, A.P. and Brito, P. (2015), Discriminant analysis of interval data: An assessment of parametric and distance-based   approaches. \emph{Journal of Classification} \bold{39}(3), 516--541.
}

\seealso{
  \code{\link{qda}},  \code{\link{MANOVA}}, \code{\link{Roblda}},  \code{\link{Robqda}}, \code{\link{snda}}, \code{\linkS4class{IData}}
}

\keyword{classes}
\keyword{interval data}

