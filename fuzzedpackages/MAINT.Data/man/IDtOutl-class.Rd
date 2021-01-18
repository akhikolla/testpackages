\name{IdtOutl-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtOutl-class}
\alias{print,IdtOutl-method}
\alias{getMahaD2}
\alias{getMahaD2,IdtOutl-method}
\alias{geteta}
\alias{geteta,IdtOutl-method}
\alias{getRefDist}
\alias{getRefDist,IdtOutl-method}
\alias{getmultiCmpCor}
\alias{getmultiCmpCor,IdtOutl-method}

\title{Class IdtOutl}
\description{
  A description of interval-valued variable outliers found by the MAINT.Data function \code{\link{getIdtOutl}}.
}

\section{Slots}{
  \describe{
    \item{\code{outliers}:}{A vector of indices of the interval data units flaged as outliers.}
    \item{\code{MD2}:}{A vector of squared robust Mahalanobis distances for all interval data units.}
    \item{eta}{Nominal size of the null hypothesis that a given observation is not an outlier.}
    \item{RefDist}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq} and \dQuote{CerioliBetaF} respectivelly for the usual Chi-squared, and the Beta and F distributions proposed by Cerioli (2010).}    
    \item{multiCmpCor}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests was performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at the \sQuote{eta} nominal level. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n)).}
    \item{NObs}{Number of original observations in the original data set.}
    \item{p}{Number of total numerical variables (MidPoints and/or LogRanges) that may be responsible for the outliers.}
    \item{h}{Size of the subsets over which the trimmed likelihood was maximized when computing the robust Mahalanobis distances.})   
    \item{boolRewind}{A logical vector indicanting which of the data units belong to the final trimmed subsetused to compute the tle estimates.})   
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "IdtOutl")}: show S4 method for the IdtOutl-class.}
    \item{plot}{\code{signature(x = "IdtOutl")}: plot S4 methods for the IdtOutl-class.}
    \item{getMahaD2}{\code{signature(x = "IdtOutl")}: retrieves the vector of squared robust Mahalanobis distances for all data units.}
    \item{geteta}{\code{signature(x = "IdtOutl")}: retrieves the nominal size of the null hypothesis used to flag observations as outliers.}
    \item{getRefDist}{\code{signature(x = "IdtOutl")}: retrieves the assumed reference distributions used to find cutoffs defining the observations assumed as outliers.}
    \item{getmultiCmpCor}{\code{signature(x = "IdtOutl")}: retrieves the multicomparison correction used when flaging observations as outliers.}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156.\cr 

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr
}


\seealso{
\code{\link{getIdtOutl}}, \code{\link{fasttle}}, \code{\link{fulltle}} 
}

\keyword{class}
\keyword{outlier detection}
\keyword{trimmed maximum likelihood estimator}
\keyword{interval data}

