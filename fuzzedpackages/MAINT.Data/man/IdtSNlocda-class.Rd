\name{IdtSNlocda-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtSNlocda-class}
\alias{predict,IdtSNlocda-method}
\alias{show,IdtSNlocda-method}
\alias{CovCase,IdtSNlocda-method}

\title{Class "IdtSNlocda"}
\description{IdtSNlocda contains the results of Discriminant Analysis for the interval data, based on a location Skew-Normal model.}

\section{Slots}{
  \describe{
    \item{\code{prior}:}{Prior probabilities of class membership; if unspecified, the class proportions for the training set are used; if present, the probabilities should be specified in the order of the factor levels.}
    \item{\code{ksi}:}{Matrix with the direct location parameter ("ksi") estimates for each group.}
    \item{\code{eta}:}{Vector with the direct scaled skewness parameter ("eta") estimates.}
    \item{\code{scaling}:}{Matrix which transforms observations to discriminant functions, normalized so that the within groups scale-association matrix ("Omega") is spherical.}
    \item{\code{mu}:}{Matrix with the centred location parameter ("mu") estimates for each group.}
    \item{\code{gamma1}:}{Vector with the centred skewness parameter ("gamma1") estimates.}
    \item{\code{N}:}{Number of observations.}
    \item{\code{CovCase}:}{Configuration case of the variance-covariance matrix: Case 1 through Case 4}
  }
}

\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "IdtSNlocda")}: Classifies interval-valued observations in conjunction with snda. }
    \item{show}{\code{signature(object = "IdtSNlocda")}: show S4 method for the IDdtlda-class}
    \item{CovCase}{\code{signature(object = "IdtSNlocda")}: Returns the configuration case of the variance-covariance matrix}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Azzalini, A. and Dalla Valle, A. (1996), The multivariate skew-normal distribution. \emph{Biometrika} \bold{83}(4), 715--726.\cr

Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Duarte Silva, A.P. and Brito, P. (2015), Discriminant analysis of interval data: An assessment of parametric and distance-based   approaches. \emph{Journal of Classification} \bold{39}(3), 516--541.
}


\seealso{
  \code{\link{snda}},  \code{\link{MANOVA}},  
\code{\linkS4class{IData}}
}

\keyword{classes}
\keyword{interval data}

