\name{IdtSNgenda-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtSNgenda-class}
\alias{predict,IdtSNgenda-method}
\alias{show,IdtSNgenda-method}
\alias{CovCase,IdtSNgenda-method}

\title{Class "IdtSNgenda"}
\description{IdtSNgenda contains the results of discriminant analysis for the interval data, based on a general Skew-Normal model.}

\section{Slots}{
  \describe{
    \item{\code{prior}:}{Prior probabilities of class membership; if unspecified, the class proportions for the training set are used; if present, the probabilities should be specified in the order of the factor levels.}
    \item{\code{ksi}:}{Matrix with the direct location parameter ("ksi") estimates for each group.}
    \item{\code{eta}:}{Matrix with the direct scaled sekwness parameter ("eta") estimates for each group.}
    \item{\code{scaling}:}{For each group g, scaling[,,g] is a matrix which transforms interval-valued observations so that in each group the scale-association matrix ("Omega") is spherical.}
    \item{\code{mu}:}{Matrix with the centred location parameter ("mu") estimates for each group.}
    \item{\code{gamma1}:}{Matrix with the centred sekwness parameter ("gamma1") estimates for each group.}
    \item{\code{ldet}:}{Vector of half log determinants of the dispersion matrix.}
    \item{\code{lev}:}{Levels of the grouping factor.}
    \item{\code{CovCase}:}{Configuration case of the variance-covariance matrix: Case 1 through Case 4}
  }
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "IdtSNgenda")}: Classifies interval-valued observations in conjunction with snda.  }
    \item{show}{\code{signature(object = "IdtSNgenda")}: show S4 method for the IdtSNgenda-class}
    \item{CovCase}{\code{signature(object = "IdtSNgenda")}: Returns the configuration case of the variance-covariance matrix}
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
  \code{\link{MANOVA}}, \code{\link{snda}}, \code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}
