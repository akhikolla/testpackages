\name{stat0072.ReadCressie}
\alias{stat0072.ReadCressie}
\alias{stat0072}

\title{The Read-Cressie test for uniformity}

\description{The Read-Cressie test \code{2nI^{lambda}} for uniformity is used

- to compute its statistic and p-value by calling function \code{\link{statcompute}};

- to compute its quantiles by calling function \code{\link{compquant}} or \code{\link{many.crit}};

- to compute its power by calling function \code{\link{powcomp.fast}} or \code{\link{powcomp.easy}}.
}

\details{If \eqn{\lambda} is not specified it assumes the default value of 1. Note that \code{n} is the sample size.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Read, Timothy R. C. and Cressie, Noel A. C. (1988), \emph{Goodness-of-fit statistics for discrete multivariate data}, 
Springer Series in Statistics. Springer-Verlag, New York.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Uniformity.tests}} for other goodness-of-fit tests for uniformity.
}
\keyword{statistic, test, uniformity, Read-Cressie}
