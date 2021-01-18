\name{stat0079.Pardo}
\alias{stat0079.Pardo}
\alias{stat0079}

\title{The Pardo test for uniformity}

\description{The Pardo test \code{E_{m,n}} for uniformity is used

- to compute its statistic and p-value by calling function \code{\link{statcompute}};

- to compute its quantiles by calling function \code{\link{compquant}} or \code{\link{many.crit}};

- to compute its power by calling function \code{\link{powcomp.fast}} or \code{\link{powcomp.easy}}.
}

\details{If \code{m} is not specified it assumes the default value of 2. Note that \eqn{m < (n/2)} where \code{n} is the sample size.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Pardo, M. C. (2003), A test for uniformity based on informational energy, 
\emph{Statistical Papers}, \bold{44}, 521--534.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Uniformity.tests}} for other goodness-of-fit tests for uniformity.
}
\keyword{statistic, test, uniformity, Pardo}
