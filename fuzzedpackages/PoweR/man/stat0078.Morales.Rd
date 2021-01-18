\name{stat0078.Morales}
\alias{stat0078.Morales}
\alias{stat0078}

\title{The Morales test for uniformity}

\description{The Morales test \code{D_{n,m}(phi_lambda)} for uniformity is used

- to compute its statistic and p-value by calling function \code{\link{statcompute}};

- to compute its quantiles by calling function \code{\link{compquant}} or \code{\link{many.crit}};

- to compute its power by calling function \code{\link{powcomp.fast}} or \code{\link{powcomp.easy}}.
}

\details{If \eqn{\lambda} and \code{m} are not specified they assume the default values of 0 and 2, respectively. 

There are 3 choices for value of \eqn{\lambda} : \eqn{\lambda} = 0, \eqn{\lambda} = -1, and \eqn{\lambda} != 0, != -1. 

Note that \code{n} is the sample size.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Morales, D., Pardo, L., Pardo, M. C. and Vajda, I. (2003), Limit laws for disparities of spacings, 
\emph{Journal of Nonparametric Statistics}, \bold{15}(3), 325--342.

M. A. Marhuenda, Y. Marhuenda, D. Morales, (2005), Uniformity tests under quantile categorization, 
\emph{Kybernetes}, \bold{34}(6), 888--901.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Uniformity.tests}} for other goodness-of-fit tests for uniformity.
}
\keyword{statistic, test, uniformity, Morales}
