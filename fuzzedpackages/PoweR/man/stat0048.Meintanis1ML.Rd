\name{stat0048.Meintanis1ML}
\alias{stat0048.Meintanis1ML}
\alias{stat0048}

\title{The 1st Meintanis test with maximum likelihood estimators for the Laplace distribution}

\description{The 1st Meintanis test \code{T_{n,a}^{(1)}} with maximum likelihood estimators test for the Laplace distribution is used

- to compute its statistic and p-value by calling function \code{\link{statcompute}};

- to compute its quantiles by calling function \code{\link{compquant}} or \code{\link{many.crit}};

- to compute its power by calling function \code{\link{powcomp.fast}} or \code{\link{powcomp.easy}}.
}

\details{If \code{a} is not specified it assumes the default value of 2.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Meintanis, S.G. (2004), A Class of Omnibus Tests for the Laplace Distribution Based on the Empirical Characteristic Function, 
\emph{Communications in Statistics - Theory and Methods}, \bold{33}(4), 925--948.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Laplace.tests}} for other goodness-of-fit tests for the Laplace distribution.
}
\keyword{statistic, test, Laplace, Meintanis}
