\name{stat0067.Kuiper}
\alias{stat0067.Kuiper}
\alias{stat0067}

\title{The Kuiper test for uniformity}

\description{The Kuiper test \code{K_n} for uniformity is used

- to compute its statistic and p-value by calling function \code{\link{statcompute}};

- to compute its quantiles by calling function \code{\link{compquant}} or \code{\link{many.crit}};

- to compute its power by calling function \code{\link{powcomp.fast}} or \code{\link{powcomp.easy}}.
}

\details{Note that \code{n} is the sample size.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Brunk, H. D. (1962), On the range of the difference between hypothetical distribution function and Pyke's modified empirical distribution function, 
\emph{Annals of Mathematical Statistics}, \bold{33}, 525--532.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Uniformity.tests}} for other goodness-of-fit tests for uniformity.
}
\keyword{statistic, test, uniformity, Kuiper}
