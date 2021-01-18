\name{MRC.1000}
\docType{data}
\alias{MRC.1000}

\title{Values of the MRC statistic for 1,000 observations (all intervals)}

\description{
Simulated values of the MRC statistic with \code{penalty="sqrt"} based on all interval lengths computed from Gaussian white noise sequences of length 1,000.

\bold{Deprecation warning:} This data set is needed for \code{\link{smuceR}} and may be removed when this function will be removed.
}

\usage{MRC.1000}

\format{A \code{\link{numeric}} vector containing 10,000 sorted values.}

\examples{
# threshold value for 95% confidence
quantile(stepR::MRC.1000, .95)
}

\keyword{datasets}
