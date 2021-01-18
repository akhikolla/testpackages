\name{MRC.asymptotic.dyadic}
\docType{data}
\alias{MRC.asymptotic.dyadic}

\title{"Asymptotic" values of the MRC statistic (dyadic intervals)}

\description{
Simulated values of the MRC statistic with \code{penalty="sqrt"} based on dyadic interval lengths computed from Gaussian white noise sequences of ("almost infinite") length 100,000.

\bold{Deprecation warning:} This data set is needed for \code{\link{smuceR}} and may be removed when this function will be removed.
}

\usage{MRC.asymptotic.dyadic}

\format{A \code{\link{numeric}} vector containing 10,000 sorted values.}

\examples{
# "asymptotic" threshold value for 95% confidence
quantile(stepR::MRC.asymptotic.dyadic, .95)
}

\keyword{datasets}
