\name{MRC.asymptotic}
\docType{data}
\alias{MRC.asymptotic}

\title{"Asymptotic" values of the MRC statistic (all intervals)}

\description{
Simulated values of the MRC statistic with \code{penalty="sqrt"} based on all interval lengths computed from Gaussian white noise sequences of ("almost infinite") length 5,000.

\bold{Deprecation warning:} This data set is needed for \code{\link{smuceR}} and may be removed when this function will be removed.
}

\usage{MRC.asymptotic}

\format{A \code{\link{numeric}} vector containing 10,000 sorted values.}

\examples{
# "asymptotic" threshold value for 95% confidence
quantile(stepR::MRC.asymptotic, .95)
}

\keyword{datasets}
