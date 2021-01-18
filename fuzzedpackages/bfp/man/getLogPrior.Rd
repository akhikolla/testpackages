\name{getLogPrior}
\alias{getLogPrior}

\title{Extract the log prior probability from a model.}
\description{
  Extract the log prior probability from the first model saved in a
  \code{\link{BayesMfp}} object.
}
\usage{
getLogPrior(x)
}

\arguments{
  \item{x}{a valid BayesMfp-Object of length 1 (otherwise only first
  element recognized)} 
}
\details{
  Note that the returned value might be only defined up to a constant,
  because some models might not be identifiable or have a dimension too
  large, and are thus implicitly assigned a zero prior probabibility.
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{getLogMargLik}}}
\keyword{regression}
\keyword{internal}

