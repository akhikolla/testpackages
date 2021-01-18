\name{stdEr--methods}
\docType{methods}
\alias{stdEr}
\alias{stdEr-methods}
\alias{stdEr,IdtNDE-method}
\alias{stdEr,IdtSNDE-method}
\alias{stdEr,IdtNandSNDE-method}

\title{Methods for function stdEr in Package `MAINT.Data'}
\description{S4 methods for function stdEr. As in the generic stdEr S3 `miscTools' method, these methods extract standard errors of the parameter estimates, for the models fitted to Interval Data.}

\keyword{methods}
\keyword{stdEr}
\keyword{Interval Data}

\usage{

\S4method{stdEr}{IdtNDE}(x, selmodel=BestModel(x), \dots)
\S4method{stdEr}{IdtSNDE}(x, selmodel=BestModel(x), \dots)
\S4method{stdEr}{IdtNandSNDE}(x, selmodel=BestModel(x), \dots)

}

\arguments{
  \item{x}{ An object representing a model fitted to interval data.}

  \item{selmodel}{Selected model from a list of candidate models saved in object x.}

  \item{\dots}{ Additional arguments for method functions.}
}

\value{ A vector of the estimated standard deviations of the parameter estimators.}

\seealso{
  \code{\link{vcov}} 
}

