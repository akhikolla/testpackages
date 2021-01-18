\name{mean--methods}
\docType{methods}
\alias{mean}
\alias{mean-methods}
\alias{mean,IdtNDE-method}
\alias{mean,IdtSNDE-method}
\alias{mean,IdtNandSNDE-method}
\alias{mean,IdtMxNDE-method}
\alias{mean,IdtMxSNDE-method}

\title{Methods for function mean in Package `MAINT.Data'}
\description{S4 methods for function mean. These methods extract estimates of mean vectors for the models fitted to Interval Data.}

\keyword{methods}
\keyword{mean}
\keyword{Interval Data}

\usage{

\S4method{mean}{IdtNDE}(x)
\S4method{mean}{IdtSNDE}(x)
\S4method{mean}{IdtNandSNDE}(x)
\S4method{mean}{IdtMxNDE}(x)
\S4method{mean}{IdtMxSNDE}(x)

}

\arguments{
  \item{x}{ An object representing a model fitted to interval data.}
}

\value{For the \code{\linkS4class{IdtNDE}}, \code{\linkS4class{IdtSNDE}} and \code{\linkS4class{IdtNandSNDE}} methods or \code{\linkS4class{IdtMxNDE}}, \code{\linkS4class{IdtMxSNDE}} methods with slot \dQuote{Hmcdt} equal to TRUE: a matrix with the estimated correlations. \cr For the \code{\linkS4class{IdtMxNDE}}, and \code{\linkS4class{IdtMxSNDE}} methods with slot \dQuote{Hmcdt} equal to FALSE: a three-dimensional array with a matrix with the estimated correlations for each group at each level of the third dimension.}

\seealso{
  \code{\link{sd}} \code{\link{var}} \code{\link{cor}} 
}

