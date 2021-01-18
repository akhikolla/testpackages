\name{vcov--methods}
\docType{methods}
\alias{vcov}
\alias{vcov-methods}
\alias{vcov,IdtNDE-method}
\alias{vcov,IdtSNDE-method}
\alias{vcov,IdtNandSNDE-method}
\alias{vcov,IdtMxNDE-method}
\alias{vcov,IdtMxSNDE-method}

\title{Methods for function vcov in Package `MAINT.Data'}
\description{S4 methods for function vcov. As in the generic vcov S3 `stats' method, these methods extract variance-covariance estimates of parameter estimators, for the models fitted to Interval Data.}

\keyword{methods}
\keyword{vcov}
\keyword{Interval Data}

\usage{

\S4method{vcov}{IdtNDE}(object, selmodel=BestModel(object), \dots)
\S4method{vcov}{IdtSNDE}(object, selmodel=BestModel(object), \dots)
\S4method{vcov}{IdtNandSNDE}(object, selmodel=BestModel(object), \dots)
\S4method{vcov}{IdtMxNDE}(object, selmodel=BestModel(object), group=NULL, \dots)
\S4method{vcov}{IdtMxSNDE}(object, selmodel=BestModel(object), group=NULL, \dots)

}

\arguments{
  \item{object}{ An object representing a model fitted to interval data.}

  \item{selmodel}{Selected model from a list of candidate models saved in object.}

  \item{group}{The group for each the estimated parameter variance-covariance will be returned. If NULL (default),
\dQuote{vcov} will return a three-dimensional array with a matrix of the estimated covariances between the parameter estimates for each group at each level of the third dimension.
Note that this argument is only used in heterocedastic models, i.e. in the \code{\linkS4class{IdtMxNDE}}, \code{\linkS4class{IdtMxSNDE}} methods when the object slot \dQuote{Hmcdt} is set to to FALSE.}

  \item{\dots}{ Additional arguments for method functions.}
}

\value{For the \code{\linkS4class{IdtNDE}}, \code{\linkS4class{IdtSNDE}} and \code{\linkS4class{IdtNandSNDE}} methods or \code{\linkS4class{IdtMxNDE}}, \code{\linkS4class{IdtMxSNDE}} methods with slot \dQuote{Hmcdt} equal to TRUE: a matrix of the estimated covariances between the parameter estimates. For the \code{\linkS4class{IdtMxNDE}}, and \code{\linkS4class{IdtMxSNDE}} methods with slot \dQuote{Hmcdt} equal to FALSE: if argument \dQuote{group} is set to NULL, a three-dimensional array with a matrix of the estimated covariances between the parameter estimates for each group at each level of the third dimension. If argument \dQuote{group} is set to an integer, the matrix with the estimated covariances between the parameter estimates, for the group chosen.}

\seealso{
  \code{\link{stdEr}} 
}

