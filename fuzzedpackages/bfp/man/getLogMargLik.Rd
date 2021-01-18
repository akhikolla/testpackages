\name{getLogMargLik}
\alias{getLogMargLik}

\title{Extract log marginal likelihood from a model.}
\description{
  Extract log marginal likelihood from a model, saved as the first
  element of a \code{\link{BayesMfp}} object.
}
\usage{
getLogMargLik(x, design=getDesignMatrix(x), nObs = nrow(design), dim = ncol(design))
}

\arguments{
  \item{x}{valid \code{\link{BayesMfp}}-Object of length 1 (otherwise only first
    element recognized)} 
  \item{design}{(centered) design matrix}
  \item{nObs}{number of observations}
  \item{dim}{number of design matrix columns}
}
\details{
  This function interfaces the C++ function \code{logMargLik}, and can
  be used to compute the marginal likelihood of a model not saved in the
  model list. But be careful to adjust the saved R^2 of the model, too,
  and not only the powers! Therefore this function is internal only...
  and is used e.g. in \code{\link{transformMfp}}.
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{getLogPrior}}}

\keyword{regression}
\keyword{internal}

