\name{getFpTransforms}
\alias{'\%bt\%'}
\alias{getFpTransforms}
\alias{getTransformName}

\title{Internal functions to handle FP transforms.}
\description{
Transform a variable according to FP transformation formula and attach
proper names to the resulting design matrix. The binary function
\code{bt} is the Box-Tidwell transform.
}
\usage{
'\%bt\%'(x, pow)
getFpTransforms(vec, powers, center=TRUE)
getTransformName(name, pow)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector}
  \item{pow}{one single power (not vectorized in this argument)}
  \item{vec}{positive (== already shifted and scaled) column vector (!) with
    proper colname} 
  \item{powers}{power vector with at least one element}
  \item{center}{should the columns be centered around zero? (default)}
  \item{name}{name of the covariate to be transformed with \code{bt}}
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{getDesignMatrix}}}
\examples{
bfp:::getFpTransforms(cbind(x=1:10), powers=c(-2, 0, 1))
}

\keyword{utilities}
\keyword{internal}

