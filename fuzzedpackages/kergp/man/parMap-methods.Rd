\name{parMap-methods}
\docType{methods}
\alias{parMap,covTS-method}
\title{
  Map the Parameters of a Structure on the Inputs and Kernel Parameters
}
\description{
  
  Map the parameters of a structure on the inputs and kernel parameters.

}
\usage{
\S4method{parMap}{covTS}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    An object with class \code{"covTS"}.

  }
  \item{\dots}{

    Not used yet.

  }
}
%%\details{}
\value{
  
  A matrix with integer values. The rows correspond to the inputs of the
  object and the columns to the \eqn{1d} kernel parameters.
  The matrix element is the number of the corresponding official
  coefficient. The same parameter of the structure can be used for
  several inputs but not (yet) for several kernel parameters. So
  each row must have different integer elements, while the same
  element can be repeated within a column.

}

\note{

  This function is for internal use only.

}
%%\seealso{}
\examples{
myCov <- covTS(d = 3, kernel = "k1Gauss",
               dep = c(range = "input"), value = c(range = 1.1))
parMap(myCov)
}
%%\keyword{ }
%%\keyword{ }% __ONLY ONE__ keyword per line
