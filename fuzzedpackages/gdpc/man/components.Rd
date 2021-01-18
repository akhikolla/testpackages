\name{components}
\alias{components}
\title{Generic Function for Getting Components From an Object}
\description{
  Generic function for getting components from an object.
}


\usage{
  components(object, which_comp)
}
\arguments{
  \item{object}{An object. Currently there is a method for objects of class \code{gdpcs}.}
  \item{which_comp}{Numeric vector indicating which components to get. Default is 1.}
}


\value{ A matrix whose columns are the desired components.}
\author{
  Daniel Peña, Ezequiel Smucler, Victor Yohai
}


\keyword{ts}