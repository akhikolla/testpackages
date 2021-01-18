\name{constant}
\alias{constant}
\alias{constant<-}
\title{
  Get or set the constant term of a spray object
}
\description{
The constant term of a spray object is the coefficient corresponding to
an index of all zeros.  These functions get or set the constant of a
spray object.
}
\usage{
constant(S,drop=FALSE)
constant(S) <- value
}
\arguments{
  \item{S}{Object of class spray}
  \item{value}{Numeric value to set the constant coefficient to}
  \item{drop}{Boolean, with default \code{FALSE} meaning to return a
    spray object and \code{TRUE} meaning to return a numeric value}
}
\value{
  In function \code{constant()}, return the coefficient, or a constant
  multivariate polynomial, depending on the value of \code{drop}.
}
\author{
Robin K. S. Hankin
}
\note{
The behaviour of \code{drop} (sort of) matches that of the spray extractor method.
}
\seealso{
\code{\link{Extract}}
}
\examples{
S <- spray(partitions::blockparts(rep(2,4),3,TRUE))

constant(S)

constant(S) <- 33

S

}

\keyword{symbolmath}

