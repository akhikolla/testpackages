\name{summary.clifford}
\alias{summary.clifford}
\alias{print.summary.clifford}
\alias{first_n_last}
\title{Summary methods for clifford objects}
\description{
Summary method for clifford objects, and a print method for summaries.
}
\usage{
\method{summary}{clifford}(object, ...)
\method{print}{summary.clifford}(x, ...)
first_n_last(x)
}
\arguments{
  \item{object,x}{Object of class clifford}
  \item{...}{Further arguments, currently ignored}
}
\details{
  Summary of a clifford object.   Note carefully that the
  \dQuote{typical terms} are implementation specific.  Function
  \code{first_n_last()} is a helper function.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{print}}}
\examples{
summary(rcliff())
}
\keyword{math}
