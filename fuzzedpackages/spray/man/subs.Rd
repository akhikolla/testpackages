\name{subs}
\alias{subs}
\alias{substitute}
\title{
Substitute values  into a spray object
}
\description{
  Substitute values  into a spray object, interpreted as a multivariate polynomial
}
\usage{
subs(S, dims, x)
}
\arguments{
  \item{S}{spray object}
  \item{dims}{Integer or logical vector with entries corresponding to
    the dimensions to be substituted}
  \item{x}{Numeric vector of values to be substituted}
}
\author{
Robin K. S. Hankin
}
\note{
It is much easier if argument \code{dims} is sorted into increasing
order.  If not, caveat emptor!
}
\seealso{\code{\link{process_dimensions}}}
\examples{

S <- spray(matrix(sample(0:3,60,replace=TRUE),nrow=12))

subs(S,c(2,5),1:2)

subs(homog(3,3),1,3)
}
\keyword{symbolmath}

