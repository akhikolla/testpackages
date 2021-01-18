\name{asum}
\alias{asum}
\alias{asum.spray}
\alias{asum_inverted}
\alias{process_dimensions}
\title{
Sum over dimension margins
}
\description{
Sum over specified dimension margins.
}
\usage{
\method{asum}{spray}(S, dims, drop=TRUE, ...)
asum_inverted(S, dims)
process_dimensions(S,dims)
}
\arguments{
  \item{S}{spray object}
  \item{dims}{Vector of strictly positive integers corresponding to
    dimensions to be summed over
}
  \item{drop}{Boolean, with default \code{TRUE} meaning to drop the
    summed dimensions, and \code{FALSE} meaning to retain them.
}
  \item{\dots}{Further arguments, currently ignored}
}
\details{

Function \code{asum.spray()} is the method for \code{asum()}.  This
takes a spray, and a vector of integers corresponding to dimensions to
be summed over.

Function \code{asum_inverted()} is the same, but takes a vector of
integers corresponding to dimensions not to sum over.  This function is
here because  there is nice \code{C++} idiom for it.

Function \code{process_dimensions()} ensures that the \code{dims}
argument is consistent with the spray \code{S} and returns a cleaned
version thereof.

}
\value{
Returns a spray object.
}
\author{
Robin K. S. Hankin
}
\examples{

S <- spray(matrix(sample(0:2,60,replace=TRUE),ncol=3),addrepeats=TRUE)
S

asum(S,1)
asum(S,1:2)

asum(S,1:2,drop=FALSE)


asum(S,c(1,3)) == asum_inverted(S,2)


}


\keyword{symbolmath}
