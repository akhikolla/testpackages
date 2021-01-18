\name{rspray}
\alias{rspray}
\title{
Random spray objects
}
\description{
Creates random spray objects as quick-and-dirty examples of multivariate polynomials
}
\usage{
rspray(n, vals = 1, arity = 3, powers = 0:2)
}
\arguments{
  \item{n}{
    Number of distinct rows (maximum); repeated rows are merged
    (argument \code{addrepeats} is \code{TRUE})}
  \item{vals}{Values to use for coefficients}
  \item{arity}{Arity of the spray; the number of columns in the index matrix}
  \item{powers}{Set from which to sample the entries of the index matrix}
}
\value{
Returns a spray object
}
\author{
Robin K. S. Hankin
}
\note{
If the index matrix contains repeated rows, the returned spray object
will  contain fewer than \code{n} entries
}
\seealso{\code{\link{spray}}}
\examples{

rspray(10)

rspray(4)*rspray(3,rnorm(3))

rspray(3,arity=7,powers=-2:2)^3

rspray(1000,vals=rnorm(1000))

}

\keyword{symbolmath}
