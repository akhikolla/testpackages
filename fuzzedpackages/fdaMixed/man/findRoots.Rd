\name{findRoots}
\alias{findRoots}
\title{Complex roots of quadratic polynomial}
\description{Find complex roots of polynomials in x that are quadratic
  polynomials in x^k}
\usage{
findRoots(coefs, k = 1)
}
\arguments{
  \item{coefs}{Coefficients \code{(c_0,c_k,c_2k)} of quadratic polynomial in
    x^k. Also accepts matrix input (J,3).} 
  \item{k}{Order of x^k}
}
\details{
  It is assumed that \code{c_2k} is non-zero, and that at least one of
  \code{c_0} and \code{c_k} are non-zero (otherwise, we have a double
  root, which is not treated by \code{\link{fdaLm}} in the present
  implementation). An error is issued if these assumptions are violated.
} 
\value{ 
  A list with components
  \item{left}{The k roots with left most real components}
  \item{right}{The k roots with right most real components}
}
\note{
  This function is intended for internal usage in \code{\link{fdaLm}} to
  find eigenvalues. If a robust and stable method of finding all the
  complex roots is a polynomial were available, then this could be used
  in \code{\link{fdaLm}} instead enhancing the scope of this function.
}
\references{
  Solved using Section 5.6 in Press et al, "Numerical Recipies in C",
  second edition.
}
\author{Bo Markussen <bomar@math.ku.dk>}
\examples{
findRoots(c(-1,0,1),1)
findRoots(c(1,-1,1),2)
}
\keyword{math}
