\name{dual}
\alias{dual}
\title{The dual of a clifford object}
\description{The dual of a clifford object \eqn{C}, written \eqn{C^*}}
\usage{dual(C, n)}
\arguments{
  \item{C}{Clifford object}
  \item{n}{Dimensionality of underlying vector space}
}
\details{
The dual of clifford object \eqn{C} is  \eqn{CI^{-1}}{CI^-1} where
\eqn{I} is the pseudoscalar.

The dual is sensitive to the signature.  Note that applying the dual
operation four times successively will return

}
\references{
C. Perwass. \dQuote{Geometric algebra with applications in
  engineering}.   Springer, 2009.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{pseudoscalar}}}
\examples{
a <- rcliff()
dual(dual(dual(dual(a,8),8),8),8) == a # should be TRUE
}