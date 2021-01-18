\name{allcliff}
\alias{allcliff}
\title{Clifford object containing all possible terms}
\description{
The Clifford algebra on basis vectors \eqn{e_1,e_2,\ldots,
  e_n}{e1,e2,...,en} has \eqn{2^n} independent multivectors.  Function
  \code{allcliff()} generates a clifford object with a nonzero
  coefficient for each multivector.
}
\usage{allcliff(n)}
\arguments{\item{n}{Integer}}
\author{Robin K. S. Hankin}
\examples{
allcliff(6)

a <- allcliff(5)
a[] <- rcliff()*100

}