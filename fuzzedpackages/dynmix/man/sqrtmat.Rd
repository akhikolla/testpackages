\name{sqrtmat}
\alias{sqrtmat}

\title{Computes the Square Root of a Matrix.}

\description{
This function computes the square root of a matrix.
}

\usage{
sqrtmat(A)
}

\arguments{
\item{A}{symmetric positive-definite {\code{\link[base]{matrix}}}}
}

\value{{\code{\link[base]{matrix}}} \eqn{B} such that \eqn{BB'=A}
}

\examples{
A <- matrix(c(5,1,1,3),2,2)
B <- sqrtmat(A)
B \%*\% t(B)
A
}

\references{
\url{https://en.wikipedia.org/wiki/Square_root_of_a_matrix}
}
