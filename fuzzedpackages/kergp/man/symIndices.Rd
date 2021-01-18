\name{symIndices}
\alias{symIndices}
\title{

  Vector of Indices Useful for Symmetric or Anti-Symmetric Matrices.

}
\usage{

symIndices(n, diag = FALSE)

}
\arguments{
  
  \item{n}{

    Size of a square matrix.

  }
  \item{diag}{

    Logical. When \code{FALSE} the diagonal is omitted in the lower and
    upper triangles.

}

}
\value{
  
A list containing the following integer vectors, each with length
\eqn{(n - 1) n / 2}. \item{i, j}{ Row and column indices for the lower
triangle to be used in a two-index style. } \item{kL}{ Indices for the
lower triangle, to be used in single-index style. The elements are
picked in column order. So if \code{X} is a square matrix with size
\code{n}, then \code{X[kL]} is the vector containing the elements of the
lower triangle of \code{X} taken in column order. } \item{kU}{ Indices
for the upper triangle, to be used in a single-index style. The
elements are picked in row order. So if \code{X} is a square matrix with
size \code{n}, then \code{X[kU]} is the vector containing the elements
of the upper triangle of \code{X} taken in row order. }

}
\description{
  
Vector of indices useful for symmetric or anti-symmetric matrices

}
\details{
  
This function is intended to provide computations which are faster than
\code{lower.tri} and \code{upper.tri}.

}
\examples{
n <- rpois(1, lambda = 10)
L <- symIndices(n)
X <- matrix(1L:(n * n), nrow = n)
max(abs(X[lower.tri(X, diag = FALSE)] - L$kL))
max(abs(t(X)[lower.tri(X, diag = FALSE)] - L$kU))
cbind(row = L$i, col = L$j)
}