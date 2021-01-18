\name{ldlt}
\alias{ldlt}

\title{Computes LDL' Matrix Decomposition.}

\description{
This function decomposes matrix \eqn{V} into \eqn{V=LDL'}, where \eqn{L} is a lower triangular matrix with unit diagonal and \eqn{D} is a diagonal matrix with non-negative terms.
}

\usage{
ldlt(A)
}

\arguments{
\item{A}{symmetric positive-definite {\code{\link[base]{matrix}}}}
}

\value{{\code{\link[base]{list}}} of 
\item{$L}{{\code{\link[base]{matrix}}} \eqn{L}}
\item{$D}{{\code{\link[base]{matrix}}} \eqn{D}}
}

\examples{
A <- matrix(c(5,1,1,3),2,2)
V <- ldlt(A)
V$L
V$D
V$L \%*\% V$D \%*\% t(V$L)
A
}

\references{
Zhuang, X., 2020, Numerical analysis. \url{https://sites.ualberta.ca/~xzhuang/Math381/Lab5.pdf}
}
