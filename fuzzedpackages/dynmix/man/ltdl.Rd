\name{ltdl}
\alias{ltdl}

\title{Computes L'DL Matrix Decomposition.}

\description{
This function decomposes matrix \eqn{V} into \eqn{V=L'DL}, where \eqn{L} is a lower triangular matrix with unit diagonal and \eqn{D} is a diagonal matrix with non-negative terms.
}

\usage{
ltdl(A)
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
V <- ltdl(A)
V$L
V$D
t(V$L) \%*\% V$D \%*\% V$L
A
}

\references{
de Jonge, P., Tiberius, C., 1996, \emph{The LAMBDA Method for Integer Ambiguity Estimation: Implementation Aspects}, Universiteitsdrukkerij TU Delft. 
}
