\name{rmvnorm} 

\alias{rmvnorm}

\title{
	Generates realizations from a multivariate normal distribution
} 

\description{
	Generates realizations from a multivariate normal distribution.
} 

\usage{
	rmvnorm(nsim = 1, mu, V, method = "eigen")
} 

\arguments{
	\item{nsim}{An integer indicating the number of realizations from the distribution.} 
	\item{mu}{A vector of length n containing the mean values of the multivariate normal distribution.} 
	\item{V}{The covariance matrix of the multivariate normal distribution.  The matrix should be symmetric and positive definite.  The size must be \eqn{n times n}.} 
	\item{method}{The method for performing a decomposition of the covariance matrix.  Possible values are "eigen", "chol", and "svd", Eigen value decomposition, Cholesky decomposition, or Singular Value Decomposoition, respectively.} 
} 

\value{
An \eqn{n \times nsim} matrix containing the \code{nsim} realizations of the multivariate normal distribution.  Each column of the matrix represents a realization of the multivariate normal distribution.
}

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{normal}
\keyword{multivariate}

\seealso{rmvnorm}

\examples{
n <- 20
mu <- 1:n
V <- exp(-dist1(matrix(rnorm(n))))
rmvnorm(nsim = 100, mu = mu, V = V, method = "eigen")
}