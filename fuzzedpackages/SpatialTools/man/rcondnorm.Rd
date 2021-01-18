\name{rcondnorm} 

\alias{rcondnorm}

\title{
	Generate from conditional normal distribution
} 

\description{
	Generates realizations from a multivariate normal distribution conditional on observed data vector
} 

\usage{
	rcondnorm(nsim = 1, y, mu, mup, V, Vp, Vop, method = "eigen")
} 

\arguments{
	\item{nsim}{An integer indicating the number of realizations from the distribution.} 
	\item{y}{A vector of length \code{n} contained the observed data.} 
	\item{mu}{The mean vector of the observed data.  Should be a vector of length \code{n}.} 
	\item{mup}{The mean vector of the responses to be generated.  Should be a vector of length \code{np}.} 
	\item{V}{The covariance matrix of the observed data.  The matrix should be symmetric and positive definite.  The size must be \eqn{n times n}.} 
	\item{Vp}{The covariance matrix of the responses to be generated.  The matrix should be symmetric and positive definite.  The size must be \eqn{np times np}.} 
	\item{Vop}{The cross-covariance matrix between the observed data and the responses to be generated.  The size must be \eqn{n times np}.} 
	\item{method}{The method for performing a decomposition of the covariance matrix.  Possible values are "eigen", "chol", and "svd", Eigen value decomposition, Cholesky decomposition, or Singular Value Decomposoition, respectively.} 
} 

\value{
An \eqn{np \times nsim} matrix containing the \code{nsim} realizations of the conditional normal distribution.  Each column of the matrix represents a realization of the multivariate normal distribution.
}

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{normal}
\keyword{multivariate}
\keyword{conditional}

\seealso{rmvnorm}

\examples{
n <- 100
np <- 100

mu <- rep(1, 100)
mup <- rep(2, 100)

coords <- matrix(runif(2 * n), ncol = 2)
pcoords <- matrix(runif(2 * np), ncol = 2)

myV <- cov.sp(coords, sp.type = "exponential", c(1, 2), error.var = 1, pcoords = pcoords)

y <- rmvnorm(1, mu = mu, V = myV$V)

rcondnorm(3, y = y, mu = mu, mup = mup, V = myV$V, Vp = myV$Vp, Vop = myV$Vop, method = "chol")
}