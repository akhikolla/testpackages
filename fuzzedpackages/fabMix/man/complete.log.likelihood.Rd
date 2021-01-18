\name{complete.log.likelihood}
\alias{complete.log.likelihood}
\title{
Complete log-likelihood function for xCx models.
}
\description{
Complete log-likelihood function for models with same error variance per component (xCx).
}
\usage{
complete.log.likelihood(x_data, w, mu, Lambda, SigmaINV, z)
}

\arguments{
  \item{x_data}{
	\eqn{n\times p} matrix containing the data
}
  \item{w}{
	a vector of length \eqn{K} containing the mixture weights
}
  \item{mu}{
	\eqn{K\times p} matrix containing the marginal means per component
}
  \item{Lambda}{
	\eqn{K\times p\times q} array of factor loadings per component
}
  \item{SigmaINV}{
	\eqn{p\times p} precision matrix (inverse covariance)
}
  \item{z}{
	A vector of length \eqn{n} containing the allocations of the \eqn{n} datapoints to the \eqn{K} mixture components
}
}
\value{
	complete log-likelihood value
}
\author{
Panagiotis Papastamoulis
}
\examples{
	library('fabMix')
	data(waveDataset1500)
	x_data <- waveDataset1500[ 1:20, -1] # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( runif(K * p), dim = c(K,p) )
	Lambda <- array( runif(K*p*q), dim = c(K,p,q) )
	SigmaINV <- array(1, dim = c(p,p))
	# compute the complete.log.likelihood
	complete.log.likelihood(x_data = x_data, w = w, mu = mu, 
		Lambda = Lambda, SigmaINV = SigmaINV, z = z)
}


