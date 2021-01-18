\name{complete.log.likelihood_q0}
\alias{complete.log.likelihood_q0}
\title{
Complete log-likelihood function for xUx models and \eqn{q=0}
}
\description{
Complete log-likelihood function for models with different error variance per component (xUx) and diagonal covariance structure per component (\eqn{q=0}. 
}
\usage{
complete.log.likelihood_q0(x_data, w, mu, SigmaINV, z)
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
  \item{SigmaINV}{
	\eqn{K\times p\times p} precision matrix (inverse covariance) per component
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
	x_data <- scale(waveDataset1500[ 1:20, -1]) # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( -0.1 + 0.2*runif(K * p), dim = c(K,p) )
	SigmaINV <- array( 1, dim = c(K,p,p))
	# compute the complete.log.likelihood ( -inf )
	complete.log.likelihood_q0(x_data = x_data, w = w, mu = mu, 
		SigmaINV = SigmaINV, z = z)
}


