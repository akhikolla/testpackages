\name{complete.log.likelihood_q0_sameSigma}
\alias{complete.log.likelihood_q0_sameSigma}
\title{
Complete log-likelihood function for xCx models and \eqn{q=0}
}
\description{
Complete log-likelihood function for models with same  error variance per component (xCx) and diagonal covariance structure per component (\eqn{q=0}.
}
\usage{
complete.log.likelihood_q0_sameSigma(x_data, w, mu, SigmaINV, z)
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
	x_data <- scale(waveDataset1500[ 1:20, -1]) # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( -0.1 + 0.2*runif(K * p), dim = c(K,p) )
	SigmaINV <- array( 1, dim = c(p,p))
	# compute the complete.log.likelihood ( -inf )
	complete.log.likelihood_q0_sameSigma(x_data = x_data, w = w, mu = mu, 
		SigmaINV = SigmaINV, z = z)
}


