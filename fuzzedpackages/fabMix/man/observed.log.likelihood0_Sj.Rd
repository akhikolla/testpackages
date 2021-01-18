\name{observed.log.likelihood0_Sj}
\alias{observed.log.likelihood0_Sj}
\title{
Log-likelihood of the mixture model
}
\description{
Log-likelihood of the mixture model evaluated only at the alive components.
}
\usage{
observed.log.likelihood0_Sj(x_data, w, mu, Lambda, Sigma, z)
}
\arguments{
  \item{x_data}{
The observed data
}
  \item{w}{
	Vector of mixture weights
}
  \item{mu}{
	Vector of marginal means
}
  \item{Lambda}{
	Factor loadings
}
  \item{Sigma}{
	\eqn{K\times p} matrix with each row containing the diagonal of the covariance matrix of the errors per cluster
}
  \item{z}{
	Allocation vector
}
}
\value{
	Log-likelihood value
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
	Sigma <- matrix(1:K, nrow = K, ncol = p)
	# compute the complete.log.likelihood
	observed.log.likelihood0_Sj(x_data = x_data, w = w, 
		mu = mu, Lambda = Lambda, Sigma = Sigma, z = z)
}

