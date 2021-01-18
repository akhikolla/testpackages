\name{compute_sufficient_statistics_given_mu}
\alias{compute_sufficient_statistics_given_mu}
\title{
Compute sufficient statistics given mu
}
\description{
Compute sufficient statistics given \eqn{y}, \eqn{z} and \eqn{\mu}.
}
\usage{
compute_sufficient_statistics_given_mu(y, z, K, x_data,mu)
}
\arguments{
  \item{y}{
	\eqn{n\times q} matrix of factors
}
  \item{z}{
	Allocation vector
}
  \item{K}{
	Number of components
}
  \item{x_data}{
	\eqn{n\times p} matrix with observed data
}
  \item{mu}{
	\eqn{K\times p} matrix with marignal means per component
}
}
\value{
A list with six entries of sufficient statistics.
	\item{cluster_size}{ Integer vector of length \eqn{K} }
	\item{sx}{ \eqn{K\times p} array }
	\item{sy}{ \eqn{K\times q} array }
	\item{sxx}{ Not used }
	\item{syy}{ \eqn{K\times q \times q} array }
	\item{sxy}{ \eqn{K\times p \times q} array }
}
\author{
Panagiotis Papastamoulis
}
\examples{
        data(waveDataset1500)
        x_data <- as.matrix(waveDataset1500[ 1:20, -1]) # data
        z <-  waveDataset1500[ 1:20, 1] # class
        p <- dim(x_data)[2]
        n <- dim(x_data)[1]
        q <- 2
        K <- length(table(z))           # 3 classes
        # give some arbitrary values to the parameters:
        set.seed(1)
        mu <- array( runif(K * p), dim = c(K,p) )
	y <- array(rnorm(n = q*n), dim = c(n,q))
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics_given_mu(y = y, 
	 z = z, K = K, x_data = x_data, mu = mu)
}

