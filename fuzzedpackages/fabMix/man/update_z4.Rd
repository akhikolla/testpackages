\name{update_z4}
\alias{update_z4}
\title{
	Collapsed Gibbs for \eqn{z}
}
\description{
	Collapsed Gibbs for \eqn{z}.
}
\usage{
update_z4(w, mu, Lambda, SigmaINV, K, x_data)
}
\arguments{
  \item{w}{
	vector with length \eqn{K} consisting of mixture weights
}
  \item{mu}{
	\eqn{K\times p} array containing the marginal means
}
  \item{Lambda}{
	\eqn{K\times p} array with factor loadings
}
  \item{SigmaINV}{
	\eqn{p\times p} precision matrix
}
  \item{K}{
	Number of components
}
  \item{x_data}{
	\eqn{n\times p} matrix containing the observed data
}
}
\value{
	A vector of length \eqn{n} with the simulated allocation of each observation among the \eqn{K} components.
}
\author{
Panagiotis Papastamoulis
}
\examples{
library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
# use the real values as input and simulate allocations
update_z4(w = syntheticDataset$weights, mu = syntheticDataset$means, 
	Lambda = syntheticDataset$factorLoadings, 
	SigmaINV = diag(1/diag(syntheticDataset$variance)), 
	K = K, x_data = syntheticDataset$data)$z

}
