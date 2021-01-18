\name{update_z_b_Sj}
\alias{update_z_b_Sj}
\title{
Gibbs sampling for \eqn{z}
}
\description{
Gibbs sampling for \eqn{z}: here the full conditional distribution is being used (that is, the distribution is also conditioned on the values of factors \eqn{y}).
}
\usage{
update_z_b_Sj(w, mu, Lambda, y, SigmaINV, K, x_data)
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
  \item{y}{
	\eqn{n\times q} Matrix of factors
}
  \item{SigmaINV}{
	\eqn{K\times p\times p} array containing the precision matrix per component
}
  \item{K}{
	Number of components
}
  \item{x_data}{
	\eqn{n\times p} matrix containing the observed data
}
}
\value{
	Allocation vector
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
SigmaINV <- array(data = 0, dim = c(K,p,p))
for(k in 1:K){
	diag(SigmaINV[k,,]) <- 1/diag(syntheticDataset$variance) + rgamma(p, shape=1, rate = 1)
}

# use the real values as input and simulate allocations
update_z_b_Sj(w = syntheticDataset$weights, mu = syntheticDataset$means, 
	Lambda = syntheticDataset$factorLoadings, 
	y = syntheticDataset$factors,
	SigmaINV = SigmaINV, 
	K = K, x_data = syntheticDataset$data)$z

}
