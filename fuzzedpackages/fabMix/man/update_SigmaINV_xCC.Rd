\name{update_SigmaINV_xCC}
\alias{update_SigmaINV_xCC}
\title{
Gibbs sampling for \eqn{\Sigma^{-1}} for xCC models
}
\description{
Gibbs sampling for \eqn{\Sigma^{-1}} for xCC models
}
\usage{
update_SigmaINV_xCC(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma)
}
\arguments{
  \item{x_data}{
	\eqn{n\times p} matrix containing the observed data
}
  \item{z}{
	Allocation vector
}
  \item{y}{
	\eqn{n\times q} matrix containing the latent factors
}
  \item{Lambda}{
	\eqn{K\times p\times q} array with factor loadings
}
  \item{mu}{
	\eqn{K\times p} array containing the marginal means
}
  \item{K}{
	Number of components
}
  \item{alpha_sigma}{
	Prior parameter \eqn{alpha}
}
  \item{beta_sigma}{
	Prior parameter \eqn{beta}
}
}
\value{
	\eqn{p\times p} matrix with the common variance of errors per component \eqn{\Sigma^{-1} = \sigma I_p}.
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

# use the real values as input and update SigmaINV
update_SigmaINV_xCC(x_data = syntheticDataset$data, 
	z = syntheticDataset$class, 
	y = syntheticDataset$factors, 
	Lambda = syntheticDataset$factorLoadings, 
	mu = syntheticDataset$means, 
	K = K, 
	alpha_sigma = 0.5, beta_sigma = 0.5)

}

