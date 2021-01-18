\name{update_SigmaINV_xUC}
\alias{update_SigmaINV_xUC}
\title{
Gibbs sampling for \eqn{\Sigma^{-1}} per component for xUC models
}
\description{
Gibbs sampling for \eqn{\Sigma^{-1}} per component for xUC models
}
\usage{
update_SigmaINV_xUC(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma)
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
	\eqn{K\times p\times p} array containing the inverse variance of errors per component under the restriction: \eqn{\Sigma^{-1}_k = \sigma_k I_p}, where \eqn{\sigma_k > 0}.
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
update_SigmaINV_xUC(x_data = syntheticDataset$data, 
	z = syntheticDataset$class, 
	y = syntheticDataset$factors, 
	Lambda = syntheticDataset$factorLoadings, 
	mu = syntheticDataset$means, 
	K = K, 
	alpha_sigma = 0.5, beta_sigma = 0.5)

}

