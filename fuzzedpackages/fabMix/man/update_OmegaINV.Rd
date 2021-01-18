\name{update_OmegaINV}
\alias{update_OmegaINV}
\title{
Gibbs sampling for \eqn{\Omega^{-1}}
}
\description{
Gibbs sampling for \eqn{\Omega^{-1}}
}
\usage{
update_OmegaINV(Lambda, K, g, h)
}
\arguments{
  \item{Lambda}{
	Factor loadings
}
  \item{K}{
	Number of components
}
  \item{g}{
	Prior parameter
}
  \item{h}{
	Prior parameter
}
}
\value{
	\eqn{q\times q} matrix \eqn{\Omega^{-1}}
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
update_OmegaINV(Lambda = syntheticDataset$factorLoadings, 
        K = K, g=0.5, h = 0.5)

}

