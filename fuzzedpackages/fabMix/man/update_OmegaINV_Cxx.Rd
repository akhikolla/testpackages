\name{update_OmegaINV_Cxx}
\alias{update_OmegaINV_Cxx}
\title{
Gibbs sampling for \eqn{\Omega^{-1}} for Cxx model
}
\description{
Gibbs sampling for \eqn{\Omega^{-1}} for Cxx model
}
\usage{
update_OmegaINV_Cxx(Lambda, K, g, h)
}
\arguments{
  \item{Lambda}{
	Factor loadings, in the form of \eqn{K\times p\times q} matrix, under the restriction that all components share the factor loadings.
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

# Use the real values as input and simulate allocations.
# Mmake sure that in this case Lambda[k,,] is the same  
# for all k = 1,..., K
update_OmegaINV_Cxx(Lambda = syntheticDataset$factorLoadings, 
        K = K, g=0.5, h = 0.5)

}

