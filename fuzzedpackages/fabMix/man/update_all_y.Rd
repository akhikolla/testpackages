\name{update_all_y}
\alias{update_all_y}
\title{
Gibbs sampling for \eqn{y} in \code{xCx} model
}
\description{
Gibbs sampling for updating the factors \eqn{y} for models with same variance of errors per component.
}
\usage{
update_all_y(x_data, mu, SigmaINV, Lambda, z)
}
\arguments{
  \item{x_data}{
	\eqn{n\times p} matrix with obseved data
}
  \item{mu}{
	\eqn{n\times p} matrix of marginal means
}
  \item{SigmaINV}{
	\eqn{p\times p} precision matrix
}
  \item{Lambda}{
	\eqn{p\times q} matrix of factor loadings
}
  \item{z}{
	Allocation vector
}
}
\value{
	A matrix with generated factors
}
\author{
Panagiotis Papastamoulis
}
\examples{
library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
# use the real values as input and simulate factors
update_all_y(x_data = syntheticDataset$data, 
		mu = syntheticDataset$means, 
		SigmaINV = diag(1/diag(syntheticDataset$variance)), 
		Lambda = syntheticDataset$factorLoadings, 
		z = syntheticDataset$class)

}
