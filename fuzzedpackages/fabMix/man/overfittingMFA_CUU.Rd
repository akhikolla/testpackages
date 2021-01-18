\name{overfittingMFA_CUU}
\alias{overfittingMFA_CUU}
\title{
	Basic MCMC sampler for the \code{CUU} model
}
\description{
	Gibbs sampling for fitting a CUU mixture model of factor analyzers.
}
\usage{
overfittingMFA_CUU(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, 
	g, h, alpha_prior, alpha_sigma, beta_sigma, 
	start_values, q, zStart, gibbs_z, lowerTriangular)
}
\arguments{
  \item{x_data}{
	normalized data
}
  \item{originalX}{
	observed raw data (only for plotting purpose)
}
  \item{outputDirectory}{
	Name of the output folder
}
  \item{Kmax}{
	Number of mixture components
}
  \item{m}{
	Number of iterations
}
  \item{thinning}{
	Thinning of chain
}
  \item{burn}{
	Burn-in period
}
  \item{g}{
	Prior parameter \eqn{g}. Default value: \eqn{g = 2}.
}
  \item{h}{
	Prior parameter \eqn{h}. Default value: \eqn{h = 1}.
}
  \item{alpha_prior}{
	Parameters of the Dirichlet prior distribution of mixture weights.
}
  \item{alpha_sigma}{
	Prior parameter \eqn{\alpha}. Default value: \eqn{\alpha = 2}.	
}
  \item{beta_sigma}{
	Prior parameter \eqn{\beta}. Default value: \eqn{\beta = 1}.
}
  \item{start_values}{
	Optional (not used)
}
  \item{q}{
	Number of factors.
}
  \item{zStart}{
	Optional (not used)
}
  \item{gibbs_z}{
	Optional
}
  \item{lowerTriangular}{
	logical value indicating whether a lower triangular parameterization should be imposed on the matrix of factor loadings (if TRUE) or not. Default: TRUE.
}
}
\value{
	Set of files written in \code{outputDirectory}.}
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
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA_CUU(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)

}


