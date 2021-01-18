\name{fabMix_parallelModels}
\alias{fabMix_parallelModels}
\title{
Function for model-level parallelization
}
\description{
This function runs multiple copies of the \code{fabMix} function in parallel. 
}
\usage{
fabMix_parallelModels(model, nChains, dirPriorAlphas, rawData, outDir, Kmax, mCycles, 
	burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize,  
	thinning, zStart, nIterPerCycle, gibbs_z, 
	warm_up_overfitting, warm_up, overfittingInitialization, 
	progressGraphs, gwar, rmDir, parallelModels, lowerTriangular)
}
\arguments{
  \item{model}{
Any subset of "UUU" "CUU" "UCU" "CCU" "UCC" "UUC" "CUC", "CCC" indicating the fitted models.  
}
  \item{nChains}{
The number of parallel heated chains. When `dirPriorAlphas` is supplied, `nChains` can be ignored.
}
  \item{dirPriorAlphas}{
vector of length \code{nChains} in the form of an increasing sequence of positive scalars. Each entry contains the (common) prior Dirichlet parameter for the corresponding chain. Default: \code{dirPriorAlphas = c(1, 1 + dN*(2:nChains - 1))/Kmax}, where \code{dN = 1}, for \code{nChains > 1}. Otherwise, \code{dirPriorAlphas = 1/Kmax}. 
}
  \item{rawData}{
The observed data as an \eqn{n\times p} matrix. Clustering is performed on the rows of the matrix.
}
  \item{outDir}{
Name of the output folder. An error is thrown if this directory already exists.
}
  \item{Kmax}{
Number of components in the overfitted mixture. Default: 20.
}
  \item{mCycles}{
Number of MCMC cycles. Each cycle consists of \code{nIterPerCycle} MCMC iterations. At the end of each cycle a swap of the state of two randomly chosen adjacent chains is attempted.
}
  \item{burnCycles}{
Number of cycles that will be discarded as burn-in period.
}
  \item{g}{
	Prior parameter \eqn{g}. Default value: \eqn{g = 0.5}.
}
  \item{h}{
	Prior parameter \eqn{h}. Default value: \eqn{h = 0.5}.
}
  \item{alpha_sigma}{
	Prior parameter \eqn{\alpha}. Default value: \eqn{\alpha = 0.5}.
}
  \item{beta_sigma}{
	Prior parameter \eqn{\beta}. Default value: \eqn{\beta = 0.5}.
}
  \item{q}{
	A vector containing the number of factors to be fitted.
}
  \item{normalize}{
	Should the observed data be normalized? Default value: TRUE. (Recommended) 
}
  \item{thinning}{
	Optional integer denoting the thinning of the keeped MCMC cycles.
}
  \item{zStart}{
	Optional starting value for the allocation vector.
}
  \item{nIterPerCycle}{
	Number of iteration per MCMC cycle. Default value: 10.
}
  \item{gibbs_z}{
	Select the gibbs sampling scheme for updating latent allocations of mixture model. Default value: 1.
}
  \item{warm_up_overfitting}{
	Number of iterations for the overfitting initialization scheme. Default value: 500.
}
  \item{warm_up}{
	Number of iterations that will be used to initialize the models before starting proposing switchings. Default value: 5000.
}
  \item{overfittingInitialization}{
	Logical value indicating whether the chains are initialized via the overfitting initialization scheme. Default: TRUE.
}
  \item{progressGraphs}{
	Logical value indicating whether to plot successive states of the chains while the sampler runs. Default: FALSE.
}
  \item{gwar}{
	Initialization parameter. Default: 0.05.
}
  \item{rmDir}{
	Logical value indicating whether to delete the \code{outDir} directory. Default: TRUE.
}
  \item{parallelModels}{
	Model-level parallelization: An optional integer specifying the number of cores that will be used in order to fit in parallel each member of \code{model}. 
}
  \item{lowerTriangular}{
	logical value indicating whether a lower triangular parameterization should be imposed on the matrix of factor loadings (if TRUE) or not. Default: TRUE.
}
}
\value{
	An object of class \code{fabMix.object} (see the \code{\link{fabMix}} function).
}
\author{
Panagiotis Papastamoulis
}
\note{
See the \code{\link{fabMix}} function for examples.
}
