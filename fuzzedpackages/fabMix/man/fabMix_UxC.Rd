\name{fabMix_UxC}
\alias{fabMix_UxC}
\title{
Function to estimate the \code{UUC} and \code{UCC} models
}
\description{
This function runs parallel chains under a prior tempering scheme of the Dirichlet prior distribution of mixture weights. 
}
\usage{
fabMix_UxC(sameSigma, dirPriorAlphas, rawData, outDir, Kmax, mCycles, 
	burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize,  
	thinning, zStart, nIterPerCycle, gibbs_z, 
	warm_up_overfitting, warm_up, overfittingInitialization, 
	progressGraphs, gwar, lowerTriangular)
}
\arguments{
  \item{sameSigma}{
Logical value denoting the parameterization of the error variance per component. If \code{TRUE}, the parameterization CCU is fitted. Otherwise, the parameterization CUU is fitted.
}
  \item{dirPriorAlphas}{
The prior Dirichlet parameters for each chain.
}
  \item{rawData}{
The observed data as an \eqn{n\times p} matrix. Clustering is performed on the rows of the matrix.
}
  \item{outDir}{
Name of the output folder.
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
	Prior parameter \eqn{g}. Default value: \eqn{g = 2}.
}
  \item{h}{
	Prior parameter \eqn{h}. Default value: \eqn{h = 1}.
}
  \item{alpha_sigma}{
	Prior parameter \eqn{\alpha}. Default value: \eqn{\alpha = 2}.
}
  \item{beta_sigma}{
	Prior parameter \eqn{\beta}. Default value: \eqn{\beta = 1}.
}
  \item{q}{
	Number of factors \eqn{q}, where \eqn{1 \leq q \leq L}. An error is thrown if the Ledermann bound (\eqn{L}) is exceeded.
}
  \item{normalize}{
	Should the observed data be normalized? Default value: TRUE. 
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
	Number of iterations for the overfitting initialization scheme. Default value: 100.
}
  \item{warm_up}{
	Number of iterations that will be used to initialize the models before starting proposing switchings. Default value: 500.
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
  \item{lowerTriangular}{
	logical value indicating whether a lower triangular parameterization should be imposed on the matrix of factor loadings (if TRUE) or not. Default: TRUE.
}
}
\value{
	List of files written to \code{outDir}
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\note{
It is recommended to always use: \code{normalize = TRUE} (default). Tuning of \code{dirPriorAlphas} may be necessary to achieve reasonable acceptance rates of chain swaps. Also note that the output is not identifiable due to label switching and the user has to subsequently call the \code{dealWithLabelSwitching} function. See the \code{\link{fabMix}} function for examples.
}
\seealso{
\code{\link{fabMix}}
}
\author{
Panagiotis Papastamoulis
}

