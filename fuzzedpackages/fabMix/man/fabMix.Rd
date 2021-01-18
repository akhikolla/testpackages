\name{fabMix}
\alias{fabMix}
\title{
Main function
}
\description{
This function runs parallel chains under a prior tempering scheme of the Dirichlet prior distribution of mixture weights. 
}
\usage{
fabMix(model, nChains, dirPriorAlphas, rawData, outDir, Kmax, mCycles, 
	burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize,  
	thinning, zStart, nIterPerCycle, gibbs_z, 
	warm_up_overfitting, warm_up, overfittingInitialization, 
	progressGraphs, gwar, rmDir, parallelModels, lowerTriangular)
}
\arguments{
  \item{model}{
Any subset of "UUU" "CUU" "UCU" "CCU" "UCC" "UUC" "CUC", "CCC" indicating the fitted models. By default, all models are fitted. 
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
Name of the output folder. An error is thrown if the directory already exists inside the current working directory. Note: it should NOT correspond to an absolute path, e.g.: \code{outDir = `fabMix_example`} is acceptable, but \code{outDir = `C:\Username\Documents\fabMix_example`} is not. 
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
	Model-level parallelization: An optional integer specifying the number of cores that will be used in order to fit in parallel each member of \code{model}. Default: NULL (no model-level parallelization).
}
  \item{lowerTriangular}{
	logical value indicating whether a lower triangular parameterization should be imposed on the matrix of factor loadings (if TRUE) or not. Default: TRUE.
}
}
\value{
	An object of class \code{fabMix.object}, that is, a list consisting of the following entries:
  \item{bic }{Bayesian Information Criterion per model and number of factors.}
  \item{class }{The estimated single best clustering of the observations according to the selected model.}
  \item{n_Clusters_per_model }{The most probable number of clusters (number of non-empty components of the overfitted mixture) per model and number of factors.}
  \item{posterior_probability }{The posterior probability of the estimated allocations according to the selected model.}
  \item{covariance_matrix }{The estimated posterior mean of the covariance matrix per cluster according to the selected model.}
  \item{mu }{The estimated posterior mean of the mean per cluster according to the selected model.}
  \item{weights }{The estimated posterior mean of the mixing proportions according to the selected model.}
  \item{selected_model }{Data frame containing the parameterization, number of clusters and factors of the selected model.}
  \item{mcmc }{A list containing the MCMC draws for the parameters of the selected model. Each entry is returned as an \code{mcmc} object, a class imported from the \code{coda} package (Plummer et al, 2006). All component-specific parameters have been reordered according to the ECR algorithm in order to undo the label switching problem. However, the output corresponding to factor scores and factor loadings is not identifiable due to sign-switching across the MCMC trace.}
  \item{data }{The observed data.}
  \item{regularizedExpression}{The regularized expressions of variable scores to each factor per cluster  (see Papastamoulis 2018, CSDA). }
  \item{Kmap_prob}{The posterior probability of the Maximum A Posteriori number of alive clusters for each parameterization and factor level.}
}
\details{
Let \eqn{\boldsymbol{X}_i}; \eqn{i =1,\ldots,n} denote independent \eqn{p}-dimensional random vectors. Let \eqn{Y_i\in  R^q} with \eqn{q < p} denote the latent factor for observation \eqn{i = 1,\ldots,n}. In the typical factor analysis model, each observation is modelled as \eqn{\boldsymbol{X}_i = \boldsymbol{\mu} + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i }, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma})}, where \eqn{\boldsymbol{\varepsilon}_i} and \eqn{Y_i} are assumed independent, \eqn{i =1,\ldots,n}. The \eqn{p\times q} matrix \eqn{\Lambda} consists of the factor loadings. Assume that there are \eqn{K} clusters and let \eqn{Z_i} denotes the latent allocation of observation \eqn{i} to one amongs the  \eqn{k} clusters, with prior probability \eqn{P(Z_i = k) = w_k}, \eqn{k = 1,\ldots,K}, independent for \eqn{i=1,\ldots,n}.  Following McNicholas et al (2008), the following parameterizations are used:

UUU model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda}_k \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma}_k)}

UCU model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda}_k \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma})}

UUC model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda}_k \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\sigma_k \boldsymbol{I}_p)}

UCC model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda}_k \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\sigma \boldsymbol{I}_p)}

CUU model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma}_k)}

CCU model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma})}

CUC model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\sigma_k \boldsymbol{I}_p)}

CCC model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\sigma \boldsymbol{I}_p)}

In all cases, \eqn{\boldsymbol{\varepsilon}_i} and \eqn{\boldsymbol{Y}_i} are assumed independent, \eqn{i =1,\ldots,n}. Note that \eqn{\boldsymbol{\Sigma}_k} and \eqn{\boldsymbol{\Sigma}} denote positive definite matrices, \eqn{\boldsymbol{I}_p} denotes the \eqn{p\times p} identity matrix and \eqn{\sigma_k}, \eqn{\sigma} denote positive scalars. 

}
\note{
It is recommended to use: \code{normalize = TRUE} (default). Tuning of \code{dirPriorAlphas} may be necessary to achieve reasonable acceptance rates of chain swaps. Note that the output is reordered in order to deal with the label switching problem, according to the ECR algorithm applied by \code{\link{dealWithLabelSwitching}} function. 

Parallelization is enabled in both the chain-level as well as the model-level. By default all heated chains (specified by the \code{nchains} argument) run in parallel using (at most) the same number of threads (if available). If \code{parallelModels = NULL} (default), then the selected parameterizations will run (serially) on the same thread. Otherwise, if \code{parallelModels = x} (where \code{x} denotes a positive integer), the algorithm will first use \code{x} threads to fit the specified model parameterizations in parallel, and furthermore will also parallelize the heated chains (according to the remaining free cores on the user's system). The user should combine \code{parallelModels} with \code{nChains} efficiently, for example: if the number of available threads equals 12 and the user wishes to run 3 heated chains per model (recall that there are 8 parameterizations in total), then, an ideal allocation would be \code{parallelModels = 4} and \code{nChains = 3} because all available threads will be constantly busy. If the user wishes to run \code{nChains = 4} heated chains per model using 12 threads, then an ideal allocation would be \code{parallelModels = 3} models running in parallel. In the case where \code{parallelModels*nChains} > \code{m}, with \code{m} denoting the available number of threads, the algorithm will first allocate min(\code{parallelModels}, \code{m}) threads to run the same number of parameterizations in parallel, and then the remaining threads (if any) will be used to process the parallel heated chains. If no other threads are available, the heated chains will be allocated to single threads. 
}
\seealso{
\code{\link{plot.fabMix.object}}
}
\references{
Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines (2006). CODA: Convergence Diagnosis and Output Analysis for MCMC, R News, vol 6, 7-11.

McNicholas, P.D. and Murphy, T.B. Stat Comput (2008) 18: 285. https://doi.org/10.1007/s11222-008-9056-0.

Papastamoulis, P. (2018). Overfitting Bayesian mixtures of factor analyzers with an unknown number of components. Computational Statistics and Data Analysis, 124: 220-234. DOI: 10.1016/j.csda.2018.03.007.

Papastamoulis, P (2019).  Clustering multivariate data using factor analytic Bayesian mixtures with an unknown number of components. Statistics and Computing, doi: 10.1007/s11222-019-09891-z.
}
\author{
Panagiotis Papastamoulis
}
\examples{
# TOY EXAMPLE (very small numbers... only for CRAN check purposes)

#################################################################
# (a) using 2 cores in parallel, each one running 2 heated chains.
#################################################################
library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2		     # true number of clusters

sINV_diag = 1/((1:p))	 # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
			sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)

# Run `fabMix` for a _small_ number of iterations 
#	but considering only the `UUU` (maximal model),
# 	using the default prior parallel heating parameters `dirPriorAlphas`.
#	NOTE: `dirPriorAlphas` may require some tuning in general.


qRange <- 2	# values for the number of factors (only the true number 
#                                                    is considered here)
Kmax <- 2	# number of components for the overfitted mixture model
nChains <- 2	# number of parallel heated chains

set.seed(1)
fm <- fabMix( model = "UUU", nChains = nChains, 
	rawData = syntheticDataset$data, outDir = "toyExample",
        Kmax = Kmax, mCycles = 4, burnCycles = 1, q = qRange,
        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
        warm_up_overfitting = 2, warm_up = 5) 

# WARNING: the following parameters: 
#  Kmax, nChains, mCycles, burnCycles, warm_up_overfitting, warm_up 
#	 should take (much) _larger_ values. E.g. a typical implementation consists of:
#        Kmax = 20, nChains >= 3, mCycles = 1100, burnCycles = 100, 
#        warm_up_overfitting = 500, warm_up = 5000. 

# You may also print and plot
# print(fm)
# plot(fm)

#################################################################
# (b) using 12 cores_____________________________________________
#_______4 models with 3 heated chains running in parallel________
#_______considering all 8 model parameterizations________________
#################################################################
\dontrun{
library('fabMix')
set.seed(99)
n = 200                # sample size
p = 30                # number of variables
q = 2                # number of factors
K = 5		     # number of clusters
sINV_diag = rep(1/20,p) 	# diagonal of inverse variance of errors
syntheticDataset <- simData(sameLambda=FALSE,K.true = K, n = n, q = q, p = p, 
			sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
qRange <- 1:3	# range of values for the number of factors
Kmax <- 20	# number of components for the overfitted mixture model
nChains <- 3	# number of parallel heated chains

# the next command takes ~ 2 hours in a Linux machine with 12 threads.

fm <- fabMix( parallelModels = 4, 
	nChains = nChains, 
	model = c("UUU","CUU","UCU","CCU","UCC","UUC","CUC","CCC"), 
	rawData = syntheticDataset$data, outDir = "toyExample_b",
        Kmax = Kmax, mCycles = 1100, burnCycles = 100, q = qRange) 

print(fm)
plot(fm, what = "BIC")
# see also
# plot(fm); summary(fm)

}


}


