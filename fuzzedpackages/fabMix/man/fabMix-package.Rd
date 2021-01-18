\name{fabMix-package}
\alias{fabMix-package}
\docType{package}
\title{
\packageTitle{fabMix}
}
\description{
\packageDescription{fabMix}

The main fuction of the package is \code{\link{fabMix}}. 
}
\author{
Panagiotis Papastamoulis

Maintainer: \packageMaintainer{fabMix}
}
\references{
Fokoue, E. and Titterington, D.M. (2003). Mixtures of Factor Analysers: Bayesian Estimation and Inference by Stochastic Simulation. Machine Learing, 50(1): 73-94.

McNicholas, P.D. and Murphy, T.B. Statistics and Computing (2008) 18: 285. https://doi.org/10.1007/s11222-008-9056-0.

Papastamoulis P. and Iliopoulos G. (2010). An artificial allocations based solution to the label switching problem in Bayesian analysis of mixtures of distributions. Journal of Computational and Graphical Statistics, 19: 313-331.

Rousseau, J. and Mengersen, K. (2011). Asymptotic behaviour of the posterior distribution in overfitted mixture models. Journal of the Royal Statistical Society, Series B (methodological), 73(5): 689-710.

van Havre, Z., White, N., Rousseau, J. and Mengersen, K. (2015). Overfitting Bayesian Mixture Models with an Unknown Number of Components. PLOS ONE, 10(7): 1-27.

Papastamoulis, P. (2016). \code{label.switching}: An R Package for Dealing with the Label Switching Problem in MCMC Outputs. Journal of Statistical Software, 69(1), 1-24.

Papastamoulis, P. (2018). Overfitting Bayesian mixtures of factor analyzers with an unknown number of components. Computational Statistics and Data Analysis, 124: 220-234. DOI: 10.1016/j.csda.2018.03.007.

Papastamoulis, P (2019).  Clustering multivariate data using factor analytic Bayesian mixtures with an unknown number of components. Statistics and Computing, doi: 10.1007/s11222-019-09891-z.
}
\keyword{ package }
\seealso{
\code{\link{fabMix}}, \code{\link{plot.fabMix.object}}
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
#	considering only `UUU` (maximal model),
# 	using the default prior parallel heating parameters `dirPriorAlphas`.
#	NOTE: `dirPriorAlphas` may require some tuning in general.


qRange <- 2	# values for the number of factors (only the true number 
#                                                    is considered here)
Kmax <- 2	# number of components for the overfitted mixture model
nChains <- 2	# number of parallel heated chains

set.seed(1)
fm <- fabMix( model = c("UUU"), nChains = nChains, 
	rawData = syntheticDataset$data, outDir = "toyExample",
        Kmax = Kmax, mCycles = 4, burnCycles = 1, q = qRange,
        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
        warm_up_overfitting = 2, warm_up = 5) 

# WARNING: the following parameters: 
#  Kmax, nChains, mCycles, burnCycles, warm_up_overfitting, warm_up 
#	 should take (much) _larger_ values. E.g. a typical implementation consists of:
#        Kmax = 20, nChains >= 3, mCycles = 1100, burnCycles = 100, 
#        warm_up_overfitting = 500, warm_up = 5000. 

# Now print a run summary and produce some plots. 
print(fm)
# you may also plot, summary the output. 

#################################################################
# (b) using 12 cores_____________________________________________
#_______4 models with 3 heated chains running in parallel________
#_______considering all 8 model parameterizations________________
#################################################################
\dontrun{
library('fabMix')
set.seed(99)
n = 100                # sample size
p = 30                # number of variables
q = 2                # number of factors
K = 5		     # number of clusters
sINV_diag = rep(1/100,p) 	# diagonal of inverse variance of errors
syntheticDataset <- simData(sameLambda=FALSE,K.true = K, n = n, q = q, p = p, 
			sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
qRange <- 1:3	# range of values for the number of factors
Kmax <- 20	# number of components for the overfitted mixture model
nChains <- 3	# number of parallel heated chains

# the next command takes ~ 1 hour in a Linux workstation with 12 threads.
fm <- fabMix( parallelModels = 4, 
	nChains = nChains, 
	model = c("UUU","CUU","UCU","CCU","UCC","UUC","CUC","CCC"), 
	rawData = syntheticDataset$data, outDir = "toyExample_b",
        Kmax = Kmax, mCycles = 600, burnCycles = 100, q = qRange,
        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
        warm_up_overfitting = 500, warm_up = 5000) 
print(fm)
plot(fm, what = "BIC")
plot(fm, what = "classification_pairs")

}

}
