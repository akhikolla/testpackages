\name{tess.pathSampling}
\alias{tess.pathSampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.pathSampling: Marginal likelihood estimation via Path-Sampling.}
\description{tess.pathSampling uses a power posterior series and path-sampling to estimate the marginal likelihood of a model. This is a very general implementation of this algorithm which can be applied basically to any model. More information can be obtained in the vignette about how to apply this method.}
\usage{
tess.pathSampling(likelihoodFunction,priorFunction,parameters,logTransforms,
                           iterations,burnin=round(iterations/3),K=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{likelihoodFunction}{The log-likelihood function which will be called internally by likelihoodFunction(parameters).}
  \item{priorFunction}{A list of functions of the log-prior-densities of each parameter.}
  \item{parameters}{The initial parameter value list.}
  \item{logTransforms}{A vector of booleans telling if log-transform for the parameters should be used (e.g. for rates).}
  \item{iterations}{The number of iterations for the MCMC.}
  \item{burnin}{The number of iterations to burn before starting the MCMC.}
  \item{K}{The number of stepping stones.}
}
\value{Returns the posterior samples for the parameters.}
\references{
Lartillot, N. and Philippe, H., 2006: Computing Bayes factors using thermodynamic integration. Systematic Biology, 55, 195

Baele et al., 2012: Improving the accuracy of demographic and molecular clock model comparison while accommodating phylogenetic uncertainty

Baele et al., 2013: Accurate Model Selection of Relaxed Molecular Clocks in Bayesian Phylogenetics 
}
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# load a test data set
data(cettiidae)
# convert the phylogeny into the branching times
times <- as.numeric( branching.times(cettiidae) )

# construct a likelihood function taking in a vector of parameters
likelihood <- function(params) {
  # We use the parameters as diversification rate and turnover rate.
  # Thus we need to transform first
  b <- params[1] + params[2]
  d <- params[2]
  
  lnl <- tess.likelihood(times,b,d,samplingProbability=1.0,log=TRUE)
  return (lnl)
}

# next, create the prior density functions
prior_diversification <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior_turnover <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priors <- c(prior_diversification,prior_turnover)

# Note, the number of iterations, the burnin
# and the number of stepping stones is too small here
# and should be adapted for real analyses
marginalLikelihood <- tess.pathSampling( likelihood,
						  priors,
						  runif(2,0,1),
						  c(TRUE,TRUE),
						  10,
						  10,
						  K=4)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line
