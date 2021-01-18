\name{tess.mcmc}
\alias{tess.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.mcmc: Markov chain Monte Carlo simulation using a general Metropolis-Hastings algorithm.}
\description{tess.mcmc constructs a Markov chain Monte Carlo sampler (MCMC) by implementing a general Metropolis-Hastings algorithm. Any model can be used where the likelihood is known and thus can be passed in as an argument. The parameters have to be continuous. Proposals are taken from a normal distribution centered around the current value. The varaince of the new proposed values is initalized with 1 but can be automatically optimized when using the option adaptive = TRUE. The algorithm creates sampels from the posterior probility distribution and returns these a CODA mcmc object. More information can be obtained in the vignette about how to apply this method.}
\usage{
tess.mcmc(likelihoodFunction,priors,parameters,logTransforms,delta,
             iterations,burnin=round(iterations/3),thinning=1,
             adaptive=TRUE,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{likelihoodFunction}{The log-likelihood function which will be called internally by likelihoodFunction(parameters).}
  \item{priors}{A list of functions of the log-prior-densities of each parameter.}
  \item{parameters}{The initial parameter value list.}
  \item{logTransforms}{A vector of booleans telling if log-transform for the parameters should be used (e.g. for rates).}
  \item{delta}{The variance of new proposed values.}
  \item{iterations}{The number of iterations for the MCMC.}
  \item{burnin}{The number of iterations to burn before starting the MCMC.}
  \item{thinning}{The frequency of taking a sample of the parameters.}
  \item{adaptive}{Should we use adaptive MCMC?}
  \item{verbose}{Do you want detailed information during the run?}
}
\value{Returns the posterior samples for the parameters.}
\references{
S. Hoehna: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374.

S. Hoehna, MR May and BR Moore: TESS: Bayesian inference of lineage diversification rates from (incompletely sampled) molecular phylogenies in R. 2015, Bioinformatics.
}
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# load in a test data set
data(cettiidae)

# convert the phylogeny into the branching times
times <- as.numeric( branching.times(cettiidae) )

# specify a likelihood function that takes in a vector of parameters
likelihood <- function(params) {
  # We use the parameters as diversification rate and turnover rate.
  # Thus we need to transform first
  b <- params[1] + params[2]
  d <- params[2]
  
  lnl <- tess.likelihood(times,b,d,samplingProbability=1.0,log=TRUE)
  return (lnl)
}

# specify a the prior functions
prior.diversification <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior.turnover <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priors <- c(prior.diversification,prior.turnover)

# Note, the number of iterations and the burnin is too small here
# and should be adapted for real analyses
samples <- tess.mcmc( likelihood,
		      priors,
		      runif(2,0,1),
		      logTransforms=c(TRUE,TRUE),
		      delta=c(0.1,0.1),
		      iterations=100,
		      burnin=20)

# now summarize and visualize the results
#plot(samples)
summary(samples)
colMeans(samples)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ tree }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line
