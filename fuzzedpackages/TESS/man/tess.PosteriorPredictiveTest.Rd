\name{tess.PosteriorPredictiveTest}
\alias{tess.PosteriorPredictiveTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.PosteriorPredictiveTest: Approximation of the posterior predictive distribution.}
\description{tess.PosteriorPredictiveTest computes the values of the statistic for the posterior predictive simulations and computes the p-value for the observed statistic.}
\usage{
tess.PosteriorPredictiveTest(samples,observation,statistic)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{Samples from the posterior predictive distribution.}
  \item{observation}{The observed value.}
  \item{statistic}{The function that computes the statistic.}
}
\value{Returns a list of the statistic for each sample.}
\references{S. Hoehna: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374 }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# We first run an MCMC to obtain samples from the posterior distribution
# and then simulate the posterior predictive distribution.

# The bird phylogeny as the test data set
data(cettiidae)
times <- as.numeric( branching.times(cettiidae) )

# The log-likelihood function
likelihood <- function(params) {
  # We use the parameters as diversification rate and turnover rate.
  # Thus we need to transform first
  b <- params[1] + params[2]
  d <- params[2]
  
  lnl <- tess.likelihood(times,b,d,samplingProbability=1.0,log=TRUE)
  return (lnl)
}

prior_diversification <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior_turnover <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priors <- c(prior_diversification,prior_turnover)

# Note, the number of iterations and the burnin is too small here
# and should be adapted for real analyses
samples <- tess.mcmc(likelihood,priors,c(1,0.1),c(TRUE,TRUE),c(0.1,0.1),10,10)

tmrca <- max(branching.times(cettiidae))
# The simulation function
sim <- function(params) {
  # We use the parameters as diversification rate and turnover rate.
  # Thus we need to transform first
  b <- params[1] + params[2]
  d <- params[2]
  
  # We need trees with at least three tips for the gamma-statistics
  repeat {
    tree <- tess.sim.age(n=1,age=tmrca,b,d,samplingProbability=1.0,MRCA=TRUE)[[1]]
    if (tree$Nnode > 1) break
  }
  return (tree)
}

# simulate trees from the posterior predictive distribution
trees <- tess.PosteriorPrediction(sim,samples)

# compute the posterior predictive test statistic
ppt <- tess.PosteriorPredictiveTest(trees,cettiidae,gammaStat)
# get the p-value of the observed test-statistic
ppt[[2]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ tree }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line
