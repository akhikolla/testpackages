\name{tess.sim.taxa}
\alias{tess.sim.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.sim.taxa.taxa: Simulate a reconstructed tree for a given number of taxa under a global, time-dependent birth-death process.}
\description{tess.sim.taxa simulates a reconstructed phylogenetic tree under a global, time-dependent birth-death process conditioned on the number of taxa sampled. The rates may be any positive function of time or a constant. The process starts at time 0 and goes forward in time, hence the rates and events should be interpreted in the time after the origin. Additionally, mass-extinction event can be provided and a uniform taxon sampling probability. It is possible to start either with the origin (1 species) or with the most recent common ancestor (2 species).}
\usage{
tess.sim.taxa(n, nTaxa, max, lambda, mu, massExtinctionTimes = c(), 
   massExtinctionSurvivalProbabilities = c(), samplingProbability = 1, 
   samplingStrategy = "uniform", SURVIVAL = TRUE, MRCA = TRUE, t_crit = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of simulations.}
  \item{nTaxa}{Number of species sampled.}
  \item{max}{Maximum time/height of the tree.}
  \item{lambda}{The speciation rate function or constant.}
  \item{mu}{The extinction rate function or constant.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{samplingStrategy}{The strategy how samples were obtained. Options are: uniform|diversified.}
  \item{SURVIVAL}{Do you want to condition on survival of the process?}
  \item{MRCA}{Does the process start with the most recent common ancestor?}
  \item{t_crit}{The critical time points when a jump in the rate function occurs. Only a help for the numerical integration routine.}
}
\value{Returns a tree in 'phylo' format.
}
\references{S. Hoehna: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374 }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

l <- function(x) { if (x > 0.5 || x < 0.3) { return (1) } else { return (2) } }
e <- function(x) { if (x > 0.5 || x < 0.3) { return (0.95) } else { return (0.5) } }

tess.sim.taxa(n=1,nTaxa=10,max=10,l,e,MRCA=TRUE)

# simulation under constant rates
tess.sim.taxa(n=1,nTaxa=10,max=10,2.0,1.0,MRCA=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ likelihood }
%\keyword{ tree }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line
