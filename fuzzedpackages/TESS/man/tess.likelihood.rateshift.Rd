\name{tess.likelihood.rateshift}
\alias{tess.likelihood.rateshift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.likelihood.rateshift: Probability density of a tree under a tree-wide time-dependent birth-death-shift process}
\description{tess.likelihood.rateshift computes the probability of a reconstructed phylogenetic tree under a rate-shift model. The rates are piecewise constant. Additionally, mass-extinction event can be provided and a uniform taxon sampling probability. You have several options for the start of the process (origin vs MRCA) and the condition of the process (time, survival or taxa; note that survival and taxa implicitly condition on the time too!). See equation (5) in the manuscript for more information. Note that constant rates lead to much faster computations.
The likelihood can be computed for incompletely sampled trees. You need to give a sampling probability != 1.0. You have three options for the sampling strategy: uniform|diversified|age. The detailed description of these can be found in the references. More information can be obtained in the vignette about how to apply this likelihood function.}
\usage{
tess.likelihood.rateshift( times,
				 lambda, 
				 mu, 
				 rateChangeTimesLambda = c(),
				 rateChangeTimesMu = c(),
				 massExtinctionTimes = c(), 
   				 massExtinctionSurvivalProbabilities = c(),
                 missingSpecies = c(),
                 timesMissingSpecies = c(),
				 samplingStrategy = "uniform", 
   				 samplingProbability = 1, 
				 MRCA = TRUE, 
   				 CONDITION = "survival", 
				 log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{The branching times of the tree.}
  \item{lambda}{The speciation rate as a vector representing the rate for each time interval.}
  \item{mu}{The extinction rate as a vector representing the rate for each time interval.}
  \item{rateChangeTimesLambda}{The times of the rate-shifts for the speciation rate.}
  \item{rateChangeTimesMu}{The times of the rate-shifts for the extinction rate.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{missingSpecies}{The number of species missed which originated in a given time interval (empirical taxon sampling).}
  \item{timesMissingSpecies}{The times intervals of the missing species (empirical taxon sampling).}
  \item{samplingStrategy}{The strategy how samples were obtained. Options are: uniform|diversified|age.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{MRCA}{Does the process start with the most recent common ancestor? If not, the tree must have a root edge!}
  \item{CONDITION}{do we condition the process on time|survival|taxa?}
  \item{log}{should the likelihood be in log-scale?}
}
\value{Returns the (log) probability of the tree, i.e., the likelihood of the parameters given the tree.}
\references{S. Hoehna: The time-dependent reconstructed evolutionary process with a key-role for mass-extinction events. 2015, Journal of Theoretical Biology, 380, 321-331.}
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# load a test data set
data(cettiidae)

# convert the phylogeny into the branching times
times <- as.numeric( branching.times(cettiidae) )

# construct speciation and extinction rate function that resemble the rate-shift
# any other function could be used too
l <- function(x) { if (x > 0.5 || x < 0.3) { return (1) } else { return (2) } }
e <- function(x) { if (x > 0.5 || x < 0.3) { return (0.95) } else { return (0.5) } }

# now compute the likelihood for the tree
tess.likelihood(times,l,e,MRCA=TRUE,log=TRUE)

# a second approach is the specific episodic birth-death process likelihood function
# we need to give the rates for each episode and the end time of the episodes
# you should see that both are equivalent in this setting
# the function approach is more general but also slower.
tess.likelihood.rateshift(times,
				lambda=c(2,1,2),
				mu=c(0.95,0.5,0.95),
				rateChangeTimesLambda=c(0.3,0.5),
				rateChangeTimesMu=c(0.3,0.5),
				MRCA=TRUE,
				log=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line
