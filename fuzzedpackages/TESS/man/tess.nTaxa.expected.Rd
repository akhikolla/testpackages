\name{tess.nTaxa.expected}
\alias{tess.nTaxa.expected}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.nTaxa.expected: The expected number of taxa at present of a tree under a global, time-dependent birth-death process (  E[ N(T) ]  )}
\description{tess.nTaxa.expected computes the expected number of taxa at the present time T (the process start at time s and times increases until the present) under time-dependent. The rates may be any positive function of time or a constant. Additionally, mass-extinction event can be provided and a uniform taxon sampling probability. You have several options for the start of the process (origin vs MRCA). One important feature is that you can compute the expected number of taxa under the reconstructed process, that is, only lineages that survive until the present.}
\usage{
tess.nTaxa.expected( begin,
			   t,
			   end,
			   lambda,
			   mu,
			   massExtinctionTimes=c(),
			   massExtinctionSurvivalProbabilities=c(),
			   samplingProbability=1.0,
			   MRCA=TRUE,
			   reconstructed=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{begin}{The time when the process starts.}
  \item{t}{The time at which we want to know the expected number of lineages (could be equal to end).}
  \item{end}{The time when the process end (e.g. the present).}
  \item{lambda}{The speciation rate function or constant.}
  \item{mu}{The extinction rate function or constant.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{MRCA}{Does the process start with the most recent common ancestor? If not, the tree must have a root edge!}
  \item{reconstructed}{Are we computing the expected number of lineage at time t in the reconstructed process?}
}
\value{Returns the expected number of taxa.}
\references{S. Hoehna: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374 }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# create the time-dependent speciation and extinction rate functions
# here we use episodic functions
l <- function(x) { if (x > 0.5 || x < 0.3) { return (1) } else { return (2) } }
e <- function(x) { if (x > 0.5 || x < 0.3) { return (0.95) } else { return (0.5) } }

# now we can compute the expected number of taxa at time t
# note that we compute here the actual diversity at time t
# if you set reconstructed=TRUE, then you get the expected
# number of lineages that will survive until the present
tess.nTaxa.expected(begin=0,t=2,end=5,l,e,MRCA=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ htest }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line
