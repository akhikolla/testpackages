\name{HACHypothetical}
\alias{HACHypothetical}

\title{Set up an object to simulate haplotype accumulation curves for a hypothetical species}

\description{Helper function which creates an object containing necessary information to run a simulation of haplotype accumulation for a hypothetical species of interest}

\usage{HACHypothetical(N, Hstar, probs, perms = 10000, p = 0.95, 
conf.level = 0.95, subsample = FALSE, prop = NULL,
progress = TRUE, num.iters = NULL, filename = NULL)}

\arguments{\item{N}{Number of individuals}
\item{Hstar}{Number of unique species' haplotypes}
\item{probs}{Haplotype frequency distribution vector}
\item{perms}{Number of permutations (replications)}
\item{p}{Proportion of haplotypes to recover}
\item{conf.level}{Desired confidence level for graphical output and interval estimation}
\item{subsample}{Is a subsample of haplotype labels desired?}
\item{prop}{If subsample = TRUE, the proportion of haplotype labels to subsample}
\item{num.iters}{Number of iterations to compute}
\item{progress}{Should iteration output be printed to the R console?}
\item{filename}{Name of file where simulation results are to be saved}
}

\value{A list object of class "HAC" with 13 elements that can be passed to \code{HAC.simrep} as follows:

  \item{input.seqs}{Should a FASTA file of aligned/trimmed DNA sequences be inputted? Default is FALSE}
  \item{subset.seqs}{Should a subsample of DNA sequences be taken? Default is FALSE}
  \item{prop.seqs}{Proportion of DNA sequences to subsample. Default is NULL}
  \item{prop.haps}{Proportion of haplotype labels to subsample. Default is NULL (can be altered by user)}
  \item{subset.haps}{Should a subsample of haplotype labels be taken? Default is NULL (can be altered by user)}
  \item{N}{Number of individuals. NA by default (provided by user)}
  \item{Hstar}{Number of unique species' haplotypes. NA by default (provided by user)}
  \item{probs}{Haplotype frequency distribution vector. NA by default (provided by user)}
  \item{p}{Proportion of haplotypes to recover. \code{p} = 0.95 by default.}
  \item{perms}{Number of permutations (replications). \code{perms} = 10000 by default.}
  \item{conf.level}{Desired confidence level for graphical output and interval estimation. \code{conf.level} = 0.95 by default.}
  \item{num.iters}{Number of iterations to compute. \code{num.iters = NULL} by default (i.e., all iterations are computed; users can specify \code{num.iters} = 1 for the first iteration.)}
  \item{progress}{Should iteration output be printed to the R console? Default is TRUE}
  \item{filename}{Name of file where simulation results are to be saved.}
}

\note{\code{N} must be greater than 1 and greater than or equal to \code{Hstar}.

\code{Hstar} must be greater than 1. 

\code{probs} must have a length equal to \code{Hstar} and its elements must sum to 1.}

\examples{
  ## Simulate hypothetical species ##
  
  N <- 100 # total number of sampled individuals
  Hstar <- 10 # total number of haplotypes
  probs <- rep(1/Hstar, Hstar) # equal haplotype frequency distribution
  
  # outputs a CSV file called "output.csv"
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, filename = "output") 
  
  ## Simulate hypothetical species - subsampling ##
  # subsamples 25\% of haplotype labels
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, 
  perms = 1000, p = 0.95, subsample = TRUE, prop = 0.25, 
  conf.level = 0.95, filename = "output") 
  
  ## Simulate hypothetical species and all paramaters changed - subsampling ##
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, 
  perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.95, 
  num.iters = 1, filename = "output")
  }
  