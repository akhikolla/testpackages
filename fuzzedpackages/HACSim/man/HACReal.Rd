\name{HACReal}
\alias{HACReal}

\title{Set up an object to simulate haplotype accumulation curves for a real species}

\description{Helper function which creates an object containing necessary information to run a simulation of haplotype accumulation for a real species of interest}

\usage{HACReal(perms = 10000, p = 0.95, conf.level = 0.95, 
subsample = FALSE, prop = NULL, progress = TRUE, num.iters = NULL, filename = NULL)}

\arguments{\item{perms}{Number of permutations (replications)}
\item{p}{Proportion of haplotypes to recover}
\item{conf.level}{Desired confidence level for graphical output and interval estimation}
\item{subsample}{Is a subsample of DNA sequences desired?}
\item{prop}{If subsample = TRUE, the proportion of DNA sequences to subsample}
\item{num.iters}{Number of iterations to compute}
\item{progress}{Should iteration output be printed to the R console?}
\item{filename}{Name of file where simulation results are to be saved}
}

\value{A list object of class "HAC" with 13 elements that can be passed to \code{HAC.simrep} as follows:

  \item{input.seqs}{Should a FASTA file of aligned/trimmed DNA sequences be inputted? Default is TRUE}
  \item{subset.seqs}{Should a subsample of DNA sequences be taken? Default is FALSE (can be altered by user)}
  \item{prop.seqs}{Proportion of DNA sequences to subsample. Default is NA (can be altered by user)}
  \item{prop.haps}{Proportion of haplotype labels to subsample. Default is NULL}
  \item{subset.haps}{Should a subsample of haplotype labels be taken? Default is NULL}
  \item{N}{Number of individuals. NA by default (computed automatically by algorithm)}
  \item{Hstar}{Number of unique species' haplotypes. NA by default (computed automatically by algorithm)}
  \item{probs}{Haplotype frequency distribution vector. NA by default (computed automatically by algorithm)}
  \item{p}{Proportion of haplotypes to recover. \code{p} = 0.95 by default.}
  \item{perms}{Number of permutations (replications). \code{perms} = 10000 by default.}
  \item{conf.level}{Desired confidence level for graphical output and interval estimation. \code{conf.level} = 0.95 by default.}
    \item{num.iters}{Number of iterations to compute. \code{num.iters = NULL} by default (i.e., all iterations are computed; users can specify \code{num.iters} = 1 for the first iteration.)}
  \item{progress}{Should iteration output be printed to the R console? Default is TRUE}
  \item{filename}{Name of file where simulation results are to be saved.}

}

\examples{
    ## Simulate real species ##
    # outputs file called "output.csv"
    HACSObj <- HACReal(filename = "output") 
    
    ## Simulate real species - subsampling ##
    # subsamples 25\% of DNA sequences
    HACSObj <- HACReal(subsample = TRUE, prop = 0.25, conf.level = 0.95, 
    filename = "output") 
    
    ## Simulate real species and all parameters changed - subsampling ##
    HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, 
    prop = 0.15, conf.level = 0.99, num.iters = 1, filename = "output")
}
