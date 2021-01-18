\name{HAC.simrep}
\alias{HAC.simrep}

\title{Run a simulation of haplotype accumulation curves for hypothetical or real species}

\description{Runs the \code{HACSim} algorithm by successively calling \code{HAC.sim} to iteratively extrapolate haplotype accumulation curves to determine likely specimen sample sizes for hypothetical or real species

The algorithm employs the following iterative methods when calculating the "Measures of Sampling Closeness":

\itemize{
\item Mean number of haplotype sampled: \eqn{H_i}
\item Mean number of haplotypes not sampled \eqn{H^* - H_i}
\item Proportion of haplotypes sampled: \eqn{\frac{H_i}{H^*}}
\item Proportion of haplotypes not sampled: \eqn{1 - \frac{H_i}{H^*}}
\item Mean value of \eqn{N*}: \eqn{\frac{N_iH^*}{H_i}}
\item Mean number of specimens not sampled: \eqn{\frac{N_iH^*}{H_i} - N_i}
}

where \eqn{H_i} is stochastically-determined through sampling from \code{probs}, the observed species' haplotype frequency distribution vector. 

As the algorithm proceeds, \eqn{H_i} will approach \eqn{H^*} asymptotically (and hence, \eqn{N_i} will converge to \eqn{N^*}), but will likely fluctuate randomly from one iteration to the next. However, estimates of \eqn{N^*} found at each iteration will be monotonically-increasing.
}

\usage{HAC.simrep(HACSObject)}

\arguments{\item{HACSObject}{object containing the desired simulation parameters}
}

\value{Iteration results are outputted to the console and graphs displayed in the plot window. Plots depict haplotype accumulation (along with shaded confidence intervals for the mean number of haplotypes found). Dashed lines correspond to the endpoint of the curve and reflect haplotype recovery for a user-defined cutoff (default \code{p} = 0.95, 95\% haplotype diversity). Output from the first iteration is useful for judging levels of haplotype diversity and recovery found in observed intraspecific sequence datasets, reflecting current sampling depth. The required sample size is displayed in the second-last iteration. All other information corresponding to the extrapolated sample size can be found in the last iteration. 
Iteration results can optionally be saved to a CSV file. Subsampled DNA sequences are automatically saved to a FASTA file.}

\note{When simulating real species via \code{HACReal(...)}, a pop-up window will appear prompting the user to select an intraspecific FASTA file of aligned/trimmed DNA sequences. The alignment must not contain missing or ambiguous nucleotides (i.e., it should only contain A, C, G or T); otherwise, haplotype diversity may be overestimated. Excluding sequences or alignment sites with missing/ambiguous data is an option.
}

\examples{

  ## Simulate hypothetical species ##
  
  N <- 100 # total number of sampled individuals
  Hstar <- 10 # total number of haplotypes
  probs <- rep(1/Hstar, Hstar) # equal haplotype frequency distribution
  
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar , probs = probs, 
  filename = "output") # outputs a CSV file called "output.csv"
  
  ## Simulate hypothetical species - subsampling ##
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, 
  perms = 1000, p = 0.95, subsample = TRUE, prop = 0.25, 
  conf.level = 0.95, filename = "output")
  
  ## Simulate hypothetical species and all paramaters changed - subsampling ##
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, 
  perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.95, 
  filename = "output")
  
  HAC.simrep(HACSObj) # runs a simulation
  
  ## Simulate real species ##
  
  \donttest{
    ## Simulate real species ##
    # outputs file called "output.csv"
    HACSObj <- HACReal(filename = "output") 
    
    ## Simulate real species - subsampling ##
    HACSObj <- HACReal(subsample = TRUE, prop = 0.15, conf.level = 0.95, 
    filename = "output")
    
    ## Simulate real species and all parameters changed - subsampling ##
    HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, 
    prop = 0.15, conf.level = 0.99, filename = "output")
    
    # user prompted to select appropriate FASTA file
    HAC.simrep(HACSObj) 
    }
}