\name{HMMsim}
\alias{HMMsim}
\title{HMMsim}
\usage{
HMMsim(n, optionalsim = list())
}
\arguments{
  \item{n}{The length of HMM to be simulated.}
  \item{optionalsim}{Optional variables as a list. Possible options include:}
  \itemize{
  \item{Ksim: Number of states of the simulated HMM. Default is \code{3}.}
  \item{P: The transition matrix of the underlying Markov chain. Default is a flat K-by-K matrix.}
  \item{mu: The mean of the observed data given each underlying state. Default is \{1, 2, ..., K\}.}
  \item{sigma: The standard deviation of the observed data given each underlying state. Default is \{0.1, 0.1, ... 0.1\}.}
  \item{pi: The distribution of the initial state. Default is an uniform distribution across all possible states.}
  \item{BoolWritetoFile: Logic variable indicating whether to write the result into file or not. Default is \code{FALSE}.}
  \item{Filenameoutput: The output file name for the simulated HMM. Default is "HMMtrace.txt".}
  }
}
\value{
It returns the sample of the simulated HMM.
}
\description{
This function simulates HMM with the observed data being conditionally Gaussian distributed given the underlying state.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}
\examples{
library(HMMmlselect)

# simulate a 25 observations HMM
obs = HMMsim ( n = 25 )$obs

# perform order selection and estimation
results = HMMmlselect ( y = obs, list(Kfits = c(2,3), boolUseMclust = FALSE) )

# visualize the results, see figure 1
PlotHMM ( y = obs, results )
}


\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

