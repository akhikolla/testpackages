\name{HMMmlselect}
\alias{HMMmlselect}
\title{HMMmlselect}
\usage{
HMMmlselect(y, optionalfit = list())
}
\arguments{
  \item{y}{The observed data.}
  \item{optionalfit}{Optional variables as a list. Possible options include:}
  \itemize{
  \item{Kfits: Possible numbers of states. The function will compute the marginal likelihood under each of these numbers, and choose the one with the highest values as the estimated number of states. Default is \{2,3,...,6\}}
  \item{RIS: Logical variable indicating whether to use the reciprocal importance sampling method in HMMml2017. Default is set to be \code{FALSE}.}
  \item{IS: Logical variable indicating whether to use the importance sampling method in HMMml2017. Default is set to be \code{TRUE}.}
  \item{RunParallel: Logical variable indicating whether to run using parallelization or not. Default is \code{TRUE}.}
  \item{boolUseMclust: Logical variable indicating whether to use the \code{Mclust} package. Default is set to be \code{FALSE}.}
  \item{priors: Lists of hyper parameters for the prior. Default is flat prior.}
  \item{boolHMM: Compute using marginal likelihood of HMM. Default is \code{TRUE}. If it is set to \code{FALSE}, then the program will use the marignal likelihood of mixture model instead.}
}}
\value{
It returns (1) the estimated number of hidden states using the marginal likelihood method,
(2) the marginal likelihood values corresponding to 2, 3, ... number of hidden states, and (3)
the fitted model parameters given the estimated number of hidden states.
}
\description{
This function computes the marginal likelihood of the HMM model
with the observed data and various number of states, and choose the one
with the highest marginal likelihood as the estimated number of states.
The method in Chen et al. (2017) is used, which we will denote it as
HMMml2017 afterward.
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

