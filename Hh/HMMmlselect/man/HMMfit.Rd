\name{HMMfit}
\alias{HMMfit}
\title{HMMfit}
\usage{
HMMfit(y, K, METHOD, optionalfit = list())
}
\arguments{
  \item{y}{The observed data.}
  \item{K}{The specified number of states of the underlying Markov chian.}
  \item{METHOD}{Integer value indicating the method of parameter estimation: (a) HMM fitting through the Expectation-Maximization al-
gorithm (METHOD = 1), (b) HMM fitting through the Markov chain Monte Carlo algorithm
(METHOD = 2), and (c) Gaussian mixture model fitting through the Markov chain Monte
Carlo algorithm (METHOD = 3)}
  \item{optionalfit}{Optional variables as a list. Possible options include:}
  \itemize{
  \item{Ngibbs: Number of samples when using MCMC. Default is 5000.}
  \item{Burnin: Length of burnin period when using MCMC. Default is 5000.}
  \item{Thin: Thinning parameter when using MCMC. Default is 10.}
  \item{Nstart: Number of starting points. Default is 50.}
  \item{verbose: Logic variable indicating pritting details or not. Default is \code{FALSE}.}
  \item{priors: Prior when using MCMC. Default is flat.}
  }
}
\value{
This functions returns the fitting parameters of the observed data given the specified number of states.
}
\description{
The following function performs (a) HMM fitting through the Expectation-Maximization al-
gorithm (METHOD = 1), (b) HMM fitting through the Markov chain Monte Carlo algorithm
(METHOD = 2), and (c) Gaussian mixture model fitting through the Markov chain Monte
Carlo algorithm (METHOD = 3).
}
\details{
See Manual.pdf in "inst/extdata" folder.
}
\examples{
library(HMMmlselect)

# Example 1: use HMMfit to inference number of states
obs = HMMsim ( n = 200 )$obs
Nest = HMMfit( y = obs, K=3, METHOD = 1)
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

