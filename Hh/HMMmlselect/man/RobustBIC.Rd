\name{RobustBIC}
\alias{RobustBIC}
\title{RobustBIC}
\usage{
RobustBIC(y, optionalbic = list())
}
\arguments{
  \item{y}{The observed data.}
  \item{optionalbic}{Optional variables as a list. Possible options include:}
  \itemize{
  \item{Kfits: Possible numbers of states. The function will compute the marginal likelihood under each of these numbers, and choose the one with the highest values as the estimated number of states. Default is \{2,3,...,6\}}
  \item{Nstart: The number of starting points for the robust BIC. Default is 50.}
  \item{verbose: Logic variable indicating pritting details or not. Default is \code{FALSE}.}
  }
}
\value{
This function returns the estimated number of hidden states through minimizing the BIC,
the BIC values of all the possible number of hidden states, and the fitted model parameters
under the estimated number of hidden states under the BIC method.
}
\description{
This function estimates the number of states of the given HMM data using robust BIC criteria.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}
\examples{
library(HMMmlselect)

# Example 1: use robust BIC to determine the order of HMM
obs = HMMsim ( n = 200 )$obs
resultsBIC = RobustBIC ( y = obs )
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2017+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

