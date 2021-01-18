\name{PlotHMM}
\alias{PlotHMM}
\title{PlotHMM}
\usage{
PlotHMM(y, results)
}
\arguments{
  \item{y}{The observed data.}
  \item{results}{The resulting output of state inference using \code{HMMmlselect}.}
}
\value{
It returns the graph with the original data and the inferenced states.
}
\description{
This function visualizes the state inference result using \code{HMMmlselect}. See the Manual.pdf under
data folder for a figure example.
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

