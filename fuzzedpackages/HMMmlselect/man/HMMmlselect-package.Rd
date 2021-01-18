\name{HMMmlselect-package}
\alias{HMMmlselect-package}
\docType{package}
\title{Determine the number of states
in hidden Markov models via marginal likelihood}
\description{
(See the Manual.pdf file in "inst/extdata"" folder for a detail description of all functions, and a detail tutorial.)

This packages provides function to estimate the number of states in a Hidden Markov model (HMM) using the  marginal likelihood method proposed by Chen, Fuh, Kao and Kou (2019+), which we would denoted as HMMml2017 afterward. \code{HMMmlselect} estimates the number of states, and \code{PlotHMM} plots. Other related functions are also provided.
}



\author{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou

Maintainer: Chu-Lan Michael Kao <chulankao@gmail.com>
}
\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
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

\details{
\tabular{ll}{
Package: \tab HMMmlselect\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2019-2-08\cr
License: \tab GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007\cr
}

}
