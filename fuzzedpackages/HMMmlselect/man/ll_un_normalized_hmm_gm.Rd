\name{ll_un_normalized_hmm_gm}
\alias{ll_un_normalized_hmm_gm}
\title{ll_un_normalized_hmm_gm}
\usage{
ll_un_normalized_hmm_gm(paras, yobs, bool_hmm, priors)
}
\arguments{
  \item{paras}{The parameter for the HMM.}
  \item{yobs}{The observed data.}
  \item{bool_hmm}{1 if it's HMM. 0 if it's Gaussian mixture.}
  \item{priors}{Prior distribution.}
}
\value{
It returns the unnormalized posterior density of HMM and Gaussian mixture.
}
\description{
Auxiliary function that computes the unnormalized posterior density of HMM and Gaussian mixture.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

