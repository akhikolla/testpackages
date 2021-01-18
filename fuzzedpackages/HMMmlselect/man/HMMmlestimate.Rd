\name{HMMmlestimate}
\alias{HMMmlestimate}
\title{HMMmlestimate}
\usage{
HMMmlestimate(y, K, optionalfit = list())
}
\arguments{
  \item{y}{The data that the marginal likelihood is computed for.}
  \item{K}{Number of states.}
  \item{optionalfit}{Optional variables.}
}
\value{
It returns the approximated marginal likelihood.
}
\description{
Auxiliary function that approximates the marginal likelihood.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

