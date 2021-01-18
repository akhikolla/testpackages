\name{find_importance_function}
\alias{find_importance_function}
\title{find_importance_function}
\usage{
find_importance_function(x, boolUseMclust)
}
\arguments{
  \item{x}{The data that the importance weight is computed on.}
  \item{boolUseMclust}{Auxiliary variable. Set to \code{TRUE}.}
}
\value{
It returns the importance weight used in approximated normalizing constant.
}
\description{
Auxiliary function that computes the importance weight used in approximate the normalizing constant.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

