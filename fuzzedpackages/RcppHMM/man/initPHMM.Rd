\name{initPHMM}
\alias{initPHMM}

\title{
  Random Initialization for a Hidden Markov Model with emissions modeled as discrete variables
}

\description{
Function used to generate a hidden Markov model with discrete observations and random parameters. This model is used when the observed data are counts that can be modelled with a mixture of Poissons. The code for the methods with categorical values or continuous data can be viewed in \code{"\link{initHMM}"} and \code{"\link{initGHMM}"}, respectively.  
}

\usage{
initPHMM(n)
}

\arguments{
  \item{n}{the number of hidden states to use.}
}

\value{
  A \code{"\link{list}"} that contains all the required values to specify the model.
  \item{Model}{it specifies that the observed values are to be modeled as a Poisson mixture model.}
  \item{StateNames}{the set of hidden state names.}
  \item{A}{the transition probabilities matrix.} 
  \item{B}{a vector with the lambda parameter for each Poisson distribution.} 
  \item{Pi}{the initial probability vector.}
}

\references{
  Cited references are listed on the \link{RcppHMM} manual page.
}

\examples{
  n <- 2
  model <- initPHMM(n)
  print(model)
}

\keyword{initialization}
