\name{conditionalShannonEntropy}
\alias{conditionalShannonEntropy}
\docType{package}
\title{
Conditional Shannon Entropy
}
\description{
Calculates the Shannon entropy of a discrete random variable \eqn{X} conditioned on a discrete (possibly multivariate) random variable \eqn{Y}.
}
\arguments{
  \item{x}{A discrete vector.}
  \item{y}{A data frame.}
}
\usage{
conditionalShannonEntropy(x, y)	
}
\examples{
set.seed(1)
x = round((runif(1000, min=0, max=5)))
Y = data.frame(round((runif(1000, min=0, max=5))), round((runif(1000, min=0, max=5))))
conditionalShannonEntropy(x=x,y=Y) ## 2.411972
}