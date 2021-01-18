\name{shannonEntropy}
\alias{shannonEntropy}
\docType{package}
\title{
Shannon Entropy
}
\description{
Calculates the Shannon entropy over data of a discrete random variable \eqn{X}.
}
\arguments{
  \item{x}{A discrete vector.}
}
\usage{
shannonEntropy(x)
}
\examples{
set.seed(1)
x = round((runif(1000, min=0, max=5)))
shannonEntropy(x=x)	## 2.522265
}