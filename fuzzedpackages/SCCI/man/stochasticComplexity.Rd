\name{stochasticComplexity}
\alias{stochasticComplexity}
\docType{package}
\title{
Stochastic Complexity for Multinomials
}
\description{
Efficient implementation of the exact computation of stochastic complexity for multinomials (Silander et al. 2008) for data over a discrete random variable \eqn{X}.
}
\arguments{
  \item{x}{A discrete vector.}
}
\usage{
stochasticComplexity(x)	
}
\examples{
set.seed(1)
x = round((runif(1000, min=0, max=5)))
stochasticComplexity(x=x)	## 2544.698
}
\references{
  Tomi Silander, Teemu Roos, Petri Kontkanen and Petri Myllymäki; Factorized Normalized Maximum Likelihood Criterion for Learning Bayesian Network Structures, Proceedings of the 4th European Workshop on Probabilistic Graphical Models, 2008
}