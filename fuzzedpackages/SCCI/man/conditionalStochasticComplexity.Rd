\name{conditionalStochasticComplexity}
\alias{conditionalStochasticComplexity}
\docType{package}
\title{
Conditional Stochastic Complexity for Multinomials
}
\description{
Calculates the Stochastic Complexity of a discrete random variable \eqn{X} conditioned on a discrete (possibly multivariate) random variable \eqn{Y}. Variants for both factorized NML (fNML, Silander et al. 2008) and quotient NML (qNML, Silander et al. 2018) are included.
}
\arguments{
  \item{x}{A discrete vector.}
  \item{y}{A discrete vector or a data frame containing multiple discrete vectors to condition X on.}
  \item{score}{Default: fNML, optionally qNML can be passed.}
}
\usage{
conditionalStochasticComplexity(x, y, score="fNML")	
}
\examples{
set.seed(1)
x = round((runif(1000, min=0, max=5)))
Y = data.frame(round((runif(1000, min=0, max=5))), round((runif(1000, min=0, max=5))))
conditionalStochasticComplexity(x=x,y=Y,score="fNML")	## 2779.477
}
\references{
  Tomi Silander, Janne Leppä-aho, Elias Jääsaari, Teemu Roos; Quotient normalized maximum likelihood criterion for learning bayesian network structures, Proceedings of the 21nd International Conference on Artificial Intelligence and Statistics (AISTATS), PMLR, 2018
  
    Tomi Silander, Teemu Roos, Petri Kontkanen and Petri Myllymäki; Factorized Normalized Maximum Likelihood Criterion for Learning Bayesian Network Structures, Proceedings of the 4th European Workshop on Probabilistic Graphical Models, 2008
}