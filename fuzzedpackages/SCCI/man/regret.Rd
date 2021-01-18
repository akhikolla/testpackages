\name{regret}
\alias{regret}
\docType{package}
\title{
Multinomial Regret Term
}
\description{
Calculates the multinomial regret term for for a discrete random variable with domain size \eqn{k} and sample size \eqn{n} (see Silander et al. 2018). Note that we use the logarithm to basis 2 to calculate the result. To compare the results to Silander et al. (2018), we need to multiply the result with \eqn{log(2)} to compare the results.
}
\arguments{
  \item{n}{Integer (sample size)}
  \item{k}{Integer (domain size)}
}
\usage{
regret(n,k)	
}
\examples{
regret(50,10)           ## 19.1
regret(50,10) * log(2)  ## 13.24 (see Silander et al. 2018)
}
\references{
  Tomi Silander, Janne Leppä-aho, Elias Jääsaari, Teemu Roos; Quotient normalized maximum likelihood criterion for learning bayesian network structures, Proceedings of the 21nd International Conference on Artificial Intelligence and Statistics (AISTATS), PMLR, 2018
}