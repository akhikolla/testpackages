\name{update_SigmaINV_faster_q0}
\alias{update_SigmaINV_faster_q0}
\title{
Gibbs sampling for \eqn{\Sigma^{-1}} per component for \eqn{q=0}
}
\description{
Gibbs sampling for \eqn{\Sigma^{-1}} per component
}
\usage{
update_SigmaINV_faster_q0( z, mu, K, alpha_sigma, beta_sigma, x_data)
}
\arguments{
  \item{z}{
	Allocation vector
}
  \item{mu}{
	Marginal means
}
  \item{K}{
	Number of components
}
  \item{alpha_sigma}{
	Prior parameter
}
  \item{beta_sigma}{
	Prior parameter
}
  \item{x_data}{
	Data
}
}
\value{
	\eqn{\Sigma^{-1}}
}
\author{
Panagiotis Papastamoulis	
}

