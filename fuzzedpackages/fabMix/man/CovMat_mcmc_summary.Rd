\name{CovMat_mcmc_summary}
\alias{CovMat_mcmc_summary}
\title{
Compute quantiles for the covariance matrix. 
}
\description{
Compute quantiles for the covariance matrix per cluster based on the MCMC output stored in a \code{fabMix.object}.
}
\usage{
CovMat_mcmc_summary(x, quantile_probs)
}
\arguments{
  \item{x}{
	An object of class \code{fabMix.object}.
}
  \item{quantile_probs}{
	Vector of probabilities for computing the corresponding quantiles.
}
}
\value{
	A list containing the quantiles for the covariance matrix per component. Each element is a \eqn{p\times p \times K} array, where \eqn{p} and \eqn{K} denote the dimension of the multivariate data and number of alive components for the selected model, respectively.

}
\author{
Panagiotis Papastamoulis
}

