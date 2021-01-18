\name{CorMat_mcmc_summary}
\alias{CorMat_mcmc_summary}
\title{
Compute quantiles for the correlation matrix. 
}
\description{
Compute quantiles for the correlation matrix per cluster based on the MCMC output stored in a \code{fabMix.object}.
}
\usage{
CorMat_mcmc_summary(x, quantile_probs)
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
	\item{quantiles }{A list containing the quantiles for the correlation matrix per component. Each element is a \eqn{p\times p \times K} array, where \eqn{p} and \eqn{K} denote the dimension of the multivariate data and number of alive components for the selected model, respectively.}
	\item{p_matrix }{ A \eqn{p\times p\times K} array, where for each \eqn{k = 1,\ldots,K} the \eqn{p\times p} matrix \code{p_matrix[,,k]} contains the posterior probability \eqn{1-2|0.5-P(\rho_{ij} > 0)|} for \eqn{i=1,\ldots,p};  \eqn{j=1,\ldots,p}.}

}
\author{
Panagiotis Papastamoulis
}

