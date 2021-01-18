

\name{summary.fabMix.object}
\alias{summary.fabMix.object}
\title{
Summary method
}
\description{
S3 method for printing a summary of a \code{fabMix.object}.
}
\usage{
\method{summary}{fabMix.object}(object, quantile_probs, ...)
}
\arguments{
  \item{object}{
An object of class \code{fabMix.object}, which is returned by the \code{fabMix} function.
}
  \item{quantile_probs}{
A vector of quantiles to evaluate for each variable.
}
  \item{...}{Ignored.}
}
\value{
        A list consisting of the following entries:
  \item{alive_cluster_labels }{The labels of the ``alive'' components of the overfitting mixture model.}
  \item{posterior_means }{Posterior means of mixing proportion, marginal means and covariance matrix per (alive) cluster. }
  \item{quantiles }{A matrix containing the quantiles for each parameter.}
}

\details{
The function prints and returns a summary of the estimated posterior means for the parameters of the selected model for a \code{fabMix.object}. In particular, the method prints the ergodic means of the mixing proportions, marginal means and covariance matrix per component, as well as the corresponding quantiles.   
}
\note{
The \code{summary} function of the \code{coda} package to the \code{mcmc} object \code{object$mcmc} is used for computing quantiles.
}
\author{
Panagiotis Papastamoulis
}

