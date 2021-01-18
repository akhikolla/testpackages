\name{logmvTdensity}
\alias{logmvTdensity}
\title{logmvTdensity}
\usage{
logmvTdensity(x, mu, sqrt_inv_sigma, lgsqrt_det_sigma, df, d, logconstT)
}
\arguments{
  \item{x}{Data to compute the log likelihood.}
  \item{mu}{Mean of the t-distribution}
  \item{sqrt_inv_sigma}{Corresponds to the standard deviation of the t-distribution.}
  \item{lgsqrt_det_sigma}{Corresponds to the determinant of the variance-covariance matrix.}
  \item{df}{Degree of freedom of the t-distribution.}
  \item{d}{Dimension of the t-distribution.}
  \item{logconstT}{Leading constant for t-distribution density function.}
}
\value{
It returns the log-likelihood of the multivariate student-t distribution.
}
\description{
Auxiliary function that computes log-likelihood based on multivariate t-distribution. Used when approximating marginal likelihood.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

