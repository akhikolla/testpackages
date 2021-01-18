\name{generate_Gaussian}
\alias{generate_Gaussian}

\title{
Generate a Gaussian distributed data set
}

\description{
This function will generate a Gaussian distributed data set with latent variables and correlated replicates.
}

\usage{
generate_Gaussian(n, R, p, l, s, sparsityA, sparsityobserved, sparsitylatent, lwb, 
upb, seed)
}

\arguments{
  \item{n}{the number of observations.}
  \item{R}{the number of replicates.}
  \item{p}{the number of observed variables.}
  \item{l}{the number of latent variables.}
  \item{s}{latent effects are generated as \eqn{s} piecewise constant across replicates. The number \eqn{s} should be a factor of \eqn{R}.}
  \item{sparsityA}{proportion of the number of zeros in the transition matrix \eqn{A}. Must be a number from 0 to 1.}
  \item{sparsityobserved}{proportion of the number of zeros in the inverse covariance of the observed variables. Must be a number from 0 to 1.}
  \item{sparsitylatent}{proportion of the number of zeros in the inverse covariances among latent variables and between observed variables and latent variables. Must be a number from 0 to 1.}
  \item{lwb}{lower bound for the elements in the inverse covariance matrix.}
  \item{upb}{upper bound for the elements in the inverse covariance matrix.}
  \item{seed}{the seed for the random number generator.}
}

\details{
This function aims to generate a Gaussian distributed data set with latent variables and correlated replicates. For each observation, the latent variables are piecewise constant across replicates, and conditioned on the latent variables, the replicates follow a one-lag vector autoregressive model, where the transition matrix \eqn{A} is sparse with non-zero elements set equal to 0.3. 
The matrix \eqn{\Sigma} is the covariance matrix for the observed variables X and the latent variables \eqn{U}, and we partition \eqn{\Sigma} into matrices that quantify the relationships among the observed variables (\eqn{\Sigma_{XX}}), between the observed variables and latent variables (\eqn{\Sigma_{XU}} or \eqn{\Sigma_{UX}}), and of the latent variables (\eqn{\Sigma_{UU}}).
In general, the data is generated with:
\deqn{
X_{i1} | U_{i1} \sim N_p(\Sigma_{XU}\Sigma^{-1}_{UU} U_{i1}, \Sigma_{XX} - \Sigma_{XU}\Sigma^{-1}_{UU}\Sigma_{UX}),
}

\deqn{
X_{it} | X_{i(t-1)},U_{it} \sim N_p(AX_{i(t-1)} + \Sigma_{XU}\Sigma^{-1}_{UU} U_{it}, \Sigma_{XX} - \Sigma_{XU}\Sigma^{-1}_{UU}\Sigma_{UX}),
}
where \eqn{1 \le i \le n} and \eqn{1 \le t \le R}. 
}

\value{
  \item{X}{the generated data, which is a list with \eqn{n} elements and each element is a matrix with \eqn{R} rows and \eqn{p} columns}
  \item{truegraph}{a matrix that encodes the conditional dependence relationships between sets of two observed variables}
}

\references{
Jin, Y., Ning, Y., and Tan, K. M. (2020), `Exponential Family Graphical Models with Correlated Replicates and Unmeasured Confounders', preprint available.
}

\examples{
data <- generate_Gaussian(n = 50, R = 20, p = 30, l = 2, s = 2, sparsityA = 0.95,
sparsityobserved = 0.9, sparsitylatent = 0.2, lwb = 0.3, upb = 0.3, seed = 1)
}
