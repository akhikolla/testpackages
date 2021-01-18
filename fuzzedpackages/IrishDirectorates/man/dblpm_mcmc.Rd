\name{dblpm_mcmc}
\alias{dblpm_mcmc}
\title{dblpm_mcmc}
\description{
Runs the Metropolis-within-Gibbs sampler on the given Dynamic Bipartite Latent Position Model (dblpm network).
}
\usage{
dblpm_mcmc(network, niter, burnin, thin, x_var, w_var, gamma_var, beta_var, verbose = T)
}
\arguments{
  \item{network}{A list identifying a dblpm network.}
  \item{niter}{Number of iterations after thinning and burnin.}
  \item{burnin}{Number of iterations to be discaded before starting the count for niter. The burnin iterations are not thinned.}
  \item{thin}{After burnin, keep one sampled observation every \code{thin} and discard the rest.}
  \item{x_var}{Proposal variance for the positions of sender nodes.}
  \item{w_var}{Proposal variance for the positions of receiver nodes.}
  \item{gamma_var}{Proposal variance for the intercept gamma.}
  \item{beta_var}{Proposal variance for the intercept beta.}
  \item{verbose}{\code{true} or \code{false} indicating whether a lengthy output should be printed out.}
}
\value{
  \item{computing_time}{Number of seconds required for the sampling process.}
  \item{samples}{Sampled values for each of the model parameters.}
  \item{tail}{dblpm network sampled in the last iteration.}
}
\examples{
data(IrishDirectoratesFit)
IrishDirectoratesFit <- dblpm_mcmc(network = IrishDirectoratesFit$tail, 
                  niter = 3, burnin = 6, thin = 3, 
                  x_var = 4.75, w_var = 0.25, gamma_var = 1.825, beta_var = 0.2175, 
                  verbose = TRUE)
# to replicate the results of the paper: niter = 10000, burnin = 500000, thin = 50
}
