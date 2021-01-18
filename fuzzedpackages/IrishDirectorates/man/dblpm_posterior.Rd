\name{dblpm_posterior}
\alias{dblpm_posterior}
\title{dblpm_posterior}
\description{
Evaluates the posterior value for a given Dynamic Bipartite Latent Position Model (dblpm network).
}
\usage{
dblpm_posterior(network)
}
\arguments{
  \item{network}{A list identifying a dblpm network.}
}
\value{
  \item{computing_time}{Number of seconds required for the evaluation.}
  \item{likelihood_value}{Likelihood value for the given network.}
  \item{posterior_value}{Posterior value for the given network.}
}
\examples{
data(IrishDirectoratesFit)
res <- dblpm_posterior(network = IrishDirectoratesFit$tail)
}


