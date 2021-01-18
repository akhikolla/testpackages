\name{slpm_elbo}
\alias{slpm_elbo}
\title{slpm_elbo}
\description{
Evaluates the evidence lower bound for a given configuration of variational parameters.
}
\usage{
slpm_elbo(X, var_pars, hyper_pars, verbose = F)
}
\arguments{
  \item{X}{Rectangular adjacency matrix with non-negative entries.}
  \item{var_pars}{A list defining the variational parameters of the model. See \emph{Details} for more specific indications.}
  \item{hyper_pars}{A list defining the hyperparameters of the model. The list should contain three vectors of length \code{K} denoted \code{delta}, \code{a_gamma} and \code{b_gamma}, where \code{K} is the number of latent dimensions.}
  \item{verbose}{\code{TRUE} or \code{FALSE} indicating whether a lengthy output should be printed out.}
}
\details{
The list \code{var_pars} must contain:
\describe{
 \item{alpha_u_tilde}{\code{M*K} matrix denoting the Gaussian means for senders.}
 \item{alpha_v_tilde}{\code{N*K} matrix denoting the Gaussian means for receivers.}
 \item{beta_u_tilde}{\code{M*K} matrix denoting the Gaussian variances for senders.}
 \item{beta_v_tilde}{\code{N*K} matrix denoting the Gaussian variances for receivers.}
 \item{lambda_tilde}{\code{M*N*K} array representing the soft clustering for the edges. This may be interpreted as the posterior probability that edge \code{ij} is determined by the \code{k}-th latent dimension.}
 \item{delta_tilde}{\code{K} dimensional vector containing the variational parameters for the mixing proportions. This may be interpreted as the importance of each latent dimension.}
 \item{a_tilde}{\code{K} dimensional vector containing the shapes of the variational Gamma distributions associated to the precisions.}
 \item{b_tilde}{\code{K} dimensional vector containing the rates of the variational Gamma distributions associated to the precisions.}
}
}
\value{
  \item{computing_time}{Number of seconds required for the evaluation.}
  \item{elbo}{Value of the ELBO for the given variational parameters.}
}
\examples{
set.seed(12345)
M <- N <- 10
K <- 2
network <- slpm_gen(M = M, N = N, K = K)
var_pars <- slpm_init(X = network$adj, K = K)
hyper_pars <- list(delta = rep(1,K), a_gamma = rep(1,K), b_gamma = rep(1,K))
slpm_elbo(X = network$adj, var_pars = var_pars, hyper_pars = hyper_pars, verbose = FALSE)
}


