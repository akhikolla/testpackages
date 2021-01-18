\name{slpm_gof}
\alias{slpm_gof}
\title{slpm_gof}
\description{
Evaluates the expected adjacency matrix for a fitted SparseLPM.
}
\usage{
slpm_gof(var_pars)
}
\arguments{
  \item{var_pars}{A list defining the variational parameters of the model. See \emph{Details} for more specific indications.}
}
\details{
The list \code{var_pars} must contain:
\describe{
 \item{alpha_u_tilde}{\code{M*K} matrix denoting the Gaussian means for senders.}
 \item{alpha_v_tilde}{\code{N*K} matrix denoting the Gaussian means for receivers.}
 \item{beta_u_tilde}{\code{M*K} matrix denoting the Gaussian variances for senders.}
 \item{beta_v_tilde}{\code{N*K} matrix denoting the Gaussian variances for receivers.}
 \item{lambda_tilde}{\code{M*N*K} array representing the soft clustering for the edges. This may be interpreted as the posterior probability that edge \code{ij} is determined by the \code{k}-th latent dimension.}
 \item{delta_tilde}{\code{K} dimensional vector containing the variational parameters for the mixing proportions. This may be interpreted as the importance of each of the latent dimensions.}
 \item{a_tilde}{\code{K} dimensional vector containing the shapes of the variational Gamma distributions associated to the precisions.}
 \item{b_tilde}{\code{K} dimensional vector containing the rates of the variational Gamma distributions associated to the precisions.}
Note that this function only uses the alphas and the lambdas. Also, to avoid numerical instability, the lambdas are automatically pre-transformed into a hard partitioning using a Maximum A Posterior method.
}
}
\value{An adjacency matrix with non-negative entries.}
\examples{
set.seed(12345)
M <- N <- 10
K <- 2
fitted_var_pars <- list()
fitted_var_pars$alpha_u_tilde = matrix(rnorm(M*K),M,K)
fitted_var_pars$alpha_v_tilde = matrix(rnorm(N*K),N,K)
fitted_var_pars$lambda_tilde = array(NA,c(M,N,K))
fitted_var_pars$lambda_tilde[,,1] = matrix(runif(M*N),M,N)
fitted_var_pars$lambda_tilde[,,2] = 1-fitted_var_pars$lambda_tilde[,,1]
expected_adj <- slpm_gof(fitted_var_pars)
}


