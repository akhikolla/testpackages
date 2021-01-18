\name{slpm_nga}
\alias{slpm_nga}
\title{slpm_nga}
\description{
Runs a Natural Gradient Ascent algorithm to maximise the variational objective for a Sparse LPM.
}
\usage{
slpm_nga(X, K, var_pars_init, hyper_pars = NULL, tol = 0.01, n_iter_max = 100000,
        natural_gradient = T, learning_rate_factor_up = 2, learning_rate_factor_down = 2,
        verbose = F)
}
\arguments{
  \item{X}{Rectangular adjacency matrix with non-negative entries.}
  \item{K}{The number of latent dimension of the model.}
  \item{var_pars_init}{List of variational parameters to be used as starting point for the optimisation. See \emph{Details} for more specific indications.}
  \item{hyper_pars}{List defining the hyperparameters of the model. The list should contain three vectors of \code{K} positive values denoted \code{delta}, \code{a_gamma} and \code{b_gamma}, respectively, where \code{K} is the number of latent dimensions. If left as null, all \code{delta} parameters are set to \code{0.001}, whereas the rest is set to \code{1}.}
  \item{tol}{Positive number setting the stop condition: the algorithm stops if one entire iteration yields an increase in the objective function smaller than this value.}
  \item{n_iter_max}{Maximum number of iterations the algorithm should be run for.}
  \item{natural_gradient}{\code{TRUE} or \code{FALSE} indicating whether the natural gradient instead of the standard gradient should be used.}
  \item{learning_rate_factor_up}{Before any natural gradient ascent update, the current step size is multiplied by this number to ensure that the algorithms tries new solutions which are relatively far from the current one.}
  \item{learning_rate_factor_down}{During any natural gradient ascent update, if a certain step size leads to a decrease in the objective function, then the step is divided by this number repeatedly until an increase is ensured.}
  \item{verbose}{\code{TRUE} or \code{FALSE} indicating whether a lengthy output should be printed out.}
}
\details{
\code{var_pars} and \code{var_pars_init} are lists with components:
\describe{
  \item{alpha_u_tilde}{\code{M*K} matrix representing the Gaussian means for the latent positions of senders.}
  \item{alpha_v_tilde}{\code{N*K} matrix representing the Gaussian means for the latent positions of receivers.}
  \item{beta_u_tilde}{\code{M*K} matrix representing the Gaussian variances for the latent positions of senders.}
  \item{beta_v_tilde}{\code{N*K} matrix representing the Gaussian variances for the latent positions of receivers.}
  \item{lambda_tilde}{\code{M*N*K} array with entries corresponding to the posterior probabilities of assigning each edge to each latent dimension.}
  \item{delta_tilde}{Vector of \code{K} positive values representing the Dirichlet parameters generating the mixing proportions.}
  \item{a_tilde}{Vector of \code{K} positive values corresponding to the shapes of the variational Gamma distribution on the precisions.}
  \item{b_tilde}{Vector of \code{K} positive values corresponding to the rates of the variational Gamma distribution on the precisions.}
}
}
\value{
A list with components:
  \item{computing_time}{Number of seconds required for the optimisation process.}
  \item{var_pars}{List containing the optimal values for the variational parameters.}
  \item{learning_rates_u}{Current step-size for the update of the variational parameters of each Gaussian distribution on the latent positions of senders.}
  \item{learning_rates_v}{Current step-size for the update of the variational parameters of each Gaussian distribution on the latent positions of receivers.}
  \item{elbo_values}{Values of the variational objective at the end of each of the iterations.}
  \item{elbo_init}{Value of the variational objective for the initial configuration.}
  \item{elbo_final}{Value of the variational objective for the optimal solution found.}
}
\references{
Rastelli, R. (2018) "The Sparse Latent Position Model for nonnegative weighted networks", \url{https://arxiv.org/abs/1808.09262}
}
\examples{
set.seed(12345)
network <- slpm_gen(M = 15, N = 10, K = 2)
K <- 6
var_pars_init <- slpm_init(X = network$adj, K = K)
res <- slpm_nga(X = network$adj, K = K, var_pars_init = var_pars_init)
}


