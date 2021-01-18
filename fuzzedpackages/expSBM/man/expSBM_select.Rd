\name{expSBM_select}
\alias{expSBM_select}
\title{expSBM_select}
\description{
Runs the variational expectation maximization algorithm for different numbers of latent groups, and selects the best overall model using the integrated completed likelihood criterion. See reference for a detailed explanation of the procedure.
}
\usage{
expSBM_select(K_max, N, edgelist, method = "SBM_gaussian", directed = F, 
              trunc = T, tol = 0.001, n_iter_max = 100, 
              init_blur_value = 1, verbose = F)
}
\arguments{
  \item{K_max}{Estimate and compare the models with number of latent groups equal to 1,2,...,\code{K_max}.}
  \item{N}{Number of nodes.}
  \item{edgelist}{A matrix with 4 columns: on the first column the sender node, on the second the receiver, on the third either a one or zero to indicate whether it is an interaction or a non-interaction respectively, on the fourth the corresponding exponential length.}
  \item{method}{Indicates the method used for the initialisation. Can be one of \code{random}, \code{SBM_binary}, \code{SBM_gaussian} or \code{spectral}. See \code{\link{expSBM_init}} for more details.}
  \item{directed}{\code{TRUE} or \code{FALSE} indicating whether interactions have an orientation or not.}
  \item{trunc}{\code{TRUE} or \code{FALSE} indicating whether the first and last interactions or non-interactions for every edge are assumed to be truncated or not.}
  \item{tol}{Stop the maximization if the relative increase in the objective function is not larger than this value.}
  \item{n_iter_max}{Stop the maximization if the number of iterations is larger than this value. This parameter can be set to zero or one for debug purposes.}
  \item{init_blur_value}{A value from zero to one, indicating if the initialized partition should be perturbed with noise. The value one means no noise, whereas the value zero has maximum noise, i.e. each node is equally likely belonging to any of the K groups.}
  \item{verbose}{\code{TRUE} or \code{FALSE} indicating whether a lengthy output should be printed out.}
}
\value{
  \item{fitted_models}{A list with the fitted values for every model considered.}
  \item{icl_values}{Integrated completed likelihood values for each model considered.}
  \item{K_star}{Optimal number of latent groups, according to the integrated completed likelihood criterion.}
  \item{best_model}{Output of the variational expectation maximization algorithm for the best overall model.}
}
\references{
R. Rastelli and M. Fop (2019) "A dynamic stochastic blockmodel for interaction lengths", \url{https://arxiv.org/abs/1901.09828}
}
\examples{
set.seed(1)
data(high_school)
\donttest{
res <- expSBM_select(K_max = 8, N = 327, edgelist = high_school$edgelist, 
     method = "random", tol = 0.01)
}
}
