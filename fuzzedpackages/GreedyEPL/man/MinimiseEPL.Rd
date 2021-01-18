% File src/library/base/man/load.Rd
\name{MinimiseEPL}
\alias{MinimiseEPL}
\title{MinimiseEPL}
\description{
  Performs greedy optimisation to find the Bayes optimal decision given a collection (sample) of partitions.
}
\usage{
MinimiseEPL(sample_of_partitions, pars = list())
}
\arguments{
  \item{sample_of_partitions}{A \code{TxN} matrix whose rows correspond to partitions of the same \code{N} data observations.}
  \item{pars}{A number of additional parameters. See 'Details'.}
}
\details{
\itemize{
\item{\code{weights}: }{a vector of size \code{T} assigning a positive weight to each of the partitions in the sample. The vector does not have to be normalised but its entries must be positive numbers. If unspecified all weights are set to 1.}
\item{\code{Kup}: }{positive integer identifying the number of groups used to generate the random initial partition. Defaults to \code{N}.}
\item{\code{decision_init}: }{a vector of size \code{N} identifying the partition used to initialise the greedy algorithm. Entries must be integers greater than one. The default is \code{sample(1:Kup, N, TRUE)}. Note that the number of groups in the initial partition also corresponds to the maximum number of groups allowed.}
\item{\code{loss_type}: }{defines the loss function to be used in the expected posterior loss minimisation. Can be one of \code{"VI"} (Variation of Information), \code{"B"} (Binder's loss), \code{"NVI"} (Normalised Variation of Information) or \code{"NID"} (Normalised Information Distance). Defaults to \code{"VI"}.}
}
}
\seealso{
  \code{\link{CollapseLabels}}.
}
\references{
Rastelli, R. and Friel, N. (2016) "Optimal Bayesian estimators for latent variable cluster models", \url{https://arxiv.org/abs/1607.02325}
}
\examples{
set.seed(123)
N <- 10
niter <- 50
Kup <- 3
sample_of_partitions <- matrix(sample(x = 1:Kup, size = niter*N, replace = TRUE),niter,N)
loss_type = "VI"
output <- MinimiseEPL(sample_of_partitions, list(Kup = Kup, loss_type = loss_type))
output$decision
output$EPL
}



