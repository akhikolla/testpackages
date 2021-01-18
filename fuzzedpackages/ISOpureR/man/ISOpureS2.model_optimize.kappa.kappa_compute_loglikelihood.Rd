\name{ISOpureS2.model_optimize.kappa.kappa_compute_loglikelihood}
\alias{ISOpureS2.model_optimize.kappa.kappa_compute_loglikelihood}
\title{Compute loglikelihood relevant to kappa for step 2}
\description{Computes the part of the loglikelihood function relevant to optimizing kappa for step 2}
\usage{
ISOpureS2.model_optimize.kappa.kappa_compute_loglikelihood(kappa, model)
}
\arguments{
	\item{kappa}{a 1xK vector strength parameter in the prior over cc given the cancer profile mm}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The part of the loglikelihood function relevant to optimizing kappa}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
