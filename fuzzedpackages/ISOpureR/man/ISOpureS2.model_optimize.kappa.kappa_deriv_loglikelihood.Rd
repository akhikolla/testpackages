\name{ISOpureS2.model_optimize.kappa.kappa_deriv_loglikelihood}
\alias{ISOpureS2.model_optimize.kappa.kappa_deriv_loglikelihood}
\title{Compute derivative of loglikelihood with respect to kappa for step 2}
\description{Computes the derivative of the part of the loglikelihood function relevant to 
	optimizing kappa for step 2. Instead of performing constrained optimization on kappa 
	directly, we optimize the log of kappa in an unconstrained fashion.}
\usage{
ISOpureS2.model_optimize.kappa.kappa_deriv_loglikelihood(log_kappa, model)
}
\arguments{
	\item{log_kappa}{the 1xD matrix log(kappa - model\$MIN_KAPPA)}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative derivative of the part of the loglikelihood function relevant to kappa
	with respect to log kappa (a Dx1 matrix).}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
