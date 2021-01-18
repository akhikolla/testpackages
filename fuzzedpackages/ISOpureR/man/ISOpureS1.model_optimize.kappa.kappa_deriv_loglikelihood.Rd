\name{ISOpureS1.model_optimize.kappa.kappa_deriv_loglikelihood}
\alias{ISOpureS1.model_optimize.kappa.kappa_deriv_loglikelihood}
\title{Compute derivative of loglikelihood with respect to kappa for step 1}
\description{Computes the derivative of the part of the loglikelihood function relevant to 
	optimizing kappa for step 1. Instead of performing constrained optimization on kappa directly, 
	we optimize the log of kappa in an unconstrained fashion.  Thus, if y=log(kappa) and L is the 
	loglikelihood function w.r.t. y, to optimize L w.r.t. y,   dL/dy = dL/dkappa * dkappa/dy, 
	where dkappa/dy = exp(y)= exp( log(kappa)). The input into the derivative function is 
	log(kappa - model\$MIN_KAPPA). }
\usage{
ISOpureS1.model_optimize.kappa.kappa_deriv_loglikelihood(log_kappa, tumordata, model)
}
\arguments{
	\item{log_kappa}{the scalar log(kappa - model\$MIN_KAPPA)}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumour samples}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative derivative of the part of the loglikelihood function relevant to kappa
	with respect to log kappa (a scalar given that for step 1 of ISOpure kappa is a scalar)}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
