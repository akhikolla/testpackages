\name{ISOpureS1.model_optimize.kappa.kappa_loglikelihood}
\alias{ISOpureS1.model_optimize.kappa.kappa_loglikelihood}
\title{Compute loglikelihood relevant to kappa for step 1}
\description{Computes the part of the loglikelihood function relevant to optimizing kappa for step 1. 
	Instead of performing constrained optimization on kappa directly, we optimize the log of kappa
   	in an unconstrained fashion.}
\usage{
ISOpureS1.model_optimize.kappa.kappa_loglikelihood(log_kappa, tumordata, model)
}
\arguments{
	\item{log_kappa}{the scalar log(kappa - model\$MIN_KAPPA)}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumour samples}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative of the loglikelihood relevant to optimizing kappa}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
