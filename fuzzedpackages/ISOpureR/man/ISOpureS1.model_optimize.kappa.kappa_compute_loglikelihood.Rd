\name{ISOpureS1.model_optimize.kappa.kappa_compute_loglikelihood}
\alias{ISOpureS1.model_optimize.kappa.kappa_compute_loglikelihood}
\title{Compute loglikelihood relevant to kappa for step 1}
\description{Computes the part of the loglikelihood function relevant to optimizing kappa for step 1}
\usage{
ISOpureS1.model_optimize.kappa.kappa_compute_loglikelihood(kappa, tumordata, model) 
}
\arguments{
	\item{kappa}{a scalar kappa, the strength parameter in the prior over the reference cancer profile}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumour samples}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The part of the loglikelihood function relevant to optimizing kappa}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
