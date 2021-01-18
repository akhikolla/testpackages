\name{ISOpureS2.model_optimize.cc.cc_deriv_loglikelihood}
\alias{ISOpureS2.model_optimize.cc.cc_deriv_loglikelihood}
\title{Compute the derivative of loglikelihood relevant to the patient cancer profiles, cc, for step 2}
\description{Computes the derivative of the part of the likelihood function relevant to optimizing 
	cc.}
\usage{
ISOpureS2.model_optimize.cc.cc_deriv_loglikelihood(ww, tumordata, dd, model)
}
\arguments{
	\item{ww}{the cc_weights for patient dd, with G entries}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{dd}{the patient number}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative derivative of the loglikelihood function relevant to optimizing cc for patient 
	dd, the cancer profile for that patient. The derivative is taken not with respect to vv but 
	with respect to unconstrained variables via a change of variables}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
