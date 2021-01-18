\name{ISOpureS2.model_optimize.theta.theta_deriv_loglikelihood}
\alias{ISOpureS2.model_optimize.theta.theta_deriv_loglikelihood}
\title{Compute the derivative of loglikelihood relevant to theta for step 2}
\description{Computes the derivative of the loglikelihood function relevant to optimizing theta, not 
	with respect to theta but with respect to unconstrained variables}
\usage{
ISOpureS2.model_optimize.theta.theta_deriv_loglikelihood(ww, tumordata, dd, model)
}
\arguments{
	\item{ww}{the theta weights corresponding to patient dd, a 1xK matrix}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{dd}{the patient number}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative derivative of the loglikelihood function relevant to optimizing theta, not 
	with respect to theta but with respect to unconstrained variables.}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
