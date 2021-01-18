\name{ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood}
\alias{ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood}
\title{Compute the derivative of loglikelihood relevant to omega for step 1}
\description{Compute the derivative of the part of the loglikelihood function relevant to omega 
	with respect to (log) omega, in step 1. Instead of performing constrained optimization on omega 
	directly, we optimize the log of omega in an unconstrained fashion.}
\usage{
ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood(ww, tumordata, model)
}
\arguments{
	\item{ww}{(K-1)x1 matrix, log(omega), where the entries in omega are constrained to add to 1
		where K-1 is the number of normal samples}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative derivative of the part of the loglikelihood function relevant to omega   
	with respect to (log) omega}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
