\name{ISOpureS1.model_optimize.mm.mm_deriv_loglikelihood}
\alias{ISOpureS1.model_optimize.mm.mm_deriv_loglikelihood}
\title{Compute the derivative of the loglikelihood relevant to mm for step 1}
\description{Computes the derivative of the loglikelihood function relevant to optimizing the 
	reference cancer profile, mm, for step 1}
\usage{
ISOpureS1.model_optimize.mm.mm_deriv_loglikelihood(ww, tumordata, model)
}
\arguments{
	\item{ww}{the mm_weights, with G entries}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative derivative the likelihood function relevant to optimizing mm.  The derivative 
	is taken not with respect to mm but with respect to unconstrained variables via a change of 
	variables.}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
