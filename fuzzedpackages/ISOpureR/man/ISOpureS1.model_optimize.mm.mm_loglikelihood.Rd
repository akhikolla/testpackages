\name{ISOpureS1.model_optimize.mm.mm_loglikelihood}
\alias{ISOpureS1.model_optimize.mm.mm_loglikelihood}
\title{Compute the loglikelihood relevant to mm for step 1}
\description{Computes the loglikelihood function relevant to optimizing the reference cancer 
	profile, mm, for step 1}
\usage{
ISOpureS1.model_optimize.mm.mm_loglikelihood(ww, tumordata, model)
}
\arguments{
	\item{ww}{the mm_weights, with G entries}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative of the likelihood function relevant to optimizing mm.}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
