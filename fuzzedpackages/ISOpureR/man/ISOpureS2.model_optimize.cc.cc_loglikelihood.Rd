\name{ISOpureS2.model_optimize.cc.cc_loglikelihood}
\alias{ISOpureS2.model_optimize.cc.cc_loglikelihood}
\title{Compute the loglikelihood relevant to the patient cancer profiles, cc, for step 2}
\description{Computes the part of the loglikelihood function relevant to optimizing cc for patient 
	dd, the cancer profile for that patient}
\usage{
ISOpureS2.model_optimize.cc.cc_loglikelihood(ww, tumordata, dd, model)
}
\arguments{
	\item{ww}{the cc_weights for patient dd, with G entries}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{dd}{the patient number}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative the part of the loglikelihood function relevant to optimizing cc for patient dd, 
	the cancer profile for that patient.}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
