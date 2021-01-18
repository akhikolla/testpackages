\name{ISOpureS2.model_core.compute_loglikelihood}
\alias{ISOpureS2.model_core.compute_loglikelihood}
\title{Compute loglikelihood given all model parameters for step 2}
\description{Computes complete loglikelihood given all model parameters for step 2}
\usage{
ISOpureS2.model_core.compute_loglikelihood(tumordata, model)
}
\arguments{
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{model}{list containing all the parameters updated in ISOpure step two iterations}
}
\value{The scalar value of the complete loglikelihood obtained given the model parameters}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
