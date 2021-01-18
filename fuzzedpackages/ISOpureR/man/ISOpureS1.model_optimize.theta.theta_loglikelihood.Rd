\name{ISOpureS1.model_optimize.theta.theta_loglikelihood}
\alias{ISOpureS1.model_optimize.theta.theta_loglikelihood}
\title{Compute the loglikelihood relevant to theta for step 1}
\description{Computes the part of the loglikelihood function relevant to optimizing theta for step 1}
\usage{
ISOpureS1.model_optimize.theta.theta_loglikelihood(ww, tumordata, dd, model)
}
\arguments{
	\item{ww}{the theta weights corresponding to patient dd, a 1xK matrix}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{dd}{the patient number}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative of the loglikelihood relevant to theta}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
