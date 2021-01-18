\name{ISOpureS1.model_optimize.omega.omega_compute_loglikelihood}
\alias{ISOpureS1.model_optimize.omega.omega_compute_loglikelihood}
\title{Compute loglikelihood relevant to omega for step 1}
\description{Computes the part of the loglikelihood function relevant to optimizing omega for step 1}
\usage{
ISOpureS1.model_optimize.omega.omega_compute_loglikelihood(omega, tumordata, model)
}
\arguments{
	\item{omega}{(K-1)x1 matrix representing the weights of the normal profiles B_i used to
         make the weighted combination that forms the mean parameter vector for the 
         Dirichlet distribution over m}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The part of the loglikelihood function relevant to optimizing omega}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
