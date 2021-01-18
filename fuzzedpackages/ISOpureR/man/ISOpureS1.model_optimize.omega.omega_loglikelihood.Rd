\name{ISOpureS1.model_optimize.omega.omega_loglikelihood}
\alias{ISOpureS1.model_optimize.omega.omega_loglikelihood}
\title{Compute the loglikelihood relevant to omega for step 1}
\description{Compute the the part of the loglikelihood function relevant to omega in step 1}
\usage{
ISOpureS1.model_optimize.omega.omega_loglikelihood(ww, tumordata, model) 
}
\arguments{
	\item{ww}{(K-1)x1 matrix, log(omega), where the entries in omega are constrained to add to 1
		where K-1 is the number of normal samples}
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative of the loglikelihood function relevant to omega}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
