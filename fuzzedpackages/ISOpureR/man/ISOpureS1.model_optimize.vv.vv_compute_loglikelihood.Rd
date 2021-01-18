\name{ISOpureS1.model_optimize.vv.vv_compute_loglikelihood}
\alias{ISOpureS1.model_optimize.vv.vv_compute_loglikelihood}
\title{Compute loglikelihood relevant to vv for step 1}
\description{Computes the part of the loglikelihood function relevant to optimizing vv for step 1.}
\usage{
ISOpureS1.model_optimize.vv.vv_compute_loglikelihood(vv, sum_log_theta, DD)
}
\arguments{
	\item{vv}{Kx1 matrix representing the weights of the normal profiles B_i used to make the 
		weighted combination that forms the mean parameter vector for the Dirichlet distribution
		over m}
	\item{sum_log_theta}{the column sums of log(theta), a 1xK matrix}
	\item{DD}{the number of patients (a scalar)}
}
\value{The negative of the loglikelihood relevant to optimizing vv}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
