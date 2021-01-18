\name{ISOpure.model_optimize.vv.vv_deriv_loglikelihood}
\alias{ISOpure.model_optimize.vv.vv_deriv_loglikelihood}
\title{Compute the derivative of the loglikelihood relevant to vv for step 1}
\description{Computes the derivative of the loglikelihood function relevant to optimizing vv for 
	step 1}
\usage{
ISOpure.model_optimize.vv.vv_deriv_loglikelihood(ww, sum_log_theta, DD)
}
\arguments{
	\item{ww}{log(vv-1), a Kx1 matrix}
	\item{sum_log_theta}{the column sums of log(theta), a 1xK matrix}
	\item{DD}{the number of patients (a scalar)}
}
\value{The negative derivative of the part of the loglikelihood function relevant to vv with 
	respect to (log) vv}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
