\name{ISOpureS2.model_optimize.vv.vv_compute_loglikelihood}
\alias{ISOpureS2.model_optimize.vv.vv_compute_loglikelihood}
\title{Compute loglikelihood relevant to vv for step 2}
\description{Computes the part of the loglikelihood function relevant to optimizing vv for step 2.}
\usage{
ISOpureS2.model_optimize.vv.vv_compute_loglikelihood(ww, sum_log_theta, D)
}
\arguments{
	\item{ww}{log(vv-1), a Kx1 matrix}
	\item{sum_log_theta}{the column sums of log(theta), a 1xK matrix}
	\item{D}{the number of patients (a scalar)}
}
\value{The negative of the loglikelihood relevant to optimizing vv}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
