\name{ISOpureS2.model_optimize.kappa.kappa_loglikelihood}
\alias{ISOpureS2.model_optimize.kappa.kappa_loglikelihood}
\title{Compute loglikelihood relevant to kappa for step 2}
\description{Computes the part of the loglikelihood function relevant to optimizing kappa for step 2. 
	Instead of performing constrained optimization on kappa directly, we optimize the log of kappa
   	in an unconstrained fashion.}
\usage{
ISOpureS2.model_optimize.kappa.kappa_loglikelihood(log_kappa, model)
}
\arguments{
	\item{log_kappa}{the 1xD matrix log(kappa - model\$MIN_KAPPA)}
	\item{model}{list containing all the parameters to be optimized}
}
\value{The negative of the loglikelihood relevant to optimizing kappa}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
