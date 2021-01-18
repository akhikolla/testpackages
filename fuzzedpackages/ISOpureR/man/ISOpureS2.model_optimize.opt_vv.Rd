\name{ISOpureS2.model_optimize.opt_vv}
\alias{ISOpureS2.model_optimize.opt_vv}
\title{Optimize vv in step 2}
\description{This function optimizes vv, the strength parameter in the prior over the reference
	cancer profile.  Note that we don't directly optimize vv because it has constraints (must be >=1
	to guarantee real-valued likelihoods).}
\usage{
ISOpureS2.model_optimize.opt_vv(
	tumordata, 
	model, 
	NUM_ITERATIONS_RMINIMIZE, 
	iter, 
	NUM_GRID_SEARCH_ITERATIONS
	) 
}
\arguments{
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumour samples}
	\item{model}{list containing all the parameters to be optimized}
	\item{NUM_ITERATIONS_RMINIMIZE}{minimum number of iteration that the minimization algorithm runs}
	\item{iter}{the iteration number}
	\item{NUM_GRID_SEARCH_ITERATIONS}{number of times to try restarting with different initial values}
}
\value{The model with the vv parameter updated}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{optimize}
\keyword{math}
