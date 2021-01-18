\name{ISOpureS1.model_optimize.opt_kappa}
\alias{ISOpureS1.model_optimize.opt_kappa}
\title{Optimize kappa in step 1}
\description{This function optimizes kappa, the strength parameter in the prior over the reference
	cancer profile.  Note that we don't directly optimize kappa because it has constraints (must be 
	greater than the minimum determined in ISOpure.step1.CPE.)}
\usage{
ISOpureS1.model_optimize.opt_kappa(
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
\value{The model with the kappa parameter updated}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{optimize}
\keyword{math}
