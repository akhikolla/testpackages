\name{ISOpureS1.model_optimize.opt_omega}
\alias{ISOpureS1.model_optimize.opt_omega}
\title{Optimize omega in step 1}
\description{This function optimizes omega, in fact the convex mixing weights that govern prior 
	over the reference cancer profile.}
\usage{
ISOpureS1.model_optimize.opt_omega(
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
\value{The model with the omega_weights and omega parameters updated}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{optimize}
\keyword{math}
