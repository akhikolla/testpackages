\name{ISOpureS2.model_optimize.opt_cc}
\alias{ISOpureS2.model_optimize.opt_cc}
\title{Optimize the tumor-specific cancer profiles in step 2}
\description{Optimize the tumor-specific cancer profiles.
	Because cc is constrained (each cc_i are parameters of multinomial/discrete distribution), we 
	don't directly optimize the likelihood function w.r.t. cc, but we perform change of variables 
	to do unconstrained optimization.  We therefore store these unconstrained variables in the
	field "cc_weights", and update these variables.}
\usage{
ISOpureS2.model_optimize.opt_cc(
	tumordata, 
	model, 
	NUM_ITERATIONS_RMINIMIZE, 
	iter, 
	NUM_GRID_SEARCH_ITERATIONS) 
}
\arguments{
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumour samples}
	\item{model}{list containing all the parameters to be optimized}
	\item{NUM_ITERATIONS_RMINIMIZE}{minimum number of iteration that the minimization algorithm runs}
	\item{iter}{the iteration number}
	\item{NUM_GRID_SEARCH_ITERATIONS}{number of times to try restarting with different initial values}
}
\value{The model with cc_weights and log_cc updated}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{optimize}
\keyword{math}
