\name{ISOpureS1.model_optimize.opt_mm}
\alias{ISOpureS1.model_optimize.opt_mm}
\title{Optimize the reference cancer profile, m, in step 1}
\description{The goal of this function is to optimize the reference cancer profile mm.
	Because mm is constrained (must be parameters of multinomial/discrete
	distribution), we don't directly optimize the likelihood function w.r.t.
	mm, but we perform change of variables to do unconstrained
	optimization.  We therefore store these unconstrained variables in the
	field "mm_weights", and update these variables.}
\usage{
ISOpureS1.model_optimize.opt_mm(
	tumordata, model, 
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
\value{The model with mm_weights updated (and log_all_rates)}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{optimize}
\keyword{math}
