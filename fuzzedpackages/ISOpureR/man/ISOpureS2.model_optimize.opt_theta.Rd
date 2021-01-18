\name{ISOpureS2.model_optimize.opt_theta}
\alias{ISOpureS2.model_optimize.opt_theta}
\title{Optimize theta in step 2}
\description{This function optimizes theta, in fact theta_weights.  Since thetas are constrained 
	(must be parameters of multinomial/discrete distribution), we don't directly optimize the 
	likelihood function w.r.t. theta, but we perform change of variables to do unconstrained
	optimization.  We therefore store these unconstrained variables in the field "theta_weights",
	and update these variables.}
\usage{
ISOpureS2.model_optimize.opt_theta(
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
\value{The model with the theta parameter updated (the first K-1 columns) corresponding to the 
	normal sample contributions}
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{optimize}
\keyword{math}
