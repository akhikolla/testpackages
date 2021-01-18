\name{ISOpureS1.model_core.optmodel}
\alias{ISOpureS1.model_core.optmodel}
\title{Optimizes the ISOpure parameters for step 1}
\description{Optimizes the ISOpure parameters for step 1 cyclically until convergence}
\usage{
ISOpureS1.model_core.optmodel(tumordata, model, NUM_ITERATIONS=35)
}
\arguments{
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{model}{list containing all the parameters to be optimized}
	\item{NUM_ITERATIONS}{(optional) minimum number of iterations of optimization algorithm, 
		default is 35}
}
\value{
	\item{model}{updated model list containing all the parameters}
    }
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
