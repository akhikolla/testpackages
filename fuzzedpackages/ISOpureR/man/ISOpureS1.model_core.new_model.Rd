\name{ISOpureS1.model_core.new_model}
\alias{ISOpureS1.model_core.new_model}
\title{Initialize a model list to hold all the parameters}
\description{Produces a list (the model) which initializes the parameters vv, log_BBtranspose, 
	PPtranspose, kappa, theta, omega, log_all_rates for step 1}
\usage{
ISOpureS1.model_core.new_model(tumordata, kappa, INITIAL_VV, PPtranspose, BBtranspose)
}
\arguments{
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{kappa}{scalar strength parameter kappa placed over the reference cancer profile mm}
	\item{INITIAL_VV}{a vector with K components, the prior over mixing proportions, theta, with 
		last entry weighed more heavily }
	\item{PPtranspose}{a (K-1)xG matrix, standardized so that all entries sum to 1, see 
		ISOpure.step1.CPE.R}
	\item{BBtranspose}{a (K-1)xG matrix of the standardized normal profiles, so that they sum to 1}
}
\value{
	\item{model}{a newly generated model list to hold all the parameters vv, log_BBtranspose, 
		PPtranspose, kappa, theta, omega, log_all_rates}
    }
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
