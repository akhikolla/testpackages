\name{ISOpureS2.model_core.new_model}
\alias{ISOpureS2.model_core.new_model}
\title{Compute loglikelihood given all model parameters for step 2}
\description{Produces a list (the model) which initializes the parameters vv, log_BBtranspose, 
	PPtranspose, kappa, theta, omega, log_all_rates for step 2}
\usage{
ISOpureS2.model_core.new_model(tumordata, kappa, INITIAL_VV, PPtranspose, BBtranspose)
}
\arguments{
	\item{tumordata}{a GxD matrix representing gene expression profiles of tumor samples}
	\item{kappa}{a 1xD matrix which represents strength parameter kappa over cc, given the reference 
		profile mm}
	\item{INITIAL_VV}{a vector with K components, the prior over mixing proportions, theta, with last 
		entry weighed more heavily}
	\item{PPtranspose}{the prior on the tumor-specific cancer profiles is just the reference cancer
        profile (1xG matrix) learned in ISOpureS1, standardized so that all entries sum to 1 }
	\item{BBtranspose}{a (K-1)xG matrix of the standardized normal profiles, so that they sum to 1}
}
\value{
	\item{model}{a newly generated model list to hold all the parameters}
    }
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
