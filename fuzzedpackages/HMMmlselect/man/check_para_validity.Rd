\name{check_para_validity}
\alias{check_para_validity}
\title{check_para_validity}
\usage{
check_para_validity(parameters_in_matrix_form, bool_hmm)
}
\arguments{
  \item{parameters_in_matrix_form}{The parameters to be checked.}
  \item{bool_hmm}{1 for HMM. 0 for Gaussian mixture.}
}
\value{
It returns the logical value of whether the parameter set is suitable for HMM/Gaussian mixture.
}
\description{
Auxiliary function that checks whether the parameter set is suitable for HMM/Gaussian mixture.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

