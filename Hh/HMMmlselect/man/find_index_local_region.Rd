\name{find_index_local_region}
\alias{find_index_local_region}
\title{find_index_local_region}
\usage{
find_index_local_region(samples, EM_result, df_t)
}
\arguments{
  \item{samples}{The sample used to consruct region.}
  \item{EM_result}{The parameters.}
  \item{df_t}{Degree of freedom for the t-distrbution.}
}
\value{
It returns the local region for locally restricted importance sampling.
}
\description{
Auxiliary function that identify the local region for locally restricted importance sampling. Used in approximate the normalizing constant.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

