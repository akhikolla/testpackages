\name{sample_mixture}
\alias{sample_mixture}
\title{sample_mixture}
\usage{
sample_mixture(N, list_paras, df_t)
}
\arguments{
  \item{N}{The number of sample to be drawn.}
  \item{list_paras}{Parameter as a dataset, including $p for probability weight, $Mu for the mean of each mixture component, and $Sigma for the standard deviation of each mixture component.}
  \item{df_t}{Degree of freedom of the t-distribution. If df_t==0, then it is treated as a normal distribution.}
}
\value{
It returns the sample of the simulated mixture normal/student-t distribution.
}
\description{
Auxiliary function that simulates mixture normal/t-distribution. Used when approximating marginal likelihood.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

