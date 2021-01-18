\name{multivariate_mixture_density}
\alias{multivariate_mixture_density}
\title{multivariate_mixture_density}
\usage{
multivariate_mixture_density(x, EM_result, df_t, logconstnormal, logconstT)
}
\arguments{
  \item{x}{Data to compute the log likelihood.}
  \item{EM_result}{Parameter as a dataset, including $p for probability weight, $Mu for the mean of each mixture component, and $Sigma for the standard deviation of each mixture component.}
  \item{df_t}{Degree of freedom of the t-distribution. If df_t==0, then it is treated as a normal distribution.}
  \item{logconstnormal}{Leading constant for normal density function. Used when df_t==0.}
  \item{logconstT}{Leading constant for t-distribution density function. Used when df_t>0.}
}
\value{
It returns the likelihood of the mixture normal/student-t distribution.
}
\description{
Auxiliary function that computes probability density of mixture normal/t-distribution. Used when approximating marginal likelihood.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

