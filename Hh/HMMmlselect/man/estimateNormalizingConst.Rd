\name{estimateNormalizingConst}
\alias{estimateNormalizingConst}
\title{estimateNormalizingConst}
\usage{
estimateNormalizingConst(SampDens, boolHMM, dft, RIS, IS, NsmpResmp, llUn, Mclust = TRUE)
}
\arguments{
  \item{SampDens}{The sample and density }
  \item{boolHMM}{Number that controls the method.}
  \item{dft}{Degree of freedom used in the t-distribution.}
  \item{RIS}{RIS in the algorithm.}
  \item{IS}{IS in the algorithm.}
  \item{NsmpResmp}{Number of resampling used in the algorithm.}
  \item{llUn}{Unnormalized log likelihood.}
  \item{Mclust}{Auxiliary variable. Set to \code{TRUE}.}
}
\value{
It returns the approximated normalizing constant.
}
\description{
Auxiliary function that approximate the normalizing constant. Used when approximating marginal likelihood.
}
\details{
See Manual.pdf in "inst/extdata" folder.
}

\references{
Yang Chen, Cheng-Der Fuh, Chu-Lan Kao, and Samuel Kou (2019+) "Determine the number of states
in hidden markov models via marginal likelihood." Submitted.
}

