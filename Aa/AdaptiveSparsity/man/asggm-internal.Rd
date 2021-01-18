\name{asggm-internal}
\alias{rCSL}
\alias{genData}
\alias{genL}
\title{
 asggm internal functions
}
\description{
These are the fitting and initialization functions used by asggm.  These should generally \emph{not} be used directly.
}
\usage{
rCSL(x, iterations = 500, init = NULL, epsilon = 1e-05, ansL = NULL)
genL(kNodes, spP)
genData(L, nSamples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
design matrix
}
  \item{iterations}{
number of iterations of the algorithm to run.
}
  \item{init}{
optional initialization, for instance, the cholesky of \code{x}.  If NULL, it defaults to the cholesky of \code{x}.
}
  \item{epsilon}{
amount to add for numerical stability.
}
 \item{ansL}{

}

  \item{kNodes}{

}
  \item{spP}{
}
  \item{L}{
  L created by genL
}
  \item{nSamples}{
number of samples.
}
}
\details{
rCSL calls the C++ code to compute the Wong EM algorithm.
genL and genData are used together to create example data.
}
\value{
rCSL returns a list with the following components:
}
\references{
\cite{Wong, Eleanor, Suyash Awate, and P. Thomas Fletcher. \dQuote{Adaptive Sparsity in Gaussian Graphical Models.}In Proceedings of the 30th International Conference on Machine Learning (ICML-13), pp. 311-319. 2013.}
}

\seealso{
\code{\link{asggm}}, which should be used directly instead of these methods
}
\keyword{ methods }
