\name{asggm}
\alias{asggm}
\alias{asggm.formula}
\alias{asggm.default}
\title{
Adaptively Sparse Gaussian Graphical Model
}
\description{
implements a parameter-free adaptively sparse Gaussian graphical model.
}
\usage{
\method{asggm}{formula}(formula, data=list(), ...)
\method{asggm}{default}(x, iterations = 100000000, init = NULL, epsilon = 0.001, ...)
}
\arguments{
  \item{formula}{
an object of class \dQuote{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
See \code{\link{lm}} Details for further information.
}
  \item{data}{
an optional data frame, list or environment containing the variables in the model.
}
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
  \item{\dots}{
further arguments
}
}
\details{
An effective approach to structure learning and parameter estimation
for Gaussian graphical models is to impose a sparsity prior, such as a
Laplace prior, on the entries of the precision matrix. We introduce a
parameter-free method for estimating a precision matrix with sparsity
that adapts to the data automatically, achieved by formulating a
hierarchical Bayesian model of the precision matrix with a
non-informative Jeffreys' hyperprior. We also naturally enforce the
symmetry and positive-definiteness constraints on the precision matrix
by parameterizing it with the Cholesky decomposition.
}
\value{
\code{asggm} returns an object of class \code{"asggm"}.

An object of class \dQuote{\code{asggm}} is a list containing at least the following components:

}
\author{
 Kristen Zygmunt, Eleanor Wong, Tom Fletcher
}

\references{
\cite{Wong, Eleanor, Suyash Awate, and P. Thomas Fletcher. \dQuote{Adaptive Sparsity in Gaussian Graphical Models.}In Proceedings of the 30th International Conference on Machine Learning (ICML-13), pp. 311-319. 2013.}
}
%%\seealso{
%%}
\examples{
A = diag(3)
asggm(A)
}
\keyword{ classes }
\keyword{ models }
\keyword{ multivariate }

