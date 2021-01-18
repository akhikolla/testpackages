\name{aslm-internal}
\alias{figEM}
\alias{fit.ols.lm}
\alias{init.ones}
\alias{init.rnorm}
\alias{init.runif}
\title{
 aslm internal functions
}
\description{
These are the fitting and initialization functions used by aslm.  These should generally \emph{not} be used.
}
\usage{
figEM(x, y, init = NULL, stopDiff = 1e-08, epsilon = 1e-06, a = 1)
fit.ols.lm(x, y)
init.ones(x, y)
init.rnorm(x, y)
init.runif(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
design matrix of dimension \code{ n * p }.
}
  \item{y}{
vector of observations of length n, or a matrix with n rows.
}
  \item{init}{
optional initialization, a list with components containing an initial estimate for \code{beta} and \code{sigma}
}
  \item{stopDiff}{
convergence criteria.  Algorithm stops once difference in beta and sigma from one iteration to the next is less than stopDiff.
}
  \item{epsilon}{
amount to add to beta for numerical stability,
}
  \item{a}{
scaling of sigmaSqr to provide numerical stability for solving steps.
}
}
\details{
figEM computes the Figueiredo EM algorithm for adaptive sparsity using Jeffreys prior.

fit.ols.lm computes an initial beta and sigma based on finding the lm.fit of the full design matrix.

init.ones computes an initial beta that is all ones and computes the associated sigmas.

init.rnorm computes an initial beta that is normally distributed with a mean of 0 and a standard deviation of 50

init.runif computes an initial beta that is uniformly distributed from 0 to 1

Currently, figEM uses fit.ols.lm to initialize beta and sigma if no init list is provided.
}
\value{
figEM returns a list with the following components:
\item{ coefficients }{ \code{p} vector (also known as \code{beta}). }
\item{ vcov }{ variance-covariance matrix. }
\item{ sigma }{ norm of the model error. }
\item{ df }{ degrees of freedom of residuals. }

fit.ols.lm and init.ones are used to initialize beta and sigma if init is not provided to figEM.
Each of these functions returns a list with the following components:
\item{ beta }{ initial \code{p} vector. }
\item{ sigma }{ initial norm of the model error based on this initial beta. }
}
\references{
\cite{Figueiredo, M.A.T.; , \dQuote{Adaptive sparseness for supervised learning}, Pattern Analysis and Machine Intelligence, IEEE Transactions on , vol.25, no.9, pp. 1150- 1159, Sept. 2003}
}

\seealso{
\code{\link{aslm}}, which should be used directly instead of these methods
}

\keyword{ methods }
