\name{aslm}
\alias{aslm}
\alias{getSparseModel}
\alias{aslm.formula}
\alias{aslm.default}
\title{
Adaptive Sparse Linear Model
}
\description{
implements the adaptive sparse linear model using Figueiredo's EM algorithm 
for adaptive sparsity (Jeffreys prior)
}
\usage{
\method{aslm}{formula}(formula, data=list(), na.action=na.omit, ...)
\method{aslm}{default}(x, y, ...)
getSparseModel(object)
}
\arguments{
  \item{formula}{
an object of class \dQuote{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
See \code{\link{lm}} Details for further information.
}
  \item{data}{
an optional data frame, list or environment containing the variables in the model.
}
  \item{na.action}{
action to use when data contains NAs.  Options include na.omit, na.exclude, na.fail
}
  \item{x}{
design matrix
}
  \item{y}{
vector of observations
}
  \item{\dots}{
further arguments
}
  \item{object}{
  an object of class \dQuote{\link{aslm}}.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\code{aslm} returns an object of class \code{c("aslm", "lm")}.

An object of class \dQuote{\code{aslm}} is a list containing at least the following components:
\item{coefficients}{ a named vector of coefficients }
\item{residuals}{ the residuals, that is response minus fitted values. }
\item{fitted.values}{ the fitted mean values }
\item{rank}{ the numeric rnak of the fitted linear model }
\item{df}{ the residual degrees of freedom }
\item{call}{ the matched call }
\item{terms}{ the terms object used }
\item{sigma}{  }

\code{getSparseModel} returns an object of class \code{"lm"} that is a model consisting of only the sparse nonzero variables from the original model.
}
\author{
 Kristen Zygmunt, Eleanor Wong, Tom Fletcher
}
\references{
\cite{Figueiredo, M.A.T.; , \dQuote{Adaptive sparseness for supervised learning}, Pattern Analysis and Machine Intelligence, IEEE Transactions on , vol.25, no.9, pp. 1150- 1159, Sept. 2003}
}
\seealso{
\code{\link{summary.aslm}}, \code{\link{logLik.aslm}}, \code{\link{print.aslm}}
}
\examples{
s = aslm(Infant.Mortality~.,data=swiss)
m = getSparseModel(s)

summary(s)
coef(m)
}
\keyword{ classes }
\keyword{ models }
\keyword{ multivariate }

