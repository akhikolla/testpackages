\name{edge.prob}
\alias{edge.prob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of posterior edge appearance probabilities in a random tree.
}
\description{
The function computes posterior edge appearance probabilities in a random tree from a (log-)weight matrix.
The (log-)weight matrix can be obtained from one of the functions \code{lweights_multinomial}, \code{lweights_gaussian} or \code{weights_gausscopula}.
The function can also account for prior edge appearance probability.
}
\usage{
edge.prob(W, log = TRUE, account.prior = FALSE, q0 = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
(log-)weight matrix
}
  \item{log}{
\code{TRUE} when using a log-weight matrix, \code{FALSE} otherwise.
}
  \item{account.prior}{
\code{FALSE} for no accounting, \code{TRUE} otherwise.
}
  \item{q0}{
Desired prior edge appearance probability.
}
}

\value{
\item{prob}{Posterior edge appearance probability matrix.}
}
\references{
This package implements the method described in the paper "Bayesian Inference of Graphical Model Structures Using Trees" by L. Schwaller, S. Robin, M. Stumpf, 2015 (submitted and availavable on arXiv).
}
\author{
Loïc Schwaller
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library('saturnin')
data(data_gaussian)

W <- lweights_gaussian(data_gaussian)
prob <- edge.prob(W, log = TRUE, account.prior = TRUE, q0 = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ posterior edge probability }

