\name{account.for.prior}
\alias{account.for.prior}

\title{
Accounting for prior edge appearance probability.
}
\description{
The function transforms the posterior edge appearance probability matrix given by \link{edge.prob} to account for prior edge appearance probability. For further
details on the transformation, we refer the reader to the paper quoted in the references section.
The function can be directly applied in \link{edge.prob} by setting \code{account.prior} to \code{TRUE}.
}
\usage{
account.for.prior(prob, q0)
}

\arguments{
  \item{prob}{
Posterior edge appearance probability matrix.
}
  \item{q0}{
Desired prior edge appearance probability.
}
}
\value{
\item{prob.q0 }{Transformed posterior edge appearance probability matrix.}
}
\references{
This package implements the method described in the paper "Bayesian Inference of Graphical Model Structures Using Trees" by L. Schwaller, S. Robin, M. Stumpf, 2015 (submitted and availavable on arXiv).
}
\author{
Loïc Schwaller
}

\examples{
library('saturnin')
data(data_gaussian)

W <- lweights_gaussian(data_gaussian)
prob <- edge.prob(W, log = TRUE)

prob.q0 <- account.for.prior(prob, q0 = 0.5)
}

\keyword{ posterior edge probability }

