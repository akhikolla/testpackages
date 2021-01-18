\name{lweights_multinomial}
\alias{lweights_multinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of the log-weight matrix in a multinomial setting.
}
\description{
The function computes the log-weights of all edges in a multinomial setting.
The result should be used in \link{edge.prob} with argument \code{log} set to \code{TRUE}.
Prior counts can be generated using the function \code{prior_unif_dirichlet}.
Computation can be parallelized by setting \code{nbcores} to more than 2. Parallelization relies on \code{parallel}.
}
\usage{
lweights_multinomial(data, prior = defaut.prior, nbcores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix containing discrete data.
}
  \item{prior}{
Prior to be used for the Dirichlet distribution.
}
  \item{nbcores}{
Number of cores to be used in parallelized computation.
}
}
\value{
\item{W}{log-weight matrix.}
}
\references{
This package implements the method described in the paper "Bayesian Inference of Graphical Model Structures Using Trees" by L. Schwaller, S. Robin, M. Stumpf, 2015 (submitted and availavable on arXiv).
}
\author{
Loïc Schwaller
}

\examples{
library('saturnin')
data(data_multinomial)

W <- lweights_multinomial(data_multinomial)

prob <- edge.prob(W, log = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multinomial }

