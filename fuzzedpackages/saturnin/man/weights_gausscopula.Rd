\name{weights_gausscopula}
\alias{weights_gausscopula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of the weight matrix in a gaussian copula setting.
}
\description{
The function computes the weights of all edges in a gaussian copula setting.
The result should be used in \link{edge.prob} with argument \code{log} set to \code{FALSE}.
The function brings the values of all variables back to \eqn{[0;1]} by computing univariate empirical cdf functions.
The prior distribution for the correlation of the bivariate gaussian copulas \code{prior} can be set to either \code{"uniform"} or \code{"beta"}. Beta prior is understood as a beta distribution with a change of variables to bring it back to \eqn{[-1;1]}.
Computation can be parallelized by setting \code{nbcores} to more than 2. Parallelization relies on \code{parallel}.
}
\usage{
weights_gausscopula(data, prior_type = "uniform", a = 1, b = 1, nbcores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix containing the data.
}
  \item{prior_type}{
Prior to be used for the correlation.
}
  \item{a}{
Shape parameter 1 for beta prior.
}
  \item{b}{
Shape parameter 2 for beta prior.
}
  \item{nbcores}{
Number of cores to be used in parallelized computation.
}
}
\value{
\item{W}{weight matrix.}
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

W <- weights_gausscopula(data_multinomial)

prob <- edge.prob(W, log = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gaussian copula }
