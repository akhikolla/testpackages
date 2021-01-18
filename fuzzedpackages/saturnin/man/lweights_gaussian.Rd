\name{lweights_gaussian}
\alias{lweights_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of the log-weight matrix in a gaussian setting.
}
\description{
The function computes the log-weights of all edges in a gaussian setting.
The result should be used in \link{edge.prob} with argument \code{log} set to \code{TRUE}.
Usual values are used as default for the prior normal-Wishart hyperparameters.
Computation can be parallelized by setting \code{nbcores} to more than 2. Parallelization relies on \code{parallel}.
}
\usage{
lweights_gaussian(data, 
      a = ncol(data), 
      mu = numeric(p), 
      au = 1, 
      T = diag(ncol(data), 
      ncol(data)), 
      nbcores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix containing continuous data.
}
  \item{a}{
Prior degree of freedom of the normal-Wishart distribution.
}
  \item{mu}{
Prior mean for the mean of the normal-Wishart distribution.
}
  \item{au}{
Prior relative precision of the normal-Wishart distribution.
}
  \item{T}{
Prior scale matrix of the normal-Wishart distribution.
}
  \item{nbcores}{
Number of cores to be used in parallelized computation.
}
}
\value{
\item{W}{log-weight matrix}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gaussian }

