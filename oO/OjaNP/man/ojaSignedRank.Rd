\name{ojaSignedRank}
\alias{ojaSignedRank}
\title{Oja Signed Ranks -- Affine Equivariant Multivariate Signed Ranks }
\description{
  The function computes the Oja signed rank of a point \code{x} w.r.t.
  a data set \code{X} or, if no
  point \code{x} is given, the Oja signed ranks of all points in \code{X}.
}
\usage{
ojaSignedRank(X, x = NULL, p = NULL, silent = FALSE, 
              na.action = na.fail)
}


\arguments{
\item{X}{numeric data.frame or matrix containing the data points as rows.}

\item{x}{\code{NULL} or a numeric vector, the point for which the Oja
  signed rank should be computed.}
  
\item{p}{\code{NULL} or a number between 0 and 1 which specifies the
  fraction of hyperplanes to be used for subsampling. If \code{p = 1},
  no subsampling is done. If \code{p = NULL}, the value of \code{p} is
  determined based on the size of the data set. See details.}
  
\item{silent}{logical, if subsampling is done or the expected
  computation time is too long, a warning message will be printed
  unless \code{silent} is \code{TRUE}. The default is \code{FALSE}.
}
  
\item{na.action}{a function which indicates what should happen
  when the data contain 'NA's. Default is to fail.}
}



\details{
The function computes the Oja signed rank of the point \code{x} w.r.t.
the data set \code{X} or, if no \code{x} is specified, the Oja signed
ranks of all data points in \code{X} w.r.t. \code{X}. For a definition
of \dfn{Oja signed rank} see Hettmansperger et al. (1997) formula (9).

The matrix \code{X} needs to have at least as many rows as columns in
order to give sensible results. The vector \code{x} has to be of length
\code{ncol(X)}. If \code{x} is specified, a vector of length
\code{ncol(X)} is returned. Otherwise the return value is a matrix of
the same dimensions as \code{X} where the \eqn{i}-th row contains the
Oja rank of the \eqn{i}-th row of \code{X}.

The function will also work for matrices \code{X} with only one column
and also vectors. Then (univariate) signed ranks are returned.

For \eqn{n = nrow(X)} data points in \eqn{R^k}, where \eqn{k = ncol(X)},
the computation of the Oja signed rank necessitates the evaluation of
\eqn{N = 2^k*choose(n,k)} 
hyperplanes in \eqn{R^k}. Thus for large data sets the function offers a
subsampling option in order to deliver (approximate) results within
reasonable time. The subsampling fraction is controlled by the parameter \code{p}: 
If \eqn{p < 1} is passed to the function, the computation is
based on a random sample of only \eqn{p N} of all possible \eqn{N}
hyperplanes. If \code{p} is not specified (which defaults to \code{p = NULL}), 
it is automatically determined based on \eqn{n} and \eqn{k} to yield a
sensible trade-off between accuracy and computing time.
If \eqn{N k^3 < 6 \cdot 10^6}{2^k*choose(n,k)*k^3 < 6e+06}, the
sample fraction \code{p} is set to 1 (no subsampling). If all Oja signed ranks of \code{X} are requested, a hyperplane sample is
drawn once and all Oja signed ranks are then computed based on this sample.

Finally, subsampling is feasible. Even for very small \code{p} useable
results can be expected, see e.g. the examples for the function \code{\link{ojaRCM}}.
}

\value{
Either a numeric vector, the Oja signed rank of \code{x}, or 
a matrix of the same dimensions as \code{X} containing the Oja signed ranks of \code{X} as rows.
}
\references{
\cite{Fischer D, Mosler K, \enc{Möttönen}{Mottonen} J, Nordhausen K, Pokotylo O and Vogel D (2020). “Computing the Oja Median in R: The Package OjaNP.” \emph{Journal of Statistical Software}, \bold{92}(8), pp. 1-36. doi: 10.18637/jss.v092.i08 (URL:
http://doi.org/10.18637/jss.v092.i08).}

\cite{Hettmansperger, T. P.,\ \enc{Möttönen}{Mottonen}, J. and Oja, H. (1997), 
  Affine invariant multivariate one-sample signed-rank tests, \emph{J.\
    Amer. Statist. Assoc.}, 
  \bold{92}, 1591--1600.}

\cite{Oja, H. (1999), Affine invariant multivariate sign and rank tests
  and corresponding estimates: A review, 
  \emph{Scand. J. Statist.}, \bold{26}, 319--343.}
}

\author{Jyrki \enc{Möttönen}{Mottonen}}




\seealso{
\code{\link{ojaSign}},
\code{\link{ojaRank}},
\code{\link{ojaRCM}},
\code{\link{hyperplane}}
}


\examples{
set.seed(123)
X <- rmvnorm(n = 30,mean = c(0,0)) # from package 'mvtnorm'
ojaSignedRank(X)
ojaSignedRank(X, x = c(0,0)) # zero
}
\keyword{multivariate}
\keyword{nonparametric}
