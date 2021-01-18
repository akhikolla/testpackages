\name{ojaRank}
\alias{ojaRank}
\title{ Oja Ranks -- Affine Equivariant Multivariate Ranks }
\description{
  The function computes the Oja rank of a point \code{x} w.r.t. a data set \code{X} or, if no
  point \code{x} is given, the Oja ranks of all points in \code{X}.
}
\usage{
ojaRank(X, x = NULL, p = NULL, silent = FALSE, na.action = na.fail)
}


\arguments{
\item{X}{numeric data.frame or matrix containing the data points as rows.}

\item{x}{\code{NULL} or a numeric vector, the point for which the Oja rank should be computed.}
  
\item{p}{\code{NULL} or a number between 0 and 1 which specifies the fraction of hyperplanes to be used for subsampling. If \code{p = 1}, no subsampling is done.
If \code{p = NULL}, the value of \code{p} is determined based on the size of the data set. See details.}
  
\item{silent}{logical, if subsampling is done or the expected computation time is too long, a warning message will be printed unless \code{silent} is \code{TRUE}. The default is \code{FALSE}.
}
  
\item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail.}
}




\details{
The function computes the Oja rank of the point \code{x} w.r.t.\ the data set \code{X} or, if no \code{x} is specified, the Oja ranks of all data points in 
\code{X} w.r.t. \code{X}. For a definition of \dfn{Oja rank} see reference below. 

The matrix \code{X} needs to have at least as many rows as columns in order to give sensible results. 
The vector \code{x} has to be of length \code{ncol(X)}. 
If \code{x} is specified, a vector of length \code{ncol(X)} is returned. Otherwise the return value is a matrix of the same dimensions as \code{X} where the \eqn{i}-th row contains the 
Oja rank of the \eqn{i}-th row of \code{X}.

The function will also work for matrices \code{X} with only one column and also vectors. Then centered and normalized (univariate) ranks are returned.


For \eqn{n = nrow(X)} data points in \eqn{R^k}, where \eqn{k = ncol(X)}, the computation of the Oja rank necessitates the evaluation of \eqn{N = choose(n,k)} hyperplanes in \eqn{R^k}. 
Thus for large data sets the function offers a subsampling option in order to deliver (approximate) results within reasonable time. The subsampling fraction is controlled by the parameter 
\code{p}: 
If \eqn{p < 1} is passed to the function, the computation is based on a random sample of only \eqn{pN} of all possible \eqn{N} hyperplanes. If \code{p} is not specified 
(which defaults to \code{p = NULL}), it is automatically determined based on \eqn{n} and \eqn{k} to yield a sensible trade-off between accuracy and computing time. If 
\eqn{N k^3 < 6 \cdot 10^6}{choose(n,k)*k^3 < 6e+06}, the sample fraction \code{p} is set to 1 (no subsampling). Otherwise \code{p} is choosen such that the computation 
(of one rank) \emph{usually} takes around 20 seconds (on a 1.66 GHz CPU and 1 GB RAM). If all Oja ranks of \code{X} are requested, a hyperplane sample is drawn once, 
all Oja ranks are then computed based on this sample.

Finally, subsampling is feasible. Even for very small \code{p} useable results can be expected, see e.g. the examples for the function \code{\link{ojaRCM}}.

Claudia \enc{Köllmann}{Kollmann} is acknowledged for bug-fixing this function.
}

\value{
either a numeric vector, the Oja rank of \code{x}, or 
a matrix of the same dimensions as \code{X} containing the Oja ranks of \code{X} as rows.
}

\author{Daniel Vogel, Jyrki \enc{Möttönen}{Mottonen}}

\references{ 
\cite{Fischer D, Mosler K, \enc{Möttönen}{Mottonen} J, Nordhausen K, Pokotylo O and Vogel D (2020). “Computing the Oja Median in R: The Package OjaNP.” \emph{Journal of Statistical Software}, \bold{92}(8), pp. 1-36. doi: 10.18637/jss.v092.i08 (URL:
http://doi.org/10.18637/jss.v092.i08).}

\cite{Oja, H. (1999), Affine invariant multivariate sign and rank tests and corresponding estimates: A review, 
\emph{Scand. J. Statist.}, \bold{26}, 319--343.}
}


\seealso{
\code{\link{ojaSign}}, \code{\link{ojaRCM}}, \code{\link{hyperplane}}, \code{\link{ojaSignedRank}}
}


\examples{
### ----<< Example 1 >>---- : 30 points in R^2
\dontrun{
set.seed(123)
X <- rmvnorm(n = 30,mean = c(0,0)) # from package 'mvtnorm'
ojaRank(X)
ojaRank(X, x = c(100,100))
ojaRank(X, x = ojaMedian(X, alg="exact")) # close to zero

# The following two return the same (only in different time)
ojaRank(X)
t(apply(X, 1, function(y){ojaRank(X,y)}))
# but the following two do not (due to different subsampling).
# 1)
set.seed(123); ojaRank(X, p = 0.9, silent = TRUE)
# 2)
set.seed(123) 
t(apply(X, 1, function(y){ojaRank(X, y, p = 0.9, silent = TRUE)}))
# In 1) one subsample for all ranks is drawn, whereas in 2)
# a different sample for each rank is drawn.
}

### ----<< Example 2 >>---- : three points in R^3: only one hyperplane
# The following commands return the same result.
\dontrun{
ojaRank(X = diag(rep(1, 3)), x = c(0,0,0))
ojaRank(X = diag(rep(1, 3)), x = c(-100,-110,-5550))
hyperplane(X = diag(rep(1,3)))[1:3]
}


### ----<< Example 3 >>---- : 300 points in R^7
# Subsampling is done.
# The following example might take a bit longer:
\dontrun{
set.seed(123)
X <- rmvnorm(n = 300, mean = rep(0, 7))
system.time(or <- ojaRank(x = 1:7, X = X))
# PLEASE NOTE: The computation of the Oja rank is based on a 
# random sub-sample of less than 1\% of all possible hyperplanes.
#
#       user      system     elapsed 
#      18.47        0.00       18.47 
print(or,d=4)
# [1]  7.733  6.613  6.839  7.383 18.237 21.851 23.700
}


### ----<< Example 4 >>---- : univariate ranks
\dontrun{
ojaRank(1:10)
ojaRank(X = 1:10, x = 5.5)
}
}

\keyword{multivariate}
\keyword{nonparametric}
