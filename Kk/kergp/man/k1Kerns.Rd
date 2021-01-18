\name{k1Exp}
\alias{k1Exp}
\alias{k1PowExp}
\alias{k1Gauss}
\alias{k1Matern3_2}
\alias{k1Matern5_2}

\docType{data}
\title{
  
  Predefined covMan Objects for 1D Kernels

}
\description{
  
   Predefined kernel Objects as \code{covMan} objects.

 }
\usage{
k1Exp
k1Matern3_2
k1Matern5_2
k1Gauss
}
\format{
  
   Objects with class \code{"covMan"}.

 }
\details{

  These objects are provided mainly as examples. They are used
  \code{\link{covTS}}.

}
%% \source{ }
%% \references{ }
\examples{
set.seed(1234)
x <- sort(runif(40))
X <- cbind(x = x)
yExp <- simulate(k1Exp, nsim = 20, X = X)
matplot(X, yExp, type = "l", col = "SpringGreen", ylab = "")
yGauss <- simulate(k1Gauss, nsim = 20, X = X)
matlines(X, yGauss, col = "orangered")
title("Simulated paths from 'k1Exp' (green) and 'k1Gauss' (orange)")

## ============================================================================
## You can build a similar object using a creator of
## 'covMan'. Although the objects 'k1Gauss' and 'myk1Gauss' differ,
## they describe the same mathematical object.
## ============================================================================

myk1Gauss <- kGauss(d = 1)


}
\keyword{datasets}
