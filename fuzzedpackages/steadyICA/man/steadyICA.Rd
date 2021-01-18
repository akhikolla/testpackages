\name{steadyICA}
\alias{steadyICA}
\title{
  Estimate independent components by minimizing distance covariance }
\description{The model is: X = S M + E, where X is n x p and has mean zero, S is n x d, M is d x p, and E is measurement error. For whitened data, we have Z = S t(W), where W is orthogonal. We find the matrix M such that S minimizes the distance covariance dependency measure.}

\usage{
steadyICA(X, n.comp = ncol(X), w.init = NULL, PIT = FALSE, bw = 'SJ', adjust = 1,
whiten = FALSE, irlba = FALSE, symmetric = FALSE, eps = 1e-08, alpha.eps = 1e-08,
maxit = 100, method = c('Cpp','R'), verbose = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The n x p data matrix, where n is the number of observations.}
  \item{n.comp}{number of components to be estimated}
  \item{w.init}{a p x d initial unmixing matrix}
  \item{PIT}{logical; if TRUE, the distribution and density of the independent components are estimated using gaussian kernel density estimates.}
  \item{bw}{Argument for bandwidth selection method; defaults to 'SJ'; see stats::density}
  \item{adjust}{adjust bandwidth selection; e.g., if observations are correlated, consider using adjust > 1; see stats::density}
  \item{whiten}{logical; if TRUE, whitens the data before applying ICA, i.e., X\%*\%whitener = Z, where Z has mean zero and empirical covariance equal to the identity matrix, and Z is then used as the input.}
  \item{irlba}{logical; when whiten=TRUE, irlbA=TRUE uses the R-package 'irlba' in the whitening, which is generally faster than base::svd though sometimes less accurate}
  \item{symmetric}{logical; if TRUE, implements the symmetric version of the ICA algorithm, which is invariant to the ordering of the columns of X but is slower}
  \item{eps}{algorithm terminates when the norm of the gradient of multidcov is less than eps}
  \item{maxit}{maximum number of iterations}
  \item{alpha.eps}{tolerance controlling the level of annealing: algorithm terminates with a warning if the learning parameter is less than alpha.eps}
  \item{method}{options 'Cpp' (default), which requires the package 'Rcpp', or 'R', which is solely written in R but is much slower}
  \item{verbose}{logical; if TRUE, prints the value of multidcov, norm of the gradient, and current value of the learning parameter.}
  }

\value{
\item{S}{the estimated independent components}
\item{W}{the estimated unmixing matrix: if whiten=TRUE, W is orthogonal and corresponds to Z W = S; if whiten=FALSE, corresponds to X ginv(M) = S}
\item{M}{Returns the estimated mixing matrix for the model X = S M, where X is not pre-whitened (although X is centered)}
\item{f}{the value of the objective function at the estimated S}
\item{Table}{summarizes algorithm status at each iteration}
\item{convergence}{1 if norm of the gradient is less than eps, 2 if the learning parameter was smaller than alpha.eps, which usually means the gradient is sufficiently small, 0 otherwise}
}

\author{
Benjamin Risk
}

\references{
Matteson, D. S. & Tsay, R. Independent component analysis via U-Statistics. 
    <http://www.stat.cornell.edu/~matteson/#ICA>
%% ~put references to the literature/web site here ~
}

\seealso{
\code{\link{multidcov}}
%%\code{\link{GPArotation::GPForth}}
%%\code{\link[GPArotation:GPForth]{GPArotation::GPForth}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(999)
nObs <- 1024
nComp <- 3
# simulate from some gamma distributions:
simS<-cbind(rgamma(nObs, shape = 1, scale = 2),
            rgamma(nObs, shape = 3, scale = 2),
            rgamma(nObs, shape = 9, scale = 0.5))

#standardize by expected value and variance:
simS[,1] = (simS[,1] - 1*2)/sqrt(1*2^2)
simS[,2] = (simS[,2] - 3*2)/sqrt(3*2^2)
simS[,3] = (simS[,3] - 9*0.5)/sqrt(9*0.5^2)

# slightly revised 'mixmat' function (from ProDenICA)
# for p>=d: uses fastICA and ProDenICA parameterization:
myMixmat <- function (p = 2, d = NULL) {
  if(is.null(d)) d = p
  a <- matrix(rnorm(d * p), d, p)
  sa <- La.svd(a)
  dL <- sort(runif(d) + 1)
  mat <- sa$u\%*\%(sa$vt*dL)
  attr(mat, "condition") <- dL[d]/dL[1]
  mat
}

simM <- myMixmat(p = 6, d = nComp)
xData <- simS\%*\%simM
xWhitened <- whitener(xData, n.comp = nComp)

#estimate mixing matrix:
est.steadyICA.v1 = steadyICA(X = xData,whiten=TRUE,n.comp=nComp,verbose = TRUE)

#Define the 'true' W:
W.true <- solve(simM\%*\%xWhitened$whitener)

frobICA(M1=est.steadyICA.v1$M,M2=simM)
frobICA(S1=est.steadyICA.v1$S,S2=simS)

\dontrun{
#now initiate from target:
est.steadyICA.v2 = steadyICA(X = xData, w.init= W.true, n.comp = nComp, whiten=TRUE, verbose=TRUE)

#estimate using PIT steadyICA such that dimension reduction is via ICA:
est.steadyICA.v3 = steadyICA(X = xData, w.init=ginv(est.steadyICA.v2$M),
PIT=TRUE, n.comp = nComp, whiten=FALSE, verbose=TRUE)

frobICA(M1=est.steadyICA.v2$M,M2=simM)
frobICA(M1=est.steadyICA.v3$M,M2=simM) 
frobICA(S1=est.steadyICA.v2$S,S2=simS)

#tends to be lower than PCA-based (i.e., whitening) methods:
frobICA(S1=est.steadyICA.v3$S,S2=simS) 

# JADE uses a different parameterization and different notation.
# Using our parameterization and notation, the arguments for 
# JADE::amari.error correspond to:
amari.error(t(W.hat), W.true)

library(JADE)

amari.error(t(est.steadyICA.v1$W), W.true) 
amari.error(t(est.steadyICA.v2$W), W.true)
##note that a square W is not estimated if PIT=TRUE and whiten=FALSE

#Compare performance to fastICA:
library(fastICA)
est.fastICA = fastICA(X = xData, n.comp = 3, tol=1e-07)
amari.error(t(est.fastICA$W), W.true)
##steadyICA usually outperforms fastICA

##Compare performance to ProDenICA:
library(ProDenICA)
est.ProDenICA = ProDenICA(x = xWhitened$Z, k = 3, maxit=40,trace=TRUE)
amari.error(t(est.ProDenICA$W), W.true)
##ProDenICA and steadyICA tend to be similar when sources
##are continuously differentiable
}
}

\keyword{distance}
\keyword{covariance}% __ONLY ONE__ keyword per line
\keyword{ICA}
\keyword{independent}
