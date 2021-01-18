\name{infomaxICA}
\alias{infomaxICA}
\title{Estimates independent components via infomax}
\description{
Estimate independent components using the infomax criteria, which is equivalent to maximum likelihood using the logistic density, exp(-S)/(1+exp(-S))^2.}
\usage{
infomaxICA(X, n.comp, W.list = NULL, whiten = FALSE, maxit = 500, eps = 1e-08, 
alpha.eps = 1e-08, verbose = FALSE, restarts=0) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the n x p data matrix}
  \item{n.comp}{number of components to be estimated}
  \item{W.list}{list of orthogonal matrices for initialization}
  \item{whiten}{Whitens the data before applying ICA, i.e., X\%*\%whitener = Z, where Z has mean zero and empirical covariance equal to the identity matrix, and Z is then used as the input.}
  \item{maxit}{maximum number of iterations}
  \item{eps}{algorithm terminates when the norm of the gradient is less than eps}
  \item{alpha.eps}{tolerance controlling the level of annealing: algorithm terminates with a warning if the learning parameter is less than alpha.eps}
  \item{verbose}{if TRUE, prints (1) the value of the infomax objective function at each iteration, (2) the norm of the gradient, and (3) current value of the learning parameter alpha.}
  \item{restarts}{An integer determining the number of initial matrices to use in estimating the ICA model. The objective function has local optima, so multiple starting values are recommended. If whiten=TRUE, then generates random orthogonal matrices. If whiten=FALSE, generate random matrices from rnorm(). See code for details.}
}
\details{
This is an R version of ICA using the infomax criteria that provides an alternative to Matlab code (\url{ftp://ftp.cnl.salk.edu/pub/tony/sep96.public}), but with a few modifications. First, we use the full data (the so-called offline algorithm) in each iteration rather than an online algorithm with batches. Secondly, we use an adaptive method to choose the step size (based upon Bernaards and Jennrich 2005), which speeds up convergence. We also omitted the bias term (intercept) included in the original formulation because we centered our data. 
}
\value{
\item{S}{the estimated independent components}
\item{W}{if whiten=TRUE, returns the orthogonal unmixing matrix; no value is returned when whiten=FALSE}
\item{M}{Returns the estimated mixing matrix for the model X = S M, where X is not pre-whitened (although X is centered)}
\item{f}{the value of the objective function at the estimated S}
\item{Table}{summarizes algorithm status at each iteration}
\item{convergence}{1 if norm of the gradient is less than eps, 2 if the learning parameter was smaller than alpha.eps, which usually means the gradient is sufficiently small, 0 otherwise}
}
\references{
Bell, A. & Sejnowski, T. An information-maximization approach to blind separation and blind deconvolution Neural computation, \emph{Neural computation}, 1995, 7, 1129-1159.

Bernaards, C. A. and Jennrich, R. I. (2005) Gradient Projection Algorithms and
  Software for Arbitrary Rotation Criteria in Factor Analysis, \emph{Educational and
  Psychological Measurement} 65, 676-696. <http://www.stat.ucla.edu/research/gpa>
}

\author{
Benjamin Risk
}
\note{
In contrast to most other ICA methods, W is not contrained to be orthogonal. 
}


%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
## Example when p > d. The MD function and amari measures 
# are not defined for M. We can compare the 
# "true W inverse", which is the mixing matrix multiplied 
# by the whitening matrix; alternatively, we can use
# multidcov::frobICA. These two approaches are
# demonstrated below:

set.seed(999)
nObs <- 1024
nComp <- 3

# simulate from gamma distributions with
# varying amounts of skewness:
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
#Define a 'true' W (uses the estimated whitening matrix):
W.true <- solve(simM\%*\%xWhitened$whitener)

estInfomax <- infomaxICA(X = xData, n.comp = nComp, whiten = TRUE, verbose = TRUE)

frobICA(estInfomax$M,simM)
library(JADE)
MD(t(estInfomax$W),t(solve(W.true)))
amari.error(t(estInfomax$W),t(solve(W.true)))
}
