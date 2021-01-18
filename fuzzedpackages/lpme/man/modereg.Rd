\name{modereg}
\alias{modereg}

\title{Nonparametric Estimators for Nonparametric Mode Regression}
\description{
This function provides the nonparametric estimators (Zhou and Huang, 2016; Zhou and Huang, 2018) for nonparametric modal regression. The corresponding estimators in the absence of measurement error are also provided. 
}
   
\usage{
modereg(Y, W, bw, xgrid=NULL, sig=NULL, nstart=4, p.order=0, maxiter = 500, 
        tol=.Machine$double.eps^0.25, mesh=NULL, PLOT=FALSE, ...)
}

\arguments{
   \item{Y}{         an n by 1 response vector. }
   
   \item{W}{         an n by 1 predictor vector.} 
   
   \item{bw}{        bandwidth.} 
   
   \item{xgrid}{     the grid values to estimate the responses.}
   
   \item{sig}{      standard deviation of the measurement error; \code{sig=NULL} returns the naive estimators ignoring measurement error.} 
   
   \item{nstart}{    the starting number of modes for each grid value.}
   
   \item{p.order}{   the order of polynomial, up to 1; \code{p.order=0} returns local constant estimators and \code{p.order=1} returns local linear estimators. }
    
    \item{maxiter}{  the maximum number of iterations performed for the mean shift algorithm if not converage.} 
    
    \item{tol}{      the deisered accurary (convergence tolerrance).}
    
    \item{mesh}{     a matrix of initial mode points, where each row corresponds a mode in \code{(x,y)} coordinate; if \code{mesh=NULL}, it will be chosen automatically according to \code{xgrid} and \code{nstart}.}
    
    \item{PLOT}{     a logical value indicating whether the estimated modes will be plotted.}
   
   \item{...}{       further arguments to be passed to or from other methods.}

}

\value{
The results include the grid points \code{xgrid} for predictor, the number of modes for each grid \code{x.num}, the initial mesh points \code{mesh}, and corresponding fitted modes \code{mode}.
}

\seealso{
\code{\link{moderegbwSIMEX},\link{moderegbw}}
}

\references{
Zhou. H. and Huang, X. (2016). Nonparametric modal regression in the presence of measurement error. \emph{Electronic Journal of Statistics}, 10: 3579-3620.

Zhou, H. and Huang, X. (2018). Bandwidth selection for nonparametric modal regression. \emph{Communications in Statistics - Simulation and Computation}, in press.
}


\author{
Haiming Zhou <zhouh@niu.edu> and Xianzheng Huang <huang@stat.sc.edu>
}

\examples{
rm(list=ls())
library(lpme)

rlaplace=function (use.n, location = 0, scale = 1) 
{
  location <- rep(location, length.out = use.n)
  scale <- rep(scale, length.out = use.n)
  rrrr <- runif(use.n)
  location - sign(rrrr - 0.5) * scale * 
    (log(2) + ifelse(rrrr < 0.5, log(rrrr), log1p(-rrrr)))
}

## sample size:
n =100;
## Function m(x) to estimate#
gofx1  = function(x){ (x+x^2) }
gofx2  = function(x){ (x+x^2)-6 }
xgrid	= seq(-2, 2, length.out=100); 
ngrid	= length(xgrid)

## Sample X
X = rnorm(n, 0, 1); sigma_x=1;
## Sample Y
Y = rep(0, n);
U = runif(n);
for(i in 1:n){
  if(U[i]<0.5){
    Y[i] = rnorm(1, gofx1(X[i]), 1);
  }else{
    Y[i] = rnorm(1, gofx2(X[i]), 1);
  }
}
## reliability ratio
lambda=0.9;
sigma_u  = sqrt(1/lambda-1)*sigma_x;
W=X+rlaplace(n,0,sigma_u/sqrt(2));

## mode estimates
hhxy = c(0.15, 1)
## Note you needs to use the following code to calculate bandwidth
## It is not run here due to the time constrain of runing examples. 
#hhxy = moderegbwSIMEX(Y, W, method="CV-density", p.order=0, 
#                       sig=sigma_u, B=5, length.h=10)$bw;
fit = modereg(Y, W, xgrid=xgrid, bw=hhxy, sig=sigma_u, p.order=0,
              PLOT=TRUE);

## Plot
plot(xgrid, gofx1(xgrid), "l", lwd="2", ylim=c(-9,7), xlim=c(-2,2));
lines(xgrid, gofx2(xgrid), "l", lwd="2");
points(rep(fit$xgrid,fit$x.num), fit$mode, col="3",lwd="2")
}

