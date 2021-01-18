\name{moderegbwSIMEX}
\alias{moderegbwSIMEX}

\title{Cross-Validation Bandwidth Selector Using SIMEX for Nonparametric Mode Regression}
\description{
This function selects the bandwidth (Zhou and Huang, 2016) for the local polynomial estimators for nonparametric modal regression in the presence of measurement error.
}

\usage{
moderegbwSIMEX(Y, W, method="CV-density", p.order=0, sig, B=5, 
               h1=NULL, h2=NULL, length.h=10, CIregion=0.95, nstart=4)
}

\arguments{
   \item{Y}{         an n by 1 response vector. }
   
   \item{W}{         an n by 1 predictor vector.} 
    
    \item{method}{   \code{method="CV-density"} is for density-based CV; \code{method="CV-mode"} is for mode-based CV; see Zhou and Huang (2016) for details. For measurement error models, \code{method="CV-density"} is recommended.}
    
    \item{p.order}{   the order of polynomial, up to 1; \code{p.order=0} returns local constant estimators and \code{p.order=1} returns local linear estimators. }
    
    \item{sig}{      standard deviation of the measurement error.}   
    
    \item{B}{        total number of cross-validation criteria to average over; defualt is 5.}
    
   \item{h1}{       bandwidth vector for h1; default is \code{NULL}, and h1 is chosen automatically. See Zhou and Huang (2017) for details. It is recommended to carefully specify a fine grid for h1. }
    
    \item{h2}{      bandwidth vector for h2; default is \code{NULL}, and h2 is chosen automatically. See Zhou and Huang (2017) for details. It is recommended to carefully specify a fine grid for h2.}
    
    \item{length.h}{ number of grid points for each of h1 and h2; default is 10.}
    
    \item{CIregion}{ used to determine the interval on which the mode is estimated.}
    
    \item{nstart}{    the starting number of modes for each grid value.}
}

\value{
  The results include the bandwidth \code{bw}.
}

\seealso{
\code{\link{modereg}, \link{moderegbw}}
}


\references{
Zhou. H. and Huang, X. (2016). Nonparametric modal regression in the presence of measurement error. \emph{Electronic Journal of Statistics}, 10: 3579-3620.
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
#                      sig=sigma_u, B=5, length.h=10)$bw;
fit = modereg(Y, W, xgrid=xgrid, bw=hhxy, sig=sigma_u, p.order=0,
              PLOT=TRUE);

## Plot
plot(xgrid, gofx1(xgrid), "l", lwd="2", ylim=c(-9,7), xlim=c(-2,2));
lines(xgrid, gofx2(xgrid), "l", lwd="2");
points(rep(fit$xgrid,fit$x.num), fit$mode, col="3",lwd="2")
}

