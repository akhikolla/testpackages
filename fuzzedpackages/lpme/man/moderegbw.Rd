\name{moderegbw}
\alias{moderegbw}

\title{Cross-Validation Bandwidth Selector for Nonparametric Mode Regression}
\description{
This function selects the bandwidth (Zhou and Huang, 2018) for the local polynomial estimators for nonparametric modal regression in the absence of measurement error.}
   
\usage{
moderegbw(Y, X, method="CV-density", p.order=0, h1=NULL, h2=NULL, nstart = 4,
          xinterval = quantile(X, probs=c(0.025, 0.975), names = FALSE),
          df=5, ncomp=5, nboot=5)
}

\arguments{
   \item{Y}{         an n by 1 response vector. }
   
   \item{X}{         an n by 1 predictor vector.} 
   
   \item{method}{    \code{method="CV-density"} is for density-based CV; \code{method="CV-mode"} is for mode-based CV; \code{metohd="bootstrap"} is for a bootstap method; see Zhou and Huang (2017) for details. For non-measurement error models, \code{method="CV-mode"} is recommended.}
   
   \item{p.order}{   the order of polynomial, up to 1; \code{p.order=0} returns local constant estimators and \code{p.order=1} returns local linear estimators. }
   
   \item{h1}{       bandwidth vector for h1; default is \code{NULL}, and h1 is chosen automatically. See Zhou and Huang (2017) for details. It is recommended to carefully specify a fine grid for h1. }
    
    \item{h2}{      bandwidth vector for h2; default is \code{NULL}, and h2 is chosen automatically. See Zhou and Huang (2017) for details. It is recommended to carefully specify a fine grid for h2.}
    
    \item{nstart}{    the starting number of modes for each grid value.}
    
    \item{xinterval}{ the interval within which the modes will be estimated.}

    \item{df}{        the degrees of freedom of splines used in the mixture normal regression for bootstrap method.} 
    
    \item{ncomp}{     the number of components used in the mixture normal regression for bootstrap method.} 
    
    \item{nboot}{     the number of bootstrap samples.}
}

\value{
The results include the bandwidth \code{bw}.
}

\seealso{
\code{\link{moderegbwSIMEX},\link{modereg}}
}

\references{
Zhou, H. and Huang, X. (2018). Bandwidth selection for nonparametric modal regression. \emph{Communications in Statistics - Simulation and Computation}, in press.
}


\author{
Haiming Zhou <zhouh@niu.edu> and Xianzheng Huang <huang@stat.sc.edu>
}

\examples{
rm(list=ls())
library(lpme)
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

## mode estimates
h1ref = c(1.06*sd(X)*n^(-0.2));
h2ref = c(1.06*sd(Y)*n^(-0.2));
## In practice moer fine grids are desired. 
hx = seq(h1ref*0.2, h1ref*1.5, length.out = 10); 
hy = seq(h2ref*0.8, h2ref, length.out = 2); 
hhxy = moderegbw(Y, X, method="CV-mode", p.order=0,
                 h1=hx, h2=hy)$bw;
fit = modereg(Y, X, xgrid=xgrid, bw=hhxy, p.order=0, PLOT=TRUE);

## Plot
plot(xgrid, gofx1(xgrid), "l", lwd="2", ylim=c(-9,7), xlim=c(-2,2));
lines(xgrid, gofx2(xgrid), "l", lwd="2");
points(rep(fit$xgrid,fit$x.num), fit$mode, col="3",lwd="2")
}

