\name{densityreg}
\alias{densityreg}

\title{Conditional Density Estimation with Covariate Measurement Error}
\description{
This function provides nonparametric estimators (Huang and Zhou, 2018) for the density of a response conditioning on an error-prone covariate. The corresponding estimators in the absence of measurement error are also provided. 
}
   
\usage{
densityreg(Y, W, bw, xgrid=NULL, ygrid = NULL, sig=NULL, K1 = "Gauss", 
           K2 = "Gauss", mean.estimate = NULL, spline.df = 5)
}

\arguments{
   \item{Y}{         an n by 1 response vector. }
   
   \item{W}{         an n by 1 predictor vector.} 
   
   \item{bw}{        bandwidth.} 
   
   \item{xgrid}{     the grid values in x-axis to estimate the conditional density.}
   
   \item{ygrid}{     the grid values in y-axis to estimate the conditional density.}
   
   \item{sig}{      standard deviation of the measurement error; \code{sig=NULL} returns the naive estimators ignoring measurement error.} 
   
   \item{K1}{      kernel function for X; choices include \code{"Gauss"} and \code{SecOrder}.} 
   
   \item{K2}{      kernel function for Y; choices include \code{"Gauss"} and \code{SecOrder}.} 
   
   \item{mean.estimate}{    method to estimate the naive mean function of Y given X; choices include \code{"spline"} and \code{"kernel"}. If \code{NULL}, the methods 1 and 3 in the reference below are considered.}
   
   \item{spline.df}{   the degrees of freedom for B-splines when \code{mean.estimate="spline"}. }
}

\value{
The results include the grid points \code{xgrid} for X and \code{ygrid} for Y, the fitted density values \code{fitxy}.
}

\seealso{
\code{\link{densityregbw}}
}

\references{
Huang, X. and Zhou, H. (2018). Conditional density estimation with covariate measurement error. Submitted.
}


\author{
Haiming Zhou <zhouh@niu.edu> and Xianzheng Huang <huang@stat.sc.edu>
}

\examples{
#############################################
## X - True covariates
## W - Observed covariates
## Y - individual response
rm(list=ls())
library(lpme)
## generate laplace
rlap=function (use.n, location = 0, scale = 1) 
{
  location <- rep(location, length.out = use.n)
  scale <- rep(scale, length.out = use.n)
  rrrr <- runif(use.n)
  location - sign(rrrr - 0.5) * scale * (log(2) + ifelse(rrrr < 0.5, log(rrrr), log1p(-rrrr)))
}
## sample size:
n =100;
## Function f(y|x) to estimate#
mofx = function(x){ x }
sofx = function(x){ 0.5 } 
fy_x=function(y,x) dnorm(y, mofx(x), sofx(x));
## Generate data
sigma_x = 1; X = rnorm(n, 0, sigma_x); 
## Sample Y
Y = rep(0, n);
for(i in 1:n){
  Y[i] = mofx(X[i]) + rnorm(1, 0, sofx(X[i]));
}

## reliability ratio
lambda=0.7;
sigma_u  = sqrt(1/lambda-1)*sigma_x;
#W=X+rnorm(n,0,sigma_u);
W=X+rlap(n,0,sigma_u/sqrt(2));

##----- naive kernel density estimate  ---- 
fitbw = densityregbw(Y, W, K1 = "Gauss", K2 = "Gauss")
fhat1 = densityreg(Y, W, bw=fitbw$bw, K1 = "Gauss", K2 = "Gauss");

###----- naive kernel density estimate with mean adjustment  ---- 
#fitbw = densityregbw(Y, W, K1 = "Gauss", K2 = "Gauss", 
#                     mean.estimate = "kernel")
#fhat2 = densityreg(Y, W, bw=fitbw$bw, K1 = "Gauss", K2 = "Gauss", 
#                   mean.estimate = "kernel");

##----- proposed method without mean adjustment  ---- 
fitbw = densityregbw(Y, W, sig=sigma_u, K1="SecOrder", K2="Gauss")
fhat3 = densityreg(Y, W, bw=fitbw$bw, sig=sigma_u, K1="SecOrder", K2="Gauss");

##----- proposed method wit mean adjustment  ---- 
#fitbw = densityregbw(Y, W, sig=sigma_u, K1="SecOrder", K2="SecOrder", 
#                     mean.estimate = "kernel")
#fhat4 = densityreg(Y, W, bw=fitbw$bw, sig=sigma_u, K1="SecOrder", K2="SecOrder", 
#                   mean.estimate = "kernel");

par(mfrow=c(2,2))
plot(W,Y, col=2)
points(X,Y)
x0 = seq(0, 1, length.out = 3)
for(i in 1:length(x0)){
  plot(fhat1$ygrid, fy_x(fhat1$ygrid, x0[i]), "l", lwd="2", xlab = "y", ylab = "density",
       main = paste("Conditional density at x=", x0[i], sep=""), ylim=c(0,1.5));
  indx = which.min(abs(fhat1$xgrid-x0[i])) ## the index of xgrid that is closest to x0[i].
  lines(fhat1$ygrid, fhat1$fitxy[indx,], lwd=3, lty=2, col=1)
  #lines(fhat2$ygrid, fhat2$fitxy[indx,], lwd=3, lty=2, col=2)
  lines(fhat3$ygrid, fhat3$fitxy[indx,], lwd=3, lty=2, col=3)
  #lines(fhat4$ygrid, fhat4$fitxy[indx,], lwd=3, lty=2, col=4)
}

}

