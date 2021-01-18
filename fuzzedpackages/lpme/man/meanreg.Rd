\name{meanreg}
\alias{meanreg}

\title{Local Polynomial Estimators for Nonparametric Mean Regression}
\description{
This function provides both the DFC (Delaigle, Fan, and Carroll, 2009) and HZ (Huang and Zhou, 2017) local linear estimators for solving the errors-in-variables problem. The local linear estimator in the absence of measurement error is also provided. 
}
   
\usage{
meanreg(Y, W, bw, xgrid=NULL, method="HZ", sig=NULL, error="laplace", FT_fu)
}

\arguments{
   \item{Y}{         an n by 1 response vector. }
   
   \item{W}{         an n by 1 predictor vector.} 
           
  \item{bw}{         bandwidth.} 
  
  \item{xgrid}{ the grid values to estimate the responses.}
    
    \item{method}{   the method to be used; \code{method="HZ"} returns the estimator proposed by Huang and Zhou (2017); \code{method="DFC"} returns the estimator proposed by Delaigle, Fan, and Carroll (2009); \code{method="naive"} returns the local linear estimator ignoring measurement error.}

    \item{sig}{      standard deviation of the measurement error; \code{sig=NULL} returns the naive estimators ignoring measurement error.}   

    \item{error}{    the distribution assumed for the measurement error; \code{error="laplace"} is for Laplace distribution; \code{error="normal"} is for Gaussian distribution; \code{error="user"} is for user-assumed distribution.} 
    
    \item{FT_fu}{    a function for the Fourier transform of density of error, which is required only when \code{error="user"}.} 
}

\value{
The results include the grid points \code{xgrid} for predictor and corresponding fitted responses \code{yhat}.
}

\seealso{
\code{\link{meanregbwSIMEX}}
}

\references{
Huang, X. and Zhou, H. (2017). An alternative local polynomial estimator for the error-in-variables problem. \emph{Journal of Nonparametric Statistics}, 29: 301-325.

Delaigle, A., Fan, J., and Carroll, R. (2009). A design-adaptive local polynomial estimator for the errors-in-variables problem. \emph{Journal of the American Statistical Association}, 104: 348-359.
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

## sample size:
n =100;
## Function gofx(x) to estimate
gofx  = function(x){ 1/4*x + x^2/4 }
xgrid = seq(-2, 2, 0.02)

## Generate data
sigma_e  = 0.5;
sigma_x = 1; X = rnorm(n, 0, sigma_x); 
## Sample Y
Y  = gofx(X) + rnorm(n, 0, sigma_e);
##------------------  method Based on X ---------------------------
ghat_X= meanreg(Y, X, 0.1, method="naive", xgrid=xgrid);

## reliability ratio
lambda=0.85;
sigma_u  = sqrt(1/lambda-1)*sigma_x;
print( sigma_x^2/(sigma_x^2 + sigma_u^2) );
W=X+rnorm(n,0,sigma_u);
#W=X+rlaplace(n,0,sigma_u/sqrt(2));
  
##------------------  method Based on W ---------------------------
ghat_W=meanreg(Y, W, 0.1, method="naive", xgrid=xgrid);

##------------------ JASA method -------------------------------------
h = 0.13; 
ghat_JASA=meanreg(Y, W, h, method="DFC", sig=sigma_u, 
                  error="laplace", xgrid=xgrid);

##------------------ Our method -------------------------
ghat_NEW=meanreg(Y, W, h, method="HZ", sig=sigma_u, 
                 error="laplace", xgrid=xgrid);

## plots
plot(xgrid, gofx(xgrid), "l", main="Individual", lwd="2")
lines(xgrid, ghat_NEW$yhat, lty="dashed", col="2",lwd="3")
lines(xgrid, ghat_JASA$yhat, lty="dotted", col="3",lwd="3")
lines(xgrid, ghat_X$yhat, lty="dashed", col="4",lwd="2")
lines(xgrid, ghat_W$yhat, lty="dashed", col="5",lwd="3")
}

