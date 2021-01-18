\name{meanregbwSIMEX}
\alias{meanregbwSIMEX}

\title{Cross-Validation Bandwidth Selector Using SIMEX for Nonparametric Mean Regression}
\description{
This function selects the bandwidth for both the DFC (Delaigle, Fan, and Carroll, 2009) and HZ (Huang and Zhou, 2017) estimators. 
}
   
\usage{
meanregbwSIMEX(Y, W, method="HZ", sig, error="laplace", k_fold=5, B=10,
               h1=NULL, h2=NULL, length.h=10, lconst=0.5, rconst=2, Wdiff=NULL)
}

\arguments{
   \item{Y}{         an n by 1 response vector. }
   
   \item{W}{         an n by 1 predictor vector.} 
    
    \item{method}{   the method to be used; \code{method="HZ"} uses the estimator proposed by Huang and Zhou (2017); \code{method="DFC"} uses the estimator proposed by Delaigle, Fan, and Carroll (2009). It currently does not support bandwidth selection for naive estimators.}

    \item{sig}{      standard deviation of the measurement error.}   

    \item{error}{    the distribution assumed for the measurement error; \code{error="laplace"} is for Laplace distribution; \code{error="normal"} is for Gaussian distribution. It currently does not support user-assumed distribution.} 
                    
    \item{k_fold}{   gives fold of cross-validation to be used; default is 2.}
    
    \item{B}{        total number of cross-validation criteria to average over; defualt is 10.}
    
    \item{h1}{       bandwidth vector for the first level error contamination; default is \code{NULL}, and h1 is chosen automatically. See Huang and Zhou (2017) for details. }
    
    \item{h2}{      bandwidth vector for the second level error contamination; defualt is \code{NULL}, and h2 is chosen automatically. See Huang and Zhou (2017) for details.}
    
    \item{length.h}{ number of grid points for each of h1 and h2; default is 10.}
    
    \item{lconst, rconst}{ used to control the searching windows for bandwidths h1 and h2. For example, \code{seq(bw1*lconst,bw1*rconst,length.out=length.h)} is used to obtain bandwidth grid points for h1, where bw1 is an initial bandwidth; see Huang and Zhou (2017) for details of finding the initial bandwith.}
    
    \item{Wdiff}{    an n by 1 vector of \code{(W1-W2)/2}, where W1, W2 are two replicated measurements; default is \code{NULL}, which indicates that the errors are generated from the assumed error distribution, otherwise, the errors are generated from Wdiff with replacement.}
}

\value{
  The results include the bandwidth \code{bw}.
}

\seealso{
\code{\link{meanreg}}
}


\references{
Huang, X. and Zhou, H. (2017). An alternative local polynomial estimator for the error-in-variables problem. \emph{Journal of Nonparametric Statistics}, 29: 301-325.

Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice in errors-in-variables problems. \emph{Journal of the American Statistical Association}, 103: 280-287.
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
location-sign(rrrr-0.5)*scale*(log(2)+ifelse(rrrr<0.5, log(rrrr), log1p(-rrrr)))
}

## sample size:
n =100;
## Function gofx(x) to estimate
gofx  = function(x){ 2*x*exp(-10*x^4/81) }

## Generate data
sigma_e  = 0.2;
sigma_x = 1; X = rnorm(n, 0, sigma_x); 
## Sample Y
Y  = gofx(X) + rnorm(n, 0, sigma_e);
## reliability ratio
lambda=0.85;
sigma_u  = sqrt(1/lambda-1)*sigma_x;
print( sigma_x^2/(sigma_x^2 + sigma_u^2) );
#W=X+rnorm(n,0,sigma_u);
W=X+rlap(n,0,sigma_u/sqrt(2));
  
#### SIMEX
#**Note: larger values for B and length.h are needed for accurate estimates.
#**e.g., k_fold=5, B=10, length.h=10 will be generally good. 
hwNEW = meanregbwSIMEX(Y, W, method="HZ", sig=sigma_u, error="laplace", k_fold=2, 
                        B=1, length.h=1)$bw
ghat_NEW = meanreg(Y, W, hwNEW , method="HZ", sig=sigma_u, error="laplace");

## plots
x = ghat_NEW$xgrid;
plot(x, gofx(x), "l", main="Individual", lwd="2")
lines(ghat_NEW$xgrid, ghat_NEW$yhat, lty="dashed", col="2",lwd="3")
}

