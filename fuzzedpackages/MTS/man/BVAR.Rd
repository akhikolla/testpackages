\name{BVAR}
\alias{BVAR}
\title{Bayesian Vector Autoregression
}
\description{Estimate a VAR(p) model using Bayesian approach, including the 
use of Minnesota prior
}
\usage{
BVAR(z,p=1,C,V0,n0=5,Phi0=NULL,include.mean=T)
}
\arguments{
  \item{z}{A matrix of vector time series, each column represents a series.
}
  \item{p}{
The AR order. Default is p=1.
}
  \item{C}{
The precision matrix of the coefficient matrix. With constant, 
the dimension of C is (kp+1)-by-(kp+1). The covariance matrix of the 
prior for the parameter vec(Beta) is Kronecker(Sigma_a,C-inverse).
}
  \item{V0}{
A k-by-k covariance matrix to be used as prior for the Sigma_a matrix
}
  \item{n0}{
The degrees of freedom used for prior of the Sigma_a matrix, the covariance matrix of the innovations. Default is n0=5.
}
  \item{Phi0}{
The prior mean for the parameters. Default is set to NULL, implying that the 
prior means are zero.
}
  \item{include.mean}{
A logical switch controls the constant term in the VAR model. Default is to include the constant term.
}
}
\details{
for a given prior, the program provide the posterior estimates of a 
VAR(p) model.
}
\value{
\item{est}{Posterior means of the parameters}
\item{Sigma}{Residual covariance matrix}
}
\references{Tsay (2014, Chapter 2).
}
\author{Ruey S. Tsay
}
\examples{
data("mts-examples",package="MTS")
z=log(qgdp[,3:5])
zt=diffM(z)*100
C=0.1*diag(rep(1,7))
V0=diag(rep(1,3))
BVAR(zt,p=2,C,V0)
}




