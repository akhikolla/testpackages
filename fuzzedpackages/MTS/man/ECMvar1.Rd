\name{ECMvar1}
\alias{ECMvar1}
\title{Error-Correction VAR Model 1
}
\description{Perform least-squares estimation of an ECM VAR(p) model with 
known co-integrating processes
}
\usage{
ECMvar1(x, p, wt, include.const = FALSE, fixed = NULL, output = TRUE)
}
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional co-integrated VAR process
}
  \item{p}{VAR order
}
  \item{wt}{A T-by-m data matrix of m-dimensional co-integrated process
}
  \item{include.const}{A logical switch to include a constant term. 
Default is no constant.
}
  \item{fixed}{A logical matrix to set zero parameter constraints
}
  \item{output}{A logical switch to control output
}
}
\value{
\item{data}{The vector time series}
\item{wt}{The co-integrated series}
\item{arorder}{VAR order}
\item{include.const}{Logical switch to include constant}
\item{coef}{Parameter estimates}
\item{aic,bic}{Information criteria of the fitted model}
\item{residuals}{The residual series}
\item{Sigma}{Residual covariance matrix}
}
\references{Tsay (2014, Chapter 5). Multivariate Time Series Analysis 
with R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{ECMvar
}
\examples{
phi=matrix(c(0.5,-0.25,-1.0,0.5),2,2); theta=matrix(c(0.2,-0.1,-0.4,0.2),2,2)
Sig=diag(2)
mm=VARMAsim(300,arlags=c(1),malags=c(1),phi=phi,theta=theta,sigma=Sig)
zt=mm$series
wt=0.5*zt[,1]+zt[,2]
m1=ECMvar1(zt,3,wt)
names(m1)
}




