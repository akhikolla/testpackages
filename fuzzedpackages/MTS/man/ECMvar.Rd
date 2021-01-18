\name{ECMvar}
\alias{ECMvar}
\title{Error-Correction VAR Models
}
\description{Performs estimation of an Error-Correction VAR(p) model using the
Quasi Maximum Likelihood Method.
}
\usage{
ECMvar(x, p, ibeta, include.const = FALSE, fixed = NULL,
                 alpha = NULL, se.alpha = NULL, se.beta = NULL, phip =
                 NULL, se.phip = NULL)
}
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional co-integrated VAR process
}
  \item{p}{VAR order
}
  \item{ibeta}{Initial estimate of the co-integrating matrix. The number of 
columns of ibeta is the number of co-integrating series
}
  \item{include.const}{A logical switch to include a constant term in the model. The default is no constant
}
  \item{fixed}{A logical matrix to set zero parameter constraints.
}
  \item{alpha}{Initial estimate of alpha, if any}
  \item{se.alpha}{Initial estimate of the standard error of alpha, if any}
  \item{se.beta}{Initial estimate of the standard error of beta, if any}
  \item{phip}{Initial estimate of the VAR coefficients, if any}
  \item{se.phip}{Initial estimate of the standard error of the VAR 
coefficients, if any}
}
\value{
\item{data}{The vector time series}
\item{ncoint}{The number of co-integrating series}
\item{arorder}{VAR order}
\item{include.const}{Logical switch to include constant}
\item{alpha,se.alpha}{Estimates and their standard errors of the alpha matrix}
\item{beta,se.beta}{Estimates and their standard errors of the beta matrix}
\item{aic,bic}{Information criteria of the fitted model}
\item{residuals}{The residual series}
\item{Sigma}{Residual covariance matrix}
\item{Phip,se.Phip}{Estimates and their standard errors of VAR coefficients}
}
\references{Tsay (2014, Chapter 5)
}
\author{Ruey S. Tsay
}
\seealso{ECMvar1
}
\examples{
phi=matrix(c(0.5,-0.25,-1.0,0.5),2,2); theta=matrix(c(0.2,-0.1,-0.4,0.2),2,2)
Sig=diag(2)
mm=VARMAsim(300,arlags=c(1),malags=c(1),phi=phi,theta=theta,sigma=Sig)
zt=mm$series[,c(2,1)]
beta=matrix(c(1,0.5),2,1)
m1=ECMvar(zt,3,ibeta=beta)
names(m1)
}




