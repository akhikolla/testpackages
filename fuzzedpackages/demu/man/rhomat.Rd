\name{rhomat}
\alias{rhomat}
\title{Calculate the correlation matrix according to the squared exponential family of models.}
\usage{
rhomat(l.d,rho,alpha=2)
}
\arguments{
\item{l.d}{Current design distance matrices in distlist format}

\item{rho}{A vector of correlation parameters taking on values in [0,1)}

\item{alpha}{Exponent parameter}
}
\value{
A list containing the constructed correlation matrix.
}
\description{
\code{rhomat()} is a helper function that constructs a correlation matrix according to the squared exponential model with parameterized by correlation parameters \code{rho} taking values in [0,1) and the exponent parameter \code{alpha}.  The default of \code{alpha=2} results in the Gaussian correlation while selecting \code{alpha=1} corresponds to the Exponential correlation model.  The design must have been already formated in distlist format using the function \code{makedistlist()}.
}
\examples{
library(demu)

design=matrix(runif(10,0,1),ncol=2,nrow=5)
rho=rep(0.01,2)
l.d=makedistlist(design)
R=rhomat(l.d,rho)$R
R
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{matern52}}
  \code{\link{wendland1}}
  \code{\link{wendland2}}
  \code{\link{generalized.wendland}}
}
