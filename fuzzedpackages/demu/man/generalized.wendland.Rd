\name{generalized.wendland}
\alias{generalized.wendland}
\title{Calculate the correlation matrix according to the generalized Wendland model.}
\usage{
generalized.wendland(l.d,theta, kap)
}
\arguments{
\item{l.d}{Current design distance matrices in distlist format}

\item{theta}{A vector of range parameters}

\item{kap}{A non-negative scalar parameter}

}
\value{
A list containing the constructed correlation matrix.
}
\description{
\code{generalized.wendland()} is a helper function that constructs a correlation matrix according to the generalized Wendland model with lengthscales given by the parameter vector \code{theta}.  When \code{kap=0} the correlation model corresponds to the Askey correlation model.  The design must have been already formated in distlist format using the function \code{makedistlist()}.
}
\examples{
library(demu)

design=matrix(runif(10,0,1),ncol=2,nrow=5)
theta=0.3
kap=3
l.d=makedistlist(design)
R=generalized.wendland(l.d,theta,kap)$R
R
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{rhomat}}
  \code{\link{matern32}}
  \code{\link{matern52}}
  \code{\link{wendland1}}
  \code{\link{wendland2}}
}
