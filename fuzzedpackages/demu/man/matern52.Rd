\name{matern52}
\alias{matern52}
\title{Calculate the correlation matrix according to the Matern model with \eqn{\nu=5/2}.}
\usage{
matern52(l.d,theta)
}
\arguments{
\item{l.d}{Current design distance matrices in distlist format}

\item{theta}{A vector of range parameters}

}
\value{
A list containing the constructed correlation matrix.
}
\description{
\code{matern52()} is a helper function that constructs a correlation matrix according to the Matern model with parameter \eqn{\nu=5/2} and lengthscales given by the parameter vector \code{theta}.  The design must have been already formated in distlist format using the function \code{makedistlist()}.
}
\examples{
library(demu)

design=matrix(runif(10,0,1),ncol=2,nrow=5)
theta=rep(0.2,2)
l.d=makedistlist(design)
R=matern52(l.d,theta)$R
R
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{rhomat}}
  \code{\link{matern32}}
  \code{\link{wendland1}}
  \code{\link{wendland2}}
  \code{\link{generalized.wendland}}
}
