\name{wendland1}
\alias{wendland1}
\title{Calculate the correlation matrix according to the Wendland1 model.}
\usage{
wendland1(l.d,theta)
}
\arguments{
\item{l.d}{Current design distance matrices in distlist format}

\item{theta}{A vector of range parameters}

}
\value{
A list containing the constructed correlation matrix.
}
\description{
\code{wendland1()} is a helper function that constructs a correlation matrix according to the Wendland 1 model with lengthscales given by the parameter vector \code{theta}.  The design must have been already formated in distlist format using the function \code{makedistlist()}.
}
\examples{
library(demu)

design=matrix(runif(10,0,1),ncol=2,nrow=5)
theta=rep(0.3,2)
l.d=makedistlist(design)
R=wendland1(l.d,theta)$R
R
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{rhomat}}
  \code{\link{matern32}}
  \code{\link{matern52}}
  \code{\link{wendland2}}
  \code{\link{generalized.wendland}}
}
