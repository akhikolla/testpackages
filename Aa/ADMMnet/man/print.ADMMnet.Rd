\name{print.ADMMnet}
\alias{print.ADMMnet}
\title{
Print a ADMMnet Object
}
\description{
Print a summary of results along the path of \code{lambda}.
}
\usage{
\method{print}{ADMMnet}(x, digits = 4, ...)
}
\arguments{
  \item{x}{fitted ADMMnet object
}
  \item{digits}{significant digits in printout
}
  \item{\dots}{additional print arguments
}
}
\details{
The performed model is printed, followed by \code{fit} and \code{fit0} (if any) from a fitted \code{ADMMnet} object.
}
\value{
The data frame above is silently returned
}
\author{
Xiang Li, Shanghong Xie, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xl2473@cumc.columbia.edu>
}

\seealso{
\code{\link{ADMMnet}}
}
\examples{
###  Linear model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]%*%beta
y=rnorm(N,xb)

fiti=ADMMnet(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
fiti
}
\keyword{Print}


