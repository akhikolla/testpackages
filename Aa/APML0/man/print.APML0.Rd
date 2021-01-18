\name{print.APML0}
\alias{print.APML0}
\title{
Print a APML0 Object
}
\description{
Print a summary of results along the path of \code{lambda}.
}
\usage{
\method{print}{APML0}(x, digits = 4, ...)
}
\arguments{
  \item{x}{fitted APML0 object
}
  \item{digits}{significant digits in printout
}
  \item{\dots}{additional print arguments
}
}
\details{
The performed model is printed, followed by \code{fit} and \code{fit0} (if any) from a fitted \code{APML0} object.
}
\value{
The data frame above is silently returned
}
\author{
Xiang Li, Shanghong Xie, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <spiritcoke@gmail.com>
}

\seealso{
\code{\link{APML0}}
}
\examples{
###  Linear model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]\%*\%beta
y=rnorm(N,xb)

fiti2=APML0(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
fiti2
}
\keyword{Print}


