\name{print.re}
\docType{methods}
\alias{print.re}
\alias{print,re-method}

\title{Print an re Object}
\description{
  Prints an \code{re} object.
}
\usage{
 \method{print}{re}(x,...)
}

\arguments{
  \item{x}{Object of class \code{re}.}
  \item{...}{Additional parameters.}
}

\details{
The print function displays in a list the significant variables for a chosen alpha and
the optimal alpha in terms of the best ratio between expected and observed significant
tests. For more details see also \code{\link[gMWT]{rejectionPlot}}.
}

\seealso{
\code{\link{rejectionPlot}}
}

\examples{
  X <- matrix(c(rnorm(500,2,1),rnorm(600,2,1),rnorm(400,2.2,1)),byrow=TRUE, ncol=10)
  colnames(X) <- letters[1:10]
  g <- c(rep(1,50),rep(2,60),rep(3,40))

  test <- gmw(X,g,test="kw",type="external")

  re1 <- getSigTests(test)
  
  re1
}



\author{Daniel Fischer}

\keyword{methods}
\keyword{print}
