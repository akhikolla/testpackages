\name{summary.re}
\docType{methods}
\alias{summary.re}
\alias{summary,re-method}

\title{Summarize an re Object}
\description{
  Summarizes and prints an \code{re} object in an informative way.
}
\usage{
 \method{summary}{re}(object, ...)
}

\arguments{
  \item{object}{Object of class \code{re}.}
  \item{...}{Additional parameters.}
}

\details{
 Summary for a \code{re} object.
}

\examples{
  X <- matrix(c(rnorm(500,2,1),rnorm(600,2,1),rnorm(400,2.2,1)),byrow=TRUE, ncol=10)
  colnames(X) <- letters[1:10]
  g <- c(rep(1,50),rep(2,60),rep(3,40))

  test <- gmw(X,g,test="kw",type="external")

  re1 <- getSigTests(test)
  
  summary(re1)
}


\author{Daniel Fischer}

\keyword{methods}
\keyword{print}
