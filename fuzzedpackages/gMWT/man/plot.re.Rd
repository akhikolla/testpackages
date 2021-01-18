\name{plot.re}
\docType{methods}
\alias{plot.re}
\alias{plot,re-method}

\title{Plot an re Object}
\description{
  Plots an \code{re} object returned from the function \code{getSigTests}.
}
\usage{
 \method{plot}{re}(x,...)
}

\arguments{
  \item{x}{Object of class \code{re}.}
  \item{...}{Additional parameters.}
}

\details{
This function plots a rejection plot for an existing \code{re} object.
}

\seealso{
\code{\link[gMWT]{getSigTests}}
}

\examples{
  X <- matrix(c(rnorm(500,2,1),rnorm(600,2,1),rnorm(400,2.2,1)),byrow=TRUE, ncol=10)
  colnames(X) <- letters[1:10]
  g <- c(rep(1,50),rep(2,60),rep(3,40))

  test <- gmw(X,g,test="kw",type="external")

  re1 <- getSigTests(test)
  
  plot(re1)
}



\author{Daniel Fischer}

\keyword{methods}
\keyword{plot}
