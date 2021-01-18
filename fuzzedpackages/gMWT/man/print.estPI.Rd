\name{print.estPI}
\docType{methods}
\alias{print.estPI}
\alias{print-method}
\alias{print,estPI-method}

\title{Print an estPI Object}
\description{
  Prints an \code{estPI} object.
}
\usage{
 \method{print}{estPI}(x,...)
}

\arguments{
  \item{x}{Object of class \code{estPI}.}
  \item{...}{Additional parameters.}
}

\details{
The print function displays the matrix of probabilistic indices.
}

\examples{
  X <- c(sample(15))
  g <- c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,4)
  res <- estPI(X,g,type="single")
  print(res)
}



\author{Daniel Fischer}

\keyword{methods}
\keyword{print}
