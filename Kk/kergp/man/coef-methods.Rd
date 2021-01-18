\name{coef-methods}
%%\alias{coef-methods}
\docType{methods}
\alias{coef,methods}
\alias{coef,covMan-method}
\alias{coef,covTS-method}

\title{
  
  Extract Coefficients of a Covariance Kernel Object as Vector, List or
  Matrix

}

\description{
  
  Extract some of or all the coefficients of a covariance kernel object as
  vector, list or matrix.

}
\usage{
\S4method{coef}{covMan}(object)

\S4method{coef}{covTS}(object, type = "all", as = "vector")
}
\arguments{
  \item{object}{

    An object representing a covariance kernel, the coefficient of which
    will be extracted.

  }
  \item{type}{
    
    Character string or vector specifying which type(s) of coefficients in the
    structure will be extracted. Can be \code{"all"} (all coefficients
    are extracted) or any parameter name(s) of the corresponding kernel.

  }
  \item{as}{

    Character string specifying the output structure to be used.  The
    default is \code{"vector"}, leading to a numeric vector. Using
    \code{"list"} one gets a list of numeric vectors, one by kernel
    parameter. Finally, using \code{"matrix"} one gets a matrix with one
    row by input (or dimension) and one column by (selected) kernel
    parameter.

  }
  %% \item{\dots}{Not used yet.}
}
\value{
  
  A numeric vector of coefficients or a structure as specified by
  \code{as} containing the coefficients selected by \code{type}.

}
\seealso{

  The \code{\link{coef<-}} replacement method which takes a vector of
  replacement values.

}
\examples{
d <- 3
myCov1 <- covTS(d = d, kernel = "k1Exp", dep = c(range = "input"),
                value = c(range = 1.1))
myCov1
## versatile 'coef' method
coef(myCov1)
coef(myCov1, as = "matrix")
coef(myCov1, as = "list")
coef(myCov1, as = "matrix", type = "range")
coef(myCov1) <- c(0.2, 0.3, 0.4, 4, 16, 25)
coef(myCov1, as = "matrix")
}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
