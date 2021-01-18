\name{RcppArmadillo-Functions}
\alias{rcpparmabasic_test}
\alias{rcpparmabasic_outerproduct}
\alias{rcpparmabasic_innerproduct}
\alias{rcpparmabasic_bothproducts}
\title{Set operation of functions in RcppArmadillo package}
\description{
  These four functions are created when
  \code{RcppArmadillo.package.skeleton()} is invoked to create a
  skeleton packages.
}
\usage{
rcpparmabasic_test()
rcpparmabasic_outerproduct(x)
rcpparmabasic_innerproduct(x)
rcpparmabasic_bothproducts(x)
}
\arguments{
  \item{x}{a numeric vector}
}
\value{
  \code{rcpparmabasic_test()} does not return a value, but displays a
  message to the console.

  \code{rcpparmabasic_outerproduct()} returns a numeric matrix computed as the
  outer (vector) product of \code{x}.

  \code{rcpparmabasic_innerproduct()} returns a double computer as the inner
  (vector) product of \code{x}.

  \code{rcpparmabasic_bothproducts()} returns a list with both the outer and
  inner products.
  
}
\details{
  These are example functions which should be largely
  self-explanatory.
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}
\examples{
  x <- sqrt(1:4)
  rcpparmabasic_innerproduct(x)
  rcpparmabasic_outerproduct(x)
}
\author{Lianfa Li}