\name{gauss.cor.test-methods}
\docType{methods}
\alias{gauss.cor.test}
\alias{gauss.cor.test-methods}
\alias{gauss.cor.test,numeric,numeric-method}
\alias{gauss.cor.test,formula,data.frame-method}
\title{Gaussian rank correlation test}
\description{
  Methods performing a Gaussian rank correlation test
}
\usage{
\S4method{gauss.cor.test}{numeric,numeric}(x, y, ...)
\S4method{gauss.cor.test}{formula,data.frame}(x, y, na.action, ...)
}
\arguments{
  \item{x}{a numeric vector or a formula; compulsory argument}
  \item{y}{compulsory argument; if \code{x} is a vector, \code{y}
    must be vector of the same length as \code{x}. If \code{x} is
    a formula, \code{y} must be a data frame.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}'s.  Defaults to \code{getOption("na.action")}.}
  \item{...}{all parameters specified are forwarded internally to the
    method \code{\link{cor.test}}, in particular, the \code{alternative}
    parameter.}
}
\details{If called for numeric vectors, \code{gauss.cor.test} performs the
  Gaussian gamma rank correlation test for \code{x} and \code{y}. This
  is done by simply performing a Pearson correlation test on the normal
  scores of the data.

  If \code{gauss.cor.test} is called for a formula \code{x} and a data
  frame \code{y}, then the method checks whether the formula \code{x} correctly
  extracts two columns from \code{y} (see examples below). If so, the
  two columns are extracted and the Gaussian gamma rank correlation test
  is applied to them according to the specified parameters.
}
\value{
  Upon successful completion, the function returns a list of class
  \code{htest} containing the results (see \code{\link{cor.test}}).
}
\author{Ulrich Bodenhofer
\email{rococo@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/rococo/}

K. Boudt, J. Cornelissen and C. Croux (2012).
The Gaussian rank correlation estimator: robustness properties.
\emph{Stat. Comput.} \bold{22}(2):471-483.}
\seealso{\code{\link{gauss.cor}}}
\examples{
## create data
f <- function(x) ifelse(x > 0.9, x - 0.9, ifelse(x < -0.9, x + 0.9, 0))
x <- rnorm(25)
y <- f(x) + rnorm(25, sd=0.1)

## perform correlation tests
gauss.cor.test(x, y, alternative="greater")

## the formula variant
require(datasets)
data(iris)
gauss.cor.test(~ Petal.Width + Petal.Length, iris,
               alternative="two.sided")
}
\keyword{htest}
\keyword{methods}
