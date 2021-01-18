\name{diabetes}
\alias{diabetes}
\docType{data}
\title{
Diabetes data
}
\description{
Subset of the diabetes data from Andrews & Herzberg (1985).
}
\usage{data(diabetes)}
\format{
  A list with components \code{x} and \code{y}.  Both components are matrices 
  with observations on different variables for the same \eqn{n = 76} persons.
  
  Component \code{x} is a matrix containing the following \eqn{p = 2} variables.
  \describe{
    \item{\code{RelativeWeight}}{relative weight.}
    \item{\code{PlasmaGlucose}}{fasting plasma glucose.}
  }
  
  Component \code{y} is a matrix containing the following \eqn{q = 3} variables.
  \describe{
    \item{\code{GlucoseIntolerance}}{glucose intolerance.}
    \item{\code{InsulinResponse}}{insulin response to oral glucose.}
    \item{\code{InsulinResistance}}{insulin resistance.}
  }
}
\source{
Andrews, D.F. and Herzberg, A.M. (1985) \emph{Data}. Springer-Verlag. Page 215.
}
\examples{
data("diabetes")
x <- diabetes$x
y <- diabetes$y

## Spearman correlation
maxCorGrid(x, y, method = "spearman")
maxCorGrid(x, y, method = "spearman", consistent = TRUE)

## Pearson correlation
maxCorGrid(x, y, method = "pearson")
}
\keyword{datasets}
