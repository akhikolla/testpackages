\name{misoreg}
\alias{misoreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply isotonic regression to each line of a matrix
}
\description{
This function calculates the isotonic regression (assuming an upward
trend) for each line of a
given matrix with a given vector of weights. It does so by using the C
implementation of the pool adjacent violators algorithm provided in the 
package <fdrtool>, looping over the lines of the matrix in compiled C++
code. This implementation is approximately 2 orders of magnitude faster
than using apply in R. 
}
\usage{
misoreg(data, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a numeric matrix, for the lines of which we want to calculate the
isotonic regression
}
  \item{weights}{
a vector of same length as the columns of data defining the weights
}
}
\details{
misoreg takes a matrix of values for the lines of which a weighted
isotonic regression is to be computed. The weights are assumed to be
equal for each line. 
}
\value{
  \item{result}{
    matrix with lines giving the isotonic regression fit for each line
    of the input matrix assuming an upward trend}
  \item{weights}{
    vector with original weights}
}
\references{
  Barlow, R. E., Bartholomew, D. J., Bremner, J. M., and Brunk, H. D. (1972)
  "Statistical inference under order restrictions"; Wiley, London.

  Robertson, T., Wright,F. T. and Dykstra, R. L. (1988)
  "Order Restricted Statistical Inference"; Wiley, New York.
}
\author{
  Florian Klinglmueller <float_at_lefant.net>
  Part of this function is C code that has been ported by Korbinian Strimmer from
  R code originally written by Kaspar Rufibach. Many thanks also to
  Romain Francois and Dirk Eddelbuettel for helping moving the code to
  the new Rcpp API.
}
\examples{
x <- matrix(rnorm(4000),nc=4)
w <- c(3,6,3,6)/18
out <- misoreg(x,w)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

