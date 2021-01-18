\name{ci}
\alias{ci,Cpp_Clmbr-method}
\alias{ci}

\title{Confidence Interval for the Changepoint}

\description{
  Confidence interval for 'theta', the changepoint's x-coordinate.
}

\usage{
\S4method{ci}{Cpp_Clmbr}(  CL =0.95, method ="CLR" )
}

\arguments{
  \item{CL}{confidence level, between 0 and 1.}
  \item{method}{"CLR" or "AF" which stand for conditional likelihood-ratio or approximate-F, see \code{\link{sl}} for details.}
}

\details{
This subroutine scans to determine the postulate values of 'theta' that have significance level greater than 1-CL.
}

\examples{
#  Data for Patient B from Smith and Cook (1980)
y <- c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
x <- 1:10
sc <- lm.br( y ~ x )
sc$ci()
sc $ ci( 0.90 )
sc $ ci( .99, 'af' )
}

