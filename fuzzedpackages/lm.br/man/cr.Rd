\name{cr}
\alias{cr,Cpp_Clmbr-method}
\alias{cr}

\title{Confidence Region for the Changepoint}

\description{
  Joint confidence region for ( theta, alpha ),  the changepoint's (x,y)-coordinates.
}

\usage{
\S4method{cr}{Cpp_Clmbr}(  CL =0.95 ,  method ="CLR",  incr,  output ="G"  )
}

\arguments{
  \item{CL}{confidence level, between 0 and 1.}
  \item{method}{"CLR" or "AF" which stand for conditional likelihood-ratio  or  approximate-F (rapid), see \code{\link{sl}} for details.}
  \item{incr}{increment of theta values for the confidence region's boundary-points.}
  \item{output}{"G", "T" or "V" which stand for graph, text  or  value.}
}

\details{
This subroutine scans to determine the postulate values of (theta, alpha) that have significance level greater than 1-CL.  It scans first along the (theta, alpha-MLE) ridge to determine the 'theta' boundary limits.
}

\value{
If 'output' is "G" or "T" then 'cr' graphs or prints-out the confidence region but does not return a value.  If 'output' is "V" then 'cr' returns an N x 3 matrix of boundary points  ( theta, min-alpha, max-alpha ).
}


\examples{
#  A quick example
y <- c( 2, 0, 2.001, 4, 6 )
x <- 1:5
t <- lm.br( y ~ x )
t $ cr()
t$cr( .9, 'af', incr = 0.1, out='t' )
}

