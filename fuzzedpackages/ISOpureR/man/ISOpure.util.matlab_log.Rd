\name{ISOpure.util.matlab_log}
\alias{ISOpure.util.matlab_log}
\title{Modified logarithm function}
\description{Logarithm function that matches Matlab behaviour on negative entries (i.e. returns a 
	complex number)}
\usage{
ISOpure.util.matlab_log(x)
}
\arguments{
  \item{x}{A numeric or complex value, vector, or matrix.}
}
\value{Returns log(x) if all entries of x > 0. For complex or negative input, x, where x = a + bi, 
	the function returns log(z) = log(abs(z)) + 1i*atan2(b,a) where atan(b,a) is on the half-closed
	interval, (-pi, pi], as for the Matlab log function.}
\author{Catalina Anghel}
\examples{
ISOpure.util.matlab_log(5)
#[1] 1.609438
ISOpure.util.matlab_log(-5)
#[1] 1.609438+3.141593i
ISOpure.util.matlab_log(complex(real=3, imaginary=4))
#[1] 1.609438+0.927295i
ISOpure.util.matlab_log(c(2,3,4,-7,1))
#[1] 0.6931472+0.000000i 1.0986123+0.000000i 1.3862944+0.000000i
#[4] 1.9459101+3.141593i 0.0000000+0.000000i
}
\keyword{arith}
\keyword{complex}
