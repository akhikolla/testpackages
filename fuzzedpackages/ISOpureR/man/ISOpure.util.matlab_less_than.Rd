\name{ISOpure.util.matlab_less_than}
\alias{ISOpure.util.matlab_less_than}
\title{Less than operator}
\description{Less than function that matches Matlab behaviour when one of the arguments is NA (i.e. returns FALSE instead of NA)}
\usage{
ISOpure.util.matlab_less_than(a, b)
}
\arguments{
  \item{a}{A numeric value (including Inf) or NA}
  \item{b}{A numeric value (including Inf) or NA}
}
\value{Logical: TRUE if a < b, FALSE if a >= b OR if one of a, b is NA or NaN}
\author{Catalina Anghel}
\examples{
ISOpure.util.matlab_less_than(5,3)
#[1] FALSE
ISOpure.util.matlab_less_than(3,5)
#[1] TRUE
ISOpure.util.matlab_less_than(5,NA)
#[1] FALSE
ISOpure.util.matlab_less_than(NA,5)
#[1] FALSE
ISOpure.util.matlab_less_than(5,Inf)
#[1] TRUE
ISOpure.util.matlab_less_than(Inf,5)
#[1] FALSE
}
\keyword{arith}
\keyword{NA}
\keyword{logic}