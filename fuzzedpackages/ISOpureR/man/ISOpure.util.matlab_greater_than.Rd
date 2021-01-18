\name{ISOpure.util.matlab_greater_than}
\alias{ISOpure.util.matlab_greater_than}
\title{Greater than operator}
\description{Greater than function that matches Matlab behaviour when one of the arguments is NA (i.e. returns FALSE instead of NA)}
\usage{
ISOpure.util.matlab_greater_than(a, b)
}
\arguments{
  \item{a}{A numeric value (including Inf) or NA}
  \item{b}{A numeric value or NA}
}
\value{Logical: TRUE if a > b, FALSE if a <= b OR if one of a, b is NA or NaN}
\author{Catalina Anghel}
\examples{
ISOpure.util.matlab_greater_than(5,3)
#[1] TRUE
ISOpure.util.matlab_greater_than(3,5)
#[1] FALSE
ISOpure.util.matlab_greater_than(5,NA)
#[1] FALSE
ISOpure.util.matlab_greater_than(NA,5)
#[1] FALSE
ISOpure.util.matlab_greater_than(5,Inf)
#[1] FALSE
ISOpure.util.matlab_greater_than(Inf,5)
#[1] TRUE
}
\keyword{arith}
\keyword{NA}
\keyword{logic}