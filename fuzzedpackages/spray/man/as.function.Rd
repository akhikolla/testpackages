\name{as.function.spray}
\alias{as.function.spray}
\title{Coerce a spray object to a function}
\description{Coerce a spray object to a function}
\usage{
\method{as.function}{spray}(x,...)
}
\arguments{
  \item{x}{spray object, interpreted as a multivariate polynomial}
  \item{...}{Further arguments, currently ignored}
}
\value{
Returns a function; this function returns a numeric vector.
}
\author{
Robin K. S. Hankin
}
\note{
Coercion is possible even if some indices are zero or negative.  The
function is not vectorized in the arity of its argument.
}
\examples{

S1 <- spray(matrix(sample(-2:2,replace=TRUE,21),ncol=3),rnorm(7),addrepeats=TRUE)
S2 <- spray(matrix(sample(-2:2,replace=TRUE,15),ncol=3),rnorm(5),addrepeats=TRUE)

f1 <- as.function(S1)
f2 <- as.function(S2)

f3 <- as.function(S1*S2)


x <- 4:6

f1(x)*f2(x)-f3(x)  # should be zero



# coercion is vectorized:

f1(matrix(1:33,ncol=3))


  
}

\keyword{symbolmath}
