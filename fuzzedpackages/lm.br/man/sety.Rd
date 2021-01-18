\name{sety}
\alias{sety,Cpp_Clmbr-method}
\alias{sety}

\title{Set y-Values}

\description{
  Reset the response values in the C++ object.  
}

\usage{
\S4method{sety}{Cpp_Clmbr}( rWy )
}

\arguments{
\item{rWy}{vector of 'y' values, pre-multiplied by the square-root of 'weights'.}
}

\details{
The 'rWy' vector is simply the y-vector if the model does not specify weights.  The square-root of a vector 'W' is the vector 'rW' of the square-roots of the elements of 'W'.  The square-root of a matrix 'W' here is the matrix 'rW' such that rW*rW = W (a stricter definition than  rW*transpose(rW) = W).  
}

\note{
The pre-multiplied vector is more convenient as input during simulation tests.  'sety' changes the y-values only for the accessor functions 'sl', 'ci', 'cr' and 'mle'.  'rW' is the inverse square-root if 'inverse' was TRUE in the 'lm.br' call.
}


\examples{
#  A simulation test
x <- c( 1.0, 1.1, 1.3, 1.7, 2.4, 3.9, 5.7, 7.6, 8.4, 8.6 )
y <- x
LLmodel <- lm.br( y ~ x )    
countCLR <- countAF <- 0
theta <- 3
\donttest{for( i in 1:10000 )  {
  y <- 0 + (-1.)*pmin(x-theta,0) + (0.5)*pmax(x-theta,0) + rnorm(10)
  LLmodel$sety( y )
  stest <- LLmodel$sl( theta, 'clr', .0001, "V" )
  if( stest > 0.05 )  countCLR <- countCLR + 1
  stest <- LLmodel$sl( theta, 'af', .0001, "V" )
  if( stest > 0.05 )  countAF <- countAF + 1
  if( floor(i/1000) - i/1000 == 0 ) cat(i, countCLR/i, countAF/i, "\n")
}}
}

