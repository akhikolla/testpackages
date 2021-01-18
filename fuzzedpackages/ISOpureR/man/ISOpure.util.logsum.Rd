\name{ISOpure.util.logsum}
\alias{ISOpure.util.logsum}
\title{Log-sum-exp}
\description{Prevents underflow/overflow using the log-sum-exp trick }
\usage{
ISOpure.util.logsum(xx, dimen);
}
\arguments{
  \item{xx}{A matrix of numerical values}
  \item{dimen}{The dimension along which the long sum is taken (1 for row, 2 for column)}
}
\value{Returns log(sum(exp(x),dimen)), the log sum of exps, summing over dimension dimen
but in a way that tries to avoid underflow/overflow.}
\author{Gerald Quon and Catalina Anghel}
\examples{
x <-  c(1, 1e20, 1e40, -1e40, -1e20, -1);
x <- as.matrix(x);

# compute log sum exp without the function
log(sum(exp(x)))
#[1] Inf

# compute log sum exp with the function
ISOpure.util.logsum(x, 1)
#[1] 1e+40 
}
\keyword{math}