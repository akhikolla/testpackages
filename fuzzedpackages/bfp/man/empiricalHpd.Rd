\name{empiricalHpd}
\alias{empiricalHpd}

\title{Construct an empirical HPD interval from samples}
\description{
Construct an empirical highest posterior density (HPD) interval from
samples which have been drawn from the distribution of a quantity of
interest. 
}
\usage{
empiricalHpd(theta, level)
}
\arguments{
  \item{theta}{the vector of samples}
  \item{level}{the credible level}
}
\value{
  A vector with the estimated lower and upper bounds of the HPD
  interval. 
}
\author{Daniel Saban\'es Bov\'e}
\examples{
## draw standard normal variates
test <- rnorm(n=1000)

## estimate the 95% HPD interval with these samples:
empiricalHpd(theta=test, level=0.95)

## compare with true HPD:
qnorm(p=c(0.025, 0.975))
}
\keyword{htest}

