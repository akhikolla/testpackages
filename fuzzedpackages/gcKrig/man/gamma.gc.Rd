
\name{gm.gc}
\alias{gm.gc}
\title{The Gamma Marginal of Class \code{\link{marginal.gc}}}
\usage{
gm.gc(shape = 1, rate = 1)
}
\arguments{
\item{shape}{a non-negative scalar, shape parameter of the Gamma distribution.}

\item{rate}{a non-negative scalar, rate parameter of the Gamma distribution.}
}
\value{
An object of class \code{\link{marginal.gc}} representing the marginal component.
}
\description{
The Gamma marginal used for simulation and computing correlation in the
trans-Gaussian random field in functions \code{simgc} and \code{corrTG} of the package \code{gcKrig}.
It cannot be used in functions \code{mlegc} nor \code{predgc} to make model inferences.
}
\details{
The Gamma distribution with parameters \code{shape = a} and \code{rate = r} has density

\deqn{\frac{r^{a}}{\Gamma(a)} x^{a-1} exp(-rx)}{r^a/\Gamma(a) x^(a-1) exp(-rx)}
for \eqn{x \ge 0}, \eqn{a > 0} and \eqn{s > 0}.
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc}}, \code{\link{beta.gc}}, \code{\link{binomial.gc}},
          \code{\link{gaussian.gc}}, \code{\link{negbin.gc}},
          \code{\link{poisson.gc}}, \code{\link{weibull.gc}}, \code{\link{zip.gc}}
}

\keyword{Marginal}
\keyword{Simulation}
