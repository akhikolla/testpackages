\name{beta.gc}
\alias{beta.gc}
\title{The Beta Marginal of Class \code{\link{marginal.gc}}}
\usage{
beta.gc(shape1 = 1, shape2 = 1)
}
\arguments{
\item{shape1}{non-negative scalar, the shape parameter of the Beta distribution.}

\item{shape2}{non-negative scalar, another shape parameter of the Beta distribution.}
}
\value{
An object of class \code{\link{marginal.gc}} representing the marginal component.
}
\description{
The Beta marginal used for simulation and computing correlation in the
trans-Gaussian random field in function \code{simgc} and \code{corrTG} of the package \code{gcKrig}.
It cannot be used in function \code{mlegc} nor \code{predgc} to make model inferences.
}

\details{
The Beta distribution with parameters \code{shape1 = a} and \code{shape2 = b} has density


\deqn{\frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)} x^{a-1} (1-x)^{b-1}}


for \eqn{a > 0, b > 0} and \eqn{0 \le x \le 1} where the boundary values at
\code{x = 0} or \code{x = 1} are defined as by continuity (as limits).
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc}}, \code{\link{binomial.gc}}, \code{\link{gm.gc}},
          \code{\link{gaussian.gc}}, \code{\link{negbin.gc}},
          \code{\link{poisson.gc}}, \code{\link{weibull.gc}}, \code{\link{zip.gc}}
}

\keyword{Marginal}
\keyword{Simulation}

