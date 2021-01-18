
\name{weibull.gc}
\alias{weibull.gc}
\title{The Weibull Marginal of Class \code{\link{marginal.gc}}}
\usage{
weibull.gc(shape = 1, scale = 1)
}
\arguments{
\item{shape}{a positive scalar of shape parameter in the Weibull distribution.}

\item{scale}{a positive scalar of scale parameter in the Weibull distribution.}
}
\value{
An object of class \code{\link{marginal.gc}} representing the marginal component.
}
\description{
The Weibull marginal used for simulation and computing correlation in the
trans-Gaussian random field in function \code{simgc} and \code{corrTG} of the package \code{gcKrig}.
It cannot be used in function \code{mlegc} nor \code{predgc} to make model inferences.
}

\details{
The Weibull distribution with \code{shape} parameter \eqn{a} and \code{scale} parameter \eqn{b}
has density given by

\deqn{(a/b) (x/b)^{a-1} exp(-(x/b)^a) }{(a/b) (x/b)^(a-1) exp(-(x/b)^a)}
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc}},  \code{\link{beta.gc}},
          \code{\link{binomial.gc}}, \code{\link{gm.gc}},
          \code{\link{gaussian.gc}}, \code{\link{negbin.gc}},
          \code{\link{poisson.gc}}, \code{\link{zip.gc}}
}

\keyword{Marginal}
\keyword{Simulation}
