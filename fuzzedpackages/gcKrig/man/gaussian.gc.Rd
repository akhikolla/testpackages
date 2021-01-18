
\name{gaussian.gc}
\alias{gaussian.gc}
\title{The Gaussian Marginal of Class \code{\link{marginal.gc}}}
\usage{
gaussian.gc(mean = 0, sd = 1)
}
\arguments{
\item{mean}{the mean of the Gaussian distribution, a scalar.}

\item{sd}{a positive scalar, the standard deviation of the Gaussian distribution.}
}
\value{
An object of class \code{\link{marginal.gc}} representing the marginal component.
}
\description{
The Gaussian marginal used for simulation and computing correlation in the
trans-Gaussian random field in function \code{simgc} and \code{corrTG} of the package \code{gcKrig}.
It cannot be used in function \code{mlegc} nor \code{predgc} to make model inferences.
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc}}, \code{\link{beta.gc}}, \code{\link{binomial.gc}},
          \code{\link{gm.gc}}, \code{\link{negbin.gc}},
          \code{\link{poisson.gc}}, \code{\link{weibull.gc}}, \code{\link{zip.gc}}
}

\keyword{Marginal}
\keyword{Simulation}
