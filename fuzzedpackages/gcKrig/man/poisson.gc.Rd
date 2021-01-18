
\name{poisson.gc}
\alias{poisson.gc}
\title{The Poisson Marginal of Class \code{\link{marginal.gc}}}
\usage{
poisson.gc(link = "log", lambda = NULL)
}
\arguments{
\item{link}{the model link function.}

\item{lambda}{a non-negative scalar of the mean parameter.}
}
\value{
An object of class \code{\link{marginal.gc}} representing the marginal component.
}
\description{
The Poisson marginal parameterized in terms of its mean.

By default, this function is used for likelihood inference and spatial prediction in function
\code{mlegc} and \code{predgc} of the package \code{gcKrig}.
When all marginal parameters are given, the function is used for simulation and computing correlation in a
trans-Gaussian random field in function \code{simgc} and \code{corrTG}.

}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc}},  \code{\link{beta.gc}},
          \code{\link{binomial.gc}}, \code{\link{gm.gc}},
          \code{\link{gaussian.gc}}, \code{\link{negbin.gc}},
          \code{\link{weibull.gc}}, \code{\link{zip.gc}}
}

\keyword{Marginal}
\keyword{Simulation}
