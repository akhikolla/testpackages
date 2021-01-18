
\name{binomial.gc}
\alias{binomial.gc}
\title{The Binomial Marginal of Class \code{\link{marginal.gc}}}
\usage{
binomial.gc(link = "logit", size = NULL, prob = NULL)
}
\arguments{
\item{link}{the model link function.}

\item{size}{number of trials (zero or more).}

\item{prob}{probability of success on each trial.}
}
\value{
An object of class \code{\link{marginal.gc}} representing the marginal component.
}
\description{
The binomial marginal parameterized in terms of its size and probability.

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
           \code{\link{gm.gc}}, \code{\link{gaussian.gc}},
          \code{\link{negbin.gc}}, \code{\link{poisson.gc}},
          \code{\link{weibull.gc}}, \code{\link{zip.gc}}
}

\keyword{Marginal}
\keyword{Simulation}
