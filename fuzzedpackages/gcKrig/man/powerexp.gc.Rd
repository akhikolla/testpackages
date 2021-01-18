
\name{powerexp.gc}
\alias{powerexp.gc}
\title{The Powered Exponential Correlation Function of Class \code{\link{corr.gc}}}
\usage{
powerexp.gc(range = NULL, kappa = 1, nugget = TRUE)
}
\arguments{
\item{range}{a non-negative scalar of the range parameter in powered
exponential correlation function.}

\item{kappa}{a scalar between 0 and 2; the value of the shape
parameter in the powered exponential correlation function.}

\item{nugget}{the nugget effect of the correlation function. If specified, it must be a scalar between 0 and 1.}
}
\value{
An object of class \code{\link{corr.gc}} representing the correlation component.
}
\description{
The powered exponential correlation function in spatial statistics.

By default, range parameter is not available, so
this function is used for likelihood inference and spatial prediction in function \code{mlegc}
and \code{predgc}.
Users need to specify the shape parameter \code{kappa} and if the correlation model includes a nugget effect
\code{nugget = TRUE} or not \code{nugget = FALSE}.

When both range and nugget parameters are given,
the function is used for simulation with function \code{simgc} in package \code{gcKrig}.


}
\details{
The powered exponential correlation function with a nugget \eqn{\tau^2} is of the form:

\deqn{\rho(h) = (1-\tau^2) exp((-h/\phi)^\kappa)}

  when \eqn{h > 0} and \eqn{\rho(h) = 1} when \eqn{h = 0}.
  Here \eqn{h} is distance, \eqn{\phi}{phi} is range parameter, \eqn{\kappa}{kappa} is the shape parameter and
  \eqn{\tau^2}{tau^2} is the nugget effect.

  When using the powered exponential correlation function, note that \eqn{0<\kappa \le 2}.
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{matern.gc}},  \code{\link{spherical.gc}}
}

\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Simulation}% __ONLY ONE__ keyword per line
