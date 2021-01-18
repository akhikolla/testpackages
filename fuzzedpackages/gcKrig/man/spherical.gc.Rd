
\name{spherical.gc}
\alias{spherical.gc}
\title{The Spherical Correlation Function of Class \code{\link{corr.gc}}}
\usage{
spherical.gc(range = NULL, nugget = TRUE)
}
\arguments{
\item{range}{a non-negative scalar of the range parameter in the spherical
correlation function.}

\item{nugget}{the nugget effect of the correlation function. If specified, it must be a scalar between 0 and 1.}
}
\value{
An object of class \code{\link{corr.gc}} representing the correlation component.
}
\description{
The spherical correlation function in spatial statistics.

By default, range parameter is not available, so
this function is used for likelihood inference and spatial prediction in function \code{mlegc}
and \code{predgc}.
Users need to specify if the correlation model includes a nugget effect
\code{nugget = TRUE} or not \code{nugget = FALSE}.

When both range and nugget parameters are given,
the function is used for simulation with function \code{simgc} in package \code{gcKrig}.

}
\details{
The spherical correlation function with a nugget \eqn{\tau^2} is of the form:

\deqn{\rho(h) = (1-\tau^2) (1 - 1.5(h/\phi) + 0.5(-h/\phi)^3) }
 when \eqn{h > 0} and \eqn{\rho(h) = 1} when \eqn{h = 0}, \eqn{h} is distance.
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{matern.gc}},  \code{\link{spherical.gc}}
}

\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Simulation}% __ONLY ONE__ keyword per line
