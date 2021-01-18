
\name{zip.gc}
\alias{zip.gc}
\title{The Zero-inflated Poisson Marginal of Class \code{\link{marginal.gc}}}
\usage{
zip.gc(link = "log", mu = NULL, od = NULL)
}
\arguments{
\item{link}{the model link function.}

\item{mu}{a non-negative scalar of the mean parameter.}

\item{od}{a non-negative scalar of the overdispersion parameter.}
}
\value{
An object of class \code{\link{marginal.gc}} representing the marginal component.
}

\description{
The zero-inflated Poisson marginal parameterized in terms of its mean and overdispersion.

By default, this function is used for likelihood inference and spatial prediction in function
\code{mlegc} and \code{predgc} of the package \code{gcKrig}.
When all marginal parameters are given, the function is used for simulation and computing correlation in a
trans-Gaussian random field in function \code{simgc} and \code{corrTG}.
}

\details{
The zero-inflated Poisson distribution with parameters \code{mu = a} and \code{od = b} has density

\deqn{b/(1+b) + exp(-(a+ab))/(1+b)}
when \eqn{y = 0}, and

\deqn{exp(-(a+ab))*(a+ab)^y/((1+b)y!)}
when \eqn{y = 1, 2, \ldots}


Under this parameterization, \eqn{var(Y)= mu + od*mu^2}, where
\eqn{mu} is the mean parameter and \eqn{od} is the overdispersion parameter.
For more details see Han and De Oliveira (2016).
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}

\references{
Han, Z. and De Oliveira, V. (2016) On the correlation structure of Gaussian
            copula models for geostatistical count data.
            \emph{Australian and New Zealand Journal of Statistics, 58:47-69}.
}

\seealso{
\code{\link{marginal.gc}},  \code{\link{beta.gc}}, \code{\link{binomial.gc}},
           \code{\link{gm.gc}}, \code{\link{gaussian.gc}},
           \code{\link{negbin.gc}}, \code{\link{poisson.gc}},
          \code{\link{weibull.gc}}
}

\keyword{Marginal}
\keyword{Simulation}
