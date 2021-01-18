
\name{negbin.gc}
\alias{negbin.gc}
\title{The Negative Binomial Marginal of Class \code{\link{marginal.gc}}}
\usage{
negbin.gc(link = "log", mu = NULL, od = NULL)
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
The negative binomial marginal parameterized in terms of its mean and overdispersion.

By default, this function is used for likelihood inference and spatial prediction in functions
\code{mlegc} and \code{predgc} of the package \code{gcKrig}.
When all marginal parameters are given, the function is used for simulation and computing correlation in a
trans-Gaussian random field in functions \code{simgc} and \code{corrTG}.

}

\details{
The negative binomial distribution with parameters \code{mu = a} and \code{od = 1/b} has density

\deqn{\frac{\Gamma(y+b)}{\Gamma(b)y!} \Big(\frac{b}{b+a}\Big)^b
        \Big(1 - \frac{b}{b+a}\Big)^y }{\Gamma(y+b)/(\Gamma(b)y!) * (b/(b+a))^b * (a/(b+a))^y
}
which is called NB2 by Cameron and Trivedi (2013).
Under this parameterization, \eqn{var(Y)= mu + od*mu^2}, where
\eqn{mu} is the mean parameter and \eqn{od} is the overdispersion parameter.
For more details see Han and De Oliveira (2016).
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}

\references{
 Cameron,A.C. and Trivedi,P.K. (2013) Regression Analysis of Count Data.
 \emph{Cambridge University Press, 2nd Edition}.

 Han, Z. and De Oliveira, V. (2016) On the correlation structure of Gaussian
            copula models for geostatistical count data.
            \emph{Australian and New Zealand Journal of Statistics, 58:47-69}.
}


\seealso{
\code{\link{marginal.gc}},  \code{\link{beta.gc}},
          \code{\link{binomial.gc}}, \code{\link{gm.gc}},
          \code{\link{gaussian.gc}}, \code{\link{poisson.gc}},
          \code{\link{weibull.gc}}, \code{\link{zip.gc}}
}

\keyword{Marginal}
\keyword{Simulation}
