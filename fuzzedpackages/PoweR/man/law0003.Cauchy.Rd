\name{law0003.Cauchy}
\alias{law0003.Cauchy}
\alias{law0003}

\title{The Cauchy Distribution}

\description{Random generation for the Cauchy distribution with parameters \code{location} and \code{scale}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{location} or \code{scale} are not specified, they assume the default values of 0 and 1 respectively.

The Cauchy distribution has density:
\deqn{ \frac{1}{\pi s(1+(\frac{x-l}{s})^2)} }{ f(x) = 1 / (\pi s (1 + ((x-l)/s)^2)) }
where \eqn{l} is the location parameter and \eqn{s} is the scale parameter, for all \eqn{x}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See function \code{\link{rcauchy}} from package stats. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(3,10000,law.pars=c(9,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Cauchy}
