\name{law0004.Logistic}
\alias{law0004.Logistic}
\alias{law0004}

\title{The Logistic Distribution}

\description{Random generation for the Logistic distribution with parameters \code{location} and \code{scale}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{location} or \code{scale} are omitted, they assume the default values of 0 and 1 respectively.

The Logistic distribution with \code{location =} \eqn{\mu} and \code{scale = s} has distribution function

\deqn{ \frac{1}{1 + exp^{-\frac{(x-\mu)}{s}}} }{ F(x) = 1 / (1 + e^(-(x - \mu)/s)) }

and density

\deqn{ \frac{exp^{-\frac{(x-\mu)}{s}}}{s(1+exp^{-\frac{(x-\mu)}{s}})^2} }{ f(x) = 1/s e^(-(x - \mu)/s) (1 + e^(-(x - \mu)/s))^-2 }

It is a long-tailed distribution with mean \eqn{\mu} and variance \eqn{(\pi^2)/3 s^2}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See function \code{\link{rlogis}} from package stats. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(4,10000,law.pars=c(9,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Logistic}
