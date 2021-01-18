\name{law0024.GeneralizedError}
\alias{law0024.GeneralizedError}
\alias{law0024}

\title{The Generalized Error Distribution}

\description{Random generation for the Generalized Error distribution with parameters \code{mu}, \code{sigma} and \code{p}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu}, \code{sigma} and \code{p} are not specified they assume the default values of 0, 1 and 1, respectively.

The Generalized Error distribution with parameters \code{mu = }\eqn{\mu}, \code{sigma = }\eqn{\sigma} and \code{p = }\eqn{p} has density:
\deqn{ \frac{1}{2p^{1/p}\Gamma(1+1/p)\sigma}\exp\left[-\frac{1}{p\sigma^p}|x-\mu|^p\right] }{ 
       f(x) = 1/(2p^(1/p) \Gamma(1 + 1/p) \sigma) e^-(1/(p \sigma^p) |x - \mu|^p) }
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(24,10000,law.pars=c(8,6,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Generalized Error}
