\name{law0026.Gumbel}
\alias{law0026.Gumbel}
\alias{law0026}

\title{The Gumbel Distribution}

\description{Random generation for the Gumbel distribution with parameters \code{mu} and \code{sigma}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu} or \code{sigma} are not specified, they assume the default values of 1.

The Gumbel distribution with parameters \code{mu = }\eqn{\mu} and \code{sigma = }\eqn{\sigma} has density:
\deqn{ \frac{1}{\sigma}\exp\left\{-\exp\left[-\left(\frac{x-\mu}{\sigma}\right)\right]-\left(\frac{x-\mu}{\sigma}\right)\right\} }{ 
       f(x) = 1/\sigma e^[ -e^(-(x - \mu)/\sigma) - (x - \mu)/\sigma ] }
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{law0028.GeneralizedExtValue}} for the
  Generalized Extreme Value distribution. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(26,10000,law.pars=c(9,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Gumbel}
