\name{law0028.GeneralizedExtValue}
\alias{law0028.GeneralizedExtValue}
\alias{law0028}

\title{The Generalized Extreme Value Distribution}

\description{Random generation for the Generalized Extreme Value distribution with parameters \code{mu}, \code{sigma} and \code{xi}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu}, \code{sigma} and \code{xi} are not specified they assume the default values of 0, 1 and 1, respectively.

The Generalized Extreme Value distribution with parameters \code{mu = }\eqn{\mu}, \code{sigma = }\eqn{\sigma} and \code{xi = }\eqn{\xi} has density:
\deqn{ [1+z]_{+}^{-\frac{1}{\xi}-1}\exp\left\{-[1+z]_{+}^{-\frac{1}{\xi}}\right\}/\sigma }{ 
       f(x) = (1 + z)_{+}^{-1/\xi - 1} {e^[(1 + z)_{+}^{-1/\xi}]}/\sigma }
for \eqn{\xi > 0} or \eqn{\xi < 0}, where \eqn{z = \xi (x - \mu)/\sigma}. If \eqn{\xi = 0}, PDF is as same as in the Gumbel distribution.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{law0026.Gumbel}} for the Gumbel
  distribution. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(28,10000,law.pars=c(8,6,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Generalized Extreme Value}
