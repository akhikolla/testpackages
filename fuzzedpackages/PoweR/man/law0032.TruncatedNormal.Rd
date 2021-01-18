\name{law0032.TruncatedNormal}
\alias{law0032.TruncatedNormal}
\alias{law0032}

\title{The Truncated Normal Distribution}

\description{Random generation for the Truncated Normal distribution with parameters \code{a} and \code{b}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}
\details{If \code{a} and \code{b} are not specified they assume the default values of 0 and 1, respectively.

The Truncated Normal distribution with parameters \code{mu = }\eqn{\mu} and \code{sigma = }\eqn{\sigma} has density:
\deqn{ \frac{\exp(-x^2/2)}{\sqrt{2\pi}(\Phi(b)-\Phi(a))} }{ 
       f(x) = [1/\sigma \phi((x - \mu)/\sigma)] / [\Phi((b - \mu)/\sigma) - \Phi((a - \mu)/\sigma)] }
for \eqn{a \le x \le b}, where \eqn{\phi(x)} is the standard normal probability density function and \eqn{\Phi(x)} is its cumulative distribution function.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{law0002.Normal}} for the Normal
  distribution. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(32,10000,law.pars=c(2,3))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Normal, Truncated Normal}
