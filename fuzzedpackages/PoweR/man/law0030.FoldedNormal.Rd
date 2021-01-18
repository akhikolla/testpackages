\name{law0030.FoldedNormal}
\alias{law0030.FoldedNormal}
\alias{law0030}

\title{The Folded Normal Distribution}

\description{Random generation for the Folded Normal distribution with parameters \code{mu} and \code{sigma}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu} and \code{sigma} are not specified they assume the default values of 0 and 1, respectively.

The Folded Normal distribution with parameters \code{mu = }\eqn{\mu} and \code{sigma = }\eqn{\sigma} has density:
\deqn{ dnorm(x,mu,sigma2)+dnorm(x,-mu,sigma2) }{ 
       f(x) = 1/(\sqrt(2\pi) \sigma) e^-((-x - \mu)^2/(2 \sigma^2)) + 1/(\sqrt(2\pi) \sigma) e^-((x - \mu)^2/(2 \sigma^2)) }
for \eqn{x \ge 0}.
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

\examples{res <- gensample(30,10000,law.pars=c(8,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Normal, Folded Normal}
