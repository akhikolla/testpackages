\name{law0035.Exponential}
\alias{law0035.Exponential}
\alias{law0035}

\title{The Exponential Distribution}

\description{Random generation for the Exponential distribution with rate \code{rate} (i.e., mean \code{1/rate}).

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.}

\details{If \code{rate} is not specified it assumes the default value of 1.

The Exponential distribution with \code{rate =} \eqn{\lambda} has density:
\deqn{ \lambda exp^{-\lambda x} }{ 
       f(x) = \lambda e^{-\lambda x} }
for \eqn{x \ge 0}. 
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(35,10000,law.pars=8)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Exponential}
