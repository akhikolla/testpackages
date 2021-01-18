\name{law0002.Normal}
\alias{law0002.Normal}
\alias{law0002}

\title{The Normal Distribution}

\description{Random generation for the Normal distribution with parameters \code{mu} and \code{sigma}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu} or \code{sigma} are not specified they assume the default values of 0 and 1, respectively.

The Normal distribution has density:
\deqn{ (\sqrt{2\pi}\sigma)^{-1}\exp^{-\frac{x^2}{2\sigma^2}} }{ f(x) = 1/(\sqrt(2\pi*\sigma)) e^-((x - \mu)^2/(2*\sigma^2)) }
where \eqn{\mu} is the mean of the distribution and \eqn{\sigma} is the standard deviation.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See function \code{link{rnorm}} from stats package. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(2,10000,law.pars=c(9,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Normal, Gaussian}
