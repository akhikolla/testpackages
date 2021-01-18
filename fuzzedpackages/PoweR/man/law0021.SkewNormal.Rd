\name{law0021.SkewNormal}
\alias{law0021.SkewNormal}
\alias{law0021}

\title{The Skew Normal Distribution}

\description{Random generation for the Skew Normal distribution with parameters \code{xi}, \code{omega^2} and \code{alpha}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{xi}, \code{omega^2} and \code{alpha} are not specified they assume the default values of 0, 1 and 0, respectively.

The Skew Normal distribution with parameters \code{xi = }\eqn{\xi}, \code{omega^2 = }\eqn{\omega^2} and \code{alpha = }\eqn{\alpha} has density:
\deqn{ \left(\frac{2}{\omega}\right)\phi\left(\frac{x-\xi}{\omega}\right)\Phi\left(\alpha\left(\frac{x-\xi}{\omega}\right)\right) }{ 
       f(x) = 2/\omega \phi((x-\xi)/\omega) \Phi(\alpha((x-\xi)/\omega)) }
where \eqn{\phi(x)} is the standard normal probability density function and \eqn{\Phi(x)} is its cumulative distribution function.
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

\examples{res <- gensample(21,10000,law.pars=c(8,6,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Normal, Skew Normal}
