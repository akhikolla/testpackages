\name{law0029.GeneralizedArcsine}
\alias{law0029.GeneralizedArcsine}
\alias{law0029}

\title{The Generalized Arcsine Distribution}

\description{Random generation for the Generalized Arcsine distribution with parameters \code{alpha}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.}

\details{If \code{alpha} is not specified it assumes the default value of 0.5.

The Generalized Arcsine distribution with parameter \code{alpha = }\eqn{\alpha} has density:
\deqn{ \frac{\sin(\pi\alpha)}{\pi}x^{-\alpha}(1-x)^{\alpha-1} }{ 
       f(x) = sin(\pi*\alpha)/\pi x^(-\alpha) (1 - x)^(\alpha - 1) }
for \eqn{0 < \alpha < 1}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(29,10000,law.pars=0.8)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Generalized Arcsine}
