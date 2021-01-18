\name{law0031.MixtureNormal}
\alias{law0031.MixtureNormal}
\alias{law0031}

\title{The Mixture Normal Distribution}

\description{Random generation for the Mixture Normal distribution with parameters \code{p}, \code{m} and \code{d}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}
\details{If \code{p}, \code{m} and \code{d} are not specified they assume the default values of 0.5, 0 and 1, respectively.

The Mixture Normal distribution has density:
\deqn{ \frac{1}{\sqrt{2\pi}}\left[pe^{-\frac{(x-m)^2}{2d^2}}+(1-p)e^{-\frac{x^2}{2}}\right] }{ 
       f(x) = 1/(\sqrt(2\pi)) (p e^-((x - m)^2)/(2 d^2) + (1 - p) e^-(x^2/2)) }
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

\examples{res <- gensample(31,10000,law.pars=c(0.9,8,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Normal, Mixture Normal}
