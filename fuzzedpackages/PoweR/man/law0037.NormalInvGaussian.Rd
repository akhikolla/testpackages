\name{law0037.NormalInvGaussian}
\alias{law0037.NormalInvGaussian}
\alias{law0037}

\title{The Normal-inverse Gaussian Distribution}

\description{Random generation for the Normal-inverse Gaussian distribution with parameters \code{shape}, \code{skewness}, \code{location} and \code{scale}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}
\details{If \code{shape}, \code{skewness}, \code{location} and \code{scale} are not specified they assume the default values of 1, 0, 0 and 1, respectively.

The Normal-inverse Gaussian distribution with parameters \code{shape = }\eqn{\alpha}, 
\code{skewness = }\eqn{\beta}, \code{location = }\eqn{\mu} and \code{scale = }\eqn{\delta} has density:
\deqn{ \frac{\alpha\delta K_1(\alpha\sqrt{\delta^2+(x-\mu)^2})}{\pi\sqrt{\delta^2+(x-\mu)^2}}e^{\delta\gamma+\beta(x-\mu)} }{ 
       f(x) = \alpha \delta K_1 (\alpha \sqrt(\delta^2 + (x - \mu)^2)) / (\pi \sqrt(\delta^2 + (x - \mu)^2)) * e^(\delta \gamma + \beta (x - \mu)) }
where \eqn{\gamma = \sqrt(\alpha^2 - \beta^2)} and 
      \eqn{K_1} denotes a modified Bessel function of the second kind.

The mean and variance of NIG are defined respectively \eqn{\mu + \beta \delta / \gamma} and  
                                                      \eqn{\delta \alpha^2 / \gamma^3}.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See package \code{fBasics}. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(37,10000,law.pars=c(3,2,1,0.5))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Normal-inverse Gaussian}
