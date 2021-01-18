\name{law0025.Stable}
\alias{law0025.Stable}
\alias{law0025}

\title{The Stable Distribution}

\description{Random generation for the Stable distribution with parameters \code{stability}, \code{skewness}, \code{scale} and \code{location}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{stability}, \code{skewness}, \code{scale} and \code{location} 
are not specified they assume the default values of 2, 0, 1 and 0, respectively.

The Stable distribution with parameters \code{stability = }\eqn{\alpha}, 
\code{skewness = }\eqn{\beta}, \code{scale = }\eqn{c} and \code{location = }\eqn{\mu}
doesn't have an analytically expressible probability density function, except for some parameter values.
The parameters have conditions : \eqn{0 < \alpha \le 2}, \eqn{-1 \le \beta \le 1} and \eqn{c > 0}.

The mean of Stable distribution is defined \eqn{\mu} when \eqn{\alpha > 1}, otherwise undefined.

The variance of Stable distribution is defined \eqn{2 c^2} when \eqn{\alpha = 2}, otherwise infinite.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(25,10000,law.pars=c(2,1,1,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Stable}
