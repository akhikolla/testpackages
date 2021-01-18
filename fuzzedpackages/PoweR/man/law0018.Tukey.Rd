\name{law0018.Tukey}
\alias{law0018.Tukey}
\alias{law0018}

\title{The Tukey Distribution}

\description{Random generation for the Tukey distribution with parameter \code{lambda}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.}

\details{If \code{lambda} is not specified it assumes the default value of 1.

The Tukey distribution with \code{lambda = }\eqn{\lambda} has \eqn{E[X] = 0} and 
  \eqn{Var[X] = 2/(\lambda^2) (1/(2\lambda+1) - \Gamma^2(\lambda+1)/\Gamma(2\lambda+2))}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{\code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(18,10000,law.pars=8)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Tukey}
