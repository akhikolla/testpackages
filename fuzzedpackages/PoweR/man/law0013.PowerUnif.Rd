\name{law0013.PowerUnif}
\alias{law0013.PowerUnif}
\alias{law0013}

\title{The Power Uniform Distribution}

\description{Random generation for the Power Uniform distribution with parameter \code{power}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.}

\details{If \code{power} is not specified it assumes the default value of 1.

The Power Uniform distribution has density:
\deqn{ \frac{1}{1+j}x^{-\frac{j}{j+1}} }{ 
	   f(x) = 1/(1+j) x^-(j/(j+1)) }
where \code{power = j}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Quesenberry and Miller (1977), Power studies of some tests for uniformity, \emph{Journal of Statistical Computation and Simulation}, \bold{5}(3), 169--191 (see p. 178)
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See \code{\link{law0007.Uniform}} for the Uniform
  distribution. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(13,10000,law.pars=8)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Power Uniform}
