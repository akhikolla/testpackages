\name{law0001.Laplace}
\alias{law0001.Laplace}
\alias{law0001}

\title{The Laplace Distribution}

\description{Random generation for the Laplace distribution with parameters \code{mu} and \code{b}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu} or \code{b} are not specified they assume the default values of 0 and 1, respectively.

The Laplace distribution has density:
\deqn{ \frac{1}{2b}\exp \left( -\frac{|x-\mu|}{b} \right) }{ f(x) = (1/2b) e^-(|x - \mu|/b) }
where \eqn{\mu} is a location parameter and \code{b} > 0, which is sometimes referred to as the diversity, is a scale parameter.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See function \code{urlaplace()} from \code{Runuran}
  package. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(1,10000,law.pars=c(9,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Laplace}
