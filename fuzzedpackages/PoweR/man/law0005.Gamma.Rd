\name{law0005.Gamma}
\alias{law0005.Gamma}
\alias{law0005}

\title{The Gamma Distribution}

\description{Random generation for the Gamma distribution with parameters \code{shape} and \code{rate}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{shape} or \code{rate} are not specified they assume the default values of 2 and 1, respectively.

The Gamma distribution has density:
\deqn{ \frac{1}{b^a\Gamma(a)}x^{a-1}exp^{-x/b} }{ f(x) = 1/(b^a \Gamma(a)) x^(a-1) e^-(x/b) }
for \eqn{x \ge 0}, a > 0 and b > 0; where \eqn{a} is the \code{shape} parameter and \eqn{b} is the \code{rate} parameter.

Here \eqn{\Gamma(a)} is the \code{gamma} function implemented by R and defined in its help.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See function \code{\link{rgamma}} from package
  stats. See. \code{\link{Distributions}} for other standard
  distributions. Type \command{help(gamma)} for additional information about the gamma function. 
}

\examples{res <- gensample(5,10000,law.pars=c(9,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Gamma}
