\name{law0023.GeneralizedPareto}
\alias{law0023.GeneralizedPareto}
\alias{law0023}

\title{The Generalized Pareto Distribution}

\description{Random generation for the Generalized Pareto distribution with parameters \code{mu}, \code{sigma} and \code{xi}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu}, \code{sigma} and \code{xi} are not specified they assume the default values of 0, 1 and 0, respectively.

The Generalized Pareto distribution with parameters \code{mu = }\eqn{\mu}, \code{sigma = }\eqn{\sigma} and \code{xi = }\eqn{\xi} has density:
\deqn{ \frac{1}{\sigma}\left(1+\frac{\xi(x-\mu)}{\sigma} \right)^{(-\frac{1}{\xi}-1)} }{ 
       f(x) = 1/\sigma (1 + \xi (x - \mu)/\sigma)^(-1/\xi - 1) }
where \eqn{x \ge \mu} if \eqn{\xi \ge 0} and \eqn{x \le \mu - \sigma/\xi} if \eqn{\xi < 0}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(23,10000,law.pars=c(8,6,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Generalized Pareto}
