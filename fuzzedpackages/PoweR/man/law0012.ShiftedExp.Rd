\name{law0012.ShiftedExp}
\alias{law0012.ShiftedExp}
\alias{law0012}

\title{The Shifted Exponential Distribution}

\description{Random generation for the Shifted Exponential distribution with parameters \code{l} and \code{rate}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{l} or \code{rate} are not specified they assume the default values of 0 and 1, respectively.

The Shifted Exponential distribution has density
\deqn{ b\exp\{-(x-l)b\} }{ 
	   f(x) = b e^(-(x-1)b) }
for \eqn{x \le 1}, where \eqn{rate = b}. 
The mean is \eqn{E(X) = l + 1/b}, and the \eqn{Var(X) = 1/(b^2)}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{\code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(12,10000,law.pars=c(8,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Shifted Exponential}
