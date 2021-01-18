\name{law0020.JohnsonSB}
\alias{law0020.JohnsonSB}
\alias{law0020}

\title{The Johnson SB Distribution}

\description{Random generation for the Johnson SB distribution with parameters \code{g} and \code{d}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{g} and \code{d} are not specified they assume the default values of 0 and 1, respectively.

The Johnson SB distribution has density:
\deqn{ \frac{d}{\sqrt{2\pi}}\frac{1}{x(1-x)}e^{-\frac{1}{2}\left(g+d\ln\frac{x}{1-x} \right)^2} }{ 
       f(x) = d/(\sqrt(2\pi)) 1/(x(1-x)) e^-(1/2(g + d ln(x/(1-x)))) }
where \eqn{d > 0}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(20,10000,law.pars=c(8,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Johnson SB}
