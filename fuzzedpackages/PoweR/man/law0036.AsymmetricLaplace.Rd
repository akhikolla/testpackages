\name{law0036.AsymmetricLaplace}
\alias{law0036.AsymmetricLaplace}
\alias{law0036}

\title{The Asymmetric Laplace Distribution}

\description{Random generation for the Asymmetric Laplace distribution with parameters \code{mu}, \code{b} and \code{k}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}
\details{If \code{mu}, \code{b} or \code{k} are not specified they assume the default values of 0, 1 and 2, respectively.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See package \code{VGAM}. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(36,10000,law.pars=c(9,2,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Asymmetric Laplace}
