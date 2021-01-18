\name{law0033.Nout}
\alias{law0033.Nout}
\alias{law0033}

\title{The Normal with outliers Distribution}

\description{Random generation for the Normal with outliers distribution with parameter \code{a} which belongs to \code{{1,2,3,4,5}}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.
}

\details{If \code{a} is not specified it assumes the default value of 1.

Five cases of standard normal distributions with outliers, hereon termed \code{Nout1} to \code{Nout5}, 
consisting of observations drawn from a standard normal distribution where some of the values 
are randomly replaced by extreme observations.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Romao, X., Delgado, R. and Costa, A. (2010), 
An empirical power comparison of univariate goodness-of-fit tests for normality,
\emph{Journal of Statistical Computation and Simulation}, \bold{80}(5), 545--591.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{law0002.Normal}} for the Normal
  distribution. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(33,10000,law.pars=4)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Normal, Normal with outliers}
