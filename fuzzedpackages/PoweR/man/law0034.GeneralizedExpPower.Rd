\name{law0034.GeneralizedExpPower}
\alias{law0034.GeneralizedExpPower}
\alias{law0034}

\title{The Generalized Exponential Power Distribution}

\description{Random generation for the Generalized Exponential Power distribution with parameters \code{t1}, \code{t2} and \code{t3}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.
}
\details{If \code{t1}, \code{t2} and \code{t3} are not specified they assume the default value of 0.5, 0 and 1, respectively.

The Generalized Exponential Power distribution has density:
\deqn{ p(x;\gamma,\delta,\alpha,\beta,z_0) \propto e^-{\delta|x|^\gamma} |x|^{-\alpha}(log|x|)^{-\beta} }{ 
       p(x;\gamma,\delta,\alpha,\beta,z_0) \propto e^-(\delta |x|^\gamma) |x|^(-\alpha) log|x|^(-\beta) }
for \eqn{x \ge z_0}, and the density equals to \eqn{p(x;\gamma,\delta,\alpha,\beta,z_0)} for \eqn{x < z_0}.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Desgagne, A., Lafaye de Micheaux, P. and Leblanc, A. (2013), 
Test of Normality Against Generalized Exponential Power Alternatives,
\emph{Communications in Statistics - Theory and Methods}, \bold{42}(1), 164--190.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(34,10000,law.pars=c(1,8,4))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Generalized Exponential Power}
