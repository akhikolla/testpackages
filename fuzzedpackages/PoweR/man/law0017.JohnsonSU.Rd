\name{law0017.JohnsonSU}
\alias{law0017.JohnsonSU}
\alias{law0017}

\title{The Johnson SU Distribution}

\description{Random generation for the Johnson SU distribution with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu}, \code{sigma}, \code{nu} and \code{tau} are not specified they assume the default values of 0, 1, 0 and 0.5, respectively.

The Johnson SU distribution with parameters \code{mu = }\eqn{\mu}, \code{sigma = }\eqn{\sigma}, \code{nu = }\eqn{\nu} and \code{tau = }\eqn{\tau} has density:
\deqn{ \frac{1}{c\sigma\tau}\frac{1}{\sqrt{z^2+1}}\frac{1}{\sqrt{2\pi}}e^{-r^2/2} }{ 
       f(x) = 1/(c*\sigma*\tau) 1/(\sqrt(z^2+1)) 1/(\sqrt(2*\pi)) e^-(r^2/2) }
where \eqn{r = -\nu + (1/\tau)sinh^-1(z)}, 
      \eqn{z = (x - (\mu + c*\sigma (\sqrt(\omega)) sinh(w)))/(c*\sigma)}, 
	  \eqn{c = ((w-1)(w cosh(2\omega)+1)/2)^-1/2},
	  \eqn{w = e^(\tau^2)} and \eqn{\omega = -\nu\tau}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(17,10000,law.pars=c(9,8,6,0.5))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Johnson SU}
