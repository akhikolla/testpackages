\name{law0038.AsymmetricPowerDistribution}
\alias{law0038.AsymmetricPowerDistribution}
\alias{law0038}

\title{The Asymmetric Power Distribution}

\description{Random generation for the Asymmetric Power Distribution with parameters \code{theta}, \code{phi}, \code{alpha} and \code{lambda}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{theta}, \code{phi}, \code{alpha} and \code{lambda} are not specified they assume the default values of 0, 1, 0.5 and 2, respectively.

The Asymmetric Power Distribution with parameters \code{theta},
\code{phi}, \code{alpha} and \code{lambda} has density:
\deqn{f(u) = \frac{1}{\phi}\frac{\delta^{1/\lambda}_{\alpha,\lambda}}{\Gamma(1+1/\lambda)}\exp\left[-\frac{\delta_{\alpha,\lambda}}{\alpha^{\lambda}}\left|\frac{u-\theta}{\phi}\right|^{\lambda}\right]}{f(u) = (1/\phi)*(\delta^{1/\lambda}_{\alpha,\lambda})/(\Gamma(1+1/\lambda))*\exp[-(\delta_{\alpha,\lambda})/(\alpha^{\lambda})*|(u-\theta)/\phi|^{\lambda}]}
if \deqn{u\leq0}{u <= 0} and
\deqn{f(u) =
  \frac{1}{\phi}\frac{\delta^{1/\lambda}_{\alpha,\lambda}}{\Gamma(1+1/\lambda)}\exp\left[-\frac{\delta_{\alpha,\lambda}}{(1-\alpha)^{\lambda}}\left|\frac{u-\theta}{\phi}\right|^{\lambda}\right]}{f(u)
  =
  (1/\phi)*(\delta^{1/\lambda}_{\alpha,\lambda})/(\Gamma(1+1/\lambda))*\exp[-(\delta_{\alpha,\lambda})/(1-\alpha)^{\lambda}*|(u-\theta)/\phi|^{\lambda}]}
if \deqn{u\leq0,}{u > 0,}where \eqn{0<\alpha<1, \lambda>0} and \eqn{\delta_{\alpha,\lambda}=\frac{2\alpha^{\lambda}(1-\alpha)^{\lambda}}{\alpha^{\lambda}+(1-\alpha)^{\lambda}}}{\delta_{\alpha,\lambda}=(2*\alpha^{\lambda}(1-\alpha)^{\lambda})/(\alpha^{\lambda}+(1-\alpha)^{\lambda})}.




The mean and variance of APD are defined respectively by
\deqn{E(U) = \theta+\phi\frac{\Gamma(2/\lambda)}{\Gamma(1/\lambda)} [1-2\alpha]\delta_{\alpha,\lambda}^{-1/\lambda}}{E(U) = \theta+\phi*\Gamma(2/\lambda)/\Gamma(1/\lambda) * [1-2*\alpha]*\delta_{\alpha,\lambda}^{-1/\lambda}} and  
\deqn{V(U) = \phi^2 
  \frac{\Gamma(3/\lambda)\Gamma(1/\lambda)[1-3\alpha+3\alpha^2]-\Gamma(2/\lambda)^2[1-2\alpha]^2}{\Gamma^2(1/\lambda)}
   \delta_{\alpha,\lambda}^{-2/\lambda}.
}{V(U) =\phi^2 *
  (\Gamma(3/\lambda)*\Gamma(1/\lambda)*[1-3\alpha+3\alpha^2]-\Gamma(2/\lambda)^2*[1-2*\alpha]^2)/(\Gamma(1/\lambda))^2
  * \delta_{\alpha,\lambda}^{-2/\lambda}.
}
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Komunjer, I. (2007), Asymmetric Power Distribution: Theory and Applications to Risk Measurement, 
\emph{Journal of Applied Econometrics}, \bold{22}, 891--921.
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See \code{\link{Distributions}} for other standard distributions.}

\examples{
res <- gensample(38,10000,law.pars=c(3,2,0.5,1))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Asymmetric Power Distribution}
