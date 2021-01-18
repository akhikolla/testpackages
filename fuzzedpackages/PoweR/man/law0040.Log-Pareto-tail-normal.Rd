\name{law0040.Log-Pareto-tail-normal}
\alias{law0040.Log-Pareto-tail-normal}
\alias{law0040}

\title{The Log-Pareto-tail-normal Distribution}

\description{Random generation for the Log-Pareto-tail-normal
  distribution with parameters \code{alpha}, \code{mu} and \code{sigma}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{alpha}, \code{mu} and \code{sigma} are not specified
  they assume the default values of 1.959964, 0.0 and 1.0 respectively.

The log-Pareto-tailed normal distribution has a symmetric and continuous density that belongs to the larger family of log-regularly varying distributions
(see Desgagne, 2015). This is essentially a normal density with log-Pareto tails. Using this distribution instead of the usual normal ensures whole robustness to outliers
in the estimation of location and scale parameters and in the estimation of parameters of a multiple linear regression. 

The density of the log-Pareto-tailed normal distribution with parameters
\code{alpha}, \code{mu} and 
\code{sigma} is given by

\deqn{g(y\mid\alpha,\mu,\sigma)=\left\{
  \begin{array}{ccc}
     \frac{1}{\sigma}\phi\left(\frac{y-\mu}{\sigma}\right) & \textrm{ if } & \mu - \alpha\sigma \le  y\le \mu + \alpha\sigma, \\
      &\\
     \phi(\alpha)\frac{\alpha}{|y-\mu|}\left(\frac{\log \alpha}{\log (|y-\mu|/\sigma)}\right)^\beta &  \textrm{ if } & |y-\mu|\ge \alpha\sigma,
   \end{array}
   \right.}{g(y|\alpha,\mu,\sigma) =
  (1/\sigma)\phi((y-\mu)(\sigma)) if \mu - \alpha \sigma <= y <= \mu +
  \alpha \sigma, and g(y|\alpha,\mu,\sigma) =  \phi(\alpha)
  (\alpha)/(|y-\mu|) ((\log \alpha)/(\log
  (|y-\mu|/\sigma)))^\beta if |y-\mu| >= \alpha \sigma.}

where \eqn{\beta = 1+2\,\phi(\alpha)\,\alpha\log(\alpha)(1-q)^{-1}} and \eqn{q=\Phi(\alpha)-\Phi(-\alpha)}. 
The functions \eqn{\phi(\alpha)=\frac{1}{\sqrt{2\pi}}\exp[-\frac{\alpha^2}{2}]} and \eqn{\Phi(\alpha)} are respectively the p.d.f. and the c.d.f. of the standard normal distribution.
The domains of the variable and the parameters are \eqn{-\infty<y<\infty}, \eqn{\alpha>1}, \eqn{-\infty<\mu<\infty} and \eqn{\sigma>0}.
 
Note that the normalizing constant \eqn{K_{(\alpha,\beta)}} (see Desgagne, 2015, Definition 3) has been set to 1. The desirable consequence is that 
the core of the density, between \eqn{\mu-\alpha\sigma}  and \eqn{\mu+\alpha\sigma}, becomes exactly the density of the \eqn{N(\mu,\sigma^2)}. This mass of the density
corresponds to \eqn{q}{q}. It follows that the parameter \eqn{\beta} is no longer free and its value depends on \eqn{\alpha} as given above. 

For example, if we set \eqn{\alpha=1.959964}, we obtain \eqn{\beta=4.083613} and \eqn{q=0.95} of the mass is comprised between \eqn{\mu-\alpha\sigma}  and \eqn{\mu+\alpha\sigma}.
Note that if one is more comfortable in choosing the central mass $q$ instead of choosing directly the parameter \eqn{\alpha}, then it suffices to use the equation \eqn{\alpha=\Phi^{-1}((1+q)/2)}, with the
contrainst \eqn{q>0.6826895\Leftrightarrow \alpha>1}{q>0.6826895 <==> \alpha>1}.


The mean and variance of Log-Pareto-tail-normal are not defined.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Desgagne, Alain. Robustness to outliers in location-scale parameter model using log-regularly varying distributions. \emph{Ann. Statist.} \bold{43} (2015), no. 4, 1568--1595. doi:10.1214/15-AOS1316. http://projecteuclid.org/euclid.aos/1434546215. 
}
\author{P. Lafaye de Micheaux}

\seealso{See \code{\link{Distributions}} for other standard distributions.}

\examples{
res <- gensample(40, 10000, law.pars = c(1.959964, 0.0, 1.0))
res$law
res$law.pars
}

\keyword{law, Log-Pareto-tail-normal}
