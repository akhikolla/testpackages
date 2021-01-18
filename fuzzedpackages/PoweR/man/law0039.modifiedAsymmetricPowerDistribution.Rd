\name{law0039.modifiedAsymmetricPowerDistribution}
\alias{law0039.modifiedAsymmetricPowerDistribution}
\alias{law0039}

\title{The modified Asymmetric Power Distribution}

\description{Random generation for the modified Asymmetric Power Distribution with parameters \code{theta}, \code{phi}, \code{alpha} and \code{lambda}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{theta}, \code{phi}, \code{alpha} and \code{lambda} are not specified they assume the default values of 0, 1, 0.5 and 2, respectively.

The modified Asymmetric Power Distribution with parameters \code{theta}, \code{phi}, \code{theta1} and \code{theta2} has density:
\deqn{f(x\mid\boldsymbol{\theta})=\frac{(\delta_{\boldsymbol{\theta}}/2)^{1/\theta_2}}{\Gamma(1+1/\theta_2)}\exp\left[-\left(\frac{2(\delta_{\boldsymbol{\theta}}/2)^{1/\theta_2}}{1+sign(x)(1-2\theta_1)}|x|\right)^{\theta_2}\right]}{f(x|\theta) = ([(\delta_{\theta}/2)^{1/\theta_2}] / [\Gamma(1+1/\theta_2)]) \exp[-(((2(\delta_{\theta}/2)^{1/\theta_2}) / (1+sign(x)(1-2\theta_1))) * |x|)^{\theta_2}]} where \eqn{\boldsymbol{\theta}=(\theta_2, \theta_1)^T}{\theta = (\theta_2, \theta_1)^T} is the vector of parameters, \eqn{\theta_2>0, 0<\theta_1<1} and \deqn{\delta_{\boldsymbol{\theta}}=\frac{2(\theta_1)^{\theta_2} (1-\theta_1)^{\theta_2}}{(\theta_1)^{\theta_2}+(1-\theta_1)^{\theta_2}}}{\delta_{\theta} = (2(\theta_1)^{\theta_2} (1-\theta_1)^{\theta_2}) / ((\theta_1)^{\theta_2}+(1-\theta_1)^{\theta_2})}.


The mean and variance of APD are defined respectively by
\deqn{E(U) = \theta + 2 ^ {1 / \theta_2}  \phi \Gamma(2 / \theta_2)  (1 -
  2 \theta_1)  \delta ^ {-1 / \theta_2} / \Gamma(1 / \theta_2)} and  
\deqn{V(U) = 2 ^ {2 / \theta_2} \phi ^ 2  \left(\Gamma(3 / \theta_2)  \Gamma(1 / \theta_2)  (1 - 3  \theta_1 + 3  \theta_1 ^ 2) - \Gamma^2(2 / \theta_2)  (1 - 2  \theta_1) ^ 2\right) \delta ^ {-2 / \theta_2} / \Gamma^2(1 / \theta_2).}
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Desgagne, A. and Lafaye de Micheaux, P. and Leblanc, A. (2016), Test of
normality based on alternate measures of skewness and kurtosis, 
\emph{}, 
}
\author{P. Lafaye de Micheaux}

\seealso{See \code{\link{Distributions}} for other standard distributions.}

\examples{
res <- gensample(39, 10000, law.pars = c(3, 2, 0.5, 1))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, modified Asymmetric Power Distribution}
