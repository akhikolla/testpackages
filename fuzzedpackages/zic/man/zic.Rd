\name{zic}

\alias{zic}

\title{Bayesian Inference for Zero-Inflated Count Models} 

\description{
\code{zic} fits zero-inflated count models via Markov chain Monte Carlo methods.
}

\usage{
zic(formula, data, a0, b0, c0, d0, e0, f0, 
    n.burnin, n.mcmc, n.thin, tune = 1.0, scale = TRUE)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit specifying the response variable and covariates.}
  \item{data}{A data frame in which to interpret the variables in \code{formula}.}
  \item{a0}{The prior variance of \eqn{\alpha}{alpha}.}
  \item{b0}{The prior variance of \eqn{\beta_j}{beta_j}.}
  \item{c0}{The prior variance of \eqn{\gamma}{gamma}.}
  \item{d0}{The prior variance of \eqn{\delta_j}{delta_j}.}
  \item{e0}{The shape parameter for the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{f0}{The inverse scale parameter the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{n.burnin}{Number of burn-in iterations of the sampler.}
  \item{n.mcmc}{Number of iterations of the sampler.}
  \item{n.thin}{Thinning interval.}
  \item{tune}{Tuning parameter of Metropolis-Hastings step.}
  \item{scale}{If true, all covariates (except binary variables) are rescaled by dividing by their respective standard errors.}
}

\value{
   A list containing the following elements:
   \item{alpha}{Posterior draws of \eqn{\alpha}{alpha} (coda mcmc object).}
   \item{beta}{Posterior draws of \eqn{\beta}{beta} (coda mcmc object) .}
   \item{gamma}{Posterior draws of \eqn{\gamma}{gamma} (coda mcmc object).}	
   \item{delta}{Posterior draws of \eqn{\delta}{delta} (coda mcmc object).}
   \item{sigma2}{Posterior draws of \eqn{\sigma^2}{sigma^2} (coda mcmc object).}
   \item{acc}{Acceptance rate of the Metropolis-Hastings step.}
}

\details{
  The considered zero-inflated count model is given by
  \deqn{y_i^* \sim \mathrm{Poisson}[\exp(\eta_i^*)],}{%
        y*_i ~ Poisson[exp(eta*_i)],}
  \deqn{\eta^*_i = \alpha + x_i'\beta + \varepsilon_i,\; \varepsilon_i \sim \mathrm{N}(0,\sigma^2),}{%
        eta*_i = x_i' * beta + epsilon_i, epsilon_i ~ N( 0, sigma^2 ),}
  \deqn{d_i^* = \gamma + x_i'\delta + \nu_i,\; \nu_i \sim \mathrm{N}(0,1),}{%
        d*_i = x_i' * delta + nu_i, nu_i ~ N(0,1),}
  \deqn{y_i = 1(d_i^*>0)y_i^*,}{%
        y_i = 1(d*_i>0) y*_i,} 
  where \eqn{y_i} and \eqn{x_i} are observed. The assumed prior distributions are
  \deqn{\alpha \sim \mathrm{N}(0,a_0),}{%
         alpha ~ N(0,a0),}
  \deqn{\beta_k \sim \mathrm{N}(0,b_0), \quad k=1,\ldots,K,}{%
         beta_k ~ N(0,b0), k=1,...,K,}
  \deqn{\gamma \sim \mathrm{N}(0,c_0),}{%
         gamma ~ N(0,c0)}
  \deqn{\delta_k \sim \mathrm{N}(0,d_0), \quad k=1,\ldots,K,}{%
         delta_k ~ N(0,d0), k=1,...,K,}
  \deqn{\sigma^2 \sim \textrm{Inv-Gamma}\left(e_0,f_0\right).}{%
        sigma^2 ~ Inv-Gamma(e0,f0).}

  The sampling algorithm described in Jochmann (2013) is used.
}

\examples{
\dontrun{
data( docvisits )
mdl <- docvisits ~ age + agesq + health + handicap + hdegree + married + schooling +
                    hhincome + children + self + civil + bluec + employed + public + addon
post <- zic( f, docvisits, 10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1000, 10000, 10, 1.0, TRUE )}
}

\references{
  Jochmann, M. (2013). ``What Belongs Where? Variable Selection for Zero-Inflated Count Models with an Application to the Demand for
  Health Care'', \emph{Computational Statistics}, 28, 1947--1964.
}
