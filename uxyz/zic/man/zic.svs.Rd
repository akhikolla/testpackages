\name{zic.svs}

\alias{zic.svs}

\title{SVS for Zero-Inflated Count Models} 

\description{
\code{zic.svs} applies SVS to zero-inflated count models
}

\usage{
zic.svs(formula, data,
        a0, g0.beta, h0.beta, nu0.beta, r0.beta, s0.beta, e0, f0, 
        c0, g0.delta, h0.delta, nu0.delta, r0.delta, s0.delta, 
        n.burnin, n.mcmc, n.thin, tune = 1.0, scale = TRUE)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit specifying the response variable and covariates.}
  \item{data}{A data frame in which to interpret the variables in \code{formula}.}
  \item{a0}{The prior variance of \eqn{\alpha}{alpha}.}
  \item{g0.beta}{The shape parameter for the inverse gamma prior on \eqn{\kappa_k^\beta}{kappa_k^beta}.}
  \item{h0.beta}{The inverse scale parameter for the inverse gamma prior on \eqn{\kappa_k^\beta}{kappa_k^beta}.}
  \item{nu0.beta}{Prior parameter for the spike of the hypervariances for the \eqn{\beta_k}{beta_k}.}
  \item{r0.beta}{Prior parameter of \eqn{\omega^\beta}{omega^beta}.}
  \item{s0.beta}{Prior parameter of \eqn{\omega^\beta}{omega^beta}.}
  \item{e0}{The shape parameter for the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{f0}{The inverse scale parameter the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{c0}{The prior variance of \eqn{\gamma}{gamma}.}	  
  \item{g0.delta}{The shape parameter for the inverse gamma prior on \eqn{\kappa_k^\delta}{kappa_k^delta}.}
  \item{h0.delta}{The inverse scale parameter for the inverse gamma prior on \eqn{\kappa_k^\delta}{kappa_k^delta}.}
  \item{nu0.delta}{Prior parameter for the spike of the hypervariances for the \eqn{\delta_k}{delta_k}.}
  \item{r0.delta}{Prior parameter of \eqn{\omega^\delta}{omega^delta}.}
  \item{s0.delta}{Prior parameter of \eqn{\omega^\delta}{omega^delta}.}
  \item{n.burnin}{Number of burn-in iterations of the sampler.}
  \item{n.mcmc}{Number of iterations of the sampler.}
  \item{n.thin}{Thinning interval.}
  \item{tune}{Tuning parameter of Metropolis-Hastings step.}
  \item{scale}{If true, all covariates (except binary variables) are rescaled by dividing by their respective standard errors.}
}

\value{
   A list containing the following elements:
   \item{alpha}{Posterior draws of \eqn{\alpha}{alpha} (coda mcmc object).}
   \item{beta}{Posterior draws of \eqn{\beta}{beta} (coda mcmc object).}
   \item{gamma}{Posterior draws of \eqn{\gamma}{gamma} (coda mcmc object).}	
   \item{delta}{Posterior draws of \eqn{\delta}{delta} (coda mcmc object).}
   \item{sigma2}{Posterior draws of \eqn{\sigma^2}{sigma^2} (coda mcmc object).}
   \item{I.beta}{Posterior draws of indicator whether \eqn{\tau_j^\beta}{tau_j^beta} is one (coda mcmc object).}
   \item{I.delta}{Posterior draws of indicator whether \eqn{\tau_j^\delta}{tau_j^delta} is one (coda mcmc object).}
   \item{omega.beta}{Posterior draws of \eqn{\omega^\beta}{omega^beta} (coda mcmc object).}
   \item{omega.delta}{Posterior draws of \eqn{\omega^\delta}{omega^delta} (coda mcmc object).}
   \item{acc}{Acceptance rate of the Metropolis-Hastings step.}
}

\details{
  The considered zero-inflated count model is given by
  \deqn{y_i^* \sim \mathrm{Poisson}[\exp(\eta^*_i)],}{%
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
  \deqn{\beta_k\sim \mathrm{N}(0,\tau^\beta_k\kappa^\beta_k),, \quad k=1,\ldots,K,}{%
        beta_k ~ N(0, tau_k^beta * kappa_k^beta), k=1,...,K}
  \deqn{\kappa^\beta_j\sim\textrm{Inv-Gamma}(g_0^\beta,h_0^\beta),}{%
        kappa_k^\beta ~ Inv-Gamma(g0^beta,h0^beta),}
  \deqn{\tau_k^\beta \sim (1-\omega^\beta)\delta_{\nu^\beta_0}+\omega^\beta\delta_1,}{%
        tau_k^\beta ~ (1-\omega^beta) delta_(nu0^\beta)+ omega^beta delta_1,}
  \deqn{\omega^\beta\sim\mathrm{Beta}(r_0^\beta,s_0^\beta),}{%
        omega^beta ~ Beta(r0^beta,s0^beta),}
  \deqn{\gamma \sim \mathrm{N}(0,c_0),}{%
         gamma ~ N(0,c0)}
  \deqn{\delta_k\sim \mathrm{N}(0,\tau^\delta_k\kappa^\delta_k), \quad k=1,\ldots,K,}{%
        delta_k ~ N(0, tau_k^delta * kappa_k^delta), k=1,...,K,}
  \deqn{\kappa^\delta_k\sim\textrm{Inv-Gamma}(g_0^\delta,h_0^\delta),}{%
        kappa_k^\delta ~ Inv-Gamma(g0^delta,h0^delta),}
  \deqn{\tau_k^\delta \sim (1-\omega^\delta)\delta_{\nu^\delta_0}+\omega^\delta\delta_1,}{%
        tau_k^\delta ~ (1-\omega^delta) delta_(nu_0^\delta)+ omega^delta delta_1,}
  \deqn{\omega^\delta\sim\mathrm{Beta}(r_0^\delta,s_0^\delta),}{%
        omega^delta ~ Beta(r0^delta,s0^delta),}
  \deqn{\sigma^2 \sim \textrm{Inv-Gamma}\left(e_0,f_0\right).}{%
        sigma^2 ~ Inv-Gamma(e0,f0).}

  The sampling algorithm described in Jochmann (2013) is used.
}

\references{
  Jochmann, M. (2013). ``What Belongs Where? Variable Selection for Zero-Inflated Count Models with an Application to the Demand for
  Health Care'', \emph{Computational Statistics}, 28, 1947--1964.
}

\examples{
\dontrun{
data( docvisits )
mdl <- docvisits ~ age + agesq + health + handicap + hdegree + married + schooling +
                    hhincome + children + self + civil + bluec + employed + public + addon
post <- zic.ssvs( mdl, docvisits,
                  10.0, 5.0, 5.0, 1.0e-04, 2.0, 2.0, 1.0, 1.0,
                  10.0, 5.0, 5.0, 1.0e-04, 2.0, 2.0,
                  1000, 10000, 10, 1.0, TRUE )}}