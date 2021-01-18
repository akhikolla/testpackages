\name{fecr_stan}
\alias{fecr_stan}
\title{Model the reduction of faecal egg count}
\usage{
fecr_stan(preFEC, postFEC, rawCounts = FALSE, preCF = 50, postCF = preCF, 
  paired = TRUE, indEfficacy = TRUE, zeroInflation = FALSE, 
  muPrior, kappaPrior, deltaPrior, phiPrior, deltakappaPrior,
  nsamples = 2000, nburnin = 1000, thinning = 1, nchain = 2, 
  ncore = 1, adaptDelta = 0.95, saveAll = FALSE, verbose = FALSE)
}
\arguments{
\item{preFEC}{numeric vector. Pre-treatment faecal egg counts.}

\item{postFEC}{numeric vector. Post-treatment faecal egg counts.}

\item{rawCounts}{logical. If TRUE, \code{preFEC} and \code{postFEC} correspond to raw counts
(as counted on equipment). Otherwise they correspond to calculated epgs (raw counts times correction factor).
Defaults to \code{FALSE}.}

\item{preCF}{a positive integer or a vector of positive integers. Pre-treatment correction factor(s).}

\item{postCF}{a positive integer or a vector of positive integers. Post-treatment correction factor(s).}

\item{paired}{logical. If TRUE, uses the model for the paired design. Otherwise uses the model for the unpaired design}

\item{indEfficacy}{logical. If TRUE, uses the paired model allowing for individual efficacy. Only use in combination with \code{paired = TRUE} and \code{zeroInflation = FALSE}.}

\item{zeroInflation}{logical. If TRUE, uses the model with zero-inflation. Otherwise uses the model without zero-inflation.}

\item{muPrior}{a named list. Prior for the group mean epg parameter \eqn{\mu}. The default prior is \code{list(priorDist = "gamma", hyperpars = c(1,0.001))}, i.e. a gamma distribution with shape 1 and rate 0.001, its 90\% probability mass lies between 51 and 2996.}

\item{kappaPrior}{a named list. Prior for the group dispersion parameter \eqn{\kappa}. The default prior is \code{list(priorDist = "gamma", hyperpars = c(1,0.7))}, i.e. a gamma distribution with shape 1 and rate 0.7, its 90\% probability mass lies between 0.1 and 4.3 with a median of 1.}

\item{deltaPrior}{a named list. Prior for the reduction parameter \eqn{\delta}. The default prior is \code{list(priorDist = "beta", hyperpars = c(1,1))}, i.e. a uniform prior between 0 and 1.}

\item{phiPrior}{a named list. Prior for the zero-inflation parameter \eqn{\phi}. The default prior is \code{list(priorDist = "beta", hyperpars = c(1,1))}, i.e. a uniform prior between 0 and 1.}

\item{deltakappaPrior}{a named list. Prior information for the shape parameter of reduction \eqn{\delta_\kappa}. The default prior is \code{list(priorDist = "normal",hyperpars=c(2,1))}. Only used if \code{indEfficacy = TRUE}.}

\item{nsamples}{a positive integer. Number of samples for each chain (including burn-in samples).}

\item{nburnin}{a positive integer. Number of burn-in samples.}

\item{thinning}{a positive integer. Thinning parameter, i.e. the period for saving samples.}

\item{nchain}{a positive integer. Number of chains.}

\item{ncore}{a positive integer. Number of cores to use when executing the chains in parallel.}

\item{adaptDelta}{numeric. The target acceptance rate, a numeric value between 0 and 1.}

\item{saveAll}{logical. If TRUE, posterior samples for all parameters are saved in the \code{stanfit} object. If FALSE, only samples for \eqn{\delta}, \eqn{\mu}, \eqn{\kappa} and \eqn{\phi} are saved. Default to FALSE.}

\item{verbose}{logical. If TRUE, prints progress and debugging information.}
}
\value{
Prints out the posterior summary of \code{FECR} as the reduction, \code{meanEPG.untreated} as the mean pre-treatment epg, and \code{meanEPG.treated} as the mean after-treatment epg. The posterior summary contains the mean, standard deviation (sd), 2.5\%,  50\% and 97.5\% percentiles, the 95\% highest posterior density interval (HPDLow95 and HPDHigh95) and the posterior mode. 

NOTE: Based on our simulation studies, we recommend to use (2.5\%, 97.5\%) as the 95\% credible interval and the median as summary statistics of reduction for the individual efficacy model. For all other models, we recommend to use the 95\% HPD interval and the mode. 

The returned value is a list that consists of:
\item{stan.samples}{an object of S4 class \code{\link[rstan]{stanfit}} representing the fitted results}
\item{posterior.summary}{a data.frame that is the same as the printed posterior summary}
}
\description{
Models the reduction in faecal egg counts with Bayesian hierarchical models. See Details for a list of model choices.
}
\details{
\subsection{List of built-in models}{
\itemize{
\item unpaired without zero-inflation: set \code{paired = FALSE}, \code{indEfficacy = FALSE}, \code{zeroInflation = FALSE}
\item unpaired with zero-inflation: set \code{paired = FALSE}, \code{indEfficacy = FALSE}, \code{zeroInflation = TRUE}
\item paired without zero-inflation: set \code{paired = TRUE}, \code{indEfficacy = FALSE}, \code{zeroInflation = FALSE}
\item paired with zero-inflation: set \code{paired = TRUE}, \code{indEfficacy = FALSE}, \code{zeroInflation = TRUE}
\item paired with individual efficacy: set \code{paired = TRUE}, \code{indEfficacy = TRUE}, \code{zeroInflation = FALSE}
}}
\subsection{Prior choice}{

Consider using non-default prior for \eqn{\delta} when,
\itemize{
\item there is on average an increase in egg counts after treatment
\item there are divergent-sample warnings
\item there are non-convergence warnings
}
Two examples of useful non-default priors include:
\enumerate{
\item \code{list(priorDist = "normal", hyperpars = c(1, 5))} for stablizing the reduction parameter without being informative. 
\item \code{list(priorDist = "beta", hyperpars = c(0, 5))} for allowing up to 4-fold increase of egg count after treatment.
}}
\subsection{Other information}{
The first time each model with non-default priors is applied, it can take up to 20 seconds to compile the model. Currently the function only support prior distributions with two parameters. For a complete list of supported priors and their parameterization, please consult the list of distributions in \href{http://mc-stan.org/documentation/}{Stan}. 

The default number of samples per chain is 2000, with 1000 burn-in samples. Normally this is sufficient in Stan. If the chains do not converge, one should tune the MCMC parameters until convergence is reached to ensure reliable results.
}}

\references{
\strong{Individual efficacy models:}
Craig Wang, Paul R. Torgerson, Ray M. Kaplan, Melissa M. George, Reinhard Furrer. (2018)
Modelling anthelmintic resistance by extending eggCounts package to allow individual efficacy,
International Journal for Parasitology: Drugs and Drug Resistance, Volume 8, Pages 386-393. \href{https://doi.org/10.1016/j.ijpddr.2018.07.003}{https://doi.org/10.1016/j.ijpddr.2018.07.003}

\strong{Zero-inflation models:}
Craig Wang, Paul R. Torgerson, Johan Hoglund, Reinhard Furrer. (2017)
Zero-inflated hierarchical models for faecal egg counts to assess anthelmintic efficacy,
Veterinary Parasitology, Volume 235, Pages 20-28. \href{http://dx.doi.org/10.1016/j.vetpar.2016.12.007}{http://dx.doi.org/10.1016/j.vetpar.2016.12.007}

\strong{Other models:}
Paul R. Torgerson, Michaela Paul, Reinhard Furrer. (2014)
Evaluating faecal egg count reduction using a specifically designed package 'eggCounts' in {R} and a user friendly web interface,
International Journal for Parasitology, Volume 44, Pages 299-303. \href{http://dx.doi.org/10.1016/j.ijpara.2014.01.005}{http://dx.doi.org/10.1016/j.ijpara.2014.01.005}
}

\author{ 
Craig Wang 
}
\seealso{
\code{\link{simData2s}} for simulating faecal egg counts data with two samples 
}
\examples{
## load sample data
data(epgs)

## apply individual efficacy model to the data vectors
model <- fecr_stan(epgs$before, epgs$after, rawCounts = FALSE, preCF = 50, 
                   paired = TRUE, indEfficacy = TRUE)
## convert to MCMC object and inspect the summary
samples <- stan2mcmc(model$stan.samples)
summary(samples)
}
\keyword{models}
