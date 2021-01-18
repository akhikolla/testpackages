\name{tess.analysis}
\alias{tess.analysis}
\alias{globalBiDe.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.analysis: Diversification rate estimation under an episodic birth-death process including mass-extinction events.}
\description{tess.analysis estimates diversification rates under an episodic birth-death process including mass-extinction events. The method uses a reversible-jump MCMC algorithm to estimate the number, timing and magnitude of rate-shifts and mass-extinction events. It is possible to fix either number of events and provide specific values that will be used. We assume a Poison process for the number of events and a lognormal distribution with fixed, but specified, hyper-parameters for the speciation and extinction rate; and an independent Poison process for the number of mass-extinction events where each survival probability follows a Beta distribution with fixed hyper-parameters.

The MCMC algorithm can be run either for a specified number of iterations, until a time limit in seconds has been reached, or until the effective sample size (ESS) has reached a given threshold. Once the first of these requirements are met TESS will stop the analysis. Internally we use scaling and sliding proposals to change the parameter values during the MCMC and a birth-move and death-move to add/remove events (rate-shifts or mass-extinction events).

The results of the MCMC run are stored within a directory that is specified by the user. Several files will be generated containing the sampled parameter values. To summarize the output see tess.process.output(...) and tess.plot.output(...).
}
\usage{
tess.analysis( tree,
               initialSpeciationRate,
               initialExtinctionRate,
               empiricalHyperPriors = TRUE,
               empiricalHyperPriorInflation = 10.0,
               empiricalHyperPriorForm = c("lognormal","normal","gamma"),
               speciationRatePriorMean = 0.0,
               speciationRatePriorStDev = 1.0,
               extinctionRatePriorMean = 0.0,
               extinctionRatePriorStDev = 1.0,
               initialSpeciationRateChangeTime = c(),
               initialExtinctionRateChangeTime = c(),
               estimateNumberRateChanges = TRUE,
               numExpectedRateChanges = 2,
               samplingProbability = 1,
               missingSpecies = c(),
               timesMissingSpecies = c(),
               tInitialMassExtinction = c(),
               pInitialMassExtinction = c(),
               pMassExtinctionPriorShape1 = 5,
               pMassExtinctionPriorShape2 = 95,
               estimateMassExtinctionTimes = TRUE,
               numExpectedMassExtinctions = 2,
               estimateNumberMassExtinctions = TRUE,
               MRCA = TRUE,
               CONDITION = "survival",
               BURNIN = 10000,
               MAX_ITERATIONS = 200000,
               THINNING = 100,
               OPTIMIZATION_FREQUENCY = 500,
               CONVERGENCE_FREQUENCY = 1000,
               MAX_TIME = Inf, MIN_ESS = 500,
               ADAPTIVE = TRUE,
               dir = "" ,
               priorOnly = FALSE,
               verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{The tree in 'phylo' format.}
  \item{initialSpeciationRate}{The initial value of the speciation rate when the MCMC is started. This can either be a single number of a vector of rates per interval.}
  \item{initialExtinctionRate}{The initial value of the extinction rate when the MCMC is started. This can either be a single number of a vector of rates per interval.}
  \item{empiricalHyperPriors}{Should we estimate the hyper-parameters empirically?}
  \item{empiricalHyperPriorInflation}{The scaling factor of the variance for the empirical hyperpriors.}
  \item{empiricalHyperPriorForm}{The possible empirical hyper prior distributions; either lognormal, normal or gamma}
  \item{speciationRatePriorMean}{The mean of the log-normal prior distribution for the speciation rate.}
  \item{speciationRatePriorStDev}{The standard deviation of the log-normal prior distribution for the speciation rate.}
  \item{extinctionRatePriorMean}{The mean of the log-normal prior distribution for the extinction rate.}
  \item{extinctionRatePriorStDev}{The standard deviation of the log-normal prior distribution for the extinction rate.}
  \item{initialSpeciationRateChangeTime}{The initial value of the time points when speciation rate-shifts occur. The number of time-shifts needs to be one smaller than the number of initial speciation rates.}
  \item{initialExtinctionRateChangeTime}{The initial value of the time points when extinction rate-shifts occur. The number of time-shifts needs to be one smaller than the number of initial extinction rates.}
  \item{estimateNumberRateChanges}{Do we estimate the number of rate shifts? Default is true.}
  \item{numExpectedRateChanges}{Expected number of rate changes which follow a Poisson process. The default gives 0.5 probability on 0 shifts.}
  \item{samplingProbability}{The extant taxa sampling probability at the present time.}
  \item{missingSpecies}{The number of species missed which originated in a given time interval (empirical taxon sampling).}
  \item{timesMissingSpecies}{The times intervals of the missing species (empirical taxon sampling).}
  \item{tInitialMassExtinction}{The initial value of the vector of times of the mass-extinction events. This is used as initial values for the MCMC.}
  \item{pInitialMassExtinction}{The initial value of the vector of survival probabilities of the mass-extinction events. This is used as initial values for the MCMC.}
  \item{pMassExtinctionPriorShape1}{The alpha (first shape) parameter of the Beta prior distribution for the survival probability of a mass-extinction event.}
  \item{pMassExtinctionPriorShape2}{The beta (second shape) parameter of the Beta prior distribution for the survival probability of a mass-extinction event.}
  \item{estimateMassExtinctionTimes}{Do we estimate the times of mass-extinction events? Default is true.}
  \item{numExpectedMassExtinctions}{Expected number of mass-extinction events which follow a Poisson process. The default gives 0.5 probability on 0 events.}
  \item{estimateNumberMassExtinctions}{Do we estimate the number of mass-extinction events? Default is true.}
  \item{MRCA}{Does the process start with the most recent common ancestor? If not, the tree must have a root edge!}
  \item{CONDITION}{do we condition the process on time|survival|taxa?}
  \item{BURNIN}{The length of the burnin period.}
  \item{MAX_ITERATIONS}{The maximum number of iteration of the MCMC. The default is 200000.}
  \item{THINNING}{The frequency how often samples are recorded during the MCMC. The default is every 100 iterations.}
  \item{OPTIMIZATION_FREQUENCY}{The frequency how often the MCMC moves are optimized. The default is every 500 iterations.}
  \item{CONVERGENCE_FREQUENCY}{The frequency how often we check for convergence? The default is every 1000 iterations.}
  \item{MAX_TIME}{The maximum time the MCMC is allowed to run in seconds. The default is Inf}
  \item{MIN_ESS}{The minimum number of effective samples (ESS) to assume convergence. The default is 500}
  \item{ADAPTIVE}{Do we use auto-tuning of the MCMC moves? The default is TRUE (recommended).}
  \item{dir}{The subdirectory in which the output will be stored. The default is the present directoy ("")}
  \item{priorOnly}{Do we sample from the prior only? The default is FALSE}
  \item{verbose}{Do you want detailed output?}
}
\value{There is no return value because all the results are stored into files.}
\references{
S. Hoehna: The time-dependent reconstructed evolutionary process with a key-role for mass-extinction events. 2015, Journal of Theoretical Biology, 380, 321-331.

S. Hoehna, MR May and BR Moore: TESS: Bayesian inference of lineage diversification rates from (incompletely sampled) molecular phylogenies in R. 2015, Bioinformatics.

MR May, S. Hoehna, and BR Moore: A Bayesian approach for detecting mass-extinction events when rates of lineage diversification vary. 2015, Systematic Biology
}
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# we load the conifers as the test data set
data(conifers)

# for the conifers we know what the total number of species is
total <- 630
# thus, we can compute what the sampling fraction is
rho <- (conifers$Nnode+1)/total


# next, we specify the prior mean and standard deviation 
# for the speciation and extinction rate
mu_lambda = 0.15
std_lambda = 0.02
mu_mu = 0.09
std_mu = 0.02

# now we can run the entire analysis.
# note that a full analyses should be run much longer
tess.analysis( tree=conifers,
               initialSpeciationRate=exp(mu_lambda),
               initialExtinctionRate=exp(mu_mu),
               empiricalHyperPriors = FALSE,
               speciationRatePriorMean = mu_lambda,
               speciationRatePriorStDev = std_lambda,
               extinctionRatePriorMean = mu_mu,
               extinctionRatePriorStDev = std_mu,
               numExpectedRateChanges = 2,
               samplingProbability = rho,
               numExpectedMassExtinctions = 2,
               BURNIN = 100,
               MAX_ITERATIONS = 200,
               THINNING = 10,
               dir = "analysis_conifer")
               
# You may want to look into the vignette for a more detailed description
# of the features for an analysis.
# also have a look at the functions tess.process.output and tess.plot.output

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line
