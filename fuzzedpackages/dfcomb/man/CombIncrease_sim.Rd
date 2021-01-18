\name{CombIncrease_sim}
\alias{CombIncrease_sim}
\alias{print.CombIncrease_sim}
\title{Combination design Simulator using Logistic model}
\description{
\code{CombIncrease_sim} is used to generate simulation replicates of phase I clinical trial for combination studies where the toxicity and efficacy of both agents is assumed to increase with the dose using the design proposed by Riviere et al. entitled "A Bayesian dose-finding design for drug combination clinical trials based on the logistic model".
}
\usage{
CombIncrease_sim(ndose_a1, ndose_a2, p_tox, target, target_min, target_max,
  prior_tox_a1, prior_tox_a2, n_cohort, cohort, tite=FALSE, time_full=0,
  poisson_rate=0, nsim, c_e=0.85, c_d=0.45, c_stop=0.95, c_t=0.5, c_over=0.25,
  cmin_overunder=2, cmin_mtd=3, cmin_recom=1, startup=1, alloc_rule=1,
  early_stop=1, nburn=2000, niter=5000, seed=14061991)
}
\arguments{
  \item{ndose_a1}{Number of dose levels for agent 1.}
  \item{ndose_a2}{Number of dose levels for agent 2.}
  \item{p_tox}{A matrix of the true toxicity probabilities associated with the combinations. True toxicity probabilities should be entered with agent 1 in row and agent 2 in column, with increasing toxicity probabilities with both row and column numbers (see examples).}
  \item{target}{Toxicity (probability) target.}
  \item{target_min}{Minimum of the targeted toxicity interval.}
  \item{target_max}{Maximum of the targeted toxicity interval.}
  \item{prior_tox_a1}{A vector of initial guesses of toxicity probabilities associated with the doses of agent 1. Must be of length \code{ndose_a1}.}
  \item{prior_tox_a2}{A vector of initial guesses of toxicity probabilities associated with the doses of agent 2. Must be of length \code{ndose_a2}.}
  \item{n_cohort}{Total number of cohorts to include in the trial.}
  \item{cohort}{Cohort size.}
  \item{tite}{A boolean indicating if the toxicity is considered as a time-to-event outcome (TRUE), or as a binary outcome (default value FALSE).}
  \item{time_full}{Full follow-up time window. This argument is used only if tite=TRUE.}
  \item{poisson_rate}{A value indicating the rate for the Poisson process used to simulate patient arrival, i.e. expected number of arrivals per observation window. This argument is used only if tite=TRUE.}
  \item{nsim}{Number of simulations.}
  \item{c_e}{Probability threshold for dose-escalation. The default value is set at 0.85.}
  \item{c_d}{Probability threshold for dose-deescalation. The default value is set at 0.45.}
  \item{c_stop}{Probability threshold for early trial termination due to over-toxicity or under-toxicity (see details). The default value is set at 0.95.}
  \item{c_t}{Probability threshold for early trial termination for finding the MTD (see details). The default value is set at 0.5.}
  \item{c_over}{Probability threshold to control over-dosing (see details).}
  \item{cmin_overunder}{Minimum number of cohorts to be included at the lowest/highest combination before possible early trial termination for over-toxicity or under-toxicity (see details). The default value is set at 2.}
  \item{cmin_mtd}{Minimum number of cohorts to be included at the recommended combination before possible early trial termination for finding the MTD (see details). The default value is set at 3.}
  \item{cmin_recom}{Minimum number of cohorts to be included at the recommended combination at the end of the trial. The default value is set at 1.}
  \item{startup}{Interger (0, 1, 2, or 3) indicating which start-up phase is used (see details). The default value is set at 1.}
  \item{alloc_rule}{Interger (1, 2, or 3) indicating which allocation rule is used (see details). The default value is set at 1.}
  \item{early_stop}{Interger (1, 2, or 3) indicating which early stopping rule is used (see details). The default value is set at 1.}
  \item{nburn}{Number of burn-in for HMC. The default value is set at 2000.}
  \item{niter}{Number of iterations for HMC. The default value is set at 5000.}
  \item{seed}{Seed of the random number generator. Default value is set at 14061991.}}
\details{
Start-up phase:
\itemize{
  \item \code{startup=0}: No startup phase.
  \item \code{startup=1} (Riviere et al 2014): Begin at the lowest combination (1,1) and increase both agent (+1, +1) until the first toxicity is observed or maximum combination is reached.
  \item \code{startup=2}: Begin at the lowest combination (1,1) and increase agent 1 (+1, 0) until a toxicity is observed or maximum dose is reached. Then begin at (1,2) and increase agent 2 (0, +1) until a toxicity is observed or maximum dose is reached.
  \item \code{startup=3}: Begin at the lowest combination (1,1) and increase alternatively each agent (+1, 0) then (0, +1) until the first toxicity is observed or maximum combination is reached.
}

Allocation rule:
\itemize{
  \item \code{alloc_rule=1} (Riviere et al 2014):
    If P(toxicity probability at combination (i,j) < \code{target}) > \code{c_e}: among combinations in the neighborhood (-1, +1), (0, +1), (+1, 0), (+1, -1), choose the combination with a higher estimated toxicity probability than the current combination and with the estimated toxicity probability closest to \code{target}.
    If P(toxicity probability at combination (i,j) > \code{target}) > 1-\code{c_d}: among neighborhood (-1, +1), (-1, 0), (0, -1), (+1, -1), choose the combination with a lower estimated toxicity probability than the current combination and with the estimated toxicity probability closest to \code{target}.
    Otherwise, remain on the same combination.
  \item \code{alloc_rule=2}: Among combinations already tested and combinations in the neighborhood (-1, 0), (-1, +1), (0, +1), (+1, 0), (+1, -1), (0, -1), (-1, -1) of a combination tested, choose the combination with the highest posterior probability to be in the targeted interval [\code{target_min}, \code{target_max}] while controling overdosing i.e. P(toxicity probability at combination (i,j) > \code{target_max}) < \code{c_over}.
  \item \code{alloc_rule=3}: Among combinations in the neighborhood (-1, 0), (-1, +1), (0, +1), (+1, 0), (+1, -1), (0, -1), (-1, -1) of the current combination, choose the combination with the highest posterior probability to be in the targeted interval [\code{target_min}, \code{target_max}] while controling overdosing i.e. P(toxicity probability at combination (i,j) > \code{target_max}) < \code{c_over}.
}

Early stopping for over-dosing:
  If the current combination is the lowest (1, 1) and at least \code{cmin_overunder} cohorts have been included at that combination and P(toxicity probability at combination (i,j) > \code{target}) >= \code{c_stop} then stop the trial and do not recommend any combination.

Early stopping for under-dosing:
  If the current combination is the highest and at least \code{cmin_overunder} cohorts have been included at that combination and P(toxicity probability at combination (i,j) < \code{target}) >= \code{c_stop} then stop the trial and do not recommend any combination.

Early stopping for identifying the MTD:
\itemize{
  \item \code{early_stop=1} (Riviere et al 2014):
  No stopping rule, include patients until maximum sample size is reached.
  \item \code{early_stop=2}:
  If the next recommended combination has been tested on at least \code{cmin_mtd} cohorts and has a posterior probability to be in the targeted interval [\code{target_min}, \code{target_max}] that is >= \code{c_t} and also control overdosing i.e. P(toxicity probability at current combination > \code{target_max}) < \code{c_over} then stop the trial and recommend this combination.
  \item \code{early_stop=3}:
  If at least \code{cmin_mtd} cohorts have been included at the next recommended combination then stop the trial and recommend this combination.
}

Stopping at the maximum sample size:
  If the maximum sample size is reached and no stopping rule is met, then the recommended combination is the one that was tested on at least \code{cmin_recom} cohorts and with the highest posterior probability to be in the targeted interval [\code{target_min}, \code{target_max}].
}
\value{
  An object of class "CombIncrease_sim" is returned, consisting of the operating characteristics of the design specified. Objects generated by \code{CombIncrease_sim} contain at least the following components:
  \item{rec_dose}{Percentage of combination selection.}
  \item{n_pat_dose}{Mean number of patients at each combination.}
  \item{n_tox_dose}{Mean number of toxicities at each combination.}
  \item{inconc}{Percentage of inclusive trials.}
  \item{early_conc}{Percentage of trials stopping with criterion for finding MTD.}
  \item{nsim}{Number of simulations (if function stopped while executed, return the current number of simulations performed with associated other outputs).}
  \item{pat_tot}{Total mean number of patients accrued.}
  \item{tab_pat}{Vector with the number of patients included for each simulation.}
}
\references{
Riviere, M-K., Yuan, Y., Dubois, F., and Zohar, S. (2014). A Bayesian dose-finding design for drug combination clinical trials based on the logistic model. Pharmaceutical Statistics.
}
\author{Jacques-Henri Jourdan and Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\seealso{
  \code{\link{CombIncrease_next}}.
}
\examples{
\donttest{
p_tox_sc1 = matrix(c(0.05,0.10,0.15,0.30,0.45,
                     0.10,0.15,0.30,0.45,0.55,
                     0.15,0.30,0.45,0.50,0.60),nrow=5,ncol=3)
prior_a1 = c(0.12, 0.2, 0.3, 0.4, 0.5)
prior_a2 = c(0.2, 0.3, 0.4)

sim1 = CombIncrease_sim(ndose_a1=5, ndose_a2=3, p_tox=p_tox_sc1, target=0.30,
  target_min=0.20, target_max=0.40, prior_tox_a1=prior_tox_a1,
  prior_tox_a2=prior_tox_a2, n_cohort=20, cohort=3, tite=FALSE, nsim=2000,
  c_over=1, cmin_overunder=3, cmin_recom=1, startup=1, alloc_rule=1,
  early_stop=1, seed=14061991)
sim1
}

# Dummy example, running quickly
useless = CombIncrease_sim(ndose_a1=3, ndose_a2=2,
  p_tox=matrix(c(0.05,0.15,0.30,0.15,0.30,0.45),nrow=3), target=0.30,
  target_min=0.20, target_max=0.40, prior_tox_a1=c(0.2,0.3,0.4),
  prior_tox_a2=c(0.2,0.3), n_cohort=2, cohort=2, nsim=1)
}
