\name{CombIncrease_next}
\alias{CombIncrease_next}
\alias{print.CombIncrease_next}
\title{Combination determination with logistic model}
\description{
\code{CombIncrease_next} is used to determine the next or recommended combination in a phase I combination clinical trial using the design proposed by Riviere et al. entitled "A Bayesian dose-finding design for drug combination clinical trials based on the logistic model".
}
\usage{
CombIncrease_next(ndose_a1, ndose_a2, target, target_min, target_max,
  prior_tox_a1, prior_tox_a2, cohort, final, pat_incl, dose_adm1,
  dose_adm2, tite=FALSE, toxicity, time_full=0, time_tox=0,
  time_follow=0, c_e=0.85, c_d=0.45, c_stop=0.95, c_t=0.5, c_over=0.25,
  cmin_overunder=2, cmin_mtd=3, cmin_recom=1, early_stop=1, alloc_rule=1,
  nburn=2000, niter=5000)
}
\arguments{
  \item{ndose_a1}{Number of dose levels for agent 1.}
  \item{ndose_a2}{Number of dose levels for agent 2.}
  \item{target}{Toxicity (probability) target.}
  \item{target_min}{Minimum of the targeted toxicity interval.}
  \item{target_max}{Maximum of the targeted toxicity interval.}
  \item{prior_tox_a1}{A vector of initial guesses of toxicity probabilities associated with the doses of agent 1. Must be of length \code{ndose_a1}.}
  \item{prior_tox_a2}{A vector of initial guesses of toxicity probabilities associated with the doses of agent 2. Must be of length \code{ndose_a2}.}
  \item{cohort}{Cohort size.}
  \item{final}{A boolean with value TRUE if the trial is finished and the recommended combination for further phases should be given, or FALSE (default value) if the combination determination is performed for the next cohort of patients.}
  \item{pat_incl}{Current number of patients included.}
  \item{dose_adm1}{A vector indicating the dose levels of agents 1 administered to each patient included in the trial. Must be of length \code{pat_incl}.}
  \item{dose_adm2}{A vector indicating the dose levels of agents 2 administered to each patient included in the trial. Must be of length \code{pat_incl}.}
  \item{tite}{A boolean indicating if the toxicity is considered as a time-to-event outcome (TRUE), or as a binary outcome (default value FALSE).}
  \item{toxicity}{A vector of observed toxicities (DLTs) for each patient included in the trial. Must be of length \code{pat_incl}. This argument is used/required only if tite=FALSE.}
  \item{time_full}{Full follow-up time window. This argument is used only if tite=TRUE.}
  \item{time_tox}{A vector of times-to-toxicity for each patient included in the trial. If no toxicity was observed for a patient, must be filled with +Inf. Must be of length \code{pat_incl}. This argument is used/required only if tite=TRUE.}
  \item{time_follow}{A vector of follow-up times for each patient included in the trial. Must be of length \code{pat_incl}. This argument is used/required only if tite=TRUE.}
  \item{c_e}{Probability threshold for dose-escalation. The default value is set at 0.85.}
  \item{c_d}{Probability threshold for dose-deescalation. The default value is set at 0.45.}
  \item{c_stop}{Probability threshold for early trial termination. The default value is set at 0.95.}
  \item{c_t}{Probability threshold for early trial termination for finding the MTD (see details). The default value is set at 0.5.}
  \item{c_over}{Probability threshold to control over-dosing (see details).}
  \item{cmin_overunder}{Minimum number of cohorts to be included at the lowest/highest combination before possible early trial termination for over-toxicity or under-toxicity (see details). The default value is set at 2.}
  \item{cmin_mtd}{Minimum number of cohorts to be included at the recommended combination before possible early trial termination for finding the MTD (see details). The default value is set at 3.}
  \item{cmin_recom}{Minimum number of cohorts to be included at the recommended combination at the end of the trial. The default value is set at 1.}
  \item{alloc_rule}{Interger (1, 2, or 3) indicating which allocation rule is used (see details). The default value is set at 1.}
  \item{early_stop}{Interger (1, 2, or 3) indicating which early stopping rule is used (see details). The default value is set at 1.}
  \item{nburn}{Number of burn-in for HMC. The default value is set at 2000.}
  \item{niter}{Number of iterations for HMC. The default value is set at 5000.}
}
\details{
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
An object of class "CombIncrease_next" is returned, consisting of determination of the next combination and estimations. Objects generated by \code{CombIncrease_next} contain at least the following components:
  \item{n_pat_comb}{Number of patients per combination.}
  \item{n_tox_comb}{Number of observed toxicities per combination.}
  \item{pi}{Estimated toxicity probabilities (if the start-up ended).}
  \item{ptox_inf}{Estimated probabilities that the toxicity probability is inferior to \code{target} (if the start-up ended).}
  \item{ptox_inf_targ}{Estimated probabilities of underdosing, i.e. to be inferior to \code{target_min} (if the start-up ended).}
  \item{ptox_targ}{Estimated probabilities to be in the targeted interval [\code{target_min},\code{target_max}] (if the start-up ended).}
  \item{ptox_sup_targ}{Estimated probabilities of overdosing, i.e. to be superior to \code{target_max} (if the start-up ended).}
  \item{(cdose1, cdose2)}{NEXT RECOMMENDED COMBINATION.}
  \item{inconc}{Boolean indicating if trial must stop for under/over dosing.}
  \item{early_conc}{Boolean indicating if trial can be stopped earlier for finding the MTD.}
}
\references{
Riviere, M-K., Yuan, Y., Dubois, F., and Zohar, S. (2014). A Bayesian dose-finding design for drug combination clinical trials based on the logistic model. Pharmaceutical Statistics.
}
\author{Jacques-Henri Jourdan and Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\seealso{
  \code{\link{CombIncrease_sim}}.
}
\examples{
prior_a1 = c(0.12, 0.2, 0.3, 0.4, 0.5)
prior_a2 = c(0.2, 0.3, 0.4)
toxicity1 = c(0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1)
dose1 = c(1,1,1,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4)
dose2 = c(1,1,1,2,2,2,3,3,3,2,2,2,1,1,1,1,1,1)
t_tox = c(rep(+Inf,8),2.9,+Inf,4.6,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,5.2)
follow = c(rep(6,15), 4.9, 3.1, 1.3)

next1 = CombIncrease_next(ndose_a1=5, ndose_a2=3, target=0.3,
  target_min=0.2, target_max=0.4, prior_tox_a1=prior_a1,
  prior_tox_a2=prior_a2, cohort=3, final=FALSE, pat_incl=18,
  dose_adm1=dose1, dose_adm2=dose2, toxicity=toxicity1, c_over=1,
  cmin_overunder=3, cmin_recom=1, early_stop=1, alloc_rule=1)
next1

next2 = CombIncrease_next(ndose_a1=5, ndose_a2=3, target=0.3,
  target_min=0.2, target_max=0.4, prior_tox_a1=prior_a1, prior_tox_a2=prior_a2,
  cohort=3, final=FALSE, pat_incl=18, dose_adm1=dose1,
  dose_adm2=dose2, tite=TRUE, time_full=6, time_tox=t_tox,
  time_follow=follow, c_over=1, cmin_overunder=3, cmin_recom=1,
  early_stop=1, alloc_rule=1)
next2
}
