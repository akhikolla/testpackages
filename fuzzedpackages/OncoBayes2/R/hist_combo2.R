#' Dataset: historical data on two single-agents to inform a combination study
#'
#' One of two datasets from the application described in Neuenschwander et al
#' (2016). The risk of DLT is to be studied as a function of dose for two drugs,
#' drug A and drug B. Historical information on the toxicity profiles of these
#' two drugs is available from single agent trials \code{trial_A} and \code{trial_B}.
#' A second dataset \code{codata_combo2} is available from this application,
#' which includes additional dose-toxicity data from \code{trial_AB} and \code{IIT} of the
#' combination of Drugs A and B.
#'
#' @format A tibble with 11 rows and 5 variables:
#' \describe{
#'   \item{group_id}{study}
#'   \item{drug_A}{dose of Drug A}
#'   \item{drug_B}{dose of Drug B}
#'   \item{num_patients}{number of patients}
#'   \item{num_toxicities}{number of DLTs}
#'   \item{cohort_time}{cohort number of patients}
#' }
#' @template ref-mac
#'
"hist_combo2"
