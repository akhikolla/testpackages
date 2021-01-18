#' Dataset: historical and concurrent data on a three-way combination
#'
#' This dataset involves a hypothetical dose-escalation study of combination
#' therapy with three treatment components. From two previous studies
#' \code{HistAgent1} and \code{HistAgent2}, historical data is available on each
#' of the treatments as single-agents, as well as two of the two-way
#' combinations. However, due to a difference in treatment schedule between the
#' \code{Combo} study and the historical studies, a stratification (through \code{stratum})
#' is made between the groups to allow differential discounting of the
#' alternate-dose data.
#'
#' @format A data frame with 18 rows and 7 variables:
#' \describe{
#'   \item{group_id}{study}
#'   \item{drug_A}{dose of Drug A}
#'   \item{drug_B}{dose of Drug B}
#'   \item{drug_C}{dose of Drug C}
#'   \item{num_patients}{number of patients}
#'   \item{num_toxicities}{number of DLTs}
#'   \item{stratum}{stratum for \code{group_id}'s used for differential discounting}
#' }
#'
#' @template start-example
#' @template example-combo3
#' @template stop-example
#'
"hist_combo3"
