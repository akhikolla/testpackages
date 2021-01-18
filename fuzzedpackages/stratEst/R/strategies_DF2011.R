#' strategies.DF2011
#'
#' List of six prisoner's dilemma strategies (Dal Bo and Frechette 2011).
#'
#' \describe{
#'  The prisoner's dilemma strategies are:
#'  \item{ALLD}{Strategy which always defects.}
#'  \item{ALLC}{Strategy which always cooperates.}
#'  \item{GRIM}{Strategy which cooperates until one player defects, then GRIM defects forever.}
#'  \item{TFT}{Strategy which cooperates unless the partner defected in the last round.}
#'  \item{WSLS}{Strategy which cooperates if both players chose the same action last round, otherwise WSLS defects.Also known as PTFT.}
#'  \item{T2}{Strategy which cooperates until either player defects, then it defects twice and returns to cooperation (regardless of the actions during the punishment phase).}
#' }
#'
#' @format Each strategy is encoded as a data.frame object. The rows of the data frame represent the states of the automaton. The first row is the start state of the automaton. Each data.frame object contains the following variables:
#' \describe{
#'   \item{\code{prob.d}}{Probability to defect.}
#'   \item{\code{prob.c}}{Probability to cooperate.}
#'   \item{\code{tremble}}{Probability of a tremble.}
#'   \item{\code{tr(cc)}}{State transition for the input cc.}
#'   \item{\code{tr(cd)}}{State transition for the input cd.}
#'   \item{\code{tr(dc)}}{State transition for the input dc.}
#'   \item{\code{tr(dd)}}{State transition for the input dd.}
#' }
#' @usage data(strategies.DF2011)
#' @examples
#' strategies <- strategies.DF2011[c("ALLD","ALLC","TFT","GRIM")]
#' @references
#' Dal Bo P, Frechette GR (2011). "The Evolution of Cooperation in Infinitely Repeated Games: Experimental Evidence." \emph{American Economic Review}, 101(1), 411-429.
#'
"strategies.DF2011"
