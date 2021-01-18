#' strategies.DFS2020
#'
#' The conformist and anticonformist strategy identifies by Dvorak, Fischbacher, and Schmelz (2020).
#'
#' \describe{
#'  The strategies are:
#'  \item{conformist}{Strategy that generally follows the own preference if the choices of the other group members are in line with the own preference and deviates from the own preference the choices of the other group members are not in line with the own preference.}
#'  \item{anticonformist}{Strategy that frequently deviates from the own preference the choices of the other group members are in line with the own preference and follows the own preference if the choices of the other group members are not in line with the own preference.}
#' }
#'
#' @format Each strategy is encoded as a stratEst.strategy object. The rows of the data frame represent the states of the strategy. The first row is the start state of the strategy. Each stratEst.strategy object contains the following variables:
#' \describe{
#'   \item{\code{prob.follow}}{Probability to follow own preference.}
#'   \item{\code{prob.deviate}}{Probability to deviate from the own preference.}
#'   \item{\code{tr(not in line)}}{State transition for the input the choices of the others are not in line with the own preference.}
#'   \item{\code{tr(in line)}}{State transition for the input the choices of the others are in line with the own preference.}
#' }
#' @usage data(strategies.DFS2020)
#' @references
#' Dvorak F, Fischbacher U, Schmelz K (2020). "Incentives for Conformity and Anticonformity." \emph{TWI Working Paper Series}.
#'
#' @examples
#' strategies <- strategies.DFS2020[c("conformist","anticonformist")]
"strategies.DFS2020"
