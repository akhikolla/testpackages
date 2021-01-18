#' strategies.RPS
#'
#' Six rock-paper-scissors strategies.
#'
#' \describe{
#'  The rock-paper-scissors strategies are:
#'  \item{rock}{Strategy which plays rock.}
#'  \item{paper}{Strategy which plays paper.}
#'  \item{scissors}{Strategy which plays scissors.}
#'  \item{nash}{Nash equilibrium strategy which plays every action with probability one-third.}
#'  \item{mixed}{Strategy which plays every action with a certain probability.}
#'  \item{imitate}{Strategy which plays a random action in the first round and subsequently imitates the last choice.}
#' }
#'
#' @format Each strategy is encoded as a stratEst.strategy object. The rows of the data frame represent the states of the strategy. The first row is the start state of the strategy. Each stratEst.strategy object contains the following variables:
#' \describe{
#'   \item{\code{prob.r}}{Probability to play rock.}
#'   \item{\code{prob.p}}{Probability to play paper.}
#'   \item{\code{prob.s}}{Probability to play scissors.}
#'   \item{\code{tremble}}{Probability of a tremble.}
#'   \item{\code{tr(r)}}{State transition for the input last choice was rock.}
#'   \item{\code{tr(p)}}{State transition for the input last choice was paper.}
#'   \item{\code{tr(s)}}{State transition for the input last choice was scissors.}
#' }
#' @usage data(strategies.RPS)
#' @examples
#' strategies <- strategies.RPS[c("nash","mixed","imitate")]
"strategies.RPS"
