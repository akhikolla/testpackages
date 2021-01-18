#' @title
#' Function for Rounding Weights
#'
#' @description
#' This function performs an optimal rounding of the provided real weights, 
#' in order to reduce a specific objective function 
#' 
#' @param weights A numerical vector of real weights to be rounded 
#' @param formula A formula to express a linear system for hitting the \code{targets}
#' @param targets A numerical vector of point-targets to hit
#' @param objective A character specifying the objective function used for calibration. By default, it is \code{"L1"}. See details for more information
#' @param tgtBnds A two-column matrix containing the bounds for the point-targets
#' @param lower A numerical vector or value defining the lower bounds of the weights
#' @param upper A numerical vector or value defining the upper bounds of the weights
#' @param sparse A logical value denoting if the linear system is sparse or not. By default, it is \code{FALSE}
#' @param scale A numerical vector of positive values 
#' @param data A \code{data.frame} or \code{matrix} object containing the data to be used for calibration
#'
#' @details 
#' The optimal rounding can be performed by considering one of the following objective functions:
#' \describe{
#'   \item{\code{"L1"}}{for the summation of absolute errors}
#'   \item{\code{"aL1"}}{for the asymmetric summation of absolute errors}
#'   \item{\code{"rL1"}}{for the summation of absolute relative errors}
#'   \item{\code{"LB1"}}{for the summation of absolute errors if outside the boundaries}
#'   \item{\code{"rB1"}}{for the summation of absolute relative errors if outside the boundaries}
#'   \item{\code{"rbLasso1"}}{for the summation of absolute relative errors if outside the boundaries plus a Lasso penalty based on the distance from the provided weights}
#'   \item{\code{"L2"}}{for the summation of square errors}
#'   \item{\code{"aL2"}}{for the asymmetric summation of square errors}
#'   \item{\code{"rL2"}}{for the summation of square relative errors}
#'   \item{\code{"LB2"}}{for the summation of square errors if outside the boundaries}
#'   \item{\code{"rB2"}}{for the summation of square relative errors if outside the boundaries}
#'   \item{\code{"rbLasso2"}}{for the summation of square relative errors if outside the boundaries plus a Lasso penalty based on the distance from the provided weights}
#' }
#' 
#' A two-column matrix must be provided to \code{tgtBnds} when \code{objective = "LB1"}, \code{objective = "rB1"}, 
#' \code{objective = "rbLasso1"}, \code{objective = "LB2"}, \code{objective = "rB2"}, and \code{objective = "rbLasso2"}.
#' 
#' The argument \code{scale} must be specified with a vector of positive reals number when \code{objective = "rL1"} 
#' or \code{objective = "rL2"}.
#' 
#' @return 
#' A vector of integer weights to be the input of the calibration algorithm
#' 
#' @examples
#' library(inca)
#' set.seed(0)
#' w <- rpois(150, 4)
#' data <- matrix(rbinom(150000, 1, .3) * rpois(150000, 4), 1000, 150)
#' y <- data %*% w
#' w <- runif(150, 0, 7.5)
#' rw <- roundWeights(w, ~. + 0, y, lower = 1, upper = 7, sparse = TRUE, data = data)
#' barplot(table(rw), main = "Rounded weigths")
#'
#' @export

"roundWeights" <- function(weights, formula, targets, objective = c("L1", "aL1", "rL1", "LB1", "rB1", "rbLasso1", 
                                                                    "L2", "aL2", "rL2", "LB2", "rB2", "rbLasso2"), 
                           tgtBnds = NULL, lower = -Inf, upper = Inf, scale = NULL, sparse = FALSE,
                           data = environment(formula)) {
  if (is.null(tgtBnds)) {
    if (objective[1] %in% c("LB1", "rB1", "rbLasso1", "LB2", "rB2", "rbLasso2"))
      stop("\"tgtBnds\" must be specified when \"objective\" is either ", paste("\"LB1\",", "\"rB1\",", "\"rbLasso1\",", "\"LB2\",", "\"rB2\",", "or \"rbLasso2\""))
    tgtBnds <- cbind(rep_len(-Inf, length(targets)),
                     rep_len(Inf, length(targets)))
  }
  else { 
    tgtBnds <- as.matrix(tgtBnds)
    if (ncol(tgtBnds) != 2) 
      stop("\"tgtBnds\" must be a data.frame or a matrix object with two columns")
  }
  if(is.null(scale)) {
    if (objective[1] %in% c("rL1", "rL2"))
      warning("The argument \"scale\" should be specified when \"objective\" is either ", paste("\"rL1\"", "or \"rL2\"\n"), 
              "In this case, \"scale\" will be set as a vector of values equal to one by default")
    scale <- rep_len(1, length(targets))
  }
  mylower <- pmax(ceiling(lower), -abs(weights) * Inf)
  myupper <- pmin(floor(upper), abs(weights) * Inf)
  wts <- adjWeights(weights, lower = mylower, upper = myupper)

  A <- model.matrix(formula, data = as.data.frame(data))
  if(sparse) {
    A <- Matrix(A, sparse = TRUE)
    wts <- .Call('sparse_round', PACKAGE = 'inca', A, targets, wts, tgtBnds, scale, objective[1])
  }
  else {
    wts <- .Call('dense_round', PACKAGE = 'inca', A, targets, wts, tgtBnds, scale, objective[1])
  }
  return(wts)
}
