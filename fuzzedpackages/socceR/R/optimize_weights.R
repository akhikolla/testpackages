#' Optimize weights from list of prediction matrices
#'
#' Computes the optimal weights to obtain the minimal loss function from a list of prediction matrices.
#'
#' @param predictionlist A list of R x T prediction matrices where each column sum to 1 and each row sums to 
#' @param outcome An integer vector listing the 
#' @param FUN The function used for optimizing the predictions. The default is top use rps for the rank probability score. Another option is logloss for log loss.
#' @return Returns a numeric vector containing an optimal vector of weights that sum to 1 and that minimizes the loss function.
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
#' @keywords manip
#' @examples
#'
#' m1 <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, .5, .5, 0, 0, .5, .5), 4)
#' m1 # Prediction where certain on the top ranks
#' m2 <- matrix(c(.5, .5, 0, 0, .5, .5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 4)
#' m2  # Prediction where the groups are okay 
#' m3 <- matrix(c(.5, .5, 0, 0, .5, .5, 0, 0, 0, 0, .5, .5, 0, 0, .5, .5), 4)
#' m3  # Prediction where no clue about anything
#' m4 <- matrix(rep(1/4, 16), 4)
#' 
#' optimize_weights(list(m1, m2, m3, m4), 1:4)
#'
#' @importFrom stats optim
#' @export
optimize_weights <- function(predictionlist, outcome, FUN=trps) {
  # Sanity checks needed:
  # Check equal dimensions of matrices
  # Check match with outcome
  
  # Start by finding their individual RPS scores
  startrps <- sapply(predictionlist, function(mat) { FUN(mat, outcome)} )
  
  # Should be possible to get much faster
  weightedrps <- function(weights) { 
    weights <- exp(weights)/sum(exp(weights))
    FUN(Reduce('+', lapply(1:length(weights), function(i){weights[i]*predictionlist[[i]]})), outcome)
  }
  
  res <- optim(exp(-startrps), weightedrps)
   
  exp(res$par)/sum(exp(res$par))
}