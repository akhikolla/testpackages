#' Calculating the normalized Gini index
#'
#' This function calculates the Gini index of a classification rule outputting probabilities. It is a classical metric in the context of Credit Scoring.
#' It is equal to 2 times the AUC (Area Under ROC Curve) minus 1.
#' @param actual The numeric binary vector of the actual labels observed.
#' @param predicted The vector of the probabilities predicted by the classification rule.
#' @concept gini index
#' @export
#' @author Adrien Ehrhardt
#' @examples
#' normalizedGini(c(1, 1, 1, 0, 0), c(0.7, 0.9, 0.5, 0.6, 0.3))
normalizedGini <- function(actual, predicted) {
  Gini <- function(a, p) {
    if (length(a) != length(p)) stop("Actual and Predicted need to be equal lengths!")
    temp.df <- data.frame(actual = a, pred = p, range = c(1:length(a)))
    temp.df <- temp.df[order(-temp.df$pred, temp.df$range), ]
    population.delta <- 1 / length(a)
    total.losses <- sum(a)
    null.losses <- rep(population.delta, length(a))
    accum.losses <- temp.df$actual / total.losses
    gini.sum <- cumsum(accum.losses - null.losses)
    sum(gini.sum) / length(a)
  }
  Gini(actual, predicted) / Gini(actual, actual)
}
