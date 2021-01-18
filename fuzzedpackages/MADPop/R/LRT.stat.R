#' Likelihood ratio test statistic for contingency tables
#'
#' Calculate the likelihood ratio test statistic for general two-way contingency tables.
#'
#' @param tab A \code{K x C} matrix (contingency table) of counts. See details.
#' @return The calculated value of the LRT statistic.
#' @details Suppose that \code{tab} consists of counts from \eqn{K} populations (rows) in \eqn{C} categories.  The likelihood ratio test statistic is computed as
#' \deqn{
#'   2 \sum_{i=1}^K \sum_{j=1}^N O_{ij} \log(p^A_{ij}/p^0_{j}),
#' }{
#'   2 \sum_ij O_ij log(p_ij/p_0j),
#' }
#' where \eqn{O_{ij}}{O_ij} is the observed number of counts in the \eqn{i}th row and \eqn{j}th column of \code{tab}, \eqn{p^A_{ij} = O_{ij}/\sum_{j=1}^C O_{ij}}{p_ij = O_ij/(\sum_j O_ij)} is the unconstrained estimate of the proportion of category \eqn{j} in population \eqn{i}, and \eqn{p^0_j = \sum_{i=1}^K O_{ij} / \sum_{i=1}^K\sum_{j=1}^C O_{ij}}{p_0j = \sum_i O_ij / \sum_ij O_ij} is the estimate of this proportion under \eqn{H_0} that the populations have indentical proportions in each category.  If any column has only zeros it is removed before calculating the LRT statistic.
#' @examples
#' # simple contingency table
#' ctab <- rbind(pop1 = c(5, 3, 0, 3),
#'                 pop2 = c(4, 10, 2, 5))
#' colnames(ctab) <- LETTERS[1:4]
#' ctab
#' LRT.stat(ctab) # likelihood ratio statistic
#' @export
LRT.stat <- function(tab) {
  K <- nrow(tab) # number of lakes
  N <- rowSums(tab) # number of fish in each lake
  p0 <- colSums(tab)/sum(N) # MLE of common probability vector
  lp0 <- matrix(rep(log(p0), each = K), nrow = K) # log term
  L0 <- sum((tab * lp0)[lp0 > -Inf]) # loglik under H0
  LA <- sum((tab * log(tab/N))[tab > 0]) # loglik under HA
  2 * (LA - L0) # LRT (>= 0)
}
