#' Equality tests for two multinomial samples
#'
#' Generate multinomial samples from a common probability vector and calculate the Chi-square and Likelihood Ratio test statistics.
#'
#' @param N1 Size of sample 1.
#' @param N2 Size of sample 2.
#' @param p0 Common probability vector from which to draw the multinomial samples.  Can also be a matrix, in which case each simulation randomly draws with replacement from the rows of p0.
#' @param nreps Number of replications of the simulation.
#' @param verbose Logical.  If \code{TRUE} prints message every \code{5000} replications.
#' @details The chi-squared and likelihood ratio test statistics are calculated from multinomial samples \eqn{(Y_1^1, Y_2^1), \ldots, (Y_1^M, Y_2^M)}{(Y_11, Y_21),\ldots,(Y_1M, Y_2M)}, where
#' \deqn{
#'   Y_k^m \stackrel{\textrm{ind}}{\sim} \textrm{Multinomial}(N_k, p_0^m),
#' }{
#'   Y_km ~ind Multinomial(N_k, p_m),
#' }
#' where \eqn{p_0^m}{p_m} is the \eqn{m}th row of \code{p0}.
#' @return An \code{nreps x 2} matrix with the simulated chi-squared and LR values.
#' @examples
#' # bootstrapped p-value calculation against equal genotype proportions
#' # in lakes Michipicoten and Simcoe
#'
#' # contingency table
#' popId <- c("Michipicoten", "Simcoe")
#' ctab <- UM.suff(fish215[fish215$Lake %in% popId,])$tab
#' ctab
#'
#' # MLE of probability vector
#' p.MLE <- colSums(ctab)/sum(ctab)
#' # sample sizes
#' N1 <- sum(ctab[1,]) # Michipicoten
#' N2 <- sum(ctab[2,]) # Simcoe
#'
#' # bootstrapped test statistics (chi^2 and LRT)
#' T.boot <- UM.eqtest(N1 = N1, N2 = N2, p0 = p.MLE, nreps = 1e3)
#'
#' # observed test statistics
#' T.obs <- c(chi2 = chi2.stat(ctab), LRT = LRT.stat(ctab))
#' # p-values
#' rowMeans(t(T.boot) > T.obs)
#' @importFrom stats rmultinom
#' @export
UM.eqtest <- function(N1, N2, p0, nreps, verbose = TRUE) {
  N <- N1+N2
  P0 <- p0
  if(is.matrix(P0)) {
    Np0 <- nrow(P0)
  } else {
    Np0 <- 1
    p0 <- P0
  }
  # output
  boot.out <- matrix(NA, nreps, 2)
  colnames(boot.out) <- c("chi2", "LRT")
  for(ii in 1:nreps) {
    if(verbose && (ii%%5e3 == 0)) message("bootstrap sample ", ii, "/", nreps)
    if(Np0 > 1) {
      p0 <- P0[sample(Np0, 1),]
    }
    tab.sim <- rbind(c(rmultinom(1, size = N1, prob = p0)),
                     c(rmultinom(1, size = N2, prob = p0)))
    boot.out[ii,1] <- chi2.stat(tab.sim)
    boot.out[ii,2] <- LRT.stat(tab.sim)
  }
  boot.out
}
