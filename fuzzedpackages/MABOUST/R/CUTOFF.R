#' Returns the superiority or futility cutoff during a MABOUST trial.
#' @param nTreat Number of active treatments in consideration, i.e. 1,...,K.
#' @param nCat Number of ordinal outcome categories, i.e. J.
#' @param Delta Value of \eqn{\Delta} to test.
#' @param gamma Length 3 vector of cutoff parameters.
#' @param n Current sample size in the trial.
#' @return The set of active treatments to continue, an optimal treatment, or a set of equally optimal treatments. Also reports posterior mean utilities and ordinal outcome probabilities as well as pairwise comparisons of utility similarity, when appropriate.
#' @references
#' [1] Chapple and Clement (2020), MABOUST: A Multi-Armed Bayesian Ordinal Outcome Utility-Based Sequential Trial. Submitted.
#' @examples
#' ###Trial parameters
#' nCat = 6
#' nTreat = 3
#' Delta=5
#' n=300
#' ###Design parameters
#' gamma= c(.5, .05, .05)
#' CUTOFF(Delta,n,nTreat,nCat,gamma)
#' @export
CUTOFF = function(Delta,n,nTreat,nCat,gamma){
  cut = gamma[1]+exp(-gamma[2]*Delta -gamma[3]*n/(nTreat*nCat) )
  return(cut)
}

