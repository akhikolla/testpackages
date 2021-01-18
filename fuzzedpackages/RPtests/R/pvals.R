#' Compute p-values for \code{RPtest} output
#'
#' Produces p-values given a list of simulated test statistics and the true test
#' statistic (which may be a vector if it is the output of multiple RP
#' functions).
#'
#' @param test_sim A list of test statisitcs, each component of which is a
#'   numeric vector.
#' @param test A numeric vector of test statistics.
#' @details In the case where the individual test statistics are vectors, the
#'   first component of test is compared against the first components of
#'   \code{test_sim[[1]]}, \code{test_sim[[2]]} etc. and the results of these
#'   multiple comparisons are combined into a single p-value (see the
#'   reference). When the lengths of the test statistics differ, the final
#'   components are first discarded to make all the test statistics have equal
#'   length.
#' @return A single p-value.
#' @references Shah, R. D., Buhlmann, P. (2016) \emph{Goodness of fit tests for
#'   high-dimensional linear models} \url{http://arxiv.org/abs/1511.03334}
#' @seealso \code{\link{RPtest}} the output of which this would typically be
#'   applied to.
#' @export
pval <- function(test, test_sim) {
  thresh <- 1e-8 # used for exact_zero function (see below)
  B <- length(test_sim)

  n_tests_max <- 1L
  combine_tests <- FALSE

  if (is.list(test_sim)) {
    n_test_sim_vec <- sapply(test_sim, length)
    n_tests_max <- max(n_test_sim_vec, length(test))
    n_tests <- min(n_test_sim_vec, length(test))
  }
  if (n_tests_max > 1L) {
    combine_tests <- TRUE
    # test is a vector of test statistics
    # finds the minimum number of test statistics. We compare from the front
    test_sim <- sapply(test_sim, function(x) x[seq_len(n_tests)])
    # test_sim is a n_tests by B matrix
    test <- test[1:n_tests]

    # Create an n_test by B matrix of means for each test and leaving out each
    # bootstrap sample in the calculation of the mean
    test_sim_sums <- rowSums(test_sim) + test # a vector of length n_tests of test statistic sums
    mean_test_sim_mat <- (test_sim_sums - test_sim) / B

    # Now the same for the variance
    test_sim_sq <- test_sim^2
    test_sim_sums_sq <- rowSums(test_sim_sq) + test^2 # a vector of length n_tests of test statistic sums
    sq_test_sim_mat <- (test_sim_sums_sq - test_sim_sq) / B

    sd_test_sim_mat <- sqrt(pmax(sq_test_sim_mat - mean_test_sim_mat^2, 0))

    test_sim_mean <- (test_sim_sums - test)/B
    test_sim_sd <- sqrt(pmax((test_sim_sums_sq - test^2) / B - test_sim_mean^2, 0))
    diff <- exact_zero(test-test_sim_mean, thresh)
    # we could have test-test_sim_mean close to machine epsilon but test_sim_sd exactly zero.
    # We would like them to both be zero and this to be discounted when taking min below
    test <- min(diff/test_sim_sd, na.rm=TRUE)
    # no NAs should be produced here It could be -Inf though

    diff <- exact_zero(test_sim - mean_test_sim_mat, thresh)
    test_sim_orig <- test_sim
    test_sim <- -colMax(-(diff) / sd_test_sim_mat, na_rm = TRUE)
    # test_sim <- as.numeric(apply(test_sim, 2, function(x) min((x-mean_test_sim)/sd_test_sim, na.rm=TRUE)))
  } else {
    if (is.list(test_sim)) test_sim <- unlist(test_sim)
  }
  return((sum(test_sim <= test)+1)/(B+1))
}


####
exact_zero <- function(mat, thresh) {
  mat[abs(mat) < thresh] <- 0
  return(mat)
}
