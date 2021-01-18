#' Title
#'
#' @param Q_alpha a list of length M composed of matrices that are the correlation structure of the CAR prior on beta.
#' @param tau2 a vector of length M that contains the CAR prior precision matrices.
#' @param use_spam a boolean that determines if the output matrix is of class "spam" (\code{use_spam = TRUE}) or of class "dgCMatrix" (\code{use_spam = FALSE}; see Matrix package for details).
#'
#' @importFrom spam bdiag.spam
#' @importFrom Matrix bdiag
#'
#' @return A sparse block diagonal matrix representing the precision matrices for all of the resolutions of the random effects.
#'
#' @examples
#' n_dims <- c(4, 8)
#' phi <- c(0.8, 0.9)
#' tau2 <- c(3, 4)
#' Q_alpha <- make_Q_alpha_2d(n_dims, phi)
#' Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2)
#' ## plot the full precision matrix structure
#' spam::display(Q_alpha_tau2)
#'
#' @export
make_Q_alpha_tau2 <- function(Q_alpha, tau2, use_spam = TRUE) {

    if (length(Q_alpha) != length(tau2))
        stop("Q_alpha must be a list of length M and tau2 must be a positive numeric vector of length M.")

    if (class(Q_alpha) != "Q_alpha")
        stop('Q_alpha must by of class "Q_alpha" which is the output of make_Q_alpha_2d()')

    if (!is_positive_numeric(tau2, length(tau2)))
        stop("tau2 must be a positive numeric vector of length M.")

    if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
        stop("use_spam must be either TRUE or FALSE.")
    }

    M <- length(Q_alpha)
    Q_alpha_tau2 <- vector(mode = "list", length = M)
    for (m in 1:M) {
        Q_alpha_tau2[[m]] <- Q_alpha[[m]] * tau2[m]
    }
    if (use_spam) {
        ## use the spam package
        Q_alpha_tau2 <- do.call(bdiag.spam, Q_alpha_tau2)
    } else {
        ## use the Matrix package
        Q_alpha_tau2 <- do.call(bdiag, Q_alpha_tau2)
    }


    return(Q_alpha_tau2)
}
