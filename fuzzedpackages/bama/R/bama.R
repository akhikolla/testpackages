#' Bayesian Mediation Analysis
#'
#' `bama` is a Bayesian inference method that uses continuous shrinkage priors
#' for high-dimensional Bayesian mediation analysis, developed by Song et al
#' (2019). \code{bama} provides estimates for the regression coefficients as
#' well as the posterior inclusion probability for ranking mediators.
#'
#' \code{bama} uses two regression models for the two conditional relationships,
#' \eqn{Y | A, M, C} and \eqn{M | A, C}. For the outcome model, \code{bama}
#' uses
#' \deqn{Y = M \beta_M  + A * \beta_A + C* \beta_C + \epsilon_Y}
#' For the mediator model, \code{bama} uses the model
#' \deqn{M = A * \alpha_A + C * \alpha_C + \epsilon_M}
#'
#' For high dimensional tractability, \code{bama} employs continuous Bayesian
#' shrinkage priors to select mediators and makes the two following assumptions:
#' First, it assumes that all the potential mediators contribute small effects
#' in mediating the exposure-outcome relationship. Second, it assumes
#' that only a small proportion of mediators exhibit large effects
#' ("active" mediators). \code{bama} uses a Metropolis-Hastings within Gibbs
#' MCMC to generate posterior samples from the model.
#'
#' @param Y Length \code{n} numeric outcome vector
#' @param A Length \code{n} numeric exposure vector
#' @param M \code{n x p} numeric matrix of mediators of Y and A
#' @param C1 \code{n x nc1} numeric matrix of extra covariates to include in the
#'     outcome model
#' @param C2 \code{n x nc2} numeric matrix of extra covariates to include in the
#'     mediator model
#' @param beta.m Length \code{p} numeric vector of initial \code{beta.m} in the
#'     outcome model
#' @param alpha.a Length \code{p} numeric vector of initial \code{alpha.a} in
#'     the mediator model
#' @param burnin number of iterations to run the MCMC before sampling
#' @param ndraws number of draws to take from MCMC after the burnin period
#' @param weights Length \code{n} numeric vector of weights
#' @param k Shape parameter prior for inverse gamma
#' @param lm0 Scale parameter prior for inverse gamma for the small normal
#'    components
#' @param lm1 Scale parameter prior for inverse gamma for the large normal
#'    components
#' @param l Scale parameter prior for the other inverse gamma distributions
#' @return
#' \code{bama} returns a object of type "bama" with 12 elements:
#' \describe{
#' \item{beta.m}{\code{ndraws x p} matrix containing outcome model mediator
#'       coefficients.
#' }
#' \item{r1}{\code{ndraws x p} matrix indicating whether or not each beta.m
#'     belongs to the larger normal component (1) or smaller normal
#'     component (0).
#' }
#' \item{alpha.a}{\code{ndraws x p} matrix containing the mediator model
#'     exposure coefficients.
#' }
#' \item{r3}{\code{ndraws x p} matrix indicating whether or not each alpha.a
#'     belongs to the larger normal component (1) or smaller normal component (0).
#' }
#' \item{beta.a}{Vector of length \code{ndraws} containing the beta.a coefficient.}
#' \item{pi.m}{Vector of length \code{ndraws} containing the proportion of
#'     non zero beta.m coefficients.
#' }
#' \item{pi.a}{Vector of length \code{ndraws} containing the proportion of
#'     non zero alpha.a coefficients.
#' }
#'   \item{sigma.m0}{Vector of length \code{ndraws} containing the standard
#'       deviation of the smaller normal component for mediator-outcome
#'       coefficients (beta.m).
#' }
#' \item{sigma.m1}{Vector of length \code{ndraws} containing standard deviation
#'     of the larger normal component for mediator-outcome coefficients (beta.m).
#' }
#' \item{sigma.ma0}{Vector of length \code{ndraws} containing standard
#'     deviation of the smaller normal component for exposure-mediator
#'     coefficients (alpha.a).
#' }
#' \item{sigma.ma1}{Vector of length \code{ndraws} containing standard deviation
#'     of the larger normal component for exposure-mediator coefficients
#'     (alpha.a).
#' }
#' \item{call}{The R call that generated the output.}
#' }
#' @examples
#' library(bama)
#'
#' Y <- bama.data$y
#' A <- bama.data$a
#'
#' # grab the mediators from the example data.frame
#' M <- as.matrix(bama.data[, paste0("m", 1:100)], nrow(bama.data))
#'
#' # We just include the intercept term in this example as we have no covariates
#' C1 <- matrix(1, 1000, 1)
#' C2 <- matrix(1, 1000, 1)
#' beta.m  <- rep(0, 100)
#' alpha.a <- rep(0, 100)
#'
#' set.seed(12345)
#' out <- bama(Y, A, M, C1, C2, beta.m, alpha.a, burnin = 1000, ndraws = 100)
#'
#' # The package includes a function to summarise output from 'bama'
#' summary <- summary(out)
#' head(summary)
#' @references
#' Song, Y, Zhou, X, Zhang, M, et al. Bayesian shrinkage estimation of high
#' dimensional causal mediation effects in omics studies. Biometrics. 2019;
#' 1-11. \href{http://doi.org/10.1111/biom.13189}{doi:10.1111/biom.13189}
#' @author Alexander Rix
#' @export
bama <- function(Y, A, M, C1, C2, beta.m, alpha.a, burnin, ndraws,
                 weights = NULL, k = 2.0, lm0 = 1e-4, lm1 = 1.0, l = 1.0)
{
    call <- match.call()

    if (!is.vector(Y) || !is.numeric(Y))
        stop("'Y' must be a numeric vector.")
    if (any(is.na(Y)))
        stop("'Y' must not have missing values.")

    n <- length(Y)

    if (!is.vector(A) || !is.numeric(A))
        stop("'A' should be a numeric vector.")
    if (any(is.na(A)))
        stop("'A' cannot have missing values.")
    if (length(A) != n)
        stop("Lengths of 'A' and 'Y' do not match.")

    if (!is.matrix(M) || !is.numeric(M))
        stop("'M' must be a numeric matrix.")
    if (any(is.na(M)))
        stop("'M' cannot have missing values.")
    if (nrow(M) != length(Y))
        stop("The number of rows in 'M' does not match the length of 'Y'.")

    if (!is.matrix(C1) || !is.numeric(C1))
        stop("'C1' must be a numeric matrix.")
    if (any(is.na(C1)))
        stop("'C1' cannot have missing values.")
    if (nrow(C1) != length(Y))
        stop("The number of rows in 'C1' does not match the length of 'Y'.")
    if (!is.matrix(C2) || !is.numeric(C2))
        stop("'C2' must be a numeric matrix.")

    if (any(is.na(C2)))
        stop("'C2' cannot have missing values.")
    if (nrow(C2) != length(Y))
        stop("The number of rows in 'C2' does not match the length of 'Y'.")

    if (!is.vector(beta.m) || !is.numeric(beta.m))
        stop("'beta.m' must be a numeric vector")
    if (any(is.na(beta.m)))
        stop("'beta.m' cannot contain missing values.")

    if (!is.vector(alpha.a) || !is.numeric(alpha.a))
        stop("'alpha.a' must be a numeric vector")
    if (any(is.na(alpha.a)))
        stop("'alpha.a' cannot contain missing values.")

    if (!is.numeric(burnin) || !is.vector(burnin) || length(burnin) != 1)
        stop("'burnin' should be a nonnegative integer.")

    if (!is.numeric(ndraws) || !is.vector(ndraws) || length(ndraws) != 1)
        stop("'ndraws' should be a nonnegative integer.")

    if (!is.numeric(k) || !is.vector(k) || length(k) != 1 || k < 0)
        stop("'k' should be a nonnegative number.")

    if (!is.numeric(lm0) || !is.vector(lm0) || length(lm0) != 1 || lm0 < 0)
        stop("'lm0' should be a nonnegative number.")

    if (!is.numeric(lm1) || !is.vector(lm1) || length(lm1) != 1 || lm1 < 0)
        stop("'lm1' should be a nonnegative number.")

    if (!is.numeric(l) || !is.vector(l) || length(l) != 1 || l < 0)
        stop("'l' should be a nonnegative number.")


    if (!is.null(weights)) {
        if (!is.numeric(weights) || !is.vector(weights) ||
            length(weights) != n || any(weights < 0))
        {
            stop("'weights' must be a length 'n' nonnegative numeric vector.")
        }

        w  <- sqrt(weights)

        Y  <- w * Y
        A  <- w * A
        M  <- apply(M, 2, function(m) m * w)
        C1 <- apply(C1, 2, function(c1) c1 * w)
        C2 <- apply(C2, 2, function(c2) c2 * w)
    }

    bama.out <- run_bama_mcmc(Y, A, M, C1, C2, beta.m, alpha.a, burnin, ndraws,
                              k, lm0, lm1, l)

    colnames(bama.out$beta.m)  <- colnames(M)
    colnames(bama.out$alpha.a) <- colnames(M)
    colnames(bama.out$r1)      <- colnames(M)
    colnames(bama.out$r3)      <- colnames(M)

    bama.out$call <- call

    structure(bama.out, class = "bama")
}

#' Summarize objects of type "bama"
#'
#' summary.bama summarizes the 'beta.m' estimates from \code{bama} and generates
#' an overall estimate, credible interval, and posterior inclusion probability.
#' @return A data.frame with 4 elements. The beta.m estimates, the estimates'
#'     *credible* interval (which by default is 95\%), and the posterior
#'      inclusion probability (pip) of each 'beta.m'.
#' @param object An object of class "bama".
#' @param rank Whether or not to rank the output by posterior inclusion
#'     probability. Default is TRUE.
#' @param ci The credible interval to calculate. \code{ci} should be a length 2
#'     numeric vector specifying the upper and lower bounds of the CI. By
#'     default, \code{ci = c(0.025, .975)}.
#' @param ... Additional optional arguments to \code{summary}
#' @export
summary.bama <- function(object, rank = F, ci = c(0.025, .975), ...)
{
    if (class(object) != "bama")
        stop("'object' is not an bama object.")

    if (!is.logical(rank) || length(rank) != 1)
        stop("'rank' should be a length 1 logical.")

    if (!is.numeric(ci) || length(ci) != 2)
        stop("'ci' should be a length 2 numeric.")
    if (ci[1] >= ci[2])
        stop("'ci[1]' should be less than 'ci[2]'.")

    pip    <- colMeans(object$r1 * object$r3)
    beta.m <- colMeans(object$beta.m)

    credible.l <- apply(object$beta.m, 2, stats::quantile, probs = ci[1])
    credible.h <- apply(object$beta.m, 2, stats::quantile, probs = ci[2])

    out <- data.frame(estimate = beta.m, ci.lower = credible.l,
                          ci.upper = credible.h, pip = pip)
    if (rank)
        out <- out[order(pip, decreasing = T), ]

    out
}

#' Printing bama objects
#'
#' Print a bama object.
#' @param x An object of class 'bama'.
#' @param ... Additional arguments to pass to print.data.frame or summary.bama
#' @export
print.bama <- function(x , ...)
{
    print(summary(x, ...), ...)
}
