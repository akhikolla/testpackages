#' Bayesian Mediation Analysis Controlling For False Discovery
#'
#' \code{fdr.bama} uses the permutation test to estimate the null PIP
#' distribution for each mediator and determines a threshold (based off of the
#' \code{fdr} parameter) for significance.
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
#' @param burnin Number of iterations to run the MCMC before sampling
#' @param ndraws Number of draws to take from MCMC after the burnin period
#' @param fdr False discovery rate. Default is 0.1
#' @param npermutations The number of permutations to generate while estimating
#'     the null pip distribution. Default is 200
#' @param weights Length \code{n} numeric vector of weights
#' @param k Shape parameter prior for inverse gamma. Default is 2.0
#' @param lm0 Scale parameter prior for inverse gamma for the small normal
#'    components. Default is 1e-4
#' @param lm1 Scale parameter prior for inverse gamma for the large normal
#'    components. Default is 1.0
#' @param l Scale parameter prior for the other inverse gamma distributions.
#'     Default is 1.0
#' @param mc.cores The number of cores to use while running \code{fdr.bama}.
#'     \code{fdr.bama} uses the \code{parallel} package for parallelization,
#'     so see that for more information. Default is 1 core
#' @param type Type of cluster to make when \code{mc.cores > 1}. See
#'     \code{makeCluster} in the \code{parallel} package for more details.
#'     Default is "PSOCK"
#' @return
#' \code{fdr.bama} returns a object of type "fdr.bama" with 5 elements:
#' \describe{
#' \item{bama.out}{Output from the \code{bama} run.}
#' \item{pip.null}{A \code{p x npermutations} matrices containing the
#'     estimated null PIP distribution for each mediator.
#' }
#' \item{threshold}{The cutoff significance threshold for each PIP controlling
#'     for the false discovery rate.
#' }
#' \item{fdr}{The false discovery rate used to calculate \code{threshold}.}
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
#' \donttest{
#' out <- fdr.bama(Y, A, M, C1, C2, beta.m, alpha.a, burnin = 1000,
#'                 ndraws = 100, npermutations = 10)
#'
#' # The package includes a function to summarise output from 'fdr.bama'
#' summary(out)
#' }
#' @references
#' Song, Y, Zhou, X, Zhang, M, et al. Bayesian shrinkage estimation of high
#' dimensional causal mediation effects in omics studies. Biometrics. 2019;
#' 1-11. \href{http://doi.org/10.1111/biom.13189}{doi:10.1111/biom.13189}
#' @author Alexander Rix
#' @export
fdr.bama <- function(Y, A, M, C1, C2, beta.m, alpha.a, burnin, ndraws,
                     weights = NULL, npermutations = 200, fdr = 0.1, k = 2.0,
                     lm0 = 1e-4, lm1 = 1.0, l = 1.0, mc.cores = 1, type = "PSOCK")
{
    call <- match.call()

    if (npermutations <= 0)
        stop("'npermutations' must be a positive integer.")

    if (fdr <= 0 || fdr >=1)
        stop("'fdr' must be in the interval (0, 1).")

    bama.out <- bama(Y, A, M, C1, C2, beta.m, alpha.a, burnin, ndraws, weights,
                     k, lm0, lm1, l)

    pi2 <- colMeans(bama.out$r1 == 0 & bama.out$r3 == 1)
    pi3 <- colMeans(bama.out$r1 == 1 & bama.out$r3 == 0)
    pi4 <- colMeans(bama.out$r1 == 0 & bama.out$r3 == 0)

    n <- length(Y)
    if (!is.null(weights)) {
        print(weights)
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

    seeds <- sample(.Machine$integer.max, npermutations + 1)
    permute.bama <- function(i)
    {
        set.seed(seeds[i])

        n <- length(Y)
        bama.r1 <- run_bama_mcmc(Y[sample(n)], A, M, C1, C2, beta.m, alpha.a,
                                 burnin, ndraws, k, lm0, lm1, l)

        bama.r3 <- run_bama_mcmc(Y, A, M[sample(n), ], C1, C2, beta.m, alpha.a,
                                 burnin, ndraws, k, lm0, lm1, l)

        bama.r1.r3 <- run_bama_mcmc(Y[sample(n)], A, M[sample(n), ], C1, C2,
                                    beta.m, alpha.a, burnin, ndraws, k, lm0,
                                    lm1, l)

        p2 <- colMeans(bama.r1$r1 == 0 & bama.r1$r3 == 1)
        p3 <- colMeans(bama.r3$r1 == 1 & bama.r3$r3 == 0)
        p4 <- colMeans(bama.r1.r3$r1 == 0 & bama.r1.r3$r3 == 0)

        # Adding small number to prevent divsion by 0 in certain cases
        denom <- pi2 + pi3 + pi4 + 1e-9
        p2 * pi2 / denom + p3 * pi3 / denom + p4 * pi4 / denom
    }

    if (mc.cores == 1) {
        pip.null <- sapply(seq(npermutations), permute.bama)
    }
    else {
        cl <- parallel::makeCluster(mc.cores, type = type)
        parallel::clusterExport(cl, list("Y", "A", "M", "C1", "C2", "beta.m",
                                         "alpha.a", "burnin", "ndraws", "k",
                                         "lm0", "lm1", "l", "pi2",  "pi3",
                                         "pi4", "seeds"),
                                envir = environment()
        )

        pip.null <- parallel::parSapply(cl, seq(npermutations), permute.bama)
        parallel::stopCluster(cl)
    }
    set.seed(seeds[npermutations + 1])

    # calculate the null pip threshold
    threshold <- apply(pip.null, 1, stats::quantile, probs = 1 - fdr)
    names(threshold) <- colnames(M)

    structure(list(bama.out = bama.out, pip.null = pip.null,
                   threshold = threshold, fdr = fdr, call = call),
              class = "fdr.bama")
}

#' Summarize objects of type "fdr.bama"
#'
#' \code{summary.fdr.bama} summarizes the \code{beta.m} estimates from
#' \code{fdr.bama} and for each mediator generates an overall estimate,
#' credible interval, posterior inclusion probability (PIP), and PIP threshold
#' for significance controlling for the specified false discovery rate (FDR).
#' @return A data.frame with 4 elements. The beta.m estimates, the estimates'
#'     *credible* interval (which by default is 95\%), and the posterior
#'      inclusion probability (pip) of each 'beta.m'.
#' @param object An object of class "bama".
#' @param rank Whether or not to rank the output by posterior inclusion
#'     probability. Default is TRUE.
#' @param ci The credible interval to calculate. \code{ci} should be a length 2
#'     numeric vector specifying the upper and lower bounds of the CI. By
#'     default, \code{ci = c(0.025, .975)}.
#' @param fdr False discovery rate. By default, it is set to whatever the
#' \code{fdr} of \code{object} is. However, it can be changed to recalculate
#' the PIP cutoff threshold.
#' @param filter Whether or not to filter out mediators with PIP less than the
#'     PIP threshold.
#' @param ... Additional optional arguments to \code{summary}
#' @export
summary.fdr.bama <- function(object, rank = F, ci = c(0.025, 0.975),
                             fdr = object$fdr, filter = T, ...)
{
    if (class(object) != "fdr.bama")
        stop("'object' is not an bama object.")

    if (!is.logical(rank) || length(rank) != 1)
    stop("'rank' should be a length 1 logical.")

    if (!is.numeric(ci) || length(ci) != 2)
        stop("'ci' should be a length 2 numeric.")
    if (ci[1] >= ci[2])
        stop("'ci[1]' should be less than 'ci[2]'.")

    pip    <- colMeans(object$bama.out$r1 * object$bama.out$r3)
    beta.m <- colMeans(object$bama.out$beta.m)

    ci.l <- apply(object$bama.out$beta.m, 2, stats::quantile, probs = ci[1])
    ci.h <- apply(object$bama.out$beta.m, 2, stats::quantile, probs = ci[2])

    out <- data.frame(estimate = beta.m, ci.lower = ci.l,
                           ci.upper = ci.h, pip = pip,
                           pip.threshold = object$threshold)
    if (rank)
        out <- out[order(pip, decreasing = T), ]

    if (filter)
        out <- out[which(out$pip > out$pip.threshold), ]

    out
}

#' Printing bama objects
#'
#' Print a bama object.
#' @param x An object of class 'bama'.
#' @param ... Additional arguments to pass to print.data.frame or summary.bama
#' @export
print.fdr.bama <- function(x , ...)
{
    print(summary(x, ...), ...)
}
