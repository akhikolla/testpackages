#' Hierarchical Integrative Group LASSO
#'
#' HiGLASSO is a regularization based selection method designed to detect
#' non-linear interactions between variables, particularly exposures in
#' environmental health studies.
#'
#'
#' There are a few things to keep in mind when using \code{higlasso}
#' \itemize{
#'     \item{\code{higlasso} uses the strong heredity principle. That is,
#'           \code{X_1} and \code{X_2} must included as main effects before the
#'           interaction \code{X_1 X_2} can be included.
#' }
#'     \item{While \code{higlasso} uses integrative weights to help with
#'           estimation, \code{higlasso} is more of a selection method.
#'           As a result, \code{higlasso} does not output coefficient estimates,
#'           only which variables are selected.
#' }
#'     \item{Simulation studies suggest that `higlasso` is a very
#'           conservative method when it comes to selecting interactions.
#'           That is, \code{higlasso} has a low false positive rate and the
#'           identification of a nonlinear interaction is a good indicator that
#'           further investigation is worthwhile.
#' }
#' \item{\code{higlasso} can be slow, so it may may be beneficial to
#'       tweak some of its settings (for example, \code{nlambda1} and
#'        \code{nlambda2}) to get a handle on how long the method will take
#'        before running the full model.
#' }}
#' @param Y A length n numeric response vector
#' @param X A n x p numeric matrix of covariates to basis expand
#' @param Z A n x m numeric matrix of non basis expanded and non
#'     regularized covariates
#' @param method Type of initialization to use. Possible choices are \code{gglasso}
#'     for group LASSO and \code{aenet} for adaptive elastic net. Default is
#'     \code{aenet}
#' @param lambda1 A numeric vector of main effect penalties on which to tune
#'     By default, \code{lambda1 = NULL} and higlasso generates a length
#'     \code{nlambda1} sequence of lambda1s based off of the data and
#'     \code{min.lambda.ratio}
#' @param lambda2 A numeric vector of interaction effects penalties on which to
#'     tune. By default, \code{lambda2 = NULL} and generates a sequence (length
#'     \code{nlambda2}) of lambda2s based off of the data and
#'     \code{min.lambda.ratio}
#' @param nlambda1 The number of lambda1 values to generate. Default is 10,
#'     minimum is 2. If \code{lambda1 != NULL}, this parameter is ignored
#' @param nlambda2 The number of lambda2 values to generate. Default is 10,
#'     minimum is 2. If \code{lambda2 != NULL}, this parameter is ignored
#' @param lambda.min.ratio Ratio that calculates min lambda from max lambda.
#'     Ignored if 'lambda1' or 'lambda2' is non NULL. Default is 0.05
#' @param sigma Scale parameter for integrative weights. Technically a third
#'     tuning parameter but defaults to 1 for computational tractability
#' @param degree Degree of \code{bs} basis expansion. Default is 2
#' @param maxit Maximum number of iterations. Default is 5000
#' @param tol Tolerance for convergence. Default is 1e-5
#' @author Alexander Rix
#' @references
#' A Hierarchical Integrative Group LASSO (HiGLASSO) Framework for Analyzing
#' Environmental Mixtures. Jonathan Boss, Alexander Rix, Yin-Hsiu Chen, Naveen N.
#' Narisetty, Zhenke Wu, Kelly K. Ferguson, Thomas F. McElrath, John D. Meeker,
#' Bhramar Mukherjee. 2020.
#' arXiv:2003.12844
#' @return
#' An object of type "higlasso" with 4 elements:
#' \describe{
#'     \item{lambda}{An \code{nlambda1 x nlambda2 x 2} array containing each
#'                   pair \code{(lambda1, lambda2)} pair.}
#'     \item{selected}{An \code{nlambda1 x nlambda2 x ncol(X)} array containing
#'                     higlasso's selections for each lambda pair.}
#'     \item{df}{The number of nonzero selections for each lambda pair.}
#'     \item{call}{The call that generated the output.}
#' }
#' @examples
#' library(higlasso)
#'
#' X <- as.matrix(higlasso.df[, paste0("V", 1:7)])
#' Y <- higlasso.df$Y
#' Z <- matrix(1, nrow(X))
#'
#' \donttest{
#' # This can take a bit of time
#' higlasso.fit <- higlasso(Y, X, Z)
#' }
#' @export
higlasso <- function(Y, X, Z, method = c("aenet", "gglasso"), lambda1 = NULL,
                     lambda2 = NULL, nlambda1 = 10, nlambda2 = 10,
                     lambda.min.ratio = .05, sigma = 1, degree = 2,
                     maxit = 5000, tol = 1e-5)
{
    call <- match.call()

    check.Y(Y)
    check.XZ(X, Y)
    check.XZ(Z, Y)

    method <- match.arg(method)
    if (!is.numeric(sigma) || sigma < 0)
        stop("'sigma' must be a nonnegative number.")
    if (length(sigma) > 1)
        stop("'sigma' must have unit length.")

    if (!is.numeric(degree) || degree < 1)
        stop("'degree' should be an integer >= 1.")

    if (!is.numeric(maxit) || maxit < 1)
        stop("'maxit' should be an integer >= 1.")

    if (!is.numeric(tol) || tol <= 0)
        stop("'tol' should be a positive number.")

    # get number of main effect variables.

    matrices <- generate_design_matrices(X, degree)

    Xm <- matrices$Xm
    Xi <- matrices$Xi
    X.xp <- matrices$X.xp
    groups <- matrices$groups
    igroups <- matrices$igroups

    # X.xp <- do.call("cbind", c(unname(Xm), Xi[j]))

    # generate lambda sequences if user does not pre-specify them
    p <- ncol(X) * degree
    YtX <- abs(Y %*% X.xp)[1,] / nrow(X.xp)
    if (!is.null(lambda1)) {
        if (!is.numeric(lambda1) || any(lambda1 <= 0))
            stop("'lambda1' must be a nonnegative numeric array.")
        nlambda1 <- length(lambda1)
    } else {
        lambda1.max <- max(YtX[1:p])
        lambda1.min <- lambda.min.ratio * lambda1.max
        lambda1 <- exp(seq(log(lambda1.max), log(lambda1.min), length.out =
                           nlambda1))
    }

    if (!is.null(lambda2)) {
        if (!is.numeric(lambda2) || any(lambda2 <= 0))
            stop("'lambda2' must be a nonnegative numeric array.")
        nlambda2 <- length(lambda2)
    } else {
        lambda2.max <- max(YtX[-(1:p)])
        lambda2.min <- lambda.min.ratio * lambda2.max
        lambda2 <- exp(seq(log(lambda2.max), log(lambda2.min), len = nlambda2))
    }

    px       <- ncol(X.xp)
    pz       <- ncol(Z)
    X.xp  <- cbind(X.xp, Z)

    fit <- higlasso.fit(Y, Xm, Xi, Z, X.xp, px, pz, method, lambda1, lambda2,
                        sigma, groups, igroups, maxit, tol, call)

    fit$coef <- NULL

    return(fit)
}

initialise_gglasso <- function(X.xp, pz, Y, tol, maxit, groups, igroups)
{
    groups <- c(groups, seq(pz) + max(groups))
    fit <- gglasso::cv.gglasso(X.xp, Y, group = groups)
    i   <- which.min(fit$cvm)
    purrr::map(igroups, ~ fit$gglasso.fit$beta[.x, i])
}


initialise_aenet <- function(X.xp, px, pz, Y, lambda2, tol, maxit, groups,
                             igroups)
{
    # penalty factor for enet. Z contains unregularized coefficents so we set
    # those weights to 0
    pf <- c(rep(1, px), rep(0, pz))

    enet <- gcdnet::cv.gcdnet(X.xp, Y, method = "ls", lambda2 = lambda2,
                              pf = pf, pf2 = pf, eps = tol,
                              maxit = max(maxit, 1e6))

    # get the best scoring lambda from gcdnet and use that to generate
    # inital weights for the adpative elastic net
    i <- which.min(enet$cvm)
    weights <- enet$gcdnet.fit$beta[1:px, i]
    weights <- 1 / abs(weights + 1 / nrow(X.xp)) ^ 2
    weights <- c(weights, rep(0, pz))
    aenet <- gcdnet::cv.gcdnet(X.xp, Y, method = "ls", lambda2 = lambda2,
                               pf = weights, pf2 = pf, eps = tol,
                               maxit = max(maxit, 1e6))

    i <- which.min(aenet$cvm)
    purrr::map(igroups, ~ aenet$gcdnet.fit$beta[.x, i])
}

higlasso.fit <- function(Y, Xm, Xi, Z, X.xp, px, pz, method, lambda1, lambda2,
                         sigma, groups, igroups, maxit, tol, call)
{
    ngroups  <- length(Xm)
    nlambda1 <- length(lambda1)
    nlambda2 <- length(lambda2)

    lambda   <- array(0, c(nlambda1, nlambda2, 2))
    coef     <- array(0, c(nlambda1, nlambda2, ncol(X.xp)))
    selected <- array(0, c(nlambda1, nlambda2, choose(ngroups, 2) + ngroups))

    nm <- names(Xm)
    k   <- purrr::map_lgl(Xi, ~ ncol(.x) > 0)
    nmi <- purrr::map_chr(purrr::cross2(nm, nm),
                          ~ paste(.x[[1]], .x[[2]], sep = "."))

    dimnames(selected) <- list(NULL, NULL, c(nm, nmi[k]))
    df <- matrix(0, nlambda1, nlambda2)


    if (method == "gglasso")
        start <- initialise_gglasso(X.xp, pz, Y, tol, maxit, groups, igroups)

    for (j in seq_along(lambda2)) {
        if (method == "aenet")
            start <- initialise_aenet(X.xp, px, pz, Y, lambda2[j], tol, maxit,
                                      groups, igroups)

        beta  <- start[1:ngroups]
        eta   <- start[-(1:ngroups)]
        for (i in seq_along(lambda1)) {

            fit <- higlasso_internal(Y, Xm, Xi, Z, beta, eta, lambda1[i],
                                     lambda2[j], sigma, maxit, tol)


            lambda[i, j,]   <- c(lambda1[i], lambda2[j])
            coef[i, j, ] <- unlist(c(fit$beta, fit$gamma[k], fit$alpha))
            selected[i, j,] <- purrr::map_lgl(c(fit$beta, fit$gamma[k]),
                                              ~ any(.x != 0))
            df[i, j]        <- sum(selected[i, j,])
        }
    }

    structure(list(lambda = lambda, coef = coef, selected = selected,
                   df = df, call = call), class = "higlasso")
}
