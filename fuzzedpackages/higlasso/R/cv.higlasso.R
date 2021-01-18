#' Cross Validated Hierarchical Integrative Group LASSO
#'
#' Does k-fold cross-validation for \code{higlasso}, and returns optimal values
#' for \code{lambda1} and \code{lambda2}.
#'
#' There are a few things to keep in mind when using \code{cv.higlasso}
#' \itemize{
#'     \item{\code{higlasso} uses the strong heredity principle. That is,
#'           \code{X_1} and \code{X_2} must included as main effects before the
#'           interaction \code{X_1 X_2} can be included.
#' }
#'     \item{While \code{higlasso} uses integrative weights to help with
#'           estimation, \code{higlasso} is more of a selection method.
#'           As a result, \code{cv.higlasso} does not output coefficient
#'           estimates, only which variables are selected.
#' }
#'     \item{Simulation studies suggest that `higlasso` is a very
#'           conservative method when it comes to selecting interactions.
#'           That is, \code{higlasso} has a low false positive rate and the
#'           identification of a nonlinear interaction is a good indicator that
#'           further investigation is worthwhile.
#' }
#' \item{\code{cv.higlasso} can be slow, so it may may be beneficial to
#'       tweak some of its settings (for example, \code{nlambda1},
#'       \code{nlambda2}, and \code{nfolds}) to get a handle on how long the
#'       method will take before running the full model.
#' }}
#' As a side effect of the conservativeness of the method, we have found that
#' using the 1 standard error rule results in overly sparse models, and that
#' \code{lambda.min} generally performs better.
#' @param Y A length n numeric response vector
#' @param X A n x p numeric matrix
#' @param Z A n x m numeric matrix
#' @param method Type of initialization to use. Possible choices are
#'     \code{gglasso} for group LASSO and \code{aenet} for adaptive elastic net.
#'     Default is \code{aenet}
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
#' @param nfolds Number of folds for cross validation. Default is 10. The
#'     minimum is 3, and while the maximum is the number of observations
#'     (ie leave one out cross validation)
#' @param foldid An optional vector of values between 1 and
#'     \code{max(foldid)} identifying what fold each observation is in. Default
#'     is NULL and \code{cv.higlasso} will automatically generate \code{foldid}
#'     based off of \code{nfolds}
#' @param sigma Scale parameter for integrative weights. Technically a third
#'     tuning parameter but defaults to 1 for computational tractability
#' @param degree Degree of \code{bs} basis expansion. Default is 2
#' @param maxit Maximum number of iterations. Default is 5000
#' @param tol Tolerance for convergence. Defaults to 1e-5
#' @author Alexander Rix
#' @references
#' A Hierarchical Integrative Group LASSO (HiGLASSO) Framework for Analyzing
#' Environmental Mixtures. Jonathan Boss, Alexander Rix, Yin-Hsiu Chen, Naveen N.
#' Narisetty, Zhenke Wu, Kelly K. Ferguson, Thomas F. McElrath, John D. Meeker,
#' Bhramar Mukherjee. 2020.
#' arXiv:2003.12844
#' @return
#' An object of type \code{cv.higlasso} with 7 elements
#' \describe{
#' \item{lambda}{An \code{nlambda1 x nlambda2 x 2} array containing each
#'     pair \code{(lambda1, lambda2)} pair.}
#' \item{lambda.min}{lambda pair with the lowest cross validation error}
#' \item{lambda.1se}{}
#' \item{cvm}{cross validation error at each lambda pair. The error is
#'     calculated from the mean square error.}
#' \item{cvse}{standard error of \code{cvm} at each lambda pair.}
#' \item{higlasso.fit}{higlasso output from fitting the whole data.}
#' \item{call}{The call that generated the output.}
#' }
#' @examples
#' library(higlasso)
#'
#' X <- as.matrix(higlasso.df[, paste0("V", 1:7)])
#' Y <- higlasso.df$Y
#' Z <- matrix(1, nrow(X))
#'
#'
#' # This can take a bit of time
#' \donttest{
#' fit <- cv.higlasso(Y, X, Z)
#'
#' print(fit)
#' }
#' @export
cv.higlasso <- function(Y, X, Z, method = c("aenet", "gglasso"), lambda1 = NULL,
                        lambda2 = NULL, nlambda1 = 10, nlambda2 = 10,
                        lambda.min.ratio = .05, nfolds = 5, foldid = NULL,
                        sigma = 1, degree = 2, maxit = 5000, tol = 1e-5)
{
    call <- match.call()
    method <- match.arg(method)
    fit <- higlasso(Y, X, Z, method, lambda1, lambda2, nlambda1, nlambda2,
                    lambda.min.ratio, sigma, degree, maxit, tol)

    lambda1 <- fit$lambda[, 1, 1]
    lambda2 <- fit$lambda[1, , 2]

    nlambda1 <- length(lambda1)
    nlambda2 <- length(lambda2)

    n <- length(Y)
    if (!is.null(foldid)) {
        if (!is.numeric(foldid) || !is.vector(foldid) || length(foldid) != n)
            stop("'foldid' must be length n numeric vector.")
        nfolds <- max(foldid)
    } else {
      r     <- n %% nfolds
      p     <- (n - r) / nfolds
      foldid <- c(rep(1:nfolds, p), seq(len = r))
      foldid <- sample(foldid, n)
    }

    matrices <- generate_design_matrices(X, degree)

    Xm <- matrices$Xm
    Xi <- matrices$Xi
    X.xp <- matrices$X.xp
    groups <- matrices$groups
    igroups <- matrices$igroups

    px       <- ncol(X.xp)
    pz       <- ncol(Z)
    X.xp  <- cbind(X.xp, Z)

    cvm <- array(0, c(nlambda1, nlambda2, nfolds))
    cvse <- matrix(0, nlambda1, nlambda2)
    for (i in 1:nfolds) {
        Y.train <- Y[foldid != i]
        Z.train <- Z[foldid != i,, drop = F]

        Xm.train <- purrr::map(Xm, ~ .x[foldid != i,, drop = F])
        Xi.train <- purrr::map(Xi, function(Xi)
            if (nrow(Xi) > 0)
                Xi[foldid != i,, drop = F]
            else
                Xi
        )
        X.xp.train <- X.xp[foldid != i,, drop = F]

        X.xp.test <- X.xp[foldid == i,, drop = F]
        Y.test  <- Y[foldid == i]

        cv.fit <- higlasso.fit(Y.train, Xm.train, Xi.train, Z.train, X.xp.train,
                               px, pz, method, lambda1,  lambda2, sigma, groups,
                               igroups, maxit, tol, call)

        # Calculate cv error for fold
        for (j in seq(nlambda2)) {
            for (k in seq(nlambda1)) {
                res <- Y.test - X.xp.test %*% cv.fit$coef[j, k, ]
                cvm[j, k, i] <- mean(res ^ 2)
            }
        }

    }

    cvse <- apply(cvm, c(1, 2), function(x) stats::sd(x) / sqrt(nfolds))
    cvm  <- apply(cvm, c(1, 2), mean)

    i <- which.min(cvm)

    lambda.min <- fit$lambda[c(i, i + nlambda1 * nlambda2)]


    j <- abs(cvm - min(cvm)) < cvse[i]
    # Inf could cause NaN if df = 0
    i <- which.min(fit$df * ifelse(j, 1, 1e9))

    lambda.1se <- fit$lambda[c(i, i + nlambda1 * nlambda2)]
    structure(list(lambda = fit$lambda, lambda.min = lambda.min,
                   lambda.1se = lambda.1se, cvm = cvm, cvse = cvse,
                   higlasso.fit = fit, call = call), class = "cv.higlasso")
}

#' Print CV HiGLASSO Objects
#'
#' \code{print.cv.higlasso} prints a fitted "cv.higlaso" object and returns it
#' invisibly.
#' @param x An object of type "cv.higlasso" to print
#' @param ... Further arguments passed to or from other methods
#' @return The original input, \code{x} (invisibly).
#' @export
print.cv.higlasso <- function(x, ...)
{
    n1 <- dim(x$lambda)[1]
    n2 <- dim(x$lambda)[2]

    dimnames(x$cvm) <- list(paste0("l1.", seq(n1)), paste0("l2.", seq(n2)))
    cat("'cv.higlaso' fit:\n")
    cat("Average cross validation error for each (lambda1, lambda2)\n")
    print(x$cvm)
    cat("Lambda min:\n")
    cat(x$lambda.min, "\n")
    cat("Lambda 1 SE:\n")
    cat(x$lambda.1se, "\n")
    invisible(x)
}
