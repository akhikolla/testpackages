
#-------------------------------------------------------------------------------
# Estimation of the overlapping group lasso for
# least squares and generalized linear models
#
# August 2015
# Jared Huling
#-------------------------------------------------------------------------------



#' Overlapping Group Lasso (OGLasso)
#'
#' @param x input matrix of dimension nobs by nvars. Each row is an observation,
#' each column corresponds to a covariate
#' @param y numeric response vector of length nobs
#' @param delta vector of length equal to the number of observations with values in 1 and 0, where a 1 indicates the
#' observed time is a death and a 0 indicates the observed time is a censoring event
#' @param group A list of length equal to the number of groups containing vectors of integers
#' indicating the variable IDs for each group. For example, \code{group = list(c(1,2), c(2,3), c(3,4,5))} specifies
#' that Group 1 contains variables 1 and 2, Group 2 contains variables 2 and 3, and Group 3 contains
#' variables 3, 4, and 5. Can also be a matrix of 0s and 1s with the number of columns equal to the
#' number of groups and the number of rows equal to the number of variables. A value of 1 in row i and
#' column j indicates that variable i is in group j and 0 indicates that variable i is not in group j.
#' @param fused matrix specifying generalized lasso penalty formulation. Each column corresponds to each variable
#'  and each row corresponds to a new penalty term, ie if row 1 has the first entry of 1 and the second
#'  entry of -1, then the penalty term lambda.fused * |beta_1 - beta_2| will be added. Not available now
#' @param family "gaussian" for least squares problems, "binomial" for binary response
#' @param nlambda The number of lambda values. Default is 100.
#' @param lambda A user-specified sequence of lambda values. Left unspecified, the a sequence of lambda values is
#' automatically computed, ranging uniformly on the log scale over the relevant range of lambda values.
#' @param lambda.min.ratio Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry
#' value (i.e. the smallest value for which all parameter estimates are zero). The default
#' depends on the sample size nobs relative to the number of variables nvars. If
#' \code{nobs > nvars}, the default is 0.0001, close to zero. If \code{nobs < nvars}, the default
#' is 0.01. A very small value of \code{lambda.min.ratio} will lead to a saturated fit in
#' the \code{nobs < nvars} case.
#' @param lambda.fused tuning parameter for fused (generalized) lasso penalty
#' @param alpha currently not used. Will be used later for fused lasso
#' @param penalty.factor vector of weights to be multiplied to the tuning parameter for the
#' group lasso penalty. A vector of length equal to the number of groups
#' @param penalty.factor.fused vector of weights to be multiplied to the tuning parameter for the
#' fused lasso penalty. A vector of length equal to the number of variables. mostly for internal usage
#' @param group.weights A vector of values representing multiplicative factors by which each group's penalty is to
#' be multiplied. Often, this is a function (such as the square root) of the number of predictors in each group.
#' The default is to use the square root of group size for the group selection methods.
#' @param adaptive.lasso Flag indicating whether or not to use adaptive lasso weights. If set to \code{TRUE} and
#' \code{group.weights} is unspecified, then this will override \code{group.weights}. If a vector is supplied to group.weights,
#' then the \code{adaptive.lasso} weights will be multiplied by the \code{group.weights} vector
#' @param adaptive.fused Flag indicating whether or not to use adaptive fused lasso weights. 
#' @param gamma power to raise the MLE estimated weights by for the adaptive lasso. defaults to 1
#' @param standardize Logical flag for \code{x} variable standardization, prior to fitting the models. 
#' The coefficients are always returned on the original scale. Default is \code{standardize = TRUE}. If 
#' variables are in the same units already, you might not wish to standardize. 
#' @param intercept Should intercept(s) be fitted (\code{default = TRUE}) or set to zero (\code{FALSE})
#' @param compute.se Should standard errors be computed? If \code{TRUE}, then models are re-fit with no penalization and the standard
#' errors are computed from the refit models. These standard errors are only theoretically valid for the
#' adaptive lasso (when \code{adaptive.lasso} is set to \code{TRUE})
#' @param rho ADMM parameter. must be a strictly positive value. By default, an appropriate value is automatically chosen
#' @param dynamic.rho \code{TRUE}/\code{FALSE} indicating whether or not the rho value should be updated throughout the course of the ADMM iterations
#' @param maxit integer. Maximum number of ADMM iterations. Default is 500.
#' @param abs.tol absolute convergence tolerance for ADMM iterations for the relative dual and primal residuals.
#' Default is \code{10^{-5}}, which is typically adequate.
#' @param rel.tol relative convergence tolerance for ADMM iterations for the relative dual and primal residuals.
#' Default is \code{10^{-5}}, which is typically adequate.
#' @param irls.maxit integer. Maximum number of IRLS iterations. Only used if \code{family != "gaussian"}. Default is 100.
#' @param irls.tol convergence tolerance for IRLS iterations. Only used if \code{family != "gaussian"}. Default is \code{10^{-5}}.
#' @return An object with S3 class "oglasso"
#' @export
#' @examples
#' library(vennLasso)
#' 
#' set.seed(123)
#' n.obs <- 1e3
#' n.vars <- 50
#'
#' true.beta <- c(rep(0,2), 1, -1, rep(0, 8), 0.5, -0.5, 1, rep(0, 35))
#'
#' x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' y <- rnorm(n.obs, sd = 3) + drop(x %*% true.beta)
#'
#' groups <- c(list(c(1,2), c(2,3), c(3,4,5), 5:10, 6:12, 7:15), lapply(16:50, function(x) x))
#' 
#' \dontrun{
#' fit <- oglasso(x = x, y = y, group = groups)
#' }
oglasso <- function(x, y,
                    delta                = NULL,
                    group,
                    fused                = NULL,
                    family               = c("gaussian", "binomial", "coxph"),
                    nlambda              = 100L,
                    lambda               = NULL,
                    lambda.min.ratio     = NULL,
                    lambda.fused         = 0,
                    alpha                = NULL, # not used currently
                    group.weights        = NULL,
                    adaptive.lasso       = FALSE,
                    adaptive.fused       = FALSE,
                    penalty.factor       = NULL,
                    penalty.factor.fused = NULL,
                    gamma                = 1,
                    standardize          = TRUE,
                    intercept            = TRUE,
                    compute.se           = FALSE,
                    rho                  = NULL,
                    dynamic.rho          = TRUE,
                    maxit                = 500L,
                    abs.tol              = 1e-5,
                    rel.tol              = 1e-5,
                    irls.tol             = 1e-5,
                    irls.maxit           = 100L)
{

    family <- match.arg(family)
    this.call = match.call()

    nvars <- ncol(x)
    nobs <- nrow(x)
    y <- drop(y)
    dimy <- dim(y)
    leny <- ifelse(is.null(dimy), length(y), dimy[1])
    stopifnot(leny == nobs)

    if (is.null(delta))
    {
        if (family == "coxph")
        {
            stop("delta must be provided for coxph")
        } else 
        {
            delta <- rep(0L, nrow(x))
        }
    }

    if(isTRUE(rho <= 0))
    {
        stop("rho should be positive")
    }

    if (missing(group)) {
        stop("Must specify group structure.")
    } else if (is.matrix(group) | inherits(group, "Matrix")) {
        if (nvars != ncol(group)) {
            stop("group matrix must have same number of rows as variables. i.e. ncol(x) = nrow(group).")
        }

        group <- as(group, "CsparseMatrix")
        group <- as(group, "dgCMatrix")
    } else if (is.list(group) & !is.data.frame(group)) {
        tmp.vc <- numeric(nvars)
        group <- as(vapply(group, function(x) {tmp.vc[x] <- 1; return(tmp.vc) }, tmp.vc), "CsparseMatrix")
        group <- as(group, "dgCMatrix")
    } else {
        stop("Supply either a list or matrix. No data frames allowed.")
    }

    # add unpenalized group if intercept wanted
    if (intercept & (family != "coxph"))
    {
        group  <- rbind(0, group)
        if (!is.null(fused))
        {
            fused <- cbind(0, fused)

        }
    }

    # Seek for variables which were not
    # included in any group and add them
    # to a final group which will be
    # unpenalized.
    rSg <- Matrix::rowSums(group) == 0
    addZeroGroup <- any(rSg)
    if (addZeroGroup) {
        group <- cbind(group, (1 * rSg))
    }
    group.idx <- as.integer(c(0, cumsum(Matrix::colSums(group))))
    ngroups <- ncol(group)
    if (!is.null(group.weights)) {
        stopifnot(length(group.weights) == ngroups)
        group.weights <- as.double(group.weights)
    } else {
        if (adaptive.lasso)
        {
            ## don't apply group size penalty by
            ## default if user asks for
            ## adaptive lasso weights
            group.weights <- sqrt(Matrix::colSums(group)) #rep(1, ncol(group))
        } else
        {
            group.weights <- sqrt(Matrix::colSums(group))
        }

    }
    if (addZeroGroup)
    {
        # Do not penalize variables in last group
        # which represents the group of variables
        # that were not put into any groups.
        # this will override any adaptive lasso
        # weighting later on (which is what we want)
        group.weights[ngroups] <- 0


        if (is.null(penalty.factor))
        {
            penalty.factor <- c(rep(1, ncol(group) - 1), 0)
        } else
        {
            penalty.factor <- c(penalty.factor, 0)
        }

        if (is.null(penalty.factor.fused))
        {
            penalty.factor.fused <- c(0, rep(1, nvars))
        } else
        {
            penalty.factor.fused <- c(0, penalty.factor.fused)
        }


    } else
    {
        if (is.null(penalty.factor))
        {
            penalty.factor <- rep(1, ncol(group))
        }

        if (is.null(penalty.factor.fused))
        {
            penalty.factor.fused <- rep(1, nvars)
        }
    }



    if (length(penalty.factor) != ncol(group))
        stop("penalty.factor must have length equal to the number of groups")

    if (length(penalty.factor.fused) != (nvars + intercept) )
        stop("penalty.factor.fused must have length equal to the number of columns of x")


    varnames <- colnames(x)
    if (is.null(varnames)) {
        varnames <- paste("V", seq(nvars), sep = "")
    }
    is.sparse = FALSE
    if (inherits(x, "sparseMatrix")) {
        is.sparse = TRUE
        x <- as(x, "CsparseMatrix")
        x <- as(x, "dgCMatrix")
    }

    if(is.null(lambda.min.ratio))
    {
        lambda.min.ratio <- ifelse(nrow(x) < ncol(x), 0.05, 0.0001)
    }

    gamma       <- as.double(gamma)
    irls.tol    <- as.double(irls.tol)
    abs.tol     <- as.double(abs.tol)
    rel.tol     <- as.double(rel.tol)
    dynamic.rho <- as.logical(dynamic.rho)
    irls.maxit  <- as.integer(irls.maxit)
    rho         <- if(is.null(rho))  -1.0  else  as.numeric(rho)

    adaptive.lasso <- as.logical(adaptive.lasso)
    adaptive.fused <- as.logical(adaptive.fused)

    if (is.null(lambda)) {
        if (lambda.min.ratio >= 1 | lambda.min.ratio <= 0) {
            stop("lambda.min.ratio should be less than 1 and greater than 0")
        }
        lambda.min.ratio <- as.double(lambda.min.ratio)
        nlambda <- as.integer(nlambda)[1]
        lambda <- numeric(0L)
    } else {
        if (any(lambda < 0)) {
            stop("lambdas should be non-negative")
        }
        lambda <- as.double(rev(sort(lambda)))
        nlambda <- as.integer(length(lambda))
    }

    opts <- list(maxit       = maxit,
                 eps_abs     = abs.tol,
                 eps_rel     = rel.tol,
                 rho         = rho,
                 dynamic_rho = dynamic.rho,
                 irls_maxit  = irls.maxit,
                 irls_tol    = irls.tol)


    fit <- oglasso.fit(family,
                       is.sparse,
                       x, y,
                       delta, group,
                       nlambda, lambda,
                       lambda.min.ratio,
                       fused, lambda.fused,
                       alpha,
                       group.weights,
                       adaptive.lasso,
                       adaptive.fused,
                       penalty.factor,
                       penalty.factor.fused,
                       gamma,
                       group.idx,
                       ngroups,
                       standardize,
                       intercept,
                       compute.se,
                       opts)
    
    fit$beta    <- fit$beta[, 1:fit$nkeep, drop = FALSE]
    fit$lambda  <- fit$lambda[1:fit$nkeep]
    fit$loss    <- fit$loss[1:fit$nkeep]
    fit$niter   <- fit$niter[1:fit$nkeep]
    fit$nlambda <- fit$nkeep

    fit$family  <- family
    fit$call    <- this.call
    fit
}



oglasso.fit <- function(family,
                        is.sparse,
                        x, y,
                        delta, group,
                        nlambda, lambda,
                        lambda.min.ratio,
                        fused, lambda.fused,
                        alpha,
                        group.weights,
                        adaptive.lasso,
                        adaptive.fused,
                        penalty.factor,
                        penalty.factor.fused,
                        gamma,
                        group.idx,
                        ngroups,
                        standardize,
                        intercept,
                        compute.se,
                        opts)
{
    if (is.null(fused))
    {
        if (is.sparse) {
            fit <- .Call("oglasso_fit_sparse",
                         x_                    = x,
                         y_                    = y,
                         group_                = group,
                         family_               = family,
                         nlambda_              = nlambda,
                         lambda_               = lambda,
                         lambda_min_ratio_     = lambda.min.ratio,
                         group_weights_        = group.weights,
                         adaptive_lasso_       = adaptive.lasso,
                         penalty_factor_       = penalty.factor,
                         gamma_                = gamma,
                         group_idx             = group.idx,
                         ngroups_              = ngroups,
                         standardize_          = standardize,
                         intercept_            = intercept,
                         compute_se_           = compute.se,
                         opts_                 = opts,
                         PACKAGE = "vennLasso")
        } else {
            fit <- .Call("admm_oglasso_dense",
                         x_                = x,
                         y_                = y,
                         delta_            = delta,
                         group_            = group,
                         family_           = family,
                         nlambda_          = nlambda,
                         lambda_           = lambda,
                         lambda_min_ratio_ = lambda.min.ratio,
                         group_weights_    = group.weights,
                         adaptive_lasso_   = adaptive.lasso,
                         penalty_factor_   = penalty.factor,
                         gamma_            = gamma,
                         group_idx         = group.idx,
                         ngroups_          = ngroups,
                         standardize_      = standardize,
                         intercept_        = intercept,
                         compute_se_       = compute.se,
                         opts_             = opts,
                         PACKAGE           = "vennLasso")
        }
    } else
    {
        if (is.sparse) {
            fit <- .Call("oglasso_fit_sparse",
                         x_                    = x,
                         y_                    = y,
                         group_                = group,
                         family_               = family,
                         nlambda_              = nlambda,
                         lambda_               = lambda,
                         lambda_min_ratio_     = lambda.min.ratio,
                         group_weights_        = group.weights,
                         adaptive_lasso_       = adaptive.lasso,
                         adaptive_fused_       = adaptive.fused,
                         penalty_factor_       = penalty.factor,
                         penalty_factor_fused_ = penalty.factor.fused,
                         gamma_                = gamma,
                         group_idx             = group.idx,
                         ngroups_              = ngroups,
                         standardize_          = standardize,
                         intercept_            = intercept,
                         compute_se_           = compute.se,
                         opts_                 = opts,
                         PACKAGE = "vennLasso")
        } else {
            fit <- .Call("admm_oglasso_fused_dense",
                         x_                    = x,
                         y_                    = y,
                         delta_                = delta,
                         group_                = group,
                         fused_                = fused,
                         family_               = family,
                         nlambda_              = nlambda,
                         lambda_               = lambda,
                         lambda_fused_         = lambda.fused,
                         lambda_min_ratio_     = lambda.min.ratio,
                         group_weights_        = group.weights,
                         adaptive_lasso_       = adaptive.lasso,
                         adaptive_fused_       = adaptive.fused,
                         penalty_factor_       = penalty.factor,
                         penalty_factor_fused_ = penalty.factor.fused,
                         gamma_                = gamma,
                         group_idx             = group.idx,
                         ngroups_              = ngroups,
                         standardize_          = standardize,
                         intercept_            = intercept,
                         compute_se_           = compute.se,
                         opts_                 = opts,
                         PACKAGE               = "vennLasso")
        }
    }
    fit$fused  <- fused
    class(fit) <- c(class(fit), "oglasso")
    fit
}
