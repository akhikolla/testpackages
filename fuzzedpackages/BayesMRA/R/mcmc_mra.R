#' Bayesian Multi-resolution Spatial Regression
#'
#' this function runs Markov Chain Monte Carlo to estimate the Bayesian multi-resolution spatial regression model.
#' @param y is a \eqn{n}{n} vector of Gaussian data.
#' @param X is a \eqn{n \times p}{n x p} matrix of fixed effects (like latitude, elevation, etc)
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param params is a list of parameter settings. The list
#' \code{params} must contain the following values:
#' * \code{n_adapt}: A positive integer number of adaptive MCMC iterations.
#' * \code{n_mcmc}: A positive integer number of total MCMC iterations
#' post adaptation.
#' * \code{n_thin}: A positive integer number of MCMC iterations per saved
#' sample.
#' * \code{n_message}: A positive integer number of frequency of iterations
#'  to output a progress message. For example, \code{n_message = 50}
#'  outputs progress messages every 50 iterations.
#' @param priors is the list of prior settings.
#' @param M The number of resolutions.
#' @param n_neighbors The expected number of neighbors for each interior basis function. This determines the basis radius parameter.
#' @param n_coarse_grid The number of basis functions in one direction (e.g. \code{n_coarse_grid = 10} results in a \eqn{10 \times 10}{10x10} course grid which is further extended by the number of additional padding basis functions given by \code{n_padding}.
#' @param n_padding The number of additional boundary points to add on each boundary. For example, n_padding = 5 will add 5 boundary knots to the both the left  and right side of the grid).
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the initial values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param verbose Should verbose output be printed? Typically this is only useful for troubleshooting.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE} (see spam and Matrix packages for details).
#' @param n_chain is the MCMC chain id. The default is 1.
#'
#' @examples
#' set.seed(111)
#' ## genereate the locations
#' locs <- matrix(runif(20), 10, 2)
#' ## generate some covariates and regression coefficients
#' X <- cbind(1, matrix(rnorm(30), 10, 3))
#' beta <- rnorm(ncol(X))
#'
#' ## simulate the MRA process
#' M <- 2
#' MRA <- mra_wendland_2d(locs, M = M, n_coarse_grid = 4)
#' W <- do.call(cbind, MRA$W)
#'
#' n_dims   <- rep(NA, length(MRA$W))
#' dims_idx <- NULL
#' for (i in 1:M) {
#'     n_dims[i] <- ncol(MRA$W[[i]])
#'     dims_idx  <- c(dims_idx, rep(i, n_dims[i]))
#' }
#' ## set up the process precision matrices
#' Q_alpha <- make_Q_alpha_2d(sqrt(n_dims), c(0.9, 0.8))
#' Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2 = c(2, 4))
#'
#' ## add in constraints so each resolution has random effects that sum to 0
#' A_constraint <- sapply(1:M, function(i){
#'     tmp = rep(0, sum(n_dims))
#'     tmp[dims_idx == i] <- 1
#'     return(tmp)
#' })
#' a_constraint <- rep(0, M)
#' alpha <- as.vector(spam::rmvnorm.prec.const(
#'     n = 1,
#'     mu = rep(0, nrow(W)),
#'     Q = Q_alpha_tau2,
#'     A = t(A_constraint),
#'     a = a_constraint))
#' ## define the data
#' y <-  as.vector(X %*% beta + W %*% alpha + rnorm(10))
#'
#' ## define the params for MCMC fitting
#' params <- list(
#'     n_mcmc    = 5,
#'     n_adapt   = 5,
#'     n_thin    = 1,
#'     n_message = 5)
#'
#' ## define the model priors
#' priors <- list(
#'     alpha_tau2   = 1,
#'     beta_tau2    = 1,
#'     alpha_sigma2 = 1,
#'     beta_sigma2  = 1,
#'     mu_beta      = rep(0, ncol(X)),
#'     Sigma_beta   = 5 * diag(ncol(X)))
#'
#' ## Fit the MRA model using MCMC
#' out     <- mcmc_mra(
#'     y             = y,
#'     X             = X,
#'     locs          = locs,
#'     params        = params,
#'     priors        = priors,
#'     M             = 2,
#'     n_coarse_grid = 4,
#'     n_cores       = 1L,
#'     verbose       = FALSE
#' )
#'
#' @export
#'
#' @importFrom mvnfast rmvn
#' @importFrom fields rdist
#' @importFrom Matrix Cholesky
#' @importFrom stats lm rgamma
#' @import spam

mcmc_mra <- function(
    y,
    X,
    locs,
    params,
    priors        = NULL,
    M             = 4,
    n_neighbors   = 68,
    n_coarse_grid = 10,
    n_padding     = 5L,
    n_cores       = 1L,
    inits         = NULL,
    config        = NULL,
    verbose       = FALSE,
    use_spam      = TRUE, ## use spam or Matrix for sparse matrix operations
    n_chain       = 1
) {

    ##
    ## Run error checks
    ##

    if (!is_numeric_vector(y, length(y)))
        stop("y must be a numeric vector of length N.")
    if (length(y) != nrow(X))
        stop("X must have the same number of rows as the length of y.")
    if (!is_numeric_matrix(X, length(y), ncol(X)))
        stop("X must be a numeric matrix with N rows.")
    if (!is_numeric_matrix(locs, length(y), 2))
        stop("locs must be a numeric matrix with N rows and 2 columns.")

    ## check the params list
    if (!is_positive_integer(params$n_adapt, 1))
        stop("params must contain a positive integer n_adapt.")
    if (!is_positive_integer(params$n_mcmc, 1))
        stop("params must contain a positive integer n_mcmc.")
    if (!is_positive_integer(params$n_thin, 1))
        stop("params must contain a positive integer n_thin.")
    if (!is_positive_integer(params$n_message, 1))
        stop("params must contain a positive integer n_message.")

    params$n_adapt   <- as.integer(params$n_adapt)
    params$n_mcmc    <- as.integer(params$n_mcmc)
    params$n_thin    <- as.integer(params$n_thin)
    params$n_message <- as.integer(params$n_message)

    ## check the priors list


    ## check mcmc input
    if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
        stop("verbose must be either TRUE or FALSE.")
    }

    if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
        stop("use_spam must be either TRUE or FALSE.")
    }

    if (!is_positive_integer(n_cores, 1)) {
        stop("n_cores must be a positive integer")
    }

    n_cores <- as.integer(n_cores)

    if (!is_positive_integer(n_chain, 1)) {
        stop("n_chain must be a positive integer")
    }

    n_chain <- as.integer(n_chain)

    ##
    ## setup config
    ##

    ## do we sample the functional relationship parameters? This is primarily
    ## used to troubleshoot model fitting using simulated data

    sample_beta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_beta']])) {
            sample_beta <- config[['sample_beta']]
            if (!is.logical(sample_beta) | is.na(sample_beta))
                stop('If specified, sample_beta must be TRUE or FALSE')
        }
    }


    # sample_rho <- TRUE
    # if (!is.null(config)) {
    #     if (!is.null(config[['sample_rho']])) {
    #         sample_rho <- config[['sample_rho']]
    #     }
    # }

    ## do we sample the MRA variance parameter? This is primarily
    ## used to troubleshoot model fitting using simulated data
    sample_tau2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_tau2']])) {
            sample_tau2 <- config[['sample_tau2']]
            if (!is.logical(sample_tau2) | is.na(sample_tau2))
                stop('If specified, sample_tau2 must be TRUE or FALSE')
        }
    }

    sample_lambda <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_lambda']])) {
            sample_lambda <- config[['sample_lambda']]
            if (!is.logical(sample_lambda) | is.na(sample_lambda))
                stop('If specified, sample_lambda must be TRUE or FALSE')
        }
    }

    ## do we sample the nugger variance parameter
    sample_sigma2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma2']])) {
            sample_sigma2 <- config[['sample_sigma2']]
            if (!is.logical(sample_sigma2) | is.na(sample_sigma2))
                stop('If specified, sample_sigma2 must be TRUE or FALSE')
        }
    }

    ## do we sample the latent intensity parameter alpha
    sample_alpha <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_alpha']])) {
            sample_alpha <- config[['sample_alpha']]
            if (!is.logical(sample_alpha) | is.na(sample_alpha))
                stop('If specified, sample_alpha must be TRUE or FALSE')
        }
    }

    ## add in a counter for the number of regularized Cholesky
    num_chol_failures <- 0

    N      <- length(y)
    # n_time <- dim(Y)[3]
    p      <- ncol(X)

    ##
    ## setup MRA spatial basis
    ##

    MRA <- mra_wendland_2d(
        locs          = locs,
        M             = M,
        n_coarse_grid = n_coarse_grid,
        n_padding     = n_padding,
        n_neighbors   = n_neighbors,
        use_spam      = use_spam
    )

    # MRA      <- mra_wendland_2d(locs, M, n_max_fine_grid = n_max_fine_grid, use_spam = use_spam)
    W_list   <- MRA$W
    tW_list  <- vector(mode = 'list', length = M)
    tWW_list <- vector(mode = 'list', length = M)
    for (m in 1:M) {
        if (use_spam) {
            tW_list[[m]] <- t(W_list[[m]])
        } else {
            stop("Only support use_spam = TRUE")
            # tW_list[[m]] <- Matrix::t(W_list[[m]])
        }
        tWW_list[[m]] <- tW_list[[m]] %*% W_list[[m]]
    }
    n_dims <- rep(NA, length(W_list))
    dims_idx <- c()
    for (i in 1:M) {
        n_dims[i] <- ncol(W_list[[i]])
        dims_idx <- c(dims_idx, rep(i, n_dims[i]))
    }
    W <- do.call(cbind, W_list)

    tW <- NULL
    if (use_spam) {
        tW <- t(W)
    } else {
        stop ('Only support use_spam = TRUE')
        # tW <- Matrix::t(W)
    }

    tWW <- tW %*% W


    ##
    ## initial values
    ##

    tX  <- t(X)
    tXX <- tX %*% X

    ##
    ## priors for beta
    ##

    mu_beta        <- rep(0, p)
    Sigma_beta     <- 10 * diag(p)

    ## check if priors for mu_beta are specified
    if (!is.null(priors[['mu_beta']])) {
        if(!is_numeric_vector(priors[['mu_beta']], p))
            stop("If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
        if (all(!is.na(priors[['mu_beta']]))) {
            mu_beta <- priors[['mu_beta']]
        }
    }

    ## check if priors for Sigma_beta are specified
    if (!is.null(priors[['Sigma_beta']])) {
        if(!is_sympd_matrix(priors[['Sigma_beta']], p))
            stop("If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
        if (all(!is.na(priors[['Sigma_beta']]))) {
            Sigma_beta <- priors[['Sigma_beta']]
        }
    }
    Sigma_beta_chol <- tryCatch(
        chol(Sigma_beta),
        error = function(e) {
            if (verbose)
                message("The Cholesky decomposition of the prior covariance Sigma_beta was ill-conditioned and mildy regularized.")
            chol(Sigma_beta + 1e-8 * diag(p))
        }
    )
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)

    ##
    ## initialize beta
    ##

    # beta   <- as.vector(rmvn(1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    beta   <- as.vector(lm (y ~ X - 1)$coeff)
    X_beta <- X %*% beta

    ##
    ## intialize a proper CAR structure to initialize the parameter alpha
    ##

    Q_alpha <- make_Q_alpha_2d(sqrt(n_dims), rep(0.999, length(n_dims)), use_spam = use_spam)

    ## prior for lambda is scale_lambda
    ## fix this so it can be set by the user
    scale_lambda <- 0.5

    lambda <- rgamma(M, 0.5, scale_lambda)
    tau2   <- rep(1, M)

    alpha_tau2 <- 1
    beta_tau2  <- 1
    ## check if priors for alpha_tau2 are specified
    if (!is.null(priors[['alpha_tau2']])) {
        if (!is_positive_numeric(priors[['alpha_tau2']], 1))
            stop("If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
        if (all(!is.na(priors[['alpha_tau2']]))) {
            alpha_tau2 <- priors[['alpha_tau2']]
        }
    }
    ## check if priors for beta_tau2 are specified
    if (!is.null(priors[['beta_tau2']])) {
        if (!is_positive_numeric(priors[['beta_tau2']], 1))
            stop("If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
        if (all(!is.na(priors[['beta_tau2']]))) {
            beta_tau2 <- priors[['beta_tau2']]
        }
    }
    #100 * pmin(pmax(1 / rgamma(M, 0.5, lambda), 1), 100)


    Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)

    ##
    ## initialize rho
    ##

    # rho      <- runif(1, -1, 1)

    ##
    ## initialize sigma2
    ##

    alpha_sigma2 <- 1
    beta_sigma2  <- 1
    ## check if priors for alpha_sigma2 are specified
    if (!is.null(priors[['alpha_sigma2']])) {
        if (!is_positive_numeric(priors[['alpha_sigma2']], 1))
            stop("If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
        if (all(!is.na(priors[['alpha_sigma2']]))) {
            alpha_sigma2 <- priors[['alpha_sigma2']]
        }
    }
    ## check if priors for beta_sigma2 are specified
    if (!is.null(priors[['beta_sigma2']])) {
        if (!is_positive_numeric(priors[['beta_sigma2']], 1))
            stop("If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
        if (all(!is.na(priors[['beta_sigma2']]))) {
            beta_sigma2 <- priors[['beta_sigma2']]
        }
    }

    sigma2  <- pmax(pmin(1 / rgamma(1, alpha_sigma2, beta_sigma2), 5), 0.1)
    sigma   <- sqrt(sigma2)

    ##
    ## initialize alpha
    ##

    alpha <- NULL
    if (use_spam) {
        A_alpha <- 1 / sigma2 * tWW + Q_alpha_tau2
        ## precalculate the sparse cholesky structure for faster Gibbs updates
        Rstruct <- chol(A_alpha)
        b_alpha <- 1 / sigma2 * tW %*% (y - X_beta)
        # alpha   <- as.vector(rmvnorm.canonical(1, b_alpha, A_alpha, Rstruct = Rstruct))
        alpha  <- rep(0, sum(n_dims))
    } else {
        stop("The only sparse matrix pacakage available is spam")
    }

    ## define the constraint matrices to ensure each resolution has random effects with mean 0
    A_constraint <- t(
        sapply(1:M, function(i){
            tmp = rep(0, sum(n_dims))
            tmp[dims_idx == i] <- 1
            return(tmp)
        })
    )
    a_constraint <- rep(0, M)


    ## intialize an ICAR structure for fitting alpha
    Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(1, length(n_dims)), use_spam = use_spam)
    Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)


    W_alpha <- W %*% alpha

    ##
    ## sampler config options -- to be added later
    ##
    #

    ## check for initial values
    ##

    ## initial values for beta
    if (!is.null(inits[['beta']])) {
        if(!is_numeric_vector(inits[['beta']], p))
            stop("initial value for beta must be a numeric vector of length p")
        if (all(!is.na(inits[['beta']]))) {
            beta <- inits[['beta']]
        }
    }

    ## intial values for sigma2
    if (!is.null(inits[['sigma2']])) {
        if(!is_positive_numeric(inits[['sigma2']], 1))
            stop("initial value for sigma2 must be a positive numeric value")
        if (all(!is.na(inits[['sigma2']]))) {
            sigma2 <- inits[['sigma2']]
        }
    }

    ## initial values for alpha
    if (!is.null(inits[['alpha']])) {
        if(!is_numeric_vector(inits[['alpha']], length(alpha)))
            stop("initial value for alpha must be positive numeric vector of length equal to the number of all grid points")
        if (all(!is.na(inits[['alpha']]))) {
            alpha <- inits[['alpha']]
        }
    }

    ## initial values for tau2
    if (!is.null(inits[['tau2']])) {
        if(!is_positive_numeric(inits[['tau2']], M) | !is.vector(inits[['tau2']]))
            stop("initial value for tau2 must be a positive numeric vector of length M")
        if (all(!is.na(inits[['tau2']]))) {
            tau2 <- inits[['tau2']]
        }

    }

    ## initial values for rho
    # if (!is.null(inits[['rho']])) {
    #     if (all(!is.na(inits[['rho']]))) {
    #         rho <- inits[['rho']]
    #     }
    # }


    ##
    ## setup save variables
    ##

    n_save       <- params$n_mcmc / params$n_thin
    beta_save    <- matrix(0, n_save, p)
    # rho_save     <- rep(0, n_save)
    tau2_save    <- matrix(0, n_save, M)
    sigma2_save  <- rep(0, n_save)
    alpha_save   <- matrix(0, n_save, sum(n_dims))
    lambda_save    <- matrix(0, n_save, M)

    ##
    ## initialize tuning
    ##

    ##
    ## tuning variables for adaptive MCMC
    ##

    # ## tuning for rho
    # rho_accept       <- 0
    # rho_accept_batch <- 0
    # rho_tune         <- 0.025

    ##
    ## Starting MCMC chain
    ##

    message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations \n")

    for (k in 1:(params$n_adapt + params$n_mcmc)) {
        if (k == params$n_adapt + 1) {
            message("Starting MCMC fitting for chain ", n_chain, ", running for ", params$n_mcmc, " iterations \n")
        }
        if (k %% params$n_message == 0) {
            if (k <= params$n_adapt) {
                message("MCMC adaptation iteration ", k, " for chain ", n_chain)
            } else {
                message("MCMC fitting iteration ", k - params$n_adapt, " for chain ", n_chain)
            }
        }

        ##
        ## sample sigma2
        ##

        if (sample_sigma2) {
            if (verbose)
                message("sample sigma2")

            devs <- y - X_beta - W_alpha
            SS        <- sum(devs^2)
            sigma2 <- 1 / rgamma(1, N / 2 + alpha_sigma2, SS / 2 + beta_sigma2)
            sigma     <- sqrt(sigma2)
        }

        ##
        ## sample beta -- double check these values
        ##

        if (sample_beta) {
            if (verbose)
                message("sample beta")

            ## only use the modern climate state to update beta
            A      <- 1 / sigma2 * tXX + Sigma_beta_inv
            b      <- 1 / sigma2 * tX %*% (y - W_alpha) + Sigma_beta_inv %*% mu_beta
            ## guarantee a symmetric matrix
            A      <- (A + t(A)) / 2
            beta   <- rmvn_arma(A, b)
            X_beta <- X %*% beta
        }

        ##
        ## sample alpha
        ##

        if (sample_alpha) {
            if (verbose)
                message("sample alpha")
            A_alpha <- 1 / sigma2 * tWW + Q_alpha_tau2
            # A_alpha_chol <- chol.spam(A_alpha, Rstruct = Rstruct)
            b_alpha <- 1 / sigma2 * tW %*% (y - X_beta)
            # alpha   <- as.vector(rmvnorm.canonical(1, b_alpha, A_alpha, Rstruct = Rstruct))
            ## sample constrained to sum to 1
            alpha   <- as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha, Rstruct = Rstruct,
                                                         A = A_constraint, a = a_constraint))
        }

        ## update W_alpha
        W_alpha <- W %*% alpha

        ##
        ## sample rho
        ##

        # if (sample_rho) {
        #     if (verbose)
        #         message("sample rho")
        #
        #     rho_vals <- rowSums(
        #         sapply(2:n_time, function(tt) {
        #             t_alpha_Q <- t(alpha[, tt-1]) %*% Q_alpha_tau2
        #             c(
        #                 t_alpha_Q %*% alpha[, tt-1],
        #                 t_alpha_Q %*% alpha[, tt]
        #             )
        #         })
        #     )
        #
        #     a_rho <- rho_vals[1]
        #     b_rho <- rho_vals[2]
        #     rho   <- rtrunc(1, "norm", a = -1, b = 1, mean = b_rho / a_rho, sd = sqrt(1 / a_rho))
        # }

        ##
        ## sample tau2
        ##

        if (sample_tau2) {
            if (verbose)
                message("sample tau2")
            for (m in 1:M) {
                devs    <- alpha[dims_idx == m]
                SS      <- as.numeric(devs %*% (Q_alpha[[m]] %*% devs))
                tau2[m] <- rgamma(1, alpha_tau2 + n_dims[m] / 2, beta_tau2 + SS / 2)
                # tau2[m] <- 1 / rgamma(1, 0.5 + n_dims[m] / 2, lambda[m] + SS / 2)
                # tau2_inv[m] <- 1 / rgamma(1, 0.5 + n_dims[m] / 2, lambda[m] + SS / 2)
                # tau2[m]     <- 1 / tau2_inv[m]
            }
        }

        Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)

        ##
        ## sample lambda
        ##

        if (sample_lambda) {
            if (verbose)
                message("sample lambda")
            for (m in 1:M) {
                lambda[m]     <- rgamma(1, 1, scale_lambda + 1 / tau2[m])
                # lambda[m]     <- rgamma(1, 1, scale_lambda + tau2[m])
            }
        }

        ##
        ## save variables
        ##

        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, ]   <- beta
                # rho_save[save_idx]       <- rho
                tau2_save[save_idx, ]   <- tau2
                sigma2_save[save_idx]   <- sigma2
                alpha_save[save_idx, ]  <- alpha
                lambda_save[save_idx, ] <- lambda
            }
        }

        ##
        ## End of MCMC loop
        ##
    }

    ## print out acceptance rates -- no tuning in this model


    if (num_chol_failures > 0)
        warning("The Cholesky decomposition for the update of alpha was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore.")

    ##
    ## return the MCMC output -- think about a better way to make this a class
    ##

    out <- list(
        beta     = beta_save,
        # rho      = rho_save,
        tau2     = tau2_save,
        sigma2   = sigma2_save,
        alpha    = alpha_save,
        lambda   = lambda_save,
        MRA      = MRA
    )

    class(out) <- "mcmc_mra"

    return(out)
}
