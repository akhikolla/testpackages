context("mcmc functions")

test_that("mcmc_mra", {
    y <- rnorm(10)
    X <- matrix(1:20, 10, 2)
    params <- list(
        n_mcmc    = 5,
        n_adapt   = 5,
        n_thin    = 1,
        n_message = 5
    )
    locs <- matrix(1:20, 10, 2)

    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )
    inits <- list(
        beta = rep(0, ncol(X)),
        sigma2 = 1,
        alpha = rep(0, 520), ## use M = 2 and n_coarse_grid = 4
        tau2  = 1:2
    )

    out <- mcmc_mra(y, X, locs, params, priors, inits = inits, M = 2, n_coarse_grid = 4)
    expect_s3_class(out, "mcmc_mra")
    expect_s3_class(out$MRA, "mra_wendland_2d")
    expect_error(mcmc_mra(y, X, locs, params, M = -3), "the number of resolutions M must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_coarse_grid = -3), "n_coarse_grid must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_neighbors = -3), "n_neighbors must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_padding = -3), "n_padding must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_cores = -3), "n_cores must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_cores = 1:4), "n_cores must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_cores = NA), "n_cores must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_cores = NULL), "n_cores must be a positive integer")

    expect_error(mcmc_mra(y, X, locs, params, n_chain = 1:4), "n_chain must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_chain = NA), "n_chain must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_chain = NULL), "n_chain must be a positive integer")

    expect_error(mcmc_mra(y, X, locs, params, use_spam = FALSE), "Only support use_spam = TRUE")

    expect_error(mcmc_mra(y, X), 'argument "locs" is missing, with no default')
    expect_error(mcmc_mra(y, X, locs), 'argument "params" is missing, with no default')
    y <- matrix(1:10, 5, 2)
    expect_error(mcmc_mra(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep(NA, 10)
    expect_error(mcmc_mra(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep("aaa", 10)
    expect_error(mcmc_mra(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep(TRUE, 10)
    expect_error(mcmc_mra(y, X, locs, params), "y must be a numeric vector of length N.")

    y <- rnorm(10)
    X <- array(0, dim=c(2, 2, 2))
    expect_error(mcmc_mra(y, X, locs, params), "X must have the same number of rows as the length of y.")
    X <- matrix(0, 5, 3)
    expect_error(mcmc_mra(y, X, locs, params), "X must have the same number of rows as the length of y.")
    X <- array(0, dim=c(10, 2, 2))
    expect_error(mcmc_mra(y, X, locs, params), "X must be a numeric matrix with N rows.")
    X <- matrix(NA, 10, 3)
    expect_error(mcmc_mra(y, X, locs, params), "X must be a numeric matrix with N rows.")

    X <- matrix(1:20, 10, 2)
    locs <- matrix(0, 10, 3)
    expect_error(mcmc_mra(y, X, locs, params), "locs must be a numeric matrix with N rows and 2 columns.")
    locs <- matrix(NA, 10, 2)
    expect_error(mcmc_mra(y, X, locs, params), "locs must be a numeric matrix with N rows and 2 columns.")
    locs <- matrix(1:20, 10, 2)

    params$n_mcmc <- -10
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- NA
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- NULL
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- "aaa"
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- 50.5
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")

    params$n_mcmc <- 500
    params$n_adapt <- -10
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- NA
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- NULL
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- "aaa"
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- 50.5
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")

    params$n_adapt <- 500
    params$n_thin <- -10
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- NA
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- NULL
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- "aaa"
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- 50.5
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")


    params$n_thin <- 1
    params$n_message <- -10
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- NA
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- NULL
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- "aaa"
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- 50.5
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")

    params$n_message = 100
    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )

    priors$mu_beta <- rep(0, ncol(X) + 1)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep(NA, ncol(X))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep("aa", ncol(X))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- matrix(rep(0, ncol(X)))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep(0, ncol(X))

    priors$Sigma_beta <- diag(ncol(X) + 1)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- matrix(1, ncol(X), ncol(X))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigam_beta <- diag(NA, ncol(X))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- matrix(rep(0, ncol(X), ncol(X)))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- diag(ncol(X))

    priors$alpha_sigma2 <- rep(1, 2)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- -1
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- "aa"
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- 0
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- NA
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- 1

    priors$beta_sigma2 <- rep(1, 2)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- -1
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- "aa"
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- 0
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- NA
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- 1

    priors$alpha_tau2 <- rep(1, 2)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- -1
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- "aa"
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- 0
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- NA
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- 1

    priors$beta_tau2 <- rep(1, 2)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- -1
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- "aa"
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- 0
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- NA
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- 1

    priors <- NULL

    expect_error(mcmc_mra(y, X, locs, params, use_spam = "aaa"), "use_spam must be either TRUE or FALSE.")
    expect_error(mcmc_mra(y, X, locs, params, use_spam = 3), "use_spam must be either TRUE or FALSE.")
    expect_error(mcmc_mra(y, X, locs, params, use_spam = NA), "use_spam must be either TRUE or FALSE.")

    expect_error(mcmc_mra(y, X, locs, params, verbose = "aaa"), "verbose must be either TRUE or FALSE.")
    expect_error(mcmc_mra(y, X, locs, params, verbose = 3), "verbose must be either TRUE or FALSE.")
    expect_error(mcmc_mra(y, X, locs, params, verbose = NA), "verbose must be either TRUE or FALSE.")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_beta = "TRUE")), "If specified, sample_beta must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_beta = NA)), "If specified, sample_beta must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_beta = 3)), "If specified, sample_beta must be TRUE or FALSE")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_tau2 = "TRUE")), "If specified, sample_tau2 must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_tau2 = NA)), "If specified, sample_tau2 must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_tau2 = 3)), "If specified, sample_tau2 must be TRUE or FALSE")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_lambda = "TRUE")), "If specified, sample_lambda must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_lambda = NA)), "If specified, sample_lambda must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_lambda = 3)), "If specified, sample_lambda must be TRUE or FALSE")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_sigma2 = "TRUE")), "If specified, sample_sigma2 must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_sigma2 = NA)), "If specified, sample_sigma2 must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_sigma2 = 3)), "If specified, sample_sigma2 must be TRUE or FALSE")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_alpha = "TRUE")), "If specified, sample_alpha must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_alpha = NA)), "If specified, sample_alpha must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_alpha = 3)), "If specified, sample_alpha must be TRUE or FALSE")


    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = 3)), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = rep("aaa", ncol(X)))), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = rep(NA, ncol(X)))), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = matrix(1:ncol(X)))), "initial value for beta must be a numeric vector of length p")

    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = 1:3)), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = matrix(1:ncol(X)))), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = NA)), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = -3)), "initial value for sigma2 must be a positive numeric value")

    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(alpha = 2:ncol(out$alpha))), "initial value for alpha must be positive numeric vector of length equal to the number of all grid points")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(alpha = matrix(1:ncol(out$alpha)))), "initial value for alpha must be positive numeric vector of length equal to the number of all grid points")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(alpha = c(2:ncol(out$alpha), NA))), "initial value for alpha must be positive numeric vector of length equal to the number of all grid points")

    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = 1:3)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = 1)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = matrix(1:ncol(X)))), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = NA)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = -3)), "initial value for tau2 must be a positive numeric vector of length M")

    # mcmc_mra(y, X, locs, params, priors,  M = 2, n_coarse_grid = 4)

    # expect_error(mcmc_mra(y, X, params), 'argument "priors" is missing, with no default')
    # priors <- default_priors_pg_lm(y, X)
    # out <- pg_lm(y, X, params, priors)
    # expect_true(class(out) == "pg_lm")

})

