#' Perform symmetry tests
#'
#' This is a generic function used to perform symmetry tests on numeric vectors
#' or objects of class lm (linear models) and objects of class fGARCH (GARCH
#' mdels fitted with the fGarch package).
#'
#' The tests are performed using bootstrap procedures or using asymptotic
#' results, where applicable. Currently, two methods of generating a bootstrap
#' sample from the null distribution are available. The "sign" method generates
#' the bootstrap sample by multiplying the existing sample by -1 or 1 at random
#' (with equal probabilities), essentially randomizing the sign of the data,
#' giving a symmetric distribution. The "reflect" method reflects the sample
#' around zero and samples length(x) elements with replacement. In practice, it
#' has been shown that the "sign" method is almost always better, thus is the
#' default.
#'
#' For numeric data, the tests can be performed around a known (parameter "mu")
#' or unknown centre. When the centre is known, the bootstrap method gives the
#' same results as a Monte Carlo simulation of the p value, for tests which are
#' distribution free. For unknown centre (when mu = NULL), bootstrap must be
#' used and the estimate of the centre used is the trimmed mean, with trim
#' parameter "trim". By default, the mean is taken (trim = 0).
#'
#' For linear models, the tests are based on a bootstrap procedure as in
#' \insertCite{Allison}{symmetry} and are used to test the symmetry of the
#' residuals around zero.
#'
#' For GARCH models (must be fitted with the 'fGarch' package), the tests are also
#' based on bootstrap and test for symmetry of the residuals around zero. An
#' approximation of the bootstrap procedure is available where the residuals are
#' treated as iid data, which is much faster and has been shown to give similar
#' results to the default bootstrap procedure (described in
#' \insertCite{Klar2012}{symmetry}).
#'
#' For a comparison of the performance of various tests of symmetry around an
#' unknown centre, see \insertCite{UNKcentre}{symmetry}).
#'
#' @param x an object of class numeric, lm or fGARCH
#' @param stat a character vector indicating the test statistic to be used (see
#'   \link[=TestStatistics]{Available Test Statistics})
#' @param mu the centre parameter around which to test symmetry
#' @param bootstrap a logical indicating whether to use bootstrap
#' @param B the number of bootstrap replications
#' @param boot_method the method of bootstrap sample generation (see Details)
#' @param trim the trim value used for estimating the centre (as used in "mean")
#' @param k the k parameter of the statistic, ignored if the test statistic
#'   doesn't depend on a parameter (see \link[=TestStatistics]{Test Statistics})
#' @param burn the number of elements to remove from the beginning of the time
#'   series for testing
#' @param approximate a logical indicating whether to use the faster approximate
#'   bootstrap method (see Details)
#' @param ... not used
#' @return An object of class "htest" containing the results of the testing.
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#'
#' # IID samples
#' x <- rnorm(50)
#' symmetry_test(x, "MOI", bootstrap = FALSE, k = 3, mu = 0)
#' symmetry_test(x, "MOI", bootstrap = TRUE, k = 3, mu = 0)
#' symmetry_test(x, "MOI", bootstrap = TRUE, k = 3)
#' x <- rsl(50, alpha = 1.5)
#' symmetry_test(x, "MOI", bootstrap = FALSE, k = 3, mu = 0)
#' symmetry_test(x, "MOI", bootstrap = TRUE, k = 3, mu = 0)
#' symmetry_test(x, "MOI", bootstrap = TRUE, k = 3)
#'
#' # Linear models
#' lin_model <- lm(dist ~ speed, cars)
#' symmetry_test(lin_model, "B1")
#'
#' # Garch models
#' library(fGarch)
#' specskew19 = fGarch::garchSpec(model = list(omega = 0.1,
#'                                     alpha = 0.3,
#'                                     beta = 0.3,
#'                                     skew = 1.9),
#'                                     cond.dist = "snorm")
#'
#' x <- fGarch::garchSim(specskew19, n = 500)
#' g <- fGarch::garchFit(~garch(1,1), x, cond.dist = "QMLE",
#'               include.mean = FALSE, trace = FALSE)
#' \donttest{symmetry_test(g, "FM", B=400, burn = 100)} # slower
#' \donttest{symmetry_test(g, "FM", B=400, burn = 100, approximate = TRUE)}
#'
#' @export
symmetry_test <- function(x, ...) {
  UseMethod("symmetry_test", x)
}
#' @rdname symmetry_test
#' @export
symmetry_test.default <- function(x, stat, mu = NULL,
                                  bootstrap = TRUE, B = 1000,
                                  boot_method = c("sign", "reflect"),
                                  trim = 0, k = 0, ...) {
  if (!is.numeric(x) && !is.logical(x)) {
    stop("Symmetry tests can be applied only to numeric vectors or objects of
         classes lm and fGARCH.")
  }

  xname <- deparse(substitute(x))
  boot_method <- match.arg(boot_method)

  stat_fun <- match.fun(stat, descend = FALSE)
  pass_k <- "k" %in% names(formals(stat))
  if (pass_k && k == 0)
    stop("Argument 'k' not specified.")

  MU <- NULL
  params <- NULL

  if (bootstrap && is.null(mu)) { # UNKNOWN MEAN BOOTSTRAP
    boot <- boot_sample(x, trim, B, boot_method, stat, k)

    MU <- c(mu = mean(x, trim = trim))
    xc <- x - MU
    tval <- if(pass_k) stat_fun(xc, k = k) else stat_fun(xc)

    pval <- mean(abs(boot) >= abs(tval))

    METHOD <- c("Symmetry test",
                "Null hypothesis: Data is symmetric")
    params <- c("B" = B)

  } else if (bootstrap && !is.null(mu)){ # KNOWN MEAN BOOTSTRAP
    boot <- boot_sample(x, mu, B, boot_method, stat, k, TRUE)

    xc <- x - mu
    tval <- if(pass_k) stat_fun(xc, k = k) else stat_fun(xc)

    pval <- mean(abs(boot) >= abs(tval))

    METHOD <- c("Symmetry test",
                paste("Null hypothesis: Data is symmetric around", mu))
    params <- c("B" = B)

  } else { # ASYMPTOTIC RESULTS
    if (!stat %in% names(asymptotic_distributions))
      stop(paste("The asymptotic distributions are available only for these statistics:",
                 paste(sort(names(asymptotic_distributions)), collapse = ", ")))
    if (is.null(mu)) {
      warning("Argument mu not provided. Using 0 by default.")
      mu <- 0
    }
    xc <- x - mu
    tval <- if(pass_k) stat_fun(xc, k = k) else stat_fun(xc)
    asymp_fun <- asymptotic_distributions[[stat]]
    pdist <- if(pass_k) asymp_fun(k) else asymp_fun
    pval <- 2 * (1 - pdist(abs(tval)))
    METHOD <- c("Symmetry test",
                paste("Null hypothesis: Data is symmetric around", mu))
  }

  names(tval) <- stat
  if(pass_k) params <- c(k=k, params)

  obj <- list(method = METHOD,
              statistic = tval,
              parameters = params,
              p.value = pval,
              estimate = MU,
              data.name = xname)
  class(obj) <- "htest"
  obj
}

#' @rdname symmetry_test
#' @export
symmetry_test.lm <- function(x, stat, B = 1000,
                             boot_method = c("sign", "reflect"), k = 0, ...) {
  boot_method <- match.arg(boot_method)
  model <- x
  stat_fun <- match.fun(stat, descend = FALSE)
  pass_k <- "k" %in% names(formals(stat))
  if (pass_k && k == 0)
    stop("Argument 'k' not specified.")

  X <- model.matrix(model)
  yfit <- fitted(model)
  res <- residuals(model)
  boot <- boot_sample_lm(X, yfit, res, B, boot_method, stat, k)
  tval <- if(pass_k) stat_fun(res, k = k) else stat_fun(res)
  names(tval) <- stat
  pval <- mean(abs(boot) >= abs(tval))

  xname <- paste("Residuals from model", deparse(substitute(x)))
  METHOD <- c("Symmetry test of linear model residuals",
              "Null hypothesis: The residuals are symmetric around 0")
  params <- c("B" = B)
  if(pass_k) params <- c(k=k, params)

  obj <- list(method = METHOD,
              statistic = tval,
              parameters = params,
              p.value = pval,
              data.name = xname)
  class(obj) <- "htest"
  obj
}

#' @rdname symmetry_test
#' @export
symmetry_test.fGARCH <- function(x, stat, B = 1000, burn = 0,
                                boot_method = c("sign", "reflect"), k = 0,
                                approximate = FALSE, ...) {
  boot_method <- match.arg(boot_method)
  iid <- approximate
  model <- x
  stat_fun <- match.fun(stat, descend = FALSE)

  pass_k <- "k" %in% names(formals(stat))
  if (pass_k && k == 0)
    stop("Argument 'k' not specified.")

  res <- fGarch::residuals(model, standardize = TRUE)

  null_sample_fun <- switch(boot_method,
                            "sign" = randomize_sign,
                            "reflect" = reflected_boot)

  coefs <- fGarch::coef(model)
  omega <- coefs["omega"]
  alpha <- coefs[grepl("alpha", names(coefs))]
  beta <- coefs[grepl("beta", names(coefs))]

  ts <- as.numeric(model@data)
  cfit <- as.numeric(fGarch::fitted(model))

  not_burned <- length(ts) - burn
  if (not_burned <= 0)
    stop("Number of points to burn is larger than the series length")

  tail_res <- tail(res, not_burned)

  tval <- if(pass_k) stat_fun(tail_res, k = k) else stat_fun(tail_res)
  names(tval) <- stat

  params <- c()

  boot <- replicate(B, {
    boot_res <- null_sample_fun(res, 0)

    if (!iid) {
      boot_y <- simulate_garch(boot_res, ts, cfit, omega, alpha, beta)

      boot_model <- fGarch::garchFit(model@formula, boot_y,
                             cond.dist = "QMLE", include.mean = FALSE,
                             trace = FALSE)
      new_res <- tail(fGarch::residuals(boot_model, standardize = TRUE), not_burned)
    } else {
      new_res <- tail(boot_res, not_burned)
    }
    if(pass_k) stat_fun(new_res, k = k) else stat_fun(new_res)
  })

  pval <- mean(abs(boot) >= abs(tval))
  params <- c("B" = B, params)
  xname <- paste("Residuals from model", deparse(substitute(x)))
  METHOD <- c("Symmetry test of GARCH model residuals",
              "Null hypothesis: The residuals are symmetric around 0")
  if(pass_k) params <- c(k=k, params)

  obj <- list(method = METHOD,
              statistic = tval,
              parameters = params,
              p.value = pval,
              data.name = xname)
  class(obj) <- "htest"
  obj
}
