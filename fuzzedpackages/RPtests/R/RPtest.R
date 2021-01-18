#' Goodness of fit tests for potentially high-dimensional linear models
#'
#' Can test for the significance of (potentially large) groups of predictors and
#' the presence of nonlinearity or heteroscedasticity in the context of both low
#' and high-dimensional linear models. Outputs a p-value. Also allows for the
#' calibration of arbitrary goodness of fit tests via specification of
#' \code{RPfunction}.
#'
#' @param x Input matrix with \code{nobs} rows, each an observation vector.
#' @param y Response vector.
#' @param resid_type Type of residuals used for the test (see details below).
#'   Use \code{Lasso} when the null model is high-dimensional; otherwise use
#'   \code{OLS}.
#' @param test Type of departure from the linear model to test for (see details
#'   below). Ignored if \code{RPfunction} is given.
#' @param x_alt If \code{test} is \code{group}, this gives the set of variables
#'   whose significance we wish to ascertain, after controlling for those in
#'   \code{x}. If \code{RPfunction} is given, it is the input matrix passed to
#'   the function \code{RPfunction}.
#' @param RPfunction A residual prediction (RP) function that must permit
#'   calling as \code{RPfunction(x_alt, resid)} where \code{resid} is a numeric
#'   vector with \code{nobs} components. The output must be either a single
#'   number or a numeric vector (in the latter case \code{RPfunction} would
#'   encode a number of RP functions).
#' @param B The number of bootstrap samples to use - note the p-value produced
#'   will always be at least 1/B.
#' @param rand_gen A function to generate the simulated errors up to an unknown
#'   scale factor. It must permit calling as \code{rand_gen(nobs*B)}. Determines
#'   the form of errors in the null model. The default \code{rnorm} equates to a
#'   null of a (sparse) Gaussian linear model. Setting \code{rand_gen=NULL}
#'   resamples residuals to generate simulated errors and approximates a null of
#'   i.i.d. errors with unknown distribution.
#' @param noise_matrix An optional matrix whose columns are the simulated errors to use.
#'   Note that \code{B} and \code{rand_gen} will be ignored if this is supplied.
#' @param mc.cores The number of cores to use. Will always be 1 in Windows.
#' @param nfolds Number of folds to use when performing cross-validation to
#'   obtain \code{beta_est}, the initial estimate of the vector of regression
#'   coefficients, via Lasso estimation.
#' @param nperms Number of permutations of the data for which \code{nfolds}
#'   cross-validation is to be performed. Thus in total prediction errors on
#'   \code{nfolds*nperms} folds are averaged over.
#' @param beta_est An optional user-supplied estimate.
#' @param resid_only If \code{TRUE} only outputs the residuals without applying
#'   an RP function.
#' @param output_all In addition to the p-value, gives further output (see Value
#'   below).
#' @param verbose Whether to print addition information.
#' @details
#'   The function works by first computing residuals from a regression of
#'   y on x. Next \code{B} sets of errors generated through \code{rand_gen} are
#'   added to a signal derived from \code{beta_est} and aritificial residuals
#'   are computed. The option \code{resid_only=TRUE} then outputs these
#'   residuals along with the original residuals, scaled to have l_2-norm
#'   squared equal to \code{nobs}. The residuals in question are OLS residuals
#'   when \code{resid_type=OLS} (case a - for use when the null hypothesis is
#'   low-dimensional so the number of columns of \code{x} is smaller than
#'   \code{nobs-1}), and square-root / scaled Lasso residuals otherwise (case
#'   b). The options for \code{test} then apply different functions to the
#'   residuals as described below.
#'   \describe{
#'     \item{\code{nonlin}}{In case (a), the test statistic is the RSS (residual
#'     sum of squares) of a \code{\link[randomForest]{randomForest}} fit from
#'     regressing the residuals on to \code{x}; case (b) is similar but the OOB
#'     error is used and the regression is carried out on the equicorrelation set
#'     rather than all of \code{x}.}
#'     \item{\code{group}}{\code{x_alt} is first residualised with
#'     respect to \code{x} by (a) OLS or (b) \code{\link{sparse_proj}}. Then the
#'     RSS from Lasso fits from regressions of the residuals on to \code{x_alt}
#'     are used.}
#'     \item{\code{hetero}}{Uses the RSS from Lasso fits from
#'     regressions of the squared residuals to the equicorrelation set (b) or all
#'     of \code{x} (a).}
#'   }
#' @return When \code{resid_only=FALSE} and \code{output_all=FALSE}, the output
#'   is a single p-value. Otherwise, a list with some of the following
#'   components is returned (\code{resid_only=FALSE} causes the last two
#'   components to be omitted):
#'   \describe{
#'     \item{\code{p-value}}{p-value}
#'     \item{\code{beta_est}}{estimated vector of regression coefficients
#'     \code{beta_est}}
#'     \item{\code{sigma_est}}{set to 1 when \code{resid_type=OLS};
#'     otherwise the normalised root-RSS derived from
#'     \code{beta_est} used in generated the simulated errors}
#'     \item{\code{resid}}{scaled residuals}
#'     \item{\code{resid_sim}}{simulated scaled residuals}
#'     \item{\code{test}}{the test statistic(s) - may be a vector if multiple RP
#'     functions are being used such as when \code{test=group}}
#'     \item{\code{test_sim}}{a list of simulated test statistics}
#'   }
#' @references Shah, R. D., Buhlmann, P. (2016) \emph{Goodness of fit tests for
#'   high-dimensional linear models} \url{http://arxiv.org/abs/1511.03334}
#' @seealso \code{\link{RPtest_single}} and \code{\link{sqrt_lasso}}
#' @examples
#' # Testing for nonlinearity
#' set.seed(1)
#' x <- scale(matrix(runif(100*200), 100, 200))
#' y <- x[, 1] + x[, 1]^4 + rnorm(nrow(x))
#' out <- RPtest(x, y, test="nonlin", B=9L, nperms=2, resid_type = "Lasso")
#'
#' # Testing significance of a group
#' y <- x[, 1:5] %*% rep(1, 5) + x[, 151] + rnorm(nrow(x))
#' (out <- RPtest(x[, 1:150], y, test="group", x_alt = x[, 151:200], B=9L, nperms=1))
#'
#' # Testing for heteroscedasticity
#' x <- scale(matrix(runif(250*100), 250, 100))
#' hetero_sig <- x[, 1] + x[, 2]
#' var_vec <- hetero_sig - min(hetero_sig) + 0.01
#' var_vec <- var_vec / mean(var_vec)
#' sd_vec <- sqrt(var_vec)
#' y <- x[, 1:5] %*% rep(1, 5) + sd_vec*rnorm(nrow(x))
#' (out <- RPtest(x, y, test="hetero", B=9L, nperms=1))
#' @export
#' @useDynLib RPtests
#' @importFrom Rcpp sourceCpp
#' @import stats
RPtest <- function(x, y, resid_type=c("Lasso", "OLS"),
                   test=c("nonlin", "group", "hetero"),
                   x_alt, RPfunction=NULL, B=49L, rand_gen = rnorm,
                   noise_matrix = NULL, mc.cores=1L,
                   nfolds=5L, nperms=2L, beta_est = NULL, resid_only=FALSE,
                   output_all=FALSE,
                   verbose = FALSE) {

  # Params not currently used
  sigma_star_fac = 1
  lam0=NULL

  resid_type <- match.arg(resid_type)

  # Check x, y and standardise
  if (!is.matrix(x)) stop("x should be a matrix with at least one column")
  np <- dim(x)
  if (is.null(np) | (np[2] < 1L))
    stop("x should be a matrix with at least one column")
  n <- as.integer(np[1])
  p <- as.integer(np[2])
  if ((p >= n-1) & (resid_type=="OLS"))
    stop("When resid_type=OLS we must have ncol(x) < nrow(x)-1; try setting resid_type=Lasso")
  y <- as.numeric(y)
  if (length(y) != n) stop("y must have nrow(x) components")

  y <- y - mean(y)
  x <- scale(x, scale=FALSE)

  # Compute simulated residuals
  if (resid_type != "OLS") {
    if (is.null(beta_est)) {
      nfolds <- as.integer(nfolds)
      if (nfolds < 2L) stop("nfolds must be an integer at least 2.")
      # Compute beta_est
      cv_glmnet_out <- glmnet::cv.glmnet(x, y, nfolds=nfolds)
      mse_vec <- cv_glmnet_out$cvm
      cv_lambda <- cv_glmnet_out$lambda
      cv_glmnet_all <- parallel::mclapply(seq_len(nperms-1L), function(i)
        glmnet::cv.glmnet(x, y, nfolds=nfolds, lambda=cv_lambda), mc.cores=mc.cores)
      for (i in seq_len(nperms-1L)) {
        mse_vec <- mse_vec[seq_along(cv_glmnet_all[[i]]$cvm)] # shorten mse_vec if necessary
        mse_vec <- cv_glmnet_all[[i]]$cvm + mse_vec
      }
      lambda_min <- cv_lambda[which.min(mse_vec)]
      beta_est <- as.numeric(glmnet::coef.glmnet(cv_glmnet_out$glmnet.fit, s=lambda_min)[-1, ])
      # coef gives an intercept term
    }

    # Comp
    sig_est <- as.numeric(x %*% beta_est)
    sigma_hat <- comp_sigma_est(x, y, beta_est)
    sigma_est <- sigma_star_fac * sigma_hat

    # Simulated errors
    if (is.null(rand_gen)) {
      rand_gen <- function(n_samp) {
        sample((y - sig_est)/sigma_hat, n_samp, replace=TRUE)
      }
    }
    error_mat <- sigma_est *
      if(is.null(noise_matrix)) matrix(rand_gen(n*B), n, B) else noise_matrix

    if(verbose){
      if(!is.null(noise_matrix)){
        cat(paste("\nNoise matrix has been provided (",
                  ncol(noise_matrix), "simulations)."))
      }
    }


    # Simulated residuals
    resid_sim <- resid_gen_lasso(x, signal = sig_est,
                                 sigma_est=sigma_est, B=B, rand_gen=rand_gen,
                                 mc.cores=mc.cores,
                                 error_mat = error_mat,
                                 lam0=NULL)

    # True residuals
    resid <- resid_lasso(x, y, lam0=NULL)

    rm(error_mat)
  } else {
    # OLS
    resid_out <- resid_ols(x, y, incl_proj=TRUE)
    resid <- resid_out$resid
    resid_sim <- resid_gen_ols(x, proj=resid_out$proj, B=B,
                               rand_gen=rand_gen,
                               error_mat = noise_matrix)
    rm(resid_out)
  }

  # Output residuals if resid_only is TRUE
  if (resid_only) return(list("call"=match.call(),
                              "beta_est"=if(exists("beta_est")) beta_est else NULL,
                              "sigma_est"=if(exists("sigma_est")) sigma_est else NULL,
                              "resid"=resid,
                              "resid_sim"=resid_sim))

  # Set RPfunction
  user_RP_func <- TRUE
  if (is.null(RPfunction)) {
    user_RP_func <- FALSE
    test <- match.arg(test)
    if (test == "nonlin") {
      x_alt <- x
      if (resid_type=="OLS") {
        RPfunction <- rf_test
      } else {
        RPfunction <- rf_test_S
      }
    } else if (test == "group") {
      if (missing(x_alt)) stop("x_alt must be specified when test = groups")
      RPfunction <- lasso_test
      if (resid_type=="OLS") {
        x_alt <- non_sparse_proj(x_null=x, x_alt=x_alt)
      } else {
        x_alt <- sparse_proj(x_null=x, x_alt=x_alt, mc.cores=mc.cores)
      }
    } else if (test == "hetero") {
      x_alt <- x
      if (resid_type == "OLS") {
        RPfunction <- sq_lasso_test
      } else {
        RPfunction <- sq_lasso_test_S
      }
    }
  } else {
    if (missing(x_alt)) stop("x_alt must be specified when the RPfunction argument is used")
  }

  # Apply RPfunction to residuals
  if (user_RP_func == TRUE) {
    # When RPfunction is user-supplied
    test <- RPfunction(x_alt, resid)
    test_sim <- parallel::mclapply(1:B, function(b) RPfunction(x_alt, resid_sim[, b]), mc.cores=mc.cores)
  } else {
    # When RPfunction is not user-supplied
    # Note the RPfunctions used when test is specified take different arguemnts
    test_out <- RPfunction(x_alt, resid, x_null=x)
    if (is.list(test_out)) {
      test <- test_out$test
      test_sim <- RPfunction(x_alt, resid_sim, test_out$test_aux, mc.cores=mc.cores, x_null=x)
    } else {
      test <- test_out
      test_sim <- parallel::mclapply(1:B, function(b) RPfunction(x_alt, resid_sim[, b], x_null=x), mc.cores=mc.cores)
    }
  }

  # Compute p-values
  pval_computed <- pval(test, test_sim)

  # Output
  if (!output_all) return(pval_computed)
  # output_all == TRUE
  return(list("call"=match.call(),
              "p-value"=pval_computed,
              "beta_est"=if(exists("beta_est")) beta_est else NULL,
              "sigma_est"=if(exists("sigma_est")) sigma_est else NULL,
              "resid"=resid,
              "resid_sim"=resid_sim,
              "test"=test,
              "test_sim"=test_sim))
}
