#' Check if a predictor is dichotomous, adopted from package \code{circGLM}
#'
#' @param x A character or numerical vector to be tested.
#'
#' @return A logical, \code{TRUE} if the \code{x} has dummy coding (0, 1),
#'   \code{FALSE} otherwise.
#'
is.dichotomous <- function(x) {
  n_unique_x <- length(unique(x))
  if (n_unique_x == 2) {
    if (all(x == 0 | x == 1)) {
      return(TRUE)
    } else {
      warning("A predictor might be dichotomous but not 0|1.")
    }
  } else if (n_unique_x > 2 & n_unique_x < 8) {
    warning(
      paste(
        "A predictor has between 3 and 7 unique values.",
        "It might be categorical with multiple categories",
        "but without dummy coding."
      )
    )
  } else if (n_unique_x == 1) {
    stop("A predictor had only a single unique value.")
  }
  FALSE
}

#' Fitting conditional binary quantile models
#'
#' The main function for running the conditional binary quantile model. The function returns a cbq \code{cbq} object that can be further investigated using standard functions such as \code{plot}, \code{print}, \code{coef}, and \code{predict}.
#'
#' The model can be passed either as a combination of a \code{formula} and a data frame \code{data}, as in \code{lm()}.
#'
#' Convergence diagnotics can be performed using either \code{print(object, "mcmc")} or \code{plot(object, "mcmc")}.
#'
#' @param formula An object of class "Formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param q The quantile value.
#' @param vi Indicating whether variantional inference should be used instead of MCMC sampling procedure.
#' @param nsim The number of iterations.
#' @param grad_samples Passed to \code{\link[rstan]{vb}} (positive integer), the number of samples for Monte Carlo estimate of gradients, defaulting to 1.
#' @param elbo_samples Passed to \code{\link[rstan]{vb}} (positive integer), the number of samples for Monte Carlo estimate of ELBO (objective function), defaulting to 100. (ELBO stands for "the evidence lower bound".)
#' @param tol_rel_obj Passed to \code{\link[rstan]{vb}} (positive double), the convergence tolerance on the relative norm of the objective, defaulting to 0.01.
#' @param output_samples Passed to \code{\link[rstan]{vb}} (positive integer), number of posterior samples to draw and save, defaults to 1000.
#' @param burnin The number of burnin iterations.
#' @param thin Thinning parameter.
#' @param CIsize The size of confidence interval.
#' @param nchain The number of parallel chains.
#' @param seeds Random seeds to replicate the results.
#' @param inverse_distr If FALSE, the ALD will not be reversed. The default is FALSE.
#' @param offset Offset values to enhance sampling stability. The default value is 1e-20.
#' @param mc_core Indicating whether the estimation will be run in multiple parallel chains. The default is TRUE.
#'
#' @return A \code{cbq} object, which can be further analyzed with its associated \code{\link{plot.cbq}}, \code{\link{coef.cbq}} and \code{\link{print.cbq}} functions.
#'
#' An object of class \code{cbq} contains the following elements
#'
#'   \describe{
#'
#'   \item{\code{Call}}{The matched call.}
#'   \item{\code{formula}}{Symbolic representation of the model.}
#'   \item{\code{q}}{The quantile value.}
#'   \item{\code{nsim}}{The number of MCMC iterations.}
#'   \item{\code{burnin}}{The number of burnin periods.}
#'   \item{\code{thin}}{Thinning.}
#'   \item{\code{seeds}}{Random seeds.}
#'   \item{\code{CIsize}}{The size of confidence interval.}
#'   \item{\code{data}}{Data used.}
#'   \item{\code{x}}{Covaraites used.}
#'   \item{\code{y}}{The dependent variable.}
#'   \item{\code{xnames}}{Names of the covariates.}
#'   \item{\code{stanfit}}{Outputs from stan.}
#'   \item{\code{sampledf}}{A matrix of posterior samples.}
#'   \item{\code{summaryout}}{A summary based on posterior samples.}
#'   \item{\code{npars}}{Number of covariates.}
#'   \item{\code{ulbs}}{Lower and upper confidence bounds.}
#'   \item{\code{means}}{Estimates at the mean.}
#'   \item{\code{vi}}{Indicating whether variational inference has been performed.}
#'   \item{\code{output_samples}}{Sample outputs.}
#'   \item{\code{fixed_var}}{Variables estimated using fixed effects.}
#'   \item{\code{random_var}}{Variables estimated using random effects.}
#'   \item{\code{xq}}{Variables indicating the choice sets.}
#'
#'
#' }
#'
#' @importFrom stats aggregate coef model.frame model.matrix quantile
#'
#' @export
#'
#' @author
#' Xiao Lu
#'
#' @references
#' Lu, Xiao. (2020). Discrete Choice Data with Unobserved Heterogeneity: A Conditional Binary Quantile Model. Political Analysis, 28(2), 147-167. https://doi.org/10.1017/pan.2019.29
#'
#'
#' @examples
#' # Simulate the data
#' x <- rnorm(50)
#' y <- ifelse(x > 0, 1, 0)
#' dat <- as.data.frame(cbind(y, x))
#'
#' # Estimate the CBQ model
#' model <- cbq(y ~ x, dat, 0.5, nchain = 1, mc_core = FALSE)
#'
#' # Show the results
#' print(model)
#' coef(model)
#' plot(model)
#'
#'
cbq <- function(formula,
                data,
                q = NULL,
                vi = FALSE,
                nsim = 1000,
                grad_samples = 1,
                elbo_samples = 100,
                tol_rel_obj = 0.01,
                output_samples = 2000,
                burnin = NULL,
                thin = 1,
                CIsize = .95,
                nchain = 1,
                seeds = 12345,
                inverse_distr = FALSE,
                offset = 1e-20,
                mc_core = TRUE) {
  if (is.null(burnin))
    burnin <- floor(nsim / 2)
  if (burnin < 0)
    stop("Burn-in must be non-negative.")
  if (thin < 1)
    stop("Thinning factor must be positive.")
  if (CIsize <= 0)
    stop("Confidence interval size 'CIsize' must be positive.")
  if (CIsize > 1)
    stop(paste0("Confidence interval size 'CIsize' ",
                "can not be larger than 1."))
  if (missing(formula) | missing(data)) {
    stop(paste0("Formula and data should be given."))
  }
  if (is.null(q)) {
    warning("Quantile is not specified. The default quantile 0.5 is used.")
    q <- 0.5
  }
  if (length(q) > 1)
    stop("Only a single quantile value is allowed.")
  if (q >= 1 |
      q <= 0)
    stop("The specified quantile is out of range. The value must be in (0,1).")
    
  old <- options(stringsAsFactors = FALSE)
  if (mc_core == TRUE){
        # activate multicore in stan
        options(mc.cores = parallel::detectCores())
    }
  f <- Formula::Formula(formula)
  data <- stats::model.frame(f, data)
  y <- c(as.matrix(stats::model.frame(f, data)[, 1]))
  if (length(f)[2] >= 2 ) {
      xq <- as.matrix(stats::model.matrix(f, data, rhs = 2)[, -1])
      nq <- dim(xq)[2]
  } else {
      xq <- NULL
      nq <- 0
  }
  if (nq == 1) {
      indx = as.integer(as.factor(c(xq)))
      tmp = stats::aggregate(y, list(indx), sum)
      if (any(tmp[, 2] != 1) == TRUE ) {
          deleteids = which(indx %in% tmp[which(tmp[,2] != 1),1] )
          data = data[-deleteids,]
          warning(
          "In each choice set, there must be only one chosen observation. Multiple 1s or all 0s are not allowed in any choice set. Problematic choice sets have been deleted."
          )
      }
  }
  if (length(f)[2] >= 2 ) {
      xq <- as.matrix(stats::model.matrix(f, data, rhs = 2)[, -1])
      nq <- dim(xq)[2]
  } else {
      xq <- NULL
      nq <- 0
  }
  if (nq == 1) {
      indx = as.integer(as.factor(c(xq)))
  }
  y <- c(as.matrix(stats::model.frame(f, data)[, 1]))
  x = stats::model.matrix(f, data)
  
  fixed_var = NULL
  random_var = NULL
  if (length(f)[2] == 3 ) {
      if (dim(as.matrix(stats::model.matrix(f, data, rhs = 3)[,-1]))[2] == 0 ) {
          fixed_var = NULL
      } else if (dim(as.matrix(stats::model.matrix(f, data, rhs = 3)[,-1]))[2] == 1 ) {
          fixed_var = c(stats::model.matrix(f, data, rhs = 3)[,-1])
      } else {
          stop("Misspecified fixed intercepts!")
      }
  }
  
  if (length(f)[2] == 4 ) {
      if (dim(as.matrix(stats::model.matrix(f, data, rhs = 3)[,-1]))[2] == 0 ) {
          fixed_var = NULL
      } else if (dim(as.matrix(stats::model.matrix(f, data, rhs = 3)[,-1]))[2] == 1 ) {
          fixed_var = c(stats::model.matrix(f, data, rhs = 3)[,-1])
      } else {
          stop("Misspecified fixed intercepts!")
      }
      if (dim(as.matrix(stats::model.matrix(f, data, rhs = 4)[,-1]))[2] == 0 ) {
          random_var = NULL
      } else if (dim(as.matrix(stats::model.matrix(f, data, rhs = 4)[,-1]))[2] == 1 ) {
          random_var = c(stats::model.matrix(f, data, rhs = 4)[,-1])
      } else {
          stop("Misspecified random intercepts!")
      }
  }

  n_covariate <- dim(x)[2]
  N <- length(y)

  if (nq > 1)
    stop("The number of index variables must be equal to or less than one.")
  if (!is.dichotomous(y))
    stop("The dependent variable must be binary with values in {0,1}.")

  
  
  
  if (is.null(fixed_var) & is.null(random_var)){
      type = "regular"
  } else if (is.null(fixed_var)) {
      random_var = as.integer(as.factor(random_var))
      if (length(unique(random_var)) < 2) {
          stop("The value of the random indicator is unique.")
      }
      type = "random"
  } else if (is.null(random_var)){
      fixed_var = as.integer(as.factor(fixed_var))
      if (length(unique(fixed_var)) < 2) {
          stop("The value of the fixed indicator is unique.")
      }
      type = "fixed"
  } else {
      random_var = as.integer(as.factor(random_var))
      fixed_var = as.integer(as.factor(fixed_var))
      if (length(unique(fixed_var)) < 2) {
          stop("The value of the fixed indicator is unique.")
      }
      if (length(unique(random_var)) < 2) {
          stop("The value of the random indicator is unique.")
      }
      type = "panel"
  }

  if (nq == 0) {
    if (inverse_distr == FALSE) {
        if (type == "regular") {
            stanmodel <- stanmodels$cbqbv
            datlist <- list(
            N = N,
            Y = c(y),
            D = n_covariate,
            X = x,
            q = q,
            offset = offset
            )
        } else if (type == "random") {
            stanmodel = stanmodels$cbqrandombv
            datlist <- list(
            N = N,
            Y = c(y),
            D = n_covariate,
            X = x,
            q = q,
            offset = offset,
            N_person = length(unique(random_var)),
            person = random_var
            )
        } else if (type == "fixed") {
            stanmodel = stanmodels$cbqfixbv
            datlist <- list(
            N = N,
            Y = c(y),
            D = n_covariate,
            X = x,
            q = q,
            offset = offset,
            N_wave = length(unique(fixed_var)),
            wave = fixed_var
            )
        } else {
            stanmodel = stanmodels$cbqpanelbv
            datlist <- list(
            N = N,
            Y = c(y),
            D = n_covariate,
            X = x,
            q = q,
            offset = offset,
            N_wave = length(unique(fixed_var)),
            wave = fixed_var,
            N_person = length(unique(random_var)),
            person = random_var
            )
        }
    } else {
        if (type == "regular") {
            stanmodel <- stanmodels$cbqb
            datlist <- list(
            N = N,
            Y = c(y),
            D = n_covariate,
            X = x,
            q = q,
            offset = offset
            )
        } else if (type == "random") {
            stanmodel = stanmodels$cbqrandomb
            datlist <- list(
            N = N,
            Y = c(y),
            D = n_covariate,
            X = x,
            q = q,
            offset = offset,
            N_person = length(unique(random_var)),
            person = random_var
            )
        } else if (type == "fixed") {
            stanmodel = stanmodels$cbqfixb
            datlist <- list(
            N = N,
            Y = c(y),
            D = n_covariate,
            X = x,
            q = q,
            offset = offset,
            N_wave = length(unique(fixed_var)),
            wave = fixed_var
            )
        } else {
            stanmodel = stanmodels$cbqpanelb
            datlist <- list(
            N = N,
            Y = c(y),
            D = n_covariate,
            X = x,
            q = q,
            offset = offset,
            N_wave = length(unique(fixed_var)),
            wave = fixed_var,
            N_person = length(unique(random_var)),
            person = random_var
            )
        }
    }
  } else {
    if (inverse_distr == FALSE) {
        if (type == "regular") {
            x <- x[order(indx, y), ]
            y <- y[order(indx, y)]
            stanmodel <- stanmodels$cbqdv
            datlist <- list(
            N = N,
            Y = c(y),
            D_common = n_covariate,
            X_common = x,
            N_indx = length(unique(indx)),
            ind = indx,
            q = q,
            offset = offset
            )
        } else if (type == "random") {
            x <- x[order(indx, y), ]
            y <- y[order(indx, y)]
            random_var = random_var[order(indx, y)]
            stanmodel <- stanmodels$cbqrandomdv
            datlist <- list(
            N = N,
            Y = c(y),
            D_common = n_covariate,
            X_common = x,
            N_indx = length(unique(indx)),
            ind = indx,
            q = q,
            offset = offset,
            N_person = length(unique(random_var)),
            person = random_var
            )
        } else if (type == "fixed") {
            x <- x[order(indx, y), ]
            y <- y[order(indx, y)]
            fixed_var = fixed_var[order(indx, y)]
            stanmodel <- stanmodels$cbqfixdv
            datlist <- list(
            N = N,
            Y = c(y),
            D_common = n_covariate,
            X_common = x,
            N_indx = length(unique(indx)),
            ind = indx,
            q = q,
            offset = offset,
            N_wave = length(unique(fixed_var)),
            wave = fixed_var
            )
        } else {
            x <- x[order(indx, y), ]
            y <- y[order(indx, y)]
            fixed_var = fixed_var[order(indx, y)]
            random_var = random_var[order(indx, y)]
            stanmodel <- stanmodels$cbqpaneldv
            datlist <- list(
            N = N,
            Y = c(y),
            D_common = n_covariate,
            X_common = x,
            N_indx = length(unique(indx)),
            ind = indx,
            q = q,
            offset = offset,
            N_person = length(unique(random_var)),
            person = random_var,
            N_wave = length(unique(fixed_var)),
            wave = fixed_var
            )
        }

    } else {
        if (type == "regular") {
            x <- x[order(indx, y), ]
            y <- y[order(indx, y)]
            stanmodel <- stanmodels$cbqd
            datlist <- list(
            N = N,
            Y = c(y),
            D_common = n_covariate,
            X_common = x,
            N_indx = length(unique(indx)),
            ind = indx,
            q = q,
            offset = offset
            )
        } else if (type == "random") {
            x <- x[order(indx, y), ]
            y <- y[order(indx, y)]
            random_var = random_var[order(indx, y)]
            stanmodel <- stanmodels$cbqrandomd
            datlist <- list(
            N = N,
            Y = c(y),
            D_common = n_covariate,
            X_common = x,
            N_indx = length(unique(indx)),
            ind = indx,
            q = q,
            offset = offset,
            N_person = length(unique(random_var)),
            person = random_var
            )
        } else if (type == "fixed") {
            x <- x[order(indx, y), ]
            y <- y[order(indx, y)]
            fixed_var = fixed_var[order(indx, y)]
            stanmodel <- stanmodels$cbqfixd
            datlist <- list(
            N = N,
            Y = c(y),
            D_common = n_covariate,
            X_common = x,
            N_indx = length(unique(indx)),
            ind = indx,
            q = q,
            offset = offset,
            N_wave = length(unique(fixed_var)),
            wave = fixed_var
            )
        } else {
            x <- x[order(indx, y), ]
            y <- y[order(indx, y)]
            fixed_var = fixed_var[order(indx, y)]
            random_var = random_var[order(indx, y)]
            stanmodel <- stanmodels$cbqpaneld
            datlist <- list(
            N = N,
            Y = c(y),
            D_common = n_covariate,
            X_common = x,
            N_indx = length(unique(indx)),
            ind = indx,
            q = q,
            offset = offset,
            N_person = length(unique(random_var)),
            person = random_var,
            N_wave = length(unique(fixed_var)),
            wave = fixed_var
            )
        }

    }
  }

  if (type == "regular") {
     pars <- "beta"
  } else if (type == "random") {
     pars <- c("beta","beta_ind","sigma_beta_ind")
  } else if (type == "fixed") {
     pars <- c("beta","beta_wave")
  } else {
     pars <- c("beta","beta_ind","sigma_beta_ind","beta_wave")
  }
  
  if (vi == FALSE) {
      stanout <- rstan::sampling(
      stanmodel,
      data = datlist,
      pars = pars,
      seed = seeds,
      iter = nsim,
      thin = thin,
      warmup = burnin,
      chains = nchain
      )
  } else {
      stanout <- rstan::vb(
      stanmodel,
      data = datlist,
      pars = pars,
      seed = seeds,
      grad_samples = grad_samples,
      elbo_samples = elbo_samples,
      tol_rel_obj = tol_rel_obj,
      output_samples = output_samples
      )
  }

  summaryout <- rstan::summary(stanout)$summary
  sampledf <- as.data.frame(stanout)[, 1:n_covariate]
  names(sampledf) <- colnames(x)

  out = list()
  class(out) <- c("cbq", class(out))
  out$Call <- match.call()
  out$formula <- formula
  out$q <- q
  out$nsim <- nsim
  out$burnin <- burnin
  out$thin <- thin
  out$seeds <- seeds
  out$CIsize  <- CIsize
  out$data   <- data
  out$x    <- x
  out$y <- y
  out$xnames <- colnames(x)
  out$stanfit <- stanout
  out$sampledf <- sampledf
  out$summaryout <- summaryout
  out$npars <- n_covariate
  out$ulbs <-
    apply(sampledf, 2, stats::quantile, probs = c((1 - CIsize) / 2, 1 - (1 -
                                                                           CIsize) / 2))
  out$means <- summaryout[1:n_covariate, 1]
  out$vi <- vi
  out$output_samples <- output_samples
  out$fixed_var <- fixed_var
  out$random_var <- random_var
  out$xq = xq
  on.exit(options(old), add = TRUE)
  return(out)


}
