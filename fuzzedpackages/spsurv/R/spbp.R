#' Semiparametric Survival Analysis Using Bernstein Polynomial
#'
#' Fits Bernstein Polynomial based Proportional regression to survival data.
#'
#' @title spbp: The BP Based Survival Analysis Function
#' @param formula a Surv object with time to event, status and explanatory terms.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains) or `rstan::optimizing`.
#' @seealso \code{\link[spsurv]{spbp.default}}
#' @examples
#'
#' library("spsurv")
#' data("veteran") ## imports from survival package
#'
#' fit_mle <- spbp(Surv(time, status) ~ karno + factor(celltype),
#'  data = veteran, model = "po")
#' summary(fit_mle)
#'
#' fit_bayes <- spbp(Surv(time, status) ~ karno + factor(celltype),
#'                   data = veteran, model = "po", approach = "bayes",
#'                    cores = 1, iter = 300, chains = 1,
#'                     priors = list(beta = c("normal(0,4)"),
#'                      gamma = "lognormal(0,4)"))
#'
#' summary(fit_bayes)
#'
#' @rdname spbp
#' @export spbp
#' @seealso  \code{\link[spsurv]{spbp.default}},  \code{\link[spsurv]{bpph}},  \code{\link[spsurv]{bppo}}, \code{\link[spsurv]{bpaft}}, \url{https://mc-stan.org/users/documentation/}
#' @return An object of class 'spbp'.

spbp <- function(formula, ...) {
  UseMethod("spbp", formula)
}
#' @title spbp: The BP Based Semiparametric Survival Analysis Function
#' @param formula a Surv object with time to event, status and explanatory terms
#' @param degree Bernstein Polynomial degree
#' @param data a data.frame object
#' @param approach Bayesian or Maximum Likelihood estimation methods, default is approach = "bayes"
#' @param model Proportional Hazards or Proportional Odds BP based regression, default is model = "ph"
#' @param priors prior settings for the Bayesian approach; `normal` or `cauchy` for beta; `gamma`, `inv_gamma` or `lognormal` for gamma (BP coefficients)
#' @param scale logical; indicates whether to center and scale the data
#' @param ... further arguments passed to or from other methods
#' @param cores number of core threads to use
#' @return An object of class \code{spbp}
#' @method spbp default
#' @export
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty
#' @importFrom MASS ginv
#' @importFrom loo waic loo
#' @importFrom coda  HPDinterval
#' @importFrom stats .getXlevels as.formula contrasts dbeta density dist formula median model.extract pbeta pchisq printCoefmat qnorm rlogis rnorm rweibull sd terms
#'

spbp.default <-
  function(formula, degree, data,
            approach = c("mle", "bayes"),
            model = c("ph", "po", "aft"),
            priors = list(beta = c("normal(0,4)"),
                         gamma = "lognormal(0,10)"),
           scale = TRUE, cores =  parallel::detectCores(),
           ...){

  # ---------------Definitions + error handling  ---------------
  ## tau degree

  if(missing(degree)){
    degree <- ceiling(sqrt(nrow(data)))
  }

  ## Call
  Call <- match.call();

  ## model
  if(length(model) == 3){model_flag = "ph"}
  else{model_flag <- model}

  model <- ifelse(match.arg(model) == "po", 0,
                  ifelse(match.arg(model) == "ph", 1, 2))

  ## approach
  approach_flag <- match.arg(approach) ### saves string input
  approach <- ifelse(approach_flag == "mle", 0, 1)

  handler1() ## error handling #1

  ## terms
  temp <- Call[c(1, aux)] # keep important args
  temp[[1L]] <- quote(stats::model.frame) # model frame call
  special <- c("frailty", "frailty.gamma", "frailty.gaussian", "frailty.t")
  temp$formula <- terms(formula, special, data = data)
  temp$formula <- terms(formula, data = data);

  ## frailty (id, distribution, column)
  handler2()

  ## Priors
  suppressMessages(handler3())

  ## stanArgs
  stanArgs <- list(...)
  handler4()

  ## Model Frame
  mf <- eval(temp, parent.frame())
  Terms <- terms(mf)
  Y <- model.extract(mf, "response") # time-to-event response
  type <- attr(Y, "type")
  handler5()

  # ---------------  Data declaration + definitions ---------------

  ## + sample size + labels
  data.n <- nrow(Y)
  labels <- attributes(temp$formula)$term.labels
  null <- 0

  if (length(labels) > 1){
    X <-  model.matrix(Terms, mf)[, -1]
  }
  else if(length(labels) == 1){
    X <- as.matrix(model.matrix(Terms, mf)[, -1], ncol = data.n)
    colnames(X) <- labels
  }
  else{
    X <- as.matrix(rep(0, data.n), ncol = data.n)
    colnames(X) <- "non-parametric"
    null <- 1
  }

  ## time + status + features
  features <- X
  attr(X, "assign") <- attr(model.matrix(Terms, mf), "assign")[-1]
  attr(X, "contrasts") <- attr(model.matrix(Terms, mf), "contrasts")
  xlevels <- .getXlevels(Terms, mf)
  contrasts <- attr(X, "contrasts")

  assign <- attrassign(X, Terms)

  q <- ncol(X)
  time <- as.vector(Y[,1])
  tau <- max(time)
  status <- as.vector(Y[,2])

  if(scale == T){
    X <- scale(X)
    ## rescaled coefficients (correction)
      means <- array(attr(X, "scaled:center"), dim = q)
      std <- array(attr(X, "scaled:scale"), dim = q)
  }
  else{
    std <- array(1, dim = q)
    means <- array(0, dim = q)
  }

  ## base calculations
  base <- bp.basis(time, degree = degree, tau = tau)

  ## priors to num
  priordist <- sapply(priordist,
                      function(x){
                        switch(x,
                               "normal" = 0,
                               "gamma" = 1,
                               "inv_gamma" = 2,
                               "lognormal" = 3)})

  priordist_beta <- sapply(priordist_beta,
                          function(x){switch(x,
                          "normal" = 0,
                          "cauchy" = 1)})
  ## Recycling the prior specs
  priordist_beta <- array(priordist_beta, dim = q)
  location_beta <- array(as.numeric(location_beta), dim = q)
  scale_beta <- array(as.numeric(scale_beta), dim = q)

  ## standata
  standata <- list(time = time,
                   tau = tau,
                   n = data.n,
                   m = base$degree,
                   q = q,
                   status = status,
                   X = X,
                   B = base$B,
                   b = base$b,
                   approach = approach,
                   M = model,
                   null = null,
                   id = rep(1, data.n),
                   dist = dist,
                   z = rep(0, data.n),
                   priordist = priordist,
                   priorpars = priorpars,
                   priordist_beta = priordist_beta,
                   location_beta = location_beta,
                   scale_beta = scale_beta,
                   std  = std,
                   means = means
  )

  # --------------- Fit  ---------------
  output <- list()
  if(approach == 0){
    spbp.mle(standata = standata, ...)
  }
  else{
    spbp.bayes(standata = standata, ...)
  }
}

spbp.mle <-
  function(standata,
           init = 0,
           hessian = TRUE,
           verbose = FALSE,
           ...){

  e <- parent.frame()
  #variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for(n in vnames) assign(n, get(n, e))

  if(!is.null(frailty_idx)){
    standata$X <- X[, -frailty_idx]
    message("Frailty ignored, change approach to `bayes` for frailty estimation.")
  }

  stanfit <-
    rstan::optimizing(stanmodels$spbp,
                             data = standata,
                             init = init,
                             hessian = hessian,
                             verbose = verbose,
                             ...)
  len <- length(stanfit$par)
  ## stanfit coefficients (beta, nu)
  aux <- stanfit$par
  coef <- aux[names(aux) %in% c(paste0("beta[", 1:q, "]"), paste0("gamma[", 1:degree, "]"))]
  ## regression estimates
  beta <- array(coef[1:q], q)
  gamma_std <- aux[names(aux) %in% paste0("gamma_std[", 1:degree, "]")]

  ## rescaled hessian matrix
  info <- - stanfit$hessian

  names(beta) <- colnames(X)
  names(coef) <- c(names(beta),
                   paste0("gamma", 1:(degree))
  )

  ## singular matrices handler
  #   if(det(-hess) == 0)
  #   stop("Optimizing hesssian matrix is singular!")

  ## rescaled fisher info
  jac <- diag(1/std, q)
  var <- jac %*% blockSolve(info, q)[1:q, 1:q] %*% t(jac)
  rownames(var) <- names(beta)
  colnames(var) <- names(beta)

  if(hessian == FALSE || null == 1){
    stanfit$hessian <- matrix(rep(NA, q^2),
                              ncol = 1:q,
                              nrow = 1:q)
  }
  nulldata <- standata
  nulldata$null <- 1
  nullfit <- rstan::optimizing(stanmodels$spbp,
                               data = nulldata,
                               init = init,
                               hessian = hessian,
                               ...)

  output <- list(coefficients = coef,
                 var = var[1:q, 1:q],
                 loglik = c(nullfit$value, stanfit$value),
                 linear.predictors = c(features %*% beta),
                 means = colMeans(features),
                 method = "optimizing",
                 n = data.n,
                 nevent = sum(status),
                 q = q,
                 terms = Terms,
                 assign = assign,
                 wald.test = coxph.wtest(var[1:q, 1:q], beta)$test,
                 y = Y,
                 formula = formula,
                 xlevels = xlevels,
                 contrasts = contrasts,
                 return_code = stanfit$return_code,
                 tau = tau,
                 call = Call)
  output$call$approach <- approach_flag
  output$call$model <- model_flag

  class(output) <- "spbp"
  message('Priors are ignored due to mle approach.')

  return(output)
}

spbp.bayes <- function(standata,
                       hessian = TRUE,
                       verbose = FALSE,
                       chains = 1,
                       ...){
  e <- parent.frame()
  #variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for(n in vnames) assign(n, get(n, e))

  # bayes
  output <- list(y = Y)

  if(dist == 0){
    output$stanfit <- rstan::sampling(stanmodels$spbp,
                                      data = standata,
                                      verbose = verbose,
                                      chains = chains, cores = cores,
                                      ...)

    samp <- rstan::extract(output$stanfit, pars = c("beta", "gamma"))
    output$pmode <- apply(X = cbind(samp[[1]], samp[[2]]), MARGIN = 2, FUN = mode)
    # output$pmode <- apply(X = samp, MARGIN = 2, FUN = mode)
  }
  else{
    standata$X <- X[, -frailty_idx]
    output$stanfit <- rstan::sampling(stanmodels$spbp_frailty,
                                      data = standata,
                                      verbose = verbose,
                                      chains = chains,
                                      ...)

    output$pmode <- apply(rstan::extract(output$stanfit, c("beta", "gamma")), 2, mode)
  }
  output$loo <- loo::loo(loo::extract_log_lik(output$stanfit), cores = cores)
  output$waic <- loo::waic(loo::extract_log_lik(output$stanfit), cores = cores)
  output$call <- Call
  output$call$approach <- approach_flag
  output$call$model <- model_flag
  class(output) <- "spbp"
  return(output)
}
