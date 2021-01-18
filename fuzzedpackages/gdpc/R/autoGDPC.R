auto.gdpc <- function(Z, crit = "LOO", normalize = 1, auto_comp = TRUE, expl_var = 0.9, num_comp = 5, tol = 1e-04, 
                      k_max = 10, niter_max = 500, ncores = 1, verbose = FALSE) {
  # Computes Generalized Dynamic Principal Components. The number of components can be supplied by the user 
  # or chosen automatically so that a given proportion of variance is explained. The number of lags is chosen
  # automatically using of the following criteria: LOO, AIC, BIC or BNG.
  # All computations are done internally using leads rather than lags. However, the final result is outputted
  # using lags.
  #INPUT
  # Z: data matrix, series by columns 
  # crit: a string: 'LOO', 'AIC', 'BIC' or 'BNG', indicating the criterion used to chose the number of lags. Default is 'LOO'
  # normalize: integer. If normalize = 1, the data is analyzed in the orig-
  # inal units, without mean and variance standarization. If normalize = 2, the data
  # is standardized to zero mean and unit variance before computing the principal
  # components, but the intercepts and the loadings are those needed to reconstruct
  # the original series. If normalize = 3 the data are standardized as in normalize = 2
  # but the intercepts and the loadings are those needed to reconstruct the standardized
  # series. The default is normalize = 1.
  # auto_comp:  logical, if TRUE compute components until the proportion of explained variance is equal to expl_var, other
  # wise use num_comp components. Default is TRUE
  # expl_var: a number between 0 and 1. Desired proportion of explained variance (only if auto_comp==TRUE).
  # default is 0.9
  # num_comp: integer, number of components to be computed (only if auto_comp==FALSE). Default is 5
  # tol: desired accuracy when computing the components. Default is 1e-4
  # k_max: maximum number of lags. Default is 10
  # niter_max : maximum number of iterations. Default is 500
  # ncores: number of cores to be used for parallel computations. Default is 1
  # verbose : logical, print progress?
  
  #OUTPUT
  # A list of length equal to the number of computed components. The i-th entry of this list is an object of class
  # gdpc, that is, a list with entries:
  # f: Coordinates of the i-th Principal Component corresponding to the periods 1,…,T
  # initial_f: Coordinates of the i-th Principal Component corresponding to the periods -k+1,…,0.
  # Only for the case k>0, else 0.
  # beta: beta matrix corresponding to f
  # alpha: alpha vector corresponding to f
  # mse: mean (in T and m) squared error of the residuals of the fit with the first i components 
  # k: number of lags used chosen using the criterion specified in crit
  # crit: the LOO, AIC, BIC, or BNG of the fitted model, according to what was specified in crit
  # expart: proportion of the variance explained by the first i components
  # call: the matched call
  # conv: Logical. Did the iterations converge?
  
  if (all(!inherits(Z, "matrix"), !inherits(Z, "mts"), !inherits(Z, "xts"), !inherits(Z, "zoo"), !inherits(Z, "data.frame"))) {
    stop("Z should belong to one of the following classes: matrix, data.frame, mts, xts, zoo")
  } else if (any(dim(Z)[2] < 2, dim(Z)[1] < 10)) {
    stop("Z should have at least ten rows and two columns")
  } else if (any(anyNA(Z), any(is.nan(Z)), any(is.infinite(Z)))) {
    stop("Z should not have missing, infinite or nan values")
  }
  if (!crit %in% c("BNG", "LOO", "BIC", "AIC")) {
    stop("crit should be LOO, AIC, BIC or BNG")
  }
  if (!inherits(auto_comp, "logical")) {
    stop("auto_comp should be logical")
  }
  if (all(!inherits(normalize, "numeric"), !inherits(normalize, "integer"))) {
    stop("normalize should be numeric")
  } else if (any(!normalize == floor(normalize), normalize <= 0, normalize > 3)) {
    stop("normalize should be either 1, 2 or 3")
  }
  if (!inherits(expl_var, "numeric")) {
    stop("expl_var should be numeric")
  } else if (!all(expl_var < 1, expl_var > 0)) {
    stop("expl_var be between 0 and 1")
  }
  if (!inherits(tol, "numeric")) {
    stop("tol should be numeric")
  } else if (!all(tol < 1, tol > 0)) {
    stop("tol be between 0 and 1")
  }
  if (!inherits(num_comp, "numeric")) {
    stop("num_comp should be numeric")
  } else if (any(!num_comp == floor(num_comp), num_comp <= 0)) {
    stop("num_comp should be a positive integer")
  }
  if (!inherits(k_max, "numeric")) {
    stop("k_max should be numeric")
  } else if (any(!k_max == floor(k_max), k_max < 0)) {
    stop("k_max should be a non-negative integer")
  }
  if (!inherits(niter_max, "numeric")) {
    stop("niter_max should be numeric")
  } else if (any(!niter_max == floor(niter_max), niter_max <= 0)) {
    stop("niter_max should be a positive integer")
  }
  if (!inherits(ncores, "numeric")) {
    stop("ncores should be numeric")
  } else if (any(!ncores == floor(ncores), ncores <= 0)) {
    stop("ncores should be a positive integer")
  }
  # Pass to matrix form. Scale and transpose data.
  if (normalize == 2 | normalize == 3) {
    V <- t(scale(as.matrix(Z)))
    mean_var_V <- 1
  } else {
    V <- t(as.matrix(Z))
    mean_var_V <- mean(apply(V, 1, var))  # Mean variance of the data
  }
  vard <- (1 - expl_var) * mean_var_V
  output <- vector("list")
  
  sel <- switch(crit, LOO = 1, AIC = 2, BIC = 3, BNG = 4)
  
  ### Set-up cluster for parallel computations
  cores <- min(detectCores(), ncores)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  ### 
  
  comp_ready <- 1
  if (verbose) {
    cat(paste("Computing component number", comp_ready, "\n"))
  }
  out <- getLeads(V, k_max, mean_var_V, tol, niter_max, sel)
  mse <- out$mse  #Mean squared error (in N and m)
  output[[comp_ready]] <- out
  V <- out$res
  
  if (auto_comp) {
    while (mse > vard) {
      comp_ready <- comp_ready + 1
      if (verbose){
        cat(paste("Computing component number", comp_ready, "\n"))
      }
      out <- getLeads(V, k_max, mean_var_V, tol, niter_max, sel)
      mse <- out$mse
      output[[comp_ready]] <- out
      V <- out$res
    }
  } else {
    while (comp_ready < num_comp) {
      comp_ready <- comp_ready + 1
      if (verbose){
        cat(paste("Computing component number", comp_ready, "\n"))
      }
      out <- getLeads(V, k_max, mean_var_V, tol, niter_max, sel)
      output[[comp_ready]] <- out
      V <- out$res
    }
  }
  
  on.exit(stopCluster(cl))
  
  fn_call <- match.call()
  fn_call$crit <- crit
  output <- construct.gdpcs(output, Z, fn_call, normalize)
  if (verbose){
    cat(paste("Total number of computed components:", comp_ready, "\n"))
  }
  
  return(output)
  
}

getLeads <- function(V, k_max, mean_var_V, tol = 1e-04, niter_max = 500, sel = 1) {
  #Auxiliary function to choose the optimal number of leads
  #INPUT
  # V : matrix of original data or residuals where each ROW is a different time series
  # k_max : maximum of numbers of leads to be considered
  # mean_var_V : mean variance of original data
  # tol : relative precision
  # niter_max: maximum number of iterations
  # sel: criterion to be used, LOO = 1, AIC = 2, BIC = 3, BNG = 4
  #OUTPUT
  # A list with entries:
  # k: optimal number of leads
  # f: dynamic component with k leads
  # beta: matrix of loadings and intercept corresponding to f. Last column are the
  # intercepts (alpha).
  # mse: mean squared error
  # crit: criterion
  # res: matrix of residuals
  # expart: proportion of the variance explained
  # conv: Logical. Did the iterations converge?
  
  
  m <- nrow(V)
  N <- ncol(V)
  
  exports <- c("gdpc_priv")
  k_lag <- NULL
  crits <- rep(0, k_max + 1)
  fits <- vector("list", k_max + 1)
  fits <- foreach(k_lag = 1:(k_max + 1), .export = exports, .packages = "gdpc") %dopar% {
    gdpc_priv(V, k = k_lag - 1, f_ini = 0, passf_ini = FALSE, tol = tol, niter_max = niter_max, sel = sel)
  }
  
  crits <- sapply(fits, function(x) { x$crit })  #Get criterion corresponding to each lead
  convs <- sapply(fits, function(x) { x$conv })
  if (!all(convs)) {
    warning("Iterations did not converge. Consider increasing niter_max.")
  }
  k_opt <- which.min(crits) - 1
  out <- fits[[k_opt + 1]]
  expart <- 1 - out$mse / mean_var_V
  out$expart <- expart
  return(out)
}

gdpc <- function(Z, k, f_ini = NULL, tol = 1e-04, niter_max = 500, crit = "LOO") {
  # A wrapper function for gdpc_priv.
  #INPUT
  # Z: data matrix each COLUMN is a different time series
  # k: number of lags used
  # f_ini: starting point for the iterations. Optional. If no argument is passed
  # the standard Principal Component completed with k leads is used.
  # tol: relative precision, stopping criterion
  # niter_max: maximum number of iterations
  # first principal component with k lags is used
  # crit: a string: "LOO", "AIC", "BIC" or "BNG"
  #OUTPUT
  # An object of class gdpc, that is, a list with entries:
  # f: coordinates of the Principal Component corresponding to the periods 1,…,T
  # initial_f: Coordinates of the Principal Component corresponding to the periods -k+1,…,0.
  # beta: beta matrix of loadings corresponding to f
  # alpha: alpha vector of intercepts corresponding to f
  # mse: mean (in T and m) squared error of the residuals of the fit
  # k: number of lags used
  # crit: the criterion of the fitted model, according to what was specified in crit
  # expart: proportion of the variance explained
  # call: the matched call
  # conv: logical. Did the iterations converge?
  
  
  if (all(!inherits(Z, "matrix"), !inherits(Z, "mts"), !inherits(Z, "xts"), !inherits(Z, "zoo"), !inherits(Z, "data.frame"))) {
    stop("Z should belong to one of the following classes: matrix, data.frame, mts, xts, zoo")
  } else if (any(dim(Z)[2] < 2, dim(Z)[1] < 10)) {
    stop("Z should have at least ten rows and two columns")
  } else if (any(anyNA(Z), any(is.nan(Z)), any(is.infinite(Z)))) {
    stop("Z should not have missing, infinite or nan values")
  }
  if (!crit %in% c("BNG", "LOO","BIC", "AIC")) {
    stop("crit should be LOO, AIC, BIC or BNG")
  }
  if (!inherits(tol, "numeric")) {
    stop("tol should be numeric")
  } else if (!all(tol < 1, tol > 0)) {
    stop("tol be between 0 and 1")
  }
  if (!inherits(k, "numeric")) {
    stop("k should be numeric")
  } else if (any(!k == floor(k), k < 0)) {
    stop("k should be a non-negative integer")
  }
  if (!inherits(niter_max, "numeric")) {
    stop("niter_max should be numeric")
  } else if (any(!niter_max == floor(niter_max), niter_max <= 0)) {
    stop("niter_max should be a positive integer")
  }
  if (!is.null(f_ini)) {
    if (all(!inherits(f_ini, "numeric"), !inherits(f_ini, "ts"), !inherits(f_ini, "xts"))) {
      stop("f_ini should belong to one of the following classes: numeric, ts, xts")
    } else if (length(f_ini) != dim(Z)[1] + k) {
      stop("f_ini should have length equal to T + k")
    }
  }
  
  
  sel <- switch(crit, LOO = 1, AIC = 2, BIC = 3, BNG = 4)
  if (is.null(f_ini)) { 
    out <- gdpc_priv(t(Z), k, 0, FALSE, tol, niter_max, sel)
  } else {
    out <- gdpc_priv(t(Z), k, f_ini, TRUE, tol, niter_max, sel)
  }
  out$expart <- 1 - out$mse/mean(apply(Z, 2, var))
  fn_call <- match.call()
  fn_call$crit <- crit
  out$call <- fn_call
  out <- construct.gdpc(out, Z)
  return(out)
}