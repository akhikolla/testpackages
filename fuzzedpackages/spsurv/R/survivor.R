#---------------------------------------------
#' Generic S3 method vcov
#'
#' @aliases survivor
#' @export
#' @param spbp an object of class spbp
#' @param ... further arguments passed to or from other methods
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{itsamp}}
#' @return estimates survival for each dataset individual (line).

#' Spbp Object Observed Survival
survivor <- function(spbp, ...) UseMethod("survivor")

#' @aliases survivor.spbp
#' @rdname survivor-methods
#' @description A method to allow survivor estimates for a spbp fit
#' @method survivor spbp
#' @export
#' @export survivor
#' @param spbp an object of the class spbp
#' @param newdata set of features referring to a specific population
#' @return Returns the probabilities that a subject will survive beyond any given times
#' @export

survivor.spbp <- function(spbp, newdata, ...){
  design <- model.matrix(spbp)
  if(missing(newdata)){
    newdata <- data.frame(t(matrix(colMeans(design))))
    colnames(newdata) <- colnames(design)
  }
  else{
    if(!is.data.frame(newdata))
      stop("newdata is not a data.frame object")
    newdata <- model.matrix(as.formula(paste("~", as.character(spbp$call$formula)[3])),
                            data =  newdata)[, -1]
  }
  res <- list()
  for(N in 1:nrow(newdata)){

    if(spbp$call$approach == "bayes"){

      beta <- rstan::extract(spbp$stanfit, pars = "beta")$beta
      if(ncol(newdata) != ncol(beta))
        stop("cols must match with `model.matrix(spbp)`")
      gamma <- rstan::extract(spbp$stanfit, pars = "gamma")$gamma
      iter <- nrow(beta)
      ####

      s <- matrix(NA, ncol = length(spbp$y[,1]), nrow = iter)
      for(i in 1:iter){
        s[i, ] <- survivor.aux(time = spbp$y[,1],
                                   arg = list(beta = beta[i,],
                                              gamma = gamma[i,]),
                                   newdata = matrix(newdata[N,], nrow = 1),
                                   model = spbp$call$model,
                                   approach = spbp$call$approach)
      }
      s <- colMeans(s)
    }
    else{
      beta <- spbp$coefficients[1:spbp$q]
      if(ncol(newdata) != length(beta))
        stop("cols must match with `model.matrix(spbp)`")
      gamma <- spbp$coefficients[(spbp$q+1):length(spbp$coefficients)]
      ####

      s <- survivor.aux(time = spbp$y[,1],
                            arg = list(beta = beta,
                                       gamma = gamma),
                            newdata = newdata[N,],
                            model = spbp$call$model,
                            approach = spbp$call$approach)
      }
    res[[N]] <- s
  }
  res <- cbind(sort(spbp$y[,1]), do.call(cbind, res))
  colnames(res) <- c("time", paste0("survival", 1:(ncol(res)-1)))
  return(data.frame(res))
}

survivor.aux <- function(time,
                         arg = list(beta = NULL, gamma = NULL),
                         newdata,
                         model = c("ph", "po", "aft"),
                         approach = c("mle", "bayes"), ...){
  e <- parent.frame()
  assign("design", get("design", envir = e))

  if(nrow(newdata)>1){
    res <- apply(newdata, 1, survivor.calc,
                 time = time, arg = arg,
                 model = model, approach = approach)
  }
  else{
    res <- survivor.calc(time = time, arg = arg,
                         newdata = as.numeric(newdata), model = model,
                         approach = approach)
  }
  return(res)
}


survivor.calc <- function(time,
                          arg = list(beta = NULL, gamma = NULL),
                          newdata,
                          model = c("ph", "po", "aft"),
                          approach = c("mle", "bayes")){
  e <- parent.frame()
  assign("design", get("design", envir = e))

  if(sum(names(arg) %in% c("beta", "gamma")) != 2)
    stop('`args` names do not match')

  ## CALL EXCEPTION HANDLING
  approach <- match.arg(approach)
  model <- match.arg(model)
  beta <- arg$beta
  gamma <- arg$gamma

  if(!is.vector(time))
    stop("time is not a vector")

  if(!is.vector(gamma))
    stop("gamma is not a vector")

  if(!is.vector(beta))
    stop("beta is not a vector")

  # if(!is.data.frame(newdata))
  #   stop("newdata is not a data.frame")

  x <- newdata
  degree <- length(gamma)
  k <- 1:degree
  y <- time[order(time)]
  tau <- max(y)
  B <- matrix(sapply(k, function(k) pbeta(y/tau, k, degree - k + 1)), ncol = degree)

  eta <- x %*% matrix(beta, ncol = 1)

  if(model == "ph"){
    H0 <- apply(B, 1, function(x){gamma %*% x})
    H <- as.vector(exp(eta)) * H0
  }
  else if(model == "po"){
    R0 <- apply(B, 1, function(x){gamma  %*% x})
    R <- as.vector(exp(eta)) * R0
    H <- log(1 + R)
  }
  else{
    y_aft <- as.matrix(y) / as.vector(exp(eta))
    tau_aft <- max(as.matrix(y)/exp(design %*% beta))
    B <- matrix(sapply(k, function(k) pbeta(y_aft / tau_aft, k, degree - k + 1)), ncol = degree)
    H <- apply(B, 1, function(x){gamma %*% x})
  }
  return(exp(-H))
}

#' @export
#' @method residuals spbp
#' @title BP based models residuals.
#' @description Residuals for a fitted \code{\link[spsurv]{spbp}} model.
#' @param object an object of class `spbp` result of a \code{\link[spsurv]{spbp}} fit.
#' @param type type of residuals, default is "cox-snell"
#' @param ... further arguments passed to or from other methods
#' @seealso \code{\link[spsurv]{spbp}}.
#' @examples
#'
#' library("spsurv")
#' data("veteran")
#'
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
#' data = veteran)
#'
#' residuals(fit)
#'

residuals.spbp <- function(object, type=c("cox-snell"), ...){
  return(-log(survivor(object)))
}
