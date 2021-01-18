##################################################################
# gtreg.sp() function                                            #
##################################################################
# This is the gtreg() function written by Boan Zhang

gtreg.sp <- function(formula, data, groupn, retest = NULL, sens = 1, 
                     spec = 1, linkf = c("logit", "probit", "cloglog"),
                     method = c("Vansteelandt", "Xie"), sens.ind = NULL, 
                     spec.ind = NULL, start = NULL, 
                     control = gtRegControl(...), ...) {
  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "groupn"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  gr <- model.extract(mf, "groupn")
  if (!is.na(pos <- match(deparse(substitute(retest)), names(data))))
    retest <- data[, pos]
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf)
  else matrix(, NROW(Y), 0)
  linkf <- match.arg(linkf)
  if ((method <- match.arg(method)) == "Vansteelandt") {
    if (!is.null(retest))
      warning("Retests cannot be used with Vansteelandt's method.")
    fit <- gtreg.fit(Y, X, gr, sens, spec, linkf, start)
  }
  else {
    if (is.null(retest)) 
      fit <- EM(Y, X, gr, sens, spec, linkf, start, control)
    else fit <-  EM.ret(Y, X, gr, retest, sens, spec, linkf,
                        sens.ind, spec.ind, start, control)
  }
  fit <- c(fit, list(call = call, formula = formula, method = method,
                     link = linkf, terms = mt))
  class(fit) <- "gt"
  fit
  
}




##################################################################
# gtreg.halving() function                                       #
##################################################################

gtreg.halving <- function(formula, data, groupn, subg, retest, sens = 1, 
                          spec = 1, linkf = c("logit", "probit", "cloglog"),
                          sens.ind = NULL, spec.ind = NULL, start = NULL, 
                          control = gtRegControl(...), ...) {
  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "groupn"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  gr <- model.extract(mf, "groupn")
  if (!is.na(pos <- match(deparse(substitute(retest)), names(data))))
    retest <- data[, pos]
  if (!is.na(pos <- match(deparse(substitute(subg)), names(data))))
    subg <- data[, pos]
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf)
  else matrix(, NROW(Y), 0)
  linkf <- match.arg(linkf)
  fit <-  EM.halving(Y, X, gr, subg, retest, sens, spec, linkf,
                     sens.ind, spec.ind, start, control)
  fit <- c(fit, list(call = call, formula = formula, method = "Xie",
                     link = linkf, terms = mt))
  class(fit) <- "gt"
  fit
  
}




##################################################################
# gtreg.mp() function                                            #
##################################################################
# "mp" refers to matrix pooling, another name for array testing

gtreg.mp <- function(formula, data, coln, rown, arrayn, retest = NULL,
                     sens = 1, spec = 1, 
                     linkf = c("logit", "probit", "cloglog"), 
                     sens.ind = NULL, spec.ind = NULL, start = NULL, 
                     control = gtRegControl(...), ...){
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "coln", "rown",
               "arrayn"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  arrayn <- model.extract(mf, "arrayn")
  rown <- model.extract(mf, "rown")
  coln <- model.extract(mf, "coln")
  if (!is.na(pos <- match(deparse(substitute(retest)), names(data))))
    retest <- data[, pos]
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf)
  else matrix(, NROW(Y), 0)
  linkf <- match.arg(linkf)
  fit <- EM.mp(Y[, 1], Y[, 2], X, coln, rown, arrayn, retest,
               sens, spec, linkf, sens.ind, spec.ind, start, control)
  fit <- c(fit, list(call = call, formula = formula, link = linkf,
                     terms = mt))
  class(fit) <- c("gt.mp", "gt")
  fit
}




# Supporting functions
##################################################################
# gtreg.fit() function                                           #
##################################################################

gtreg.fit <- function (Y, X, groupn, sens, spec, linkf, start = NULL){
  z <- tapply(Y, groupn, tail, n = 1)
  num.g <- max(groupn)
  K <- ncol(X)
  sam <- length(Y)
  if (is.null(start)) {
    if (K == 1) {
      cova.mean <- as.matrix(tapply(X, groupn, mean)) 
      optim.meth <- "BFGS"
    } 
    else {
      temp <- by(X, groupn, colMeans)
      cova.mean <- do.call(rbind, temp)
      optim.meth <- "Nelder-Mead"
    }
    beta.group <- glm.fit(cova.mean, as.vector(z), 
                          family = binomial(link = linkf))$coefficients
  } 
  else {
    beta.group <- start
    names(beta.group) <- dimnames(X)[[2]]
    optim.meth <- ifelse(K == 1, "BFGS", "Nelder-Mead")
  }   
  logL <- function(beta) {
    pijk <- switch(linkf, logit = plogis(X %*% beta), probit = pnorm(X %*%
                                                                       beta), cloglog = 1 - exp(-exp(X %*% beta)))
    prodp <- tapply(1 - pijk, groupn, prod)
    -sum(z * log(sens + (1 - sens - spec) * prodp) +
           (1 - z) * log(1 - sens - (1 - sens - spec) * prodp))
  }
  mod.fit <- optim(par = beta.group, fn = logL, method = optim.meth,
                   control = list(trace = 0, maxit = 1000), hessian = TRUE)
  if (det(mod.fit$hessian) == 0)
    mod.fit <- optim(par = beta.group, fn = logL, method = "SANN", hessian = TRUE)
  logL0 <- function(beta) {
    inter <- rep(beta, sam)
    pijk <- switch(linkf, logit = plogis(inter), probit = pnorm(inter),
                   cloglog = 1 - exp(-exp(inter)))
    prodp <- tapply(1 - pijk, groupn, prod)
    -sum(z * log(sens + (1 - sens - spec) * prodp) +
           (1 - z) * log(1 - sens - (1 - sens - spec) * prodp))
  }
  mod.fit0 <- optim(par = binomial()$linkfun(mean(z)), 
                    fn = logL0, method = "BFGS", control = list(trace = 0, maxit = 1000))
  nulld <- 2 * mod.fit0$value
  residd <- 2 * mod.fit$value
  xib <- X %*% mod.fit$par
  pijk <- switch(linkf, logit = plogis(xib), probit = pnorm(xib),
                 cloglog = 1 - exp(-exp(xib)))
  prodp <- tapply(1 - pijk, groupn, prod)
  zhat <- sens + (1 - sens - spec) * prodp
  residual <- z - zhat
  aic <- residd + 2 * K
  if (mod.fit$convergence == 0)
    counts <- mod.fit$counts[[1]]
  else warning("Maximum number of iterations exceeded.")
  list(coefficients = mod.fit$par, hessian = mod.fit$hessian,
       fitted.values = zhat, deviance = residd, df.residual = num.g - K,
       null.deviance = nulld, df.null = num.g - 1, aic = aic, counts = counts,
       residuals = residual, z = z)
}




##################################################################
# gtRegControl() function                                        #
##################################################################

#' @title Auxiliary for controlling group testing regression
#' 
#' @description Auxiliary function to control fitting parameters 
#' of the EM algorithm used internally in \code{\link{gtReg}} 
#' for simple pooling (\kbd{type="sp"}) with \kbd{method="Xie"} 
#' or for array testing (\kbd{type="array"}). 
#' 
#' @param tol convergence criterion
#' @param n.gibbs the Gibbs sample size to be used in each E step 
#' of the EM algorithm, for array testing. The default is 1000.
#' @param n.burnin the number of samples in the burn-in period, 
#' for array testing. The default is 20.
#' @param maxit maximum number of iterations in the EM algorithm.
#' @param trace a logical value indicating whether the output should 
#' be printed for each iteration. The default is \kbd{FALSE}. 
#' @param time a logical value indicating whether the length of time 
#' for the model fitting should be printed. The default is \kbd{TRUE}.
#'   
#' @return A list with components named as the input arguments.
#' 
#' @author This function was originally written as the \code{gt.control} 
#' function for the binGroup package. Minor modifications have been 
#' made for inclusion in the binGroup2 package.
#' 
#' @examples 
#' # The default settings:
#' gtRegControl()

gtRegControl <- function (tol = 0.0001, n.gibbs = 1000, n.burnin = 20, 
                        maxit = 500, trace = FALSE, time = TRUE){
  if (!is.numeric(tol) || tol <= 0) 
    stop("value of 'tol' must be > 0")
  if (round(n.gibbs) != n.gibbs || n.gibbs <= 0) 
    stop("value of 'n.gibbs' must be a positive integer")
  if (round(n.burnin) != n.burnin || n.burnin <= 0) 
    stop("value of 'n.burnin' must be a positive integer")    
  if (!is.numeric(maxit) || maxit <= 0) 
    stop("maximum number of iterations must be > 0")
  list(tol = tol, n.gibbs = n.gibbs, n.burnin = n.burnin, maxit = maxit, 
       trace = trace, time = time)
}




##################################################################
# EM() function                                                  #
##################################################################
EM <- function (Y, X, groupn, sens, spec, linkf, start = NULL, control = gtRegControl())
{
  if (control$time)
    start.time <- proc.time()
  z <- tapply(Y, groupn, tail, n = 1)
  num.g <- max(groupn)
  K <- ncol(X)
  if (is.null(start)) {
    if (K == 1)
      cova.mean <- as.matrix(tapply(X, groupn, mean))
    else {
      temp <- by(X, groupn, colMeans)
      cova.mean <- do.call(rbind, temp)
    }
    beta.old <- lm.fit(cova.mean, z)$coefficients
  } 
  else beta.old <- start
  sam <- length(Y)
  vec <- 1:sam
  group.sizes <- tapply(Y, groupn, length)
  diff <- 1
  counts <- 1
  extra.loop <- FALSE
  next.loop <- TRUE
  while (next.loop) {
    xib <- X %*% beta.old
    pijk <- switch(linkf, logit = plogis(xib),
                   probit = pnorm(xib), cloglog = 1 - exp(-exp(xib)))
    prodp <- tapply(1 - pijk, groupn, prod)
    den <- rep((1 - spec) * prodp + sens * (1 - prodp), group.sizes)
    den2 <- rep(spec * prodp + (1 - sens) * (1 - prodp),
                group.sizes)
    expect <- rep(NA, times = sam)
    for (i in vec) {
      if (Y[i] == 0)
        expect[i] <- (1 - sens) * pijk[i]/den2[i]
      else expect[i] <- sens * pijk[i]/den[i]
    }
    if (!extra.loop) {
      suppress <- function(w) 
        if(any(grepl("non-integer #successes in a binomial glm", w))) 
          invokeRestart("muffleWarning")
      mod.fit <- withCallingHandlers(glm.fit(X, expect, 
                                             family = binomial(link = linkf)), warning = suppress)
      diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
      beta.old <- mod.fit$coefficients
      if (control$trace)
        cat("beta is", beta.old, "\tdiff is", diff, "\n")
      counts <- counts + 1
      if (diff <= control$tol || counts > control$maxit) 
        extra.loop <- TRUE
    } 
    else next.loop <- FALSE
  }    
  erf <- 2 * pijk - 1
  pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
                probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
                                                                      erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
  pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
                                              erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
                                                                                                      pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
                cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
                              exp(xib))/(exp(-exp(xib)) - 1)^2)
  nm <- pt1 + expect * pt2
  sign1 <- as.vector(sign(nm))
  nn <- as.vector(sqrt(abs(nm)))
  x2 <- X * nn
  m <- (t(x2) %*% (sign1 * x2))
  b <- array(NA, c(K, K, sum(group.sizes^2)))
  p <- 1
  for (i in vec) for (j in vec[groupn == groupn[i]]) {
    wii <- ifelse(i == j, expect[i] - expect[i]^2, expect[i] *
                    (pijk[j] - expect[j]))
    coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i]^2 +
                                                         xib[j]^2)/2)/((1 - erf[i]^2) * (1 - erf[j]^2) * pi),
                  cloglog = exp(xib[i] + xib[j])/((exp(-exp(xib[i])) -
                                                     1) * (exp(-exp(xib[j])) - 1)))
    b[, , p] <- wii * coe * X[i, ] %*% t(X[j, ])
    p <- p + 1
  }
  m1 <- apply(b, c(1, 2), sum)
  H <- -(m + m1)
  zhat <- sens + (1 - sens - spec) * prodp
  residual <- z - zhat
  residd <- -2 * sum(z * log(zhat) + (1 - z) * log(1 - zhat))
  logL0 <- function(beta) {
    inter <- rep(beta, sam)
    pijk <- switch(linkf, logit = plogis(inter), probit = pnorm(inter),
                   cloglog = 1 - exp(-exp(inter)))
    prodp <- tapply(1 - pijk, groupn, prod)
    -sum(z * log(sens + (1 - sens - spec) * prodp) +
           (1 - z) * log(1 - sens - (1 - sens - spec) * prodp))
  }
  mod.fit0 <- optim(par = binomial()$linkfun(mean(z)), fn = logL0,
                    method = "BFGS", control = list(trace = 0, maxit = 1000))
  nulld <- 2 * mod.fit0$value
  aic <- residd + 2 * K
  if (diff > control$tol && counts > control$maxit)
    warning("EM algorithm did not converge.")
  if (control$time) {
    end.time <- proc.time()
    save.time <- end.time - start.time
    cat("\n Number of minutes running:", round(save.time[3]/60, 2), "\n \n")
  }
  list(coefficients = beta.old, hessian = H, fitted.values = zhat,
       deviance = residd, df.residual = num.g - K, null.deviance = nulld,
       df.null = num.g - 1, aic = aic, counts = counts - 1, residuals = residual,
       z = z)
}




##################################################################
# EM.ret() function                                              #
##################################################################

EM.ret <- function (Y, X, groupn, ret, sens, spec, linkf,
                    sens.ind, spec.ind,
                    start = NULL, control = gtRegControl())
{
  if (control$time)
    start.time <- proc.time()
  if (is.null(sens.ind))
    sens.ind <- sens
  if (is.null(spec.ind))
    spec.ind <- spec
  z <- tapply(Y, groupn, tail, n = 1)
  num.g <- max(groupn)
  K <- ncol(X)
  if (is.null(start)) {
    if (K == 1)
      cova.mean <- as.matrix(tapply(X, groupn, mean))
    else {
      temp <- by(X, groupn, colMeans)
      cova.mean <- do.call(rbind, temp)
    }
    beta.old <- lm.fit(cova.mean, z)$coefficients
  }
  else beta.old <- start
  sam <- length(Y)
  vec <- 1:sam
  group.sizes <- tapply(Y, groupn, length)
  diff <- 1
  counts <- 1
  extra.loop <- FALSE
  next.loop <- TRUE
  a0 <- ifelse(ret == 1, sens.ind, 1 - sens.ind)
  a1 <- ifelse(ret == 0, spec.ind, 1 - spec.ind)
  while (next.loop) {
    xib <- X %*% beta.old
    pijk <- switch(linkf, logit = plogis(xib),
                   probit = pnorm(xib), cloglog = 1 -
                     exp(-exp(xib)))
    erf <- 2 * pijk - 1
    prodp <- tapply(1 - pijk, groupn, prod)
    den2 <- rep(spec * prodp + (1 - sens) * (1 - prodp),
                group.sizes)
    expect <- rep(NA, times = sam)
    i <- 1
    while (i <= sam) {
      if (Y[i] == 0)
        expect[i] <- (1 - sens) * pijk[i]/den2[i]
      else {
        vec1 <- vec[groupn == groupn[i]]
        mb2 <- 1
        for (l in vec1) {
          temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
          mb2 <- mb2 * temp
        }
        null <- 1
        for (l in vec1) {
          temp <- a1[l] * (1 - pijk[l])
          null <- null * temp
        }                                                       
        den <- mb2 * sens + null * (1 - sens - spec)
        for (l1 in vec1) {                         
          temp <- a0[l1] * pijk[l1] + a1[l1] * (1 - pijk[l1])
          num <- mb2/temp * a0[l1] * pijk[l1] * sens
          expect[l1] <- num/den
        }
        i <- l1
      }
      i <- i + 1
    }
    expect[expect > 1] <- 1
    expect[expect < 0] <- 0
    if (!extra.loop) {
      suppress <- function(w)
        if (any(grepl("non-integer #successes in a binomial glm", w)))
          invokeRestart("muffleWarning")
      mod.fit <- withCallingHandlers(glm.fit(X, expect,
                                             family = binomial(link = linkf)), warning = suppress)
      diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
      beta.old <- mod.fit$coefficients
      if (control$trace)
        cat("beta is", beta.old, "\tdiff is", diff, "\n")
      counts <- counts + 1
      if (diff <= control$tol || counts > control$maxit) 
        extra.loop <- TRUE
    } 
    else next.loop <- FALSE
  }
  pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
                probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
                                                                      erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
  pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
                                              erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
                                                                                                      pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
                cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
                              exp(xib))/(exp(-exp(xib)) - 1)^2)
  nm <- pt1 + expect * pt2
  sign1 <- as.vector(sign(nm))
  nn <- as.vector(sqrt(abs(nm)))
  x2 <- X * nn
  m <- (t(x2) %*% (sign1 * x2))
  m1 <- 0    
  for (i in vec) {
    vec1 <- vec[groupn == groupn[i]]
    if (Y[i] == 0) {
      for (j in vec1) {
        coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i]^2 +
                                                             xib[j]^2)/2)/((1 - erf[i]^2) * (1 - erf[j]^2) * pi),
                      cloglog = exp(xib[i] + xib[j])/((exp(-exp(xib[i])) -
                                                         1) * (exp(-exp(xib[j])) - 1)))            
        wii <- ifelse(i == j, expect[i] - expect[i]^2, expect[i] *
                        (pijk[j] - expect[j]))
        tim <- wii * coe * X[i, ] %*% t(X[j, ])
        m1 <- m1 + tim
      }                
    }        
    else {         
      for (j in vec1) {      
        temp <- a0[j] * pijk[j] + a1[j] * (1 - pijk[j])
        eii <- expect[i]/temp * a0[j] * pijk[j]
        wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
        coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i]^2 +
                                                             xib[j]^2)/2)/((1 - erf[i]^2) * (1 - erf[j]^2) * pi),
                      cloglog = exp(xib[i] + xib[j])/((exp(-exp(xib[i])) -
                                                         1) * (exp(-exp(xib[j])) - 1)))
        tim <- wii * coe * X[i, ] %*% t(X[j, ])
        m1 <- m1 + tim
      }
    }
  }
  H <- -(m + m1)
  zhat <- sens + (1 - sens - spec) * prodp
  residual <- z - zhat
  logl <- 0
  for (grn in 1:num.g) {
    if (z[grn] == 1) {
      vec1 <- vec[groupn == grn]
      mb2 <- 1
      for (l in vec1) {
        temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
        mb2 <- mb2 * temp
      }
      null <- 1
      for (l in vec1) {
        temp <- a1[l] * (1 - pijk[l])
        null <- null * temp
      }                                                       
      prob1 <- mb2 * sens + null * (1 - sens - spec)
    } else prob1 <- 1 - zhat[grn]
    logl <- logl - log(prob1)
  }
  aic <- 2 * logl + 2 * K
  if (diff > control$tol && counts > control$maxit)
    warning("EM algorithm did not converge.")
  if (control$time) {
    end.time <- proc.time()
    save.time <- end.time - start.time
    cat("\n Number of minutes running:", round(save.time[3]/60, 2), "\n \n")
  }
  list(coefficients = beta.old, hessian = H, fitted.values = zhat,
       deviance = 2 * logl, aic = aic, counts = counts - 1, residuals = residual, 
       z = z)
}




##################################################################
# EM.mp() function                                               #
##################################################################

EM.mp <- function (col.resp, row.resp, X, coln, rown, sqn, ret, sens,
                   spec, linkf, sens.ind, spec.ind, start = NULL, control = gtRegControl())
{
  if (control$time)
    start.time <- proc.time()
  if (is.null(sens.ind))
    sens.ind <- sens
  if (is.null(spec.ind))
    spec.ind <- spec
  len <- max(sqn)
  diff <- 1
  counts <- 1
  sam <- length(sqn)
  col.groupn <- coln[sqn == 1]
  if (len > 1) {
    for (i in 2:len) {
      temp <- max(col.groupn) + coln[sqn == i]
      col.groupn <- c(col.groupn, temp)
    }
  }
  if (is.null(start)) {
    mod.fit <- try(gtreg.fit(col.resp, X, col.groupn,
                             sens, spec, linkf))
    if (class(mod.fit) == "try-error") {
      row.groupn <- rown[sqn == 1]
      if (len > 1) {
        for (i in 2:len) {
          temp <- max(row.groupn) + rown[sqn == i]
          row.groupn <- c(row.groupn, temp)
        }
      }
      mod.fit <- gtreg.fit(row.resp, X, row.groupn,
                           sens, spec, linkf)
    }
    beta.old <- mod.fit$coefficients
  }
  else beta.old <- start
  extra.loop <- FALSE
  next.loop <- TRUE
  while (next.loop) {
    xib <- X %*% beta.old
    pijk.all <- switch(linkf, logit = plogis(xib),
                       probit = pnorm(xib), cloglog = 1 - exp(-exp(xib)))
    expect.all <- numeric(0)
    mat2 <- index <- 0        
    erf <- 2 * pijk.all - 1
    for (arrayn in 1:len) {
      index.r <- index.c <- vector("logical", length = sam)
      for (i in 1:sam) {
        if (rown[i] == 1 && sqn[i] == arrayn)
          index.c[i] <- TRUE
        else index.c[i] <- FALSE
        if (coln[i] == 1 && sqn[i] == arrayn)
          index.r[i] <- TRUE
        else index.r[i] <- FALSE
      }
      n.row <- max(rown[index.r])
      n.col <- max(coln[index.c])
      rowresp <- row.resp[index.r]
      colresp <- col.resp[index.c]
      index <- max(index) + 1:(n.row * n.col)
      if (!is.null(ret)) {
        re.ind <- na.omit(cbind(coln[sqn == arrayn],
                                rown[sqn == arrayn], ret[sqn == arrayn]))
        re <- ifelse(re.ind[, 3] == 1, sens.ind, 1 -
                       sens.ind)
        re1 <- ifelse(re.ind[, 3] == 0, spec.ind, 1 -
                        spec.ind)
      }
      pijk <- matrix(pijk.all[sqn == arrayn], nrow = n.row)
      a <- ifelse(rowresp == 1, sens, 1 - sens)
      b <- ifelse(colresp == 1, sens, 1 - sens)
      a1 <- ifelse(rowresp == 0, spec, 1 - spec)
      b1 <- ifelse(colresp == 0, spec, 1 - spec)
      mat <- array(NA, c(n.row, n.col, control$n.gibbs))
      y <- matrix(0, nrow = n.row, ncol = n.col)
      for (k in 1:(control$n.gibbs + control$n.burnin)) {
        l <- 1
        for (j in 1:n.col) for (i in 1:n.row) {
          num <- a[i] * b[j] * pijk[i, j]
          den.r <- ifelse(sum(y[i, ]) - y[i, j] > 0,
                          a[i], a1[i])
          den.c <- ifelse(sum(y[, j]) - y[i, j] > 0,
                          b[j], b1[j])
          den2 <- den.r * den.c * (1 - pijk[i, j])
          if (!is.null(ret)) {
            if (l <= length(re) && j == re.ind[l, 1] &&
                i == re.ind[l, 2]) {
              num <- num * re[l]
              den2 <- den2 * re1[l]
              l <- l + 1
            }
          }
          den <- num + den2
          if (den != 0) {
            cond.p <- num/den
            y[i, j] <- rbinom(1, 1, cond.p)
          }
          else y[i, j] <- 0
        }
        if (k > control$n.burnin) {
          mat[, , k - control$n.burnin] <- y
          vec <- as.vector(y)
          if (extra.loop)
            for (i1 in index[vec == 1]) for (j1 in index[vec ==
                                                         1]) {
              bq <- switch(linkf, logit = 1, probit = 8 *
                             exp(-(xib[i1]^2 + xib[j1]^2)/2)/((1 - erf[i1]^2) *
                                                                (1 - erf[j1]^2) * pi), cloglog = exp(xib[i1] +
                                                                                                       xib[j1])/((exp(-exp(xib[i1])) - 1) * (exp(-exp(xib[j1])) -
                                                                                                                                               1))) * X[i1, ] %*% t(X[j1, ])
              mat2 <- mat2 + bq
            }
        }
      }
      expect.m <- apply(mat, c(1, 2), mean)
      expect <- as.vector(expect.m)
      expect.all <- c(expect.all, expect)
    }
    if (!extra.loop) {
      suppress <- function(w) 
        if(any(grepl("non-integer #successes in a binomial glm", w))) 
          invokeRestart("muffleWarning")
      mod.fit <- withCallingHandlers(glm.fit(X, expect.all, 
                                             family = binomial(link = linkf)), warning = suppress)
      diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
      beta.old <- mod.fit$coefficients
      if (control$trace)
        cat("beta is", beta.old, "\tdiff is", diff, "\n")
      counts <- counts + 1
      if (diff <= control$tol || counts > control$maxit) 
        extra.loop <- TRUE
    } 
    else next.loop <- FALSE
  }
  index <- 0
  first <- mat2/control$n.gibbs
  second <- 0
  for (arrayn in 1:len) {
    n.row <- max(rown[sqn == arrayn])
    n.col <- max(coln[sqn == arrayn])
    index <- max(index) + 1:(n.row * n.col)
    expect <- expect.all[index]
    for (i1 in index) for (j1 in index) {
      coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i1]^2 +
                                                           xib[j1]^2)/2)/((1 - erf[i1]^2) * (1 - erf[j1]^2) *
                                                                            pi), cloglog = exp(xib[i1] + xib[j1])/((exp(-exp(xib[i1])) -
                                                                                                                      1) * (exp(-exp(xib[j1])) - 1)))
      tim <- expect.all[i1] * expect.all[j1] * coe * X[i1,
                                                       ] %*% t(X[j1, ])
      second <- second + tim
    }
  }
  m1 <- first - second
  pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
                probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
                                                                      erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
  pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
                                              erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
                                                                                                      pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
                cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
                              exp(xib))/(exp(-exp(xib)) - 1)^2)
  nm <- pt1 + expect.all * pt2
  sign1 <- as.vector(sign(nm))
  nn <- as.vector(sqrt(abs(nm)))
  x2 <- X * nn
  m <- (t(x2) %*% (sign1 * x2))
  H <- -(m + m1)
  if (diff > control$tol && counts > control$maxit)
    warning("EM algorithm did not converge.")
  if (control$time) {
    end.time <- proc.time()
    save.time <- end.time - start.time
    cat("\n Number of minutes running:", round(save.time[3]/60, 2), "\n \n")
  }
  list(coefficients = beta.old, hessian = H, Gibbs.sample.size = control$n.gibbs,
       counts = counts - 1)
}



##################################################################
# EM.halving() function                                          #
##################################################################

EM.halving <- function (Y, X, groupn, subg, ret, sens, spec, linkf,
                        sens.ind, spec.ind,
                        start = NULL, control = gtRegControl())
{
  if (control$time)
    start.time <- proc.time()
  if (is.null(sens.ind))
    sens.ind <- sens
  if (is.null(spec.ind))
    spec.ind <- spec
  z <- tapply(Y, groupn, tail, n = 1)
  num.g <- max(groupn)
  K <- ncol(X)
  if (is.null(start)) {
    if (K == 1)
      cova.mean <- as.matrix(tapply(X, groupn, mean))
    else {
      temp <- by(X, groupn, colMeans)
      cova.mean <- do.call(rbind, temp)
    }
    beta.old <- lm.fit(cova.mean, z)$coefficients
  } else beta.old <- start
  sam <- length(Y)
  vec <- 1:sam
  group.sizes <- tapply(Y, groupn, length)
  diff <- 1
  counts <- 1
  extra.loop <- FALSE
  next.loop <- TRUE
  a0 <- ifelse(ret == 1, sens.ind, 1 - sens.ind)
  a1 <- ifelse(ret == 0, spec.ind, 1 - spec.ind)
  while (next.loop) {
    xib <- X %*% beta.old
    pijk <- switch(linkf, logit = plogis(xib),
                   probit = pnorm(xib), cloglog = 1 -
                     exp(-exp(xib)))
    erf <- 2 * pijk - 1
    prodp <- tapply(1 - pijk, groupn, prod)
    den2 <- rep(spec * prodp + (1 - sens) * (1 - prodp),
                group.sizes)
    expect <- rep(NA, times = sam)
    i <- 1
    while (i <= sam) {
      if (Y[i] == 0)
        expect[i] <- (1 - sens) * pijk[i]/den2[i]
      else {
        if (subg[i] == 0) {
          vec1 <- vec[groupn == groupn[i]]
          gs <- length(vec1)
          sub1 <- vec1[1:ceiling(gs/2)]
          sub2 <- vec1[(ceiling(gs/2) + 1):gs]
          if (subg[vec1[gs]] == 0) {
            den <- (1-spec)*spec^2*prod(1-pijk[sub1])*prod(1-pijk[sub2])+
              spec*(1-sens)*sens*prod(1-pijk[sub1])*(1-prod(1-pijk[sub2]))+
              spec*(1-sens)*sens*prod(1-pijk[sub2])*(1-prod(1-pijk[sub1]))+
              (1-sens)^2*sens*(1-prod(1-pijk[sub1]))*(1-prod(1-pijk[sub2]))
            ab1 <- (1-sens)*sens*(spec*prod(1-pijk[sub2])+
                                    (1-sens)*(1-prod(1-pijk[sub2])))
            ab2 <- (1-sens)*sens*(spec*prod(1-pijk[sub1])+
                                    (1-sens)*(1-prod(1-pijk[sub1])))
            for (l1 in sub1) {
              expect[l1]<-ab1*pijk[l1]/den
            }
            for (l1 in sub2) {
              expect[l1]<-ab2*pijk[l1]/den
            }
          }
          if (subg[vec1[gs]] == 1) {
            mb2 <- 1
            for (l in sub2) {
              temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
              mb2 <- mb2 * temp
            }
            null <- 1
            for (l in sub2) {
              temp <- a1[l] * (1 - pijk[l])
              null <- null * temp
            }
            den <- (1-spec)^2*spec*null*prod(1-pijk[sub1])+
              (1-spec)*(1-sens)*sens*null*(1-prod(1-pijk[sub1]))+
              spec*sens^2*(mb2-null)*prod(1-pijk[sub1])+
              (1-sens)*sens^2*(mb2-null)*(1-prod(1-pijk[sub1]))
            ab1 <- (1-sens)*sens*(mb2*sens+null*(1-sens-spec))
            for (l1 in sub1) {
              expect[l1]<-ab1*pijk[l1]/den
            }
            for (l1 in sub2) {
              temp <- a0[l1] * pijk[l1] + a1[l1] * (1 - pijk[l1])
              num <- mb2/temp * a0[l1] * pijk[l1] * sens^2*(spec*prod(1-pijk[sub1])+
                                                              (1-sens)*(1-prod(1-pijk[sub1])))
              expect[l1]<-num/den
            }
          }
          i <- l1
        } else {
          vec1 <- vec[groupn == groupn[i]]
          gs <- length(vec1)
          sub1 <- vec1[1:ceiling(gs/2)]
          sub2 <- vec1[(ceiling(gs/2) + 1):gs]
          if (subg[vec1[gs]] == 0) {
            mb2 <- 1
            for (l in sub1) {
              temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
              mb2 <- mb2 * temp
            }
            null <- 1
            for (l in sub1) {
              temp <- a1[l] * (1 - pijk[l])
              null <- null * temp
            }
            den <- (1-spec)^2*spec*null*prod(1-pijk[sub2])+
              (1-spec)*(1-sens)*sens*null*(1-prod(1-pijk[sub2]))+
              spec*sens^2*(mb2-null)*prod(1-pijk[sub2])+
              (1-sens)*sens^2*(mb2-null)*(1-prod(1-pijk[sub2]))
            ab1 <- (1-sens)*sens*(mb2*sens+null*(1-sens-spec))
            for (l1 in sub1) {
              temp <- a0[l1]*pijk[l1]+a1[l1]*(1-pijk[l1])
              num <- mb2/temp*a0[l1]*pijk[l1]*sens^2*(spec*prod(1-pijk[sub2])+
                                                        (1-sens)*(1-prod(1-pijk[sub2])))
              expect[l1]<-num/den
              
            }
            for (l1 in sub2) {
              expect[l1]<-ab1*pijk[l1]/den
            }
          }
          if (subg[vec1[gs]] == 1) {
            mb2 <- 1
            for (l in sub1) {
              temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
              mb2 <- mb2 * temp
            }
            null <- 1
            for (l in sub1) {
              temp <- a1[l] * (1 - pijk[l])
              null <- null * temp
            }
            mb2a <- 1
            for (l in sub2) {
              temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
              mb2a <- mb2a * temp
            }
            nulla <- 1
            for (l in sub2) {
              temp <- a1[l] * (1 - pijk[l])
              nulla <- nulla * temp
            }
            den <- (1-spec)^3*null*nulla+
              (1-spec)*sens^2*null*(mb2a-nulla)+
              (1-spec)*sens^2*(mb2-null)*nulla+
              sens^3*(mb2-null)*(mb2a-nulla)
            for (l1 in sub1) {
              temp <- a0[l1]*pijk[l1]+a1[l1]*(1-pijk[l1])
              num <- mb2/temp*a0[l1]*pijk[l1]*sens^2*(mb2a*sens+nulla*(1-sens-spec))
              expect[l1]<-num/den
            }
            for (l1 in sub2) {
              temp <- a0[l1]*pijk[l1]+a1[l1]*(1-pijk[l1])
              num <- mb2a/temp*a0[l1]*pijk[l1]*sens^2*(mb2*sens+null*(1-sens-spec))
              expect[l1]<-num/den
            }
          }
          i <- l1
        }
      }
      i <- i + 1
    }
    expect[expect > 1] <- 1
    expect[expect < 0] <- 0
    if (!extra.loop) {
      suppress <- function(w)
        if (any(grepl("non-integer #successes in a binomial glm", w)))
          invokeRestart("muffleWarning")
      mod.fit <- withCallingHandlers(glm.fit(X, expect,
                                             family = binomial(link = linkf)), warning = suppress)
      diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
      beta.old <- mod.fit$coefficients
      if (control$trace)
        cat("beta is", beta.old, "\tdiff is", diff, "\n")
      counts <- counts + 1
      if (diff <= control$tol || counts > control$maxit)
        extra.loop <- TRUE
    }
    else next.loop <- FALSE
  }
  pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
                probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
                                                                      erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
  pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
                                              erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
                                                                                                      pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
                cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
                              exp(xib))/(exp(-exp(xib)) - 1)^2)
  nm <- pt1 + expect * pt2
  sign1 <- as.vector(sign(nm))
  nn <- as.vector(sqrt(abs(nm)))
  x2 <- X * nn
  m <- (t(x2) %*% (sign1 * x2))
  m1 <- 0
  i <- 1
  while (i <= sam) {
    vec1 <- vec[groupn == groupn[i]]
    gs <- length(vec1)
    if (Y[i] == 0) {
      for (j in vec1) {
        wii <- ifelse(i == j, expect[i] - expect[i]^2, expect[i] *
                        (pijk[j] - expect[j]))
        tim <- wii * X[i, ] %*% t(X[j, ])
        m1 <- m1 + tim
      }
    } else {
      sub1 <- vec1[1:ceiling(gs/2)]
      sub2 <- vec1[(ceiling(gs/2) + 1):gs]
      for (i in sub1) {
        for (j in sub1) { 
          if (subg[j] == 0) { 
            eii <- expect[i] * pijk[j]
          } else { 
            temp <- a0[j] * pijk[j] + a1[j] * (1 - pijk[j])
            eii <- expect[i]/temp * a0[j] * pijk[j]
          }
          wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
          tim <- wii * X[i, ] %*% t(X[j, ])
          m1 <- m1 + tim
        }
        for (j in sub2) { 
          if (subg[j] == 0) { 
            temp<-spec*prod(1-pijk[sub2])+(1-sens)*(1-prod(1-pijk[sub2]))
            eii <- expect[i]*(1-sens)*pijk[j]/temp
          } else { 
            mb2a <- 1
            for (l in sub2) {
              temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
              mb2a <- mb2a * temp
            }
            nulla <- 1
            for (l in sub2) {
              temp <- a1[l] * (1 - pijk[l])
              nulla <- nulla * temp
            }
            temp <- a0[j]*pijk[j]+a1[j]*(1-pijk[j])
            tempa <- mb2a * sens + nulla * (1 - sens - spec)
            eii <- expect[i]/tempa*sens*a0[j]*pijk[j]*mb2a/temp
          }
          wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
          tim <- wii * X[i, ] %*% t(X[j, ])
          m1 <- m1 + tim
        }
      }
      for (i in sub2) {
        for (j in sub1) { 
          if (subg[j] == 0) { 
            temp<-spec*prod(1-pijk[sub1])+(1-sens)*(1-prod(1-pijk[sub1]))
            eii <- expect[i] * (1-sens)* pijk[j]/temp
          } else { 
            mb2 <- 1
            for (l in sub1) {
              temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
              mb2 <- mb2 * temp
            }
            null <- 1
            for (l in sub1) {
              temp <- a1[l] * (1 - pijk[l])
              null <- null * temp
            }
            temp <- a0[j]*pijk[j]+a1[j]*(1-pijk[j])
            tempa <- mb2*sens+null*(1-sens-spec)
            eii <- expect[i]/tempa*sens*a0[j]*pijk[j]*mb2/temp
          }
          wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
          tim <- wii * X[i, ] %*% t(X[j, ])
          m1 <- m1 + tim
        }
        for (j in sub2) { 
          if (subg[j] == 0) { 
            eii <- expect[i] * pijk[j]
          } else { 
            temp <- a0[j] * pijk[j] + a1[j] * (1 - pijk[j])
            eii <- expect[i]/temp * a0[j] * pijk[j]
          }
          wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
          tim <- wii * X[i, ] %*% t(X[j, ])
          m1 <- m1 + tim
        }
      }
    }
    i <- i + 1
  }
  H <- -(m + m1)
  zhat <- sens + (1 - sens - spec) * prodp
  residual <- z - zhat
  logl <- 0
  for (grn in 1:num.g) {
    if (z[grn] == 1) {
      vec1 <- vec[groupn == grn]
      gs <- length(vec1)
      sub1 <- vec1[1:ceiling(gs/2)]
      sub2 <- vec1[(ceiling(gs/2) + 1):gs]
      if (subg[vec1[1]] == 0) {
        if (subg[vec1[gs]] == 0) {
          prob1 <- (1-spec)*spec^2*prod(1-pijk[sub1])*prod(1-pijk[sub2])+
            spec*(1-sens)*sens*prod(1-pijk[sub1])*(1-prod(1-pijk[sub2]))+
            spec*(1-sens)*sens*prod(1-pijk[sub2])*(1-prod(1-pijk[sub1]))+
            (1-sens)^2*sens*(1-prod(1-pijk[sub1]))*(1-prod(1-pijk[sub2]))
        }
        if (subg[vec1[gs]] == 1) {
          mb2 <- 1
          for (l in sub2) {
            temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
            mb2 <- mb2 * temp
          }
          null <- 1
          for (l in sub2) {
            temp <- a1[l] * (1 - pijk[l])
            null <- null * temp
          }
          prob1 <- (1-spec)^2*spec*null*prod(1-pijk[sub1])+
            (1-spec)*(1-sens)*sens*null*(1-prod(1-pijk[sub1]))+
            spec*sens^2*(mb2-null)*prod(1-pijk[sub1])+
            (1-sens)*sens^2*(mb2-null)*(1-prod(1-pijk[sub1]))
        }
      } else {
        if (subg[vec1[gs]] == 0) {
          mb2 <- 1
          for (l in sub1) {
            temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
            mb2 <- mb2 * temp
          }
          null <- 1
          for (l in sub1) {
            temp <- a1[l] * (1 - pijk[l])
            null <- null * temp
          }
          prob1 <- (1-spec)^2*spec*null*prod(1-pijk[sub2])+
            (1-spec)*(1-sens)*sens*null*(1-prod(1-pijk[sub2]))+
            spec*sens^2*(mb2-null)*prod(1-pijk[sub2])+
            (1-sens)*sens^2*(mb2-null)*(1-prod(1-pijk[sub2]))
        }
        if (subg[vec1[gs]] == 1) {
          mb2 <- 1
          for (l in sub1) {
            temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
            mb2 <- mb2 * temp
          }
          null <- 1
          for (l in sub1) {
            temp <- a1[l] * (1 - pijk[l])
            null <- null * temp
          }
          mb2a <- 1
          for (l in sub2) {
            temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
            mb2a <- mb2a * temp
          }
          nulla <- 1
          for (l in sub2) {
            temp <- a1[l] * (1 - pijk[l])
            nulla <- nulla * temp
          }
          prob1 <- (1-spec)^3*null*nulla+
            (1-spec)*sens^2*null*(mb2a-nulla)+
            (1-spec)*sens^2*(mb2-null)*nulla+
            sens^3*(mb2-null)*(mb2a-nulla)
        }
      }
    } else prob1 <- 1 - zhat[grn]
    logl <- logl - log(prob1)
  }
  aic <- 2 * logl + 2 * K
  if (diff > control$tol && counts > control$maxit)
    warning("EM algorithm did not converge.")
  if (control$time) {
    end.time <- proc.time()
    save.time <- end.time - start.time
    cat("\n Number of minutes running:", round(save.time[3]/60, 2), "\n \n")
  }
  list(coefficients = beta.old, hessian = H, fitted.values = zhat,
       deviance = 2 * logl, aic = aic,
       counts = counts - 1, residuals = residual, z = z)
}




##################################################################
# gtReg() function                                               #
##################################################################

#' @title Fitting group testing regression models
#' 
#' @description Fits the group testing regression model specified 
#' through a symbolic description of the linear predictor and 
#' descriptions of the group testing setting. This function allows 
#' for fitting regression models with simple pooling, halving, or array 
#' testing data. 
#' 
#' @param type \kbd{"sp"} for simple pooling, \kbd{"halving"} for halving 
#' protocol, or \kbd{"array"} for array testing. See 'Details' for 
#' descriptions of the group testing algorithms.
#' @param formula an object of class "formula" (or one that 
#' can be coerced to that class); a symbolic description of 
#' the model to be fitted. The details of model specification 
#' are under 'Details'.
#' @param data an optional data frame, list, or environment 
#' (or object coercible by \kbd{as.data.frame} to a data frame) 
#' containing the variables in the model. If not found in data, 
#' the variables are taken from \kbd{environment(formula)}, 
#' typically the environment from which \code{gtReg} is called.
#' @param groupn a vector, list, or data frame of the group 
#' numbers that designates individuals to groups (for use with 
#' simple pooling, \kbd{type="sp"}, or the halving protocol, 
#' \kbd{type="halving"}).
#' @param subg a vector, list, or data frame of the group numbers 
#' that designates individuals to subgroups (for use with the 
#' halving protocol, \kbd{type="halving"}).
#' @param coln a vector, list, or data frame that specifies the
#' column group number for each sample (for use with array 
#' testing, \kbd{type="array"}).
#' @param rown a vector, list, or data frame that specifies the
#' row group number for each sample (for use with array testing, 
#' \kbd{type="array"}).
#' @param arrayn a vector, list, or data frame that specifies the
#' array number for each sample (for use with array testing, 
#' \kbd{type="array"}).
#' @param retest a vector, list, or data frame of individual 
#' retest results. Default value is \kbd{NULL} for no retests. 
#' See 'Details' for details on how to specify \kbd{retest}.
#' @param sens sensitivity of the test. Default value is set 
#' to 1.
#' @param spec specificity of the test. Default value is set 
#' to 1.
#' @param linkf a character string specifying one of the three 
#' link functions for a binomial model: \kbd{"logit"} (default), 
#' \kbd{"probit"}, or \kbd{"cloglog"}.
#' @param method the method to fit the regression model. 
#' Options include \kbd{"Vansteelandt"} (default) or \kbd{"Xie"}. 
#' The \kbd{"Vansteelandt"} option finds estimates by directly 
#' maximizing the likelihood function based on the group responses, 
#' while the \kbd{"Xie"} option uses the EM algorithm to 
#' maximize the likelihood function in terms of the unobserved 
#' individual responses.
#' @param sens.ind sensitivity of the individual retests. If NULL, 
#' set to be equal to sens.
#' @param spec.ind specificity of the individual retests. If NULL, 
#' set to be equal to spec. 
#' @param start starting values for the parameters in the linear 
#' predictor.
#' @param control a list of parameters for controlling the fitting
#' process in method \kbd{"Xie"}. These parameters will be passed 
#' to the \code{\link{gtRegControl}} function for use.
#' @param ... arguments to be passed to \code{\link{gtRegControl}} by 
#' default. See argument \kbd{control}. 
#' 
#' @details With simple pooling and halving, a typical predictor 
#' has the form \kbd{groupresp ~ covariates} where \kbd{groupresp} 
#' is the (numeric) group response vector. With array testing, 
#' individual samples are placed in a matrix-like grid where 
#' samples are pooled within each row and within each column. 
#' This leads to two kinds of group responses: row and column 
#' group responses. Thus, a typical predictor has the form 
#' \kbd{cbind(col.resp, row.resp) ~ covariates}, where 
#' \kbd{col.resp} is the (numeric) column group response vector 
#' and \kbd{row.resp} is the (numeric) row group response vector. 
#' For all methods, \kbd{covariates} is a series of terms which 
#' specifies a linear predictor for individual responses.
#' Note that it is actually the unobserved individual responses, 
#' not the observed group responses, which are modeled by the 
#' covariates. When denoting group responses (\kbd{groupresp}, 
#' \kbd{col.resp}, and \kbd{row.resp}), a 0 denotes a negative 
#' response and a 1 denotes a positive response, where the 
#' probability of an individual positive response is being 
#' modeled directly. 
#' 
#' A terms specification of the form  
#' \kbd{first + second} indicates all the terms in \kbd{first} 
#' together with all the terms in \kbd{second} with duplicates 
#' removed. A specification of the form \kbd{first:second} 
#' indicates the set of terms obtained by taking the interactions 
#' of all terms in \kbd{first} with all terms in \kbd{second}. 
#' The specification \kbd{first*second} indicates the cross of 
#' \kbd{first} and \kbd{second}. This is the same as \kbd{first + 
#' second + first:second}. The terms in the formula will be 
#' re-ordered so that main effects come first, followed by the 
#' interactions, all second-order, all third-order, and so on; 
#' to avoid this, pass a terms object as the formula.
#' 
#' For simple pooling (\kbd{type="sp"}), the functions \kbd{gtreg.fit}, 
#' \kbd{EM}, and \kbd{EM.ret}, where the first corresponds to Vansteelandt's 
#' method described in Vansteelandt et al. (2000) and the last two correspond 
#' to Xie's method described in Xie (2001), are called to carry out the 
#' model fitting. The \kbd{gtreg.fit} function uses the \kbd{optim} 
#' function with default method \kbd{"Nelder-Mead"} to maximize
#' the likelihood function of the observed group responses. 
#' If this optimization method produces a Hessian matrix of all 
#' zero elements, the \kbd{"SANN"} method in \kbd{optim} is 
#' employed to find the coefficients and Hessian matrix. For 
#' the \kbd{"SANN"} method, the number of iterations in \kbd{optim} 
#' is set to be 10000. For the background on the use of \kbd{optim}, 
#' see \kbd{help(optim)}.
#' 
#' The \kbd{EM} and \kbd{EM.ret} functions apply Xie's EM 
#' algorithm to the likelihood function written in terms of the 
#' unobserved individual responses; the functions use \kbd{glm.fit} 
#' to update the parameter estimates within each M step. The 
#' \kbd{EM} function is used when there are no retests and 
#' \kbd{EM.ret} is used when individual retests are available. 
#' Thus, within the \kbd{retest} argument, individual observations 
#' in observed positive groups are 0 (negative) or 1 (positive); 
#' the remaining individual observations are \kbd{NA}s, meaning 
#' that no retest is performed for them. Retests cannot be used 
#' with Vansteelandt's method; a warning message will be given 
#' in this case, and the individual retests will be ignored in 
#' the model fitting. There could be slight differences in the 
#' estimates between Vansteelandt's and Xie's methods (when 
#' retests are not available) due to different convergence criteria.
#' 
#' With simple pooling (i.e., Dorfman testing, two-stage hierarchical 
#' testing), each individual appears in exactly one pool. When only the 
#' group responses are observed, the null degrees of freedom are the number 
#' of groups minus 1 and the residual degrees of freedom are the number of 
#' groups minus the number of parameters. When individual retests are 
#' observed too, it is an open research question for what the degrees of 
#' freedom and the deviance for the null model should be; therefore, the 
#' degrees of freedom and \kbd{null.deviance} will not be displayed.
#' 
#' Under the halving protocol, the \kbd{EM.halving} function 
#' applies Xie's EM algorithm to the 
#' likelihood function written in terms of the unobserved 
#' individual responses; the functions use \kbd{glm.fit} to update 
#' the parameter estimates within each M step. In the halving 
#' protocol, if the initial group tests positive, it is split 
#' into two subgroups. The two subgroups are subsequently tested 
#' and if either subgroup tests positive, the third and final 
#' step is to test all individuals within the subgroup. Thus, 
#' within \kbd{subg}, subgroup responses in observed positive 
#' groups are 0 (negative) or 1 (positive); the remaining 
#' subgroup responses are \kbd{NA}s, meaning that no tests are 
#' performed for them. The individual retests are similarly coded.
#' 
#' With array testing (also known as matrix pooling), the 
#' \kbd{EM.mp} function applies Xie's 
#' EM algorithm to the likelihood function written in terms of the 
#' unobserved individual responses. In each E step, the Gibbs 
#' sampling technique is used to estimate the conditional 
#' probabilities. Because of the large number of Gibbs samples 
#' needed to achieve convergence, the model fitting process could 
#' be quite slow, especially when multiple positive rows and 
#' columns are observed. In this case, we can either increase the 
#' Gibbs sample size to help achieve convergence or loosen the 
#' convergence criteria by increasing \kbd{tol} at the expense 
#' of perhaps poorer estimates. If follow-up retests are performed, 
#' the retest results going into the model will help achieve 
#' convergence faster with the same Gibbs sample size and 
#' convergence criteria. In each M step, we use \kbd{glm.fit} to 
#' update the parameter estimates. 
#' 
#' For simple pooling, \kbd{retest} provides individual retest 
#' results for Dorfman's retesting procedure. Under the halving 
#' protocol, \kbd{retest} provides individual retest results 
#' within a subgroup that tests positive. The \kbd{retest} 
#' argument provides individual retest results, where a 0 
#' denotes negative and 1 denotes positive status. A \kbd{NA} 
#' denotes that no retest is performed for that individual. 
#' The default value is \kbd{NULL} for no retests.
#' 
#' For simple pooling, \kbd{control} provides parameters for 
#' controlling the fitting process in the \kbd{"Xie"} method only.
#' 
#' \kbd{gtReg} returns an object of class \kbd{"gtReg"}. 
#' The function \kbd{summary} (i.e., \code{\link{summary.gtReg}} 
#' is used to obtain or print a summary of the results. 
#' The group testing function \kbd{predict} (i.e., 
#' \code{\link{predict.gtReg}}) is used to make predictions 
#' on \kbd{"gtReg"} objects.
#' 
#' @return An object of class \kbd{"gtReg"}, a list which may include:
#' \item{coefficients}{a named vector of coefficients.}
#' \item{hessian}{estimated Hessian matrix of the negative 
#' log-likelihood function. This serves as an estimate of the 
#' information matrix.}
#' \item{residuals}{the response residuals. This is the difference 
#' of the observed group responses and the fitted group 
#' responses. Not included for array testing.}
#' \item{fitted.values}{the fitted mean values of group responses.
#' Not included for array testing.}
#' \item{deviance}{the deviance between the fitted model and the 
#' saturated model. Not included for array testing.}
#' \item{aic}{Akaike's Information Criterion. This is minus twice 
#' the maximized log-likelihood plus twice the number of 
#' coefficients. Not included for array testing.}
#' \item{null.deviance}{the deviance for the null model,
#' comparable with \kbd{deviance}. The null model will 
#' include only the intercept, if there is one in the model.
#' Provided for simple pooling, \kbd{type="sp"}, only.}
#' \item{counts}{the number of iterations in \kbd{optim} 
#' (Vansteelandt's method) or the number of iterations in the 
#' EM algorithm (Xie's method, halving, and array testing).}
#' \item{Gibbs.sample.size}{the number of Gibbs samples 
#' generated in each E step. Provided for array testing, 
#' \kbd{type="array"}, only.}
#' \item{df.residual}{the residual degrees of freedom.
#' Provided for simple pooling, \kbd{type="sp"}, only.}
#' \item{df.null}{the residual degrees of freedom for the null model.
#' Provided for simple pooling, \kbd{type="sp"}, only.}
#' \item{z}{the vector of group responses. Not included for array testing.}
#' \item{call}{the matched call.}
#' \item{formula}{the formula supplied.}
#' \item{terms}{the terms object used.}
#' \item{method}{the method (\kbd{"Vansteelandt"} or \kbd{"Xie"}) 
#' used to fit the model. For the halving protocol, the 
#' \kbd{"Xie"} method is used. Not included for array testing.}
#' \item{link}{the link function used in the model.}
#' 
#' @author The majority of this function was originally written as 
#' \kbd{gtreg.sp}, \kbd{gtreg.halving}, and \kbd{gtreg.mp} by Boan Zhang 
#' for the \code{binGroup} package. Minor modifications have been made for 
#' inclusion of the functions in the \code{binGroup2} package.
#' 
#' @references 
#' \insertRef{Vansteelandt2000}{binGroup2}
#' 
#' \insertRef{Xie2001}{binGroup2}
#'
#' @seealso \code{\link{gtSim}} for simulation of data in the 
#' group testing form to be used by \kbd{gtReg}, 
#' \code{\link{summary.gtReg}} and \code{\link{predict.gtReg}} 
#' for \kbd{gtreg} methods. 
#' 
#' @examples 
#' # Estimated running time for all examples was calculated 
#' #   using a computer with 16 GB of RAM and one core of 
#' #   an Intel i7-6500U processor. Please take this into 
#' #   account when interpreting the run times given.
#' 
#' data(hivsurv)
#' fit1 <- gtReg(type="sp", formula = groupres ~ AGE + EDUC., 
#'               data = hivsurv, groupn = gnum, sens = 0.9, 
#'               spec = 0.9, method = "Xie")
#' fit1
#' 
#' set.seed(46)
#' gt.data <- gtSim(type="sp", par=c(-12, 0.2), 
#'                  size1=700, size2=5)
#' fit2 <- gtReg(type="sp", formula=gres~x, data=gt.data, 
#'               groupn=groupn)
#' fit2
#' 
#' set.seed(21)
#' gt.data <- gtSim(type="sp", par=c(-12, 0.2), 
#'                  size1=700, size2=6, sens=0.95, spec=0.95, 
#'                  sens.ind=0.98, spec.ind=0.98)
#' fit3 <- gtReg(type="sp", formula=gres~x, data=gt.data, 
#'               groupn=groupn, retest=retest, method="Xie", 
#'               sens=0.95, spec=0.95, sens.ind=0.98, 
#'               spec.ind=0.98, trace=TRUE)
#' summary(fit3)
#' 
#' set.seed(46)
#' gt.data <- gtSim(type="halving", par=c(-6, 0.1), gshape=17, 
#'                  gscale=1.4, size1=5000, size2=5, 
#'                  sens=0.95, spec=0.95)
#' fit4 <- gtReg(type="halving", formula=gres~x, 
#'               data=gt.data, groupn=groupn, subg=subgroup, 
#'               retest=retest, sens=0.95, spec=0.95, 
#'               start=c(-6, 0.1), trace=TRUE)
#' summary(fit4)
#' 
#' # This example takes approximately 15 seconds to run.
#' # 5x6 and 4x5 array
#' set.seed(9128)
#' sa1a <- gtSim(type="array", par=c(-7, 0.1), size1=c(5,4), 
#'               size2=c(6,5), sens=0.95, spec=0.95)
#' sa1 <- sa1a$dframe
#' \donttest{
#' fit5 <- gtReg(type="array", 
#'               formula=cbind(col.resp, row.resp)~x, 
#'               data=sa1, coln=coln, rown=rown, 
#'               arrayn=arrayn, sens=0.95, spec=0.95, 
#'               tol=0.005, n.gibbs=2000, trace=TRUE)
#' fit5
#' summary(fit5)}
#' 
#' # The example below shows how long this fitting process 
#' #   may take. It takes approximately 1.5 minutes to achieve 
#' #   convergence.
#' set.seed(9012)
#' sa2a <- gtSim(type="array", par=c(-7, 0.1), 
#'               size1=rep(10, 4), size2=rep(10, 4), 
#'               sens=0.95, spec=0.95)
#' sa2 <- sa2a$dframe
#' \donttest{
#' fit6 <- gtReg(type="array", 
#'               formula=cbind(col.resp, row.resp)~x, 
#'               data=sa2, coln=coln, rown=rown, 
#'               arrayn=arrayn, retest=retest, 
#'               sens=0.95, spec=0.95, 
#'               start=c(-7, 0.1), tol=0.005)
#' fit6
#' summary(fit6)}

# Brianna Hitt - 01-07-2020
gtReg <- function(type="sp", formula, data, groupn=NULL, 
                  subg=NULL, coln=NULL, rown=NULL, arrayn=NULL, 
                  retest=NULL, sens=1, spec=1, 
                  linkf=c("logit", "probit", "cloglog"),
                  method=c("Vansteelandt", "Xie"), 
                  sens.ind=NULL, spec.ind=NULL, start=NULL, 
                  control=gtRegControl(...), ...){
  
  call <- match.call()
  mf <- match.call(expand.dots=FALSE)
  
  if(type %in% c("sp", "halving")){
    m <- match(c("formula", "data", "groupn"), names(mf), 0)
  } else if(type=="array"){
    m <- match(c("formula", "data", "coln", "rown", "arrayn"), 
               names(mf), 0)
  }
  
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  if(type %in% c("sp", "halving")){
    gr <- model.extract(mf, "groupn")
  } else if(type=="array"){
    arrayn <- model.extract(mf, "arrayn")
    rown <- model.extract(mf, "rown")
    coln <- model.extract(mf, "coln")
  }
  
  if (!is.na(pos <- match(deparse(substitute(retest)), names(data))))
    retest <- data[, pos]
  if(type=="halving"){
    if (!is.na(pos <- match(deparse(substitute(subg)), names(data))))
      subg <- data[, pos]
  }
  
  Y <- model.response(mf, "any")
  if(length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm))
      names(Y) <- nm
  }
  X <- if(!is.empty.model(mt))
    model.matrix(mt, mf)
  else matrix(, NROW(Y), 0)
  linkf <- match.arg(linkf)
  
  if(type=="sp"){
    if ((method <- match.arg(method)) == "Vansteelandt") {
      if (!is.null(retest))
        warning("Retests cannot be used with Vansteelandt's method.")
      fit <- gtreg.fit(Y, X, gr, sens, spec, linkf, start)
    }
    else {
      if (is.null(retest)) 
        fit <- EM(Y, X, gr, sens, spec, linkf, start, control)
      else fit <-  EM.ret(Y, X, gr, retest, sens, spec, linkf,
                          sens.ind, spec.ind, start, control)
    }
    fit <- c(fit, list(call = call, formula = formula, method = method,
                       link = linkf, terms = mt))
    class(fit) <- "gt"
  } else if(type=="halving"){
    fit <-  EM.halving(Y, X, gr, subg, retest, sens, spec, linkf,
                       sens.ind, spec.ind, start, control)
    fit <- c(fit, list(call = call, formula = formula, method = "Xie",
                       link = linkf, terms = mt))
    class(fit) <- "gt"
  } else if(type=="array"){
    fit <- EM.mp(Y[, 1], Y[, 2], X, coln, rown, arrayn, retest,
                 sens, spec, linkf, sens.ind, spec.ind, start, control)
    fit <- c(fit, list(call = call, formula = formula, link = linkf,
                       terms = mt))
    class(fit) <- c("gt.mp", "gt")
  }
  
  class(fit) <- c(class(fit), "gtReg")
  fit
}




# Summary function
##################################################################
# summary.gtReg() function                                       #
##################################################################

#' @title Summary method for group testing regression model fits
#' 
#' @description Produce a summary list for objects of class 
#' \kbd{"gtReg"} returned by \code{\link{gtReg}}. 
#' 
#' @param object a fitted object of class \kbd{"gtReg"}.
#' @param ... currently not used.
#' 
#' @details The \kbd{coefficients} component of the results gives 
#' the estimated coefficients and their estimated standard errors, 
#' together with their ratio. This third column is labeled 
#' \kbd{z ratio} using Wald tests. A fourth column gives the 
#' two-tailed p-value corresponding to the z-ratio based on a 
#' Wald test. Note that it is possible that there are no residual 
#' degrees of freedom from which to estimate, in which case the 
#' estimate is \kbd{NaN}. 
#' 
#' @return \kbd{summary.gtReg} returns an object of class 
#' \kbd{"summary.gtReg"}, a list containing:
#' \item{call}{the component from \kbd{object}.}
#' \item{link}{the component from \kbd{object}.}
#' \item{deviance}{the component from \kbd{object}, 
#' for simple pooling (\kbd{type="sp"} in \code{\link{gtReg}}) only.}
#' \item{aic}{the component from \kbd{object}, 
#' for simple pooling (\kbd{type="sp"} in \code{\link{gtReg}}) only.}
#' \item{df.residual}{the component from \kbd{object}, 
#' for simple pooling (\kbd{type="sp"} in \code{\link{gtReg}}) only.}
#' \item{null.deviance}{the component from \kbd{object}, 
#' for simple pooling (\kbd{type="sp"} in \code{\link{gtReg}}) only.}
#' \item{df.null}{the component from \kbd{object}, 
#' for simple pooling (\kbd{type="sp"} in \code{\link{gtReg}}) only.}
#' \item{deviance.resid}{the deviance residuals, 
#' for simple pooling (\kbd{type="sp"} in \code{\link{gtReg}}) only.}
#' \item{coefficients}{the matrix of coefficients, standard errors, 
#' z-values, and p-values. Aliased coefficients are omitted.}
#' \item{counts}{the component from \kbd{object}.}
#' \item{method}{the component from \kbd{object}, 
#' for simple pooling (\kbd{type="sp"} in \code{\link{gtReg}}) only.}
#' \item{Gibbs.sample.size}{the component from \kbd{object}, 
#' for array testing (\kbd{type="array"} in \code{\link{gtReg}}) only.}
#' \item{cov.mat}{the estimated covariance matrix of the estimated 
#' coefficients.}
#' 
#' @author The majority of this function was originally written as 
#' \code{summary.gt} and \code{summary.gt.mp} by Boan Zhang for the 
#' \code{binGroup} package. Minor modifications were made to the function 
#' for inclusion in the \code{binGroup2} package.
#' 
#' @seealso \code{\link{gtReg}} for creating an object of class 
#' \kbd{"gtReg"}.
#' 
#' @examples 
#' data(hivsurv)
#' fit1 <- gtReg(type="sp", formula=groupres ~ AGE + EDUC., 
#'               data=hivsurv, groupn=gnum, sens=0.9, spec=0.9, 
#'               method="Xie")
#' summary(fit1)
#' 
#' # This examples takes approximately 5 seconds to run.
#' # 5x6 and 4x5 array
#' set.seed(9128)
#' sa2a <- gtSim(type="array", par=c(-7,0.1), size1=c(5,4), 
#'               size2=c(6,5), sens=0.95, spec=0.95)
#' sa2 <- sa2a$dframe
#' \donttest{
#' fit2 <- gtReg(type="array", formula=cbind(col.resp, row.resp) ~ x, 
#'               data=sa2, coln=coln, rown=rown, arrayn=arrayn, 
#'               sens=0.95, spec=0.95, linkf="logit", 
#'               n.gibbs=1000, tol=0.005)
#' summary(fit2)}

summary.gtReg <- function(object, ...)
{
  coef.p <- object$coefficients
  cov.mat <- solve(object$hessian)
  dimnames(cov.mat) <- list(names(coef.p), names(coef.p))
  var.cf <- diag(cov.mat)
  s.err <- sqrt(var.cf)
  zvalue <- coef.p/s.err
  dn <- c("Estimate", "Std. Error")
  pvalue <- 2 * pnorm(-abs(zvalue))
  coef.table <- cbind(coef.p, s.err, zvalue, pvalue)
  dimnames(coef.table) <- list(names(coef.p), c(dn, "z value",
                                                "Pr(>|z|)"))
  
  if ("gt.mp" %in% class(object)) {
    keep <- match(c("call", "link", "Gibbs.sample.size", "counts"), 
                  names(object), 0)
    ans <- c(object[keep], list(coefficients = coef.table, 
                                cov.mat = cov.mat))
  } else{
    keep <- match(c("call", "link", "aic", "deviance", "df.residual",
                    "null.deviance", "df.null", "counts", "method", "z"),
                  names(object), 0)
    ans <- c(object[keep], list(coefficients = coef.table, 
                                deviance.resid = residuals(object, type = "deviance"), 
                                cov.mat = cov.mat))
  }
  class(ans) <- "summary.gtReg"
  ans
}




# Print method for summary function
##################################################################
# print.summary.gtReg() function                                 #
##################################################################

#' @title Print method for \kbd{summary.gtReg}
#' 
#' @description Print method for objects obtained by 
#' \code{\link{summary.gtReg}}.
#' 
#' @param x An object of class "summary.gtReg" created by 
#' \code{\link{summary.gtReg}}.
#' @param digits digits for rounding.
#' @param signif.stars a logical value indicating whether significance 
#' stars should be shown.
#' @param ... Additional arguments to be passed to \code{printCoefmat}.
#' 
#' @return A print out of the function call, deviance residuals (for 
#' simple pooling and halving only), coefficients, null and 
#' residual deviance and degrees of freedom (for simple pooling only), 
#' AIC (for simple pooling and halving only), number of 
#' Gibbs samples (for array testing only), and the number of 
#' iterations.
#' 
#' @author This function combines code from 
#' \code{print.summary.gt} and 
#' \code{print.summary.gt.mp}, written by Boan Zhang 
#' for the \code{binGroup} package. Minor modifications were 
#' made for inclusion in the \code{binGroup2} package.

print.summary.gtReg <- 
  function(x, digits = max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...) {
    
    obj <- x
    cat("\nCall:\n")
    cat(paste(deparse(obj$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    
    # for simple pooling and array testing only
    if ("deviance.resid" %in% names(obj)) {
      cat("Deviance Residuals: \n")
      if (length(obj$z) > 5) {
        obj$deviance.resid <- quantile(obj$deviance.resid, na.rm = TRUE)
        names(obj$deviance.resid) <- c("Min", "1Q", "Median",
                                       "3Q", "Max")
      }
      print.default(obj$deviance.resid, digits = digits, na.print = "",
                    print.gap = 2)
    }
    
    cat("\nCoefficients:\n")
    coefs <- obj$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    
    if (!is.null(unlist(obj["df.null"]))) {
      cat("\n", 
          apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), 
                            "deviance:"), 
                      format(unlist(obj[c("null.deviance", "deviance")]), 
                             digits = 4), 
                      " on", format(unlist(obj[c("df.null", "df.residual")])), 
                      " degrees of freedom\n"), 1, paste,
                collapse = " "), sep = "")
    }

    if ("method" %in% names(obj)) {
      if (obj$method == "Vansteelandt") {
        cat("AIC: ", format(obj$aic, digits = 4), 
            "\n\n", "Number of iterations in optim(): ",
            obj$counts, "\n", sep = "")
      } else {
        cat("AIC: ", format(obj$aic, digits = 4), 
            "\n\n", "Number of iterations in EM: ",
            obj$counts, "\n", sep = "")
      }
    }

    if ("Gibbs.sample.size" %in% names(obj)) {
      cat("\nNumber of Gibbs samples generated in each E step: ", 
          obj$Gibbs.sample.size, "\n", 
          "Number of iterations in EM algorithm: ", 
          obj$counts, "\n", sep = "")
    }
    
    cat("\n")
    invisible(obj)
    
  }
    



# Other regression-related functions
##################################################################
# predict.gtReg() function                                       #
##################################################################

#' @title Predict method for group testing regression model fits
#' 
#' @description Obtains predictions for individual observations and 
#' optionally estimates standard errors of those predictions from 
#' objects of class \kbd{"gtReg"} returned by \code{\link{gtReg}}.
#' 
#' @param object a fitted object of class \kbd{"gtReg"}. 
#' @param newdata an optional data frame in which to look for 
#' variables with which to predict. If omitted, the fitted linear 
#' predictors are used.
#' @param type the type of prediction required. The \kbd{"link"} 
#' option is on the scale of the linear predictors. The \kbd{"response"} 
#' option is on the scale of the response variable. Thus, for the 
#' logit model, the \kbd{"link"} predictions are of log-odds 
#' (probabilities on the logit scale) and  \kbd{type="response"} 
#' gives the predicted probabilities.
#' @param se.fit a logical value indicating whether standard errors 
#' are required.
#' @param conf.level the confidence level of the interval for the 
#' predicted values.
#' @param na.action a function determining what should be done with 
#' missing values in \kbd{newdata}. The default is to predict \kbd{NA}. 
#' @param ... currently not used.
#' 
#' @details If \kbd{newdata} is omitted, the predictions are based 
#' on the data used for the fit. When \kbd{newdata} is present and 
#' contains missing values, how the missing values will be dealt with 
#' is determined by the \kbd{na.action} argument. In this case, if 
#' \kbd{na.action=na.omit}, omitted cases will not appear, whereas 
#' if \kbd{na.action=na.exclude}, omitted cases will appear (in 
#' predictions and standard errors) with value \kbd{NA}.
#' 
#' @return If \kbd{se=FALSE}, a vector or matrix of predictions. If 
#' \kbd{se=TRUE}, a list containing:
#' \item{fit}{predictions.}
#' \item{se.fit}{estimated standard errors.}
#' \item{lower}{the lower bound of the confidence interval, 
#' if calculated.}
#' \item{upper}{the upper bound of the confidence interval, 
#' if calculated.}
#' 
#' @author Boan Zhang
#' 
#' @examples 
#' data(hivsurv)
#' fit1 <- gtReg(formula = groupres ~ AGE + EDUC., data = hivsurv,  
#'               groupn = gnum, sens = 0.9, spec = 0.9, 
#'               linkf = "logit", method = "V")
#' pred.data <- data.frame(AGE = c(15, 25, 30), EDUC. = c(1, 3, 2))
#' predict(object = fit1, newdata = pred.data, type = "link", 
#'         se.fit = TRUE)
#' predict(object = fit1, newdata = pred.data, type = "response", 
#'         se.fit = TRUE, conf.level = 0.9)
#' predict(object = fit1, type = "response", se.fit = TRUE, 
#'         conf.level = 0.9)

predict.gtReg <- function (object, newdata, type = c("link", "response"), 
                        se.fit = FALSE, conf.level = NULL, 
                        na.action = na.pass, ...) 
{
  tt <- terms(object)
  Terms <- delete.response(tt)
  if (missing(newdata) || is.null(newdata)) {
    m <- model.frame(object)
    newd <- model.matrix(Terms, m)
  }
  else {
    m <- model.frame(Terms, newdata, na.action = na.action)
    newd <- model.matrix(Terms, m)
  }
  type <- match.arg(type)
  lin.pred <- as.vector(newd %*% object$coefficients)
  link <- object$link
  res <- switch(link, logit = plogis(lin.pred), probit = pnorm(lin.pred), 
                cloglog = 1 - exp(-exp(lin.pred)))
  if (type == "response") 
    pred <- res
  else pred <- lin.pred
  if (se.fit) {
    cov <- solve(object$hessian)
    var.lin.pred <- diag(newd %*% cov %*% t(newd))
    var.res <- switch(link, logit = exp(2 * lin.pred)/(1 + exp(lin.pred))^4, 
                      probit = dnorm(lin.pred)^2, 
                      cloglog = (exp(-exp(lin.pred))*exp(lin.pred))^2)*var.lin.pred
    if (type == "response") 
      se <- sqrt(var.res)
    else se <- sqrt(var.lin.pred)
    if (!is.null(conf.level)) {
      alpha <- 1 - conf.level
      lower <- lin.pred - qnorm(1 - alpha/2) * sqrt(var.lin.pred)
      upper <- lin.pred + qnorm(1 - alpha/2) * sqrt(var.lin.pred)
      res.lower <- switch(link, logit = plogis(lower), 
                          probit = pnorm(lower), cloglog = 1 - exp(-exp(lower)))
      res.upper <- switch(link, logit = plogis(upper), 
                          probit = pnorm(upper), cloglog = 1 - exp(-exp(upper)))
      if (type == "response") {
        lwr <- res.lower
        upr <- res.upper
      }
      else {
        lwr <- lower
        upr <- upper
      }
    }
  }
  names(pred) <- 1:length(lin.pred)
  if (!is.null(conf.level)) {
    list(fit = pred, se.fit = se, lower = lwr, upper = upr)
  }
  else if (se.fit) 
    list(fit = pred, se.fit = se)
  else pred
}

#