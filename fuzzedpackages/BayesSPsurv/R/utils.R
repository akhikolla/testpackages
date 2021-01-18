#' @useDynLib BayesSPsurv
#' @importFrom stats dgamma runif rgamma dnorm model.frame as.formula model.matrix model.response na.omit
#' @importFrom MCMCpack riwish
#' @importFrom coda mcmc
#' @import grDevices
#' @import graphics
#' @importFrom FastGP rcpp_rmvnorm rcpp_log_dmvnorm
#' @importFrom dplyr distinct
#' @importFrom reshape2 acast
#' @export

.onUnload <- function (libpath) {
  library.dynam.unload("BayesSPsurv", libpath)
}


nameob <- c('capdist', 'numa', 'numb', 'kmdist', 'midist')

if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(nameob))
  utils::suppressForeignCheck(c(nameob))
}



formcall <- function(duration,
                     immune,
                     Y0,
                     LY,
                     S = NULL,
                     A = NULL,
                     data,
                     N,
                     burn,
                     thin,
                     w = c(1, 1, 1),
                     m = 10,
                     form,
                     prop.var = NULL,
                     model = character())
{

  formula1 <- as.formula(duration)
  formula2 <- as.formula(immune)
  variable <- unique(c(all.vars(formula1), all.vars(formula2)))
  mf1 <- model.frame(formula = duration, data = data)
  mf2 <- model.frame(formula = immune,   data = data)
  X <- model.matrix(attr(mf1, "terms"), data = mf1)
  Z <- model.matrix(attr(mf2, "terms"), data = mf2)
  Y <- as.matrix(model.response(mf1))
  C <- as.matrix(model.response(mf2))
  Y0 <- data[,Y0]
  LY <- data[,LY]
  burn <-  burn
  if (is.null(w)) w <- c(1,1,1) else w <- w
  if (is.null(m)) m <- 10 else m <- m
  form <-  form
  cnx <- colnames(X)
  cnz <- colnames(Z)

  if(model == "SPsurv"){

    dataset <- data.frame(cbind(Y, Y0, C, LY, X, Z))
    dataset <- na.omit(dataset)
    Y  <- as.matrix(dataset[,1])
    Y0 <- as.matrix(dataset[,2])
    C  <- as.matrix(dataset[,3])
    LY <- as.matrix(dataset[,4])
    X  <- as.matrix(dataset[,5:(4+ncol(X))])
    Z  <- as.matrix(dataset[,(5+ncol(X)):ncol(dataset)])
    colnames(X) <- cnx
    colnames(Z) <- cnz
    fm <- list(Y = Y, Y0 = Y0, C = C, LY = LY, X = X, Z = Z, N = N, burn = burn,
               thin = thin, w = w, m = m, form = form)

  } else {

    prop.var <-  prop.var
    S  <- data[, S]
    dataset <- data.frame(cbind(Y, Y0, C, LY, S, X, Z))
    dataset <- na.omit(dataset)
    Y  <- as.matrix(dataset[,1])
    Y0 <- as.matrix(dataset[,2])
    C  <- as.matrix(dataset[,3])
    LY <- as.matrix(dataset[,4])
    S  <- as.matrix(dataset[,5])
    X  <- as.matrix(dataset[,6:(5+ncol(X))])
    Z  <- as.matrix(dataset[,(6+ncol(X)):ncol(dataset)])
    colnames(X) <- cnx
    colnames(Z) <- cnz
    fm <- list(Y = Y, Y0 = Y0, C = C, LY = LY, X = X, Z = Z, S = S, N = N, burn = burn,
               thin = thin, w = w, m = m, form = form, prop.var = prop.var)

  }

  if(!is.null(A)) fm$A <- as.matrix(A)

  return(fm)

}



# @title betas.slice.sampling
# @description slice sampling for betas
#
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling. A vector of values for beta, gamma, rho.
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

betas.slice.sampling = function(Sigma.b,
                                Y,
                                Y0,
                                X,
                                W,
                                betas,
                                delta,
                                C,
                                LY,
                                rho,
                                w,
                                m,
                                form) {
  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    betas[p] = univ.betas.slice.sampling(betas[p], p, Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, w, m, form = form)
  }
  return(betas)
}

# @title betas.slice.sampling2
# @description slice sampling for betas
#
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
betas.slice.sampling2 = function(Sigma.b, Y, Y0, X, W, betas, delta, C,  LY, rho, w, m, form) {
  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    betas[p] = univ.betas.slice.sampling2(betas[p], p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w, m, form = form)
  }
  return(betas)
}


# @title univ.betas.slice.sampling
# @description univariate slice sampling for betas.p
#
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#

univ.betas.slice.sampling = function(betas.p,
                                     p,
                                     Sigma.b,
                                     Y,
                                     Y0,
                                     X,
                                     W,
                                     betas,
                                     delta,
                                     C,
                                     LY,
                                     rho,
                                     w,
                                     m,
                                     lower = -Inf,
                                     upper = +Inf,
                                     form) {
  b0 = betas.p
  b.post0 = betas.post(b0, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form)

  e.b.post0 = exp(b.post0)
  if (is.infinite(e.b.post0)){e.b.post0 = exp(700)}

  if (exp(b.post0) > 0) { b.post0 = log(runif(1, 0, e.b.post0))}

  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y,Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y,Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    b1 = runif(1, L, R)
    b.post1 = betas.post(b1, p, Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, form)

    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)
}



# @title univ.betas.slice.sampling2
# @description univariate slice sampling for betas.p
#
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#
univ.betas.slice.sampling2 = function(betas.p, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w, m, lower = -Inf, upper = +Inf, form) {
  b0 = betas.p
  b.post0 = betas.post2(b0, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form)

  e.b.post0 = exp(b.post0)
  if (is.infinite(e.b.post0)){e.b.post0 = exp(700)}

  if (exp(b.post0) > 0) { b.post0 = log(runif(1, 0, e.b.post0))}

  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post2(L, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post2(R, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, is.numeric(m)))
    K = (is.numeric(m) - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (betas.post2(L, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (betas.post2(R, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    b1 = runif(1, L, R)
    b.post1 = betas.post2(b1, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form)
    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)
}




# @title gammas.slice.sampling
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

gammas.slice.sampling = function(Sigma.g,
                                 Y,
                                 Y0,
                                 eXB,
                                 Z,
                                 gammas,
                                 C,
                                 LY,
                                 rho,
                                 w,
                                 m,
                                 form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling(gammas[p], p, Sigma.g, Y,Y0, eXB, Z, gammas, C, LY, rho, w, m, form = form)
  }
  return(gammas)
}


# @title gammas.slice.sampling2
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling. A vector of values for beta, gamma, rho.
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

gammas.slice.sampling2 = function(Sigma.g,
                                  Y,
                                  Y0,
                                  eXB,
                                  Z,
                                  V,
                                  gammas,
                                  C,
                                  LY,
                                  rho,
                                  w,
                                  m,
                                  form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling2(gammas[p], p, Sigma.g, Y,Y0, eXB, Z, V, gammas, C, LY, rho, w, m, form = form)
  }
  return(gammas)
}

# @title gammas.slice.sampling3
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#

gammas.slice.sampling3 = function(Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, w, m, form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling3(gammas[p], p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, w, m, form = form)
  }
  return(gammas)
}
# @title gammas.slice.sampling4
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#

gammas.slice.sampling4 = function(Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, w, m, form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling4(gammas[p], p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, w, m, form = form)
  }
  return(gammas)
}



# @title univ.gammas.slice.sampling2
# @description univariate slice sampling for gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling. A vector of values for beta, gamma, rho.
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

univ.gammas.slice.sampling2 = function(gammas.p,
                                       p, Sigma.g,
                                       Y,
                                       Y0,
                                       eXB,
                                       Z,
                                       V,
                                       gammas,
                                       C,
                                       LY,
                                       rho,
                                       w,
                                       m,
                                       lower = -Inf,
                                       upper = +Inf,
                                       form) {
  g0 = gammas.p
  g.post0 = gammas.post2(g0, p, Sigma.g, Y,Y0, eXB, Z, V, gammas, C, LY, rho, form)

  e.g.post0 = exp(g.post0)

  if (is.infinite(e.g.post0)){e.g.post0 = exp(700)}

  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, e.g.post0))}

  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post2(L, p, Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post2(R, p, Sigma.g, Y,Y0, eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post2(L, p, Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post2(R, p, Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post2(g1, p, Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

# @title univ.gammas.slice.sampling
# @description univariate slice sampling for gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

univ.gammas.slice.sampling = function(gammas.p,
                                      p,
                                      Sigma.g,
                                      Y,
                                      Y0,
                                      eXB,
                                      Z,
                                      gammas,
                                      C,
                                      LY,
                                      rho,
                                      w,
                                      m,
                                      lower = -Inf,
                                      upper = +Inf,
                                      form) {
  g0 = gammas.p
  g.post0 = gammas.post(g0, p, Sigma.g, Y, Y0,eXB, Z,gammas, C, LY, rho, form)

  e.g.post0 = exp(g.post0)

  if (is.infinite(e.g.post0)){e.g.post0 = exp(700)}

  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, e.g.post0))}


  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y,Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post(g1, p, Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}


# @title univ.gammas.slice.sampling3
# @description univariate slice sampling for gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
univ.gammas.slice.sampling3 = function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, w, m, lower = -Inf, upper = +Inf, form) {
  g0 = gammas.p
  g.post0 = gammas.post3(g0, p, Sigma.g, Y, Y0, eXB, Z,gammas, C, LY, rho, form)

  e.g.post0 = exp(g.post0)

  if (is.infinite(e.g.post0)){e.g.post0 = exp(700)}

  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, e.g.post0))}


  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post3(L, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post3(R, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post3(L, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post3(R, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post3(g1, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

# @title univ.gammas.slice.sampling4
# @description univariate slice sampling for gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#
univ.gammas.slice.sampling4 = function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, w, m, lower = -Inf, upper = +Inf, form) {
  g0 = gammas.p
  g.post0 = gammas.post4(g0, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form)

  e.g.post0 = exp(g.post0)

  if (is.infinite(e.g.post0)){e.g.post0 = exp(700)}

  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, e.g.post0))}

  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post4(L, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post4(R, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, is.numeric(m)))
    K = (is.numeric(m) - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post4(L, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C,  LY, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post4(R, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post4(g1, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}



# @title rho.slice.sampling
# @description univariate slice sampling for rho
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling. A vector of values for beta, gamma, rho.
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
#
# @return One sample update using slice sampling
#

rho.slice.sampling = function(Y,
                              Y0,
                              eXB,
                              delta,
                              C,
                              LY,
                              rho,
                              w,
                              m,
                              lower = 0.01,
                              upper = +Inf) {
  l0 = rho
  l.post0 = rho.post(Y, Y0,eXB, delta, C, LY, l0)

  e.l.post0 = exp(l.post0)
  if (is.infinite(e.l.post0)){e.l.post0 = exp(700)}

  if (exp(l.post0) > 0) { l.post0 = log(runif(1, 0, e.l.post0))}


  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (rho.post(Y, Y0,eXB, delta, C, LY, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (is.na(rho.post(Y, Y0,eXB, delta, C, LY, R))) #browser()
      if (rho.post(Y, Y0,eXB, delta, C, LY, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (rho.post(Y, Y0,eXB, delta, C, LY, L) <= l.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (is.na(rho.post(Y, Y0,eXB, delta, C, LY, R))) #browser()
      if (rho.post(Y, Y0,eXB, delta, C, LY, R) <= l.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    l1 = runif(1, L, R)
    l.post1 = rho.post(Y, Y0,eXB, delta, C, LY, l1)

    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)
}


# @title rho.slice.sampling2
# @description univariate slice sampling for rho
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
#
# @return One sample update using slice sampling
#
rho.slice.sampling2 = function(Y, Y0, eXB, delta, C, LY, rho, w, m, lower = 0.01, upper = +Inf) {
  l0 = rho
  l.post0 = rho.post2(Y, Y0, eXB, delta, C, LY, l0)

  e.l.post0 = exp(l.post0)
  if (is.infinite(e.l.post0)){e.l.post0 = exp(700)}

  if (exp(l.post0) > 0) { l.post0 = log(runif(1, 0, e.l.post0))}


  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (rho.post2(Y, Y0, eXB, delta, C, LY, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (is.na(rho.post2(Y, Y0, eXB, delta, C, LY, R))) #browser()
        if (rho.post2(Y, Y0, eXB, delta, C, LY, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (rho.post2(Y, Y0, eXB, delta, C, LY, L) <= l.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (is.na(rho.post2(Y, Y0, eXB, delta, C, LY, R))) #browser()
        if (rho.post2(Y, Y0, eXB, delta, C, LY, R) <= l.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    l1 = runif(1, L, R)
    l.post1 = rho.post2(Y, Y0, eXB, delta, C, LY, l1)

    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)
}

# @title betas.post
# @description log-posterior distribution of betas with pth element fixed as betas.p
#
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#

betas.post = function(betas.p,
                      p, Sigma.b,
                      Y,
                      Y0,
                      X,
                      W,
                      betas,
                      delta,
                      C,
                      LY,
                      rho,
                      form) {
  betas[p] = betas.p
  eXB = exp(-(X %*% betas) + W)
  lprior = rcpp_log_dmvnorm(Sigma.b, rep(0, length(betas)), betas, FALSE)
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}
# @title betas.post
# @description log-posterior distribution of betas with pth element fixed as betas.p
#
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Log Lokelihood)
#
# @return log- posterior density of betas
#
betas.post2 = function(betas.p, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) {
  betas[p] = betas.p
  eXB = exp(-(X %*% betas) + W)
  lprior = rcpp_log_dmvnorm(Sigma.b, rep(0, length(betas)), betas, FALSE)
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}
# @title gammas.post
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#

gammas.post = function(gammas.p,
                       p,
                       Sigma.g,
                       Y,
                       Y0,
                       eXB,
                       Z,
                       gammas,
                       C,
                       LY,
                       rho,
                       form) {
  gammas[p] = gammas.p
  num = exp(Z %*% gammas)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikWeibull(Y, Y0,eXB, delta, C, LY, rho) + lprior
  return(lpost)
}


# @title gammas.post2
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#

gammas.post2 = function(gammas.p,
                        p,
                        Sigma.g,
                        Y,
                        Y0,
                        eXB, Z, V, gammas, C, LY, rho, form) {
  gammas[p] = gammas.p
  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikWeibull(Y, Y0,eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title gammas.post3
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Log Likelihood)
#
# @return log- posterior density of betas
#

gammas.post3 = function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) {
  gammas[p] = gammas.p
  num = exp(Z %*% gammas)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title gammas.post4
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#

gammas.post4 = function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form) {
  gammas[p] = gammas.p
  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title rho.post
# @description log-posterior distribution of rho
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#

rho.post = function(Y,
                    Y0,
                    eXB,
                    delta,
                    C,
                    LY,
                    rho,
                    a = 1,
                    b = 1) {
  lprior = dgamma(rho, a, b, log = TRUE)
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title rho.post2
# @description log-posterior distribution of rho
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#


rho.post2 = function(Y, Y0, eXB, delta, C, LY, rho, a = 1, b = 1) {
  lprior = dgamma(rho, a, b, log = TRUE)
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title W.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of W (Weibull or exponential)
#

W.post = function(S,
                  A,
                  lambda,
                  Y,
                  Y0,
                  X,
                  W,
                  betas,
                  delta,
                  C,
                  LY,
                  rho) {
  eXB = exp(-(X %*% betas) + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
  adj = which(A[s,] == 1)
  m_j = length(adj)

  if(m_j == 0){
    m_j <- exp(-700)
  } else {
    m_j <- m_j
  }

  wj = S_uniq[which(S_uniq[,1] %in% adj),2]

  if (length(wj) == 0){
    W_j_bar = 0
  } else {
    W_j_bar = mean(wj)
  }

  lprior = lprior + dnorm(S_uniq[s, 2], W_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title W.post2
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of W (Log likelihood)
#
W.post2 = function(S, A, lambda, Y, Y0, X, W, betas, delta, C, LY, rho) {
  eXB = exp(-(X %*% betas) + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
    adj = which(A[s,] == 1)
    m_j = length(adj)
    if(m_j == 0){
      m_j <- exp(-700)
    } else {
      m_j <- m_j
    }

    wj = S_uniq[which(S_uniq[,1] %in% adj),2]
    if (length(wj) == 0){
      W_j_bar = 0
    } else {
      W_j_bar = mean(wj)
    }

    lprior = lprior + dnorm(S_uniq[s, 2], W_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}



# @title W.MH.sampling
# @description MH Sampling for W
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#

W.MH.sampling = function(S,
                         A,
                         lambda,
                         Y,
                         Y0,
                         X,
                         W,
                         betas,
                         delta,
                         C,
                         LY,
                         rho,
                         prop.var) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.var * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  w1 = W.post(S, A, lambda, Y, Y0,X, W_new[S], betas, delta, C, LY, rho)
  w1[which(w1==-Inf)] <- -740
  w2 = W.post(S, A, lambda, Y,Y0, X, W_old[S], betas, delta, C, LY, rho)
  w2[which(w2==-Inf)] <- -740
  temp = w1- w2
  alpha = min(0, temp)
  if (u <= alpha) {
  	W = W_new[S]
  } else {
  	W = W_old[S]
  }
  return(W)
}


# @title W.MH.sampling2
# @description slice sampling for W
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#
W.MH.sampling2 = function(S, A, lambda, Y, Y0, X, W, betas, delta, C, LY, rho, prop.var) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.var * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  w1 = W.post2(S, A, lambda, Y, Y0, X, W_new[S], betas, delta, C, LY, rho)
  #w1[which(w1==-Inf)] <- -740
  #w1[which(w1==Inf)] <- exp(700)
  w2 = W.post2(S, A, lambda, Y, Y0, X, W_old[S], betas, delta, C, LY, rho)
  #w2[which(w2==Inf)] <- exp(700)
  temp = w1- w2
  alpha = min(0, temp)
  if (u <= alpha) {
    W = W_new[S]
  } else {
    W = W_old[S]
  }
  return(W)
}

# @title W.F.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of W
#

W.F.post = function(Sigma.w,
                    S,
                    Y,
                    Y0,
                    X,
                    W,
                    betas,
                    delta,
                    C,
                    LY,
                    rho) {
  eXB = exp(-(X %*% betas) + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:length(unique(S))){
    lprior = lprior + dnorm(S_uniq[s, 2], 0, Sigma.w[s,s])
  }
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title W.F.post2
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of W
#
W.F.post2 = function(Sigma.w, S, Y, Y0, X, W, betas, delta, C, LY, rho) {
  eXB = exp(-(X %*% betas) + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:length(unique(S))){
    lprior = lprior + dnorm(S_uniq[s, 2], 0, Sigma.w[s,s])

  }
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title W.F.MH.sampling (Cure Model with Frailties)
# @description MH sampling for W
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#

W.F.MH.sampling = function(Sigma.w,
                           S,
                           Y,
                           Y0,
                           X,
                           W,
                           betas,
                           delta,
                           C,
                           LY,
                           rho,
                           prop.var) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.var * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  w1 = W.F.post(Sigma.w, S, Y, Y0, X, W_new[S], betas, delta, C, LY, rho)
  w2 = W.F.post(Sigma.w, S, Y, Y0, X, W_old[S], betas, delta, C, LY, rho)
  temp = w1- w2
  alpha = min(0, temp)
  if (u <= alpha) {
    W = W_new[S]
  } else {
    W = W_old[S]
  }
  return(W)
}

# @title W.F.MH.sampling2 (Cure Model with Frailties)
# @description MH sampling for W
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#
W.F.MH.sampling2 = function(Sigma.w, S, Y, Y0, X, W, betas, delta, C, LY, rho, prop.var) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.var * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  w1 = W.F.post2(Sigma.w, S, Y, Y0, X, W_new[S], betas, delta, C, LY, rho)
  w1[which(w1==-Inf)] <- -740
  w2 = W.F.post2(Sigma.w, S,  Y, Y0, X, W_old[S], betas, delta, C, LY, rho)
  w2[which(w2==-Inf)] <- -740
  temp = w1- w2
  alpha = min(0, temp)
  if (u <= alpha) {
    W = W_new[S]
  } else {
    W = W_old[S]
  }
  return(W)
}



# @title V.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas
#

V.post = function(S,
                  A,
                  lambda,
                  Y,
                  Y0,
                  eXB,
                  Z,
                  V,
                  gammas,
                  C,
                  LY,
                  rho) {

  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
  adj = which(A[s,] == 1)

  m_j = length(adj)
  if(m_j == 0){
    m_j <- exp(-700)
  } else {
    m_j <- m_j
  }

  vj = S_uniq[which(S_uniq[,1] %in% adj),2]
  if (length(vj) == 0){
    V_j_bar = 0
  } else {
    V_j_bar = mean(vj)
  }

  lprior = lprior + dnorm(S_uniq[s, 2], V_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title V.post2
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas (log likelihood)
#
V.post2 = function(S, A, lambda, Y, Y0, eXB, Z, V, gammas, C, LY, rho) {

  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  #delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
    adj = which(A[s,] == 1)

    m_j = length(adj)
    if(m_j == 0){
      m_j <- exp(-700)
    } else {
      m_j <- m_j
    }

    vj = S_uniq[which(S_uniq[,1] %in% adj),2]
    if (length(vj) == 0){
      V_j_bar = 0
    } else {
      V_j_bar = mean(vj)
    }

    lprior = lprior + dnorm(S_uniq[s, 2], V_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}





# @title V.MH.sampling
# @description MH sampling for rcpp_log_dmvnorm
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#

V.MH.sampling = function(S,
                         A,
                         lambda,
                         Y,
                         Y0,
                         eXB,
                         Z,
                         V,
                         gammas,
                         C,
                         LY,
                         rho,
                         prop.var) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.var * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  v1 = V.post(S, A, lambda, Y,Y0, eXB, Z, V_new[S], gammas, C, LY, rho)
  v1[which(v1==-Inf)] <- -740
  v2 = V.post(S, A, lambda, Y, Y0,eXB, Z, V_old[S], gammas, C, LY, rho)
  v2[which(v2==-Inf)] <- -740
  temp = v1- v2
  alpha = min(0, temp)
  if (u <= alpha) {
  	V = V_new[S]
  } else {
  	V = V_old[S]
  }
  return(V)
}

# @title V.MH.sampling2
# @description slice sampling for rcpp_log_dmvnorm
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#
V.MH.sampling2 = function(S, A, lambda, Y, Y0, eXB, Z, V, gammas, C, LY, rho, prop.var) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.var * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  v1 = V.post2(S, A, lambda, Y, Y0, eXB, Z, V_new[S], gammas, C, LY, rho)
  #v1[which(v1==-Inf)] <- -740
  #V1[which(V1==Inf)] <- exp(700)
  v2 = V.post2(S, A, lambda, Y, Y0, eXB, Z, V_old[S], gammas, C, LY, rho)
  #v2[which(v2==-Inf)] <- -740
  #v2[which(v2==Inf)] <- exp(700)
  temp = v1- v2
  alpha = min(0, temp)
  if (u <= alpha) {
    V = V_new[S]
  } else {
    V = V_old[S]
  }
  return(V)
}



# @title V.F.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas (Weibull and Exponential)
#

V.F.post = function(Sigma.v,
                    S,
                    Y,
                    Y0,
                    eXB,
                    Z,
                    V,
                    gammas,
                    C,
                    LY,
                    rho) {

  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:length(unique(S))){
    lprior = lprior + dnorm(S_uniq[s, 2], 0, Sigma.v[s,s])
  }
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title V.F.post2
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas (log likelihood)
#

V.F.post2 = function(Sigma.v, S, Y, Y0, eXB, Z, V, gammas, C, LY, rho) {

  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  #delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:length(unique(S))){
    lprior = lprior + dnorm(S_uniq[s, 2], 0, Sigma.v[s,s])
  }
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}


# @title V.F.MH.sampling (Cure Model with non-spatial Frailties)
# @description MH sampling for rcpp_log_dmvnorm
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#

V.F.MH.sampling = function(Sigma.v,
                           S,
                           Y,
                           Y0,
                           eXB,
                           Z,
                           V,
                           gammas,
                           C,
                           LY,
                           rho,
                           prop.var) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.var * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  v1 = V.F.post(Sigma.v,S, Y,Y0, eXB, Z, V_new[S], gammas, C, LY, rho)
  v1[which(v1==-Inf)] <- -740
  v2 = V.F.post(Sigma.v,S, Y, Y0,eXB, Z, V_old[S], gammas, C, LY, rho)
  v2[which(v2==-Inf)] <- -740
  temp = v1- v2
  alpha = min(0, temp)
  if (u <= alpha) {
    V = V_new[S]
  } else {
    V = V_old[S]
  }
  return(V)
}

# @title V.F.MH.sampling2 (Cure Model with non-spatial Frailties)
# @description MH sampling for rcpp_log_dmvnorm
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling (log likelihood)
#
V.F.MH.sampling2 = function(Sigma.v, S,  Y, Y0, eXB, Z, V, gammas, C, LY, rho, prop.var) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.var * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  v1 = V.F.post2(Sigma.v, S, Y, Y0, eXB, Z, V_new[S], gammas, C, LY, rho)
  #v1[which(v1==-Inf)] <- -740
  v2 = V.F.post2(Sigma.v, S, Y, Y0, eXB, Z, V_old[S], gammas, C, LY, rho)
  #v2[which(v2==-Inf)] <- -740
  temp = v1- v2
  alpha = min(0, temp)
  if (u <= alpha) {
    V = V_new[S]
  } else {
    V = V_old[S]
  }
  return(V)
}


# @title lambda.gibbs.sampling2
# @description log-posterior distribution of rho
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param W spatial random effects
# @param V spatial random effects
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#

lambda.gibbs.sampling2 <- function(S,
                                   A,
                                   W,
                                   V,
                                   a = 1,    # prior
                                   b = 1) {  # prior
  S_uniq = unique(cbind(S, W, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  J = nrow(S_uniq)
  sums = 0
  for (j in 1:J) {
    adj = which(A[j,]==1)
    m_j = length(adj)

    if(m_j == 0){
      m_j <- exp(-700)
    } else {
      m_j <- m_j
    }

    wj = S_uniq[which(S_uniq[,1] %in% adj),2]
    if (length(wj) == 0){
      W_j_bar = 0
    } else {
      W_j_bar = mean(wj)
    }

    vj = S_uniq[which(S_uniq[,1] %in% adj),3]
    if (length(vj) == 0){
      V_j_bar = 0
    } else {
      V_j_bar = mean(vj)
    }

    W_j = S_uniq[j, 2]
    V_j = S_uniq[j, 3]

    sums = sums + m_j/2 * ((W_j-W_j_bar)^2 + (V_j-V_j_bar)^2)

  }
  lambda = rgamma(1, J + a, sums + b)
  return(lambda)
}



# @title mcmcSP
# @description Markov Chain Monte Carlo (MCMC) to run Bayesian split population survival model with no frailties
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling. A vector of values for beta, gamma, rho.
# @param form type of parametric model (Exponential or Weibull)
#
# @return chain of the variables of interest
#

mcmcSP <- function(Y,
                   Y0,
                   C,
                   LY,
                   X,
                   Z,
                   N,
                   burn,
                   thin,
                   w = c(1, 1, 1),
                   m = 10,
                   form,
                   propvar) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  rho = 1
  lambda = 1
  W = rep(0, length(Y))
  V = rep(0, length(Y))
  delta = exp(Z %*% gammas)/ (1 + exp(Z %*% gammas))
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  for (iter in 1:N) {
    # if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    betas = betas.slice.sampling(Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp(-(X %*% betas))
    gammas = gammas.slice.sampling(Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, Y0,eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      delta.samp[(iter - burn) / thin] = mean(delta)

    }
  }
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas = betas.samp, gammas = gammas.samp, rho = rho.samp, delta = delta.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, form = form)))
}

# @title mcmcSPlog
# @description Markov Chain Monte Carlo (MCMC) to run Bayesian cure model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Log Likelihood)
#
# @return chain of the variables of interest
#
mcmcSPlog <- function(Y, C, Y0, X, LY, Z, N, burn, thin, w = c(1, 1, 1), m, form) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  rho = 1
  lambda = 1
  W = rep(0, length(Y))
  V = rep(0, length(Y))
  delta = exp(Z %*% gammas)/ (1 + exp(Z %*% gammas))
  Sigma.b = 5 * p1 * diag(p1)
  Sigma.g = 5 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  for (iter in 1:N) {
    #if (iter %% 1000 == 0) print(iter) #F #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    betas = betas.slice.sampling2(Sigma.b, Y, Y0, X, W, betas, delta, C,  LY, rho, w[1], m, form = form)
    eXB = exp((X %*% betas))
    gammas = gammas.slice.sampling3(Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    if (form %in% "loglog") {
      rho = rho.slice.sampling2(Y, Y0, eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
    }
  }
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas = betas.samp, gammas = gammas.samp, rho = rho.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, form = form)))
}



# @title mcmcspatialSP
# @description Markov Chain Monte Carlo (MCMC) routine for Bayesian spatial split population survival model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param LY last observation year
# @param X covariates for betas
# @param Z covariates for gammas
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling. A vector of values for beta, gamma, rho.
# @param form type of parametric model (Exponential or Weibull)
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return chain of the variables of interest
#

mcmcspatialSP <- function(Y,
                          Y0,
                          C,
                          LY,
                          X,
                          Z,
                          S,
                          A,
                          N,
                          burn,
                          thin, w = c(1, 1, 1),
                          m = 10,
                          form,
                          prop.var,
                          id_WV) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  rho = 1
  lambda = 1
  W = rep(0, length(Y))
  V = rep(0, length(Y))
  delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))
  V.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))
  for (iter in 1:N) {
    #if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    #CAR model
    lambda = lambda.gibbs.sampling2(S, A, W, V)
    W = W.MH.sampling(S, A, lambda, Y, Y0,X, W, betas, delta, C, LY, rho, prop.var)
    betas = betas.slice.sampling(Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp(-(X %*% betas) + W)
    V = V.MH.sampling(S, A, lambda, Y,Y0, eXB, Z, V, gammas, C, LY, rho, prop.var)
 	  gammas = gammas.slice.sampling2(Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, w[2], m, form = form)
 	  num = exp(Z %*% gammas + V)
 	  num[which(is.infinite(num))] <- exp(700)
 	  denom = (1 + exp(Z %*% gammas + V))
 	  num[which(is.infinite(denom))] <- exp(700)
 	  delta = num/denom


    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, Y0,eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      delta.samp[(iter - burn) / thin] = mean(delta)
      S_uniq = unique(cbind(S, W, V))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
      V.samp[(iter - burn) / thin, ] = S_uniq[,3]
    }
  }

  colnames(V.samp) <- id_WV #colnames(A)           #ADC
  colnames(W.samp) <- id_WV #colnames(A)           #ADC
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas = betas.samp, gammas = gammas.samp, rho = rho.samp, lambda = lambda.samp,
              delta = delta.samp, W = W.samp, V = V.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, S = S, form = form)))
}


# @title mcmcSpatialLog
# @description Markov Chain Monte Carlo (MCMC) to run Bayesian spatial cure model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return chain of the variables of interest
#

mcmcSpatialLog <- function(Y,
                           Y0,
                           C,
                           LY,
                           X,
                           Z,
                           S,
                           A,
                           N,
                           burn,
                           thin,
                           w = c(1, 1, 1),
                           m,
                           form,
                           prop.var = .00001,
                           id_WV) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  rho = 1
  lambda = 1
  W = rep(0, length(Y))
  V = rep(0, length(Y))
  delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  #delta = 1/(1 + exp(-Z %*% gammas + V))
  Sigma.b = 5 * p1 * diag(p1)
  Sigma.g = 5 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))
  V.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))
  for (iter in 1:N) {
    #if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    #CAR model
    lambda = lambda.gibbs.sampling2(S, A, W, V)
    W = W.MH.sampling2(S, A, lambda, Y, Y0, X, W, betas, delta, C, LY, rho, prop.var)
    betas = betas.slice.sampling2(Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp((X %*% betas) + W)
    V = V.MH.sampling2(S, A, lambda, Y, Y0, eXB, Z, V, gammas, C,  LY, rho, prop.var)
    gammas = gammas.slice.sampling4(Sigma.g, Y, Y0, eXB, Z, V, gammas, C,  LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas + V)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas + V))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    #delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))

    if (form == "loglog") {
      rho = rho.slice.sampling2(Y, Y0, eXB, delta, C, LY, rho, w[3], m)
    }
    # print(delta) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      delta.samp[(iter - burn) / thin] = mean(delta)
      S_uniq = unique(cbind(S, W, V))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
      V.samp[(iter - burn) / thin, ] = S_uniq[,3]
      # print(100) #### ***** #### ***** #### ***** #### ***** ####
    }
  }
  colnames(V.samp) <- id_WV #colnames(A)           #ADC
  colnames(W.samp) <- id_WV #colnames(A)           #ADC
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas = betas.samp, gammas = gammas.samp, rho = rho.samp, delta= delta.samp, lambda = lambda.samp, W = W.samp, V = V.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, S = S, form = form)))
}



# @title mcmcfrailtySP
# @description Markov Chain Monte Carlo (MCMC) routine to run Bayesian non-spatial frailties split population survival model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param LY last observation year
# @param X covariates for betas
# @param Z covariates for gammas
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling. A vector of values for beta, gamma, rho.
# @param form type of parametric model (Exponential or Weibull)
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return chain of the variables of interest
#

mcmcfrailtySP <- function(Y,
                          Y0,
                          C,
                          LY,
                          X,
                          Z,
                          S,
                          N,
                          burn,
                          thin, w = c(1, 1, 1),
                          m,
                          form,
                          prop.var,
                          id_WV) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  p3 = length(unique(S))
  p4 = length(unique(S))
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  rho = 1
  lambda = 1
  W = rep(0, length(Y))
  V = rep(0, length(Y))
  WS = rep(0, p3)
  VS = rep(0, p4)
  delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  Sigma.w = 10 * p3 * diag(p3)
  Sigma.v = 10 * p4 * diag(p4)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = length(unique(S)))
  V.samp = matrix(NA, nrow = (N - burn) / thin, ncol = length(unique(S)))
  for (iter in 1:N) {
    #if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
      Sigma.w = riwish(1 + p3, WS %*% t(WS) + p3 * diag(p3))
      Sigma.v = riwish(1 + p4, VS %*% t(VS) + p4 * diag(p4))
    }
    #non-spatial Frailty model
    W = W.F.MH.sampling(Sigma.w, S,Y, Y0,X, W, betas, delta, C, LY, rho, prop.var)
    betas = betas.slice.sampling(Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp(-(X %*% betas) + W)
    V = V.F.MH.sampling(Sigma.v, S, Y,Y0, eXB, Z, V, gammas, C, LY, rho, prop.var)
    gammas = gammas.slice.sampling2(Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas + V)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas + V))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, Y0,eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      delta.samp[(iter - burn) / thin] = mean(delta)
      S_uniq = unique(cbind(S, W, V))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
      V.samp[(iter - burn) / thin, ] = S_uniq[,3]
    }
  }


  colnames(V.samp) <- id_WV #colnames(A)           #ADC
  colnames(W.samp) <- id_WV #colnames(A)           #ADC
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas = betas.samp, gammas = gammas.samp, rho = rho.samp, lambda = lambda.samp, delta = delta.samp, W = W.samp, V = V.samp,
         spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, S = S, form = form)))
}

# @title mcmc Cure with Non-spatial frailties
# @description Markov Chain Monte Carlo (MCMC) to run Bayesian spatial cure model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param C censoring indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param S spatial information (e.g. district)
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return chain of the variables of interest
#

mcmcfrailtySPlog <- function(Y,
                             Y0,
                             C,
                             LY,
                             X,
                             Z,
                             S,
                             N,
                             burn,
                             thin,
                             w = c(1, 1, 1),
                             m = 10,
                             form,
                             prop.var= .00001,
                             id_WV) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  p3 = length(unique(S))
  p4 = length(unique(S))
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  rho = 1
  lambda = 1
  W = rep(0, length(Y))
  V = rep(0, length(Y))
  WS = rep(0, p3)
  VS = rep(0, p4)
  delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  prop.var = prop.var
  #delta = 1/(1 + exp(-Z %*% gammas + V))
  Sigma.b = 5 * p1 * diag(p1)
  Sigma.g = 5 * p2 * diag(p2)
  Sigma.w = 5 * p3 * diag(p3)
  Sigma.v = 5 * p4 * diag(p4)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = length(unique(S)))
  V.samp = matrix(NA, nrow = (N - burn) / thin, ncol = length(unique(S)))
  for (iter in 1:N) {
    #if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
      Sigma.w = riwish(1 + p3, WS %*% t(WS) + p3 * diag(p3))
      Sigma.v = riwish(1 + p4, VS %*% t(VS) + p4 * diag(p4))
    }
    #non-spatial Frailty model.
    W = W.F.MH.sampling2(Sigma.w, S,Y, Y0, X, W, betas, delta, C, LY, rho, prop.var)
    betas = betas.slice.sampling2(Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp(-(X %*% betas) + W)
    V = V.F.MH.sampling2(Sigma.v, S, Y, Y0, eXB, Z, V, gammas, C, LY, rho, prop.var)
    gammas = gammas.slice.sampling4(Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas + V)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas + V))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    if (form == "loglog") {
      rho = rho.slice.sampling2(Y, Y0, eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      delta.samp[(iter - burn) / thin] = mean(delta)
      S_uniq = unique(cbind(S, W, V))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
      V.samp[(iter - burn) / thin, ] = S_uniq[,3]
    }
  }


  colnames(V.samp) <- id_WV #colnames(A)           #ADC
  colnames(W.samp) <- id_WV #colnames(A)           #ADC
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas = betas.samp, gammas = gammas.samp, rho = rho.samp, lambda = lambda.samp, delta = delta.samp, W = W.samp, V = V.samp,
         spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, S = S, form = form)))
}




# @title llFun
# @description Log-likelihood function for exchnageable and spatial frailty split population models
#
# @param est parameter estimates from the model
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C immunity indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param W survival stage (spatial) frailty parameter estimates
# @param V split stage (spatial) frailty parameter estimates
# @param data data that contains the variables of interest
#
# @return Log-likelihood of the exchangeable or the spatial split population survial model
#

llFun <- function(est,
                  Y,
                  Y0,
                  C,
                  X,
                  Z,
                  W,
                  V,
                  data,
                  form){		#Note the extra variable Y0 passed to the time varying MF Weibull

  n     <- nrow(data)
  llik  <- matrix(0, nrow = n, ncol = 1)
  gamma <- est[1:ncol(Z)]
  beta  <- est[(ncol(Z)+1):(ncol(Z)+ncol(X))]
  p     <- est[length(est)]
  XB    <- X %*% beta
  ZG    <- Z %*% gamma
  phi   <- exp(-ZG + V)/(1+exp(-ZG + V))
  eXB   <- exp(-XB + W)
  if(form == "loglog"){
    llik <- C*(log((1-phi)*eXB*p*((eXB*Y)^(p-1))*(1+((eXB*Y0)^(p)))/(1+exp(-(eXB*Y)^(p)))^2)+(1-C)*(log(phi+(1-phi)*(1+((eXB*Y0)^p)))^2/(1+((eXB*Y)^p)))^2)
  } else {
    llik <- C*(log((1-phi)*eXB*p*((eXB*Y)^(p-1))*exp(-(eXB*Y))^p/exp(-(eXB*Y0))^p))+(1-C)*(log(phi+(1-phi)*((exp(-eXB*Y))^p)/((exp(-eXB*Y0))^p)))
  }
  one   <- nrow(llik)
  llik  <- subset(llik, is.finite(llik))
  two   <- nrow(llik)
  llik  <- subset(llik, llik[,1] > -1000)
  three <- nrow(llik)
  llik  <- -1*sum(llik)
  list(llik = llik, one = one, two = two, three = three)

}

# @title rllFun
# @description Log-likelihood function for pooled split population models
#
# @param est parameter estimates from the model
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C immunity indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param data data that contains the variables of interest
#
# @return Log-likelihood of the poold slit population survival model
#

rllFun <- function(est,
                   Y,
                   Y0,
                   C,
                   X,
                   Z,
                   data,
                   form){		#Note the extra variable Y0 passed to the time varying MF Weibull

  n     <- nrow(data)
  llik  <- matrix(0, nrow = n, ncol = 1)
  gamma <- est[1:ncol(Z)]
  beta  <- est[(ncol(Z)+1):(ncol(Z)+ncol(X))]
  p     <- est[length(est)]
  XB    <- X %*% beta
  ZG    <- Z %*% gamma
  phi   <- exp(-ZG )/(1+exp(-ZG ))
  eXB   <- exp(-XB )
  if(form == "loglog"){
    llik <- C*(log((1-phi)*eXB*p*((eXB*Y)^(p-1))*(1+((eXB*Y0)^(p)))/(1+exp(-(eXB*Y)^(p)))^2)+(1-C)*(log(phi+(1-phi)*(1+((eXB*Y0)^p)))^2/(1+((eXB*Y)^p)))^2)
  } else {
    llik <- C*(log((1-phi)*eXB*p*((eXB*Y)^(p-1))*exp(-(eXB*Y))^p/exp(-(eXB*Y0))^p))+(1-C)*(log(phi+(1-phi)*((exp(-eXB*Y))^p)/((exp(-eXB*Y0))^p)))
  }
  one   <- nrow(llik)
  llik  <- subset(llik, is.finite(llik))
  two   <- nrow(llik)
  llik  <- subset(llik, llik[,1] > -1000)
  three <- nrow(llik)
  llik  <- -1*sum(llik)
  list(llik = llik, one = one, two = two, three = three)

}


JointCount<- function(X, W) {
  PB<- mean (X)
  PW<- 1- mean(X)
  PBB<- PB*PB
  PWW<- PW*PW
  PBW<- 2*PB*(1-PB)
  n<-length(X)
  nB<-sum(X)
  EBB<- 0.5*sum(W)*PBB
  EWW<- 0.5*sum(W)*PWW
  EBW<- 0.5*sum(W)*PBW
  s.2j<-rep(NA,ncol(W))
  s.2i<-rep(NA,nrow(W))
  s.3j<-rep(NA,ncol(W))
  s.3i<-rep(NA,nrow(W))
  s_2<-matrix(NA,nrow(W),ncol(W))
  bb<-matrix(NA,nrow(W),ncol(W))
  ww<-matrix(NA,nrow(W),ncol(W))
  bw<-matrix(NA,nrow(W),ncol(W))
  for (i in 1:nrow(W)){
    for (j in 1:ncol(W)){
      s_2[i,j]<-(W[i,j]+W[j,i])^2
      bb[i,j]<-W[i,j]*X[i]*X[j]
      ww[i,j]<-W[i,j]*(1-X[i])*(1-X[j])
      bw[i,j]<-W[i,j]*((X[i]-X[j])^2)
    }
    s.3j[i]<-(sum(W[i,])+ sum(W[,i]))^2
  }
  BB<-0.5*sum(bb)
  WW<-0.5*sum(ww)
  BW<-0.5*sum(bw)
  s2<-sum(s_2)*0.5
  s3<-sum(s.3j)
  s1<-sum(W)
  varBW<-0.25*(((2*s2*nB*(n-nB))/(n*(n-1)))+(((s3-s1)*nB*(n-nB))/(n*(n-1)))+((4*(s1^2+s2-s3)*nB*(nB-1)*(n-nB)*(n-nB-1))/(n*(n-1)*(n-2)*(n-3)) ))-(EBW^2)
  ZBW<-(BW-EBW)/sqrt(varBW)
  out<-data.frame(BW,EBW,sqrt(varBW),ZBW)
  colnames(out)<-c("Obs","Exp","sd","Z")
  return(out)
}

data_plots <- function(data, var_id = character(), var_time = character(), n = 0){

  dat  <- lapply(split(data, data[,var_time]), function(x){x[!duplicated(x[, var_id]), ]})
  dat  <- dat[sapply(dat, nrow) >= n]
  mats <- lapply(dat, function(x){BayesSPsurv::spatial_SA(x, var_ccode = var_id)$matrixA})
  w2 <- list()
  for(i in 1:length(dat)){w2[[i]] <- dat[[i]][order(dat[[i]][, var_id]),]}
  list(mats = mats, w2 = w2)

}




