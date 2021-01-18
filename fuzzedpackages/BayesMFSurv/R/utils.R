

# @title betas.slice.sampling2
# @description slice sampling for betas
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param X covariates for betas
# @param betas current value of betas
# @param alpha probability of true censoring
# @param C censoring indicator
# @param lambda current value of lambda
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
# @return One sample update using slice sampling

betas.slice.sampling2 <- function(Sigma.b, Y, Y0, X, betas, alpha, C, lambda, w, m, form){

  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    betas[p] = univ.betas.slice.sampling2(betas[p], p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, w, m, form = form)
  }
  return(betas)

}

# @title univ.betas.slice.sampling2
# @description univariate slice sampling for betas.p
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param X covariates for betas
# @param betas current value of betas
# @param alpha probability of true censoring
# @param C censoring indicator
# @param lambda current value of lambda
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
# @return One sample update using slice sampling

univ.betas.slice.sampling2 <- function(betas.p, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, w, m, lower = -Inf, upper = +Inf, form){

  b0 = betas.p
  b.post0 = betas.post2(b0, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, form)
  if (exp(b.post0) > 0) { b.post0 = log(runif(1, 0, exp(b.post0)))}

  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post2(L, p, Sigma.b, Y, Y0,X, betas, alpha, C, lambda, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post2(R, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (betas.post2(L, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, form) <= b.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (betas.post2(R, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, form) <= b.post0) break
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
    b.post1 = betas.post2(b1, p, Sigma.b, Y,  Y0, X, betas, alpha, C, lambda, form)

    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)

}

# @title gammas.slice.sampling2
# @description slice sampling for gammas
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param lambda current value of lambda
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
# @return One sample update using slice sampling

gammas.slice.sampling2 <- function(Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, w, m, form){

  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling2(gammas[p], p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, w, m, form = form)
  }
  return(gammas)

}

# @title univ.gammas.slice.sampling2
# @description univariate slice sampling for gammas.p
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param lambda current value of lambda
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
# @return One sample update using slice sampling


univ.gammas.slice.sampling2 <- function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, w, m, lower = -Inf, upper = +Inf, form){

  g0 = gammas.p
  g.post0 = gammas.post2(g0, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form)
  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, exp(g.post0)))}

  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post2(L, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post2(R, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post2(L, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post2(R, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form) <= g.post0) break
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
    g.post1 = gammas.post2(g1, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)

}

# @title lambda.slice.sampling2
# @description univariate slice sampling for lambda
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param eXB exponentiated vector of covariates times betas
# @param alpha probability of true censoring
# @param C censoring indicator
# @param lambda current value of lambda
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @return One sample update using slice sampling

lambda.slice.sampling2 <- function(Y, Y0, eXB, alpha, C, lambda, w, m, lower = 0.01, upper = +Inf){

  l0 = lambda
  l.post0 = lambda.post2(Y, Y0, eXB, alpha, C, l0)
  if (exp(l.post0) > 0) { l.post0 = log(runif(1, 0, exp(l.post0)))}

  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (lambda.post2(Y, Y0, eXB, alpha, C, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (lambda.post2(Y, Y0, eXB, alpha, C, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (lambda.post2(Y, Y0, eXB, alpha, C, L) <= l.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (lambda.post2(Y, Y0, eXB, alpha, C, R) <= l.post0) break
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
    l.post1 = lambda.post2(Y, Y0, eXB, alpha, C, l1)
    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)

}


# @title betas.post2
# @description log-posterior distribution of betas with pth element fixed as betas.p
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param X covariates for betas
# @param betas current value of betas
# @param alpha probability of true censoring
# @param C censoring indicator
# @param lambda current value of lambda
# @param form type of parametric model (Exponential or Weibull)
# @return log- posterior density of betas

betas.post2 <- function(betas.p, p, Sigma.b, Y, Y0,  X, betas, alpha, C, lambda, form){

  betas[p] = betas.p
  if (form %in% "Weibull") {
    eXB = exp(X %*% betas)
  } else {
    eXB = exp(X %*% betas)
  }
  lprior = dmvnorm(betas, rep(0, length(betas)), Sigma.b, log = TRUE)
  lpost = llikWeibull2(Y, Y0, eXB, alpha, C, lambda) + lprior
  return(lpost)

}

# @title gammas.post2
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param lambda current value of lambda
# @param form type of parametric model (Exponential or Weibull)
# @return log- posterior density of betas


gammas.post2 <- function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form){

  gammas[p] = gammas.p
  if (form %in% "Weibull") {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  } else {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  }
  lprior = dmvnorm(gammas, rep(0, length(gammas)), Sigma.g, log = TRUE)
  lpost = llikWeibull2(Y, Y0, eXB, alpha, C, lambda) + lprior
  return(lpost)

}

# @title lambda.post2
# @description log-posterior distribution of lambda
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param eXB exponentiated vector of covariates times betas
# @param alpha probability of true censoring
# @param C censoring indicator
# @param lambda current value of lambda
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
# @return log- posterior density of betas


lambda.post2 <- function(Y, Y0, eXB, alpha, C, lambda, a = 1, b = 1){

  lprior = dgamma(lambda, a, b, log = TRUE)
  lpost = llikWeibull2(Y, Y0, eXB, alpha, C, lambda) + lprior
  return(lpost)

}

jointpost = function(Y, Y0, X, Z, betas, Sigma.b, gammas, Sigma.g, alpha, C, lambda, a = 1, b = 1, form){

  if (form %in% "Weibull") {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  } else {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  }
  if (form %in% "Weibull") {
    eXB = exp(X %*% betas)
  } else {
    eXB = exp(X %*% betas)
  }
  lprior0 = dmvnorm(betas, rep(0, length(betas)), Sigma.b, log = TRUE)
  lprior1 = dmvnorm(gammas, rep(0, length(gammas)), Sigma.g, log = TRUE)
  lprior2 = dgamma(lambda, a, b, log = TRUE)
  lpost = llikWeibull2(Y, Y0, eXB, alpha, C, lambda) + lprior0 + lprior1 + lprior2
  return(lpost)

}


# @title bayes.mfsurv.est
# @description Raw form of Markov Chain Monte Carlo (MCMC) to run Bayesian parametric MF model. For user-friendly formula-oriented command, use \code{ \link{mfsurv}}.
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param C the censoring (or failure) dependent variable for the misclassification stage.
# @param X covariates for the survival stage.
# @param Z covariates for the
# @param N number of MCMC iterations.
# @param burn burn-ins to be discarded.
# @param thin thinning to prevent autocorrelation of chain of samples by only taking the n-th values.
# @param w size of the slice in the slice sampling for (betas, gammas, lambda). The default is c(1,1,1). This value may be changed by the user to meet one's needs.
# @param m limit on steps in the slice sampling. The default is 10. This value may be changed by the user to meet one's needs.
# @param form type of parametric model distribution to be used. Options are "Exponential" or "Weibull". The default is "Weibull".
# @param na.action a function indicating what should happen when NAs are included in the data. Options are "na.omit" or "na.fail". The default is "na.omit".

bayes.mfsurv.est <- function(Y, Y0,  C, X, Z, N, burn, thin, w, m, form, na.action){

  # na.action
  na.action <- na.action

  p1 = dim(X)[2]
  p2 = dim(Z)[2]

  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  lambda = 1
  if (form %in% "Weibull") {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  } else {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  }
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  lambda.samp = rep(NA, (N - burn) / thin)
  jointpost.samp = rep(NA, (N - burn) / thin)
  for (iter in 1:N) {
    if (iter %% 1000 == 0) print(iter)
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    betas = betas.slice.sampling2(Sigma.b, Y, Y0, X, betas, alpha, C, lambda, w[1], m, form = form)

    eXB = exp(X %*% betas)

    gammas = gammas.slice.sampling2(Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, w[2], m, form = form)
    if (form %in% "Weibull") {
      alpha = 1 / (1 + exp(-Z %*% gammas))
    } else {
      alpha = 1 / (1 + exp(-Z %*% gammas))
    }
    if (form %in% "Weibull") {
      lambda = lambda.slice.sampling2(Y, Y0, eXB, alpha, C, lambda, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      lambda.samp[(iter - burn) / thin] = lambda
      jointpost.samp[(iter - burn) / thin] = jointpost(Y, Y0, X, Z, betas, Sigma.b, gammas, Sigma.g, alpha, C, lambda, form = form)
    }
  }
  betas = as.data.frame(betas.samp)
  gammas = as.data.frame(gammas.samp)

  lambda = lambda.samp
  post = jointpost.samp
  names(betas) <- colnames(X)
  names(gammas) <- colnames(Z)

  Y <- as.matrix(Y)
  Y0 <- as.matrix(Y0)
  C <- as.matrix(C)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  return(list(Y = Y, Y0 = Y0, C = C, X = X, Z = Z, betas = betas, gammas = gammas, lambda = lambda, post = post, iterations = N, burn_in = burn,
              thinning = thin, betan = nrow(betas), gamman = nrow(gammas), distribution = form))

}


# @title bayes.mfsurv.default
# @description Fit a parametric Bayesian MF model via Markov Chain Monte Carlo (MCMC) to estimate the misclassification in the first stage
#  and the hazard in the second stage.
# @param Y the time (duration) dependent variable for the survival stage (t).
# @param Y0 the elapsed time since inception until the beginning of time period (t-1).
# @param C the censoring (or failure) dependent variable for the misclassification stage (t-1).
# @param X covariates for the survival stage.
# @param Z covariates for the misclassification stage.
# @param N number of MCMC iterations.
# @param burn burn-ins to be discarded.
# @param thin thinning to prevent autocorrelation of chain of samples by only taking the n-th values.
# @param w size of the slice in the slice sampling for (betas, gammas, lambda). The default is c(1,1,1). This value may be changed by the user to meet one's needs.
# @param m limit on steps in the slice sampling. The default is 10. This value may be changed by the user to meet one's needs.
# @param form type of parametric model distribution to be used. Options are "Exponential" or "Weibull". The default is "Weibull".
# @param na.action a function indicating what should happen when NAs are included in the data. Options are "na.omit" or "na.fail". The default is "na.omit".

bayes.mfsurv.default<-function(Y, Y0, C, X, Z, N, burn, thin, w, m, form, na.action){

  Y <- as.numeric(Y)
  Y0 <- as.numeric(Y0)
  C <- as.numeric(C)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  N <- as.numeric(N)
  burn <- as.numeric(burn)
  thin <- as.numeric(thin)
  m <- as.numeric(m)
  w <- as.vector(w)
  form <- as.character(form)
  na.action <- na.action
  est <- bayes.mfsurv.est(Y, Y0, C, X, Z, N, burn, thin, w, m, form, na.action)
  est$call <- match.call()
  class(est) <- "mfsurv"
  est

}



# @title llFun
# @description A function to calculate the log-likelihood of the data.
# @param est a matrix of posterior distribution sample such as posterior mean or the full chain of posteior samples.
# @param Y a matrix the time (duration) dependent variable for the survival stage (t).
# @param Y0 a matrix of the elapsed time since inception until the beginning of time period (t-1).
# @param C a matrix of the censoring (or failure) dependent variable for the misclassification stage.
# @param X a matrix of covariates for the survival stage.
# @param Z a matrix of covariates for the misclassification stage.
# @param data a data frame that contains the Y, Y0, C, X, and Z variables.


llFun <- function(est,Y,Y0,C,X,Z,data){		#Note the extra variable Y0 passed to the time varying MF Weibull

  n <- nrow(data)
  llik <- matrix(0, nrow = n, ncol = 1)
  gamma <- est[1:ncol(Z)]
  beta <- est[(ncol(Z)+1):(length(est)-1)]
  p <- est[length(est)]
  p <- exp(p)
  XB <- X%*%beta
  ZG <- Z%*%gamma
  phi <- 1/(1+exp(-(ZG/p)))
  llik <- C*(log((1-phi)+phi*exp(XB/p)*p*((exp(XB/p)*Y)^(p-1))*exp(-(exp(XB/p)*Y)^p))/exp(-(exp(XB/p)*Y0)^p))+(1-C)*(log(phi)+-(exp(XB/p)*Y)^p--(exp(XB/p)*Y0)^p)
  one <- nrow(llik)
  llik <- subset(llik, is.finite(llik))
  two <- nrow(llik)
  llik <- subset(llik, llik[,1] > -1000)
  three <- nrow(llik)
  llik <- -1*sum(llik)
  list(llik = llik, one = one, two = two, three = three)

}
