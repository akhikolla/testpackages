if(!require(HiddenMarkov)) {
    install.packages("HiddenMarkov", repos = "http://cran.us.r-project.org")
}
if(!require(mclust)) {
    install.packages("mclust", repos = "http://cran.us.r-project.org")
}
if(!require(mvtnorm)) {
    install.packages("mvtnorm", repos = "http://cran.us.r-project.org")
}
if(!require(MCMCpack)) {
  install.packages("MCMCpack", repos = "http://cran.us.r-project.org")
}
if(!require(Rcpp)) {
    install.packages("Rcpp", repos = "http://cran.us.r-project.org")
}

library(HiddenMarkov)
library(mclust)
library(mvtnorm)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")

# identify local region for locally restricted importance sampling

find_index_local_region <- function(samples, EM_result, df_t){

  p <- EM_result$pi
  Mu <- EM_result$Mu
  Sigmas <- EM_result$Sigmas
  K <- length(p)
  d <- length(Mu[1, ])
  inv_Sigmas <- EM_result$inv_Sigmas

  idx_within_region <- array(0, c(dim(samples)[1], K))
  for(k in 1:K){
    idx_within_region[, k] <- apply(samples, 1, function(x) {t(x - Mu[k, ]) %*% inv_Sigmas[[k]] %*% (x - Mu[k, ]) < d})
  }
  idx <- apply(idx_within_region, 1, function(x) {max(x)})
  idx <- which(idx > 0)
  if(df_t == 0){
    log_alpha <- log(pchisq(d, df = d))
  }
  if(df_t > 0){
    log_alpha <- log(pf(d, df1 = 1, df2 = df_t))
  }
  return(list(idx = idx, log_alpha = log_alpha))

}

# multivariate t #

logmvTdensity <- function(x, mu, sqrt_inv_sigma, lgsqrt_det_sigma, df, d, logconstT){

  # logconstT = lgamma((df + d) / 2) - lgamma(df / 2) - log(PI * df) * d / 2

  y <- sqrt_inv_sigma %*% (x - mu)
  ll = - log(sum(y^2) / df + 1) * (df + d) / 2 + logconstT - lgsqrt_det_sigma

  return(ll)

}

# multivariate normal #

logmvNormdensity <- function(x, mu, sqrt_inv_sigma, lgsqrt_det_sigma, d, logconstnormal){

  # logconstnormal = log(sqrt(2 * PI))

  y <- sqrt_inv_sigma %*% (x - mu)
  ll = - sum(y^2) / 2 - d * logconstnormal - lgsqrt_det_sigma

  return(ll)

}

# multivariate density, df_t = 0 corresponds to multivariate Gaussian mixture #

multivariate_mixture_density <- function(x, EM_result, df_t, logconstnormal, logconstT){

  p <- EM_result$pi
  Mu <- EM_result$Mu
  Sigmas <- EM_result$Sigmas
  sqrt_inv_sigma <- EM_result$sqrt_inv_sigma
  lgsqrt_det_sigma <- EM_result$lgsqrt_det_sigma
  K <- length(p)
  d <- dim(Mu)[2]

  if(df_t > 0){
    log_comp <- unlist(sapply(1:K, function(k) {logmvTdensity(x, Mu[k, ], sqrt_inv_sigma[[k]], lgsqrt_det_sigma[k], df_t, d, logconstT)}))
  }
  if(df_t == 0){
    log_comp <- unlist(sapply(1:K, function(k) {logmvNormdensity(x, Mu[k, ], sqrt_inv_sigma[[k]], lgsqrt_det_sigma[k], d, logconstnormal)}))
  }

  log_t <- log(sum(p * exp(log_comp - max(log_comp)))) + max(log_comp)

  return(log_t)
}

# sample from a mixture distribution #

sample_mixture <- function(N, list_paras, df_t){
  p <- list_paras$pi
  Mu <- list_paras$Mu
  Sigmas <- list_paras$Sigmas
  K <- length(p)
  d <- length(Mu[1, ])
  idx <- sample(1:K, N, replace = TRUE, prob = p)
  if(df_t == 0){
    samples <- t(sapply(1:N, function(k) rmvnorm(1, Mu[idx[k], ], Sigmas[[idx[k]]])))
  }
  if(df_t > 0){
    samples <- t(sapply(1:N, function(k) rmvt(1, sigma = Sigmas[[idx[k]]], df = df_t, delta = Mu[idx[k], ])))
  }
  return(samples)
}


# check the validity of parameters #

check_para_validity <- function(parameters_in_matrix_form, bool_hmm){

  # bool_hmm = 1: hidden Markov model
  # bool_hmm = 0: Gaussian mixture model
  # bool_hmm = -1: other models without constraint on parameters

  if(1 - is.matrix(parameters_in_matrix_form)){
    stop("parameters must be in matrix for for validity check\n")
  }
  Npara <- dim(parameters_in_matrix_form)[2]
  if(bool_hmm == 1){
    Nparas <- (1:20)^2 + 1:20
    K <- which(Npara == Nparas)
    bool_index_valid_para <- apply(parameters_in_matrix_form, 1, function(x) return(prod(x[(K + 1):(2 * K)] > 0)))
    bool_index_valid_para <- bool_index_valid_para * apply(parameters_in_matrix_form, 1, function(x) return(prod(x[(2 * K + 1):(K^2 + K)] > 0) * prod(x[(2 * K + 1):(K^2 + K)] < 1)))
    for(j in 1:K){
      bool_index_valid_para <- bool_index_valid_para * apply(parameters_in_matrix_form, 1, function(x) return(sum(x[(K + 2 + j * (K - 1)):(2 * K + j * (K - 1))]) < 1))
    }
  }
  if(bool_hmm == 0){
    Nparas <- (1:20) * 3 - 1
    K <- which(Npara == Nparas)
    bool_index_valid_para <- apply(parameters_in_matrix_form, 1, function(x) return(prod(x[(K + 1):(2 * K)] > 0)))
    bool_index_valid_para <- bool_index_valid_para * apply(parameters_in_matrix_form, 1, function(x) return(prod(x[(2 * K + 1):(3 * K - 1)] > 0) * prod(x[(2 * K + 1):(3 * K - 1)] < 1) * (sum(x[(2 * K + 1):(3 * K - 1)]) < 1)))
  }
  if(bool_hmm == -1){
    bool_index_valid_para <- rep(1, dim(parameters_in_matrix_form)[1])
  }
  return(which(bool_index_valid_para == 1))

}

# unnormalized posterior density of HMM and Gaussian mixture #

ll_un_normalized_hmm_gm <- function(paras, yobs, bool_hmm, priors){

  if(bool_hmm){
    Npara <- (1:50)^2 + (1:50)
  }
  if(1 - bool_hmm){
    Npara <- 3 * (1:50) - 1
  }
  paras <- as.numeric(paras)
  K <- which(Npara == length(paras))
  mu <- as.numeric(paras[1:K])
  sigmas <- as.numeric(paras[(K + 1):(2 * K)])
  if(bool_hmm){
    transmat <- diag(K)
    u <- 1
    for(i in 1:K){
      for(j in 1:K){
        if(j != i){
          transmat[i, j] = paras[2 * K + u]
          u <- u + 1
        }
      }
    }
    for(i in 1:K){
      transmat[i, i] <- 1 - sum(transmat[i, -i])
    }

    obj <- dthmm(yobs, transmat, rep(1/K, K), "norm", list(mean = mu, sd = sqrt(sigmas)))

    ll <- logLik(obj)
    for(j in 1:K){
      ll <- ll + dnorm(mu[j], priors$mu_prior_mean[j], priors$mu_prior_sd, log = TRUE)
      ll <- ll + log(priors$s2[j] * priors$nu[j] / 2.0) * priors$nu[j] / 2.0 - lgamma(priors$nu[j] / 2.0) - priors$s2[j] * priors$nu[j] / (2.0 * sigmas[j]) - log(sigmas[j]) * (priors$nu[j] / 2.0 + 1.0)
      ll <- ll + log(ddirichlet(transmat[j, ], priors$P_prior[j, ]))
    }

  }
  if(1 - bool_hmm){
    pis <- as.numeric(paras[(2 * K + 1):(3 * K - 1)])
    pis <- c(1 - sum(pis), pis)
    ll <- 0
    for(i in 1:length(yobs)){
      ll <- ll + log(sum(sapply(1:K, function(k) pis[k] * dnorm(yobs[i], mu[k], sqrt(sigmas[k])))))
    }
    for(j in 1:K){
      ll <- ll + dnorm(mu[j], priors$mu_prior_mean[j], priors$mu_prior_sd, log = TRUE)
      ll <- ll + log(priors$s2[j] * priors$nu[j] / 2.0) * priors$nu[j] / 2.0 - lgamma(priors$nu[j] / 2.0) - priors$s2[j] * priors$nu[j] / (2.0 * sigmas[j]) - log(sigmas[j]) * (priors$nu[j] / 2.0 + 1.0)
      ll <- ll + (priors$pi_prior[j] - 1.0) * log(pis[j] + 0.000000001) - lgamma(priors$pi_prior[j])
    }
    ll <- ll + lgamma(sum(priors$pi_prior))

  }
  return(ll)

}

find_importance_function <- function(x, boolUseMclust){

  d <- dim(x)[2]

  if (!boolUseMclust){

    Mu <- array(apply(x, 2, mean), c(1, d))
    Sigmas <- list(var(x))
    if (rcond(Sigmas[[1]]) <= 1e-6){
      EM_result = NaN
    }
    if (rcond(Sigmas[[1]]) > 1e-6){
      inv_Sigmas <- list(solve(Sigmas[[1]]))
      sqrt_inv_sigma <- list(chol(inv_Sigmas[[1]]))
      lgsqrt_det_sigma <- log(det(Sigmas[[1]])) / 2

      EM_result = list(pi = 1, Mu = Mu, Sigmas = Sigmas, lgsqrt_det_sigma = lgsqrt_det_sigma, inv_Sigmas = inv_Sigmas, sqrt_inv_sigma = sqrt_inv_sigma)
    }
  }

  if (boolUseMclust){
    mclustresult = Mclust(data = x)
    K = mclustresult$G
    Mu = t(mclustresult$parameters$mean)
    Sigmas <- list()
    inv_Sigmas <- list()
    sqrt_inv_sigma <- list()
    lgsqrt_det_sigma <- rep(0, K)
    for(k in 1:K){
      Sigmas[[k]] <- mclustresult$parameters$variance$sigma[, , k]
      inv_Sigmas[[k]] <- solve(Sigmas[[k]])
      sqrt_inv_sigma[[k]] <- chol(inv_Sigmas[[k]])
      lgsqrt_det_sigma[k] <- log(det(Sigmas[[k]])) / 2
    }
    EM_result = list(pi = mclustresult$parameters$pro, Mu = Mu, Sigmas = Sigmas, lgsqrt_det_sigma = lgsqrt_det_sigma, inv_Sigmas = inv_Sigmas, sqrt_inv_sigma = sqrt_inv_sigma)
  }

  return(EM_result)
}


### estiamte normalizing constant based on posterior samples ###
# samples_and_density: matrix, last column being the log-likelihood
# df_t = 0 is Gaussian tail
# RIS: reciprocal importance sampling
# IS: importance Sampling
# Nsamples_resample: number of re-samples
# ll_un_normalized: log posterior

estimateNormalizingConst <- function(SampDens, boolHMM, dft, RIS, IS, NsmpResmp, llUn, Mclust = TRUE){

  # first fit gaussian mixture

  print("Finding proper importance function.")

  ptm <- proc.time()
  d <- dim(SampDens)[2] - 1
  x <- as.matrix(SampDens[, 1:d])
  ly <- unlist(SampDens[, d + 1])

  EM_result <- find_importance_function (x, Mclust)
  if (is.na(EM_result)[1]){
    return(-Inf)
  }

  print("Finished finding importance function.")

  logconstnormal = log(2 * 3.1415926) / 2
  logconstT = 0
  if (dft > 0){
    logconstT = lgamma((dft + d) / 2) - lgamma(dft / 2) - d / 2 * log(dft * 3.1415826)
  }
  print(proc.time() - ptm)

  # assign the local region

  print("Assigning proper local region")

  ptm <- proc.time()
  temp <- find_index_local_region(x, EM_result, dft)
  idx_local_region_gen_t <- intersect(temp$idx, check_para_validity(x, boolHMM))
  log_alpha_t <- temp$log_alpha
  n_omega_gen_t <- length(idx_local_region_gen_t)
  lg_P_omega_t <- log(n_omega_gen_t) - log(dim(x)[1])
  print(proc.time() - ptm)

  if (RIS){

    # reciprocal importance sampling with t tail

    print("Estimating marginal likelihood with reciprocal importance sampling.")
    ptm <- proc.time()
    lg <- unlist(apply(x[idx_local_region_gen_t, ], 1, function(u) {multivariate_mixture_density(u, EM_result, dft, logconstnormal, logconstT)}))
    temp <- lg - ly[idx_local_region_gen_t]
    C <- log_alpha_t - log(sum(exp(temp - max(temp)))) - max(temp) + log(dim(x)[1])
    print(proc.time() - ptm)

  }

  if (IS){

    # reversely, sampling from g

    print("Estimating marginal likelihood with importance sampling.")
    ptm <- proc.time()
    Kfit <- length(EM_result$pi)
    samples_t <- sample_mixture(NsmpResmp, EM_result, dft)

    # find which samples are in the local region
    print("Locate local region.")
    idx_samp_t <- intersect(find_index_local_region(samples_t, EM_result, dft)$idx, check_para_validity(samples_t, boolHMM))

    # importance sampling with gaussian and t tail
    print("Importance sampling")
    ly_t <- apply(samples_t[idx_samp_t, ], 1, function(x) llUn(x))
    lg_t <- apply(samples_t[idx_samp_t, ], 1, function(x) multivariate_mixture_density(x, EM_result, dft, logconstnormal, logconstT))
    temp_t <- ly_t - lg_t
    temp_t <- temp_t[is.na(temp_t) == 0]
    print("Estimating marginal likelihood.")
    C <- log(sum(exp(temp_t - max(temp_t)))) + max(temp_t) - log(NsmpResmp) - lg_P_omega_t
    print(proc.time() - ptm)
  }

  return(C)

}



# Output: resultHMM (EM result or posterior samples)
# METHOD = 1: EM
# METHOD = 2: HMM mcmc
# METHOD = 3: GM mcmc

HMMfit <- function(y, K, METHOD, optionalfit = list()){

    Ngibbs <- ifelse (exists('Ngibbs', optionalfit), optionalfit$Ngibbs, 5000)
    Burnin <- ifelse (exists('Burnin', optionalfit), optionalfit$Burnin, 5000)
    Thin <- ifelse (exists('Thin', optionalfit), optionalfit$Thin, 10)
    Nstart <- ifelse (exists('Nstart', optionalfit), optionalfit$Nstart, 0)
    verbose <- ifelse (exists('verbose', optionalfit), optionalfit$verbose, 0)
    SigmaInitMethod <- ifelse(exists('SigmaInitMethod', optionalfit), optionalfit$SigmaInitMethod, 1)
    L = length(y)
    MuPriorMean = as.numeric(quantile(y, probs = seq(0.05, 0.95, length.out = K)))
    MuPriorSd = 100.0
    MuPriorVar = rep(MuPriorSd^2, K)
    AlphaPrior = rep(1, K^2)
    alphaPi = rep(1, K)
    AInit = rep(1/K, K^2)
    PiInit = rep(1/K, K)
    MuInit = MuPriorMean
    if (SigmaInitMethod == 1){
        Sigma2Init = rep(sum(sapply(y, function(x) min((x - MuInit)^2))) / L, K)
        #print(Sigma2Init)
    }
    if (SigmaInitMethod == 2){
        Sigma2Init = rep((as.numeric(quantile(y, 0.75) - quantile(y, 0.25))) / (2 * K), K)^2
        #print(Sigma2Init)
    }
    nu = rep(3.0, K)
    s2 = Sigma2Init
    if (exists('priors', optionalfit)){
        priors = optionalfit$priors
    }
    if (!exists('priors', optionalfit)){
        priors <- list(mu_prior_mean = MuPriorMean, mu_prior_sd = MuPriorSd, nu = nu, s2 = s2, P_prior = matrix(AlphaPrior, nrow = K, ncol = K, byrow = TRUE), pi_prior = alphaPi)
    }
    tuningparameters <- list(Method = METHOD, Y = y, Kfit = K, Nstart = Nstart, Ngibbs = Ngibbs, Burnin = Burnin, Thin = Thin, MuInit = MuInit, Sigma2Init = Sigma2Init, PiInit = PiInit, nu = nu, s2 = s2, alphaPi = alphaPi, MuPriorMean = MuPriorMean, MuPriorVar = MuPriorVar, AInit = AInit, AlphaPrior = AlphaPrior, updates2EM = 0, verbose = verbose)

    resultHMM = HMMfitting(tuningparameters)

    if (METHOD == 1){
        resultHMM = list(Pi = resultHMM[1, ], mu = resultHMM[2, ], sigma = sqrt(resultHMM[3, ]), P = resultHMM[4:(3 + K), ], loglike = resultHMM[4 + K, 1])
    }

    return(list(resultHMM = resultHMM, priors = priors, tuningparameters = tuningparameters))
}


# Return list contains: obs, hidden

HMMsim <- function(n, optionalsim = list()){

    Ksim = ifelse (exists('Ksim', optionalsim), optionalsim$Ksim, NA)
    transmat = NA
    mu = NA
    sigma = NA
    pivec = NA

    if (exists('P', optionalsim)){
        transmat = optionalsim$P
        Ksim = dim(transmat)[1]
    }
    if (exists('mu', optionalsim)){
        mu = optionalsim$mu
        Ksim = length(mu)
    }
    if (exists('sigma', optionalsim)){
        sigma = optionalsim$sigma
        Ksim = length(sigma)
    }
    if (exists('pi', optionalsim)){
        pivec = optionalsim$pi
        Ksim = length(pivec)
    }
    if (is.na(Ksim)){
        Ksim = 3
    }
    if (is.na(transmat[1])){
        transmat = matrix(1/Ksim, Ksim, Ksim)
    }
    if (is.na(mu[1])){
        mu = 1:Ksim
    }
    if (is.na(sigma[1])){
        temp = mean(mu[2:length(mu)] - mu[1:(length(mu) - 1)])
        sigma = rep(temp/Ksim, Ksim)
    }
    if (is.na(pivec[1])){
        pivec = rep(1/Ksim, Ksim)
    }

    tuningparameters = list(T = n, Mu = mu, Sigma2 = sigma^2, Pi = pivec, A = as.vector(t(transmat)))

    result = HMMsimulate(tuningparameters)

    if (exists('BoolWritetoFile', optionalsim)){
        if (optionalsim$BoolWritetoFile){
            filename = 'HMMtrace.txt'
            if (exists('Filenameoutput', optionalsim)){
                filename = optionalsim$Filenameoutput
            }
            x = matrix(cbind(result$obs, result$hidden), ncol = 2)
            write.table(x, filename, row.names = FALSE, col.names = c('obs', 'hidden'))
        }
    }

    return(result)

}

# returns estimated number of hidden states, marginal likelihoods, posterior mean estimate of model parameters under the estimated number of states, and corresponding posterior samples

HMMmlselect <- function(y, optionalfit = list()){

    Kfits = 2:6
    if (exists('Kfits', optionalfit)){
        Kfits = optionalfit$Kfits
    }

    METHOD = 2
    boolHMM = 1
    if(exists('boolHMM', optionalfit)){
        boolHMM = optionalfit$boolHMM
        METHOD = ifelse(boolHMM, 2, 3)
    }

    temp = lapply(Kfits, function(K){
        HMMmlestimate(y, K, optionalfit = optionalfit)
    })

    NCs = Kfits
    for (k in 1:length(Kfits)){
        NCs[k] = temp[[k]]$C
    }
    kchoice = which.max(NCs)
    Kest = Kfits[kchoice]
    postsamples = temp[[kchoice]]$posteriorsamplesHMM
    if(boolHMM){
        postmeantrans = apply(postsamples[, (2 * Kest + 1):(Kest^2 + Kest)], 2, mean)
        postmeantransmat = matrix(0, Kest, Kest)
        i = 1
        for (j in 1:Kest){
            for (k in 1:Kest){
                if(k != j){
                    postmeantransmat[j, k] = postmeantrans[i]
                    i = i + 1
                }
            }
        }
        for (k in 1:Kest){
            postmeantransmat[k, k] = 1 - sum(postmeantransmat[k, ])
        }
        postmean = list(mu = apply(postsamples[, 1:Kest], 2, mean), sigma = sqrt(apply(postsamples[, (Kest + 1):(2 * Kest)], 2, mean)), P = postmeantransmat)
    }
    if(!boolHMM){
        if (Kest > 2){
            postmeanpi = apply(postsamples[, (2 * Kest + 1):(3 * Kest - 1)], 2, mean)
            postmeanpi = c(1 - sum(postmeanpi), postmeanpi)
        }
        if (Kest == 2){
            postmeanpi = mean(postsamples[, 5])
            postmeanpi = c(1 - postmeanpi, postmeanpi)
        }
        postmean = list(mu = apply(postsamples[, 1:Kest], 2, mean), sigma = sqrt(apply(postsamples[, (Kest + 1):(2 * Kest)], 2, mean)), Pi = postmeanpi)
    }

    return(list(Kest = Kest, MLs = NCs, Parameters = postmean, postsamples = postsamples))

}



### estimate normalizing constant ###

HMMmlestimate <- function(y, K, optionalfit = list()){

    df_t <- ifelse (exists('df_t', optionalfit), optionalfit$df_t, 2)
    RIS <- ifelse (exists('RIS', optionalfit), optionalfit$RIS, 1)
    IS <- ifelse (exists('IS', optionalfit), optionalfit$IS, 0)
    Nsamples_resample <- ifelse (exists('Nsamples_resample', optionalfit), optionalfit$Nsamples_resample, 5000)
    boolUseMclust <- ifelse(exists('boolUseMclust', optionalfit), optionalfit$boolUseMclust, FALSE)

    METHOD = 2
    boolHMM = 1
    if(exists('boolHMM', optionalfit)){
        boolHMM = optionalfit$boolHMM
        METHOD = ifelse(boolHMM, 2, 3)
    }

    ptm <- proc.time()
    print(paste("begin", K, " state model fitting."))
    temp = HMMfit (y, K, METHOD, optionalfit = optionalfit)
    posteriorsamplesHMM = temp$resultHMM
    priors = temp$priors

    print(paste("finished", K, " state model fitting."))
    timeSampling <- proc.time() - ptm
    print("Time for sampling from posterior distribution:")
    print(timeSampling)

    ll_un_normalized <- function(paras){
        ll_un_normalized_hmm_gm(paras, y, boolHMM, priors)
    }

    ptm <- proc.time()
    print(paste("begin estimating marginal likelihood (K = ", K, ")"))
    C = estimateNormalizingConst (posteriorsamplesHMM, boolHMM, df_t, RIS = RIS, IS = IS, Nsamples_resample, ll_un_normalized, boolUseMclust)
    print(paste("finished estimating marginal likelihood (K = ", K, ")"))
    print("Total time for estimating marginal likelihood:")
    timeConstant <- proc.time() - ptm
    print(timeConstant)

    return(list(posteriorsamplesHMM = posteriorsamplesHMM, C = C,
    time = list(timeSampling = timeSampling, timeConstant = timeConstant)))


}

# Estimate BICs using multiple starting points

RobustBIC <- function(y, optionalbic = list()){


    Kfits = 2:6
    if (exists('Kfits', optionalbic)){
        Kfits = optionalbic$Kfits
    }

    temp = lapply(Kfits, function(K){
        HMMfit (y, K, 1, optionalfit = optionalbic)})

    BICs = Kfits
    for (k in 1:length(Kfits)){
        LL = temp[[k]]$resultHMM$loglike
        K = Kfits[k]
        BICs[k] = - 2 * LL + log(length(y)) * (K^2 + 2 * K - 1)
    }
    kchoice = which.max(-BICs)

    return(list(Kest = Kfits[kchoice], BICs = BICs, Parameters = temp[[kchoice]]$resultHMM))


}



# make figures to visualize results

PlotHMM <- function(y, results){

    mu = results$Parameters$mu
    sigma = results$Parameters$sigma
    K = length(mu)
    Pi = rep(1/K, K)
    P = matrix(1/K, K, K)
    if (exists('Pi', results$Parameters)){
        Pi = results$Parameters$Pi
    }
    if (exists('P', results$Parameters)){
        P = results$Parameters$P
    }

    obj <- dthmm(y, P, Pi, "norm", list(mean = mu, sd = sigma))
    x = Viterbi (obj)

    plot(y, type = 'l', col = 'blue', xlab = '', ylab = '', main = paste('K = ', K, ', Observed (Blue) and Fitted (Black) HMM'))
    lines(mu[x], col = 'black')

}













