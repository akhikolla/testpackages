.groveEB.all.random <- function(W, 
                                formula, 
                                X, 
                                alpha, 
                                nu, 
                                eta.rho,
                                eta.kappa,
                                gamma.rho,
                                gamma.kappa,
                                n.samples = 500, 
                                verbose = FALSE,
                                transition.mode = "Markov",
                                method = "Nelder-Mead") {
  
  input <- .ReshapeInput(W,
                         formula,
                         transition.mode,
                         X)
  
  kNuInf <- .GetNuInf()
  beta <- .InitBetaPar()
  sigma <- .InitSigmaPar(W)
  tau <- .InitTauPar(W, sigma, input$p_len)
  
  fr <- function(data) { 
    tau.internal <- exp(data[1 : input$p_len])
    eta.rho.internal <- exp(data[input$p_len + 1])
    eta.kap.internal <- rep(exp(data[input$p_len + 2]), input$p_len - 1)
    gamma.rho.internal <- 0.5 * exp(data[input$p_len + 3]) / 
      (1 + exp(data[input$p_len + 3]))      
    gamma.kap.internal <- rep(0.5 * exp(data[input$p_len + 4]) / 
                                (1 + exp(data[input$p_len + 4])), 
                              input$p_len - 1)      
    sigma.internal <- exp(data[input$p_len + 5])
    alpha.internal <- data[input$p_len + 6]
    # If the nu_par is set to infinity, fit the model with stationary noise;
    # otherwise take in the initial nu value for EB.
    nu.internal <- ifelse(nu != Inf, exp(data[input$p_len + 7]), kNuInf)
    
    output <- fitGroveML(W$D, 
                         input$XX, 
                         p = input$p, 
                         tau_par = tau.internal, 
                         eta_par = c(eta.rho.internal, eta.kap.internal),
                         gamma_par = c(gamma.rho.internal, gamma.kap.internal),
                         init_state = input$init_state,
                         nu = nu.internal, 
                         sigma0 = sigma.internal, 
                         alpha = alpha.internal, 
                         beta = beta, 
                         transition_mode = input$transition.mode)
    if (verbose == TRUE) {
      print(-output$marginal_likelihood)
    }    
    return(-output$marginal_likelihood)
  } 
  # Organize starting values in a vector.
  eb.par <- c(log(tau),
              log(eta.rho), 
              log(eta.kappa),
              log(gamma.rho / (0.5 - gamma.rho)), 
              log(gamma.kappa / (0.5 - gamma.kappa)),
              log(sigma),
              alpha)
  # Run EB.
  empirical_bayes <- ifelse(nu != Inf,
                            optim(par = c(eb.par, log(nu)), 
                                  fn = fr,  
                                  method = method),
                            optim(par = eb.par, 
                                  fn = fr,  
                                  method = method))
  # Define final values for EB.
  tau_par_final <- exp(empirical_bayes[[1]][1 : input$p_len])
  eta_par_final <- c(exp(empirical_bayes[[1]][input$p_len + 1]), 
                     rep(exp(empirical_bayes[[1]][input$p_len + 2]), 
                         input$p_len - 1))
  gamma_par_final <- c(0.5 * exp(empirical_bayes[[1]][input$p_len + 3]) / 
                         (1 + exp(empirical_bayes[[1]][input$p_len + 3])),
                       rep(0.5 * exp(empirical_bayes[[1]][input$p_len + 4]) /
                             (1 + exp(empirical_bayes[[1]][input$p_len + 4])), 
                           input$p_len - 1))
  sigma_par_final <- exp(empirical_bayes[[1]][input$p_len + 5])
  alpha_par_final <- empirical_bayes[[1]][input$p_len + 6]
  beta_par_final <- beta
  nu_par_final  <- ifelse(nu != Inf, 
                          exp(empirical_bayes[[1]][input$p_len + 7]),
                          kNuInf)
  # Run model based on optimal hyperparameters.
  ans <- .DenoiseFinalRun(W, 
                          input$XX, 
                          X,
                          formula,
                          input$C_hat,
                          input$p,
                          tau_par_final,
                          eta_par_final,
                          gamma_par_final,
                          input$init_state,
                          nu_par_final,
                          sigma_par_final,
                          alpha_par_final,
                          beta_par_final,
                          n.samples,
                          input$transition.mode)
  return(ans)  
}  
  
.groveEB.fixed.kappa <- function(W, 
                                 formula, 
                                 X, 
                                 alpha, 
                                 nu, 
                                 eta.rho,
                                 eta.kappa,
                                 gamma.rho,
                                 gamma.kappa,
                                 n.samples = 500, 
                                 verbose = FALSE,
                                 transition.mode = "Markov",
                                 method = "Nelder-Mead") {
  
  input <- .ReshapeInput(W,
                         formula,
                         transition.mode,
                         X)
  
  kNuInf <- .GetNuInf()
  beta <- .InitBetaPar()
  sigma <- .InitSigmaPar(W)
  tau <- .InitTauPar(W, sigma, input$p_len)
  
  fr <- function(data) { 
    tau.internal <- exp(data[1 : input$p_len])
    eta.rho.internal <- exp(data[input$p_len + 1])
    gamma.rho.internal <- 0.5 * exp(data[input$p_len + 2]) / 
      (1 + exp(data[input$p_len + 2]))      
    sigma.internal <- exp(data[input$p_len + 3])
    alpha.internal <- data[input$p_len + 4]
    # If the nu_par is set to infinity, fit the model with stationary noise;
    # otherwise take in the initial nu value for EB.
    nu.internal <- ifelse(nu != Inf, exp(data[input$p_len + 5]), kNuInf)
    
    output <- fitGroveML(W$D, 
                         input$XX, 
                         p = input$p, 
                         tau_par = tau.internal, 
                         eta_par = c(eta.rho.internal, 
                                     rep(eta.kappa, input$p_len - 1)),
                         gamma_par = c(gamma.rho.internal, 
                                       rep(gamma.kappa, input$p_len - 1)),
                         init_state = input$init_state,
                         nu = nu.internal, 
                         sigma0 = sigma.internal, 
                         alpha = alpha.internal, 
                         beta = beta, 
                         transition_mode = input$transition.mode)
    if (verbose == TRUE) {
      print(- output$marginal_likelihood)
    }    
    return(- output$marginal_likelihood)
  } 
  # Organize starting values in a vector.
  eb.par <- c(log(tau),
              log(eta.rho), 
              log(gamma.rho / (0.5 - gamma.rho)), 
              log(sigma),
              alpha)
  # Run EB.
  empirical_bayes <- ifelse(nu != Inf,
                            optim(par = c(eb.par, log(nu)), 
                                  fn = fr,  
                                  method = method),
                            optim(par = eb.par, 
                                  fn = fr,  
                                  method = method))
  # Define final values for EB.
  tau_par_final <- exp(empirical_bayes[[1]][1 : input$p_len])
  eta_par_final <- c(exp(empirical_bayes[[1]][input$p_len + 1]), 
                     rep(eta.kappa, input$p_len - 1))
  gamma_par_final <- c(0.5 * exp(empirical_bayes[[1]][input$p_len + 2]) / 
                         (1 + exp(empirical_bayes[[1]][input$p_len + 2])),
                       rep(gamma.kappa, input$p_len - 1))
  sigma_par_final <- exp(empirical_bayes[[1]][input$p_len + 3])
  alpha_par_final <- empirical_bayes[[1]][input$p_len + 4]
  beta_par_final <- beta
  nu_par_final  <- ifelse(nu != Inf, 
                          exp(empirical_bayes[[1]][input$p_len + 5]),
                          kNuInf)
  # Run model based on optimal hyperparameters.
  ans <- .DenoiseFinalRun(W, 
                          input$XX, 
                          X,
                          formula,
                          input$C_hat,
                          input$p,
                          tau_par_final,
                          eta_par_final,
                          gamma_par_final,
                          input$init_state,
                          nu_par_final,
                          sigma_par_final,
                          alpha_par_final,
                          beta_par_final,
                          n.samples,
                          input$transition.mode)
  return(ans)  
}  