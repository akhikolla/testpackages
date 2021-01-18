.ReshapeInput <- function(W,
                          formula,
                          transition.mode,
                          X) {
  # Reshape input parameters for functions Denoise and FAnova. 

  if (class(W) != "DWT")  {
    stop("W should be a DWT object.")
  }
  
  if (strsplit(as.character(formula)[2], split = "")[[1]][1] != 1)  {
    stop("formula should include an intercept.")
  }
  
  if (transition.mode == "Markov") {
    transition.mode <- 1
  } else if (transition.mode == "Independent") {
    transition.mode <- 0
  } else {
    transition.mode <- 1
    print("WARNING: Unrecognized transition mode. Default to Markov.")
  }
  
  XX <- model.matrix(formula, X)
  p <- table(attr(XX, "assign"))  
  p_len <- length(p)
  
  frml <- paste("C", paste(formula, sep = "", collapse = " "), sep = " ")
  C_hat <- lm(frml, data = cbind(C = W$C, X))$coefficients
  init_state <- rep(0, p_len)

  return(list(XX= XX,
              p = p,
              p_len = p_len,
              C_hat = C_hat,
              init_state = init_state,
              transition.mode = transition.mode))
}

.InitSigmaPar <- function(W) {
  return(mad(tail(t(W$D), n = 10)))
}

.InitEtaRhoPar <- function() {
  return(2)
}

.InitEtaKappaPar <- function(n) {
  return(0.1)
}

.InitTauPar <- function(W, sigma, n) {
  return(rep((sd(head(t(W$D), n = 10)) / sigma) ^ 2, n))
}

.InitBetaPar <- function() {
  return(1)
}

.InitGammaRhoPar <- function() {
  return(0.3)
}

.InitAlphaPar <- function() {
  return(0.5)
}

.InitGammaKappaPar <- function(n) {
  return(0.3)
}

.GetNuInf <- function() {
  return(999999999)
}

.DenoiseFinalRun <- function(W, 
                             XX, 
                             X,
                             formula,
                             C_hat,
                             p,
                             tau_par,
                             eta_par,
                             gamma_par,
                             init_state,
                             nu_par,
                             sigma_par,
                             alpha_par,
                             beta_par,
                             n.samples,
                             transition.mode) {
  ans <- fitGrove(W$D, 
                  XX, 
                  p = p, 
                  tau_par = tau_par, 
                  eta_par = eta_par, 
                  gamma_par = gamma_par, 
                  init_state = init_state,
                  nu = nu_par, 
                  sigma0 = sigma_par, 
                  alpha = alpha_par, 
                  beta = beta_par, 
                  n_samp = n.samples, 
                  transition_mode = transition.mode)
  ans$data$formula <- formula
  ans$data$X <- X
  ans$data$W <- W
  colnames(ans$data$design_matrix) <- colnames(XX)
  ans$C_hat <- C_hat
  ans$hyperparameters <- list(tau = tau_par, 
                              eta = eta_par, 
                              gamma = gamma_par, 
                              sigma = sigma_par, 
                              nu = nu_par, 
                              alpha = alpha_par)
  class(ans) <- "grove"
  return( ans )
}