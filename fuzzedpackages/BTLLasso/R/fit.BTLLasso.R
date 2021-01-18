fit.BTLLasso <- function(response, design, penalty, lambda, k, 
  m, control, trace) {
  
  adaptive <- control$adaptive
  norm <- control$norm
  epsilon <- control$epsilon
  lambda2 <- control$lambda2
  c <- control$c
  
  #### initialize for estimation
  coefs <- matrix(0, nrow = length(lambda), ncol = ncol(design))
  colnames(coefs) <- colnames(design)
  df <- c()
  start <- NULL
  
  ## calculate adaptive if needed
  if (adaptive) {
    if (k > 2) {
      m0 <- cum.fit.Cpp(response, design, kat = k, epsilon = epsilon, 
        start = start, penalty = penalty, lambda = 0, 
        max.iter = 100, norm = norm, adaptive = NULL, 
        control = list(c = c, gama = 20, index = 1), 
        m = m, hat.matrix = FALSE, lambda2 = lambda2)
    } else {
      m0 <- bin.fit.Cpp(response, design, kat = k, epsilon = epsilon, 
        start = start, penalty = penalty, lambda = 0, 
        max.iter = 100, norm = norm, adaptive = NULL, 
        control = list(c = c, gama = 20, index = 1), 
        m = m, hat.matrix = FALSE, lambda2 = lambda2)
    }
    adaptive <- m0$coef
    if (any(is.nan(adaptive))) {
      stop("Unpenalized parameters for adaptive weights can not be estimated! Please increase lambda2 in 
       the control argument or set adaptive = FALSE!")
    }
  } else {
    adaptive <- NULL
  }
  
  
  
  ## start estimation
  for (i in seq_along(lambda)) {
    if (trace) {
      cat("lambda =", lambda[i], "\n")
    }
    if (k > 2) {
      m1 <- cum.fit.Cpp(response, design, kat = k, epsilon = epsilon, 
        start = start, penalty = penalty, lambda = lambda[i], 
        max.iter = 100, norm = norm, adaptive = adaptive, 
        control = list(c = c, gama = 20, index = 1), 
        m = m, hat.matrix = FALSE, lambda2 = lambda2)
    } else {
      m1 <- bin.fit.Cpp(response, design, kat = k, epsilon = epsilon, 
        start = start, penalty = penalty, lambda = lambda[i], 
        max.iter = 100, norm = norm, adaptive = adaptive, 
        control = list(c = c, gama = 20, index = 1), 
        m = m, hat.matrix = FALSE, lambda2 = lambda2)
    }
    coefs[i, ] <- m1$coef
    start <- m1$coef
    df[i] <- m1$df
    
  }
  
  return(list(coefs = coefs, df = df))
  
  
  
}