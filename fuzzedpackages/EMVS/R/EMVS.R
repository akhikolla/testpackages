    
EMVS <- function(Y,
                 X,
                 v0,
                 v1,
                 type = c("betabinomial", "fixed"),
                 independent = TRUE,
                 beta_init,
                 sigma_init = 1,
                 epsilon = 10^(-5),
                 temperature,
                 theta = 0.5,
                 a = 1,
                 b = 1,
                 v1_g,
                 direction = c("backward", "forward", "null"),
                 standardize = TRUE,
                 log_v0 = FALSE) {
  
  type <- match.arg(type)
  direction <- match.arg(direction)
  
  if (missing(beta_init)) beta_init = numeric(ncol(X))
  if (missing(v1_g)) {v1_g = v1}
  if (missing(temperature)){
    temperature = 1
  } else { 
    temperature = 1/temperature
  }
  
  
  if (!is(Y, "numeric"))  { 
    tmp <- try(Y <- as.numeric(Y), silent = TRUE)
    if (is(tmp, "try-error")) 
      stop("Y must numeric or able to be coerced to numeric")
  }  
  
  if (!is(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0 + ., data = X), silent = TRUE)
    if (is(tmp, "try-error")) 
      stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (any(is.na(Y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing Y and X to EMVS")
  
  if(is.unsorted(v0)) stop("Ladder of v0 values must be increasing")
  if (sum(v1 < v0) > 0) stop("v1 has to be larger than v0")
  
  if (!is(beta_init, "numeric")) {
    tmp <- try(beta_init <- as.numeric(beta_init), silent = TRUE)
    if (is(tmp, "try-error")) stop("beta_init must numeric or able to be coerced to numeric")
  }
  
  if (length(beta_init)!=ncol(X))  stop("beta_init has to be of length ncol(X)")
  if (length(a) > 1) stop("a has to be a scalar")
  if (length(b) > 1) stop("b has to be a scalar")
  if (length(v1) > 1) stop("v1 has to be a scalar")
  if (length(theta) > 1) stop("theta has to be a scalar")
  if (length(v1_g) > 1) stop("v1_g has to be a scalar")
  if (length(sigma_init) > 1) stop("sigma_init has to be a scalar")
  if (length(temperature) >1) stop("temperature has to be a scalar")
  

  ## Set up XX, YY
  if (standardize) {
    XX <- t(X)
    center <- apply(XX, 1, mean)
    XX <- XX - center
    scale <- sqrt(apply(XX, 1, function(x) sum(x^2)) / nrow(X))
    XX <- 1/scale * XX
    XX <- t(XX)
    nz <- which(scale > 1e-6)
    if (length(nz) != ncol(XX)) XX <- XX[ ,nz, drop=FALSE]
  } else {
    XX <- X
  }
  
  p <- ncol(XX)
  YY <- Y - mean(Y)
  n <- nrow(XX)

  if(n!=length(Y)) stop("length of Y has to be nrow(X)")

  if(independent == T) {
    res <- .ind_EMVS(YY, XX, v0, v1, type, beta_init, sigma_init, epsilon, temperature, theta, a, b, direction)
    } else {
      res <- .conj_EMVS(YY, XX, v0, v1, type, beta_init, sigma_init, epsilon, temperature, theta, a, b, v1_g, direction)
    }

  ## Unstandardize
  if (standardize) {
    betas <- res$betas
    betas <- t(res$betas)
    bbb <- betas/scale[nz]
    res$betas <- t(bbb)
  }
  
  if(independent == T) {
    res$independent = T
  } else {
    res$independent = F
  }
  
  if(log_v0 == T){
    res$log_v0 = T
  } else {
    res$log_v0 = F
  }
  
  return(res)
         
}


