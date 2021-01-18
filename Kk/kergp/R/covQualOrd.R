setClass("covOrd",   	
         representation(
           covLevels = "function",
           covLevMat = "matrix",
           hasGrad = "logical",           
           acceptLowerSQRT = "logical",    
           label = "character",
           d = "integer",
           inputNames = "character",     
           nlevels = "integer",                
           levels = "list",          
           parLower = "numeric", 
           parUpper = "numeric",       
           par = "numeric",             
           parN = "integer",              
           kernParNames  = "character",
           k1Fun1 = "function",           ## NEW
           warpFun = "list",              ## NEW
           cov = "integer",               ## NEW   0 : corr, 1 : homo
           parNk1 = "integer",            ## NEW   number of par in kern1fun, usually 0
           parNwarp = "integer",           ## NEW   number of par in warpFun
           k1ParNames = "character",      ## NEW
           warpParNames  = "character",    ## NEW
           ordered = "logical"
         ),
         validity = function(object){
           if (object@parN != object@parNwarp + object@parNk1 + (object@cov == 1)) {
             return("Incorrect number of parameters")
           } else return(TRUE)
         },
         contains = "covQual"
)


covOrd <- function(ordered,
                   k1Fun1 = k1Fun1Matern5_2,
                   warpFun = c("norm", "unorm", "power", "spline1", "spline2"), 
                   cov = c("corr", "homo"),
                   hasGrad = TRUE,
                   inputs = "u",
                   par = NULL,
                   parLower = NULL,
                   parUpper = NULL,
                   label = "Ordinal kernel",
                   intAsChar = TRUE, 
                   ...) {
  
  ## if (intAsChar) {
  ##   warning("With this object, an input of class \"integer\" will be coerced ",
  ##               "into \"character\", not into \"factor\". Use `intAsChar = FALSE` ",
  ##               " to change this behaviour")
  ## }
  
  ordered <- as.ordered(ordered)
  L <- nlevels(ordered)
  cov <- match.arg(cov)
  cov <- switch(cov, corr = 0, homo = 1)
  eps <- 1e-10

  warpFunName <- match.arg(warpFun)
  warpFunName <- switch(warpFunName, 
                         norm = "warpNorm", 
                         unorm = "warpUnorm", 
                         power = "warpPower", 
                         spline1 = "warpSpline1", 
                         spline2 = "warpSpline2")
  warpFun <- get(warpFunName, mode = "list")
  
  warpParNames <- warpFun$parNames
  par <- warpFun$parDefault
  parUpper <- warpFun$parUpper
  parLower <- warpFun$parLower
  hasGrad <- warpFun$hasGrad
  parNk1 <- as.integer(warpFun$isCdf)
  
  if (warpFunName == "warpSpline1") {
    warpParNames <- paste(warpParNames, 1:(L - 1), sep="_")
    par <- rep(par, L - 1)
    parUpper <- rep(parUpper, L - 1)
    parLower <- rep(parLower, L - 1)
  }
  
  if (warpFunName == "warpSpline2") {
    warpParNames <- paste(warpParNames, 1:L, sep="_")
    par <- rep(par, L)
    parUpper <- rep(parUpper, L)
    parLower <- rep(parLower, L)
  }
  
  kernParNames <- warpParNames
  parNwarp <- length(warpParNames)
  
  if (warpFun$isCdf){
    par <- c(par, theta = 1)
    parUpper <- c(parUpper, theta = Inf)
    parLower <- c(parLower, theta = eps)
    kernParNames <- c(kernParNames, "theta")
  }
  
  if (cov == 1){
    par <- c(par, sigma2 = 1)
    parUpper <- c(parUpper, sigma2 = Inf)
    parLower <- c(parLower, sigma2 = eps)
    kernParNames <- c(kernParNames, "sigma2")
  }
  
  parN <- length(par)
  
  if (warpFun$isCdf){
    thisCovLevel <- function(par, lowerSQRT = FALSE, compGrad = FALSE){
      z <- seq(from = 0, to = 1, length.out = L)
      x <- do.call(warpFun$fun, list(z, par[1:parNwarp], L = L))
      K <- matrix(1, nrow = L, ncol = L)
      sI <- symIndices(L)
      Hsym <- (x[sI$i] - x[sI$j]) / par[parNwarp + 1]
      # faster than: outer(x, x, "-")[sI$kL] / par[parNwarp + 1]
      Ksym <- k1Fun1(Hsym)
      K[sI$kL] <- K[sI$kU] <- Ksym
      if (cov == 1) K <- K * par[parNwarp + 2]    # if covariance kernel
      rownames(K) <- colnames(K) <- levels(ordered)
      
      if (compGrad){
        theta <- par[parNwarp + 1]
        dK <- array(0, dim = c(L, L, parN))
        dKsymdx <- attr(Ksym, "der")[, 1]
        for (i in 1:parNwarp){
          dF <- attr(x, "gradient")[, i]
          dFsym <- dF[sI$i] - dF[sI$j]   # faster than: outer(dF, dF, "-")[sI$kL]
          dK[, , i][sI$kL] <- dKsymdx * dFsym / theta
          dK[, , i][sI$kU] <- dK[, , i][sI$kL]  
        }
        dK[, , parNwarp + 1][sI$kL] <- - dKsymdx * Hsym / theta 
        dK[, , parNwarp + 1][sI$kU] <- dK[, , parNwarp + 1][sI$kL]   
        if (cov == 1) {
          dK[, , 1:(parNwarp + 1)] <- dK[, , 1:(parNwarp + 1)] * par[parNwarp + 2]
          dK[, , parNwarp + 2] <- K
        }
        attr(K, "gradient") <- dK
      }
      attr(K, "der") <- NULL
      
      return(K)
    }
  } else {
    thisCovLevel <- function(par, lowerSQRT = FALSE, compGrad = FALSE){
      z <- seq(from = 0, to = 1, length.out = L)
      x <- do.call(warpFun$fun, list(z, par[1:parNwarp], L = L))
      K <- matrix(1, nrow = L, ncol = L)
      sI <- symIndices(L)
      Hsym <- (x[sI$i] - x[sI$j]) 
      Ksym <- k1Fun1(Hsym)
      K[sI$kL] <- K[sI$kU] <- Ksym
      if (cov == 1) K <- K * par[parNwarp + 1]    # if covariance kernel
      rownames(K) <- colnames(K) <- levels(ordered)
      
      if (compGrad){
        dK <- array(0, dim = c(L, L, parN))
        dKsymdx <- attr(Ksym, "der")[, 1]
        for (i in 1:parNwarp){
          dF <- attr(x, "gradient")[, i]
          dFsym <- dF[sI$i] - dF[sI$j]   # faster than: outer(dF, dF, "-")[sI$kL]
          dK[, , i][sI$kL] <- dKsymdx * dFsym 
          dK[, , i][sI$kU] <- dK[, , i][sI$kL]  
        }
        if (cov == 1) {
          dK[, , parNwarp + 1] <- K
        }
        attr(K, "gradient") <- dK
      }
      attr(K, "der") <- NULL
      
      return(K)
    }
  }
  
  new("covOrd", 
      covLevels = thisCovLevel,
      covLevMat = thisCovLevel(par),
      hasGrad = hasGrad,
      acceptLowerSQRT = FALSE,    
      label = "ordKernel",
      d = 1L,
      inputNames = inputs,   
      nlevels = L,            
      levels = list(levels(ordered)),      
      parLower = parLower,          
      parUpper = parUpper, 
      par  = par,              
      parN = parN,  
      kernParNames = kernParNames,
      k1Fun1 = k1Fun1,     
      warpFun = warpFun,
      cov = as.integer(cov),               
      parNk1 = parNk1,            
      parNwarp = parNwarp,     
      k1ParNames = character(0), 
      warpParNames  = warpParNames,
      ordered = TRUE,
      intAsChar = intAsChar,
      ...
  )
  
}

