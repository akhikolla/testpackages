checkArgs1 <- function(lower, upper, mean, sigma, uselog2, nu = NULL)
{
  if(!is.numeric(lower) || !is.vector(lower))
    stop('lower should be numeric vector')
  if(!is.numeric(upper) || !is.vector(upper))
    stop('upper should be numeric vector')
  if(!is.numeric(mean) || !is.vector(mean))
    stop('mean should be numeric vector')
  rec <- cbind(lower, upper, mean) # <--> recycling to same length
  lower <- rec[,"lower"]
  upper <- rec[,"upper"]
  if (!all(lower <= upper))
    stop("at least one element of ", sQuote("lower"), " is larger than ",
         sQuote("upper"))
  mean <- rec[,"mean"]
  if (anyNA(mean))
    stop("mean contains NA")
  n <- length(lower)
  buildSigma <- is.null(sigma)
  if(!buildSigma) # check sigma
  {
    if(ncol(sigma) != n || nrow(sigma) != n)
      stop("Row number or column number of sigma is incorrect")
    if(anyNA(sigma))
      stop("sigma cannot contain NA")
  }
  if(!is.logical(uselog2))
    stop("uselog2 should be logical")
  if(!is.null(nu)) # check nu
  {
    if(!is.numeric(nu))
      stop("nu should be numeric")
    if(is.na(nu))
      stop("nu cannot be NA")
    if(length(nu) != 1)
      stop("nu should be of length 1")
    if(nu <= 0)
      stop("nu should be positive")
  }
  list(lower = lower, upper = upper, mean = mean, sigma = sigma, 
       uselog2 = uselog2, nu = nu, buildSigma = buildSigma)
}

checkArgs2 <- function(n, geom, kernelType, para)
{
  if(is.null(geom))
    stop("geom cannot be NULL when sigma is not given")
  if(nrow(geom) != n)
    stop("Row number of geom is incorrect")
  if(ncol(geom) != 2 || any(geom > 1) || any(geom < 0))
    stop("Currently geom is required to be within the 2D unit square")
  if(is.null(kernelType))
    stop("kernelType cannot be NULL when sigma is not given")
  if(is.null(para))
    stop("para cannot be NULL when sigma is not given")
  if(!is.character(kernelType))
    stop("kernelType should be a character string")
  if(tolower(kernelType) == "matern")
  {
    if(!is.numeric(para) || !is.vector(para))
      stop("para should be a numeric vector")
    if(length(para) != 4)
      stop("para for the matern kernel should be of length 4, denoting 
           scale, range, smoothness, and nugget")
    if(anyNA(para))
      stop("para contains NA")
    if(any(para[1 : 3] <= 0) || para[4] < 0)
      stop("para[1:3] for the matern kernel should all be positive and 
           para[4] (nugget) should be non-negative")
  }
  # else if <- Other cov kernels
  else
  {
    stop("Unsupported kernelType")
  }
}

checkN <- function(N)
{
  if(length(N) != 1)
    stop("N should be of length 1")
  if(!is.numeric(N))
    stop("N should be numeric")
  if(N < 100)
    stop("Monte Carlo sample size N should be at least 100")
}

checkTLR <- function(n, m, epsl)
{
  if(length(m) != 1)
    stop("m should be of length 1")
  if(!is.numeric(m))
    stop("m should be numeric")
  if(m < 4)
    stop("Block size m should be at least 4")
  if(m >= n)
    stop("Block size m should be smaller than the problem dimension")
  
  if(length(epsl) != 1)
    stop("Truncation level epsl should be of length 1")
  if(!is.numeric(epsl))
    stop("Truncation level epsl should be numeric")
  if(epsl <= 0)
    stop("Truncation level epsl should be positive")
}

GenzBretz <- function(N = 499) {
  structure(list(N = N), class = "GenzBretz")
}

TLRQMC <- function(N = 499, m = 64, epsl = 1e-4) {
  structure(list(N = N, m = m, epsl = epsl), class = "TLRQMC")
}

pmvn <- function(lower = -Inf, upper = Inf, mean = 0, sigma = NULL, uselog2 = FALSE, 
                 algorithm = GenzBretz(), ...)
{
  # check arguments
  cArgs <- checkArgs1(lower, upper, mean, sigma, uselog2)
  if(cArgs$buildSigma)
  {
    addArgs <- list(...)
    n <- length(lower)
    checkArgs2(n, addArgs$geom, addArgs$kernelType, addArgs$para)
    if(tolower(addArgs$kernelType) == "matern")
      addArgs$kernelType <- 1
    #else if
  }
  checkN(algorithm$N) # check number of samples
  if(class(algorithm) == "TLRQMC")
  {
    n <- length(lower)
    checkTLR(n, algorithm$m, algorithm$epsl)
  }
  # pass to the internal function
  cArgs$lower <- cArgs$lower - cArgs$mean
  cArgs$upper <- cArgs$upper - cArgs$mean
  if(class(algorithm) == "GenzBretz")
  {
    if(cArgs$buildSigma)
    {
      RET <- mvn_internal2(cArgs$lower, cArgs$upper, addArgs$geom, 
                    addArgs$kernelType, addArgs$para[1:3], 
                    addArgs$para[4], cArgs$uselog2, algorithm$N)
    }else
    {
      RET <- mvn_internal(cArgs$lower, cArgs$upper, cArgs$sigma, 
                          cArgs$uselog2, algorithm$N)
    }
  }else if(class(algorithm) == "TLRQMC")
  {
    if(cArgs$buildSigma)
    {
      RET <- tlrmvn_internal2(cArgs$lower, cArgs$upper, addArgs$geom, 
                       addArgs$kernelType, addArgs$para[1:3], 
                       addArgs$para[4], cArgs$uselog2, algorithm$m, 
                       algorithm$epsl, algorithm$N)
    }else
    {
      RET <- tlrmvn_internal(cArgs$lower, cArgs$upper, cArgs$sigma, 
                      cArgs$uselog2, algorithm$m, algorithm$epsl, 
                      algorithm$N)
    }
  }else
    stop("Invalid algorithm class")
  if(!is.null(RET$Error))
    attr(RET$Estimation, "error") <- RET$Error
  return(RET$Estimation)
}

pmvt <- function(lower = -Inf, upper = Inf, delta = 0, df = 1, 
                 sigma = NULL, uselog2 = FALSE, 
                 algorithm = GenzBretz(), 
                 type = "Kshirsagar", ...)
{
  # check arguments
  cArgs <- checkArgs1(lower, upper, delta, sigma, uselog2, df)
  if(cArgs$buildSigma)
  {
    addArgs <- list(...)
    n <- length(lower)
    checkArgs2(n, addArgs$geom, addArgs$kernelType, addArgs$para)
    if(tolower(addArgs$kernelType) == "matern")
      addArgs$kernelType <- 1
    #else if
  }
  checkN(algorithm$N) # check number of samples
  if(class(algorithm) == "TLRQMC")
  {
    n <- length(lower)
    checkTLR(n, algorithm$m, algorithm$epsl)
  }
  # check if type is "shifted"
  if(tolower(type) == "shifted")
  {
    cArgs$lower <- cArgs$lower - cArgs$mean
    cArgs$upper <- cArgs$upper - cArgs$mean
    cArgs$mean <- rep(0, length(cArgs$lower))
  }
  if(class(algorithm) == "GenzBretz")
  {
    if(cArgs$buildSigma)
      RET <- mvt_internal2(cArgs$lower, cArgs$upper, cArgs$mean, cArgs$nu, 
                    addArgs$geom, addArgs$kernelType, addArgs$para[1:3], 
                    addArgs$para[4], cArgs$uselog2, algorithm$N)
    else
      RET <- mvt_internal(cArgs$lower, cArgs$upper, cArgs$mean, cArgs$nu,
                   cArgs$sigma, cArgs$uselog2, algorithm$N)
  }else if(class(algorithm) == "TLRQMC")
  {
    if(cArgs$buildSigma)
      RET <- tlrmvt_internal2(cArgs$lower, cArgs$upper, cArgs$nu, cArgs$mean,
                       addArgs$geom, addArgs$kernelType, addArgs$para[1:3], 
                       addArgs$para[4], cArgs$uselog2, algorithm$m, 
                       algorithm$epsl, algorithm$N)
    else
      RET <- tlrmvt_internal(cArgs$lower, cArgs$upper, cArgs$nu, cArgs$mean, 
                      cArgs$sigma, cArgs$uselog2, algorithm$m, 
                      algorithm$epsl, algorithm$N)
    
  }else
    stop("Invalid algorithm class")
  if(!is.null(RET$Error))
    attr(RET$Estimation, "error") <- RET$Error
  return(RET$Estimation)
}
