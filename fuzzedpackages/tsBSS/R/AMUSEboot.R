AMUSE_boot_teststatistic <- function(X, k, tau) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  Yt <- Y[1:(n - tau), ]
  Yti <- Y[(1 + tau):n, ]
  Ri <- crossprod(Yt, Yti)/nrow(Yt)
  R <- (Ri + t(Ri))/2
 
  EVDvalues <- as.vector(.Call("EIGEN", R, PACKAGE = "tsBSS")$values)
  ds <- EVDvalues^2
  ord <- order(ds, decreasing = TRUE)
  DS <- ds[ord]
  TEST.STATISTIC.X <- mean(DS[(k + 1):p])
  TEST.STATISTIC.X
}

AMUSE_boot_par <- function(Z1, Winv, MEAN, k, tau) {
  n <- nrow(Z1)
  p <- length(MEAN)
  Z2tilde <- matrix(rnorm(n*(p - k)), nrow = n)
  Zstar <- cbind(Z1, Z2tilde)
  Xstar <- tcrossprod(Zstar, Winv)
  Xstar <- sweep(Xstar, 2, MEAN, "+")
  
  TEST.STATISTIC.Xstar <- AMUSE_boot_teststatistic(Xstar, k, tau)
  TEST.STATISTIC.Xstar
}


AMUSE_boot_nonpar2 <- function(Z1, Z2, Winv, MEAN, k, tau) {
  n <- nrow(Z1)
  p <- length(MEAN)
  Z2tilde <- apply(Z2, 2, sample, replace = TRUE)
  Zstar <- cbind(Z1, Z2tilde)
  Xstar <- tcrossprod(Zstar, Winv)
  Xstar <- sweep(Xstar, 2, MEAN, "+")
  TEST.STATISTIC.Xstar <- AMUSE_boot_teststatistic(Xstar, k, tau)
  TEST.STATISTIC.Xstar
}

AMUSE_boot_nonpar1 <- function(Z1, Z2, Winv, MEAN, k, tau) {
  n <- nrow(Z1)
  Z2tilde <- matrix(sample(Z2, replace = TRUE), nrow = n)
  Zstar <- cbind(Z1, Z2tilde)
  Xstar <- tcrossprod(Zstar, Winv)
  Xstar <- sweep(Xstar, 2, MEAN, "+")
  TEST.STATISTIC.Xstar <- AMUSE_boot_teststatistic(Xstar, k, tau)
  TEST.STATISTIC.Xstar
}

AMUSE_boot_nonpar3 <- function(Z1, Z2, Winv, MEAN, k, tau) {
  n <- nrow(Z1)
  p <- length(MEAN)
  ind <- sample(1:n, replace = TRUE)
  Z2tilde <- Z2[ind, ]
  Zstar <- cbind(Z1, Z2tilde)
  Xstar <- tcrossprod(Zstar, Winv)
  Xstar <- sweep(Xstar, 2, MEAN, "+")
  TEST.STATISTIC.Xstar <- AMUSE_boot_teststatistic(Xstar, k, tau)
  TEST.STATISTIC.Xstar
}

AMUSE_boot_p <- function(X, k, tau, n.boot = 200, ncores = NULL, iseed = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  MEAN <- prep$MEAN

  Yt <- Y[1:(n - tau), ]
  Yti <- Y[(1 + tau):n, ]
  Ri <- crossprod(Yt, Yti)/nrow(Yt)
  R <- (Ri + t(Ri))/2
  
  EVD <- .Call("EIGEN", R, PACKAGE = "tsBSS")
  EVDvalues <- as.vector(EVD$values)
  ds <- EVDvalues^2
  ord <- order(ds, decreasing = TRUE)
  DS <- ds[ord]
  D <- EVDvalues[ord]
  W <- crossprod(EVD$vectors, prep$COV.sqrt.i)[ord, ]
  TEST.STATISTIC.X <- mean(DS[(k + 1):p])
  names(TEST.STATISTIC.X) = "T"
  PARAMETER <- n.boot
  names(PARAMETER) <- c("replications")
  Z <- tcrossprod(prep$X.C, W)
  Z1 <- Z[, 0:k, drop = FALSE]
  Winv <- solve(W)
  
  if(!is.null(ncores) && ncores > 1){
  
  if(is.null(iseed)){
    if (exists(".Random.seed", envir = globalenv())){
      oldseed <- get(".Random.seed", envir = globalenv())
      rm(.Random.seed, envir = globalenv())
      on.exit(assign(".Random.seed", oldseed, envir = globalenv()))
    }
  }
  
  type <- "PSOCK"
  cl <- makeCluster(ncores, type = type)
  clusterExport(cl, c("n.boot", "Z1", "MEAN", "k", "tau", "iseed", "AMUSE_boot_par", "AMUSE_boot_teststatistic"), envir = environment())
  clusterSetRNGStream(cl = cl, iseed = iseed)
  TEST.STATISTICS.Xstar <- parSapply(cl, 1:n.boot, function(i, ...) {
    t(AMUSE_boot_par(Z1, Winv, MEAN, k, tau))})
  stopCluster(cl)
  } else {
  TEST.STATISTICS.Xstar <- t(replicate(n.boot, AMUSE_boot_par(Z1, Winv, MEAN, k, tau)))
}
  
  Z <- ts(Z)
  colnames(Z) <- paste0("Series", 1:p)
  PVAL <- (sum(TEST.STATISTIC.X < TEST.STATISTICS.Xstar) + 1)/(n.boot + 1)
  ALTERNATIVE <- paste0("the last ", p - k, " components are not white noise")
  RES <- list(statistic = n * TEST.STATISTIC.X, p.value = PVAL, 
              parameter = PARAMETER, alternative = ALTERNATIVE, k = k, 
              W = W, S = Z, D = D, MU = MEAN, tau = tau)
  RES
}


AMUSE_boot_np2 <- function(X, k, tau, n.boot = 200, ncores = NULL, iseed = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  MEAN <- prep$MEAN
  
  Yt <- Y[1:(n - tau), ]
  Yti <- Y[(1 + tau):n, ]
  Ri <- crossprod(Yt, Yti)/nrow(Yt)
  R <- (Ri + t(Ri))/2
  
  EVD <- .Call("EIGEN", R, PACKAGE = "tsBSS")
  EVDvalues <- as.vector(EVD$values)
  ds <- EVDvalues^2
  ord <- order(ds, decreasing = TRUE)
  DS <- ds[ord]
  D <- EVDvalues[ord]
  W <- crossprod(EVD$vectors, prep$COV.sqrt.i)[ord, ]
  TEST.STATISTIC.X <- mean(DS[(k + 1):p])
  names(TEST.STATISTIC.X) = "T"
  PARAMETER <- n.boot
  names(PARAMETER) <- c("replications")
  Z <- tcrossprod(prep$X.C, W)
  Z1 <- Z[, 0:k, drop = FALSE]
  Z2 <- Z[, (k + 1):p, drop = FALSE]
  Winv <- solve(W)
  
  
  if (!is.null(ncores) && ncores > 1) {
    
    if (is.null(iseed)) {
      if (exists(".Random.seed", envir = globalenv())) {
        oldseed <- get(".Random.seed", envir = globalenv())
        rm(.Random.seed, envir = globalenv())
        on.exit(assign(".Random.seed", oldseed, envir = globalenv()))
      }
    }
    
    type <- "PSOCK"
    cl <- makeCluster(ncores, type = type)
    clusterExport(cl, c("n.boot", "Z1", "MEAN", "Z2", "Winv", "k", "tau", "iseed", "AMUSE_boot_nonpar2", "AMUSE_boot_teststatistic"), envir = environment())
    clusterSetRNGStream(cl = cl, iseed = iseed)
    TEST.STATISTICS.Xstar <- parSapply(cl, 1:n.boot, function(i, ...) {
      t(AMUSE_boot_nonpar2(Z1, Z2, Winv, MEAN, k, tau))})
    stopCluster(cl)
  } else {
    TEST.STATISTICS.Xstar <- t(replicate(n.boot, AMUSE_boot_nonpar2(Z1, Z2, Winv, MEAN, k, tau)))
  }
  Z <- ts(Z)
  colnames(Z) <- paste0("Series", 1:p)
  PVAL <- (sum(TEST.STATISTIC.X < TEST.STATISTICS.Xstar) + 1)/(n.boot + 1)
  ALTERNATIVE <- paste0("the last ", p - k, " components are not white noise")
  RES <- list(statistic = n * TEST.STATISTIC.X, p.value = PVAL, 
              parameter = PARAMETER, alternative = ALTERNATIVE, k = k, 
              W = W, S = Z, D = D, MU = MEAN, tau = tau)
  RES
}

AMUSE_boot_np1 <- function(X, k, tau, n.boot = 200, ncores = NULL, iseed = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  MEAN <- prep$MEAN
  
  Yt <- Y[1:(n - tau), ]
  Yti <- Y[(1 + tau):n, ]
  Ri <- crossprod(Yt, Yti)/nrow(Yt)
  R <- (Ri + t(Ri))/2
  
  EVD <- .Call("EIGEN", R, PACKAGE = "tsBSS")
  EVDvalues <- as.vector(EVD$values)
  ds <- EVDvalues^2
  ord <- order(ds, decreasing = TRUE)
  DS <- ds[ord]
  D <- EVDvalues[ord]
  W <- crossprod(EVD$vectors, prep$COV.sqrt.i)[ord, ]
  TEST.STATISTIC.X <- mean(DS[(k + 1):p])
  names(TEST.STATISTIC.X) = "T"
  PARAMETER <- n.boot
  names(PARAMETER) <- c("replications")
  Z <- tcrossprod(prep$X.C, W)
  Z1 <- Z[, 0:k, drop = FALSE]
  Z2 <- Z[, (k + 1):p, drop = FALSE]
  Winv <- solve(W)
  
  if (!is.null(ncores) && ncores > 1) {
    if (is.null(iseed)) {
      if (exists(".Random.seed", envir = globalenv())) {
      oldseed <- get(".Random.seed", envir = globalenv())
      rm(.Random.seed, envir = globalenv())
      on.exit(assign(".Random.seed", oldseed, envir = globalenv()))
      }
    }
    type <- "PSOCK"
    cl <- makeCluster(ncores, type = type)
    clusterExport(cl, c("n.boot", "Z1", "MEAN", "Z2", "Winv", "k", "tau", "iseed", "AMUSE_boot_nonpar1", "AMUSE_boot_teststatistic"), envir = environment())
    clusterSetRNGStream(cl = cl, iseed = iseed)
    TEST.STATISTICS.Xstar <- parSapply(cl, 1:n.boot, function(i, ...) {
      t(AMUSE_boot_nonpar1(Z1, Z2, Winv, MEAN, k, tau))})
    stopCluster(cl)
  } else {
    TEST.STATISTICS.Xstar <- t(replicate(n.boot, AMUSE_boot_nonpar1(Z1, Z2, Winv, MEAN, k, tau)))
  }
  
  Z <- ts(Z)
  colnames(Z) <- paste0("Series", 1:p)
  PVAL <- (sum(TEST.STATISTIC.X < TEST.STATISTICS.Xstar) + 1)/(n.boot + 1)
  ALTERNATIVE <- paste0("the last ", p - k, " components are not white noise")
  RES <- list(statistic = n * TEST.STATISTIC.X, p.value = PVAL, 
              parameter = PARAMETER, alternative = ALTERNATIVE, k = k, 
              W = W, S = Z, D = D, MU = MEAN, tau = tau)
  RES
}



AMUSE_boot_np3 <- function (X, k, tau, n.boot = 200, ncores = NULL, iseed = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  MEAN <- prep$MEAN
  
  Yt <- Y[1:(n - tau), ]
  Yti <- Y[(1 + tau):n, ]
  Ri <- crossprod(Yt, Yti)/nrow(Yt)
  R <- (Ri + t(Ri))/2

  EVD <- .Call("EIGEN", R, PACKAGE = "tsBSS")
  EVDvalues <- as.vector(EVD$values)
  ds <- EVDvalues^2
  ord <- order(ds, decreasing = TRUE)
  DS <- ds[ord]
  D <- EVDvalues[ord]
  W <- crossprod(EVD$vectors, prep$COV.sqrt.i)[ord, ]
  TEST.STATISTIC.X <- mean(DS[(k + 1):p])
  names(TEST.STATISTIC.X) = "T"
  PARAMETER <- n.boot
  names(PARAMETER) <- c("replications")
  Z <- tcrossprod(prep$X.C, W)
  Z1 <- Z[, 0:k, drop = FALSE]
  Z2 <- Z[, (k+1):p, drop = FALSE]
  Winv <- solve(W)
  
  if(!is.null(ncores) && ncores > 1){
    if(is.null(iseed)){
      if (exists(".Random.seed", envir = globalenv())){
        oldseed <- get(".Random.seed", envir = globalenv())
        rm(.Random.seed, envir = globalenv())
        on.exit(assign(".Random.seed", oldseed, envir = globalenv()))
      }
    }
    
    type <- "PSOCK"
    cl <- makeCluster(ncores, type = type)
    clusterExport(cl, c("n.boot", "Z1", "MEAN", "Z2", "Winv", "k", "tau", "iseed", "AMUSE_boot_nonpar3", "AMUSE_boot_teststatistic"),envir=environment())
    clusterSetRNGStream(cl = cl, iseed = iseed)
    TEST.STATISTICS.Xstar <- parSapply(cl, 1:n.boot, function(i, ...) {
      t(AMUSE_boot_nonpar3(Z1, Z2, Winv, MEAN, k, tau))})
    stopCluster(cl)
  } else {
    TEST.STATISTICS.Xstar <- t(replicate(n.boot, AMUSE_boot_nonpar3(Z1, Z2, Winv, MEAN, k, tau)))
  }
  
  Z <- ts(Z)
  colnames(Z) <- paste0("Series", 1:p)
  PVAL <- (sum(TEST.STATISTIC.X < TEST.STATISTICS.Xstar) + 1)/(n.boot + 1)
  ALTERNATIVE <- paste0("the last ", p - k, " components are not white noise")
  RES <- list(statistic = n * TEST.STATISTIC.X, p.value = PVAL, 
              parameter = PARAMETER, alternative = ALTERNATIVE, k = k, 
              W = W, S = Z, D = D, MU = MEAN, tau = tau)
  RES
}

# Method AMUSEboot
AMUSEboot <- function(X, ...) UseMethod("AMUSEboot")

AMUSEboot.default <- function (X, k, tau = 1, n.boot = 200, s.boot = c("p", "np1", "np2", "np3"),
                               ncores = NULL, iseed = NULL, ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  DNAME <- deparse(substitute(X))
  k <- abs(as.integer(k))
  tau <- abs(as.integer(tau))
  n.boot <- abs(as.integer(n.boot))
  X <- as.matrix(X)
  s.boot <- match.arg(s.boot)
  
  switch(s.boot,
         p = {
           RES <- AMUSE_boot_p(X = X, k = k, tau = tau, n.boot = n.boot, ncores = ncores, iseed = iseed)
           METHOD <- "ICA subwhite noise bootstrapping test using AMUSE and strategy p"  
         },
         np1 = {
           RES <- AMUSE_boot_np1(X = X, k = k, tau = tau, n.boot = n.boot, ncores = ncores, iseed = iseed)
           METHOD <- "ICA subwhite noise bootstrapping test using AMUSE and strategy np1"
         },
         np2 = {
           RES <- AMUSE_boot_np2(X = X, k = k, tau = tau, n.boot = n.boot, ncores = ncores, iseed = iseed)
           METHOD <- "ICA subwhite noise bootstrapping test using AMUSE and strategy np2"
         },
         np3 = {
           RES <- AMUSE_boot_np3(X = X, k = k, tau = tau, n.boot = n.boot, ncores = ncores, iseed = iseed)
           METHOD <- "ICA subwhite noise bootstrapping test using AMUSE and strategy np3"
         }
         )
  RES <- c(RES, method = METHOD, data.name = DNAME, s.boot = s.boot)
  class(RES) <- c("ictest", "htest")
  RES
}

AMUSEboot.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEboot.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

AMUSEboot.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEboot.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

AMUSEboot.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEboot.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}
