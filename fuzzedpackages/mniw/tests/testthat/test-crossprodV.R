#require(testthat)
## library(mniw)
source("mniw-testfunctions.R")
context("Cross-Products")

tol <- 1e-6

test_that("Cross-product X'VX is same in C++ as R", {
  case.par <- expand.grid(p = 1:4, q = 1:4,
                          X = c("single", "multi"),
                          V = c("single", "multi"),
                          inverse = c(TRUE, FALSE), stringsAsFactors = FALSE)
  ncases <- nrow(case.par)
  n <- 10
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    args <- list(X = list(p = p, q = q, rtype = cp$X, vtype = "matrix"),
                 V = list(p = p, rtype = cp$V, vtype = "matrix"))
    args <- get_args(n = n, args = args, drop = FALSE)
    # R test
    ipR <- array(NA, dim = c(q,q,n))
    for(jj in 1:n) {
      if(cp$inverse) {
        VX <- solve(a = args$R$V[[jj]], b = args$R$X[[jj]])
      } else {
        VX <- args$R$V[[jj]] %*% args$R$X[[jj]]
      }
      ipR[,,jj] <- crossprod(args$R$X[[jj]], VX)
    }
    # C++ test
    ipcpp <- do.call(crossprodV, args = c(args$cpp, list(inverse = cp$inverse)))
    if(all_single(cp)) {
      ipcpp <- array(ipcpp, c(q,q,n))
    }
    expect_equal(ipR, ipcpp, tolerance = tol)
  }
})

test_that("Cross-product X'VY is same in C++ as R", {
  case.par <- expand.grid(p = 1:4, q = 1:4, r = 1:4,
                          X = c("single", "multi"),
                          V = c("single", "multi"),
                          Y = c("single", "multi"),
                          inverse = c(TRUE, FALSE), stringsAsFactors = FALSE)
  ncases <- nrow(case.par)
  n <- 10
  for(ii in 1:nrow(case.par)) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    r <- cp$r
    args <- list(X = list(p = p, q = q, rtype = cp$X, vtype = "matrix"),
                 V = list(p = p, rtype = cp$V, vtype = "matrix"),
                 Y = list(p = p, q = r, rtype = cp$Y, vtype = "matrix"))
    args <- get_args(n = n, args = args, drop = FALSE)
    # R test
    ipR <- array(NA, dim = c(q,r,n))
    for(jj in 1:n) {
      if(cp$inverse) {
        VY <- solve(a = args$R$V[[jj]], b = args$R$Y[[jj]])
      } else {
        VY <- args$R$V[[jj]] %*% args$R$Y[[jj]]
      }
      ipR[,,jj] <- crossprod(args$R$X[[jj]], VY)
    }
    # C++ test
    ipcpp <- do.call(crossprodV, args = c(args$cpp, list(inverse = cp$inverse)))
    if(all_single(cp)) {
      ipcpp <- array(ipcpp, c(q,r,n))
    }
    expect_equal(ipR, ipcpp, tolerance = tol)
  }
})
