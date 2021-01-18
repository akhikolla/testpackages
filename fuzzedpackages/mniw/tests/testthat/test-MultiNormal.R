## library(mniw)
source("mniw-testfunctions.R")
context("Multivariate-Normal Distribution")

tol <- 1e-6

test_that("Multivariate Normal simulation is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(q = c(1,2,4),
                          mu = c("none", "single", "multi"),
                          Sigma = c("none", "single", "multi"),
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  # remove cases where dimensions can't be identified
  case.par <- case.par[!with(case.par, {
    mu == "none" & Sigma == "none"}),]
  ncases <- nrow(case.par)
  n <- 12 # number of random draws
  test.seed <- sample(1e6, ncases)
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    q <- cp$q
    args <- list(mu = list(p = 1, q = q, rtype = cp$mu, vtype = "vector"),
                 Sigma = list(q = q, rtype = cp$Sigma, vtype = "matrix"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    yR <- matrix(NA, n, q)
    set.seed(test.seed[ii])
    for(jj in 1:n) {
      yR[jj,] <- rmNormR(mu = args$R$mu[[jj]],
                         V = args$R$Sigma[[jj]])
    }
    # C++ test
    set.seed(test.seed[ii])
    ycpp <- do.call(rmNorm, args = c(args$cpp, list(n = n)))
    mx <- arDiff(yR, ycpp)
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      ## expect_equal(mx, 0, tolerance = tol)
      expect_Rcpp_equal("rmNorm", ii, mx, tolerance = tol)
    }
  }
})

test_that("Matrix Normal density is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(q = c(1,2,4),
                          x = c("single", "multi"),
                          mu = c("none", "single", "multi"),
                          Sigma = c("none", "single", "multi"),
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  ncases <- nrow(case.par)
  n <- 12 # number of random draws
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    q <- cp$q
    args <- list(x = list(p = 1, q = q, rtype = cp$x, vtype = "vector"),
                 mu = list(p = 1, q = q, rtype = cp$mu, vtype = "vector"),
                 Sigma = list(q = q, rtype = cp$Sigma, vtype = "matrix"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    llR <- rep(NA, n)
    for(jj in 1:n) {
      llR[jj] <- dmNormR(x = args$R$x[[jj]],
                         mu = args$R$mu[[jj]],
                         V = args$R$Sigma[[jj]], log = TRUE)
    }
    # C++ test
    llcpp <- do.call(dmNorm, args = c(args$cpp, list(log = TRUE)))
    if(all_single(cp)) llcpp <- rep(llcpp, n)
    mx <- arDiff(llR, llcpp)
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      ## expect_equal(mx, 0, tolerance = tol)
      expect_Rcpp_equal("dmNorm", ii, mx, tolerance = tol)
    }
  }
})
