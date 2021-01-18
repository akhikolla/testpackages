#require(testthat)
## library(mniw)
source("mniw-testfunctions.R")
context("Matrix Normal Distribution")

tol <- 1e-6

test_that("Matrix Normal density is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(1,2,4), q = c(1,2,3),
                          X = c("single", "multi"),
                          Lambda = c("none", "single", "multi"),
                          SigmaR = c("none", "single", "multi"),
                          SigmaC = c("none", "single", "multi"),
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  ncases <- nrow(case.par)
  n <- 10
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    args <- list(X = list(p = p, q = q, rtype = cp$X, vtype = "matrix"),
                 Lambda = list(p = p, q = q, rtype = cp$Lambda, vtype = "matrix"),
                 SigmaR = list(p = p, rtype = cp$SigmaR, vtype = "matrix"),
                 SigmaC = list(q = q, rtype = cp$SigmaC, vtype = "matrix"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    llR <- rep(NA, n)
    for(jj in 1:n) {
      llR[jj] <- dMNormR(X = args$R$X[[jj]],
                         Lambda = args$R$Lambda[[jj]],
                         SigmaU = args$R$SigmaR[[jj]],
                         SigmaV = args$R$SigmaC[[jj]],
                         log = TRUE)
    }
    # C++ test
    llcpp <- do.call(dMNorm, args = c(args$cpp, list(log = TRUE)))
    # if all inputs to c++ are single it returns only one value
    if(all_single(cp)) {
      llcpp <- rep(llcpp, n)
    }
    mx <- abs(llR-llcpp)
    mx <- min(max(mx), max(mx/abs(llR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      ## expect_equal(mx, 0, tolerance = tol)
      expect_Rcpp_equal("dMNorm", ii, mx, tolerance = tol)
    }
  }
})

test_that("Matrix Normal simulation is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(1,2,4), q = c(1,2,3),
                          Lambda = c("none", "single", "multi"),
                          SigmaR = c("none", "single", "multi"),
                          SigmaC = c("none", "single", "multi"),
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  # remove cases where dimensions can't be identified
  case.par <- case.par[!with(case.par, {
    Lambda == "none" & ((SigmaR == "none") | (SigmaC == "none"))}),]
  ncases <- nrow(case.par)
  rownames(case.par) <- 1:ncases
  n <- 10
  TestSeed <- sample(1e6, ncases)
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:nrow(case.par)) {
    set.seed(TestSeed[ii]) # seed
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    args <- list(Lambda = list(p = p, q = q, rtype = cp$Lambda, vtype = "matrix"),
                 SigmaR = list(p = p, rtype = cp$SigmaR, vtype = "matrix"),
                 SigmaC = list(q = q, rtype = cp$SigmaC, vtype = "matrix"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    set.seed(TestSeed[ii]) # seed
    XR <- array(NA, dim = c(p,q,n))
    for(jj in 1:n) {
      XR[,,jj] <- rMNormR(Lambda = args$R$Lambda[[jj]],
                          SigmaU = args$R$SigmaR[[jj]],
                          SigmaV = args$R$SigmaC[[jj]])
    }
    # C++ test
    set.seed(TestSeed[ii]) # seed
    Xcpp <- do.call(rMNorm, args = c(args$cpp, list(n = n)))
    mx <- abs(range(XR - Xcpp))
    mx <- min(max(mx), max(mx/abs(XR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      ## expect_equal(mx, 0, tolerance = tol)
      expect_Rcpp_equal("rMNorm", ii, mx, tolerance = tol)
    }
  }
})
