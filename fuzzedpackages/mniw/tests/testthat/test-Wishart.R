#require(testthat)
## library(mniw)
source("mniw-testfunctions.R")
context("Wishart and Inverse Wishart Distributions")

tol <- 1e-6

test_that("Wishart density is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(q = c(1,2,4),
                          X = c("single", "multi"),
                          Psi = c("none", "single", "multi"),
                          nu = c("single", "multi"),
                          inverse = c(TRUE, FALSE),
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  ncases <- nrow(case.par)
  if(calc.diff) MaxDiff <- rep(NA, ncases)
  n <- 5 # tests per case
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    q <- cp$q
    args <- list(X = list(q = q, rtype = cp$X, vtype = "matrix"),
                 Psi = list(q = q, rtype = cp$Psi, vtype = "matrix"),
                 nu = list(q = q, rtype = cp$nu, vtype = "scalar"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    llR <- rep(NA, n)
    for(jj in 1:n) {
      llR[jj] <- dwishR(X = args$R$X[[jj]],
                        Psi = args$R$Psi[[jj]],
                        nu = args$R$nu[[jj]],
                        inverse = cp$inverse, log = TRUE)
    }
    # C++ test
    llcpp <- do.call(dwishart,
                     args = c(args$cpp,
                              list(inverse = cp$inverse, log = TRUE)))
    # C++ produces single output if all inputs are single
    if(all_single(cp)) llcpp <- rep(llcpp, n)
    mx <- abs(llR - llcpp)
    mx <- min(max(mx), max(mx/abs(llR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      expect_Rcpp_equal("dwishart", ii, mx, tolerance = tol)
      ## expect_equal(mx, 0, tolerance = tol)
    }
  }
})

test_that("Wishart sampling is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(q = c(1,2,4),
                          Psi = c("single", "multi"),
                          nu = c("single", "multi"),
                          inverse = c(TRUE, FALSE),
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  ncases <- nrow(case.par)
  n <- 10
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    q <- cp$q
    args <- list(Psi = list(q = q, rtype = cp$Psi, vtype = "matrix"),
                 nu = list(q = q, rtype = cp$nu, vtype = "scalar"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    test.seed <- sample(1e6, 1)
    set.seed(test.seed)
    llR <- array(NA, dim = c(q,q,n))
    for(jj in 1:n) {
      llR[,,jj] <- rwishR(Psi = args$R$Psi[[jj]],
                          nu = args$R$nu[[jj]],
                          inverse = cp$inverse)
    }
    # C++ test
    set.seed(test.seed)
    llcpp <- do.call(rwishart,
                     args = c(args$cpp, list(n = n, inverse = cp$inverse)))
    mx <- abs(llR - llcpp)
    mx <- min(max(mx), max(mx/abs(llR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      expect_Rcpp_equal("rwishart", ii, mx, tolerance = tol)
    }
  }
})
