## library(mniw)
source("mniw-testfunctions.R")
context("Random-Effects Normal Distribution")

tol <- 1e-6

test_that("Random-Effects Normal sampling is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(q = c(1,2,4),
                          x = c("single", "multi"),
                          V = c("none", "single", "multi"),
                          lambda = c("none", "single", "multi"),
                          Sigma = c("none", "single", "multi"),
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  ncases <- nrow(case.par)
  n <- 5 # number of random draws
  test.seed <- sample(1e6, ncases)
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    args <- list(x = list(p = 1, q = q, rtype = cp$x, vtype = "vector"),
                 V = list(q = q, rtype = cp$V, vtype = "matrix"),
                 lambda = list(p = 1, q = q, rtype = cp$lambda, vtype = "vector"),
                 Sigma = list(q = q, rtype = cp$Sigma, vtype = "matrix"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    muR <- matrix(NA, n, q)
    set.seed(test.seed[ii])
    for(jj in 1:n) {
      muR[jj,] <- rmNormRER(y = args$R$x[[jj]],
                            V = args$R$V[[jj]],
                            lambda = args$R$lambda[[jj]],
                            A = args$R$Sigma[[jj]])
    }
    # C++ test
    set.seed(test.seed[ii])
    mucpp <- do.call(rRxNorm, args = c(args$cpp, list(n = n)))
    mx <- arDiff(muR, mucpp)
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      ## expect_equal(mx, 0, tolerance = tol)
      expect_Rcpp_equal("rRxNorm", ii, mx, tolerance = tol)
    }
  }
})
