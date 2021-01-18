## library(mniw)
source("mniw-testfunctions.R")
context("Matrix-Normal Inverse-Wishart Distribution")

tol <- 1e-6

test_that("MNIW sampling is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(1,2,4), q = c(1,2,3),
                          Lambda = c("none", "single", "multi"),
                          Sigma = c("none", "single", "multi"),
                          Psi = c("none", "single", "multi"),
                          nu = c("single", "multi"),
                          prec = c(TRUE, FALSE),
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  case.par <- case.par[with(case.par, {
    !(Lambda == "none" & ((Sigma == "none") | (Psi == "none")))
  }),]
  ncases <- nrow(case.par)
  rownames(case.par) <- 1:ncases
  n <- 5
  TestSeed <- sample(1e6, ncases)
  if(calc.diff) {
    MaxDiff <- matrix(NA, ncases, 2)
  }
  for(ii in 1:ncases) {
    set.seed(TestSeed[ii])
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    args <- list(Lambda = list(p = p, q = q, rtype = cp$Lambda, vtype = "matrix"),
                 Sigma = list(p = p, rtype = cp$Sigma, vtype = "matrix"),
                 Psi = list(q = q, rtype = cp$Psi, vtype = "matrix"),
                 nu = list(q = q, rtype = cp$nu, vtype = "scalar"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    XR <- array(NA, dim = c(p,q,n))
    VR <- array(NA, dim = c(q,q,n))
    set.seed(TestSeed[ii]) # seed
    for(jj in 1:n) {
      XVR <- rmniwR(Lambda = args$R$Lambda[[jj]],
                    Sigma = args$R$Sigma[[jj]],
                    Psi = args$R$Psi[[jj]],
                    nu = args$R$nu[[jj]], prec = cp$prec)
      XR[,,jj] <- XVR$X
      VR[,,jj] <- XVR$V
    }
    # C++ test
    set.seed(TestSeed[ii]) # seed
    XVcpp <- do.call(rMNIW,
                     args = c(args$cpp, list(n = n, prec = cp$prec)))
    mx <- c(arDiff(XR, XVcpp$X), arDiff(VR, XVcpp$V))
    if(calc.diff) {
      MaxDiff[ii,] <- mx
    } else {
      ## expect_equal(mx, c(0,0), tolerance = tol)
      expect_Rcpp_equal("rMNIW", ii, mx, tolerance = tol)
    }
  }
})


## range(MaxDiff)

## cp <- case.par[ii,]
## p <- cp$p
## q <- cp$q
## Lambda <- if(cp$Lambda == "none") matrix(0,p,q) else rMnorm(p,q)
## Sigma <- if(cp$Sigma == "none") diag(p) else crossprod(rMnorm(p))
## Psi <- if(cp$Psi == "none") crossprod(rMnorm(q)) else diag(q)
## nu <- runif(1, q, 2*q)
## prec <- cp$prec

## n <- 1
## set.seed(test.seed)
## XVR2 <- rmniwR(Lambda = Lambda,
##                Sigma = Sigma,
##                Psi = Psi,
##                nu = nu,
##                prec = prec)
## set.seed(test.seed)
## XVcpp2 <- rMNIW(n = 1,
##                 Lambda = Lambda,
##                 Sigma = Sigma,
##                 Psi = Psi,
##                 nu = nu,
##                 prec = prec)
## tmp <- expect_equal(XVR2, XVcpp2, tolerance = tol)
## tmp
