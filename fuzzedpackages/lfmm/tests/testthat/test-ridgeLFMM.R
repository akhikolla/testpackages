library(testthat)
context("RidgeLFMM")

#### test qui marche (Basile 09/08/19)

test_that("RidgeLFMM_main", {

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  lambda <- 1e-5
  P.list <- compute_P(X = dat$X, lambda = lambda)

  m <- ridgeLFMM(K = K,
                 lambda = lambda)

  res <- ridgeLFMM_main(m, dat, P.list)

  svd.res <- svd(P.list$sqrt.P %*% dat$Y)
  expect_lte(mean(abs(svd.res$v[,1:K] - res$V)), 0.06)
  ## RMK: error very high, it is because for this data first K singular values
  ## are quite the same.
  ## svd.res$d

})


