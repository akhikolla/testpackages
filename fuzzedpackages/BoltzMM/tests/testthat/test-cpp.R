context("Check Exported C++ Functions")
library(BoltzMM)
tol_C<-1.0e-3

test_that("Test fitfvbm ", {
  bvec<-c(0,0.5,0.25)
  Mmat<- matrix(0.1,3,3) - diag(0.1,3,3)
  set.seed(1)
  data <- matrix(sample(c(-1,1),30, replace =TRUE), ncol=3)
  fitfvbmResult<-fitfvbm(data, bvec, Mmat, delta_crit = 0.001, max_it = 1000L)

  tmp2 <- "./ffitfvbm"

  expect_is(fitfvbmResult, "list")
  expect_is(fitfvbmResult[[1]], "numeric")
  expect_is(fitfvbmResult[[2]], "numeric")
  expect_is(fitfvbmResult[[3]], "matrix")
  expect_is(fitfvbmResult[[4]], "integer")
  expect_equal(length(fitfvbmResult), 4)
  expect_equal(length(fitfvbmResult[[1]]), 1)
  expect_equal(length(fitfvbmResult[[2]]), 3)
  expect_equal(dim(fitfvbmResult[[3]]), c(3,3))
  expect_equal(length(fitfvbmResult[[4]]), 1)

  # The first run always succeeds, but warns
  ##expect_known_output(fitfvbmResult, tmp2, tolerance=tol_C, print = TRUE, update=FALSE)

})


test_that("Test pfvbm", {
  xval1 <- c(-1,1,-1)
  xval2 <- c(1,1,1)
  xval3 <- c(-1,-1,-1)
  xval4 <- c(0,-1,-1)
  xval5 <- c(NA,-1,-1)
  xval6 <- c(NA,0,NA)
  xval7 <- c(-1,1)

  bvec <- c(0,0.5,0.25)
  Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)

  expect_equal(pfvbm(xval1,bvec,Mmat), 0.1213876, tolerance=tol_C )
  expect_equal(pfvbm(xval2,bvec,Mmat), 0.2985652, tolerance=tol_C )
  expect_equal(pfvbm(xval3,bvec,Mmat), 0.0666189, tolerance=tol_C )
  expect_equal(pfvbm(xval4,bvec,Mmat), 0.05454294, tolerance=tol_C )
  expect_true(is.na(pfvbm(xval5,bvec,Mmat)))
  expect_true(is.na(pfvbm(xval6,bvec,Mmat)))
  expect_equal(pfvbm(xval7,bvec,Mmat), 0)

})


test_that("Test fvbmpartiald", {

  model<-list()
  model$pll<-NA
  model$bvec<-c(0,0.5,0.25)
  model$Mmat<- matrix(0.1,3,3) - diag(0.1,3,3)
  model$itt<-NA
  set.seed(1)
  data <- matrix(sample(c(-1,1),300, replace =TRUE), ncol=3)

  fvbmpartialdResult<-fvbmpartiald(data, model)
  tmp2 <- "./fvbmpartiald"

  expect_is(fvbmpartialdResult, "list")
  expect_is(fvbmpartialdResult[[1]], "numeric")
  expect_is(fvbmpartialdResult[[2]], "matrix")
  expect_equal(length(fvbmpartialdResult), 2)
  expect_equal(length(fvbmpartialdResult[[1]]), 3)
  expect_equal(dim(fvbmpartialdResult[[2]]), c(3,3))
  ##expect_known_output(fvbmpartialdResult, tmp2,tolerance=tol_C, print = TRUE, update=FALSE)

})

test_that(" Test allpfvbm", {
  bvec <- c(0,0.5,0.25)
  Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)


  res1<-matrix(c(0.0666189, 0.04465599, 0.1213876, 0.1213876,
          0.07362527, 0.07362527, 0.2001342,
          0.2985652), nrow=1)

  bvec2 <- c(0,0.5)

  expect_length(allpfvbm(bvec,Mmat), 8)
  expect_equal(allpfvbm(bvec,Mmat), res1,tolerance=tol_C)
  expect_equal(allpfvbm(bvec2,Mmat),matrix(c(0,0,0,0),1,4))

})


test_that("Test fvbmcov", {

  model<-list()
  model$pll<-NA
  model$bvec<-c(0,0.5,0.25)
  model$Mmat<- matrix(0.1,3,3) - diag(0.1,3,3)
  model$itt<-NA
  set.seed(1)
  data <- matrix(sample(c(-1,1),300, replace =TRUE), ncol=3)
  fvbmcovResult<-fvbmcov(data, model, fvbmHess)
  tmp2 <- "./fvbmcov"

  # The first run alwafvbmcovs succeeds, but warns
  expect_is(fvbmcovResult, "matrix")
  expect_equal(dim(fvbmcovResult), c(6,6))
  #expect_known_output(fvbmcovResult, tmp2, tolerance=tol_C, print = TRUE)

})

#test issues due to RNG from R in C++ and back to R doenst always work
# test_that("Test rfvbm", {
#   num <- 1000
#   bvec <- c(0,0.5,0.25)
#   Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
#   set.seed(1)
#   rfvbmResult<-rfvbm(num, bvec,Mmat)
#   tmp2 <- "./rfvbm"
#
#   # The first run alwarfvbms succeeds, but warns
#   #expect_known_output(rfvbmResult, tmp2, print = TRUE)
#
# })
