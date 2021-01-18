context("Check native R Functions")
library(BoltzMM)

test_that("Check that fvbmHess calculates Hesssian Correctly.",{
  model<-list()
  model$pll<-NA
  model$bvec<-c(0,0.5,0.25)
  model$Mmat<- matrix(0.1,3,3) - diag(0.1,3,3)
  model$itt<-NA
  set.seed(1)
  data <- matrix(sample(c(-1,1),300, replace =TRUE), ncol=3)
  HessResult<-fvbmHess(data, model)

  tmp1 <- "./fvbmHess"

  expect_equal(dim(HessResult), c(6,6))
  expect_is(HessResult, "matrix")
  #expect_known_output(HessResult, tmp1, print = TRUE, update=FALSE)

})

test_that("Check that fvbmstderr calculates stderr Correctly.",{
  set.seed(1)
  data <- matrix(sample(c(-1,1),300, replace =TRUE), ncol=3)
  bvec <- c(0,0.5,0.25)
  Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
  covarmat <- diag(length(bvec)[1]*2)
  stderrResult<-fvbmstderr(data, covarmat)
  tmp2 <- "./fvbmstderr"


  expect_is(stderrResult, "list")
  expect_is(stderrResult[[1]], "numeric")
  expect_is(stderrResult[[2]], "matrix")
  expect_equal(length(stderrResult), 2)
  expect_equal(length(stderrResult[[1]]), 3)
  expect_equal(dim(stderrResult[[2]]), c(3,3))
  #expect_known_output(stderrResult, tmp2, print = TRUE, update=FALSE)

})

test_that("Check that marginpfvbm calculates marginal probailities Correctly.",{
  bvec <- c(0,0.5,0.25)
  Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
  marginResult<-marginpfvbm(bvec, Mmat)
  tmp2 <- "./marginpfvbm"

  expect_equal(length(marginResult), 3)
  # The first run always succeeds, but warns
  #expect_known_output(marginResult, tmp2, print = TRUE, update=FALSE)

})

test_that("Check that fvbmtests calculates scores and p-values correctly.",{
  set.seed(1)
  bvec <- c(0,0.5,0.25)
  Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
  num <- 1000
  data <- rfvbm(num,bvec,Mmat)
  model <- fitfvbm(data,bvec,Mmat)
  nullmodel <- list(bvec = c(0,0,0), Mmat = matrix(0,3,3))
  testResult<-fvbmtests(data,model,nullmodel)
  tmp2 <- "./fvbmtests"

  expect_is(testResult, "list")
  expect_is(testResult[[1]], "numeric")
  expect_is(testResult[[3]], "matrix")
  expect_is(testResult[[2]], "numeric")
  expect_is(testResult[[4]], "matrix")
  expect_equal(length(testResult), 4)
  expect_equal(length(testResult[[1]]), 3)
  expect_equal(length(testResult[[2]]), 3)
  expect_equal(dim(testResult[[3]]), c(3,3))
  expect_equal(dim(testResult[[4]]), c(3,3))
  #expect_known_output(testResult, tmp2, print = TRUE, update=FALSE)

})
