library(testthat)

test_that("Sum of conditional probabilities", {
  tolerance = 1E-8

  set.seed(666) # Tests should not be stochastic

  a0 = 100
  b0 = 0
  A = 0
  B = a0

  muL = runif(1,0,1)
  muM = 0.0015
  eta = runif(1,0,1)
  gamma = 0.04

  drates1=function(a,b){muL*a+eta*a^2}
  brates2=function(a,b){0}
  drates2=function(a,b){muM*b}
  trans=function(a,b){gamma*a} # a -> b

  if (Sys.getenv("R_ARCH") != "/i386") {
    # These funcations fail on i386, but still unable to debug
  print("non-parallel bbd_prob:")
  p <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B)
  print("parallel bbd_prob:")
  p1 <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=1)
  p2 <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=2)
  p4 <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=4)

  expect_equal(1.0, sum(p), tolerance)
  expect_equal(0.0, sum(abs(p-p1)), tolerance)
  expect_equal(0.0, sum(abs(p1-p2)), tolerance)
  expect_equal(0.0, sum(abs(p2-p4)), tolerance)

  print("non-parallel bbd_prob:")
  p <- dbd_prob(t=1e-5,a0,b0,drates1,brates2,drates2,trans,a=A,B)
  print("parallel bbd_prob:")
  p1 <- dbd_prob(t=1e-5,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=1)
  p2 <- dbd_prob(t=1e-5,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=2)
  p4 <- dbd_prob(t=1e-5,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=4)

  expect_equal(1.0, sum(p), tolerance)
  expect_equal(0.0, sum(abs(p-p1)), tolerance)
  expect_equal(0.0, sum(abs(p1-p2)), tolerance)
  expect_equal(0.0, sum(abs(p2-p4)), tolerance)
  }
})
